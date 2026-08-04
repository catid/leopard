#!/usr/bin/env python3
"""Reproduce Leopard2 tiny direct-encoder crossover measurements.

Screening mode consumes already-built ``bench_leopard2_direct_encode``
binaries.  The legacy generic ``pinned`` command is deliberately fail-closed:
hashing an externally supplied executable and its neighboring build artifacts
does not prove a clean source-to-object-to-archive-to-link closure.
``historical-avx2`` instead creates or strictly resumes a retained fresh
Release/explicit-AVX2 build while holding the canonical authoritative lock,
freezes that executable, and only then begins isolated ABBA measurement.  Every
cell is a resumable JSON job with deterministic input data and retained logs.

Typical use::

    tools/leopard2_direct_encode_crossover.py screen --build-root build
    tools/leopard2_direct_encode_crossover.py historical-avx2 \
        --cpu 16 --sibling 80
    tools/leopard2_direct_encode_crossover.py analyze \
        --result-dir results/leopard2/direct-encode-crossover/historical-avx2

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
import ctypes
import errno
import fcntl
import hashlib
import importlib.util
import inspect
import json
import math
import os
import platform
import re
import resource
import select
import shutil
import signal
import stat
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import types
import unicodedata
from pathlib import Path


SCHEMA = "leopard2-direct-encode-crossover/v7"
JOB_SCHEMA = "leopard2-direct-encode-crossover-job/v7"
ANALYSIS_SCHEMA = "leopard2-direct-encode-crossover-analysis/v4"
BENCHMARK_SCHEMA = "leopard2-direct-encode-benchmark-v2"
KNOWN_BACKENDS = ("scalar", "ssse3", "avx2", "avx512")
FROZEN_EXECUTABLE_SCHEMA = "leopard2-frozen-executable/v3"
AUTHORITATIVE_COMMANDS = ("historical-avx2",)
RUN_COMMANDS = ("screen",) + AUTHORITATIVE_COMMANDS
UNPROVED_PINNED_COMMAND = "pinned"
AUTHORITATIVE_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
DIRECT_COMMAND_SUPERVISOR_MODE = "--internal-direct-command-supervisor"
DIRECT_COMMAND_OWNER_MODE = "--internal-direct-command-owner"
MAX_SUPERVISOR_CONTROL_BYTES = 64 * 1024
MAX_RAW_JSON_BYTES = 64 * 1024 * 1024
MAX_RETAINED_LOG_BYTES = 128 * 1024 * 1024
SUPERVISOR_REAP_GRACE_SECONDS = 12.0
RAW_OUTPUT_DESCRIPTOR = 198
EXECUTABLE_DESCRIPTOR = 199
GIT_EXECUTABLE_DESCRIPTOR = 197
TASKSET_EXECUTABLE_DESCRIPTOR = 196
CONTROLLED_CMAKE_DESCRIPTOR = 195
CONTROLLED_NINJA_DESCRIPTOR = 194
CANONICAL_GIT = Path("/usr/bin/git")
CONTROLLED_BUILD_SCHEMA = "leopard2-direct-controlled-build/v7"
BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2 = (
    "leopard2-benchmark-build-configuration-attestation/v2"
)
BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V3 = (
    "leopard2-benchmark-build-configuration-attestation/v3"
)
BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V4 = (
    "leopard2-benchmark-build-configuration-attestation/v4"
)
BUILD_CONFIGURATION_ATTESTATION_SCHEMA = (
    "leopard2-benchmark-build-configuration-attestation/v5"
)
BUILD_CONFIGURATION_FILE_SCHEMA_V2 = (
    "leopard2-benchmark-build-configuration/v2"
)
BUILD_CONFIGURATION_FILE_SCHEMA_V3 = (
    "leopard2-benchmark-build-configuration/v3"
)
BUILD_CONFIGURATION_FILE_SCHEMA_V6 = (
    "leopard2-benchmark-build-configuration/v6"
)
BUILD_CONFIGURATION_FILE_SCHEMA = (
    "leopard2-benchmark-build-configuration/v7"
)
BUILD_CONFIGURATION_RELATIVE_PATH = (
    "generated/leopard2-benchmark-attestation/"
    "leopard2_benchmark_build_configuration.txt"
)
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
BUILD_CONFIGURATION_VARIABLES_V3 = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:-1],
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
BUILD_CONFIGURATION_VARIABLES_V6 = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:16],
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
    *BUILD_CONFIGURATION_VARIABLES_V2[16:-1],
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING",
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED",
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
BUILD_CONFIGURATION_VARIABLES = BUILD_CONFIGURATION_VARIABLES_V6
BUILD_CONFIGURATION_EXPERIMENT_SELECTORS = (
    *BUILD_CONFIGURATION_VARIABLES[16:],
)
BUILD_CONFIGURATION_EXPERIMENT_SELECTORS_V6 = \
    BUILD_CONFIGURATION_EXPERIMENT_SELECTORS
BUILD_CONFIGURATION_CANONICAL_SELECTORS_V2 = {
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
}
BUILD_CONFIGURATION_CANONICAL_SELECTORS_V3 = {
    **BUILD_CONFIGURATION_CANONICAL_SELECTORS_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
}
BUILD_CONFIGURATION_CANONICAL_SELECTORS_V6 = {
    **BUILD_CONFIGURATION_CANONICAL_SELECTORS_V2,
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
}
BUILD_CONFIGURATION_CANONICAL_SELECTORS = \
    BUILD_CONFIGURATION_CANONICAL_SELECTORS_V6
CMAKE_CACHE_ENTRY_TYPES = frozenset((
    "BOOL", "FILEPATH", "INTERNAL", "PATH", "STATIC", "STRING",
    "UNINITIALIZED",
))
CMAKE_CACHE_REQUIRED_ENTRY_TYPES = {
    # These bindings have canonical CMake-declared types.  Reject a
    # hand-edited cache that preserves text while changing CMake semantics.
    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
        frozenset(("INTERNAL",)),
    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
        frozenset(("INTERNAL",)),
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK":
        frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED":
        frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": frozenset(("STRING",)),
}
BENCHMARK_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
GIT_ENVIRONMENT = dict(BENCHMARK_ENVIRONMENT)
GIT_ENVIRONMENT.update({
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_SYSTEM": "/dev/null",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
})
SOURCE_IDENTITY_LOCK = threading.Lock()
RUNTIME_LAUNCHER_NAMES = (
    "python", "script", "build_provenance", "git_capture",
    "link_common", "containment",
)
RUNTIME_PYTHON_FLAGS = ("-I", "-S", "-B")


class CrossoverError(Exception):
    """An actionable configuration, execution, or result error."""


class AuthoritativeGuardError(CrossoverError):
    """A held canonical or physical-pair guard was lost or replaced."""


def attested_text(value, description):
    if not isinstance(value, str):
        raise CrossoverError("{} is not text".format(description))
    try:
        value.encode("utf-8")
    except UnicodeError as error:
        raise CrossoverError(
            "{} is not valid UTF-8 text: {}".format(description, error)
        )
    forbidden = sorted(set(value) & {"\0", "\n", "\r"})
    if forbidden:
        raise CrossoverError(
            "{} contains a forbidden record delimiter".format(description)
        )
    forbidden_categories = sorted(set(
        unicodedata.category(character)
        for character in value
        if unicodedata.category(character) in (
            "Cc", "Cf", "Zl", "Zp",
        )
    ))
    if forbidden_categories:
        raise CrossoverError(
            "{} contains forbidden Unicode categories {}".format(
                description, ",".join(forbidden_categories)
            )
        )
    return value


def checked_path(value, description):
    try:
        path_value = os.fspath(value)
    except (OSError, TypeError, ValueError) as error:
        raise CrossoverError("{} is not a valid path: {}".format(
            description, error
        ))
    if isinstance(path_value, bytes):
        try:
            path_value = os.fsdecode(path_value)
        except (OSError, UnicodeError, ValueError) as error:
            raise CrossoverError("{} is not a valid path: {}".format(
                description, error
            ))
    if not isinstance(path_value, str) or "\0" in path_value:
        raise CrossoverError(
            "{} is not a valid NUL-free path".format(description)
        )
    try:
        return Path(path_value)
    except (OSError, TypeError, ValueError) as error:
        raise CrossoverError("{} is not a valid path: {}".format(
            description, error
        ))


def checked_absolute_path(value, description):
    path = checked_path(value, description)
    try:
        absolute = path.is_absolute()
    except (OSError, ValueError) as error:
        raise CrossoverError("{} cannot be checked: {}".format(
            description, error
        ))
    if not absolute:
        raise CrossoverError("{} is not absolute".format(description))
    return path


def checked_resolve(value, description, strict=False):
    path = checked_path(value, description)
    try:
        return path.resolve(strict=strict)
    except (OSError, RuntimeError, ValueError) as error:
        raise CrossoverError("{} cannot be resolved: {}".format(
            description, error
        ))


def load_isolation_support(source):
    path = (
        Path(source) / "experiments" / "leopard2" /
        "main_compare" / "run_abba.py"
    ).resolve()
    if not path.is_file():
        raise CrossoverError(
            "canonical isolation support is missing: {}".format(path)
        )
    name = "leopard2_direct_crossover_isolation_support"
    specification = importlib.util.spec_from_file_location(name, str(path))
    if specification is None or specification.loader is None:
        raise CrossoverError(
            "cannot load canonical isolation support: {}".format(path)
        )
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    try:
        specification.loader.exec_module(module)
    except Exception as error:
        raise CrossoverError(
            "cannot initialize canonical isolation support: {}".format(error)
        )
    for function in (
            "cpu_stat_snapshot", "isolation_record", "PairLease",
            "validate_isolation", "validate_pair_lease_identity",
            "validate_topology"):
        if not callable(getattr(module, function, None)):
            raise CrossoverError(
                "canonical isolation support omits {}".format(function)
            )
    if getattr(module, "MAX_SIBLING_NONIDLE_JIFFIES", None) != 0:
        raise CrossoverError(
            "canonical isolation support no longer requires an idle sibling"
        )
    return module


class CanonicalAuthoritativeLock:
    def __init__(self, path=AUTHORITATIVE_LOCK):
        self.path = Path(path)
        self.descriptor = None
        self.identity = None

    def validate_current(self):
        if self.descriptor is None or self.identity is None:
            raise AuthoritativeGuardError(
                "canonical authoritative lock is not held"
            )
        try:
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(str(self.path))
        except OSError as error:
            raise AuthoritativeGuardError(
                "canonical authoritative lock revalidation failed: {}".format(
                    error
                )
            )
        current = {
            "device": descriptor.st_dev,
            "inode": descriptor.st_ino,
            "lock": "exclusive",
            "mode": stat.S_IMODE(descriptor.st_mode),
            "path": str(self.path),
            "uid": descriptor.st_uid,
        }
        if (not stat.S_ISREG(descriptor.st_mode) or
                descriptor.st_nlink != 1 or
                (descriptor.st_dev, descriptor.st_ino) !=
                (path.st_dev, path.st_ino) or current != self.identity):
            raise AuthoritativeGuardError(
                "canonical authoritative lock path was replaced or changed"
            )
        return self.identity

    def __enter__(self):
        if self.path != AUTHORITATIVE_LOCK:
            raise CrossoverError(
                "authoritative runs require canonical lock {}".format(
                    AUTHORITATIVE_LOCK
                )
            )
        flags = os.O_RDWR | os.O_CREAT
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            self.descriptor = os.open(str(self.path), flags, 0o600)
            metadata = os.fstat(self.descriptor)
            path_metadata = os.lstat(str(self.path))
            if (not stat.S_ISREG(metadata.st_mode) or
                    metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                    stat.S_IMODE(metadata.st_mode) != 0o600 or
                    (metadata.st_dev, metadata.st_ino) !=
                    (path_metadata.st_dev, path_metadata.st_ino)):
                raise CrossoverError(
                    "canonical authoritative lock has unsafe metadata"
                )
            fcntl.flock(self.descriptor, fcntl.LOCK_EX)
            locked = os.fstat(self.descriptor)
            self.identity = {
                "device": locked.st_dev,
                "inode": locked.st_ino,
                "lock": "exclusive",
                "mode": stat.S_IMODE(locked.st_mode),
                "path": str(self.path),
                "uid": locked.st_uid,
            }
            self.validate_current()
            return self
        except BaseException:
            if self.descriptor is not None:
                os.close(self.descriptor)
                self.descriptor = None
                self.identity = None
            raise

    def __exit__(self, exc_type, exc, traceback):
        validation_error = None
        try:
            self.validate_current()
        except Exception as error:
            validation_error = error
        descriptor = self.descriptor
        self.descriptor = None
        self.identity = None
        # flock(LOCK_UN) would release the lock for the shared open-file
        # description even while a timed child still owns an inherited
        # descriptor.  Close-only lifetime lets the kernel release the lock
        # only after the coordinator and every descendant have closed it.
        if descriptor is not None:
            os.close(descriptor)
        if validation_error is not None:
            raise validation_error


def canonical_authoritative_lock(path=AUTHORITATIVE_LOCK):
    return CanonicalAuthoritativeLock(path)


class AuthoritativePairGuard:
    """Add fail-closed exit validation to the shared pair-lease helper."""

    def __init__(self, guard):
        self.guard = guard
        self.identity = None

    def validate_current(self):
        try:
            self.guard.validate_current()
        except Exception as error:
            raise AuthoritativeGuardError(
                "physical-pair guard revalidation failed: {}".format(error)
            )
        return self.identity

    def __enter__(self):
        self.identity = self.guard.__enter__()
        try:
            self.validate_current()
        except Exception:
            try:
                self.guard.__exit__(None, None, None)
            finally:
                self.identity = None
            raise
        return self.identity

    def __exit__(self, exc_type, exc, traceback):
        validation_error = None
        try:
            self.validate_current()
        except Exception as error:
            validation_error = error
        try:
            result = self.guard.__exit__(exc_type, exc, traceback)
        finally:
            self.identity = None
        if validation_error is not None:
            raise validation_error
        return result


def validate_authoritative_guards(lock_guard, pair_guard):
    try:
        lock_guard.validate_current()
        pair_guard.validate_current()
    except AuthoritativeGuardError:
        raise
    except Exception as error:
        raise AuthoritativeGuardError(
            "authoritative guard revalidation failed: {}".format(error)
        )


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False
    ).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical_bytes(value))


def metadata_signature(metadata):
    """Return all mutable regular-file metadata relevant to evidence."""
    return (
        metadata.st_dev, metadata.st_ino, metadata.st_mode,
        metadata.st_uid, metadata.st_nlink, metadata.st_size,
        metadata.st_mtime_ns, metadata.st_ctime_ns,
    )


def read_descriptor_bytes(descriptor, maximum_bytes, description):
    """Double-read one retained inode and bind bytes to stable metadata."""
    if (type(descriptor) is not int or descriptor < 0 or
            type(maximum_bytes) is not int or maximum_bytes < 0):
        raise CrossoverError(
            "{} descriptor-read contract is invalid".format(description)
        )
    try:
        before = os.fstat(descriptor)
        values = []
        for unused_pass in range(2):
            os.lseek(descriptor, 0, os.SEEK_SET)
            chunks = []
            total = 0
            while total <= maximum_bytes:
                block = os.read(
                    descriptor,
                    min(1024 * 1024, maximum_bytes + 1 - total),
                )
                if not block:
                    break
                chunks.append(block)
                total += len(block)
            if total > maximum_bytes:
                raise CrossoverError(
                    "{} exceeds {} bytes".format(
                        description, maximum_bytes
                    )
                )
            values.append(b"".join(chunks))
        after = os.fstat(descriptor)
    except OSError as error:
        raise CrossoverError(
            "cannot read {}: {}".format(description, error)
        )
    if (not stat.S_ISREG(before.st_mode) or
            values[0] != values[1] or
            metadata_signature(before) != metadata_signature(after) or
            len(values[0]) != before.st_size):
        raise CrossoverError(
            "{} changed while its exact inode was read".format(description)
        )
    return values[0]


def linux_memfd_create(name, flags):
    """Call Linux memfd_create when the Python wrapper is unavailable."""
    wrapper = getattr(os, "memfd_create", None)
    if callable(wrapper):
        return wrapper(name, flags)
    if sys.platform != "linux":
        raise CrossoverError("Linux memfd sealing is unavailable")
    try:
        function = ctypes.CDLL(None, use_errno=True).memfd_create
    except (AttributeError, OSError) as error:
        raise CrossoverError(
            "Linux memfd sealing is unavailable: {}".format(error)
        )
    function.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    descriptor = function(name.encode("utf-8"), ctypes.c_uint(flags))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    raise CrossoverError(
        "Linux memfd_create failed: {}".format(
            os.strerror(number or errno.EPERM)
        )
    )


def create_sealed_executable_snapshot(value, description):
    """Create a byte-immutable Linux memfd suitable for exact execution."""
    if not isinstance(value, bytes) or not value:
        raise CrossoverError(
            "{} snapshot payload is invalid".format(description)
        )
    if sys.platform != "linux":
        raise CrossoverError(
            "{} requires Linux memfd sealing".format(description)
        )
    descriptor = None
    transferred = False
    try:
        descriptor = linux_memfd_create(
            "leopard2-executable-snapshot",
            getattr(os, "MFD_CLOEXEC", 0x0001) |
            getattr(os, "MFD_ALLOW_SEALING", 0x0002),
        )
        write_descriptor_all(descriptor, value)
        os.fchmod(descriptor, 0o555)
        os.fsync(descriptor)
        add_seals = getattr(fcntl, "F_ADD_SEALS", 1033)
        get_seals = getattr(fcntl, "F_GET_SEALS", 1034)
        required_seals = (
            getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
            getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
            getattr(fcntl, "F_SEAL_GROW", 0x0004) |
            getattr(fcntl, "F_SEAL_WRITE", 0x0008)
        )
        fcntl.fcntl(descriptor, add_seals, required_seals)
        if fcntl.fcntl(descriptor, get_seals) != required_seals:
            raise CrossoverError(
                "{} did not acquire the exact immutable seal set".format(
                    description
                )
            )
        metadata = os.fstat(descriptor)
        retained = read_descriptor_bytes(
            descriptor, len(value), description + " sealed snapshot"
        )
        after = os.fstat(descriptor)
        if (retained != value or not stat.S_ISREG(metadata.st_mode) or
                stat.S_IMODE(metadata.st_mode) != 0o555 or
                metadata.st_size != len(value) or
                metadata_signature(metadata) != metadata_signature(after)):
            raise CrossoverError(
                "{} sealed executable snapshot changed".format(description)
            )
        result = {
            "descriptor": descriptor,
            "sha256": digest_bytes(value),
            "size": len(value),
        }
        transferred = True
        return result
    except OSError as error:
        transferred = False
        raise CrossoverError(
            "cannot seal {}: {}".format(description, error)
        )
    except BaseException:
        transferred = False
        raise
    finally:
        # A successful return transfers the retained descriptor.  If a
        # BaseException is delivered while evaluating that return, ownership
        # remains here and the descriptor must still be closed.
        if descriptor is not None and not transferred:
            os.close(descriptor)


def validate_sealed_executable_snapshot(snapshot, expected_sha256,
                                        description):
    if (not isinstance(snapshot, dict) or
            set(snapshot) not in (
                {"descriptor", "sha256", "size"},
                {"descriptor", "path", "sha256", "size"},
            ) or
            type(snapshot.get("descriptor")) is not int or
            type(snapshot.get("size")) is not int or
            snapshot["size"] <= 0 or
            snapshot.get("sha256") != expected_sha256):
        raise CrossoverError(
            "{} sealed snapshot contract is invalid".format(description)
        )
    descriptor = snapshot["descriptor"]
    try:
        seals = fcntl.fcntl(
            descriptor, getattr(fcntl, "F_GET_SEALS", 1034)
        )
        required_seals = (
            getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
            getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
            getattr(fcntl, "F_SEAL_GROW", 0x0004) |
            getattr(fcntl, "F_SEAL_WRITE", 0x0008)
        )
        value = read_descriptor_bytes(
            descriptor, snapshot["size"], description
        )
        metadata = os.fstat(descriptor)
    except OSError as error:
        raise CrossoverError(
            "cannot validate {}: {}".format(description, error)
        )
    if (seals != required_seals or
            stat.S_IMODE(metadata.st_mode) != 0o555 or
            metadata.st_size != snapshot["size"] or
            digest_bytes(value) != expected_sha256):
        raise CrossoverError(
            "{} sealed snapshot identity is stale".format(description)
        )
    return value


def close_sealed_snapshot(snapshot):
    if snapshot is not None and snapshot.get("descriptor") is not None:
        try:
            os.close(snapshot["descriptor"])
        finally:
            snapshot["descriptor"] = None


def open_exact_executable_snapshot(
        path, description, require_executable=True, allow_writable=False):
    """Seal bytes read from one stable no-follow external file inode."""
    path = checked_absolute_path(path, description)
    descriptor = None
    snapshot = None
    snapshot_transferred = False
    try:
        descriptor = os.open(
            str(path),
            os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
        )
        before = os.fstat(descriptor)
        current = os.lstat(str(path))
        value = read_descriptor_bytes(
            descriptor, MAX_RAW_JSON_BYTES, description
        )
        after = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                (not allow_writable and
                 stat.S_IMODE(before.st_mode) & 0o022) or
                (require_executable and not before.st_mode & 0o111) or
                metadata_signature(before) != metadata_signature(after) or
                (before.st_dev, before.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "{} is not one stable permitted regular-file inode".format(
                    description
                )
            )
        snapshot = create_sealed_executable_snapshot(value, description)
        snapshot["path"] = str(path)
        snapshot_transferred = True
        return snapshot
    except OSError as error:
        snapshot_transferred = False
        raise CrossoverError(
            "cannot retain {}: {}".format(description, error)
        )
    except BaseException:
        snapshot_transferred = False
        raise
    finally:
        try:
            if descriptor is not None:
                os.close(descriptor)
        except BaseException:
            snapshot_transferred = False
            if snapshot is not None:
                close_sealed_snapshot(snapshot)
            raise
        if snapshot is not None and not snapshot_transferred:
            close_sealed_snapshot(snapshot)


def open_current_interpreter_snapshot(python_path=None):
    """Seal the running executable and bind a byte-identical prefix argv[0]."""
    descriptor = None
    snapshot = None
    candidate_snapshot = None
    transferred = False
    try:
        descriptor = os.open(
            "/proc/self/exe",
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NONBLOCK", 0),
        )
        before = os.fstat(descriptor)
        value = read_descriptor_bytes(
            descriptor, MAX_RAW_JSON_BYTES,
            "current Python interpreter executable",
        )
        after = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                not before.st_mode & 0o111 or
                metadata_signature(before) != metadata_signature(after)):
            raise CrossoverError(
                "current Python interpreter executable is not stable"
            )
        snapshot = create_sealed_executable_snapshot(
            value, "current Python interpreter executable"
        )

        if python_path is None:
            prefix = Path(os.path.abspath(sys.base_prefix))
            name = "python{}.{}{}".format(
                sys.version_info.major, sys.version_info.minor, sys.abiflags
            )
            candidate = prefix / "bin" / name
        else:
            candidate = Path(python_path)
        candidate = candidate.resolve(strict=True)
        if not candidate.is_absolute():
            raise CrossoverError(
                "Python prefix-discovery executable is not absolute"
            )
        candidate_snapshot = open_exact_executable_snapshot(
            candidate, "Python prefix-discovery executable",
            require_executable=True, allow_writable=True,
        )
        if (candidate_snapshot["size"] != snapshot["size"] or
                candidate_snapshot["sha256"] != snapshot["sha256"]):
            raise CrossoverError(
                "Python prefix-discovery executable differs from the "
                "running interpreter inode"
            )
        snapshot["path"] = str(candidate)
        transferred = True
        return snapshot
    except OSError as error:
        transferred = False
        raise CrossoverError(
            "cannot retain current Python interpreter: {}".format(error)
        ) from error
    except BaseException:
        transferred = False
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)
        close_sealed_snapshot(candidate_snapshot)
        if snapshot is not None and not transferred:
            close_sealed_snapshot(snapshot)


def sealed_snapshot_identity(snapshot):
    if not isinstance(snapshot, dict):
        raise CrossoverError("sealed snapshot is not an object")
    return {
        "path": snapshot["path"],
        "sha256": snapshot["sha256"],
        "size": snapshot["size"],
    }


def validate_sealed_snapshot_identity(value, description):
    if not isinstance(value, dict) or set(value) != {
            "path", "sha256", "size"}:
        raise CrossoverError(
            "{} sealed identity is invalid".format(description)
        )
    path = attested_text(value.get("path"), description + " path")
    try:
        path_is_absolute = Path(path).is_absolute()
    except (OSError, TypeError, ValueError):
        path_is_absolute = False
    if (not path_is_absolute or
            not isinstance(value.get("sha256"), str) or
            not re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) or
            type(value.get("size")) is not int or value["size"] <= 0):
        raise CrossoverError(
            "{} sealed identity is invalid".format(description)
        )
    return value


def open_runtime_launcher_snapshots(
        python_path=None, script_path=None):
    script_path = (
        Path(__file__).resolve()
        if script_path is None else Path(script_path).resolve()
    )
    repository = Path(__file__).resolve().parents[1]
    paths = {
        # The current interpreter may come from a user-managed installation
        # with group-writable mode bits.  It is never executed by pathname:
        # exact bytes are read once and then executed only from a sealed memfd.
        "python": (python_path, True, True),
        "script": (script_path, False, True),
        "build_provenance": (
            repository / "tools" / "leopard2_build_provenance.py",
            False, True,
        ),
        "git_capture": (
            repository / "experiments" / "leopard2" /
            "main_compare" / "git_capture.py",
            False, True,
        ),
        "link_common": (
            repository / "experiments" / "leopard2" /
            "decoder_dispatch" / "balanced_evidence_common.py",
            False, True,
        ),
        "containment": (
            repository / "experiments" / "leopard2" /
            "main_compare" / "run_abba.py",
            False, True,
        ),
    }
    snapshots = {}
    try:
        snapshots["python"] = open_current_interpreter_snapshot(python_path)
        for name in RUNTIME_LAUNCHER_NAMES:
            if name == "python":
                continue
            path, require_executable, allow_writable = paths[name]
            snapshots[name] = open_exact_executable_snapshot(
                path, "command-owner {} launcher input".format(name),
                require_executable, allow_writable,
            )
        snapshots["identity"] = {
            name: sealed_snapshot_identity(snapshots[name])
            for name in RUNTIME_LAUNCHER_NAMES
        }
        return snapshots
    except BaseException as primary:
        cleanup_errors = []
        for name in reversed(RUNTIME_LAUNCHER_NAMES):
            if name in snapshots:
                try:
                    close_sealed_snapshot(snapshots[name])
                except BaseException as error:
                    cleanup_errors.append((name, error))
        if cleanup_errors:
            raise CrossoverError(
                "runtime launcher setup failed: {}; cleanup also failed: "
                "{}".format(
                    primary,
                    "; ".join(
                        "{} {}: {}".format(
                            name, type(error).__name__, error
                        ) for name, error in cleanup_errors
                    ),
                )
            ) from primary
        raise


def validate_runtime_launcher_snapshots(launcher):
    if (not isinstance(launcher, dict) or
            set(launcher) != {"identity"} | set(RUNTIME_LAUNCHER_NAMES) or
            not isinstance(launcher.get("identity"), dict) or
            set(launcher["identity"]) != set(RUNTIME_LAUNCHER_NAMES)):
        raise CrossoverError("runtime launcher contract is invalid")
    for name in RUNTIME_LAUNCHER_NAMES:
        identity = validate_sealed_snapshot_identity(
            launcher["identity"].get(name),
            "runtime launcher " + name,
        )
        snapshot = launcher.get(name)
        if (not isinstance(snapshot, dict) or
                sealed_snapshot_identity(snapshot) != identity):
            raise CrossoverError(
                "runtime launcher {} snapshot is stale".format(name)
            )
        validate_sealed_executable_snapshot(
            snapshot, identity["sha256"],
            "runtime launcher " + name,
        )
    return launcher["identity"]


def runtime_launcher_contract(launcher):
    identity = validate_runtime_launcher_snapshots(launcher)
    result = {"identity": identity}
    result.update({
        name + "_descriptor": launcher[name]["descriptor"]
        for name in RUNTIME_LAUNCHER_NAMES
    })
    return result


def runtime_launcher_from_contract(value):
    descriptor_keys = {
        name + "_descriptor" for name in RUNTIME_LAUNCHER_NAMES
    }
    if (not isinstance(value, dict) or
            set(value) != {"identity"} | descriptor_keys or
            not isinstance(value.get("identity"), dict) or
            set(value["identity"]) != set(RUNTIME_LAUNCHER_NAMES) or
            any(type(value.get(key)) is not int or value[key] < 3
                for key in descriptor_keys) or
            len({value[key] for key in descriptor_keys}) !=
            len(descriptor_keys)):
        raise CrossoverError("runtime launcher transport is invalid")
    launcher = {"identity": value["identity"]}
    for name in RUNTIME_LAUNCHER_NAMES:
        identity = validate_sealed_snapshot_identity(
            value["identity"].get(name), "runtime launcher " + name
        )
        launcher[name] = {
            "descriptor": value[name + "_descriptor"],
            "path": identity["path"],
            "sha256": identity["sha256"],
            "size": identity["size"],
        }
    validate_runtime_launcher_snapshots(launcher)
    return launcher


def close_runtime_launcher_snapshots(launcher):
    errors = []
    if launcher is None:
        return
    for name in reversed(RUNTIME_LAUNCHER_NAMES):
        try:
            close_sealed_snapshot(launcher.get(name))
        except BaseException as error:
            errors.append(error)
    if errors:
        raise CrossoverError(
            "cannot close runtime launcher snapshots: {}".format(
                "; ".join(str(error) for error in errors)
            )
        )


def current_runtime_launcher_identity():
    launcher = None
    try:
        launcher = open_runtime_launcher_snapshots()
        return validate_runtime_launcher_snapshots(launcher)
    finally:
        if launcher is not None:
            close_runtime_launcher_snapshots(launcher)


def validate_current_executable_identity(identity, description):
    identity = validate_sealed_snapshot_identity(identity, description)
    snapshot = None
    try:
        snapshot = open_exact_executable_snapshot(
            identity["path"], description, True
        )
        if sealed_snapshot_identity(snapshot) != identity:
            raise CrossoverError(
                "{} executable identity changed".format(description)
            )
        return identity
    finally:
        if snapshot is not None:
            close_sealed_snapshot(snapshot)


def load_sealed_python_module(snapshot, module_name, description):
    if (not isinstance(module_name, str) or
            not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", module_name)):
        raise CrossoverError(
            "{} module name is invalid".format(description)
        )
    value = validate_sealed_executable_snapshot(
        snapshot, snapshot.get("sha256"), description
    )
    module = types.ModuleType(module_name)
    module.__file__ = snapshot["path"]
    module.__package__ = ""
    previous = sys.modules.get(module_name)
    sys.modules[module_name] = module
    try:
        code = compile(value, snapshot["path"], "exec")
        exec(code, module.__dict__)  # pylint: disable=exec-used
    except BaseException:
        if previous is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    return module


def evidence_contract(frozen_executable_required):
    return {
        "benchmark_schema": BENCHMARK_SCHEMA,
        "build_configuration_attestation_schema":
            BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
        "build_configuration_file_schema":
            BUILD_CONFIGURATION_FILE_SCHEMA,
        "execution_mad_path": (
            "metrics.encode_execution.mad_us_per_batch_call"
        ),
        "execution_median_path": (
            "metrics.encode_execution.median_us_per_batch_call"
        ),
        "frozen_executable_required": bool(frozen_executable_required),
        "numeric_policy": "finite JSON numbers; median > 0 and MAD >= 0",
        "parity_identity_fields": [
            "algorithm", "digest", "hashed_bytes",
            "requested_parity_indices",
        ],
        "parity_oracles": [
            "selected_transform_reference_parity_match",
            "direct_transform_parity_match on force-direct invocation",
            "same non-vacuous parity identity across direct/transform ABBA",
            "unrequested_outputs_untouched",
        ],
    }


def decode_json_bytes(value, description):
    def reject_duplicate_keys(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise CrossoverError(
                    "{} contains duplicate JSON key {!r}".format(description, key)
                )
            result[key] = item
        return result

    def reject_nonstandard_constant(constant):
        raise CrossoverError(
            "{} contains non-standard JSON number {}".format(
                description, constant
            )
        )

    def parse_finite_float(text):
        try:
            number = float(text)
        except (OverflowError, ValueError):
            number = float("nan")
        if not math.isfinite(number):
            raise CrossoverError(
                "{} contains a non-finite JSON number {}".format(
                    description, text
                )
            )
        return number

    def parse_bounded_integer(text):
        try:
            number = int(text)
        except ValueError:
            raise CrossoverError(
                "{} contains an invalid JSON integer".format(description)
            )
        if abs(number) > (1 << 64) - 1:
            raise CrossoverError(
                "{} contains a JSON integer outside the 64-bit evidence "
                "domain".format(description)
            )
        return number

    try:
        text = value.decode("utf-8")
    except (AttributeError, UnicodeError) as error:
        raise CrossoverError("cannot decode {}: {}".format(description, error))
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonstandard_constant,
            parse_float=parse_finite_float,
            parse_int=parse_bounded_integer,
        )
    except CrossoverError:
        raise
    except (ValueError, OverflowError, RecursionError) as error:
        raise CrossoverError("cannot parse {}: {}".format(description, error))


def normalized_output(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def _open_absolute_directory_components(
        path, description, create, required_mode=0o700):
    """Open an absolute directory one no-follow component at a time.

    ``O_NOFOLLOW`` on an absolute pathname protects only its final component.
    Starting from a retained root dirfd closes the resolve/open race for every
    intermediate component as well.
    """
    if required_mode not in (0o555, 0o700):
        raise CrossoverError(
            "{} directory mode contract is invalid".format(description)
        )
    if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
        raise CrossoverError(
            "{} requires O_NOFOLLOW and O_DIRECTORY".format(description)
        )
    lexical = Path(os.path.abspath(os.fspath(path)))
    if not lexical.is_absolute():
        raise CrossoverError("{} is not absolute".format(description))
    flags = (
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
        getattr(os, "O_CLOEXEC", 0)
    )
    descriptor = None
    try:
        components = list(lexical.parts[1:])
        if (len(components) >= 4 and components[:3] ==
                ["proc", "self", "fd"] and components[3].isdigit()):
            retained_descriptor = int(components[3])
            descriptor = os.dup(retained_descriptor)
            retained_metadata = os.fstat(descriptor)
            if not stat.S_ISDIR(retained_metadata.st_mode):
                raise CrossoverError(
                    "{} retained procfd is not a directory".format(
                        description
                    )
                )
            current = Path("/proc/self/fd") / components[3]
            components = components[4:]
        else:
            descriptor = os.open("/", flags)
            current = Path("/")
        for component in components:
            if component in ("", ".", "..") or "/" in component:
                raise CrossoverError(
                    "{} has an unsafe component".format(description)
                )
            try:
                child = os.open(component, flags, dir_fd=descriptor)
            except FileNotFoundError:
                if not create:
                    raise
                os.mkdir(component, 0o700, dir_fd=descriptor)
                child = os.open(component, flags, dir_fd=descriptor)
            except FileExistsError:
                child = os.open(component, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
            current = current / component
        metadata = os.fstat(descriptor)
        if (not stat.S_ISDIR(metadata.st_mode) or
                metadata.st_uid != os.getuid() or
                stat.S_IMODE(metadata.st_mode) != required_mode):
            raise CrossoverError(
                "{} is not an exact mode-{:04o} owner directory: {}".format(
                    description, required_mode, current
                )
            )
        return {
            "descriptor": descriptor,
            "identity": (metadata.st_dev, metadata.st_ino),
            "path": lexical,
            "required_mode": required_mode,
        }
    except FileNotFoundError:
        if descriptor is not None:
            os.close(descriptor)
        raise
    except OSError as error:
        if descriptor is not None:
            os.close(descriptor)
        raise CrossoverError(
            "cannot open {} without following links: {}".format(
                description, error
            )
        )
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        raise


def owned_canonical_directory(path, description):
    """Create/open one owner-controlled directory without following links."""
    return _open_absolute_directory_components(path, description, True)


def validate_owned_directory(directory, description):
    current = None
    try:
        try:
            metadata = os.fstat(directory["descriptor"])
            current = _open_absolute_directory_components(
                directory["path"], description + " current path", False,
                directory.get("required_mode", 0o700),
            )
            current_metadata = os.fstat(current["descriptor"])
        except OSError as error:
            raise CrossoverError(
                "cannot revalidate {}: {}".format(description, error)
            )
        required_mode = directory.get("required_mode", 0o700)
        if (not stat.S_ISDIR(metadata.st_mode) or
                metadata.st_uid != os.getuid() or
                stat.S_IMODE(metadata.st_mode) != required_mode or
                (metadata.st_dev, metadata.st_ino) != directory["identity"] or
                (current_metadata.st_dev, current_metadata.st_ino) !=
                directory["identity"] or
                stat.S_IMODE(current_metadata.st_mode) != required_mode):
            raise CrossoverError(
                "{} was replaced or became unsafe".format(description)
            )
    finally:
        if current is not None:
            close_owned_directory(current)


def close_owned_directory(directory):
    descriptor = directory.get("descriptor")
    if descriptor is not None:
        os.close(descriptor)
        directory["descriptor"] = None


def duplicate_owned_directory(directory, description):
    validate_owned_directory(directory, description)
    return {
        "descriptor": os.dup(directory["descriptor"]),
        "identity": directory["identity"],
        "path": directory["path"],
        "required_mode": directory.get("required_mode", 0o700),
    }


def result_relative_parts(relative, description):
    try:
        path = Path(relative)
    except (TypeError, ValueError) as error:
        raise CrossoverError("{} is invalid: {}".format(description, error))
    if path.is_absolute():
        raise CrossoverError("{} must be relative".format(description))
    parts = path.parts
    if (not parts or any(
            part in ("", ".", "..") or "/" in part or "\0" in part
            for part in parts)):
        raise CrossoverError("{} is unsafe".format(description))
    return parts


def open_result_directory(
        root, relative, description, create=False, required_mode=0o700):
    """Traverse retained result directories through O_NOFOLLOW dirfds."""
    if required_mode not in (0o555, 0o700):
        raise CrossoverError(
            "{} directory mode contract is invalid".format(description)
        )
    parts = result_relative_parts(relative, description)
    validate_owned_directory(root, "canonical result directory")
    descriptor = os.dup(root["descriptor"])
    current_path = root["path"]
    flags = (
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
        getattr(os, "O_CLOEXEC", 0)
    )
    try:
        for component_index, component in enumerate(parts):
            try:
                child = os.open(component, flags, dir_fd=descriptor)
            except FileNotFoundError:
                if not create:
                    raise
                os.mkdir(component, 0o700, dir_fd=descriptor)
                child = os.open(component, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
            current_path = current_path / component
            metadata = os.fstat(descriptor)
            component_mode = (
                required_mode
                if component_index == len(parts) - 1 else 0o700
            )
            if (not stat.S_ISDIR(metadata.st_mode) or
                    metadata.st_uid != os.getuid() or
                    stat.S_IMODE(metadata.st_mode) != component_mode):
                raise CrossoverError(
                    "{} is not an exact mode-{:04o} owner directory".format(
                        description, component_mode
                    )
                )
        metadata = os.fstat(descriptor)
        result = {
            "descriptor": descriptor,
            "identity": (metadata.st_dev, metadata.st_ino),
            "path": current_path,
            "required_mode": required_mode,
        }
        validate_owned_directory(result, description)
        validate_owned_directory(root, "canonical result directory")
        return result
    except BaseException:
        os.close(descriptor)
        raise


def open_existing_result_directory(
        root, relative, description, required_mode=0o700):
    return open_result_directory(
        root, relative, description, False, required_mode
    )


def read_result_regular(
        root, relative, maximum_bytes, description, mutation_hook=None,
        required_mode=0o600):
    """Read one retained artifact via a stable no-follow result-root walk."""
    parts = result_relative_parts(relative, description)
    if len(parts) == 1:
        directory = {
            "descriptor": os.dup(root["descriptor"]),
            "identity": root["identity"],
            "path": root["path"],
        }
    else:
        directory = open_existing_result_directory(
            root, Path(*parts[:-1]), description + " directory"
        )
    try:
        value = read_owned_regular(
            directory, parts[-1], maximum_bytes, description,
            mutation_hook=mutation_hook,
            required_mode=required_mode,
        )
        validate_owned_directory(root, "canonical result directory")
        return value
    finally:
        close_owned_directory(directory)


def list_result_regular_names(root, relative, suffix, description):
    """Enumerate a retained directory and reject every non-regular entry."""
    directory = open_existing_result_directory(root, relative, description)
    try:
        before = os.fstat(directory["descriptor"])
        before_identity = (
            before.st_dev, before.st_ino, before.st_size,
            before.st_mtime_ns, before.st_ctime_ns,
        )
        names = sorted(os.listdir(directory["descriptor"]))
        result = []
        for name in names:
            if (not isinstance(name, str) or name in ("", ".", "..") or
                    "/" in name or "\0" in name):
                raise CrossoverError(
                    "{} contains an unsafe name".format(description)
                )
            metadata = os.stat(
                name, dir_fd=directory["descriptor"],
                follow_symlinks=False,
            )
            if (not stat.S_ISREG(metadata.st_mode) or
                    metadata.st_uid != os.getuid() or
                    metadata.st_nlink != 1 or
                    stat.S_IMODE(metadata.st_mode) != 0o600):
                raise CrossoverError(
                    "{} contains unsafe entry {!r}".format(
                        description, name
                    )
                )
            if suffix is None or name.endswith(suffix):
                result.append(name)
            else:
                raise CrossoverError(
                    "{} contains unexpected entry {!r}".format(
                        description, name
                    )
                )
        after_names = sorted(os.listdir(directory["descriptor"]))
        after = os.fstat(directory["descriptor"])
        after_identity = (
            after.st_dev, after.st_ino, after.st_size,
            after.st_mtime_ns, after.st_ctime_ns,
        )
        if names != after_names or before_identity != after_identity:
            raise CrossoverError(
                "{} changed while it was being enumerated".format(description)
            )
        validate_owned_directory(directory, description)
        validate_owned_directory(root, "canonical result directory")
        return result
    finally:
        close_owned_directory(directory)


def validate_atomic_destination(directory_descriptor, name, description):
    try:
        metadata = os.stat(
            name, dir_fd=directory_descriptor, follow_symlinks=False
        )
    except FileNotFoundError:
        return
    except OSError as error:
        raise CrossoverError(
            "cannot inspect {}: {}".format(description, error)
        )
    if (not stat.S_ISREG(metadata.st_mode) or
            metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
            stat.S_IMODE(metadata.st_mode) != 0o600):
        raise CrossoverError(
            "refusing unsafe existing {}".format(description)
        )
    return metadata


def link_descriptor_noreplace(
        source_descriptor, destination_directory_descriptor,
        destination_name, description):
    """Link one exact anonymous inode to a new name without replacement."""
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = libc.linkat
    linkat.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_int,
    )
    linkat.restype = ctypes.c_int
    result = linkat(
            source_descriptor, ctypes.c_char_p(b""),
            destination_directory_descriptor,
            os.fsencode(destination_name), 0x1000)  # AT_EMPTY_PATH
    if result != 0:
        error_number = ctypes.get_errno()
        if error_number in (errno.ENOENT, errno.EPERM):
            # Unprivileged O_TMPFILE publication uses the kernel's exact fd
            # symlink with AT_SYMLINK_FOLLOW.
            result = linkat(
                -100, os.fsencode("/proc/self/fd/{}".format(
                    source_descriptor
                )),
                destination_directory_descriptor,
                os.fsencode(destination_name), 0x400,
            )
            if result == 0:
                return
            error_number = ctypes.get_errno()
        raise CrossoverError(
            "cannot publish {} without replacement: {}".format(
                description, os.strerror(error_number)
            )
        )


def rename_noreplace(
        directory_descriptor, source_name, destination_name, description):
    """Atomically rename one staged name only if its destination is absent."""
    libc = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = libc.renameat2
    except AttributeError:
        raise CrossoverError(
            "{} requires renameat2(RENAME_NOREPLACE)".format(description)
        )
    renameat2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        directory_descriptor, os.fsencode(source_name),
        directory_descriptor, os.fsencode(destination_name),
        0x1,  # RENAME_NOREPLACE
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise CrossoverError(
            "cannot atomically publish {} without replacement: {}".format(
                description, os.strerror(error_number)
            )
        )


def rename_exchange(
        directory_descriptor, first_name, second_name, description):
    """Atomically exchange two existing names in one retained directory."""
    libc = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = libc.renameat2
    except AttributeError:
        raise CrossoverError(
            "{} requires renameat2(RENAME_EXCHANGE)".format(description)
        )
    renameat2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        directory_descriptor, os.fsencode(first_name),
        directory_descriptor, os.fsencode(second_name),
        0x2,  # RENAME_EXCHANGE
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise CrossoverError(
            "cannot atomically exchange {}: {}".format(
                description, os.strerror(error_number)
            )
        )


def named_inode_identity(directory_descriptor, name, description):
    """Return one no-follow name identity, or None when it is absent."""
    try:
        metadata = os.stat(
            name, dir_fd=directory_descriptor, follow_symlinks=False
        )
    except FileNotFoundError:
        return None
    except OSError as error:
        raise CrossoverError(
            "cannot reconcile {}: {}".format(description, error)
        )
    return metadata.st_dev, metadata.st_ino


def reconcile_atomic_publication(
        directory_descriptor, destination_name, staged_name,
        temporary_descriptor, old_descriptor, description):
    """Classify state after an interrupted renameat2 publication syscall.

    Python can deliver a BaseException after renameat2 committed but before
    the caller executes its next bytecode.  Bind both names to the retained
    new/prior file descriptors so that this boundary is either recognized as
    committed and rolled back, or rejected as ambiguous without unlinking an
    unrelated inode.
    """
    new_metadata = os.fstat(temporary_descriptor)
    new_identity = (new_metadata.st_dev, new_metadata.st_ino)
    old_identity = None
    if old_descriptor is not None:
        old_metadata = os.fstat(old_descriptor)
        old_identity = (old_metadata.st_dev, old_metadata.st_ino)
    destination_identity = named_inode_identity(
        directory_descriptor, destination_name,
        description + " destination",
    )
    staged_identity = named_inode_identity(
        directory_descriptor, staged_name, description + " staging name"
    )
    if old_identity is None:
        if (destination_identity is None and
                staged_identity == new_identity):
            return "uncommitted"
        if (destination_identity == new_identity and
                staged_identity is None):
            return "new"
    else:
        if (destination_identity == old_identity and
                staged_identity == new_identity):
            return "uncommitted"
        if (destination_identity == new_identity and
                staged_identity == old_identity):
            return "replacement"
    raise CrossoverError(
        "{} left an ambiguous destination/staging inode state".format(
            description
        )
    )


def write_result_bytes(
        root, relative, encoded, description="atomic JSON destination",
        replace=False, precommit_hook=None, postlink_hook=None,
        postcommit_hook=None):
    """Publish bytes relative to a retained root from one exact O_TMPFILE inode.

    This binds the bytes and inode through the atomic commit and its immediate
    retained-FD verification.  Mode 0600 intentionally leaves evidence
    owner-writable for resumable campaigns, so no claim is made against a
    malicious same-UID writer that mutates it after this function returns;
    readers therefore repeat exact-inode byte/metadata validation.
    """
    if not isinstance(encoded, bytes):
        raise CrossoverError("{} payload is not bytes".format(description))
    parts = result_relative_parts(relative, description)
    if len(parts) == 1:
        directory = duplicate_owned_directory(
            root, description + " directory"
        )
    else:
        directory = open_result_directory(
            root, Path(*parts[:-1]), description + " directory", True
        )
    name = parts[-1]
    temporary_descriptor = None
    old_descriptor = None
    old_payload = None
    staged_name = None
    committed_state = "uncommitted"
    try:
        existing = validate_atomic_destination(
            directory["descriptor"], name, description
        )
        if existing is not None and not replace:
            raise CrossoverError(
                "{} already exists; refusing to replace it".format(description)
            )
        if existing is not None:
            old_descriptor = os.open(
                name,
                os.O_RDONLY | os.O_NOFOLLOW |
                getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NONBLOCK", 0),
                dir_fd=directory["descriptor"],
            )
            old_metadata = os.fstat(old_descriptor)
            if metadata_signature(old_metadata) != metadata_signature(
                    existing):
                raise CrossoverError(
                    "{} changed while retaining prior evidence".format(
                        description
                    )
                )
            old_payload = read_descriptor_bytes(
                old_descriptor, old_metadata.st_size,
                description + " prior evidence",
            )
        if not hasattr(os, "O_TMPFILE"):
            raise CrossoverError(
                "{} requires O_TMPFILE exact-inode publication".format(
                    description
                )
            )
        temporary_descriptor = os.open(
            ".", os.O_RDWR | os.O_TMPFILE | getattr(os, "O_CLOEXEC", 0),
            0o600, dir_fd=directory["descriptor"],
        )
        os.fchmod(temporary_descriptor, 0o600)
        offset = 0
        while offset < len(encoded):
            try:
                written = os.write(
                    temporary_descriptor, encoded[offset:]
                )
            except InterruptedError:
                continue
            if written <= 0:
                raise CrossoverError(
                    "atomic JSON write made no progress"
                )
            offset += written
        os.fsync(temporary_descriptor)
        metadata = os.fstat(temporary_descriptor)
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or
                metadata.st_nlink != 0 or
                stat.S_IMODE(metadata.st_mode) != 0o600 or
                metadata.st_size != len(encoded)):
            raise CrossoverError(
                "atomic JSON temporary metadata changed"
            )
        expected_sha256 = digest_bytes(encoded)
        if read_descriptor_bytes(
                temporary_descriptor, len(encoded),
                description + " staged payload") != encoded:
            raise CrossoverError(
                "{} staged payload differs from requested bytes".format(
                    description
                )
            )
        validate_owned_directory(
            directory, description + " directory"
        )
        current = validate_atomic_destination(
            directory["descriptor"], name, description
        )
        if ((existing is None) != (current is None) or
                (existing is not None and
                 (existing.st_dev, existing.st_ino) !=
                 (current.st_dev, current.st_ino))):
            raise CrossoverError(
                "{} changed before publication".format(description)
            )
        for unused_attempt in range(32):
            candidate = ".leo2-stage-{}-{}-{}".format(
                os.getpid(), threading.get_ident(),
                os.urandom(12).hex(),
            )
            try:
                link_descriptor_noreplace(
                    temporary_descriptor, directory["descriptor"],
                    candidate, description + " staged inode",
                )
                staged_name = candidate
                break
            except CrossoverError as error:
                if "File exists" not in str(error):
                    raise
        if staged_name is None:
            raise CrossoverError(
                "{} could not allocate a unique staging name".format(
                    description
                )
            )
        if postlink_hook is not None:
            postlink_hook(
                temporary_descriptor, directory["descriptor"], staged_name
            )
        staged = os.stat(
            staged_name, dir_fd=directory["descriptor"],
            follow_symlinks=False,
        )
        retained_before_commit = os.fstat(temporary_descriptor)
        retained_payload = read_descriptor_bytes(
            temporary_descriptor, len(encoded),
            description + " linked staged payload",
        )
        retained_after_read = os.fstat(temporary_descriptor)
        if (retained_payload != encoded or
                digest_bytes(retained_payload) != expected_sha256 or
                (staged.st_dev, staged.st_ino) !=
                (retained_before_commit.st_dev,
                 retained_before_commit.st_ino) or
                retained_before_commit.st_nlink != 1 or
                retained_after_read.st_nlink != 1 or
                metadata_signature(retained_before_commit) !=
                metadata_signature(retained_after_read) or
                stat.S_IMODE(staged.st_mode) != 0o600 or
                staged.st_size != len(encoded)):
            raise CrossoverError(
                "{} linked staged inode or payload changed".format(
                    description
                )
            )
        if precommit_hook is not None:
            precommit_hook()
        current = validate_atomic_destination(
            directory["descriptor"], name, description
        )
        if ((existing is None) != (current is None) or
                (existing is not None and
                 metadata_signature(existing) !=
                 metadata_signature(current))):
            raise CrossoverError(
                "{} changed immediately before publication".format(
                    description
                )
            )
        publication_staged_name = staged_name
        commit_error = None
        try:
            if current is None:
                rename_noreplace(
                    directory["descriptor"], staged_name, name, description
                )
                committed_state = "new"
                staged_name = None
            else:
                rename_exchange(
                    directory["descriptor"], staged_name, name, description
                )
                committed_state = "replacement"
        except BaseException as error:
            try:
                reconciled_state = reconcile_atomic_publication(
                    directory["descriptor"], name,
                    publication_staged_name, temporary_descriptor,
                    old_descriptor, description + " interrupted commit",
                )
            except BaseException as reconciliation_error:
                raise CrossoverError(
                    "{} commit failed: {}; state reconciliation also "
                    "failed: {}: {}".format(
                        description, error,
                        type(reconciliation_error).__name__,
                        reconciliation_error,
                    )
                ) from error
            if reconciled_state == "uncommitted":
                raise
            committed_state = reconciled_state
            staged_name = (
                None if reconciled_state == "new"
                else publication_staged_name
            )
            commit_error = error
        try:
            if commit_error is not None:
                raise commit_error
            if postcommit_hook is not None:
                postcommit_hook(
                    temporary_descriptor, directory["descriptor"], name
                )
            final = os.stat(
                name, dir_fd=directory["descriptor"],
                follow_symlinks=False,
            )
            retained = os.fstat(temporary_descriptor)
            final_payload = read_descriptor_bytes(
                temporary_descriptor, len(encoded),
                description + " published payload",
            )
            retained_after_final_read = os.fstat(temporary_descriptor)
            if (not stat.S_ISREG(final.st_mode) or
                    final.st_uid != os.getuid() or final.st_nlink != 1 or
                    stat.S_IMODE(final.st_mode) != 0o600 or
                    final.st_size != len(encoded) or
                    (final.st_dev, final.st_ino) !=
                    (retained.st_dev, retained.st_ino) or
                    final_payload != encoded or
                    digest_bytes(final_payload) != expected_sha256 or
                    metadata_signature(retained) !=
                    metadata_signature(retained_after_final_read)):
                raise CrossoverError(
                    "{} final inode or payload differs from its retained "
                    "temporary".format(description)
                )
            os.fsync(directory["descriptor"])
            validate_owned_directory(
                directory, description + " directory"
            )
            if committed_state == "replacement":
                staged_old = os.stat(
                    staged_name, dir_fd=directory["descriptor"],
                    follow_symlinks=False,
                )
                retained_old = os.fstat(old_descriptor)
                if ((staged_old.st_dev, staged_old.st_ino) !=
                        (retained_old.st_dev, retained_old.st_ino) or
                        read_descriptor_bytes(
                            old_descriptor, len(old_payload),
                            description + " retained prior evidence",
                        ) != old_payload):
                    raise CrossoverError(
                        "{} prior-evidence rollback inode changed".format(
                            description
                        )
                    )
                try:
                    os.unlink(
                        staged_name, dir_fd=directory["descriptor"]
                    )
                except BaseException:
                    # A signal may be delivered after unlinkat(2) committed
                    # but before Python reports its return.  If the exact
                    # retained old-evidence name is already absent, rollback
                    # is no longer possible and reporting failure would
                    # falsely describe the already durable replacement.
                    try:
                        os.stat(
                            staged_name,
                            dir_fd=directory["descriptor"],
                            follow_symlinks=False,
                        )
                    except FileNotFoundError:
                        staged_name = None
                        committed_state = "success"
                    else:
                        raise
                else:
                    staged_name = None
            committed_state = "success"
        except BaseException as postcommit_error:
            rollback_errors = []
            try:
                if committed_state == "replacement":
                    try:
                        rename_exchange(
                            directory["descriptor"], staged_name, name,
                            description + " rollback",
                        )
                    except BaseException:
                        rollback_state = reconcile_atomic_publication(
                            directory["descriptor"], name, staged_name,
                            temporary_descriptor, old_descriptor,
                            description + " interrupted rollback",
                        )
                        if rollback_state != "uncommitted":
                            raise
                    committed_state = "rolled-back"
                    restored = os.stat(
                        name, dir_fd=directory["descriptor"],
                        follow_symlinks=False,
                    )
                    retained_old = os.fstat(old_descriptor)
                    if ((restored.st_dev, restored.st_ino) !=
                            (retained_old.st_dev, retained_old.st_ino) or
                            stat.S_IMODE(restored.st_mode) != 0o600 or
                            read_descriptor_bytes(
                                old_descriptor, len(old_payload),
                                description + " restored prior evidence",
                            ) != old_payload):
                        raise CrossoverError(
                            "{} rollback did not restore prior evidence".format(
                                description
                            )
                        )
                    os.fsync(directory["descriptor"])
                elif committed_state == "new":
                    retained = os.fstat(temporary_descriptor)
                    published_identity = named_inode_identity(
                        directory["descriptor"], name,
                        description + " new rollback destination",
                    )
                    if (published_identity !=
                            (retained.st_dev, retained.st_ino)):
                        raise CrossoverError(
                            "{} new publication changed before rollback".format(
                                description
                            )
                        )
                    try:
                        os.unlink(name, dir_fd=directory["descriptor"])
                    except BaseException:
                        if named_inode_identity(
                                directory["descriptor"], name,
                                description + " interrupted new rollback",
                        ) is not None:
                            raise
                    os.fsync(directory["descriptor"])
                    committed_state = "rolled-back"
            except BaseException as rollback_error:
                rollback_errors.append(rollback_error)
            if rollback_errors:
                raise CrossoverError(
                    "{} post-commit failure could not be rolled back: {}; "
                    "rollback: {}".format(
                        description, postcommit_error,
                        "; ".join(
                            "{}: {}".format(type(error).__name__, error)
                            for error in rollback_errors
                        ),
                    )
                ) from postcommit_error
            raise
        # The new name and its containing directory were fully verified and
        # fsynced before prior evidence was unlinked.  Cleanup durability is
        # best-effort from this point: reporting failure would be false,
        # because the replacement is already the committed evidence.
        try:
            os.fsync(directory["descriptor"])
        except OSError:
            pass
    except OSError as error:
        raise CrossoverError(
            "cannot publish {}: {}".format(description, error)
        )
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors = []
        if staged_name is not None:
            try:
                staged = os.stat(
                    staged_name, dir_fd=directory["descriptor"],
                    follow_symlinks=False,
                )
                retained = (
                    os.fstat(temporary_descriptor)
                    if temporary_descriptor is not None else None
                )
                if (retained is not None and
                        (staged.st_dev, staged.st_ino) ==
                        (retained.st_dev, retained.st_ino)):
                    os.unlink(
                        staged_name, dir_fd=directory["descriptor"]
                    )
                    os.fsync(directory["descriptor"])
            except FileNotFoundError:
                pass
            except BaseException as error:
                # A same-UID actor can mutate owner-writable evidence after
                # any validation returns.  Never unlink an ambiguous name.
                cleanup_errors.append(error)
        for label, descriptor in (
                ("temporary", temporary_descriptor),
                ("prior evidence", old_descriptor)):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except BaseException as error:
                    cleanup_errors.append(CrossoverError(
                        "{} descriptor close failed: {}".format(label, error)
                    ))
        try:
            close_owned_directory(directory)
        except BaseException as error:
            cleanup_errors.append(CrossoverError(
                "publication directory close failed: {}".format(error)
            ))
        if cleanup_errors and committed_state != "success":
            detail = "; ".join(
                "{}: {}".format(type(error).__name__, error)
                for error in cleanup_errors
            )
            if active_error is not None:
                raise CrossoverError(
                    "{} failed: {}; cleanup also failed: {}".format(
                        description, active_error, detail
                    )
                ) from active_error
            raise CrossoverError(
                "{} cleanup failed: {}".format(description, detail)
            )


def write_bytes_atomic(path, encoded):
    path = checked_path(path, "atomic JSON destination")
    if path.name in ("", ".", ".."):
        raise CrossoverError("atomic JSON destination has an unsafe name")
    directory = owned_canonical_directory(
        path.parent, "atomic JSON destination directory"
    )
    try:
        write_result_bytes(
            directory, path.name, encoded,
            "atomic JSON destination", replace=True,
        )
    finally:
        close_owned_directory(directory)


def canonical_json_bytes(value):
    return (
        json.dumps(
            value, indent=2, sort_keys=True, ensure_ascii=True,
            allow_nan=False
        ) + "\n"
    ).encode("utf-8")


def atomic_write_json(path, value):
    write_bytes_atomic(Path(path), canonical_json_bytes(value))


def write_result_json(root, relative, value, description, replace=False):
    write_result_bytes(
        root, relative, canonical_json_bytes(value), description, replace
    )


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


def open_canonical_git_snapshot():
    """Retain and seal the exact Git bytes used for one source snapshot."""
    descriptor = None
    snapshot = None
    snapshot_transferred = False
    try:
        descriptor = os.open(
            str(CANONICAL_GIT),
            os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
        )
        before = os.fstat(descriptor)
        current = os.lstat(str(CANONICAL_GIT))
        value = read_descriptor_bytes(
            descriptor, MAX_RAW_JSON_BYTES, "canonical Git executable"
        )
        after = os.fstat(descriptor)
        if (CANONICAL_GIT.resolve() != CANONICAL_GIT or
                not stat.S_ISREG(before.st_mode) or
                not before.st_mode & 0o111 or
                stat.S_IMODE(before.st_mode) & 0o022 or
                metadata_signature(before) != metadata_signature(after) or
                (before.st_dev, before.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "canonical Git executable is not one stable executable inode"
            )
        snapshot = create_sealed_executable_snapshot(
            value, "canonical Git executable"
        )
        snapshot["path"] = str(CANONICAL_GIT)
        snapshot_transferred = True
        return snapshot
    except OSError as error:
        snapshot_transferred = False
        raise CrossoverError(
            "canonical Git executable is unavailable: {}".format(error)
        )
    except BaseException:
        snapshot_transferred = False
        raise
    finally:
        try:
            if descriptor is not None:
                os.close(descriptor)
        except BaseException:
            snapshot_transferred = False
            if snapshot is not None:
                close_sealed_snapshot(snapshot)
            raise
        if snapshot is not None and not snapshot_transferred:
            close_sealed_snapshot(snapshot)


def close_canonical_git_snapshot(snapshot):
    if snapshot is not None:
        snapshot.pop("path", None)
        close_sealed_snapshot(snapshot)


def git_command_bytes(source, arguments, description, git_tool=None):
    owns_git_tool = git_tool is None
    if git_tool is None:
        git_tool = open_canonical_git_snapshot()
    elif (not isinstance(git_tool, dict) or
            set(git_tool) != {"descriptor", "path", "sha256", "size"}):
        raise CrossoverError(
            "{} exact Git tool contract is invalid".format(description)
        )
    try:
        with tempfile.TemporaryDirectory(
                prefix="leopard2-git-command-") as temporary:
            log_root = Path(temporary)
            log_root.chmod(0o700)
            stdout_path = log_root / "stdout.log"
            stderr_path = log_root / "stderr.log"
            validate_sealed_executable_snapshot(
                git_tool, git_tool["sha256"], "canonical Git executable"
            )
            completed = run_command(
                [
                    str(Path("/proc/self/fd") /
                        str(GIT_EXECUTABLE_DESCRIPTOR))
                ] + list(arguments),
                source,
                stdout_path, stderr_path, 30, dict(GIT_ENVIRONMENT),
                inherited_descriptors=(git_tool["descriptor"],),
                descriptor_remaps={
                    GIT_EXECUTABLE_DESCRIPTOR: git_tool["descriptor"],
                },
            )
            validate_sealed_executable_snapshot(
                git_tool, git_tool["sha256"], "canonical Git executable"
            )
            stdout = read_path_owned_regular(
                stdout_path, MAX_RETAINED_LOG_BYTES,
                description + " stdout",
            )
            stderr = read_path_owned_regular(
                stderr_path, MAX_RETAINED_LOG_BYTES,
                description + " stderr",
            )
    except (OSError, subprocess.SubprocessError) as error:
        raise CrossoverError("cannot {}: {}".format(description, error))
    finally:
        if owns_git_tool:
            close_canonical_git_snapshot(git_tool)
    if completed["returncode"] != 0 or completed["timed_out"]:
        raise CrossoverError(
            "cannot {}: {}".format(
                description,
                normalized_output(stderr).decode(
                    "utf-8", errors="replace"
                ).strip()
            )
        )
    return stdout


def canonical_git_identity(git_tool=None):
    owns_git_tool = git_tool is None
    if git_tool is None:
        git_tool = open_canonical_git_snapshot()
    try:
        validate_sealed_executable_snapshot(
            git_tool, git_tool["sha256"], "canonical Git executable"
        )
        sha256 = git_tool["sha256"]
        path = git_tool["path"]
    finally:
        if owns_git_tool:
            close_canonical_git_snapshot(git_tool)
    return {
        "environment": dict(GIT_ENVIRONMENT),
        "path": path,
        "sha256": sha256,
    }


def parse_stage_zero_entries(value, description):
    entries = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            header, encoded_path = encoded.split(b"\t", 1)
            mode, object_id, stage = header.decode("ascii").split(" ")
            relative = encoded_path.decode("utf-8")
        except (UnicodeError, ValueError) as error:
            raise CrossoverError(
                "{} contains an invalid entry: {}".format(description, error)
            )
        if (stage != "0" or not re.fullmatch(r"[0-9a-f]{40}", object_id) or
                mode not in ("100644", "100755", "120000", "160000") or
                not relative or relative in entries):
            raise CrossoverError(
                "{} contains a conflicted, unsupported, or duplicate entry "
                "{!r}".format(description, relative)
            )
        entries[relative] = (mode, object_id)
    return entries


def parse_head_tree_entries(value, description):
    entries = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            header, encoded_path = encoded.split(b"\t", 1)
            mode, object_type, object_id = header.decode("ascii").split(" ")
            relative = encoded_path.decode("utf-8")
        except (UnicodeError, ValueError) as error:
            raise CrossoverError(
                "{} contains an invalid entry: {}".format(description, error)
            )
        if (object_type not in ("blob", "commit") or
                not re.fullmatch(r"[0-9a-f]{40}", object_id) or
                mode not in ("100644", "100755", "120000", "160000") or
                not relative or relative in entries or
                (mode == "160000") != (object_type == "commit")):
            raise CrossoverError(
                "{} contains an unsupported or duplicate entry {!r}".format(
                    description, relative
                )
            )
        entries[relative] = (mode, object_id)
    return entries


def parse_index_tags(value, description):
    tags = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            tag = encoded[:1].decode("ascii")
            separator = encoded[1:2]
            relative = encoded[2:].decode("utf-8")
        except UnicodeError as error:
            raise CrossoverError(
                "{} contains an invalid tag: {}".format(description, error)
            )
        if (separator != b" " or not relative or relative in tags or
                tag != "H"):
            raise CrossoverError(
                "{} found non-default index flag tag {!r} for {!r}".format(
                    description, tag, relative
                )
            )
        tags[relative] = tag
    return tags


def repository_source_controls(
        source, stage_entries, expected_git, git_tool):
    replacements = normalized_output(git_command_bytes(
        source, (
            "for-each-ref", "--format=%(refname)", "refs/replace/",
        ), "enumerate Git replacement refs", git_tool
    )).decode("utf-8", errors="strict").rstrip("\n")
    replacement_refs = replacements.splitlines() if replacements else []
    if replacement_refs:
        raise CrossoverError(
            "Git replacement refs are forbidden: {}".format(
                ", ".join(replacement_refs)
            )
        )

    graft_name = normalized_output(git_command_bytes(
        source, ("rev-parse", "--git-path", "info/grafts"),
        "locate legacy Git grafts", git_tool
    )).decode("utf-8", errors="strict").strip()
    graft_path = Path(graft_name)
    if not graft_path.is_absolute():
        graft_path = (source / graft_path).resolve()
    if graft_path.exists():
        raise CrossoverError(
            "legacy Git grafts are forbidden: {}".format(graft_path)
        )

    visible_tags = parse_index_tags(git_command_bytes(
        source, ("ls-files", "-v", "-z"),
        "enumerate assume-unchanged/skip-worktree index flags", git_tool
    ), "git ls-files -v")
    fsmonitor_tags = parse_index_tags(git_command_bytes(
        source, ("ls-files", "-f", "-z"),
        "enumerate fsmonitor-valid index flags", git_tool
    ), "git ls-files -f")
    if (set(visible_tags) != set(stage_entries) or
            set(fsmonitor_tags) != set(stage_entries)):
        raise CrossoverError(
            "Git index flag enumeration differs from the stage-0 closure"
        )

    head_entries = parse_head_tree_entries(git_command_bytes(
        source, ("ls-tree", "-r", "-z", expected_git["head"]),
        "enumerate the exact captured HEAD tree", git_tool
    ), "git ls-tree HEAD")
    return {
        "head_index_match": head_entries == stage_entries,
        "head_tree_sha256": digest_value({
            name: list(head_entries[name]) for name in sorted(head_entries)
        }),
        "index_entry_count": len(stage_entries),
        "index_flags": "all ls-files -v/-f tags are canonical H",
        "index_sha256": digest_value({
            name: list(stage_entries[name]) for name in sorted(stage_entries)
        }),
        "legacy_grafts_absent": True,
        "replace_refs": [],
    }


def git_blob_object_id(value):
    header = b"blob " + str(len(value)).encode("ascii") + b"\0"
    return hashlib.sha1(header + value).hexdigest()


def _tracked_git_closure_bound(source, git_tool, expected_git):
    source = Path(source).resolve()
    output = git_command_bytes(
        source, ("ls-files", "-s", "-z"),
        "enumerate tracked source closure", git_tool
    )
    stage_entries = parse_stage_zero_entries(
        output, "Git stage-0 source closure"
    )
    controls = repository_source_controls(
        source, stage_entries, expected_git, git_tool
    )
    files = {}
    for relative, (mode, object_id) in stage_entries.items():
        path = source / relative
        try:
            if mode == "160000":
                resolved = path.resolve(strict=True)
                common = os.path.commonpath((str(source), str(resolved)))
                if common != str(source) or not resolved.is_dir():
                    raise CrossoverError(
                        "tracked gitlink escapes or is not initialized: {}".format(
                            relative
                        )
                    )
                submodule_git = git_source_identity(
                    resolved, False, git_tool
                )
                head = submodule_git["head"]
                status_lines = submodule_git["status"]
                if head != object_id:
                    raise CrossoverError(
                        "initialized submodule {} HEAD {} differs from index "
                        "{}".format(relative, head, object_id)
                    )
                closure = _tracked_git_closure_bound(
                    resolved, git_tool, submodule_git
                )
                record = {
                    "index_mode": mode,
                    "index_object": object_id,
                    "submodule": {
                        "digest": closure["digest"],
                        "files": closure["files"],
                        "head": head,
                        "repository_controls":
                            closure["repository_controls"],
                        "status": (
                            status_lines
                        ),
                        "worktree_clean": not bool(status_lines),
                    },
                }
            elif mode in ("100644", "100755", "120000"):
                resolved = path.resolve(strict=True)
                common = os.path.commonpath((str(source), str(resolved)))
                if common != str(source):
                    raise CrossoverError(
                        "tracked path escapes the source root: {}".format(relative)
                    )
                if mode == "120000":
                    if not path.is_symlink():
                        raise CrossoverError(
                            "tracked symlink changed type: {}".format(relative)
                        )
                    before = path.lstat()
                    if not stat.S_ISLNK(before.st_mode):
                        raise CrossoverError(
                            "tracked symlink changed type while hashing: "
                            "{}".format(relative)
                        )
                    value = os.readlink(str(path)).encode("utf-8")
                    after = path.lstat()
                    worktree_mode = "120000"
                else:
                    if not path.is_file() or path.is_symlink():
                        raise CrossoverError(
                            "tracked file changed type: {}".format(relative)
                        )
                    before = path.lstat()
                    if not stat.S_ISREG(before.st_mode):
                        raise CrossoverError(
                            "tracked file changed type while hashing: "
                            "{}".format(relative)
                        )
                    value = path.read_bytes()
                    after = path.lstat()
                    worktree_mode = (
                        "100755" if before.st_mode & 0o111 else "100644"
                    )
                if (
                    before.st_dev, before.st_ino, before.st_mode,
                    before.st_size, before.st_mtime_ns, before.st_ctime_ns
                ) != (
                    after.st_dev, after.st_ino, after.st_mode,
                    after.st_size, after.st_mtime_ns, after.st_ctime_ns
                ):
                    raise CrossoverError(
                        "tracked source changed while hashing: {}".format(
                            relative
                        )
                    )
                record = {
                    "index_mode": mode,
                    "index_object": object_id,
                    "worktree_git_object": git_blob_object_id(value),
                    "worktree_mode": worktree_mode,
                    "worktree_sha256": digest_bytes(value),
                }
            else:
                raise CrossoverError(
                    "unsupported Git index mode {} for {}".format(mode, relative)
                )
        except CrossoverError:
            raise
        except (OSError, UnicodeError) as error:
            raise CrossoverError(
                "cannot hash tracked source {}: {}".format(relative, error)
            )
        files[relative] = record
    if not files:
        raise CrossoverError("Git tracked source closure is empty")
    ordered = {name: files[name] for name in sorted(files)}
    material = {
        "files": ordered,
        "repository_controls": controls,
    }
    final_stage_entries = parse_stage_zero_entries(
        git_command_bytes(
            source, ("ls-files", "-s", "-z"),
            "re-enumerate tracked source closure", git_tool,
        ),
        "final Git stage-0 source closure",
    )
    final_controls = repository_source_controls(
        source, final_stage_entries, expected_git, git_tool
    )
    if (final_stage_entries != stage_entries or
            canonical_bytes(final_controls) != canonical_bytes(controls)):
        raise CrossoverError(
            "Git index or repository controls changed across the tracked "
            "closure"
        )
    final_git = git_source_identity(source, False, git_tool)
    if canonical_bytes(final_git) != canonical_bytes(expected_git):
        raise CrossoverError(
            "Git source identity changed across the tracked closure"
        )
    return {
        "digest": digest_value(material),
        "files": ordered,
        "repository_controls": controls,
    }


def tracked_git_closure(source, git_tool=None, expected_git=None):
    owns_git_tool = git_tool is None
    if git_tool is None:
        git_tool = open_canonical_git_snapshot()
    try:
        if expected_git is None:
            expected_git = git_source_identity(source, False, git_tool)
        return _tracked_git_closure_bound(
            source, git_tool, expected_git
        )
    finally:
        if owns_git_tool:
            close_canonical_git_snapshot(git_tool)


def source_fingerprint(source):
    with SOURCE_IDENTITY_LOCK:
        return tracked_git_closure(Path(source).resolve())


def git_source_identity(source, require_clean, git_tool=None):
    def git(arguments):
        return normalized_output(git_command_bytes(
            source, arguments, "inspect Git source identity", git_tool
        )).decode(
            "utf-8", errors="strict"
        ).strip()

    revision_lines = git((
        "rev-parse", "HEAD", "HEAD^{tree}",
    )).splitlines()
    if len(revision_lines) != 2:
        raise CrossoverError(
            "Git HEAD/tree query did not return exactly two identities"
        )
    head, tree = revision_lines
    branch = git(("branch", "--show-current"))
    status = normalized_output(git_command_bytes(
        source, (
            "status", "--porcelain=v1", "--untracked-files=all",
            "--ignore-submodules=none",
        ), "read Git worktree status", git_tool
    )).decode("utf-8", errors="strict").rstrip("\n")
    for label, value in (("HEAD", head), ("tree", tree)):
        if not re.fullmatch(r"[0-9a-f]{40}", value):
            raise CrossoverError("Git {} is not a full lowercase SHA-1".format(label))
    if require_clean and status:
        raise CrossoverError(
            "authoritative source tree must be clean in Git: {}".format(
                status.replace("\n", "; ")
            )
        )
    return {
        "branch": branch or None,
        "head": head,
        "status": status.splitlines() if status else [],
        "worktree_clean": not bool(status),
        "tree": tree,
    }


def source_identity(source, require_clean):
    """Capture one cross-bound Git/index/worktree state with one Git image.

    The process lock serializes this runner's own capture threads, and exact
    pre/post Git identities reject a state transition that overlaps capture.
    Like any unprivileged filesystem observation, this is not a claim that a
    malicious same-UID actor cannot perform a complete A-B-A between probes or
    mutate the repository after return; authoritative callers recapture and
    compare the identity at campaign boundaries.
    """
    source = Path(source).resolve()
    with SOURCE_IDENTITY_LOCK:
        git_tool = open_canonical_git_snapshot()
        try:
            git_before = git_source_identity(
                source, require_clean, git_tool
            )
            result = _tracked_git_closure_bound(
                source, git_tool, git_before
            )
            git_after = git_source_identity(
                source, require_clean, git_tool
            )
            if canonical_bytes(git_after) != canonical_bytes(git_before):
                raise CrossoverError(
                    "Git source identity changed across source capture"
                )
            result["git"] = git_before
            result["git_tool"] = canonical_git_identity(git_tool)
            return validate_source_state(
                result, "captured source identity", require_clean
            )
        finally:
            close_canonical_git_snapshot(git_tool)


def validate_source_state(value, description, require_clean):
    if (not isinstance(value, dict) or
            set(value) != {
                "digest", "files", "git", "git_tool",
                "repository_controls"}):
        raise CrossoverError("{} has an invalid source identity".format(description))
    files = value.get("files")
    git = value.get("git")
    git_tool = value.get("git_tool")
    controls = value.get("repository_controls")
    if (not isinstance(files, dict) or not files or
            value.get("digest") != digest_value({
                "files": files,
                "repository_controls": controls,
            })):
        raise CrossoverError(
            "{} has an invalid tracked-source closure".format(description)
        )
    submodules_clean = True
    tracked_files_clean = True

    def validate_controls(control_value, entries, label):
        expected_keys = {
            "head_index_match", "head_tree_sha256", "index_entry_count",
            "index_flags", "index_sha256", "legacy_grafts_absent",
            "replace_refs",
        }
        index_material = {
            name: [entries[name]["index_mode"], entries[name]["index_object"]]
            for name in sorted(entries)
        }
        if (not isinstance(control_value, dict) or
                set(control_value) != expected_keys or
                type(control_value.get("head_index_match")) is not bool or
                type(control_value.get("legacy_grafts_absent")) is not bool or
                control_value["legacy_grafts_absent"] is not True or
                control_value.get("replace_refs") != [] or
                control_value.get("index_flags") !=
                "all ls-files -v/-f tags are canonical H" or
                type(control_value.get("index_entry_count")) is not int or
                control_value["index_entry_count"] != len(entries) or
                control_value.get("index_sha256") !=
                digest_value(index_material) or
                not isinstance(control_value.get("head_tree_sha256"), str) or
                not re.fullmatch(
                    r"[0-9a-f]{64}", control_value["head_tree_sha256"]) or
                (control_value["head_index_match"] and
                 control_value["head_tree_sha256"] !=
                 control_value["index_sha256"])):
            raise CrossoverError(
                "{} has invalid repository controls for {}".format(
                    description, label
                )
            )
        return control_value["head_index_match"]

    def validate_files(entries, prefix):
        nonlocal submodules_clean, tracked_files_clean
        if not isinstance(entries, dict) or not entries:
            raise CrossoverError(
                "{} {} closure is empty".format(description, prefix)
            )
        if list(entries) != sorted(entries):
            raise CrossoverError(
                "{} {} closure is not canonically ordered".format(
                    description, prefix
                )
            )
        for name, record in entries.items():
            label = "{}{}".format(prefix, name)
            if (not isinstance(name, str) or not name or
                    not isinstance(record, dict) or
                    record.get("index_mode") not in (
                        "100644", "100755", "120000", "160000") or
                    not isinstance(record.get("index_object"), str) or
                    not re.fullmatch(
                        r"[0-9a-f]{40}", record["index_object"])):
                raise CrossoverError(
                    "{} has invalid tracked entry {}".format(description, label)
                )
            if record["index_mode"] == "160000":
                if set(record) != {
                        "index_mode", "index_object", "submodule"}:
                    raise CrossoverError(
                        "{} has invalid gitlink {}".format(description, label)
                    )
                submodule = record.get("submodule")
                if (not isinstance(submodule, dict) or set(submodule) != {
                        "digest", "files", "head", "repository_controls",
                        "status", "worktree_clean"} or
                        submodule.get("head") != record["index_object"] or
                        not isinstance(submodule.get("status"), list) or
                        any(not isinstance(line, str) or not line
                            for line in submodule["status"]) or
                        type(submodule.get("worktree_clean")) is not bool or
                        submodule["worktree_clean"] !=
                        (not bool(submodule["status"])) or
                        submodule.get("digest") != digest_value({
                            "files": submodule.get("files"),
                            "repository_controls":
                                submodule.get("repository_controls"),
                        })):
                    raise CrossoverError(
                        "{} has invalid submodule state {}".format(
                            description, label
                        )
                    )
                submodule_controls_clean = validate_controls(
                    submodule["repository_controls"],
                    submodule["files"], label
                )
                submodules_clean = (
                    submodules_clean and submodule["worktree_clean"] and
                    submodule_controls_clean
                )
                validate_files(submodule["files"], label + "/")
            else:
                if (set(record) != {
                        "index_mode", "index_object", "worktree_git_object",
                        "worktree_mode", "worktree_sha256"} or
                        not isinstance(record.get("worktree_sha256"), str) or
                        not re.fullmatch(
                            r"[0-9a-f]{64}", record["worktree_sha256"]) or
                        not isinstance(
                            record.get("worktree_git_object"), str) or
                        not re.fullmatch(
                            r"[0-9a-f]{40}",
                            record["worktree_git_object"]) or
                        record.get("worktree_mode") not in (
                            "100644", "100755", "120000")):
                    raise CrossoverError(
                        "{} has invalid file state {}".format(
                            description, label
                        )
                    )
                tracked_files_clean = (
                    tracked_files_clean and
                    record["worktree_git_object"] ==
                        record["index_object"] and
                    record["worktree_mode"] == record["index_mode"]
                )

    controls_clean = validate_controls(controls, files, "top-level source")
    validate_files(files, "")
    if (not isinstance(git, dict) or set(git) != {
            "branch", "head", "status", "tree", "worktree_clean"} or
            not isinstance(git.get("status"), list) or
            any(not isinstance(line, str) or not line for line in git["status"]) or
            type(git.get("worktree_clean")) is not bool or
            git["worktree_clean"] != (not bool(git["status"])) or
            not isinstance(git.get("branch"), (str, type(None))) or
            not all(isinstance(git.get(name), str) and
                    re.fullmatch(r"[0-9a-f]{40}", git[name])
                    for name in ("head", "tree"))):
        raise CrossoverError(
            "{} has an invalid exact Git identity".format(description)
        )
    if require_clean and (
            git["worktree_clean"] is not True or not submodules_clean or
            not tracked_files_clean or not controls_clean):
        raise CrossoverError("{} is not a clean Git source".format(description))
    if canonical_bytes(git_tool) != canonical_bytes(canonical_git_identity()):
        raise CrossoverError(
            "{} used a different canonical Git executable".format(description)
        )
    return value


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
    matches = []
    seen = set()
    for candidate in candidates:
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            resolved = candidate.resolve()
            if str(resolved) not in seen:
                seen.add(str(resolved))
                matches.append(resolved)
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise CrossoverError(
            "benchmark executable lookup for {} is ambiguous: {}".format(
                backend, ", ".join(str(item) for item in matches)
            )
        )
    raise CrossoverError(
        "benchmark executable for {} was not found; checked {}".format(
            backend, ", ".join(str(item) for item in candidates)
        )
    )


def build_configuration_contract(attestation_schema):
    if attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA:
        return (
            BUILD_CONFIGURATION_FILE_SCHEMA,
            BUILD_CONFIGURATION_VARIABLES,
            BUILD_CONFIGURATION_EXPERIMENT_SELECTORS,
        )
    if attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V4:
        return (
            BUILD_CONFIGURATION_FILE_SCHEMA_V6,
            BUILD_CONFIGURATION_VARIABLES_V6,
            BUILD_CONFIGURATION_EXPERIMENT_SELECTORS_V6,
        )
    if attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2:
        return (
            BUILD_CONFIGURATION_FILE_SCHEMA_V2,
            BUILD_CONFIGURATION_VARIABLES_V2,
            tuple(
                variable
                for variable in BUILD_CONFIGURATION_EXPERIMENT_SELECTORS
                if variable in BUILD_CONFIGURATION_VARIABLES_V2
            ),
        )
    if attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V3:
        return (
            BUILD_CONFIGURATION_FILE_SCHEMA_V3,
            BUILD_CONFIGURATION_VARIABLES_V3,
            tuple(
                variable
                for variable in BUILD_CONFIGURATION_EXPERIMENT_SELECTORS
                if variable in BUILD_CONFIGURATION_VARIABLES_V3
            ),
        )
    raise CrossoverError(
        "benchmark effective-configuration attestation schema is invalid"
    )


def build_configuration_digest(
        entries, variables=BUILD_CONFIGURATION_VARIABLES):
    if not isinstance(entries, dict):
        raise CrossoverError(
            "effective CMake configuration is not an object"
        )
    if tuple(variables) not in (
            BUILD_CONFIGURATION_VARIABLES,
            BUILD_CONFIGURATION_VARIABLES_V6,
            BUILD_CONFIGURATION_VARIABLES_V3,
            BUILD_CONFIGURATION_VARIABLES_V2):
        raise CrossoverError(
            "effective CMake configuration contract is invalid"
        )
    if set(entries) != set(variables):
        raise CrossoverError(
            "effective CMake configuration has unexpected variables"
        )
    for variable in variables:
        attested_text(
            entries[variable],
            "effective CMake configuration value {}".format(variable)
        )
    if entries["LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] not in (
            "0", "1", "2"):
        raise CrossoverError(
            "effective CMake configuration has an invalid "
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"
        )
    material = "".join(
        "{}={}\n".format(variable, entries[variable])
        for variable in variables
    )
    return digest_bytes(material.encode("utf-8"))


def read_build_configuration_attestation(path):
    configuration_path = checked_path(
        path, "effective CMake configuration path"
    )
    try:
        encoded = configuration_path.read_bytes()
        text = encoded.decode("utf-8")
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError(
            "cannot read effective CMake configuration {}: {}".format(
                path, error
            )
        )
    if "\0" in text or "\r" in text:
        raise CrossoverError(
            "effective CMake configuration contains a forbidden delimiter"
        )
    if not text.endswith("\n"):
        raise CrossoverError(
            "effective CMake configuration is not newline-terminated"
        )
    # The framing delimiter is exactly LF.  Values are subsequently passed
    # through attested_text(), which also rejects Unicode line/paragraph
    # separators and invisible/control format categories.
    lines = text[:-1].split("\n")
    if not lines:
        raise CrossoverError(
            "effective CMake configuration has an invalid schema"
        )
    schema_prefix = "schema="
    file_schema = (
        lines[0][len(schema_prefix):]
        if lines[0].startswith(schema_prefix) else ""
    )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA:
        attestation_schema = BUILD_CONFIGURATION_ATTESTATION_SCHEMA
        variables = BUILD_CONFIGURATION_VARIABLES
    elif file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V6:
        attestation_schema = BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V4
        variables = BUILD_CONFIGURATION_VARIABLES_V6
    elif file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V3:
        attestation_schema = BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V3
        variables = BUILD_CONFIGURATION_VARIABLES_V3
    elif file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V2:
        attestation_schema = BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2
        variables = BUILD_CONFIGURATION_VARIABLES_V2
    else:
        raise CrossoverError(
            "effective CMake configuration has an invalid schema"
        )
    expected_count = 2 + len(variables)
    if len(lines) != expected_count:
        raise CrossoverError(
            "effective CMake configuration has {} lines, expected {}".format(
                len(lines), expected_count
            )
        )
    if not lines[1].startswith("sha256="):
        raise CrossoverError(
            "effective CMake configuration omits its digest"
        )
    declared_digest = lines[1][len("sha256="):]
    entries = {}
    for variable, line in zip(
            variables, lines[2:]):
        prefix = "{}=".format(variable)
        if not line.startswith(prefix):
            raise CrossoverError(
                "effective CMake configuration expected {}, found {!r}".format(
                    variable, line
                )
            )
        entries[variable] = line[len(prefix):]
    actual_digest = build_configuration_digest(entries, variables)
    if (not re.fullmatch(r"[0-9a-f]{64}", declared_digest) or
            declared_digest != actual_digest):
        raise CrossoverError(
            "effective CMake configuration digest is invalid"
        )
    return {
        "entries": entries,
        "path": str(checked_resolve(
            configuration_path, "effective CMake configuration path"
        )),
        "schema": attestation_schema,
        "sha256": actual_digest,
    }


def validate_build_configuration_attestation(
        value, expected_path=None,
        expected_schema=BUILD_CONFIGURATION_ATTESTATION_SCHEMA):
    expected_keys = {"entries", "path", "schema", "sha256"}
    if not isinstance(value, dict) or set(value) != expected_keys:
        raise CrossoverError(
            "benchmark effective-configuration attestation is invalid"
        )
    if not isinstance(value.get("path"), str):
        raise CrossoverError(
            "benchmark effective-configuration attestation is invalid"
        )
    try:
        attested_text(value.get("schema"), "attestation schema")
        if value.get("schema") != expected_schema:
            raise CrossoverError(
                "benchmark effective-configuration attestation schema "
                "does not match its enclosing contract"
            )
        attested_path = checked_absolute_path(
            value.get("path"), "attestation path"
        )
        attested_text(value.get("sha256"), "attestation digest")
        file_schema, variables, selectors = \
            build_configuration_contract(value.get("schema"))
        del file_schema
        if (value.get("schema") ==
                BUILD_CONFIGURATION_ATTESTATION_SCHEMA):
            canonical_selectors = BUILD_CONFIGURATION_CANONICAL_SELECTORS
        elif (value.get("schema") ==
              BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V4):
            canonical_selectors = BUILD_CONFIGURATION_CANONICAL_SELECTORS_V6
        elif (value.get("schema") ==
                BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V3):
            canonical_selectors = BUILD_CONFIGURATION_CANONICAL_SELECTORS_V3
        else:
            canonical_selectors = BUILD_CONFIGURATION_CANONICAL_SELECTORS_V2
        valid = (
            re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None and
            value["sha256"] == build_configuration_digest(
                value.get("entries"), variables
            ) and
            set(canonical_selectors) == set(selectors) and
            all(value["entries"].get(name) == expected
                for name, expected in canonical_selectors.items())
        )
    except CrossoverError:
        raise
    if not valid:
        raise CrossoverError(
            "benchmark effective-configuration attestation is invalid"
        )
    if (expected_path is not None and checked_resolve(
            attested_path, "attestation path") != checked_resolve(
                expected_path, "expected attestation path")):
        raise CrossoverError(
            "benchmark effective-configuration path is unexpected"
        )
    return value["sha256"]


def parse_selected_cmake_cache(cache_bytes, prefixes, description):
    try:
        cache_lines = cache_bytes.decode("utf-8").split("\n")
    except (AttributeError, UnicodeError) as error:
        raise CrossoverError("cannot read {}: {}".format(description, error))
    entries = {}
    for line in cache_lines:
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        typed_key, value = line.split("=", 1)
        key = typed_key.split(":", 1)[0]
        if not any(
                key == prefix or key.startswith(prefix + "_")
                for prefix in prefixes):
            continue
        if typed_key.count(":") != 1:
            raise CrossoverError(
                "{} contains malformed selected cache key {!r}".format(
                    description, typed_key
                )
            )
        key, entry_type = typed_key.split(":", 1)
        if (not key or not re.fullmatch(r"[A-Za-z0-9_.+-]+", key) or
                entry_type not in CMAKE_CACHE_ENTRY_TYPES or
                key in entries):
            raise CrossoverError(
                "{} contains malformed or duplicate selected cache key "
                "{!r}".format(description, key)
            )
        allowed_types = CMAKE_CACHE_REQUIRED_ENTRY_TYPES.get(key)
        if allowed_types is not None and entry_type not in allowed_types:
            raise CrossoverError(
                "{} cache key {} has type {}, expected one of {}".format(
                    description, key, entry_type,
                    sorted(allowed_types)
                )
            )
        entries[key] = attested_text(
            value, "CMake cache value {}".format(key)
        )
    return entries


def cmake_build_metadata(executable):
    """Fingerprint the CMake inputs nearest an already-built benchmark."""
    executable = checked_path(executable, "benchmark executable")
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
    except (OSError, UnicodeError) as error:
        raise CrossoverError("cannot read {}: {}".format(cache, error))
    prefixes = (
        "CMAKE_BINARY_DIR", "CMAKE_BUILD_TYPE", "CMAKE_CONFIGURATION_TYPES",
        "CMAKE_C_COMPILER", "CMAKE_C_FLAGS", "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_FLAGS", "CMAKE_GENERATOR",
        "CMAKE_HOME_DIRECTORY",
        "ENABLE_OPENMP", "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
        "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS",
        "LEO2_BENCHMARK_GIT_EXECUTABLE",
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION",
        "LEO2_BUILD_FUZZERS", "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA",
        *BUILD_CONFIGURATION_EXPERIMENT_SELECTORS,
    )
    entries = parse_selected_cmake_cache(
        cache_bytes, prefixes, str(cache)
    )
    build_root = cache.parent
    build_configuration_path = (
        build_root / BUILD_CONFIGURATION_RELATIVE_PATH
    )
    build_configuration_attestation = (
        read_build_configuration_attestation(build_configuration_path)
    )
    validate_build_configuration_attestation(
        build_configuration_attestation,
        build_configuration_path,
        BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
    )
    (build_configuration_file_schema, unused_variables,
     build_configuration_selectors) = build_configuration_contract(
        build_configuration_attestation["schema"])
    del unused_variables
    if (entries.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") !=
            build_configuration_file_schema or
            entries.get(
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"
            ) != build_configuration_attestation["sha256"] or
            any(entries.get(variable) !=
                build_configuration_attestation["entries"].get(variable)
                for variable in
                build_configuration_selectors)):
        raise CrossoverError(
            "effective CMake configuration differs from its cache binding"
        )
    extra_files = {}
    for relative in (
            "compile_commands.json", "CMakeFiles/CMakeConfigureLog.yaml",
            "CMakeFiles/CMakeOutput.log", "build.ninja",
            BUILD_CONFIGURATION_RELATIVE_PATH,
            "generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h",
            "CMakeFiles/leopard.dir/flags.make",
            "CMakeFiles/bench_leopard2_direct_encode.dir/DependInfo.cmake",
            "CMakeFiles/bench_leopard2_direct_encode.dir/flags.make",
            "CMakeFiles/bench_leopard2_direct_encode.dir/link.txt",
            "CMakeFiles/bench_leopard2_direct_encode.dir/"
            "bench/leopard2/direct_encode_benchmark.cpp.o",
            "CMakeFiles/leopard_test_hooks.dir/flags.make",
            "CMakeFiles/leopard_test_hooks.dir/link.txt",
            "CMakeFiles/leopard2_backend_avx2_test_hooks.dir/flags.make",
            "libleopard.a", "libleopard_test_hooks.a"):
        path = build_root / relative
        if path.is_file():
            extra_files[relative] = digest_bytes(path.read_bytes())
    for directory in (
            "CMakeFiles/leopard_test_hooks.dir",
            "CMakeFiles/leopard2_backend_ssse3_test_hooks.dir",
            "CMakeFiles/leopard2_backend_avx2_test_hooks.dir",
            "CMakeFiles/leopard2_backend_avx512_test_hooks.dir"):
        root = build_root / directory
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.o")):
            relative = str(path.relative_to(build_root))
            extra_files[relative] = digest_bytes(path.read_bytes())
    executable_stat = executable.stat()
    return {
        "binding_scope": (
            "exact executable, effective-configuration sidecar, CMake cache, "
            "available compile database/generator graph, direct benchmark "
            "object/link recipe, and present test-hook archive/object hashes; "
            "embedded clean Git commit/tree attestation and full tracked-"
            "worktree hashes are validated separately"
        ),
        "build_root": str(checked_resolve(
            build_root, "CMake build root"
        )),
        "cmake_cache": str(checked_resolve(cache, "CMake cache")),
        "cmake_cache_sha256": digest_bytes(cache_bytes),
        "effective_configuration_attestation":
            build_configuration_attestation,
        "entries": entries,
        "executable": {
            "mtime_ns": executable_stat.st_mtime_ns,
            "path": str(checked_resolve(
                executable, "benchmark executable"
            )),
            "sha256": digest_bytes(executable.read_bytes()),
            "size": executable_stat.st_size,
        },
        "extra_file_sha256": extra_files,
    }


def cmake_configuration_types(entries):
    if not isinstance(entries, dict):
        raise CrossoverError(
            "effective CMake configuration is not an object"
        )
    build_type = attested_text(
        entries.get("CMAKE_BUILD_TYPE"),
        "effective CMAKE_BUILD_TYPE"
    )
    encoded_types = attested_text(
        entries.get("CMAKE_CONFIGURATION_TYPES"),
        "effective CMAKE_CONFIGURATION_TYPES"
    )
    if not encoded_types:
        return build_type, ()
    configuration_types = tuple(encoded_types.split(";"))
    if (any(not value for value in configuration_types) or
            len(set(configuration_types)) != len(configuration_types)):
        raise CrossoverError(
            "effective CMAKE_CONFIGURATION_TYPES is malformed"
        )
    return build_type, configuration_types


def cmake_generator_is_multi_config(entries):
    if not isinstance(entries, dict):
        raise CrossoverError(
            "effective CMake configuration is not an object"
        )
    generator = attested_text(
        entries.get("CMAKE_GENERATOR"), "effective CMake generator"
    )
    if not generator:
        raise CrossoverError("effective CMake generator is empty")
    return (
        generator == "Xcode" or
        generator.startswith("Visual Studio") or
        "Multi-Config" in generator
    )


def validate_embedded_build_type(
        entries, embedded_build_type, authoritative):
    embedded_build_type = attested_text(
        embedded_build_type, "embedded benchmark build type"
    )
    build_type, configuration_types = cmake_configuration_types(entries)
    if cmake_generator_is_multi_config(entries):
        if not configuration_types:
            raise CrossoverError(
                "multi-config CMake generator omits configuration types"
            )
        if (not embedded_build_type or
                embedded_build_type not in configuration_types):
            raise CrossoverError(
                "benchmark multi-config build type {!r} is not one of "
                "{!r}".format(embedded_build_type, configuration_types)
            )
    elif embedded_build_type != build_type:
        raise CrossoverError(
            "benchmark single-config build type {!r} differs from {!r}".format(
                embedded_build_type, build_type
            )
        )
    if authoritative and embedded_build_type != "Release":
        raise CrossoverError(
            "authoritative benchmark executable must be the Release "
            "configuration"
        )
    return embedded_build_type


def validate_build_source_binding(
        metadata, source, source_state, backend, require_fresh):
    if not isinstance(metadata, dict):
        raise CrossoverError("CMake build provenance is not an object")
    entries = metadata.get("entries")
    extra = metadata.get("extra_file_sha256")
    executable = metadata.get("executable")
    configuration_attestation = metadata.get(
        "effective_configuration_attestation"
    )
    if not all(isinstance(value, dict) for value in (
            entries, extra, executable, configuration_attestation)):
        raise CrossoverError("CMake build provenance is incomplete")
    build_root = metadata.get("build_root")
    if not isinstance(build_root, str):
        raise CrossoverError("CMake build provenance omits its build root")
    try:
        build_root_path = checked_absolute_path(
            build_root, "CMake build provenance root"
        )
    except CrossoverError:
        raise CrossoverError("CMake build provenance omits its build root")
    configuration_sha256 = validate_build_configuration_attestation(
        configuration_attestation,
        build_root_path / BUILD_CONFIGURATION_RELATIVE_PATH
    )
    (configuration_file_schema, unused_variables,
     configuration_selectors) = build_configuration_contract(
        configuration_attestation["schema"])
    del unused_variables
    configuration_entries = configuration_attestation["entries"]
    if (entries.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") !=
            configuration_file_schema or
            entries.get(
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"
            ) != configuration_sha256 or
            any(entries.get(variable) != configuration_entries.get(variable)
                for variable in
                configuration_selectors)):
        raise CrossoverError(
            "effective CMake configuration differs from its cache binding"
        )
    home_directory = entries.get("CMAKE_HOME_DIRECTORY")
    if not isinstance(home_directory, str) or not home_directory:
        raise CrossoverError("CMake build source root differs from --source")
    try:
        home_path = checked_absolute_path(
            home_directory, "CMake build source root"
        )
        source_path = checked_resolve(source, "--source")
    except CrossoverError:
        raise CrossoverError("CMake build source root differs from --source")
    if checked_resolve(
            home_path, "CMake build source root") != source_path:
        raise CrossoverError("CMake build source root differs from --source")
    if configuration_entries.get(
            "LEO2_BACKEND_VARIANT", "").lower() != backend:
        raise CrossoverError(
            "CMake LEO2_BACKEND_VARIANT is {!r}, expected {!r}".format(
                configuration_entries.get("LEO2_BACKEND_VARIANT"), backend
            )
        )
    if configuration_entries.get(
            "LEO2_BUILD_BENCHMARKS", "").upper() not in ("1", "ON", "TRUE"):
        raise CrossoverError("CMake build did not explicitly enable benchmarks")
    if not require_fresh:
        return
    git = source_state.get("git", {})
    if git.get("worktree_clean") is not True:
        raise CrossoverError("authoritative build requires a clean Git source tree")
    configured_build_type, configured_types = cmake_configuration_types(
        configuration_entries
    )
    is_multi_config = cmake_generator_is_multi_config(
        configuration_entries
    )
    if ((is_multi_config and "Release" not in configured_types) or
            (not is_multi_config and configured_build_type != "Release")):
        raise CrossoverError(
            "authoritative build must provide the Release configuration"
        )
    if configuration_entries.get(
            "LEO2_BUILD_TESTS", "").upper() not in ("1", "ON", "TRUE"):
        raise CrossoverError("authoritative build did not explicitly enable tests")
    embedded_git = configuration_entries.get(
        "LEO2_BENCHMARK_GIT_EXECUTABLE"
    )
    if not isinstance(embedded_git, str) or not embedded_git:
        raise CrossoverError(
            "authoritative benchmark attestation did not use canonical "
            "{}".format(CANONICAL_GIT)
        )
    try:
        embedded_git_path = checked_absolute_path(
            embedded_git, "attested Git executable"
        )
    except CrossoverError:
        raise CrossoverError(
            "authoritative benchmark attestation did not use canonical "
            "{}".format(CANONICAL_GIT)
        )
    if checked_resolve(
            embedded_git_path, "attested Git executable") != CANONICAL_GIT:
        raise CrossoverError(
            "authoritative benchmark attestation did not use canonical "
            "{}".format(CANONICAL_GIT)
        )
    required = {
        "compile_commands.json",
        "CMakeFiles/bench_leopard2_direct_encode.dir/"
        "bench/leopard2/direct_encode_benchmark.cpp.o",
        BUILD_CONFIGURATION_RELATIVE_PATH,
        "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h",
    }
    missing = sorted(required - set(extra))
    if missing:
        raise CrossoverError(
            "authoritative build provenance omits {}".format(", ".join(missing))
        )
    if not any(name in extra for name in (
            "build.ninja",
            "CMakeFiles/bench_leopard2_direct_encode.dir/link.txt")):
        raise CrossoverError("authoritative build provenance omits the link graph")
    if not any(name in extra for name in (
            "libleopard.a", "libleopard_test_hooks.a")):
        raise CrossoverError("authoritative build provenance omits the linked archive")
    def source_mtimes(root, entries):
        for relative, record in entries.items():
            path = root / relative
            yield path.lstat().st_mtime_ns
            if record["index_mode"] == "160000":
                yield from source_mtimes(
                    path.resolve(), record["submodule"]["files"]
                )

    try:
        newest_source_mtime = max(
            source_mtimes(source, source_state["files"])
        )
    except OSError as error:
        raise CrossoverError(
            "cannot establish authoritative source freshness: {}".format(error)
        )
    if executable.get("mtime_ns", 0) < newest_source_mtime:
        raise CrossoverError(
            "authoritative executable predates the captured source closure; "
            "use a fresh dedicated build directory"
        )


def frozen_executable_identity(
        backend, origin_executable, executable_sha256, build_metadata,
        source_state):
    return {
        "backend": backend,
        "build_metadata": build_metadata,
        "executable_sha256": executable_sha256,
        "origin_executable": str(Path(origin_executable).resolve()),
        "source_fingerprint": source_state,
    }


def validate_frozen_executable(
        artifact, build_metadata=None, source_state=None, result_root=None):
    expected_keys = {
        "artifact_id", "backend", "directory", "executable",
        "executable_sha256", "provenance", "provenance_sha256", "schema",
    }
    if (not isinstance(artifact, dict) or set(artifact) != expected_keys or
            artifact.get("schema") != FROZEN_EXECUTABLE_SCHEMA):
        raise CrossoverError("frozen executable artifact has an invalid schema")
    artifact_id = artifact.get("artifact_id")
    if (not isinstance(artifact_id, str) or
            not re.fullmatch(r"[0-9a-f]{64}", artifact_id)):
        raise CrossoverError("frozen executable artifact ID is invalid")
    if result_root is None:
        directory = Path(artifact["directory"]).resolve()
        executable = Path(artifact["executable"]).resolve()
        provenance_path = Path(artifact["provenance"]).resolve()
    else:
        directory = Path(os.path.abspath(artifact["directory"]))
        executable = Path(os.path.abspath(artifact["executable"]))
        provenance_path = Path(os.path.abspath(artifact["provenance"]))
    if executable.parent != directory or provenance_path.parent != directory:
        raise CrossoverError("frozen executable paths escape the artifact directory")
    executable_held = None
    provenance_held = None
    artifact_directory = None
    try:
        if result_root is None:
            executable_bytes = executable.read_bytes()
            provenance_bytes = provenance_path.read_bytes()
            directory_mode = directory.stat().st_mode
            executable_mode = executable.stat().st_mode
            provenance_mode = provenance_path.stat().st_mode
        else:
            try:
                directory_relative = Path(os.path.abspath(
                    os.fspath(directory)
                )).relative_to(result_root["path"])
                executable_relative = Path(os.path.abspath(
                    os.fspath(executable)
                )).relative_to(result_root["path"])
                provenance_relative = Path(os.path.abspath(
                    os.fspath(provenance_path)
                )).relative_to(result_root["path"])
            except ValueError:
                raise CrossoverError(
                    "frozen executable artifact escapes the held result root"
                )
            artifact_directory = open_existing_result_directory(
                result_root, directory_relative,
                "frozen executable artifact directory",
                required_mode=0o555,
            )
            directory_mode = os.fstat(
                artifact_directory["descriptor"]
            ).st_mode
            executable_held = open_result_regular_held(
                result_root, executable_relative, MAX_RAW_JSON_BYTES,
                "frozen executable", 0o555,
                directory_required_mode=0o555,
            )
            provenance_held = open_result_regular_held(
                result_root, provenance_relative, MAX_RAW_JSON_BYTES,
                "frozen executable provenance", 0o444,
                directory_required_mode=0o555,
            )
            executable_bytes = read_held_regular(
                executable_held, MAX_RAW_JSON_BYTES, "frozen executable"
            )
            provenance_bytes = read_held_regular(
                provenance_held, MAX_RAW_JSON_BYTES,
                "frozen executable provenance",
            )
            executable_mode = os.fstat(executable_held["descriptor"]).st_mode
            provenance_mode = os.fstat(provenance_held["descriptor"]).st_mode
    except OSError as error:
        raise CrossoverError("cannot read frozen executable artifact: {}".format(error))
    finally:
        if provenance_held is not None:
            close_log(provenance_held)
        if executable_held is not None:
            close_log(executable_held)
        if artifact_directory is not None:
            close_owned_directory(artifact_directory)
    if (stat.S_IMODE(directory_mode) != 0o555 or
            stat.S_IMODE(executable_mode) != 0o555 or
            stat.S_IMODE(provenance_mode) != 0o444):
        raise CrossoverError(
            "frozen executable artifact modes are not exactly "
            "0555/0555/0444"
        )
    executable_sha256 = digest_bytes(executable_bytes)
    if (executable_sha256 != artifact.get("executable_sha256") or
            not executable_mode & 0o111):
        raise CrossoverError("frozen executable hash or mode is invalid")
    if digest_bytes(provenance_bytes) != artifact.get("provenance_sha256"):
        raise CrossoverError("frozen executable provenance hash is invalid")
    provenance = decode_json_bytes(provenance_bytes, str(provenance_path))
    provenance_keys = {
        "artifact_id", "backend", "build_metadata", "executable",
        "origin_executable", "schema", "source_fingerprint",
    }
    if (not isinstance(provenance, dict) or set(provenance) != provenance_keys or
            provenance.get("schema") != FROZEN_EXECUTABLE_SCHEMA or
            provenance.get("artifact_id") != artifact_id or
            provenance.get("backend") != artifact.get("backend") or
            provenance.get("executable") != {
                "name": executable.name, "sha256": executable_sha256,
            }):
        raise CrossoverError("frozen executable provenance is inconsistent")
    origin = provenance.get("origin_executable")
    if (not isinstance(origin, dict) or set(origin) != {"path", "sha256"} or
            origin.get("sha256") != executable_sha256):
        raise CrossoverError("frozen executable origin provenance is invalid")
    identity = frozen_executable_identity(
        artifact["backend"], origin["path"], executable_sha256,
        provenance["build_metadata"], provenance["source_fingerprint"]
    )
    if digest_value(identity) != artifact_id:
        raise CrossoverError("frozen executable artifact ID does not match provenance")
    if (build_metadata is not None and
            canonical_bytes(provenance["build_metadata"]) !=
            canonical_bytes(build_metadata)):
        raise CrossoverError("frozen executable build provenance is stale")
    if (source_state is not None and
            canonical_bytes(provenance["source_fingerprint"]) !=
            canonical_bytes(source_state)):
        raise CrossoverError("frozen executable source provenance is stale")
    return artifact


def freeze_executable(
        result_dir, backend, executable, build_metadata, source_state,
        result_root=None):
    executable = Path(executable).resolve()
    try:
        executable_bytes = executable.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot read executable to freeze {}: {}".format(executable, error)
        )
    current_build_metadata = cmake_build_metadata(executable)
    if canonical_bytes(current_build_metadata) != canonical_bytes(build_metadata):
        raise CrossoverError(
            "origin executable or CMake metadata changed before freezing"
        )
    executable_sha256 = digest_bytes(executable_bytes)
    identity = frozen_executable_identity(
        backend, executable, executable_sha256, build_metadata, source_state
    )
    artifact_id = digest_value(identity)
    frozen_root = Path(result_dir).resolve() / "frozen-executables"
    artifact_dir = frozen_root / "{}-{}".format(backend, artifact_id[:24])
    frozen_name = "bench_leopard2_direct_encode"
    frozen_path = artifact_dir / frozen_name
    provenance_path = artifact_dir / "provenance.json"
    provenance = {
        "artifact_id": artifact_id,
        "backend": backend,
        "build_metadata": build_metadata,
        "executable": {
            "name": frozen_name,
            "sha256": executable_sha256,
        },
        "origin_executable": {
            "path": str(executable),
            "sha256": executable_sha256,
        },
        "schema": FROZEN_EXECUTABLE_SCHEMA,
        "source_fingerprint": source_state,
    }

    if result_root is None:
        result_root = owned_canonical_directory(
            result_dir, "frozen executable result directory"
        )
        owns_result_root = True
    else:
        owns_result_root = False
    frozen_directory = open_result_directory(
        result_root, "frozen-executables",
        "frozen executable parent", True,
    )
    artifact_relative = Path("frozen-executables") / artifact_dir.name
    artifact_directory = None
    try:
        try:
            artifact_directory = open_existing_result_directory(
                result_root, artifact_relative,
                "frozen executable artifact", required_mode=0o555,
            )
        except FileNotFoundError:
            try:
                os.mkdir(
                    artifact_dir.name, 0o700,
                    dir_fd=frozen_directory["descriptor"],
                )
            except FileExistsError:
                pass
            artifact_directory = open_existing_result_directory(
                result_root, artifact_relative,
                "frozen executable artifact",
            )
            try:
                executable_file = open_exclusive_owned_file(
                    artifact_directory, frozen_name,
                    "frozen executable",
                )
                try:
                    write_descriptor_all(
                        executable_file["descriptor"], executable_bytes
                    )
                    os.fsync(executable_file["descriptor"])
                    os.fchmod(executable_file["descriptor"], 0o555)
                finally:
                    close_log(executable_file)
                write_result_json(
                    artifact_directory, "provenance.json", provenance,
                    "frozen executable provenance", replace=False,
                )
                os.chmod(
                    "provenance.json", 0o444,
                    dir_fd=artifact_directory["descriptor"],
                    follow_symlinks=False,
                )
                os.fchmod(artifact_directory["descriptor"], 0o555)
                os.fsync(artifact_directory["descriptor"])
                os.fsync(frozen_directory["descriptor"])
            except BaseException:
                # Leave a fail-closed partial artifact for explicit inspection.
                # An existing artifact was also treated this way historically.
                raise
        finally:
            if artifact_directory is not None:
                close_owned_directory(artifact_directory)
    finally:
        close_owned_directory(frozen_directory)
        if owns_result_root:
            close_owned_directory(result_root)

    try:
        origin_bytes_after = executable.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot re-read origin executable after freezing: {}".format(error)
        )
    if (origin_bytes_after != executable_bytes or
            canonical_bytes(cmake_build_metadata(executable)) !=
            canonical_bytes(build_metadata)):
        raise CrossoverError(
            "origin executable or CMake metadata changed while freezing"
        )

    provenance_sha256 = digest_bytes(canonical_json_bytes(provenance))
    artifact = {
        "artifact_id": artifact_id,
        "backend": backend,
        "directory": str(artifact_dir),
        "executable": str(frozen_path),
        "executable_sha256": executable_sha256,
        "provenance": str(provenance_path),
        "provenance_sha256": provenance_sha256,
        "schema": FROZEN_EXECUTABLE_SCHEMA,
    }
    return validate_frozen_executable(
        artifact, build_metadata, source_state,
        None if owns_result_root else result_root,
    )


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


def historical_avx2_grid():
    """Corrected high-profile AVX2 cells, including the omitted 17.88x cell."""
    values = [
        # Exact historical cell omitted by the invalid 2026-07-28 screen.
        cell(
            "candidate_historical_exact", "avx2", 2, 16, "high",
            "gf8", 4096, 1, "0"
        ),
        # The nine retained screen cells, now interpreted only as campaign shape.
        cell("excluded_historical_single_side_control", "avx2", 2, 1, "high", "gf8", 1024, 1, "0"),
        cell("excluded_historical_single_side_control", "avx2", 3, 1, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 4, 4, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 8, 8, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 8, 8, "high", "gf8", 65536, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 4, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 65536, 1, "0"),
        # One-coordinate neighbors make K, R, byte size, mask position, and Q
        # independently visible around the historical cell.
        cell("historical_neighbor_k", "avx2", 3, 16, "high", "gf8", 4096, 1, "0"),
        cell("historical_neighbor_r", "avx2", 2, 15, "high", "gf8", 4096, 1, "0"),
        cell(
            "historical_neighbor_bytes_lower", "avx2", 2, 16, "high", "gf8",
            4032, 1, "0"
        ),
        cell(
            "historical_neighbor_bytes_upper", "avx2", 2, 16, "high", "gf8",
            4160, 1, "0"
        ),
        cell(
            "historical_neighbor_sparse_output", "avx2", 2, 16, "high",
            "gf8", 4096, 1, "15"
        ),
        cell(
            "historical_neighbor_q2", "avx2", 2, 16, "high", "gf8",
            4096, 2, "0-1"
        ),
    ]
    identities = {digest_value(value) for value in values}
    if len(identities) != len(values):
        raise CrossoverError("historical AVX2 campaign contains duplicate cells")
    return values


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
    if mode in AUTHORITATIVE_COMMANDS:
        return ("direct", "transform", "transform", "direct") * abba_rounds
    # Alternating the two-cell screening order prevents one path from always
    # receiving the colder process launch while keeping the order deterministic.
    return (("direct", "transform") if int(job_id[:8], 16) % 2 == 0
            else ("transform", "direct"))


def job_identity(
        cell_value, executable, executable_artifact, executable_sha256,
        build_metadata, source_state, machine, settings):
    return {
        "cell": cell_value,
        "build_metadata": build_metadata,
        "executable": str(executable),
        "executable_artifact": executable_artifact,
        "executable_sha256": executable_sha256,
        "machine": machine,
        "settings": settings,
        "source_identity": source_state,
    }


def make_jobs(
        cells, executables, build_metadata, source_state, machine, settings,
        executable_artifacts=None):
    executable_artifacts = executable_artifacts or {}
    jobs = []
    for cell_value in cells:
        executable = executables[cell_value["backend"]]
        executable_sha256 = digest_bytes(executable.read_bytes())
        executable_artifact = executable_artifacts.get(cell_value["backend"])
        identity = job_identity(
            cell_value, executable, executable_artifact, executable_sha256,
            build_metadata[cell_value["backend"]], source_state, machine,
            settings
        )
        configuration_id = digest_value(identity)
        job_id = configuration_id[:24]
        jobs.append({
            "cell": cell_value,
            "build_metadata": build_metadata[cell_value["backend"]],
            "configuration_id": configuration_id,
            "executable": str(executable),
            "executable_artifact": executable_artifact,
            "executable_sha256": executable_sha256,
            "invocation_order": list(invocation_order(
                settings["mode"], job_id, settings["abba_rounds"])),
            "job_id": job_id,
            "seed": stable_seed(cell_value),
            "source_identity": source_state,
        })
    return jobs


def benchmark_argv(job, timed_mode, raw_path, settings):
    item = job["cell"]
    benchmark = settings["benchmark"]
    argv = [
        (
            str(Path("/proc/self/fd") / str(EXECUTABLE_DESCRIPTOR))
            if job.get("executable_artifact") is not None
            else job["executable"]
        ),
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
    if settings["mode"] in AUTHORITATIVE_COMMANDS:
        argv = [
            str(Path("/proc/self/fd") / str(
                TASKSET_EXECUTABLE_DESCRIPTOR
            )),
            "-c", str(settings["pin_cpu"]),
        ] + argv
    return argv


_SUBREAPER_LOCK = threading.Lock()
_SUBREAPER_ENABLED = False
_PR_SET_CHILD_SUBREAPER = 36
_PR_GET_CHILD_SUBREAPER = 37


def ensure_child_subreaper():
    """Make detached command descendants reapable by this coordinator.

    Authoritative execution is Linux-only in practice, but screening remains
    useful on other POSIX systems.  A new session still contains ordinary
    descendants there; Linux additionally adopts children that deliberately
    leave that session so the timeout path can reap them rather than merely
    turning them into zombies owned by PID 1.
    """
    global _SUBREAPER_ENABLED
    if sys.platform != "linux":
        return False
    with _SUBREAPER_LOCK:
        if _SUBREAPER_ENABLED:
            return True
        libc = ctypes.CDLL(None, use_errno=True)
        state = ctypes.c_int(-1)
        if libc.prctl(
                _PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) != 0:
            error_number = ctypes.get_errno()
            raise CrossoverError(
                "cannot enable command descendant reaping: {}".format(
                    os.strerror(error_number)
                )
            )
        if libc.prctl(
                _PR_GET_CHILD_SUBREAPER, ctypes.byref(state),
                0, 0, 0) != 0 or state.value != 1:
            error_number = ctypes.get_errno()
            raise CrossoverError(
                "cannot verify command descendant reaping: {}".format(
                    os.strerror(error_number) if error_number else
                    "kernel did not retain child-subreaper state"
                )
            )
        _SUBREAPER_ENABLED = True
        return True


def parse_proc_stat(value):
    """Return PID/session identity fields from one Linux /proc stat record."""
    try:
        open_parenthesis = value.index("(")
        close_parenthesis = value.rindex(")")
        pid = int(value[:open_parenthesis].strip())
        fields = value[close_parenthesis + 1:].split()
        if len(fields) < 20:
            return None
        return {
            "pid": pid,
            "ppid": int(fields[1]),
            "session": int(fields[3]),
            "starttime_ticks": int(fields[19]),
            "state": fields[0],
        }
    except (IndexError, ValueError):
        return None


def open_proc_record(pid):
    """Return one process record retained by its exact /proc task inode."""
    if sys.platform != "linux":
        return None
    proc_descriptor = None
    task_descriptor = None
    stat_descriptor = None
    transferred_record = None
    try:
        flags = (
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
            getattr(os, "O_CLOEXEC", 0)
        )
        proc_descriptor = os.open("/proc", flags)
        task_descriptor = os.open(
            str(pid), flags, dir_fd=proc_descriptor
        )
        task_metadata = os.fstat(task_descriptor)
        current = os.stat(
            str(pid), dir_fd=proc_descriptor, follow_symlinks=False
        )
        if (not stat.S_ISDIR(task_metadata.st_mode) or
                (task_metadata.st_dev, task_metadata.st_ino) !=
                (current.st_dev, current.st_ino)):
            return None
        stat_descriptor = os.open(
            "stat", os.O_RDONLY | os.O_NOFOLLOW |
            getattr(os, "O_CLOEXEC", 0),
            dir_fd=task_descriptor,
        )
        chunks = []
        total = 0
        while total <= 65536:
            block = os.read(stat_descriptor, 65537 - total)
            if not block:
                break
            chunks.append(block)
            total += len(block)
        if total > 65536:
            return None
        value = b"".join(chunks).decode("ascii")
        record = parse_proc_stat(value)
        if record is None or record["pid"] != pid:
            return None
        record["proc_identity"] = (
            task_metadata.st_dev, task_metadata.st_ino
        )
        record["proc_descriptor"] = task_descriptor
        transferred_record = record
        return record
    except (OSError, UnicodeError):
        return None
    except BaseException:
        transferred_record = None
        raise
    finally:
        for descriptor in (stat_descriptor, proc_descriptor):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
        if (task_descriptor is not None and
                transferred_record is None):
            try:
                os.close(task_descriptor)
            except OSError:
                pass


def proc_record(pid):
    record = open_proc_record(pid)
    if record is None:
        return None
    descriptor = record.pop("proc_descriptor")
    os.close(descriptor)
    return record


def process_has_descriptor_identity(pid, identity):
    if sys.platform != "linux":
        return False
    directory = Path("/proc") / str(pid) / "fd"
    try:
        entries = list(directory.iterdir())
    except OSError:
        return False
    for entry in entries:
        try:
            metadata = entry.stat()
        except OSError:
            continue
        if (metadata.st_dev, metadata.st_ino) == identity:
            return True
    return False


def direct_child_identities():
    if sys.platform != "linux":
        return set()
    result = set()
    try:
        entries = list(Path("/proc").iterdir())
    except OSError:
        return result
    for entry in entries:
        if not entry.name.isdigit():
            continue
        record = proc_record(int(entry.name))
        if record is not None and record["ppid"] == os.getpid():
            result.add((
                record["pid"], record["starttime_ticks"],
                record.get("proc_identity"),
            ))
    return result


def contained_processes(
        session, containment_identity, adopted_baseline=None):
    """Inventory this launch by session, ancestry, and its opaque pipe."""
    if sys.platform != "linux":
        return []
    records = {}
    try:
        entries = list(Path("/proc").iterdir())
    except OSError:
        return []
    for entry in entries:
        if not entry.name.isdigit():
            continue
        record = proc_record(int(entry.name))
        if record is None or record["pid"] == os.getpid():
            continue
        records[record["pid"]] = record

    def descends_from_leader(record):
        seen = set()
        parent = record["ppid"]
        while parent > 0 and parent not in seen:
            if parent == session:
                return True
            seen.add(parent)
            ancestor = records.get(parent)
            if ancestor is None:
                return False
            parent = ancestor["ppid"]
        return False

    result = []
    for record in records.values():
        identity = (
            record["pid"], record["starttime_ticks"],
            record.get("proc_identity"),
        )
        adopted = (
            adopted_baseline is not None and
            record["ppid"] == os.getpid() and
            identity not in adopted_baseline
        )
        if (record["session"] == session or
                descends_from_leader(record) or adopted or
                process_has_descriptor_identity(
                    record["pid"], containment_identity
                )):
            result.append(record)
    return result


def linux_pidfd_open(pid):
    wrapper = getattr(os, "pidfd_open", None)
    if callable(wrapper):
        return wrapper(pid, 0)
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise CrossoverError(
            "Linux pidfd support is unavailable: {}".format(error)
        )
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


def linux_pidfd_send_signal(descriptor, signal_number):
    wrapper = getattr(signal, "pidfd_send_signal", None)
    if callable(wrapper):
        wrapper(descriptor, signal_number, None, 0)
        return
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise CrossoverError(
            "Linux pidfd signaling is unavailable: {}".format(error)
        )
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


def validate_pidfd_runtime(description):
    """Prove pidfd open and signal operations work before any child starts."""
    descriptor = None
    try:
        descriptor = linux_pidfd_open(os.getpid())
        linux_pidfd_send_signal(descriptor, 0)
    except (CrossoverError, OSError) as error:
        raise CrossoverError(
            "{} requires working Linux pidfd operations: {}".format(
                description, error
            )
        ) from error
    finally:
        if descriptor is not None:
            os.close(descriptor)


def await_launch_gate(descriptor, description):
    """Block a freshly execed helper until its pidfd identity is retained."""
    if type(descriptor) is not int or descriptor < 3:
        raise CrossoverError("{} launch gate is invalid".format(description))
    try:
        value = os.read(descriptor, 2)
    except OSError as error:
        raise CrossoverError(
            "cannot read {} launch gate: {}".format(description, error)
        ) from error
    finally:
        os.close(descriptor)
    if value != b"G":
        raise CrossoverError(
            "{} launch gate was not released".format(description)
        )


def retain_process_identity(record, description):
    """Open a pidfd only while one retained /proc task inode stays current."""
    if sys.platform != "linux":
        raise CrossoverError(
            "{} requires Linux pidfd support".format(description)
        )
    before = open_proc_record(record["pid"])
    if (before is None or before["starttime_ticks"] !=
            record["starttime_ticks"] or
            (record.get("proc_identity") is not None and
             before["proc_identity"] != record["proc_identity"])):
        if before is not None:
            os.close(before["proc_descriptor"])
        return None
    descriptor = None
    after = None
    transferred = None
    try:
        descriptor = linux_pidfd_open(record["pid"])
        after = open_proc_record(record["pid"])
        if (after is None or
                after["starttime_ticks"] != before["starttime_ticks"] or
                after["starttime_ticks"] != record["starttime_ticks"] or
                after["proc_identity"] != before["proc_identity"]):
            os.close(descriptor)
            descriptor = None
            if after is not None:
                os.close(after["proc_descriptor"])
                after = None
            return None
        retained = dict(after)
        retained["pidfd"] = descriptor
        transferred = retained
        return retained
    except ProcessLookupError:
        if descriptor is not None:
            os.close(descriptor)
            descriptor = None
        return None
    except OSError as error:
        if descriptor is not None:
            os.close(descriptor)
            descriptor = None
        raise CrossoverError(
            "cannot retain {} PID {}: {}".format(
                description, record["pid"], error
            )
        )
    except BaseException:
        transferred = None
        raise
    finally:
        os.close(before["proc_descriptor"])
        if descriptor is not None and transferred is None:
            os.close(descriptor)
        if after is not None and transferred is None:
            os.close(after["proc_descriptor"])


def retain_launched_process(process, description):
    record = proc_record(process.pid)
    if record is None:
        raise CrossoverError(
            "{} exited before its identity could be retained".format(
                description
            )
        )
    retained = retain_process_identity(record, description)
    if retained is None:
        raise CrossoverError(
            "{} identity changed while opening its pidfd".format(description)
        )
    return retained


def close_retained_process(record):
    for key in ("pidfd", "proc_descriptor"):
        descriptor = record.get(key)
        if descriptor is not None:
            try:
                os.close(descriptor)
            finally:
                record[key] = None


def retained_process_alive(record):
    descriptor = record.get("pidfd")
    if descriptor is None:
        return False
    try:
        poller = select.poll()
        poller.register(descriptor, select.POLLIN | select.POLLHUP |
                       select.POLLERR)
        return not poller.poll(0)
    except OSError:
        return False


def wait_retained_process_exit(record, timeout):
    """Wait for exact pidfd readiness without reaping/releasing the PID."""
    if (not isinstance(timeout, (int, float)) or isinstance(timeout, bool) or
            not math.isfinite(float(timeout)) or timeout < 0):
        raise CrossoverError("retained-process timeout is invalid")
    descriptor = record.get("pidfd")
    if descriptor is None:
        return True
    deadline = time.monotonic() + float(timeout)
    poller = select.poll()
    poller.register(descriptor, select.POLLIN | select.POLLHUP |
                    select.POLLERR)
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return not retained_process_alive(record)
        try:
            if poller.poll(max(1, int(math.ceil(remaining * 1000.0)))):
                return True
        except InterruptedError:
            continue


def signal_retained_process(record, signal_number):
    descriptor = record.get("pidfd")
    if descriptor is None:
        return
    try:
        linux_pidfd_send_signal(descriptor, signal_number)
    except (ProcessLookupError, PermissionError):
        pass


def reap_retained_process(record):
    """Reap an adopted exact child through its retained pidfd."""
    descriptor = record.get("pidfd")
    if descriptor is None or not hasattr(os, "waitid"):
        return
    try:
        os.waitid(
            getattr(os, "P_PIDFD", 3), descriptor,
            os.WEXITED | os.WNOHANG)
    except (ChildProcessError, ProcessLookupError):
        pass


def retain_discovered_processes(
        session, containment_identity, adopted_baseline, retained):
    """Retain every discovered identity before it can be acted upon."""
    for snapshot in contained_processes(
            session, containment_identity, adopted_baseline):
        key = (
            snapshot["pid"], snapshot["starttime_ticks"],
            snapshot.get("proc_identity"),
        )
        if key in retained:
            continue
        record = retain_process_identity(snapshot, "command descendant")
        if record is not None:
            retained[key] = record


def open_containment_descriptor():
    flags = getattr(os, "O_CLOEXEC", 0)
    try:
        if hasattr(os, "pipe2"):
            read_descriptor, write_descriptor = os.pipe2(flags)
        else:
            read_descriptor, write_descriptor = os.pipe()
            os.set_inheritable(read_descriptor, False)
            os.set_inheritable(write_descriptor, False)
        metadata = os.fstat(read_descriptor)
        if not stat.S_ISFIFO(metadata.st_mode):
            raise CrossoverError(
                "command containment descriptor is not a pipe"
            )
        return (
            read_descriptor, write_descriptor,
            (metadata.st_dev, metadata.st_ino)
        )
    except BaseException:
        for descriptor in (
                locals().get("read_descriptor"),
                locals().get("write_descriptor")):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
        raise


def open_owned_directory(path, description):
    path = checked_path(path, description)
    if not path.is_absolute():
        raise CrossoverError("{} is not absolute".format(description))
    directory = _open_absolute_directory_components(path, description, False)
    transferred = False
    try:
        transferred = True
        return directory["descriptor"]
    except BaseException:
        transferred = False
        raise
    finally:
        if not transferred:
            close_owned_directory(directory)


def open_existing_canonical_directory(path, description):
    descriptor = None
    transferred = False
    try:
        descriptor = open_owned_directory(
            Path(os.path.abspath(os.fspath(path))), description
        )
        metadata = os.fstat(descriptor)
        result = {
            "descriptor": descriptor,
            "identity": (metadata.st_dev, metadata.st_ino),
            "path": Path(os.path.abspath(os.fspath(path))),
            "required_mode": 0o700,
        }
        transferred = True
        return result
    except BaseException:
        transferred = False
        raise
    finally:
        if descriptor is not None and not transferred:
            os.close(descriptor)


def open_exclusive_log(path, description):
    path = checked_path(path, description)
    if path.name in ("", ".", ".."):
        raise CrossoverError("{} has an unsafe name".format(description))
    directory_descriptor = None
    descriptor = None
    identity = None
    try:
        directory_descriptor = open_owned_directory(
            path.parent, description + " directory"
        )
        flags = (
            os.O_RDWR | os.O_CREAT | os.O_EXCL |
            os.O_CLOEXEC | os.O_NOFOLLOW
        )
        descriptor = os.open(
            path.name, flags, 0o600, dir_fd=directory_descriptor
        )
        os.fchmod(descriptor, 0o600)
        metadata = os.fstat(descriptor)
        current = os.stat(
            path.name, dir_fd=directory_descriptor, follow_symlinks=False
        )
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != 0o600 or
                (metadata.st_dev, metadata.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "{} is not a stable owner-only regular file".format(
                    description
                )
            )
        identity = (metadata.st_dev, metadata.st_ino)
        result = {
            "descriptor": descriptor,
            "directory_identity": (
                os.fstat(directory_descriptor).st_dev,
                os.fstat(directory_descriptor).st_ino,
            ),
            "directory_descriptor": directory_descriptor,
            "directory_path": path.parent,
            "identity": identity,
            "name": path.name,
            "path": path,
        }
        return result
    except FileExistsError as error:
        primary = CrossoverError(
            "{} already exists; refusing to replace or follow it".format(
                description
            )
        )
        primary.__cause__ = error
    except BaseException as error:
        primary = error
    cleanup_errors = []
    if descriptor is not None:
        try:
            descriptor_metadata = os.fstat(descriptor)
            current = os.stat(
                path.name, dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
            descriptor_identity = (
                descriptor_metadata.st_dev, descriptor_metadata.st_ino
            )
            if (descriptor_identity !=
                    (current.st_dev, current.st_ino) or
                    (identity is not None and
                     descriptor_identity != identity)):
                raise CrossoverError(
                    "{} rollback name changed".format(description)
                )
            os.unlink(path.name, dir_fd=directory_descriptor)
            os.fsync(directory_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
        try:
            os.close(descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if directory_descriptor is not None:
        try:
            os.close(directory_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if cleanup_errors:
        raise CrossoverError(
            "{} failed and cleanup also failed: {}; cleanup: {}".format(
                description, primary,
                "; ".join(
                    "{}: {}".format(type(error).__name__, error)
                    for error in cleanup_errors
                ),
            )
        ) from primary
    raise primary


def open_exclusive_owned_file(directory, name, description):
    """Create one owner-only file relative to an already retained directory."""
    if (not isinstance(name, str) or name in ("", ".", "..") or
            "/" in name):
        raise CrossoverError("{} has an unsafe name".format(description))
    directory_descriptor = None
    descriptor = None
    identity = None
    try:
        validate_owned_directory(directory, description + " directory")
        directory_descriptor = os.dup(directory["descriptor"])
        flags = (
            os.O_RDWR | os.O_CREAT | os.O_EXCL |
            os.O_CLOEXEC | os.O_NOFOLLOW
        )
        descriptor = os.open(
            name, flags, 0o600, dir_fd=directory_descriptor
        )
        os.fchmod(descriptor, 0o600)
        metadata = os.fstat(descriptor)
        current = os.stat(
            name, dir_fd=directory_descriptor, follow_symlinks=False
        )
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != 0o600 or
                metadata.st_size != 0 or
                (metadata.st_dev, metadata.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "{} is not a new stable owner-only regular file".format(
                    description
                )
            )
        validate_owned_directory(directory, description + " directory")
        identity = (metadata.st_dev, metadata.st_ino)
        result = {
            "descriptor": descriptor,
            "directory_identity": directory["identity"],
            "directory_descriptor": directory_descriptor,
            "directory_path": directory["path"],
            "identity": identity,
            "name": name,
            "path": directory["path"] / name,
        }
        return result
    except FileExistsError as error:
        primary = CrossoverError(
            "{} already exists; refusing to replace or follow it".format(
                description
            )
        )
        primary.__cause__ = error
    except BaseException as error:
        primary = error
    cleanup_errors = []
    if descriptor is not None:
        try:
            retained = os.fstat(descriptor)
            current = os.stat(
                name, dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
            retained_identity = (retained.st_dev, retained.st_ino)
            if (retained_identity !=
                    (current.st_dev, current.st_ino) or
                    (identity is not None and
                     retained_identity != identity)):
                raise CrossoverError(
                    "{} rollback name changed".format(description)
                )
            os.unlink(name, dir_fd=directory_descriptor)
            os.fsync(directory_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
        try:
            os.close(descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if directory_descriptor is not None:
        try:
            os.close(directory_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if cleanup_errors:
        raise CrossoverError(
            "{} failed and cleanup also failed: {}; cleanup: {}".format(
                description, primary,
                "; ".join(
                    "{}: {}".format(type(error).__name__, error)
                    for error in cleanup_errors
                ),
            )
        ) from primary
    raise primary


def close_log(log):
    errors = []
    for key in ("descriptor", "directory_descriptor"):
        descriptor = log.get(key)
        if descriptor is not None:
            try:
                os.close(descriptor)
            except BaseException as error:
                errors.append(error)
            finally:
                log[key] = None
    if errors:
        raise CrossoverError(
            "cannot close retained file descriptors: {}".format(
                "; ".join(
                    "{}: {}".format(type(error).__name__, error)
                    for error in errors
                )
            )
        )


def validate_log_identity(log, description):
    directory_current = None
    try:
        metadata = os.fstat(log["descriptor"])
        directory_metadata = os.fstat(log["directory_descriptor"])
        directory_required_mode = log.get(
            "directory_required_mode", 0o700
        )
        directory_current = _open_absolute_directory_components(
            log["directory_path"], description + " directory", False,
            directory_required_mode,
        )
        current = os.stat(
            log["name"], dir_fd=log["directory_descriptor"],
            follow_symlinks=False
        )
        directory_current_metadata = os.fstat(
            directory_current["descriptor"]
        )
        required_mode = log.get("required_mode", 0o600)
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != required_mode or
                (directory_metadata.st_dev, directory_metadata.st_ino) !=
                log["directory_identity"] or
                (directory_current_metadata.st_dev,
                 directory_current_metadata.st_ino) !=
                log["directory_identity"] or
                stat.S_IMODE(directory_metadata.st_mode) !=
                directory_required_mode or
                (metadata.st_dev, metadata.st_ino) != log["identity"] or
                (current.st_dev, current.st_ino) != log["identity"]):
            raise CrossoverError(
                "{} changed during command execution".format(description)
            )
    finally:
        if directory_current is not None:
            close_owned_directory(directory_current)


def open_result_regular_held(
        root, relative, maximum_bytes, description, required_mode,
        directory_required_mode=0o700):
    """Retain one exact regular-file inode below an already held result root."""
    parts = result_relative_parts(relative, description)
    if len(parts) == 1:
        directory = duplicate_owned_directory(
            root, description + " directory"
        )
    else:
        directory = open_existing_result_directory(
            root, Path(*parts[:-1]), description + " directory",
            required_mode=directory_required_mode,
        )
    descriptor = None
    transferred = False
    try:
        descriptor = os.open(
            parts[-1],
            os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NONBLOCK", 0),
            dir_fd=directory["descriptor"],
        )
        metadata = os.fstat(descriptor)
        current = os.stat(
            parts[-1], dir_fd=directory["descriptor"],
            follow_symlinks=False,
        )
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != required_mode or
                metadata.st_size > maximum_bytes or
                (metadata.st_dev, metadata.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "{} is not the expected bounded owner-controlled file".format(
                    description
                )
            )
        held = {
            "descriptor": descriptor,
            "directory_identity": directory["identity"],
            "directory_descriptor": directory["descriptor"],
            "directory_path": directory["path"],
            "directory_required_mode":
                directory.get("required_mode", 0o700),
            "identity": (metadata.st_dev, metadata.st_ino),
            "name": parts[-1],
            "path": directory["path"] / parts[-1],
            "required_mode": required_mode,
        }
        validate_log_identity(held, description)
        validate_owned_directory(root, "canonical result directory")
        transferred = True
        return held
    except BaseException:
        transferred = False
        raise
    finally:
        if descriptor is not None and not transferred:
            os.close(descriptor)
        if not transferred:
            close_owned_directory(directory)


def unlink_new_log(log):
    errors = []
    try:
        metadata = os.fstat(log["descriptor"])
        current = os.stat(
            log["name"], dir_fd=log["directory_descriptor"],
            follow_symlinks=False,
        )
        if ((metadata.st_dev, metadata.st_ino) != log["identity"] or
                (current.st_dev, current.st_ino) != log["identity"]):
            raise CrossoverError(
                "new command log changed before rollback"
            )
        os.unlink(log["name"], dir_fd=log["directory_descriptor"])
        os.fsync(log["directory_descriptor"])
    except BaseException as error:
        errors.append(error)
    try:
        close_log(log)
    except BaseException as error:
        errors.append(error)
    if errors:
        raise CrossoverError(
            "cannot roll back new command log: {}".format(
                "; ".join(
                    "{}: {}".format(type(error).__name__, error)
                    for error in errors
                )
            )
        )


def open_command_logs(stdout_path, stderr_path):
    stdout_log = None
    stderr_log = None
    try:
        stdout_log = open_exclusive_log(stdout_path, "stdout log")
        stderr_log = open_exclusive_log(stderr_path, "stderr log")
        return stdout_log, stderr_log
    except BaseException as primary:
        cleanup_errors = []
        for label, log in (
                ("stderr", stderr_log), ("stdout", stdout_log)):
            if log is None:
                continue
            try:
                unlink_new_log(log)
            except BaseException as cleanup:
                cleanup_errors.append((label, cleanup))
        if cleanup_errors:
            raise CrossoverError(
                "command-log reservation failed: {}; rollback failed: "
                "{}".format(
                    primary,
                    "; ".join(
                        "{} {}: {}".format(
                            label, type(cleanup).__name__, cleanup
                        ) for label, cleanup in cleanup_errors
                    ),
                )
            ) from primary
        raise


def append_log(log, value):
    validate_log_identity(log, "command log")
    os.lseek(log["descriptor"], 0, os.SEEK_END)
    remaining = value
    while remaining:
        written = os.write(log["descriptor"], remaining)
        if written <= 0:
            raise CrossoverError("cannot append command log")
        remaining = remaining[written:]


def read_normalized_log(log, description):
    validate_log_identity(log, description)
    os.lseek(log["descriptor"], 0, os.SEEK_SET)
    chunks = []
    while True:
        chunk = os.read(log["descriptor"], 1024 * 1024)
        if not chunk:
            break
        chunks.append(chunk)
    raw = b"".join(chunks)
    value = normalized_output(raw)
    if value != raw:
        os.lseek(log["descriptor"], 0, os.SEEK_SET)
        os.ftruncate(log["descriptor"], 0)
        remaining = value
        while remaining:
            written = os.write(log["descriptor"], remaining)
            if written <= 0:
                raise CrossoverError("cannot normalize {}".format(description))
            remaining = remaining[written:]
    os.fsync(log["descriptor"])
    validate_log_identity(log, description)
    return value


def duplicate_lock_descriptor(descriptor):
    if descriptor is None:
        return None
    if type(descriptor) is not int or descriptor < 0:
        raise CrossoverError("inherited lock descriptor is invalid")
    try:
        metadata = os.fstat(descriptor)
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1):
            raise CrossoverError(
                "inherited lock descriptor is not an owned regular file"
            )
        return os.dup(descriptor)
    except OSError as error:
        raise CrossoverError(
            "cannot duplicate inherited lock descriptor: {}".format(error)
        )


def write_descriptor_all(descriptor, value):
    offset = 0
    while offset < len(value):
        try:
            written = os.write(descriptor, value[offset:])
        except InterruptedError:
            continue
        if written <= 0:
            raise CrossoverError(
                "direct command supervisor write made no progress"
            )
        offset += written


class BoundedControlDrain(object):
    """Drain a control pipe concurrently while retaining only a bounded prefix."""

    def __init__(self, descriptor, maximum_bytes, description):
        if (type(descriptor) is not int or descriptor < 0 or
                type(maximum_bytes) is not int or maximum_bytes <= 0):
            raise CrossoverError(
                "{} drain contract is invalid".format(description)
            )
        self._descriptor = descriptor
        self._maximum_bytes = maximum_bytes
        self._description = description
        self._payload = bytearray()
        self._oversize = False
        self._error = None
        self._stop = threading.Event()
        try:
            os.set_blocking(self._descriptor, False)
        except (AttributeError, OSError) as error:
            os.close(self._descriptor)
            self._descriptor = None
            raise CrossoverError(
                "cannot make {} nonblocking: {}".format(description, error)
            )
        self._thread = threading.Thread(
            target=self._drain,
            name="leopard2-bounded-control-drain",
        )
        self._thread.daemon = False
        self._thread.start()

    def _drain(self):
        try:
            poller = select.poll()
            poller.register(
                self._descriptor,
                select.POLLIN | select.POLLHUP | select.POLLERR,
            )
            while True:
                if self._stop.is_set():
                    break
                try:
                    events = poller.poll(100)
                except InterruptedError:
                    continue
                if not events:
                    continue
                while True:
                    try:
                        block = os.read(self._descriptor, 65536)
                    except BlockingIOError:
                        break
                    except InterruptedError:
                        continue
                    if not block:
                        return
                    room = self._maximum_bytes + 1 - len(self._payload)
                    if room > 0:
                        self._payload.extend(block[:room])
                    if (len(block) > room or
                            len(self._payload) > self._maximum_bytes):
                        self._oversize = True
                    # Continue draining after the cap so the writer cannot block.
        except BaseException as error:
            self._error = error
        finally:
            try:
                os.close(self._descriptor)
            except OSError:
                pass
            self._descriptor = None

    def finish(self, timeout=2.0):
        if (not isinstance(timeout, (int, float)) or
                isinstance(timeout, bool) or
                not math.isfinite(float(timeout)) or timeout <= 0):
            raise CrossoverError(
                "{} drain timeout is invalid".format(self._description)
            )
        self._thread.join(float(timeout))
        if self._thread.is_alive():
            self._stop.set()
            self._thread.join(0.5)
            if self._thread.is_alive():
                raise CrossoverError(
                    "{} drain did not stop after its deadline".format(
                        self._description
                    )
                )
            raise CrossoverError(
                "{} writer remained open after its deadline".format(
                    self._description
                )
            )
        if self._error is not None:
            raise CrossoverError(
                "cannot drain {}: {}".format(
                    self._description, self._error
                )
            )
        if self._oversize or len(self._payload) > self._maximum_bytes:
            raise CrossoverError(
                "{} exceeds {} bytes".format(
                    self._description, self._maximum_bytes
                )
            )
        return bytes(self._payload)


def remap_inherited_descriptors(inherited, remaps, forbidden):
    """Apply an arbitrary FD remap graph without clobbering any source.

    Every source is snapshotted above the complete source/target/internal-FD
    namespace before the first ``dup2``.  An inherited descriptor that is
    overwritten but is not itself a remap source is retained at a backup FD;
    this preserves containment and lock identities even under collisions.
    """
    inherited_set = set(inherited)
    targets = set(remaps)
    sources = set(remaps.values())
    occupied = inherited_set | targets | sources | set(forbidden)
    minimum = max(occupied | {2}) + 1
    source_copies = {}
    retained_backups = {}
    transferred = False
    try:
        for source in sorted(sources):
            copied = fcntl.fcntl(
                source, fcntl.F_DUPFD_CLOEXEC, minimum
            )
            source_copies[source] = copied
            occupied.add(copied)
            minimum = max(minimum, copied + 1)
        for target in sorted(targets):
            if target in inherited_set and target not in sources:
                copied = fcntl.fcntl(
                    target, fcntl.F_DUPFD_CLOEXEC, minimum
                )
                retained_backups[target] = copied
                occupied.add(copied)
                minimum = max(minimum, copied + 1)
        for target, source in sorted(remaps.items()):
            os.dup2(source_copies[source], target, inheritable=True)
        result = (
            (inherited_set - sources - targets) |
            set(targets) | set(retained_backups.values())
        )
        transferred = True
        return tuple(sorted(result))
    except OSError as error:
        transferred = False
        raise CrossoverError(
            "cannot apply collision-safe descriptor remaps: {}".format(
                error
            )
        )
    except BaseException:
        transferred = False
        raise
    finally:
        for descriptor in source_copies.values():
            try:
                os.close(descriptor)
            except OSError:
                pass
        # Successful backups belong to the returned descriptor set.  On
        # failure, close them here because no child will inherit them.
        if not transferred:
            for descriptor in retained_backups.values():
                try:
                    os.close(descriptor)
                except OSError:
                    pass


def relocate_control_descriptor(control_descriptor, inherited, remaps):
    """Move a supervisor control FD away from arbitrary benchmark targets."""
    if control_descriptor not in remaps:
        return control_descriptor
    relocation_floor = max(
        set(inherited) | set(remaps) | set(remaps.values()) |
        {control_descriptor}
    ) + 1
    try:
        relocated = fcntl.fcntl(
            control_descriptor, fcntl.F_DUPFD_CLOEXEC,
            relocation_floor,
        )
        os.close(control_descriptor)
        return relocated
    except OSError as error:
        raise CrossoverError(
            "cannot relocate colliding supervisor control descriptor: "
            "{}".format(error)
        )


def direct_command_supervisor(arguments):
    """Single-threaded boundary for run_abba's audited Linux subreaper."""
    if len(arguments) != 9:
        raise CrossoverError(
            "direct command supervisor argument count changed"
        )
    try:
        launch_gate_descriptor = int(arguments[8])
    except (TypeError, ValueError, OverflowError) as error:
        raise CrossoverError(
            "direct command supervisor launch gate is invalid"
        ) from error
    await_launch_gate(
        launch_gate_descriptor, "direct command supervisor")
    argv = decode_json_bytes(
        arguments[0].encode("utf-8"),
        "direct command supervisor argv"
    )
    environment = decode_json_bytes(
        arguments[2].encode("utf-8"),
        "direct command supervisor environment"
    )
    inherited = decode_json_bytes(
        arguments[4].encode("utf-8"),
        "direct command supervisor descriptors"
    )
    remaps = decode_json_bytes(
        arguments[5].encode("utf-8"),
        "direct command supervisor descriptor remaps"
    )
    launcher_value = decode_json_bytes(
        arguments[7].encode("utf-8"),
        "direct command supervisor launcher"
    )
    launcher = runtime_launcher_from_contract(launcher_value)
    try:
        timeout = float(arguments[3])
        control_descriptor = int(arguments[6])
    except (TypeError, ValueError, OverflowError) as error:
        raise CrossoverError(
            "direct command supervisor scalar argument is invalid: {}".format(
                error
            )
        )
    if (not isinstance(argv, list) or not argv or
            not all(isinstance(item, str) and item and "\0" not in item
                    for item in argv)):
        raise CrossoverError(
            "direct command supervisor argv is invalid"
        )
    if (not isinstance(environment, dict) or
            not all(isinstance(key, str) and isinstance(value, str)
                    for key, value in environment.items())):
        raise CrossoverError(
            "direct command supervisor environment is invalid"
        )
    if (not isinstance(inherited, list) or
            not all(type(descriptor) is int and descriptor >= 0
                    for descriptor in inherited) or
            not isinstance(remaps, dict) or
            type(control_descriptor) is not int or control_descriptor < 0 or
            not math.isfinite(timeout) or timeout <= 0 or timeout > 3600):
        raise CrossoverError(
            "direct command supervisor execution contract is invalid"
        )
    inherited = tuple(sorted(set(inherited)))
    parsed_remaps = {}
    for target_text, source in remaps.items():
        try:
            target = int(target_text)
        except (TypeError, ValueError) as error:
            raise CrossoverError(
                "direct command supervisor remap target is invalid"
            ) from error
        if (str(target) != target_text or target < 3 or
                type(source) is not int or source not in inherited):
            raise CrossoverError(
                "direct command supervisor descriptor remap is unsafe"
            )
        parsed_remaps[target] = source
    launcher_descriptors = tuple(
        launcher[name]["descriptor"] for name in RUNTIME_LAUNCHER_NAMES
    )
    if (set(launcher_descriptors) & (
            set(inherited) | set(remaps) | set(remaps.values()) |
            {control_descriptor})):
        raise CrossoverError(
            "runtime launcher descriptors collide with command transport"
        )
    for descriptor in inherited + launcher_descriptors + (
            control_descriptor,):
        try:
            os.fstat(descriptor)
        except OSError as error:
            raise CrossoverError(
                "direct command supervisor descriptor {} is invalid: "
                "{}".format(descriptor, error)
            )
    control_descriptor = relocate_control_descriptor(
        control_descriptor, inherited, parsed_remaps
    )
    inherited = remap_inherited_descriptors(
        inherited, parsed_remaps, (control_descriptor,)
    )
    # Benchmark-created JSON must remain owner-only regardless of the
    # coordinator's inherited shell umask.  The kernel additionally enforces
    # the evidence byte ceiling on the exact inherited output inode.
    os.umask(0o077)
    try:
        resource.setrlimit(
            resource.RLIMIT_FSIZE,
            (MAX_RAW_JSON_BYTES, MAX_RAW_JSON_BYTES),
        )
    except (OSError, ValueError) as error:
        raise CrossoverError(
            "cannot enforce benchmark output ceiling: {}".format(error)
        )

    load_sealed_python_module(
        launcher["build_provenance"], "leopard2_build_provenance",
        "sealed build-provenance module",
    )
    load_sealed_python_module(
        launcher["git_capture"], "git_capture",
        "sealed Git-capture module",
    )
    load_sealed_python_module(
        launcher["link_common"], "balanced_evidence_common",
        "sealed link-evidence module",
    )
    containment_module = load_sealed_python_module(
        launcher["containment"], "run_abba",
        "sealed command-containment module",
    )
    try:
        ContainmentError = containment_module.EvidenceError
        audited_run_process_bounded = (
            containment_module.run_process_bounded
        )
    except AttributeError as error:
        raise CrossoverError(
            "sealed command-containment module omits its audited API"
        ) from error

    status = None
    try:
        completed = audited_run_process_bounded(
            argv, cwd=Path(arguments[1]), environment=environment,
            timeout=timeout, inherited_descriptors=inherited
        )
        write_descriptor_all(1, completed.stdout)
        write_descriptor_all(2, completed.stderr)
        status = {
            "returncode": completed.returncode,
            "status": "ok",
        }
    except ContainmentError as error:
        message = str(error)
        status = {
            "message": message,
            "status": (
                "timeout" if message.startswith("command exceeded ")
                else "error"
            ),
        }
    except BaseException as error:
        status = {
            "message": "{}: {}".format(type(error).__name__, error),
            "status": "error",
        }
    encoded = canonical_bytes(status)
    if not encoded or len(encoded) > MAX_SUPERVISOR_CONTROL_BYTES:
        raise CrossoverError(
            "direct command supervisor control record is invalid"
        )
    try:
        write_descriptor_all(control_descriptor, encoded)
    finally:
        os.close(control_descriptor)
    return 0


def decode_direct_supervisor_control(payload):
    if not payload or len(payload) > MAX_SUPERVISOR_CONTROL_BYTES:
        raise CrossoverError(
            "direct command supervisor returned no valid control record"
        )
    control = decode_json_bytes(
        payload, "direct command supervisor control"
    )
    if not isinstance(control, dict):
        raise CrossoverError(
            "direct command supervisor control is not an object"
        )
    return control


def remove_stale_owned_file(directory, name, description):
    """Remove only an exact owner-controlled regular file for safe resume."""
    if not isinstance(name, str) or name in ("", ".", "..") or "/" in name:
        raise CrossoverError("{} has an unsafe name".format(description))
    validate_owned_directory(directory, description + " directory")
    try:
        metadata = os.stat(
            name, dir_fd=directory["descriptor"], follow_symlinks=False
        )
    except FileNotFoundError:
        return
    except OSError as error:
        raise CrossoverError(
            "cannot inspect stale {}: {}".format(description, error)
        )
    if (not stat.S_ISREG(metadata.st_mode) or
            metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
            stat.S_IMODE(metadata.st_mode) != 0o600):
        raise CrossoverError(
            "refusing to remove unsafe stale {}".format(description)
        )
    os.unlink(name, dir_fd=directory["descriptor"])
    os.fsync(directory["descriptor"])
    validate_owned_directory(directory, description + " directory")


def read_owned_regular(
        directory, name, maximum_bytes, description,
        mutation_hook=None, required_mode=0o600):
    """Read a bounded exact file from a held, revalidated directory."""
    if (not isinstance(name, str) or name in ("", ".", "..") or
            "/" in name or type(maximum_bytes) is not int or
            maximum_bytes <= 0 or required_mode not in (0o600, 0o755)):
        raise CrossoverError("{} read contract is invalid".format(description))
    validate_owned_directory(directory, description + " directory")
    flags = (
        os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NONBLOCK", 0)
    )
    descriptor = None
    try:
        descriptor = os.open(
            name, flags, dir_fd=directory["descriptor"]
        )
        metadata = os.fstat(descriptor)
        current = os.stat(
            name, dir_fd=directory["descriptor"], follow_symlinks=False
        )
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != required_mode or
                metadata.st_size > maximum_bytes or
                (metadata.st_dev, metadata.st_ino) !=
                (current.st_dev, current.st_ino)):
            raise CrossoverError(
                "{} is not a bounded owner-controlled regular file".format(
                    description
                )
            )
        initial_identity = (
            metadata.st_dev, metadata.st_ino, metadata.st_size,
            metadata.st_mtime_ns, metadata.st_ctime_ns,
        )
        chunks = []
        total = 0
        hook_called = False
        while total <= maximum_bytes:
            block = os.read(
                descriptor, min(1024 * 1024, maximum_bytes + 1 - total)
            )
            if not block:
                break
            chunks.append(block)
            total += len(block)
            if mutation_hook is not None and not hook_called:
                hook_called = True
                mutation_hook()
        if total > maximum_bytes:
            raise CrossoverError(
                "{} exceeds {} bytes".format(description, maximum_bytes)
            )
        os.lseek(descriptor, 0, os.SEEK_SET)
        verification = bytearray()
        while len(verification) <= maximum_bytes:
            block = os.read(
                descriptor,
                min(1024 * 1024,
                    maximum_bytes + 1 - len(verification))
            )
            if not block:
                break
            verification.extend(block)
        if (len(verification) > maximum_bytes or
                bytes(verification) != b"".join(chunks)):
            raise CrossoverError(
                "{} changed while it was being read".format(description)
            )
        final_metadata = os.fstat(descriptor)
        final_current = os.stat(
            name, dir_fd=directory["descriptor"], follow_symlinks=False
        )
        final_identity = (
            final_metadata.st_dev, final_metadata.st_ino,
            final_metadata.st_size, final_metadata.st_mtime_ns,
            final_metadata.st_ctime_ns,
        )
        if (total != metadata.st_size or
                final_identity != initial_identity or
                (final_current.st_dev, final_current.st_ino) !=
                (metadata.st_dev, metadata.st_ino) or
                final_current.st_size != metadata.st_size or
                final_current.st_mtime_ns != metadata.st_mtime_ns or
                final_current.st_ctime_ns != metadata.st_ctime_ns):
            raise CrossoverError(
                "{} changed while it was being read".format(description)
            )
        validate_owned_directory(directory, description + " directory")
        return b"".join(chunks)
    except OSError as error:
        raise CrossoverError(
            "cannot read {}: {}".format(description, error)
        )
    finally:
        if descriptor is not None:
            os.close(descriptor)


def read_held_regular(
        held, maximum_bytes, description, mutation_hook=None):
    """Read/hash/parse bytes through the exact pre-opened regular-file inode."""
    if (type(maximum_bytes) is not int or maximum_bytes <= 0 or
            not isinstance(held, dict) or held.get("descriptor") is None):
        raise CrossoverError("{} held-read contract is invalid".format(
            description
        ))
    validate_log_identity(held, description)
    descriptor = held["descriptor"]
    metadata = os.fstat(descriptor)
    if metadata.st_size > maximum_bytes:
        raise CrossoverError(
            "{} exceeds {} bytes".format(description, maximum_bytes)
        )
    initial = (
        metadata.st_dev, metadata.st_ino, metadata.st_size,
        metadata.st_mtime_ns, metadata.st_ctime_ns,
    )
    values = []
    for pass_index in range(2):
        os.lseek(descriptor, 0, os.SEEK_SET)
        chunks = []
        total = 0
        hook_called = False
        while total <= maximum_bytes:
            block = os.read(
                descriptor,
                min(1024 * 1024, maximum_bytes + 1 - total),
            )
            if not block:
                break
            chunks.append(block)
            total += len(block)
            if (mutation_hook is not None and pass_index == 0 and
                    not hook_called):
                hook_called = True
                mutation_hook()
        if total > maximum_bytes:
            raise CrossoverError(
                "{} exceeds {} bytes".format(description, maximum_bytes)
            )
        values.append(b"".join(chunks))
    final = os.fstat(descriptor)
    final_identity = (
        final.st_dev, final.st_ino, final.st_size,
        final.st_mtime_ns, final.st_ctime_ns,
    )
    validate_log_identity(held, description)
    if (values[0] != values[1] or initial != final_identity or
            len(values[0]) != metadata.st_size):
        raise CrossoverError(
            "{} changed while it was being read".format(description)
        )
    return values[0]


def read_path_owned_regular(path, maximum_bytes, description):
    path = checked_path(path, description)
    try:
        resolved_parent = path.parent.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise CrossoverError(
            "cannot resolve {} directory: {}".format(description, error)
        )
    if resolved_parent != path.parent:
        raise CrossoverError(
            "{} directory contains a symbolic-link component".format(
                description
            )
        )
    directory = owned_canonical_directory(
        path.parent, description + " directory"
    )
    try:
        return read_owned_regular(
            directory, path.name, maximum_bytes, description
        )
    finally:
        close_owned_directory(directory)


def procfd_exact_path(held, description,
                      child_descriptor=RAW_OUTPUT_DESCRIPTOR):
    validate_log_identity(held, description)
    if type(child_descriptor) is not int or child_descriptor < 3:
        raise CrossoverError(
            "{} child descriptor is invalid".format(description)
        )
    return Path("/proc/self/fd") / str(child_descriptor)


def cleanup_command_tree(
        process, containment_identity, timed_out,
        adopted_baseline=None, leader_identity=None):
    """Terminate and reap every observable member of one command launch."""
    if leader_identity is None:
        raise CrossoverError(
            "command cleanup requires a retained launch pidfd"
        )
    retained = {
        (
            leader_identity["pid"], leader_identity["starttime_ticks"],
            leader_identity.get("proc_identity"),
        ):
            leader_identity
    }
    try:
        retain_discovered_processes(
            process.pid, containment_identity, adopted_baseline, retained
        )
        residual = sorted(
            record["pid"] for record in retained.values()
            if record is not leader_identity
        )
        if not timed_out and not residual:
            process.wait(timeout=0)
            return residual

        # Never signal a numeric PID, process group, or session here.  A pidfd
        # names the exact task retained between two matching /proc starttime
        # reads and remains safe after the numeric identifier can be reused.
        for record in retained.values():
            signal_retained_process(record, signal.SIGTERM)

        term_deadline = time.monotonic() + 0.35
        quiescent_scans = 0
        while time.monotonic() < term_deadline:
            count_before = len(retained)
            retain_discovered_processes(
                process.pid, containment_identity, adopted_baseline, retained
            )
            for record in retained.values():
                if record is not leader_identity:
                    reap_retained_process(record)
            live = [
                record for record in retained.values()
                if retained_process_alive(record)
            ]
            if not live:
                quiescent_scans = (
                    quiescent_scans + 1
                    if len(retained) == count_before else 0
                )
                if quiescent_scans >= 3:
                    process.wait(timeout=0)
                    return residual
            else:
                quiescent_scans = 0
            time.sleep(0.01)

        for record in retained.values():
            signal_retained_process(record, signal.SIGKILL)

        kill_deadline = time.monotonic() + 2.0
        quiescent_scans = 0
        while time.monotonic() < kill_deadline:
            count_before = len(retained)
            retain_discovered_processes(
                process.pid, containment_identity, adopted_baseline, retained
            )
            for record in retained.values():
                signal_retained_process(record, signal.SIGKILL)
                if record is not leader_identity:
                    reap_retained_process(record)
            remaining = [
                record for record in retained.values()
                if retained_process_alive(record)
            ]
            if not remaining:
                quiescent_scans = (
                    quiescent_scans + 1
                    if len(retained) == count_before else 0
                )
                if quiescent_scans >= 3:
                    process.wait(timeout=0)
                    return residual
            else:
                quiescent_scans = 0
            time.sleep(0.01)
        raise CrossoverError(
            "command descendant cleanup failed for retained PIDs {}".format(
                ",".join(str(record["pid"]) for record in remaining)
            )
        )
    finally:
        for record in retained.values():
            if record is not leader_identity:
                close_retained_process(record)


def _run_command_owned(
        argv, cwd, stdout_path, stderr_path, timeout, environment,
        inherited_lock_descriptor=None, inherited_descriptors=(),
        descriptor_remaps=None, runtime_launcher=None):
    validate_pidfd_runtime("direct command supervisor")
    ensure_child_subreaper()
    adopted_baseline = direct_child_identities()
    if adopted_baseline:
        raise CrossoverError(
            "isolated command owner has pre-existing child processes"
        )
    if (not isinstance(inherited_descriptors, (tuple, list)) or
            not all(type(descriptor) is int and descriptor >= 0
                    for descriptor in inherited_descriptors)):
        raise CrossoverError(
            "command inherited descriptor set is invalid"
        )
    for descriptor in inherited_descriptors:
        try:
            os.fstat(descriptor)
        except OSError as error:
            raise CrossoverError(
                "command inherited descriptor {} is invalid: {}".format(
                    descriptor, error
                )
            )
    if descriptor_remaps is None:
        descriptor_remaps = {}
    if (not isinstance(descriptor_remaps, dict) or
            not all(type(target) is int and target >= 3 and
                    type(source) is int and source in inherited_descriptors
                    for target, source in descriptor_remaps.items())):
        raise CrossoverError(
            "command inherited descriptor remap is invalid"
        )
    launcher_identity = validate_runtime_launcher_snapshots(
        runtime_launcher
    )
    launcher_contract = runtime_launcher_contract(runtime_launcher)
    launcher_descriptors = tuple(
        runtime_launcher[name]["descriptor"]
        for name in RUNTIME_LAUNCHER_NAMES
    )
    if (set(launcher_descriptors) & (
            set(inherited_descriptors) | set(descriptor_remaps) |
            set(descriptor_remaps.values()))):
        raise CrossoverError(
            "runtime launcher descriptors collide with command transport"
        )
    stdout_log, stderr_log = open_command_logs(
        stdout_path, stderr_path
    )
    containment_descriptor = None
    containment_writer = None
    lock_descriptor = None
    process = None
    leader_identity = None
    control_read = None
    control_write = None
    control_drain = None
    launch_gate_read = None
    launch_gate_write = None
    timed_out = False
    residual = []
    try:
        (containment_descriptor, containment_writer,
         containment_identity) = open_containment_descriptor()
        lock_descriptor = duplicate_lock_descriptor(
            inherited_lock_descriptor
        )
        control_read, control_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        launch_gate_read, launch_gate_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        benchmark_descriptors = list(inherited_descriptors)
        benchmark_descriptors.append(containment_descriptor)
        if lock_descriptor is not None:
            benchmark_descriptors.append(lock_descriptor)
        benchmark_descriptors = sorted(set(benchmark_descriptors))
        pass_descriptors = (
            benchmark_descriptors + list(launcher_descriptors) +
            [control_write, launch_gate_read]
        )
        process = subprocess.Popen(
            [
                runtime_launcher["identity"]["python"]["path"],
                *RUNTIME_PYTHON_FLAGS,
                str(Path("/proc/self/fd") / str(
                    runtime_launcher["script"]["descriptor"]
                )),
                DIRECT_COMMAND_SUPERVISOR_MODE,
                json.dumps(
                    [str(item) for item in argv],
                    ensure_ascii=True, separators=(",", ":")
                ),
                str(cwd),
                json.dumps(
                    dict(environment), ensure_ascii=True,
                    sort_keys=True, separators=(",", ":")
                ),
                repr(float(timeout)),
                json.dumps(benchmark_descriptors, separators=(",", ":")),
                json.dumps(
                    {str(target): source for target, source
                     in sorted(descriptor_remaps.items())},
                    sort_keys=True, separators=(",", ":")
                ),
                str(control_write),
                json.dumps(
                    launcher_contract, ensure_ascii=True,
                    sort_keys=True, separators=(",", ":")
                ),
                str(launch_gate_read),
            ],
            cwd=str(cwd), env=environment,
            stdin=subprocess.DEVNULL,
            stdout=stdout_log["descriptor"],
            stderr=stderr_log["descriptor"],
            start_new_session=True,
            pass_fds=tuple(pass_descriptors),
            executable=str(Path("/proc/self/fd") / str(
                runtime_launcher["python"]["descriptor"]
            )),
        )
        os.close(launch_gate_read)
        launch_gate_read = None
        leader_identity = retain_launched_process(
            process, "direct command supervisor"
        )
        write_descriptor_all(launch_gate_write, b"G")
        os.close(launch_gate_write)
        launch_gate_write = None
        os.close(control_write)
        control_write = None
        control_drain = BoundedControlDrain(
            control_read, MAX_SUPERVISOR_CONTROL_BYTES,
            "direct command supervisor control",
        )
        control_read = None
        timed_out = not wait_retained_process_exit(
            leader_identity,
            float(timeout) + SUPERVISOR_REAP_GRACE_SECONDS,
        )
        residual = cleanup_command_tree(
            process, containment_identity, timed_out, adopted_baseline,
            leader_identity,
        )
        if timed_out:
            raise CrossoverError(
                "direct command supervisor exceeded its containment bound"
            )
        control = decode_direct_supervisor_control(control_drain.finish())
        control_drain = None
        if process.returncode != 0:
            raise CrossoverError(
                "direct command supervisor exited with status {}".format(
                    process.returncode
                )
            )
        status = control.get("status")
        if status == "timeout":
            if (set(control) != {"message", "status"} or
                    not isinstance(control.get("message"), str)):
                raise CrossoverError(
                    "direct command supervisor timeout record is malformed"
                )
            timed_out = True
            returncode = 124
        elif status == "error":
            if (set(control) != {"message", "status"} or
                    not isinstance(control.get("message"), str)):
                raise CrossoverError(
                    "direct command supervisor error record is malformed"
                )
            returncode = 125
            append_log(
                stderr_log,
                (
                    "leopard2 runner: command containment failed: " +
                    control["message"] + "\n"
                ).encode("utf-8", errors="replace")
            )
        elif (status == "ok" and
                set(control) == {"returncode", "status"} and
                type(control.get("returncode")) is int):
            returncode = control["returncode"]
        else:
            raise CrossoverError(
                "direct command supervisor success record is malformed"
            )
        if residual:
            returncode = 125
            append_log(
                stderr_log,
                b"leopard2 runner: residual command descendants were "
                b"terminated\n"
            )
        stdout = read_normalized_log(stdout_log, "stdout log")
        stderr = read_normalized_log(stderr_log, "stderr log")
        return {
            "argv": [str(item) for item in argv],
            "cwd": str(cwd),
            "launcher_identity": launcher_identity,
            "returncode": returncode,
            "stderr_log": checked_path(stderr_path, "stderr log").name,
            "stderr_sha256": digest_bytes(stderr),
            "stdout_log": checked_path(stdout_path, "stdout log").name,
            "stdout_sha256": digest_bytes(stdout),
            "timed_out": timed_out,
        }
    except BaseException:
        if process is not None and process.poll() is None:
            try:
                if leader_identity is not None:
                    cleanup_command_tree(
                        process, containment_identity, True,
                        adopted_baseline, leader_identity,
                    )
                else:
                    # The child is still an unreaped launch owned by this
                    # Popen object, so its numeric PID cannot yet be reused.
                    process.kill()
                    process.wait(timeout=1.0)
            except Exception:
                pass
        raise
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors = []
        if control_drain is not None:
            try:
                control_drain.finish()
            except BaseException as error:
                cleanup_errors.append(("control drain", error))
        if leader_identity is not None:
            try:
                close_retained_process(leader_identity)
            except BaseException as error:
                cleanup_errors.append(("retained supervisor", error))
        # The lock duplicate is deliberately closed only after process-tree
        # cleanup.  Never issue LOCK_UN: descendants inherit the same
        # open-file description and must keep it locked if this coordinator
        # is killed.
        for descriptor in (
                lock_descriptor, containment_descriptor,
                containment_writer, control_read, control_write,
                launch_gate_read, launch_gate_write):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except BaseException as error:
                    cleanup_errors.append(("transport descriptor", error))
        for label, log in (
                ("stderr log", stderr_log), ("stdout log", stdout_log)):
            try:
                close_log(log)
            except BaseException as error:
                cleanup_errors.append((label, error))
        if cleanup_errors:
            detail = "; ".join(
                "{} {}: {}".format(
                    label, type(error).__name__, error
                ) for label, error in cleanup_errors
            )
            if active_error is not None:
                raise CrossoverError(
                    "owned command failed: {}; cleanup also failed: {}".format(
                        active_error, detail
                    )
                ) from active_error
            raise CrossoverError(
                "owned command cleanup failed: {}".format(detail)
            )


def direct_command_owner(arguments):
    """Run exactly one command inside an isolated subreaper process."""
    if len(arguments) != 12:
        raise CrossoverError(
            "direct command owner argument count changed"
        )
    try:
        launch_gate_descriptor = int(arguments[11])
    except (TypeError, ValueError, OverflowError) as error:
        raise CrossoverError(
            "direct command owner launch gate is invalid"
        ) from error
    await_launch_gate(launch_gate_descriptor, "direct command owner")
    argv = decode_json_bytes(
        arguments[0].encode("utf-8"), "direct command owner argv"
    )
    environment = decode_json_bytes(
        arguments[5].encode("utf-8"),
        "direct command owner environment"
    )
    inherited = decode_json_bytes(
        arguments[7].encode("utf-8"),
        "direct command owner descriptors"
    )
    remaps_value = decode_json_bytes(
        arguments[8].encode("utf-8"),
        "direct command owner descriptor remaps"
    )
    launcher_value = decode_json_bytes(
        arguments[10].encode("utf-8"),
        "direct command owner launcher"
    )
    launcher = runtime_launcher_from_contract(launcher_value)
    try:
        timeout = float(arguments[4])
        lock_descriptor_value = int(arguments[6])
        control_descriptor = int(arguments[9])
    except (TypeError, ValueError, OverflowError) as error:
        raise CrossoverError(
            "direct command owner scalar argument is invalid: {}".format(
                error
            )
        )
    if (not isinstance(argv, list) or not argv or
            not all(isinstance(item, str) and item and "\0" not in item
                    for item in argv) or
            not isinstance(environment, dict) or
            not all(isinstance(key, str) and isinstance(value, str)
                    for key, value in environment.items()) or
            not isinstance(inherited, list) or
            not all(type(descriptor) is int and descriptor >= 0
                    for descriptor in inherited) or
            not isinstance(remaps_value, dict) or
            not math.isfinite(timeout) or timeout <= 0 or timeout > 3600 or
            lock_descriptor_value < -1 or
            type(control_descriptor) is not int or control_descriptor < 0):
        raise CrossoverError(
            "direct command owner execution contract is invalid"
        )
    remaps = {}
    for target_text, source in remaps_value.items():
        try:
            target = int(target_text)
        except (TypeError, ValueError) as error:
            raise CrossoverError(
                "direct command owner remap target is invalid"
            ) from error
        if (str(target) != target_text or target < 3 or
                type(source) is not int or source not in inherited):
            raise CrossoverError(
                "direct command owner remap is invalid"
            )
        remaps[target] = source
    lock_descriptor = (
        None if lock_descriptor_value < 0 else lock_descriptor_value
    )
    status = None
    try:
        record = _run_command_owned(
            argv, Path(arguments[1]), Path(arguments[2]),
            Path(arguments[3]), timeout, environment,
            inherited_lock_descriptor=lock_descriptor,
            inherited_descriptors=tuple(inherited),
            descriptor_remaps=remaps,
            runtime_launcher=launcher,
        )
        status = {"record": record, "status": "ok"}
    except BaseException as error:
        status = {
            "message": "{}: {}".format(type(error).__name__, error),
            "status": "error",
        }
    encoded = canonical_bytes(status)
    if not encoded or len(encoded) > MAX_SUPERVISOR_CONTROL_BYTES:
        raise CrossoverError(
            "direct command owner control record is invalid"
        )
    try:
        write_descriptor_all(control_descriptor, encoded)
    finally:
        os.close(control_descriptor)
    return 0


def run_command(
        argv, cwd, stdout_path, stderr_path, timeout, environment,
        inherited_lock_descriptor=None, inherited_descriptors=(),
        descriptor_remaps=None):
    """Execute through a per-command owner so concurrent jobs never mix."""
    validate_pidfd_runtime("direct command owner")
    if descriptor_remaps is None:
        descriptor_remaps = {}
    if (not isinstance(argv, (tuple, list)) or not argv or
            not all(isinstance(item, (str, Path)) and
                    str(item) and "\0" not in str(item)
                    for item in argv) or
            not isinstance(environment, dict) or
            not all(isinstance(key, str) and isinstance(value, str)
                    for key, value in environment.items()) or
            not isinstance(timeout, (int, float)) or
            isinstance(timeout, bool) or not math.isfinite(float(timeout)) or
            timeout <= 0 or timeout > 3600 or
            not isinstance(inherited_descriptors, (tuple, list)) or
            not all(type(descriptor) is int and descriptor >= 0
                    for descriptor in inherited_descriptors) or
            not isinstance(descriptor_remaps, dict) or
            not all(type(target) is int and target >= 3 and
                    type(source) is int and source in inherited_descriptors
                    for target, source in descriptor_remaps.items())):
        raise CrossoverError(
            "command owner execution contract is invalid"
        )
    for descriptor in inherited_descriptors:
        try:
            os.fstat(descriptor)
        except OSError as error:
            raise CrossoverError(
                "command owner inherited descriptor {} is invalid: "
                "{}".format(descriptor, error)
            )
    pass_descriptors = list(inherited_descriptors)
    if inherited_lock_descriptor is not None:
        if (type(inherited_lock_descriptor) is not int or
                inherited_lock_descriptor < 0):
            raise CrossoverError(
                "command owner lock descriptor is invalid"
            )
        pass_descriptors.append(inherited_lock_descriptor)
    containment_descriptor = None
    containment_writer = None
    control_read = None
    control_write = None
    process = None
    leader_identity = None
    control_drain = None
    launch_gate_read = None
    launch_gate_write = None
    launcher = None
    try:
        launcher = open_runtime_launcher_snapshots()
        launcher_identity = validate_runtime_launcher_snapshots(launcher)
        launcher_contract = runtime_launcher_contract(launcher)
        launcher_descriptors = tuple(
            launcher[name]["descriptor"]
            for name in RUNTIME_LAUNCHER_NAMES
        )
        if (set(launcher_descriptors) & (
                set(inherited_descriptors) | set(descriptor_remaps) |
                set(descriptor_remaps.values()) |
                ({inherited_lock_descriptor}
                 if inherited_lock_descriptor is not None else set()))):
            raise CrossoverError(
                "runtime launcher descriptors collide with command transport"
            )
        (containment_descriptor, containment_writer,
         containment_identity) = open_containment_descriptor()
        forwarded = sorted(set(
            list(inherited_descriptors) + [containment_descriptor]
        ))
        control_read, control_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        launch_gate_read, launch_gate_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        pass_descriptors.extend(
            launcher_descriptors +
            (containment_descriptor, control_write, launch_gate_read)
        )
        process = subprocess.Popen(
            [
                launcher["identity"]["python"]["path"],
                *RUNTIME_PYTHON_FLAGS,
                str(Path("/proc/self/fd") / str(
                    launcher["script"]["descriptor"]
                )),
                DIRECT_COMMAND_OWNER_MODE,
                json.dumps(
                    [str(item) for item in argv],
                    ensure_ascii=True, separators=(",", ":")
                ),
                str(cwd), str(stdout_path), str(stderr_path),
                repr(float(timeout)),
                json.dumps(
                    dict(environment), ensure_ascii=True,
                    sort_keys=True, separators=(",", ":")
                ),
                str(
                    -1 if inherited_lock_descriptor is None
                    else inherited_lock_descriptor
                ),
                json.dumps(forwarded, separators=(",", ":")),
                json.dumps(
                    {str(target): source for target, source
                     in sorted(descriptor_remaps.items())},
                    sort_keys=True, separators=(",", ":")
                ),
                str(control_write),
                json.dumps(
                    launcher_contract, ensure_ascii=True,
                    sort_keys=True, separators=(",", ":")
                ),
                str(launch_gate_read),
            ],
            cwd=str(cwd), env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True,
            pass_fds=tuple(sorted(set(pass_descriptors))),
            executable=str(Path("/proc/self/fd") / str(
                launcher["python"]["descriptor"]
            )),
        )
        os.close(launch_gate_read)
        launch_gate_read = None
        leader_identity = retain_launched_process(
            process, "direct command owner"
        )
        write_descriptor_all(launch_gate_write, b"G")
        os.close(launch_gate_write)
        launch_gate_write = None
        os.close(control_write)
        control_write = None
        control_drain = BoundedControlDrain(
            control_read, MAX_SUPERVISOR_CONTROL_BYTES,
            "direct command owner control",
        )
        control_read = None
        timed_out = False
        timed_out = not wait_retained_process_exit(
            leader_identity,
            float(timeout) + 2.0 * SUPERVISOR_REAP_GRACE_SECONDS,
        )
        cleanup_command_tree(
            process, containment_identity, timed_out,
            leader_identity=leader_identity,
        )
        if timed_out:
            raise CrossoverError(
                "direct command owner exceeded its containment bound"
            )
        control = decode_direct_supervisor_control(control_drain.finish())
        control_drain = None
        if process.returncode != 0:
            raise CrossoverError(
                "direct command owner exited with status {}".format(
                    process.returncode
                )
            )
        if (control.get("status") == "error" and
                set(control) == {"message", "status"} and
                isinstance(control.get("message"), str)):
            raise CrossoverError(
                "direct command owner rejected execution: " +
                control["message"]
            )
        if (control.get("status") != "ok" or
                set(control) != {"record", "status"} or
                not isinstance(control.get("record"), dict)):
            raise CrossoverError(
                "direct command owner control is malformed"
            )
        record = control["record"]
        if canonical_bytes(record.get("launcher_identity")) != canonical_bytes(
                launcher_identity):
            raise CrossoverError(
                "direct command owner launcher identity changed"
            )
        return record
    except BaseException:
        if process is not None and process.poll() is None:
            try:
                if leader_identity is not None:
                    cleanup_command_tree(
                        process, containment_identity, True,
                        leader_identity=leader_identity,
                    )
                else:
                    process.kill()
                    process.wait(timeout=1.0)
            except Exception:
                pass
        raise
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors = []
        if control_drain is not None:
            try:
                control_drain.finish()
            except BaseException as error:
                cleanup_errors.append(("owner control drain", error))
        if leader_identity is not None:
            try:
                close_retained_process(leader_identity)
            except BaseException as error:
                cleanup_errors.append(("retained owner", error))
        for descriptor in (
                containment_descriptor, containment_writer,
                control_read, control_write,
                launch_gate_read, launch_gate_write):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except BaseException as error:
                    cleanup_errors.append(("owner transport descriptor", error))
        if launcher is not None:
            try:
                close_runtime_launcher_snapshots(launcher)
            except BaseException as error:
                cleanup_errors.append(("runtime launcher", error))
        if cleanup_errors:
            detail = "; ".join(
                "{} {}: {}".format(
                    label, type(error).__name__, error
                ) for label, error in cleanup_errors
            )
            if active_error is not None:
                raise CrossoverError(
                    "command owner failed: {}; cleanup also failed: {}".format(
                        active_error, detail
                    )
                ) from active_error
            raise CrossoverError(
                "command owner cleanup failed: {}".format(detail)
            )


def controlled_avx2_configure_argv(cmake, ninja, source, build_root):
    """Return the stable historical mode-0 configure command.

    The default small-direct mode deliberately remains implicit.  Adding an
    explicit ``-D...=0`` would make otherwise equivalent historical campaign
    commands differ even though production compilation is unchanged.
    """
    return [
        str(cmake),
        "-S", str(source),
        "-B", str(build_root),
        "-G", "Ninja",
        "-DCMAKE_MAKE_PROGRAM={}".format(ninja),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DLEO2_BACKEND_VARIANT=avx2",
        "-DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git",
        "-DLEO2_BUILD_BENCHMARKS=ON",
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
        "-DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF",
        "-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF",
        "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF",
        "-DLEOPARD_ENABLE_GF8=ON",
        "-DLEOPARD_ENABLE_GF16=ON",
    ]


def controlled_avx2_build(
        source, result_dir, source_state, parallel, validate_guard,
        guard_identity, inherited_lock_descriptor, result_root):
    if type(parallel) is not int or parallel <= 0 or parallel > 128:
        raise CrossoverError("controlled build parallelism is outside [1,128]")
    if not callable(validate_guard):
        raise CrossoverError("controlled build requires a live guard validator")
    if (not isinstance(guard_identity, dict) or set(guard_identity) != {
            "canonical_lock", "pair_lease"}):
        raise CrossoverError(
            "controlled build requires exact campaign guard identities"
        )
    validate_guard()
    validate_owned_directory(result_root, "canonical result directory")
    build_root = result_dir / "controlled-build-avx2"
    record_path = result_dir / "controlled-build.json"
    try:
        build_metadata_entry = os.stat(
            "controlled-build-avx2", dir_fd=result_root["descriptor"],
            follow_symlinks=False,
        )
    except FileNotFoundError:
        build_metadata_entry = None
    try:
        record_metadata_entry = os.stat(
            "controlled-build.json", dir_fd=result_root["descriptor"],
            follow_symlinks=False,
        )
    except FileNotFoundError:
        record_metadata_entry = None
    if build_metadata_entry is not None or record_metadata_entry is not None:
        if (build_metadata_entry is None or record_metadata_entry is None or
                not stat.S_ISDIR(build_metadata_entry.st_mode) or
                not stat.S_ISREG(record_metadata_entry.st_mode)):
            raise CrossoverError(
                "controlled build is only partially retained; use a new "
                "--result-dir"
            )
        try:
            record_bytes = read_result_regular(
                result_root, "controlled-build.json",
                MAX_RAW_JSON_BYTES, "controlled-build record",
            )
            record = decode_json_bytes(record_bytes, str(record_path))
            executable = Path(record["executable"]).resolve()
            try:
                executable_relative = executable.relative_to(result_dir)
            except ValueError:
                raise CrossoverError(
                    "controlled-build executable escapes its result directory"
                )
            executable_bytes = read_result_regular(
                result_root, executable_relative,
                MAX_RAW_JSON_BYTES, "controlled-build executable",
                required_mode=0o755,
            )
            metadata = cmake_build_metadata(executable)
        except (CrossoverError, KeyError, OSError, TypeError) as error:
            raise CrossoverError(
                "cannot resume controlled AVX2 build: {}".format(error)
            )
        if (not isinstance(record, dict) or set(record) != {
                "backend", "build_metadata", "build_root", "build_tools",
                "commands",
                "executable", "executable_sha256", "guard_identity",
                "schema", "source_identity"} or
                record.get("schema") != CONTROLLED_BUILD_SCHEMA or
                record.get("backend") != "avx2" or
                record.get("build_root") != str(build_root.resolve()) or
                canonical_bytes(record.get("source_identity")) !=
                canonical_bytes(source_state) or
                canonical_bytes(record.get("guard_identity")) !=
                canonical_bytes(guard_identity) or
                canonical_bytes(record.get("build_metadata")) !=
                canonical_bytes(metadata) or
                record.get("executable_sha256") !=
                digest_bytes(executable_bytes)):
            raise CrossoverError(
                "retained controlled AVX2 build differs from this campaign"
            )
        build_tools = record.get("build_tools")
        if (not isinstance(build_tools, dict) or
                set(build_tools) != {"cmake", "ninja"}):
            raise CrossoverError(
                "retained controlled AVX2 build tool identity is invalid"
            )
        validate_current_executable_identity(
            build_tools["cmake"], "retained controlled CMake"
        )
        validate_current_executable_identity(
            build_tools["ninja"], "retained controlled Ninja"
        )
        for command_index, command in enumerate(record["commands"]):
            for stream in ("stdout", "stderr"):
                relative = command.get(stream + "_log")
                value = read_result_regular(
                    result_root, relative, MAX_RETAINED_LOG_BYTES,
                    "controlled-build {} {} log".format(
                        command_index, stream
                    ),
                )
                if digest_bytes(value) != command.get(stream + "_sha256"):
                    raise CrossoverError(
                        "retained controlled-build log hash changed"
                    )
        validate_build_source_binding(
            metadata, source, source_state, "avx2", require_fresh=True
        )
        validate_guard()
        return executable, {
            "path": str(record_path.relative_to(result_dir)),
            "schema": CONTROLLED_BUILD_SCHEMA,
            "sha256": digest_bytes(record_bytes),
        }
    cmake = shutil.which("cmake", path=BENCHMARK_ENVIRONMENT["PATH"])
    ninja = shutil.which("ninja", path=BENCHMARK_ENVIRONMENT["PATH"])
    if not cmake or not ninja:
        raise CrossoverError(
            "controlled authoritative build requires /usr/bin cmake and ninja"
        )
    # These same-UID launchers are snapshotted because they orchestrate every
    # build command.  The compiler/linker paths and their effective flags are
    # captured by the CMake/object/link attestation below; the campaign threat
    # model treats their root-owned installation as the OS/toolchain trust
    # boundary rather than recursively snapshotting the dynamic loader,
    # shared libraries, assembler, and linker dependency closure.
    log_dir = result_dir / "build-logs"
    build_log_directory = open_result_directory(
        result_root, "build-logs", "controlled-build log directory", True
    )
    close_owned_directory(build_log_directory)
    cmake_procfd = str(
        Path("/proc/self/fd") / str(CONTROLLED_CMAKE_DESCRIPTOR)
    )
    ninja_procfd = str(
        Path("/proc/self/fd") / str(CONTROLLED_NINJA_DESCRIPTOR)
    )
    configure = controlled_avx2_configure_argv(
        cmake_procfd, ninja_procfd, source, build_root
    )
    build = [
        cmake_procfd, "--build", str(build_root),
        "--target", "bench_leopard2_direct_encode",
        "--parallel", str(parallel),
    ]
    commands = []
    cmake_snapshot = None
    ninja_snapshot = None
    try:
        cmake_snapshot = open_exact_executable_snapshot(
            Path(cmake).resolve(), "controlled CMake executable", True
        )
        ninja_snapshot = open_exact_executable_snapshot(
            Path(ninja).resolve(), "controlled Ninja executable", True
        )
        build_tools = {
            "cmake": sealed_snapshot_identity(cmake_snapshot),
            "ninja": sealed_snapshot_identity(ninja_snapshot),
        }
    except BaseException as primary:
        cleanup_errors = []
        for snapshot in (ninja_snapshot, cmake_snapshot):
            if snapshot is not None:
                try:
                    close_sealed_snapshot(snapshot)
                except BaseException as error:
                    cleanup_errors.append(error)
        if cleanup_errors:
            raise CrossoverError(
                "controlled build tool snapshot failed: {}; cleanup also "
                "failed: {}".format(
                    primary,
                    "; ".join(
                        "{}: {}".format(type(error).__name__, error)
                        for error in cleanup_errors
                    ),
                )
            ) from primary
        raise
    build_primary = None
    try:
        for label, argv in (("configure", configure), ("build", build)):
            validate_guard()
            stdout_path = log_dir / (label + ".stdout.log")
            stderr_path = log_dir / (label + ".stderr.log")
            bound_stdout_path = (
                Path("/proc/self/fd") / str(result_root["descriptor"]) /
                stdout_path.relative_to(result_dir)
            )
            bound_stderr_path = (
                Path("/proc/self/fd") / str(result_root["descriptor"]) /
                stderr_path.relative_to(result_dir)
            )
            command = run_command(
                argv, source, bound_stdout_path, bound_stderr_path, 1800,
                dict(GIT_ENVIRONMENT),
                inherited_lock_descriptor=inherited_lock_descriptor,
                inherited_descriptors=(
                    result_root["descriptor"],
                    cmake_snapshot["descriptor"],
                    ninja_snapshot["descriptor"],
                ),
                descriptor_remaps={
                    CONTROLLED_CMAKE_DESCRIPTOR:
                        cmake_snapshot["descriptor"],
                    CONTROLLED_NINJA_DESCRIPTOR:
                        ninja_snapshot["descriptor"],
                },
            )
            command["environment"] = dict(GIT_ENVIRONMENT)
            command["stdout_log"] = str(stdout_path.relative_to(result_dir))
            command["stderr_log"] = str(stderr_path.relative_to(result_dir))
            commands.append(command)
            validate_guard()
            validate_sealed_executable_snapshot(
                cmake_snapshot, build_tools["cmake"]["sha256"],
                "executed controlled CMake snapshot",
            )
            validate_sealed_executable_snapshot(
                ninja_snapshot, build_tools["ninja"]["sha256"],
                "executed controlled Ninja snapshot",
            )
            if command["returncode"] != 0 or command["timed_out"]:
                raise CrossoverError(
                    "controlled AVX2 {} failed; see {} and {}".format(
                        label, stdout_path, stderr_path
                    )
                )
    except BaseException as error:
        build_primary = error
    close_errors = []
    for label, snapshot in (
            ("Ninja", ninja_snapshot), ("CMake", cmake_snapshot)):
        try:
            close_sealed_snapshot(snapshot)
        except BaseException as error:
            close_errors.append((label, error))
    if build_primary is not None:
        if close_errors:
            raise CrossoverError(
                "controlled build failed: {}; tool teardown failed: {}".format(
                    build_primary,
                    "; ".join(
                        "{} {}: {}".format(
                            label, type(error).__name__, error
                        ) for label, error in close_errors
                    ),
                )
            ) from build_primary
        raise build_primary
    if close_errors:
        raise CrossoverError(
            "controlled build tool teardown failed: {}".format(
                "; ".join(
                    "{} {}: {}".format(
                        label, type(error).__name__, error
                    ) for label, error in close_errors
                )
            )
        )
    executable = (build_root / "bench_leopard2_direct_encode").resolve()
    if not executable.is_file() or not os.access(str(executable), os.X_OK):
        raise CrossoverError(
            "controlled AVX2 build omitted {}".format(executable)
        )
    metadata = cmake_build_metadata(executable)
    validate_build_source_binding(
        metadata, source, source_state, "avx2", require_fresh=True
    )
    validate_guard()
    executable_bytes = read_result_regular(
        result_root, executable.relative_to(result_dir),
        MAX_RAW_JSON_BYTES, "controlled-build executable",
        required_mode=0o755,
    )
    record = {
        "backend": "avx2",
        "build_metadata": metadata,
        "build_root": str(build_root.resolve()),
        "build_tools": build_tools,
        "commands": commands,
        "executable": str(executable),
        "executable_sha256": digest_bytes(executable_bytes),
        "guard_identity": guard_identity,
        "schema": CONTROLLED_BUILD_SCHEMA,
        "source_identity": source_state,
    }
    validate_owned_directory(result_root, "canonical result directory")
    write_result_json(
        result_root, "controlled-build.json", record,
        "controlled-build record", replace=False,
    )
    validate_owned_directory(result_root, "canonical result directory")
    try:
        validate_guard()
    except Exception:
        try:
            record_path.unlink()
        except FileNotFoundError:
            pass
        raise
    record_bytes = read_result_regular(
        result_root, "controlled-build.json", MAX_RAW_JSON_BYTES,
        "controlled-build record",
    )
    return executable, {
        "path": str(record_path.relative_to(result_dir)),
        "schema": CONTROLLED_BUILD_SCHEMA,
        "sha256": digest_bytes(record_bytes),
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


def required_mapping(value, path):
    if not isinstance(value, dict):
        raise CrossoverError("benchmark {} must be a JSON object".format(path))
    return value


def require_exact_keys(value, expected, path):
    mapping = required_mapping(value, path)
    actual = set(mapping)
    expected = set(expected)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise CrossoverError(
            "benchmark {} has invalid v2 fields (missing={}, extra={})".format(
                path, missing, extra
            )
        )
    return mapping


def required_finite_number(value, path):
    number = None
    if not isinstance(value, bool) and isinstance(value, (int, float)):
        try:
            number = float(value)
        except (OverflowError, TypeError, ValueError):
            number = None
    if number is None or not math.isfinite(number):
        raise CrossoverError(
            "benchmark {} must be a finite JSON number".format(path)
        )
    return number


def required_finite_metric(value, path, allow_zero):
    number = required_finite_number(value, path)
    if (number < 0) if allow_zero else (number <= 0):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise CrossoverError(
            "benchmark {} must be a finite {} JSON number".format(path, qualifier)
        )
    return number


def validate_raw(raw, job, timed_mode, settings):
    raw = require_exact_keys(raw, {
        "build", "correctness", "memory", "methodology", "metrics",
        "operation_model", "parameters", "resolved", "schema",
    }, "document")
    if raw.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    build = require_exact_keys(raw.get("build"), {
        "backend_variant", "build_configuration_sha256", "build_type",
        "compiler", "compiler_version", "cplusplus", "source_commit",
        "source_tracked_dirty", "source_tree", "test_hooks",
    }, "build")
    parameters = require_exact_keys(raw.get("parameters"), {
        "K", "Q", "R", "batch", "forced_mode", "iterations",
        "requested_field", "requested_parity_indices", "requested_profile",
        "reuse", "seed", "shard_bytes", "thread_count", "warmups",
    }, "parameters")
    resolved = require_exact_keys(raw.get("resolved"), {
        "backend", "direct_capable", "field", "padded_side", "parent_count",
        "profile", "thread_count", "timed_path_is_direct",
    }, "resolved")
    correctness = require_exact_keys(raw.get("correctness"), {
        "direct_transform_parity_match", "parity_checksum_fnv1a64",
        "selected_transform_reference_parity_match",
        "unrequested_outputs_untouched",
    }, "correctness")
    operation_model = require_exact_keys(raw.get("operation_model"), {
        "direct_model_applies_to_timed_path",
        "direct_output_accumulations", "direct_output_initializations",
        "direct_row_terms",
        "fixed_coefficient_symbol_terms_before_unit_specialization",
        "hardware_counters", "model_scope", "modeled_output_bytes_read",
        "modeled_output_bytes_written", "modeled_source_bytes_read",
        "transform_operation_counts", "xor_accumulation_symbols",
    }, "operation_model")
    required_mapping(raw.get("memory"), "memory")
    required_mapping(raw.get("methodology"), "methodology")
    metrics = required_mapping(raw.get("metrics"), "metrics")
    execution = require_exact_keys(metrics.get("encode_execution"), {
        "logical_input_GB_per_s", "logical_input_plus_output_GB_per_s",
        "mad_us_per_batch_call", "maximum_us_per_batch_call",
        "median_us_per_batch_call", "minimum_us_per_batch_call",
        "requested_parity_output_GB_per_s",
    }, "metrics.encode_execution")
    item = job["cell"]
    try:
        git_identity = job["source_identity"]["git"]
        expected_build_identity = {
            "source_commit": git_identity["head"],
            "source_tree": git_identity["tree"],
            "source_tracked_dirty": not git_identity["worktree_clean"],
        }
    except (KeyError, TypeError):
        raise CrossoverError("job omits its exact Git source identity")
    for key, value in expected_build_identity.items():
        if canonical_bytes(build.get(key)) != canonical_bytes(value):
            raise CrossoverError(
                "benchmark build {} is {!r}, expected {!r}".format(
                    key, build.get(key), value
                )
            )
    if build.get("test_hooks") is not True:
        raise CrossoverError("benchmark build omits required test hooks")
    if build.get("backend_variant") != item["backend"]:
        raise CrossoverError(
            "benchmark executable variant {!r} differs from cell {!r}".format(
                build.get("backend_variant"), item["backend"]
            )
        )
    try:
        entries = job["build_metadata"]["entries"]
        build_root = job["build_metadata"]["build_root"]
        configuration_attestation = job["build_metadata"][
            "effective_configuration_attestation"
        ]
    except (KeyError, TypeError):
        raise CrossoverError(
            "job omits CMake effective-configuration attestation"
        )
    if not isinstance(build_root, str):
        raise CrossoverError("job omits its absolute CMake build directory")
    try:
        build_root_path = checked_absolute_path(
            build_root, "job CMake build directory"
        )
    except CrossoverError:
        raise CrossoverError("job omits its absolute CMake build directory")
    home_directory = entries.get("CMAKE_HOME_DIRECTORY")
    if not isinstance(home_directory, str) or not home_directory:
        raise CrossoverError(
            "job omits its absolute CMake source directory"
        )
    try:
        checked_absolute_path(
            home_directory, "job CMake source directory"
        )
    except CrossoverError:
        raise CrossoverError(
            "job omits its absolute CMake source directory"
        )
    attested_configuration_sha256 = (
        validate_build_configuration_attestation(
            configuration_attestation,
            build_root_path / BUILD_CONFIGURATION_RELATIVE_PATH
        )
    )
    (configuration_file_schema, unused_variables,
     configuration_selectors) = build_configuration_contract(
        configuration_attestation["schema"])
    del unused_variables
    if (entries.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") !=
            configuration_file_schema or
            entries.get(
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"
            ) != attested_configuration_sha256 or
            any(entries.get(variable) !=
                configuration_attestation["entries"].get(variable)
                for variable in
                configuration_selectors)):
        raise CrossoverError(
            "job effective CMake configuration differs from its cache binding"
        )
    validate_embedded_build_type(
        configuration_attestation["entries"], build.get("build_type"),
        settings.get("mode") in AUTHORITATIVE_COMMANDS
    )
    if (build.get("build_configuration_sha256") !=
            attested_configuration_sha256):
        raise CrossoverError(
            "benchmark executable build configuration differs from CMake"
        )
    if (not isinstance(build.get("compiler"), str) or
            not build["compiler"] or
            not isinstance(build.get("compiler_version"), str) or
            not build["compiler_version"] or
            type(build.get("cplusplus")) is not int or
            build["cplusplus"] < 201103):
        raise CrossoverError("benchmark build metadata has invalid exact types")
    parity_indices = requested_indices(item["mask"], item["r"])
    if not parity_indices or len(parity_indices) != item["q"]:
        raise CrossoverError(
            "generated parity mask must select exactly Q nonempty outputs"
        )
    expected = {
        "K": item["k"], "R": item["r"], "Q": item["q"],
        "forced_mode": "force_" + timed_mode,
        "shard_bytes": item["shard_bytes"],
        "seed": job["seed"],
        "requested_parity_indices": parity_indices,
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
        if canonical_bytes(parameters.get(key)) != canonical_bytes(value):
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
    if correctness.get("selected_transform_reference_parity_match") is not True:
        raise CrossoverError("benchmark selected/reference parity check failed")
    if (timed_mode == "direct" and
            correctness.get("direct_transform_parity_match") is not True):
        raise CrossoverError("benchmark direct/transform correctness check failed")
    if (timed_mode == "transform" and
            correctness.get("direct_transform_parity_match") is not None):
        raise CrossoverError(
            "force-transform must not claim a direct/transform comparison"
        )
    if correctness.get("unrequested_outputs_untouched") is not True:
        raise CrossoverError("benchmark unrequested-output guard check failed")
    checksum = correctness.get("parity_checksum_fnv1a64")
    if not isinstance(checksum, str) or not re.fullmatch(r"0x[0-9a-fA-F]{16}", checksum):
        raise CrossoverError("benchmark omitted a well-formed parity checksum")
    checksum = checksum.lower()
    if (resolved.get("direct_capable") is not True or
            type(resolved.get("thread_count")) is not int or
            resolved.get("thread_count") != 1 or
            type(resolved.get("parent_count")) is not int or
            resolved.get("parent_count") <= 0 or
            type(resolved.get("padded_side")) is not int or
            resolved.get("padded_side") <= 0):
        raise CrossoverError("benchmark resolved unexpected capability or thread count")
    expected_direct = timed_mode == "direct"
    if resolved.get("timed_path_is_direct") is not expected_direct:
        raise CrossoverError("benchmark timed a path other than the forced path")

    batch = settings["benchmark"]["batch"]
    direct_terms = item["k"] * item["q"] * batch
    direct_initializations = item["q"] * batch
    direct_accumulations = (item["k"] - 1) * item["q"] * batch
    symbol_bytes = 2 if item["field"] == "gf16" else 1
    if item["shard_bytes"] % symbol_bytes != 0:
        raise CrossoverError("benchmark cell contains a partial field symbol")
    symbols_per_shard = item["shard_bytes"] // symbol_bytes
    expected_operation_model = {
        "direct_model_applies_to_timed_path": expected_direct,
        "direct_output_accumulations": direct_accumulations,
        "direct_output_initializations": direct_initializations,
        "direct_row_terms": direct_terms,
        "fixed_coefficient_symbol_terms_before_unit_specialization":
            direct_terms * symbols_per_shard,
        "hardware_counters": None,
        "model_scope": (
            "direct streaming kernels before cache effects; unit "
            "coefficients specialize to copy/XOR"
        ),
        "modeled_output_bytes_read":
            direct_accumulations * item["shard_bytes"],
        "modeled_output_bytes_written":
            direct_terms * item["shard_bytes"],
        "modeled_source_bytes_read":
            direct_terms * item["shard_bytes"],
        "transform_operation_counts": None,
        "xor_accumulation_symbols":
            direct_accumulations * symbols_per_shard,
    }
    if canonical_bytes(operation_model) != canonical_bytes(
            expected_operation_model):
        raise CrossoverError(
            "benchmark operation_model differs from the exact v2 model"
        )

    median_us = required_finite_metric(
        execution.get("median_us_per_batch_call"),
        "metrics.encode_execution.median_us_per_batch_call", False
    )
    mad_us = required_finite_metric(
        execution.get("mad_us_per_batch_call"),
        "metrics.encode_execution.mad_us_per_batch_call", True
    )
    minimum_us = required_finite_metric(
        execution.get("minimum_us_per_batch_call"),
        "metrics.encode_execution.minimum_us_per_batch_call", False
    )
    maximum_us = required_finite_metric(
        execution.get("maximum_us_per_batch_call"),
        "metrics.encode_execution.maximum_us_per_batch_call", False
    )
    if minimum_us > median_us or median_us > maximum_us:
        raise CrossoverError(
            "benchmark encode-execution minimum/median/maximum are inconsistent"
        )
    for key in (
            "logical_input_GB_per_s",
            "logical_input_plus_output_GB_per_s",
            "requested_parity_output_GB_per_s"):
        required_finite_metric(
            execution.get(key), "metrics.encode_execution." + key, True
        )
    hashed_bytes = (
        settings["benchmark"]["batch"] * item["q"] * item["shard_bytes"]
    )
    if hashed_bytes <= 0:
        raise CrossoverError("benchmark parity identity covers no output bytes")
    parity_identity = {
        "algorithm": "fnv1a64",
        "digest": checksum,
        "hashed_bytes": hashed_bytes,
        "requested_parity_indices": parity_indices,
    }
    return median_us, mad_us, parity_identity


def validate_execution_inputs(job, result_root=None):
    executable = Path(job["executable"])
    executable_artifact = job.get("executable_artifact")
    if executable_artifact is not None:
        validate_frozen_executable(
            executable_artifact, build_metadata=job["build_metadata"],
            source_state=job["source_identity"], result_root=result_root,
        )
        if ((Path(os.path.abspath(os.fspath(executable))) !=
             Path(os.path.abspath(executable_artifact["executable"]))) or
                job.get("executable_sha256") !=
                executable_artifact["executable_sha256"]):
            raise CrossoverError(
                "job executable differs from its frozen artifact"
            )
        return
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


def validate_job_artifacts(
        result_dir, result, expected_job, settings, result_root):
    validate_owned_directory(result_root, "canonical result directory")
    result_dir = result_root["path"]
    if not isinstance(result, dict) or result.get("schema") != JOB_SCHEMA:
        raise CrossoverError("job result has a legacy or unknown schema")
    common_keys = {
        "build_metadata", "cell", "commands", "configuration_id",
        "executable", "executable_artifact", "executable_sha256", "job_id",
        "isolation", "measurements", "reason", "resumed", "schema", "seed",
        "source_identity", "status",
    }
    if result.get("status") != "passed":
        if set(result) != common_keys or result.get("status") != "failed":
            raise CrossoverError("failed job result has an invalid v6 shape")
        return
    if set(result) != common_keys | {"parity_identity", "summary"}:
        raise CrossoverError("passed job result has an invalid v6 shape")
    if result.get("reason") != "" or not isinstance(result.get("resumed"), bool):
        raise CrossoverError("passed job contains inconsistent status metadata")
    if settings.get("mode") in AUTHORITATIVE_COMMANDS:
        isolation_settings = settings["isolation"]
        try:
            support = load_isolation_support(
                expected_job["build_metadata"]["entries"][
                    "CMAKE_HOME_DIRECTORY"
                ]
            )
            support.validate_isolation(
                result.get("isolation"),
                isolation_settings["benchmark_cpu"],
                isolation_settings["reserved_sibling"],
                require_accepted=True,
            )
            support.validate_pair_lease_identity(
                isolation_settings["pair_lease"],
                isolation_settings["benchmark_cpu"],
                isolation_settings["reserved_sibling"],
            )
            if canonical_bytes(
                    result["isolation"]["pair_lease"]) != canonical_bytes(
                    isolation_settings["pair_lease"]):
                raise CrossoverError(
                    "job isolation pair lease differs from campaign settings"
                )
        except Exception as error:
            raise CrossoverError(
                "passed job isolation evidence is invalid: {}".format(error)
            )
    elif result.get("isolation") is not None:
        raise CrossoverError("screening job unexpectedly claims isolation")
    for key in (
            "build_metadata", "cell", "configuration_id", "executable",
            "executable_artifact", "executable_sha256", "job_id", "seed",
            "source_identity"):
        if canonical_bytes(result.get(key)) != canonical_bytes(expected_job.get(key)):
            raise CrossoverError(
                "passed job {} has stale {} metadata".format(
                    expected_job["job_id"], key
                )
            )
    validate_execution_inputs(expected_job, result_root)
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
    parity_identities = set()
    log_root = result_dir / "logs" / expected_job["job_id"]
    try:
        source = Path(
            expected_job["build_metadata"]["entries"]["CMAKE_HOME_DIRECTORY"]
        ).resolve()
    except (KeyError, TypeError, ValueError):
        raise CrossoverError("job build provenance omits the source root")
    command_keys = {
        "argv", "cwd", "environment", "launcher_identity", "returncode",
        "stderr_log",
        "stderr_sha256", "stdout_log", "stdout_sha256", "timed_out",
    }
    expected_launcher_identity = current_runtime_launcher_identity()
    measurement_keys = {
        "benchmark_json", "benchmark_json_sha256", "mad_us", "median_us",
        "parity_identity", "sequence_index", "timed_mode",
    }
    for index, timed_mode in enumerate(order):
        command = commands[index]
        measurement = measurements[index]
        if (not isinstance(command, dict) or set(command) != command_keys or
                not isinstance(measurement, dict) or
                set(measurement) != measurement_keys):
            raise CrossoverError(
                "passed job {} has a malformed sequence entry {}".format(
                    expected_job["job_id"], index
                )
            )
        label = "{:02d}-{}".format(index, timed_mode)
        raw_path = (
            result_dir / "raw" / expected_job["job_id"] / (label + ".json")
        ).resolve()
        child_raw_path = (
            Path("/proc/self/fd") / str(RAW_OUTPUT_DESCRIPTOR)
        )
        expected_relative_raw = str(raw_path.relative_to(result_dir.resolve()))
        expected_environment = dict(BENCHMARK_ENVIRONMENT)
        expected_environment["LEO2_EXPECT_BACKEND"] = (
            expected_job["cell"]["backend"]
        )
        if (type(command.get("returncode")) is not int or
                command.get("returncode") != 0 or
                command.get("timed_out") is not False or
                type(measurement.get("sequence_index")) is not int or
                measurement.get("sequence_index") != index or
                measurement.get("timed_mode") != timed_mode or
                command.get("argv") != benchmark_argv(
                    expected_job, timed_mode, child_raw_path, settings
                ) or canonical_bytes(
                    command.get("launcher_identity")
                ) != canonical_bytes(expected_launcher_identity) or
                command.get("cwd") != str(source) or
                command.get("environment") != expected_environment or
                command.get("stdout_log") != label + ".stdout.log" or
                command.get("stderr_log") != label + ".stderr.log" or
                measurement.get("benchmark_json") != expected_relative_raw):
            raise CrossoverError(
                "passed job {} has an invalid sequence entry {}".format(
                    expected_job["job_id"], index
                )
            )
        for stream in ("stdout", "stderr"):
            name = command.get(stream + "_log")
            if not isinstance(name, str) or Path(name).name != name:
                raise CrossoverError("job log name is missing or unsafe")
            path = log_root / name
            value = read_result_regular(
                result_root, path.relative_to(result_dir),
                MAX_RETAINED_LOG_BYTES, stream + " log"
            )
            if digest_bytes(value) != command.get(stream + "_sha256"):
                raise CrossoverError("{} hash does not match the job record".format(path))
        raw_bytes = read_result_regular(
            result_root, raw_path.relative_to(result_dir),
            MAX_RAW_JSON_BYTES, "benchmark JSON"
        )
        raw = decode_json_bytes(raw_bytes, str(raw_path))
        if digest_bytes(raw_bytes) != measurement.get("benchmark_json_sha256"):
            raise CrossoverError("{} hash does not match the job record".format(raw_path))
        median_us, mad_us, parity_identity = validate_raw(
            raw, expected_job, timed_mode, settings
        )
        recorded_median = required_finite_metric(
            measurement.get("median_us"),
            "job.measurements[{}].median_us".format(index), False
        )
        recorded_mad = required_finite_metric(
            measurement.get("mad_us"),
            "job.measurements[{}].mad_us".format(index), True
        )
        if recorded_median != median_us or recorded_mad != mad_us:
            raise CrossoverError("raw metrics differ from the passed job summary")
        if canonical_bytes(measurement.get("parity_identity")) != canonical_bytes(
                parity_identity):
            raise CrossoverError(
                "raw parity identity differs from the passed job summary"
            )
        parity_identities.add(canonical_bytes(parity_identity))
    if len(parity_identities) != 1:
        raise CrossoverError(
            "passed job raw artifacts contain different parity identities"
        )
    parity_identity = json.loads(next(iter(parity_identities)).decode("utf-8"))
    if canonical_bytes(result.get("parity_identity")) != canonical_bytes(
            parity_identity):
        raise CrossoverError("passed job omits its non-vacuous parity identity")
    recomputed = summarize_measurements(measurements)
    if canonical_bytes(recomputed) != canonical_bytes(result.get("summary")):
        raise CrossoverError("passed job aggregate does not match its raw measurements")


def median(values):
    try:
        result = float(statistics.median(values))
    except (OverflowError, TypeError, ValueError, statistics.StatisticsError) as error:
        raise CrossoverError("cannot compute finite median: {}".format(error))
    if not math.isfinite(result):
        raise CrossoverError("cannot compute finite median")
    return result


def relative_gain_percent(candidate, reference):
    try:
        result = (reference / candidate - 1.0) * 100.0
    except (OverflowError, ZeroDivisionError):
        result = float("nan")
    if not math.isfinite(result):
        raise CrossoverError("cannot compute a finite encoder gain")
    return result


def paired_log_inference(rounds):
    if len(rounds) != 3:
        raise CrossoverError(
            "authoritative ABBA inference requires exactly three rounds"
        )
    contrasts = []
    for index, round_value in enumerate(rounds):
        round_value = require_exact_keys(round_value, {
            "direct_geometric_mean_us", "gain_percent", "log_contrast",
            "transform_geometric_mean_us",
        }, "rounds[{}]".format(index))
        direct = required_finite_metric(
            round_value.get("direct_geometric_mean_us"),
            "rounds[{}].direct_geometric_mean_us".format(index), False
        )
        transform = required_finite_metric(
            round_value.get("transform_geometric_mean_us"),
            "rounds[{}].transform_geometric_mean_us".format(index), False
        )
        contrast = required_finite_number(
            round_value.get("log_contrast"),
            "rounds[{}].log_contrast".format(index)
        )
        expected_contrast = math.log(transform) - math.log(direct)
        recorded_gain = required_finite_number(
            round_value.get("gain_percent"),
            "rounds[{}].gain_percent".format(index)
        )
        expected_gain = relative_gain_percent(direct, transform)
        if not math.isclose(
                contrast, expected_contrast, rel_tol=1e-15, abs_tol=1e-15):
            raise CrossoverError(
                "round {} log contrast does not match its geometric means".format(
                    index
                )
            )
        if not math.isclose(
                recorded_gain, expected_gain, rel_tol=1e-12, abs_tol=1e-12):
            raise CrossoverError(
                "round {} gain does not match its geometric means".format(index)
            )
        contrasts.append(contrast)
    mean = statistics.mean(contrasts)
    standard_error = statistics.stdev(contrasts) / math.sqrt(3.0)
    # Two-sided 95% Student-t critical value at df=2.
    margin = 4.302652729911275 * standard_error
    try:
        lower = math.exp(mean - margin)
        upper = math.exp(mean + margin)
        speedup = math.exp(mean)
    except OverflowError:
        raise CrossoverError(
            "paired-log confidence interval overflowed finite arithmetic"
        )
    if not all(math.isfinite(value) and value > 0 for value in (
            lower, speedup, upper)):
        raise CrossoverError("paired-log confidence interval is not finite")
    return {
        "confidence": 0.95,
        "degrees_of_freedom": 2,
        "gain_percent": (speedup - 1.0) * 100.0,
        "gain_percent_student_t_interval": [
            (lower - 1.0) * 100.0,
            (upper - 1.0) * 100.0,
        ],
        "log_contrasts": contrasts,
        "speedup_geometric_mean": speedup,
        "speedup_student_t_interval": [lower, upper],
    }


def summarize_measurements(measurements):
    by_mode = {"direct": [], "transform": []}
    if not isinstance(measurements, list):
        raise CrossoverError("job measurements must be a JSON array")
    for index, measurement in enumerate(measurements):
        measurement = required_mapping(
            measurement, "job.measurements[{}]".format(index)
        )
        timed_mode = measurement.get("timed_mode")
        if timed_mode not in by_mode:
            raise CrossoverError("job measurement has an unknown timed mode")
        by_mode[timed_mode].append(required_finite_metric(
            measurement.get("median_us"),
            "job.measurements[{}].median_us".format(index), False
        ))
    if not by_mode["direct"] or not by_mode["transform"]:
        raise CrossoverError("job did not measure both encoder paths")
    direct = median(by_mode["direct"])
    transform = median(by_mode["transform"])
    gain = relative_gain_percent(direct, transform)
    rounds = []
    if len(measurements) % 4 == 0 and len(measurements) >= 4:
        for offset in range(0, len(measurements), 4):
            group = measurements[offset:offset + 4]
            if [item["timed_mode"] for item in group] != [
                    "direct", "transform", "transform", "direct"]:
                rounds = []
                break
            direct_logs = [
                math.log(required_finite_metric(
                    group[position]["median_us"],
                    "ABBA direct invocation", False
                )) for position in (0, 3)
            ]
            transform_logs = [
                math.log(required_finite_metric(
                    group[position]["median_us"],
                    "ABBA transform invocation", False
                )) for position in (1, 2)
            ]
            round_direct = math.exp(statistics.mean(direct_logs))
            round_transform = math.exp(statistics.mean(transform_logs))
            log_contrast = (
                statistics.mean(transform_logs) -
                statistics.mean(direct_logs)
            )
            rounds.append({
                "direct_geometric_mean_us": round_direct,
                "gain_percent": relative_gain_percent(
                    round_direct, round_transform
                ),
                "log_contrast": log_contrast,
                "transform_geometric_mean_us": round_transform,
            })
    inference = paired_log_inference(rounds) if len(rounds) == 3 else None
    return {
        "direct_invocation_medians_us": by_mode["direct"],
        "direct_median_us": direct,
        "gain_percent": gain,
        "paired_log_inference": inference,
        "rounds": rounds,
        "transform_invocation_medians_us": by_mode["transform"],
        "transform_median_us": transform,
    }


def run_job(job, context):
    isolation_context = context.get("isolation")
    campaign_result_root = context.get("result_root")
    if campaign_result_root is None:
        result_root = owned_canonical_directory(
            context["result_dir"], "job canonical result directory"
        )
    else:
        result_root = duplicate_owned_directory(
            campaign_result_root, "campaign canonical result directory"
        )

    def validate_guards():
        if isolation_context is not None:
            validate_authoritative_guards(
                isolation_context["canonical_guard"],
                isolation_context["pair_guard"],
            )

    result_path = context["result_dir"] / "jobs" / (job["job_id"] + ".json")
    if context["resume"]:
        try:
            validate_guards()
            try:
                previous = decode_json_bytes(
                    read_result_regular(
                        result_root,
                        result_path.relative_to(context["result_dir"]),
                        MAX_RAW_JSON_BYTES, "resume job JSON"
                    ),
                    str(result_path)
                )
                if (not isinstance(previous, dict) or
                        previous.get("schema") != JOB_SCHEMA):
                    raise CrossoverError(
                        "resume job has a legacy or unknown schema"
                    )
                if (previous.get("configuration_id") ==
                        job["configuration_id"] and
                        previous.get("status") == "passed"):
                    validate_job_artifacts(
                        context["result_dir"], previous, job,
                        context["settings"], result_root,
                    )
                    validate_guards()
                    close_owned_directory(result_root)
                    return previous
            except (CrossoverError, OSError):
                validate_guards()
        except BaseException as primary:
            try:
                close_owned_directory(result_root)
            except BaseException as cleanup:
                raise CrossoverError(
                    "resume validation failed: {}; result-root cleanup "
                    "failed: {}: {}".format(
                        primary, type(cleanup).__name__, cleanup
                    )
                ) from primary
            raise

    log_relative = Path("logs") / job["job_id"]
    raw_relative = Path("raw") / job["job_id"]
    log_dir = context["result_dir"] / log_relative
    raw_dir = context["result_dir"] / raw_relative
    result = {
        "build_metadata": job["build_metadata"],
        "cell": job["cell"],
        "commands": [],
        "configuration_id": job["configuration_id"],
        "executable": job["executable"],
        "executable_artifact": job["executable_artifact"],
        "executable_sha256": job["executable_sha256"],
        "job_id": job["job_id"],
        "isolation": None,
        "measurements": [],
        "reason": "",
        "resumed": False,
        "schema": JOB_SCHEMA,
        "seed": job["seed"],
        "source_identity": job["source_identity"],
        "status": "failed",
    }
    environment = dict(BENCHMARK_ENVIRONMENT)
    environment["LEO2_EXPECT_BACKEND"] = job["cell"]["backend"]
    parity_identities = set()
    isolation_before = None
    execution_file = None
    execution_snapshot = None
    taskset_snapshot = None
    log_directory = None
    raw_directory = None
    try:
        log_directory = open_result_directory(
            result_root, log_relative, "job log directory", True
        )
        raw_directory = open_result_directory(
            result_root, raw_relative, "job raw-output directory", True
        )
    except BaseException as primary:
        cleanup_errors = []
        for label, directory in (
                ("raw directory", raw_directory),
                ("log directory", log_directory),
                ("result root", result_root)):
            if directory is None:
                continue
            try:
                close_owned_directory(directory)
            except BaseException as cleanup:
                cleanup_errors.append((label, cleanup))
        if cleanup_errors:
            raise CrossoverError(
                "job directory setup failed: {}; cleanup: {}".format(
                    primary,
                    "; ".join(
                        "{} {}: {}".format(
                            label, type(error).__name__, error
                        )
                        for label, error in cleanup_errors
                    ),
                )
            ) from primary
        raise
    unhandled_error = None
    try:
        if context["settings"].get("mode") in AUTHORITATIVE_COMMANDS:
            taskset_snapshot = open_exact_executable_snapshot(
                context["settings"]["taskset"],
                "authoritative taskset executable", True,
            )
            if (taskset_snapshot["sha256"] !=
                    context["settings"]["taskset_sha256"]):
                raise CrossoverError(
                    "authoritative taskset snapshot differs from settings"
                )
        if isolation_context is not None:
            support = isolation_context["support"]
            validate_guards()
            isolation_before = {
                "monotonic_ns": time.monotonic_ns(),
                "benchmark_cpu": support.cpu_stat_snapshot(
                    isolation_context["cpu"]
                ),
                "reserved_sibling": support.cpu_stat_snapshot(
                    isolation_context["sibling"]
                ),
            }
        if job.get("executable_artifact") is not None:
            try:
                executable_relative = Path(os.path.abspath(
                    job["executable_artifact"]["executable"]
                )).relative_to(result_root["path"])
            except (KeyError, TypeError, ValueError):
                raise CrossoverError(
                    "frozen job executable escapes the held result root"
                )
            execution_file = open_result_regular_held(
                result_root, executable_relative, MAX_RAW_JSON_BYTES,
                "frozen job executable", 0o555,
                directory_required_mode=0o555,
            )
            if digest_bytes(read_held_regular(
                    execution_file, MAX_RAW_JSON_BYTES,
                    "frozen job executable")) != job["executable_sha256"]:
                raise CrossoverError(
                    "frozen job executable descriptor hash is stale"
                )
            execution_snapshot = create_sealed_executable_snapshot(
                read_held_regular(
                    execution_file, MAX_RAW_JSON_BYTES,
                    "frozen job executable",
                ),
                "frozen job executable",
            )
            validate_sealed_executable_snapshot(
                execution_snapshot, job["executable_sha256"],
                "frozen job executable snapshot",
            )
        # Failed jobs own these deterministic names.  Remove only verified
        # regular files before retrying so partial evidence cannot make
        # --resume permanently fail on exclusive log creation.
        for index, timed_mode in enumerate(job["invocation_order"]):
            label = "{:02d}-{}".format(index, timed_mode)
            remove_stale_owned_file(
                raw_directory, label + ".json", "benchmark JSON"
            )
            remove_stale_owned_file(
                log_directory, label + ".stdout.log", "stdout log"
            )
            remove_stale_owned_file(
                log_directory, label + ".stderr.log", "stderr log"
            )
        for index, timed_mode in enumerate(job["invocation_order"]):
            validate_guards()
            validate_execution_inputs(job, result_root)
            label = "{:02d}-{}".format(index, timed_mode)
            raw_path = raw_dir / (label + ".json")
            raw_file = open_exclusive_owned_file(
                raw_directory, raw_path.name, "benchmark JSON"
            )
            stdout_path = log_dir / (label + ".stdout.log")
            stderr_path = log_dir / (label + ".stderr.log")
            try:
                child_raw_path = procfd_exact_path(
                    raw_file, "benchmark JSON"
                )
                argv = benchmark_argv(
                    job, timed_mode, child_raw_path, context["settings"]
                )
                inherited_descriptors = [raw_file["descriptor"]]
                descriptor_remaps = {
                    RAW_OUTPUT_DESCRIPTOR: raw_file["descriptor"],
                }
                if execution_snapshot is not None:
                    inherited_descriptors.append(
                        execution_snapshot["descriptor"]
                    )
                    descriptor_remaps[EXECUTABLE_DESCRIPTOR] = (
                        execution_snapshot["descriptor"]
                    )
                if taskset_snapshot is not None:
                    inherited_descriptors.append(
                        taskset_snapshot["descriptor"]
                    )
                    descriptor_remaps[TASKSET_EXECUTABLE_DESCRIPTOR] = (
                        taskset_snapshot["descriptor"]
                    )
                bound_stdout_path = (
                    Path("/proc/self/fd") /
                    str(result_root["descriptor"]) /
                    stdout_path.relative_to(context["result_dir"])
                )
                bound_stderr_path = (
                    Path("/proc/self/fd") /
                    str(result_root["descriptor"]) /
                    stderr_path.relative_to(context["result_dir"])
                )
                inherited_descriptors.append(result_root["descriptor"])
                command = run_command(
                    argv, context["source"], bound_stdout_path,
                    bound_stderr_path,
                    context["timeout"], environment,
                    inherited_lock_descriptor=(
                        isolation_context["canonical_guard"].descriptor
                        if isolation_context is not None else None
                    ),
                    inherited_descriptors=tuple(inherited_descriptors),
                    descriptor_remaps=descriptor_remaps,
                )
                command["environment"] = dict(environment)
                result["commands"].append(command)
                validate_guards()
                validate_execution_inputs(job, result_root)
                if execution_snapshot is not None:
                    validate_sealed_executable_snapshot(
                        execution_snapshot, job["executable_sha256"],
                        "executed frozen executable snapshot",
                    )
                if taskset_snapshot is not None:
                    validate_sealed_executable_snapshot(
                        taskset_snapshot,
                        context["settings"]["taskset_sha256"],
                        "executed taskset snapshot",
                    )
                if command["returncode"] != 0:
                    raise CrossoverError(
                        "{} exited with status {}".format(
                            label, command["returncode"]
                        )
                    )
                raw_bytes = read_held_regular(
                    raw_file, MAX_RAW_JSON_BYTES, "benchmark JSON"
                )
            finally:
                primary = sys.exc_info()[1]
                try:
                    close_log(raw_file)
                except BaseException as cleanup:
                    if primary is not None:
                        raise CrossoverError(
                            "job invocation failed: {}; raw-file cleanup "
                            "also failed: {}: {}".format(
                                primary, type(cleanup).__name__, cleanup
                            )
                        ) from primary
                    raise
            raw = decode_json_bytes(raw_bytes, str(raw_path))
            median_us, mad_us, parity_identity = validate_raw(
                raw, job, timed_mode, context["settings"]
            )
            parity_identities.add(canonical_bytes(parity_identity))
            result["measurements"].append({
                "benchmark_json": str(raw_path.relative_to(context["result_dir"])),
                "benchmark_json_sha256": digest_bytes(raw_bytes),
                "mad_us": mad_us,
                "median_us": median_us,
                "parity_identity": parity_identity,
                "sequence_index": index,
                "timed_mode": timed_mode,
            })
        if len(parity_identities) != 1:
            raise CrossoverError("forced paths emitted different parity identities")
        result["parity_identity"] = json.loads(
            next(iter(parity_identities)).decode("utf-8")
        )
        result["summary"] = summarize_measurements(result["measurements"])
        result["status"] = "passed"
    except (CrossoverError, OSError, ValueError) as error:
        result["reason"] = str(error)
    except BaseException as error:
        unhandled_error = error

    teardown_errors = []
    if isolation_context is not None and isolation_before is not None:
        try:
            support = isolation_context["support"]
            after_cpu = support.cpu_stat_snapshot(
                isolation_context["cpu"]
            )
            after_sibling = support.cpu_stat_snapshot(
                isolation_context["sibling"]
            )
            after_monotonic = time.monotonic_ns()
            validate_guards()
            result["isolation"] = support.isolation_record(
                isolation_context["cpu"],
                isolation_context["sibling"],
                isolation_context["pair_lease"],
                isolation_before["monotonic_ns"],
                after_monotonic,
                isolation_before["benchmark_cpu"],
                after_cpu,
                isolation_before["reserved_sibling"],
                after_sibling,
            )
            support.validate_isolation(
                result["isolation"],
                isolation_context["cpu"],
                isolation_context["sibling"],
                require_accepted=True,
            )
        except BaseException as error:
            if isinstance(error, Exception):
                result["status"] = "failed"
                result["reason"] = (
                    result["reason"] + "; " if result["reason"] else ""
                ) + "isolation validation failed: {}".format(error)
                result.pop("parity_identity", None)
                result.pop("summary", None)
            else:
                teardown_errors.append(("isolation validation", error))
    for label, cleanup in (
            ("taskset snapshot",
             (lambda: close_sealed_snapshot(taskset_snapshot))
             if taskset_snapshot is not None else None),
            ("executable snapshot",
             (lambda: close_sealed_snapshot(execution_snapshot))
             if execution_snapshot is not None else None),
            ("frozen executable",
             (lambda: close_log(execution_file))
             if execution_file is not None else None),
            ("raw directory",
             (lambda: close_owned_directory(raw_directory))),
            ("log directory",
             (lambda: close_owned_directory(log_directory)))):
        if cleanup is None:
            continue
        try:
            cleanup()
        except BaseException as error:
            teardown_errors.append((label, error))
    if unhandled_error is not None or teardown_errors:
        try:
            close_owned_directory(result_root)
        except BaseException as error:
            teardown_errors.append(("result root", error))
    if unhandled_error is not None:
        if teardown_errors:
            raise CrossoverError(
                "job execution raised {}: {}; teardown also failed: {}".format(
                    type(unhandled_error).__name__, unhandled_error,
                    "; ".join(
                        "{} {}: {}".format(
                            label, type(error).__name__, error
                        )
                        for label, error in teardown_errors
                    ),
                )
            ) from unhandled_error
        raise unhandled_error
    if teardown_errors:
        raise CrossoverError(
            "job teardown failed{}: {}".format(
                " after " + result["reason"] if result["reason"] else "",
                "; ".join(
                    "{} {}: {}".format(
                        label, type(error).__name__, error
                    )
                    for label, error in teardown_errors
                ),
            )
        )
    publication_error = None
    try:
        validate_owned_directory(result_root, "canonical result directory")
        write_result_json(
            result_root, result_path.relative_to(context["result_dir"]),
            result, "job result JSON", replace=True,
        )
        validate_owned_directory(result_root, "canonical result directory")
    except BaseException as error:
        publication_error = error
    close_error = None
    try:
        close_owned_directory(result_root)
    except BaseException as error:
        close_error = error
    if publication_error is not None:
        if close_error is not None:
            raise CrossoverError(
                "job-result publication failed: {}; result-root close also "
                "failed: {}: {}".format(
                    publication_error, type(close_error).__name__, close_error
                )
            ) from publication_error
        raise publication_error
    if close_error is not None:
        raise CrossoverError(
            "job-result publication succeeded but result-root close failed: "
            "{}".format(close_error)
        ) from close_error
    return result


def summarize_region(results, promotion_percent):
    passed = [item for item in results if item.get("status") == "passed"]
    gains = []
    round_gains = []
    confidence_lower_gains = []
    confidence_missing_count = 0
    for index, item in enumerate(passed):
        summary = required_mapping(
            item.get("summary"), "jobs[{}].summary".format(index)
        )
        gains.append(required_finite_number(
            summary.get("gain_percent"),
            "jobs[{}].summary.gain_percent".format(index)
        ))
        rounds = summary.get("rounds")
        if not isinstance(rounds, list):
            raise CrossoverError("passed job summary rounds must be a JSON array")
        for round_index, round_value in enumerate(rounds):
            round_value = required_mapping(
                round_value,
                "jobs[{}].summary.rounds[{}]".format(index, round_index)
            )
            round_gains.append(required_finite_number(
                round_value.get("gain_percent"),
                "jobs[{}].summary.rounds[{}].gain_percent".format(
                    index, round_index
                )
            ))
        inference = summary.get("paired_log_inference")
        if inference is None:
            confidence_missing_count += 1
        else:
            inference = required_mapping(
                inference,
                "jobs[{}].summary.paired_log_inference".format(index)
            )
            interval = inference.get("gain_percent_student_t_interval")
            if (not isinstance(interval, list) or len(interval) != 2):
                raise CrossoverError("passed job has an invalid confidence interval")
            lower = required_finite_number(
                interval[0],
                "jobs[{}].summary.paired_log_inference.lower".format(index)
            )
            upper = required_finite_number(
                interval[1],
                "jobs[{}].summary.paired_log_inference.upper".format(index)
            )
            if lower > upper:
                raise CrossoverError("passed job confidence interval is reversed")
            confidence_lower_gains.append(lower)
    gains.sort()
    round_gains.sort()
    return {
        "cell_count": len(results),
        "confidence_interval_missing_count": confidence_missing_count,
        "confidence_promotion_count": sum(
            value >= promotion_percent for value in confidence_lower_gains
        ),
        "confidence_worst_lower_gain_percent": (
            min(confidence_lower_gains) if confidence_lower_gains else None
        ),
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
    analysis["candidate"][
        "all_cells_confidently_meet_promotion_threshold"
    ] = bool(candidate) and (
        candidate_summary["failed_count"] == 0 and
        candidate_summary["confidence_interval_missing_count"] == 0 and
        candidate_summary["confidence_promotion_count"] ==
        candidate_summary["cell_count"]
    )
    return analysis


def manifest_identity(manifest):
    return {
        "evidence_contract": manifest["evidence_contract"],
        "executables": manifest["executables"],
        "jobs": manifest["jobs"],
        "machine": manifest["machine"],
        "settings": manifest["settings"],
        "source_fingerprint": manifest["source_fingerprint"],
    }


def validate_cell_value(value, description):
    expected_keys = {
        "backend", "field", "k", "mask", "profile", "q", "r", "region",
        "shard_bytes",
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        raise CrossoverError("{} has an invalid cell shape".format(description))
    if (value["backend"] not in KNOWN_BACKENDS or
            value["field"] not in ("gf8", "gf16") or
            value["profile"] not in ("low", "high") or
            not isinstance(value["region"], str) or not value["region"] or
            not isinstance(value["mask"], str)):
        raise CrossoverError("{} has invalid cell labels".format(description))
    for key in ("k", "r", "q", "shard_bytes"):
        if type(value[key]) is not int or value[key] <= 0:
            raise CrossoverError(
                "{} cell {} must be a positive integer".format(description, key)
            )
    if value["k"] > 16 or value["r"] > 16 or value["q"] > value["r"]:
        raise CrossoverError("{} is outside the bounded direct domain".format(
            description
        ))
    if len(requested_indices(value["mask"], value["r"])) != value["q"]:
        raise CrossoverError("{} mask cardinality differs from Q".format(description))


def validate_manifest_settings(settings, jobs, path):
    core_keys = {
        "abba_rounds", "benchmark", "frozen_executable_required", "mode",
        "pin_cpu", "placement_policy", "taskset", "taskset_sha256",
        "timeout_seconds_per_invocation", "workers",
    }
    mode = settings.get("mode")
    validate_supported_run_mode(mode, str(path))
    expected_keys = set(core_keys)
    if mode in AUTHORITATIVE_COMMANDS:
        expected_keys.add("isolation")
    if mode == "historical-avx2":
        expected_keys.update(("campaign", "controlled_build"))
    if set(settings) != expected_keys:
        raise CrossoverError("{} has invalid v6 settings fields".format(path))
    benchmark = settings.get("benchmark")
    if (not isinstance(benchmark, dict) or
            set(benchmark) != {"batch", "iterations", "reuse", "warmups"}):
        raise CrossoverError("{} has invalid benchmark settings".format(path))
    for key in ("batch", "iterations", "reuse", "warmups"):
        if type(benchmark[key]) is not int or benchmark[key] <= 0:
            raise CrossoverError("{} benchmark {} is invalid".format(path, key))
    for key in ("abba_rounds", "timeout_seconds_per_invocation", "workers"):
        if type(settings.get(key)) is not int or settings[key] < 0:
            raise CrossoverError("{} setting {} is invalid".format(path, key))
    if settings["timeout_seconds_per_invocation"] == 0 or settings["workers"] == 0:
        raise CrossoverError("{} has a zero execution setting".format(path))
    if not isinstance(settings.get("placement_policy"), str):
        raise CrossoverError("{} has invalid placement policy".format(path))
    frozen = settings.get("frozen_executable_required")
    if type(frozen) is not bool or frozen != (mode in AUTHORITATIVE_COMMANDS):
        raise CrossoverError("{} has inconsistent frozen policy".format(path))
    if mode in AUTHORITATIVE_COMMANDS:
        if (settings["abba_rounds"] != 3 or settings["workers"] != 1 or
                type(settings.get("pin_cpu")) is not int or
                not isinstance(settings.get("taskset"), str) or
                Path(settings["taskset"]) != Path("/usr/bin/taskset") or
                not isinstance(settings.get("taskset_sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", settings["taskset_sha256"])):
            raise CrossoverError("{} has invalid authoritative ABBA settings".format(
                path
            ))
        try:
            taskset_sha256 = digest_bytes(Path("/usr/bin/taskset").read_bytes())
        except OSError as error:
            raise CrossoverError(
                "{} cannot revalidate /usr/bin/taskset: {}".format(path, error)
            )
        if taskset_sha256 != settings["taskset_sha256"]:
            raise CrossoverError("{} taskset executable changed".format(path))
        isolation = settings.get("isolation")
        isolation_keys = {
            "allowed_cpu_set_at_launch", "benchmark_cpu", "canonical_lock",
            "child_environment", "housekeeping_cpu_set", "pair_lease",
            "reserved_sibling",
        }
        if not isinstance(isolation, dict) or set(isolation) != isolation_keys:
            raise CrossoverError(
                "{} has incomplete authoritative isolation settings".format(path)
            )
        cpu = isolation.get("benchmark_cpu")
        sibling = isolation.get("reserved_sibling")
        allowed = isolation.get("allowed_cpu_set_at_launch")
        housekeeping = isolation.get("housekeeping_cpu_set")
        if (type(cpu) is not int or type(sibling) is not int or cpu == sibling or
                cpu != settings["pin_cpu"] or
                not isinstance(allowed, list) or
                not isinstance(housekeeping, list) or
                any(type(value) is not int or value < 0
                    for value in allowed + housekeeping) or
                allowed != sorted(set(allowed)) or
                housekeeping != sorted(set(housekeeping)) or
                set(housekeeping) != set(allowed) - {cpu, sibling} or
                not housekeeping):
            raise CrossoverError(
                "{} has invalid authoritative CPU-pair settings".format(path)
            )
        if isolation.get("child_environment") != BENCHMARK_ENVIRONMENT:
            raise CrossoverError(
                "{} has an unsanitized benchmark environment".format(path)
            )
        lock = isolation.get("canonical_lock")
        if (not isinstance(lock, dict) or set(lock) != {
                "device", "inode", "lock", "mode", "path", "uid"} or
                lock.get("path") != str(AUTHORITATIVE_LOCK) or
                lock.get("lock") != "exclusive" or lock.get("mode") != 0o600 or
                lock.get("uid") != os.getuid() or
                any(type(lock.get(key)) is not int or lock[key] < 0
                    for key in ("device", "inode"))):
            raise CrossoverError(
                "{} has invalid canonical-lock identity".format(path)
            )
        pair = isolation.get("pair_lease")
        payload = pair.get("payload") if isinstance(pair, dict) else None
        if (not isinstance(pair, dict) or set(pair) != {
                "device", "directory_device", "directory_inode", "inode",
                "lock", "path", "payload", "sha256"} or
                not isinstance(payload, dict) or
                payload.get("cpus") != sorted((cpu, sibling)) or
                payload.get("schema") != "leopard2-cpu-pair-lease/v1"):
            raise CrossoverError(
                "{} has invalid physical-pair lease identity".format(path)
            )
        try:
            source_root = jobs[0]["build_metadata"]["entries"][
                "CMAKE_HOME_DIRECTORY"
            ]
            support = load_isolation_support(source_root)
            support.validate_pair_lease_identity(pair, cpu, sibling)
        except Exception as error:
            raise CrossoverError(
                "{} physical-pair lease failed canonical validation: {}".format(
                    path, error
                )
            )
    elif (settings["abba_rounds"] != 0 or settings.get("pin_cpu") is not None or
          settings.get("taskset") is not None or
          settings.get("taskset_sha256") is not None):
        raise CrossoverError("{} has invalid screening settings".format(path))
    if mode == "historical-avx2":
        cells = [job.get("cell") for job in jobs]
        expected_cells = historical_avx2_grid()
        expected_campaign = {
            "cell_count": len(expected_cells),
            "cells_sha256": digest_value(expected_cells),
            "name": "corrected-high-profile-explicit-avx2",
        }
        if (canonical_bytes(cells) != canonical_bytes(expected_cells) or
                canonical_bytes(settings.get("campaign")) !=
                canonical_bytes(expected_campaign) or benchmark != {
                    "batch": 1, "iterations": 15, "reuse": 64, "warmups": 4,
                }):
            raise CrossoverError(
                "{} is not the exact historical AVX2 campaign".format(path)
            )
        controlled = settings.get("controlled_build")
        if (not isinstance(controlled, dict) or set(controlled) != {
                "path", "schema", "sha256"} or
                controlled.get("schema") != CONTROLLED_BUILD_SCHEMA or
                controlled.get("path") != "controlled-build.json" or
                not isinstance(controlled.get("sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", controlled["sha256"])):
            raise CrossoverError(
                "{} has invalid controlled-build descriptor".format(path)
            )


def validate_controlled_build_held(
        settings, jobs, manifest_path, result_root):
    if settings.get("mode") != "historical-avx2":
        return
    descriptor = settings["controlled_build"]
    result_dir = result_root["path"]
    record_path = result_dir / descriptor["path"]
    record_bytes = read_result_regular(
        result_root, descriptor["path"], MAX_RAW_JSON_BYTES,
        "controlled build record",
    )
    if digest_bytes(record_bytes) != descriptor["sha256"]:
        raise CrossoverError("controlled build record hash is stale")
    record = decode_json_bytes(record_bytes, str(record_path))
    expected_keys = {
        "backend", "build_metadata", "build_root", "build_tools", "commands",
        "executable",
        "executable_sha256", "guard_identity", "schema", "source_identity",
    }
    if (not isinstance(record, dict) or set(record) != expected_keys or
            record.get("schema") != CONTROLLED_BUILD_SCHEMA or
            record.get("backend") != "avx2"):
        raise CrossoverError("controlled build record has an invalid schema")
    build_root = (result_dir / "controlled-build-avx2").resolve()
    executable = (build_root / "bench_leopard2_direct_encode").resolve()
    if (record.get("build_root") != str(build_root) or
            record.get("executable") != str(executable)):
        raise CrossoverError("controlled build paths differ from the campaign")
    executable_bytes = read_result_regular(
        result_root,
        executable.relative_to(result_dir),
        MAX_RAW_JSON_BYTES,
        "controlled build executable",
        required_mode=0o755,
    )
    if (digest_bytes(executable_bytes) != record.get("executable_sha256") or
            not os.access(str(executable), os.X_OK)):
        raise CrossoverError("controlled build executable changed")
    representative = jobs[0]
    if (canonical_bytes(record.get("build_metadata")) != canonical_bytes(
            representative.get("build_metadata")) or
            canonical_bytes(record.get("source_identity")) != canonical_bytes(
                representative.get("source_identity"))):
        raise CrossoverError(
            "controlled build identity differs from manifest jobs"
        )
    expected_guard_identity = {
        "canonical_lock": settings["isolation"]["canonical_lock"],
        "pair_lease": settings["isolation"]["pair_lease"],
    }
    if canonical_bytes(record.get("guard_identity")) != canonical_bytes(
            expected_guard_identity):
        raise CrossoverError(
            "controlled build guard identity differs from the campaign"
        )
    build_tools = record.get("build_tools")
    if (not isinstance(build_tools, dict) or
            set(build_tools) != {"cmake", "ninja"}):
        raise CrossoverError(
            "controlled build tool identity is malformed"
        )
    validate_current_executable_identity(
        build_tools["cmake"], "controlled CMake"
    )
    validate_current_executable_identity(
        build_tools["ninja"], "controlled Ninja"
    )
    commands = record.get("commands")
    if not isinstance(commands, list) or len(commands) != 2:
        raise CrossoverError("controlled build command sequence is incomplete")
    command_keys = {
        "argv", "cwd", "environment", "launcher_identity", "returncode",
        "stderr_log",
        "stderr_sha256", "stdout_log", "stdout_sha256", "timed_out",
    }
    if any(not isinstance(command, dict) or set(command) != command_keys or
           not isinstance(command.get("argv"), list) or
           not command["argv"] or
           any(not isinstance(item, str) or not item
               for item in command["argv"])
           for command in commands):
        raise CrossoverError("controlled build command shape is invalid")
    try:
        source = Path(
            record["build_metadata"]["entries"]["CMAKE_HOME_DIRECTORY"]
        ).resolve()
    except (KeyError, TypeError, ValueError):
        raise CrossoverError("controlled build omits its source directory")
    current_metadata = cmake_build_metadata(executable)
    if canonical_bytes(current_metadata) != canonical_bytes(
            record["build_metadata"]):
        raise CrossoverError(
            "controlled build CMake/object/link provenance changed"
        )
    validate_build_source_binding(
        current_metadata, source, record["source_identity"], "avx2",
        require_fresh=True
    )
    cmake = str(
        Path("/proc/self/fd") / str(CONTROLLED_CMAKE_DESCRIPTOR)
    )
    ninja = str(
        Path("/proc/self/fd") / str(CONTROLLED_NINJA_DESCRIPTOR)
    )
    if commands[0]["argv"][0] != cmake or commands[1]["argv"][0] != cmake:
        raise CrossoverError(
            "controlled build did not use its sealed CMake descriptor"
        )
    expected_argv = (
        controlled_avx2_configure_argv(cmake, ninja, source, build_root),
        [
            cmake, "--build", str(build_root), "--target",
            "bench_leopard2_direct_encode", "--parallel",
            str(min(len(settings["isolation"]["housekeeping_cpu_set"]), 128)),
        ],
    )
    expected_launcher_identity = current_runtime_launcher_identity()
    for index, command in enumerate(commands):
        if (command.get("argv") != expected_argv[index] or
                canonical_bytes(
                    command.get("launcher_identity")
                ) != canonical_bytes(expected_launcher_identity) or
                command.get("cwd") != str(source) or
                command.get("environment") != GIT_ENVIRONMENT or
                type(command.get("returncode")) is not int or
                command.get("returncode") != 0 or
                command.get("timed_out") is not False):
            raise CrossoverError(
                "controlled build command {} is invalid".format(index)
            )
        for stream in ("stdout", "stderr"):
            log_path = result_dir / command[stream + "_log"]
            value = read_result_regular(
                result_root, command[stream + "_log"],
                MAX_RETAINED_LOG_BYTES,
                "controlled build {} log".format(stream),
            )
            if digest_bytes(value) != command[stream + "_sha256"]:
                raise CrossoverError("controlled build log hash changed")


def validate_controlled_build(
        settings, jobs, manifest_path, result_root=None):
    if settings.get("mode") != "historical-avx2":
        return
    owns_root = result_root is None
    if result_root is None:
        result_root = open_existing_canonical_directory(
            Path(manifest_path).parent,
            "controlled-build validation result directory",
        )
    try:
        return validate_controlled_build_held(
            settings, jobs, manifest_path, result_root
        )
    finally:
        if owns_root:
            close_owned_directory(result_root)


def validate_manifest(manifest, path, result_root=None):
    expected_keys = {
        "configuration_id", "evidence_contract", "executables", "jobs",
        "machine", "schema", "settings", "source_fingerprint",
    }
    if (not isinstance(manifest, dict) or set(manifest) != expected_keys or
            manifest.get("schema") != SCHEMA):
        raise CrossoverError(
            "{} has an unknown, legacy, or incomplete schema".format(path)
        )
    settings = manifest.get("settings")
    jobs = manifest.get("jobs")
    executables = manifest.get("executables")
    source_state = manifest.get("source_fingerprint")
    if (not isinstance(settings, dict) or not isinstance(jobs, list) or
            not jobs or not isinstance(executables, dict) or not executables or
            not isinstance(source_state, dict)):
        raise CrossoverError("{} omits required v5 manifest fields".format(path))
    expected_job_keys = {
        "build_metadata", "cell", "configuration_id", "executable",
        "executable_artifact", "executable_sha256", "invocation_order",
        "job_id", "seed", "source_identity",
    }
    if any(not isinstance(job, dict) or set(job) != expected_job_keys
           for job in jobs):
        raise CrossoverError("{} contains an invalid v6 job".format(path))
    validate_manifest_settings(settings, jobs, path)
    validate_controlled_build(settings, jobs, path, result_root)
    contract = evidence_contract(settings.get("frozen_executable_required") is True)
    if canonical_bytes(manifest.get("evidence_contract")) != canonical_bytes(contract):
        raise CrossoverError("{} has a stale or relabeled evidence contract".format(path))
    validate_source_state(
        source_state, str(path) + " source_fingerprint",
        settings.get("mode") in AUTHORITATIVE_COMMANDS
    )
    if manifest.get("configuration_id") != digest_value(manifest_identity(manifest)):
        raise CrossoverError("{} configuration ID does not match its content".format(path))

    seen = set()
    expected_executables = {}
    for job in jobs:
        if not isinstance(job, dict) or set(job) != expected_job_keys:
            raise CrossoverError("{} contains an invalid v6 job".format(path))
        validate_cell_value(job.get("cell"), "manifest job")
        identity = job_identity(
            job["cell"], job["executable"], job["executable_artifact"],
            job["executable_sha256"], job["build_metadata"],
            job["source_identity"],
            manifest["machine"], settings
        )
        configuration_id = digest_value(identity)
        job_id = configuration_id[:24]
        if (job.get("configuration_id") != configuration_id or
                job.get("job_id") != job_id or job_id in seen or
                canonical_bytes(job.get("source_identity")) !=
                canonical_bytes(source_state) or
                job.get("seed") != stable_seed(job["cell"]) or
                job.get("invocation_order") != list(invocation_order(
                    settings.get("mode"), job_id,
                    settings.get("abba_rounds", 0)
                ))):
            raise CrossoverError("{} contains a stale or duplicate job".format(path))
        frozen_required = settings.get("frozen_executable_required") is True
        if frozen_required != (job.get("executable_artifact") is not None):
            raise CrossoverError("{} has inconsistent frozen job policy".format(path))
        if frozen_required:
            validate_frozen_executable(
                job["executable_artifact"], job["build_metadata"], source_state,
                result_root,
            )
        backend = job["cell"].get("backend")
        if backend not in expected_executables:
            expected_executables[backend] = {
                "build_metadata": job["build_metadata"],
                "execution_artifact": job["executable_artifact"],
                "origin_path": job["build_metadata"]["executable"]["path"],
                "origin_sha256": job["build_metadata"]["executable"]["sha256"],
            }
        elif canonical_bytes(expected_executables[backend]) != canonical_bytes({
                "build_metadata": job["build_metadata"],
                "execution_artifact": job["executable_artifact"],
                "origin_path": job["build_metadata"]["executable"]["path"],
                "origin_sha256": job["build_metadata"]["executable"]["sha256"],
        }):
            raise CrossoverError("{} has inconsistent backend jobs".format(path))
        seen.add(job_id)
    if canonical_bytes(executables) != canonical_bytes(expected_executables):
        raise CrossoverError("{} executable map differs from its jobs".format(path))
    return manifest


def load_manifest(result_dir, result_root=None):
    owns_root = result_root is None
    if result_root is None:
        result_root = open_existing_canonical_directory(
            result_dir, "manifest canonical result directory"
        )
    path = result_root["path"] / "manifest.json"
    try:
        manifest_bytes = read_result_regular(
            result_root, "manifest.json", MAX_RAW_JSON_BYTES,
            "campaign manifest",
        )
        manifest = decode_json_bytes(manifest_bytes, str(path))
        return validate_manifest(manifest, path, result_root)
    finally:
        if owns_root:
            close_owned_directory(result_root)


def require_compatible_result_dir(
        result_dir, manifest, result_root=None):
    owns_root = result_root is None
    if result_root is None:
        result_root = owned_canonical_directory(
            result_dir, "resume canonical result directory"
        )
    try:
        try:
            metadata = os.stat(
                "manifest.json", dir_fd=result_root["descriptor"],
                follow_symlinks=False,
            )
        except FileNotFoundError:
            try:
                names = list_result_regular_names(
                    result_root, "jobs", ".json", "resume job directory"
                )
            except FileNotFoundError:
                names = []
            if names:
                raise CrossoverError(
                    "result directory has jobs but no manifest: {}".format(
                        result_dir
                    )
                )
            return
        if (not stat.S_ISREG(metadata.st_mode) or
                metadata.st_uid != os.getuid() or
                metadata.st_nlink != 1 or
                stat.S_IMODE(metadata.st_mode) != 0o600):
            raise CrossoverError("retained manifest is unsafe")
        previous = load_manifest(result_dir, result_root)
        if previous.get("configuration_id") != manifest.get(
                "configuration_id"):
            raise CrossoverError(
                "result directory belongs to configuration {}; new "
                "configuration is {}; select a new --result-dir rather than "
                "mixing jobs".format(
                    previous.get("configuration_id"),
                    manifest.get("configuration_id")
                )
            )
    finally:
        if owns_root:
            close_owned_directory(result_root)


def load_job_results(result_dir, manifest, result_root=None):
    owns_root = result_root is None
    if result_root is None:
        result_root = open_existing_canonical_directory(
            result_dir, "job-results canonical result directory"
        )
    job_dir = result_root["path"] / "jobs"
    expected = {}
    for job in manifest["jobs"]:
        job_id = job.get("job_id")
        configuration_id = job.get("configuration_id")
        if not job_id or not configuration_id or job_id in expected:
            raise CrossoverError("manifest contains a duplicate or incomplete job")
        expected[job_id] = job
    try:
        actual_names = list_result_regular_names(
            result_root, "jobs", ".json", "job result directory"
        )
    except FileNotFoundError:
        if owns_root:
            close_owned_directory(result_root)
        raise CrossoverError(
            "job directory does not exist: {}".format(job_dir)
        )
    actual_paths = {
        Path(name).stem: job_dir / name for name in actual_names
    }
    missing = sorted(set(expected) - set(actual_paths))
    extra = sorted(set(actual_paths) - set(expected))
    if missing or extra:
        raise CrossoverError(
            "job set does not match manifest (missing: {}; stale/extra: {})".format(
                ",".join(missing) or "none", ",".join(extra) or "none"
            )
        )
    try:
        results = []
        for job_id in sorted(expected):
            path = actual_paths[job_id]
            item_bytes = read_result_regular(
                result_root, Path("jobs") / path.name,
                MAX_RAW_JSON_BYTES, "job result JSON",
            )
            item = decode_json_bytes(item_bytes, str(path))
            if not isinstance(item, dict) or item.get("schema") != JOB_SCHEMA:
                raise CrossoverError(
                    "{} has a legacy or unknown job schema".format(path)
                )
            if item.get("job_id") != job_id:
                raise CrossoverError(
                    "{} contains the wrong job ID".format(path)
                )
            if item.get("configuration_id") != expected[job_id][
                    "configuration_id"]:
                raise CrossoverError(
                    "{} belongs to a stale configuration".format(path)
                )
            validate_job_artifacts(
                result_dir, item, expected[job_id],
                manifest["settings"], result_root,
            )
            results.append(item)
        final_names = list_result_regular_names(
            result_root, "jobs", ".json", "job result directory"
        )
        if final_names != actual_names:
            raise CrossoverError(
                "job result inventory changed while artifacts were validated"
            )
        return results
    finally:
        if owns_root:
            close_owned_directory(result_root)


def write_merged(
        result_dir, manifest, results, source_end, promotion_percent,
        execution_input_error="", result_root=None):
    ordered = sorted(results, key=lambda item: item["job_id"])
    source_changed = source_end is not None and (
        canonical_bytes(source_end) !=
        canonical_bytes(manifest["source_fingerprint"])
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
    if result_root is None:
        atomic_write_json(result_dir / "matrix.json", merged)
        atomic_write_json(result_dir / "analysis.json", analysis)
    else:
        write_result_json(
            result_root, "matrix.json", merged, "merged matrix", replace=True
        )
        write_result_json(
            result_root, "analysis.json", analysis,
            "merged analysis", replace=True,
        )
    return merged


def invalidate_authoritative_result_dir_held(result_root, reason):
    """Ensure a failed authoritative guard can never leave retained PASS."""
    result_dir = result_root["path"]
    job_dir = result_dir / "jobs"
    retained = {}
    try:
        job_names = list_result_regular_names(
            result_root, "jobs", ".json", "invalidation job directory"
        )
    except FileNotFoundError:
        job_names = []
    for name in job_names:
        path = job_dir / name
        value = decode_json_bytes(
            read_result_regular(
                result_root, Path("jobs") / name,
                MAX_RAW_JSON_BYTES, "invalidation job JSON",
            ),
            str(path),
        )
        if (not isinstance(value, dict) or
                value.get("schema") != JOB_SCHEMA):
            continue
        if value.get("status") == "passed":
            value.pop("parity_identity", None)
            value.pop("summary", None)
            value["status"] = "failed"
            value["reason"] = (
                "authoritative guard validation failed: {}".format(reason)
            )
            write_result_json(
                result_root, Path("jobs") / name, value,
                "invalidated job result", replace=True,
            )
            validate_owned_directory(
                result_root, "canonical result directory"
            )
        retained[value.get("job_id")] = value

    matrix_path = result_dir / "matrix.json"
    manifest_path = result_dir / "manifest.json"
    try:
        matrix_metadata = os.stat(
            "matrix.json", dir_fd=result_root["descriptor"],
            follow_symlinks=False,
        )
    except FileNotFoundError:
        return
    if not stat.S_ISREG(matrix_metadata.st_mode):
        raise CrossoverError("invalidation matrix is not a regular file")
    matrix_bytes = read_result_regular(
        result_root, "matrix.json", MAX_RAW_JSON_BYTES,
        "invalidation matrix",
    )
    matrix = decode_json_bytes(matrix_bytes, str(matrix_path))
    try:
        manifest = decode_json_bytes(
            read_result_regular(
                result_root, "manifest.json", MAX_RAW_JSON_BYTES,
                "invalidation manifest",
            ),
            str(manifest_path),
        )
        ordered = [retained[job["job_id"]] for job in manifest["jobs"]]
        promotion = required_finite_number(
            matrix["analysis"]["promotion_percent"],
            "retained failure promotion_percent"
        )
        write_merged(
            result_dir, manifest, ordered,
            matrix.get("source_fingerprint_after"), promotion,
            matrix.get("execution_input_error") or "", result_root,
        )
        validate_owned_directory(
            result_root, "canonical result directory"
        )
    except Exception:
        # Even if older/partial artifacts cannot be recomputed, never leave a
        # human- or machine-visible PASS after this process observed guard loss.
        if isinstance(matrix, dict):
            matrix["status"] = "failed"
            analysis = matrix.get("analysis")
            if isinstance(analysis, dict):
                analysis["run_status"] = "failed"
            write_result_json(
                result_root, "matrix.json", matrix,
                "invalidated matrix", replace=True,
            )
            if isinstance(analysis, dict):
                write_result_json(
                    result_root, "analysis.json", analysis,
                    "invalidated analysis", replace=True,
                )


def invalidate_authoritative_result_dir(result_dir, reason):
    result_root = open_existing_canonical_directory(
        result_dir, "invalidation canonical result directory"
    )
    try:
        return invalidate_authoritative_result_dir_held(result_root, reason)
    finally:
        close_owned_directory(result_root)


def print_analysis(analysis):
    candidate = analysis["candidate"]
    print(
        "candidate cells: {}/{} passed; gain min={}; median={}; "
        "regressions={}; >= {:.1f}%={}; 95% lower min={}; "
        "confident >= threshold={}/{}; promotion gate={}".format(
            candidate["cell_count"] - candidate["failed_count"],
            candidate["cell_count"],
            ("n/a" if candidate["gain_min_percent"] is None else
             "{:.2f}%".format(candidate["gain_min_percent"])),
            ("n/a" if candidate["gain_median_percent"] is None else
             "{:.2f}%".format(candidate["gain_median_percent"])),
            candidate["regression_count"], analysis["promotion_percent"],
            candidate["promotion_count"],
            ("n/a" if candidate["confidence_worst_lower_gain_percent"] is None
             else "{:.2f}%".format(
                 candidate["confidence_worst_lower_gain_percent"])),
            candidate["confidence_promotion_count"],
            candidate["cell_count"],
            ("passed" if candidate[
                "all_cells_confidently_meet_promotion_threshold"]
             else "not passed"),
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


def validate_supported_run_mode(command, context):
    """Reject modes that cannot prove their claimed build closure."""
    if command == UNPROVED_PINNED_COMMAND:
        raise CrossoverError(
            "{}: generic pinned mode is disabled because an externally built "
            "bench_leopard2_direct_encode has no runner-owned clean-rebuild "
            "source/object/archive/link proof; use screen for diagnostics or "
            "historical-avx2 for authoritative evidence".format(context)
        )
    if command not in RUN_COMMANDS:
        raise CrossoverError("{} has an unknown runner mode".format(context))


def run_matrix(arguments):
    # This must precede topology, lock, source, and artifact handling.  In
    # particular, touching an executable linked from stale objects must never
    # turn the old generic pinned path into authoritative evidence.
    validate_supported_run_mode(arguments.command, "runner")
    if arguments.command not in AUTHORITATIVE_COMMANDS:
        return run_matrix_body(arguments, None)
    if arguments.cpu is None or arguments.sibling is None:
        raise CrossoverError(
            "authoritative v2 requires explicit --cpu and --sibling"
        )
    source = Path(arguments.source).resolve()
    support = load_isolation_support(source)
    try:
        allowed, housekeeping = support.validate_topology(
            arguments.cpu, arguments.sibling
        )
    except Exception as error:
        raise CrossoverError(
            "authoritative CPU topology validation failed: {}".format(error)
        )
    original_affinity = set(os.sched_getaffinity(0))
    if original_affinity != set(allowed):
        raise CrossoverError(
            "launch affinity changed during topology validation"
        )
    result_dir = resolve_path(source, arguments.result_dir)
    campaign_state = {"evidence_started": False}
    try:
        with canonical_authoritative_lock() as lock_guard:
            pair_guard = AuthoritativePairGuard(
                support.PairLease(arguments.cpu, arguments.sibling)
            )
            try:
                with pair_guard as pair_identity:
                    validate_authoritative_guards(lock_guard, pair_guard)
                    lock_identity = lock_guard.identity
                    os.sched_setaffinity(0, housekeeping)
                    isolation = {
                        "allowed": sorted(allowed),
                        "canonical_lock": lock_identity,
                        "canonical_guard": lock_guard,
                        "campaign_state": campaign_state,
                        "cpu": arguments.cpu,
                        "housekeeping": sorted(housekeeping),
                        "pair_guard": pair_guard,
                        "pair_lease": pair_identity,
                        "sibling": arguments.sibling,
                        "support": support,
                    }
                    result = run_matrix_body(arguments, isolation)
                    validate_authoritative_guards(lock_guard, pair_guard)
                    return result
            finally:
                os.sched_setaffinity(0, original_affinity)
    except Exception as error:
        if (isinstance(error, AuthoritativeGuardError) and
                campaign_state["evidence_started"]):
            try:
                invalidate_authoritative_result_dir(result_dir, error)
            except Exception as invalidation_error:
                raise CrossoverError(
                    "authoritative isolation failed: {}; retained evidence "
                    "invalidation also failed: {}".format(
                        error, invalidation_error
                    )
                )
        if isinstance(error, CrossoverError):
            raise
        raise CrossoverError(
            "authoritative isolation failed: {}".format(error)
        )


def run_matrix_body_held(arguments, isolation, result_root):
    # Defense in depth for tests or callers that invoke the body directly.
    validate_supported_run_mode(arguments.command, "runner")
    source = Path(arguments.source).resolve()
    executable_root = resolve_path(
        source, arguments.executable_root or arguments.build_root
    )
    result_dir = result_root["path"]
    backends = parse_backends(arguments.backends)
    authoritative = arguments.command in AUTHORITATIVE_COMMANDS

    def validate_campaign_guards():
        if isolation is not None:
            validate_authoritative_guards(
                isolation["canonical_guard"], isolation["pair_guard"]
            )

    validate_campaign_guards()
    if arguments.command == "historical-avx2" and backends != ["avx2"]:
        raise CrossoverError(
            "historical-avx2 requires exactly --backends avx2"
        )
    cpus = (
        list(isolation["allowed"]) if isolation is not None
        else allowed_cpus()
    )
    pin_cpu = None
    taskset = None
    if authoritative:
        pin_cpu = arguments.cpu
        if pin_cpu not in cpus:
            raise CrossoverError(
                "pinned CPU {} is outside allowed affinity {}".format(
                    pin_cpu, compact_cpu_list(cpus)
                )
            )
        taskset = shutil.which(
            arguments.taskset, path=BENCHMARK_ENVIRONMENT["PATH"]
        )
        if not taskset or Path(taskset).resolve() != Path("/usr/bin/taskset"):
            raise CrossoverError(
                "authoritative mode requires canonical /usr/bin/taskset"
            )
        if arguments.workers != 1:
            raise CrossoverError("pinned ABBA measurements require --workers 1")
    source_state = source_identity(source, require_clean=authoritative)
    machine = machine_identity(cpus)
    if pin_cpu is not None:
        machine["pinned_cpu_topology"] = cpu_topology(pin_cpu)
        machine["reserved_sibling_topology"] = cpu_topology(arguments.sibling)
    controlled_build = None
    if arguments.command == "historical-avx2":
        validate_campaign_guards()
        executable, controlled_build = controlled_avx2_build(
            source, result_dir, source_state,
            min(len(isolation["housekeeping"]), 128),
            validate_campaign_guards,
            {
                "canonical_lock": isolation["canonical_lock"],
                "pair_lease": isolation["pair_lease"],
            },
            isolation["canonical_guard"].descriptor,
            result_root,
        )
        validate_campaign_guards()
        executables = {"avx2": executable}
    else:
        executables = {
            backend: find_executable(executable_root, backend)
            for backend in backends
        }
    build_metadata = {
        backend: cmake_build_metadata(executable)
        for backend, executable in executables.items()
    }
    for backend in backends:
        validate_build_source_binding(
            build_metadata[backend], source, source_state, backend,
            require_fresh=authoritative
        )
    settings = {
        "abba_rounds": arguments.abba_rounds if authoritative else 0,
        "benchmark": {
            "batch": arguments.batch,
            "iterations": arguments.iterations,
            "reuse": arguments.reuse,
            "warmups": arguments.warmups,
        },
        "mode": arguments.command,
        "frozen_executable_required": authoritative,
        "pin_cpu": pin_cpu,
        "placement_policy": (
            "external taskset single-CPU pin"
            if authoritative else
            "unpinned discovery jobs inherit allowed affinity; OS schedules workers"
        ),
        "taskset": str(Path(taskset).resolve()) if taskset else None,
        "taskset_sha256": (
            digest_bytes(Path(taskset).read_bytes()) if taskset else None
        ),
        "timeout_seconds_per_invocation": arguments.timeout,
        "workers": arguments.workers,
    }
    if authoritative:
        settings["isolation"] = {
            "allowed_cpu_set_at_launch": isolation["allowed"],
            "benchmark_cpu": isolation["cpu"],
            "canonical_lock": isolation["canonical_lock"],
            "child_environment": dict(BENCHMARK_ENVIRONMENT),
            "housekeeping_cpu_set": isolation["housekeeping"],
            "pair_lease": isolation["pair_lease"],
            "reserved_sibling": isolation["sibling"],
        }
    if arguments.command == "historical-avx2":
        cells = historical_avx2_grid()
        settings["campaign"] = {
            "cell_count": len(cells),
            "cells_sha256": digest_value(cells),
            "name": "corrected-high-profile-explicit-avx2",
        }
        settings["controlled_build"] = controlled_build
    else:
        r_values = parse_r_values(arguments.r)
        cells = sorted_grid(backends, r_values, arguments.full)

    executable_artifacts = {}
    execution_executables = dict(executables)
    if authoritative:
        for backend in backends:
            validate_campaign_guards()
            artifact = freeze_executable(
                result_dir, backend, executables[backend],
                build_metadata[backend], source_state, result_root
            )
            validate_campaign_guards()
            executable_artifacts[backend] = artifact
            execution_executables[backend] = Path(artifact["executable"])
        source_after_freeze = source_identity(source, require_clean=True)
        if canonical_bytes(source_after_freeze) != canonical_bytes(source_state):
            raise CrossoverError("source identity changed while freezing executables")

    jobs = make_jobs(
        cells, execution_executables, build_metadata, source_state, machine,
        settings, executable_artifacts
    )
    contract = evidence_contract(authoritative)
    executable_manifest = {
        backend: {
            "build_metadata": build_metadata[backend],
            "execution_artifact": executable_artifacts.get(backend),
            "origin_path": str(path),
            "origin_sha256": digest_bytes(path.read_bytes()),
        } for backend, path in sorted(executables.items())
    }
    manifest = {
        "configuration_id": None,
        "evidence_contract": contract,
        "executables": executable_manifest,
        "jobs": jobs,
        "machine": machine,
        "schema": SCHEMA,
        "settings": settings,
        "source_fingerprint": source_state,
    }
    manifest["configuration_id"] = digest_value(manifest_identity(manifest))
    validate_manifest(
        manifest, result_dir / "manifest.json", result_root
    )
    require_compatible_result_dir(result_dir, manifest, result_root)
    validate_campaign_guards()
    validate_owned_directory(result_root, "canonical result directory")
    try:
        os.stat(
            "manifest.json", dir_fd=result_root["descriptor"],
            follow_symlinks=False,
        )
    except FileNotFoundError:
        write_result_json(
            result_root, "manifest.json", manifest,
            "campaign manifest", replace=False,
        )
    else:
        # require_compatible_result_dir already parsed and validated this exact
        # retained manifest.  Keeping it avoids a crash window in which jobs
        # could outlive a needlessly replaced manifest.
        load_manifest(result_dir, result_root)
    validate_owned_directory(result_root, "canonical result directory")
    if isolation is not None:
        isolation["campaign_state"]["evidence_started"] = True
    validate_campaign_guards()
    print(
        "direct-encode crossover: mode={} cells={} backends={} cpus={}{}".format(
            arguments.command, len(jobs), ",".join(backends),
            machine["allowed_cpu_list"],
            " pinned={}".format(pin_cpu) if pin_cpu is not None else "",
        ), flush=True
    )
    context = {
        "isolation": isolation,
        "result_dir": result_dir,
        "result_root": result_root,
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
    source_end = source_identity(source, require_clean=False)
    validate_campaign_guards()
    execution_input_error = ""
    try:
        checked = set()
        for job in jobs:
            key = (job["executable"], job["executable_sha256"])
            if key not in checked:
                checked.add(key)
                validate_execution_inputs(job, result_root)
    except CrossoverError as error:
        execution_input_error = str(error)
    validate_campaign_guards()
    merged = write_merged(
        result_dir, manifest, results, source_end, arguments.promotion_percent,
        execution_input_error, result_root
    )
    validate_owned_directory(result_root, "canonical result directory")
    validate_campaign_guards()
    print_analysis(merged["analysis"])
    promotion_passed = merged["analysis"]["candidate"][
        "all_cells_confidently_meet_promotion_threshold"
    ]
    print(
        "matrix integrity: {}; promotion gate: {} ({})".format(
            merged["status"], "passed" if promotion_passed else "not passed",
            result_dir / "matrix.json"
        )
    )
    if merged["status"] != "passed":
        return 1
    if authoritative and not promotion_passed:
        return 2
    return 0


def run_matrix_body(arguments, isolation):
    source = Path(arguments.source).resolve()
    result_dir = resolve_path(source, arguments.result_dir)
    result_root = owned_canonical_directory(
        result_dir, "run canonical result directory"
    )
    try:
        return run_matrix_body_held(arguments, isolation, result_root)
    finally:
        close_owned_directory(result_root)


def analyze_command_held(arguments, result_root):
    result_dir = result_root["path"]
    manifest = load_manifest(result_dir, result_root)
    results = load_job_results(result_dir, manifest, result_root)
    matrix_path = result_dir / "matrix.json"
    try:
        matrix = decode_json_bytes(
            read_result_regular(
                result_root, "matrix.json", MAX_RAW_JSON_BYTES,
                "completed matrix",
            ),
            str(matrix_path),
        )
    except (CrossoverError, OSError) as error:
        raise CrossoverError("cannot read completed matrix {}: {}".format(
            matrix_path, error
        ))
    matrix_keys = {
        "analysis", "execution_input_error", "jobs",
        "manifest_configuration_id", "schema", "source_changed_during_run",
        "source_fingerprint", "source_fingerprint_after", "status",
    }
    if (not isinstance(matrix, dict) or set(matrix) != matrix_keys or
            matrix.get("schema") != SCHEMA or
            matrix.get("manifest_configuration_id") != manifest["configuration_id"]):
        raise CrossoverError("matrix does not match the current manifest")
    if canonical_bytes(matrix.get("source_fingerprint")) != canonical_bytes(
            manifest.get("source_fingerprint")):
        raise CrossoverError("matrix source fingerprint does not match the manifest")
    source_after = validate_source_state(
        matrix.get("source_fingerprint_after"),
        str(matrix_path) + " source_fingerprint_after", False
    )
    matrix_jobs_value = matrix.get("jobs")
    if (not isinstance(matrix_jobs_value, list) or
            any(not isinstance(item, dict) for item in matrix_jobs_value)):
        raise CrossoverError("matrix jobs snapshot has an invalid shape")
    matrix_jobs = sorted(
        matrix_jobs_value, key=lambda item: item.get("job_id", "")
    )
    ordered_results = sorted(results, key=lambda item: item.get("job_id", ""))
    if canonical_bytes(matrix_jobs) != canonical_bytes(ordered_results):
        raise CrossoverError("matrix job snapshot differs from the validated job files")
    source_changed = canonical_bytes(source_after) != canonical_bytes(
        manifest["source_fingerprint"]
    )
    if matrix.get("source_changed_during_run") is not source_changed:
        raise CrossoverError("matrix source-change flag was not derived")
    execution_input_error = ""
    try:
        checked = set()
        for job in manifest["jobs"]:
            key = (job["executable"], job["executable_sha256"])
            if key not in checked:
                checked.add(key)
                validate_execution_inputs(job, result_root)
    except CrossoverError as error:
        execution_input_error = str(error)
    retained_input_error = matrix.get("execution_input_error")
    if retained_input_error != (execution_input_error or None):
        raise CrossoverError("matrix execution-input error was not derived")
    failed = any(item.get("status") != "passed" for item in results)
    run_status = (
        "failed"
        if failed or source_changed or execution_input_error
        else "passed"
    )
    if matrix.get("status") != run_status:
        raise CrossoverError("matrix status was not derived from retained evidence")
    retained_analysis = required_mapping(matrix.get("analysis"), "matrix.analysis")
    retained_promotion = required_finite_number(
        retained_analysis.get("promotion_percent"),
        "matrix.analysis.promotion_percent"
    )
    if retained_promotion < 0:
        raise CrossoverError("matrix retained a negative promotion threshold")
    recomputed_retained_analysis = analyze_results(
        results, retained_promotion, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    if canonical_bytes(retained_analysis) != canonical_bytes(
            recomputed_retained_analysis):
        raise CrossoverError("matrix analysis was not derived from retained jobs")
    analysis = analyze_results(
        results, arguments.promotion_percent, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    output = Path(arguments.output).resolve() if arguments.output else result_dir / "analysis.json"
    validate_owned_directory(result_root, "canonical result directory")
    try:
        output_relative = output.relative_to(result_dir)
    except ValueError:
        atomic_write_json(output, analysis)
    else:
        write_result_json(
            result_root, output_relative, analysis,
            "reanalyzed output", replace=True,
        )
    validate_owned_directory(result_root, "canonical result directory")
    print_analysis(analysis)
    print("analysis: {}".format(output))
    if analysis["jobs_failed"] != 0 or run_status != "passed":
        return 1
    if (manifest["settings"]["mode"] in AUTHORITATIVE_COMMANDS and
            not analysis["candidate"][
                "all_cells_confidently_meet_promotion_threshold"]):
        return 2
    return 0


def analyze_command(arguments):
    result_root = open_existing_canonical_directory(
        Path(arguments.result_dir),
        "analyze canonical result directory",
    )
    try:
        return analyze_command_held(arguments, result_root)
    finally:
        close_owned_directory(result_root)


def self_test():
    def check(condition, message):
        if not condition:
            raise CrossoverError("self-test failed: {}".format(message))

    class InvalidFilesystemPath(object):
        def __init__(self, error):
            self.error = error

        def __fspath__(self):
            raise self.error

    for label, invalid_path in (
            ("NUL", "/self-test/\0path"),
            ("ValueError", InvalidFilesystemPath(ValueError("invalid"))),
            ("OSError", InvalidFilesystemPath(OSError("unavailable")))):
        try:
            checked_resolve(invalid_path, "self-test malformed path")
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} malformed path was accepted".format(
                    label
                )
            )

    check(compact_cpu_list([3, 2, 1, 7, 9, 8]) == "1-3,7-9",
          "CPU-list compaction")
    check(digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2}),
          "canonical JSON digest")
    check(stable_seed({"cell": 1}) == stable_seed({"cell": 1}),
          "stable cell seed")
    check(parse_backends("avx512") == ["avx512"], "AVX-512 backend parsing")
    check(parser().parse_args(["screen"]).backends == "scalar,ssse3,avx2",
          "default backend list")
    repository = Path(__file__).resolve().parents[1]
    repository_state = source_identity(repository, require_clean=False)
    validate_source_state(repository_state, "self-test repository", False)
    with tempfile.TemporaryDirectory(
            prefix="leo2-crossover-path-shim-") as shim_directory:
        shim = Path(shim_directory) / "git"
        shim.write_text("#!/bin/sh\nexit 97\n", encoding="utf-8")
        shim.chmod(0o755)
        old_path = os.environ.get("PATH")
        os.environ["PATH"] = shim_directory
        try:
            shim_state = source_identity(repository, require_clean=False)
        finally:
            if old_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = old_path
        check(
            canonical_bytes(shim_state["git_tool"]) ==
            canonical_bytes(repository_state["git_tool"]) and
            shim_state["git_tool"]["path"] == str(CANONICAL_GIT),
            "source identity ignored the canonical absolute Git executable"
        )
    repository_gitlink = repository_state["files"].get("sse2neon")
    check(
        isinstance(repository_gitlink, dict) and
        repository_gitlink.get("index_mode") == "160000" and
        repository_gitlink.get("submodule", {}).get("head") ==
            repository_gitlink.get("index_object"),
        "real tracked sse2neon gitlink closure"
    )
    grid = sorted_grid(["scalar", "ssse3", "avx2"], [16], False)
    candidates = [item for item in grid if item["region"] == "candidate"]
    check(bool(candidates), "compact grid has candidate cells")
    for item in candidates:
        check(item["profile"] == "low" and item["q"] == 1,
              "candidate profile and Q")
        check(item["shard_bytes"] >= 1024 and item["shard_bytes"] % 64 == 0,
              "candidate byte region")
        check(item["k"] >= threshold_k(item["backend"]),
              "candidate backend K threshold")
    regions = {item["region"] for item in grid}
    check({
        "excluded_scalar_k2", "excluded_q2", "excluded_high_profile",
        "excluded_ragged_tail",
    }.issubset(regions), "excluded-neighbor regions")
    check({item["field"] for item in candidates} == {"gf8", "gf16"},
          "candidate field coverage")
    check(any(item["region"] == "candidate_sparse_output" for item in grid),
          "sparse-output candidate")
    check(any(item["region"] == "excluded_q2_holey" for item in grid),
          "holey-Q2 neighbor")
    r1 = sorted_grid(["scalar"], [1], False)
    check(all(item["r"] == 1 and item["q"] == 1 for item in r1),
          "R=1 grid bounds Q")
    check(parse_r_values("1,16,1") == [1, 16], "R-list normalization")
    check(parse_r_values("all") == list(range(1, 17)), "all-R expansion")
    check(requested_indices("0,2-3", 4) == [0, 2, 3],
          "parity-mask expansion")
    check(invocation_order("historical-avx2", "00000000", 2) == (
        "direct", "transform", "transform", "direct",
        "direct", "transform", "transform", "direct",
    ), "historical authoritative ABBA order")
    with tempfile.TemporaryDirectory(
            prefix="leo2-crossover-stale-pinned-") as stale_directory:
        stale_root = Path(stale_directory) / "avx2"
        stale_object = (
            stale_root / "CMakeFiles" / "leopard_test_hooks.dir" /
            "leopard2.cpp.o"
        )
        stale_archive = stale_root / "libleopard_test_hooks.a"
        touched_executable = stale_root / "bench_leopard2_direct_encode"
        stale_object.parent.mkdir(parents=True)
        stale_object.write_bytes(b"stale object")
        stale_archive.write_bytes(b"archive containing stale object")
        touched_executable.write_bytes(b"recently relinked-looking executable")
        touched_executable.chmod(0o755)
        now_ns = time.time_ns()
        stale_ns = now_ns - 10_000_000_000
        os.utime(stale_object, ns=(stale_ns, stale_ns))
        os.utime(stale_archive, ns=(stale_ns + 1, stale_ns + 1))
        os.utime(touched_executable, ns=(now_ns, now_ns))
        check(
            touched_executable.stat().st_mtime_ns >
                stale_archive.stat().st_mtime_ns,
            "adversarial executable appears newer than stale link inputs"
        )
        stale_pinned_arguments = parser().parse_args([
            "pinned",
            "--source", str(repository),
            "--executable-root", str(Path(stale_directory)),
            "--backends", "avx2",
            "--cpu", "0",
            "--sibling", "1",
        ])
        try:
            run_matrix(stale_pinned_arguments)
        except CrossoverError as error:
            rejection = str(error)
            check(
                "generic pinned mode is disabled" in rejection and
                "source/object/archive/link proof" in rejection,
                "stale/touched pinned artifact is rejected for unproved closure"
            )
        else:
            raise CrossoverError(
                "self-test failed: stale/touched pinned artifact was accepted"
            )
        try:
            validate_manifest_settings(
                {"mode": "pinned"}, [], Path(stale_directory) / "manifest.json"
            )
        except CrossoverError as error:
            check(
                "generic pinned mode is disabled" in str(error),
                "retained pinned manifest is rejected for unproved closure"
            )
        else:
            raise CrossoverError(
                "self-test failed: retained pinned manifest was accepted"
            )
    mode_zero_configure = controlled_avx2_configure_argv(
        "/usr/bin/cmake", "/usr/bin/ninja",
        "/self-test/source", "/self-test/build"
    )
    check(
        not any(
            argument.startswith(
                "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE="
            )
            for argument in mode_zero_configure
        ),
        "historical mode-0 configure command remains implicit"
    )
    check(
        mode_zero_configure[0] == "/usr/bin/cmake" and
        "-DCMAKE_MAKE_PROGRAM=/usr/bin/ninja" in mode_zero_configure,
        "controlled configure binds the selected CMake and Ninja launchers"
    )

    historical = historical_avx2_grid()
    exact = [
        item for item in historical
        if item["region"] == "candidate_historical_exact"
    ]
    check(len(exact) == 1, "one exact historical cell")
    check(exact[0] == cell(
        "candidate_historical_exact", "avx2", 2, 16, "high",
        "gf8", 4096, 1, "0"
    ), "exact K=2,R=16,Q=1,4096-byte historical cell")
    check({
        "historical_neighbor_k", "historical_neighbor_r",
        "historical_neighbor_bytes_lower", "historical_neighbor_bytes_upper",
        "historical_neighbor_sparse_output", "historical_neighbor_q2",
    }.issubset({item["region"] for item in historical}),
          "historical one-coordinate neighbors")

    self_test_build_root = Path("/self-test/build")
    self_test_effective_entries = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
        "CMAKE_GENERATOR": "Ninja",
        "CMAKE_CXX_COMPILER": "/self-test/c++",
        # These deliberately differ from the stale cache-level values below.
        "CMAKE_CXX_FLAGS": " -Wall -Wextra -fopenmp",
        "CMAKE_CXX_FLAGS_DEBUG": "-g -O0",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
        "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
        "ENABLE_OPENMP": "ON",
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
        "LEO2_BACKEND_VARIANT": "avx2",
        "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_TESTS": "ON",
        **BUILD_CONFIGURATION_CANONICAL_SELECTORS,
    }
    self_test_digest = build_configuration_digest(
        self_test_effective_entries
    )
    for selector in (
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL"):
        changed_entries = dict(self_test_effective_entries)
        changed_entries[selector] = (
            "OFF" if changed_entries[selector] == "ON" else "ON")
        check(
            build_configuration_digest(changed_entries) != self_test_digest,
            "{} changes the effective-configuration digest".format(selector)
        )
    with tempfile.TemporaryDirectory(
            prefix="leo2-build-config-self-test-") as configuration_dir:
        configuration_path = Path(configuration_dir) / "configuration.txt"
        configuration_material = "".join(
            "{}={}\n".format(
                variable, self_test_effective_entries[variable]
            )
            for variable in BUILD_CONFIGURATION_VARIABLES
        )
        valid_configuration = (
            "schema={}\nsha256={}\n{}".format(
                BUILD_CONFIGURATION_FILE_SCHEMA,
                self_test_digest,
                configuration_material
            )
        )
        configuration_path.write_text(
            valid_configuration, encoding="utf-8"
        )
        configuration_attestation = (
            read_build_configuration_attestation(configuration_path)
        )
        check(
            configuration_attestation["sha256"] == self_test_digest and
            configuration_attestation["entries"] ==
                self_test_effective_entries,
            "effective configuration comes from the exact CMake sidecar"
        )
        mode_digests = {}
        for mode in ("0", "1", "2"):
            mode_entries = dict(self_test_effective_entries)
            mode_entries["LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] = mode
            mode_digest = build_configuration_digest(mode_entries)
            mode_material = "".join(
                "{}={}\n".format(variable, mode_entries[variable])
                for variable in BUILD_CONFIGURATION_VARIABLES
            )
            configuration_path.write_text(
                "schema={}\nsha256={}\n{}".format(
                    BUILD_CONFIGURATION_FILE_SCHEMA,
                    mode_digest,
                    mode_material
                ),
                encoding="utf-8"
            )
            mode_attestation = read_build_configuration_attestation(
                configuration_path
            )
            check(
                mode_attestation["entries"][
                    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] == mode and
                mode_attestation["sha256"] == mode_digest,
                "small-direct mode {} sidecar round-trip".format(mode)
            )
            mode_digests[mode] = mode_digest
        check(
            len(set(mode_digests.values())) == 3,
            "small-direct mode sidecars have distinct digests"
        )
        historical_entries = {
            variable: self_test_effective_entries[variable]
            for variable in BUILD_CONFIGURATION_VARIABLES_V2
        }
        historical_digest = build_configuration_digest(
            historical_entries, BUILD_CONFIGURATION_VARIABLES_V2)
        historical_material = "".join(
            "{}={}\n".format(variable, historical_entries[variable])
            for variable in BUILD_CONFIGURATION_VARIABLES_V2)
        configuration_path.write_text(
            "schema={}\nsha256={}\n{}".format(
                BUILD_CONFIGURATION_FILE_SCHEMA_V2,
                historical_digest,
                historical_material),
            encoding="utf-8")
        historical_attestation = read_build_configuration_attestation(
            configuration_path)
        check(
            historical_attestation == {
                "entries": historical_entries,
                "path": str(configuration_path.resolve()),
                "schema": BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
                "sha256": historical_digest,
            },
            "historical v2 sidecar retains its exact selector set")
        try:
            validate_build_configuration_attestation(
                historical_attestation, configuration_path)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: current evidence accepted a historical "
                "v2 effective-configuration attestation")
        check(
            validate_build_configuration_attestation(
                historical_attestation,
                configuration_path,
                BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
            ) == historical_digest,
            "historical v2 attestation requires an explicit old contract")
        historical_v3_entries = {
            variable: self_test_effective_entries[variable]
            for variable in BUILD_CONFIGURATION_VARIABLES_V3
        }
        historical_v3_digest = build_configuration_digest(
            historical_v3_entries, BUILD_CONFIGURATION_VARIABLES_V3)
        historical_v3_material = "".join(
            "{}={}\n".format(variable, historical_v3_entries[variable])
            for variable in BUILD_CONFIGURATION_VARIABLES_V3)
        configuration_path.write_text(
            "schema={}\nsha256={}\n{}".format(
                BUILD_CONFIGURATION_FILE_SCHEMA_V3,
                historical_v3_digest,
                historical_v3_material),
            encoding="utf-8")
        historical_v3_attestation = (
            read_build_configuration_attestation(configuration_path))
        check(
            validate_build_configuration_attestation(
                historical_v3_attestation,
                configuration_path,
                BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V3,
            ) == historical_v3_digest,
            "historical v3 attestation retains its General=OFF contract")
        try:
            validate_build_configuration_attestation(
                historical_v3_attestation, configuration_path)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: current evidence accepted a historical "
                "v3 effective-configuration attestation")
        current_enabled = {
            "entries": dict(self_test_effective_entries),
            "path": str(configuration_path.resolve()),
            "schema": BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
            "sha256": "",
        }
        current_enabled["entries"][
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = "ON"
        current_enabled["sha256"] = build_configuration_digest(
            current_enabled["entries"])
        try:
            validate_build_configuration_attestation(
                current_enabled, configuration_path)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: current effective configuration accepted "
                "GENERAL_ONE_LOSS_DIRECT=ON")
        historical_enabled = {
            "entries": dict(historical_entries),
            "path": str(configuration_path.resolve()),
            "schema": BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
            "sha256": "",
        }
        historical_enabled["entries"][
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE"] = "ON"
        historical_enabled["sha256"] = build_configuration_digest(
            historical_enabled["entries"],
            BUILD_CONFIGURATION_VARIABLES_V2)
        try:
            validate_build_configuration_attestation(
                historical_enabled, configuration_path,
                BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: historical effective configuration "
                "accepted HIGH_DIRECT_ENCODE=ON")
        separator_entries = dict(self_test_effective_entries)
        separator_entries["CMAKE_CXX_FLAGS"] += "\u2028preserved"
        try:
            build_configuration_digest(separator_entries)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: U+2028 line separator was blessed in "
                "attested text"
            )
        for category, character in (
                ("Cc", "\t"), ("Cf", "\u2066"),
                ("Zl", "\u2028"), ("Zp", "\u2029")):
            try:
                attested_text(
                    "prefix" + character + "suffix",
                    "self-test Unicode category " + category,
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: Unicode category {} was accepted in "
                    "attested text".format(category)
                )

        def reject_configuration_file(name, contents):
            configuration_path.write_text(contents, encoding="utf-8")
            try:
                read_build_configuration_attestation(configuration_path)
            except CrossoverError:
                return
            raise CrossoverError(
                "self-test failed: configuration sidecar mutation was "
                "accepted: {}".format(name)
            )

        reject_configuration_file(
            "wrong schema",
            valid_configuration.replace(
                BUILD_CONFIGURATION_FILE_SCHEMA, "unknown", 1
            )
        )
        reject_configuration_file(
            "malformed digest",
            valid_configuration.replace(self_test_digest, "not-a-digest", 1)
        )
        reject_configuration_file(
            "material mismatch",
            valid_configuration.replace(
                "CMAKE_CXX_FLAGS= -Wall",
                "CMAKE_CXX_FLAGS= -Werror", 1
            )
        )
        reject_configuration_file(
            "missing variable",
            valid_configuration.replace(
                "LEO2_BUILD_TESTS=ON\n", "", 1
            )
        )
        reject_configuration_file(
            "missing small-direct mode",
            valid_configuration.replace(
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\n", "", 1
            )
        )
        reject_configuration_file(
            "missing current general one-loss selector",
            valid_configuration.replace(
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\n", "", 1
            )
        )
        reject_configuration_file(
            "duplicate small-direct mode",
            valid_configuration +
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\n"
        )
        reject_configuration_file(
            "forged small-direct mode",
            valid_configuration.replace(
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\n",
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=3\n", 1
            )
        )
        reject_configuration_file(
            "extra variable",
            valid_configuration + "UNVERSIONED=1\n"
        )
        reject_configuration_file(
            "missing final newline",
            valid_configuration[:-1]
        )
        reject_configuration_file(
            "embedded NUL",
            valid_configuration.replace(
                "CMAKE_CXX_FLAGS= -Wall",
                "CMAKE_CXX_FLAGS= -Wall\0", 1
            )
        )
        reject_configuration_file(
            "embedded multiline value",
            valid_configuration.replace(
                "CMAKE_CXX_FLAGS= -Wall",
                "CMAKE_CXX_FLAGS= -Wall\ninjected", 1
            )
        )
        reject_configuration_file(
            "embedded carriage return",
            valid_configuration.replace(
                "CMAKE_CXX_FLAGS= -Wall",
                "CMAKE_CXX_FLAGS= -Wall\rinjected", 1
            )
        )

        for label, delimiter in (
                ("NUL", "\0"), ("LF", "\n"), ("CR", "\r"),
                ("lone surrogate", "\ud800")):
            invalid_entries = dict(self_test_effective_entries)
            invalid_entries["CMAKE_CXX_FLAGS"] += delimiter + "injected"
            try:
                build_configuration_digest(invalid_entries)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} in-memory attestation value was "
                    "accepted".format(label)
                )
        for invalid_mode in ("", "3", "01", "ON", "1\n2", "1\r2"):
            invalid_entries = dict(self_test_effective_entries)
            invalid_entries[
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] = invalid_mode
            try:
                build_configuration_digest(invalid_entries)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: forged small-direct mode {!r} was "
                    "accepted".format(invalid_mode)
                )

        for label, invalid_path in (
                ("NUL sidecar path", "/self-test/\0configuration.txt"),
                ("NUL expected path", "/self-test/\0expected.txt")):
            try:
                if label == "NUL sidecar path":
                    read_build_configuration_attestation(invalid_path)
                else:
                    validate_build_configuration_attestation(
                        configuration_attestation, invalid_path
                    )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} was accepted".format(label)
                )

        invalid_path_attestation = dict(configuration_attestation)
        invalid_path_attestation["path"] = "/self-test/\0configuration.txt"
        try:
            validate_build_configuration_attestation(
                invalid_path_attestation
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: NUL attestation path was accepted"
            )

    cache_prefixes = (
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER",
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
        "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
    )
    valid_cache = (
        "CMAKE_BUILD_TYPE:STRING=Release\n"
        "CMAKE_CXX_COMPILER:STRING=/usr/bin/clang++-18\n"
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN:BOOL=OFF\n"
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE:BOOL=OFF\n"
        "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT:BOOL=OFF\n"
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE:STRING=0\n"
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:INTERNAL={}\n"
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:INTERNAL={}\n"
    ).format(BUILD_CONFIGURATION_FILE_SCHEMA, self_test_digest).encode("utf-8")
    parsed_cache = parse_selected_cmake_cache(
        valid_cache, cache_prefixes, "self-test CMake cache"
    )
    check(
        parsed_cache["CMAKE_BUILD_TYPE"] == "Release" and
        parsed_cache["CMAKE_CXX_COMPILER"] == "/usr/bin/clang++-18" and
        parsed_cache["LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN"] == "OFF" and
        parsed_cache["LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE"] == "OFF" and
        parsed_cache["LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] == "OFF" and
        parsed_cache["LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] == "0" and
        parsed_cache[
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"
        ] == self_test_digest,
        "selected CMake cache parsing"
    )
    for compiler_type in ("FILEPATH", "UNINITIALIZED"):
        variant_cache = valid_cache.replace(
            b"CMAKE_CXX_COMPILER:STRING=",
            ("CMAKE_CXX_COMPILER:{}=".format(
                compiler_type
            )).encode("ascii"),
        )
        check(
            parse_selected_cmake_cache(
                variant_cache, cache_prefixes, "self-test CMake cache"
            )["CMAKE_CXX_COMPILER"] == "/usr/bin/clang++-18",
            "caller-selected compiler cache type {}".format(compiler_type)
        )

    for label, mutated_cache in (
            (
                "duplicate selected CMake cache key",
                valid_cache + b"CMAKE_BUILD_TYPE:STRING=Debug\n",
            ),
            (
                "missing selected CMake cache entry type",
                valid_cache.replace(
                    b"CMAKE_BUILD_TYPE:STRING=",
                    b"CMAKE_BUILD_TYPE=",
                ),
            ),
            (
                "wrong attestation-binding cache entry type",
                valid_cache.replace(
                    b"LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:INTERNAL=",
                    b"LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:STRING=",
                ),
            ),
            (
                "unknown selected CMake cache entry type",
                valid_cache.replace(
                    b"CMAKE_BUILD_TYPE:STRING=",
                    b"CMAKE_BUILD_TYPE:UNVERSIONED=",
                ),
            ),
            (
                "NUL selected CMake cache value",
                valid_cache.replace(b"Release\n", b"Release\0injected\n"),
            )):
        try:
            parse_selected_cmake_cache(
                mutated_cache, cache_prefixes, "self-test CMake cache"
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} was accepted".format(label)
            )

    single_configuration = dict(self_test_effective_entries)
    multi_configuration = dict(self_test_effective_entries)
    multi_configuration["CMAKE_BUILD_TYPE"] = "Release"
    multi_configuration["CMAKE_GENERATOR"] = "Ninja Multi-Config"
    validate_embedded_build_type(
        single_configuration, "Release", authoritative=False
    )
    validate_embedded_build_type(
        multi_configuration, "Debug", authoritative=False
    )
    validate_embedded_build_type(
        multi_configuration, "Release", authoritative=False
    )
    validate_embedded_build_type(
        multi_configuration, "Release", authoritative=True
    )
    for label, configuration, build_type, authoritative in (
            (
                "single-config Debug",
                single_configuration, "Debug", False,
            ),
            (
                "multi-config unknown configuration",
                multi_configuration, "RelWithDebInfo", False,
            ),
            (
                "authoritative multi-config Debug",
                multi_configuration, "Debug", True,
            ),
            (
                "authoritative single-config Debug",
                dict(single_configuration, CMAKE_BUILD_TYPE="Debug"),
                "Debug", True,
            ),
            (
                "multi-config without configuration types",
                dict(
                    multi_configuration,
                    CMAKE_CONFIGURATION_TYPES="",
                ),
                "Release", False,
            )):
        try:
            validate_embedded_build_type(
                configuration, build_type, authoritative
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} build type was accepted".format(label)
            )

    raw_job = {
        "build_metadata": {
            "build_root": str(self_test_build_root),
            "entries": {
                "CMAKE_BUILD_TYPE": "",
                "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
                "CMAKE_CXX_COMPILER": "/self-test/c++",
                "CMAKE_CXX_FLAGS": "",
                "CMAKE_CXX_FLAGS_DEBUG": "-g",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
                "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
                "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
                "CMAKE_GENERATOR": "Ninja",
                "ENABLE_OPENMP": "ON",
                "LEOPARD_ENABLE_GF8": "ON",
                "LEOPARD_ENABLE_GF16": "ON",
                "LEO2_BACKEND_VARIANT": "avx2",
                "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
                "LEO2_BUILD_BENCHMARKS": "ON",
                "LEO2_BUILD_TESTS": "ON",
                **BUILD_CONFIGURATION_CANONICAL_SELECTORS,
                "CMAKE_HOME_DIRECTORY": "/self-test/source",
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    BUILD_CONFIGURATION_FILE_SCHEMA,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                    self_test_digest,
            },
            "effective_configuration_attestation": {
                "entries": self_test_effective_entries,
                "path": str(
                    self_test_build_root /
                    BUILD_CONFIGURATION_RELATIVE_PATH
                ),
                "schema": BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
                "sha256": self_test_digest,
            },
        },
        "cell": exact[0],
        "seed": 12345,
        "source_identity": {
            "git": {
                "head": "1" * 40, "tree": "2" * 40,
                "worktree_clean": True,
            },
        },
    }
    raw_settings = {
        "mode": "screen",
        "benchmark": {
            "batch": 1, "iterations": 15, "reuse": 64, "warmups": 4,
        },
    }
    raw_fixture = {
        "schema": BENCHMARK_SCHEMA,
        "build": {
            "backend_variant": "avx2",
            "build_type": "Release",
            "build_configuration_sha256": "",
            "compiler": "self-test", "compiler_version": "1",
            "cplusplus": 201103, "test_hooks": True,
            "source_commit": "1" * 40,
            "source_tree": "2" * 40,
            "source_tracked_dirty": False,
        },
        "parameters": {
            "K": 2, "R": 16, "Q": 1,
            "forced_mode": "force_direct",
            "shard_bytes": 4096,
            "seed": 12345,
            "requested_parity_indices": [0],
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8",
            "batch": 1, "reuse": 64, "iterations": 15, "warmups": 4,
            "thread_count": 1,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": "avx2",
            "direct_capable": True, "thread_count": 1,
            "parent_count": 32, "padded_side": 16,
            "timed_path_is_direct": True,
        },
        "correctness": {
            "selected_transform_reference_parity_match": True,
            "direct_transform_parity_match": True,
            "unrequested_outputs_untouched": True,
            "parity_checksum_fnv1a64": "0x0123456789abcdef",
        },
        "memory": {},
        "operation_model": {
            "direct_row_terms": 2,
            "direct_output_initializations": 1,
            "direct_output_accumulations": 1,
            "fixed_coefficient_symbol_terms_before_unit_specialization":
                8192,
            "xor_accumulation_symbols": 4096,
            "modeled_source_bytes_read": 8192,
            "modeled_output_bytes_read": 4096,
            "modeled_output_bytes_written": 8192,
            "model_scope": (
                "direct streaming kernels before cache effects; unit "
                "coefficients specialize to copy/XOR"
            ),
            "direct_model_applies_to_timed_path": True,
            "transform_operation_counts": None,
            "hardware_counters": None,
        },
        "metrics": {
            # This deliberate distractor is the value selected by the invalid
            # recursive report-20 extractor.
            "codec_setup": {"median_us": 987654.0, "mad_us": 0.0},
            "encode_execution": {
                "median_us_per_batch_call": 7.25,
                "mad_us_per_batch_call": 0.5,
                "minimum_us_per_batch_call": 7.0,
                "maximum_us_per_batch_call": 8.0,
                "logical_input_GB_per_s": 1.0,
                "requested_parity_output_GB_per_s": 1.0,
                "logical_input_plus_output_GB_per_s": 2.0,
            },
        },
        "methodology": {},
    }
    raw_fixture["build"]["build_configuration_sha256"] = self_test_digest
    median_us, mad_us, parity_identity = validate_raw(
        raw_fixture, raw_job, "direct", raw_settings
    )
    check(median_us == 7.25 and mad_us == 0.5,
          "exact encode-execution metric path")
    check(parity_identity == {
        "algorithm": "fnv1a64",
        "digest": "0x0123456789abcdef",
        "hashed_bytes": 4096,
        "requested_parity_indices": [0],
    }, "non-vacuous parity identity")
    for label, encoded in (
            ("duplicate key", b'{"metric": 1, "metric": 2}'),
            ("non-standard NaN", b'{"metric": NaN}'),
            ("non-standard infinity", b'{"metric": Infinity}'),
            ("overflowing exponent", b'{"metric": 1e9999}'),
            ("huge integer", b'{"metric": 100000000000000000000}')):
        try:
            decode_json_bytes(encoded, "self-test {}".format(label))
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} JSON was accepted".format(label)
            )
    try:
        decode_json_bytes(
            b"[" * 10000 + b"0" + b"]" * 10000,
            "self-test deeply nested JSON"
        )
    except CrossoverError:
        pass
    else:
        raise CrossoverError(
            "self-test failed: deeply nested JSON was accepted"
        )

    def reject_raw_mutation(name, mutator):
        changed = json.loads(json.dumps(raw_fixture))
        mutator(changed)
        try:
            validate_raw(changed, raw_job, "direct", raw_settings)
        except CrossoverError:
            return
        raise CrossoverError(
            "self-test failed: raw mutation was accepted: {}".format(name)
        )

    def reject_job_mutation(name, mutator):
        changed = json.loads(json.dumps(raw_job))
        mutator(changed)
        try:
            validate_raw(
                raw_fixture, changed, "direct", raw_settings
            )
        except CrossoverError:
            return
        raise CrossoverError(
            "self-test failed: job mutation was accepted: {}".format(name)
        )

    reject_raw_mutation(
        "wrong schema", lambda value: value.update({"schema": "unknown"})
    )
    reject_raw_mutation(
        "mismatched effective-configuration digest",
        lambda value: value["build"].update(
            {"build_configuration_sha256": "b" * 64}
        )
    )
    reject_job_mutation(
        "missing effective-configuration attestation",
        lambda value: value["build_metadata"].pop(
            "effective_configuration_attestation"
        )
    )
    reject_job_mutation(
        "malformed effective-configuration digest",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ].update({"sha256": "not-a-digest"})
    )
    reject_job_mutation(
        "mutated cache schema binding",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA": "unknown"
        })
    )
    reject_job_mutation(
        "mutated cache digest binding",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "b" * 64
        })
    )
    reject_job_mutation(
        "small-direct mode cache mismatch",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "1"
        })
    )
    reject_job_mutation(
        "direct-source-plan cache mismatch",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "ON"
        })
    )
    reject_job_mutation(
        "high-direct-encode cache mismatch",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "ON"
        })
    )
    reject_job_mutation(
        "general-one-loss-direct cache mismatch",
        lambda value: value["build_metadata"]["entries"].update({
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON"
        })
    )
    reject_job_mutation(
        "unexpected effective-configuration path",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ].update({"path": "/self-test/other.txt"})
    )
    reject_job_mutation(
        "NUL effective-configuration path",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ].update({"path": "/self-test/\0configuration.txt"})
    )
    reject_job_mutation(
        "NUL CMake build root",
        lambda value: value["build_metadata"].update({
            "build_root": "/self-test/\0build"
        })
    )
    reject_job_mutation(
        "NUL CMake source root",
        lambda value: value["build_metadata"]["entries"].update({
            "CMAKE_HOME_DIRECTORY": "/self-test/\0source"
        })
    )
    reject_job_mutation(
        "multiline effective-configuration entry",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ]["entries"].update({"CMAKE_CXX_FLAGS": "-Wall\ninjected"})
    )
    reject_job_mutation(
        "non-UTF-8 effective-configuration entry",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ]["entries"].update({"CMAKE_CXX_FLAGS": "\ud800"})
    )
    reject_job_mutation(
        "extra effective-configuration field",
        lambda value: value["build_metadata"][
            "effective_configuration_attestation"
        ].update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "missing exact metric object",
        lambda value: value["metrics"].pop("encode_execution")
    )
    reject_raw_mutation(
        "codec-setup substitution",
        lambda value: value["metrics"]["encode_execution"].pop(
            "median_us_per_batch_call"
        )
    )
    reject_raw_mutation(
        "numeric string",
        lambda value: value["metrics"]["encode_execution"].update(
            {"median_us_per_batch_call": "7.25"}
        )
    )
    reject_raw_mutation(
        "boolean median",
        lambda value: value["metrics"]["encode_execution"].update(
            {"median_us_per_batch_call": True}
        )
    )
    for label, invalid in (
            ("NaN median", float("nan")),
            ("positive-infinity median", float("inf")),
            ("negative-infinity median", float("-inf")),
            ("overflowing integer median", 10 ** 10000)):
        reject_raw_mutation(
            label,
            lambda value, invalid=invalid:
                value["metrics"]["encode_execution"].update(
                    {"median_us_per_batch_call": invalid}
                )
        )
    for label, invalid in (
            ("negative MAD", -1.0),
            ("NaN MAD", float("nan")),
            ("infinite MAD", float("inf"))):
        reject_raw_mutation(
            label,
            lambda value, invalid=invalid:
                value["metrics"]["encode_execution"].update(
                    {"mad_us_per_batch_call": invalid}
                )
        )
    reject_raw_mutation(
        "empty parity digest string",
        lambda value: value["correctness"].update(
            {"parity_checksum_fnv1a64": ""}
        )
    )
    reject_raw_mutation(
        "empty parity digest object",
        lambda value: value["correctness"].update(
            {"parity_checksum_fnv1a64": {}}
        )
    )
    reject_raw_mutation(
        "empty correctness object",
        lambda value: value.update({"correctness": {}})
    )
    reject_raw_mutation(
        "extra top-level field",
        lambda value: value.update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra build field",
        lambda value: value["build"].update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra parameter field",
        lambda value: value["parameters"].update({"unversioned_extra": 1})
    )
    reject_raw_mutation(
        "extra resolved field",
        lambda value: value["resolved"].update({"unversioned_extra": 1})
    )
    reject_raw_mutation(
        "extra correctness field",
        lambda value: value["correctness"].update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra operation-model field",
        lambda value: value["operation_model"].update(
            {"unversioned_extra": 1}
        )
    )
    reject_raw_mutation(
        "extra encode-execution field",
        lambda value: value["metrics"]["encode_execution"].update(
            {"unversioned_extra": 1}
        )
    )
    reject_raw_mutation(
        "missing operation model",
        lambda value: value.pop("operation_model")
    )
    reject_raw_mutation(
        "flipped operation-model applicability",
        lambda value: value["operation_model"].update(
            {"direct_model_applies_to_timed_path": False}
        )
    )
    reject_raw_mutation(
        "wrong direct-row term count",
        lambda value: value["operation_model"].update(
            {"direct_row_terms": 3}
        )
    )
    reject_raw_mutation(
        "boolean K parameter",
        lambda value: value["parameters"].update({"K": True})
    )
    reject_raw_mutation(
        "wrong embedded backend variant",
        lambda value: value["build"].update({"backend_variant": "auto"})
    )
    reject_raw_mutation(
        "wrong embedded source commit",
        lambda value: value["build"].update({"source_commit": "3" * 40})
    )

    transform_fixture = json.loads(json.dumps(raw_fixture))
    transform_fixture["parameters"]["forced_mode"] = "force_transform"
    transform_fixture["resolved"]["timed_path_is_direct"] = False
    transform_fixture["correctness"]["direct_transform_parity_match"] = None
    transform_fixture["operation_model"][
        "direct_model_applies_to_timed_path"
    ] = False
    validate_raw(transform_fixture, raw_job, "transform", raw_settings)
    transform_fixture["correctness"]["direct_transform_parity_match"] = True
    try:
        validate_raw(transform_fixture, raw_job, "transform", raw_settings)
    except CrossoverError:
        pass
    else:
        raise CrossoverError(
            "self-test failed: transform invocation claimed direct parity"
        )

    measurements = [
        {"timed_mode": "direct", "median_us": 8.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "direct", "median_us": 8.0},
    ]
    summary = summarize_measurements(measurements)
    check(summary["gain_percent"] == 25.0, "gain summary")
    check(math.isclose(
        summary["rounds"][0]["gain_percent"], 25.0,
        rel_tol=1e-12, abs_tol=1e-12
    ), "ABBA round summary")
    inference_summary = summarize_measurements(measurements * 3)
    inference = inference_summary["paired_log_inference"]
    check(
        inference["degrees_of_freedom"] == 2 and
        abs(inference["gain_percent"] - 25.0) < 1e-12 and
        abs(inference["gain_percent_student_t_interval"][0] - 25.0) < 1e-12,
        "three-round paired-log Student-t inference"
    )
    asymmetric = summarize_measurements([
        {"timed_mode": "direct", "median_us": 4.0},
        {"timed_mode": "transform", "median_us": 9.0},
        {"timed_mode": "transform", "median_us": 25.0},
        {"timed_mode": "direct", "median_us": 16.0},
    ])
    check(math.isclose(
        asymmetric["rounds"][0]["direct_geometric_mean_us"], 8.0,
        rel_tol=1e-15, abs_tol=1e-15
    ) and math.isclose(
        asymmetric["rounds"][0]["transform_geometric_mean_us"], 15.0,
        rel_tol=1e-15, abs_tol=1e-15
    ), "ABBA round uses geometric representatives")
    for label, invalid in (
            ("string", "8.0"), ("boolean", True),
            ("NaN", float("nan")), ("infinity", float("inf"))):
        changed = json.loads(json.dumps(measurements))
        changed[0]["median_us"] = invalid
        try:
            summarize_measurements(changed)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} summary metric was accepted".format(label)
            )
    fake = {
        "cell": candidates[0], "job_id": "a", "status": "passed",
        "summary": summary,
    }
    analysis = analyze_results([fake], 5.0)
    check(analysis["candidate"]["gain_min_percent"] == 25.0,
          "candidate minimum gain")
    check(analysis["candidate"]["regression_count"] == 0,
          "candidate regression count")
    with tempfile.TemporaryDirectory(prefix="leo2-crossover-self-test-") as directory:
        root = Path(directory)
        path = root / "stable.json"
        atomic_write_json(path, {"z": [3, 2, 1], "a": "stable"})
        check(json.loads(path.read_text(encoding="utf-8"))["a"] == "stable",
              "atomic JSON output")
        atomic_victim = root / "atomic-victim"
        atomic_victim.mkdir()
        redirected_parent = root / "redirected-parent"
        redirected_parent.symlink_to(atomic_victim.name)
        try:
            atomic_write_json(
                redirected_parent / "manifest.json", {"unsafe": True}
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: atomic JSON followed a parent symlink"
            )
        check(
            not (atomic_victim / "manifest.json").exists(),
            "atomic JSON parent symlink victim is untouched"
        )
        race_parent = root / "raced-parent"
        race_parent.mkdir(mode=0o700)
        race_leaf = race_parent / "leaf"
        race_leaf.mkdir(mode=0o700)
        race_target = root / "raced-target"
        race_target.mkdir(mode=0o700)
        (race_target / "leaf").mkdir(mode=0o700)
        race_saved = root / "raced-parent-saved"
        original_os_open = os.open
        raced_component = [False]

        def race_intermediate_component(path_value, flags, *values, **options):
            if (not raced_component[0] and path_value == race_parent.name and
                    options.get("dir_fd") is not None):
                raced_component[0] = True
                race_parent.rename(race_saved)
                race_parent.symlink_to(
                    race_target.name, target_is_directory=True
                )
            return original_os_open(path_value, flags, *values, **options)

        try:
            os.open = race_intermediate_component
            try:
                owned_canonical_directory(
                    race_leaf, "self-test raced directory"
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: raced intermediate symlink was followed"
                )
        finally:
            os.open = original_os_open
        check(
            raced_component[0],
            "componentwise directory race fixture executed",
        )

        class InjectedAbort(BaseException):
            pass

        def descriptor_closed(descriptor):
            try:
                os.fstat(descriptor)
            except OSError as error:
                return error.errno == errno.EBADF
            return False

        def source_line_number(function, fragment):
            lines, first_line = inspect.getsourcelines(function)
            matches = [
                first_line + index for index, line in enumerate(lines)
                if line.strip() == fragment
            ]
            if len(matches) != 1:
                raise CrossoverError(
                    "self-test cannot locate unique {!r} in {}".format(
                        fragment, function.__name__
                    )
                )
            return matches[0]

        def run_trace_abort(function, predicate, invocation, description):
            triggered = [False]

            def abort_trace(frame, event, unused_argument):
                if (not triggered[0] and event == "line" and
                        frame.f_code is function.__code__ and
                        predicate(frame)):
                    triggered[0] = True
                    raise InjectedAbort(description)
                return abort_trace

            previous_trace = sys.gettrace()
            try:
                sys.settrace(abort_trace)
                try:
                    invocation()
                except InjectedAbort:
                    pass
                else:
                    raise CrossoverError(
                        "self-test failed: {} was ignored".format(description)
                    )
            finally:
                sys.settrace(previous_trace)
            check(
                triggered[0],
                "{} trace fixture executed".format(description),
            )

        sealed_transfer_descriptors = []
        sealed_return_line = source_line_number(
            create_sealed_executable_snapshot, "return result"
        )

        def abort_sealed_return(frame):
            if frame.f_lineno != sealed_return_line:
                return False
            sealed_transfer_descriptors.append(
                frame.f_locals["result"]["descriptor"]
            )
            return True

        run_trace_abort(
            create_sealed_executable_snapshot, abort_sealed_return,
            lambda: create_sealed_executable_snapshot(
                b"#!/bin/sh\nexit 0\n",
                "self-test interrupted sealed snapshot",
            ),
            "sealed-snapshot return BaseException",
        )
        check(
            len(sealed_transfer_descriptors) == 1 and
            descriptor_closed(sealed_transfer_descriptors[0]),
            "sealed-snapshot return interruption closes its retained FD",
        )

        caller_context_snapshot = None
        try:
            raise RuntimeError("self-test active caller exception")
        except RuntimeError:
            caller_context_snapshot = create_sealed_executable_snapshot(
                b"#!/bin/sh\nexit 0\n",
                "self-test caller-exception snapshot",
            )
        try:
            check(
                not descriptor_closed(
                    caller_context_snapshot["descriptor"]
                ) and
                validate_sealed_executable_snapshot(
                    caller_context_snapshot,
                    caller_context_snapshot["sha256"],
                    "self-test caller-exception snapshot",
                ) == b"#!/bin/sh\nexit 0\n",
                "a caller's actively handled exception does not revoke a "
                "successful descriptor transfer",
            )
        finally:
            close_sealed_snapshot(caller_context_snapshot)

        exact_executable_path = root / "return-interrupt-executable"
        exact_executable_path.write_bytes(b"#!/bin/sh\nexit 0\n")
        exact_executable_path.chmod(0o755)
        exact_snapshot_descriptors = []
        exact_return_line = source_line_number(
            open_exact_executable_snapshot, "return snapshot"
        )

        def abort_exact_snapshot_return(frame):
            if frame.f_lineno != exact_return_line:
                return False
            exact_snapshot_descriptors.append(
                frame.f_locals["snapshot"]["descriptor"]
            )
            return True

        run_trace_abort(
            open_exact_executable_snapshot, abort_exact_snapshot_return,
            lambda: open_exact_executable_snapshot(
                exact_executable_path,
                "self-test interrupted exact executable",
            ),
            "exact-executable return BaseException",
        )
        check(
            len(exact_snapshot_descriptors) == 1 and
            descriptor_closed(exact_snapshot_descriptors[0]),
            "exact-executable return interruption closes its sealed FD",
        )

        fd_fixture = owned_canonical_directory(
            root, "self-test descriptor-cleanup directory"
        )
        original_open_components = globals()[
            "_open_absolute_directory_components"
        ]
        original_fstat = os.fstat
        opened_current = []

        def open_current_then_abort(
                unused_path, unused_description, unused_create,
                required_mode=0o700):
            descriptor = os.dup(fd_fixture["descriptor"])
            opened_current.append(descriptor)
            metadata = original_fstat(descriptor)
            return {
                "descriptor": descriptor,
                "identity": (metadata.st_dev, metadata.st_ino),
                "path": root,
                "required_mode": required_mode,
            }

        def abort_current_fstat(descriptor):
            if opened_current and descriptor == opened_current[-1]:
                raise InjectedAbort("injected current-directory fstat abort")
            return original_fstat(descriptor)

        try:
            globals()["_open_absolute_directory_components"] = (
                open_current_then_abort
            )
            os.fstat = abort_current_fstat
            try:
                validate_owned_directory(
                    fd_fixture, "self-test descriptor cleanup"
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: directory validation abort was ignored"
                )
        finally:
            os.fstat = original_fstat
            globals()["_open_absolute_directory_components"] = (
                original_open_components
            )
        check(
            len(opened_current) == 1 and
            descriptor_closed(opened_current[0]),
            "directory revalidation closes a newly opened FD on "
            "BaseException",
        )

        held_log_path = root / "fd-validation.log"
        held_log = open_exclusive_log(
            held_log_path, "self-test held validation log"
        )
        original_stat = os.stat
        log_current = []

        def open_log_current(
                unused_path, unused_description, unused_create,
                required_mode=0o700):
            descriptor = os.dup(held_log["directory_descriptor"])
            log_current.append(descriptor)
            metadata = original_fstat(descriptor)
            return {
                "descriptor": descriptor,
                "identity": (metadata.st_dev, metadata.st_ino),
                "path": root,
                "required_mode": required_mode,
            }

        def abort_log_stat(*unused_values, **unused_options):
            raise InjectedAbort("injected log stat abort")

        try:
            globals()["_open_absolute_directory_components"] = (
                open_log_current
            )
            os.stat = abort_log_stat
            try:
                validate_log_identity(
                    held_log, "self-test log descriptor cleanup"
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: log validation abort was ignored"
                )
        finally:
            os.stat = original_stat
            globals()["_open_absolute_directory_components"] = (
                original_open_components
            )
        check(
            len(log_current) == 1 and descriptor_closed(log_current[0]),
            "log revalidation closes a newly opened directory FD on "
            "BaseException",
        )
        close_log(held_log)

        original_open_owned = globals()["open_owned_directory"]
        transferred = []

        def open_owned_then_abort(unused_path, unused_description):
            descriptor = os.dup(fd_fixture["descriptor"])
            transferred.append(descriptor)
            return descriptor

        def abort_transferred_fstat(descriptor):
            if transferred and descriptor == transferred[-1]:
                raise InjectedAbort("injected canonical-directory abort")
            return original_fstat(descriptor)

        try:
            globals()["open_owned_directory"] = open_owned_then_abort
            os.fstat = abort_transferred_fstat
            try:
                open_existing_canonical_directory(
                    root, "self-test canonical transfer cleanup"
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: canonical transfer abort was ignored"
                )
        finally:
            os.fstat = original_fstat
            globals()["open_owned_directory"] = original_open_owned
        check(
            len(transferred) == 1 and descriptor_closed(transferred[0]),
            "canonical directory opener closes its FD before ownership "
            "transfer on BaseException",
        )

        canonical_return_descriptors = []
        canonical_return_line = source_line_number(
            open_existing_canonical_directory, "return result"
        )

        def abort_canonical_return(frame):
            if frame.f_lineno != canonical_return_line:
                return False
            canonical_return_descriptors.append(
                frame.f_locals["result"]["descriptor"]
            )
            return True

        run_trace_abort(
            open_existing_canonical_directory, abort_canonical_return,
            lambda: open_existing_canonical_directory(
                root, "self-test canonical return cleanup"
            ),
            "canonical-directory return BaseException",
        )
        check(
            len(canonical_return_descriptors) == 1 and
            descriptor_closed(canonical_return_descriptors[0]),
            "canonical directory return interruption closes its FD",
        )
        close_owned_directory(fd_fixture)

        abort_log_path = root / "abort-open.log"
        original_fchmod = os.fchmod
        abort_fchmod_once = [True]

        def abort_created_log_fchmod(descriptor, mode):
            if abort_fchmod_once[0]:
                abort_fchmod_once[0] = False
                raise InjectedAbort("injected exclusive-log abort")
            return original_fchmod(descriptor, mode)

        try:
            os.fchmod = abort_created_log_fchmod
            try:
                open_exclusive_log(
                    abort_log_path, "self-test aborting exclusive log"
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: exclusive-log abort was ignored"
                )
        finally:
            os.fchmod = original_fchmod
        check(
            not abort_log_path.exists(),
            "exclusive-log BaseException removes the exact newly created "
            "inode",
        )

        return_abort_log_path = root / "abort-log-return.log"
        returned_log_descriptors = []
        log_return_line = source_line_number(
            open_exclusive_log, "return result"
        )

        def abort_log_return(frame):
            if frame.f_lineno != log_return_line:
                return False
            returned_log_descriptors.extend((
                frame.f_locals["result"]["descriptor"],
                frame.f_locals["result"]["directory_descriptor"],
            ))
            return True

        run_trace_abort(
            open_exclusive_log, abort_log_return,
            lambda: open_exclusive_log(
                return_abort_log_path,
                "self-test interrupted log return",
            ),
            "exclusive-log return BaseException",
        )
        check(
            len(returned_log_descriptors) == 2 and
            all(descriptor_closed(value)
                for value in returned_log_descriptors) and
            not return_abort_log_path.exists(),
            "exclusive-log return interruption closes both FDs and removes "
            "its inode",
        )

        abort_raw_directory = owned_canonical_directory(
            root, "self-test aborting raw-file directory"
        )
        abort_fchmod_once[0] = True
        try:
            os.fchmod = abort_created_log_fchmod
            try:
                open_exclusive_owned_file(
                    abort_raw_directory, "abort-raw.json",
                    "self-test aborting raw file",
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: raw-file abort was ignored"
                )
        finally:
            os.fchmod = original_fchmod
            close_owned_directory(abort_raw_directory)
        check(
            not (root / "abort-raw.json").exists(),
            "raw-file BaseException removes the exact newly created inode",
        )

        raw_directory = owned_canonical_directory(
            root, "self-test raw-directory transfer"
        )
        original_dup = os.dup
        raw_directory_duplicates = []

        def capture_raw_directory_dup(descriptor):
            duplicate = original_dup(descriptor)
            if descriptor == raw_directory["descriptor"]:
                raw_directory_duplicates.append(duplicate)
            return duplicate

        def abort_after_raw_directory_dup(frame):
            return (
                frame.f_locals.get("directory_descriptor") is not None and
                frame.f_locals.get("descriptor") is None
            )

        try:
            os.dup = capture_raw_directory_dup
            run_trace_abort(
                open_exclusive_owned_file,
                abort_after_raw_directory_dup,
                lambda: open_exclusive_owned_file(
                    raw_directory, "pre-file-abort.json",
                    "self-test interrupted raw-directory transfer",
                ),
                "raw-directory duplicate BaseException",
            )
        finally:
            os.dup = original_dup
            close_owned_directory(raw_directory)
        check(
            len(raw_directory_duplicates) == 1 and
            descriptor_closed(raw_directory_duplicates[0]) and
            not (root / "pre-file-abort.json").exists(),
            "raw-file setup interruption closes its directory duplicate",
        )

        command_stdout = root / "abort-command.stdout.log"
        command_stderr = root / "abort-command.stderr.log"
        original_open_exclusive = globals()["open_exclusive_log"]
        log_open_count = [0]

        def abort_second_log(path_value, description):
            log_open_count[0] += 1
            if log_open_count[0] == 2:
                raise InjectedAbort("injected stderr reservation abort")
            return original_open_exclusive(path_value, description)

        try:
            globals()["open_exclusive_log"] = abort_second_log
            try:
                open_command_logs(command_stdout, command_stderr)
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: command-log abort was ignored"
                )
        finally:
            globals()["open_exclusive_log"] = original_open_exclusive
        check(
            not command_stdout.exists() and not command_stderr.exists(),
            "stderr reservation BaseException rolls back the exact stdout "
            "inode",
        )

        transfer_stdout = root / "transfer-command.stdout.log"
        transfer_stderr = root / "transfer-command.stderr.log"

        def abort_after_stdout_reservation(frame):
            return (
                frame.f_locals.get("stdout_log") is not None and
                frame.f_locals.get("stderr_log") is None
            )

        run_trace_abort(
            open_command_logs, abort_after_stdout_reservation,
            lambda: open_command_logs(transfer_stdout, transfer_stderr),
            "stdout-log ownership-transfer BaseException",
        )
        check(
            not transfer_stdout.exists() and not transfer_stderr.exists(),
            "stdout ownership-transfer interruption rolls back both log "
            "names",
        )

        publication_path = root / "held-publication"
        publication_path.mkdir(mode=0o700)
        publication_root = owned_canonical_directory(
            publication_path, "self-test held publication root"
        )
        original_link_descriptor = globals()["link_descriptor_noreplace"]
        publication_saved = root / "held-publication-saved"
        publication_replacement = root / "held-publication-replacement"

        def publication_aba(
                source_descriptor, destination_descriptor,
                destination_name, description):
            publication_path.rename(publication_saved)
            publication_path.mkdir(mode=0o700)
            original_link_descriptor(
                source_descriptor, destination_descriptor,
                destination_name, description,
            )
            publication_path.rename(publication_replacement)
            publication_saved.rename(publication_path)

        try:
            globals()["link_descriptor_noreplace"] = publication_aba
            write_result_json(
                publication_root, "manifest.json", {"held": True},
                "self-test held publication", replace=False,
            )
        finally:
            globals()["link_descriptor_noreplace"] = original_link_descriptor
        check(
            json.loads((publication_path / "manifest.json").read_text(
                encoding="utf-8"
            )) == {"held": True} and
            not (publication_replacement / "manifest.json").exists(),
            "publication remains bound to the held result-root inode",
        )

        held_return_descriptors = []
        held_return_line = source_line_number(
            open_result_regular_held, "return held"
        )

        def abort_held_return(frame):
            if frame.f_lineno != held_return_line:
                return False
            held_return_descriptors.extend((
                frame.f_locals["held"]["descriptor"],
                frame.f_locals["held"]["directory_descriptor"],
            ))
            return True

        run_trace_abort(
            open_result_regular_held, abort_held_return,
            lambda: open_result_regular_held(
                publication_root, "manifest.json", MAX_RAW_JSON_BYTES,
                "self-test interrupted held-file return", 0o600,
            ),
            "held-file return BaseException",
        )
        check(
            len(held_return_descriptors) == 2 and
            all(descriptor_closed(value)
                for value in held_return_descriptors),
            "held-file return interruption closes file and directory FDs",
        )

        collision_path = publication_path / "collision.json"

        def collide_exact_publication():
            collision_descriptor = os.open(
                collision_path.name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600, dir_fd=publication_root["descriptor"],
            )
            try:
                write_descriptor_all(collision_descriptor, b"collision\n")
                os.fsync(collision_descriptor)
            finally:
                os.close(collision_descriptor)

        try:
            try:
                write_result_bytes(
                    publication_root, collision_path.name,
                    canonical_json_bytes({"must_not_publish": True}),
                    "self-test collided publication", replace=False,
                    precommit_hook=collide_exact_publication,
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: no-replace publication collision "
                    "was accepted"
                )
        finally:
            pass
        check(
            collision_path.read_bytes() == b"collision\n",
            "no-replace publication leaves a raced destination untouched",
        )

        replacement_path = publication_path / "replace.json"
        write_result_bytes(
            publication_root, replacement_path.name, b"old!",
            "self-test replace baseline", replace=False,
        )

        interrupted_new_path = publication_path / "interrupted-new.json"
        original_rename_noreplace = globals()["rename_noreplace"]
        original_unlink = os.unlink
        interrupted_new_unlink = [False]

        def commit_new_then_abort(*values, **options):
            original_rename_noreplace(*values, **options)
            raise InjectedAbort("injected post-rename new-publication abort")

        def unlink_new_then_abort(path_value, *values, **options):
            result = original_unlink(path_value, *values, **options)
            if (path_value == interrupted_new_path.name and
                    not interrupted_new_unlink[0]):
                interrupted_new_unlink[0] = True
                raise InjectedAbort(
                    "injected post-unlink new-rollback abort"
                )
            return result

        try:
            globals()["rename_noreplace"] = commit_new_then_abort
            os.unlink = unlink_new_then_abort
            try:
                write_result_bytes(
                    publication_root, interrupted_new_path.name, b"new!",
                    "self-test interrupted new publication", replace=False,
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: post-rename new-publication abort "
                    "was ignored"
                )
        finally:
            os.unlink = original_unlink
            globals()["rename_noreplace"] = original_rename_noreplace
        check(
            interrupted_new_unlink[0] and
            not interrupted_new_path.exists() and
            not any(
                path.name.startswith(".leo2-stage-")
                for path in publication_path.iterdir()
            ),
            "post-rename/post-unlink BaseExceptions restore an absent "
            "new-publication destination",
        )

        original_rename_exchange = globals()["rename_exchange"]
        interrupted_exchange_count = [0]

        def exchange_then_abort(*values, **options):
            original_rename_exchange(*values, **options)
            interrupted_exchange_count[0] += 1
            raise InjectedAbort(
                "injected post-exchange replacement abort"
            )

        try:
            globals()["rename_exchange"] = exchange_then_abort
            try:
                write_result_bytes(
                    publication_root, replacement_path.name, b"new!",
                    "self-test interrupted replacement", replace=True,
                )
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: post-exchange replacement abort "
                    "was ignored"
                )
        finally:
            globals()["rename_exchange"] = original_rename_exchange
        check(
            interrupted_exchange_count[0] == 2 and
            replacement_path.read_bytes() == b"old!" and
            not any(
                path.name.startswith(".leo2-stage-")
                for path in publication_path.iterdir()
            ),
            "post-exchange BaseExceptions reconcile both commit and rollback "
            "without changing prior evidence",
        )

        def fail_before_publication():
            raise CrossoverError("injected precommit failure")

        try:
            write_result_bytes(
                publication_root, replacement_path.name, b"new!",
                "self-test precommit replacement", replace=True,
                precommit_hook=fail_before_publication,
            )
        except CrossoverError as error:
            check(
                "injected precommit failure" in str(error),
                "precommit publication failure is surfaced",
            )
        else:
            raise CrossoverError(
                "self-test failed: injected precommit failure was ignored"
            )
        check(
            replacement_path.read_bytes() == b"old!",
            "precommit replacement failure preserves old evidence",
        )

        original_fsync = os.fsync
        post_replace_fsync_target = [None]

        def arm_post_replace_fsync(
                unused_descriptor, directory_descriptor, unused_name):
            post_replace_fsync_target[0] = directory_descriptor

        def fail_one_post_replace_fsync(descriptor):
            if descriptor == post_replace_fsync_target[0]:
                post_replace_fsync_target[0] = None
                raise OSError(
                    errno.EIO, "injected post-replace fsync failure"
                )
            return original_fsync(descriptor)

        try:
            os.fsync = fail_one_post_replace_fsync
            try:
                write_result_bytes(
                    publication_root, replacement_path.name, b"new!",
                    "self-test post-replace fsync rollback", replace=True,
                    postcommit_hook=arm_post_replace_fsync,
                )
            except CrossoverError as error:
                check(
                    "injected post-replace fsync failure" in str(error),
                    "post-replace fsync failure is surfaced after rollback",
                )
            else:
                raise CrossoverError(
                    "self-test failed: post-replace fsync failure was ignored"
                )
        finally:
            os.fsync = original_fsync
        check(
            replacement_path.read_bytes() == b"old!" and
            not any(
                path.name.startswith(".leo2-stage-")
                for path in publication_path.iterdir()
            ),
            "post-replace fsync failure restores old evidence and removes "
            "staging",
        )

        def mutate_linked_staging(
                descriptor, unused_directory, unused_name):
            os.pwrite(descriptor, b"evil", 0)
            os.fsync(descriptor)

        try:
            write_result_bytes(
                publication_root, replacement_path.name, b"new!",
                "self-test mutated staged replacement", replace=True,
                postlink_hook=mutate_linked_staging,
            )
        except CrossoverError as error:
            check(
                "staged inode or payload changed" in str(error),
                "same-size post-link staging mutation is rejected",
            )
        else:
            raise CrossoverError(
                "self-test failed: same-size staged mutation was accepted"
            )
        check(
            replacement_path.read_bytes() == b"old!",
            "mutated staged replacement preserves old evidence",
        )

        def mutate_after_replace(
                descriptor, unused_directory, unused_name):
            os.pwrite(descriptor, b"evil", 0)
            os.fsync(descriptor)

        try:
            write_result_bytes(
                publication_root, replacement_path.name, b"new!",
                "self-test post-replace verification rollback",
                replace=True, postcommit_hook=mutate_after_replace,
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: post-replace payload mutation was "
                "accepted"
            )
        check(
            replacement_path.read_bytes() == b"old!",
            "post-replace verification failure restores old evidence",
        )
        write_result_bytes(
            publication_root, replacement_path.name, b"new!",
            "self-test successful atomic replacement", replace=True,
        )
        check(
            replacement_path.read_bytes() == b"new!" and
            stat.S_IMODE(replacement_path.stat().st_mode) == 0o600,
            "replacement atomically publishes exact mode-0600 evidence",
        )
        final_fsync_count = [0]
        final_fsync_descriptor = [None]

        def arm_cleanup_fsync(
                unused_descriptor, directory_descriptor, unused_name):
            final_fsync_descriptor[0] = directory_descriptor

        def fail_cleanup_fsync_after_commit(descriptor):
            if descriptor == final_fsync_descriptor[0]:
                final_fsync_count[0] += 1
                if final_fsync_count[0] == 2:
                    raise OSError(
                        errno.EIO, "injected irreversible cleanup fsync"
                    )
            return original_fsync(descriptor)

        try:
            os.fsync = fail_cleanup_fsync_after_commit
            write_result_bytes(
                publication_root, replacement_path.name, b"done",
                "self-test irreversible replacement", replace=True,
                postcommit_hook=arm_cleanup_fsync,
            )
        finally:
            os.fsync = original_fsync
        check(
            final_fsync_count[0] == 2 and
            replacement_path.read_bytes() == b"done",
            "post-cleanup fsync failure never reports a false replacement "
            "failure after old evidence is irreversibly removed",
        )
        replacement_path.chmod(0o640)
        try:
            read_result_regular(
                publication_root, replacement_path.name, 16,
                "self-test weak-mode evidence",
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: non-0600 evidence mode was accepted"
            )
        replacement_path.chmod(0o600)
        close_owned_directory(publication_root)

        weak_directory = root / "weak-mode-directory"
        weak_directory.mkdir(mode=0o700)
        weak_directory.chmod(0o750)
        try:
            owned_canonical_directory(
                weak_directory, "self-test weak-mode directory"
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: non-0700 result directory was accepted"
            )

        retained_path = root / "retained-result"
        retained_path.mkdir(mode=0o700)
        retained_root = owned_canonical_directory(
            retained_path, "self-test retained result"
        )
        try:
            def write_retained_fixture(name, value):
                fixture = retained_path / name
                fixture.write_bytes(value)
                fixture.chmod(0o600)
                return fixture

            write_retained_fixture("stable.bin", b"stable")
            check(
                read_result_regular(
                    retained_root, "stable.bin", 16,
                    "self-test stable retained file",
                ) == b"stable",
                "retained file stable snapshot",
            )
            write_retained_fixture("victim.bin", b"victim")
            (retained_path / "linked.bin").symlink_to("victim.bin")
            os.mkfifo(retained_path / "fifo.bin", 0o600)
            os.link(
                retained_path / "victim.bin",
                retained_path / "hardlink.bin",
            )
            write_retained_fixture("oversize.bin", b"x" * 17)
            for name, maximum in (
                    ("linked.bin", 32), ("fifo.bin", 32),
                    ("hardlink.bin", 32), ("oversize.bin", 16)):
                try:
                    read_result_regular(
                        retained_root, name, maximum,
                        "self-test unsafe retained file",
                    )
                except CrossoverError:
                    pass
                else:
                    raise CrossoverError(
                        "self-test failed: unsafe retained {} was "
                        "accepted".format(name)
                    )

            write_retained_fixture(
                "mutating.bin",
                b"A" * (2 * 1024 * 1024)
            )

            def mutate_retained_file():
                with open(
                        retained_path / "mutating.bin", "r+b",
                        buffering=0) as stream:
                    stream.write(b"B")
                    os.fsync(stream.fileno())

            try:
                read_result_regular(
                    retained_root, "mutating.bin", 3 * 1024 * 1024,
                    "self-test mutating retained file",
                    mutation_hook=mutate_retained_file,
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: in-place retained mutation was accepted"
                )

            inventory_path = retained_path / "inventory"
            inventory_path.mkdir(mode=0o700)
            inventory_expected = inventory_path / "expected.json"
            inventory_expected.write_text("{}", encoding="ascii")
            inventory_expected.chmod(0o600)
            original_listdir = os.listdir
            inventory_raced = [False]

            def race_inventory(descriptor):
                names = original_listdir(descriptor)
                if not inventory_raced[0]:
                    inventory_raced[0] = True
                    extra_descriptor = os.open(
                        "extra.json",
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                        0o600, dir_fd=descriptor,
                    )
                    os.close(extra_descriptor)
                return names

            try:
                os.listdir = race_inventory
                try:
                    list_result_regular_names(
                        retained_root, "inventory", ".json",
                        "self-test raced inventory",
                    )
                except CrossoverError:
                    pass
                else:
                    raise CrossoverError(
                        "self-test failed: raced result inventory was accepted"
                    )
            finally:
                os.listdir = original_listdir
            check(
                inventory_raced[0],
                "result inventory race fixture executed",
            )

            write_retained_fixture(
                "parent-swap.bin",
                b"P" * (2 * 1024 * 1024)
            )
            swapped_path = root / "retained-result-swapped"

            def swap_retained_parent():
                retained_path.rename(swapped_path)
                retained_path.mkdir(mode=0o700)
                replacement = retained_path / "parent-swap.bin"
                replacement.write_bytes(b"replacement")
                replacement.chmod(0o600)

            try:
                read_result_regular(
                    retained_root, "parent-swap.bin", 3 * 1024 * 1024,
                    "self-test parent-swapped retained file",
                    mutation_hook=swap_retained_parent,
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: retained parent swap was accepted"
                )
        finally:
            close_owned_directory(retained_root)
        pidfd_probe_marker = root / "pidfd-probe-must-not-run"
        original_pidfd_open = globals()["linux_pidfd_open"]

        def unavailable_pidfd(unused_pid):
            raise OSError(errno.ENOSYS, "pidfd unavailable")

        try:
            globals()["linux_pidfd_open"] = unavailable_pidfd
            try:
                run_command(
                    [
                        sys.executable, "-c",
                        "from pathlib import Path;"
                        "Path({!r}).write_text('bad')".format(
                            str(pidfd_probe_marker)
                        ),
                    ],
                    root, root / "pidfd-probe.stdout.log",
                    root / "pidfd-probe.stderr.log", 5,
                    os.environ.copy(),
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: unavailable pidfd operations were "
                    "accepted"
                )
        finally:
            globals()["linux_pidfd_open"] = original_pidfd_open
        check(
            not pidfd_probe_marker.exists(),
            "pidfd operations are probed before command launch",
        )

        launch_gate_marker = root / "launch-gate-must-not-run"
        original_retain_launch = globals()["retain_launched_process"]

        def reject_launch_retention(unused_process, unused_description):
            raise CrossoverError("injected launch-retention failure")

        try:
            globals()["retain_launched_process"] = reject_launch_retention
            try:
                run_command(
                    [
                        sys.executable, "-c",
                        "from pathlib import Path;"
                        "Path({!r}).write_text('bad')".format(
                            str(launch_gate_marker)
                        ),
                    ],
                    root, root / "launch-gate.stdout.log",
                    root / "launch-gate.stderr.log", 5,
                    os.environ.copy(),
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: unretained command owner was accepted"
                )
        finally:
            globals()["retain_launched_process"] = original_retain_launch
        check(
            not launch_gate_marker.exists(),
            "command owner cannot run before pidfd identity retention",
        )

        record = run_command(
            [sys.executable, "-c", "print('ok')"], root,
            root / "stdout.log", root / "stderr.log", 5, os.environ.copy()
        )
        check(record["returncode"] == 0, "self-test subprocess exit")
        check((root / "stdout.log").read_bytes() == b"ok\n",
              "self-test subprocess output")

        small_read, small_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        try:
            fcntl.fcntl(small_read, fcntl.F_SETPIPE_SZ, 4096)
            drain = BoundedControlDrain(
                small_read, 16 * 1024, "self-test 4KiB control pipe"
            )
            small_read = None
            payload = b"c" * 8192
            writer = threading.Thread(
                target=write_descriptor_all,
                args=(small_write, payload),
            )
            writer.start()
            writer.join(timeout=5)
            check(not writer.is_alive(), "4KiB control writer did not deadlock")
            os.close(small_write)
            small_write = None
            check(
                drain.finish() == payload,
                "4KiB control pipe was drained concurrently",
            )
        finally:
            for descriptor in (small_read, small_write):
                if descriptor is not None:
                    os.close(descriptor)

        bounded_read, bounded_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0)
        )
        bounded_drain = BoundedControlDrain(
            bounded_read, 1024, "self-test bounded open-writer drain"
        )
        bounded_started = time.monotonic()
        try:
            try:
                bounded_drain.finish(timeout=0.1)
            except CrossoverError as error:
                check(
                    "writer remained open" in str(error),
                    "open-writer drain has an actionable bounded failure",
                )
            else:
                raise CrossoverError(
                    "self-test failed: open-writer drain unexpectedly finished"
                )
            check(
                time.monotonic() - bounded_started < 1.0,
                "control drain finish respects its deadline",
            )
        finally:
            os.close(bounded_write)

        proc_return_descriptors = []
        proc_return_line = source_line_number(
            open_proc_record, "return record"
        )

        def abort_proc_record_return(frame):
            if frame.f_lineno != proc_return_line:
                return False
            proc_return_descriptors.append(
                frame.f_locals["record"]["proc_descriptor"]
            )
            return True

        run_trace_abort(
            open_proc_record, abort_proc_record_return,
            lambda: open_proc_record(os.getpid()),
            "process-record return BaseException",
        )
        check(
            len(proc_return_descriptors) == 1 and
            descriptor_closed(proc_return_descriptors[0]),
            "process-record return interruption closes its retained /proc FD",
        )

        process_snapshot = proc_record(os.getpid())
        check(
            process_snapshot is not None,
            "self-test current process identity is readable",
        )
        retained_return_descriptors = []
        retained_return_line = source_line_number(
            retain_process_identity, "return retained"
        )

        def abort_retained_process_return(frame):
            if frame.f_lineno != retained_return_line:
                return False
            retained_return_descriptors.extend((
                frame.f_locals["retained"]["pidfd"],
                frame.f_locals["retained"]["proc_descriptor"],
            ))
            return True

        run_trace_abort(
            retain_process_identity, abort_retained_process_return,
            lambda: retain_process_identity(
                process_snapshot,
                "self-test interrupted retained process",
            ),
            "retained-process return BaseException",
        )
        check(
            len(retained_return_descriptors) == 2 and
            all(descriptor_closed(value)
                for value in retained_return_descriptors),
            "retained-process return interruption closes pidfd and /proc FD",
        )

        # Deterministically model PID reuse between pidfd_open and the second
        # /proc read.  The transient pidfd must be closed and no identity
        # returned for signaling.
        original_open_proc_record = globals()["open_proc_record"]
        original_pidfd_open = globals()["linux_pidfd_open"]
        reuse_read = os.open("/dev/null", os.O_RDONLY)
        before_proc = os.open("/dev/null", os.O_RDONLY)
        after_proc = os.open("/dev/null", os.O_RDONLY)
        reuse_reads = iter((
            {
                "pid": 4242, "starttime_ticks": 11,
                "proc_identity": (1, 101),
                "proc_descriptor": before_proc,
            },
            {
                "pid": 4242, "starttime_ticks": 11,
                "proc_identity": (1, 102),
                "proc_descriptor": after_proc,
            },
        ))
        try:
            globals()["open_proc_record"] = lambda unused_pid: next(reuse_reads)
            globals()["linux_pidfd_open"] = lambda unused_pid: reuse_read
            check(
                retain_process_identity(
                    {
                        "pid": 4242, "starttime_ticks": 11,
                        "proc_identity": (1, 101),
                    },
                    "self-test reused process",
                ) is None,
                "same-jiffy PID reuse is rejected around pidfd_open",
            )
            try:
                os.fstat(reuse_read)
            except OSError as error:
                check(error.errno == errno.EBADF,
                      "reused PID transient pidfd is closed")
            else:
                raise CrossoverError(
                    "self-test failed: reused PID pidfd leaked"
                )
            reuse_read = None
            before_proc = None
            after_proc = None
        finally:
            globals()["open_proc_record"] = original_open_proc_record
            globals()["linux_pidfd_open"] = original_pidfd_open
            for descriptor in (reuse_read, before_proc, after_proc):
                if descriptor is not None:
                    try:
                        os.close(descriptor)
                    except OSError:
                        pass

        remap_targets = (210, 211, 212, 213, 214)
        saved_remap_targets = {}
        remap_seed_descriptors = []
        remapped_descriptors = ()
        relocated_control = None
        remap_lock_descriptor = None
        try:
            for target in remap_targets:
                try:
                    saved_remap_targets[target] = os.dup(target)
                except OSError as error:
                    if error.errno != errno.EBADF:
                        raise
            seed_identities = {}
            for label, target in zip(("a", "b", "c", "control"),
                                     remap_targets):
                descriptor = linux_memfd_create(
                    "leo2-remap-" + label,
                    getattr(os, "MFD_CLOEXEC", 0x0001),
                )
                remap_seed_descriptors.append(descriptor)
                write_descriptor_all(
                    descriptor, label.encode("ascii")
                )
                os.dup2(descriptor, target, inheritable=True)
                metadata = os.fstat(target)
                seed_identities[label] = (
                    metadata.st_dev, metadata.st_ino
                )
            for descriptor in remap_seed_descriptors:
                os.close(descriptor)
            remap_seed_descriptors = []
            remaps = {
                210: 211,
                211: 210,
                212: 210,
                213: 211,
            }
            relocated_control = relocate_control_descriptor(
                213, (210, 211, 212), remaps
            )
            check(
                relocated_control != 213 and
                (os.fstat(relocated_control).st_dev,
                 os.fstat(relocated_control).st_ino) ==
                seed_identities["control"],
                "colliding supervisor control FD is retained elsewhere",
            )
            try:
                raise RuntimeError(
                    "self-test active remap caller exception"
                )
            except RuntimeError:
                remapped_descriptors = remap_inherited_descriptors(
                    (210, 211, 212), remaps, (relocated_control,)
                )
            remapped_identities = {
                target: (
                    os.fstat(target).st_dev, os.fstat(target).st_ino
                ) for target in remaps
            }
            backup_identities = {
                (
                    os.fstat(descriptor).st_dev,
                    os.fstat(descriptor).st_ino,
                )
                for descriptor in remapped_descriptors
                if descriptor not in remap_targets
            }
            check(
                remapped_identities == {
                    210: seed_identities["b"],
                    211: seed_identities["a"],
                    212: seed_identities["a"],
                    213: seed_identities["b"],
                } and seed_identities["c"] in backup_identities,
                "descriptor remap cycles/fan-out preserve every source and "
                "colliding inherited identity",
            )
            remap_lock_path = root / "remap-graph.lock"
            remap_lock_descriptor = os.open(
                remap_lock_path,
                os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600,
            )
            os.fchmod(remap_lock_descriptor, 0o600)
            os.dup2(remap_lock_descriptor, 214, inheritable=True)
            supervisor_remap_record = run_command(
                [
                    sys.executable, "-c",
                    "import os;"
                    "print(os.pread(210,1,0).decode(),"
                    "os.pread(211,1,0).decode(),"
                    "os.pread(214,1,0).decode())",
                ],
                root, root / "remap-graph.stdout.log",
                root / "remap-graph.stderr.log", 5,
                os.environ.copy(),
                inherited_lock_descriptor=214,
                inherited_descriptors=(210, 211, 214),
                descriptor_remaps={
                    210: 211,
                    211: 210,
                    214: 210,
                },
            )
            check(
                supervisor_remap_record["returncode"] == 0 and
                (root / "remap-graph.stdout.log").read_bytes() ==
                b"a b b\n",
                "direct supervisor applies collision-safe cyclic remaps",
            )
        finally:
            for descriptor in remapped_descriptors:
                if descriptor not in remap_targets:
                    try:
                        os.close(descriptor)
                    except OSError:
                        pass
            if relocated_control is not None:
                try:
                    os.close(relocated_control)
                except OSError:
                    pass
            if remap_lock_descriptor is not None:
                try:
                    os.close(remap_lock_descriptor)
                except OSError:
                    pass
            for descriptor in remap_seed_descriptors:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            for target in remap_targets:
                try:
                    os.close(target)
                except OSError:
                    pass
                saved = saved_remap_targets.get(target)
                if saved is not None:
                    os.dup2(saved, target)
                    os.close(saved)

        raw_directory = owned_canonical_directory(
            root / "raw-descriptor", "self-test raw directory"
        )
        try:
            raw_name = "output.json"
            raw_file = open_exclusive_owned_file(
                raw_directory, raw_name, "self-test raw output"
            )
            try:
                raw_child_path = procfd_exact_path(
                    raw_file, "self-test raw output"
                )
                raw_record = run_command(
                    [
                        sys.executable, "-c",
                        "import pathlib,sys;"
                        "pathlib.Path(sys.argv[1]).write_bytes("
                        "b'{\"ok\":true}')",
                        str(raw_child_path),
                    ],
                    root, root / "raw.stdout.log", root / "raw.stderr.log",
                    5, os.environ.copy(),
                    inherited_descriptors=(raw_file["descriptor"],),
                    descriptor_remaps={
                        RAW_OUTPUT_DESCRIPTOR: raw_file["descriptor"],
                    },
                )
                check(
                    raw_record["returncode"] == 0 and
                    read_held_regular(
                        raw_file, MAX_RAW_JSON_BYTES,
                        "self-test raw output"
                    ) == b'{"ok":true}',
                    "raw output is bound to the exact held file descriptor"
                )
            finally:
                close_log(raw_file)

            raw_victim = root / "raw-symlink-victim"
            raw_victim.write_bytes(b"victim")
            os.symlink(
                str(raw_victim), "raw-symlink.json",
                dir_fd=raw_directory["descriptor"],
            )
            try:
                open_exclusive_owned_file(
                    raw_directory, "raw-symlink.json",
                    "self-test raw symlink",
                )
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: raw symlink was accepted"
                )
            check(
                raw_victim.read_bytes() == b"victim",
                "raw symlink victim is untouched",
            )
            saved_target = None
            collision_read = None
            collision_write = None
            collision_file = None
            try:
                try:
                    saved_target = os.dup(RAW_OUTPUT_DESCRIPTOR)
                except OSError as error:
                    if error.errno != errno.EBADF:
                        raise
                collision_read, collision_write = os.pipe()
                os.dup2(collision_read, RAW_OUTPUT_DESCRIPTOR)
                collision_name = "descriptor-collision.json"
                collision_file = open_exclusive_owned_file(
                    raw_directory, collision_name,
                    "self-test colliding raw output"
                )
                collision_child_path = procfd_exact_path(
                    collision_file, "self-test colliding raw output"
                )
                collision_record = run_command(
                    [
                        sys.executable, "-c",
                        "import pathlib,sys;"
                        "pathlib.Path(sys.argv[1]).write_bytes(b'collision')",
                        str(collision_child_path),
                    ],
                    root, root / "collision-raw.stdout.log",
                    root / "collision-raw.stderr.log",
                    5, os.environ.copy(),
                    inherited_descriptors=(
                        collision_file["descriptor"],
                        RAW_OUTPUT_DESCRIPTOR,
                    ),
                    descriptor_remaps={
                        RAW_OUTPUT_DESCRIPTOR:
                            collision_file["descriptor"],
                    },
                )
                check(
                    collision_record["returncode"] == 0 and
                    read_held_regular(
                        collision_file,
                        MAX_RAW_JSON_BYTES,
                        "self-test colliding raw output"
                    ) == b"collision",
                    "raw descriptor remap preserves a colliding inherited FD"
                )
            finally:
                if collision_file is not None:
                    close_log(collision_file)
                os.close(RAW_OUTPUT_DESCRIPTOR)
                if saved_target is not None:
                    os.dup2(saved_target, RAW_OUTPUT_DESCRIPTOR)
                    os.close(saved_target)
                if (collision_read is not None and
                        collision_read != RAW_OUTPUT_DESCRIPTOR):
                    os.close(collision_read)
                if collision_write is not None:
                    os.close(collision_write)

            capped_file = open_exclusive_owned_file(
                raw_directory, "over-cap.json",
                "self-test capped raw output",
            )
            try:
                capped_record = run_command(
                    [
                        sys.executable, "-c",
                        "import os;"
                        "block=b'x'*(1024*1024);"
                        "fd=int(os.environ['LEO2_RAW_FD']);"
                        "[os.write(fd,block) for _ in range(65)]",
                    ],
                    root, root / "capped-raw.stdout.log",
                    root / "capped-raw.stderr.log",
                    10,
                    dict(os.environ, LEO2_RAW_FD=str(RAW_OUTPUT_DESCRIPTOR)),
                    inherited_descriptors=(capped_file["descriptor"],),
                    descriptor_remaps={
                        RAW_OUTPUT_DESCRIPTOR: capped_file["descriptor"],
                    },
                )
                check(
                    capped_record["returncode"] != 0 and
                    os.fstat(capped_file["descriptor"]).st_size <=
                    MAX_RAW_JSON_BYTES,
                    "kernel-enforced raw JSON write ceiling",
                )
            finally:
                close_log(capped_file)
            mutation_name = "mutating.json"
            mutation_descriptor = os.open(
                mutation_name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600, dir_fd=raw_directory["descriptor"]
            )
            try:
                os.fchmod(mutation_descriptor, 0o600)
                write_descriptor_all(
                    mutation_descriptor, b"A" * (2 * 1024 * 1024)
                )
                os.fsync(mutation_descriptor)
            finally:
                os.close(mutation_descriptor)

            def mutate_raw_during_read():
                descriptor = os.open(
                    mutation_name, os.O_WRONLY | os.O_NOFOLLOW,
                    dir_fd=raw_directory["descriptor"]
                )
                try:
                    os.pwrite(descriptor, b"B", 0)
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)

            try:
                read_owned_regular(
                    raw_directory, mutation_name, MAX_RAW_JSON_BYTES,
                    "self-test mutating raw output",
                    mutation_hook=mutate_raw_during_read
                )
            except CrossoverError as error:
                check(
                    "changed while it was being read" in str(error),
                    "raw mutation has an actionable rejection"
                )
            else:
                raise CrossoverError(
                    "self-test failed: in-place raw mutation was accepted"
                )
        finally:
            close_owned_directory(raw_directory)

        executable_result_root = owned_canonical_directory(
            root, "self-test executable result root"
        )
        executable_directory = root / "exact-executable"
        executable_directory.mkdir(mode=0o700)
        exact_executable_path = executable_directory / "codec"
        exact_executable_path.write_bytes(
            b"#!/bin/sh\nprintf 'held-executable\\n'\n"
        )
        exact_executable_path.chmod(0o555)
        exact_executable = open_result_regular_held(
            executable_result_root, Path("exact-executable") / "codec",
            MAX_RAW_JSON_BYTES, "self-test exact executable", 0o555,
        )
        saved_executable_path = executable_directory / "codec.saved"
        try:
            validate_log_identity(
                exact_executable, "self-test exact executable"
            )
            exact_executable_path.rename(saved_executable_path)
            exact_executable_path.write_bytes(
                b"#!/bin/sh\nprintf 'substituted-executable\\n'\n"
            )
            exact_executable_path.chmod(0o555)
            exact_record = run_command(
                [str(Path("/proc/self/fd") / str(EXECUTABLE_DESCRIPTOR))],
                root, root / "exact-executable.stdout.log",
                root / "exact-executable.stderr.log", 5,
                os.environ.copy(),
                inherited_descriptors=(exact_executable["descriptor"],),
                descriptor_remaps={
                    EXECUTABLE_DESCRIPTOR: exact_executable["descriptor"],
                },
            )
            check(
                exact_record["returncode"] == 0 and
                (root / "exact-executable.stdout.log").read_bytes() ==
                b"held-executable\n",
                "execution uses the exact retained executable inode",
            )
            exact_executable_path.unlink()
            saved_executable_path.rename(exact_executable_path)
            validate_log_identity(
                exact_executable, "self-test exact executable"
            )

            original_executable_bytes = read_held_regular(
                exact_executable, MAX_RAW_JSON_BYTES,
                "self-test exact executable",
            )
            execution_snapshot = create_sealed_executable_snapshot(
                original_executable_bytes,
                "self-test exact executable",
            )
            source_identity_before = exact_executable["identity"]
            try:
                exact_executable_path.chmod(0o755)
                mutation_descriptor = os.open(
                    exact_executable_path,
                    os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW,
                )
                try:
                    write_descriptor_all(
                        mutation_descriptor,
                        b"#!/bin/sh\nprintf 'evil-executable\\n'\n",
                    )
                    os.fsync(mutation_descriptor)
                finally:
                    os.close(mutation_descriptor)
                exact_executable_path.chmod(0o555)
                snapshot_record = run_command(
                    [
                        str(Path("/proc/self/fd") /
                            str(EXECUTABLE_DESCRIPTOR))
                    ],
                    root, root / "sealed-executable.stdout.log",
                    root / "sealed-executable.stderr.log", 5,
                    os.environ.copy(),
                    inherited_descriptors=(
                        execution_snapshot["descriptor"],
                    ),
                    descriptor_remaps={
                        EXECUTABLE_DESCRIPTOR:
                            execution_snapshot["descriptor"],
                    },
                )
                exact_executable_path.chmod(0o755)
                restoration_descriptor = os.open(
                    exact_executable_path,
                    os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW,
                )
                try:
                    write_descriptor_all(
                        restoration_descriptor,
                        original_executable_bytes,
                    )
                    os.fsync(restoration_descriptor)
                finally:
                    os.close(restoration_descriptor)
                exact_executable_path.chmod(0o555)
                validate_sealed_executable_snapshot(
                    execution_snapshot,
                    digest_bytes(original_executable_bytes),
                    "self-test exact executable snapshot",
                )
                check(
                    snapshot_record["returncode"] == 0 and
                    (root / "sealed-executable.stdout.log").read_bytes() ==
                    b"held-executable\n" and
                    exact_executable["identity"] ==
                    source_identity_before and
                    read_held_regular(
                        exact_executable, MAX_RAW_JSON_BYTES,
                        "self-test restored exact executable",
                    ) == original_executable_bytes,
                    "sealed execution resists same-inode A-B-A source "
                    "mutation at exec",
                )
            finally:
                if exact_executable_path.exists():
                    exact_executable_path.chmod(0o755)
                    restoration_descriptor = os.open(
                        exact_executable_path,
                        os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW,
                    )
                    try:
                        write_descriptor_all(
                            restoration_descriptor,
                            original_executable_bytes,
                        )
                    finally:
                        os.close(restoration_descriptor)
                    exact_executable_path.chmod(0o555)
                close_sealed_snapshot(execution_snapshot)
        finally:
            if saved_executable_path.exists():
                if exact_executable_path.exists():
                    exact_executable_path.unlink()
                saved_executable_path.rename(exact_executable_path)
            close_log(exact_executable)
            close_owned_directory(executable_result_root)

        retained_git_path = root / "retained-git"
        retained_git_a = b"#!/bin/sh\nprintf 'retained-git-a\\n'\n"
        retained_git_b = b"#!/bin/sh\nprintf 'retained-git-b\\n'\n"
        retained_git_path.write_bytes(retained_git_a)
        retained_git_path.chmod(0o755)
        original_canonical_git = globals()["CANONICAL_GIT"]
        retained_git_snapshot = None
        try:
            globals()["CANONICAL_GIT"] = retained_git_path
            retained_git_snapshot = open_canonical_git_snapshot()
            retained_git_path.write_bytes(retained_git_b)
            retained_git_path.chmod(0o755)
            retained_git_output = git_command_bytes(
                root, ("--version",),
                "self-test exact retained Git", retained_git_snapshot,
            )
            retained_git_path.write_bytes(retained_git_a)
            retained_git_path.chmod(0o755)
            check(
                retained_git_output == b"retained-git-a\n" and
                validate_sealed_executable_snapshot(
                    retained_git_snapshot,
                    digest_bytes(retained_git_a),
                    "self-test retained Git snapshot",
                ) == retained_git_a,
                "Git commands use one sealed exact-tool snapshot across "
                "same-inode A-B-A mutation",
            )
        finally:
            globals()["CANONICAL_GIT"] = original_canonical_git
            if retained_git_path.exists():
                retained_git_path.write_bytes(retained_git_a)
                retained_git_path.chmod(0o755)
            close_canonical_git_snapshot(retained_git_snapshot)

        launcher_script = root / "launcher-script.py"
        launcher_a = b"print('launcher-a')\n"
        launcher_b = b"print('launcher-b')\n"
        launcher_script.write_bytes(launcher_a)
        launcher_script.chmod(0o755)
        unrelated_python = root / "unrelated-python"
        unrelated_python.write_bytes(b"#!/bin/sh\nexit 0\n")
        unrelated_python.chmod(0o755)
        try:
            open_runtime_launcher_snapshots(
                python_path=unrelated_python, script_path=launcher_script
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: unrelated Python launcher inode was "
                "accepted"
            )
        launcher_snapshots = open_runtime_launcher_snapshots(
            python_path=sys.executable, script_path=launcher_script
        )
        try:
            launcher_script.write_bytes(launcher_b)
            launcher_script.chmod(0o755)
            launcher_process = subprocess.Popen(
                [
                    launcher_snapshots["identity"]["python"]["path"],
                    *RUNTIME_PYTHON_FLAGS,
                    str(Path("/proc/self/fd") / str(
                        launcher_snapshots["script"]["descriptor"]
                    )),
                ],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                pass_fds=(
                    launcher_snapshots["python"]["descriptor"],
                    launcher_snapshots["script"]["descriptor"],
                ),
                executable=str(Path("/proc/self/fd") / str(
                    launcher_snapshots["python"]["descriptor"]
                )),
            )
            launcher_stdout, launcher_stderr = launcher_process.communicate(
                timeout=10
            )
            launcher_script.write_bytes(launcher_a)
            launcher_script.chmod(0o755)
            check(
                launcher_process.returncode == 0 and
                launcher_stdout == b"launcher-a\n" and
                launcher_stderr == b"" and
                validate_runtime_launcher_snapshots(
                    launcher_snapshots
                ) == launcher_snapshots["identity"],
                "Python/runner re-exec uses sealed snapshots across a "
                "same-inode A-B-A mutation",
            )
        finally:
            if launcher_script.exists():
                launcher_script.write_bytes(launcher_a)
                launcher_script.chmod(0o755)
            close_runtime_launcher_snapshots(launcher_snapshots)

        sealed_tools = []
        try:
            for label, descriptor_target in (
                    ("taskset", TASKSET_EXECUTABLE_DESCRIPTOR),
                    ("cmake", CONTROLLED_CMAKE_DESCRIPTOR),
                    ("ninja", CONTROLLED_NINJA_DESCRIPTOR)):
                tool_path = root / ("sealed-" + label)
                tool_a = (
                    "#!/bin/sh\nprintf '{}-a\\\\n'\n".format(label)
                ).encode("ascii")
                tool_b = (
                    "#!/bin/sh\nprintf '{}-b\\\\n'\n".format(label)
                ).encode("ascii")
                tool_path.write_bytes(tool_a)
                tool_path.chmod(0o755)
                snapshot = open_exact_executable_snapshot(
                    tool_path, "self-test sealed " + label, True
                )
                sealed_tools.append(snapshot)
                tool_path.write_bytes(tool_b)
                tool_path.chmod(0o755)
                record = run_command(
                    [
                        str(Path("/proc/self/fd") /
                            str(descriptor_target))
                    ],
                    root,
                    root / ("sealed-" + label + ".stdout.log"),
                    root / ("sealed-" + label + ".stderr.log"),
                    5, os.environ.copy(),
                    inherited_descriptors=(snapshot["descriptor"],),
                    descriptor_remaps={
                        descriptor_target: snapshot["descriptor"],
                    },
                )
                tool_path.write_bytes(tool_a)
                tool_path.chmod(0o755)
                check(
                    record["returncode"] == 0 and
                    (root / (
                        "sealed-" + label + ".stdout.log"
                    )).read_bytes() ==
                    ("{}-a\\n".format(label)).encode("ascii") and
                    validate_sealed_executable_snapshot(
                        snapshot, digest_bytes(tool_a),
                        "self-test sealed " + label,
                    ) == tool_a,
                    "{} execution uses a sealed descriptor across a "
                    "same-inode A-B-A mutation".format(label),
                )
        finally:
            for snapshot in reversed(sealed_tools):
                close_sealed_snapshot(snapshot)

        fake_git = root / "fake-git"
        fake_git_child_pid = root / "fake-git-child.pid"
        fake_git.write_text(
            "#!/usr/bin/python3\n"
            "import os,pathlib,signal,time\n"
            "pid=os.fork()\n"
            "if pid==0:\n"
            " os.setsid()\n"
            " target=pathlib.Path({!r})\n"
            " tmp=target.with_suffix('.tmp')\n"
            " tmp.write_text(str(os.getpid()),encoding='ascii')\n"
            " os.replace(str(tmp),str(target))\n"
            " signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
            " time.sleep(30)\n"
            " os._exit(0)\n"
            "deadline=time.monotonic()+5\n"
            "while not pathlib.Path({!r}).exists() and "
            "time.monotonic()<deadline: time.sleep(.01)\n"
            "raise SystemExit(0)\n".format(
                str(fake_git_child_pid), str(fake_git_child_pid)
            ),
            encoding="utf-8",
        )
        fake_git.chmod(0o700)
        original_canonical_git = globals()["CANONICAL_GIT"]
        try:
            globals()["CANONICAL_GIT"] = fake_git
            git_command_bytes(root, ("status",), "self-test bounded Git")
        finally:
            globals()["CANONICAL_GIT"] = original_canonical_git
        fake_git_pid = int(fake_git_child_pid.read_text(encoding="ascii"))
        check(
            proc_record(fake_git_pid) is None,
            "Git provenance helper descendant is killed and reaped",
        )

        resume_root = root / "resumable-job"
        resume_log_root = resume_root / "logs" / "resume-fixture"
        resume_raw_root = resume_root / "raw" / "resume-fixture"
        resume_job_root = resume_root / "jobs"
        for directory_path in (
                resume_log_root, resume_raw_root, resume_job_root):
            directory_path.mkdir(parents=True, mode=0o700)
            directory_path.chmod(0o700)
        resume_root.chmod(0o700)
        (resume_root / "logs").chmod(0o700)
        (resume_root / "raw").chmod(0o700)
        resume_order = ["direct", "transform"]
        for index, mode in enumerate(resume_order):
            label = "{:02d}-{}".format(index, mode)
            (resume_log_root / (label + ".stdout.log")).write_bytes(
                b"stale stdout"
            )
            (resume_log_root / (label + ".stderr.log")).write_bytes(
                b"stale stderr"
            )
            (resume_raw_root / (label + ".json")).write_bytes(
                b'{"stale":true}'
            )
            for stale_path in (
                    resume_log_root / (label + ".stdout.log"),
                    resume_log_root / (label + ".stderr.log"),
                    resume_raw_root / (label + ".json")):
                stale_path.chmod(0o600)
        atomic_write_json(
            resume_job_root / "resume-fixture.json",
            {
                "configuration_id": "resume-configuration",
                "schema": JOB_SCHEMA,
                "status": "failed",
            }
        )
        resume_job = {
            "build_metadata": {},
            "cell": {"backend": "scalar"},
            "configuration_id": "resume-configuration",
            "executable": "/self-test/executable",
            "executable_artifact": None,
            "executable_sha256": "a" * 64,
            "invocation_order": resume_order,
            "job_id": "resume-fixture",
            "seed": 7,
            "source_identity": {},
        }
        resume_context = {
            "isolation": None,
            "result_dir": resume_root,
            "resume": True,
            "settings": {},
            "source": root,
            "timeout": 5,
        }
        original_benchmark_argv = globals()["benchmark_argv"]
        original_validate_execution_inputs = globals()[
            "validate_execution_inputs"
        ]
        original_validate_raw = globals()["validate_raw"]
        original_summarize_measurements = globals()[
            "summarize_measurements"
        ]

        def resume_benchmark_argv(
                unused_job, unused_mode, output_path, unused_settings):
            return [
                sys.executable, "-c",
                "import pathlib,sys;"
                "pathlib.Path(sys.argv[1]).write_bytes(b'{}')",
                str(output_path),
            ]

        def resume_validate_raw(
                unused_raw, unused_job, unused_mode, unused_settings):
            return 1.0, 0.0, {
                "algorithm": "fixture",
                "digest": "stable",
                "hashed_bytes": 1,
                "requested_parity_indices": [0],
            }

        try:
            globals()["benchmark_argv"] = resume_benchmark_argv
            globals()["validate_execution_inputs"] = (
                lambda unused_job, unused_root=None: None
            )
            globals()["validate_raw"] = resume_validate_raw
            globals()["summarize_measurements"] = (
                lambda unused_measurements: {"fixture": "passed"}
            )
            resumed_result = run_job(resume_job, resume_context)
        finally:
            globals()["benchmark_argv"] = original_benchmark_argv
            globals()["validate_execution_inputs"] = (
                original_validate_execution_inputs
            )
            globals()["validate_raw"] = original_validate_raw
            globals()["summarize_measurements"] = (
                original_summarize_measurements
            )
        check(
            resumed_result["status"] == "passed" and
            all(
                (resume_raw_root / "{:02d}-{}.json".format(
                    index, mode
                )).read_bytes() == b"{}"
                for index, mode in enumerate(resume_order)
            ) and
            all(
                (resume_log_root / "{:02d}-{}.stdout.log".format(
                    index, mode
                )).read_bytes() == b""
                for index, mode in enumerate(resume_order)
            ),
            "partial job artifacts are safely cleaned and resumable: "
            "{}".format(resumed_result.get("reason"))
        )

        def live_descriptor_count():
            count = 0
            for name in os.listdir("/proc/self/fd"):
                try:
                    descriptor = int(name)
                    os.fstat(descriptor)
                except (OSError, ValueError):
                    continue
                count += 1
            return count

        class AbortingIsolationSupport(object):
            @staticmethod
            def cpu_stat_snapshot(unused_cpu):
                raise InjectedAbort("injected isolation snapshot abort")

        abort_job = dict(resume_job)
        abort_job["configuration_id"] = "abort-configuration"
        abort_job["job_id"] = "abort-fixture"
        abort_context = dict(resume_context)
        abort_context["resume"] = False
        abort_context["settings"] = {"mode": "screen"}
        abort_context["isolation"] = {
            "canonical_guard": object(),
            "cpu": 0,
            "pair_guard": object(),
            "pair_lease": {},
            "sibling": 1,
            "support": AbortingIsolationSupport(),
        }
        original_validate_authoritative_guards = globals()[
            "validate_authoritative_guards"
        ]
        descriptors_before_abort = live_descriptor_count()
        try:
            globals()["validate_authoritative_guards"] = (
                lambda unused_canonical, unused_pair: None
            )
            try:
                run_job(abort_job, abort_context)
            except InjectedAbort:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: isolation BaseException was ignored"
                )
        finally:
            globals()["validate_authoritative_guards"] = (
                original_validate_authoritative_guards
            )
        check(
            live_descriptor_count() == descriptors_before_abort,
            "run_job closes result/log/raw directories after a "
            "BaseException during isolation",
        )

        escaped_pid_path = root / "closed-fd-escape.pid"
        escape_program = (
            "import os,pathlib,signal,sys,time\n"
            "first=os.fork()\n"
            "if first==0:\n"
            " os.setsid()\n"
            " daemon=os.fork()\n"
            " if daemon!=0: os._exit(0)\n"
            " null=os.open('/dev/null',os.O_RDWR)\n"
            " os.dup2(null,0);os.dup2(null,1);os.dup2(null,2)\n"
            " if null>2: os.close(null)\n"
            " for fd in range(3,256):\n"
            "  try: os.close(fd)\n"
            "  except OSError: pass\n"
            " pathlib.Path(sys.argv[1]).write_text("
            "str(os.getpid()),encoding='ascii')\n"
            " signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
            " time.sleep(30)\n"
            "deadline=time.monotonic()+5\n"
            "while not os.path.exists(sys.argv[1]) and "
            "time.monotonic()<deadline: time.sleep(.01)\n"
            "sys.exit(0 if os.path.exists(sys.argv[1]) else 93)\n"
        )
        escaped_record = run_command(
            [
                sys.executable, "-c", escape_program,
                str(escaped_pid_path),
            ],
            root, root / "escape.stdout.log", root / "escape.stderr.log",
            5, os.environ.copy()
        )
        escaped_pid = int(escaped_pid_path.read_text(encoding="ascii"))
        check(
            escaped_record["returncode"] == 0 and
            proc_record(escaped_pid) is None,
            "closed-descriptor double-fork escape is killed and reaped"
        )

        killed_supervisor_pid_path = root / "killed-supervisor-escape.pid"
        killed_supervisor_program = (
            "import os,pathlib,signal,sys,time\n"
            "first=os.fork()\n"
            "if first==0:\n"
            " os.setsid()\n"
            " daemon=os.fork()\n"
            " if daemon!=0: os._exit(0)\n"
            " null=os.open('/dev/null',os.O_RDWR)\n"
            " os.dup2(null,0);os.dup2(null,1);os.dup2(null,2)\n"
            " if null>2: os.close(null)\n"
            " for fd in range(3,256):\n"
            "  try: os.close(fd)\n"
            "  except OSError: pass\n"
            " pathlib.Path(sys.argv[1]).write_text("
            "str(os.getpid()),encoding='ascii')\n"
            " signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
            " time.sleep(30)\n"
            "deadline=time.monotonic()+5\n"
            "while not os.path.exists(sys.argv[1]) and "
            "time.monotonic()<deadline: time.sleep(.01)\n"
            "if not os.path.exists(sys.argv[1]): sys.exit(93)\n"
            "os.kill(os.getppid(),signal.SIGKILL)\n"
            "time.sleep(.05)\n"
            "os._exit(0)\n"
        )
        try:
            run_command(
                [
                    sys.executable, "-c", killed_supervisor_program,
                    str(killed_supervisor_pid_path),
                ],
                root,
                root / "killed-supervisor.stdout.log",
                root / "killed-supervisor.stderr.log",
                5, os.environ.copy()
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: killed command supervisor was accepted"
            )
        killed_supervisor_escape_pid = int(
            killed_supervisor_pid_path.read_text(encoding="ascii")
        )
        check(
            proc_record(killed_supervisor_escape_pid) is None,
            "outer owner reaps a closed-FD escape after supervisor SIGKILL"
        )

        timeout_pid_path = root / "timeout-tree.pids"
        timeout_child_program = (
            "import signal,time; "
            "signal.signal(signal.SIGTERM,signal.SIG_IGN); "
            "time.sleep(30)"
        )
        timeout_program = (
            "import os,pathlib,signal,subprocess,sys,time\n"
            "child=subprocess.Popen("
            "[sys.executable,'-c',{!r}],start_new_session=True)\n"
            "pathlib.Path({!r}).write_text("
            "str(os.getpid())+' '+str(child.pid),encoding='ascii')\n"
            "signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
            "time.sleep(30)\n"
        ).format(timeout_child_program, str(timeout_pid_path))
        timeout_started = time.monotonic()
        timeout_record = run_command(
            [sys.executable, "-c", timeout_program], root,
            root / "timeout.stdout.log", root / "timeout.stderr.log",
            2, os.environ.copy()
        )
        timeout_elapsed = time.monotonic() - timeout_started
        check(
            timeout_record["timed_out"] is True and
            timeout_record["returncode"] == 124 and
            timeout_elapsed < 6,
            "timed-out process tree is terminated in bounded time"
        )
        check(
            timeout_pid_path.is_file(),
            "timed-out process tree recorded its identities"
        )
        timeout_pids = [
            int(value) for value in timeout_pid_path.read_text(
                encoding="ascii"
            ).split()
        ]
        check(
            len(timeout_pids) == 2 and
            all(proc_record(pid) is None for pid in timeout_pids),
            "timed-out leader and detached SIGTERM-ignoring child are reaped"
        )

        collision_stdout = root / "collision.stdout.log"
        collision_stderr = root / "collision.stderr.log"
        collision_stdout.write_bytes(b"pre-existing stdout\n")
        try:
            run_command(
                [sys.executable, "-c", "print('must not run')"], root,
                collision_stdout, collision_stderr, 5, os.environ.copy()
            )
        except CrossoverError as error:
            check(
                "already exists" in str(error),
                "pre-existing stdout collision has an actionable rejection"
            )
        else:
            raise CrossoverError(
                "self-test failed: pre-existing stdout log was replaced"
            )
        check(
            collision_stdout.read_bytes() == b"pre-existing stdout\n" and
            not collision_stderr.exists(),
            "pre-existing stdout collision is untouched"
        )

        collision_stdout.unlink()
        collision_stderr.write_bytes(b"pre-existing stderr\n")
        try:
            run_command(
                [sys.executable, "-c", "print('must not run')"], root,
                collision_stdout, collision_stderr, 5, os.environ.copy()
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: pre-existing stderr log was replaced"
            )
        check(
            not collision_stdout.exists() and
            collision_stderr.read_bytes() == b"pre-existing stderr\n",
            "stderr collision rolls back its newly reserved stdout log"
        )

        collision_stderr.unlink()
        collision_victim = root / "collision-victim.txt"
        collision_victim.write_bytes(b"symlink victim\n")
        collision_stdout.symlink_to(collision_victim)
        try:
            run_command(
                [sys.executable, "-c", "print('must not run')"], root,
                collision_stdout, collision_stderr, 5, os.environ.copy()
            )
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: symbolic-link stdout log was followed"
            )
        check(
            collision_stdout.is_symlink() and
            collision_victim.read_bytes() == b"symlink victim\n" and
            not collision_stderr.exists(),
            "symbolic-link log collision is untouched"
        )

        inherited_lock_path = root / "coordinator-exit.lock"
        coordinator_root = root / "coordinator-exit"
        coordinator_root.mkdir(mode=0o700)
        coordinator_root.chmod(0o700)
        child_started_path = coordinator_root / "child-started.pid"
        child_finished_path = coordinator_root / "child-finished.txt"
        lock_child_program = (
            "import os,pathlib,time; "
            "started=pathlib.Path({!r}); tmp=started.with_suffix('.tmp'); "
            "tmp.write_text(str(os.getpid()),encoding='ascii'); "
            "os.replace(str(tmp),str(started)); time.sleep(1.25); "
            "finished=pathlib.Path({!r}); tmp=finished.with_suffix('.tmp'); "
            "tmp.write_text('finished',encoding='ascii'); "
            "os.replace(str(tmp),str(finished))"
        ).format(str(child_started_path), str(child_finished_path))
        module_path = str(Path(__file__).resolve())
        coordinator_program = (
            "import fcntl,importlib.util,os,sys; "
            "spec=importlib.util.spec_from_file_location("
            "'leo2_lock_lifetime_fixture',{!r}); "
            "module=importlib.util.module_from_spec(spec); "
            "sys.modules[spec.name]=module; spec.loader.exec_module(module); "
            "lock=os.open({!r},os.O_RDWR|os.O_CREAT,0o600); "
            "fcntl.flock(lock,fcntl.LOCK_EX); "
            "module.run_command([sys.executable,'-c',{!r}],{!r},"
            "{!r},{!r},10,os.environ.copy(),"
            "inherited_lock_descriptor=lock)"
        ).format(
            module_path, str(inherited_lock_path), lock_child_program,
            str(coordinator_root),
            str(coordinator_root / "stdout.log"),
            str(coordinator_root / "stderr.log"),
        )
        coordinator = None
        try:
            coordinator = subprocess.Popen(
                [sys.executable, "-c", coordinator_program],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, start_new_session=True,
            )
            inherited_child_pid = None
            child_started_deadline = time.monotonic() + 5
            while time.monotonic() < child_started_deadline:
                if coordinator.poll() is not None:
                    break
                try:
                    candidate_pid = int(child_started_path.read_text(
                        encoding="ascii"
                    ))
                    if candidate_pid > 0:
                        inherited_child_pid = candidate_pid
                        break
                except (FileNotFoundError, OSError, UnicodeError, ValueError):
                    pass
                time.sleep(0.01)
            check(
                inherited_child_pid is not None and
                coordinator.poll() is None,
                "coordinator-exit fixture launched its timed child"
            )
            coordinator.kill()
            coordinator.wait(timeout=5)
        finally:
            if coordinator is not None:
                if coordinator.poll() is None:
                    coordinator.kill()
                try:
                    coordinator.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    coordinator.kill()
                    coordinator.wait(timeout=5)

        def lock_can_be_acquired():
            descriptor = os.open(
                str(inherited_lock_path), os.O_RDWR | os.O_CLOEXEC
            )
            try:
                try:
                    fcntl.flock(
                        descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB
                    )
                except BlockingIOError:
                    return False
                return True
            finally:
                # Close-only lifetime is the contract under test.
                os.close(descriptor)

        check(
            not lock_can_be_acquired(),
            "coordinator SIGKILL does not release a timed child's lock"
        )
        child_finished_deadline = time.monotonic() + 5
        while (not child_finished_path.is_file() and
               time.monotonic() < child_finished_deadline):
            time.sleep(0.01)
        check(
            child_finished_path.is_file(),
            "inherited-lock child completed after coordinator SIGKILL"
        )
        lock_release_deadline = time.monotonic() + 5
        while time.monotonic() < lock_release_deadline:
            if lock_can_be_acquired():
                break
            time.sleep(0.01)
        else:
            raise CrossoverError(
                "self-test failed: inherited lock was not released after "
                "the timed child closed its descriptor"
            )
        try:
            os.waitpid(inherited_child_pid, 0)
        except ChildProcessError:
            pass

        result_root = root / "results"
        manifest_path = result_root / "manifest.json"
        for label, invalid_manifest in (
                ("v1 manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v1",
                    "settings": {},
                }),
                ("v2 manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v2",
                    "settings": {},
                }),
                ("v3 manifest without small-direct mode binding", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v3",
                    "settings": {},
                }),
                ("v4 directory-fd raw-output manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v4",
                    "settings": {},
                }),
                ("v5 mutable-launcher manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v5",
                    "settings": {},
                }),
                ("relabeled v1 manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema": SCHEMA, "settings": {},
                })):
            atomic_write_json(manifest_path, invalid_manifest)
            try:
                load_manifest(result_root)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} was accepted".format(label)
                )
        submodule = root / "submodule"
        superproject = root / "superproject"

        def fixture_git(cwd, arguments):
            completed = subprocess.run(
                [str(CANONICAL_GIT)] + list(arguments), cwd=str(cwd),
                env=GIT_ENVIRONMENT,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=20, check=False
            )
            if completed.returncode != 0:
                raise CrossoverError(
                    "self-test git fixture failed: {}".format(
                        normalized_output(completed.stderr).decode(
                            "utf-8", errors="replace"
                        ).strip()
                    )
                )
            return normalized_output(completed.stdout)

        submodule.mkdir()
        fixture_git(submodule, ("init", "-q"))
        fixture_git(submodule, ("config", "user.name", "Self Test"))
        fixture_git(submodule, ("config", "user.email", "self@test.invalid"))
        (submodule / "tracked.txt").write_text("clean\n", encoding="utf-8")
        fixture_git(submodule, ("add", "tracked.txt"))
        fixture_git(submodule, ("commit", "-q", "-m", "initial"))
        (submodule / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )
        fixture_git(submodule, ("commit", "-q", "-am", "current"))
        superproject.mkdir()
        fixture_git(superproject, ("init", "-q"))
        fixture_git(superproject, ("config", "user.name", "Self Test"))
        fixture_git(superproject, ("config", "user.email", "self@test.invalid"))
        fixture_git(superproject, (
            "-c", "protocol.file.allow=always", "submodule", "add", "-q",
            str(submodule.resolve()), "module",
        ))
        fixture_git(superproject, ("commit", "-q", "-am", "gitlink"))
        (superproject / "module" / "tracked.txt").write_text(
            "worktree-only change\n", encoding="utf-8"
        )
        fixture_state = source_identity(superproject, require_clean=False)
        validate_source_state(fixture_state, "self-test gitlink fixture", False)
        fixture_submodule = fixture_state["files"]["module"]["submodule"]
        check(
            fixture_submodule["worktree_clean"] is False and
            fixture_submodule["status"] == [" M tracked.txt"],
            "gitlink retains exact worktree-only porcelain status"
        )
        module = superproject / "module"
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )
        clean_fixture = source_identity(
            superproject, require_clean=True
        )

        def reject_hidden_source(label):
            try:
                source_identity(superproject, require_clean=True)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} was accepted".format(label)
                )

        fixture_git(module, (
            "update-index", "--assume-unchanged", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "assume-unchanged mutation\n", encoding="utf-8"
        )
        reject_hidden_source("assume-unchanged dirty submodule")
        fixture_git(module, (
            "update-index", "--no-assume-unchanged", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )

        fixture_git(module, (
            "update-index", "--skip-worktree", "tracked.txt",
        ))
        (module / "tracked.txt").unlink()
        reject_hidden_source("skip-worktree missing submodule file")
        fixture_git(module, (
            "update-index", "--no-skip-worktree", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )

        fixture_git(module, ("config", "core.filemode", "false"))
        (module / "tracked.txt").chmod(0o755)
        reject_hidden_source("status-hidden executable-bit mutation")
        (module / "tracked.txt").chmod(0o644)

        fixture_git(module, ("replace", "HEAD", "HEAD^"))
        reject_hidden_source("Git replacement ref")
        fixture_git(module, ("replace", "-d", "HEAD"))

        hostile_config = root / "hostile.gitconfig"
        hostile_config.write_text(
            "[core]\n\tbare = true\n", encoding="utf-8"
        )
        inherited_names = ("GIT_CONFIG_GLOBAL", "HOME", "PATH")
        inherited_before = {
            name: os.environ.get(name) for name in inherited_names
        }
        os.environ["GIT_CONFIG_GLOBAL"] = str(hostile_config)
        os.environ["HOME"] = str(root / "hostile-home")
        os.environ["PATH"] = str(root)
        try:
            hostile_fixture = source_identity(
                superproject, require_clean=True
            )
        finally:
            for name, old_value in inherited_before.items():
                if old_value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = old_value
        check(
            canonical_bytes(hostile_fixture) ==
                canonical_bytes(clean_fixture),
            "hostile inherited Git config or PATH changed source identity"
        )

        race_repository = root / "snapshot-race-repository"
        race_repository.mkdir()
        fixture_git(race_repository, ("init", "-q"))
        fixture_git(
            race_repository, ("config", "user.name", "Self Test")
        )
        fixture_git(
            race_repository,
            ("config", "user.email", "self@test.invalid"),
        )
        (race_repository / "tracked.txt").write_text(
            "snapshot a\n", encoding="utf-8"
        )
        fixture_git(race_repository, ("add", "tracked.txt"))
        fixture_git(
            race_repository, ("commit", "-q", "-m", "snapshot a")
        )
        race_commit_a = fixture_git(
            race_repository, ("rev-parse", "HEAD")
        ).decode("ascii").strip()
        (race_repository / "tracked.txt").write_text(
            "snapshot b\n", encoding="utf-8"
        )
        fixture_git(
            race_repository, ("commit", "-q", "-am", "snapshot b")
        )
        race_commit_b = fixture_git(
            race_repository, ("rev-parse", "HEAD")
        ).decode("ascii").strip()
        fixture_git(
            race_repository, ("reset", "--hard", "-q", race_commit_a)
        )
        race_worktree = root / "snapshot-race-worktree"
        fixture_git(race_repository, (
            "worktree", "add", "-q", "-b", "snapshot-race-worktree",
            str(race_worktree), race_commit_a,
        ))

        def reject_cross_bound_snapshot(
                capture_root, mutation_root, trigger_root,
                baseline_commit, raced_commit, label):
            original_git_command = globals()["git_command_bytes"]
            triggered = [False]

            def racing_git_command(
                    source_value, arguments, description, git_tool=None):
                if (not triggered[0] and
                        Path(source_value).resolve() ==
                        Path(trigger_root).resolve() and
                        tuple(arguments) ==
                        ("ls-files", "-s", "-z")):
                    triggered[0] = True
                    fixture_git(
                        mutation_root,
                        ("reset", "--hard", "-q", raced_commit),
                    )
                return original_git_command(
                    source_value, arguments, description, git_tool
                )

            try:
                globals()["git_command_bytes"] = racing_git_command
                try:
                    source_identity(capture_root, require_clean=False)
                except CrossoverError:
                    pass
                else:
                    raise CrossoverError(
                        "self-test failed: {} source-snapshot race was "
                        "accepted".format(label)
                    )
                check(
                    triggered[0],
                    "{} source-snapshot race fixture executed".format(label),
                )
            finally:
                globals()["git_command_bytes"] = original_git_command
                fixture_git(
                    mutation_root,
                    ("reset", "--hard", "-q", baseline_commit),
                )

        reject_cross_bound_snapshot(
            race_repository, race_repository, race_repository,
            race_commit_a, race_commit_b, "ordinary repository",
        )
        reject_cross_bound_snapshot(
            race_worktree, race_worktree, race_worktree,
            race_commit_a, race_commit_b, "linked worktree",
        )
        module_baseline = fixture_git(
            module, ("rev-parse", "HEAD")
        ).decode("ascii").strip()
        module_prior = fixture_git(
            module, ("rev-parse", "HEAD^")
        ).decode("ascii").strip()
        reject_cross_bound_snapshot(
            superproject, module, module,
            module_baseline, module_prior, "submodule",
        )

        fixture_git(module, ("config", "user.name", "Self Test"))
        fixture_git(module, ("config", "user.email", "self@test.invalid"))
        (module / "tracked.txt").write_text(
            "different submodule commit\n", encoding="utf-8"
        )
        fixture_git(module, ("commit", "-q", "-am", "different head"))
        try:
            source_identity(superproject, require_clean=False)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: mismatched gitlink HEAD was accepted"
            )

        fixture_git(superproject, (
            "submodule", "deinit", "-q", "-f", "--", "module",
        ))
        try:
            source_identity(superproject, require_clean=False)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: uninitialized gitlink was accepted"
            )
    for invalid in ("1", "4"):
        try:
            exactly_three_abba_rounds(invalid)
        except argparse.ArgumentTypeError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: authoritative ABBA rounds {} accepted".format(
                    invalid
                )
            )
    print(
        "leopard2 direct-encode crossover self-test passed "
        "(strict JSON/schema, provenance, parity, metric, summary, "
        "Git-closure, ABBA mutations, process-tree cleanup, inherited-lock "
        "lifetime, and exclusive-log collisions; no codec required)"
    )
    return 0


def finite_percentage(text):
    try:
        value = float(text)
    except (OverflowError, TypeError, ValueError):
        raise argparse.ArgumentTypeError("must be a finite nonnegative number")
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be a finite nonnegative number")
    return value


def exactly_three_abba_rounds(text):
    try:
        value = int(text)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            "authoritative v2 requires exactly 3 ABBA rounds"
        )
    if value != 3:
        raise argparse.ArgumentTypeError(
            "authoritative v2 requires exactly 3 ABBA rounds"
        )
    return value


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
    parser.add_argument("--promotion-percent", type=finite_percentage, default=5.0)
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
    pinned_help = (
        "Disabled legacy command: external artifacts cannot prove clean "
        "source/object/archive/link closure. Use screen for diagnostics or "
        "historical-avx2 for authoritative evidence."
    )
    pinned = subparsers.add_parser(
        "pinned", help=pinned_help, description=pinned_help,
    )
    add_run_arguments(
        pinned, "results/leopard2/direct-encode-crossover/pinned"
    )
    pinned.set_defaults(iterations=15, warmups=4, reuse=64, workers=1)
    pinned.add_argument("--cpu", type=int, default=None)
    pinned.add_argument("--sibling", type=int, default=None)
    pinned.add_argument("--taskset", default="taskset")
    pinned.add_argument(
        "--abba-rounds", type=exactly_three_abba_rounds, default=3
    )
    historical = subparsers.add_parser(
        "historical-avx2",
        help=(
            "run the corrected high-profile cells with a frozen explicit-AVX2 "
            "binary and pinned ABBA order"
        ),
    )
    add_run_arguments(
        historical,
        "results/leopard2/direct-encode-crossover/historical-avx2",
    )
    historical.set_defaults(
        backends="avx2", iterations=15, warmups=4, reuse=64, workers=1
    )
    historical.add_argument("--cpu", type=int, default=None)
    historical.add_argument("--sibling", type=int, default=None)
    historical.add_argument("--taskset", default="taskset")
    historical.add_argument(
        "--abba-rounds", type=exactly_three_abba_rounds, default=3
    )
    analyze = subparsers.add_parser(
        "analyze", help="deterministically reanalyze completed job JSON"
    )
    analyze.add_argument(
        "--result-dir",
        default="results/leopard2/direct-encode-crossover/historical-avx2",
    )
    analyze.add_argument(
        "--promotion-percent", type=finite_percentage, default=5.0
    )
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
        if arguments.command in ("screen", "pinned", "historical-avx2"):
            numeric = (
                arguments.batch, arguments.reuse, arguments.iterations,
                arguments.warmups, arguments.timeout, arguments.workers,
            )
            if any(value <= 0 for value in numeric):
                raise CrossoverError(
                    "batch, reuse, iterations, warmups, timeout, and workers "
                    "must be positive"
                )
            required_finite_number(
                arguments.promotion_percent, "promotion_percent"
            )
            if (arguments.command in AUTHORITATIVE_COMMANDS and
                    arguments.abba_rounds != 3):
                raise CrossoverError(
                    "authoritative paired-log inference requires --abba-rounds 3"
                )
            return run_matrix(arguments)
        parser().print_help()
        return 2
    except (CrossoverError, OSError, subprocess.SubprocessError) as error:
        print("direct-encode crossover error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    if (len(sys.argv) >= 2 and
            sys.argv[1] == DIRECT_COMMAND_OWNER_MODE):
        sys.exit(direct_command_owner(sys.argv[2:]))
    if (len(sys.argv) >= 2 and
            sys.argv[1] == DIRECT_COMMAND_SUPERVISOR_MODE):
        sys.exit(direct_command_supervisor(sys.argv[2:]))
    sys.exit(main())
