#!/usr/bin/env python3
"""Acquisition primitives for the sealed exact-main Leopard1 baseline.

The evidence schemas live in :mod:`leopard2_exact_main_baseline_record` and the
read-only consumer lives in :mod:`leopard2_exact_main_baseline_verifier`.  This
module is the producer-side boundary.  In particular, it deliberately does not
import the verifier: a completed lane must be checked by launching that separate
program after the owner-only seal is complete.

The first public layer is intentionally host-independent.  It validates the
eight acquisition roots and publishes an already-constructed authority or
failure record through an fd-anchored, exclusive, crash-conservative lane
writer.  The build/acquisition state machine is layered on these primitives so
that fault-injection tests and the eventual real acquisition use exactly the
same sealing implementation.
"""

from __future__ import annotations

import base64
import copy
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import importlib.util
import math
import os
from pathlib import Path
import re
import resource
import secrets
import selectors
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Mapping, NoReturn, Sequence


def _load_local_contract(module_name: str, filename: str) -> Any:
    """Load one sibling contract by exact path under isolated Python."""
    expected = Path(__file__).resolve().with_name(filename)
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        observed = Path(getattr(loaded, "__file__", "")).resolve()
        if observed != expected:
            raise RuntimeError(
                f"{module_name} was loaded from an unexpected path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact-main contract {expected}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    specification.loader.exec_module(module)
    observed = Path(getattr(module, "__file__", "")).resolve()
    if observed != expected:
        raise RuntimeError(f"{module_name} resolved to an unexpected path")
    return module


identity_contract = _load_local_contract(
    "leopard2_exact_main_baseline", "leopard2_exact_main_baseline.py")
record_contract = _load_local_contract(
    "leopard2_exact_main_baseline_record",
    "leopard2_exact_main_baseline_record.py",
)

ExactMainBaselineError = identity_contract.ExactMainBaselineError
canonical_json_bytes = identity_contract.canonical_json_bytes


TREE_METADATA_SCHEMA = "leopard2-exact-main-baseline-tree-metadata/v1"
BUILD_CLOSURE_SCHEMA = record_contract.BUILD_CLOSURE_SCHEMA
LANE_FILE_MODE = 0o400
LANE_DIRECTORY_MODE = 0o500
WRITABLE_FILE_MODE = 0o600
WRITABLE_DIRECTORY_MODE = 0o700
MAX_PATH_LENGTH = 4096
MAX_TREE_NODES = 1024
MAX_TREE_DEPTH = 32
MAX_SEALED_FILE_BYTES = 2 * 1024 * 1024 * 1024
MAX_SEALED_TOTAL_BYTES = 8 * 1024 * 1024 * 1024
READ_CHUNK = 1024 * 1024
MAX_COMMAND_ARGUMENTS = 512
MAX_COMMAND_ARGUMENT_BYTES = 65536
MAX_COMMAND_OUTPUT_BYTES = 64 * 1024 * 1024
MAX_COMMAND_TIMEOUT_SECONDS = 6 * 60 * 60
COMMAND_TERMINATION_SECONDS = 5.0
MAX_CONTAINED_PROCESSES = 4096
MIN_CONTAINED_PROCESSES = 64
MAX_CPUINFO_BYTES = record_contract.MAX_CPU_COUNT * 4096
CANONICAL_LOCK_PATH = "/tmp/leopard-gf8-authoritative.lock"
SOURCE_ACQUISITION_TIMEOUT_SECONDS = 30 * 60
TOOL_COMMAND_TIMEOUT_SECONDS = 60
MAX_SOURCE_LOG_BYTES = 16 * 1024 * 1024
BUILD_TIMEOUT_SECONDS = 2 * 60 * 60
ATTESTATION_TIMEOUT_SECONDS = 30 * 60
VERIFICATION_TIMEOUT_SECONDS = 30 * 60
BUILD_LOG_SCHEMA = "leopard2-exact-main-build-command-log/v1"
BUILD_STAGE_LOG_SCHEMA = "leopard2-exact-main-build-stage-log/v1"
ASSEMBLY_LOG_SCHEMA = "leopard2-exact-main-authority-assembly-log/v1"
VERIFICATION_LOG_SCHEMA = \
    "leopard2-exact-main-independent-verification-log/v1"
VERIFICATION_FAILURE_LOG_SCHEMA = \
    "leopard2-exact-main-independent-verification-failure-log/v1"
SEAL_LOG_SCHEMA = "leopard2-exact-main-seal-log/v1"
ACQUISITION_REPORT_SCHEMA = "leopard2-exact-main-acquisition-report/v1"
_TREE_METADATA_PATH = "TREE-METADATA.json"
_SHA256SUMS_PATH = "SHA256SUMS"
_TERMINAL_PATHS = frozenset(("baseline-authority.json", "FAILED.json"))
_LDD_ADDRESS = r"\(0x[0-9A-Fa-f]+\)"
_LDD_MAPPED = re.compile(
    rf"^[ \t]*(?P<soname>[^ /\t=]+)[ \t]+=>[ \t]+"
    rf"(?P<path>/[^ \t]+)[ \t]+{_LDD_ADDRESS}[ \t]*$")
_LDD_DIRECT = re.compile(
    rf"^[ \t]*(?P<path>/[^ \t]+)[ \t]+{_LDD_ADDRESS}[ \t]*$")
_LDD_VIRTUAL = re.compile(
    rf"^[ \t]*(?P<soname>linux-(?:vdso(?:32|64)?|gate)\.so\.[0-9]+)[ \t]+"
    rf"{_LDD_ADDRESS}[ \t]*$")


class AcquisitionError(ExactMainBaselineError):
    """The exact-main producer cannot create one trustworthy lane."""


class CommandExecutionError(AcquisitionError):
    """One exact acquisition command returned a producer-visible failure."""

    def __init__(self, message: str, exit_status: int):
        super().__init__(message)
        self.exit_status = exit_status


class SourceStageError(AcquisitionError):
    """A source-stage failure carrying its complete command transcript."""

    def __init__(self, error: BaseException, log: bytes):
        super().__init__(str(error) or type(error).__name__)
        self.log = log
        self.error_kind = "command_error" if isinstance(
            error, CommandExecutionError) else "acquisition_error"
        self.exit_status = error.exit_status if isinstance(
            error, CommandExecutionError) else 1


class BuildStageError(AcquisitionError):
    """A per-role build-stage failure with its complete stage transcript."""

    def __init__(self, role: str, error: BaseException, log: bytes):
        super().__init__(str(error) or type(error).__name__)
        self.role = role
        self.stage = role + "_build"
        self.log = log
        self.error_kind = "command_error" if isinstance(
            error, CommandExecutionError) else "acquisition_error"
        self.exit_status = error.exit_status if isinstance(
            error, CommandExecutionError) else 1


class VerificationError(AcquisitionError):
    """Independent verification failed with retained exact child streams."""

    def __init__(
        self, message: str, *, exit_status: int,
        stdout: bytes, stderr: bytes,
        process_result_observed: bool = True,
    ) -> None:
        super().__init__(message)
        _require(type(process_result_observed) is bool,
                 "verification process-result state is not boolean")
        self.observed_exit_status = exit_status if type(exit_status) is int \
            else None
        self.exit_status = exit_status if (
            type(exit_status) is int and 1 <= exit_status <= 255) else 1
        self.stdout = stdout
        self.stderr = stderr
        self.process_result_observed = process_result_observed


@dataclass(frozen=True)
class AcquisitionOutcome:
    """One sealed terminal outcome from the complete acquisition driver."""

    outcome: str
    lane_root: str
    terminal: str
    seal: dict[str, Any]
    verdict: dict[str, Any] | None
    authority_lane_root: str | None
    authority_record_sha256: str | None
    verification_program: dict[str, Any] | None
    verification_exit_status: int | None
    message: str | None


def _fail(message: str) -> NoReturn:
    raise AcquisitionError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _portable_absolute_path(value: Any, label: str) -> str:
    _require(type(value) is str and 1 < len(value) <= MAX_PATH_LENGTH,
             f"{label} is not a bounded absolute path")
    _require(all(0x21 <= ord(character) <= 0x7E for character in value),
             f"{label} is not printable ASCII")
    _require(value.startswith("/") and value != "/" and "//" not in value,
             f"{label} is not a canonical absolute POSIX path")
    _require(all(part not in ("", ".", "..")
                 for part in value.split("/")[1:]),
             f"{label} is not a canonical absolute POSIX path")
    return value


def _safe_relative_path(value: Any, label: str) -> str:
    _require(type(value) is str and 0 < len(value) <= MAX_PATH_LENGTH,
             f"{label} is not a bounded relative path")
    _require(all(0x21 <= ord(character) <= 0x7E for character in value),
             f"{label} is not printable ASCII")
    _require(not value.startswith("/") and not value.endswith("/") and
             "//" not in value and
             all(part not in ("", ".", "..")
                 for part in value.split("/")),
             f"{label} is not canonical")
    return value


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _write_all(descriptor: int, content: bytes, label: str) -> None:
    view = memoryview(content)
    cursor = 0
    while cursor < len(view):
        written = os.write(descriptor, view[cursor:])
        _require(written > 0, f"{label} stopped while being written")
        cursor += written


def _hash_fd(descriptor: int, size: int, label: str) -> str:
    _require(0 <= size <= MAX_SEALED_FILE_BYTES,
             f"{label} exceeds the retained byte bound")
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    remaining = size
    while remaining:
        content = os.read(descriptor, min(READ_CHUNK, remaining))
        _require(bool(content), f"{label} ended while being hashed")
        digest.update(content)
        remaining -= len(content)
    _require(os.read(descriptor, 1) == b"", f"{label} grew while hashed")
    return digest.hexdigest()


def _read_bounded_system_file(path: str, maximum: int, label: str) -> bytes:
    _require(type(maximum) is int and 0 < maximum <= MAX_SEALED_FILE_BYTES,
             f"{label} byte bound is invalid")
    descriptor = -1
    try:
        descriptor = os.open(
            path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(descriptor, min(READ_CHUNK, maximum + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            _require(total <= maximum, f"{label} is oversized")
        content = b"".join(chunks)
        _require(bool(content), f"{label} is empty")
        return content
    except OSError as error:
        raise AcquisitionError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _parse_linux_cpu_list(content: bytes) -> list[int]:
    _require(type(content) is bytes and 0 < len(content) <= 65536,
             "Linux online CPU list is not bounded bytes")
    try:
        text = content.decode("ascii")
    except UnicodeDecodeError as error:
        raise AcquisitionError(
            "Linux online CPU list is not strict ASCII") from error
    _require(text.endswith("\n") and text.count("\n") == 1,
             "Linux online CPU list is not one LF-framed line")
    body = text[:-1]
    _require(bool(body) and all(character in "0123456789,-" for character in body),
             "Linux online CPU list has invalid characters")
    cpus: list[int] = []
    for index, field in enumerate(body.split(",")):
        _require(bool(field),
                 f"Linux online CPU field {index} is empty")
        if "-" in field:
            parts = field.split("-")
            _require(len(parts) == 2 and all(parts),
                     f"Linux online CPU range {index} is invalid")
            try:
                first, last = (int(part, 10) for part in parts)
            except ValueError as error:
                raise AcquisitionError(
                    f"Linux online CPU range {index} is not numeric") from error
            _require(0 <= first <= last <= (1 << 31) - 1,
                     f"Linux online CPU range {index} is invalid")
            count = last - first + 1
            _require(len(cpus) + count <= record_contract.MAX_CPU_COUNT,
                     "Linux online CPU list exceeds the contract bound")
            cpus.extend(range(first, last + 1))
        else:
            try:
                cpu = int(field, 10)
            except ValueError as error:
                raise AcquisitionError(
                    f"Linux online CPU field {index} is not numeric") from error
            _require(0 <= cpu <= (1 << 31) - 1,
                     f"Linux online CPU field {index} is invalid")
            _require(len(cpus) < record_contract.MAX_CPU_COUNT,
                     "Linux online CPU list exceeds the contract bound")
            cpus.append(cpu)
    _require(cpus == sorted(set(cpus)),
             "Linux online CPU list is not sorted and unique")
    return cpus


def _parse_linux_cpu_models(content: bytes) -> dict[int, str]:
    _require(type(content) is bytes and 0 < len(content) <= MAX_CPUINFO_BYTES,
             "host CPU identity is not bounded bytes")
    try:
        text = content.decode("ascii")
    except UnicodeDecodeError as error:
        raise AcquisitionError(
            "host CPU identity is not strict ASCII") from error
    records: dict[int, str] = {}
    for block_index, block in enumerate(text.strip().split("\n\n")):
        fields: dict[str, str] = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            name, value = line.split(":", 1)
            key = name.strip()
            retained = value.strip()
            if key in ("processor", "model name"):
                _require(key not in fields,
                         f"host CPU block {block_index} repeats {key!r}")
                fields[key] = retained
        _require(set(fields) == {"processor", "model name"},
                 f"host CPU block {block_index} lacks its x86 identity")
        try:
            processor = int(fields["processor"], 10)
        except ValueError as error:
            raise AcquisitionError(
                f"host CPU block {block_index} has a nonnumeric processor") \
                from error
        model = fields["model name"]
        _require(0 <= processor <= (1 << 31) - 1 and
                 processor not in records and model and
                 len(model) <= record_contract.MAX_TEXT_LENGTH and
                 all(character not in "\0\r\n" for character in model),
                 f"host CPU block {block_index} has an invalid identity")
        records[processor] = model
        _require(len(records) <= record_contract.MAX_CPU_COUNT,
                 "host CPU identity exceeds the contract bound")
    _require(bool(records), "host CPU identity has no processor records")
    return records


@dataclass(frozen=True)
class LanePlan:
    """The nine immutable roots and attempt identity for one acquisition."""

    lane_root: str
    failure_lane_root: str
    attempt: int
    repository: str
    verifier: str
    canonical_adapter_root: str
    canonical_baseline_root: str
    canonical_build_root: str
    variant_adapter_root: str
    variant_baseline_root: str
    variant_build_root: str
    scratch_root: str


_PLAN_PATH_FIELDS = (
    "lane_root",
    "failure_lane_root",
    "repository",
    "verifier",
    "canonical_adapter_root",
    "canonical_baseline_root",
    "canonical_build_root",
    "variant_adapter_root",
    "variant_baseline_root",
    "variant_build_root",
    "scratch_root",
)


def validate_lane_plan(value: Any) -> LanePlan:
    """Validate and detach one acquisition plan.

    The lane, six build/source roots, and scratch root are mutually
    non-containing both component-wise and as raw UTF-8 byte strings.  The
    latter is deliberately stricter than the record contract because the ELF
    census searches for exact path bytes and must never attribute one root's
    occurrence to another root.
    """
    _require(isinstance(value, LanePlan), "lane plan has the wrong type")
    paths = {
        field: _portable_absolute_path(getattr(value, field), f"lane {field}")
        for field in _PLAN_PATH_FIELDS
    }
    _require(type(value.attempt) is int and 1 <= value.attempt <= 3,
             "lane attempt is outside the frozen three-attempt budget")
    root_fields = (
        "lane_root", "failure_lane_root",
        "canonical_adapter_root", "canonical_baseline_root",
        "canonical_build_root", "variant_adapter_root",
        "variant_baseline_root", "variant_build_root", "scratch_root",
    )
    roots = [paths[field] for field in root_fields]
    _require(len(roots) == len(set(roots)), "lane roots contain duplicates")
    for index, left in enumerate(roots):
        for right in roots[index + 1:]:
            _require(not (left + "/").startswith(right + "/") and
                     not (right + "/").startswith(left + "/"),
                     "lane roots overlap by path containment")
            left_bytes = left.encode("ascii")
            right_bytes = right.encode("ascii")
            _require(left_bytes not in right_bytes and
                     right_bytes not in left_bytes,
                     "lane roots overlap as census byte strings")
    repository = paths["repository"]
    verifier = paths["verifier"]
    _require(verifier.startswith(repository + "/"),
             "lane verifier is not beneath its repository")
    for root in roots:
        _require(not (root + "/").startswith(repository + "/") and
                 not (repository + "/").startswith(root + "/"),
                 "lane root overlaps its controller repository")
    return LanePlan(
        lane_root=paths["lane_root"],
        failure_lane_root=paths["failure_lane_root"],
        attempt=value.attempt,
        repository=paths["repository"],
        verifier=paths["verifier"],
        canonical_adapter_root=paths["canonical_adapter_root"],
        canonical_baseline_root=paths["canonical_baseline_root"],
        canonical_build_root=paths["canonical_build_root"],
        variant_adapter_root=paths["variant_adapter_root"],
        variant_baseline_root=paths["variant_baseline_root"],
        variant_build_root=paths["variant_build_root"],
        scratch_root=paths["scratch_root"],
    )


_GIT_CAPTURE_MODULE_NAME = "_leopard2_exact_main_git_capture"
_GIT_CAPTURE_RELATIVE_PATH = \
    "experiments/leopard2/main_compare/git_capture.py"
_PREPARED_ROOT_FIELDS = (
    "canonical_adapter_root", "canonical_baseline_root",
    "canonical_build_root", "variant_adapter_root",
    "variant_baseline_root", "variant_build_root", "scratch_root",
)


def _path_is_absent(path: str) -> bool:
    try:
        os.lstat(path)
    except FileNotFoundError:
        return True
    except OSError as error:
        raise AcquisitionError(
            f"cannot inspect acquisition path {path!r}: {error}") from error
    return False


def _load_repository_module(
    repository: str, module_name: str, relative_path: str,
) -> Any:
    """Load one repository controller by an exact, anchored source path."""
    root = Path(_portable_absolute_path(
        repository, "repository module root"))
    _require(module_name == _GIT_CAPTURE_MODULE_NAME and
             relative_path == _GIT_CAPTURE_RELATIVE_PATH,
             "repository module request is not allowlisted")
    try:
        resolved_root = root.resolve(strict=True)
        expected = (root / relative_path).resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise AcquisitionError(
            f"cannot resolve repository module: {error}") from error
    _require(resolved_root == root and
             expected == root / relative_path and expected.is_file(),
             "repository module path is not canonical")
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        observed = Path(getattr(loaded, "__file__", "")).resolve()
        _require(observed == expected,
                 "repository module was loaded from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(
        module_name, expected)
    _require(specification is not None and
             specification.loader is not None,
             "repository module could not be loaded")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    observed = Path(getattr(module, "__file__", "")).resolve()
    _require(observed == expected,
             "repository module resolved to another path")
    return module


def frozen_child_environment() -> dict[str, str]:
    """Return exactly the child environment bound by the frozen profile."""
    profile = record_contract.exact_main_build_profile()
    values = profile["environment"]
    _require(type(values) is list and values and all(
        type(item) is dict and set(item) == {"name", "value"} and
        type(item["name"]) is str and type(item["value"]) is str
        for item in values), "frozen child environment is invalid")
    result = {item["name"]: item["value"] for item in values}
    _require(len(result) == len(values),
             "frozen child environment repeats a name")
    return _command_environment(result)


class PreparedAcquisitionRoots:
    """Own and revalidate all writable roots other than the sealed lane."""

    def __init__(self, plan: LanePlan):
        self.plan = validate_lane_plan(plan)
        self.descriptors: dict[str, int] = {}
        self.identities: dict[str, tuple[int, int, int, int, int]] = {}
        self.created: list[str] = []
        self.entered = False

    @staticmethod
    def _identity(status: os.stat_result) -> tuple[int, int, int, int, int]:
        _require(stat.S_ISDIR(status.st_mode) and
                 status.st_uid == os.geteuid() and
                 status.st_gid == os.getegid() and
                 stat.S_IMODE(status.st_mode) == WRITABLE_DIRECTORY_MODE,
                 "acquisition root has unsafe ownership or mode")
        return (status.st_dev, status.st_ino, status.st_uid, status.st_gid,
                stat.S_IMODE(status.st_mode))

    def __enter__(self) -> "PreparedAcquisitionRoots":
        _require(not self.entered and not self.descriptors,
                 "acquisition roots are already prepared")
        _require(getattr(shutil.rmtree, "avoids_symlink_attacks", False),
                 "safe acquisition-root cleanup is unavailable")
        try:
            repository = Path(self.plan.repository)
            verifier = Path(self.plan.verifier)
            _require(repository.resolve(strict=True) == repository and
                     repository.is_dir(),
                     "controller repository is not a canonical directory")
            _require(verifier.resolve(strict=True) == verifier and
                     verifier.is_file(),
                     "independent verifier is not a canonical file")
            _require(_path_is_absent(self.plan.lane_root) and
                     _path_is_absent(self.plan.failure_lane_root),
                     "lane or failure root already exists and cannot be reused")
            for field in _PREPARED_ROOT_FIELDS:
                path = getattr(self.plan, field)
                _require(_path_is_absent(path),
                         f"acquisition root {field} already exists")
                try:
                    os.mkdir(path, WRITABLE_DIRECTORY_MODE)
                except OSError as error:
                    raise AcquisitionError(
                        f"cannot create acquisition root {field}: {error}") \
                        from error
                self.created.append(field)
                descriptor = -1
                try:
                    descriptor = os.open(
                        path,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                        os.O_CLOEXEC)
                    status = os.fstat(descriptor)
                    pathname = os.lstat(path)
                    identity = self._identity(status)
                    _require(os.path.realpath(path) == path and
                             (status.st_dev, status.st_ino) ==
                             (pathname.st_dev, pathname.st_ino),
                             f"acquisition root {field} is not canonical")
                    self.descriptors[field] = descriptor
                    self.identities[field] = identity
                    descriptor = -1
                except OSError as error:
                    raise AcquisitionError(
                        f"cannot retain acquisition root {field}: {error}") \
                        from error
                finally:
                    if descriptor >= 0:
                        try:
                            os.close(descriptor)
                        except OSError as error:
                            raise AcquisitionError(
                                f"cannot release acquisition root {field}: "
                                f"{error}") from error
            self.entered = True
            self.validate_current()
            return self
        except BaseException:
            self._cleanup()
            raise

    def validate_current(self) -> None:
        _require(self.entered and
                 set(self.descriptors) == set(_PREPARED_ROOT_FIELDS),
                 "acquisition roots are not prepared")
        for field in _PREPARED_ROOT_FIELDS:
            path = getattr(self.plan, field)
            descriptor = os.fstat(self.descriptors[field])
            try:
                pathname = os.lstat(path)
            except OSError as error:
                raise AcquisitionError(
                    f"acquisition root {field} disappeared: {error}") \
                    from error
            _require(self._identity(descriptor) == self.identities[field] and
                     (descriptor.st_dev, descriptor.st_ino) ==
                     (pathname.st_dev, pathname.st_ino) and
                     os.path.realpath(path) == path,
                     f"acquisition root {field} was replaced")

    def reset_root(self, field: str) -> None:
        """Empty one build root in place without changing its identity."""
        _require(field in ("canonical_build_root", "variant_build_root") and
                 type(field) is str,
                 "only an exact-main build root may be reset")
        self.validate_current()
        descriptor = self.descriptors[field]
        try:
            names = sorted(os.listdir(descriptor))
            _require(len(names) <= MAX_TREE_NODES,
                     "build root exceeds the reset node bound")
            for name in names:
                _require(type(name) is str and name not in ("", ".", "..") and
                         "/" not in name and "\0" not in name,
                         "build root contains an unsafe entry name")
                status = os.stat(
                    name, dir_fd=descriptor, follow_symlinks=False)
                if stat.S_ISDIR(status.st_mode):
                    # ``shutil.rmtree(dir_fd=...)`` is only available from
                    # Python 3.12, while the repository supports Python 3.10.
                    # The prepared-root identity and the no-follow stat above
                    # establish the same canonical child for this producer
                    # boundary before using the portable pathname form.
                    shutil.rmtree(getattr(self.plan, field) + "/" + name)
                else:
                    os.unlink(name, dir_fd=descriptor)
            os.fsync(descriptor)
            _require(os.listdir(descriptor) == [],
                     "build root is not empty after reset")
            self.validate_current()
        except OSError as error:
            raise AcquisitionError(
                f"cannot reset acquisition root {field}: {error}") from error

    def _cleanup(self) -> None:
        errors: list[str] = []
        for field in reversed(self.created):
            path = getattr(self.plan, field)
            descriptor = self.descriptors.pop(field, -1)
            identity = self.identities.pop(field, None)
            try:
                if descriptor >= 0:
                    current = os.fstat(descriptor)
                    pathname = os.lstat(path)
                    _require(identity is not None and
                             self._identity(current) == identity and
                             (current.st_dev, current.st_ino) ==
                             (pathname.st_dev, pathname.st_ino),
                             f"acquisition root {field} changed before cleanup")
                shutil.rmtree(path)
                _require(_path_is_absent(path),
                         f"acquisition root {field} survived cleanup")
            except BaseException as error:
                errors.append(f"{field}: {type(error).__name__}: {error}")
            finally:
                if descriptor >= 0:
                    try:
                        os.close(descriptor)
                    except OSError as error:
                        errors.append(f"{field}: close: {error}")
        self.created.clear()
        self.entered = False
        if errors:
            raise AcquisitionError(
                "acquisition-root cleanup failed: " + "; ".join(errors))

    def __exit__(self, _kind: object, value: object,
                 _traceback: object) -> None:
        try:
            self._cleanup()
        except BaseException as cleanup_error:
            if value is None:
                raise
            raise AcquisitionError(
                f"acquisition-root cleanup failed after "
                f"{type(value).__name__}: {value}; {cleanup_error}") \
                from cleanup_error


def prepare_acquisition_roots(plan: LanePlan) -> PreparedAcquisitionRoots:
    """Return the exclusive context that owns one plan's writable roots."""
    return PreparedAcquisitionRoots(plan)


def canonical_ldd_text(rows: Sequence[Mapping[str, Any]]) -> bytes:
    """Serialize the inverse of the frozen canonical-ldd-C-v1 parser."""
    _require(type(rows) in (list, tuple) and
             0 < len(rows) <= record_contract.MAX_DEPENDENCIES,
             "canonical ldd rows are not a bounded non-empty sequence")
    lines: list[str] = []
    canonical_rows: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        _require(type(row) is dict and set(row) == {"soname", "kind", "path"},
                 f"canonical ldd row {index} has an unexpected key set")
        soname = row["soname"]
        _require(type(soname) is str and
                 0 < len(soname) <= record_contract.MAX_TEXT_LENGTH and
                 all(0x21 <= ord(character) <= 0x7E for character in soname) and
                 "/" not in soname and "\t" not in soname,
                 f"canonical ldd row {index} has an invalid soname")
        if row["kind"] == "virtual":
            _require(row["path"] is None,
                     f"canonical ldd row {index} gives a virtual path")
            lines.append(f"{soname}\tvirtual\n")
            canonical_rows.append({
                "soname": soname, "kind": "virtual", "path": None})
        else:
            _require(row["kind"] == "file",
                     f"canonical ldd row {index} has an invalid kind")
            path = _portable_absolute_path(
                row["path"], f"canonical ldd row {index} path")
            lines.append(f"{soname}\tfile\t{path}\n")
            canonical_rows.append({
                "soname": soname, "kind": "file", "path": path})
    sonames = [row["soname"] for row in rows]
    _require(sonames == sorted(set(sonames)),
             "canonical ldd rows are not sorted and unique")
    content = "".join(lines).encode("ascii")
    try:
        parsed = record_contract.parse_canonical_ldd_output(content)
    except ExactMainBaselineError as error:
        raise AcquisitionError("canonical ldd rows are invalid") from error
    _require(identity_contract.exact_json_equal(list(parsed), canonical_rows),
             "canonical ldd serialization did not round trip")
    return content


def normalize_ldd_output(raw: bytes) -> tuple[dict[str, Any], ...]:
    """Normalize one successful ``LC_ALL=C ldd`` stdout transcript.

    Load addresses and presentation whitespace are discarded.  File-backed
    dependencies retain the exact absolute path printed by the frozen ldd
    executable; the separately retained dependency identities bind those
    paths to bytes.  Unexpected diagnostics, missing libraries, static-link
    messages, duplicate SONAMEs, and duplicate file paths fail closed.
    """
    _require(type(raw) is bytes and
             0 < len(raw) <= record_contract.MAX_CANONICAL_LDD_BYTES,
             "raw ldd output is not bounded non-empty bytes")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as error:
        raise AcquisitionError("raw ldd output is not strict ASCII") from error
    _require(text.endswith("\n") and "\r" not in text and "\0" not in text,
             "raw ldd output is not LF-framed")
    lines = text[:-1].split("\n")
    _require(0 < len(lines) <= record_contract.MAX_DEPENDENCIES and
             all(line.strip(" \t") for line in lines),
             "raw ldd output row count is invalid")
    rows: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        match = _LDD_MAPPED.fullmatch(line)
        if match is not None:
            soname = match.group("soname")
            path = _portable_absolute_path(
                match.group("path"), f"raw ldd row {index} path")
            rows.append({"soname": soname, "kind": "file", "path": path})
            continue
        match = _LDD_DIRECT.fullmatch(line)
        if match is not None:
            path = _portable_absolute_path(
                match.group("path"), f"raw ldd row {index} path")
            rows.append({
                "soname": path.rsplit("/", 1)[-1],
                "kind": "file",
                "path": path,
            })
            continue
        match = _LDD_VIRTUAL.fullmatch(line)
        if match is not None:
            rows.append({
                "soname": match.group("soname"),
                "kind": "virtual",
                "path": None,
            })
            continue
        _fail(f"raw ldd row {index} is not a supported dependency record")
    rows.sort(key=lambda item: item["soname"])
    sonames = [item["soname"] for item in rows]
    file_paths = [item["path"] for item in rows if item["kind"] == "file"]
    _require(sonames == sorted(set(sonames)),
             "raw ldd output contains duplicate SONAMEs")
    _require(len(file_paths) == len(set(file_paths)),
             "raw ldd output contains duplicate dependency paths")
    canonical = canonical_ldd_text(rows)
    try:
        parsed = record_contract.parse_canonical_ldd_output(canonical)
    except ExactMainBaselineError as error:
        raise AcquisitionError("normalized ldd output is invalid") from error
    return tuple(copy.deepcopy(parsed))


@dataclass(frozen=True)
class CommandResult:
    """One bounded child-process result with byte-exact output streams."""

    exit_status: int
    stdout: bytes
    stderr: bytes


@dataclass(frozen=True)
class StreamedCommandResult:
    """One child result whose stdout was streamed to an owned file."""

    exit_status: int
    stdout_size: int
    stdout_sha256: str
    stderr: bytes


@dataclass(frozen=True)
class SourceStageResult:
    """Detached source/tool evidence plus its retained byte/path projection."""

    source: dict[str, Any]
    adapter: dict[str, Any]
    toolchain: dict[str, Any]
    retained_bytes: dict[str, bytes]
    retained_paths: dict[str, str]
    log: bytes


@dataclass(frozen=True)
class BuildStageResult:
    """One complete role's build, runtime, attestation, and ELF evidence."""

    role: str
    build: dict[str, Any]
    runtime: dict[str, Any]
    attestation: dict[str, Any]
    identity: dict[str, Any]
    retained_bytes: dict[str, bytes]
    retained_paths: dict[str, str]
    log: bytes


def _command_arguments(value: Any) -> list[str]:
    _require(type(value) in (list, tuple) and
             1 <= len(value) <= MAX_COMMAND_ARGUMENTS,
             "command argument vector is invalid")
    arguments: list[str] = []
    for index, argument in enumerate(value):
        _require(type(argument) is str and argument and "\0" not in argument and
                 len(argument.encode("utf-8")) <= MAX_COMMAND_ARGUMENT_BYTES,
                 f"command argument {index} is invalid")
        arguments.append(argument)
    return arguments


def _command_environment(value: Any) -> dict[str, str]:
    _require(type(value) is dict and len(value) <= 128,
             "command environment is invalid")
    _require(all(type(name) is str for name in value),
             "command environment contains a non-string name")
    result: dict[str, str] = {}
    for name in sorted(value):
        retained = value[name]
        _require(type(name) is str and name and "=" not in name and
                 "\0" not in name and len(name.encode("utf-8")) <= 1024,
                 "command environment name is invalid")
        _require(type(retained) is str and "\0" not in retained and
                 len(retained.encode("utf-8")) <= MAX_COMMAND_ARGUMENT_BYTES,
                 f"command environment value {name!r} is invalid")
        result[name] = retained
    return result


def _linux_subreaper_state() -> int:
    state = ctypes.c_int(-1)
    function = ctypes.CDLL(None, use_errno=True).prctl
    ctypes.set_errno(0)
    if function(37, ctypes.byref(state), 0, 0, 0) != 0:
        number = ctypes.get_errno()
        raise AcquisitionError(
            "cannot read Linux child-subreaper state: " +
            os.strerror(number or errno.EPERM))
    _require(state.value in (0, 1),
             "Linux child-subreaper state is invalid")
    return state.value


def _set_linux_subreaper(state: int) -> None:
    _require(type(state) is int and state in (0, 1),
             "Linux child-subreaper state is invalid")
    function = ctypes.CDLL(None, use_errno=True).prctl
    ctypes.set_errno(0)
    if function(36, ctypes.c_ulong(state), 0, 0, 0) != 0:
        number = ctypes.get_errno()
        raise AcquisitionError(
            "cannot set Linux child-subreaper state: " +
            os.strerror(number or errno.EPERM))
    _require(_linux_subreaper_state() == state,
             "Linux child-subreaper state did not change")


def _linux_process_identity(pid: int) -> tuple[int, int] | None:
    """Return ``(parent_pid, start_time_ticks)`` from Linux procfs."""
    try:
        raw = Path(f"/proc/{pid}/stat").read_bytes()
    except (FileNotFoundError, ProcessLookupError):
        return None
    except PermissionError as error:
        try:
            owner = os.stat(f"/proc/{pid}", follow_symlinks=False).st_uid
        except (FileNotFoundError, ProcessLookupError):
            return None
        except OSError as owner_error:
            raise AcquisitionError(
                f"cannot identify hidden Linux process {pid}: "
                f"{owner_error}") from owner_error
        if owner != os.getuid():
            return None
        raise AcquisitionError(
            f"cannot inspect owned Linux process {pid}: {error}") from error
    except OSError as error:
        if error.errno in (errno.EACCES, errno.EPERM):
            try:
                owner = os.stat(
                    f"/proc/{pid}", follow_symlinks=False).st_uid
            except (FileNotFoundError, ProcessLookupError):
                return None
            except OSError as owner_error:
                raise AcquisitionError(
                    f"cannot identify hidden Linux process {pid}: "
                    f"{owner_error}") from owner_error
            if owner != os.getuid():
                return None
            raise AcquisitionError(
                f"cannot inspect owned Linux process {pid}: {error}") from error
        raise AcquisitionError(
            f"cannot inspect Linux process {pid}: {error}") from error
    _require(0 < len(raw) <= 1024 * 1024 and b"\0" not in raw,
             f"Linux process {pid} stat record is invalid")
    closing = raw.rfind(b")")
    _require(closing > 1 and closing + 2 < len(raw) and
             raw[closing + 1:closing + 2] == b" ",
             f"Linux process {pid} stat record is malformed")
    fields = raw[closing + 2:].split()
    # fields[0] is proc field 3 (state), fields[1] is ppid, and fields[19]
    # is field 22 (starttime).  Starttime closes the PID-reuse race.
    _require(len(fields) >= 20,
             f"Linux process {pid} stat record is truncated")
    try:
        parent = int(fields[1], 10)
        started = int(fields[19], 10)
    except ValueError as error:
        raise AcquisitionError(
            f"Linux process {pid} stat record is not numeric") from error
    _require(parent >= 0 and started >= 0,
             f"Linux process {pid} stat identity is invalid")
    return parent, started


def _linux_process_snapshot() -> dict[int, tuple[int, int]]:
    try:
        names = os.listdir("/proc")
    except OSError as error:
        raise AcquisitionError(
            f"cannot enumerate Linux procfs: {error}") from error
    result: dict[int, tuple[int, int]] = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name, 10)
        identity = _linux_process_identity(pid)
        if identity is not None:
            result[pid] = identity
    self_identity = _linux_process_identity(os.getpid())
    observed_self = result.get(os.getpid())
    _require(self_identity is not None and observed_self is not None and
             observed_self[1] == self_identity[1],
             "Linux procfs does not expose the acquisition process")
    return result


def _linux_libc_function(name: str) -> Any | None:
    """Return one Linux libc entry point when CPython omits its wrapper."""
    try:
        return getattr(ctypes.CDLL(None, use_errno=True), name, None)
    except OSError:
        return None


def _linux_pidfd_open(pid: int, flags: int = 0) -> int:
    """Open a pidfd through CPython or the equivalent Linux libc entry."""
    native = getattr(os, "pidfd_open", None)
    if callable(native):
        return native(pid, flags)
    _require(sys.platform.startswith("linux"),
             "bounded command containment requires Linux pidfd_open")
    function = _linux_libc_function("pidfd_open")
    _require(function is not None,
             "bounded command containment requires Linux pidfd_open")
    function.argtypes = (ctypes.c_int, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(flags))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno() or errno.EPERM
    if number == errno.ESRCH:
        raise ProcessLookupError(number, os.strerror(number), str(pid))
    raise OSError(number, os.strerror(number), str(pid))


def _linux_pidfd_send_signal(
    descriptor: int, signal_number: int, flags: int = 0,
) -> None:
    """Signal through a pidfd through CPython or Linux libc."""
    native = getattr(signal, "pidfd_send_signal", None)
    if callable(native):
        native(descriptor, signal_number, None, flags)
        return
    _require(sys.platform.startswith("linux"),
             "bounded command containment requires Linux pidfd_send_signal")
    function = _linux_libc_function("pidfd_send_signal")
    _require(function is not None,
             "bounded command containment requires Linux pidfd_send_signal")
    function.argtypes = (
        ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    status = function(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number), None,
        ctypes.c_uint(flags))
    if status == 0:
        return
    number = ctypes.get_errno() or errno.EPERM
    if number == errno.ESRCH:
        raise ProcessLookupError(
            number, os.strerror(number), str(descriptor))
    raise OSError(number, os.strerror(number), str(descriptor))


def _validate_linux_pidfd_runtime() -> None:
    """Prove both pidfd operations work before launching any command."""
    descriptor = -1
    try:
        descriptor = _linux_pidfd_open(os.getpid())
        _linux_pidfd_send_signal(descriptor, 0)
    except (AcquisitionError, OSError) as error:
        raise AcquisitionError(
            f"bounded command containment requires working Linux pidfds: "
            f"{error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


class _LinuxChildContainment:
    """Contain one command tree, including setsid and double-fork escapes."""

    def __init__(self) -> None:
        self.owner_pid = os.getpid()
        self.previous_subreaper: int | None = None
        self.baseline_children: set[tuple[int, int]] = set()
        self.process: subprocess.Popen[bytes] | None = None
        self.leader_identity: tuple[int, int] | None = None
        self.handles: dict[tuple[int, int], int] = {}
        self.reserve_descriptor = -1
        self.containment_bound_reached = False
        self.pidfd_resource_exhausted = False
        self.process_bound = MAX_CONTAINED_PROCESSES
        self.active = False

    @staticmethod
    def _children(
        snapshot: Mapping[int, tuple[int, int]], parent: int,
    ) -> set[tuple[int, int]]:
        return {
            (pid, identity[1]) for pid, identity in snapshot.items()
            if identity[0] == parent
        }

    def __enter__(self) -> "_LinuxChildContainment":
        _require(sys.platform.startswith("linux") and Path("/proc/self/task").is_dir(),
                 "bounded command containment requires Linux procfs")
        try:
            threads = [name for name in os.listdir("/proc/self/task")
                       if name.isascii() and name.isdigit()]
        except OSError as error:
            raise AcquisitionError(
                f"cannot inspect acquisition process threads: {error}") from error
        _require(len(threads) == 1,
                 "bounded command containment requires one acquisition thread")
        _validate_linux_pidfd_runtime()
        try:
            open_descriptors = len(os.listdir("/proc/self/fd"))
            descriptor_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
        except (OSError, ValueError) as error:
            raise AcquisitionError(
                f"cannot inspect command descriptor capacity: {error}") from error
        available = MAX_CONTAINED_PROCESSES \
            if descriptor_limit == resource.RLIM_INFINITY else \
            max(0, descriptor_limit - open_descriptors - 16)
        self.process_bound = min(MAX_CONTAINED_PROCESSES, available)
        _require(self.process_bound >= MIN_CONTAINED_PROCESSES,
                 "bounded command containment lacks pidfd capacity; "
                 "raise the process soft file-descriptor limit")
        self.previous_subreaper = _linux_subreaper_state()
        try:
            self.reserve_descriptor = os.open(
                "/dev/null", os.O_RDONLY | os.O_CLOEXEC)
            _set_linux_subreaper(1)
            snapshot = _linux_process_snapshot()
            self.baseline_children = self._children(snapshot, self.owner_pid)
            _require(not self.baseline_children,
                     "bounded command containment found a pre-existing child")
            self.active = True
            return self
        except BaseException:
            if self.reserve_descriptor >= 0:
                os.close(self.reserve_descriptor)
                self.reserve_descriptor = -1
            _set_linux_subreaper(self.previous_subreaper)
            self.previous_subreaper = None
            raise

    def _retain(self, pid: int, started: int) -> None:
        identity = (pid, started)
        if identity in self.handles:
            return
        if len(self.handles) >= self.process_bound:
            self.containment_bound_reached = True
            return
        try:
            descriptor = _linux_pidfd_open(pid)
        except ProcessLookupError:
            return
        except OSError as error:
            if (error.errno in (errno.EMFILE, errno.ENFILE) and
                    self.reserve_descriptor >= 0):
                os.close(self.reserve_descriptor)
                self.reserve_descriptor = -1
                self.pidfd_resource_exhausted = True
                try:
                    descriptor = _linux_pidfd_open(pid)
                except ProcessLookupError:
                    return
                except OSError as retry_error:
                    if retry_error.errno in (errno.EMFILE, errno.ENFILE):
                        return
                    raise AcquisitionError(
                        f"cannot retain Linux process {pid}: "
                        f"{retry_error}") from retry_error
            elif error.errno in (errno.EMFILE, errno.ENFILE):
                self.pidfd_resource_exhausted = True
                return
            else:
                raise AcquisitionError(
                    f"cannot retain Linux process {pid}: {error}") from error
        try:
            observed = _linux_process_identity(pid)
        except BaseException:
            os.close(descriptor)
            raise
        if observed is None or observed[1] != started:
            os.close(descriptor)
            return
        self.handles[identity] = descriptor

    def _prune_handles(self) -> None:
        for identity, descriptor in list(self.handles.items()):
            observed = _linux_process_identity(identity[0])
            if observed is not None and observed[1] == identity[1]:
                continue
            os.close(descriptor)
            del self.handles[identity]

    def attach(self, process: subprocess.Popen[bytes]) -> None:
        _require(self.active and self.process is None and process.pid > 0,
                 "bounded command containment attachment is invalid")
        self.process = process
        observed = _linux_process_identity(process.pid)
        _require(observed is not None and observed[0] == self.owner_pid,
                 "spawned command is not an owned direct child")
        self.leader_identity = (process.pid, observed[1])
        self._retain(process.pid, observed[1])
        _require(self.leader_identity in self.handles,
                 "spawned command could not be retained by pidfd")

    def _discover(self) -> set[tuple[int, int]]:
        snapshot = _linux_process_snapshot()
        targets: set[tuple[int, int]] = set()
        if self.leader_identity is not None:
            leader = snapshot.get(self.leader_identity[0])
            if leader is not None and leader[1] == self.leader_identity[1]:
                targets.add(self.leader_identity)
        targets.update(
            identity for identity in self._children(snapshot, self.owner_pid)
            if identity not in self.baseline_children)
        changed = True
        while changed:
            changed = False
            parent_pids = {identity[0] for identity in targets}
            for pid, record in snapshot.items():
                identity = (pid, record[1])
                if record[0] in parent_pids and identity not in targets:
                    targets.add(identity)
                    changed = True
        for pid, started in sorted(targets):
            self._retain(pid, started)
        return targets

    def _signal_all(self) -> None:
        for identity, descriptor in sorted(
                self.handles.items(), reverse=True):
            observed = _linux_process_identity(identity[0])
            if observed is None or observed[1] != identity[1]:
                continue
            try:
                _linux_pidfd_send_signal(descriptor, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except OSError as error:
                if error.errno != errno.ESRCH:
                    raise AcquisitionError(
                        f"cannot terminate Linux process {identity[0]}: "
                        f"{error}") from error

    def _reap_adopted(self) -> None:
        leader_pid = self.process.pid if self.process is not None else -1
        for identity in sorted(self.handles):
            pid = identity[0]
            if identity == self.leader_identity or pid == leader_pid:
                continue
            try:
                waited, _status = os.waitpid(pid, os.WNOHANG)
            except ChildProcessError:
                continue
            _require(waited in (0, pid),
                     "an unexpected adopted command child was reaped")

    def _cleanup(self) -> None:
        process = self.process
        if process is None:
            snapshot = _linux_process_snapshot()
            _require(self._children(snapshot, self.owner_pid) ==
                     self.baseline_children,
                     "an unattached child appeared during command containment")
            return
        deadline = time.monotonic() + COMMAND_TERMINATION_SECONDS
        empty_scans = 0
        while True:
            self._prune_handles()
            targets = self._discover()
            self._signal_all()
            if process.poll() is None:
                remaining = deadline - time.monotonic()
                _require(remaining > 0,
                         "command leader could not be reaped after SIGKILL")
                try:
                    process.wait(timeout=min(0.05, remaining))
                except subprocess.TimeoutExpired:
                    pass
            self._reap_adopted()
            self._prune_handles()
            snapshot = _linux_process_snapshot()
            current = {
                identity for identity in self._children(snapshot, self.owner_pid)
                if identity not in self.baseline_children
            }
            if not targets and not current and process.poll() is not None:
                empty_scans += 1
                if empty_scans >= 2:
                    _require(not self.containment_bound_reached,
                             "command exceeded the contained-process bound")
                    _require(not self.pidfd_resource_exhausted,
                             "command exhausted the reserved pidfd resource")
                    return
            else:
                empty_scans = 0
            _require(time.monotonic() < deadline,
                     "command descendants remained after bounded SIGKILL")
            time.sleep(0.01)

    def __exit__(self, kind: object, value: object, _traceback: object) -> None:
        if not self.active:
            return
        cleanup_error: BaseException | None = None
        restore_error: BaseException | None = None
        close_error: BaseException | None = None
        try:
            try:
                self._cleanup()
            except BaseException as error:
                cleanup_error = error
        finally:
            for descriptor in self.handles.values():
                try:
                    os.close(descriptor)
                except OSError as error:
                    if close_error is None:
                        close_error = error
            self.handles.clear()
            if self.reserve_descriptor >= 0:
                try:
                    os.close(self.reserve_descriptor)
                except OSError as error:
                    if close_error is None:
                        close_error = error
                self.reserve_descriptor = -1
            previous = self.previous_subreaper
            self.previous_subreaper = None
            self.active = False
            try:
                _require(previous is not None,
                         "previous child-subreaper state was lost")
                _set_linux_subreaper(previous)
            except BaseException as error:
                restore_error = error
        if cleanup_error is not None or restore_error is not None or \
                close_error is not None:
            details = []
            if cleanup_error is not None:
                details.append(
                    f"cleanup {type(cleanup_error).__name__}: {cleanup_error}")
            if restore_error is not None:
                details.append(
                    f"restore {type(restore_error).__name__}: {restore_error}")
            if close_error is not None:
                details.append(
                    f"close {type(close_error).__name__}: {close_error}")
            if value is not None:
                details.append(f"primary {type(value).__name__}: {value}")
            raise AcquisitionError(
                "bounded command containment failed: " + "; ".join(details)) \
                from (cleanup_error or restore_error or close_error)


class HostEnvironment:
    """Injectable Linux host edge for acquisition orchestration.

    The class deliberately has no retained mutable checkout state.  Tests may
    override these methods with scripted facts, while the production default
    uses exact argv, exact child environments, and byte-bounded pipe capture.
    """

    def __init__(
        self,
        *,
        anchor_path: str | None = None,
        canonical_lock_path: str = CANONICAL_LOCK_PATH,
    ) -> None:
        runtime = f"/run/user/{os.getuid()}" \
            if anchor_path is None else anchor_path
        self.anchor_path = _portable_absolute_path(
            runtime, "global evidence lease anchor")
        self.canonical_lock_path = _portable_absolute_path(
            canonical_lock_path, "canonical evidence lock")

    def independent_verifier_argv(
        self, *, interpreter: str, verifier: str, lane_root: str,
    ) -> list[str]:
        """Return the exact fresh-process independent-verifier invocation."""
        return [interpreter, "-I", "-S", "-B", verifier, lane_root]

    def run(
        self,
        argv: Sequence[str],
        *,
        cwd: str,
        env: Mapping[str, str],
        timeout: float,
        maximum_bytes: int,
    ) -> CommandResult:
        """Run one new-session child with bounded stdout+stderr capture."""
        arguments = _command_arguments(argv)
        working_directory = _portable_absolute_path(cwd, "command cwd")
        environment = _command_environment(env)
        _require(type(timeout) in (int, float) and type(timeout) is not bool and
                 math.isfinite(float(timeout)) and
                 0 < float(timeout) <= MAX_COMMAND_TIMEOUT_SECONDS,
                 "command timeout is invalid")
        _require(type(maximum_bytes) is int and
                 0 <= maximum_bytes <= MAX_COMMAND_OUTPUT_BYTES,
                 "command output bound is invalid")
        process: subprocess.Popen[bytes] | None = None
        selector = selectors.DefaultSelector()
        outputs: dict[int, bytearray] = {}
        failure: str | None = None
        try:
            with _LinuxChildContainment() as containment:
                try:
                    process = subprocess.Popen(
                        arguments,
                        cwd=working_directory,
                        env=environment,
                        stdin=subprocess.DEVNULL,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        start_new_session=True,
                        close_fds=True,
                    )
                except OSError as error:
                    raise AcquisitionError(
                        f"command could not be started: {error}") from error
                containment.attach(process)
                _require(process.stdout is not None and
                         process.stderr is not None,
                         "command output pipes were not created")
                for stream in (process.stdout, process.stderr):
                    descriptor = stream.fileno()
                    os.set_blocking(descriptor, False)
                    selector.register(stream, selectors.EVENT_READ)
                    outputs[descriptor] = bytearray()
                deadline = time.monotonic() + float(timeout)
                while selector.get_map():
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        failure = "command exceeded its timeout"
                        break
                    events = selector.select(min(remaining, 0.1))
                    for key, _event in events:
                        descriptor = key.fileobj.fileno()
                        try:
                            chunk = os.read(descriptor, READ_CHUNK)
                        except BlockingIOError:
                            continue
                        except OSError as error:
                            raise AcquisitionError(
                                f"cannot read command output: {error}") from error
                        if not chunk:
                            selector.unregister(key.fileobj)
                            continue
                        outputs[descriptor].extend(chunk)
                        if sum(len(value) for value in outputs.values()) > \
                                maximum_bytes:
                            failure = \
                                "command exceeded its combined output bound"
                            break
                    if failure is not None:
                        break
                if failure is not None:
                    _fail(failure)
                remaining = deadline - time.monotonic()
                _require(remaining > 0,
                         "command exceeded its timeout after closing streams")
                try:
                    exit_status = process.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    _fail("command exceeded its timeout after closing streams")
                stdout_descriptor = process.stdout.fileno()
                stderr_descriptor = process.stderr.fileno()
                result = CommandResult(
                    exit_status=exit_status,
                    stdout=bytes(outputs[stdout_descriptor]),
                    stderr=bytes(outputs[stderr_descriptor]),
                )
            return result
        finally:
            selector.close()
            if process is not None:
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()

    def run_to_path(
        self,
        argv: Sequence[str],
        *,
        cwd: str,
        env: Mapping[str, str],
        timeout: float,
        destination_fd: int,
        maximum_bytes: int,
    ) -> StreamedCommandResult:
        """Run one child while streaming bounded stdout to an owned file."""
        arguments = _command_arguments(argv)
        working_directory = _portable_absolute_path(cwd, "command cwd")
        environment = _command_environment(env)
        _require(type(timeout) in (int, float) and type(timeout) is not bool and
                 math.isfinite(float(timeout)) and
                 0 < float(timeout) <= MAX_COMMAND_TIMEOUT_SECONDS,
                 "command timeout is invalid")
        _require(type(destination_fd) is int and destination_fd >= 0,
                 "streamed command destination descriptor is invalid")
        _require(type(maximum_bytes) is int and
                 0 < maximum_bytes <= MAX_SEALED_FILE_BYTES,
                 "streamed command output bound is invalid")
        before = os.fstat(destination_fd)
        _require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                 before.st_uid == os.geteuid() and
                 before.st_gid == os.getegid() and
                 stat.S_IMODE(before.st_mode) == WRITABLE_FILE_MODE and
                 before.st_size == 0 and
                 os.lseek(destination_fd, 0, os.SEEK_CUR) == 0,
                 "streamed command destination is not an empty owned file")
        process: subprocess.Popen[bytes] | None = None
        selector = selectors.DefaultSelector()
        stderr = bytearray()
        stdout_size = 0
        stdout_digest = hashlib.sha256()
        failure: str | None = None
        try:
            with _LinuxChildContainment() as containment:
                try:
                    process = subprocess.Popen(
                        arguments,
                        cwd=working_directory,
                        env=environment,
                        stdin=subprocess.DEVNULL,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        start_new_session=True,
                        close_fds=True,
                    )
                except OSError as error:
                    raise AcquisitionError(
                        f"streamed command could not be started: {error}") \
                        from error
                containment.attach(process)
                _require(process.stdout is not None and
                         process.stderr is not None,
                         "streamed command pipes were not created")
                stdout_descriptor = process.stdout.fileno()
                stderr_descriptor = process.stderr.fileno()
                for stream in (process.stdout, process.stderr):
                    os.set_blocking(stream.fileno(), False)
                    selector.register(stream, selectors.EVENT_READ)
                deadline = time.monotonic() + float(timeout)
                while selector.get_map():
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        failure = "streamed command exceeded its timeout"
                        break
                    events = selector.select(min(remaining, 0.1))
                    for key, _event in events:
                        descriptor = key.fileobj.fileno()
                        try:
                            chunk = os.read(descriptor, READ_CHUNK)
                        except BlockingIOError:
                            continue
                        except OSError as error:
                            raise AcquisitionError(
                                f"cannot read streamed command output: "
                                f"{error}") from error
                        if not chunk:
                            selector.unregister(key.fileobj)
                            continue
                        if descriptor == stdout_descriptor:
                            if stdout_size + len(chunk) > maximum_bytes:
                                failure = \
                                    "streamed command exceeded its stdout bound"
                                break
                            _write_all(
                                destination_fd, chunk,
                                "streamed command stdout")
                            stdout_digest.update(chunk)
                            stdout_size += len(chunk)
                        else:
                            _require(descriptor == stderr_descriptor,
                                     "streamed command returned an unknown pipe")
                            stderr.extend(chunk)
                            if len(stderr) > MAX_COMMAND_OUTPUT_BYTES:
                                failure = \
                                    "streamed command exceeded its stderr bound"
                                break
                    if failure is not None:
                        break
                if failure is not None:
                    _fail(failure)
                remaining = deadline - time.monotonic()
                _require(remaining > 0,
                         "streamed command exceeded its timeout after output")
                try:
                    exit_status = process.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    _fail("streamed command exceeded its timeout after output")
                os.fsync(destination_fd)
                after = os.fstat(destination_fd)
                _require((after.st_dev, after.st_ino, after.st_uid,
                          after.st_gid, after.st_nlink,
                          stat.S_IMODE(after.st_mode)) ==
                         (before.st_dev, before.st_ino, before.st_uid,
                          before.st_gid, before.st_nlink,
                          stat.S_IMODE(before.st_mode)) and
                         after.st_size == stdout_size,
                         "streamed command destination changed")
                result = StreamedCommandResult(
                    exit_status=exit_status,
                    stdout_size=stdout_size,
                    stdout_sha256=stdout_digest.hexdigest(),
                    stderr=bytes(stderr),
                )
            return result
        finally:
            selector.close()
            if process is not None:
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()

    def capture_git_identity(
        self,
        controller_repository: str,
        source_root: str,
        commit: str,
        *,
        require_detached: bool,
    ) -> dict[str, Any]:
        """Delegate guarded Git capture to the source-bound controller."""
        repository = _portable_absolute_path(
            controller_repository, "Git capture controller repository")
        root = _portable_absolute_path(source_root, "Git capture source root")
        _require(type(commit) is str and len(commit) == 40 and all(
            character in "0123456789abcdef" for character in commit),
            "Git capture commit is invalid")
        _require(type(require_detached) is bool,
                 "Git capture detached policy is invalid")
        module = _load_repository_module(
            repository, "_leopard2_exact_main_git_capture",
            "experiments/leopard2/main_compare/git_capture.py")
        try:
            result = module.capture_git_identity(
                root, commit, require_detached=require_detached)
        except module.GitCaptureError as error:
            raise AcquisitionError(
                f"cannot capture Git source identity: {error}") from error
        _require(type(result) is dict,
                 "Git capture controller returned a non-object")
        return copy.deepcopy(result)

    def now_utc(self) -> str:
        return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def host_facts(self) -> dict[str, Any]:
        """Capture the exact portable host fields owned by the record schema."""
        identity = os.uname()
        cpuinfo = _read_bounded_system_file(
            "/proc/cpuinfo", MAX_CPUINFO_BYTES, "host CPU identity")
        models = _parse_linux_cpu_models(cpuinfo)
        cpus = _parse_linux_cpu_list(_read_bounded_system_file(
            "/sys/devices/system/cpu/online", 65536,
            "Linux online CPU list"))
        _require(sorted(models) == cpus and
                 len(set(models.values())) == 1,
                 "host online CPU identities are absent or inconsistent")
        try:
            ticks = os.sysconf("SC_CLK_TCK")
        except (AttributeError, OSError, ValueError) as error:
            raise AcquisitionError(
                f"cannot capture host CPU facts: {error}") from error
        facts = {
            "hostname": identity.nodename,
            "kernel": identity.release,
            "architecture": identity.machine,
            "cpu_model": models[cpus[0]],
            "online_cpus": cpus,
            "sc_clk_tck": ticks,
        }
        # The authority constructor is the ultimate schema validator.  These
        # local checks give callers a deterministic failure before any build.
        _require(identity.machine == "x86_64" and
                 type(ticks) is int and ticks > 0 and
                 cpus == sorted(set(cpus)) and bool(cpus),
                 "host facts are outside the exact-main contract")
        for label in ("hostname", "kernel", "cpu_model"):
            value = facts[label]
            _require(type(value) is str and value and
                     len(value) <= record_contract.MAX_TEXT_LENGTH and
                     all(
                         ord(character) >= 0x20 and
                         not (0x7F <= ord(character) <= 0x9F) and
                         not (0xD800 <= ord(character) <= 0xDFFF)
                         for character in value),
                     f"host {label} is invalid")
        return copy.deepcopy(facts)


class StableLeaseAnchor:
    """Nonblocking flock over a stable, owner-safe per-UID directory."""

    def __init__(self, path: str):
        self.path = _portable_absolute_path(path, "global evidence lease anchor")
        self.descriptor = -1
        self.identity: tuple[int, int, int, int] | None = None

    @staticmethod
    def _metadata(status: os.stat_result) -> tuple[int, int, int, int]:
        mode = stat.S_IMODE(status.st_mode)
        _require(stat.S_ISDIR(status.st_mode) and
                 status.st_uid == os.getuid() and
                 status.st_gid == os.getgid() and mode & 0o022 == 0,
                 "global evidence lease anchor has unsafe ownership or mode")
        return status.st_dev, status.st_ino, status.st_uid, mode

    def validate_current(self) -> None:
        _require(self.descriptor >= 0 and self.identity is not None,
                 "global evidence lease anchor is not held")
        descriptor = os.fstat(self.descriptor)
        try:
            pathname = os.lstat(self.path)
        except OSError as error:
            raise AcquisitionError(
                "global evidence lease anchor path disappeared") from error
        _require(self._metadata(descriptor) == self.identity and
                 (descriptor.st_dev, descriptor.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 "global evidence lease anchor path was replaced")

    def __enter__(self) -> "StableLeaseAnchor":
        _require(self.descriptor < 0,
                 "global evidence lease anchor is already held")
        try:
            self.descriptor = os.open(
                self.path,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
            try:
                fcntl.flock(
                    self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise AcquisitionError(
                    "global evidence lease is already held") from error
            descriptor = os.fstat(self.descriptor)
            pathname = os.lstat(self.path)
            self.identity = self._metadata(descriptor)
            _require((descriptor.st_dev, descriptor.st_ino) ==
                     (pathname.st_dev, pathname.st_ino),
                     "global evidence lease anchor changed during acquisition")
            self.validate_current()
            return self
        except BaseException as error:
            self.__exit__(None, None, None)
            if isinstance(error, OSError):
                raise AcquisitionError(
                    f"cannot acquire global evidence lease: {error}") from error
            raise

    def __exit__(self, _kind: object, _value: object,
                 _traceback: object) -> None:
        if self.descriptor >= 0:
            descriptor = self.descriptor
            self.descriptor = -1
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)
        self.identity = None


class CanonicalFileLock:
    """Exclusive canonical build/test lock with stable inode validation."""

    def __init__(self, path: str, *, blocking: bool = True):
        self.path = _portable_absolute_path(path, "canonical evidence lock")
        _require(type(blocking) is bool,
                 "canonical evidence lock policy is not boolean")
        self.blocking = blocking
        self.descriptor = -1
        self.identity: tuple[int, int, int, int, int] | None = None

    @staticmethod
    def _metadata(status: os.stat_result) -> tuple[int, int, int, int, int]:
        mode = stat.S_IMODE(status.st_mode)
        _require(stat.S_ISREG(status.st_mode) and
                 status.st_uid == os.getuid() and
                 status.st_gid == os.getgid() and
                 status.st_nlink == 1 and mode == 0o600,
                 "canonical evidence lock has unsafe metadata")
        return (status.st_dev, status.st_ino, status.st_uid,
                status.st_gid, mode)

    def validate_current(self) -> None:
        _require(self.descriptor >= 0 and self.identity is not None,
                 "canonical evidence lock is not held")
        descriptor = os.fstat(self.descriptor)
        try:
            pathname = os.lstat(self.path)
        except OSError as error:
            raise AcquisitionError(
                "canonical evidence lock path disappeared") from error
        _require(self._metadata(descriptor) == self.identity and
                 (descriptor.st_dev, descriptor.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 "canonical evidence lock path was replaced")

    def __enter__(self) -> "CanonicalFileLock":
        _require(self.descriptor < 0,
                 "canonical evidence lock is already held")
        try:
            self.descriptor = os.open(
                self.path,
                os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | os.O_CLOEXEC,
                0o600)
            operation = fcntl.LOCK_EX
            if not self.blocking:
                operation |= fcntl.LOCK_NB
            try:
                fcntl.flock(self.descriptor, operation)
            except BlockingIOError as error:
                raise AcquisitionError(
                    "canonical evidence lock is already held") from error
            descriptor = os.fstat(self.descriptor)
            pathname = os.lstat(self.path)
            self.identity = self._metadata(descriptor)
            _require((descriptor.st_dev, descriptor.st_ino) ==
                     (pathname.st_dev, pathname.st_ino),
                     "canonical evidence lock changed during acquisition")
            self.validate_current()
            return self
        except BaseException as error:
            self.__exit__(None, None, None)
            if isinstance(error, OSError):
                raise AcquisitionError(
                    f"cannot acquire canonical evidence lock: {error}") from error
            raise

    def __exit__(self, _kind: object, _value: object,
                 _traceback: object) -> None:
        if self.descriptor >= 0:
            descriptor = self.descriptor
            self.descriptor = -1
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)
        self.identity = None


class AcquisitionLocks:
    """Hold and jointly revalidate both acquisition serialization anchors."""

    def __init__(self, environment: HostEnvironment, *, blocking: bool = True):
        _require(isinstance(environment, HostEnvironment),
                 "acquisition lock environment is invalid")
        self.anchor = StableLeaseAnchor(environment.anchor_path)
        self.canonical = CanonicalFileLock(
            environment.canonical_lock_path, blocking=blocking)
        self._entered = False

    def __enter__(self) -> "AcquisitionLocks":
        _require(not self._entered, "acquisition locks are already held")
        try:
            self.anchor.__enter__()
            self.canonical.__enter__()
            self._entered = True
            self.validate_current()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def validate_current(self) -> None:
        _require(self._entered,
                 "acquisition locks are not both held")
        self.anchor.validate_current()
        self.canonical.validate_current()

    def __exit__(self, _kind: object, _value: object,
                 _traceback: object) -> None:
        try:
            self.canonical.__exit__(None, None, None)
        finally:
            self.anchor.__exit__(None, None, None)
            self._entered = False


def _bytes_identity(relative_path: str, content: bytes) -> dict[str, Any]:
    relative = _safe_relative_path(relative_path, "retained file path")
    _require(type(content) is bytes and len(content) <= MAX_SEALED_FILE_BYTES,
             f"retained file {relative!r} is not bounded bytes")
    return {
        "relative_path": relative,
        "size": len(content),
        "sha256": _sha256(content),
    }


def _owned_file_identity(
    path: str, label: str, *, maximum_bytes: int = MAX_SEALED_FILE_BYTES,
) -> dict[str, Any]:
    absolute = _portable_absolute_path(path, label)
    _require(type(maximum_bytes) is int and
             0 < maximum_bytes <= MAX_SEALED_FILE_BYTES,
             f"{label} bound is invalid")
    descriptor = -1
    try:
        descriptor = os.open(
            absolute, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.lstat(absolute)
        _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1 and
                 status.st_uid == os.geteuid() and
                 status.st_gid == os.getegid() and
                 stat.S_IMODE(status.st_mode) == WRITABLE_FILE_MODE and
                 0 < status.st_size <= maximum_bytes and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{label} is not one owned regular file")
        digest = _hash_fd(descriptor, status.st_size, label)
        after = os.fstat(descriptor)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) ==
                 (status.st_dev, status.st_ino, status.st_size,
                  status.st_mtime_ns, status.st_ctime_ns),
                 f"{label} changed while hashed")
        return {"path": absolute, "size": status.st_size, "sha256": digest}
    except OSError as error:
        raise AcquisitionError(f"cannot inspect {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _host_file_identity(
    path: str, label: str, *, maximum_bytes: int = MAX_SEALED_FILE_BYTES,
    minimum_size: int = 1,
) -> dict[str, Any]:
    """Hash one stable regular host file without owner or mode assumptions."""
    absolute = _portable_absolute_path(path, label)
    _require(type(maximum_bytes) is int and
             0 < maximum_bytes <= MAX_SEALED_FILE_BYTES and
             minimum_size in (0, 1),
             f"{label} bound is invalid")
    descriptor = -1
    try:
        descriptor = os.open(
            absolute, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.lstat(absolute)
        _require(stat.S_ISREG(status.st_mode) and
                 minimum_size <= status.st_size <= maximum_bytes and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{label} is not one stable regular file")
        digest = _hash_fd(descriptor, status.st_size, label)
        after = os.fstat(descriptor)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) ==
                 (status.st_dev, status.st_ino, status.st_size,
                  status.st_mtime_ns, status.st_ctime_ns),
                 f"{label} changed while hashed")
        return {"path": absolute, "size": status.st_size, "sha256": digest}
    except OSError as error:
        raise AcquisitionError(f"cannot inspect {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _runtime_dependency_identity(
    path: str, label: str, *, maximum_bytes: int = MAX_SEALED_FILE_BYTES,
) -> dict[str, Any]:
    """Hash the target named by one exact path printed by frozen ``ldd``.

    Dynamic-linker SONAME paths are normally symlinks.  The record must keep
    the path printed by ``ldd`` while binding the regular target bytes, so this
    narrowly scoped reader follows the final link and then proves that the
    pathname still resolves to the opened file before and after hashing.
    Other host-file readers remain no-follow boundaries.
    """
    absolute = _portable_absolute_path(path, label)
    _require(type(maximum_bytes) is int and
             0 < maximum_bytes <= MAX_SEALED_FILE_BYTES,
             f"{label} bound is invalid")
    descriptor = -1
    try:
        descriptor = os.open(absolute, os.O_RDONLY | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.stat(absolute, follow_symlinks=True)
        _require(stat.S_ISREG(status.st_mode) and
                 0 < status.st_size <= maximum_bytes and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{label} does not resolve to one stable regular file")
        digest = _hash_fd(descriptor, status.st_size, label)
        after = os.fstat(descriptor)
        pathname_after = os.stat(absolute, follow_symlinks=True)
        fingerprint = (status.st_dev, status.st_ino, status.st_size,
                       status.st_mtime_ns, status.st_ctime_ns)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) == fingerprint and
                 (pathname_after.st_dev, pathname_after.st_ino,
                  pathname_after.st_size, pathname_after.st_mtime_ns,
                  pathname_after.st_ctime_ns) == fingerprint,
                 f"{label} changed while hashed")
        return {"path": absolute, "size": status.st_size, "sha256": digest}
    except OSError as error:
        raise AcquisitionError(f"cannot inspect {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _read_host_file(
    path: str, label: str, *, maximum_bytes: int,
) -> bytes:
    identity = _host_file_identity(path, label, maximum_bytes=maximum_bytes)
    descriptor = -1
    try:
        descriptor = os.open(
            identity["path"], os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        content = bytearray()
        while len(content) < identity["size"]:
            chunk = os.read(
                descriptor, min(READ_CHUNK, identity["size"] - len(content)))
            _require(bool(chunk), f"{label} ended while read")
            content.extend(chunk)
        _require(os.read(descriptor, 1) == b"" and
                 _sha256(bytes(content)) == identity["sha256"],
                 f"{label} changed after hashing")
        return bytes(content)
    except OSError as error:
        raise AcquisitionError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def stage_build_output(
    source_path: str,
    destination_path: str,
    label: str,
    *,
    maximum_bytes: int = MAX_SEALED_FILE_BYTES,
) -> dict[str, Any]:
    """Copy one stable build artifact into an owner-only scratch file."""
    source = _portable_absolute_path(source_path, f"{label} source")
    destination = _portable_absolute_path(
        destination_path, f"{label} destination")
    _require(type(maximum_bytes) is int and
             0 < maximum_bytes <= MAX_SEALED_FILE_BYTES,
             f"{label} bound is invalid")
    source_descriptor = -1
    destination_descriptor = -1
    created = False
    published = False
    try:
        source_descriptor = os.open(
            source, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        before = os.fstat(source_descriptor)
        pathname = os.lstat(source)
        _require(stat.S_ISREG(before.st_mode) and before.st_nlink >= 1 and
                 0 < before.st_size <= maximum_bytes and
                 (before.st_dev, before.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{label} source is not one stable regular file")
        destination_descriptor = _exclusive_output_file(
            destination, f"{label} staged output")
        created = True
        digest = hashlib.sha256()
        remaining = before.st_size
        while remaining:
            chunk = os.read(
                source_descriptor, min(READ_CHUNK, remaining))
            _require(bool(chunk), f"{label} source ended while copied")
            _write_all(destination_descriptor, chunk, f"{label} staged output")
            digest.update(chunk)
            remaining -= len(chunk)
        _require(os.read(source_descriptor, 1) == b"",
                 f"{label} source grew while copied")
        after = os.fstat(source_descriptor)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) ==
                 (before.st_dev, before.st_ino, before.st_size,
                  before.st_mtime_ns, before.st_ctime_ns),
                 f"{label} source changed while copied")
        os.fsync(destination_descriptor)
        os.close(destination_descriptor)
        destination_descriptor = -1
        identity = _owned_file_identity(
            destination, f"{label} staged output", maximum_bytes=maximum_bytes)
        _require(identity["size"] == before.st_size and
                 identity["sha256"] == digest.hexdigest(),
                 f"{label} staged output changed after copy")
        published = True
        return identity
    except OSError as error:
        raise AcquisitionError(f"cannot stage {label}: {error}") from error
    finally:
        if source_descriptor >= 0:
            os.close(source_descriptor)
        if destination_descriptor >= 0:
            os.close(destination_descriptor)
        if created and not published:
            try:
                os.unlink(destination)
            except FileNotFoundError:
                pass
            except OSError as error:
                raise AcquisitionError(
                    f"cannot remove failed {label} staging output: {error}") \
                    from error


def build_root_census(
    build_root: str, role: str,
) -> tuple[dict[str, Any], bytes]:
    """Derive the bounded exact regular-file census of one build tree."""
    root = _portable_absolute_path(build_root, "build census root")
    _require(role in record_contract.BUILD_ROLES and type(role) is str,
             "build census role is invalid")
    try:
        root_status = os.lstat(root)
    except OSError as error:
        raise AcquisitionError(f"cannot inspect build census root: {error}") \
            from error
    _require(stat.S_ISDIR(root_status.st_mode) and
             not stat.S_ISLNK(root_status.st_mode),
             "build census root is not a directory")
    files: list[dict[str, Any]] = []
    node_count = 0
    total_bytes = 0

    def visit(directory: str, relative_directory: str, depth: int) -> None:
        nonlocal node_count, total_bytes
        _require(depth <= MAX_TREE_DEPTH,
                 "build closure exceeds the depth bound")
        try:
            with os.scandir(directory) as iterator:
                entries = sorted(iterator, key=lambda item: item.name)
        except OSError as error:
            raise AcquisitionError(
                f"cannot scan build closure directory: {error}") from error
        for entry in entries:
            node_count += 1
            _require(node_count <= record_contract.MAX_CLOSURE_FILES,
                     "build closure exceeds the node bound")
            relative = entry.name if not relative_directory else \
                relative_directory + "/" + entry.name
            _safe_relative_path(relative, "build closure path")
            status = entry.stat(follow_symlinks=False)
            _require(not stat.S_ISLNK(status.st_mode),
                     "build closure contains a symbolic link")
            if stat.S_ISDIR(status.st_mode):
                visit(entry.path, relative, depth + 1)
                continue
            _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1,
                     "build closure contains a non-regular or linked file")
            identity = _host_file_identity(
                entry.path, f"build closure file {relative}", minimum_size=0)
            total_bytes += identity["size"]
            _require(total_bytes <= MAX_SEALED_TOTAL_BYTES,
                     "build closure exceeds the total-byte bound")
            files.append({
                "relative_path": relative,
                "size": identity["size"],
                "sha256": identity["sha256"],
            })

    visit(root, "", 0)
    _require(0 < len(files) <= record_contract.MAX_CLOSURE_FILES,
             "build closure contains no regular files")
    files.sort(key=lambda item: item["relative_path"])
    closure = build_closure_document(role, root, files)
    content = canonical_json_bytes(closure)
    _require(0 < len(content) <= MAX_SEALED_FILE_BYTES,
             "build closure document exceeds its byte bound")
    return closure, content


def _read_owned_file(
    path: str, label: str, *, maximum_bytes: int,
) -> bytes:
    identity = _owned_file_identity(
        path, label, maximum_bytes=maximum_bytes)
    descriptor = -1
    try:
        descriptor = os.open(
            identity["path"], os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        content = bytearray()
        while len(content) < identity["size"]:
            chunk = os.read(
                descriptor, min(READ_CHUNK, identity["size"] - len(content)))
            _require(bool(chunk), f"{label} ended while read")
            content.extend(chunk)
        _require(os.read(descriptor, 1) == b"" and
                 _sha256(bytes(content)) == identity["sha256"],
                 f"{label} changed after hashing")
        return bytes(content)
    except OSError as error:
        raise AcquisitionError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _owned_files_identical(left: str, right: str, label: str) -> bool:
    left_identity = _owned_file_identity(left, f"{label} first file")
    right_identity = _owned_file_identity(right, f"{label} second file")
    if (left_identity["size"], left_identity["sha256"]) != \
            (right_identity["size"], right_identity["sha256"]):
        return False
    left_descriptor = -1
    right_descriptor = -1
    try:
        left_descriptor = os.open(
            left, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        right_descriptor = os.open(
            right, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        remaining = left_identity["size"]
        while remaining:
            count = min(READ_CHUNK, remaining)
            left_chunk = os.read(left_descriptor, count)
            right_chunk = os.read(right_descriptor, count)
            _require(len(left_chunk) == count and len(right_chunk) == count,
                     f"{label} changed while compared")
            if left_chunk != right_chunk:
                return False
            remaining -= count
        _require(os.read(left_descriptor, 1) == b"" and
                 os.read(right_descriptor, 1) == b"",
                 f"{label} grew while compared")
        return True
    except OSError as error:
        raise AcquisitionError(f"cannot compare {label}: {error}") from error
    finally:
        if left_descriptor >= 0:
            os.close(left_descriptor)
        if right_descriptor >= 0:
            os.close(right_descriptor)


def _one_line_output(content: bytes, label: str) -> str:
    _require(type(content) is bytes and content.endswith(b"\n") and
             content.count(b"\n") == 1 and b"\0" not in content and
             b"\r" not in content,
             f"{label} is not one LF-terminated line")
    try:
        value = content[:-1].decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise AcquisitionError(f"{label} is not strict UTF-8") from error
    _require(value and len(value) <= record_contract.MAX_TEXT_LENGTH and
             all(ord(character) >= 0x20 and
                 not (0x7F <= ord(character) <= 0x9F) and
                 not (0xD800 <= ord(character) <= 0xDFFF)
                 for character in value),
             f"{label} contains unsafe text")
    return value


def _command_log_record(
    argv: Sequence[str], cwd: str, result: CommandResult | StreamedCommandResult,
) -> dict[str, Any]:
    if isinstance(result, StreamedCommandResult):
        stdout_size = result.stdout_size
        stdout_sha256 = result.stdout_sha256
    else:
        stdout_size = len(result.stdout)
        stdout_sha256 = _sha256(result.stdout)
    return {
        "argv": _command_arguments(argv),
        "cwd": _portable_absolute_path(cwd, "logged command cwd"),
        "exit_status": result.exit_status,
        "stdout": {"size": stdout_size, "sha256": stdout_sha256},
        "stderr": {"size": len(result.stderr),
                   "sha256": _sha256(result.stderr)},
    }


def _command_transcript_bytes(
    argv: Sequence[str], cwd: str, result: CommandResult,
) -> bytes:
    """Serialize exact command streams without timing fields."""
    _require(isinstance(result, CommandResult),
             "command transcript result has the wrong type")
    arguments = _command_arguments(argv)
    working_directory = _portable_absolute_path(cwd, "command transcript cwd")
    _require(type(result.exit_status) is int,
             "command transcript exit status is not integral")

    def stream(content: bytes) -> dict[str, Any]:
        _require(type(content) is bytes and
                 len(content) <= MAX_COMMAND_OUTPUT_BYTES,
                 "command transcript stream exceeds its bound")
        return {
            "size": len(content),
            "sha256": _sha256(content),
            "base64": base64.b64encode(content).decode("ascii"),
        }

    content = canonical_json_bytes({
        "schema": BUILD_LOG_SCHEMA,
        "argv": arguments,
        "cwd": working_directory,
        "exit_status": result.exit_status,
        "stdout": stream(result.stdout),
        "stderr": stream(result.stderr),
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "command transcript exceeds its byte bound")
    return content


def _build_stage_log_bytes(
    role: str,
    status: str,
    commands: Sequence[Mapping[str, Any]],
    *,
    controller_sha256: str,
    error: BaseException | None,
) -> bytes:
    _require(role in record_contract.BUILD_ROLES and type(role) is str,
             "build stage role is invalid")
    _require(status in ("complete", "failed"),
             "build stage log status is invalid")
    _require(type(controller_sha256) is str and
             len(controller_sha256) == 64 and all(
                 character in "0123456789abcdef"
                 for character in controller_sha256),
             "build stage controller SHA-256 is invalid")
    content = canonical_json_bytes({
        "schema": BUILD_STAGE_LOG_SCHEMA,
        "role": role,
        "stage": role + "_build",
        "status": status,
        "command_count": len(commands),
        "commands": copy.deepcopy(list(commands)),
        "controller_sha256": controller_sha256,
        "error": None if error is None else _safe_failure_message(error),
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "build stage log exceeds its byte bound")
    return content


def _checked_run(
    environment: HostEnvironment,
    argv: Sequence[str],
    *,
    cwd: str,
    child_environment: Mapping[str, str],
    log: list[dict[str, Any]],
    label: str,
    expected_statuses: tuple[int, ...] = (0,),
    timeout: float = SOURCE_ACQUISITION_TIMEOUT_SECONDS,
) -> CommandResult:
    result = environment.run(
        argv, cwd=cwd, env=child_environment, timeout=timeout,
        maximum_bytes=MAX_COMMAND_OUTPUT_BYTES)
    log.append(_command_log_record(argv, cwd, result))
    _require(type(result.exit_status) is int,
             f"{label} returned a non-integral status")
    if result.exit_status not in expected_statuses:
        detail = result.stderr or result.stdout
        flattened = " ".join(
            detail.decode("utf-8", errors="replace").replace("\r", " ").
            replace("\n", " ").split())[:1024]
        raise CommandExecutionError(
            f"{label} failed with status {result.exit_status}" +
            (f": {flattened}" if flattened else ""),
            _normalized_exit_status(result.exit_status))
    return result


def _normalized_exit_status(status: int) -> int:
    _require(type(status) is int, "command exit status is not integral")
    if status < 0:
        return min(255, 128 + abs(status))
    return min(255, max(1, status))


def _repository_source_identity(
    environment: HostEnvironment,
    plan: LanePlan,
    child_environment: Mapping[str, str],
    log: list[dict[str, Any]],
) -> tuple[str, str, str]:
    git = "/usr/bin/git"
    commit_result = _checked_run(
        environment,
        [git, "-C", plan.repository, "rev-parse", "--verify",
         "HEAD^{commit}"],
        cwd=plan.repository, child_environment=child_environment, log=log,
        label="adapter repository commit query")
    tree_result = _checked_run(
        environment,
        [git, "-C", plan.repository, "rev-parse", "--verify",
         "HEAD^{tree}"],
        cwd=plan.repository, child_environment=child_environment, log=log,
        label="adapter repository tree query")
    status_result = _checked_run(
        environment,
        [git, "-C", plan.repository, "status", "--porcelain=v1",
         "--untracked-files=normal", "--ignore-submodules=none"],
        cwd=plan.repository, child_environment=child_environment, log=log,
        label="adapter repository status query")
    _require(status_result.stdout == b"" and status_result.stderr == b"",
             "adapter repository is not clean")
    commit = _one_line_output(
        commit_result.stdout, "adapter repository commit")
    tree = _one_line_output(tree_result.stdout, "adapter repository tree")
    _require(len(commit) == 40 and len(tree) == 40 and all(
        character in "0123456789abcdef" for character in commit + tree),
        "adapter repository identities are not lowercase Git IDs")
    submodule_result = _checked_run(
        environment,
        [git, "-C", plan.repository, "ls-tree", commit, "--", "sse2neon"],
        cwd=plan.repository, child_environment=child_environment, log=log,
        label="adapter repository submodule query")
    _require(submodule_result.stdout.endswith(b"\n") and
             submodule_result.stdout.count(b"\n") == 1 and
             b"\0" not in submodule_result.stdout and
             b"\r" not in submodule_result.stdout,
             "adapter repository submodule is not one LF record")
    try:
        row = submodule_result.stdout[:-1].decode("ascii", errors="strict")
    except UnicodeDecodeError as error:
        raise AcquisitionError(
            "adapter repository submodule is not strict ASCII") from error
    match = re.fullmatch(
        r"160000 commit ([0-9a-f]{40})\tsse2neon", row)
    _require(match is not None,
             "adapter repository sse2neon gitlink changed")
    return commit, tree, match.group(1)


def stage_detached_source(
    environment: HostEnvironment,
    *,
    source_repository: str,
    submodule_repository: str,
    destination: str,
    commit: str,
    tree: str,
    submodule_commit: str,
    child_environment: Mapping[str, str],
    log: list[dict[str, Any]],
) -> None:
    """Clone and prove one clean detached source tree without network access."""
    source = _portable_absolute_path(
        source_repository, "detached source repository")
    submodule_source = _portable_absolute_path(
        submodule_repository, "detached submodule repository")
    root = _portable_absolute_path(destination, "detached source root")
    for value, label in ((commit, "commit"), (tree, "tree"),
                         (submodule_commit, "submodule commit")):
        _require(type(value) is str and len(value) == 40 and all(
            character in "0123456789abcdef" for character in value),
            f"detached source {label} is invalid")
    try:
        with os.scandir(root) as entries:
            _require(next(entries, None) is None,
                     "detached source destination is not empty")
    except OSError as error:
        raise AcquisitionError(
            f"cannot inspect detached source destination: {error}") from error
    git = "/usr/bin/git"
    _checked_run(
        environment,
        [git, "clone", "--no-hardlinks", "--no-checkout", source, root],
        cwd=str(Path(root).parent), child_environment=child_environment,
        log=log, label="detached source clone")
    _checked_run(
        environment,
        [git, "-C", root, "checkout", "--detach", commit],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source checkout")
    submodule_root = root + "/sse2neon"
    if not _path_is_absent(submodule_root):
        try:
            submodule_status = os.lstat(submodule_root)
            with os.scandir(submodule_root) as entries:
                submodule_empty = next(entries, None) is None
        except OSError as error:
            raise AcquisitionError(
                f"cannot inspect detached submodule destination: {error}") \
                from error
        _require(stat.S_ISDIR(submodule_status.st_mode) and
                 not stat.S_ISLNK(submodule_status.st_mode) and
                 submodule_status.st_uid == os.geteuid() and
                 submodule_status.st_gid == os.getegid() and submodule_empty,
                 "detached source submodule path is not one empty directory")
    _checked_run(
        environment,
        [git, "clone", "--no-hardlinks", "--no-checkout",
         submodule_source, submodule_root],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source submodule clone")
    _checked_run(
        environment,
        [git, "-C", submodule_root, "checkout", "--detach",
         submodule_commit],
        cwd=submodule_root, child_environment=child_environment, log=log,
        label="detached source submodule checkout")
    observed_commit = _one_line_output(_checked_run(
        environment,
        [git, "-C", root, "rev-parse", "--verify", "HEAD^{commit}"],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source commit verification").stdout,
        "detached source observed commit")
    observed_tree = _one_line_output(_checked_run(
        environment,
        [git, "-C", root, "rev-parse", "--verify", "HEAD^{tree}"],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source tree verification").stdout,
        "detached source observed tree")
    observed_submodule = _one_line_output(_checked_run(
        environment,
        [git, "-C", root + "/sse2neon", "rev-parse", "--verify",
         "HEAD^{commit}"],
        cwd=root + "/sse2neon", child_environment=child_environment,
        log=log, label="detached source submodule verification").stdout,
        "detached source observed submodule")
    _require((observed_commit, observed_tree, observed_submodule) ==
             (commit, tree, submodule_commit),
             "detached source identity differs from its frozen request")
    symbolic = _checked_run(
        environment,
        [git, "-C", root, "symbolic-ref", "-q", "HEAD"],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source HEAD verification", expected_statuses=(1,))
    _require(symbolic.stdout == b"" and symbolic.stderr == b"",
             "detached source HEAD verification emitted output")
    submodule_symbolic = _checked_run(
        environment,
        [git, "-C", submodule_root, "symbolic-ref", "-q", "HEAD"],
        cwd=submodule_root, child_environment=child_environment, log=log,
        label="detached source submodule HEAD verification",
        expected_statuses=(1,))
    _require(submodule_symbolic.stdout == b"" and
             submodule_symbolic.stderr == b"",
             "detached source submodule HEAD verification emitted output")
    status_result = _checked_run(
        environment,
        [git, "-C", root, "status", "--porcelain=v1",
         "--untracked-files=normal", "--ignore-submodules=none"],
        cwd=root, child_environment=child_environment, log=log,
        label="detached source cleanliness verification")
    _require(status_result.stdout == b"" and status_result.stderr == b"",
             "detached source is not clean")


def _exclusive_output_file(path: str, label: str) -> int:
    absolute = _portable_absolute_path(path, label)
    _require(_path_is_absent(absolute), f"{label} already exists")
    try:
        return os.open(
            absolute,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW |
            os.O_CLOEXEC,
            WRITABLE_FILE_MODE,
        )
    except OSError as error:
        raise AcquisitionError(f"cannot create {label}: {error}") from error


def canonical_git_archive(
    environment: HostEnvironment,
    *,
    source_repository: str,
    commit: str,
    prefix: str,
    scratch_root: str,
    destination_name: str,
    child_environment: Mapping[str, str],
    log: list[dict[str, Any]],
) -> dict[str, Any]:
    """Write one config-isolated canonical Git archive into scratch."""
    source = _portable_absolute_path(
        source_repository, "canonical archive source repository")
    scratch = _portable_absolute_path(
        scratch_root, "canonical archive scratch root")
    _require(type(commit) is str and len(commit) == 40 and all(
        character in "0123456789abcdef" for character in commit),
        "canonical archive commit is invalid")
    _require(prefix in (
        "leopard-main-source/", "leopard2-adapter-source/",
        "sse2neon-source/"), "canonical archive prefix is invalid")
    name = _safe_relative_path(
        destination_name, "canonical archive destination name")
    _require("/" not in name and name.endswith(".tar"),
             "canonical archive destination is not one tar basename")
    destination = scratch + "/" + name
    archive_root = ""
    descriptor = -1
    try:
        archive_root = tempfile.mkdtemp(
            prefix=".leopard-canonical-git-archive.", dir=scratch)
        _require(os.path.realpath(archive_root) == archive_root,
                 "canonical archive work root is not canonical")
        empty_template = archive_root + "/empty-template"
        os.mkdir(empty_template, WRITABLE_DIRECTORY_MODE)
        git = "/usr/bin/git"
        common_result = _checked_run(
            environment,
            [git, "-C", source, "rev-parse", "--path-format=absolute",
             "--git-common-dir"],
            cwd=source, child_environment={
                **child_environment, "GIT_ATTR_NOSYSTEM": "1"},
            log=log, label="canonical archive object-store query")
        common = _one_line_output(
            common_result.stdout, "canonical archive common Git directory")
        common = os.path.realpath(common)
        objects = os.path.realpath(common + "/objects")
        _require(Path(objects).is_dir() and ":" not in objects,
                 "canonical archive object store is invalid")
        bare = archive_root + "/repository.git"
        _checked_run(
            environment,
            [git, "init", "--bare", "--quiet",
             "--template=" + empty_template, bare],
            cwd=archive_root, child_environment={
                **child_environment, "GIT_ATTR_NOSYSTEM": "1"},
            log=log, label="canonical archive bare repository creation")
        _require(_path_is_absent(bare + "/info/attributes"),
                 "canonical archive bare repository gained attributes")
        archive_environment = {
            **child_environment,
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": objects,
            "GIT_ATTR_NOSYSTEM": "1",
            "GIT_ATTR_SOURCE": commit,
            "GIT_DIR": bare,
        }
        observed = _one_line_output(_checked_run(
            environment,
            [git, "-c", "core.attributesFile=/dev/null",
             "-c", "tar.umask=0002", "rev-parse", "--verify",
             commit + "^{commit}"],
            cwd=archive_root, child_environment=archive_environment,
            log=log, label="canonical archive commit verification").stdout,
            "canonical archive observed commit")
        _require(observed == commit,
                 "canonical archive resolved another commit")
        descriptor = _exclusive_output_file(
            destination, "canonical archive output")
        argv = [
            git, "-c", "core.attributesFile=/dev/null",
            "-c", "tar.umask=0002", "archive", "--format=tar",
            "--prefix=" + prefix, commit,
        ]
        result = environment.run_to_path(
            argv, cwd=archive_root, env=archive_environment,
            timeout=SOURCE_ACQUISITION_TIMEOUT_SECONDS,
            destination_fd=descriptor,
            maximum_bytes=MAX_SEALED_FILE_BYTES)
        log.append(_command_log_record(argv, archive_root, result))
        if result.exit_status != 0:
            raise CommandExecutionError(
                f"canonical archive failed with status {result.exit_status}",
                _normalized_exit_status(result.exit_status))
        _require(result.stdout_size > 0 and result.stderr == b"",
                 "canonical archive produced empty bytes or diagnostics")
        os.close(descriptor)
        descriptor = -1
        identity = _owned_file_identity(
            destination, "canonical archive output")
        _require(identity["size"] == result.stdout_size and
                 identity["sha256"] == result.stdout_sha256,
                 "canonical archive output changed after capture")
        return identity
    except OSError as error:
        raise AcquisitionError(
            f"cannot produce canonical Git archive: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if archive_root:
            try:
                shutil.rmtree(archive_root)
            except FileNotFoundError:
                pass
            except OSError as error:
                raise AcquisitionError(
                    f"cannot remove canonical archive work root: {error}") \
                    from error


def _read_source_regular(path: str, label: str) -> bytes:
    absolute = _portable_absolute_path(path, label)
    descriptor = -1
    try:
        descriptor = os.open(
            absolute, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.lstat(absolute)
        _require(stat.S_ISREG(status.st_mode) and status.st_nlink >= 1 and
                 0 < status.st_size <= MAX_SEALED_FILE_BYTES and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{label} is not one stable regular file")
        chunks: list[bytes] = []
        total = 0
        while total < status.st_size:
            chunk = os.read(
                descriptor, min(READ_CHUNK, status.st_size - total))
            _require(bool(chunk), f"{label} ended while read")
            chunks.append(chunk)
            total += len(chunk)
        _require(os.read(descriptor, 1) == b"",
                 f"{label} grew while read")
        after = os.fstat(descriptor)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) ==
                 (status.st_dev, status.st_ino, status.st_size,
                  status.st_mtime_ns, status.st_ctime_ns),
                 f"{label} changed while read")
        return b"".join(chunks)
    except OSError as error:
        raise AcquisitionError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _git_blob_sha1(content: bytes) -> str:
    digest = hashlib.sha1(usedforsecurity=False)
    digest.update(f"blob {len(content)}\0".encode("ascii"))
    digest.update(content)
    return digest.hexdigest()


def adapter_inventory(
    adapter_root: str,
    capture: Mapping[str, Any],
) -> tuple[dict[str, Any], bytes]:
    """Bind the three frozen adapter files to worktree and capture bytes."""
    root = _portable_absolute_path(adapter_root, "adapter inventory root")
    _require(type(capture) is dict and
             type(capture.get("tracked_files")) is list,
             "adapter Git capture lacks its tracked inventory")
    tracked: dict[str, Mapping[str, Any]] = {}
    for item in capture["tracked_files"]:
        _require(type(item) is dict and type(item.get("path")) is str and
                 item["path"] not in tracked,
                 "adapter Git capture has an invalid tracked path")
        tracked[item["path"]] = item
    files: list[dict[str, Any]] = []
    controller = b""
    for path in record_contract.ADAPTER_PATHS:
        _safe_relative_path(path, "adapter inventory path")
        _require(path in tracked and tracked[path].get("kind") == "regular" and
                 tracked[path].get("git_mode") in ("100644", "100755") and
                 type(tracked[path].get("object_id")) is str,
                 f"adapter Git capture does not bind {path!r}")
        content = _read_source_regular(
            root + "/" + path, f"adapter file {path}")
        blob = _git_blob_sha1(content)
        _require(blob == tracked[path]["object_id"],
                 f"adapter file {path!r} differs from its Git blob")
        record = {
            "path": path,
            "git_blob_sha1": blob,
            "size": len(content),
            "sha256": _sha256(content),
        }
        files.append(record)
        if path == record_contract.ADAPTER_PATHS[2]:
            controller = content
    _require(bool(controller), "adapter attestation controller is empty")
    return ({
        "minimum_harness_commit": record_contract.MINIMUM_HARNESS_COMMIT,
        "files": copy.deepcopy(files),
        "files_sha256": identity_contract.canonical_sha256(files),
    }, controller)


def _archive_record(
    relative_path: str,
    prefix: str,
    original: Mapping[str, Any],
    replay: Mapping[str, Any],
) -> dict[str, Any]:
    _require(type(original) is dict and type(replay) is dict and
             set(original) == {"path", "size", "sha256"} and
             set(replay) == {"path", "size", "sha256"} and
             original["size"] == replay["size"] and
             original["sha256"] == replay["sha256"] and
             _owned_files_identical(
                 original["path"], replay["path"],
                 f"{relative_path} canonical archive replay"),
             f"{relative_path} was not replayed byte-identically")
    return {
        "relative_path": _safe_relative_path(
            relative_path, "source archive retained path"),
        "prefix": prefix,
        "size": original["size"],
        "sha256": original["sha256"],
        "replay_sha256": replay["sha256"],
        "replay_identical": True,
    }


def capture_source_evidence(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    adapter_commit: str,
    adapter_tree: str,
    log: list[dict[str, Any]],
    child_environment: Mapping[str, str],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, bytes],
           dict[str, str], tuple[dict[str, Any], dict[str, Any]]]:
    """Capture canonical source archives, Git identities, and adapter files."""
    baseline_original = canonical_git_archive(
        environment, source_repository=plan.repository,
        commit=record_contract.BASELINE_COMMIT,
        prefix="leopard-main-source/", scratch_root=plan.scratch_root,
        destination_name="leopard-main-source.original.tar",
        child_environment=child_environment, log=log)
    baseline_replay = canonical_git_archive(
        environment, source_repository=plan.canonical_baseline_root,
        commit=record_contract.BASELINE_COMMIT,
        prefix="leopard-main-source/", scratch_root=plan.scratch_root,
        destination_name="leopard-main-source.replay.tar",
        child_environment=child_environment, log=log)
    adapter_original = canonical_git_archive(
        environment, source_repository=plan.repository,
        commit=adapter_commit, prefix="leopard2-adapter-source/",
        scratch_root=plan.scratch_root,
        destination_name="leopard2-adapter-source.original.tar",
        child_environment=child_environment, log=log)
    adapter_replay = canonical_git_archive(
        environment, source_repository=plan.canonical_adapter_root,
        commit=adapter_commit, prefix="leopard2-adapter-source/",
        scratch_root=plan.scratch_root,
        destination_name="leopard2-adapter-source.replay.tar",
        child_environment=child_environment, log=log)
    baseline_capture = environment.capture_git_identity(
        plan.repository, plan.canonical_baseline_root,
        record_contract.BASELINE_COMMIT, require_detached=True)
    adapter_capture = environment.capture_git_identity(
        plan.repository, plan.canonical_adapter_root,
        adapter_commit, require_detached=True)
    _require(baseline_capture.get("head") == record_contract.BASELINE_COMMIT and
             baseline_capture.get("tree") == record_contract.BASELINE_TREE and
             baseline_capture.get("path") == plan.canonical_baseline_root and
             baseline_capture.get("detached") is True and
             baseline_capture.get("tracked_status") == "clean",
             "Leopard1 Git capture differs from its frozen source")
    baseline_submodules = baseline_capture.get("submodules")
    _require(type(baseline_submodules) is list and
             len(baseline_submodules) == 1 and
             baseline_submodules[0].get("path") == "sse2neon" and
             baseline_submodules[0].get("object_id") ==
             record_contract.BASELINE_SSE2NEON_COMMIT,
             "Leopard1 Git capture has another sse2neon identity")
    _require(adapter_capture.get("head") == adapter_commit and
             adapter_capture.get("tree") == adapter_tree and
             adapter_capture.get("path") == plan.canonical_adapter_root and
             adapter_capture.get("detached") is True and
             adapter_capture.get("tracked_status") == "clean",
             "adapter Git capture differs from its frozen source")
    baseline_capture_bytes = canonical_json_bytes(baseline_capture)
    adapter_capture_bytes = canonical_json_bytes(adapter_capture)
    adapter, controller = adapter_inventory(
        plan.canonical_adapter_root, adapter_capture)
    baseline_archive = _archive_record(
        "source/leopard-main-source.tar", "leopard-main-source/",
        baseline_original, baseline_replay)
    adapter_archive = _archive_record(
        "source/leopard2-adapter-source.tar", "leopard2-adapter-source/",
        adapter_original, adapter_replay)
    source = {
        "baseline": {
            "commit": record_contract.BASELINE_COMMIT,
            "tree": record_contract.BASELINE_TREE,
            "submodule": {
                "path": "sse2neon",
                "commit": record_contract.BASELINE_SSE2NEON_COMMIT,
            },
            "git_capture": _bytes_identity(
                "source/leopard-main-git-capture.json",
                baseline_capture_bytes),
            "archive": baseline_archive,
        },
        "adapter_repository": {
            "commit": adapter_commit,
            "tree": adapter_tree,
            "clean": True,
            "git_capture": _bytes_identity(
                "source/adapter-git-capture.json", adapter_capture_bytes),
            "archive": adapter_archive,
        },
    }
    retained_bytes = {
        "source/leopard-main-git-capture.json": baseline_capture_bytes,
        "source/adapter-git-capture.json": adapter_capture_bytes,
        "controllers/test_legacy_main_benchmark.py": controller,
    }
    retained_paths = {
        baseline_archive["relative_path"]: baseline_original["path"],
        adapter_archive["relative_path"]: adapter_original["path"],
    }
    return (copy.deepcopy(source), copy.deepcopy(adapter), retained_bytes,
            retained_paths, (baseline_capture, adapter_capture))


def _tool_requested_paths() -> dict[str, str]:
    return {
        "archiver": "/usr/bin/ar",
        "cmake": "/usr/bin/cmake",
        "compiler": "/usr/bin/c++",
        "ctest": "/usr/bin/ctest",
        "git": "/usr/bin/git",
        "ldd": "/usr/bin/ldd",
        "linker": "/usr/bin/ld",
        "make": "/usr/bin/make",
        "python": _portable_absolute_path(
            sys.executable, "acquisition Python executable"),
        "ranlib": "/usr/bin/ranlib",
    }


def _tool_record(role: str, path: str) -> dict[str, Any]:
    requested = _portable_absolute_path(path, f"{role} tool path")
    resolved = os.path.realpath(requested)
    _portable_absolute_path(resolved, f"{role} resolved tool path")
    descriptor = -1
    try:
        descriptor = os.open(
            resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.stat(requested)
        _require(stat.S_ISREG(status.st_mode) and status.st_size > 0 and
                 status.st_size <= MAX_SEALED_FILE_BYTES and
                 stat.S_IMODE(status.st_mode) & 0o111 != 0 and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{role} tool is not one executable regular file")
        digest = _hash_fd(descriptor, status.st_size, f"{role} tool")
        return {
            "role": role,
            "path": requested,
            "resolved_path": resolved,
            "size": status.st_size,
            "mode": stat.S_IMODE(status.st_mode),
            "sha256": digest,
        }
    except OSError as error:
        raise AcquisitionError(f"cannot capture {role} tool: {error}") \
            from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def resolve_toolchain(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    child_environment: Mapping[str, str],
    log: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, bytes]]:
    """Resolve, hash, version, and retain the frozen build tool inventory."""
    requested = _tool_requested_paths()
    _require(tuple(requested) == record_contract.TOOL_ROLES,
             "producer tool role order differs from the record contract")
    tools = [_tool_record(role, requested[role])
             for role in record_contract.TOOL_ROLES]
    compiler = next(item for item in tools if item["role"] == "compiler")
    subtool_programs = {
        "assembler": "as", "cc1plus": "cc1plus", "collect2": "collect2"}
    subtools: list[dict[str, Any]] = []
    for role in record_contract.SUBTOOL_ROLES:
        program = subtool_programs[role]
        result = _checked_run(
            environment,
            [compiler["resolved_path"], "-print-prog-name=" + program],
            cwd=plan.repository, child_environment=child_environment,
            log=log, label=f"{role} compiler-subtool query",
            timeout=TOOL_COMMAND_TIMEOUT_SECONDS)
        observed = _one_line_output(
            result.stdout, f"{role} compiler-subtool path")
        _require(result.stderr == b"" and
                 ("/" not in observed or observed.startswith("/")),
                 f"{role} compiler-subtool query returned diagnostics")
        if not observed.startswith("/"):
            located = shutil.which(observed, path="/usr/bin:/bin")
            _require(located is not None,
                     f"{role} compiler subtool is not in the frozen path")
            observed = located
        subtools.append(_tool_record(role, observed))
    role_records = {
        item["role"]: item for item in tools + subtools}
    _require(tuple(sorted(role_records)) == record_contract.VERSION_ROLES,
             "producer version role inventory differs from the contract")
    retained: dict[str, bytes] = {}
    versions: list[dict[str, Any]] = []
    for role in record_contract.VERSION_ROLES:
        argv = [role_records[role]["resolved_path"], "--version"]
        result = _checked_run(
            environment, argv, cwd=plan.repository,
            child_environment=child_environment, log=log,
            label=f"{role} version query",
            timeout=TOOL_COMMAND_TIMEOUT_SECONDS)
        stdout_path = f"toolchain/versions/{role}.stdout"
        stderr_path = f"toolchain/versions/{role}.stderr"
        retained[stdout_path] = result.stdout
        retained[stderr_path] = result.stderr
        versions.append({
            "role": role,
            "argv": argv,
            "stdout": _bytes_identity(stdout_path, result.stdout),
            "stderr": _bytes_identity(stderr_path, result.stderr),
            "exit_status": 0,
        })
    return ({
        "tools": copy.deepcopy(tools),
        "subtools": copy.deepcopy(subtools),
        "versions": versions,
    }, retained)


def _capture_git_tool_binding(
    captures: Sequence[Mapping[str, Any]], toolchain: Mapping[str, Any],
) -> None:
    git = next(item for item in toolchain["tools"] if item["role"] == "git")
    for capture in captures:
        executable = capture.get("git_executable")
        _require(type(executable) is dict and
                 type(executable.get("source")) is dict,
                 "Git capture omits its executable identity")
        source = executable["source"]
        _require(source.get("path") == git["resolved_path"] and
                 source.get("size") == git["size"] and
                 source.get("mode") == git["mode"] and
                 source.get("sha256") == git["sha256"],
                 "Git capture and toolchain name different executables")


def _source_stage_log_bytes(
    status: str,
    commands: Sequence[Mapping[str, Any]],
    *,
    adapter_commit: str | None,
    adapter_tree: str | None,
    retained_byte_paths: Sequence[str],
    retained_path_sources: Sequence[str],
    error: BaseException | None,
) -> bytes:
    _require(status in ("complete", "failed"),
             "source acquisition log status is invalid")
    value = canonical_json_bytes({
        "schema": "leopard2-exact-main-source-acquisition-log/v1",
        "status": status,
        "command_count": len(commands),
        "commands": copy.deepcopy(list(commands)),
        "adapter_commit": adapter_commit,
        "adapter_tree": adapter_tree,
        "retained_byte_paths": sorted(retained_byte_paths),
        "retained_path_sources": sorted(retained_path_sources),
        "error": None if error is None else _safe_failure_message(error),
    })
    _require(0 < len(value) <= MAX_SOURCE_LOG_BYTES,
             "source acquisition log exceeds its byte bound")
    return value


def acquire_source_stage(
    environment: HostEnvironment,
    plan: LanePlan,
    prepared: PreparedAcquisitionRoots,
) -> SourceStageResult:
    """Run the complete no-build source/tool acquisition stage."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment) and
             isinstance(prepared, PreparedAcquisitionRoots) and
             prepared.plan == canonical_plan,
             "source acquisition dependencies are invalid")
    prepared.validate_current()
    child_environment = frozen_child_environment()
    log: list[dict[str, Any]] = []
    adapter_commit: str | None = None
    adapter_tree: str | None = None
    retained_bytes: dict[str, bytes] = {}
    retained_paths: dict[str, str] = {}
    try:
        adapter_commit, adapter_tree, adapter_submodule = \
            _repository_source_identity(
                environment, canonical_plan, child_environment, log)
        submodule_repository = canonical_plan.repository + "/sse2neon"
        _require(Path(submodule_repository).resolve(strict=True) ==
                 Path(submodule_repository) and
                 Path(submodule_repository).is_dir(),
                 "controller sse2neon source is not canonical")
        for destination, commit, tree, submodule in (
            (canonical_plan.canonical_adapter_root, adapter_commit,
             adapter_tree, adapter_submodule),
            (canonical_plan.canonical_baseline_root,
             record_contract.BASELINE_COMMIT, record_contract.BASELINE_TREE,
             record_contract.BASELINE_SSE2NEON_COMMIT),
            (canonical_plan.variant_adapter_root, adapter_commit,
             adapter_tree, adapter_submodule),
            (canonical_plan.variant_baseline_root,
             record_contract.BASELINE_COMMIT, record_contract.BASELINE_TREE,
             record_contract.BASELINE_SSE2NEON_COMMIT),
        ):
            stage_detached_source(
                environment,
                source_repository=canonical_plan.repository,
                submodule_repository=submodule_repository,
                destination=destination, commit=commit, tree=tree,
                submodule_commit=submodule,
                child_environment=child_environment, log=log)
            prepared.validate_current()
        source, adapter, retained_bytes, retained_paths, captures = \
            capture_source_evidence(
                environment, canonical_plan, adapter_commit=adapter_commit,
                adapter_tree=adapter_tree, log=log,
                child_environment=child_environment)
        toolchain, version_bytes = resolve_toolchain(
            environment, canonical_plan,
            child_environment=child_environment, log=log)
        _capture_git_tool_binding(captures, toolchain)
        retained_bytes.update(version_bytes)
        log_bytes = _source_stage_log_bytes(
            "complete", log, adapter_commit=adapter_commit,
            adapter_tree=adapter_tree,
            retained_byte_paths=tuple(retained_bytes),
            retained_path_sources=tuple(retained_paths), error=None)
        prepared.validate_current()
        return SourceStageResult(
            source=copy.deepcopy(source), adapter=copy.deepcopy(adapter),
            toolchain=copy.deepcopy(toolchain),
            retained_bytes=copy.deepcopy(retained_bytes),
            retained_paths=copy.deepcopy(retained_paths), log=log_bytes)
    except SourceStageError:
        raise
    except (ExactMainBaselineError, OSError) as error:
        log_bytes = _source_stage_log_bytes(
            "failed", log, adapter_commit=adapter_commit,
            adapter_tree=adapter_tree,
            retained_byte_paths=tuple(retained_bytes),
            retained_path_sources=tuple(retained_paths), error=error)
        raise SourceStageError(error, log_bytes) from error


def _safe_failure_message(error: BaseException) -> str:
    raw = str(error) or type(error).__name__
    flattened = " ".join(raw.replace("\r", " ").replace("\n", " ").split())
    safe = "".join(
        character if 0x20 <= ord(character) <= 0x7E else "?"
        for character in flattened)
    return (safe or "source acquisition failed")[:4096]


def seal_stage_failure(
    environment: HostEnvironment,
    plan: LanePlan,
    error: BaseException,
    *,
    stage: str,
    stage_logs: Mapping[str, bytes],
    retained_bytes: Mapping[str, bytes] | None = None,
    retained_paths: Mapping[str, str] | None = None,
    diagnostics: Mapping[str, bytes] | None = None,
    lane_root: str | None = None,
    authority_record_sha256: str | None = None,
) -> dict[str, Any]:
    """Seal one immutable acquisition or independent-verification failure."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment) and
             isinstance(error, BaseException),
             "stage failure inputs are invalid")
    _require(type(stage) is str and
             stage in record_contract.FAILURE_STAGES and
             ((stage == "independent_verification") ==
              (authority_record_sha256 is not None)),
             "failure stage and authority binding are inconsistent")
    target_root = canonical_plan.lane_root if lane_root is None else \
        _portable_absolute_path(lane_root, "failure lane root")
    _require(target_root in (
                 canonical_plan.lane_root, canonical_plan.failure_lane_root),
             "failure lane root is outside the acquisition plan")
    if isinstance(error, SourceStageError):
        _require(stage == "source_acquisition",
                 "source-stage error was assigned to another stage")
    if isinstance(error, BuildStageError):
        _require(error.stage == stage,
                 "build-stage error was assigned to another stage")
    stage_index = record_contract.STAGE_NAMES.index(stage)
    expected_stage_names = record_contract.STAGE_NAMES[:stage_index + 1]
    _require(type(stage_logs) is dict and
             list(stage_logs) == list(expected_stage_names),
             "acquisition failure stage-log prefix changed")
    message = _safe_failure_message(error)
    retained: dict[str, bytes] = {}
    path_sources = {} if retained_paths is None else dict(retained_paths)
    byte_sources = {} if retained_bytes is None else dict(retained_bytes)
    _require(type(byte_sources) is dict and type(path_sources) is dict,
             "acquisition failure retained sources are invalid")
    _require(set(byte_sources).isdisjoint(path_sources),
             "acquisition failure retained sources overlap")
    for path in sorted(byte_sources):
        relative = _safe_relative_path(path, "failure retained byte path")
        content = byte_sources[path]
        _require(not relative.startswith("logs/") and
                 type(content) is bytes and
                 len(content) <= MAX_SEALED_FILE_BYTES,
                 f"failure retained byte {relative!r} is invalid")
        retained[relative] = content
    canonical_paths: dict[str, str] = {}
    for path in sorted(path_sources):
        relative = _safe_relative_path(path, "failure retained path source")
        source = path_sources[path]
        _require(not relative.startswith("logs/") and type(source) is str,
                 f"failure retained path {relative!r} is invalid")
        canonical_paths[relative] = source
    diagnostic_values = {} if diagnostics is None else diagnostics
    _require(type(diagnostic_values) is dict,
             "failure diagnostics are not a mapping")
    for path in sorted(diagnostic_values):
        relative = _safe_relative_path(path, "failure diagnostic path")
        content = diagnostic_values[path]
        _require(relative.startswith("diagnostics/") and
                 type(content) is bytes and content and
                 len(content) <= MAX_SEALED_FILE_BYTES and
                 relative not in retained and relative not in canonical_paths,
                 f"failure diagnostic {relative!r} is invalid")
        retained[relative] = content

    stages: list[dict[str, Any]] = []
    for index, name in enumerate(expected_stage_names):
        content = stage_logs[name]
        _require(type(content) is bytes and
                 0 < len(content) <= MAX_SOURCE_LOG_BYTES,
                 f"{name} stage failure log is invalid")
        path = f"logs/{index:02d}-{name}.log"
        _require(path not in retained and path not in canonical_paths,
                 "failure stage log collides with retained evidence")
        retained[path] = content
        stages.append({
            "name": name,
            "status": "failed" if index == stage_index else "complete",
            "log": _bytes_identity(path, content),
        })

    claims = [_bytes_identity(path, retained[path])
              for path in sorted(retained)
              if not path.startswith("logs/")]
    for path in sorted(canonical_paths):
        identity = _owned_file_identity(
            canonical_paths[path], f"failure retained path {path!r}")
        claims.append({
            "relative_path": path,
            "size": identity["size"],
            "sha256": identity["sha256"],
        })
    claims.sort(key=lambda item: item["relative_path"])
    _require(len(claims) <= record_contract.MAX_RETAINED_FILES,
             "failure retained file inventory exceeds its bound")
    lane = {
        "root": target_root,
        "attempt": canonical_plan.attempt,
        "attempt_budget": 3,
        "record_relative_path": "FAILED.json",
        "seal_protocol": record_contract.SEAL_PROTOCOL,
        "stages": stages,
    }
    is_verification_error = isinstance(error, VerificationError)
    is_command_error = isinstance(error, CommandExecutionError) or (
        isinstance(error, (SourceStageError, BuildStageError)) and
        error.error_kind == "command_error")
    exit_status = error.exit_status if (
        is_command_error or is_verification_error) else 1
    error_record = {
        "kind": ("verification_error" if is_verification_error else
                 "command_error" if is_command_error else
                 "acquisition_error"),
        "message": message,
        "exit_status": _normalized_exit_status(exit_status),
    }
    if authority_record_sha256 is None:
        record = record_contract.baseline_acquisition_failure_record(
            created_utc=environment.now_utc(), lane=lane,
            stage=stage, error=error_record, retained_files=claims)
    else:
        record = record_contract.baseline_verification_failure_record(
            created_utc=environment.now_utc(), lane=lane,
            stage=stage, error=error_record, retained_files=claims,
            authority_record_sha256=authority_record_sha256)
    with LaneWriter(target_root) as writer:
        return writer.seal_record(
            record, retained, retained_paths=canonical_paths)


def seal_verification_failure(
    environment: HostEnvironment,
    plan: LanePlan,
    error: BaseException,
    *,
    stage_logs: Mapping[str, bytes],
    authority_record_sha256: str,
    retained_bytes: Mapping[str, bytes] | None = None,
    diagnostics: Mapping[str, bytes] | None = None,
) -> dict[str, Any]:
    """Seal a stage-4 failure beside an already claimed authority lane."""
    return seal_stage_failure(
        environment, plan, error, stage="independent_verification",
        stage_logs=stage_logs, retained_bytes=retained_bytes,
        diagnostics=diagnostics, lane_root=plan.failure_lane_root,
        authority_record_sha256=authority_record_sha256)


def seal_source_acquisition_failure(
    environment: HostEnvironment,
    plan: LanePlan,
    error: BaseException,
    *,
    log: bytes | None = None,
    diagnostics: Mapping[str, bytes] | None = None,
) -> dict[str, Any]:
    """Seal one immutable stage-0 failure that the offline verifier accepts."""
    inherited_log = error.log if isinstance(error, SourceStageError) else None
    stage_bytes = log if log is not None else inherited_log
    if stage_bytes is None:
        stage_bytes = _source_stage_log_bytes(
            "failed", (), adapter_commit=None, adapter_tree=None,
            retained_byte_paths=(), retained_path_sources=(), error=error)
    return seal_stage_failure(
        environment, plan, error, stage="source_acquisition",
        stage_logs={"source_acquisition": stage_bytes},
        diagnostics=diagnostics)


def build_closure_document(
    role: str,
    build_root: str,
    files: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Construct the exact build-closure document replayed by the verifier."""
    _require(role in record_contract.BUILD_ROLES and type(role) is str,
             "build closure role is invalid")
    root = _portable_absolute_path(build_root, "build closure root")
    _require(type(files) in (list, tuple) and 0 < len(files) <=
             record_contract.MAX_CLOSURE_FILES,
             "build closure file inventory is invalid")
    canonical: list[dict[str, Any]] = []
    total = 0
    for index, raw in enumerate(files):
        _require(type(raw) is dict and
                 set(raw) == {"relative_path", "size", "sha256"},
                 f"build closure file {index} has an unexpected key set")
        relative = _safe_relative_path(
            raw["relative_path"], f"build closure file {index} path")
        size = raw["size"]
        digest = raw["sha256"]
        _require(type(size) is int and 0 <= size <= MAX_SEALED_FILE_BYTES,
                 f"build closure file {index} size is invalid")
        _require(type(digest) is str and len(digest) == 64 and
                 all(character in "0123456789abcdef" for character in digest),
                 f"build closure file {index} digest is invalid")
        total += size
        _require(total <= MAX_SEALED_TOTAL_BYTES,
                 "build closure total bytes exceed the producer bound")
        canonical.append({
            "relative_path": relative,
            "size": size,
            "sha256": digest,
        })
    _require([item["relative_path"] for item in canonical] ==
             sorted({item["relative_path"] for item in canonical}),
             "build closure paths are not sorted and unique")
    return {
        "schema": BUILD_CLOSURE_SCHEMA,
        "role": role,
        "build_root": root,
        "files": copy.deepcopy(canonical),
        "file_count": len(canonical),
    }


def build_role_roots(plan: LanePlan, role: str) -> dict[str, str]:
    """Resolve one oriented build role to its frozen source/build roots."""
    canonical_plan = validate_lane_plan(plan)
    _require(role in record_contract.BUILD_ROLES and type(role) is str,
             "exact-main build role is invalid")
    if role == "path_variant":
        roots = {
            "adapter_source_root": canonical_plan.variant_adapter_root,
            "baseline_source_root": canonical_plan.variant_baseline_root,
            "build_root": canonical_plan.variant_build_root,
        }
    else:
        roots = {
            "adapter_source_root": canonical_plan.canonical_adapter_root,
            "baseline_source_root": canonical_plan.canonical_baseline_root,
            "build_root": canonical_plan.canonical_build_root,
        }
    record_contract.exact_main_build_cache_requirements(roots)
    return roots


def _tool_paths_from_toolchain(toolchain: Mapping[str, Any]) -> dict[str, str]:
    _require(type(toolchain) is dict and type(toolchain.get("tools")) is list,
             "exact-main toolchain is missing its tool inventory")
    result: dict[str, str] = {}
    for index, item in enumerate(toolchain["tools"]):
        _require(type(item) is dict and type(item.get("role")) is str and
                 type(item.get("resolved_path")) is str,
                 f"exact-main tool {index} is invalid")
        role = item["role"]
        _require(role in record_contract.TOOL_ROLES and role not in result,
                 "exact-main tool inventory is not unique")
        result[role] = _portable_absolute_path(
            item["resolved_path"], f"exact-main {role} resolved path")
    _require(set(result) == set(record_contract.TOOL_ROLES),
             "exact-main tool inventory changed")
    return result


def parse_ctest_success_summary(content: bytes) -> tuple[int, int]:
    """Validate only CTest's fixed pass/fail summary; ignore timing lines."""
    _require(type(content) is bytes and
             0 < len(content) <= MAX_COMMAND_OUTPUT_BYTES,
             "CTest stdout is not bounded non-empty bytes")
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise AcquisitionError("CTest stdout is not strict UTF-8") from error
    _require(text.endswith("\n") and "\r" not in text and "\0" not in text,
             "CTest stdout is not canonical LF text")
    summaries = [
        line for line in text[:-1].split("\n")
        if re.fullmatch(
            r"[0-9]+% tests passed, [0-9]+ tests failed out of [0-9]+",
            line)
    ]
    _require(summaries == [record_contract.CTEST_SUMMARY_LINE],
             "CTest stdout lacks the exact one-test success summary")
    return (1, 0)


def _identity_census_is_zero(identity: Mapping[str, Any]) -> bool:
    return all(
        root["occurrences"] == 0 and all(
            row["occurrences"] == 0 for row in root["sections"])
        for root in identity["path_string_census"]["roots"]
    )


def acquire_build_stage(
    environment: HostEnvironment,
    plan: LanePlan,
    prepared: PreparedAcquisitionRoots,
    *,
    role: str,
    toolchain: Mapping[str, Any],
    controller_sha256: str,
) -> BuildStageResult:
    """Acquire one complete build/runtime/correctness evidence role."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment) and
             isinstance(prepared, PreparedAcquisitionRoots) and
             prepared.plan == canonical_plan and
             role in record_contract.BUILD_ROLES and type(role) is str,
             "build-stage dependencies are invalid")
    _require(type(controller_sha256) is str and
             len(controller_sha256) == 64 and all(
                 character in "0123456789abcdef"
                 for character in controller_sha256),
             "build-stage controller SHA-256 is invalid")
    prepared.validate_current()
    roots = build_role_roots(canonical_plan, role)
    tools = _tool_paths_from_toolchain(toolchain)
    role_path = role.replace("_", "-")
    build_field = "variant_build_root" if role == "path_variant" else \
        "canonical_build_root"
    _require(os.listdir(prepared.descriptors[build_field]) == [],
             f"{role} build root is not empty")
    child_environment = frozen_child_environment()
    commands: list[dict[str, Any]] = []
    retained_bytes: dict[str, bytes] = {}
    retained_paths: dict[str, str] = {}
    try:
        configure_argv = record_contract.exact_main_configure_argv(
            cmake=tools["cmake"], compiler=tools["compiler"], roots=roots)
        configure_result = _checked_run(
            environment, configure_argv, cwd=canonical_plan.scratch_root,
            child_environment=child_environment, log=commands,
            label=f"{role} CMake configure", timeout=BUILD_TIMEOUT_SECONDS)
        retained_bytes[f"builds/{role_path}/configure.log"] = \
            _command_transcript_bytes(
                configure_argv, canonical_plan.scratch_root, configure_result)
        prepared.validate_current()

        build_argv = record_contract.exact_main_build_argv(
            cmake=tools["cmake"], roots=roots)
        build_result = _checked_run(
            environment, build_argv, cwd=canonical_plan.scratch_root,
            child_environment=child_environment, log=commands,
            label=f"{role} serial build", timeout=BUILD_TIMEOUT_SECONDS)
        retained_bytes[f"builds/{role_path}/build.log"] = \
            _command_transcript_bytes(
                build_argv, canonical_plan.scratch_root, build_result)
        prepared.validate_current()

        cache_path = roots["build_root"] + "/CMakeCache.txt"
        cache_bytes = _read_host_file(
            cache_path, f"{role} CMake cache",
            maximum_bytes=MAX_COMMAND_OUTPUT_BYTES)
        record_contract.validate_cmake_cache(cache_bytes, roots)
        retained_bytes[f"builds/{role_path}/CMakeCache.txt"] = cache_bytes
        prepared.validate_current()

        compile_path = roots["build_root"] + "/compile_commands.json"
        compile_bytes = _read_host_file(
            compile_path, f"{role} compile commands",
            maximum_bytes=MAX_COMMAND_OUTPUT_BYTES)
        compile_value = identity_contract.strict_json_loads(
            compile_bytes, f"{role} compile commands JSON")
        record_contract.validate_compile_commands(
            compile_value, roots=roots, compiler=tools["compiler"],
            profile=record_contract.exact_main_build_profile())
        retained_bytes[f"builds/{role_path}/compile_commands.json"] = \
            compile_bytes
        prepared.validate_current()

        closure, closure_bytes = build_root_census(roots["build_root"], role)
        closure_relative = f"builds/{role_path}/build-closure.json"
        retained_bytes[closure_relative] = closure_bytes
        prepared.validate_current()

        executable_source = roots["build_root"] + "/leopard_main_benchmark"
        archive_source = roots["build_root"] + "/libleopard_main_exact.a"
        executable_stage = canonical_plan.scratch_root + \
            f"/{role_path}-leopard_main_benchmark"
        archive_stage = canonical_plan.scratch_root + \
            f"/{role_path}-libleopard_main_exact.a"
        executable_identity = stage_build_output(
            executable_source, executable_stage, f"{role} executable",
            maximum_bytes=identity_contract.MAX_ELF_INPUT_BYTES)
        archive_identity = stage_build_output(
            archive_source, archive_stage, f"{role} archive")
        executable_relative = \
            f"artifacts/{role_path}/leopard_main_benchmark"
        archive_relative = f"artifacts/{role_path}/libleopard_main_exact.a"
        retained_paths[executable_relative] = executable_stage
        retained_paths[archive_relative] = archive_stage
        prepared.validate_current()

        ldd_argv = [tools["ldd"], executable_source]
        ldd_result = _checked_run(
            environment, ldd_argv, cwd=canonical_plan.scratch_root,
            child_environment=child_environment, log=commands,
            label=f"{role} runtime dependency discovery",
            timeout=ATTESTATION_TIMEOUT_SECONDS)
        _require(ldd_result.stderr == b"",
                 f"{role} runtime dependency discovery emitted diagnostics")
        ldd_rows = normalize_ldd_output(ldd_result.stdout)
        ldd_bytes = canonical_ldd_text(ldd_rows)
        ldd_relative = f"runtime/{role_path}/ldd.txt"
        retained_bytes[ldd_relative] = ldd_bytes
        dependencies: list[dict[str, Any]] = []
        for row in ldd_rows:
            if row["kind"] == "virtual":
                dependencies.append({
                    "soname": row["soname"], "kind": "virtual",
                    "path": None, "size": None, "sha256": None,
                })
            else:
                dependency = _runtime_dependency_identity(
                    row["path"], f"{role} runtime dependency {row['soname']}")
                dependencies.append({
                    "soname": row["soname"], "kind": "file",
                    "path": dependency["path"],
                    "size": dependency["size"],
                    "sha256": dependency["sha256"],
                })
        prepared.validate_current()

        benchmark_argv = record_contract.exact_main_benchmark_argv(
            executable_path=executable_source)
        benchmark_result = _checked_run(
            environment, benchmark_argv, cwd=canonical_plan.scratch_root,
            child_environment=child_environment, log=commands,
            label=f"{role} benchmark correctness attestation",
            timeout=ATTESTATION_TIMEOUT_SECONDS)
        benchmark_value = identity_contract.strict_json_loads(
            benchmark_result.stdout, f"{role} benchmark attestation JSON")
        record_contract.validate_attestation_stdout(
            benchmark_value, argv=benchmark_argv,
            reported_schema=record_contract.BENCHMARK_SCHEMA)
        benchmark_stdout_relative = \
            f"attestations/{role_path}/benchmark.stdout.json"
        benchmark_stderr_relative = \
            f"attestations/{role_path}/benchmark.stderr"
        retained_bytes[benchmark_stdout_relative] = benchmark_result.stdout
        retained_bytes[benchmark_stderr_relative] = benchmark_result.stderr
        prepared.validate_current()

        ctest_argv = record_contract.exact_main_ctest_argv(
            ctest=tools["ctest"], build_root=roots["build_root"])
        ctest_result = _checked_run(
            environment, ctest_argv, cwd=canonical_plan.scratch_root,
            child_environment=child_environment, log=commands,
            label=f"{role} CTest correctness attestation",
            timeout=ATTESTATION_TIMEOUT_SECONDS)
        passed, failed = parse_ctest_success_summary(ctest_result.stdout)
        ctest_stdout_relative = f"attestations/{role_path}/ctest.stdout.log"
        ctest_stderr_relative = f"attestations/{role_path}/ctest.stderr.log"
        retained_bytes[ctest_stdout_relative] = ctest_result.stdout
        retained_bytes[ctest_stderr_relative] = ctest_result.stderr
        prepared.validate_current()

        attested_executable = _host_file_identity(
            executable_source, f"{role} attested executable",
            maximum_bytes=identity_contract.MAX_ELF_INPUT_BYTES)
        _require(attested_executable["size"] == executable_identity["size"] and
                 attested_executable["sha256"] ==
                 executable_identity["sha256"],
                 f"{role} attested executable changed during attestation")

        executable_bytes = _read_owned_file(
            executable_stage, f"{role} staged executable",
            maximum_bytes=identity_contract.MAX_ELF_INPUT_BYTES)
        identity = identity_contract.normalized_code_identity_from_elf_bytes(
            executable_bytes, roots=roots)
        _require(identity["artifact"]["size"] == executable_identity["size"] and
                 identity["artifact"]["sha256"] ==
                 executable_identity["sha256"],
                 f"{role} normalized identity changed its artifact")
        _require(_identity_census_is_zero(identity),
                 f"{role} selected ELF sections retain an acquisition root")

        build = {
            "role": role,
            "roots": copy.deepcopy(roots),
            "configure_argv": configure_argv,
            "build_argv": build_argv,
            "configure_log": _bytes_identity(
                f"builds/{role_path}/configure.log",
                retained_bytes[f"builds/{role_path}/configure.log"]),
            "build_log": _bytes_identity(
                f"builds/{role_path}/build.log",
                retained_bytes[f"builds/{role_path}/build.log"]),
            "cmake_cache": _bytes_identity(
                f"builds/{role_path}/CMakeCache.txt", cache_bytes),
            "compile_commands": _bytes_identity(
                f"builds/{role_path}/compile_commands.json", compile_bytes),
            "executable": {
                "name": "leopard_main_benchmark",
                "build_relative_path": "leopard_main_benchmark",
                "retained_relative_path": executable_relative,
                "size": executable_identity["size"],
                "sha256": executable_identity["sha256"],
            },
            "archive": {
                "name": "libleopard_main_exact.a",
                "build_relative_path": "libleopard_main_exact.a",
                "retained_relative_path": archive_relative,
                "size": archive_identity["size"],
                "sha256": archive_identity["sha256"],
            },
            "closure": {
                **_bytes_identity(closure_relative, closure_bytes),
                "file_count": closure["file_count"],
            },
        }
        build = record_contract.validate_exact_main_build(
            build, role=role, tools=tools)
        record_contract.validate_build_closure(
            closure, role=role, build=build)

        runtime = {
            "role": role,
            "executable_sha256": executable_identity["sha256"],
            "canonical_ldd_output": _bytes_identity(
                ldd_relative, ldd_bytes),
            "dependencies": dependencies,
        }
        attestation = {
            "role": role,
            "argv": benchmark_argv,
            "stdout": _bytes_identity(
                benchmark_stdout_relative, benchmark_result.stdout),
            "stderr": _bytes_identity(
                benchmark_stderr_relative, benchmark_result.stderr),
            "exit_status": 0,
            "reported_schema": record_contract.BENCHMARK_SCHEMA,
            "main_source_commit": record_contract.BASELINE_COMMIT,
            "pure_avx2": True,
            "round_trip": True,
            "ctest": {
                "argv": ctest_argv,
                "stdout": _bytes_identity(
                    ctest_stdout_relative, ctest_result.stdout),
                "stderr": _bytes_identity(
                    ctest_stderr_relative, ctest_result.stderr),
                "exit_status": 0,
                "passed": passed,
                "failed": failed,
            },
        }
        log_bytes = _build_stage_log_bytes(
            role, "complete", commands,
            controller_sha256=controller_sha256, error=None)
        prepared.validate_current()
        return BuildStageResult(
            role=role, build=copy.deepcopy(build),
            runtime=copy.deepcopy(runtime),
            attestation=copy.deepcopy(attestation),
            identity=copy.deepcopy(identity),
            retained_bytes=copy.deepcopy(retained_bytes),
            retained_paths=copy.deepcopy(retained_paths), log=log_bytes)
    except BuildStageError:
        raise
    except (ExactMainBaselineError, OSError) as error:
        log_bytes = _build_stage_log_bytes(
            role, "failed", commands,
            controller_sha256=controller_sha256, error=error)
        raise BuildStageError(role, error, log_bytes) from error


def expected_sha256sums(digests: Mapping[str, str]) -> bytes:
    """Return the exact outer checksum ledger bytes."""
    _require(type(digests) is dict and _SHA256SUMS_PATH not in digests,
             "checksum input is not an exact path-to-digest mapping")
    rows: list[bytes] = []
    for path in sorted(digests):
        relative = _safe_relative_path(path, "checksum path")
        digest = digests[path]
        _require(type(digest) is str and len(digest) == 64 and
                 all(character in "0123456789abcdef" for character in digest),
                 f"checksum digest for {relative!r} is invalid")
        rows.append(f"{digest}  ./{relative}\n".encode("ascii"))
    return b"".join(rows)


def _derived_directories(paths: Sequence[str]) -> set[str]:
    directories = {"."}
    for path in paths:
        components = path.split("/")[:-1]
        for index in range(1, len(components) + 1):
            directories.add("/".join(components[:index]))
    return directories


def _metadata_entry(path: str, status: os.stat_result,
                    node_type: str) -> dict[str, Any]:
    return {
        "gid": status.st_gid,
        "mode": f"{stat.S_IMODE(status.st_mode):04o}",
        "nlink": status.st_nlink,
        "path": path,
        "type": node_type,
        "uid": status.st_uid,
    }


def expected_tree_metadata(
    nodes: Mapping[str, tuple[str, os.stat_result]],
) -> dict[str, Any]:
    """Return the exact metadata document for one final-mode tree snapshot."""
    _require(type(nodes) is dict and "." in nodes,
             "tree metadata snapshot omits its root")
    root_type, root_status = nodes["."]
    _require(root_type == "directory", "tree metadata root is not a directory")
    entries = [
        _metadata_entry(path, nodes[path][1], nodes[path][0])
        for path in sorted(nodes)
        if path != _TREE_METADATA_PATH
    ]
    return {
        "entries": entries,
        "excluded_paths": [_TREE_METADATA_PATH],
        "final_mode_policy": "observed mode with all write bits removed",
        "root": ".",
        "schema": TREE_METADATA_SCHEMA,
        "self_policy": {
            "gid": root_status.st_gid,
            "mode": "0400",
            "nlink": 1,
            "sha256_binding":
                "exactly one ./TREE-METADATA.json checksum entry",
            "type": "file",
            "uid": root_status.st_uid,
        },
        "uid_gid_policy": {
            "gid": root_status.st_gid,
            "rule": "every retained node has the invoking effective uid and gid",
            "uid": root_status.st_uid,
        },
    }


class LaneWriter:
    """Fd-anchored exclusive publisher for one never-reused evidence lane."""

    def __init__(self, root: str):
        self.root_path = _portable_absolute_path(root, "lane root")
        self._root_descriptor = -1
        self._closed = False
        self._sealed = False
        self._published: dict[str, dict[str, Any]] = {}
        self._directories: set[str] = {"."}
        try:
            os.mkdir(self.root_path, WRITABLE_DIRECTORY_MODE)
            self._root_descriptor = os.open(
                self.root_path,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
            )
            self._validate_directory_fd(
                self._root_descriptor, ".", WRITABLE_DIRECTORY_MODE)
        except FileExistsError as error:
            raise AcquisitionError("lane root already exists and cannot be reused") \
                from error
        except BaseException:
            if self._root_descriptor >= 0:
                os.close(self._root_descriptor)
                self._root_descriptor = -1
            raise

    def __enter__(self) -> "LaneWriter":
        return self

    def __exit__(self, _kind: object, _value: object,
                 _traceback: object) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        if self._root_descriptor >= 0:
            os.close(self._root_descriptor)
        self._root_descriptor = -1
        self._closed = True

    def _require_open(self) -> None:
        _require(not self._closed and self._root_descriptor >= 0,
                 "lane writer is closed")
        _require(not self._sealed, "sealed lane cannot be changed")

    def _validate_directory_fd(self, descriptor: int, relative: str,
                               expected_mode: int) -> os.stat_result:
        status = os.fstat(descriptor)
        _require(stat.S_ISDIR(status.st_mode),
                 f"lane directory {relative!r} is not a directory")
        _require(status.st_dev == os.fstat(self._root_descriptor).st_dev,
                 f"lane directory {relative!r} crosses a device boundary")
        _require(status.st_uid == os.geteuid() and
                 status.st_gid == os.getegid(),
                 f"lane directory {relative!r} has another owner")
        _require(stat.S_IMODE(status.st_mode) == expected_mode,
                 f"lane directory {relative!r} has an unsafe mode")
        return status

    def _open_directory(
        self,
        relative: str,
        *,
        expected_mode: int = WRITABLE_DIRECTORY_MODE,
    ) -> int:
        if relative == ".":
            descriptor = os.dup(self._root_descriptor)
            try:
                self._validate_directory_fd(descriptor, ".", expected_mode)
                return descriptor
            except BaseException:
                os.close(descriptor)
                raise
        descriptor = os.dup(self._root_descriptor)
        prefix: list[str] = []
        try:
            for component in relative.split("/"):
                prefix.append(component)
                child = os.open(
                    component,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                    os.O_CLOEXEC,
                    dir_fd=descriptor,
                )
                os.close(descriptor)
                descriptor = child
                self._validate_directory_fd(
                    descriptor, "/".join(prefix), expected_mode)
            return descriptor
        except BaseException:
            os.close(descriptor)
            raise

    def _ensure_parent(self, relative_path: str) -> tuple[int, str]:
        components = relative_path.split("/")
        basename = components[-1]
        descriptor = os.dup(self._root_descriptor)
        prefix: list[str] = []
        try:
            for component in components[:-1]:
                prefix.append(component)
                relative = "/".join(prefix)
                if relative not in self._directories:
                    try:
                        os.mkdir(component, WRITABLE_DIRECTORY_MODE,
                                 dir_fd=descriptor)
                    except FileExistsError as error:
                        raise AcquisitionError(
                            f"lane directory {relative!r} appeared concurrently") \
                            from error
                    os.fsync(descriptor)
                    self._directories.add(relative)
                child = os.open(
                    component,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                    os.O_CLOEXEC,
                    dir_fd=descriptor,
                )
                os.close(descriptor)
                descriptor = child
                self._validate_directory_fd(
                    descriptor, relative, WRITABLE_DIRECTORY_MODE)
            return descriptor, basename
        except BaseException:
            os.close(descriptor)
            raise

    def publish_bytes(self, relative_path: str, content: bytes) -> dict[str, Any]:
        """Publish one file without ever replacing an existing destination."""
        self._require_open()
        relative = _safe_relative_path(relative_path, "lane file path")
        _require(relative not in self._published,
                 f"lane file {relative!r} was already published")
        _require(type(content) is bytes and len(content) <= MAX_SEALED_FILE_BYTES,
                 f"lane file {relative!r} is not bounded bytes")
        content_digest = _sha256(content)
        parent, basename = self._ensure_parent(relative)
        staging = ".leopard-stage-" + secrets.token_hex(16)
        descriptor = -1
        linked = False
        try:
            descriptor = os.open(
                staging,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW |
                os.O_CLOEXEC,
                WRITABLE_FILE_MODE,
                dir_fd=parent,
            )
            _write_all(descriptor, content, relative)
            os.fsync(descriptor)
            status = os.fstat(descriptor)
            _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1 and
                     status.st_uid == os.geteuid() and
                     status.st_gid == os.getegid() and
                     stat.S_IMODE(status.st_mode) == WRITABLE_FILE_MODE and
                     status.st_size == len(content),
                     f"staged lane file {relative!r} changed")
            os.close(descriptor)
            descriptor = -1
            try:
                os.link(staging, basename, src_dir_fd=parent,
                        dst_dir_fd=parent, follow_symlinks=False)
            except FileExistsError as error:
                raise AcquisitionError(
                    f"lane file {relative!r} already exists") from error
            linked = True
            os.unlink(staging, dir_fd=parent)
            linked = False
            os.fsync(parent)
            final_descriptor = os.open(
                basename, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                dir_fd=parent)
            try:
                final = os.fstat(final_descriptor)
                _require(stat.S_ISREG(final.st_mode) and final.st_nlink == 1 and
                         final.st_dev == os.fstat(self._root_descriptor).st_dev and
                         final.st_uid == os.geteuid() and
                         final.st_gid == os.getegid() and
                         stat.S_IMODE(final.st_mode) == WRITABLE_FILE_MODE and
                         final.st_size == len(content) and
                         _hash_fd(final_descriptor, final.st_size, relative) ==
                         content_digest,
                         f"published lane file {relative!r} changed")
            finally:
                os.close(final_descriptor)
        except BaseException:
            if descriptor >= 0:
                os.close(descriptor)
            try:
                os.unlink(staging, dir_fd=parent)
            except FileNotFoundError:
                pass
            if linked:
                # If unlinking the staging name failed after the final link, the
                # final inode has nlink=2 and is intentionally left fail-closed.
                pass
            raise
        finally:
            os.close(parent)
        identity = {
            "relative_path": relative,
            "size": len(content),
            "sha256": content_digest,
        }
        self._published[relative] = identity
        return copy.deepcopy(identity)

    def publish_path(
        self, relative_path: str, source_path: str,
    ) -> dict[str, Any]:
        """Stream one owner-only source file into the lane without buffering."""
        self._require_open()
        relative = _safe_relative_path(relative_path, "lane file path")
        source = _portable_absolute_path(
            source_path, f"lane source for {relative!r}")
        _require(relative not in self._published,
                 f"lane file {relative!r} was already published")
        parent, basename = self._ensure_parent(relative)
        staging = ".leopard-stage-" + secrets.token_hex(16)
        source_descriptor = -1
        destination_descriptor = -1
        linked = False
        try:
            source_descriptor = os.open(
                source, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
            source_before = os.fstat(source_descriptor)
            source_pathname = os.lstat(source)
            _require(stat.S_ISREG(source_before.st_mode) and
                     source_before.st_nlink == 1 and
                     source_before.st_uid == os.geteuid() and
                     source_before.st_gid == os.getegid() and
                     stat.S_IMODE(source_before.st_mode) ==
                     WRITABLE_FILE_MODE and
                     0 < source_before.st_size <= MAX_SEALED_FILE_BYTES and
                     (source_before.st_dev, source_before.st_ino) ==
                     (source_pathname.st_dev, source_pathname.st_ino) and
                     os.path.realpath(source) == source,
                     f"lane source for {relative!r} is not one owned file")
            destination_descriptor = os.open(
                staging,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW |
                os.O_CLOEXEC,
                WRITABLE_FILE_MODE,
                dir_fd=parent,
            )
            digest = hashlib.sha256()
            copied = 0
            while copied < source_before.st_size:
                chunk = os.read(
                    source_descriptor,
                    min(READ_CHUNK, source_before.st_size - copied))
                _require(bool(chunk),
                         f"lane source for {relative!r} ended while copied")
                _write_all(destination_descriptor, chunk, relative)
                digest.update(chunk)
                copied += len(chunk)
            _require(os.read(source_descriptor, 1) == b"",
                     f"lane source for {relative!r} grew while copied")
            source_after = os.fstat(source_descriptor)
            _require((source_after.st_dev, source_after.st_ino,
                      source_after.st_size, source_after.st_mtime_ns,
                      source_after.st_ctime_ns) ==
                     (source_before.st_dev, source_before.st_ino,
                      source_before.st_size, source_before.st_mtime_ns,
                      source_before.st_ctime_ns),
                     f"lane source for {relative!r} changed while copied")
            content_digest = digest.hexdigest()
            os.fsync(destination_descriptor)
            staged = os.fstat(destination_descriptor)
            _require(stat.S_ISREG(staged.st_mode) and staged.st_nlink == 1 and
                     staged.st_uid == os.geteuid() and
                     staged.st_gid == os.getegid() and
                     stat.S_IMODE(staged.st_mode) == WRITABLE_FILE_MODE and
                     staged.st_size == copied,
                     f"staged lane file {relative!r} changed")
            os.close(destination_descriptor)
            destination_descriptor = -1
            os.close(source_descriptor)
            source_descriptor = -1
            try:
                os.link(staging, basename, src_dir_fd=parent,
                        dst_dir_fd=parent, follow_symlinks=False)
            except FileExistsError as error:
                raise AcquisitionError(
                    f"lane file {relative!r} already exists") from error
            linked = True
            os.unlink(staging, dir_fd=parent)
            linked = False
            os.fsync(parent)
            final_descriptor = os.open(
                basename, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                dir_fd=parent)
            try:
                final = os.fstat(final_descriptor)
                _require(stat.S_ISREG(final.st_mode) and final.st_nlink == 1 and
                         final.st_dev == os.fstat(
                             self._root_descriptor).st_dev and
                         final.st_uid == os.geteuid() and
                         final.st_gid == os.getegid() and
                         stat.S_IMODE(final.st_mode) == WRITABLE_FILE_MODE and
                         final.st_size == copied and
                         _hash_fd(final_descriptor, final.st_size, relative) ==
                         content_digest,
                         f"published lane file {relative!r} changed")
            finally:
                os.close(final_descriptor)
        except BaseException as error:
            if destination_descriptor >= 0:
                os.close(destination_descriptor)
            if source_descriptor >= 0:
                os.close(source_descriptor)
            try:
                os.unlink(staging, dir_fd=parent)
            except FileNotFoundError:
                pass
            if linked:
                # Preserve a fail-closed final inode if staging unlink failed.
                pass
            if isinstance(error, OSError):
                raise AcquisitionError(
                    f"cannot publish lane source {relative!r}: {error}") \
                    from error
            raise
        finally:
            os.close(parent)
        identity = {
            "relative_path": relative,
            "size": copied,
            "sha256": content_digest,
        }
        self._published[relative] = identity
        return copy.deepcopy(identity)

    def _open_placeholder(self, relative: str) -> int:
        parent_path, basename = relative.rsplit("/", 1) \
            if "/" in relative else (".", relative)
        parent = self._open_directory(parent_path)
        descriptor = -1
        try:
            descriptor = os.open(
                basename, os.O_RDWR | os.O_NOFOLLOW | os.O_CLOEXEC,
                dir_fd=parent)
            status = os.fstat(descriptor)
            _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1 and
                     status.st_uid == os.geteuid() and
                     status.st_gid == os.getegid() and
                     stat.S_IMODE(status.st_mode) == WRITABLE_FILE_MODE and
                     status.st_size == 0,
                     f"seal placeholder {relative!r} changed")
            return descriptor
        except BaseException:
            if descriptor >= 0:
                os.close(descriptor)
            raise
        finally:
            os.close(parent)

    def _replace_placeholder(
        self,
        relative: str,
        descriptor: int,
        content: bytes,
    ) -> None:
        try:
            before = os.fstat(descriptor)
            _require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                     before.st_uid == os.geteuid() and
                     before.st_gid == os.getegid() and
                     stat.S_IMODE(before.st_mode) == LANE_FILE_MODE,
                     f"seal placeholder {relative!r} changed")
            os.ftruncate(descriptor, 0)
            os.lseek(descriptor, 0, os.SEEK_SET)
            _write_all(descriptor, content, relative)
            os.fsync(descriptor)
            after = os.fstat(descriptor)
            _require(after.st_nlink == 1 and after.st_size == len(content) and
                     stat.S_IMODE(after.st_mode) == LANE_FILE_MODE,
                     f"seal placeholder {relative!r} did not close")
        except OSError as error:
            raise AcquisitionError(
                f"seal placeholder {relative!r} could not be written") \
                from error
        self._published[relative] = {
            "relative_path": relative,
            "size": len(content),
            "sha256": _sha256(content),
        }

    def _scan_tree(self, expected_files: set[str],
                   expected_directories: set[str],
                   *,
                   hash_files: bool = True,
                   ) -> tuple[dict[str, tuple[str, os.stat_result]],
                              dict[str, str]]:
        _require(type(hash_files) is bool,
                 "lane scan hash policy is not boolean")
        root_status = self._validate_directory_fd(
            self._root_descriptor, ".", LANE_DIRECTORY_MODE)
        root_device = root_status.st_dev
        nodes: dict[str, tuple[str, os.stat_result]] = {
            ".": ("directory", root_status)}
        digests: dict[str, str] = {}
        total_bytes = 0

        def walk(directory: int, prefix: str, depth: int) -> None:
            nonlocal total_bytes
            _require(depth <= MAX_TREE_DEPTH,
                     "lane exceeds its directory-depth bound")
            entries = []
            with os.scandir(directory) as iterator:
                for entry in iterator:
                    entries.append(entry)
                    _require(len(nodes) + len(entries) <= MAX_TREE_NODES,
                             "lane contains too many nodes")
            entries.sort(key=lambda item: item.name)
            for entry in entries:
                name = _safe_relative_path(entry.name, "lane node name")
                relative = name if not prefix else prefix + "/" + name
                status = os.stat(name, dir_fd=directory, follow_symlinks=False)
                _require(status.st_dev == root_device and
                         status.st_uid == os.geteuid() and
                         status.st_gid == os.getegid(),
                         f"lane node {relative!r} has an invalid identity")
                if stat.S_ISDIR(status.st_mode):
                    _require(stat.S_IMODE(status.st_mode) ==
                             LANE_DIRECTORY_MODE,
                             f"lane directory {relative!r} is not sealed")
                    child = os.open(
                        name,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                        os.O_CLOEXEC,
                        dir_fd=directory,
                    )
                    try:
                        current = os.fstat(child)
                        _require((current.st_dev, current.st_ino) ==
                                 (status.st_dev, status.st_ino),
                                 f"lane directory {relative!r} was replaced")
                        nodes[relative] = ("directory", current)
                        walk(child, relative, depth + 1)
                    finally:
                        os.close(child)
                    continue
                _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1 and
                         stat.S_IMODE(status.st_mode) == LANE_FILE_MODE,
                         f"lane file {relative!r} is not a sealed regular file")
                _require(status.st_size <= MAX_SEALED_FILE_BYTES,
                         f"lane file {relative!r} is oversized")
                total_bytes += status.st_size
                _require(total_bytes <= MAX_SEALED_TOTAL_BYTES,
                         "lane total bytes exceed the producer bound")
                descriptor = os.open(
                    name, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                    dir_fd=directory)
                try:
                    current = os.fstat(descriptor)
                    _require((current.st_dev, current.st_ino, current.st_size) ==
                             (status.st_dev, status.st_ino, status.st_size),
                             f"lane file {relative!r} was replaced")
                    if hash_files:
                        digests[relative] = _hash_fd(
                            descriptor, current.st_size, relative)
                    nodes[relative] = ("file", current)
                finally:
                    os.close(descriptor)

        walk(self._root_descriptor, "", 1)
        observed_files = {
            path for path, (kind, _status) in nodes.items() if kind == "file"}
        observed_directories = {
            path for path, (kind, _status) in nodes.items()
            if kind == "directory"}
        _require(observed_files == expected_files,
                 "lane file set differs from the terminal inventory")
        _require(observed_directories == expected_directories,
                 "lane directory set differs from the terminal inventory")
        return nodes, digests

    def _apply_final_modes(self, files: set[str], directories: set[str]) -> None:
        for relative in sorted(files):
            parent_path, basename = relative.rsplit("/", 1) \
                if "/" in relative else (".", relative)
            parent = self._open_directory(parent_path)
            descriptor = -1
            try:
                descriptor = os.open(
                    basename, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                    dir_fd=parent)
                status = os.fstat(descriptor)
                _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1,
                         f"lane file {relative!r} changed before sealing")
                os.fchmod(descriptor, LANE_FILE_MODE)
                os.fsync(descriptor)
            finally:
                if descriptor >= 0:
                    os.close(descriptor)
                os.close(parent)
        for relative in sorted(
                (path for path in directories if path != "."),
                key=lambda path: (-path.count("/"), path),
        ):
            descriptor = self._open_directory(relative)
            try:
                os.fchmod(descriptor, LANE_DIRECTORY_MODE)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        os.fchmod(self._root_descriptor, LANE_DIRECTORY_MODE)
        os.fsync(self._root_descriptor)

    def seal_record(
        self,
        record: Mapping[str, Any],
        retained_files: Mapping[str, bytes],
        *,
        retained_paths: Mapping[str, str] | None = None,
    ) -> dict[str, Any]:
        """Publish and seal bytes/path sources plus their exact terminal."""
        self._require_open()
        _require(type(record) is dict, "lane terminal is not a JSON object")
        schema = record.get("schema")
        try:
            if schema == record_contract.AUTHORITY_SCHEMA:
                validated = record_contract.validate_baseline_authority_record(
                    record)
                inventory = record_contract.authority_retained_inventory(
                    validated)
                terminal = "baseline-authority.json"
            else:
                validated = record_contract.validate_baseline_failure_record(
                    record)
                inventory = record_contract.failure_retained_inventory(validated)
                terminal = "FAILED.json"
        except ExactMainBaselineError as error:
            raise AcquisitionError("lane terminal record is invalid") from error
        _require(validated["lane"]["root"] == self.root_path,
                 "lane terminal names another root")
        _require(validated["lane"]["record_relative_path"] == terminal,
                 "lane terminal path changed")
        _require(type(retained_files) is dict,
                 "retained lane files are not a path-to-bytes mapping")
        path_files = {} if retained_paths is None else retained_paths
        _require(type(path_files) is dict and
                 set(retained_files).isdisjoint(path_files),
                 "retained path-backed lane files are invalid")
        expected = {item["relative_path"]: item for item in inventory}
        _require(set(expected) ==
                 set(retained_files) | set(path_files) | {terminal},
                 "retained lane files differ from the terminal inventory")
        for path, content in retained_files.items():
            relative = _safe_relative_path(path, "retained lane path")
            _require(relative not in _TERMINAL_PATHS,
                     "retained lane mapping supplied a terminal")
            _require(type(content) is bytes,
                     f"retained lane file {relative!r} is not bytes")
            claim = expected[relative]
            _require(claim["size"] == len(content) and
                     claim["sha256"] == _sha256(content),
                     f"retained lane file {relative!r} differs from its claim")
        for path, source in path_files.items():
            relative = _safe_relative_path(path, "retained lane path")
            _require(relative not in _TERMINAL_PATHS and type(source) is str,
                     "retained path-backed lane file is invalid")
            identity = _owned_file_identity(
                source, f"retained path-backed lane file {relative!r}")
            claim = expected[relative]
            _require(claim["size"] == identity["size"] and
                     claim["sha256"] == identity["sha256"],
                     f"retained lane file {relative!r} differs from its claim")
        terminal_content = canonical_json_bytes(validated)
        for path in sorted(retained_files):
            self.publish_bytes(path, retained_files[path])
        for path in sorted(path_files):
            published = self.publish_path(path, path_files[path])
            claim = expected[path]
            _require(published["size"] == claim["size"] and
                     published["sha256"] == claim["sha256"],
                     f"retained lane file {path!r} changed before publication")
        self.publish_bytes(terminal, terminal_content)
        self.publish_bytes(_TREE_METADATA_PATH, b"")
        self.publish_bytes(_SHA256SUMS_PATH, b"")
        metadata_descriptor = self._open_placeholder(_TREE_METADATA_PATH)
        ledger_descriptor = -1
        try:
            ledger_descriptor = self._open_placeholder(_SHA256SUMS_PATH)
        except BaseException:
            os.close(metadata_descriptor)
            raise
        expected_files = set(expected) | {
            _TREE_METADATA_PATH, _SHA256SUMS_PATH}
        expected_directories = _derived_directories(sorted(expected_files))
        try:
            _require(self._directories == expected_directories,
                     "lane contains an unexpected derived directory")
            self._apply_final_modes(expected_files, expected_directories)
            nodes, _digests = self._scan_tree(
                expected_files, expected_directories, hash_files=False)
            metadata_content = canonical_json_bytes(
                expected_tree_metadata(nodes))
            self._replace_placeholder(
                _TREE_METADATA_PATH, metadata_descriptor, metadata_content)
            _nodes, digests = self._scan_tree(
                expected_files, expected_directories)
            ledger_content = expected_sha256sums({
                path: digest for path, digest in digests.items()
                if path != _SHA256SUMS_PATH
            })
            self._replace_placeholder(
                _SHA256SUMS_PATH, ledger_descriptor, ledger_content)
            _nodes, digests = self._scan_tree(
                expected_files, expected_directories)
            _require(all(
                digests[path] == identity["sha256"]
                for path, identity in self._published.items()),
                "lane files changed after publication")
        finally:
            os.close(metadata_descriptor)
            os.close(ledger_descriptor)
        self._sealed = True
        return {
            "root": self.root_path,
            "terminal": terminal,
            "terminal_record_sha256": validated["record_sha256"],
            "file_count": len(expected_files),
            "directory_count": len(expected_directories),
            "tree_metadata_sha256": digests[_TREE_METADATA_PATH],
            "sha256sums_sha256": digests[_SHA256SUMS_PATH],
        }


def _verification_file_identity(
    role: str, path: str, *, executable: bool,
) -> dict[str, Any]:
    """Hash one canonical verifier trust input with stable pathname binding."""
    _require(type(role) is str and role in (
                 "verifier", "identity_contract", "record_contract", "python"),
             "verification trust role is invalid")
    absolute = _portable_absolute_path(path, f"{role} verification program")
    _require(os.path.realpath(absolute) == absolute,
             f"{role} verification program is not canonical")
    descriptor = -1
    try:
        descriptor = os.open(
            absolute, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        before = os.fstat(descriptor)
        pathname = os.lstat(absolute)
        mode = stat.S_IMODE(before.st_mode)
        _require(stat.S_ISREG(before.st_mode) and
                 0 < before.st_size <= MAX_SEALED_FILE_BYTES and
                 (not executable or mode & 0o111 != 0) and
                 (before.st_dev, before.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 f"{role} verification program is not one stable file")
        digest = _hash_fd(descriptor, before.st_size, role)
        after = os.fstat(descriptor)
        pathname_after = os.lstat(absolute)
        fingerprint = (before.st_dev, before.st_ino, before.st_size,
                       before.st_mtime_ns, before.st_ctime_ns)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) == fingerprint and
                 (pathname_after.st_dev, pathname_after.st_ino,
                  pathname_after.st_size, pathname_after.st_mtime_ns,
                  pathname_after.st_ctime_ns) == fingerprint,
                 f"{role} verification program changed while hashed")
        return {
            "role": role,
            "path": absolute,
            "size": before.st_size,
            "mode": mode,
            "sha256": digest,
        }
    except OSError as error:
        raise AcquisitionError(
            f"cannot inspect {role} verification program: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def verification_program_identity(
    plan: LanePlan,
    toolchain: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind the verifier, its two pure contracts, and its interpreter."""
    canonical_plan = validate_lane_plan(plan)
    tools_directory = Path(__file__).resolve().parent
    verifier_path = str(
        tools_directory / "leopard2_exact_main_baseline_verifier.py")
    identity_path = str(
        tools_directory / "leopard2_exact_main_baseline.py")
    record_path = str(
        tools_directory / "leopard2_exact_main_baseline_record.py")
    _require(canonical_plan.verifier == verifier_path,
             "lane verifier is not the source-bound independent verifier")
    python_record: Mapping[str, Any] | None = None
    interpreter = os.path.realpath(sys.executable)
    if toolchain is not None:
        _require(type(toolchain) is dict and
                 type(toolchain.get("tools")) is list,
                 "verification toolchain is invalid")
        matches = [item for item in toolchain["tools"]
                   if type(item) is dict and item.get("role") == "python"]
        _require(len(matches) == 1,
                 "verification toolchain has no unique Python interpreter")
        python_record = matches[0]
        interpreter = _portable_absolute_path(
            python_record.get("resolved_path"),
            "verification Python interpreter")
    files = [
        _verification_file_identity(
            "verifier", verifier_path, executable=False),
        _verification_file_identity(
            "identity_contract", identity_path, executable=False),
        _verification_file_identity(
            "record_contract", record_path, executable=False),
        _verification_file_identity(
            "python", interpreter, executable=True),
    ]
    if python_record is not None:
        observed = files[-1]
        _require(
            observed["path"] == python_record.get("resolved_path") and
            observed["size"] == python_record.get("size") and
            observed["mode"] == python_record.get("mode") and
            observed["sha256"] == python_record.get("sha256"),
            "verification interpreter differs from the acquired toolchain")
    return {
        "schema": "leopard2-exact-main-verification-program/v1",
        "files": files,
    }


def verification_argv(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    interpreter: str,
    lane_root: str | None = None,
) -> list[str]:
    """Build and constrain one injectable independent-verifier argv."""
    canonical_plan = validate_lane_plan(plan)
    target = canonical_plan.lane_root if lane_root is None else \
        _portable_absolute_path(lane_root, "verification lane root")
    _require(target in (
                 canonical_plan.lane_root, canonical_plan.failure_lane_root),
             "verification lane root is outside the acquisition plan")
    python = _portable_absolute_path(interpreter, "verification interpreter")
    argv = _command_arguments(environment.independent_verifier_argv(
        interpreter=python, verifier=canonical_plan.verifier,
        lane_root=target))
    _require(argv[0] == python and argv[-1] == target,
             "independent verifier argv changed its interpreter or lane")
    return argv


def verification_plan_log_bytes(
    *,
    argv: Sequence[str],
    program: Mapping[str, Any],
    cwd: str,
    child_environment: Mapping[str, str],
) -> bytes:
    """Freeze the stage-4 policy before the authority lane is sealed."""
    arguments = _command_arguments(argv)
    working_directory = _portable_absolute_path(cwd, "verification plan cwd")
    values = _command_environment(child_environment)
    content = canonical_json_bytes({
        "schema": VERIFICATION_LOG_SCHEMA,
        "status": "planned",
        "program": copy.deepcopy(dict(program)),
        "argv": arguments,
        "cwd": working_directory,
        "environment": [
            {"name": name, "value": values[name]} for name in sorted(values)
        ],
        "policy": {
            "expected_exit_status": 0,
            "expected_verdict_schema":
                record_contract.VERIFIER_VERDICT_SCHEMA,
            "expected_outcome": "promoted_authority",
            "expected_promoted": True,
            "bind_record_sha256": True,
            "bind_seal_counts_and_digests": True,
            "canonical_stdout_required": True,
            "trust_inputs_rehashed_after_run": True,
        },
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "independent verification plan log exceeds its byte bound")
    return content


def seal_plan_log_bytes(
    plan: LanePlan,
    *,
    terminal: str,
    retained_byte_paths: Sequence[str],
    retained_path_sources: Sequence[str],
) -> bytes:
    """Freeze the stage-5 path and final-mode policy without circular hashes."""
    canonical_plan = validate_lane_plan(plan)
    _require(terminal == "baseline-authority.json",
             "seal plan terminal is invalid")
    byte_paths = sorted(
        _safe_relative_path(path, "seal plan retained byte path")
        for path in retained_byte_paths)
    path_sources = sorted(
        _safe_relative_path(path, "seal plan retained source path")
        for path in retained_path_sources)
    _require(len(byte_paths) == len(set(byte_paths)) and
             len(path_sources) == len(set(path_sources)) and
             set(byte_paths).isdisjoint(path_sources),
             "seal plan retained paths overlap")
    content = canonical_json_bytes({
        "schema": SEAL_LOG_SCHEMA,
        "status": "planned",
        "seal_protocol": record_contract.SEAL_PROTOCOL,
        "lane_root": canonical_plan.lane_root,
        "attempt": canonical_plan.attempt,
        "terminal_relative_path": terminal,
        "retained_byte_paths": byte_paths,
        "retained_path_sources": path_sources,
        "retained_byte_count": len(byte_paths),
        "retained_path_source_count": len(path_sources),
        "file_mode": f"{LANE_FILE_MODE:04o}",
        "directory_mode": f"{LANE_DIRECTORY_MODE:04o}",
        "final_mode_policy": "owner-only tree with all write bits removed",
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "seal plan log exceeds its byte bound")
    return content


def assembly_failure_log_bytes(
    build_stage_log: bytes, error: BaseException,
    *,
    authority_lane_root: str | None = None,
    authority_record_sha256: str | None = None,
) -> bytes:
    """Represent a post-build gate/seal failure without a contradictory log."""
    _require(type(build_stage_log) is bytes and
             0 < len(build_stage_log) <= MAX_SOURCE_LOG_BYTES,
             "completed path-variant build log is invalid")
    _require((authority_lane_root is None) ==
             (authority_record_sha256 is None),
             "assembly failure authority binding is incomplete")
    authority = None
    if authority_lane_root is not None:
        root = _portable_absolute_path(
            authority_lane_root, "assembly failure authority lane")
        _require(type(authority_record_sha256) is str and
                 re.fullmatch(
                     r"[0-9a-f]{64}", authority_record_sha256) is not None,
                 "assembly failure authority SHA-256 is invalid")
        authority = {
            "lane_root": root,
            "record_sha256": authority_record_sha256,
        }
    content = canonical_json_bytes({
        "schema": ASSEMBLY_LOG_SCHEMA,
        "status": "failed",
        "authority": authority,
        "completed_build_log": _bytes_identity(
            "diagnostics/completed-path-variant-build.log", build_stage_log),
        "error": _safe_failure_message(error),
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "authority assembly failure log exceeds its byte bound")
    return content


def verification_failure_log_bytes(
    *,
    planned_log: bytes,
    argv: Sequence[str],
    program: Mapping[str, Any],
    error: VerificationError,
    authority_record_sha256: str,
    authority_seal: Mapping[str, Any],
) -> bytes:
    """Record a stage-4 failure with exact retained child-stream claims."""
    _require(type(planned_log) is bytes and planned_log and
             isinstance(error, VerificationError),
             "verification failure log inputs are invalid")
    content = canonical_json_bytes({
        "schema": VERIFICATION_FAILURE_LOG_SCHEMA,
        "status": "failed",
        "planned_log": _bytes_identity(
            "verification/planned-independent-verification.log",
            planned_log),
        "program": copy.deepcopy(dict(program)),
        "argv": _command_arguments(argv),
        "authority": {
            "lane_root": authority_seal["root"],
            "record_sha256": authority_record_sha256,
            "tree_metadata_sha256":
                authority_seal["tree_metadata_sha256"],
            "sha256sums_sha256": authority_seal["sha256sums_sha256"],
        },
        "result": {
            "exit_status": error.observed_exit_status,
            "stdout": _bytes_identity(
                "verification/verifier.stdout", error.stdout),
            "stderr": _bytes_identity(
                "verification/verifier.stderr", error.stderr),
        },
        "error": _safe_failure_message(error),
    })
    _require(0 < len(content) <= MAX_SOURCE_LOG_BYTES,
             "verification failure log exceeds its byte bound")
    return content


def assemble_authority_record(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    host: Mapping[str, Any],
    source_stage: SourceStageResult,
    results: Sequence[BuildStageResult],
    stage_logs: Mapping[str, bytes],
) -> dict[str, Any]:
    """Derive every cross-build gate and call the sole authority constructor."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment) and
             isinstance(source_stage, SourceStageResult) and
             type(stage_logs) is dict and
             list(stage_logs) == list(record_contract.STAGE_NAMES) and
             type(results) in (list, tuple) and len(results) == 3,
             "authority assembly inputs are invalid")
    role_results = {result.role: result for result in results
                    if isinstance(result, BuildStageResult)}
    _require(tuple(role_results) == record_contract.BUILD_ROLES,
             "authority build result order changed")
    first = role_results["canonical_first"]
    second = role_results["canonical_second"]
    variant = role_results["path_variant"]
    canonical_identity_equal = identity_contract.exact_json_equal(
        first.identity, second.identity)
    canonical_archive_equal = (
        first.build["archive"]["size"] == second.build["archive"]["size"] and
        first.build["archive"]["sha256"] ==
        second.build["archive"]["sha256"])
    variant_raw_differs = (
        variant.identity["artifact"]["sha256"] !=
        first.identity["artifact"]["sha256"])
    combined = first.identity["combined_sha256"]
    normalized_match = all(
        result.identity["combined_sha256"] == combined for result in results)
    census_zero = all(_identity_census_is_zero(result.identity)
                      for result in results)
    failures = []
    if not canonical_identity_equal:
        failures.append("same-path normalized identities differ")
    if not canonical_archive_equal:
        failures.append("same-path archives differ")
    if not variant_raw_differs:
        failures.append("path-variant raw executable did not differ")
    if not normalized_match:
        failures.append("path-variant normalized identity changed")
    if not census_zero:
        failures.append("selected ELF sections retain an acquisition root")
    _require(not failures,
             "authority promotion gates failed: " + "; ".join(failures))
    controller_path = "controllers/test_legacy_main_benchmark.py"
    _require(controller_path in source_stage.retained_bytes,
             "source stage omitted the attestation controller")
    controller = _bytes_identity(
        controller_path, source_stage.retained_bytes[controller_path])
    lane = {
        "root": canonical_plan.lane_root,
        "attempt": canonical_plan.attempt,
        "attempt_budget": 3,
        "record_relative_path": "baseline-authority.json",
        "seal_protocol": record_contract.SEAL_PROTOCOL,
        "stages": [{
            "name": name,
            "status": "complete",
            "log": _bytes_identity(
                f"logs/{index:02d}-{name}.log", stage_logs[name]),
        } for index, name in enumerate(record_contract.STAGE_NAMES)],
    }
    identities = {role: role_results[role].identity
                  for role in record_contract.BUILD_ROLES}
    return record_contract.baseline_authority_record(
        created_utc=environment.now_utc(), host=copy.deepcopy(dict(host)),
        lane=lane, source=source_stage.source, adapter=source_stage.adapter,
        toolchain=source_stage.toolchain,
        builds={role: role_results[role].build
                for role in record_contract.BUILD_ROLES},
        runtime_closure={
            "schema": record_contract.RUNTIME_CLOSURE_SCHEMA,
            "normalization": record_contract.CANONICAL_LDD_NORMALIZATION,
            "records": [role_results[role].runtime
                        for role in record_contract.BUILD_ROLES],
        },
        attestation={
            "schema": record_contract.ATTESTATION_SCHEMA,
            "test_controller": controller,
            "records": [role_results[role].attestation
                        for role in record_contract.BUILD_ROLES],
        },
        identity={
            **identities,
            "combined_sha256": combined,
            "canonical_raw_bytes_identical": canonical_identity_equal,
            "path_variant_raw_bytes_differ": variant_raw_differs,
            "normalized_match": normalized_match,
        },
    )


def _verification_program_interpreter(
    value: Mapping[str, Any],
) -> str:
    _require(type(value) is dict and set(value) == {"schema", "files"} and
             value["schema"] == "leopard2-exact-main-verification-program/v1" and
             type(value["files"]) is list and len(value["files"]) == 4,
             "independent verification program identity is invalid")
    expected_roles = (
        "verifier", "identity_contract", "record_contract", "python")
    _require(tuple(item.get("role") if type(item) is dict else None
                   for item in value["files"]) == expected_roles,
             "independent verification program role order changed")
    for index, item in enumerate(value["files"]):
        _require(type(item) is dict and set(item) == {
                     "role", "path", "size", "mode", "sha256"} and
                 type(item["size"]) is int and item["size"] > 0 and
                 type(item["mode"]) is int and 0 <= item["mode"] <= 0o7777 and
                 type(item["sha256"]) is str and
                 re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is not None,
                 f"independent verification program file {index} is invalid")
        _portable_absolute_path(
            item["path"], f"independent verification program file {index}")
    return value["files"][-1]["path"]


def _load_sealed_terminal_record(
    lane_root: str, terminal: str,
) -> dict[str, Any]:
    """Read one just-sealed terminal without relaxing its independent replay."""
    root = _portable_absolute_path(lane_root, "sealed terminal lane")
    _require(terminal in _TERMINAL_PATHS,
             "sealed terminal name is invalid")
    path = root + "/" + terminal
    descriptor = -1
    try:
        descriptor = os.open(
            path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        status = os.fstat(descriptor)
        pathname = os.lstat(path)
        _require(stat.S_ISREG(status.st_mode) and status.st_nlink == 1 and
                 status.st_uid == os.geteuid() and
                 status.st_gid == os.getegid() and
                 stat.S_IMODE(status.st_mode) == LANE_FILE_MODE and
                 0 < status.st_size <= MAX_SOURCE_LOG_BYTES and
                 (status.st_dev, status.st_ino) ==
                 (pathname.st_dev, pathname.st_ino),
                 "sealed terminal is not one owner-only regular file")
        chunks: list[bytes] = []
        remaining = status.st_size
        while remaining:
            chunk = os.read(descriptor, min(READ_CHUNK, remaining))
            _require(bool(chunk), "sealed terminal ended while read")
            chunks.append(chunk)
            remaining -= len(chunk)
        _require(os.read(descriptor, 1) == b"",
                 "sealed terminal grew while read")
        content = b"".join(chunks)
    except OSError as error:
        raise AcquisitionError(
            f"cannot read sealed terminal {terminal!r}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    try:
        if terminal == "baseline-authority.json":
            return record_contract.load_baseline_authority_record(content)
        return record_contract.load_baseline_failure_record(content)
    except ExactMainBaselineError as error:
        raise AcquisitionError("sealed terminal record is invalid") from error


def run_independent_verification(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    record: Mapping[str, Any],
    seal: Mapping[str, Any] | None,
    toolchain: Mapping[str, Any] | None = None,
    program: Mapping[str, Any] | None = None,
    lane_root: str | None = None,
) -> dict[str, Any]:
    """Launch and fully cross-bind the source-bound verifier process."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment) and
             type(record) is dict and (seal is None or type(seal) is dict),
             "independent verification inputs are invalid")
    try:
        if record.get("schema") == record_contract.AUTHORITY_SCHEMA:
            validated = record_contract.validate_baseline_authority_record(
                record)
            expected_terminal = "baseline-authority.json"
            expected_status = 0
            expected_outcome = "promoted_authority"
            expected_promoted = True
        else:
            validated = record_contract.validate_baseline_failure_record(
                record)
            expected_terminal = "FAILED.json"
            expected_status = 3
            expected_outcome = "verified_failure"
            expected_promoted = False
    except ExactMainBaselineError as error:
        raise AcquisitionError(
            "independent verification record is invalid") from error
    target = canonical_plan.lane_root if lane_root is None else \
        _portable_absolute_path(lane_root, "independent verification lane")
    _require(target in (
                 canonical_plan.lane_root, canonical_plan.failure_lane_root) and
             validated["lane"]["root"] == target,
             "independent verification target differs from its record")
    if seal is not None:
        _require(set(seal) == {
                     "root", "terminal", "terminal_record_sha256",
                     "file_count", "directory_count",
                     "tree_metadata_sha256", "sha256sums_sha256"} and
                 seal["root"] == target and
                 seal["terminal"] == expected_terminal and
                 seal["terminal_record_sha256"] ==
                    validated["record_sha256"],
                 "independent verification target differs from its seal")
    before = verification_program_identity(canonical_plan, toolchain) \
        if program is None else copy.deepcopy(dict(program))
    interpreter = _verification_program_interpreter(before)
    argv = verification_argv(
        environment, canonical_plan, interpreter=interpreter,
        lane_root=target)
    child_environment = frozen_child_environment()
    try:
        result = environment.run(
            argv, cwd=canonical_plan.scratch_root, env=child_environment,
            timeout=VERIFICATION_TIMEOUT_SECONDS,
            maximum_bytes=MAX_COMMAND_OUTPUT_BYTES)
    except BaseException as error:
        raise VerificationError(
            "independent verifier could not run: " +
            _safe_failure_message(error), exit_status=1,
            stdout=b"", stderr=b"", process_result_observed=False) from error
    if not isinstance(result, CommandResult):
        raise VerificationError(
            "independent verifier returned an invalid process result",
            exit_status=1, stdout=b"", stderr=b"",
            process_result_observed=False)
    if type(result.exit_status) is not int or not (
            0 <= result.exit_status <= 255) or \
            type(result.stdout) is not bytes or \
            type(result.stderr) is not bytes or \
            len(result.stdout) + len(result.stderr) > MAX_COMMAND_OUTPUT_BYTES:
        raise VerificationError(
            "independent verifier returned malformed process metadata",
            exit_status=(result.exit_status
                         if type(result.exit_status) is int else 1),
            stdout=(result.stdout if type(result.stdout) is bytes else b""),
            stderr=(result.stderr if type(result.stderr) is bytes else b""),
            process_result_observed=False)

    def reject(message: str) -> NoReturn:
        raise VerificationError(
            message, exit_status=result.exit_status,
            stdout=result.stdout, stderr=result.stderr)

    try:
        if result.exit_status != expected_status:
            reject(
                f"independent verifier returned status {result.exit_status}; "
                f"expected {expected_status}")
        if result.stderr != b"":
            reject("independent verifier wrote stderr")
        verdict = identity_contract.strict_json_loads(
            result.stdout, "independent verifier verdict JSON")
        if canonical_json_bytes(verdict) != result.stdout:
            reject("independent verifier verdict is not canonical JSON")
        if type(verdict) is not dict or set(verdict) != {
                "schema", "outcome", "promoted", "record", "seal",
                "recomputed", "producer_attested", "verdict_sha256"}:
            reject("independent verifier verdict has an unexpected key set")
        if verdict["schema"] != record_contract.VERIFIER_VERDICT_SCHEMA or \
                verdict["outcome"] != expected_outcome or \
                verdict["promoted"] is not expected_promoted:
            reject("independent verifier returned the wrong terminal outcome")
        record_view = verdict["record"]
        if type(record_view) is not dict or set(record_view) != {
                "relative_path", "schema", "record_sha256",
                "canonical_bytes_identical"} or \
                record_view["relative_path"] != expected_terminal or \
                record_view["schema"] != validated["schema"] or \
                record_view["record_sha256"] != validated["record_sha256"] or \
                record_view["canonical_bytes_identical"] is not True:
            reject("independent verifier did not bind the terminal record")
        seal_view = verdict["seal"]
        if type(seal_view) is not dict or set(seal_view) != {
                "protocol", "tree_metadata_schema", "tree_metadata_sha256",
                "sha256sums_sha256", "directory_count", "file_count",
                "checksum_line_count", "metadata_entry_count"} or \
                seal_view["protocol"] != record_contract.SEAL_PROTOCOL or \
                seal_view["tree_metadata_schema"] != TREE_METADATA_SCHEMA or \
                type(seal_view["file_count"]) is not int or \
                seal_view["file_count"] <= 0 or \
                type(seal_view["directory_count"]) is not int or \
                seal_view["directory_count"] <= 0 or \
                type(seal_view["checksum_line_count"]) is not int or \
                seal_view["checksum_line_count"] != \
                    seal_view["file_count"] - 1 or \
                type(seal_view["metadata_entry_count"]) is not int or \
                seal_view["metadata_entry_count"] != \
                    (seal_view["file_count"] +
                     seal_view["directory_count"] - 1) or \
                any(type(seal_view[name]) is not str or
                    re.fullmatch(r"[0-9a-f]{64}", seal_view[name]) is None
                    for name in (
                        "tree_metadata_sha256", "sha256sums_sha256")):
            reject("independent verifier did not bind the sealed tree")
        if seal is not None and (
                seal_view["tree_metadata_sha256"] !=
                    seal["tree_metadata_sha256"] or
                seal_view["sha256sums_sha256"] !=
                    seal["sha256sums_sha256"] or
                seal_view["file_count"] != seal["file_count"] or
                seal_view["directory_count"] != seal["directory_count"]):
            reject("independent verifier did not bind the producer seal")
        expected_verdict_sha256 = identity_contract.canonical_sha256({
            key: copy.deepcopy(value) for key, value in verdict.items()
            if key != "verdict_sha256"
        })
        if verdict["verdict_sha256"] != expected_verdict_sha256:
            reject("independent verifier verdict SHA-256 is inconsistent")
        after = verification_program_identity(canonical_plan, toolchain)
        if not identity_contract.exact_json_equal(before, after):
            reject("independent verification program changed during replay")
    except VerificationError:
        raise
    except BaseException as error:
        reject("independent verifier output is invalid: " +
               _safe_failure_message(error))
    return copy.deepcopy(verdict)


def _seal_result_from_verdict(
    lane_root: str,
    terminal: str,
    record: Mapping[str, Any],
    verdict: Mapping[str, Any],
) -> dict[str, Any]:
    """Recover the producer seal view only after trusted independent replay."""
    root = _portable_absolute_path(lane_root, "recovered authority lane")
    _require(terminal in _TERMINAL_PATHS and type(record) is dict and
             type(verdict) is dict and type(verdict.get("seal")) is dict and
             type(record.get("record_sha256")) is str,
             "recovered authority seal inputs are invalid")
    seal = verdict["seal"]
    return {
        "root": root,
        "terminal": terminal,
        "terminal_record_sha256": record["record_sha256"],
        "file_count": seal["file_count"],
        "directory_count": seal["directory_count"],
        "tree_metadata_sha256": seal["tree_metadata_sha256"],
        "sha256sums_sha256": seal["sha256sums_sha256"],
    }


def _merge_retained_evidence(
    source_stage: SourceStageResult | None,
    results: Sequence[BuildStageResult],
) -> tuple[dict[str, bytes], dict[str, str]]:
    retained_bytes: dict[str, bytes] = {}
    retained_paths: dict[str, str] = {}

    def merge_bytes(values: Mapping[str, bytes]) -> None:
        _require(type(values) is dict,
                 "retained byte evidence is not a mapping")
        for path, content in values.items():
            relative = _safe_relative_path(path, "retained byte evidence path")
            _require(relative not in retained_bytes and
                     relative not in retained_paths and type(content) is bytes,
                     f"retained evidence path {relative!r} is duplicated")
            retained_bytes[relative] = content

    def merge_paths(values: Mapping[str, str]) -> None:
        _require(type(values) is dict,
                 "retained path evidence is not a mapping")
        for path, source in values.items():
            relative = _safe_relative_path(path, "retained path evidence path")
            _require(relative not in retained_bytes and
                     relative not in retained_paths and type(source) is str,
                     f"retained evidence path {relative!r} is duplicated")
            retained_paths[relative] = source

    if source_stage is not None:
        _require(isinstance(source_stage, SourceStageResult),
                 "source-stage retained evidence is invalid")
        merge_bytes(source_stage.retained_bytes)
        merge_paths(source_stage.retained_paths)
    for result in results:
        _require(isinstance(result, BuildStageResult),
                 "build-stage retained evidence is invalid")
        merge_bytes(result.retained_bytes)
        merge_paths(result.retained_paths)
    return retained_bytes, retained_paths


def _seal_and_verify_failure(
    environment: HostEnvironment,
    plan: LanePlan,
    error: BaseException,
    *,
    stage: str,
    stage_logs: Mapping[str, bytes],
    retained_bytes: Mapping[str, bytes] | None = None,
    retained_paths: Mapping[str, str] | None = None,
    diagnostics: Mapping[str, bytes] | None = None,
    lane_root: str | None = None,
    authority_record_sha256: str | None = None,
    validate_state: Any = None,
) -> AcquisitionOutcome:
    """Seal one failure and require a second-process nonpromotion verdict."""
    canonical_plan = validate_lane_plan(plan)
    target = canonical_plan.lane_root if lane_root is None else \
        _portable_absolute_path(lane_root, "failed acquisition lane")
    seal = seal_stage_failure(
        environment, canonical_plan, error, stage=stage,
        stage_logs=stage_logs, retained_bytes=retained_bytes,
        retained_paths=retained_paths, diagnostics=diagnostics,
        lane_root=target,
        authority_record_sha256=authority_record_sha256)
    _require(validate_state is None or callable(validate_state),
             "failure-state validator is not callable")
    if validate_state is not None:
        validate_state()
    record = _load_sealed_terminal_record(target, "FAILED.json")
    if validate_state is not None:
        validate_state()
    program = verification_program_identity(canonical_plan, None)
    verdict = run_independent_verification(
        environment, canonical_plan, record=record, seal=seal,
        program=program, lane_root=target)
    if validate_state is not None:
        validate_state()
    return AcquisitionOutcome(
        outcome=("verification_failure" if authority_record_sha256 is not None
                 else "acquisition_failure"),
        lane_root=target, terminal="FAILED.json",
        seal=copy.deepcopy(seal), verdict=copy.deepcopy(verdict),
        authority_lane_root=(canonical_plan.lane_root
                             if authority_record_sha256 is not None else None),
        authority_record_sha256=authority_record_sha256,
        verification_program=copy.deepcopy(program),
        verification_exit_status=3,
        message=_safe_failure_message(error),
    )


def acquire_exact_main_baseline(
    environment: HostEnvironment,
    plan: LanePlan,
    *,
    mode: str,
) -> AcquisitionOutcome:
    """Run the exact no-performance-timing acquisition state machine."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(environment, HostEnvironment),
             "acquisition environment is invalid")
    _require(type(mode) is str and mode in ("rehearsal", "authoritative"),
             "acquisition mode is invalid")

    with AcquisitionLocks(environment, blocking=False) as locks:
        with prepare_acquisition_roots(canonical_plan) as prepared:
            def checkpoint() -> None:
                locks.validate_current()
                prepared.validate_current()

            checkpoint()
            try:
                host = environment.host_facts()
            except BaseException as error:
                source_log = _source_stage_log_bytes(
                    "failed", (), adapter_commit=None, adapter_tree=None,
                    retained_byte_paths=(), retained_path_sources=(),
                    error=error)
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, error,
                    stage="source_acquisition",
                    stage_logs={"source_acquisition": source_log},
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()

            try:
                source_stage = acquire_source_stage(
                    environment, canonical_plan, prepared)
            except SourceStageError as error:
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, error,
                    stage="source_acquisition",
                    stage_logs={"source_acquisition": error.log},
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            except BaseException as error:
                wrapped = SourceStageError(
                    error, _source_stage_log_bytes(
                        "failed", (), adapter_commit=None, adapter_tree=None,
                        retained_byte_paths=(), retained_path_sources=(),
                        error=error))
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, wrapped,
                    stage="source_acquisition",
                    stage_logs={"source_acquisition": wrapped.log},
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()

            controller_path = "controllers/test_legacy_main_benchmark.py"
            _require(controller_path in source_stage.retained_bytes,
                     "source stage omitted the attestation controller")
            controller_sha256 = _sha256(
                source_stage.retained_bytes[controller_path])
            stage_logs: dict[str, bytes] = {
                "source_acquisition": source_stage.log,
            }
            results: list[BuildStageResult] = []

            def build_role(role: str) -> AcquisitionOutcome | None:
                stage = role + "_build"
                try:
                    result = acquire_build_stage(
                        environment, canonical_plan, prepared, role=role,
                        toolchain=source_stage.toolchain,
                        controller_sha256=controller_sha256)
                except BuildStageError as error:
                    stage_logs[stage] = error.log
                    retained, paths = _merge_retained_evidence(
                        source_stage, results)
                    checkpoint()
                    return _seal_and_verify_failure(
                        environment, canonical_plan, error, stage=stage,
                        stage_logs=stage_logs, retained_bytes=retained,
                        retained_paths=paths, validate_state=checkpoint)
                except BaseException as error:
                    failed_log = _build_stage_log_bytes(
                        role, "failed", (),
                        controller_sha256=controller_sha256, error=error)
                    wrapped = BuildStageError(role, error, failed_log)
                    stage_logs[stage] = failed_log
                    retained, paths = _merge_retained_evidence(
                        source_stage, results)
                    checkpoint()
                    return _seal_and_verify_failure(
                        environment, canonical_plan, wrapped, stage=stage,
                        stage_logs=stage_logs, retained_bytes=retained,
                        retained_paths=paths, validate_state=checkpoint)
                results.append(result)
                stage_logs[stage] = result.log
                checkpoint()
                return None

            outcome = build_role("canonical_first")
            if outcome is not None:
                checkpoint()
                return outcome

            try:
                prepared.reset_root("canonical_build_root")
            except BaseException as error:
                failed_log = _build_stage_log_bytes(
                    "canonical_second", "failed", (),
                    controller_sha256=controller_sha256, error=error)
                wrapped = BuildStageError(
                    "canonical_second", error, failed_log)
                stage_logs["canonical_second_build"] = failed_log
                retained, paths = _merge_retained_evidence(
                    source_stage, results)
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, wrapped,
                    stage="canonical_second_build", stage_logs=stage_logs,
                    retained_bytes=retained, retained_paths=paths,
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()

            outcome = build_role("canonical_second")
            if outcome is not None:
                checkpoint()
                return outcome
            outcome = build_role("path_variant")
            if outcome is not None:
                checkpoint()
                return outcome

            retained_bytes, retained_paths = _merge_retained_evidence(
                source_stage, results)
            completed_path_variant_log = stage_logs["path_variant_build"]
            try:
                program = verification_program_identity(
                    canonical_plan, source_stage.toolchain)
                interpreter = _verification_program_interpreter(program)
                verifier_arguments = verification_argv(
                    environment, canonical_plan, interpreter=interpreter)
                child_environment = frozen_child_environment()
                stage_logs["independent_verification"] = \
                    verification_plan_log_bytes(
                        argv=verifier_arguments, program=program,
                        cwd=canonical_plan.scratch_root,
                        child_environment=child_environment)
                all_log_paths = [
                    f"logs/{index:02d}-{name}.log"
                    for index, name in enumerate(record_contract.STAGE_NAMES)
                ]
                stage_logs["seal"] = seal_plan_log_bytes(
                    canonical_plan, terminal="baseline-authority.json",
                    retained_byte_paths=(
                        list(retained_bytes) + all_log_paths),
                    retained_path_sources=list(retained_paths))
                authority = assemble_authority_record(
                    environment, canonical_plan, host=host,
                    source_stage=source_stage, results=results,
                    stage_logs=stage_logs)
                for index, name in enumerate(record_contract.STAGE_NAMES):
                    path = f"logs/{index:02d}-{name}.log"
                    _require(path not in retained_bytes and
                             path not in retained_paths,
                             "authority stage log collides with evidence")
                    retained_bytes[path] = stage_logs[name]
            except BaseException as error:
                failed_log = assembly_failure_log_bytes(
                    completed_path_variant_log, error)
                failure_stage_logs = {
                    name: stage_logs[name]
                    for name in record_contract.STAGE_NAMES[:3]
                }
                failure_stage_logs["path_variant_build"] = failed_log
                diagnostics = {
                    "diagnostics/completed-path-variant-build.log":
                        completed_path_variant_log,
                }
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, error,
                    stage="path_variant_build",
                    stage_logs=failure_stage_logs,
                    retained_bytes=retained_bytes,
                    retained_paths=retained_paths,
                    diagnostics=diagnostics, validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()

            try:
                with LaneWriter(canonical_plan.lane_root) as writer:
                    authority_seal = writer.seal_record(
                        authority, retained_bytes,
                        retained_paths=retained_paths)
            except BaseException as error:
                primary_verification_error: VerificationError | None = None
                primary_terminal = canonical_plan.lane_root + \
                    "/baseline-authority.json"
                if not _path_is_absent(primary_terminal):
                    # A writer exception can occur after every byte and final
                    # mode is already correct.  Let the source-bound verifier
                    # classify that state before emitting any contrary lane.
                    checkpoint()
                    try:
                        recovered_verdict = run_independent_verification(
                            environment, canonical_plan, record=authority,
                            seal=None, toolchain=source_stage.toolchain,
                            program=program)
                    except VerificationError as recovery_error:
                        if not recovery_error.process_result_observed:
                            raise AcquisitionError(
                                "authority seal failed and the primary lane "
                                "could not be independently classified: " +
                                _safe_failure_message(error) + "; " +
                                _safe_failure_message(recovery_error)) \
                                from recovery_error
                        primary_verification_error = recovery_error
                    except BaseException as recovery_error:
                        raise AcquisitionError(
                            "authority seal failed and the primary lane could "
                            "not be independently classified: " +
                            _safe_failure_message(error) + "; " +
                            _safe_failure_message(recovery_error)) \
                            from recovery_error
                    else:
                        recovered_seal = _seal_result_from_verdict(
                            canonical_plan.lane_root,
                            "baseline-authority.json", authority,
                            recovered_verdict)
                        checkpoint()
                        return AcquisitionOutcome(
                            outcome="promoted_authority",
                            lane_root=canonical_plan.lane_root,
                            terminal="baseline-authority.json",
                            seal=copy.deepcopy(recovered_seal),
                            verdict=copy.deepcopy(recovered_verdict),
                            authority_lane_root=canonical_plan.lane_root,
                            authority_record_sha256=
                                authority["record_sha256"],
                            verification_program=copy.deepcopy(program),
                            verification_exit_status=0, message=None)
                failed_log = assembly_failure_log_bytes(
                    completed_path_variant_log, error,
                    authority_lane_root=canonical_plan.lane_root,
                    authority_record_sha256=authority["record_sha256"])
                failure_stage_logs = {
                    name: stage_logs[name]
                    for name in record_contract.STAGE_NAMES[:3]
                }
                failure_stage_logs["path_variant_build"] = failed_log
                diagnostics = {
                    "diagnostics/completed-path-variant-build.log":
                        completed_path_variant_log,
                }
                # The primary root has been claimed and is never reused, even
                # when independent replay rejected its terminal.
                failure_retained = {
                    path: content for path, content in retained_bytes.items()
                    if not path.startswith("logs/")
                }
                if primary_verification_error is not None:
                    failure_retained.update({
                        "verification/primary-verifier.stdout":
                            primary_verification_error.stdout,
                        "verification/primary-verifier.stderr":
                            primary_verification_error.stderr,
                    })
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, error,
                    stage="path_variant_build",
                    stage_logs=failure_stage_logs,
                    retained_bytes=failure_retained,
                    retained_paths=retained_paths,
                    diagnostics=diagnostics,
                    lane_root=canonical_plan.failure_lane_root,
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()

            try:
                verdict = run_independent_verification(
                    environment, canonical_plan, record=authority,
                    seal=authority_seal, toolchain=source_stage.toolchain,
                    program=program)
            except BaseException as observed_error:
                error = observed_error if isinstance(
                    observed_error, VerificationError) else VerificationError(
                        _safe_failure_message(observed_error), exit_status=1,
                        stdout=b"", stderr=b"")
                failed_log = verification_failure_log_bytes(
                    planned_log=stage_logs["independent_verification"],
                    argv=verifier_arguments, program=program, error=error,
                    authority_record_sha256=authority["record_sha256"],
                    authority_seal=authority_seal)
                failure_stage_logs = {
                    name: stage_logs[name]
                    for name in record_contract.STAGE_NAMES[:4]
                }
                failure_stage_logs["independent_verification"] = failed_log
                failure_retained = {
                    "verification/planned-independent-verification.log":
                        stage_logs["independent_verification"],
                    "verification/verifier.stdout": error.stdout,
                    "verification/verifier.stderr": error.stderr,
                }
                checkpoint()
                outcome = _seal_and_verify_failure(
                    environment, canonical_plan, error,
                    stage="independent_verification",
                    stage_logs=failure_stage_logs,
                    retained_bytes=failure_retained,
                    lane_root=canonical_plan.failure_lane_root,
                    authority_record_sha256=authority["record_sha256"],
                    validate_state=checkpoint)
                checkpoint()
                return outcome
            checkpoint()
            return AcquisitionOutcome(
                outcome="promoted_authority",
                lane_root=canonical_plan.lane_root,
                terminal="baseline-authority.json",
                seal=copy.deepcopy(authority_seal),
                verdict=copy.deepcopy(verdict),
                authority_lane_root=canonical_plan.lane_root,
                authority_record_sha256=authority["record_sha256"],
                verification_program=copy.deepcopy(program),
                verification_exit_status=0, message=None)


_CLI_PATH_OPTIONS = (
    "lane-root", "failure-lane-root", "repository", "verifier",
    "canonical-adapter-root", "canonical-baseline-root",
    "canonical-build-root", "variant-adapter-root",
    "variant-baseline-root", "variant-build-root", "scratch-root",
)
_CLI_OPTIONS = ("mode", "attempt") + _CLI_PATH_OPTIONS
_USAGE = (
    "usage: leopard2_exact_main_baseline_acquire.py "
    "--mode {rehearsal|authoritative} --attempt {1|2|3} "
    "--lane-root P --failure-lane-root P --repository P --verifier P "
    "--canonical-adapter-root P --canonical-baseline-root P "
    "--canonical-build-root P --variant-adapter-root P "
    "--variant-baseline-root P --variant-build-root P --scratch-root P"
)


def _parse_acquisition_cli(
    arguments: Sequence[str],
) -> tuple[str, LanePlan]:
    _require(type(arguments) in (list, tuple) and
             len(arguments) == 2 * len(_CLI_OPTIONS),
             "acquisition CLI requires every option exactly once")
    values: dict[str, str] = {}
    for index in range(0, len(arguments), 2):
        option = arguments[index]
        value = arguments[index + 1]
        _require(type(option) is str and option.startswith("--") and
                 option[2:] in _CLI_OPTIONS and option[2:] not in values and
                 type(value) is str and bool(value),
                 "acquisition CLI has an unknown, duplicate, or empty option")
        values[option[2:]] = value
    _require(set(values) == set(_CLI_OPTIONS),
             "acquisition CLI option set is incomplete")
    _require(values["mode"] in ("rehearsal", "authoritative"),
             "acquisition mode must be rehearsal or authoritative")
    _require(values["attempt"] in ("1", "2", "3"),
             "acquisition attempt must be 1, 2, or 3")
    plan = LanePlan(
        lane_root=values["lane-root"],
        failure_lane_root=values["failure-lane-root"],
        attempt=int(values["attempt"], 10),
        repository=values["repository"], verifier=values["verifier"],
        canonical_adapter_root=values["canonical-adapter-root"],
        canonical_baseline_root=values["canonical-baseline-root"],
        canonical_build_root=values["canonical-build-root"],
        variant_adapter_root=values["variant-adapter-root"],
        variant_baseline_root=values["variant-baseline-root"],
        variant_build_root=values["variant-build-root"],
        scratch_root=values["scratch-root"],
    )
    return values["mode"], validate_lane_plan(plan)


def acquisition_report(
    plan: LanePlan,
    outcome: AcquisitionOutcome,
    *,
    mode: str,
) -> dict[str, Any]:
    """Build the canonical no-performance-timing operator report."""
    canonical_plan = validate_lane_plan(plan)
    _require(isinstance(outcome, AcquisitionOutcome) and
             type(mode) is str and mode in ("rehearsal", "authoritative") and
             outcome.outcome in (
                 "promoted_authority", "acquisition_failure",
                 "verification_failure") and
             outcome.lane_root in (
                 canonical_plan.lane_root, canonical_plan.failure_lane_root) and
             outcome.terminal in _TERMINAL_PATHS and
             type(outcome.seal) is dict and
             outcome.seal.get("root") == outcome.lane_root and
             outcome.seal.get("terminal") == outcome.terminal and
             type(outcome.verdict) is dict and
             type(outcome.verification_program) is dict and
             outcome.verification_exit_status in (0, 3),
             "acquisition outcome cannot be reported")
    _verification_program_interpreter(outcome.verification_program)
    verdict_sha256 = outcome.verdict.get("verdict_sha256")
    _require(type(verdict_sha256) is str and
             re.fullmatch(r"[0-9a-f]{64}", verdict_sha256) is not None,
             "acquisition verdict SHA-256 is invalid")
    return {
        "schema": ACQUISITION_REPORT_SCHEMA,
        "mode": mode,
        "authoritative": mode == "authoritative",
        "outcome": outcome.outcome,
        "lane_roots": {
            "authority": canonical_plan.lane_root,
            "failure": canonical_plan.failure_lane_root,
        },
        "terminal": {
            "lane_root": outcome.lane_root,
            "relative_path": outcome.terminal,
            "record_sha256": outcome.seal["terminal_record_sha256"],
        },
        "authority_record_sha256": outcome.authority_record_sha256,
        "seal": {
            "protocol": record_contract.SEAL_PROTOCOL,
            "file_count": outcome.seal["file_count"],
            "directory_count": outcome.seal["directory_count"],
            "tree_metadata_sha256":
                outcome.seal["tree_metadata_sha256"],
            "sha256sums_sha256": outcome.seal["sha256sums_sha256"],
        },
        "verification": {
            "exit_status": outcome.verification_exit_status,
            "verdict_sha256": verdict_sha256,
            "program": copy.deepcopy(outcome.verification_program),
        },
        "message": outcome.message,
    }


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        mode, plan = _parse_acquisition_cli(arguments)
    except BaseException as error:
        sys.stderr.write(_USAGE + "\n" + _safe_failure_message(error) + "\n")
        return 2
    try:
        outcome = acquire_exact_main_baseline(
            HostEnvironment(), plan, mode=mode)
        report = acquisition_report(plan, outcome, mode=mode)
        sys.stdout.buffer.write(canonical_json_bytes(report))
        sys.stdout.buffer.flush()
    except BaseException as error:
        sys.stderr.write(
            "exact-main baseline acquisition failed: " +
            _safe_failure_message(error) + "\n")
        return 1
    if outcome.outcome != "promoted_authority":
        sys.stderr.write((outcome.message or "acquisition did not promote") +
                         "\n")
        return 3
    return 0


__all__ = (
    "AcquisitionLocks",
    "AcquisitionError",
    "AcquisitionOutcome",
    "ACQUISITION_REPORT_SCHEMA",
    "ASSEMBLY_LOG_SCHEMA",
    "ATTESTATION_TIMEOUT_SECONDS",
    "BUILD_CLOSURE_SCHEMA",
    "BUILD_LOG_SCHEMA",
    "BUILD_STAGE_LOG_SCHEMA",
    "BUILD_TIMEOUT_SECONDS",
    "BuildStageError",
    "BuildStageResult",
    "CANONICAL_LOCK_PATH",
    "CanonicalFileLock",
    "CommandExecutionError",
    "CommandResult",
    "HostEnvironment",
    "LANE_DIRECTORY_MODE",
    "LANE_FILE_MODE",
    "LanePlan",
    "LaneWriter",
    "PreparedAcquisitionRoots",
    "SourceStageError",
    "SourceStageResult",
    "StableLeaseAnchor",
    "StreamedCommandResult",
    "TREE_METADATA_SCHEMA",
    "VERIFICATION_FAILURE_LOG_SCHEMA",
    "VERIFICATION_LOG_SCHEMA",
    "VERIFICATION_TIMEOUT_SECONDS",
    "VerificationError",
    "SEAL_LOG_SCHEMA",
    "acquire_exact_main_baseline",
    "acquire_build_stage",
    "acquire_source_stage",
    "adapter_inventory",
    "acquisition_report",
    "assemble_authority_record",
    "assembly_failure_log_bytes",
    "build_closure_document",
    "build_role_roots",
    "build_root_census",
    "canonical_git_archive",
    "canonical_ldd_text",
    "expected_sha256sums",
    "expected_tree_metadata",
    "frozen_child_environment",
    "normalize_ldd_output",
    "parse_ctest_success_summary",
    "prepare_acquisition_roots",
    "resolve_toolchain",
    "seal_stage_failure",
    "seal_source_acquisition_failure",
    "seal_verification_failure",
    "seal_plan_log_bytes",
    "stage_build_output",
    "stage_detached_source",
    "validate_lane_plan",
    "verification_argv",
    "verification_failure_log_bytes",
    "verification_plan_log_bytes",
    "verification_program_identity",
    "run_independent_verification",
)


if __name__ == "__main__":
    raise SystemExit(main())
