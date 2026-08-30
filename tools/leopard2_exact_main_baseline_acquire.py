#!/usr/bin/env python3
"""Acquisition primitives for the sealed exact-main Leopard1 baseline.

The evidence schemas live in :mod:`leopard2_exact_main_baseline_record` and the
read-only consumer lives in :mod:`leopard2_exact_main_baseline_verifier`.  This
module is the producer-side boundary.  In particular, it deliberately does not
import the verifier: a completed lane must be checked by launching that separate
program after the owner-only seal is complete.

The first public layer is intentionally host-independent.  It validates the
seven acquisition roots and publishes an already-constructed authority or
failure record through an fd-anchored, exclusive, crash-conservative lane
writer.  The build/acquisition state machine is layered on these primitives so
that fault-injection tests and the eventual real acquisition use exactly the
same sealing implementation.
"""

from __future__ import annotations

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
import signal
import stat
import subprocess
import sys
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
BUILD_CLOSURE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-build-closure/v1"
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
    """The seven immutable roots and attempt identity for one acquisition."""

    lane_root: str
    attempt: int
    repository: str
    verifier: str
    canonical_adapter_root: str
    canonical_baseline_root: str
    canonical_build_root: str
    variant_adapter_root: str
    variant_baseline_root: str
    variant_build_root: str


_PLAN_PATH_FIELDS = (
    "lane_root",
    "repository",
    "verifier",
    "canonical_adapter_root",
    "canonical_baseline_root",
    "canonical_build_root",
    "variant_adapter_root",
    "variant_baseline_root",
    "variant_build_root",
)


def validate_lane_plan(value: Any) -> LanePlan:
    """Validate and detach one acquisition plan.

    The lane plus the six acquisition roots are mutually non-containing both
    component-wise and as raw UTF-8 byte strings.  The latter is deliberately
    stricter than the record contract because the ELF census searches for exact
    path bytes and must never attribute one root's occurrence to another root.
    """
    _require(isinstance(value, LanePlan), "lane plan has the wrong type")
    paths = {
        field: _portable_absolute_path(getattr(value, field), f"lane {field}")
        for field in _PLAN_PATH_FIELDS
    }
    _require(type(value.attempt) is int and 1 <= value.attempt <= 3,
             "lane attempt is outside the frozen three-attempt budget")
    root_fields = (
        "lane_root",
        "canonical_adapter_root", "canonical_baseline_root",
        "canonical_build_root", "variant_adapter_root",
        "variant_baseline_root", "variant_build_root",
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
        attempt=value.attempt,
        repository=paths["repository"],
        verifier=paths["verifier"],
        canonical_adapter_root=paths["canonical_adapter_root"],
        canonical_baseline_root=paths["canonical_baseline_root"],
        canonical_build_root=paths["canonical_build_root"],
        variant_adapter_root=paths["variant_adapter_root"],
        variant_baseline_root=paths["variant_baseline_root"],
        variant_build_root=paths["variant_build_root"],
    )


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
    ) -> dict[str, Any]:
        """Publish and seal one validated terminal plus its exact inventory."""
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
        expected = {item["relative_path"]: item for item in inventory}
        _require(set(expected) == set(retained_files) | {terminal},
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
        terminal_content = canonical_json_bytes(validated)
        for path in sorted(retained_files):
            self.publish_bytes(path, retained_files[path])
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


__all__ = (
    "AcquisitionLocks",
    "AcquisitionError",
    "BUILD_CLOSURE_SCHEMA",
    "CANONICAL_LOCK_PATH",
    "CanonicalFileLock",
    "CommandResult",
    "HostEnvironment",
    "LANE_DIRECTORY_MODE",
    "LANE_FILE_MODE",
    "LanePlan",
    "LaneWriter",
    "StableLeaseAnchor",
    "TREE_METADATA_SCHEMA",
    "build_closure_document",
    "canonical_ldd_text",
    "expected_sha256sums",
    "expected_tree_metadata",
    "normalize_ldd_output",
    "validate_lane_plan",
)
