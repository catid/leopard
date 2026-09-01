#!/usr/bin/env python3
"""Acquire and seal one K65 generation-3 candidate build authority.

The controller captures a clean detached source, canonical Git archive, pure
AVX2 production build provenance, clean reproducible rebuild proof, and the
current runtime-file closure.  It invokes only bounded loader discovery for the
ELF and never enters ``bench_leopard2`` or any timing, qualification, or
campaign mode.  A lane is published only after a fresh isolated verifier has
replayed its complete sealed tree.
"""

from __future__ import annotations

import argparse
import copy
import ctypes
import errno
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import secrets
import shutil
import stat
import struct
import sys
import tempfile
from typing import Any, Mapping, Sequence


_REPOSITORY_PYCACHE_PREFIX = "/dev/null"
if sys.pycache_prefix not in (None, _REPOSITORY_PYCACHE_PREFIX):
    raise RuntimeError("candidate producer Python cache prefix is unsafe")
sys.pycache_prefix = _REPOSITORY_PYCACHE_PREFIX
sys.dont_write_bytecode = True


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
TOOLS = REPO_ROOT / "tools"
MAIN_COMPARE = HERE.parent / "main_compare"


def _inside_repository(path: Path) -> bool:
    try:
        path.relative_to(REPO_ROOT)
        return True
    except ValueError:
        return False


_initial_search_path = list(sys.path)
_safe_search_path: list[str] = []
for _entry in sys.path:
    try:
        _resolved_entry = Path(_entry or os.getcwd()).resolve(strict=True)
    except OSError:
        _safe_search_path.append(_entry)
        continue
    if not _inside_repository(_resolved_entry):
        _safe_search_path.append(_entry)
sys.path[:] = _safe_search_path
del _entry, _resolved_entry, _safe_search_path


def _load_repository_module(
    module_name: str, expected_path: Path, label: str,
) -> Any:
    """Execute one repository module from its exact canonical pathname."""
    expected = expected_path.resolve(strict=True)
    if expected != expected_path or not expected.is_file():
        raise RuntimeError(f"{label} path is not one canonical file")
    specification = importlib.util.spec_from_file_location(
        module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {label} from its exact path")
    module = importlib.util.module_from_spec(specification)
    previous = sys.modules.get(module_name)
    original_search_path = list(sys.path)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
        observed = Path(
            str(getattr(module, "__file__", ""))).resolve(strict=True)
        if observed != expected:
            raise RuntimeError(f"{label} resolved outside this source tree")
    except BaseException:
        if previous is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    finally:
        sys.path[:] = original_search_path
    return module


authority = _load_repository_module(
    "leopard2_k65_gen3_candidate_authority_for_acquisition",
    HERE / "k65_gen3_candidate_authority.py",
    "candidate-authority module")
# The authority loader has already installed and retained exact objects for
# these shared dependencies.  Reuse them so producer and verifier-side logic
# do not diverge merely because the same source file was executed twice.
build_provenance = authority.build_provenance
git_capture = authority.git_capture
sealed_acquire = _load_repository_module(
    "leopard2_k65_gen3_exact_main_acquire",
    TOOLS / "leopard2_exact_main_baseline_acquire.py",
    "sealed-lane acquisition module")
if __name__ != "__main__":
    sys.path[:] = _initial_search_path
del _initial_search_path


contract = authority.contract
RECEIPT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-candidate-receipt/v2"
VERIFICATION_PROGRAM_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-verification-program/v1"
DEFAULT_CANONICAL_LOCK = "/tmp/leopard-gf8-authoritative.lock"
MAX_VERIFIER_OUTPUT_BYTES = 4 * 1024 * 1024
MAX_LDD_OUTPUT_BYTES = 4 * 1024 * 1024
MAX_PROGRAM_FILE_BYTES = 64 * 1024 * 1024
VERIFIER_TIMEOUT_SECONDS = 900
LDD_TIMEOUT_SECONDS = 120
_HEX40 = re.compile(r"[0-9a-f]{40}")
_IN_MODIFY = 0x00000002
_IN_ATTRIB = 0x00000004
_IN_CLOSE_WRITE = 0x00000008
_IN_MOVED_FROM = 0x00000040
_IN_MOVED_TO = 0x00000080
_IN_CREATE = 0x00000100
_IN_DELETE = 0x00000200
_IN_DELETE_SELF = 0x00000400
_IN_MOVE_SELF = 0x00000800
_IN_UNMOUNT = 0x00002000
_IN_Q_OVERFLOW = 0x00004000
_IN_IGNORED = 0x00008000
_IN_MASK_ADD = 0x20000000
_PROGRAM_MUTATION_MASK = (
    _IN_MODIFY | _IN_ATTRIB | _IN_CLOSE_WRITE | _IN_MOVED_FROM |
    _IN_MOVED_TO | _IN_CREATE | _IN_DELETE | _IN_DELETE_SELF |
    _IN_MOVE_SELF | _IN_UNMOUNT | _IN_Q_OVERFLOW)


class CandidateAcquisitionError(RuntimeError):
    """Candidate authority acquisition failed closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CandidateAcquisitionError(message)


def _canonical_directory(
    value: Path | str, label: str, *, private: bool = False,
) -> Path:
    lexical = Path(os.path.abspath(os.fspath(value)))
    try:
        resolved = lexical.resolve(strict=True)
        status = os.lstat(lexical)
    except OSError as error:
        raise CandidateAcquisitionError(f"{label} is unavailable") from error
    mode = stat.S_IMODE(status.st_mode)
    _require(
        lexical == resolved and stat.S_ISDIR(status.st_mode) and
        not stat.S_ISLNK(status.st_mode) and
        status.st_uid == os.geteuid() and status.st_gid == os.getegid() and
        (mode == 0o700 if private else mode & 0o022 == 0),
        f"{label} is not one safe owned canonical directory")
    return resolved


def _absent_output(value: Path | str) -> tuple[Path, Path]:
    output = Path(os.path.abspath(os.fspath(value)))
    _require(output.name not in ("", ".", "..") and
             output.parent != output,
             "candidate authority output path is invalid")
    parent = _canonical_directory(
        output.parent, "candidate authority output parent", private=True)
    _require(output.parent == parent and not output.exists() and
             not output.is_symlink(),
             "candidate authority output lane is not absent")
    return output, parent


def _path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _hash_descriptor(descriptor: int, size: int, label: str) -> str:
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        chunk = os.pread(descriptor, min(1 << 20, size - offset), offset)
        _require(bool(chunk), f"{label} ended while hashed")
        digest.update(chunk)
        offset += len(chunk)
    return digest.hexdigest()


def _program_file_identity(
    role: str, path: Path | str, *, executable: bool = False,
) -> dict[str, Any]:
    """Hash one canonical acquisition/verifier trust input by descriptor."""
    _require(type(role) is str and bool(role),
             "verification program role is invalid")
    lexical = Path(os.path.abspath(os.fspath(path)))
    descriptor = -1
    try:
        resolved = lexical.resolve(strict=True)
        _require(lexical == resolved,
                 f"verification program {role} path is not canonical")
        descriptor = os.open(
            resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        before = os.fstat(descriptor)
        pathname = os.stat(resolved, follow_symlinks=False)
        mode = stat.S_IMODE(before.st_mode)
        fingerprint = (
            before.st_dev, before.st_ino, before.st_mode, before.st_uid,
            before.st_gid, before.st_nlink, before.st_size,
            before.st_mtime_ns, before.st_ctime_ns,
        )
        _require(
            stat.S_ISREG(before.st_mode) and
            0 < before.st_size <= MAX_PROGRAM_FILE_BYTES and
            (not executable or mode & 0o111 != 0) and
            (pathname.st_dev, pathname.st_ino, pathname.st_mode,
             pathname.st_uid, pathname.st_gid, pathname.st_nlink,
             pathname.st_size, pathname.st_mtime_ns,
             pathname.st_ctime_ns) == fingerprint,
            f"verification program {role} is not one stable file")
        digest = _hash_descriptor(descriptor, before.st_size, role)
        after = os.fstat(descriptor)
        pathname_after = os.stat(resolved, follow_symlinks=False)
        _require(
            (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
             after.st_gid, after.st_nlink, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns) == fingerprint and
            (pathname_after.st_dev, pathname_after.st_ino,
             pathname_after.st_mode, pathname_after.st_uid,
             pathname_after.st_gid, pathname_after.st_nlink,
             pathname_after.st_size, pathname_after.st_mtime_ns,
             pathname_after.st_ctime_ns) == fingerprint,
            f"verification program {role} changed while hashed")
        return {
            "role": role, "path": str(resolved), "size": before.st_size,
            "mode": mode, "sha256": digest,
        }
    except CandidateAcquisitionError:
        raise
    except OSError as error:
        raise CandidateAcquisitionError(
            f"verification program {role} could not be captured") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def verification_program_identity() -> dict[str, Any]:
    """Bind the producer plus the isolated verifier's complete repo closure."""
    paths = {
        "acquisition_controller":
            HERE / "k65_gen3_candidate_authority_acquire.py",
        "authority_verifier": HERE / "k65_gen3_candidate_authority.py",
        "preregistration_contract": HERE / "k65_gen3_preregistration.py",
        "pair_contract": MAIN_COMPARE / "pair_qualification_contract.py",
        "pair_bridge_contract":
            MAIN_COMPARE / "pair_qualification_bridge_contract.py",
        "git_capture": MAIN_COMPARE / "git_capture.py",
        "build_provenance": TOOLS / "leopard2_build_provenance.py",
        "sealed_lane_acquisition":
            TOOLS / "leopard2_exact_main_baseline_acquire.py",
        "sealed_tree_verifier":
            TOOLS / "leopard2_exact_main_baseline_verifier.py",
        "exact_main_identity_contract":
            TOOLS / "leopard2_exact_main_baseline.py",
        "exact_main_record_contract":
            TOOLS / "leopard2_exact_main_baseline_record.py",
        "python": Path(os.path.realpath(sys.executable)),
    }
    files = [
        _program_file_identity(
            role, paths[role], executable=(role == "python"))
        for role in sorted(paths)
    ]
    _require(len({item["role"] for item in files}) == len(files) and
             len({item["path"] for item in files}) == len(files),
             "verification program closure is not unique")
    return {"schema": VERIFICATION_PROGRAM_SCHEMA, "files": files}


def _program_path(program: Mapping[str, Any], role: str) -> str:
    _require(type(program) is dict and
             program.get("schema") == VERIFICATION_PROGRAM_SCHEMA and
             type(program.get("files")) is list,
             "verification program identity is invalid")
    matches = [
        item for item in program["files"]
        if type(item) is dict and item.get("role") == role
    ]
    _require(len(matches) == 1 and type(matches[0].get("path")) is str,
             f"verification program has no unique {role}")
    return matches[0]["path"]


class VerificationProgramGuard:
    """Retain every verifier input and detect even restored path mutation."""

    def __init__(self, program: Mapping[str, Any]) -> None:
        self.program = copy.deepcopy(dict(program))
        self.inotify_descriptor = -1
        self.files: dict[str, dict[str, Any]] = {}
        self.directory_targets: dict[int, set[str]] = {}
        self.file_watches: set[int] = set()

    def _add_watch(self, path: Path) -> int:
        add = ctypes.CDLL(None, use_errno=True).inotify_add_watch
        add.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32)
        add.restype = ctypes.c_int
        watch = add(
            self.inotify_descriptor, os.fsencode(path),
            _PROGRAM_MUTATION_MASK | _IN_MASK_ADD)
        if watch < 0:
            number = ctypes.get_errno() or errno.EPERM
            raise OSError(number, os.strerror(number), str(path))
        return watch

    def _watch_path(self, path: Path) -> None:
        parent = Path("/")
        for component in path.parts[1:]:
            watch = self._add_watch(parent)
            self.directory_targets.setdefault(watch, set()).add(component)
            parent = parent / component
        self.file_watches.add(self._add_watch(path))

    def _retain(self, record: Mapping[str, Any]) -> None:
        _require(type(record) is dict and set(record) == {
                     "role", "path", "size", "mode", "sha256"},
                 "verification program file record is invalid")
        role = record["role"]
        path = Path(str(record["path"]))
        _require(type(role) is str and role not in self.files and
                 path.is_absolute() and path.resolve(strict=True) == path,
                 "verification program retained path is invalid")
        self._watch_path(path)
        descriptor = -1
        try:
            descriptor = os.open(
                path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
            before = os.fstat(descriptor)
            pathname = os.stat(path, follow_symlinks=False)
            fingerprint = (
                before.st_dev, before.st_ino, before.st_mode, before.st_uid,
                before.st_gid, before.st_nlink, before.st_size,
                before.st_mtime_ns, before.st_ctime_ns,
            )
            _require(
                stat.S_ISREG(before.st_mode) and
                before.st_size == record["size"] and
                stat.S_IMODE(before.st_mode) == record["mode"] and
                _hash_descriptor(descriptor, before.st_size, role) ==
                    record["sha256"] and
                (pathname.st_dev, pathname.st_ino, pathname.st_mode,
                 pathname.st_uid, pathname.st_gid, pathname.st_nlink,
                 pathname.st_size, pathname.st_mtime_ns,
                 pathname.st_ctime_ns) == fingerprint,
                f"verification program {role} differs while retained")
            self.files[role] = {
                "descriptor": descriptor, "path": path,
                "fingerprint": fingerprint,
            }
            descriptor = -1
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    def __enter__(self) -> "VerificationProgramGuard":
        _require(
            self.inotify_descriptor < 0 and not self.files and
            self.program.get("schema") == VERIFICATION_PROGRAM_SCHEMA and
            type(self.program.get("files")) is list,
            "verification program guard input is invalid")
        initialize = ctypes.CDLL(None, use_errno=True).inotify_init1
        initialize.argtypes = (ctypes.c_int,)
        initialize.restype = ctypes.c_int
        self.inotify_descriptor = initialize(os.O_NONBLOCK | os.O_CLOEXEC)
        if self.inotify_descriptor < 0:
            number = ctypes.get_errno() or errno.EPERM
            raise CandidateAcquisitionError(
                "verification program mutation guard could not start") \
                from OSError(number, os.strerror(number))
        try:
            for record in self.program["files"]:
                self._retain(record)
            _require(len(self.files) == len(self.program["files"]),
                     "verification program retained closure is not unique")
            self.validate_current()
            return self
        except BaseException:
            self.close()
            raise

    def descriptor(self, role: str) -> int:
        _require(role in self.files,
                 f"verification program does not retain {role}")
        return self.files[role]["descriptor"]

    def _reject_events(self) -> None:
        while True:
            try:
                data = os.read(self.inotify_descriptor, 64 * 1024)
            except BlockingIOError:
                return
            except OSError as error:
                raise CandidateAcquisitionError(
                    "verification program mutation queue could not be read") \
                    from error
            _require(bool(data),
                     "verification program mutation queue closed")
            offset = 0
            while offset < len(data):
                _require(offset + 16 <= len(data),
                         "verification program mutation event is truncated")
                watch, mask, unused_cookie, name_size = struct.unpack_from(
                    "iIII", data, offset)
                end = offset + 16 + name_size
                _require(end <= len(data),
                         "verification program mutation name is truncated")
                name = data[offset + 16:end].split(b"\0", 1)[0]
                offset = end
                _require(mask & (_IN_Q_OVERFLOW | _IN_IGNORED) == 0,
                         "verification program mutation watch lost authority")
                if watch in self.file_watches:
                    raise CandidateAcquisitionError(
                        "retained verification program file was mutated")
                targets = self.directory_targets.get(watch)
                _require(targets is not None,
                         "verification program mutation watch is unknown")
                decoded = os.fsdecode(name)
                if not decoded or decoded in targets:
                    raise CandidateAcquisitionError(
                        "retained verification program path was mutated")

    def validate_current(self) -> None:
        _require(self.inotify_descriptor >= 0 and bool(self.files),
                 "verification program guard is closed")
        for role, retained in self.files.items():
            descriptor = retained["descriptor"]
            path = retained["path"]
            expected = retained["fingerprint"]
            try:
                current = os.fstat(descriptor)
                pathname = os.stat(path, follow_symlinks=False)
                resolved = path.resolve(strict=True)
            except OSError as error:
                raise CandidateAcquisitionError(
                    f"verification program {role} is unavailable") from error
            _require(
                (current.st_dev, current.st_ino, current.st_mode,
                 current.st_uid, current.st_gid, current.st_nlink,
                 current.st_size, current.st_mtime_ns,
                 current.st_ctime_ns) == expected and
                (pathname.st_dev, pathname.st_ino, pathname.st_mode,
                 pathname.st_uid, pathname.st_gid, pathname.st_nlink,
                 pathname.st_size, pathname.st_mtime_ns,
                 pathname.st_ctime_ns) == expected and resolved == path,
                f"verification program {role} identity changed")
        self._reject_events()

    def close(self) -> None:
        for retained in self.files.values():
            try:
                os.close(retained["descriptor"])
            except OSError:
                pass
        self.files.clear()
        if self.inotify_descriptor >= 0:
            descriptor = self.inotify_descriptor
            self.inotify_descriptor = -1
            try:
                os.close(descriptor)
            except OSError:
                pass

    def __exit__(self, unused_kind: object, unused_value: object,
                 unused_traceback: object) -> None:
        self.close()


def _runtime_file_identity(path: str, role: str) -> dict[str, Any]:
    """Bind a raw loader/launcher path to its resolved current regular file."""
    raw = Path(path)
    _require(raw.is_absolute() and "\0" not in path and
             type(role) is str and bool(role),
             "runtime-file request is malformed")
    descriptor = -1
    try:
        resolved = raw.resolve(strict=True)
        descriptor = os.open(
            resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        before = os.fstat(descriptor)
        pathname = os.stat(resolved, follow_symlinks=False)
        _require(
            stat.S_ISREG(before.st_mode) and
            0 < before.st_size <= authority.MAX_EXECUTABLE_BYTES and
            (before.st_dev, before.st_ino) ==
                (pathname.st_dev, pathname.st_ino),
            f"runtime {role} {path!r} is not one bounded regular file")
        digest = _hash_descriptor(descriptor, before.st_size, path)
        after = os.fstat(descriptor)
        pathname_after = os.stat(resolved, follow_symlinks=False)
        fingerprint = (
            before.st_dev, before.st_ino, before.st_mode, before.st_uid,
            before.st_gid, before.st_nlink, before.st_size,
            before.st_mtime_ns, before.st_ctime_ns,
        )
        _require(
            (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
             after.st_gid, after.st_nlink, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns) ==
                fingerprint and
            (pathname_after.st_dev, pathname_after.st_ino,
             pathname_after.st_mode, pathname_after.st_uid,
             pathname_after.st_gid, pathname_after.st_nlink,
             pathname_after.st_size, pathname_after.st_mtime_ns,
             pathname_after.st_ctime_ns) == fingerprint and
            raw.resolve(strict=True) == resolved,
            f"runtime {role} {path!r} changed while captured")
        return {
            "path": path,
            "sha256": digest,
            "size": before.st_size,
            "mode": stat.S_IMODE(before.st_mode),
            "role": role,
        }
    except CandidateAcquisitionError:
        raise
    except OSError as error:
        raise CandidateAcquisitionError(
            f"runtime {role} {path!r} could not be captured") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def capture_runtime_closure(
    environment: sealed_acquire.HostEnvironment,
    executable: Path,
    scratch_root: Path,
) -> dict[str, Any]:
    """Discover the exact file set twice and bind every current target."""
    child_environment = sealed_acquire.frozen_child_environment()
    ldd_before = _runtime_file_identity(authority.LDD_PATH, "discovery")
    interpreter_before = _runtime_file_identity(
        authority.LDD_INTERPRETER_PATH, authority.LDD_INTERPRETER_ROLE)

    def discover() -> tuple[dict[str, Any], ...]:
        result = environment.run(
            [authority.LDD_PATH, str(executable)], cwd=str(scratch_root),
            env=child_environment, timeout=LDD_TIMEOUT_SECONDS,
            maximum_bytes=MAX_LDD_OUTPUT_BYTES)
        _require(result.exit_status == 0 and result.stderr == b"",
                 "candidate runtime discovery failed or emitted diagnostics")
        try:
            return authority.normalize_runtime_discovery(result.stdout)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate runtime discovery output is invalid") from error

    first_rows = discover()
    second_rows = discover()
    _require(contract.exact_json_equal(list(first_rows), list(second_rows)),
             "candidate runtime dependency set changed between discoveries")
    discovery_sha256 = authority.runtime_discovery_sha256(first_rows)
    dependency_paths = sorted(
        row["path"] for row in first_rows if row["kind"] == "file")
    _require(len(dependency_paths) == len(set(dependency_paths)),
             "candidate runtime dependency paths are duplicated")
    dependencies = [
        _runtime_file_identity(path, "dependency")
        for path in dependency_paths
    ]
    ldd_after = _runtime_file_identity(authority.LDD_PATH, "discovery")
    interpreter_after = _runtime_file_identity(
        authority.LDD_INTERPRETER_PATH, authority.LDD_INTERPRETER_ROLE)
    _require(
        contract.exact_json_equal(ldd_before, ldd_after) and
        contract.exact_json_equal(interpreter_before, interpreter_after),
        "candidate runtime discovery program changed during capture")
    ldd_before["role"] = authority.LDD_DISCOVERY_ROLE_PREFIX + discovery_sha256
    dependencies.extend((interpreter_before, ldd_before))
    dependencies.sort(key=lambda item: (item["role"], item["path"]))
    launchers = [
        _runtime_file_identity(path, "launcher")
        for path in ("/usr/bin/prlimit", "/usr/bin/taskset")
    ]
    # Establish that no current runtime target changed between the first and
    # final byte bindings used by the sealed record.
    replayed = [
        _runtime_file_identity(item["path"], item["role"])
        for item in dependencies + launchers
    ]
    _require(contract.exact_json_equal(dependencies + launchers, replayed),
             "candidate runtime files changed during closure capture")
    return authority.validate_runtime_closure({
        "schema": authority.RUNTIME_CLOSURE_SCHEMA,
        "dependencies": dependencies,
        "launchers": launchers,
    })


def _owned_payload_identity(
    path: Path, label: str, *, maximum_bytes: int,
) -> dict[str, Any]:
    try:
        identity = sealed_acquire._owned_file_identity(
            str(path), label, maximum_bytes=maximum_bytes)
    except Exception as error:
        raise CandidateAcquisitionError(
            f"{label} is not one stable owner-only payload") from error
    return {"size": identity["size"], "sha256": identity["sha256"]}


def _canonical_bytes(value: Any, label: str) -> bytes:
    try:
        data = contract.canonical_json_bytes(value)
    except Exception as error:
        raise CandidateAcquisitionError(f"{label} is not canonical JSON") \
            from error
    _require(0 < len(data) <= authority.MAX_AUTHORITY_PAYLOAD_BYTES,
             f"{label} exceeds its payload bound")
    return data


def _canonical_archive_pair(
    environment: sealed_acquire.HostEnvironment,
    source_root: Path,
    commit: str,
    scratch_root: Path,
) -> dict[str, Any]:
    child_environment = sealed_acquire.frozen_child_environment()
    commands: list[dict[str, Any]] = []
    try:
        first = sealed_acquire.canonical_git_archive(
            environment, source_repository=str(source_root), commit=commit,
            prefix="candidate-source/", scratch_root=str(scratch_root),
            destination_name="candidate-source-first.tar",
            child_environment=child_environment, log=commands)
        second = sealed_acquire.canonical_git_archive(
            environment, source_repository=str(source_root), commit=commit,
            prefix="candidate-source/", scratch_root=str(scratch_root),
            destination_name="candidate-source-second.tar",
            child_environment=child_environment, log=commands)
        _require(
            (first["size"], first["sha256"]) ==
                (second["size"], second["sha256"]) and
            sealed_acquire._owned_files_identical(
                first["path"], second["path"],
                "candidate canonical source archive"),
            "candidate canonical source archive did not reproduce")
        _require(first["size"] <= authority.MAX_AUTHORITY_PAYLOAD_BYTES,
                 "candidate canonical source archive exceeds its bound")
        return first
    except CandidateAcquisitionError:
        raise
    except Exception as error:
        raise CandidateAcquisitionError(
            "candidate canonical source archive could not be produced") \
            from error


def _rename_directory_noreplace(staging: Path, output: Path) -> None:
    _require(staging.parent == output.parent,
             "candidate lane publication crosses a directory boundary")
    parent_descriptor = -1
    try:
        before = os.lstat(staging)
        _require(stat.S_ISDIR(before.st_mode) and not output.exists() and
                 not output.is_symlink(),
                 "candidate lane publication paths differ")
        parent_descriptor = os.open(
            staging.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
        libc = ctypes.CDLL(None, use_errno=True)
        renameat2 = getattr(libc, "renameat2", None)
        _require(renameat2 is not None,
                 "candidate lane publication requires Linux renameat2")
        renameat2.argtypes = (
            ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_uint)
        renameat2.restype = ctypes.c_int
        ctypes.set_errno(0)
        status = renameat2(
            parent_descriptor, os.fsencode(staging.name), parent_descriptor,
            os.fsencode(output.name), 1)
        if status != 0:
            number = ctypes.get_errno() or errno.EPERM
            raise OSError(number, os.strerror(number), str(output))
        after = os.lstat(output)
        _require(not staging.exists() and not staging.is_symlink() and
                 (after.st_dev, after.st_ino) ==
                    (before.st_dev, before.st_ino),
                 "candidate lane publication identity differs")
        os.fsync(parent_descriptor)
    except CandidateAcquisitionError:
        raise
    except OSError as error:
        raise CandidateAcquisitionError(
            "candidate lane could not be published without replacement") \
            from error
    finally:
        if parent_descriptor >= 0:
            os.close(parent_descriptor)


def _fresh_verdict(
    environment: sealed_acquire.HostEnvironment,
    *,
    lane_root: Path,
    source_root: Path,
    controller_sha256: str,
    record: Mapping[str, Any],
    record_sha256: str,
    ledger_sha256: str,
    expected_program: Mapping[str, Any] | None = None,
    expected_runtime: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    program = verification_program_identity() if expected_program is None \
        else copy.deepcopy(dict(expected_program))
    with VerificationProgramGuard(program) as guard:
        before = verification_program_identity()
        _require(contract.exact_json_equal(before, program),
                 "verification program differs before isolated replay")
        guard.validate_current()
        interpreter_descriptor = guard.descriptor("python")
        verifier_descriptor = guard.descriptor("authority_verifier")
        argv = [
            f"/proc/self/fd/{interpreter_descriptor}",
            "-I", "-S", "-B", "-X", "pycache_prefix=/dev/null",
            f"/proc/self/fd/{verifier_descriptor}",
            str(lane_root), "--source-root", str(source_root),
            "--controller-bindings-sha256", controller_sha256,
        ]
        try:
            result = environment.run(
                argv, cwd=str(source_root),
                env=sealed_acquire.frozen_child_environment(),
                timeout=VERIFIER_TIMEOUT_SECONDS,
                maximum_bytes=MAX_VERIFIER_OUTPUT_BYTES,
                inherited_descriptors=(
                    interpreter_descriptor, verifier_descriptor),
                executable_descriptor=interpreter_descriptor)
        except BaseException as error:
            raise CandidateAcquisitionError(
                "independent candidate verifier could not run") from error
        finally:
            after = verification_program_identity()
            guard.validate_current()
            _require(contract.exact_json_equal(after, program),
                     "verification program changed during isolated replay")
    _require(result.exit_status == 0 and result.stderr == b"" and
             bool(result.stdout),
             "independent candidate verifier failed or emitted diagnostics")
    try:
        value = contract.strict_json_loads(
            result.stdout, "independent candidate verdict")
        _require(type(value) is dict and
                 result.stdout == contract.canonical_json_bytes(value),
                 "independent candidate verdict bytes are not canonical")
        verdict = authority.validate_candidate_authority_verdict(
            value, expected_root=lane_root,
            expected_record_sha256=record_sha256,
            expected_ledger_sha256=ledger_sha256,
            expected_controller_bindings_sha256=controller_sha256)
    except CandidateAcquisitionError:
        raise
    except Exception as error:
        raise CandidateAcquisitionError(
            "independent candidate verdict is invalid") from error
    _require(
        contract.exact_json_equal(verdict["source"], {
            "commit": record["source"]["commit"],
            "tree": record["source"]["tree"], "detached": True}) and
        contract.exact_json_equal(verdict["candidate"], record["candidate"]) and
        contract.exact_json_equal(
            verdict["build_closure"], record["build_closure"]),
        "independent candidate verdict differs from its terminal record")
    if expected_runtime is not None:
        try:
            runtime = authority.validate_runtime_closure(expected_runtime)
        except Exception as error:
            raise CandidateAcquisitionError(
                "expected candidate runtime closure is invalid") from error
        discovery = next(
            item for item in runtime["dependencies"]
            if item["role"].startswith(authority.LDD_DISCOVERY_ROLE_PREFIX))
        interpreter_entry = next(
            item for item in runtime["dependencies"]
            if item["role"] == authority.LDD_INTERPRETER_ROLE)
        view = verdict["runtime_discovery"]
        _require(
            view["transcript_sha256"] == discovery["role"][
                len(authority.LDD_DISCOVERY_ROLE_PREFIX):] and
            view["ldd_sha256"] == discovery["sha256"] and
            view["interpreter_sha256"] == interpreter_entry["sha256"],
            "independent runtime discovery differs from the sealed closure")
    return verdict


def _candidate_receipt(
    verdict: Mapping[str, Any], program: Mapping[str, Any],
) -> dict[str, Any]:
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "root": verdict["root"],
        "record_sha256": verdict["record"]["sha256"],
        "authority_ledger_sha256": verdict["seal"]["sha256sums_sha256"],
        "source": copy.deepcopy(verdict["source"]),
        "candidate": copy.deepcopy(verdict["candidate"]),
        "build_closure": copy.deepcopy(verdict["build_closure"]),
        "controller_bindings_sha256":
            verdict["controller_bindings_sha256"],
        "verification_program": copy.deepcopy(dict(program)),
        "verification_program_sha256": contract.canonical_sha256(program),
        "runtime_discovery": copy.deepcopy(verdict["runtime_discovery"]),
        "verifier": copy.deepcopy(dict(verdict)),
    }
    encoded = _canonical_bytes(receipt, "candidate acquisition receipt")
    _require(len(encoded) <= MAX_VERIFIER_OUTPUT_BYTES,
             "candidate acquisition receipt exceeds its output bound")
    return receipt


def _publish_verified_receipt(
    environment: sealed_acquire.HostEnvironment,
    *,
    staging: Path,
    output_lane: Path,
    source_root: Path,
    controller_sha256: str,
    record: Mapping[str, Any],
    record_sha256: str,
    ledger_sha256: str,
    program: Mapping[str, Any],
    runtime: Mapping[str, Any],
    lock: sealed_acquire.CanonicalFileLock,
) -> dict[str, Any]:
    """Commit, replay, and either return a receipt or restore quarantine."""
    lock.validate_current()
    try:
        _rename_directory_noreplace(staging, output_lane)
        verdict = _fresh_verdict(
            environment, lane_root=output_lane, source_root=source_root,
            controller_sha256=controller_sha256, record=record,
            record_sha256=record_sha256, ledger_sha256=ledger_sha256,
            expected_program=program, expected_runtime=runtime)
        lock.validate_current()
        return _candidate_receipt(verdict, program)
    except BaseException as error:
        staging_present = os.path.lexists(staging)
        output_present = os.path.lexists(output_lane)
        if output_present and not staging_present:
            try:
                _rename_directory_noreplace(output_lane, staging)
            except BaseException as rollback_error:
                # renameat2 may have committed even when the following fsync
                # failed.  Reinspect before declaring the official path live.
                staging_present = os.path.lexists(staging)
                output_present = os.path.lexists(output_lane)
                if staging_present and not output_present:
                    parent_descriptor = -1
                    try:
                        parent_descriptor = os.open(
                            staging.parent,
                            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                            os.O_CLOEXEC)
                        os.fsync(parent_descriptor)
                    except OSError as sync_error:
                        raise CandidateAcquisitionError(
                            f"candidate authority official lane is absent but "
                            f"rollback durability failed: {sync_error}") \
                            from error
                    finally:
                        if parent_descriptor >= 0:
                            os.close(parent_descriptor)
                else:
                    raise CandidateAcquisitionError(
                        f"candidate authority publication failed after "
                        f"commit; official lane {output_lane} may remain "
                        f"occupied because rollback failed: "
                        f"{rollback_error}") from error
        elif not (staging_present and not output_present):
            raise CandidateAcquisitionError(
                f"candidate authority publication state is ambiguous: "
                f"staging={staging_present}, official={output_present}") \
                from error
        raise CandidateAcquisitionError(
            f"candidate authority publication was rolled back to retained "
            f"staging lane {staging}") from error


def _acquire_locked(
    *,
    source_root: Path,
    build_root: Path,
    output_lane: Path,
    output_parent: Path,
    expected_commit: str,
    expected_tree: str,
    jobs: int,
    environment: sealed_acquire.HostEnvironment,
    lock: sealed_acquire.CanonicalFileLock,
) -> dict[str, Any]:
    scratch = Path(tempfile.mkdtemp(
        prefix=".k65-gen3-candidate-work-", dir=output_parent))
    os.chmod(scratch, 0o700)
    staging = output_parent / (
        "." + output_lane.name + ".staging-" + secrets.token_hex(16))
    staging_created = False
    body_failed = False
    scratch_cleaned = False
    try:
        inherited = (lock.descriptor,)
        lock.validate_current()
        try:
            first_git = git_capture.capture_git_identity(
                source_root, expected_commit, require_detached=True,
                inherited_descriptors=inherited)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate source Git identity could not be captured") \
                from error
        _require(first_git["tree"] == expected_tree and
                 first_git["tracked_status"] == "clean",
                 "candidate source commit/tree/cleanliness differs")
        controller_sha256 = authority.controller_bindings_sha256(source_root)
        program = verification_program_identity()
        lock.validate_current()

        executable = build_root / "bench_leopard2"
        try:
            provenance = build_provenance.candidate_build_provenance(
                build_root, source_root, executable, "bench_leopard2",
                inherited_descriptors=inherited)
            proof = build_provenance.verify_reproducible_candidate_build(
                provenance, jobs=jobs, inherited_descriptors=inherited)
            build_provenance.validate_reproducible_build_proof(
                proof, provenance, label="K65 candidate")
            core = build_provenance.compare_reproducible_builds(
                provenance, provenance)
            header = build_provenance.\
                _canonical_replay_attestation_header_bytes(provenance)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate pure-AVX2 reproducible build did not validate") \
                from error
        provenance_data = _canonical_bytes(
            provenance, "candidate build provenance")
        proof_data = _canonical_bytes(
            proof, "candidate reproducible-build proof")
        core_data = _canonical_bytes(
            core, "candidate reproducible-build core")
        _require(type(header) is bytes and 0 < len(header) <=
                 authority.MAX_AUTHORITY_PAYLOAD_BYTES,
                 "candidate attestation header is invalid")
        staged_executable = scratch / "bench_leopard2"
        try:
            staged_identity = sealed_acquire.stage_build_output(
                str(executable), str(staged_executable),
                "K65 candidate executable",
                maximum_bytes=authority.MAX_EXECUTABLE_BYTES)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate executable could not be staged") from error
        _require(
            staged_identity["size"] == provenance["executable"]["size"] and
            staged_identity["sha256"] ==
                provenance["executable"]["sha256"],
            "staged candidate executable differs from build provenance")

        archive = _canonical_archive_pair(
            environment, source_root, expected_commit, scratch)
        runtime = capture_runtime_closure(
            environment, staged_executable, scratch)
        runtime_data = _canonical_bytes(runtime, "candidate runtime closure")
        git_data = _canonical_bytes(first_git, "candidate Git capture")
        lock.validate_current()
        try:
            final_git = git_capture.capture_git_identity(
                source_root, expected_commit, require_detached=True,
                inherited_descriptors=inherited)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate source Git identity could not be recaptured") \
                from error
        _require(contract.exact_json_equal(first_git, final_git),
                 "candidate source changed during authority acquisition")

        retained_files = {
            "build/benchmark-source-attestation.h": header,
            "build/build-provenance.json": provenance_data,
            "build/reproducible-build-core.json": core_data,
            "build/reproducible-build-proof.json": proof_data,
            "runtime/runtime-closure.json": runtime_data,
            "source/git-capture.json": git_data,
        }
        retained_paths = {
            "artifacts/bench_leopard2": str(staged_executable),
            "source/candidate-source.tar": archive["path"],
        }
        identities: dict[str, dict[str, Any]] = {
            path: authority.payload_identity(data)
            for path, data in retained_files.items()
        }
        identities["artifacts/bench_leopard2"] = _owned_payload_identity(
            staged_executable, "staged candidate executable",
            maximum_bytes=authority.MAX_EXECUTABLE_BYTES)
        identities["source/candidate-source.tar"] = _owned_payload_identity(
            Path(archive["path"]), "candidate source archive",
            maximum_bytes=authority.MAX_AUTHORITY_PAYLOAD_BYTES)
        record = authority.candidate_authority_record(
            source_commit=expected_commit, source_tree=expected_tree,
            controller_bindings_sha256=controller_sha256,
            payload_identities=identities)
        record_data = _canonical_bytes(record, "candidate authority record")
        record_sha256 = hashlib.sha256(record_data).hexdigest()

        lock.validate_current()
        try:
            with sealed_acquire.LaneWriter(str(staging)) as writer:
                staging_created = True
                seal = writer.seal_payload(
                    terminal=authority.TERMINAL_PATH,
                    terminal_content=record_data,
                    retained_files=retained_files,
                    retained_paths=retained_paths)
        except Exception as error:
            raise CandidateAcquisitionError(
                "candidate authority lane could not be sealed") from error
        _require(seal["terminal_sha256"] == record_sha256,
                 "candidate authority sealer returned another terminal hash")
        _fresh_verdict(
            environment, lane_root=staging, source_root=source_root,
            controller_sha256=controller_sha256, record=record,
            record_sha256=record_sha256,
            ledger_sha256=seal["sha256sums_sha256"],
            expected_program=program, expected_runtime=runtime)
        lock.validate_current()
        try:
            shutil.rmtree(scratch)
        except OSError as error:
            raise CandidateAcquisitionError(
                "candidate acquisition scratch cleanup failed before "
                "publication") from error
        scratch_cleaned = True
        lock.validate_current()
        return _publish_verified_receipt(
            environment, staging=staging, output_lane=output_lane,
            source_root=source_root, controller_sha256=controller_sha256,
            record=record, record_sha256=record_sha256,
            ledger_sha256=seal["sha256sums_sha256"], program=program,
            runtime=runtime, lock=lock)
    except BaseException as error:
        body_failed = True
        if staging_created:
            if output_lane.exists() and not staging.exists():
                message = (f"{error}; publication may have committed at "
                           f"{output_lane} and requires independent replay")
            else:
                message = (f"{error}; fail-closed staging lane retained at "
                           f"{staging}")
            raise CandidateAcquisitionError(message) from error
        raise
    finally:
        if not scratch_cleaned:
            try:
                shutil.rmtree(scratch)
            except FileNotFoundError:
                pass
            except OSError as error:
                if not body_failed:
                    raise CandidateAcquisitionError(
                        "candidate acquisition scratch cleanup failed") \
                        from error


def acquire_candidate_authority(
    *,
    source_root: Path | str,
    build_root: Path | str,
    output_lane: Path | str,
    expected_commit: str,
    expected_tree: str,
    jobs: int,
    canonical_lock_path: Path | str = DEFAULT_CANONICAL_LOCK,
    environment: sealed_acquire.HostEnvironment | None = None,
) -> dict[str, Any]:
    """Run the build-only acquisition and return its verified final receipt."""
    source = _canonical_directory(source_root, "candidate source root")
    build = _canonical_directory(build_root, "candidate build root")
    output, parent = _absent_output(output_lane)
    _require(source == REPO_ROOT and
             not _path_is_within(build, source) and
             not _path_is_within(source, build) and
             not _path_is_within(parent, source),
             "candidate source/build/output topology is unsafe")
    _require(type(expected_commit) is str and
             _HEX40.fullmatch(expected_commit) is not None and
             type(expected_tree) is str and
             _HEX40.fullmatch(expected_tree) is not None,
             "candidate expected commit/tree is invalid")
    _require(type(jobs) is int and 1 <= jobs <= 128,
             "candidate reproducible-build jobs must be in 1..128")
    lock_path = Path(os.path.abspath(os.fspath(canonical_lock_path)))
    host = environment or sealed_acquire.HostEnvironment(
        canonical_lock_path=str(lock_path))
    _require(isinstance(host, sealed_acquire.HostEnvironment),
             "candidate acquisition environment is invalid")
    with sealed_acquire.CanonicalFileLock(str(lock_path)) as lock:
        return _acquire_locked(
            source_root=source, build_root=build, output_lane=output,
            output_parent=parent, expected_commit=expected_commit,
            expected_tree=expected_tree, jobs=jobs, environment=host,
            lock=lock)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="acquire a sealed build-only K65 Gen3 candidate authority")
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--build-root", required=True, type=Path)
    parser.add_argument("--output-lane", required=True, type=Path)
    parser.add_argument("--expected-commit", required=True)
    parser.add_argument("--expected-tree", required=True)
    parser.add_argument("--jobs", required=True, type=int)
    parser.add_argument(
        "--canonical-lock", type=Path,
        default=Path(DEFAULT_CANONICAL_LOCK))
    arguments = parser.parse_args(argv)
    try:
        receipt = acquire_candidate_authority(
            source_root=arguments.source_root,
            build_root=arguments.build_root,
            output_lane=arguments.output_lane,
            expected_commit=arguments.expected_commit,
            expected_tree=arguments.expected_tree,
            jobs=arguments.jobs,
            canonical_lock_path=arguments.canonical_lock)
        sys.stdout.buffer.write(contract.canonical_json_bytes(receipt))
        return 0
    except Exception as error:
        message = " ".join(str(error).replace("\r", " ").replace(
            "\n", " ").split())[:4096]
        sys.stderr.write((message or type(error).__name__) + "\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
