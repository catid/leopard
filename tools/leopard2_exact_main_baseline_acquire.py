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
from dataclasses import dataclass
import hashlib
import importlib.util
import os
from pathlib import Path
import secrets
import stat
import sys
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
_TREE_METADATA_PATH = "TREE-METADATA.json"
_SHA256SUMS_PATH = "SHA256SUMS"
_TERMINAL_PATHS = frozenset(("baseline-authority.json", "FAILED.json"))


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
    "AcquisitionError",
    "BUILD_CLOSURE_SCHEMA",
    "LANE_DIRECTORY_MODE",
    "LANE_FILE_MODE",
    "LanePlan",
    "LaneWriter",
    "TREE_METADATA_SCHEMA",
    "build_closure_document",
    "canonical_ldd_text",
    "expected_sha256sums",
    "expected_tree_metadata",
    "validate_lane_plan",
)
