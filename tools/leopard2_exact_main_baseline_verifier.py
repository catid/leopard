#!/usr/bin/env python3
"""Independent read-only verifier for sealed exact-main baseline lanes.

The verifier consumes only retained bytes plus the fixed pure contracts in
``leopard2_exact_main_baseline`` and
``leopard2_exact_main_baseline_record``.  It never imports acquisition code,
executes retained code, launches a subprocess, or writes to the lane.  A lane
may be relocated, but it must remain owned by the invoking effective uid/gid
and sealed with exact owner-only modes (files 0400, directories 0500).
"""

from __future__ import annotations

import base64
import binascii
import copy
import hashlib
import importlib.util
import os
import stat
import sys
from typing import Any, Mapping, NoReturn, Sequence

_TOOLS_DIRECTORY = os.path.dirname(os.path.realpath(__file__))


def _load_local_contract(module_name: str, filename: str) -> Any:
    path = os.path.join(_TOOLS_DIRECTORY, filename)
    existing = sys.modules.get(module_name)
    if existing is not None and os.path.realpath(
            getattr(existing, "__file__", "")) == path:
        return existing
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact-main verifier dependency {path}")
    module = importlib.util.module_from_spec(specification)
    previous = sys.modules.get(module_name)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        if previous is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    return module


_identity_contract = _load_local_contract(
    "leopard2_exact_main_baseline", "leopard2_exact_main_baseline.py")
_record_contract = _load_local_contract(
    "leopard2_exact_main_baseline_record",
    "leopard2_exact_main_baseline_record.py")

from leopard2_exact_main_baseline import (
    ExactMainBaselineError,
    MAX_ELF_INPUT_BYTES,
    canonical_json_bytes,
    canonical_sha256,
    exact_json_equal,
    strict_json_loads,
    verify_normalized_code_identity_against_elf_bytes,
)
from leopard2_exact_main_baseline_record import (
    ACQUISITION_FAILURE_SCHEMA,
    ADAPTER_PATHS,
    AUTHORITY_SCHEMA,
    BASELINE_COMMIT,
    BASELINE_CPP_SOURCES,
    BASELINE_SSE2NEON_COMMIT,
    BASELINE_TREE,
    BUILD_CLOSURE_SCHEMA,
    BUILD_ROLES,
    CANONICAL_LDD_NORMALIZATION,
    SEAL_PROTOCOL,
    VERIFICATION_FAILURE_SCHEMA,
    authority_retained_inventory,
    failure_retained_inventory,
    load_baseline_authority_record,
    load_baseline_failure_record,
    parse_canonical_ldd_output,
    validate_attestation_stdout,
    validate_build_closure,
    validate_cmake_cache,
    validate_compile_commands,
)


VERIFIER_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline-verification/v1"
TREE_METADATA_SCHEMA = "leopard2-exact-main-baseline-tree-metadata/v1"
GIT_CAPTURE_SCHEMA = "leopard2-git-source-capture/v2"
LANE_FILE_MODE = 0o400
LANE_DIRECTORY_MODE = 0o500

MAX_TREE_NODES = 1024
MAX_TREE_DEPTH = 32
MAX_SEALED_FILE_BYTES = 2 * 1024 * 1024 * 1024
MAX_SEALED_TOTAL_BYTES = 8 * 1024 * 1024 * 1024
MAX_PARSED_JSON_BYTES = 64 * 1024 * 1024
MAX_TAR_MEMBERS = 1 << 16
MAX_TAR_PAX_BYTES = 1 << 20
MAX_GIT_COMMIT_BYTES = 16 * 1024 * 1024
MAX_GIT_OBJECT_BYTES = 64 * 1024 * 1024
MAX_GIT_TRACKED_FILE_BYTES = 64 * 1024 * 1024
MAX_GIT_TREE_OBJECTS = 1 << 16
MAX_GIT_TREE_TOTAL_BYTES = 64 * 1024 * 1024
MAX_GIT_TRACKED_FILES = 16 * 1024
READ_CHUNK = 1024 * 1024

_SHA256SUMS_PATH = "SHA256SUMS"
_TREE_METADATA_PATH = "TREE-METADATA.json"
_AUTHORITY_PATH = "baseline-authority.json"
_FAILURE_PATH = "FAILED.json"
_BASELINE_CPP_SOURCES = BASELINE_CPP_SOURCES

_PRODUCER_ATTESTED = tuple(sorted((
    "/created_utc",
    "/attestation/records/*/{exit_status,stderr/content}",
    "/attestation/records/*/stdout/uninterpreted-fields",
    "/attestation/records/*/ctest/"
        "{exit_status,failed,passed,stdout/content,stderr/content}",
    "/builds/*/{configure_argv,build_argv,configure_log/content,"
        "build_log/content}",
    "/builds/*/closure/files/<non-artifact>/{size,sha256}",
    "/host",
    "/lane/{attempt,attempt_budget,root}",
    "/lane/stages/*/log/content",
    "/runtime_closure/records/*/dependencies/*/{sha256,size}",
    "/source/*/archive/{replay_identical,replay_sha256}",
    "/source/*/git_capture/{config,git_executable,git_metadata,"
        "superproject_worktree,worktree_guard_policy}",
    "/source/*/git_capture/{detached,head_ref,tracked_status}",
    "/source/*/git_capture/submodules/*/identity/interior",
    "/toolchain/subtools",
    "/toolchain/tools",
    "/toolchain/versions",
)))

_GIT_CAPTURE_KEYS = frozenset((
    "schema", "path", "head", "tree", "detached", "head_ref",
    "superproject_worktree", "tracked_tree_listing_sha256",
    "tracked_status", "commit_object", "tree_objects", "git_executable",
    "git_metadata", "worktree_guard_policy", "config", "index",
    "tracked_files", "tracked_files_sha256", "submodules",
))


class SealedLaneError(ExactMainBaselineError):
    """The sealed lane is not exact, authentic, or semantically replayable."""


def _fail(message: str) -> NoReturn:
    raise SealedLaneError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _safe_relative_path(path: Any, label: str) -> str:
    _require(type(path) is str and 0 < len(path) <= 4096,
             f"{label} is not a bounded relative path")
    _require(all(0x21 <= ord(character) <= 0x7E for character in path),
             f"{label} is not printable ASCII")
    _require(not path.startswith("/") and not path.endswith("/") and
             "//" not in path and
             all(part not in ("", ".", "..") for part in path.split("/")),
             f"{label} is not canonical")
    return path


def _fingerprint(status: os.stat_result) -> tuple[int, ...]:
    return (
        status.st_dev, status.st_ino, status.st_mode, status.st_nlink,
        status.st_uid, status.st_gid, status.st_size,
        status.st_mtime_ns, status.st_ctime_ns,
    )


def _read_exact_fd(descriptor: int, size: int, label: str) -> bytes:
    _require(0 <= size <= MAX_SEALED_FILE_BYTES,
             f"{label} exceeds its retained byte bound")
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        chunk = os.read(descriptor, min(READ_CHUNK, remaining))
        _require(bool(chunk), f"{label} ended before its recorded size")
        chunks.append(chunk)
        remaining -= len(chunk)
    _require(os.read(descriptor, 1) == b"", f"{label} grew while read")
    return b"".join(chunks)


def _hash_fd(descriptor: int, size: int, label: str) -> str:
    _require(0 <= size <= MAX_SEALED_FILE_BYTES,
             f"{label} exceeds its retained byte bound")
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    remaining = size
    while remaining:
        chunk = os.read(descriptor, min(READ_CHUNK, remaining))
        _require(bool(chunk), f"{label} ended before its recorded size")
        digest.update(chunk)
        remaining -= len(chunk)
    _require(os.read(descriptor, 1) == b"", f"{label} grew while hashed")
    return digest.hexdigest()


def _node_record(path: str, status: os.stat_result, node_type: str) -> dict[str, Any]:
    return {
        "path": path,
        "type": node_type,
        "mode": stat.S_IMODE(status.st_mode),
        "nlink": status.st_nlink,
        "uid": status.st_uid,
        "gid": status.st_gid,
        "size": status.st_size,
        "fingerprint": _fingerprint(status),
    }


def _observe_tree(
    root_descriptor: int,
    *,
    retain_descriptors: bool,
) -> tuple[dict[str, dict[str, Any]], dict[str, str], dict[str, int]]:
    root_status = os.fstat(root_descriptor)
    _require(stat.S_ISDIR(root_status.st_mode), "sealed lane root is not a directory")
    _require(stat.S_IMODE(root_status.st_mode) == LANE_DIRECTORY_MODE,
             "sealed lane root mode is not owner-only 0500")
    _require(root_status.st_uid == os.geteuid() and
             root_status.st_gid == os.getegid(),
             "sealed lane root owner differs from the invoking identity")
    root_device = root_status.st_dev
    nodes: dict[str, dict[str, Any]] = {
        ".": _node_record(".", root_status, "directory")}
    digests: dict[str, str] = {}
    descriptors: dict[str, int] = {}
    total_bytes = 0

    def walk(directory_descriptor: int, prefix: str, depth: int) -> None:
        nonlocal total_bytes
        _require(depth <= MAX_TREE_DEPTH,
                 "sealed lane exceeds its directory-depth bound")
        entries = []
        with os.scandir(directory_descriptor) as iterator:
            for entry in iterator:
                entries.append(entry)
                _require(len(nodes) + len(entries) <= MAX_TREE_NODES,
                         "sealed lane contains too many nodes")
        entries.sort(key=lambda entry: entry.name)
        for entry in entries:
            name = _safe_relative_path(entry.name, "sealed lane node name")
            relative = name if not prefix else prefix + "/" + name
            _safe_relative_path(relative, "sealed lane relative path")
            _require(relative not in nodes,
                     "sealed lane contains a duplicate relative path")
            status = os.stat(name, dir_fd=directory_descriptor,
                             follow_symlinks=False)
            _require(status.st_dev == root_device,
                     f"sealed lane node {relative!r} crosses a device boundary")
            _require(status.st_uid == root_status.st_uid and
                     status.st_gid == root_status.st_gid,
                     f"sealed lane node {relative!r} has another owner")
            if stat.S_ISLNK(status.st_mode):
                _fail(f"sealed lane node {relative!r} is a symlink")
            if stat.S_ISDIR(status.st_mode):
                _require(stat.S_IMODE(status.st_mode) == LANE_DIRECTORY_MODE,
                         f"sealed lane directory {relative!r} is not mode 0500")
                child = os.open(
                    name,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
                    dir_fd=directory_descriptor,
                )
                try:
                    current = os.fstat(child)
                    _require(_fingerprint(current) == _fingerprint(status),
                             f"sealed lane directory {relative!r} was replaced")
                    nodes[relative] = _node_record(relative, current, "directory")
                    _require(len(nodes) <= MAX_TREE_NODES,
                             "sealed lane contains too many nodes")
                    walk(child, relative, depth + 1)
                    _require(_fingerprint(os.fstat(child)) == _fingerprint(current),
                             f"sealed lane directory {relative!r} changed")
                finally:
                    os.close(child)
                continue
            _require(stat.S_ISREG(status.st_mode),
                     f"sealed lane node {relative!r} is not a regular file")
            _require(stat.S_IMODE(status.st_mode) == LANE_FILE_MODE,
                     f"sealed lane file {relative!r} is not mode 0400")
            _require(status.st_nlink == 1,
                     f"sealed lane file {relative!r} is hard-linked")
            _require(0 <= status.st_size <= MAX_SEALED_FILE_BYTES,
                     f"sealed lane file {relative!r} is oversized")
            total_bytes += status.st_size
            _require(total_bytes <= MAX_SEALED_TOTAL_BYTES,
                     "sealed lane total bytes exceed the verifier bound")
            descriptor = os.open(
                name, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                dir_fd=directory_descriptor)
            keep = False
            try:
                current = os.fstat(descriptor)
                _require(_fingerprint(current) == _fingerprint(status),
                         f"sealed lane file {relative!r} was replaced")
                digest = _hash_fd(descriptor, current.st_size, relative)
                _require(_fingerprint(os.fstat(descriptor)) ==
                         _fingerprint(current),
                         f"sealed lane file {relative!r} changed while hashed")
                nodes[relative] = _node_record(relative, current, "file")
                digests[relative] = digest
                _require(len(nodes) <= MAX_TREE_NODES,
                         "sealed lane contains too many nodes")
                if retain_descriptors:
                    descriptors[relative] = descriptor
                    keep = True
            finally:
                if not keep:
                    os.close(descriptor)

    try:
        walk(root_descriptor, "", 1)
    except BaseException:
        for descriptor in descriptors.values():
            os.close(descriptor)
        raise
    return nodes, digests, descriptors


class SealedTree:
    """One fd-anchored, revalidatable snapshot of an owner-only lane."""

    def __init__(self, root: os.PathLike[str] | str):
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
        self.root_path = os.fspath(root)
        self._root_descriptor = os.open(self.root_path, flags)
        self._closed = False
        try:
            self.nodes, self.digests, self._descriptors = _observe_tree(
                self._root_descriptor, retain_descriptors=True)
        except BaseException:
            os.close(self._root_descriptor)
            self._closed = True
            raise

    def __enter__(self) -> "SealedTree":
        return self

    def __exit__(self, _kind: object, _value: object, _traceback: object) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        for descriptor in self._descriptors.values():
            os.close(descriptor)
        os.close(self._root_descriptor)
        self._closed = True

    @property
    def files(self) -> tuple[str, ...]:
        return tuple(sorted(self._descriptors))

    @property
    def directories(self) -> tuple[str, ...]:
        return tuple(sorted(
            path for path, node in self.nodes.items()
            if node["type"] == "directory"))

    def size(self, relative_path: str) -> int:
        _safe_relative_path(relative_path, "retained file path")
        _require(relative_path in self._descriptors,
                 f"retained file {relative_path!r} is absent")
        return int(self.nodes[relative_path]["size"])

    def read_file(self, relative_path: str, *, maximum_bytes: int) -> bytes:
        _safe_relative_path(relative_path, "retained file path")
        _require(relative_path in self._descriptors,
                 f"retained file {relative_path!r} is absent")
        size = self.size(relative_path)
        _require(size <= maximum_bytes,
                 f"retained file {relative_path!r} exceeds its semantic bound")
        descriptor = self._descriptors[relative_path]
        before = os.fstat(descriptor)
        _require(_fingerprint(before) == self.nodes[relative_path]["fingerprint"],
                 f"retained file {relative_path!r} changed before replay")
        content = _read_exact_fd(descriptor, size, relative_path)
        _require(hashlib.sha256(content).hexdigest() ==
                 self.digests[relative_path] and
                 _fingerprint(os.fstat(descriptor)) == _fingerprint(before),
                 f"retained file {relative_path!r} changed during replay")
        return content

    def pread(self, relative_path: str, offset: int, size: int) -> bytes:
        _safe_relative_path(relative_path, "retained file path")
        _require(relative_path in self._descriptors and
                 type(offset) is int and type(size) is int and
                 0 <= offset <= self.size(relative_path) and
                 0 <= size <= self.size(relative_path) - offset,
                 "retained range request is invalid")
        descriptor = self._descriptors[relative_path]
        chunks: list[bytes] = []
        cursor = offset
        remaining = size
        while remaining:
            chunk = os.pread(descriptor, min(READ_CHUNK, remaining), cursor)
            _require(bool(chunk),
                     f"retained file {relative_path!r} ended during range read")
            chunks.append(chunk)
            cursor += len(chunk)
            remaining -= len(chunk)
        return b"".join(chunks)

    def reverify(self) -> None:
        current_nodes, current_digests, current_descriptors = _observe_tree(
            self._root_descriptor, retain_descriptors=False)
        _require(not current_descriptors,
                 "internal verifier descriptor accounting failed")
        _require(exact_json_equal(current_nodes, self.nodes) and
                 exact_json_equal(current_digests, self.digests),
                 "sealed lane changed during verification")


def read_sealed_tree(root: os.PathLike[str] | str) -> SealedTree:
    """Open and snapshot a sealed lane without following links."""
    try:
        return SealedTree(root)
    except SealedLaneError:
        raise
    except OSError as error:
        raise SealedLaneError(f"cannot open sealed lane: {error}") from error


def _metadata_entry(node: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "gid": node["gid"],
        "mode": f'{node["mode"]:04o}',
        "nlink": node["nlink"],
        "path": node["path"],
        "type": node["type"],
        "uid": node["uid"],
    }


def _expected_tree_metadata(tree: SealedTree) -> dict[str, Any]:
    root = tree.nodes["."]
    entries = [
        _metadata_entry(tree.nodes[path])
        for path in sorted(tree.nodes)
        if path != _TREE_METADATA_PATH
    ]
    return {
        "entries": entries,
        "excluded_paths": [_TREE_METADATA_PATH],
        "final_mode_policy": "observed mode with all write bits removed",
        "root": ".",
        "schema": TREE_METADATA_SCHEMA,
        "self_policy": {
            "gid": root["gid"],
            "mode": "0400",
            "nlink": 1,
            "sha256_binding": (
                "exactly one ./TREE-METADATA.json checksum entry"),
            "type": "file",
            "uid": root["uid"],
        },
        "uid_gid_policy": {
            "gid": root["gid"],
            "rule": "every retained node has the invoking effective uid and gid",
            "uid": root["uid"],
        },
    }


def verify_tree_metadata(tree: SealedTree) -> dict[str, Any]:
    content = tree.read_file(
        _TREE_METADATA_PATH, maximum_bytes=MAX_PARSED_JSON_BYTES)
    try:
        observed = strict_json_loads(content, "exact-main tree metadata JSON")
    except ExactMainBaselineError as error:
        raise SealedLaneError("tree metadata is not strict JSON") from error
    expected = _expected_tree_metadata(tree)
    _require(exact_json_equal(observed, expected),
             "tree metadata differs from the observed sealed tree")
    _require(content == canonical_json_bytes(expected),
             "tree metadata bytes are not canonical")
    return copy.deepcopy(expected)


def _expected_sha256sums(tree: SealedTree) -> bytes:
    return b"".join(
        f"{tree.digests[path]}  ./{path}\n".encode("ascii")
        for path in sorted(tree.digests)
        if path != _SHA256SUMS_PATH
    )


def verify_sha256sums(tree: SealedTree) -> dict[str, Any]:
    content = tree.read_file(
        _SHA256SUMS_PATH, maximum_bytes=MAX_PARSED_JSON_BYTES)
    expected = _expected_sha256sums(tree)
    _require(content == expected,
             "SHA256SUMS is not the exact recomputed checksum ledger")
    return {
        "sha256": hashlib.sha256(content).hexdigest(),
        "line_count": len(tree.digests) - 1,
    }


def _derived_directories(paths: Sequence[str]) -> set[str]:
    directories = {"."}
    for path in paths:
        components = path.split("/")[:-1]
        for index in range(1, len(components) + 1):
            directories.add("/".join(components[:index]))
    return directories


def _verify_inventory(
    tree: SealedTree,
    inventory: Sequence[Mapping[str, Any]],
) -> None:
    expected_files = {
        item["relative_path"] for item in inventory
    } | {_SHA256SUMS_PATH, _TREE_METADATA_PATH}
    _require(set(tree.files) == expected_files,
             "sealed lane file set differs from the terminal inventory")
    _require(set(tree.directories) == _derived_directories(expected_files),
             "sealed lane directory set differs from the terminal inventory")
    for item in inventory:
        path = item["relative_path"]
        if item["size"] is None:
            _require(item["sha256"] is None,
                     "terminal inventory entry has a partial byte identity")
            continue
        _require(tree.size(path) == item["size"] and
                 tree.digests[path] == item["sha256"],
                 f"retained file {path!r} differs from its record claim")


def _canonical_terminal(
    tree: SealedTree,
) -> tuple[str, dict[str, Any], tuple[dict[str, Any], ...], bytes]:
    authority_present = _AUTHORITY_PATH in tree.files
    failure_present = _FAILURE_PATH in tree.files
    _require(authority_present != failure_present,
             "sealed lane does not contain exactly one terminal record")
    path = _AUTHORITY_PATH if authority_present else _FAILURE_PATH
    content = tree.read_file(path, maximum_bytes=MAX_PARSED_JSON_BYTES)
    try:
        if authority_present:
            record = load_baseline_authority_record(content)
            inventory = authority_retained_inventory(record)
        else:
            record = load_baseline_failure_record(content)
            inventory = failure_retained_inventory(record)
    except ExactMainBaselineError as error:
        raise SealedLaneError("terminal record is invalid") from error
    _require(content == canonical_json_bytes(record),
             "terminal record bytes are not canonical")
    _require(record["lane"]["seal_protocol"] == SEAL_PROTOCOL,
             "terminal record names another seal protocol")
    _verify_inventory(tree, inventory)
    return path, record, inventory, content


def _strict_json_file(tree: SealedTree, path: str, label: str) -> Any:
    content = tree.read_file(path, maximum_bytes=MAX_PARSED_JSON_BYTES)
    try:
        return strict_json_loads(content, label)
    except ExactMainBaselineError as error:
        raise SealedLaneError(f"{label} is not strict JSON") from error


def _exact_object(value: Any, keys: set[str] | frozenset[str], label: str
                  ) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == set(keys),
             f"{label} has an unexpected object shape")
    return value


def _git_blob_sha1(content: bytes) -> str:
    return _git_object_sha1("blob", content)


def _git_object_sha1(kind: str, content: bytes) -> str:
    digest = hashlib.sha1(usedforsecurity=False)
    digest.update(f"{kind} {len(content)}\0".encode("ascii"))
    digest.update(content)
    return digest.hexdigest()


def _git_capture_digest(value: Any) -> str:
    content = canonical_json_bytes(value)
    _require(content.endswith(b"\n"),
             "internal canonical Git-capture encoding lost its delimiter")
    return hashlib.sha256(content[:-1]).hexdigest()


def _is_lower_hex(value: Any, length: int) -> bool:
    return type(value) is str and len(value) == length and all(
        character in "0123456789abcdef" for character in value)


def _byte_identity(content: bytes) -> dict[str, Any]:
    return {"size": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def _validate_byte_identity(value: Any, label: str) -> dict[str, Any]:
    identity = _exact_object(value, {"size", "sha256"}, label)
    _require(type(identity["size"]) is int and 0 <= identity["size"] and
             _is_lower_hex(identity["sha256"], 64),
             f"{label} byte identity is invalid")
    return identity


def _validate_git_object_identity(
    value: Any,
    *,
    kind: str,
    expected_id: str,
    label: str,
    maximum_bytes: int,
) -> bytes:
    record = _exact_object(
        value, {"encoding", "size", "sha256", "object_id", "base64"}, label)
    _require(record["encoding"] == "base64" and
             record["object_id"] == expected_id and
             _is_lower_hex(expected_id, 40) and
             type(record["size"]) is int and
             0 <= record["size"] <= maximum_bytes and
             _is_lower_hex(record["sha256"], 64) and
             type(record["base64"]) is str and
             len(record["base64"]) <= ((maximum_bytes + 2) // 3) * 4,
             f"{label} object identity is invalid")
    try:
        content = base64.b64decode(record["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise SealedLaneError(f"{label} object is not canonical base64") \
            from error
    _require(base64.b64encode(content).decode("ascii") == record["base64"] and
             len(content) == record["size"] and
             hashlib.sha256(content).hexdigest() == record["sha256"] and
             _git_object_sha1(kind, content) == expected_id,
             f"{label} object bytes differ")
    return content


def _parse_git_tree_object(content: bytes, label: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    offset = 0
    while offset < len(content):
        space = content.find(b" ", offset)
        nul = content.find(b"\0", space + 1 if space >= 0 else offset)
        _require(space > offset and nul > space + 1 and nul + 21 <= len(content),
                 f"{label} contains a truncated tree entry")
        mode_raw = content[offset:space]
        name_raw = content[space + 1:nul]
        object_raw = content[nul + 1:nul + 21]
        offset = nul + 21
        try:
            mode = mode_raw.decode("ascii", errors="strict")
            name = name_raw.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise SealedLaneError(f"{label} contains non-canonical text") \
                from error
        kind = (
            "tree" if mode == "40000" else
            "blob" if mode in ("100644", "100755", "120000") else
            "commit" if mode == "160000" else None)
        _require(kind is not None and "/" not in name and
                 _safe_relative_path(name, f"{label} entry name") == name,
                 f"{label} contains an unsupported entry")
        entries.append({
            "name": name,
            "git_mode": mode,
            "git_type": kind,
            "object_id": object_raw.hex(),
        })
        _require(len(entries) <= MAX_GIT_TRACKED_FILES + MAX_GIT_TREE_OBJECTS,
                 f"{label} contains too many entries")
    sort_keys = [
        entry["name"].encode("utf-8") +
        (b"/" if entry["git_type"] == "tree" else b"")
        for entry in entries
    ]
    _require(sort_keys == sorted(sort_keys) and
             len(sort_keys) == len(set(sort_keys)),
             f"{label} entries are duplicated or not in Git tree order")
    return entries


def _flatten_git_tree_objects(
    root_tree: str,
    objects: Mapping[str, bytes],
) -> list[dict[str, str]]:
    _require(root_tree in objects,
             "retained Git tree-object closure omits its root")
    result: list[dict[str, str]] = []
    reachable: set[str] = set()
    stack: list[tuple[str, Any, str, frozenset[str]]] = [
        ("tree", root_tree, "", frozenset())]
    expansions = 0
    while stack:
        kind, payload, prefix, active = stack.pop()
        expansions += 1
        _require(expansions <= MAX_GIT_TRACKED_FILES + MAX_GIT_TREE_OBJECTS,
                 "retained Git recursive tree expansion exceeds its bound")
        if kind == "leaf":
            entry = payload
            _require(type(entry) is dict,
                     "retained Git recursive tree leaf is invalid")
            result.append({
                "path": prefix,
                "git_mode": entry["git_mode"],
                "git_type": entry["git_type"],
                "object_id": entry["object_id"],
            })
            _require(len(result) <= MAX_GIT_TRACKED_FILES,
                     "retained Git recursive tree has too many leaves")
            continue
        object_id = payload
        _require(type(object_id) is str and object_id not in active and
                 object_id in objects,
                 "retained Git tree-object closure is cyclic or incomplete")
        reachable.add(object_id)
        next_active = active | {object_id}
        entries = _parse_git_tree_object(
            objects[object_id], f"retained Git tree object {object_id}")
        for entry in reversed(entries):
            relative = f"{prefix}/{entry['name']}" if prefix else entry["name"]
            _safe_relative_path(relative, "retained Git recursive tree path")
            if entry["git_type"] == "tree":
                stack.append(("tree", entry["object_id"], relative, next_active))
            else:
                stack.append(("leaf", entry, relative, next_active))
    _require(reachable == set(objects),
             "retained Git tree-object closure contains unreachable objects")
    paths = [entry["path"] for entry in result]
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "retained Git recursive tree inventory is non-canonical")
    return result


def _tar_octal(field: bytes, label: str) -> int:
    stripped = field.rstrip(b"\0 ").lstrip(b" ")
    _require(bool(stripped) and all(0x30 <= byte <= 0x37 for byte in stripped),
             f"{label} is not canonical octal")
    return int(stripped, 8)


def _tar_text(field: bytes, label: str) -> str:
    raw = field.split(b"\0", 1)[0]
    _require(b"\0" not in raw and bool(raw), f"{label} is empty")
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise SealedLaneError(f"{label} is not strict UTF-8") from error
    _require(all(ord(character) >= 0x20 and
                 not (0x7F <= ord(character) <= 0x9F) and
                 not (0xD800 <= ord(character) <= 0xDFFF)
                 for character in text),
             f"{label} contains unsafe text")
    return text


def _safe_tar_name(name: str, expected_prefix: str, *, directory: bool) -> str:
    _require(not name.startswith("/") and "//" not in name,
             "source archive contains a non-canonical path")
    stripped = name[:-1] if directory and name.endswith("/") else name
    _require(bool(stripped) and
             all(part not in ("", ".", "..") for part in stripped.split("/")),
             "source archive contains an unsafe path component")
    prefix_root = expected_prefix[:-1]
    _require(stripped == prefix_root or stripped.startswith(expected_prefix),
             "source archive member escaped its frozen prefix")
    return stripped


def _parse_pax(content: bytes, label: str) -> dict[str, str]:
    _require(len(content) <= MAX_TAR_PAX_BYTES,
             f"{label} exceeds its byte bound")
    result: dict[str, str] = {}
    offset = 0
    while offset < len(content):
        space = content.find(b" ", offset)
        _require(space > offset, f"{label} has a malformed record length")
        length_raw = content[offset:space]
        _require(length_raw.isdigit() and not length_raw.startswith(b"0"),
                 f"{label} has a non-canonical record length")
        length = int(length_raw)
        end = offset + length
        _require(end <= len(content) and content[end - 1:end] == b"\n",
                 f"{label} contains a truncated record")
        payload = content[space + 1:end - 1]
        key_raw, separator, value_raw = payload.partition(b"=")
        _require(bool(separator), f"{label} record has no value")
        try:
            key = key_raw.decode("ascii")
            value = value_raw.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise SealedLaneError(f"{label} record is not canonical text") \
                from error
        _require(key in ("comment", "path", "linkpath") and key not in result,
                 f"{label} contains an unsupported or duplicate key")
        _require(value and "\0" not in value and "\n" not in value and
                 "\r" not in value,
                 f"{label} contains an unsafe value")
        result[key] = value
        offset = end
    _require(offset == len(content), f"{label} has trailing bytes")
    return result


def _tar_member_bytes(
    tree: SealedTree,
    archive_path: str,
    offset: int,
    size: int,
) -> bytes:
    _require(size <= MAX_SEALED_FILE_BYTES,
             "source archive member exceeds its byte bound")
    return tree.pread(archive_path, offset, size)


def _tar_member_identity(
    tree: SealedTree,
    archive_path: str,
    offset: int,
    size: int,
) -> dict[str, Any]:
    _require(0 <= size <= MAX_SEALED_FILE_BYTES,
             "source archive member exceeds its byte bound")
    sha256 = hashlib.sha256()
    git_sha1 = hashlib.sha1(usedforsecurity=False)
    git_sha1.update(f"blob {size}\0".encode("ascii"))
    cursor = offset
    remaining = size
    while remaining:
        chunk = tree.pread(
            archive_path, cursor, min(READ_CHUNK, remaining))
        sha256.update(chunk)
        git_sha1.update(chunk)
        cursor += len(chunk)
        remaining -= len(chunk)
    return {
        "size": size,
        "sha256": sha256.hexdigest(),
        "git_blob_sha1": git_sha1.hexdigest(),
    }


def _verify_source_archive(
    tree: SealedTree,
    archive_path: str,
    *,
    expected_prefix: str,
    expected_commit: str,
) -> dict[str, dict[str, Any]]:
    total_size = tree.size(archive_path)
    _require(total_size >= 1024 and total_size % 512 == 0,
             "source archive is not a block-aligned tar stream")
    observed: dict[str, dict[str, Any]] = {}
    global_pax: dict[str, str] = {}
    pending_pax: dict[str, str] = {}
    member_count = 0
    offset = 0
    zero_blocks = 0
    while offset < total_size:
        block = tree.pread(archive_path, offset, 512)
        if block == b"\0" * 512:
            zero_blocks += 1
            offset += 512
            if zero_blocks < 2:
                continue
            while offset < total_size:
                chunk_size = min(READ_CHUNK, total_size - offset)
                _require(tree.pread(archive_path, offset, chunk_size) ==
                         b"\0" * chunk_size,
                         "source archive has bytes after its terminator")
                offset += chunk_size
            break
        _require(zero_blocks == 0,
                 "source archive has a lone zero header block")
        member_count += 1
        _require(member_count <= MAX_TAR_MEMBERS,
                 "source archive contains too many members")
        recorded_checksum = _tar_octal(block[148:156], "tar header checksum")
        checksum_block = block[:148] + b" " * 8 + block[156:]
        _require(sum(checksum_block) == recorded_checksum,
                 "source archive header checksum is invalid")
        _require(block[257:263] == b"ustar\0" and
                 block[263:265] in (b"00", b"\0\0"),
                 "source archive is not canonical ustar")
        size = _tar_octal(block[124:136], "tar member size")
        mode = _tar_octal(block[100:108], "tar member mode")
        type_flag = block[156:157]
        name = _tar_text(block[0:100], "tar member name")
        prefix_raw = block[345:500].split(b"\0", 1)[0]
        if prefix_raw:
            prefix = _tar_text(block[345:500], "tar header prefix")
            name = prefix + "/" + name
        data_offset = offset + 512
        padded_size = (size + 511) // 512 * 512
        _require(data_offset + padded_size <= total_size,
                 "source archive member extends beyond the archive")
        if type_flag in (b"g", b"x"):
            _require(size <= MAX_TAR_PAX_BYTES,
                     "source archive pax header exceeds its byte bound")
            pax = _parse_pax(
                _tar_member_bytes(tree, archive_path, data_offset, size),
                "tar pax header")
            if type_flag == b"g":
                _require(set(pax) == {"comment"} and not global_pax,
                         "source archive global pax header changed")
                global_pax = pax
            else:
                _require(not pending_pax,
                         "source archive stacked extended pax headers")
                pending_pax = pax
            offset = data_offset + padded_size
            continue
        if "path" in pending_pax:
            name = pending_pax["path"]
        _require("linkpath" not in pending_pax,
                 "source archive retained a link target")
        pending_pax = {}
        if type_flag in (b"0", b"\0"):
            canonical_name = _safe_tar_name(
                name, expected_prefix, directory=False)
            _require(mode in (0o664, 0o775),
                     "source archive file mode changed")
            relative = canonical_name[len(expected_prefix):]
            _require(relative not in observed,
                     "source archive duplicated a file member")
            observed[relative] = {
                **_tar_member_identity(
                    tree, archive_path, data_offset, size),
                "git_mode": "100755" if mode == 0o775 else "100644",
            }
        elif type_flag == b"5":
            _safe_tar_name(name, expected_prefix, directory=True)
            _require(mode == 0o775 and size == 0,
                     "source archive directory metadata changed")
        else:
            _fail("source archive contains a link or special member")
        offset = data_offset + padded_size
    _require(zero_blocks >= 2 and not pending_pax,
             "source archive has no complete terminator")
    _require(global_pax.get("comment") == expected_commit,
             "source archive global commit identity changed")
    _require(bool(observed), "source archive contains no regular files")
    return observed


def _verify_archive_git_binding(
    archive_files: Mapping[str, Mapping[str, Any]],
    tracked_files: Mapping[str, Mapping[str, Any]],
    label: str,
) -> None:
    _require(not any(item["kind"] == "symlink"
                     for item in tracked_files.values()),
             f"{label} Git capture contains an unsupported symlink")
    expected = {
        path: item for path, item in tracked_files.items()
        if item["kind"] == "regular"
    }
    _require(set(archive_files) == set(expected),
             f"{label} archive file set differs from its Git tree")
    for path, identity in archive_files.items():
        tracked = expected[path]
        _require(identity["git_mode"] == tracked["git_mode"] and
                 identity["git_blob_sha1"] == tracked["object_id"],
                 f"{label} archive member {path!r} differs from its Git blob")


def _verify_git_capture(
    value: Any,
    *,
    expected_path: str,
    expected_head: str,
    expected_tree: str,
    expected_submodule: tuple[str, str] | None,
) -> dict[str, dict[str, Any]]:
    capture = _exact_object(value, _GIT_CAPTURE_KEYS, "Git source capture")
    _require(capture["schema"] == GIT_CAPTURE_SCHEMA and
             capture["path"] == expected_path and
             capture["head"] == expected_head and
             capture["tree"] == expected_tree and
             _is_lower_hex(expected_head, 40) and
             _is_lower_hex(expected_tree, 40) and
             type(capture["detached"]) is bool and
             ((capture["detached"] and capture["head_ref"] is None) or
              (not capture["detached"] and
               type(capture["head_ref"]) is str and
               capture["head_ref"].startswith("refs/"))) and
             capture["tracked_status"] == "clean" and
             _is_lower_hex(capture["tracked_tree_listing_sha256"], 64) and
             type(capture["tree_objects"]) is list and
             type(capture["tracked_files"]) is list and
             type(capture["submodules"]) is list,
             "Git source capture identity changed")

    commit = _validate_git_object_identity(
        capture["commit_object"], kind="commit", expected_id=expected_head,
        label="retained Git commit", maximum_bytes=MAX_GIT_COMMIT_BYTES)
    _require(b"\n\n" in commit,
             "retained Git commit has no header/message boundary")
    headers = commit.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + expected_tree.encode("ascii")
    _require(headers and headers[0] == tree_line and
             [line for line in headers if line.startswith(b"tree ")] ==
             [tree_line],
             "retained Git commit names a different tree")

    raw_tree_objects = capture["tree_objects"]
    _require(0 < len(raw_tree_objects) <= MAX_GIT_TREE_OBJECTS,
             "retained Git tree-object closure is invalid")
    tree_objects: dict[str, bytes] = {}
    tree_total = 0
    previous_object_id: str | None = None
    for index, raw in enumerate(raw_tree_objects):
        object_id = raw.get("object_id") if type(raw) is dict else None
        _require(_is_lower_hex(object_id, 40) and
                 (previous_object_id is None or
                  previous_object_id < object_id) and
                 object_id not in tree_objects,
                 "retained Git tree-object closure is not canonical")
        content = _validate_git_object_identity(
            raw, kind="tree", expected_id=object_id,
            label=f"retained Git tree object {index}",
            maximum_bytes=MAX_GIT_OBJECT_BYTES)
        tree_total += len(content)
        _require(tree_total <= MAX_GIT_TREE_TOTAL_BYTES,
                 "retained Git tree-object closure exceeds its byte bound")
        tree_objects[object_id] = content
        previous_object_id = object_id
    tree_entries = _flatten_git_tree_objects(expected_tree, tree_objects)

    index = _exact_object(
        capture["index"], {"entry_count", "stage", "flags_v", "flags_f"},
        "retained Git index identity")
    records = capture["tracked_files"]
    _require(type(index["entry_count"]) is int and
             0 <= index["entry_count"] <= MAX_GIT_TRACKED_FILES and
             len(records) == index["entry_count"] and
             len(records) <= MAX_GIT_TRACKED_FILES and
             capture["tracked_files_sha256"] == _git_capture_digest(records),
             "retained tracked-file inventory is invalid")
    paths: list[str] = []
    normalized_records: list[dict[str, Any]] = []
    submodule_records: dict[str, dict[str, Any]] = {}
    for item_index, raw in enumerate(records):
        _require(type(raw) is dict and raw.get("kind") in (
                     "regular", "symlink", "submodule"),
                 f"retained tracked-file record {item_index} is invalid")
        kind = raw["kind"]
        common = {"path", "git_mode", "git_type", "object_id", "kind"}
        if kind == "regular":
            expected_keys = common
            valid_mode = raw.get("git_mode") in ("100644", "100755")
        elif kind == "symlink":
            expected_keys = common | {
                "size", "sha256", "target_encoding", "target_base64"}
            valid_mode = raw.get("git_mode") == "120000"
        else:
            expected_keys = common | {"identity_sha256"}
            valid_mode = raw.get("git_mode") == "160000"
        _require(set(raw) == expected_keys and valid_mode and
                 raw.get("git_type") ==
                    ("commit" if kind == "submodule" else "blob") and
                 _is_lower_hex(raw.get("object_id"), 40),
                 f"retained tracked-file record {item_index} shape differs")
        path = _safe_relative_path(
            raw.get("path"), f"retained tracked-file record {item_index} path")
        paths.append(path)
        if kind == "symlink":
            _require(raw.get("target_encoding") == "base64" and
                     type(raw.get("target_base64")) is str and
                     type(raw.get("size")) is int and
                     0 <= raw["size"] <= MAX_GIT_TRACKED_FILE_BYTES and
                     _is_lower_hex(raw.get("sha256"), 64) and
                     len(raw["target_base64"]) <=
                        ((MAX_GIT_TRACKED_FILE_BYTES + 2) // 3) * 4,
                     "retained tracked symlink identity is invalid")
            try:
                target = base64.b64decode(raw["target_base64"], validate=True)
            except (ValueError, binascii.Error) as error:
                raise SealedLaneError(
                    "retained tracked symlink target is not canonical base64") \
                    from error
            _require(base64.b64encode(target).decode("ascii") ==
                     raw["target_base64"] and len(target) == raw["size"] and
                     hashlib.sha256(target).hexdigest() == raw["sha256"] and
                     _git_blob_sha1(target) == raw["object_id"],
                     "retained tracked symlink bytes differ")
        if kind == "submodule":
            _require(_is_lower_hex(raw.get("identity_sha256"), 64),
                     "retained submodule digest is invalid")
            submodule_records[path] = raw
        normalized_records.append(raw)
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "retained tracked-file inventory is non-canonical")
    _require([{
        "path": item["path"],
        "git_mode": item["git_mode"],
        "git_type": item["git_type"],
        "object_id": item["object_id"],
    } for item in normalized_records] == tree_entries,
             "retained tracked-file inventory differs from recursive trees")

    tree_listing = b"".join(
        (f"{item['git_mode']} {item['git_type']} {item['object_id']}\t"
         f"{item['path']}\0").encode("utf-8")
        for item in normalized_records)
    index_stage = b"".join(
        (f"{item['git_mode']} {item['object_id']} 0\t{item['path']}\0").
        encode("utf-8") for item in normalized_records)
    default_flags = b"".join(
        f"H {item['path']}\0".encode("utf-8")
        for item in normalized_records)
    _require(capture["tracked_tree_listing_sha256"] ==
             hashlib.sha256(tree_listing).hexdigest() and
             _validate_byte_identity(index["stage"], "Git index stage") ==
             _byte_identity(index_stage) and
             _validate_byte_identity(index["flags_v"], "Git index flags-v") ==
             _byte_identity(default_flags) and
             _validate_byte_identity(index["flags_f"], "Git index flags-f") ==
             _byte_identity(default_flags),
             "retained Git tree/index transcripts differ")
    _validate_byte_identity(capture["config"], "retained Git config")

    submodules = capture["submodules"]
    _require(len(submodules) == len(submodule_records),
             "retained submodule inventory differs from Git links")
    observed_submodules: list[str] = []
    for raw in submodules:
        item = _exact_object(
            raw, {"path", "object_id", "identity_sha256", "identity"},
            "retained submodule record")
        path = _safe_relative_path(item["path"], "retained submodule path")
        nested = item["identity"]
        _require(path in submodule_records and
                 item["object_id"] == submodule_records[path]["object_id"] and
                 item["identity_sha256"] ==
                    submodule_records[path]["identity_sha256"] and
                 item["identity_sha256"] == _git_capture_digest(nested) and
                 type(nested) is dict and
                 nested.get("schema") == GIT_CAPTURE_SCHEMA and
                 nested.get("head") == item["object_id"] and
                 nested.get("tracked_status") == "clean",
                 "retained submodule record is inconsistent")
        observed_submodules.append(path)
    _require(observed_submodules == sorted(observed_submodules) and
             len(observed_submodules) == len(set(observed_submodules)),
             "retained submodule inventory is non-canonical")
    if expected_submodule is not None:
        submodule_path, submodule_commit = expected_submodule
        _require(len(submodules) == 1 and
                 submodules[0]["path"] == submodule_path and
                 submodules[0]["object_id"] == submodule_commit,
                 "Leopard1 submodule capture identity changed")
    _require(type(capture["git_executable"]) is dict and
             type(capture["git_metadata"]) is dict and
             type(capture["worktree_guard_policy"]) is str,
             "retained producer-attested Git metadata shape differs")
    return {item["path"]: item for item in normalized_records}


def _verify_build_closure(
    value: Any,
    role: str,
    build: Mapping[str, Any],
) -> None:
    try:
        validate_build_closure(value, role=role, build=build)
    except ExactMainBaselineError as error:
        raise SealedLaneError(str(error)) from error


def _verify_compile_commands(
    value: Any,
    build: Mapping[str, Any],
    compiler_path: str,
) -> None:
    try:
        validate_compile_commands(
            value, roots=build["roots"], compiler=compiler_path,
            profile=_record_contract.exact_main_build_profile())
    except ExactMainBaselineError as error:
        raise SealedLaneError(str(error)) from error


def _verify_attestation_stdout(value: Any, record: Mapping[str, Any]) -> None:
    try:
        validate_attestation_stdout(
            value, argv=record["argv"],
            reported_schema=record["reported_schema"])
    except ExactMainBaselineError as error:
        raise SealedLaneError(str(error)) from error


def _verify_authority_semantics(
    tree: SealedTree,
    record: Mapping[str, Any],
) -> dict[str, Any]:
    rederived: dict[str, dict[str, Any]] = {}
    executable_digests: dict[str, str] = {}
    archive_digests: dict[str, str] = {}
    for role in BUILD_ROLES:
        build = record["builds"][role]
        executable_path = build["executable"]["retained_relative_path"]
        executable = tree.read_file(
            executable_path, maximum_bytes=MAX_ELF_INPUT_BYTES)
        executable_digests[role] = tree.digests[executable_path]
        try:
            rederived[role] = \
                verify_normalized_code_identity_against_elf_bytes(
                    executable, record["identity"][role],
                    roots=build["roots"])
        except ExactMainBaselineError as error:
            raise SealedLaneError(
                f"{role} normalized ELF identity replay failed") from error
        archive_path = build["archive"]["retained_relative_path"]
        archive_digests[role] = tree.digests[archive_path]

    canonical_executable_equal = (
        executable_digests["canonical_first"] ==
        executable_digests["canonical_second"])
    canonical_archive_equal = (
        archive_digests["canonical_first"] ==
        archive_digests["canonical_second"])
    variant_raw_differs = (
        executable_digests["path_variant"] !=
        executable_digests["canonical_first"])
    combined = [rederived[role]["combined_sha256"] for role in BUILD_ROLES]
    normalized_match = len(set(combined)) == 1
    census_zero = all(
        root["occurrences"] == 0 and
        all(row["occurrences"] == 0 for row in root["sections"])
        for role in BUILD_ROLES
        for root in rederived[role]["path_string_census"]["roots"]
    )
    _require(canonical_executable_equal,
             "same-path canonical executable bytes differ")
    _require(canonical_archive_equal,
             "same-path canonical archive bytes differ")
    _require(variant_raw_differs,
             "path-variant executable did not differ in raw bytes")
    _require(normalized_match,
             "path-variant normalized ELF identity changed")
    _require(census_zero,
             "selected allocatable ELF sections retain an acquisition root")

    canonical_build = record["builds"]["canonical_first"]
    baseline_capture = _strict_json_file(
        tree, record["source"]["baseline"]["git_capture"]["relative_path"],
        "Leopard1 Git capture JSON")
    baseline_tracked_files = _verify_git_capture(
        baseline_capture,
        expected_path=canonical_build["roots"]["baseline_source_root"],
        expected_head=BASELINE_COMMIT,
        expected_tree=BASELINE_TREE,
        expected_submodule=("sse2neon", BASELINE_SSE2NEON_COMMIT),
    )
    adapter_capture = _strict_json_file(
        tree,
        record["source"]["adapter_repository"]["git_capture"]["relative_path"],
        "adapter Git capture JSON")
    adapter_tracked_files = _verify_git_capture(
        adapter_capture,
        expected_path=canonical_build["roots"]["adapter_source_root"],
        expected_head=record["source"]["adapter_repository"]["commit"],
        expected_tree=record["source"]["adapter_repository"]["tree"],
        expected_submodule=None,
    )

    adapter_archive = record["source"]["adapter_repository"]["archive"]
    adapter_members = _verify_source_archive(
        tree, adapter_archive["relative_path"],
        expected_prefix=adapter_archive["prefix"],
        expected_commit=record["source"]["adapter_repository"]["commit"],
    )
    _verify_archive_git_binding(
        adapter_members, adapter_tracked_files, "adapter source")
    for expected, path in zip(record["adapter"]["files"], ADAPTER_PATHS):
        _require(path in adapter_members and
                 adapter_members[path]["size"] == expected["size"] and
                 adapter_members[path]["sha256"] == expected["sha256"] and
                 adapter_members[path]["git_blob_sha1"] ==
                    expected["git_blob_sha1"],
                 f"adapter source archive member {path!r} changed")

    baseline_archive = record["source"]["baseline"]["archive"]
    baseline_members = _verify_source_archive(
        tree, baseline_archive["relative_path"],
        expected_prefix=baseline_archive["prefix"],
        expected_commit=BASELINE_COMMIT,
    )
    _verify_archive_git_binding(
        baseline_members, baseline_tracked_files, "Leopard1 source")
    _require(all(path in baseline_members for path in _BASELINE_CPP_SOURCES),
             "Leopard1 source archive omits a compiled translation unit")

    controller_path = record["attestation"]["test_controller"]["relative_path"]
    controller = tree.read_file(
        controller_path, maximum_bytes=MAX_PARSED_JSON_BYTES)
    adapter_controller = record["adapter"]["files"][2]
    _require(len(controller) == adapter_controller["size"] and
             hashlib.sha256(controller).hexdigest() ==
             adapter_controller["sha256"] and
             _git_blob_sha1(controller) == adapter_controller["git_blob_sha1"],
             "retained attestation controller differs from adapter source")

    compiler_path = next(
        item["resolved_path"] for item in record["toolchain"]["tools"]
        if item["role"] == "compiler")
    for role in BUILD_ROLES:
        build = record["builds"][role]
        cache_content = tree.read_file(
            build["cmake_cache"]["relative_path"],
            maximum_bytes=MAX_PARSED_JSON_BYTES)
        try:
            validate_cmake_cache(cache_content, build["roots"])
        except ExactMainBaselineError as error:
            raise SealedLaneError(str(error)) from error
        compile_commands = _strict_json_file(
            tree, build["compile_commands"]["relative_path"],
            f"{role} compile commands JSON")
        _verify_compile_commands(
            compile_commands, build, compiler_path)
        closure = _strict_json_file(
            tree, build["closure"]["relative_path"],
            f"{role} build closure JSON")
        _verify_build_closure(closure, role, build)

    expected_ldd_bytes: bytes | None = None
    for role, runtime_record, attestation_record in zip(
            BUILD_ROLES, record["runtime_closure"]["records"],
            record["attestation"]["records"]):
        _require(record["runtime_closure"]["normalization"] ==
                 CANONICAL_LDD_NORMALIZATION,
                 "runtime closure normalization changed")
        ldd_path = runtime_record["canonical_ldd_output"]["relative_path"]
        ldd_bytes = tree.read_file(
            ldd_path, maximum_bytes=MAX_PARSED_JSON_BYTES)
        try:
            parsed_ldd = parse_canonical_ldd_output(ldd_bytes)
        except ExactMainBaselineError as error:
            raise SealedLaneError(
                f"{role} canonical runtime closure is invalid") from error
        recorded_projection = tuple({
            "soname": dependency["soname"],
            "kind": dependency["kind"],
            "path": dependency["path"],
        } for dependency in runtime_record["dependencies"])
        _require(exact_json_equal(parsed_ldd, recorded_projection),
                 f"{role} canonical runtime dependency projection changed")
        if expected_ldd_bytes is None:
            expected_ldd_bytes = ldd_bytes
        else:
            _require(ldd_bytes == expected_ldd_bytes,
                     "runtime dependency closure changed between builds")

        stdout = _strict_json_file(
            tree, attestation_record["stdout"]["relative_path"],
            f"{role} benchmark attestation JSON")
        _verify_attestation_stdout(stdout, attestation_record)

    promotion = record["promotion"]
    recomputed_promotion = {
        "same_path_executable_bytes_identical": canonical_executable_equal,
        "same_path_archive_bytes_identical": canonical_archive_equal,
        "path_variant_raw_executable_differs": variant_raw_differs,
        "path_variant_normalized_match": normalized_match,
        "selected_section_census_zero": census_zero,
        "pure_avx2_attested": all(
            item["pure_avx2"] is True and item["round_trip"] is True
            for item in record["attestation"]["records"]),
        "historical_references_non_authoritative": all(
            item["authority"] is False
            for item in record["superseded_references"]),
    }
    _require(all(recomputed_promotion.values()) and
             promotion["promoted"] is True and
             all(promotion[key] is value
                 for key, value in recomputed_promotion.items()),
             "authority promotion gate differs from independent replay")
    return {
        "adapter_archive_members_verified": len(ADAPTER_PATHS),
        "attestation_records_verified": len(BUILD_ROLES),
        "build_closures_verified": len(BUILD_ROLES),
        "canonical_archive_bytes_identical": canonical_archive_equal,
        "canonical_executable_bytes_identical": canonical_executable_equal,
        "git_capture_object_closures_verified": 2,
        "normalized_identity_records_verified": len(BUILD_ROLES),
        "path_variant_normalized_match": normalized_match,
        "path_variant_raw_executable_differs": variant_raw_differs,
        "runtime_closures_verified": len(BUILD_ROLES),
        "selected_section_census_zero": census_zero,
        "source_archive_git_bindings_verified": 2,
    }


def _verdict_sha256(value: Mapping[str, Any]) -> str:
    return canonical_sha256({
        key: copy.deepcopy(field) for key, field in value.items()
        if key != "verdict_sha256"
    })


def _seal_verdict(
    tree: SealedTree,
    metadata: Mapping[str, Any],
    ledger: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "protocol": SEAL_PROTOCOL,
        "tree_metadata_schema": TREE_METADATA_SCHEMA,
        "tree_metadata_sha256": tree.digests[_TREE_METADATA_PATH],
        "sha256sums_sha256": ledger["sha256"],
        "directory_count": len(tree.directories),
        "file_count": len(tree.files),
        "checksum_line_count": ledger["line_count"],
        "metadata_entry_count": len(metadata["entries"]),
    }


def _build_verdict(
    *,
    outcome: str,
    promoted: bool,
    path: str,
    record: Mapping[str, Any],
    content: bytes,
    tree: SealedTree,
    metadata: Mapping[str, Any],
    ledger: Mapping[str, Any],
    recomputed: Mapping[str, Any],
) -> dict[str, Any]:
    value = {
        "schema": VERIFIER_SCHEMA,
        "outcome": outcome,
        "promoted": promoted,
        "record": {
            "relative_path": path,
            "schema": record["schema"],
            "record_sha256": record["record_sha256"],
            "canonical_bytes_identical": content == canonical_json_bytes(record),
        },
        "seal": _seal_verdict(tree, metadata, ledger),
        "recomputed": copy.deepcopy(dict(recomputed)),
        "producer_attested": list(_PRODUCER_ATTESTED),
        "verdict_sha256": "0" * 64,
    }
    value["verdict_sha256"] = _verdict_sha256(value)
    return value


def verify_authority_lane(tree: SealedTree) -> dict[str, Any]:
    metadata = verify_tree_metadata(tree)
    ledger = verify_sha256sums(tree)
    path, record, _inventory, content = _canonical_terminal(tree)
    _require(path == _AUTHORITY_PATH and record["schema"] == AUTHORITY_SCHEMA,
             "sealed authority lane has another terminal schema")
    recomputed = _verify_authority_semantics(tree, record)
    tree.reverify()
    return _build_verdict(
        outcome="promoted_authority", promoted=True, path=path,
        record=record, content=content, tree=tree, metadata=metadata,
        ledger=ledger, recomputed=recomputed)


def verify_failure_lane(tree: SealedTree) -> dict[str, Any]:
    metadata = verify_tree_metadata(tree)
    ledger = verify_sha256sums(tree)
    path, record, inventory, content = _canonical_terminal(tree)
    _require(path == _FAILURE_PATH and record["schema"] in (
                 ACQUISITION_FAILURE_SCHEMA, VERIFICATION_FAILURE_SCHEMA) and
             record["promoted"] is False,
             "sealed failure lane was relabeled as authority")
    recomputed = {
        "failure_record_verified": True,
        "retained_file_claims_verified": len(inventory) - 1,
        "stage_log_count": len(record["lane"]["stages"]),
        "terminal_nonpromotion_verified": True,
    }
    tree.reverify()
    return _build_verdict(
        outcome="verified_failure", promoted=False, path=path,
        record=record, content=content, tree=tree, metadata=metadata,
        ledger=ledger, recomputed=recomputed)


def verify_sealed_lane(root: os.PathLike[str] | str) -> dict[str, Any]:
    """Verify one complete sealed lane and return a canonical verdict."""
    with read_sealed_tree(root) as tree:
        authority_present = _AUTHORITY_PATH in tree.files
        failure_present = _FAILURE_PATH in tree.files
        _require(authority_present != failure_present,
                 "sealed lane does not contain exactly one terminal record")
        if authority_present:
            return verify_authority_lane(tree)
        return verify_failure_lane(tree)


def _safe_error(error: BaseException) -> str:
    raw = str(error) or error.__class__.__name__
    flattened = " ".join(raw.replace("\r", " ").replace("\n", " ").split())
    safe = "".join(
        character if 0x20 <= ord(character) <= 0x7E else "?"
        for character in flattened)
    return (safe or "invalid evidence")[:4096]


def _silence_broken_stdout() -> None:
    try:
        descriptor = os.open(os.devnull, os.O_WRONLY | os.O_CLOEXEC)
        try:
            os.dup2(descriptor, sys.stdout.fileno())
        finally:
            os.close(descriptor)
    except (AttributeError, OSError, ValueError):
        pass


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 1 or arguments[0] in ("-h", "--help"):
        sys.stderr.write(
            "usage: leopard2_exact_main_baseline_verifier.py <lane-root>\n")
        return 2
    try:
        verdict = verify_sealed_lane(arguments[0])
        output = canonical_json_bytes(verdict)
        sys.stdout.buffer.write(output)
        sys.stdout.buffer.flush()
    except Exception as error:  # Every malformed-evidence exception fails closed.
        if isinstance(error, BrokenPipeError):
            _silence_broken_stdout()
        sys.stderr.write(
            "invalid sealed exact-main baseline: " + _safe_error(error) + "\n")
        return 1
    return 0 if verdict["promoted"] is True else 3


__all__ = (
    "BUILD_CLOSURE_SCHEMA",
    "LANE_DIRECTORY_MODE",
    "LANE_FILE_MODE",
    "SEAL_PROTOCOL",
    "TREE_METADATA_SCHEMA",
    "VERIFIER_SCHEMA",
    "SealedLaneError",
    "SealedTree",
    "main",
    "read_sealed_tree",
    "verify_authority_lane",
    "verify_failure_lane",
    "verify_sealed_lane",
    "verify_sha256sums",
    "verify_tree_metadata",
)


if __name__ == "__main__":
    raise SystemExit(main())
