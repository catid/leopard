#!/usr/bin/env python3
"""Simultaneous, fail-closed Git source capture for benchmark evidence.

The benchmark runners must not assemble a source identity from unrelated Git
processes.  This module retains one exact Git executable, executes a sealed
copy of it for every query, holds the worktree and all ordinary/linked-worktree
metadata throughout capture, and verifies the complete query series twice.
Tracked bytes are opened through retained directory descriptors and remain
open until the Git state and pathnames have been revalidated.
"""

from __future__ import annotations

import base64
import binascii
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import struct
import sys
from typing import Any, Mapping, Sequence


TOOLS_DIRECTORY = Path(__file__).resolve().parents[3] / "tools"
if str(TOOLS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIRECTORY))
import leopard2_build_provenance as _build_provenance  # noqa: E402
_BUILD_PROVENANCE_PATH = (
    TOOLS_DIRECTORY / "leopard2_build_provenance.py").resolve(strict=True)
if Path(_build_provenance.__file__).resolve(strict=True) != \
        _BUILD_PROVENANCE_PATH:
    raise RuntimeError(
        "Git capture provenance primitives resolved outside this source tree")
from leopard2_build_provenance import (  # noqa: E402
    BuildProvenanceError,
    MAX_METADATA_BYTES,
    MAX_TRACKED_SOURCE_FILES,
    MAX_TRACKED_SOURCE_FILE_BYTES,
    MAX_TRACKED_SOURCE_PATH_BYTES,
    MAX_TRACKED_SOURCE_TOTAL_BYTES,
    _InotifyMutationGuard,
    _OpenDirectoryTree,
    _OwnedDescriptor,
    _RetainedFileSnapshot,
    _RetainedGitMetadata,
    _read_bounded_descriptor,
    _require_safe_unicode,
    _run,
    _stable_fields,
)


SCHEMA_V1 = "leopard2-git-source-capture/v1"
SCHEMA = "leopard2-git-source-capture/v2"
DIRECT_SCHEMA = "leopard2-direct-git-source-capture/v1"
GIT_EXECUTABLE_SOURCE_PATH = "/usr/bin/git"
GIT_EXECUTABLE_PROTOCOL = "linux-sealed-git-executable-memfd/v1"
DIRECT_INDEX_PROTOCOL = "linux-sealed-git-index-memfd/v1"
METADATA_GUARD_POLICY = \
    "retained-gitdir-commondir-recursive-inotify/v1"
WORKTREE_GUARD_POLICY = \
    "retained-root-and-tracked-path-components-inotify/v1"
DIRECT_WORKTREE_GUARD_POLICY = \
    "retained-recursive-topology-ignore-and-head-index-union-inotify/v3"
MAX_GIT_COMMIT_BYTES = 16 * 1024 * 1024
MAX_GIT_TREE_BYTES = 64 * 1024 * 1024
MAX_GIT_TREE_OBJECTS = 4 * MAX_TRACKED_SOURCE_FILES
MAX_GIT_TREE_TOTAL_BYTES = MAX_GIT_TREE_BYTES
MAX_GIT_LISTING_BYTES = 64 * 1024 * 1024
MAX_GIT_CONFIG_BYTES = 16 * 1024 * 1024
MAX_SUBMODULE_DEPTH = 16
# A productive expanded tree visit can be charged to one component of a
# retained leaf path.  Empty subtrees have no such leaf, so allow at most one
# additional visit per retained tree object before treating repeated empty-DAG
# expansion as an unsupported resource-amplification dialect.
MAX_GIT_TREE_EXPANSION_VISITS = (
    MAX_GIT_TREE_OBJECTS +
    MAX_TRACKED_SOURCE_FILES * ((MAX_TRACKED_SOURCE_PATH_BYTES + 1) // 2)
)
# Materializing recursive pathnames necessarily copies their bytes.  Admit the
# complete ordinary-tree envelope (one path per non-root tree plus every leaf),
# but reject compact shared DAGs that multiply those occurrences far beyond the
# declared object and leaf inventories.
MAX_GIT_TREE_EXPANDED_PATH_BYTES = (
    (MAX_GIT_TREE_OBJECTS + MAX_TRACKED_SOURCE_FILES) *
    MAX_TRACKED_SOURCE_PATH_BYTES
)
HEX40 = re.compile(r"[0-9a-f]{40}")
HEX256 = re.compile(r"[0-9a-f]{64}")
REQUIRED_SEALS = 0x0001 | 0x0002 | 0x0004 | 0x0008


class GitCaptureError(RuntimeError):
    """The repository could not be captured as one coherent source state."""


def _raise_walk_error(error: OSError) -> None:
    """Make metadata enumeration fail closed on ``scandir`` errors."""
    raise error


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GitCaptureError(message)


def _canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError, RecursionError) as error:
        raise GitCaptureError(
            f"Git capture identity is not canonical JSON: {error}") from error


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _sha1_object(kind: str, content: bytes) -> str:
    return hashlib.sha1(
        f"{kind} {len(content)}\0".encode("ascii") + content,
        usedforsecurity=False).hexdigest()


def _portable_file_identity(snapshot: _RetainedFileSnapshot) -> dict[str, Any]:
    identity = snapshot.identity
    return {
        "path": str(snapshot.resolved),
        "size": identity["size"],
        "mode": stat.S_IMODE(identity["mode"]),
        "sha256": identity["sha256"],
    }


def _byte_identity(content: bytes) -> dict[str, Any]:
    return {
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def _byte_payload(content: bytes) -> dict[str, Any]:
    """Retain exact bounded command output for offline source validation."""
    return {
        "encoding": "base64",
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "base64": base64.b64encode(content).decode("ascii"),
    }


def _object_identity(kind: str, content: bytes) -> dict[str, Any]:
    return {
        "encoding": "base64",
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "object_id": _sha1_object(kind, content),
        "base64": base64.b64encode(content).decode("ascii"),
    }


def _decode_text(raw: bytes, label: str, *, strip_lf: bool = True) -> str:
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise GitCaptureError(f"{label} is not strict UTF-8") from error
    if strip_lf:
        _require(text.endswith("\n") and text.count("\n") == 1,
                 f"{label} is not one canonical LF-terminated record")
        text = text[:-1]
    _require_safe_unicode(text, label)
    return text


def _safe_path(raw: bytes, label: str) -> str:
    _require(0 < len(raw) <= MAX_TRACKED_SOURCE_PATH_BYTES,
             f"{label} is empty or exceeds its byte bound")
    try:
        path = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise GitCaptureError(f"{label} is not strict UTF-8") from error
    _require_safe_unicode(path, label)
    _require(
        not path.startswith("/") and
        all(component not in ("", ".", "..") for component in path.split("/")),
        f"{label} is not a safe relative pathname")
    return path


def _parse_tree_listing(raw: bytes) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    if raw:
        _require(raw.endswith(b"\0"),
                 "Git tree listing is not NUL terminated")
    encoded_records = [] if not raw else raw[:-1].split(b"\0")
    _require(all(encoded_records),
             "Git tree listing contains an empty record")
    for encoded in encoded_records:
        header, separator, path_raw = encoded.partition(b"\t")
        _require(bool(separator), "Git tree listing record has no pathname")
        fields = header.split(b" ")
        _require(len(fields) == 3,
                 "Git tree listing record has an invalid header")
        try:
            mode = fields[0].decode("ascii")
            kind = fields[1].decode("ascii")
            object_id = fields[2].decode("ascii")
        except UnicodeDecodeError as error:
            raise GitCaptureError(
                "Git tree listing header is not ASCII") from error
        expected_kind = (
            "blob" if mode in ("100644", "100755", "120000") else
            "commit" if mode == "160000" else None)
        _require(expected_kind is not None and kind == expected_kind and
                 HEX40.fullmatch(object_id) is not None,
                 "Git tree listing contains an unsupported mode, type, or ID")
        entries.append({
            "path": _safe_path(path_raw, "tracked source path"),
            "git_mode": mode,
            "git_type": kind,
            "object_id": object_id,
        })
    paths = [entry["path"] for entry in entries]
    _require(len(entries) <= MAX_TRACKED_SOURCE_FILES and
             len(paths) == len(set(paths)),
             "Git tree listing is oversized or contains duplicate paths")
    _require(paths == sorted(paths),
             "Git tree listing is not in canonical pathname order")
    return entries


def _parse_tree_object(raw: bytes, label: str) -> list[dict[str, str]]:
    """Parse one canonical Git tree object, including subtree entries."""
    entries: list[dict[str, str]] = []
    offset = 0
    while offset < len(raw):
        space = raw.find(b" ", offset)
        nul = raw.find(b"\0", space + 1 if space >= 0 else offset)
        _require(space > offset and nul > space + 1 and
                 nul + 21 <= len(raw),
                 f"{label} contains a truncated tree entry")
        mode_raw = raw[offset:space]
        name_raw = raw[space + 1:nul]
        object_raw = raw[nul + 1:nul + 21]
        offset = nul + 21
        try:
            mode = mode_raw.decode("ascii", errors="strict")
        except UnicodeDecodeError as error:
            raise GitCaptureError(f"{label} mode is not ASCII") from error
        expected_kind = (
            "tree" if mode == "40000" else
            "blob" if mode in ("100644", "100755", "120000") else
            "commit" if mode == "160000" else None)
        _require(expected_kind is not None and
                 b"/" not in name_raw and name_raw not in (b".", b".."),
                 f"{label} contains an unsupported mode or pathname")
        name = _safe_path(name_raw, f"{label} entry name")
        entries.append({
            "name": name,
            "git_mode": mode,
            "git_type": expected_kind,
            "object_id": object_raw.hex(),
        })
        _require(
            len(entries) <=
                MAX_TRACKED_SOURCE_FILES + MAX_GIT_TREE_OBJECTS,
            f"{label} contains too many entries")
    names = [entry["name"] for entry in entries]
    _require(len(names) == len(set(names)),
             f"{label} contains duplicate entry names")
    sort_keys = [
        entry["name"].encode("utf-8") +
        (b"/" if entry["git_type"] == "tree" else b"")
        for entry in entries
    ]
    _require(sort_keys == sorted(sort_keys) and
             len(sort_keys) == len(set(sort_keys)),
             f"{label} entries are duplicated or not in Git tree order")
    return entries


def _capture_tree_objects(
    git_snapshot: _RetainedFileSnapshot,
    metadata: _RetainedGitMetadata,
    root_descriptor: int,
    root_tree: str,
    root_content: bytes,
    inherited_descriptors: Sequence[int],
) -> list[dict[str, Any]]:
    """Retain the complete reachable tree-object closure exactly once."""
    _require(_sha1_object("tree", root_content) == root_tree,
             "retained Git root-tree bytes differ")
    contents: dict[str, bytes] = {root_tree: root_content}
    pending = [root_tree]
    queued = {root_tree}
    total = 0
    while pending:
        object_id = pending.pop()
        content = contents.get(object_id)
        if content is None:
            content = _invoke_git(
                git_snapshot, metadata, root_descriptor,
                ("cat-file", "tree", object_id),
                f"Git tree object {object_id}",
                maximum_bytes=MAX_GIT_TREE_BYTES,
                inherited_descriptors=inherited_descriptors)
            _require(_sha1_object("tree", content) == object_id,
                     "retained Git subtree bytes differ from their object ID")
            contents[object_id] = content
        total += len(content)
        _require(len(contents) <= MAX_GIT_TREE_OBJECTS and
                 total <= MAX_GIT_TREE_TOTAL_BYTES,
                 "retained Git tree-object closure exceeds its bound")
        for entry in _parse_tree_object(
                content, f"Git tree object {object_id}"):
            if entry["git_type"] != "tree":
                continue
            child = entry["object_id"]
            if child not in queued:
                queued.add(child)
                _require(len(queued) <= MAX_GIT_TREE_OBJECTS,
                         "retained Git tree-object closure has too many trees")
                pending.append(child)
    return [
        _object_identity("tree", contents[object_id])
        for object_id in sorted(contents)
    ]


def _flatten_tree_objects(
    root_tree: str, objects: Mapping[str, bytes],
) -> list[dict[str, str]]:
    """Derive the exact recursive ls-tree leaf inventory from retained bytes."""
    _require(root_tree in objects,
             "retained Git tree-object closure omits its root")

    # Analyze each distinct object once before expanding path occurrences.
    # A path-local cycle set alone is insufficient here: an acyclic object DAG
    # can name the same empty child twice at every level and cause exponential
    # work without ever producing a leaf for the existing leaf-count bound.
    parsed: dict[str, list[dict[str, str]]] = {}
    reachable: set[str] = set()
    active: set[str] = set()
    # object -> (expanded leaf count, expanded tree visits, longest suffix,
    #            materialized path count, materialized path bytes)
    analysis: dict[str, tuple[int, int, int, int, int]] = {}
    pending: list[tuple[str, bool]] = [(root_tree, False)]
    while pending:
        object_id, exiting = pending.pop()
        if exiting:
            entries = parsed[object_id]
            leaf_count = 0
            expansion_visits = 1
            longest_suffix_bytes = 0
            expanded_path_count = 0
            expanded_path_bytes = 0
            for entry in entries:
                name_bytes = len(entry["name"].encode("utf-8"))
                if entry["git_type"] == "tree":
                    child = analysis.get(entry["object_id"])
                    _require(child is not None,
                             "retained Git tree-object analysis is incomplete")
                    (child_leaves, child_visits, child_suffix_bytes,
                     child_path_count, child_path_bytes) = child
                    leaf_count += child_leaves
                    expansion_visits = min(
                        MAX_GIT_TREE_EXPANSION_VISITS + 1,
                        expansion_visits + child_visits)
                    suffix_bytes = (
                        name_bytes if child_suffix_bytes == 0 else
                        name_bytes + 1 + child_suffix_bytes)
                    if child_leaves:
                        expanded_path_count = min(
                            MAX_GIT_TREE_EXPANSION_VISITS +
                                MAX_TRACKED_SOURCE_FILES + 1,
                            expanded_path_count + 1 + child_path_count)
                        expanded_path_bytes = min(
                            MAX_GIT_TREE_EXPANDED_PATH_BYTES + 1,
                            expanded_path_bytes + name_bytes +
                                child_path_bytes +
                                child_path_count * (name_bytes + 1))
                else:
                    leaf_count += 1
                    suffix_bytes = name_bytes
                    expanded_path_count += 1
                    expanded_path_bytes = min(
                        MAX_GIT_TREE_EXPANDED_PATH_BYTES + 1,
                        expanded_path_bytes + name_bytes)
                _require(leaf_count <= MAX_TRACKED_SOURCE_FILES,
                         "retained Git recursive tree has too many leaves")
                _require(
                    suffix_bytes <= MAX_TRACKED_SOURCE_PATH_BYTES,
                    "retained Git recursive tree path exceeds its byte bound")
                longest_suffix_bytes = max(
                    longest_suffix_bytes, suffix_bytes)
            analysis[object_id] = (
                leaf_count, expansion_visits, longest_suffix_bytes,
                expanded_path_count, expanded_path_bytes)
            active.remove(object_id)
            continue

        if object_id in analysis:
            continue
        _require(object_id not in active and object_id in objects,
                 "retained Git tree-object closure is cyclic or incomplete")
        active.add(object_id)
        reachable.add(object_id)
        entries = _parse_tree_object(
            objects[object_id], f"retained Git tree object {object_id}")
        parsed[object_id] = entries
        pending.append((object_id, True))
        scheduled_children: set[str] = set()
        for entry in reversed(entries):
            if entry["git_type"] != "tree":
                continue
            child = entry["object_id"]
            _require(child not in active and child in objects,
                     "retained Git tree-object closure is cyclic or incomplete")
            if child not in analysis and child not in scheduled_children:
                scheduled_children.add(child)
                pending.append((child, False))

    _require(reachable == set(objects),
             "retained Git tree-object closure contains unreachable objects")
    (root_leaves, root_visits, _root_suffix_bytes,
     root_path_count, root_path_bytes) = analysis[root_tree]
    _require(
        root_visits <= MAX_GIT_TREE_EXPANSION_VISITS,
        "retained Git tree-object closure expands beyond its traversal bound")
    _require(
        root_path_bytes <= MAX_GIT_TREE_EXPANDED_PATH_BYTES,
        "retained Git tree-object closure exceeds its expanded path-byte bound")

    result: list[dict[str, str]] = []
    expanded_visits = 0
    materialized_path_count = 0
    materialized_path_bytes = 0
    stack: list[tuple[str, object, str, int]] = [
        ("tree", root_tree, "", 0)]
    while stack:
        kind, payload, prefix, prefix_bytes = stack.pop()
        if kind == "leaf":
            entry = payload
            _require(isinstance(entry, dict),
                     "retained Git recursive tree leaf is invalid")
            result.append({
                "path": prefix,
                "git_mode": entry["git_mode"],
                "git_type": entry["git_type"],
                "object_id": entry["object_id"],
            })
            _require(len(result) <= MAX_TRACKED_SOURCE_FILES,
                     "retained Git recursive tree has too many leaves")
            continue
        object_id = payload
        _require(isinstance(object_id, str) and object_id in parsed,
                 "retained Git recursive tree traversal is invalid")
        expanded_visits += 1
        _require(
            expanded_visits <= MAX_GIT_TREE_EXPANSION_VISITS,
            "retained Git tree-object closure expands beyond its traversal "
            "bound")
        for entry in reversed(parsed[object_id]):
            if entry["git_type"] == "tree":
                # Empty subtrees contribute no recursive ls-tree leaf.  Their
                # complete reachability, path length, cycles, and algebraic
                # expansion were already checked once in the analysis pass.
                if analysis[entry["object_id"]][0] == 0:
                    continue
            name_bytes = len(entry["name"].encode("utf-8"))
            relative_bytes = (
                prefix_bytes + 1 + name_bytes if prefix else name_bytes)
            _require(
                relative_bytes <= MAX_TRACKED_SOURCE_PATH_BYTES,
                "retained Git recursive tree path exceeds its byte bound")
            # Each component was already decoded, category-checked, and
            # rejected if empty, '.', '..', or slash-bearing by
            # _parse_tree_object.  Concatenating those components with one '/'
            # cannot introduce a new unsafe character or component, so only
            # the O(1) cumulative byte-length check is needed here.
            relative = (
                f"{prefix}/{entry['name']}" if prefix else entry["name"])
            materialized_path_count += 1
            materialized_path_bytes += relative_bytes
            _require(
                materialized_path_bytes <=
                    MAX_GIT_TREE_EXPANDED_PATH_BYTES,
                "retained Git tree-object closure exceeds its expanded "
                "path-byte bound")
            if entry["git_type"] == "tree":
                stack.append((
                    "tree", entry["object_id"], relative, relative_bytes))
            else:
                stack.append(("leaf", entry, relative, relative_bytes))
    _require(len(result) == root_leaves and
             materialized_path_count == root_path_count and
             materialized_path_bytes == root_path_bytes,
             "retained Git recursive tree analysis differs")
    paths = [entry["path"] for entry in result]
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "retained Git recursive tree inventory is non-canonical")
    return result


def _parse_index(raw: bytes) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    if raw:
        _require(raw.endswith(b"\0"),
                 "Git index-stage output is not NUL terminated")
    encoded_records = [] if not raw else raw[:-1].split(b"\0")
    _require(all(encoded_records),
             "Git index-stage output contains an empty record")
    for encoded in encoded_records:
        header, separator, path_raw = encoded.partition(b"\t")
        fields = header.split(b" ")
        _require(bool(separator) and len(fields) == 3,
                 "Git index-stage record is malformed")
        try:
            mode = fields[0].decode("ascii")
            object_id = fields[1].decode("ascii")
            stage = fields[2].decode("ascii")
        except UnicodeDecodeError as error:
            raise GitCaptureError("Git index-stage record is not ASCII") \
                from error
        _require(mode in ("100644", "100755", "120000", "160000") and
                 HEX40.fullmatch(object_id) is not None and stage == "0",
                 "Git index contains an unsupported mode, ID, or stage")
        entries.append({
            "path": _safe_path(path_raw, "Git index path"),
            "git_mode": mode,
            "object_id": object_id,
        })
    paths = [entry["path"] for entry in entries]
    _require(len(entries) <= MAX_TRACKED_SOURCE_FILES and
             len(paths) == len(set(paths)) and paths == sorted(paths),
             "Git index is oversized, duplicated, or non-canonical")
    return entries


def _parse_raw_index(raw: bytes) -> list[dict[str, str]]:
    """Parse a bounded SHA-1 Git index and return its stage-zero material.

    Direct evidence retains the exact index bytes used through
    ``GIT_INDEX_FILE``.  Parsing the ordinary v2/v3 on-disk format here binds
    that immutable file to the separately retained ``ls-files`` transcripts
    without invoking Git during offline validation.  Split, sparse, extended,
    conflicted, or otherwise non-default entries are rejected by this dialect.
    """
    _require(32 <= len(raw) <= MAX_GIT_LISTING_BYTES,
             "raw direct Git index is empty or exceeds its byte bound")
    _require(raw[:4] == b"DIRC",
             "raw direct Git index has an invalid signature")
    version, count = struct.unpack(">II", raw[4:12])
    _require(version in (2, 3) and count <= MAX_TRACKED_SOURCE_FILES,
             "raw direct Git index version or entry count is unsupported")
    payload_end = len(raw) - 20
    _require(hashlib.sha1(
                 raw[:payload_end], usedforsecurity=False).digest() ==
             raw[payload_end:],
             "raw direct Git index checksum differs")
    entries: list[dict[str, str]] = []
    offset = 12
    for _unused in range(count):
        entry_start = offset
        _require(offset + 62 <= payload_end,
                 "raw direct Git index entry is truncated")
        mode = struct.unpack(">I", raw[offset + 24:offset + 28])[0]
        object_id = raw[offset + 40:offset + 60].hex()
        flags = struct.unpack(">H", raw[offset + 60:offset + 62])[0]
        offset += 62
        _require(not flags & 0x8000 and not flags & 0x4000,
                 "raw direct Git index uses non-default entry flags")
        stage = (flags >> 12) & 0x3
        declared_name_length = flags & 0x0fff
        nul = raw.find(b"\0", offset, payload_end)
        _require(nul >= offset,
                 "raw direct Git index pathname is unterminated")
        path_raw = raw[offset:nul]
        _require(
            stage == 0 and
            (declared_name_length == 0x0fff or
             declared_name_length == len(path_raw)),
            "raw direct Git index stage or pathname length differs")
        mode_text = format(mode, "o")
        _require(mode_text in ("100644", "100755", "120000", "160000") and
                 HEX40.fullmatch(object_id) is not None,
                 "raw direct Git index mode or object ID is unsupported")
        path = _safe_path(path_raw, "raw direct Git index path")
        entries.append({
            "path": path,
            "git_mode": mode_text,
            "object_id": object_id,
        })
        used = nul + 1 - entry_start
        padded = (used + 7) & ~7
        offset = entry_start + padded
        _require(offset <= payload_end and
                 raw[nul + 1:offset] == b"\0" * (offset - nul - 1),
                 "raw direct Git index padding differs")

    # Optional extensions have an uppercase first byte.  Lowercase extensions
    # affect required semantics (including split/sparse index) and are outside
    # this proof dialect.  Every extension remains authenticated by the raw
    # index checksum and byte payload even when its optional contents are not
    # otherwise interpreted here.
    while offset < payload_end:
        _require(offset + 8 <= payload_end,
                 "raw direct Git index extension is truncated")
        signature = raw[offset:offset + 4]
        extension_size = struct.unpack(">I", raw[offset + 4:offset + 8])[0]
        _require(signature[:1].isalpha() and signature[:1].isupper(),
                 "raw direct Git index has an unsupported mandatory extension")
        offset += 8
        _require(extension_size <= payload_end - offset,
                 "raw direct Git index extension exceeds the payload")
        offset += extension_size
    _require(offset == payload_end,
             "raw direct Git index payload has trailing bytes")
    paths = [entry["path"] for entry in entries]
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "raw direct Git index paths are non-canonical")
    return entries


def _parse_default_flags(raw: bytes, label: str) -> list[str]:
    paths: list[str] = []
    if raw:
        _require(raw.endswith(b"\0"), f"{label} is not NUL terminated")
    encoded_records = [] if not raw else raw[:-1].split(b"\0")
    _require(all(encoded_records), f"{label} contains an empty record")
    for encoded in encoded_records:
        _require(encoded.startswith(b"H "),
                 "Git index uses assume-unchanged, skip-worktree, "
                 "fsmonitor-valid, or another non-default flag")
        paths.append(_safe_path(encoded[2:], label))
    _require(len(paths) == len(set(paths)) and paths == sorted(paths),
             f"{label} is duplicated or non-canonical")
    return paths


def _parse_direct_refs(raw: bytes) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    if raw:
        _require(raw.endswith(b"\n") and b"\r" not in raw,
                 "direct Git reference inventory is not LF terminated")
    lines = [] if not raw else raw[:-1].split(b"\n")
    _require(all(lines),
             "direct Git reference inventory contains an empty record")
    for line in lines:
        fields = line.split(b" ")
        _require(len(fields) == 2,
                 "direct Git reference record is malformed")
        try:
            name = fields[0].decode("utf-8", errors="strict")
            object_id = fields[1].decode("ascii", errors="strict")
        except UnicodeDecodeError as error:
            raise GitCaptureError(
                "direct Git reference record is not UTF-8/ASCII") from error
        _require(name.startswith("refs/") and
                 HEX40.fullmatch(object_id) is not None,
                 "direct Git reference name or object ID is invalid")
        _require_safe_unicode(name, "direct Git reference")
        records.append({"name": name, "object_id": object_id})
    names = [record["name"] for record in records]
    _require(names == sorted(names) and len(names) == len(set(names)),
             "direct Git reference inventory is non-canonical")
    return records


def _validate_local_config(raw: bytes) -> None:
    if raw:
        _require(raw.endswith(b"\0"),
                 "local Git configuration is not NUL terminated")
    records = [] if not raw else raw[:-1].split(b"\0")
    _require(all(records),
             "local Git configuration contains an empty record")
    for record in records:
        key, separator, _value = record.partition(b"\n")
        _require(bool(separator) and key,
                 "local Git configuration record is malformed")
        try:
            name = key.decode("ascii", errors="strict").lower()
        except UnicodeDecodeError as error:
            raise GitCaptureError(
                "local Git configuration key is not ASCII") from error
        _require(not name.startswith(("include.", "includeif.")),
                 "local Git configuration contains an external include")


def _invoke_git(
    git_snapshot: _RetainedFileSnapshot,
    metadata: _RetainedGitMetadata,
    root_descriptor: int,
    arguments: Sequence[str],
    label: str,
    *,
    maximum_bytes: int = MAX_METADATA_BYTES,
    inherited_descriptors: Sequence[int] = (),
    index_snapshot: _RetainedFileSnapshot | None = None,
) -> bytes:
    """Run one patchable query through the exact sealed Git executable."""
    _require(git_snapshot.resolved is not None,
             "retained Git executable pathname was lost")
    metadata.verify()
    git_snapshot.verify()
    # Git otherwise consults the invoking account's default global excludes
    # file even when the global/system configuration files are disabled.
    # /dev/null is a stable non-directory, so neither HOME nor XDG can name an
    # external ignore file outside the retained repository snapshot.
    environment_overrides = {
        "HOME": "/dev/null",
        "XDG_CONFIG_HOME": "/dev/null",
    }
    index_descriptor = None
    if index_snapshot is not None:
        index_snapshot.verify()
        index_descriptor = index_snapshot.executable_descriptor
        environment_overrides["GIT_INDEX_FILE"] = \
            f"/proc/self/fd/{index_descriptor}"
    inherited = (
        *inherited_descriptors, root_descriptor, metadata.descriptor,
        *((index_descriptor,) if index_descriptor is not None else ()))
    output = _run(
        (
            str(git_snapshot.resolved),
            "-c", "core.fsmonitor=false",
            "-c", "core.untrackedCache=false",
            f"--git-dir=/proc/self/fd/{metadata.descriptor}",
            f"--work-tree=/proc/self/fd/{root_descriptor}",
            *arguments,
        ),
        label,
        maximum_bytes=maximum_bytes,
        inherited_descriptors=inherited,
        executable_descriptor=git_snapshot.executable_descriptor,
        environment_overrides=environment_overrides,
    )
    metadata.verify()
    git_snapshot.verify()
    if index_snapshot is not None:
        index_snapshot.verify()
    return output


def _guard_metadata_files(
    metadata: _RetainedGitMetadata,
) -> list[str]:
    """Watch every metadata-file inode, including writes through hardlinks."""
    gitdir = metadata.gitdir
    common = metadata.common
    _require(gitdir is not None and gitdir.resolved is not None and
             gitdir.guard is not None and common is not None and
             common.resolved is not None and common.guard is not None,
             "retained Git metadata tree is unavailable")
    roots: list[tuple[str, Path, _InotifyMutationGuard]] = [
        ("gitdir", gitdir.resolved, gitdir.guard),
    ]
    if common is not gitdir:
        roots.append(("commondir", common.resolved, common.guard))
    records: list[str] = []
    for role, root, guard in roots:
        try:
            for directory, child_directories, files in os.walk(
                    root, topdown=True, followlinks=False,
                    onerror=_raise_walk_error):
                child_directories.sort()
                files.sort()
                parent = Path(directory)
                for name in files:
                    path = parent / name
                    status = path.lstat()
                    _require(stat.S_ISREG(status.st_mode),
                             f"Git metadata leaf {path} is not a regular file")
                    guard.add_file_path(path)
                    relative = path.relative_to(root).as_posix()
                    _require_safe_unicode(
                        relative, "retained Git metadata path")
                    records.append(f"{role}:{relative}")
        except OSError as error:
            raise GitCaptureError(
                f"cannot enumerate retained Git metadata files: {error}") \
                from error
    records.sort()
    _require(len(records) == len(set(records)),
             "retained Git metadata file inventory is non-canonical")
    metadata.verify()
    return records


def _require_direct_grafts_absent(metadata: _RetainedGitMetadata) -> None:
    """Reject legacy graft files under every retained metadata root."""
    gitdir = metadata.gitdir
    common = metadata.common
    _require(gitdir is not None and gitdir.resolved is not None and
             common is not None and common.resolved is not None,
             "retained direct Git metadata root was lost")
    for root in {gitdir.resolved, common.resolved}:
        grafts = root / "info" / "grafts"
        _require(not grafts.exists() and not grafts.is_symlink(),
                 f"legacy Git grafts are forbidden: {grafts}")
    metadata.verify()


_SERIES = (
    ("top", ("rev-parse", "--show-toplevel"), MAX_TRACKED_SOURCE_PATH_BYTES),
    ("replace_refs",
     ("for-each-ref", "--format=%(refname)", "refs/replace"),
     MAX_METADATA_BYTES),
    ("head", ("rev-parse", "--verify", "HEAD"), 256),
    ("tree", ("rev-parse", "--verify", "HEAD^{tree}"), 256),
    ("head_ref", ("rev-parse", "--symbolic-full-name", "HEAD"), 65536),
    ("superproject",
     ("rev-parse", "--show-superproject-working-tree"),
     MAX_TRACKED_SOURCE_PATH_BYTES),
    ("commit_object", ("cat-file", "commit", "{head}"),
     MAX_GIT_COMMIT_BYTES),
    ("tree_object", ("cat-file", "tree", "{tree}"), MAX_GIT_TREE_BYTES),
    ("tree_listing", ("ls-tree", "-r", "-z", "--full-tree", "HEAD"),
     MAX_GIT_LISTING_BYTES),
    ("index_stage", ("ls-files", "--stage", "-z"), MAX_GIT_LISTING_BYTES),
    ("index_v", ("ls-files", "-v", "-z"), MAX_GIT_LISTING_BYTES),
    ("index_f", ("ls-files", "-f", "-z"), MAX_GIT_LISTING_BYTES),
    ("status",
     ("status", "--porcelain=v1", "-z", "--untracked-files=normal",
      "--ignore-submodules=none"),
     MAX_GIT_LISTING_BYTES),
    ("config",
     # Global and system configuration are disabled by _run's fixed
     # environment.  Omitting --local deliberately includes config.worktree
     # for linked worktrees, while --no-includes keeps every consumed byte
     # inside the retained gitdir/commondir closure.
     ("config", "--no-includes", "--null", "--list"),
     MAX_GIT_CONFIG_BYTES),
)


def _query_series(
    git_snapshot: _RetainedFileSnapshot,
    metadata: _RetainedGitMetadata,
    root_descriptor: int,
    inherited_descriptors: Sequence[int],
) -> dict[str, bytes]:
    result: dict[str, bytes] = {}
    for name, template, maximum in _SERIES:
        arguments = tuple(
            result["head"].decode("ascii").strip() if item == "{head}" else
            result["tree"].decode("ascii").strip() if item == "{tree}" else
            item
            for item in template)
        result[name] = _invoke_git(
            git_snapshot, metadata, root_descriptor, arguments,
            f"Git source {name}", maximum_bytes=maximum,
            inherited_descriptors=inherited_descriptors)
    return result


def _query_direct_series(
    git_snapshot: _RetainedFileSnapshot,
    metadata: _RetainedGitMetadata,
    index_snapshot: _RetainedFileSnapshot,
    root_descriptor: int,
    inherited_descriptors: Sequence[int],
) -> dict[str, bytes]:
    """Query the direct runner's clean-or-dirty source-state dialect."""
    result: dict[str, bytes] = {}
    for name, template, maximum in _SERIES:
        if name == "status":
            # Direct dirtiness is derived from exact retained layers below.
            # A read-only sealed index can preserve stale stat-cache fields,
            # making porcelain report a false modification even when the
            # retained blob and mode are byte-for-byte identical.
            continue
        arguments = tuple(
            result["head"].decode("ascii").strip() if item == "{head}" else
            result["tree"].decode("ascii").strip() if item == "{tree}" else
            item
            for item in template)
        result[name] = _invoke_git(
            git_snapshot, metadata, root_descriptor, arguments,
            f"direct Git source {name}", maximum_bytes=maximum,
            inherited_descriptors=inherited_descriptors,
            index_snapshot=index_snapshot)
    result["shared_index"] = _invoke_git(
        git_snapshot, metadata, root_descriptor,
        ("rev-parse", "--shared-index-path"),
        "direct Git shared-index path", maximum_bytes=MAX_METADATA_BYTES,
        inherited_descriptors=inherited_descriptors,
        index_snapshot=index_snapshot)
    result["refs"] = _invoke_git(
        git_snapshot, metadata, root_descriptor,
        ("for-each-ref", "--sort=refname",
         "--format=%(refname) %(objectname)"),
        "direct Git reference inventory", maximum_bytes=MAX_METADATA_BYTES,
        inherited_descriptors=inherited_descriptors,
        index_snapshot=index_snapshot)
    result["untracked_paths"] = _invoke_git(
        git_snapshot, metadata, root_descriptor,
        ("ls-files", "--others", "--exclude-standard", "-z"),
        "direct Git untracked path inventory",
        maximum_bytes=MAX_GIT_LISTING_BYTES,
        inherited_descriptors=inherited_descriptors,
        index_snapshot=index_snapshot)
    return result


def _verify_query_series(
    root: Path, requested_commit: str, series: Mapping[str, bytes],
    preliminary_listing: bytes,
) -> tuple[str, str, bool, str | None, str | None,
           list[dict[str, str]]]:
    top = Path(_decode_text(series["top"], "Git top level"))
    try:
        top = top.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise GitCaptureError(f"Git top level cannot be resolved: {error}") \
            from error
    _require(top == root,
             f"source root is not the Git top level: {root} != {top}")
    _require(series["replace_refs"] == b"",
             "source repository contains Git replace refs")
    head = _decode_text(series["head"], "Git HEAD")
    tree_id = _decode_text(series["tree"], "Git root tree")
    _require(HEX40.fullmatch(head) is not None and
             HEX40.fullmatch(tree_id) is not None,
             "Git HEAD or root tree is not one lowercase SHA-1 ID")
    _require(head == requested_commit,
             f"source HEAD mismatch: requested {requested_commit}, got {head}")

    head_name = _decode_text(series["head_ref"], "Git symbolic HEAD")
    detached = head_name == "HEAD"
    head_ref = None if detached else head_name
    if head_ref is not None:
        _require(head_ref.startswith("refs/"),
                 "Git symbolic HEAD is not a canonical reference")
        _require_safe_unicode(head_ref, "Git symbolic HEAD")

    raw_superproject = series["superproject"]
    if raw_superproject:
        superproject = _decode_text(
            raw_superproject, "Git superproject worktree")
        _require(Path(superproject).is_absolute(),
                 "Git superproject worktree is not absolute")
    else:
        superproject = None

    commit_object = series["commit_object"]
    _require(_sha1_object("commit", commit_object) == head,
             "Git commit bytes do not match HEAD")
    _require(b"\n\n" in commit_object,
             "Git commit object has no header/message boundary")
    headers = commit_object.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + tree_id.encode("ascii")
    tree_headers = [line for line in headers if line.startswith(b"tree ")]
    _require(headers and headers[0] == tree_line and
             tree_headers == [tree_line],
             "Git commit object names a different root tree")
    _require(_sha1_object("tree", series["tree_object"]) == tree_id,
             "Git root-tree bytes do not match the retained tree ID")
    _require(series["tree_listing"] == preliminary_listing,
             "Git tree listing changed before the guarded query series")
    tree_entries = _parse_tree_listing(series["tree_listing"])
    index_entries = _parse_index(series["index_stage"])
    _require(
        [(entry["path"], entry["git_mode"], entry["object_id"])
         for entry in index_entries] ==
        [(entry["path"], entry["git_mode"], entry["object_id"])
         for entry in tree_entries],
        "Git index entries do not exactly match the committed tree")
    expected_paths = [entry["path"] for entry in tree_entries]
    _require(_parse_default_flags(series["index_v"], "Git -v index path") ==
             expected_paths and
             _parse_default_flags(series["index_f"], "Git -f index path") ==
             expected_paths,
             "Git index flag inventories do not match the committed tree")
    _require(series["status"] == b"",
             "source has tracked or untracked modifications, or submodule "
             "modifications")
    _validate_local_config(series["config"])
    return head, tree_id, detached, head_ref, superproject, tree_entries


def _verify_direct_query_series(
    root: Path, series: Mapping[str, bytes],
) -> tuple[str, str, bool, str | None, str | None,
           list[dict[str, str]], list[dict[str, str]], bytes]:
    """Cross-bind a direct runner query series without requiring clean state."""
    top = Path(_decode_text(series["top"], "direct Git top level"))
    try:
        top = top.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise GitCaptureError(
            f"direct Git top level cannot be resolved: {error}") from error
    _require(top == root,
             f"source root is not the direct Git top level: {root} != {top}")
    _require(series["replace_refs"] == b"",
             "direct source repository contains Git replace refs")
    _require(series.get("shared_index") == b"",
             "direct source repository uses an unsupported split index")
    head = _decode_text(series["head"], "direct Git HEAD")
    tree_id = _decode_text(series["tree"], "direct Git root tree")
    _require(HEX40.fullmatch(head) is not None and
             HEX40.fullmatch(tree_id) is not None,
             "direct Git HEAD or root tree is not one lowercase SHA-1 ID")

    head_name = _decode_text(
        series["head_ref"], "direct Git symbolic HEAD")
    detached = head_name == "HEAD"
    head_ref = None if detached else head_name
    if head_ref is not None:
        _require(head_ref.startswith("refs/"),
                 "direct Git symbolic HEAD is not a canonical reference")
        _require_safe_unicode(head_ref, "direct Git symbolic HEAD")
    refs = _parse_direct_refs(series.get("refs", b""))
    if head_ref is not None:
        matching = [
            record for record in refs if record["name"] == head_ref]
        _require(len(matching) == 1 and
                 matching[0]["object_id"] == head,
                 "direct Git symbolic HEAD differs from its retained ref")

    raw_superproject = series["superproject"]
    if raw_superproject:
        superproject = _decode_text(
            raw_superproject, "direct Git superproject worktree")
        _require(Path(superproject).is_absolute(),
                 "direct Git superproject worktree is not absolute")
    else:
        superproject = None

    commit_object = series["commit_object"]
    _require(_sha1_object("commit", commit_object) == head,
             "direct Git commit bytes do not match HEAD")
    _require(b"\n\n" in commit_object,
             "direct Git commit object has no header/message boundary")
    headers = commit_object.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + tree_id.encode("ascii")
    tree_headers = [line for line in headers if line.startswith(b"tree ")]
    _require(headers and headers[0] == tree_line and
             tree_headers == [tree_line],
             "direct Git commit object names a different root tree")
    _require(_sha1_object("tree", series["tree_object"]) == tree_id,
             "direct Git root-tree bytes do not match the retained tree ID")
    tree_entries = _parse_tree_listing(series["tree_listing"])
    index_entries = _parse_index(series["index_stage"])
    expected_index_paths = [entry["path"] for entry in index_entries]
    _require(
        _parse_default_flags(
            series["index_v"], "direct Git -v index path") ==
        expected_index_paths and
        _parse_default_flags(
            series["index_f"], "direct Git -f index path") ==
        expected_index_paths,
        "direct Git index flag inventories differ from the staged index")
    _validate_local_config(series["config"])
    return (
        head, tree_id, detached, head_ref, superproject,
        tree_entries, index_entries, series["untracked_paths"])


def _path_status(tree: _OpenDirectoryTree, relative: str) -> os.stat_result:
    parts = tuple(relative.split("/"))
    parent = tree.directory(parts[:-1])
    try:
        return os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
    except OSError as error:
        raise GitCaptureError(
            f"tracked source path {relative!r} cannot be inspected: {error}") \
            from error


def _open_symlink(
    tree: _OpenDirectoryTree, relative: str, stack: ExitStack,
) -> tuple[int, os.stat_result, bytes]:
    parts = tuple(relative.split("/"))
    parent = tree.directory(parts[:-1])
    flags = (
        getattr(os, "O_PATH", os.O_RDONLY) |
        getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NOFOLLOW", 0))
    owner = _OwnedDescriptor()

    # Register ownership before opening the descriptor.  A line-trace,
    # signal, or other asynchronous exception at any later Python boundary
    # will therefore unwind through this callback rather than stranding the
    # raw descriptor between this function's return and caller registration.
    stack.callback(owner.close)
    try:
        descriptor = owner.open(parts[-1], flags, dir_fd=parent)
        status = os.fstat(descriptor)
        target = os.readlink(parts[-1], dir_fd=parent)
        _require(stat.S_ISLNK(status.st_mode),
                 f"tracked path {relative!r} is not a symlink")
    except BaseException as error:
        try:
            owner.close()
        except BaseException as cleanup_error:
            if hasattr(error, "add_note"):
                error.add_note(
                    "tracked symlink descriptor cleanup also failed: "
                    f"{cleanup_error}")
        if isinstance(error, OSError):
            raise GitCaptureError(
                f"tracked symlink {relative!r} cannot be retained: {error}") \
                from error
        raise
    return descriptor, status, os.fsencode(target)


def _capture_repository(
    root: Path,
    requested_commit: str,
    git_snapshot: _RetainedFileSnapshot,
    stack: ExitStack,
    *,
    depth: int,
    inherited_descriptors: Sequence[int],
) -> dict[str, Any]:
    _require(depth <= MAX_SUBMODULE_DEPTH,
             "Git submodule nesting exceeds its recursion bound")
    metadata = stack.enter_context(_RetainedGitMetadata(root))
    metadata_files = _guard_metadata_files(metadata)
    worktree_guard = _InotifyMutationGuard(
        f"tracked worktree {root}")
    stack.callback(worktree_guard.close)
    worktree_guard.add_directory_path(root)
    tree = _OpenDirectoryTree(root, stack, f"tracked worktree {root}")

    preliminary_listing = _invoke_git(
        git_snapshot, metadata, tree.root_descriptor,
        ("ls-tree", "-r", "-z", "--full-tree", "HEAD"),
        "preliminary Git source tree", maximum_bytes=MAX_GIT_LISTING_BYTES,
        inherited_descriptors=inherited_descriptors)
    preliminary_entries = _parse_tree_listing(preliminary_listing)
    for entry in preliminary_entries:
        path = root / entry["path"]
        if entry["git_mode"] == "160000":
            worktree_guard.add_directory_path(path)
        else:
            worktree_guard.add_file_path(path)
    worktree_guard.verify()

    first = _query_series(
        git_snapshot, metadata, tree.root_descriptor, inherited_descriptors)
    (head, tree_id, detached, head_ref, superproject,
     entries) = _verify_query_series(
        root, requested_commit, first, preliminary_listing)
    _require(entries == preliminary_entries,
             "Git tree inventory changed during guarded capture")
    tree_objects = _capture_tree_objects(
        git_snapshot, metadata, tree.root_descriptor, tree_id,
        first["tree_object"], inherited_descriptors)
    tree_object_bytes = {
        record["object_id"]: base64.b64decode(
            record["base64"], validate=True)
        for record in tree_objects
    }
    _require(_flatten_tree_objects(tree_id, tree_object_bytes) == entries,
             "recursive Git tree objects differ from the flattened listing")

    retained_regular: list[
        tuple[str, int, tuple[int, ...], bytes]
    ] = []
    retained_symlinks: list[
        tuple[str, int, tuple[int, ...], bytes]
    ] = []
    records: list[dict[str, Any]] = []
    submodules: list[dict[str, Any]] = []
    total = 0
    for entry in entries:
        relative = entry["path"]
        mode = entry["git_mode"]
        object_id = entry["object_id"]
        if mode in ("100644", "100755"):
            descriptor, before = tree.open_regular(relative)
            _require(before.st_size <= MAX_TRACKED_SOURCE_FILE_BYTES,
                     f"tracked file {relative!r} exceeds its byte bound")
            total += before.st_size
            _require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
                     "tracked source bytes exceed their total bound")
            content = _read_bounded_descriptor(
                descriptor, before.st_size, f"tracked file {relative!r}",
                MAX_TRACKED_SOURCE_FILE_BYTES)
            after = os.fstat(descriptor)
            _require(_stable_fields(before) == _stable_fields(after),
                     f"tracked file {relative!r} changed while read")
            is_executable = bool(stat.S_IMODE(after.st_mode) & 0o111)
            _require(is_executable == (mode == "100755") and
                     _sha1_object("blob", content) == object_id,
                     f"tracked file {relative!r} differs from its Git blob")
            retained_regular.append((
                relative, descriptor, _stable_fields(after), content))
            records.append({
                **entry,
                "kind": "regular",
            })
        elif mode == "120000":
            descriptor, before, target = _open_symlink(
                tree, relative, stack)
            total += len(target)
            _require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES and
                     len(target) <= MAX_TRACKED_SOURCE_FILE_BYTES,
                     f"tracked symlink {relative!r} exceeds its byte bound")
            _require(_sha1_object("blob", target) == object_id,
                     f"tracked symlink {relative!r} differs from its Git blob")
            retained_symlinks.append((
                relative, descriptor, _stable_fields(before), target))
            records.append({
                **entry,
                "kind": "symlink",
                "size": len(target),
                "sha256": hashlib.sha256(target).hexdigest(),
                "target_encoding": "base64",
                "target_base64": base64.b64encode(target).decode("ascii"),
            })
        else:
            subroot = root / relative
            try:
                sub_status = subroot.lstat()
            except OSError as error:
                raise GitCaptureError(
                    f"submodule {relative!r} is not initialized: {error}") \
                    from error
            _require(stat.S_ISDIR(sub_status.st_mode) and
                     not subroot.is_symlink(),
                     f"submodule {relative!r} is not an initialized directory")
            nested = _capture_repository(
                subroot.resolve(strict=True), object_id, git_snapshot, stack,
                depth=depth + 1,
                inherited_descriptors=inherited_descriptors)
            nested_digest = _digest(nested)
            records.append({
                **entry,
                "kind": "submodule",
                "identity_sha256": nested_digest,
            })
            submodules.append({
                "path": relative,
                "object_id": object_id,
                "identity_sha256": nested_digest,
                "identity": nested,
            })

    second = _query_series(
        git_snapshot, metadata, tree.root_descriptor, inherited_descriptors)
    _require(second == first,
             "Git command results changed during source capture")
    _verify_query_series(root, requested_commit, second, preliminary_listing)

    for relative, descriptor, fields, content in retained_regular:
        current = os.fstat(descriptor)
        _require(_stable_fields(current) == fields and
                 _read_bounded_descriptor(
                     descriptor, current.st_size,
                     f"retained tracked file {relative!r}",
                     MAX_TRACKED_SOURCE_FILE_BYTES) == content and
                 _stable_fields(_path_status(tree, relative)) == fields,
                 f"tracked file {relative!r} changed or was replaced")
    for relative, descriptor, fields, target in retained_symlinks:
        current = os.fstat(descriptor)
        parts = tuple(relative.split("/"))
        parent = tree.directory(parts[:-1])
        try:
            current_target = os.fsencode(
                os.readlink(parts[-1], dir_fd=parent))
        except OSError as error:
            raise GitCaptureError(
                f"tracked symlink {relative!r} cannot be re-read: {error}") \
                from error
        _require(_stable_fields(current) == fields and
                 _stable_fields(_path_status(tree, relative)) == fields and
                 current_target == target,
                 f"tracked symlink {relative!r} changed or was replaced")

    metadata.verify()
    git_snapshot.verify()
    worktree_guard.verify()
    gitdir = metadata.gitdir
    common = metadata.common
    _require(gitdir is not None and gitdir.resolved is not None and
             common is not None and common.resolved is not None,
             "retained Git metadata pathname was lost")
    metadata_paths = sorted({
        str(gitdir.resolved), str(common.resolved)})
    executable_status = os.fstat(git_snapshot.executable_descriptor)
    result = {
        "schema": SCHEMA,
        "path": str(root),
        "head": head,
        "tree": tree_id,
        "detached": detached,
        "head_ref": head_ref,
        "superproject_worktree": superproject,
        "tracked_tree_listing_sha256":
            hashlib.sha256(first["tree_listing"]).hexdigest(),
        "tracked_status": "clean",
        "commit_object": _object_identity("commit", first["commit_object"]),
        "tree_objects": tree_objects,
        "git_executable": {
            "source": _portable_file_identity(git_snapshot),
            "sealed": {
                "protocol": GIT_EXECUTABLE_PROTOCOL,
                "size": git_snapshot.executable_identity["size"],
                "mode": stat.S_IMODE(executable_status.st_mode),
                "sha256": git_snapshot.executable_identity["sha256"],
                "seals": git_snapshot.executable_identity["seals"],
                "source_sha256":
                    git_snapshot.executable_identity["source_sha256"],
            },
        },
        "git_metadata": {
            "layout": (
                "linked-worktree" if metadata.entry_file is not None
                else "ordinary"),
            "gitdir": str(gitdir.resolved),
            "commondir": str(common.resolved),
            "guarded_components": metadata_paths,
            "guard_policy": METADATA_GUARD_POLICY,
            "guarded_file_count": len(metadata_files),
            "guarded_files_sha256": _digest(metadata_files),
        },
        "worktree_guard_policy": WORKTREE_GUARD_POLICY,
        "config": _byte_identity(first["config"]),
        "index": {
            "entry_count": len(entries),
            "stage": _byte_identity(first["index_stage"]),
            "flags_v": _byte_identity(first["index_v"]),
            "flags_f": _byte_identity(first["index_f"]),
        },
        "tracked_files": records,
        "tracked_files_sha256": _digest(records),
        "submodules": submodules,
    }
    validate_git_capture(
        result, str(root), requested_commit, require_detached=False)
    return result


def _sealed_file_snapshot_identity(
    snapshot: _RetainedFileSnapshot, protocol: str,
) -> dict[str, Any]:
    """Describe one retained source inode and its immutable memfd copy."""
    sealed_descriptor = snapshot.executable_descriptor
    sealed_status = os.fstat(sealed_descriptor)
    return {
        "source": _portable_file_identity(snapshot),
        "sealed": {
            "protocol": protocol,
            "size": snapshot.executable_identity["size"],
            "mode": stat.S_IMODE(sealed_status.st_mode),
            "sha256": snapshot.executable_identity["sha256"],
            "seals": snapshot.executable_identity["seals"],
            "source_sha256":
                snapshot.executable_identity["source_sha256"],
        },
    }


def _guard_direct_worktree_path(
    root: Path, relative: str, guard: _InotifyMutationGuard,
) -> None:
    """Guard one tracked union path, including a currently missing leaf."""
    parts = tuple(relative.split("/"))
    _require(parts and all(part not in ("", ".", "..") for part in parts),
             f"direct tracked path {relative!r} is unsafe")
    current = root
    for offset, component in enumerate(parts):
        candidate = current / component
        try:
            status = candidate.lstat()
        except FileNotFoundError:
            guard.add_exact_directory_entries(current, (component,))
            return
        except OSError as error:
            raise GitCaptureError(
                f"cannot guard direct tracked path {relative!r}: {error}") \
                from error
        if offset + 1 < len(parts):
            _require(stat.S_ISDIR(status.st_mode) and
                     not stat.S_ISLNK(status.st_mode),
                     f"direct tracked path {relative!r} has an unsafe "
                     "directory component")
            current = candidate
            continue
        if stat.S_ISDIR(status.st_mode):
            guard.add_directory_path(candidate)
        elif stat.S_ISLNK(status.st_mode):
            guard.add_exact_directory_entries(current, (component,))
        else:
            guard.add_file_path(candidate)
            _require(
                _stable_fields(candidate.lstat()) ==
                    _stable_fields(status),
                f"direct tracked path {relative!r} changed while its inode "
                "watch was armed")


def _direct_worktree_is_clean(
    tree_entries: Sequence[Mapping[str, str]],
    index_entries: Sequence[Mapping[str, str]],
    records: Sequence[Mapping[str, Any]], status: bytes,
) -> bool:
    head_layer = [
        (entry["path"], entry["git_mode"], entry["object_id"])
        for entry in tree_entries
    ]
    index_layer = [
        (entry["path"], entry["git_mode"], entry["object_id"])
        for entry in index_entries
    ]
    if status or head_layer != index_layer:
        return False
    for record in records:
        index_mode = record["index_mode"]
        index_object = record["index_object"]
        worktree = record["worktree"]
        kind = worktree["kind"]
        if index_mode is None or index_object is None:
            return False
        if kind in ("regular", "symlink"):
            if (worktree["mode"] != index_mode or
                    worktree["git_object"] != index_object):
                return False
        elif kind == "submodule":
            identity = worktree["identity"]
            if (index_mode != "160000" or identity["head"] != index_object or
                    identity["clean"] is not True):
                return False
        else:
            return False
    return True


def _prepare_direct_guard_repository(
    root: Path,
    git_snapshot: _RetainedFileSnapshot,
    stack: ExitStack,
    *,
    depth: int,
    inherited_descriptors: Sequence[int],
) -> None:
    """Arm a descendant repository's ABA guards before parent series one."""
    _require(depth <= MAX_SUBMODULE_DEPTH,
             "direct Git guard preparation exceeds its recursion bound")
    metadata = stack.enter_context(_RetainedGitMetadata(root))
    _guard_metadata_files(metadata)
    _require_direct_grafts_absent(metadata)
    gitdir = metadata.gitdir
    _require(gitdir is not None and gitdir.resolved is not None,
             "prepared direct Git directory was lost")
    index_snapshot = stack.enter_context(_RetainedFileSnapshot(
        gitdir.resolved / "index", "prepared direct Git index",
        maximum_bytes=MAX_GIT_LISTING_BYTES))
    index_snapshot.executable_descriptor
    worktree_guard = _InotifyMutationGuard(
        f"prepared direct tracked worktree {root}")
    stack.callback(worktree_guard.close)
    worktree_guard.add_topology_tree(
        root, excluded_root_entries=(".git",),
        guarded_leaf_names=(".gitignore",))
    tree = _OpenDirectoryTree(
        root, stack, f"prepared direct tracked worktree {root}")
    preliminary = _query_direct_series(
        git_snapshot, metadata, index_snapshot, tree.root_descriptor,
        inherited_descriptors)
    (_, _, _, _, _, tree_entries, index_entries,
     _) = _verify_direct_query_series(root, preliminary)
    head_by_path = {entry["path"]: entry for entry in tree_entries}
    index_by_path = {entry["path"]: entry for entry in index_entries}
    union_paths = sorted(set(head_by_path) | set(index_by_path))
    _require(union_paths and
             len(union_paths) <= 2 * MAX_TRACKED_SOURCE_FILES,
             "prepared direct tracked source union is empty or oversized")
    for relative in union_paths:
        _guard_direct_worktree_path(root, relative, worktree_guard)
    worktree_guard.verify()
    armed = _query_direct_series(
        git_snapshot, metadata, index_snapshot, tree.root_descriptor,
        inherited_descriptors)
    _require(armed == preliminary,
             "direct descendant changed while arming ABA guards")
    for relative in union_paths:
        head_entry = head_by_path.get(relative)
        index_entry = index_by_path.get(relative)
        modes = {
            entry["git_mode"] for entry in (head_entry, index_entry)
            if entry is not None
        }
        if "160000" not in modes:
            continue
        path = root / relative
        try:
            status = path.lstat()
        except FileNotFoundError:
            continue
        except OSError as error:
            raise GitCaptureError(
                f"cannot prepare direct submodule {relative!r}: {error}") \
                from error
        if stat.S_ISDIR(status.st_mode) and not stat.S_ISLNK(status.st_mode):
            _prepare_direct_guard_repository(
                path.resolve(strict=True), git_snapshot, stack,
                depth=depth + 1,
                inherited_descriptors=inherited_descriptors)
    metadata.verify()
    index_snapshot.verify()
    git_snapshot.verify()
    worktree_guard.verify()


def _capture_direct_repository(
    root: Path,
    git_snapshot: _RetainedFileSnapshot,
    stack: ExitStack,
    *,
    depth: int,
    inline_paths: frozenset[str],
    parent_worktree: str | None,
    inherited_descriptors: Sequence[int],
) -> dict[str, Any]:
    """Capture HEAD, sealed index, and worktree as three guarded layers."""
    _require(depth <= MAX_SUBMODULE_DEPTH,
             "direct Git submodule nesting exceeds its recursion bound")
    metadata = stack.enter_context(_RetainedGitMetadata(root))
    metadata_files = _guard_metadata_files(metadata)
    _require_direct_grafts_absent(metadata)
    gitdir = metadata.gitdir
    _require(gitdir is not None and gitdir.resolved is not None,
             "retained direct Git directory was lost")
    index_snapshot = stack.enter_context(_RetainedFileSnapshot(
        gitdir.resolved / "index", "direct Git index",
        maximum_bytes=MAX_GIT_LISTING_BYTES))
    # Materialize the sealed index before any Git command.  Every index-aware
    # query below receives this exact memfd through GIT_INDEX_FILE.
    index_snapshot.executable_descriptor

    worktree_guard = _InotifyMutationGuard(
        f"direct tracked worktree {root}")
    stack.callback(worktree_guard.close)
    worktree_guard.add_topology_tree(
        root, excluded_root_entries=(".git",),
        guarded_leaf_names=(".gitignore",))
    tree = _OpenDirectoryTree(root, stack, f"direct tracked worktree {root}")

    preliminary = _query_direct_series(
        git_snapshot, metadata, index_snapshot, tree.root_descriptor,
        inherited_descriptors)
    (_, _, _, _, _, preliminary_tree, preliminary_index,
     _) = _verify_direct_query_series(root, preliminary)
    union_paths = sorted(
        {entry["path"] for entry in preliminary_tree} |
        {entry["path"] for entry in preliminary_index})
    _require(union_paths and
             len(union_paths) <= 2 * MAX_TRACKED_SOURCE_FILES,
             "direct tracked source union is empty or oversized")
    for relative in union_paths:
        _guard_direct_worktree_path(root, relative, worktree_guard)
    worktree_guard.verify()

    preliminary_head = {
        entry["path"]: entry for entry in preliminary_tree}
    preliminary_staged = {
        entry["path"]: entry for entry in preliminary_index}
    for relative in union_paths:
        modes = {
            entry["git_mode"] for entry in (
                preliminary_head.get(relative),
                preliminary_staged.get(relative))
            if entry is not None
        }
        if "160000" not in modes:
            continue
        path = root / relative
        try:
            status = path.lstat()
        except FileNotFoundError:
            continue
        except OSError as error:
            raise GitCaptureError(
                f"cannot prepare direct submodule {relative!r}: {error}") \
                from error
        if stat.S_ISDIR(status.st_mode) and not stat.S_ISLNK(status.st_mode):
            _prepare_direct_guard_repository(
                path.resolve(strict=True), git_snapshot, stack,
                depth=depth + 1,
                inherited_descriptors=inherited_descriptors)

    first = _query_direct_series(
        git_snapshot, metadata, index_snapshot, tree.root_descriptor,
        inherited_descriptors)
    _require(first == preliminary,
             "direct Git command results changed while arming worktree "
             "guards")
    (head, tree_id, detached, head_ref, superproject,
     tree_entries, index_entries, untracked_output) = \
        _verify_direct_query_series(root, first)
    tree_objects = _capture_tree_objects(
        git_snapshot, metadata, tree.root_descriptor, tree_id,
        first["tree_object"], inherited_descriptors)
    tree_object_bytes = {
        record["object_id"]: base64.b64decode(
            record["base64"], validate=True)
        for record in tree_objects
    }
    _require(_flatten_tree_objects(tree_id, tree_object_bytes) == tree_entries,
             "direct recursive Git tree objects differ from their listing")

    head_by_path = {entry["path"]: entry for entry in tree_entries}
    index_by_path = {entry["path"]: entry for entry in index_entries}
    retained_regular: list[
        tuple[str, int, tuple[int, ...], bytes]
    ] = []
    retained_symlinks: list[
        tuple[str, int, tuple[int, ...], bytes]
    ] = []
    retained_missing: list[str] = []
    retained_submodules: list[tuple[str, tuple[int, ...]]] = []
    records: list[dict[str, Any]] = []
    total = 0
    for relative in union_paths:
        head_entry = head_by_path.get(relative)
        index_entry = index_by_path.get(relative)
        record: dict[str, Any] = {
            "path": relative,
            "head_mode": (
                head_entry["git_mode"] if head_entry is not None else None),
            "head_object": (
                head_entry["object_id"] if head_entry is not None else None),
            "index_mode": (
                index_entry["git_mode"]
                if index_entry is not None else None),
            "index_object": (
                index_entry["object_id"]
                if index_entry is not None else None),
        }
        path = root / relative
        try:
            path_status = path.lstat()
        except FileNotFoundError:
            _require(record["index_mode"] != "160000",
                     f"direct submodule {relative!r} is not initialized")
            record["worktree"] = {"kind": "missing"}
            retained_missing.append(relative)
            records.append(record)
            continue
        except OSError as error:
            raise GitCaptureError(
                f"cannot inspect direct tracked path {relative!r}: {error}") \
                from error

        if stat.S_ISREG(path_status.st_mode):
            descriptor, before = tree.open_regular(relative)
            _require(before.st_size <= MAX_TRACKED_SOURCE_FILE_BYTES,
                     f"direct tracked file {relative!r} exceeds its bound")
            total += before.st_size
            _require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
                     "direct tracked source bytes exceed their total bound")
            content = _read_bounded_descriptor(
                descriptor, before.st_size,
                f"direct tracked file {relative!r}",
                MAX_TRACKED_SOURCE_FILE_BYTES)
            after = os.fstat(descriptor)
            _require(_stable_fields(before) == _stable_fields(after),
                     f"direct tracked file {relative!r} changed while read")
            retained_regular.append((
                relative, descriptor, _stable_fields(after), content))
            mode = "100755" if stat.S_IMODE(after.st_mode) & 0o111 \
                else "100644"
            worktree_record = {
                "kind": "regular",
                "mode": mode,
                "git_object": _sha1_object("blob", content),
            }
            if relative in inline_paths:
                worktree_record["payload"] = _byte_payload(content)
            record["worktree"] = worktree_record
        elif stat.S_ISLNK(path_status.st_mode):
            descriptor, before, target = _open_symlink(
                tree, relative, stack)
            total += len(target)
            _require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES and
                     len(target) <= MAX_TRACKED_SOURCE_FILE_BYTES,
                     f"direct tracked symlink {relative!r} exceeds its bound")
            retained_symlinks.append((
                relative, descriptor, _stable_fields(before), target))
            record["worktree"] = {
                "kind": "symlink",
                "mode": "120000",
                "size": len(target),
                "sha256": hashlib.sha256(target).hexdigest(),
                "git_object": _sha1_object("blob", target),
                "target_encoding": "base64",
                "target_base64": base64.b64encode(target).decode("ascii"),
            }
        elif stat.S_ISDIR(path_status.st_mode):
            expected_modes = {
                mode for mode in (
                    record["head_mode"], record["index_mode"])
                if mode is not None
            }
            _require(expected_modes <= {"160000"},
                     f"direct tracked blob {relative!r} became a directory")
            nested_root = path.resolve(strict=True)
            _require(nested_root.is_dir(),
                     f"direct submodule {relative!r} is not a directory")
            nested = _capture_direct_repository(
                nested_root, git_snapshot, stack, depth=depth + 1,
                inline_paths=frozenset(),
                parent_worktree=str(root),
                inherited_descriptors=inherited_descriptors)
            nested_digest = _digest(nested)
            retained_submodules.append(
                (relative, _stable_fields(path_status)))
            record["worktree"] = {
                "kind": "submodule",
                "identity_sha256": nested_digest,
                "identity": nested,
            }
        else:
            raise GitCaptureError(
                f"direct tracked path {relative!r} has an unsupported type")
        records.append(record)

    retained_inline_paths = {
        record["path"] for record in records
        if isinstance(record.get("worktree"), dict) and
        "payload" in record["worktree"]
    }
    _require(retained_inline_paths == inline_paths,
             "direct inline source paths are missing or are not regular")

    second = _query_direct_series(
        git_snapshot, metadata, index_snapshot, tree.root_descriptor,
        inherited_descriptors)
    _require(second == first,
             "direct Git command results changed during source capture")
    _verify_direct_query_series(root, second)

    for relative, descriptor, fields, content in retained_regular:
        current = os.fstat(descriptor)
        _require(_stable_fields(current) == fields and
                 _read_bounded_descriptor(
                     descriptor, current.st_size,
                     f"retained direct file {relative!r}",
                     MAX_TRACKED_SOURCE_FILE_BYTES) == content and
                 _stable_fields(_path_status(tree, relative)) == fields,
                 f"direct tracked file {relative!r} changed or was replaced")
    for relative, descriptor, fields, target in retained_symlinks:
        current = os.fstat(descriptor)
        parts = tuple(relative.split("/"))
        parent = tree.directory(parts[:-1])
        try:
            current_target = os.fsencode(
                os.readlink(parts[-1], dir_fd=parent))
        except OSError as error:
            raise GitCaptureError(
                f"direct tracked symlink {relative!r} cannot be re-read: "
                f"{error}") from error
        _require(_stable_fields(current) == fields and
                 _stable_fields(_path_status(tree, relative)) == fields and
                 current_target == target,
                 f"direct tracked symlink {relative!r} changed or was replaced")
    for relative in retained_missing:
        _require(not (root / relative).exists() and
                 not (root / relative).is_symlink(),
                 f"missing direct tracked path {relative!r} appeared")
    for relative, fields in retained_submodules:
        _require(_stable_fields((root / relative).lstat()) == fields,
                 f"direct submodule path {relative!r} changed or was replaced")

    metadata.verify()
    index_snapshot.verify()
    git_snapshot.verify()
    worktree_guard.verify()
    common = metadata.common
    _require(gitdir.resolved is not None and common is not None and
             common.resolved is not None,
             "retained direct Git metadata pathname was lost")
    metadata_paths = sorted({str(gitdir.resolved), str(common.resolved)})
    status_output = _encode_direct_status(
        _derive_direct_status_records(
            {record["path"]: record for record in records},
            _parse_direct_path_list(
                untracked_output, "direct Git untracked path inventory"),
            core_filemode=_direct_boolean_config(
                first["config"], "core.filemode", default=True),
        )
    )
    clean = _direct_worktree_is_clean(
        tree_entries, index_entries, records, status_output)
    result = {
        "schema": DIRECT_SCHEMA,
        "path": str(root),
        "head": head,
        "tree": tree_id,
        "detached": detached,
        "head_ref": head_ref,
        # Explicit --git-dir/--work-tree queries cannot rediscover their
        # superproject relationship.  Recursive capture supplies the exact
        # containing proof path instead, which the offline validator binds.
        "superproject_worktree": parent_worktree,
        "refs": _byte_payload(first["refs"]),
        "commit_object": _object_identity("commit", first["commit_object"]),
        "tree_objects": tree_objects,
        "head_tree_listing": _byte_payload(first["tree_listing"]),
        "git_executable": _sealed_file_snapshot_identity(
            git_snapshot, GIT_EXECUTABLE_PROTOCOL),
        "git_metadata": {
            "layout": (
                "linked-worktree" if metadata.entry_file is not None
                else "ordinary"),
            "gitdir": str(gitdir.resolved),
            "commondir": str(common.resolved),
            "git_entry": (
                _byte_payload(metadata.entry_file.content)
                if metadata.entry_file is not None else None),
            "commondir_file": (
                _byte_payload(metadata.common_file.content)
                if metadata.common_file is not None else None),
            "guarded_components": metadata_paths,
            "guard_policy": METADATA_GUARD_POLICY,
            "guarded_files": metadata_files,
            "guarded_file_count": len(metadata_files),
            "guarded_files_sha256": _digest(metadata_files),
            "legacy_grafts_absent": True,
        },
        "worktree_guard_policy": DIRECT_WORKTREE_GUARD_POLICY,
        "config": _byte_payload(first["config"]),
        "index": {
            **_sealed_file_snapshot_identity(
                index_snapshot, DIRECT_INDEX_PROTOCOL),
            "entry_count": len(index_entries),
            "raw": _byte_payload(index_snapshot.content),
            "stage": _byte_payload(first["index_stage"]),
            "flags_v": _byte_payload(first["index_v"]),
            "flags_f": _byte_payload(first["index_f"]),
        },
        "status": _byte_payload(status_output),
        "untracked_paths": _byte_payload(first["untracked_paths"]),
        "inline_paths": sorted(inline_paths),
        "paths": records,
        "paths_sha256": _digest(records),
        "clean": clean,
    }
    result["identity_sha256"] = _digest(result)
    validate_direct_git_capture(
        result, str(root), require_clean=False,
        _expected_superproject=parent_worktree)
    return result


def capture_git_identities(
    requests: Sequence[tuple[Path | str, str, bool]],
    *, inherited_descriptors: Sequence[int] = (),
) -> list[dict[str, Any]]:
    """Capture repositories while one Git snapshot and every guard stay live."""
    try:
        _require(isinstance(requests, Sequence) and
                 not isinstance(requests, (str, bytes)) and requests,
                 "Git source capture request set is empty or invalid")
        _require(isinstance(inherited_descriptors, Sequence) and
                 not isinstance(inherited_descriptors, (str, bytes)) and
                 all(type(descriptor) is int and descriptor >= 0
                     for descriptor in inherited_descriptors),
                 "Git source inherited descriptor set is invalid")
        for descriptor in set(inherited_descriptors):
            try:
                os.fstat(descriptor)
            except OSError as error:
                raise GitCaptureError(
                    f"Git source inherited descriptor {descriptor} is "
                    f"invalid: {error}") from error
        normalized: list[tuple[Path, str, bool]] = []
        for request in requests:
            _require(isinstance(request, Sequence) and
                     not isinstance(request, (str, bytes)) and
                     len(request) == 3,
                     "Git source capture request is malformed")
            source_root, requested_commit, require_detached = request
            _require(isinstance(source_root, (str, os.PathLike)),
                     "Git source root request is invalid")
            _require(isinstance(requested_commit, str) and
                     HEX40.fullmatch(requested_commit) is not None,
                     "requested source commit must be one lowercase SHA-1 ID")
            _require(type(require_detached) is bool,
                     "Git source detached requirement is invalid")
            root = Path(source_root).resolve(strict=True)
            _require_safe_unicode(str(root), "Git source root")
            _require(root.is_dir(), f"source root is not a directory: {root}")
            try:
                git_entry = (root / ".git").lstat()
            except FileNotFoundError:
                raise GitCaptureError(
                    f"source root is not the Git top level: {root}")
            _require(stat.S_ISDIR(git_entry.st_mode) or
                     stat.S_ISREG(git_entry.st_mode),
                     f"source root is not the Git top level: {root}")
            normalized.append((root, requested_commit, require_detached))
        with ExitStack() as stack:
            git_snapshot = stack.enter_context(
                _RetainedFileSnapshot(
                    Path(GIT_EXECUTABLE_SOURCE_PATH),
                    "benchmark source Git executable"))
            # Materialize and attest the immutable executable exactly once.
            git_snapshot.executable_descriptor
            results = []
            for root, requested_commit, require_detached in normalized:
                result = _capture_repository(
                    root, requested_commit, git_snapshot, stack, depth=0,
                    inherited_descriptors=inherited_descriptors)
                _require(not require_detached or result["detached"] is True,
                         f"exact-main source {root} is not detached")
                results.append(result)
            git_snapshot.verify()
            return results
    except GitCaptureError:
        raise
    except BuildProvenanceError as error:
        raise GitCaptureError(str(error)) from error
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        raise GitCaptureError(f"cannot capture Git source identity: {error}") \
            from error


def capture_git_identity(
    source_root: Path | str, requested_commit: str, *,
    require_detached: bool = False,
    inherited_descriptors: Sequence[int] = (),
) -> dict[str, Any]:
    """Capture one source tree under a single retained/sealed Git identity."""
    return capture_git_identities(
        ((source_root, requested_commit, require_detached),),
        inherited_descriptors=inherited_descriptors)[0]


def capture_direct_git_identity(
    source_root: Path | str, *, require_clean: bool,
    expected_head: str | None = None,
    inline_paths: Sequence[str] = (),
    inherited_descriptors: Sequence[int] = (),
) -> dict[str, Any]:
    """Capture one direct-runner source with clean or dirty state retained."""
    try:
        _require(isinstance(source_root, (str, os.PathLike)),
                 "direct Git source root request is invalid")
        _require(type(require_clean) is bool,
                 "direct Git clean-source requirement is invalid")
        _require(expected_head is None or
                 (isinstance(expected_head, str) and
                  HEX40.fullmatch(expected_head) is not None),
                 "expected direct Git HEAD must be one lowercase SHA-1 ID")
        _require(isinstance(inline_paths, Sequence) and
                 not isinstance(inline_paths, (str, bytes)),
                 "direct Git inline-path request is invalid")
        normalized_inline_paths = frozenset(
            _safe_path(
                path.encode("utf-8"), "direct Git inline source path")
            for path in inline_paths
            if isinstance(path, str)
        )
        _require(len(normalized_inline_paths) == len(inline_paths),
                 "direct Git inline source paths are duplicated or invalid")
        _require(isinstance(inherited_descriptors, Sequence) and
                 not isinstance(inherited_descriptors, (str, bytes)) and
                 all(type(descriptor) is int and descriptor >= 0
                     for descriptor in inherited_descriptors),
                 "direct Git inherited descriptor set is invalid")
        for descriptor in set(inherited_descriptors):
            try:
                os.fstat(descriptor)
            except OSError as error:
                raise GitCaptureError(
                    f"direct Git inherited descriptor {descriptor} is "
                    f"invalid: {error}") from error
        root = Path(source_root).resolve(strict=True)
        _require_safe_unicode(str(root), "direct Git source root")
        _require(root.is_dir(),
                 f"direct Git source root is not a directory: {root}")
        try:
            git_entry = (root / ".git").lstat()
        except FileNotFoundError:
            raise GitCaptureError(
                f"direct source root is not the Git top level: {root}")
        _require(stat.S_ISDIR(git_entry.st_mode) or
                 stat.S_ISREG(git_entry.st_mode),
                 f"direct source root is not the Git top level: {root}")
        with ExitStack() as stack:
            git_snapshot = stack.enter_context(_RetainedFileSnapshot(
                Path(GIT_EXECUTABLE_SOURCE_PATH),
                "direct source Git executable"))
            git_snapshot.executable_descriptor
            result = _capture_direct_repository(
                root, git_snapshot, stack, depth=0,
                inline_paths=normalized_inline_paths,
                parent_worktree=None,
                inherited_descriptors=inherited_descriptors)
            validate_direct_git_capture(
                result, str(root), require_clean=require_clean,
                expected_head=expected_head)
            git_snapshot.verify()
            return result
    except GitCaptureError:
        raise
    except BuildProvenanceError as error:
        raise GitCaptureError(str(error)) from error
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        raise GitCaptureError(
            f"cannot capture direct Git source identity: {error}") from error


def _parse_direct_path_list(raw: bytes, label: str) -> list[str]:
    """Parse one canonical NUL-terminated Git pathname inventory."""
    if raw:
        _require(raw.endswith(b"\0"), f"{label} is not NUL terminated")
    encoded = [] if not raw else raw[:-1].split(b"\0")
    _require(all(encoded), f"{label} contains an empty pathname")
    paths = [_safe_path(path, label) for path in encoded]
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             f"{label} is duplicated or non-canonical")
    return paths


def _parse_direct_status(raw: bytes) -> list[dict[str, str]]:
    """Parse the direct dialect's no-renames porcelain-v1 status."""
    if not raw:
        return []
    _require(raw.endswith(b"\0"),
             "direct Git status is not NUL terminated")
    records = raw[:-1].split(b"\0")
    _require(records and all(records),
             "direct Git status contains an empty record")
    parsed: list[dict[str, str]] = []
    tracked_codes = {
        b" M", b" T", b" D",
        b"M ", b"MM", b"MT", b"MD",
        b"T ", b"TM", b"TT", b"TD",
        b"A ", b"AM", b"AT", b"AD",
        b"D ",
    }
    for record in records:
        _require(len(record) >= 4 and record[2:3] == b" ",
                 "direct Git status record has an invalid prefix")
        status_code = record[:2]
        _require(status_code == b"??" or status_code in tracked_codes,
                 "direct Git status code is invalid or uses rename/copy "
                 "detection")
        path = _safe_path(
            record[3:], "direct Git status pathname")
        parsed.append({
            "code": status_code.decode("ascii"),
            "path": path,
        })
    return parsed


def _direct_status_lines(raw: bytes) -> list[str]:
    """Project canonical direct status records to compatibility strings."""
    return [
        record["code"] + " " + record["path"]
        for record in _parse_direct_status(raw)
    ]


def _direct_mode_kind(mode: str | None) -> str | None:
    if mode in {"100644", "100755"}:
        return "regular"
    if mode == "120000":
        return "symlink"
    if mode == "160000":
        return "submodule"
    return None


def _direct_record_status(record: Mapping[str, Any]) -> tuple[str, str]:
    """Derive claimed staged/worktree changes from retained exact layers."""
    head_mode = record["head_mode"]
    head_object = record["head_object"]
    index_mode = record["index_mode"]
    index_object = record["index_object"]
    if head_mode is None:
        staged = "A"
    elif index_mode is None:
        staged = "D"
    elif head_mode == index_mode and head_object == index_object:
        staged = " "
    elif _direct_mode_kind(head_mode) != _direct_mode_kind(index_mode):
        staged = "T"
    else:
        staged = "M"

    if index_mode is None:
        return staged, " "
    worktree = record["worktree"]
    kind = worktree["kind"]
    if kind == "missing":
        return staged, "D"
    if kind == "submodule":
        nested = worktree["identity"]
        worktree_status = (
            " " if index_mode == "160000" and
            nested["head"] == index_object and nested["clean"] is True
            else "M" if index_mode == "160000" else "T"
        )
        return staged, worktree_status
    if _direct_mode_kind(index_mode) != kind:
        return staged, "T"
    if (worktree["mode"] == index_mode and
            worktree["git_object"] == index_object):
        return staged, " "
    return staged, "M"


def _derive_direct_status_records(
    records_by_path: Mapping[str, Mapping[str, Any]],
    untracked_paths: Sequence[str],
    *, core_filemode: bool,
) -> list[dict[str, str]]:
    """Derive canonical no-renames porcelain from exact retained layers."""
    result: list[dict[str, str]] = []
    for path in sorted(records_by_path):
        record = records_by_path[path]
        _staged, code = _direct_record_status(record)
        worktree = record["worktree"]
        mode_only_hidden = (
            not core_filemode and code == "M" and
            record["index_mode"] in {"100644", "100755"} and
            worktree["kind"] == "regular" and
            worktree["git_object"] == record["index_object"] and
            worktree["mode"] != record["index_mode"]
        )
        if mode_only_hidden:
            code = " "
        staged, _unused = _direct_record_status(record)
        combined = staged + code
        if combined != "  ":
            result.append({"code": combined, "path": path})

    for path in untracked_paths:
        record = records_by_path.get(path)
        if record is not None:
            _require(record["index_mode"] is None and
                     record["worktree"]["kind"] != "missing",
                     "direct Git untracked status contradicts its tracked "
                     "source layers")
        result.append({"code": "??", "path": path})
    return result


def _encode_direct_status(
    records: Sequence[Mapping[str, str]],
) -> bytes:
    return b"".join(
        record["code"].encode("ascii") + b" " +
        record["path"].encode("utf-8") + b"\0"
        for record in records
    )


def direct_source_projection(
    value: Mapping[str, Any], *, _expected_superproject: str | None = None,
) -> dict[str, Any]:
    """Derive the direct runner's compatibility view from one rich proof."""
    _require(isinstance(value, Mapping),
             "direct Git projection input is invalid")
    retained = dict(value)
    validate_direct_git_capture(
        retained, retained.get("path", ""), require_clean=False,
        _expected_superproject=_expected_superproject)
    tree_objects = {
        record["object_id"]: base64.b64decode(
            record["base64"], validate=True)
        for record in retained["tree_objects"]
    }
    tree_entries = _flatten_tree_objects(retained["tree"], tree_objects)
    index_entries = _parse_index(_validate_byte_payload(
        retained["index"]["stage"], "projected direct Git index"))
    head_material = {
        entry["path"]: [entry["git_mode"], entry["object_id"]]
        for entry in tree_entries
    }
    index_material = {
        entry["path"]: [entry["git_mode"], entry["object_id"]]
        for entry in index_entries
    }
    records = {record["path"]: record for record in retained["paths"]}
    files: dict[str, Any] = {}
    for entry in index_entries:
        relative = entry["path"]
        record = records[relative]
        worktree = record["worktree"]
        projected: dict[str, Any] = {
            "index_mode": entry["git_mode"],
            "index_object": entry["object_id"],
        }
        if worktree["kind"] in {"regular", "symlink"}:
            projected.update({
                "worktree_git_object": worktree["git_object"],
                "worktree_mode": worktree["mode"],
            })
            if worktree["kind"] == "symlink":
                projected["worktree_sha256"] = worktree["sha256"]
            elif "payload" in worktree:
                projected["worktree_sha256"] = \
                    worktree["payload"]["sha256"]
        elif worktree["kind"] == "missing":
            projected["worktree_missing"] = True
        else:
            nested = worktree["identity"]
            nested_projection = direct_source_projection(
                nested, _expected_superproject=retained["path"])
            nested_status = _direct_status_lines(_validate_byte_payload(
                nested["status"], "projected direct submodule status"))
            projected["submodule"] = {
                "digest": nested_projection["digest"],
                "files": nested_projection["files"],
                "head": nested["head"],
                "repository_controls":
                    nested_projection["repository_controls"],
                "status": nested_status,
                "worktree_clean": nested["clean"],
            }
        files[relative] = projected
    controls = {
        "head_index_match": head_material == index_material,
        "head_tree_sha256": _digest(head_material),
        "index_entry_count": len(index_entries),
        "index_flags": "all ls-files -v/-f tags are canonical H",
        "index_sha256": _digest(index_material),
        "legacy_grafts_absent": True,
        "replace_refs": [],
    }
    status_lines = _direct_status_lines(_validate_byte_payload(
        retained["status"], "projected direct Git status"))
    head_ref = retained["head_ref"]
    branch = (
        head_ref[len("refs/heads/"):]
        if isinstance(head_ref, str) and head_ref.startswith("refs/heads/")
        else None
    )
    material = {
        "files": files,
        "repository_controls": controls,
    }
    return {
        "digest": _digest(material),
        "files": files,
        "git": {
            "branch": branch,
            "head": retained["head"],
            "status": status_lines,
            "worktree_clean": retained["clean"],
            "tree": retained["tree"],
        },
        "git_tool": {
            "environment": dict(_build_provenance.GIT_ENVIRONMENT),
            "path": retained["git_executable"]["source"]["path"],
            "sha256": retained["git_executable"]["source"]["sha256"],
        },
        "repository_controls": controls,
    }


def legacy_projection(
    identity: Mapping[str, Any], *, include_commit_object: bool,
) -> dict[str, Any]:
    """Project a current capture into the immutable historical runner shape."""
    result = {
        key: identity[key]
        for key in (
            "path", "head", "tree", "detached",
            "tracked_tree_listing_sha256", "tracked_status")
    }
    if include_commit_object:
        result["commit_object"] = identity["commit_object"]
    return result


def _validate_byte_identity(value: object, label: str) -> dict[str, Any]:
    _require(isinstance(value, dict) and set(value) == {"size", "sha256"} and
             type(value.get("size")) is int and value["size"] >= 0 and
             isinstance(value.get("sha256"), str) and
             HEX256.fullmatch(value["sha256"]) is not None,
             f"{label} byte identity is invalid")
    return value


def _validate_byte_payload(value: object, label: str) -> bytes:
    _require(isinstance(value, dict) and set(value) == {
                 "encoding", "size", "sha256", "base64"} and
             value.get("encoding") == "base64" and
             type(value.get("size")) is int and value["size"] >= 0 and
             isinstance(value.get("sha256"), str) and
             HEX256.fullmatch(value["sha256"]) is not None and
             isinstance(value.get("base64"), str),
             f"{label} byte payload is invalid")
    try:
        content = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise GitCaptureError(f"{label} is not canonical base64") from error
    _require(base64.b64encode(content).decode("ascii") == value["base64"] and
             len(content) == value["size"] and
             hashlib.sha256(content).hexdigest() == value["sha256"],
             f"{label} retained bytes differ")
    return content


def _validate_git_path_reference(
    value: object, base: Path, label: str, *, prefix: bytes = b"",
) -> str:
    """Derive one retained gitfile/commondir target without filesystem I/O."""
    content = _validate_byte_payload(value, label)
    _require(content.endswith(b"\n") and content.count(b"\n") == 1,
             f"{label} is not one canonical LF-terminated record")
    record = content[:-1]
    _require(record.startswith(prefix) and len(record) > len(prefix),
             f"{label} has an invalid prefix")
    try:
        path_text = record[len(prefix):].decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise GitCaptureError(f"{label} is not strict UTF-8") from error
    _require_safe_unicode(path_text, label)
    _require("\0" not in path_text and
             len(os.fsencode(path_text)) <= MAX_TRACKED_SOURCE_PATH_BYTES,
             f"{label} path is invalid or oversized")
    path = Path(path_text)
    if not path.is_absolute():
        path = base / path
    normalized = str(Path(os.path.abspath(os.fspath(path))))
    _require(Path(normalized).is_absolute(),
             f"{label} target is not absolute")
    _require_safe_unicode(normalized, f"{label} normalized target")
    return normalized


def _validate_sealed_file_proof(
    value: object, protocol: str, label: str, *,
    require_source_executable: bool, expected_source_path: str | None = None,
) -> dict[str, Any]:
    _require(isinstance(value, dict) and set(value) == {"source", "sealed"},
             f"{label} snapshot shape differs")
    source = value["source"]
    sealed = value["sealed"]
    _require(isinstance(source, dict) and set(source) == {
                 "path", "size", "mode", "sha256"} and
             isinstance(source.get("path"), str) and
             Path(source["path"]).is_absolute() and
             (expected_source_path is None or
              source["path"] == expected_source_path) and
             type(source.get("size")) is int and source["size"] >= 0 and
             type(source.get("mode")) is int and
             0 <= source["mode"] <= 0o7777 and
             (not require_source_executable or source["mode"] & 0o111) and
             isinstance(source.get("sha256"), str) and
             HEX256.fullmatch(source["sha256"]) is not None,
             f"{label} source identity is invalid")
    _require_safe_unicode(source["path"], f"{label} source path")
    _require(isinstance(sealed, dict) and set(sealed) == {
                 "protocol", "size", "mode", "sha256", "seals",
                 "source_sha256"} and
             sealed.get("protocol") == protocol and
             sealed.get("size") == source["size"] and
             sealed.get("sha256") == source["sha256"] and
             sealed.get("source_sha256") == source["sha256"] and
             type(sealed.get("mode")) is int and
             sealed["mode"] == source["mode"] and
             (not require_source_executable or sealed["mode"] & 0o111) and
             type(sealed.get("seals")) is int and
             sealed["seals"] & REQUIRED_SEALS == REQUIRED_SEALS,
             f"{label} sealed identity is invalid")
    return value


def _validate_direct_config(raw: bytes) -> None:
    _validate_local_config(raw)
    unsupported = {"core.sparsecheckout", "index.sparse"}
    external_path_settings = {"core.excludesfile"}
    truthy = {b"1", b"true", b"yes", b"on"}
    for record in raw.split(b"\0"):
        if not record:
            continue
        key, separator, value = record.partition(b"\n")
        _require(bool(separator),
                 "direct Git configuration record is malformed")
        try:
            name = key.decode("ascii", errors="strict").lower()
        except UnicodeDecodeError as error:
            raise GitCaptureError(
                "direct Git configuration key is not ASCII") from error
        _require(not (name in unsupported and value.strip().lower() in truthy),
                 "direct Git sparse index/worktree mode is unsupported")
        _require(name != "extensions.objectformat",
                 "direct Git non-SHA-1 object format is unsupported")
        _require(name not in external_path_settings,
                 "direct Git external excludes file is unsupported")


def _direct_boolean_config(
    raw: bytes, name: str, *, default: bool,
) -> bool:
    """Read one behavior-bearing local boolean from validated config bytes."""
    values: list[bytes] = []
    for record in raw[:-1].split(b"\0") if raw else ():
        key, _separator, value = record.partition(b"\n")
        if key.decode("ascii", errors="strict").lower() == name:
            values.append(value.strip().lower())
    _require(len(values) <= 1,
             f"direct Git configuration repeats {name}")
    if not values:
        return default
    truthy = {b"1", b"true", b"yes", b"on"}
    falsy = {b"0", b"false", b"no", b"off", b""}
    _require(values[0] in truthy | falsy,
             f"direct Git configuration has an invalid {name} boolean")
    return values[0] in truthy


def _validate_object_identity(
    value: object, kind: str, expected_id: str, label: str,
) -> bytes:
    _require(isinstance(value, dict) and set(value) == {
                 "encoding", "size", "sha256", "object_id", "base64"} and
             value.get("encoding") == "base64" and
             value.get("object_id") == expected_id and
             isinstance(value.get("base64"), str),
             f"{label} object identity shape differs")
    try:
        content = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise GitCaptureError(f"{label} object is not canonical base64") \
            from error
    _require(base64.b64encode(content).decode("ascii") == value["base64"] and
             type(value.get("size")) is int and value["size"] == len(content) and
             value.get("sha256") == hashlib.sha256(content).hexdigest() and
             _sha1_object(kind, content) == expected_id,
             f"{label} object bytes differ")
    return content


def validate_direct_git_capture(
    value: object, expected_path: str, *, require_clean: bool,
    expected_head: str | None = None, _depth: int = 0,
    _expected_git_executable: Mapping[str, Any] | None = None,
    _expected_superproject: str | None = None,
) -> dict[str, Any]:
    """Validate a portable direct-runner capture without filesystem I/O."""
    _require(_depth <= MAX_SUBMODULE_DEPTH,
             "retained direct Git submodule identity exceeds its bound")
    _require(type(require_clean) is bool,
             "retained direct Git clean-source requirement is invalid")
    expected_keys = {
        "schema", "identity_sha256", "path", "head", "tree",
        "detached", "head_ref", "superproject_worktree", "refs",
        "commit_object", "tree_objects", "head_tree_listing",
        "git_executable", "git_metadata", "worktree_guard_policy",
        "config", "index", "status", "untracked_paths", "inline_paths",
        "paths", "paths_sha256", "clean",
    }
    _require(isinstance(value, dict) and set(value) == expected_keys and
             value.get("schema") == DIRECT_SCHEMA,
             "direct Git source capture shape or schema differs")
    digest_material = {
        key: item for key, item in value.items()
        if key != "identity_sha256"
    }
    _require(isinstance(value.get("identity_sha256"), str) and
             HEX256.fullmatch(value["identity_sha256"]) is not None and
             value["identity_sha256"] == _digest(digest_material),
             "direct Git source capture digest differs")
    _require(isinstance(expected_path, str) and
             Path(expected_path).is_absolute() and
             value.get("path") == expected_path and
             isinstance(value.get("head"), str) and
             HEX40.fullmatch(value["head"]) is not None and
             (expected_head is None or value["head"] == expected_head) and
             isinstance(value.get("tree"), str) and
             HEX40.fullmatch(value["tree"]) is not None and
             type(value.get("detached")) is bool and
             type(value.get("clean")) is bool,
             "direct Git source capture identity is invalid")
    _require_safe_unicode(expected_path, "retained direct Git source root")
    head_ref = value.get("head_ref")
    _require((value["detached"] and head_ref is None) or
             (not value["detached"] and isinstance(head_ref, str) and
              head_ref.startswith("refs/")),
             "direct Git symbolic HEAD is inconsistent")
    if head_ref is not None:
        _require_safe_unicode(head_ref, "retained direct Git symbolic HEAD")
    superproject = value.get("superproject_worktree")
    _require(superproject is None or
             (isinstance(superproject, str) and
              Path(superproject).is_absolute()),
             "direct Git superproject worktree is invalid")
    if superproject is not None:
        _require_safe_unicode(
            superproject, "retained direct Git superproject worktree")
    if _expected_superproject is None:
        _require(superproject is None,
                 "top-level direct capture unexpectedly names a superproject")
    else:
        _require(superproject == _expected_superproject,
                 "retained direct submodule names a different superproject")

    refs_raw = _validate_byte_payload(
        value.get("refs"), "retained direct Git references")
    refs = _parse_direct_refs(refs_raw)
    _require(not any(
                 record["name"].startswith("refs/replace/")
                 for record in refs),
             "retained direct Git source contains replacement refs")
    if head_ref is not None:
        matching = [record for record in refs if record["name"] == head_ref]
        _require(len(matching) == 1 and
                 matching[0]["object_id"] == value["head"],
                 "retained direct Git HEAD differs from its reference")

    commit = _validate_object_identity(
        value.get("commit_object"), "commit", value["head"],
        "retained direct Git commit")
    _require(b"\n\n" in commit,
             "retained direct Git commit has no header/message boundary")
    headers = commit.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + value["tree"].encode("ascii")
    _require(headers and headers[0] == tree_line and
             [line for line in headers if line.startswith(b"tree ")] ==
             [tree_line],
             "retained direct Git commit names a different tree")
    raw_tree_objects = value.get("tree_objects")
    _require(isinstance(raw_tree_objects, list) and
             0 < len(raw_tree_objects) <= MAX_GIT_TREE_OBJECTS,
             "retained direct Git tree-object closure is invalid")
    tree_objects: dict[str, bytes] = {}
    tree_total = 0
    previous_object_id: str | None = None
    for tree_record in raw_tree_objects:
        object_id = (
            tree_record.get("object_id")
            if isinstance(tree_record, dict) else None)
        _require(isinstance(object_id, str) and
                 HEX40.fullmatch(object_id) is not None and
                 (previous_object_id is None or
                  previous_object_id < object_id) and
                 object_id not in tree_objects,
                 "retained direct Git tree-object closure is non-canonical")
        content = _validate_object_identity(
            tree_record, "tree", object_id,
            "retained direct Git tree object")
        tree_total += len(content)
        _require(tree_total <= MAX_GIT_TREE_TOTAL_BYTES,
                 "retained direct Git tree-object closure exceeds its bound")
        tree_objects[object_id] = content
        previous_object_id = object_id
    tree_entries = _flatten_tree_objects(value["tree"], tree_objects)
    expected_tree_listing = b"".join(
        (
            f"{entry['git_mode']} {entry['git_type']} "
            f"{entry['object_id']}\t{entry['path']}\0"
        ).encode("utf-8")
        for entry in tree_entries)
    _require(_validate_byte_payload(
                 value.get("head_tree_listing"),
                 "retained direct Git tree listing") ==
             expected_tree_listing,
             "retained direct Git listing differs from its tree objects")

    git_executable = _validate_sealed_file_proof(
        value.get("git_executable"), GIT_EXECUTABLE_PROTOCOL,
        "retained direct Git executable", require_source_executable=True,
        expected_source_path=GIT_EXECUTABLE_SOURCE_PATH)
    if _expected_git_executable is not None:
        _require(_canonical_bytes(git_executable) ==
                 _canonical_bytes(_expected_git_executable),
                 "retained direct submodule used a different Git executable")
    else:
        _expected_git_executable = git_executable
    metadata = value.get("git_metadata")
    _require(isinstance(metadata, dict) and set(metadata) == {
                 "layout", "gitdir", "commondir", "guarded_components",
                 "git_entry", "commondir_file", "guard_policy",
                 "guarded_files", "guarded_file_count",
                 "guarded_files_sha256", "legacy_grafts_absent"} and
             metadata.get("layout") in {"ordinary", "linked-worktree"} and
             isinstance(metadata.get("gitdir"), str) and
             Path(metadata["gitdir"]).is_absolute() and
             isinstance(metadata.get("commondir"), str) and
             Path(metadata["commondir"]).is_absolute() and
             isinstance(metadata.get("guarded_components"), list) and
             metadata["guarded_components"] ==
                 sorted(set(metadata["guarded_components"])) and
             set(metadata["guarded_components"]) ==
                 {metadata["gitdir"], metadata["commondir"]} and
             metadata.get("guard_policy") == METADATA_GUARD_POLICY and
             type(metadata.get("guarded_file_count")) is int and
             metadata["guarded_file_count"] >= 0 and
             isinstance(metadata.get("guarded_files_sha256"), str) and
             HEX256.fullmatch(metadata["guarded_files_sha256"]) is not None and
             metadata.get("legacy_grafts_absent") is True and
             value.get("worktree_guard_policy") ==
                 DIRECT_WORKTREE_GUARD_POLICY,
             "retained direct Git metadata identity is invalid")
    for metadata_path in metadata["guarded_components"]:
        _require(isinstance(metadata_path, str) and
                 Path(metadata_path).is_absolute(),
                 "retained direct Git metadata path is invalid")
        _require_safe_unicode(
            metadata_path, "retained direct Git metadata path")
    git_entry = metadata.get("git_entry")
    if git_entry is None:
        _require(metadata["layout"] == "ordinary" and
                 metadata["gitdir"] == str(Path(value["path"]) / ".git"),
                 "retained ordinary Git layout differs from its .git entry")
    else:
        _require(metadata["layout"] == "linked-worktree" and
                 _validate_git_path_reference(
                     git_entry, Path(value["path"]),
                     "retained direct .git file", prefix=b"gitdir: ") ==
                     metadata["gitdir"],
                 "retained linked Git layout differs from its .git file")
    commondir_file = metadata.get("commondir_file")
    if commondir_file is None:
        _require(metadata["commondir"] == metadata["gitdir"],
                 "retained direct Git commondir is not derived")
    else:
        _require(_validate_git_path_reference(
                     commondir_file, Path(metadata["gitdir"]),
                     "retained direct Git commondir file") ==
                     metadata["commondir"],
                 "retained direct Git commondir file differs")
    guarded_files = metadata.get("guarded_files")
    _require(isinstance(guarded_files, list) and
             guarded_files == sorted(set(guarded_files)),
             "retained direct Git metadata inventory is non-canonical")
    for guarded_file in guarded_files:
        _require(isinstance(guarded_file, str),
                 "retained direct Git metadata inventory entry is invalid")
        role, separator, relative = guarded_file.partition(":")
        _require(bool(separator) and role in {"gitdir", "commondir"} and
                 _safe_path(
                     relative.encode("utf-8"),
                     "retained direct Git metadata inventory path") ==
                     relative,
                 "retained direct Git metadata inventory entry differs")
    _require(metadata["guarded_file_count"] == len(guarded_files) and
             metadata["guarded_files_sha256"] == _digest(guarded_files),
             "retained direct Git metadata inventory summary differs")

    config_raw = _validate_byte_payload(
        value.get("config"), "retained direct Git configuration")
    _validate_direct_config(config_raw)
    index = value.get("index")
    _require(isinstance(index, dict) and set(index) == {
                 "source", "sealed", "entry_count", "stage", "flags_v",
                 "flags_f", "raw"} and
             type(index.get("entry_count")) is int and
             0 <= index["entry_count"] <= MAX_TRACKED_SOURCE_FILES,
             "retained direct Git index identity is invalid")
    expected_index_path = str(Path(metadata["gitdir"]) / "index")
    _validate_sealed_file_proof(
        {"source": index["source"], "sealed": index["sealed"]},
        DIRECT_INDEX_PROTOCOL, "retained direct Git index",
        require_source_executable=False,
        expected_source_path=expected_index_path)
    raw_index = _validate_byte_payload(
        index.get("raw"), "retained raw direct Git index")
    _require(len(raw_index) == index["source"]["size"] and
             hashlib.sha256(raw_index).hexdigest() ==
                 index["source"]["sha256"],
             "retained raw direct Git index differs from its sealed file")
    raw_index_entries = _parse_raw_index(raw_index)
    stage_raw = _validate_byte_payload(
        index.get("stage"), "retained direct Git staged index")
    flags_v_raw = _validate_byte_payload(
        index.get("flags_v"), "retained direct Git -v index")
    flags_f_raw = _validate_byte_payload(
        index.get("flags_f"), "retained direct Git -f index")
    index_entries = _parse_index(stage_raw)
    index_paths = [entry["path"] for entry in index_entries]
    _require(index_entries == raw_index_entries and
             len(index_entries) == index["entry_count"] and
             _parse_default_flags(
                 flags_v_raw, "retained direct Git -v index") == index_paths and
             _parse_default_flags(
                 flags_f_raw, "retained direct Git -f index") == index_paths,
             "retained direct Git index transcripts are inconsistent")
    status_raw = _validate_byte_payload(
        value.get("status"), "retained direct Git status")
    status_records = _parse_direct_status(status_raw)
    untracked_paths = _parse_direct_path_list(
        _validate_byte_payload(
            value.get("untracked_paths"),
            "retained direct Git untracked path inventory"),
        "retained direct Git untracked path inventory")

    inline_paths = value.get("inline_paths")
    _require(isinstance(inline_paths, list) and
             inline_paths == sorted(set(inline_paths)) and
             all(isinstance(path, str) and
                 _safe_path(
                     path.encode("utf-8"),
                     "retained direct inline source path") == path
                 for path in inline_paths),
             "retained direct inline source paths are invalid")
    inline_path_set = set(inline_paths)
    if _expected_superproject is not None:
        _require(not inline_path_set,
                 "retained direct submodule has unexpected inline payloads")
    records = value.get("paths")
    _require(isinstance(records, list) and records and
             len(records) <= 2 * MAX_TRACKED_SOURCE_FILES and
             value.get("paths_sha256") == _digest(records),
             "retained direct Git path inventory is invalid")
    head_by_path = {entry["path"]: entry for entry in tree_entries}
    index_by_path = {entry["path"]: entry for entry in index_entries}
    expected_paths = sorted(set(head_by_path) | set(index_by_path))
    observed_paths: list[str] = []
    observed_inline_payloads: set[str] = set()
    total = 0
    for record in records:
        _require(isinstance(record, dict) and set(record) == {
                     "path", "head_mode", "head_object", "index_mode",
                     "index_object", "worktree"} and
                 isinstance(record.get("path"), str),
                 "retained direct Git path record shape differs")
        path = _safe_path(
            record["path"].encode("utf-8"),
            "retained direct Git path")
        observed_paths.append(path)
        expected_head_entry = head_by_path.get(path)
        expected_index_entry = index_by_path.get(path)
        _require(record.get("head_mode") == (
                     expected_head_entry["git_mode"]
                     if expected_head_entry is not None else None) and
                 record.get("head_object") == (
                     expected_head_entry["object_id"]
                     if expected_head_entry is not None else None) and
                 record.get("index_mode") == (
                     expected_index_entry["git_mode"]
                     if expected_index_entry is not None else None) and
                 record.get("index_object") == (
                     expected_index_entry["object_id"]
                     if expected_index_entry is not None else None),
                 "retained direct Git path differs from HEAD or index")
        worktree = record.get("worktree")
        _require(isinstance(worktree, dict) and
                 worktree.get("kind") in {
                     "regular", "symlink", "missing", "submodule"},
                 "retained direct Git worktree record is invalid")
        kind = worktree["kind"]
        if kind == "missing":
            _require(set(worktree) == {"kind"} and
                     record.get("index_mode") != "160000",
                     "retained missing worktree record differs")
        elif kind == "regular":
            expected_worktree_keys = {"kind", "mode", "git_object"}
            if path in inline_path_set:
                expected_worktree_keys.add("payload")
            _require(set(worktree) == expected_worktree_keys and
                     worktree.get("mode") in {"100644", "100755"} and
                     isinstance(worktree.get("git_object"), str) and
                     HEX40.fullmatch(worktree["git_object"]) is not None,
                     "retained direct regular-file identity is invalid")
            if path in inline_path_set:
                content = _validate_byte_payload(
                    worktree["payload"],
                    f"retained inline direct file {path!r}")
                _require(len(content) <= MAX_TRACKED_SOURCE_FILE_BYTES and
                         _sha1_object("blob", content) ==
                             worktree["git_object"],
                         "retained inline direct file bytes differ")
                total += len(content)
                observed_inline_payloads.add(path)
        elif kind == "symlink":
            _require(set(worktree) == {
                         "kind", "mode", "size", "sha256", "git_object",
                         "target_encoding", "target_base64"} and
                     worktree.get("mode") == "120000" and
                     type(worktree.get("size")) is int and
                     0 <= worktree["size"] <=
                         MAX_TRACKED_SOURCE_FILE_BYTES and
                     isinstance(worktree.get("sha256"), str) and
                     HEX256.fullmatch(worktree["sha256"]) is not None and
                     isinstance(worktree.get("git_object"), str) and
                     HEX40.fullmatch(worktree["git_object"]) is not None and
                     worktree.get("target_encoding") == "base64" and
                     isinstance(worktree.get("target_base64"), str),
                     "retained direct symlink identity is invalid")
            try:
                target = base64.b64decode(
                    worktree["target_base64"], validate=True)
            except (ValueError, binascii.Error) as error:
                raise GitCaptureError(
                    "retained direct symlink target is invalid base64") \
                    from error
            _require(base64.b64encode(target).decode("ascii") ==
                     worktree["target_base64"] and
                     len(target) == worktree["size"] and
                     hashlib.sha256(target).hexdigest() ==
                         worktree["sha256"] and
                     _sha1_object("blob", target) ==
                         worktree["git_object"],
                     "retained direct symlink bytes differ")
            total += worktree["size"]
        else:
            _require(set(worktree) == {
                         "kind", "identity_sha256", "identity"} and
                     record.get("head_mode") in {None, "160000"} and
                     record.get("index_mode") in {None, "160000"} and
                     isinstance(worktree.get("identity_sha256"), str) and
                     HEX256.fullmatch(
                         worktree["identity_sha256"]) is not None and
                     isinstance(worktree.get("identity"), dict) and
                     worktree["identity_sha256"] ==
                         _digest(worktree["identity"]),
                     "retained direct submodule identity is invalid")
            nested_path = os.path.abspath(os.path.join(value["path"], path))
            validate_direct_git_capture(
                worktree["identity"], nested_path, require_clean=False,
                _depth=_depth + 1,
                _expected_git_executable=_expected_git_executable,
                _expected_superproject=value["path"])
    _require(inline_path_set == observed_inline_payloads,
             "retained direct inline source declaration differs from its "
             "regular payloads")
    _require(observed_paths == expected_paths and
             observed_paths == sorted(set(observed_paths)),
             "retained direct Git path inventory is non-canonical")
    expected_status = _derive_direct_status_records(
        {record["path"]: record for record in records},
        untracked_paths,
        core_filemode=_direct_boolean_config(
            config_raw, "core.filemode", default=True),
    )
    _require(status_records == expected_status and
             status_raw == _encode_direct_status(expected_status),
             "retained direct Git status is not canonically derived from "
             "its exact layers")
    _require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
             "retained direct worktree bytes exceed their total bound")
    derived_clean = _direct_worktree_is_clean(
        tree_entries, index_entries, records, status_raw)
    _require(value["clean"] == derived_clean,
             "retained direct Git cleanliness is not derived from its layers")
    _require(not require_clean or derived_clean,
             "retained direct Git source is not clean")
    return value


def validate_git_capture(
    value: object, expected_path: str, expected_head: str, *,
    require_detached: bool,
    _depth: int = 0,
    _expected_git_executable: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate the portable, current-schema capture without filesystem I/O."""
    _require(_depth <= MAX_SUBMODULE_DEPTH,
             "retained Git submodule identity exceeds its recursion bound")
    expected_keys = {
        "schema", "path", "head", "tree", "detached", "head_ref",
        "superproject_worktree", "tracked_tree_listing_sha256",
        "tracked_status", "commit_object", "tree_objects", "git_executable",
        "git_metadata", "worktree_guard_policy", "config", "index",
        "tracked_files", "tracked_files_sha256", "submodules",
    }
    _require(isinstance(value, dict) and set(value) == expected_keys and
             value.get("schema") == SCHEMA,
             "Git source capture shape or schema differs")
    _require(isinstance(expected_path, str) and Path(expected_path).is_absolute()
             and value.get("path") == expected_path and
             isinstance(expected_head, str) and
             HEX40.fullmatch(expected_head) is not None and
             value.get("head") == expected_head and
             isinstance(value.get("tree"), str) and
             HEX40.fullmatch(value["tree"]) is not None and
             type(value.get("detached")) is bool and
             (not require_detached or value["detached"] is True) and
             value.get("tracked_status") == "clean" and
             isinstance(value.get("tracked_tree_listing_sha256"), str) and
             HEX256.fullmatch(value["tracked_tree_listing_sha256"]) is not None,
             "Git source capture identity is invalid")
    _require_safe_unicode(expected_path, "retained Git source root")
    head_ref = value.get("head_ref")
    _require((value["detached"] and head_ref is None) or
             (not value["detached"] and isinstance(head_ref, str) and
              head_ref.startswith("refs/")),
             "Git source capture symbolic HEAD is inconsistent")
    superproject = value.get("superproject_worktree")
    _require(superproject is None or
             (isinstance(superproject, str) and
              Path(superproject).is_absolute()),
             "Git source capture superproject worktree is invalid")
    if superproject is not None:
        _require_safe_unicode(
            superproject, "retained Git superproject worktree")

    commit = _validate_object_identity(
        value.get("commit_object"), "commit", expected_head,
        "retained Git commit")
    _require(b"\n\n" in commit,
             "retained Git commit has no header/message boundary")
    headers = commit.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + value["tree"].encode("ascii")
    _require(headers and headers[0] == tree_line and
             [line for line in headers if line.startswith(b"tree ")] ==
             [tree_line],
             "retained Git commit names a different tree")
    raw_tree_objects = value.get("tree_objects")
    _require(isinstance(raw_tree_objects, list) and
             0 < len(raw_tree_objects) <= MAX_GIT_TREE_OBJECTS,
             "retained Git tree-object closure is invalid")
    tree_objects: dict[str, bytes] = {}
    tree_total = 0
    previous_object_id: str | None = None
    for record in raw_tree_objects:
        object_id = record.get("object_id") if isinstance(record, dict) else None
        _require(isinstance(object_id, str) and
                 HEX40.fullmatch(object_id) is not None and
                 (previous_object_id is None or
                  previous_object_id < object_id) and
                 object_id not in tree_objects,
                 "retained Git tree-object closure is not canonical")
        content = _validate_object_identity(
            record, "tree", object_id, "retained Git tree object")
        tree_total += len(content)
        _require(tree_total <= MAX_GIT_TREE_TOTAL_BYTES,
                 "retained Git tree-object closure exceeds its byte bound")
        tree_objects[object_id] = content
        previous_object_id = object_id
    tree_entries = _flatten_tree_objects(value["tree"], tree_objects)

    executable = value.get("git_executable")
    _require(isinstance(executable, dict) and
             set(executable) == {"source", "sealed"},
             "retained Git executable identity shape differs")
    source = executable["source"]
    sealed = executable["sealed"]
    _require(isinstance(source, dict) and set(source) == {
                 "path", "size", "mode", "sha256"} and
             source.get("path") == GIT_EXECUTABLE_SOURCE_PATH and
             type(source.get("size")) is int and source["size"] > 0 and
             type(source.get("mode")) is int and
             0 <= source["mode"] <= 0o7777 and source["mode"] & 0o111 and
             isinstance(source.get("sha256"), str) and
             HEX256.fullmatch(source["sha256"]) is not None,
             "retained Git source executable identity is invalid")
    _require_safe_unicode(source["path"], "retained Git executable path")
    _require(isinstance(sealed, dict) and set(sealed) == {
                 "protocol", "size", "mode", "sha256", "seals",
                 "source_sha256"} and
             sealed.get("protocol") == GIT_EXECUTABLE_PROTOCOL and
             sealed.get("size") == source["size"] and
             sealed.get("sha256") == source["sha256"] and
             sealed.get("source_sha256") == source["sha256"] and
             type(sealed.get("mode")) is int and sealed["mode"] & 0o111 and
             type(sealed.get("seals")) is int and
             sealed["seals"] & REQUIRED_SEALS == REQUIRED_SEALS,
             "sealed Git executable identity is invalid")
    if _expected_git_executable is None:
        _expected_git_executable = executable
    else:
        _require(
            _canonical_bytes(executable) ==
                _canonical_bytes(_expected_git_executable),
            "retained submodule used a different Git executable")

    metadata = value.get("git_metadata")
    _require(isinstance(metadata, dict) and set(metadata) == {
                 "layout", "gitdir", "commondir", "guarded_components",
                 "guard_policy", "guarded_file_count",
                 "guarded_files_sha256"} and
             metadata.get("layout") in {"ordinary", "linked-worktree"} and
             isinstance(metadata.get("gitdir"), str) and
             Path(metadata["gitdir"]).is_absolute() and
             isinstance(metadata.get("commondir"), str) and
             Path(metadata["commondir"]).is_absolute() and
             isinstance(metadata.get("guarded_components"), list) and
             all(isinstance(component, str)
                 for component in metadata["guarded_components"]) and
             metadata["guarded_components"] ==
                 sorted(set(metadata["guarded_components"])) and
             set(metadata["guarded_components"]) ==
                 {metadata["gitdir"], metadata["commondir"]} and
             metadata.get("guard_policy") == METADATA_GUARD_POLICY and
             type(metadata.get("guarded_file_count")) is int and
             metadata["guarded_file_count"] >= 0 and
             isinstance(metadata.get("guarded_files_sha256"), str) and
             HEX256.fullmatch(metadata["guarded_files_sha256"]) is not None and
             value.get("worktree_guard_policy") == WORKTREE_GUARD_POLICY,
             "retained Git metadata identity is invalid")
    ordinary_gitdir = str(Path(value["path"]) / ".git")
    _require(
        (metadata["layout"] == "ordinary" and
         metadata["gitdir"] == ordinary_gitdir) or
        (metadata["layout"] == "linked-worktree" and
         metadata["gitdir"] != ordinary_gitdir),
        "retained Git metadata layout differs from its source root")
    for metadata_path in metadata["guarded_components"]:
        _require(isinstance(metadata_path, str) and
                 Path(metadata_path).is_absolute(),
                 "retained Git metadata component path is invalid")
        _require_safe_unicode(
            metadata_path, "retained Git metadata component path")
    _validate_byte_identity(value.get("config"), "retained Git config")
    index = value.get("index")
    _require(isinstance(index, dict) and set(index) == {
                 "entry_count", "stage", "flags_v", "flags_f"} and
             type(index.get("entry_count")) is int and
             0 <= index["entry_count"] <= MAX_TRACKED_SOURCE_FILES,
             "retained Git index identity is invalid")
    for name in ("stage", "flags_v", "flags_f"):
        _validate_byte_identity(index.get(name), f"retained Git index {name}")

    records = value.get("tracked_files")
    _require(isinstance(records, list) and
             len(records) == index["entry_count"] and
             len(records) <= MAX_TRACKED_SOURCE_FILES and
             value.get("tracked_files_sha256") == _digest(records),
             "retained tracked-file inventory is invalid")
    paths: list[str] = []
    submodule_records: dict[str, dict[str, Any]] = {}
    for record in records:
        _require(isinstance(record, dict) and
                 record.get("kind") in {
                     "regular", "symlink", "submodule"},
                 "retained tracked-file record is invalid")
        kind = record["kind"]
        common_keys = {
            "path", "git_mode", "git_type", "object_id", "kind"}
        if kind == "regular":
            # The object ID is authenticated by the retained recursive tree
            # closure.  Do not retain an independently editable size/SHA-256
            # transcript here: without the blob bytes those fields cannot be
            # derived offline and would create a second, unauthenticated truth.
            expected_record_keys = common_keys
            valid_mode = record.get("git_mode") in {"100644", "100755"}
        elif kind == "symlink":
            expected_record_keys = common_keys | {
                "size", "sha256", "target_encoding", "target_base64"}
            valid_mode = record.get("git_mode") == "120000"
        else:
            expected_record_keys = common_keys | {"identity_sha256"}
            valid_mode = record.get("git_mode") == "160000"
        _require(set(record) == expected_record_keys and valid_mode and
                 record.get("git_type") ==
                    ("commit" if kind == "submodule" else "blob") and
                 isinstance(record.get("object_id"), str) and
                 HEX40.fullmatch(record["object_id"]) is not None and
                 isinstance(record.get("path"), str),
                 "retained tracked-file record shape differs")
        path = _safe_path(
            record["path"].encode("utf-8"), "retained tracked path")
        paths.append(path)
        if kind == "symlink":
            _require(type(record.get("size")) is int and
                     0 <= record["size"] <= MAX_TRACKED_SOURCE_FILE_BYTES and
                     isinstance(record.get("sha256"), str) and
                     HEX256.fullmatch(record["sha256"]) is not None,
                     "retained tracked byte identity is invalid")
        if kind == "symlink":
            _require(record.get("target_encoding") == "base64" and
                     isinstance(record.get("target_base64"), str),
                     "retained tracked symlink identity is invalid")
            try:
                target = base64.b64decode(
                    record["target_base64"], validate=True)
            except (ValueError, binascii.Error) as error:
                raise GitCaptureError(
                    "retained tracked symlink target is invalid base64") \
                    from error
            _require(base64.b64encode(target).decode("ascii") ==
                     record["target_base64"] and
                     len(target) == record["size"] and
                     hashlib.sha256(target).hexdigest() == record["sha256"] and
                     _sha1_object("blob", target) == record["object_id"],
                     "retained tracked symlink bytes differ")
        if kind == "submodule":
            _require(isinstance(record.get("identity_sha256"), str) and
                     HEX256.fullmatch(record["identity_sha256"]) is not None,
                     "retained submodule digest is invalid")
            submodule_records[path] = record
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "retained tracked-file inventory is non-canonical")
    _require([
        {
            "path": record["path"],
            "git_mode": record["git_mode"],
            "git_type": record["git_type"],
            "object_id": record["object_id"],
        }
        for record in records
    ] == tree_entries,
        "retained tracked-file inventory differs from recursive tree objects")
    tree_listing = b"".join(
        (
            f"{record['git_mode']} {record['git_type']} "
            f"{record['object_id']}\t{record['path']}\0"
        ).encode("utf-8")
        for record in records)
    index_stage = b"".join(
        (
            f"{record['git_mode']} {record['object_id']} 0\t"
            f"{record['path']}\0"
        ).encode("utf-8")
        for record in records)
    default_flags = b"".join(
        f"H {record['path']}\0".encode("utf-8")
        for record in records)
    _require(
        value["tracked_tree_listing_sha256"] ==
            hashlib.sha256(tree_listing).hexdigest() and
        index["stage"] == _byte_identity(index_stage) and
        index["flags_v"] == _byte_identity(default_flags) and
        index["flags_f"] == _byte_identity(default_flags),
        "retained Git tree/index transcripts do not match tracked records")

    submodules = value.get("submodules")
    _require(isinstance(submodules, list) and
             len(submodules) == len(submodule_records),
             "retained submodule inventory differs from Git links")
    observed_submodule_paths: list[str] = []
    for record in submodules:
        _require(isinstance(record, dict) and set(record) == {
                     "path", "object_id", "identity_sha256", "identity"} and
                 isinstance(record.get("path"), str) and
                 record["path"] in submodule_records and
                 record.get("object_id") ==
                    submodule_records[record["path"]]["object_id"] and
                 record.get("identity_sha256") ==
                    submodule_records[record["path"]]["identity_sha256"] and
                 record["identity_sha256"] == _digest(record.get("identity")),
                 "retained submodule record is inconsistent")
        nested_path = os.path.abspath(
            os.path.join(value["path"], record["path"]))
        validate_git_capture(
            record["identity"], nested_path, record["object_id"],
            require_detached=False, _depth=_depth + 1,
            _expected_git_executable=_expected_git_executable)
        observed_submodule_paths.append(record["path"])
    _require(observed_submodule_paths == sorted(observed_submodule_paths) and
             len(observed_submodule_paths) ==
                 len(set(observed_submodule_paths)),
             "retained submodule inventory is non-canonical")
    return value
