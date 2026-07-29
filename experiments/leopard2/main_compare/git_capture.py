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
GIT_EXECUTABLE_PROTOCOL = "linux-sealed-git-executable-memfd/v1"
METADATA_GUARD_POLICY = \
    "retained-gitdir-commondir-recursive-inotify/v1"
WORKTREE_GUARD_POLICY = \
    "retained-root-and-tracked-path-components-inotify/v1"
MAX_GIT_COMMIT_BYTES = 16 * 1024 * 1024
MAX_GIT_TREE_BYTES = 64 * 1024 * 1024
MAX_GIT_TREE_OBJECTS = 4 * MAX_TRACKED_SOURCE_FILES
MAX_GIT_TREE_TOTAL_BYTES = MAX_GIT_TREE_BYTES
MAX_GIT_LISTING_BYTES = 64 * 1024 * 1024
MAX_GIT_CONFIG_BYTES = 16 * 1024 * 1024
MAX_SUBMODULE_DEPTH = 16
HEX40 = re.compile(r"[0-9a-f]{40}")
HEX256 = re.compile(r"[0-9a-f]{64}")
REQUIRED_SEALS = 0x0001 | 0x0002 | 0x0004 | 0x0008


class GitCaptureError(RuntimeError):
    """The repository could not be captured as one coherent source state."""


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
    for encoded in raw.split(b"\0"):
        if not encoded:
            continue
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
    result: list[dict[str, str]] = []
    reachable: set[str] = set()
    stack: list[tuple[str, object, str, frozenset[str]]] = [
        ("tree", root_tree, "", frozenset())]
    while stack:
        kind, payload, prefix, active = stack.pop()
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
        _require(isinstance(object_id, str) and object_id not in active and
                 object_id in objects,
                 "retained Git tree-object closure is cyclic or incomplete")
        reachable.add(object_id)
        next_active = active | {object_id}
        entries = _parse_tree_object(
            objects[object_id], f"retained Git tree object {object_id}")
        for entry in reversed(entries):
            relative = (
                f"{prefix}/{entry['name']}" if prefix else entry["name"])
            _safe_path(
                relative.encode("utf-8"), "retained Git recursive tree path")
            if entry["git_type"] == "tree":
                stack.append((
                    "tree", entry["object_id"], relative, next_active))
            else:
                stack.append(("leaf", entry, relative, next_active))
    _require(reachable == set(objects),
             "retained Git tree-object closure contains unreachable objects")
    paths = [entry["path"] for entry in result]
    _require(paths == sorted(paths) and len(paths) == len(set(paths)),
             "retained Git recursive tree inventory is non-canonical")
    return result


def _parse_index(raw: bytes) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for encoded in raw.split(b"\0"):
        if not encoded:
            continue
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


def _parse_default_flags(raw: bytes, label: str) -> list[str]:
    paths: list[str] = []
    for encoded in raw.split(b"\0"):
        if not encoded:
            continue
        _require(encoded.startswith(b"H "),
                 "Git index uses assume-unchanged, skip-worktree, "
                 "fsmonitor-valid, or another non-default flag")
        paths.append(_safe_path(encoded[2:], label))
    _require(len(paths) == len(set(paths)) and paths == sorted(paths),
             f"{label} is duplicated or non-canonical")
    return paths


def _validate_local_config(raw: bytes) -> None:
    for record in raw.split(b"\0"):
        if not record:
            continue
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
) -> bytes:
    """Run one patchable query through the exact sealed Git executable."""
    _require(git_snapshot.resolved is not None,
             "retained Git executable pathname was lost")
    metadata.verify()
    git_snapshot.verify()
    inherited = (
        *inherited_descriptors, root_descriptor, metadata.descriptor)
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
    )
    metadata.verify()
    git_snapshot.verify()
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
                    root, topdown=True, followlinks=False):
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
                    Path("/usr/bin/git"),
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


def validate_git_capture(
    value: object, expected_path: str, expected_head: str, *,
    require_detached: bool,
    _depth: int = 0,
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
             isinstance(source.get("path"), str) and
             Path(source["path"]).is_absolute() and
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
            require_detached=False, _depth=_depth + 1)
        observed_submodule_paths.append(record["path"])
    _require(observed_submodule_paths == sorted(observed_submodule_paths) and
             len(observed_submodule_paths) ==
                 len(set(observed_submodule_paths)),
             "retained submodule inventory is non-canonical")
    return value
