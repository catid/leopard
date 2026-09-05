#!/usr/bin/python3
"""Bounded source-file lifetime owner; leopard-79h.38.5.4.8.2.2.2.1.

This owns already staged file bytes, not their Git authenticity or creation
history. The caller must hold the host/lock owner, validate exact source/Git
identities while this context is held, and separately own build/runtime stages.
No Git, compiler, workload, acquisition, or arming entry point exists here.
"""
from __future__ import annotations

from contextlib import ExitStack
import hashlib
import importlib.util
import os
from pathlib import Path
import stat

HERE = Path(__file__).resolve().parent
_dependency = HERE / "v19_retained_preflight.py"
if _dependency.resolve(strict=True) != _dependency:
    raise RuntimeError("streamed source dependency is not canonical")
_spec = importlib.util.spec_from_file_location("v19_streamed_preflight", _dependency)
preflight = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(preflight)
owners = preflight.owners
host, provenance = owners.host, owners.provenance
require = host.require

SOURCE_ROOTS = ("candidate-source", "leopard1-source")
HASH_BLOCK_BYTES = 64 << 10
MAX_FILE_BYTES = 64 << 20
MAX_TOTAL_BYTES = 256 << 20
MAX_FILES = 4096
MAX_DIRECTORIES = 1024
MAX_DEPTH = 32
MAX_PATH_BYTES = 4096


def _directory_identity(value):
    return value.st_dev, value.st_ino, value.st_mode, value.st_uid, value.st_gid


class StreamingSourceOwner:
    """Retain complete staged source/.git inventories with bounded hash buffers.

    Both fixed source roots must be direct children of an owner-only workspace.
    Their group-writable modes are allowed inside that private boundary; world
    writes, links, special files and shared file inodes are rejected. Initial
    ordinary Git directories are required, but their authenticity is a separate
    gate. Sibling build directories may be created without changing sources.

    Checks are sequential boundary observations, not an atomic snapshot or
    protection against arbitrary same-process manipulation. Inotify records
    ordinary mutations; current-byte hashing also detects persistent mmap drift.
    """

    def __init__(self, workspace: Path):
        self.workspace = Path(workspace)
        host.canonical_path(str(self.workspace), nonroot=True)
        self._stack = ExitStack()
        self._state = "new"
        self._pid = os.getpid()
        self._workspace_fd = -1
        self._workspace_identity = None
        self._guard = None
        self._directories = {}
        self._files = {}
        self._inodes = set()
        self._total_bytes = 0

    def _own_fd(self, descriptor):
        try:
            self._stack.callback(os.close, descriptor)
        except BaseException:
            os.close(descriptor)
            raise
        return descriptor

    @staticmethod
    def _owned(value):
        require(value.st_uid == os.getuid() and value.st_gid == os.getgid() and
                not stat.S_IMODE(value.st_mode) & 0o7002,
                "source entry has unsafe owner or mode")

    def _name(self, relative, name):
        provenance._require_safe_unicode(name, "streamed source entry")
        require(name not in ("", ".", "..") and "/" not in name and "\\" not in name and
                len(os.fsencode(name)) <= 255 and
                len(os.fsencode(relative)) <= MAX_PATH_BYTES and
                len(relative.split("/")) <= MAX_DEPTH,
                "source entry path exceeds its canonical bounds")

    def _check_workspace(self):
        require(os.getpid() == self._pid and not os.get_inheritable(self._workspace_fd),
                "source owner process or descriptor inheritance changed")
        require(_directory_identity(os.fstat(self._workspace_fd)) == self._workspace_identity,
                "private workspace descriptor changed")
        descriptor = host.LinuxReader.open_directory(str(self.workspace))
        try:
            require(_directory_identity(os.fstat(descriptor)) == self._workspace_identity,
                    "private workspace pathname changed")
        finally:
            os.close(descriptor)

    def _metadata(self, entry):
        descriptor, parent, name, fields = entry[:4]
        require(not os.get_inheritable(descriptor) and
                provenance._stable_fields(os.fstat(descriptor)) == fields and
                provenance._stable_fields(os.stat(name, dir_fd=parent, follow_symlinks=False)) == fields,
                "held source metadata or pathname changed")

    def _hash(self, entry):
        self._metadata(entry)
        descriptor, _parent, _name, fields = entry[:4]
        size = fields[6]
        require(0 <= size <= MAX_FILE_BYTES, "held source size exceeds bound")
        digest = hashlib.sha256()
        offset = 0
        while offset < size:
            part = os.pread(descriptor, min(HASH_BLOCK_BYTES, size - offset), offset)
            require(part, "held source ended while hashed")
            digest.update(part)
            offset += len(part)
        require(not os.pread(descriptor, 1, size), "held source grew while hashed")
        self._metadata(entry)
        return digest.hexdigest()

    def _scan(self, relative, parent, name):
        self._name(relative, name)
        before = os.stat(name, dir_fd=parent, follow_symlinks=False)
        self._owned(before)
        require(stat.S_ISREG(before.st_mode) or stat.S_ISDIR(before.st_mode),
                "source entry is a link or special file")
        require(before.st_dev == self._workspace_identity[0], "source crosses a filesystem boundary")
        inode = before.st_dev, before.st_ino
        require(inode not in self._inodes, "source entries share an inode")
        self._inodes.add(inode)
        directory = stat.S_ISDIR(before.st_mode)
        if directory:
            require(len(self._directories) < MAX_DIRECTORIES, "source directory count exceeds bound")
        else:
            require(len(self._files) < MAX_FILES and before.st_nlink == 1 and
                    0 <= before.st_size <= MAX_FILE_BYTES,
                    "source file count, links or size exceeds bound")
            self._total_bytes += before.st_size
            require(self._total_bytes <= MAX_TOTAL_BYTES, "source total exceeds byte bound")
        descriptor = self._own_fd(os.open(name, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW |
            (os.O_DIRECTORY if directory else os.O_NONBLOCK), dir_fd=parent))
        fields = provenance._stable_fields(before)
        entry = descriptor, parent, name, fields
        self._metadata(entry)
        path = self.workspace / relative
        if directory:
            self._guard._add_watch(path, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                provenance.IN_DONT_FOLLOW | provenance.IN_ATTRIB, None)
            self._directories[relative] = entry + (set(),)
            children = self._directories[relative][4]
            with os.scandir(descriptor) as entries:
                for child in entries:
                    require(len(children) < MAX_FILES + MAX_DIRECTORIES,
                            "source directory inventory exceeds bound")
                    children.add(child.name)
            for child_name in sorted(children):
                self._scan(relative + "/" + child_name, descriptor, child_name)
        else:
            self._guard._add_watch(path, self._guard._FILE_MASK, None)
            self._files[relative] = entry + (self._hash(entry),)
        self._metadata(entry)
        self._guard.verify()

    def __enter__(self):
        require(self._state == "new", "source owner cannot be reused")
        self._state = "entering"
        try:
            self._workspace_fd = self._own_fd(host.LinuxReader.open_directory(str(self.workspace)))
            value = os.fstat(self._workspace_fd)
            require(value.st_uid == os.getuid() and value.st_gid == os.getgid() and
                    stat.S_IMODE(value.st_mode) == 0o700, "source workspace is not owner-only")
            self._workspace_identity = _directory_identity(value)
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 streamed sources"))
            self._guard.add_directory_path(self.workspace)
            self._guard._add_watch(self.workspace, provenance.IN_ATTRIB | provenance.IN_ONLYDIR, set())
            for name in SOURCE_ROOTS:
                self._guard.add_directory_path(self.workspace / name)
                self._scan(name, self._workspace_fd, name)
                require(name + "/.git" in self._directories, "source lacks an ordinary Git directory")
            require(self._files, "staged source inventory is empty")
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "closed"
            self._stack.close()
            raise

    def _check_directories(self):
        for entry in self._directories.values():
            self._metadata(entry)
            observed = set()
            with os.scandir(entry[0]) as entries:
                for child in entries:
                    require(child.name in entry[4] and len(observed) < MAX_FILES + MAX_DIRECTORIES,
                            "source directory gained an entry")
                    observed.add(child.name)
            require(observed == entry[4], "source directory lost an entry")
            self._metadata(entry)

    def validate_current(self, *, evict_cache: bool = False):
        require(self._state == "held", "source owner is not usable")
        try:
            require(type(evict_cache) is bool, "cache-eviction selector is not boolean")
            self._check_workspace()
            self._guard.verify()
            self._check_directories()
            for entry in self._files.values():
                require(self._hash(entry) == entry[4], "current source bytes changed")
                if evict_cache:
                    # These exact held nlink=1 descriptors were captured beneath
                    # the private workspace; never reopen a pathname for eviction.
                    os.fsync(entry[0])
                    os.posix_fadvise(entry[0], 0, 0, os.POSIX_FADV_DONTNEED)
                    self._metadata(entry)
            self._check_directories()
            self._guard.verify()
            self._check_workspace()
        except BaseException:
            self._state = "failed"
            raise

    def record(self, *, evict_cache: bool = False):
        self.validate_current(evict_cache=evict_cache)
        files = [{"path": path, "size": entry[3][6], "mode": stat.S_IMODE(entry[3][2]),
                  "sha256": entry[4]} for path, entry in sorted(self._files.items())]
        return {"schema": "leopard2-v19-streamed-source-files/v1", "workspace": str(self.workspace),
            "files": files, "file_manifest_sha256": hashlib.sha256(host.canonical_bytes(files)).hexdigest(),
            "file_count": len(files), "directory_count": len(self._directories),
            "total_bytes": self._total_bytes, "hash_block_bytes": HASH_BLOCK_BYTES,
            "staged_file_ownership_verified": True, "cache_eviction_requested": evict_cache,
            "git_identity_verified": False, "source_creation_verified": False,
            "runtime_closure_verified": False, "atomic_snapshot": False,
            "live_acquisition_armed": False}

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held":
                self.validate_current()
            elif value is None:
                raise host.PreflightError("failed source owner cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
