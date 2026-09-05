#!/usr/bin/python3
"""Dormant v19 lock and artifact owners (leopard-79h.38.5.4.8.2.2.2).

No CLI, build, qualification or benchmark dispatch exists here. The future
orchestrator must first verify retained lineage and own the source/build stages.
These owners bind output bytes to the existing preregistration, not a claimed
new build history. They do not establish runtime-library closure or authorize
execution. Host captures are boundary observations, not continuous authority.
"""

from __future__ import annotations

from contextlib import ExitStack
import copy
import fcntl
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import stat


HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path):
    if path.resolve(strict=True) != path:
        raise RuntimeError("v19 owner dependency is not canonical")
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


host = _load("v19_owned_artifacts_host", HERE / "v19_host_preflight.py")
provenance = _load("v19_owned_artifacts_provenance",
                   HERE.parents[2] / "tools/leopard2_build_provenance.py")
require = host.require
LOCK_PATH = "/tmp/leopard-gf8-authoritative.lock"
MAX_ARTIFACT_BYTES = 16 << 20
MAX_TOTAL_BYTES = 32 << 20
ARTIFACTS = {
    "candidate_executable": ("candidate", "bench_leopard2", "candidate_binary_sha256", 0o500),
    "candidate_archive": ("candidate", "libleopard.a", "candidate_archive_sha256", 0o400),
    "baseline_executable": ("baseline", "leopard_main_benchmark", "baseline_binary_sha256", 0o500),
    "baseline_archive": ("baseline", "libleopard_main_exact.a", "baseline_archive_sha256", 0o400),
}


def owned_directory(path: Path, *, private: bool = False) -> tuple[int, int, int]:
    """Check every path component without following directory symlinks."""
    descriptor = host.LinuxReader.open_directory(str(path))
    try:
        value = os.fstat(descriptor)
        mode = stat.S_IMODE(value.st_mode)
        require(value.st_uid == os.getuid() and value.st_gid == os.getgid() and
                (mode == 0o700 if private else mode & 0o022 == 0),
                "artifact directory has unsafe owner or permissions")
        return value.st_dev, value.st_ino, mode
    finally:
        os.close(descriptor)


def verify_current_bytes(snapshot, expected: str) -> None:
    """Rehash the held file, not its cached bytes (mmap need not notify inotify)."""
    snapshot._path_guard.verify()
    snapshot.verify()
    size = snapshot.identity["size"]
    require(0 < size <= MAX_ARTIFACT_BYTES, "retained artifact exceeds byte bound")
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        part = os.pread(snapshot.descriptor, min(1 << 20, size - offset), offset)
        require(part, "retained artifact ended while hashed")
        digest.update(part)
        offset += len(part)
    require(digest.hexdigest() == expected, "current artifact bytes differ from preregistration")
    snapshot.verify()
    snapshot._path_guard.verify()


def validate_flock_record(text: str, status: os.stat_result, pid: int) -> None:
    """Require the exclusive FLOCK attached to this open file description.

    Linux /proc/PID/fdinfo exposes lock rows in /proc/locks format. Do not
    reacquire the flock here: doing so would hide an earlier release.
    """
    rows = [line.split() for line in text.splitlines() if line.startswith("lock:")]
    require(len(rows) == 1, "canonical descriptor lacks one exclusive lock")
    fields = rows[0]
    require(len(fields) == 9 and fields[0] == "lock:" and
            re.fullmatch(r"[0-9]+:", fields[1]) is not None and
            fields[2:6] == ["FLOCK", "ADVISORY", "WRITE", str(pid)] and
            fields[7:] == ["0", "EOF"], "canonical descriptor lock mode or owner differs")
    match = re.fullmatch(r"([0-9a-f]+):([0-9a-f]+):([0-9]+)", fields[6])
    require(match is not None and
            (int(match[1], 16), int(match[2], 16), int(match[3])) ==
            (os.major(status.st_dev), os.minor(status.st_dev), status.st_ino),
            "canonical descriptor lock inode differs")


class _LockMutationGuard(provenance._InotifyMutationGuard):
    # A competing cooperative worker may open O_RDWR, fail flock, and close
    # without writing. That produces CLOSE_WRITE but is not a mutation. Actual
    # writes, metadata changes, rename/unlink, and event overflow still reject.
    _DIRECTORY_MASK = provenance._InotifyMutationGuard._DIRECTORY_MASK & ~provenance.IN_CLOSE_WRITE
    _FILE_MASK = provenance._InotifyMutationGuard._FILE_MASK & ~provenance.IN_CLOSE_WRITE


class CanonicalLock:
    """Own a nonblocking flock and retain its pathname mutation history.

    The descriptor is private, close-on-exec, and never lent to children.
    Revalidation detects an observed release/downgrade or pathname ABA; it
    cannot detect a malicious same-process unlock/relock between observations.
    """

    def __init__(self, *, _path: str = LOCK_PATH):
        host.canonical_path(_path, nonroot=True)
        self.path = Path(_path)
        self._stack = ExitStack()
        self._file = provenance._OwnedDescriptor()
        self._guard = None
        self._identity = None
        self._pid = os.getpid()
        self._state = "new"

    @staticmethod
    def _metadata(value: os.stat_result) -> tuple:
        require(stat.S_ISREG(value.st_mode) and value.st_nlink == 1 and
                value.st_uid == os.getuid() and value.st_gid == os.getgid() and
                stat.S_IMODE(value.st_mode) == 0o600,
                "canonical lock has unsafe metadata")
        return provenance._stable_fields(value)

    def __enter__(self):
        require(self._state == "new", "canonical lock owner cannot be reused")
        self._state = "entering"
        try:
            self._stack.enter_context(self._file)
            directory = host.LinuxReader.open_directory(str(self.path.parent))
            try:
                descriptor = self._file.open(
                    self.path.name, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW |
                    os.O_CLOEXEC | os.O_NONBLOCK, 0o600, dir_fd=directory)
            finally:
                os.close(directory)
            self._identity = self._metadata(os.fstat(descriptor))
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self._guard = self._stack.enter_context(
                _LockMutationGuard("v19 canonical lock"))
            self._guard.add_file_path(self.path)
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "closed"
            self._stack.close()
            raise

    def validate_current(self) -> None:
        require(self._state == "held" and os.getpid() == self._pid,
                "canonical lock is not held by this owner")
        self._guard.verify()
        status = os.fstat(self._file.descriptor)
        require(self._metadata(status) == self._identity and
                self._metadata(os.lstat(self.path)) == self._identity and
                not os.get_inheritable(self._file.descriptor),
                "canonical lock identity changed")
        text = host.LinuxReader().read_text(
            f"/proc/{self._pid}/fdinfo/{self._file.descriptor}")
        validate_flock_record(text, status, self._pid)
        self._guard.verify()

    def __exit__(self, kind, value, traceback):
        try:
            self.validate_current()
        finally:
            self._state = "closed"
            # Closing our sole description releases the flock. Never unlock
            # an integer that might now name somebody else's description.
            self._stack.__exit__(kind, value, traceback)


class BuildArtifactLease:
    """Hold the build-stage lock, refresh host observations, and freeze bytes.

    Private injection arguments exist only for offline tests. No retained JSON
    capture can be supplied as live host authority. Keep this context open
    through every consumer and terminal sealing; returned FDs are borrowed.
    """

    def __init__(self, preregistration_bytes: bytes, *,
                 _reader_factory=None, _lock_path: str = LOCK_PATH):
        self._preregistration_bytes = preregistration_bytes
        self._preregistration = host.load_preregistration(preregistration_bytes)
        require(self._preregistration["artifact_execution"]["canonical_lock"] == LOCK_PATH,
                "canonical lock differs from preregistration")
        self._reader_factory = host.LinuxReader if _reader_factory is None else _reader_factory
        self._lock = CanonicalLock(_path=_lock_path)
        self._stack = ExitStack()
        self._state = "new"
        self._initial_host = None
        self._lane = None
        self._lane_identity = None
        self._artifact_root = None
        self._artifact_identity = None
        self._files = {}
        self._artifact_guard = None

    def _capture(self) -> dict:
        return host.capture_host(self._preregistration_bytes, self._reader_factory())

    def __enter__(self):
        require(self._state == "new", "artifact owner cannot be reused")
        self._state = "entering"
        try:
            self._initial_host = self._capture()  # Before lock or output creation.
            self._stack.enter_context(self._lock)
            self._state = "locked"
            self.validate_current()
            return self
        except BaseException:
            self._state = "closed"
            self._stack.close()
            raise

    def validate_current(self) -> None:
        require(self._state in ("locked", "frozen"), "artifact lease is not usable")
        self._lock.validate_current()
        current = self._capture()
        require(current["observations"][0] == self._initial_host["observations"][0] and
                current["file_identity_sha256"] == self._initial_host["file_identity_sha256"],
                "host/resource identity differs from the current lease")
        if self._state == "frozen":
            require(owned_directory(self._lane, private=True) == self._lane_identity and
                    owned_directory(self._artifact_root) == self._artifact_identity,
                    "lane or frozen artifact directory changed")
            self._artifact_guard.verify()
            for role, snapshot in self._files.items():
                expected = self._preregistration["build_preflight"][ARTIFACTS[role][2]]
                verify_current_bytes(snapshot, expected)
            self._artifact_guard.verify()
        self._lock.validate_current()

    def freeze(self, lane: Path, build_roots: dict[str, Path]) -> dict:
        """Copy fixed build outputs once; partial failures are retained, not reused.

        The caller must separately own/validate the detached sources, build
        recipes, child containment, and prior-attempt lineage. This method is
        an output-byte gate, not a substitute for those acquisition stages.
        """
        require(self._state == "locked", "artifact freeze is duplicate or out of order")
        self.validate_current()
        require(type(build_roots) is dict and set(build_roots) == {"candidate", "baseline"},
                "build roots differ")
        lane = Path(lane)
        roots = {role: Path(path) for role, path in build_roots.items()}
        self._lane_identity = owned_directory(lane, private=True)
        for root in roots.values():
            owned_directory(root)
        paths = [lane, *roots.values()]
        require(all(a != b and a not in b.parents and b not in a.parents
                    for i, a in enumerate(paths) for b in paths[i + 1:]),
                "lane and build roots alias or overlap")
        destination = lane / "v19-artifacts"
        self._lane, self._artifact_root = lane, destination
        try:
            with ExitStack() as inputs:
                snapshots = {}
                total = 0
                for role, (build_role, name, key, _mode) in ARTIFACTS.items():
                    path = roots[build_role] / name
                    metadata = os.lstat(path)
                    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                            metadata.st_uid == os.getuid() and metadata.st_gid == os.getgid() and
                            not stat.S_IMODE(metadata.st_mode) & 0o7022 and
                            0 < metadata.st_size <= MAX_ARTIFACT_BYTES,
                            "build output has unsafe type, owner, links, mode or size")
                    total += metadata.st_size
                    require(total <= MAX_TOTAL_BYTES, "build output total exceeds byte bound")
                    snapshot = inputs.enter_context(provenance._RetainedFileSnapshot(
                        path, "v19 input " + role, maximum_bytes=MAX_ARTIFACT_BYTES))
                    require(provenance._stable_fields(metadata) ==
                            provenance._stable_fields(os.fstat(snapshot.descriptor)) and
                            snapshot.resolved == path and
                            snapshot.identity["sha256"] == self._preregistration["build_preflight"][key],
                            "build output differs from pinned preflight bytes")
                    magic = b"\x7fELF" if role.endswith("executable") else b"!<arch>\n"
                    require(snapshot.content.startswith(magic), "build output format differs")
                    snapshots[role] = snapshot
                self.validate_current()
                require(owned_directory(lane, private=True) == self._lane_identity,
                        "lane changed before artifact creation")
                parent = host.LinuxReader.open_directory(str(lane))
                try:
                    os.mkdir("v19-artifacts", 0o700, dir_fd=parent)
                    with provenance._OwnedDescriptor() as directory:
                        directory.open("v19-artifacts", os.O_RDONLY | os.O_DIRECTORY |
                                       os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=parent)
                        for role, snapshot in snapshots.items():
                            with provenance._OwnedDescriptor() as output:
                                descriptor = output.open(role, os.O_WRONLY | os.O_CREAT |
                                    os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
                                    0o600, dir_fd=directory.descriptor)
                                remaining = memoryview(snapshot.content)
                                while remaining:
                                    count = os.write(descriptor, remaining[:1 << 20])
                                    require(count > 0, "artifact write stalled")
                                    remaining = remaining[count:]
                                os.fchmod(descriptor, ARTIFACTS[role][3])
                                os.fsync(descriptor)
                        os.fchmod(directory.descriptor, 0o500)
                        os.fsync(directory.descriptor)
                finally:
                    os.close(parent)
                for role, snapshot in snapshots.items():
                    verify_current_bytes(snapshot,
                        self._preregistration["build_preflight"][ARTIFACTS[role][2]])
            self._artifact_identity = owned_directory(destination)
            require(self._artifact_identity[2] == 0o500,
                    "artifact directory is not sealed")
            self._artifact_guard = self._stack.enter_context(
                provenance._InotifyMutationGuard("v19 frozen artifacts"))
            self._artifact_guard.add_tree(destination)
            descriptor = host.LinuxReader.open_directory(str(destination))
            try:
                names = set()
                with os.scandir(descriptor) as entries:
                    for entry in entries:
                        require(entry.name in ARTIFACTS, "unexpected artifact entry")
                        names.add(entry.name)
                require(names == set(ARTIFACTS), "frozen artifact inventory differs")
            finally:
                os.close(descriptor)
            for role in ARTIFACTS:
                snapshot = self._stack.enter_context(provenance._RetainedFileSnapshot(
                    destination / role, "v19 frozen " + role, maximum_bytes=MAX_ARTIFACT_BYTES))
                metadata = os.fstat(snapshot.descriptor)
                require(metadata.st_nlink == 1 and metadata.st_uid == os.getuid() and
                        metadata.st_gid == os.getgid() and
                        stat.S_IMODE(metadata.st_mode) == ARTIFACTS[role][3] and
                        snapshot.identity["sha256"] ==
                        self._preregistration["build_preflight"][ARTIFACTS[role][2]],
                        "frozen artifact identity differs")
                if role.endswith("executable"):
                    snapshot.executable_descriptor  # Allocate a sealed immutable copy now.
                self._files[role] = snapshot
            self._state = "frozen"
            self.validate_current()
            return self.record()
        except BaseException:
            self._state = "failed"
            raise

    def record(self) -> dict:
        require(self._state == "frozen", "artifacts have not been frozen")
        self.validate_current()
        preflight = self._preregistration["build_preflight"]
        return copy.deepcopy({
            "schema": "leopard2-v19-owned-build-artifacts/v1",
            "preregistration_sha256": host.PREREGISTRATION_SHA256,
            "preregistered_source": {key: preflight[key] for key in
                ("source_commit", "source_tree", "baseline_commit", "baseline_tree")},
            "artifacts": {role: {"path": str(self._artifact_root / role),
                "sha256": snapshot.identity["sha256"], "size": len(snapshot.content)}
                for role, snapshot in self._files.items()},
            "live_acquisition_armed": False,
            "source_build_history_proved": False,
            "runtime_closure_proved": False,
            "resource_lifetime_proved": False,
        })

    def executable_descriptor(self, role: str) -> int:
        require(self._state == "frozen" and role in ("candidate", "baseline"),
                "executable request is premature or invalid")
        self.validate_current()
        return self._files[role + "_executable"].executable_descriptor

    def __exit__(self, kind, value, traceback):
        try:
            if self._state in ("locked", "frozen"):
                self.validate_current()
            elif value is None:
                raise host.PreflightError("failed artifact lease cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
