#!/usr/bin/python3
"""Physical v18 ancestry custody; leopard-79h.38.5.4.8.2.2.2.2.

Authenticate every retained archive byte without executing historical code or
caching the archive bodies. Disclosure labels are not physical path authority.
This verifies custody of pinned historical evidence, not a new failure replay,
runtime closure, host/lock ownership, or permission to acquire timings.
"""
from __future__ import annotations

from contextlib import ExitStack
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat

HERE = Path(__file__).resolve().parent
_dependency = HERE / "v19_source_identity.py"
if _dependency.resolve(strict=True) != _dependency:
    raise RuntimeError("lineage dependency is not canonical")
_spec = importlib.util.spec_from_file_location("v19_lineage_identity", _dependency)
identity = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(identity)
streamed, preflight, owners = identity.streamed, identity.preflight, identity.owners
host, provenance, require = identity.host, identity.provenance, identity.require

SOURCE_COMMIT = "c8f825d0a033d31d220b0ebce9cc8871e8c2fc6d"
SOURCE_TREE = "2c17a0a7bcea20274d2593cb204442c4c817e464"
ARCHIVES = (
    (1, "c8f825d-v18-passive-main-a1", "ce65c3a49ef1c1d89ba51ea03d0af4742d6790e6f2ea2662917d9ef9a9d945d7"),
    (2, "c8f825d-v18-passive-main-a2", "a1bf0eda157c251f33f7260ebd76931d88054d460bd07a97bcba2811384b2c10"),
    (3, "c8f825d-v18-passive-main-a3", "fe5b40cc98753cbd794ee019cb0e2643d0ccee0aca4c5fd7b2e0b27df8a86139"),
)
MAX_FILES, MAX_DIRECTORIES = 6144, 2048
MAX_FILE_BYTES, MAX_TOTAL_BYTES, MAX_JSON_BYTES = 128 << 20, 512 << 20, 4 << 20
BLOCK_BYTES = 64 << 10


def lineage_record():
    return {"schema": "leopard2-main-compare-v18-failure-lineage/v1", "source_commit": SOURCE_COMMIT,
            "source_tree": SOURCE_TREE, "attempts": [{"attempt": attempt,
                "envelope": ".research/leopard-79h/" + name,
                "envelope_sha256sums_sha256": digest} for attempt, name, digest in ARCHIVES]}


def checksums(data):
    require(type(data) is bytes and 0 < len(data) <= 512 << 10 and data.endswith(b"\n"),
            "lineage checksum manifest is empty, oversized or unterminated")
    result = {}
    for line in data.decode("ascii").split("\n")[:-1]:
        match = re.fullmatch(r"([0-9a-f]{64})  \./([A-Za-z0-9_./+-]{1,4096})", line)
        require(match is not None and len(result) < MAX_FILES, "invalid lineage checksum row")
        digest, path = match.groups()
        require(len(path.split("/")) <= 24 and all(part not in ("", ".", "..") for part in path.split("/"))
                and path not in result and path != "SHA256SUMS", "unsafe or duplicate lineage checksum path")
        result[path] = digest
    require(list(result) == sorted(result), "lineage checksums are not sorted")
    return result


def _json(data):
    require(type(data) is bytes and 0 < len(data) <= MAX_JSON_BYTES, "lineage JSON exceeds bound")
    value = json.loads(data, object_pairs_hook=preflight._unique_object)
    require(type(value) is dict, "lineage JSON is not an object")
    return value


class PinnedV18Lineage:
    """Hold the three pinned archive trees beneath a caller-selected parent.

    No input tree is modified, copied, evicted, or executed. All nodes remain
    open and no-follow anchored; source bytes are rehashed at boundaries to
    cover persistent mmap writes that need not notify inotify. Observations
    are sequential, not atomic or a defense against arbitrary same-process
    manipulation. The caller separately owns the host and canonical lock.
    """
    def __init__(self, preregistration_bytes: bytes, archive_parent: Path):
        self.preregistration = host.load_preregistration(preregistration_bytes)
        self.parent = Path(archive_parent)
        host.canonical_path(str(self.parent), nonroot=True)
        self.lineage = lineage_record()
        self.lineage_sha256 = hashlib.sha256(host.canonical_bytes(self.lineage) + b"\n").hexdigest()
        require(self.lineage_sha256 == self.preregistration["attempt_contract"]["failure_lineage_sha256"],
                "physical lineage inventory differs from v19 preregistration")
        self._stack = ExitStack()
        self._state = "new"
        self._pid = os.getpid()
        self._files, self._directories, self._digests = {}, {}, {}
        self._inodes = set()
        self._total_bytes = 0
        self._guard = None
        self._parent_fd = -1
        self._parent_identity = None

    def _own(self, descriptor):
        try: self._stack.callback(os.close, descriptor)
        except BaseException:
            os.close(descriptor)
            raise
        return descriptor

    @staticmethod
    def _safe(value, directory):
        require(value.st_uid == os.geteuid() and value.st_gid == os.getegid() and
                not stat.S_IMODE(value.st_mode) & 0o7222 and
                (stat.S_ISDIR(value.st_mode) if directory else
                 stat.S_ISREG(value.st_mode) and value.st_nlink == 1),
                "lineage node is not a sealed, owned, ordinary single-link file/directory")

    @staticmethod
    def _tree_entry(relative, value, directory):
        return {"path": relative or ".", "type": "directory" if directory else "file",
                "mode": format(stat.S_IMODE(value.st_mode), "04o"), "nlink": value.st_nlink,
                "uid": value.st_uid, "gid": value.st_gid}

    def _metadata(self, entry):
        descriptor, parent, name, fields, _record = entry
        require(not os.get_inheritable(descriptor) and provenance._stable_fields(os.fstat(descriptor)) == fields ==
                provenance._stable_fields(os.stat(name, dir_fd=parent, follow_symlinks=False)),
                "retained lineage node metadata or pathname changed")

    def _hash(self, entry):
        self._metadata(entry)
        size = os.fstat(entry[0]).st_size
        digest = hashlib.sha256()
        for offset in range(0, size, BLOCK_BYTES):
            part = os.pread(entry[0], min(BLOCK_BYTES, size - offset), offset)
            require(len(part) == min(BLOCK_BYTES, size - offset), "retained lineage file truncated")
            digest.update(part)
        self._metadata(entry)
        return digest.hexdigest()

    def _bytes(self, key):
        entry = self._files[key]
        self._metadata(entry)
        size = os.fstat(entry[0]).st_size
        require(size <= MAX_JSON_BYTES, "lineage semantic input exceeds bound")
        # Only small semantic records are materialized; syscall reads retain
        # the same bound as archive hashing (including the ~2.2 MiB failure).
        parts = []
        for offset in range(0, size, BLOCK_BYTES):
            part = os.pread(entry[0], min(BLOCK_BYTES, size - offset), offset)
            require(len(part) == min(BLOCK_BYTES, size - offset), "lineage semantic input truncated")
            parts.append(part)
        value = b"".join(parts)
        require(len(value) == size and hashlib.sha256(value).hexdigest() == self._digests[key],
                "lineage semantic bytes differ")
        self._metadata(entry)
        return value

    def _visit(self, archive, relative, parent, name, depth=0):
        require(depth <= 24 and len(relative.encode()) <= 4096, "lineage tree exceeds path bound")
        before = os.stat(name, dir_fd=parent, follow_symlinks=False)
        directory = stat.S_ISDIR(before.st_mode)
        self._safe(before, directory)
        require(len(self._directories) < MAX_DIRECTORIES if directory else len(self._files) < MAX_FILES,
                "lineage node count exceeds bound")
        descriptor = self._own(os.open(name, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_NONBLOCK |
                                      (os.O_DIRECTORY if directory else 0), dir_fd=parent))
        opened = os.fstat(descriptor)
        require(provenance._stable_fields(opened) == provenance._stable_fields(before),
                "lineage node changed while opened")
        key = archive + ("/" + relative if relative else "")
        path = self.parent / key
        entry = (descriptor, parent, name, provenance._stable_fields(opened),
                 self._tree_entry(relative, opened, directory))
        if directory:
            self._guard.add_directory_path(path)
            self._guard._add_watch(path, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, None)
            self._directories[key] = entry
            names = []
            with os.scandir(descriptor) as entries:
                for child in entries:
                    require(len(names) < MAX_FILES + MAX_DIRECTORIES and
                            re.fullmatch(r"[A-Za-z0-9_.+-]{1,255}", child.name) is not None and
                            child.name not in (".", ".."), "unsafe or excessive lineage directory entries")
                    names.append(child.name)
            for child in sorted(names):
                self._visit(archive, relative + "/" + child if relative else child, descriptor, child, depth + 1)
        else:
            require(0 <= opened.st_size <= MAX_FILE_BYTES, "lineage file exceeds bound")
            self._total_bytes += opened.st_size
            require(self._total_bytes <= MAX_TOTAL_BYTES, "lineage total byte count exceeds bound")
            inode = opened.st_dev, opened.st_ino
            require(inode not in self._inodes, "lineage files share an inode")
            self._inodes.add(inode)
            self._guard.add_file_path(path)
            self._files[key] = entry
            self._digests[key] = self._hash(entry)
        self._metadata(entry)

    def _semantics(self, attempt, name, expected_outer):
        prefix = name + "/"
        require(all(prefix + relative in self._files for relative in (
            "SHA256SUMS", "TREE-METADATA.json", "FAILED.json", "core/SHA256SUMS",
            "core/manifest.json", "core/attempt-lineage.json", "core/campaign/failure.json")),
            "v18 archive lacks required evidence files")
        require(self._digests[prefix + "SHA256SUMS"] == expected_outer, "v18 outer manifest differs")
        outer = checksums(self._bytes(prefix + "SHA256SUMS"))
        actual = {key[len(prefix):]: digest for key, digest in self._digests.items() if key.startswith(prefix)}
        require(actual == {**outer, "SHA256SUMS": expected_outer}, "v18 physical file inventory or hashes differ")
        metadata_raw = self._bytes(prefix + "TREE-METADATA.json")
        metadata = _json(metadata_raw)
        require(metadata_raw == host.canonical_bytes(metadata) + b"\n", "v18 tree metadata is not canonical")
        policy = {"uid": os.geteuid(), "gid": os.getegid(),
                  "rule": "every retained node has the invoking effective uid and gid"}
        expected_entries = [entry[4] for key, entry in {**self._directories, **self._files}.items()
                            if (key == name or key.startswith(prefix)) and key != prefix + "TREE-METADATA.json"]
        expected_metadata = {"schema": "leopard2-authoritative-tree-metadata/v1", "root": ".",
            "excluded_paths": ["TREE-METADATA.json"], "final_mode_policy": "observed mode with all write bits removed",
            "uid_gid_policy": policy, "self_policy": {"uid": policy["uid"], "gid": policy["gid"],
                "mode": "0400", "nlink": 1, "type": "file",
                "sha256_binding": "exactly one ./TREE-METADATA.json checksum entry"},
            "entries": sorted(expected_entries, key=lambda row: row["path"])}
        require(host.canonical_bytes(metadata) == host.canonical_bytes(expected_metadata) and
                self._files[prefix + "TREE-METADATA.json"][4]["mode"] == "0400", "v18 physical tree metadata differs")
        core = checksums(self._bytes(prefix + "core/SHA256SUMS"))
        # Historical core publication used find ! -name SHA256SUMS, excluding
        # nested checksum files too. The outer envelope includes and binds all
        # of them; reproduce that exact core dialect, not a new inventory rule.
        require(core == {key[5:]: digest for key, digest in outer.items()
                         if key.startswith("core/") and key.rsplit("/", 1)[-1] != "SHA256SUMS"},
                "v18 nested core checksums differ")
        common = {"status": "failed", "acquisition_generation": "passive-v2", "attempt": attempt,
                  "attempt_budget": 3, "source_commit": SOURCE_COMMIT, "source_tree": SOURCE_TREE,
                  "promotion_passed": False, "campaign_exit_status": 1, "failure_verified": True,
                  "attempt_lineage_sha256": outer["core/attempt-lineage.json"]}
        terminal = _json(self._bytes(prefix + "FAILED.json"))
        expected_terminal = {**common, "schema": "leopard2-v18-gfni-main-failed-envelope/v1",
                             "core_sha256sums_sha256": outer["core/SHA256SUMS"]}
        require(host.canonical_bytes(terminal) == host.canonical_bytes(expected_terminal), "v18 failure terminal differs")
        manifest = _json(self._bytes(prefix + "core/manifest.json"))
        preflight._subset(manifest, {**common, "schema": "leopard2-v18-gfni-main-passive-failed-core-manifest/v1",
            "failure_verify_status": 0, "baseline_commit": self.preregistration["build_preflight"]["baseline_commit"],
            "failure_sha256": outer["core/campaign/failure.json"], "canonical_lock": owners.LOCK_PATH,
            "cpu": 52, "sibling": 116}, "v18 core terminal")
        lineage = _json(self._bytes(prefix + "core/attempt-lineage.json"))
        expected_prior = [{"attempt": number, "envelope": "/home/catid/leopard/.research/leopard-79h/" + prior_name,
                           "envelope_sha256sums_sha256": digest, "terminal": "FAILED.json",
                           "terminal_schema": "leopard2-v18-gfni-main-failed-envelope/v1"}
                          for number, prior_name, digest in ARCHIVES if number < attempt]
        require(host.canonical_bytes(lineage) == host.canonical_bytes({"schema": "leopard2-v18-gfni-main-attempt-lineage/v1",
                "acquisition_generation": "passive-v2", "attempt": attempt, "attempt_budget": 3,
                "source_commit": SOURCE_COMMIT, "source_tree": SOURCE_TREE, "prior_attempts": expected_prior}),
                "v18 prior-attempt chain differs")
        failure = _json(self._bytes(prefix + "core/campaign/failure.json"))
        preflight._subset(failure, {"schema": "leopard2-main-compare-failure/v18", "status": "failed",
                                   "valid": False, "error_type": "EvidenceError"}, "v18 retained failure")

    def __enter__(self):
        require(self._state == "new", "lineage owner cannot be reused")
        self._state = "entering"
        try:
            self._parent_fd = self._own(host.LinuxReader.open_directory(str(self.parent)))
            value = os.fstat(self._parent_fd)
            require(value.st_uid == os.geteuid() and value.st_gid == os.getegid() and
                    not stat.S_IMODE(value.st_mode) & 0o022, "lineage parent has unsafe ownership or permissions")
            self._parent_identity = streamed._directory_identity(value)
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 retained v18 ancestry"))
            self._guard.add_directory_path(self.parent)
            self._guard._add_watch(self.parent, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, {b""})
            for attempt, name, digest in ARCHIVES:
                self._visit(name, "", self._parent_fd, name)
                self._semantics(attempt, name, digest)
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def validate_current(self):
        require(self._state == "held" and os.getpid() == self._pid, "lineage owner is not live")
        try:
            self._guard.verify()
            require(not os.get_inheritable(self._parent_fd) and
                    streamed._directory_identity(os.fstat(self._parent_fd)) == self._parent_identity ==
                    streamed._directory_identity(self.parent.lstat()), "lineage parent identity changed")
            for entry in self._directories.values(): self._metadata(entry)
            for key, entry in self._files.items():
                require(self._hash(entry) == self._digests[key], "current v18 archive bytes differ")
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-retained-v18-lineage/v1", "physical_parent": str(self.parent),
            "lineage": self.lineage, "lineage_sha256": self.lineage_sha256, "file_count": len(self._files),
            "directory_count": len(self._directories), "total_bytes": self._total_bytes,
            "physical_archives_verified": True, "historical_failures_replayed": False,
            "runtime_closure_verified": False, "atomic_snapshot": False, "live_acquisition_armed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise host.PreflightError("failed lineage owner cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
