#!/usr/bin/python3
"""Retain the pinned v19 build preflight; leopard-79h.38.5.4.8.2.2.2.

This is a read-only historical-input owner, not a builder or an execution
authority. It never runs a retained script, Git, a compiler, or a benchmark.
The orchestrator must hold its live host/lock owner around this context and
separately verify v18 lineage, detached sources, fresh builds and runtime files.
"""

from __future__ import annotations

from contextlib import ExitStack
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat


HERE = Path(__file__).resolve().parent
_owner_path = HERE / "v19_owned_artifacts.py"
if _owner_path.resolve(strict=True) != _owner_path:
    raise RuntimeError("retained preflight owner dependency is not canonical")
_spec = importlib.util.spec_from_file_location(
    "v19_retained_preflight_owners", _owner_path)
owners = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(owners)
host, provenance = owners.host, owners.provenance
require = host.require
MAX_FILES = 64
MAX_FILE_BYTES = 16 << 20
MAX_TOTAL_BYTES = 32 << 20
MAX_JSON_BYTES = 1 << 20
REPOSITORY = Path("/home/catid/leopard")


def checksum_entries(data: bytes) -> dict[str, str]:
    require(type(data) is bytes and 0 < len(data) <= 16384 and data.endswith(b"\n"),
            "preflight checksum manifest is empty, oversized or unterminated")
    result = {}
    for line in data.decode("ascii", errors="strict").split("\n")[:-1]:
        require(len(result) < MAX_FILES, "preflight checksum inventory exceeds bound")
        match = re.fullmatch(r"([0-9a-f]{64})  \./([A-Za-z0-9_./-]{1,256})", line)
        require(match is not None, "preflight checksum row is not canonical")
        digest, path = match.groups()
        components = path.split("/")
        require(len(components) <= 8 and all(part not in ("", ".", "..") for part in components)
                and path != "SHA256SUMS" and path not in result,
                "preflight checksum path is unsafe or duplicate")
        result[path] = digest
    require(list(result) == sorted(result), "preflight checksum entries are not sorted")
    require(not any(parent in result for path in result for parent in
                    (str(value) for value in PurePosixPath(path).parents if str(value) != ".")),
            "preflight checksum files overlap directories")
    return result


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "preflight JSON has duplicate keys")
        result[key] = value
    return result


def _json(data: bytes) -> dict:
    require(0 < len(data) <= MAX_JSON_BYTES, "preflight JSON exceeds byte bound")
    value = json.loads(data, object_pairs_hook=_unique_object)
    require(type(value) is dict, "preflight JSON is not an object")
    return value


def _subset(value, expected: dict, label: str) -> None:
    require(type(value) is dict and all(key in value for key in expected) and
            host.canonical_bytes({key: value[key] for key in expected}) == host.canonical_bytes(expected),
            label + " differs from frozen preflight")


class PinnedPreflight:
    """Hold every manifest-listed file and rehash it through the retained fd.

    The original preflight directories are owner-only, not write-protected;
    their exact inventories and mutation histories remain guarded throughout.
    Private root injection is for synthetic filesystem tests, not a CLI option.
    """

    def __init__(self, preregistration_bytes: bytes, *, _root: Path | None = None):
        self._preregistration = host.load_preregistration(preregistration_bytes)
        self._contract = self._preregistration["build_preflight"]
        self.root = REPOSITORY / self._contract["retained_lane"] if _root is None else Path(_root)
        host.canonical_path(str(self.root), nonroot=True)
        self._stack = ExitStack()
        self._files = {}
        self._directories = {}
        self._digests = {}
        self._guard = None
        self._state = "new"
        self._total_bytes = 0
        self._run_root = None

    @staticmethod
    def _directory(path: Path) -> tuple:
        descriptor = host.LinuxReader.open_directory(str(path))
        try:
            value = os.fstat(descriptor)
            require(value.st_uid == os.getuid() and value.st_gid == os.getgid() and
                    stat.S_IMODE(value.st_mode) in (0o700, 0o500),
                    "preflight directory is not owner-only")
            return value.st_dev, value.st_ino, value.st_mode, value.st_uid, value.st_gid
        finally:
            os.close(descriptor)

    def _retain(self, relative: str):
        path = self.root / relative
        descriptor = host.LinuxReader.open_directory(str(path.parent))
        try:
            value = os.stat(path.name, dir_fd=descriptor, follow_symlinks=False)
        finally:
            os.close(descriptor)
        require(stat.S_ISREG(value.st_mode) and value.st_uid == os.getuid() and
                value.st_gid == os.getgid() and value.st_nlink == 1 and
                stat.S_IMODE(value.st_mode) == 0o444 and
                0 <= value.st_size <= MAX_FILE_BYTES,
                "preflight file has unsafe type, owner, mode, links or size")
        self._total_bytes += value.st_size
        require(self._total_bytes <= MAX_TOTAL_BYTES, "preflight total exceeds byte bound")
        snapshot = self._stack.enter_context(provenance._RetainedFileSnapshot(
            path, "v19 retained " + relative, maximum_bytes=MAX_FILE_BYTES))
        require(snapshot.resolved == path and provenance._stable_fields(value) ==
                provenance._stable_fields(os.fstat(snapshot.descriptor)),
                "preflight file changed while opened")
        require(snapshot.identity["sha256"] == self._digests[relative],
                "preflight file hash differs: " + relative)
        self._files[relative] = snapshot
        return snapshot

    def _inventory(self, *, retain: bool = False) -> None:
        children = {"": set()}
        for relative in self._digests:
            parts = relative.split("/")
            for index, part in enumerate(parts):
                children.setdefault("/".join(parts[:index]), set()).add(part)
        for relative, expected in children.items():
            path = self.root / relative
            if retain:
                self._guard._add_watch(path, self._guard._DIRECTORY_MASK |
                    provenance.IN_ONLYDIR | provenance.IN_ATTRIB, None)
            identity = self._directory(path)
            if retain:
                self._directories[relative] = identity
            else:
                require(identity == self._directories[relative], "preflight directory changed")
            descriptor = host.LinuxReader.open_directory(str(path))
            try:
                observed = set()
                with os.scandir(descriptor) as entries:
                    for entry in entries:
                        require(len(observed) < MAX_FILES + 1 and entry.name in expected,
                                "unexpected preflight tree entry")
                        observed.add(entry.name)
                require(observed == expected, "preflight tree inventory differs")
            finally:
                os.close(descriptor)

    def __enter__(self):
        require(self._state == "new", "preflight owner cannot be reused")
        self._state = "entering"
        try:
            initial_root = self._directory(self.root)
            self._digests["SHA256SUMS"] = self._contract["outer_sha256sums_sha256"]
            manifest = self._retain("SHA256SUMS")
            self._digests.update(checksum_entries(manifest.content))
            self._guard = self._stack.enter_context(
                provenance._InotifyMutationGuard("v19 retained preflight tree"))
            self._guard.add_file_path(self.root / "SHA256SUMS")
            self._inventory(retain=True)
            require(self._directories[""] == initial_root, "preflight root changed while opened")
            for relative in self._digests:
                if relative != "SHA256SUMS":
                    self._retain(relative)
            self._validate_semantics()
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "closed"
            self._stack.close()
            raise

    def _bytes(self, relative: str) -> bytes:
        require(relative in self._files, "required preflight file is missing: " + relative)
        return self._files[relative].content

    def _validate_semantics(self) -> None:
        preflight = self._contract
        source = {"candidate_commit": preflight["source_commit"],
                  "candidate_tree": preflight["source_tree"],
                  "baseline_commit": preflight["baseline_commit"],
                  "baseline_tree": preflight["baseline_tree"]}
        identities = "".join(f"{key}={value}\n" for key, value in source.items()).encode("ascii")
        require(self._bytes("source-identities.txt") == identities,
                "retained source identity file differs")
        require(self._bytes("terminal.env") == b"exit_status=0\nstage=complete\n" + identities +
                b"memory_max=536870912\nmemory_swap_max=0\nbenchmark_invoked=false\n" and
                self._bytes("build-stage.txt") == b"complete\n" and
                self._bytes("provenance-validator-nofile65536.status") == b"0\n",
                "retained preflight did not complete successfully")
        run_root = self._bytes("run-root.txt").decode("ascii", errors="strict")
        require(re.fullmatch(r"/home/catid/leopard-v19-build-preflight\.[A-Za-z0-9]{8}\n", run_root)
                is not None, "retained canonical build root differs")
        self._run_root = run_root[:-1]
        result = _json(self._bytes("RESULT.json"))
        _subset(result, {"schema": "leopard2-v19-gf16-build-preflight-result/v1",
                         "status": "passed", "canonical_lock_held": True,
                         "v19_attempt_consumed": False, "source": source}, "preflight result")
        _subset(result.get("host"), {"hostname": "ripper", "shared_host": True, "swap_rows": 0},
                "preflight host")
        _subset(result.get("build"), {"configuration": "Release", "generator": "Unix Makefiles",
            "compiler": "c++ 13.3.0", "parallelism": 1, "same_path_replay": True,
            "candidate_target": "bench_leopard2", "baseline_target": "leopard_main_benchmark",
            **{role + "_byte_identical": True for role in owners.ARTIFACTS}}, "preflight build")
        events = self._bytes("memory.events")
        host.memory_events(events.decode("ascii"))
        require(self._bytes("memory.events.local") == events,
                "retained memory event inventories differ")
        peak = host.decimal(self._bytes("memory.peak").decode("ascii"), "retained memory peak")
        current = host.decimal(self._bytes("memory.current").decode("ascii"), "retained memory usage")
        require(0 <= current <= peak == preflight["observed_memory_peak_bytes"] <= host.MEMORY_MAX,
                "retained memory peak or current usage differs")
        _subset(result.get("resource_envelope"), {"benchmark_invoked": False,
            "memory_max_bytes": host.MEMORY_MAX, "memory_swap_max_bytes": 0,
            "memory_peak_bytes": peak, "oom": 0, "oom_kill": 0,
            "memory_events_sha256": hashlib.sha256(events).hexdigest()}, "preflight resources")
        _subset(result.get("provenance"), {"candidate_record_sha256": preflight["provenance_record_sha256"],
            "schema": "leopard2-production-build-closure/v2", "validated": True,
            "validated_nofile": 65536}, "preflight provenance result")
        artifact_digests = {}
        for role, (build_role, _filename, key, _mode) in owners.ARTIFACTS.items():
            expected = preflight[key]
            artifact_digests[role + "_sha256"] = expected
            name = build_role + ("" if role.endswith("executable") else "-archive")
            for generation in ("first", "replay"):
                relative = f"artifacts/{name}-{generation}"
                require(self._digests.get(relative) == expected and
                        self._files[relative].identity["size"] > 0,
                        "first/replay artifact differs from pinned preflight")
        _subset(result.get("artifacts"), artifact_digests, "preflight artifacts")
        require(self._digests.get("candidate-build-provenance.json") == preflight["provenance_record_sha256"],
                "candidate provenance record hash differs")
        record = _json(self._bytes("candidate-build-provenance.json"))
        _subset(record, {"schema": "leopard2-production-build-closure/v2",
            "source_root": self._run_root + "/candidate-source",
            "physical_source_root": self._run_root + "/candidate-source",
            "build_root": self._run_root + "/candidate-build", "executable_target": "bench_leopard2"},
            "retained candidate build paths")
        manifest = record.get("tracked_source_manifest")
        require(type(manifest) is dict, "retained source manifest is missing")
        _subset(manifest.get("git"), {"commit": source["candidate_commit"], "tree": source["candidate_tree"],
            "dirty": False, "status_sha256": hashlib.sha256(b"").hexdigest()}, "retained source manifest")
        for kind, filename, key in (("executable", "bench_leopard2", "candidate_binary_sha256"),
                                    ("archive", "libleopard.a", "candidate_archive_sha256")):
            _subset(record.get(kind), {"path": self._run_root + "/candidate-build/" + filename,
                                      "sha256": preflight[key]}, "retained candidate " + kind)

    def validate_current(self) -> None:
        require(self._state == "held", "retained preflight is not held")
        self._guard.verify()
        self._inventory()
        for relative, snapshot in self._files.items():
            snapshot._path_guard.verify()
            snapshot.verify()
            current = provenance._read_bounded_descriptor(snapshot.descriptor,
                snapshot.identity["size"], "current preflight " + relative, MAX_FILE_BYTES)
            require(hashlib.sha256(current).hexdigest() == self._digests[relative],
                    "current retained preflight hash differs: " + relative)
            snapshot.verify()
            snapshot._path_guard.verify()
        self._guard.verify()

    def record(self) -> dict:
        self.validate_current()
        return copy.deepcopy({
            "schema": "leopard2-v19-retained-build-preflight-capture/v1",
            "root": str(self.root), "preregistration_sha256": host.PREREGISTRATION_SHA256,
            "outer_sha256sums_sha256": self._contract["outer_sha256sums_sha256"],
            "provenance_record_sha256": self._contract["provenance_record_sha256"],
            "file_count": len(self._files), "total_bytes": self._total_bytes,
            "canonical_run_root": self._run_root,
            "source_pins": {key: self._contract[key] for key in
                ("source_commit", "source_tree", "baseline_commit", "baseline_tree")},
            "retained_preflight_verified": True, "new_build_verified": False,
            "live_acquisition_armed": False, "v18_archives_verified": False,
        })

    def __exit__(self, kind, value, traceback):
        try:
            self.validate_current()
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
