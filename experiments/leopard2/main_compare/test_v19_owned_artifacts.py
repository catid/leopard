#!/usr/bin/python3
"""Offline real-descriptor tests for leopard-79h.38.5.4.8.2.2.2.

Host data and pinned executable bytes are synthetic; locks, temporary files,
inotify history and memfd seals use the real kernel. Nothing is executed.
"""

from __future__ import annotations

import copy
import fcntl
import hashlib
import importlib.util
import mmap
import os
from pathlib import Path
import stat
import subprocess
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def load(name, filename):
    specification = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


owner = load("v19_owner_tested", "v19_owned_artifacts.py")
fixtures = load("v19_owner_host_fixtures", "test_v19_host_preflight.py")
FAILURES = (owner.host.PreflightError, owner.provenance.BuildProvenanceError, OSError)


class ArtifactOwnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preregistration = subprocess.run(
            [str(HERE / "run_authoritative_v17_gfni_main_compare.sh"),
             "--print-conditioned-v19-preregistration"], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        cls.real_contract = owner.host.load_preregistration(cls.preregistration)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-owner-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.lane = self.root / "lane"
        self.lane.mkdir(mode=0o700)
        self.roots = {role: self.root / role for role in ("candidate", "baseline")}
        for path in self.roots.values():
            path.mkdir(mode=0o700)
        self.contract = copy.deepcopy(self.real_contract)
        for role, (build_role, name, key, _mode) in owner.ARTIFACTS.items():
            magic = b"\x7fELF" if role.endswith("executable") else b"!<arch>\n"
            data = magic + role.encode() + b"\0" * 64
            (self.roots[build_role] / name).write_bytes(data)
            (self.roots[build_role] / name).chmod(0o644)
            self.contract["build_preflight"][key] = hashlib.sha256(data).hexdigest()
        self.reader = fixtures.FakeReader()
        # Only the output hashes are synthetic. The test never exposes an
        # alternate production preregistration or host-authorization switch.
        real_loader = owner.host.load_preregistration
        def fixture_contract(data):
            real_loader(data)
            return copy.deepcopy(self.contract)
        patcher = mock.patch.object(owner.host, "load_preregistration", side_effect=fixture_contract)
        patcher.start()
        self.addCleanup(patcher.stop)

    def lease(self):
        return owner.BuildArtifactLease(self.preregistration,
            _reader_factory=lambda: copy.deepcopy(self.reader),
            _lock_path=str(self.root / "canonical.lock"))

    def test_real_lock_contention_and_release(self):
        path = str(self.root / "lock")
        with owner.CanonicalLock(_path=path) as first:
            first.validate_current()
            with self.assertRaises(BlockingIOError):
                with owner.CanonicalLock(_path=path):
                    self.fail("second description acquired held lock")
            descriptor = os.open(path, os.O_RDWR | os.O_CLOEXEC)
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                os.close(descriptor)
            first.validate_current()
            self.assertFalse(os.get_inheritable(first._file.descriptor))
        with owner.CanonicalLock(_path=path) as second:
            second.validate_current()
        with self.assertRaises(owner.host.PreflightError):
            first.validate_current()
        with self.assertRaises(owner.host.PreflightError):
            first.__enter__()

    def test_same_inode_unlock_or_downgrade_is_rejected(self):
        for operation in (fcntl.LOCK_UN, fcntl.LOCK_SH):
            with self.subTest(operation=operation), self.assertRaises(FAILURES):
                with owner.CanonicalLock(_path=str(self.root / "lock")) as lock:
                    fcntl.flock(lock._file.descriptor, operation)
                    lock.validate_current()

    def test_flock_record_rejects_splices(self):
        with owner.CanonicalLock(_path=str(self.root / "lock")) as lock:
            descriptor = lock._file.descriptor
            status = os.fstat(descriptor)
            text = owner.host.LinuxReader().read_text(f"/proc/{os.getpid()}/fdinfo/{descriptor}")
            line = next(line for line in text.splitlines() if line.startswith("lock:"))
            for changed in ("", line + "\n" + line, line.replace("FLOCK", "POSIX"),
                            line.replace("WRITE", "READ"), line.replace("ADVISORY", "MANDATORY"),
                            line.replace(str(os.getpid()) + " ", "-1 "),
                            line.replace(":" + str(status.st_ino) + " ", ":1 ")):
                with self.subTest(changed=changed), self.assertRaises(owner.host.PreflightError):
                    owner.validate_flock_record(changed, status, os.getpid())

    def test_actual_lock_write_and_inherited_descriptor_are_rejected(self):
        for mutation in ("write", "inheritable"):
            with self.subTest(mutation=mutation), self.assertRaises(FAILURES):
                with owner.CanonicalLock(_path=str(self.root / "lock")) as lock:
                    if mutation == "write":
                        with lock.path.open("r+b") as stream:
                            stream.write(b"changed")
                    else:
                        os.set_inheritable(lock._file.descriptor, True)
                    lock.validate_current()

    def test_lock_path_aba_is_rejected(self):
        path = self.root / "lock"
        with self.assertRaises(FAILURES):
            with owner.CanonicalLock(_path=str(path)) as lock:
                path.rename(self.root / "retained-lock")
                (self.root / "retained-lock").rename(path)
                lock.validate_current()

    def test_unsafe_lock_file_rejected_without_truncation(self):
        path = self.root / "lock"
        path.write_bytes(b"retained diagnostic")
        path.chmod(0o666)
        with self.assertRaises(FAILURES):
            with owner.CanonicalLock(_path=str(path)):
                self.fail("unsafe lock acquired")
        self.assertEqual(path.read_bytes(), b"retained diagnostic")
        path.chmod(0o600)
        link = self.root / "lock-link"
        link.symlink_to(path)
        with self.assertRaises(FAILURES):
            with owner.CanonicalLock(_path=str(link)):
                self.fail("symlink lock acquired")

    def test_bad_preregistration_or_wrong_host_precedes_lock_creation(self):
        with self.assertRaises(owner.host.PreflightError):
            owner.BuildArtifactLease(self.preregistration + b"\n")
        self.reader.process_value["hostname"] = "foureyes"
        with self.assertRaises(FAILURES):
            with self.lease():
                self.fail("wrong host entered")
        self.assertFalse((self.root / "canonical.lock").exists())
        self.assertEqual(list(self.lane.iterdir()), [])

    def test_host_drift_is_rechecked_before_freeze(self):
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                self.reader.process_value["namespaces"]["cgroup"]["inode"] += 1
                lease.freeze(self.lane, self.roots)
        self.assertEqual(list(self.lane.iterdir()), [])

    def test_real_frozen_files_and_sealed_descriptors(self):
        with self.lease() as lease:
            record = lease.freeze(self.lane, self.roots)
            self.assertEqual(set(record["artifacts"]), set(owner.ARTIFACTS))
            self.assertFalse(record["live_acquisition_armed"])
            self.assertFalse(record["source_build_history_proved"])
            self.assertFalse(record["runtime_closure_proved"])
            self.assertFalse(record["resource_lifetime_proved"])
            self.assertEqual(stat.S_IMODE((self.lane / "v19-artifacts").stat().st_mode), 0o500)
            for role in ("candidate", "baseline"):
                descriptor = lease.executable_descriptor(role)
                seals = fcntl.fcntl(descriptor, getattr(fcntl, "F_GET_SEALS", 1034))
                self.assertEqual(seals & 15, 15)
                self.assertFalse(os.get_inheritable(descriptor))
                with self.assertRaises(PermissionError):
                    os.pwrite(descriptor, b"mutate", 0)
                frozen = record["artifacts"][role + "_executable"]
                self.assertEqual(hashlib.sha256(os.pread(descriptor, frozen["size"], 0)).hexdigest(),
                                 frozen["sha256"])
                original = self.roots[role] / owner.ARTIFACTS[role + "_executable"][1]
                self.assertNotEqual(original.stat().st_ino, Path(frozen["path"]).stat().st_ino)
                # Once frozen, subsequent build output edits cannot replace
                # the lane-owned executable or immutable descriptor.
                original.write_bytes(b"a later unrelated build")
            lease.validate_current()
            record["artifacts"].clear()
            self.assertEqual(len(lease.record()["artifacts"]), 4)
        with self.assertRaises(owner.host.PreflightError):
            lease.executable_descriptor("candidate")

    def test_wrong_hash_fails_before_any_artifact_creation(self):
        (self.roots["candidate"] / "bench_leopard2").write_bytes(b"wrong build")
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                lease.freeze(self.lane, self.roots)
        self.assertEqual(list(self.lane.iterdir()), [])

    def test_unsafe_input_rejected(self):
        path = self.roots["candidate"] / "bench_leopard2"
        for kind in ("symlink", "hardlink", "mode", "fifo", "oversized"):
            retained = path.read_bytes()
            path.unlink()
            extra = self.root / "extra"
            try:
                if kind in ("symlink", "hardlink"):
                    extra.write_bytes(retained)
                    if kind == "symlink": path.symlink_to(extra)
                    else: os.link(extra, path)
                elif kind == "fifo": os.mkfifo(path)
                else:
                    path.write_bytes(retained)
                    if kind == "mode": path.chmod(0o666)
                    if kind == "oversized":
                        with path.open("r+b") as stream:
                            stream.truncate(owner.MAX_ARTIFACT_BYTES + 1)
                with self.subTest(kind=kind), self.assertRaises(FAILURES):
                    with self.lease() as lease:
                        lease.freeze(self.lane, self.roots)
                self.assertEqual(list(self.lane.iterdir()), [])
            finally:
                path.unlink()
                if extra.exists(): extra.unlink()
                path.write_bytes(retained)
                path.chmod(0o644)

    def test_build_roots_cannot_alias_lane_or_each_other(self):
        for roots in ({"candidate": self.lane, "baseline": self.roots["baseline"]},
                      {"candidate": self.roots["candidate"], "baseline": self.roots["candidate"]}):
            with self.subTest(roots=roots), self.lease() as lease:
                with self.assertRaises(FAILURES):
                    lease.freeze(self.lane, roots)
        alias = self.root / "build-alias"
        alias.symlink_to(self.roots["candidate"], target_is_directory=True)
        with self.lease() as lease, self.assertRaises(FAILURES):
            lease.freeze(self.lane, dict(self.roots, candidate=alias))

    def test_frozen_content_write_restore_is_rejected(self):
        reached_frozen = False
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                lease.freeze(self.lane, self.roots)
                reached_frozen = True
                path = self.lane / "v19-artifacts/candidate_executable"
                original = path.read_bytes()
                path.chmod(0o700)
                path.write_bytes(original)
                path.chmod(0o500)
                lease.validate_current()
        self.assertTrue(reached_frozen)

    def test_preexisting_writable_mapping_cannot_bypass_fresh_hash(self):
        target = self.lane / "v19-artifacts/candidate_executable"
        real_snapshot = owner.provenance._RetainedFileSnapshot
        mappings = []
        def map_before_capture(path, *arguments, **keywords):
            if Path(path) == target:
                target.chmod(0o700)
                try:
                    with target.open("r+b") as stream:
                        mapping = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_WRITE)
                        mapping[0] = mapping[0]  # Fault writable page before baseline metadata.
                        mappings.append(mapping)
                finally:
                    target.chmod(0o500)
            return real_snapshot(path, *arguments, **keywords)
        mutated = False
        try:
            with self.assertRaises(FAILURES):
                with self.lease() as lease:
                    with mock.patch.object(owner.provenance, "_RetainedFileSnapshot",
                                           side_effect=map_before_capture):
                        lease.freeze(self.lane, self.roots)
                    mappings[0][4] ^= 1
                    mutated = True
                    lease.validate_current()
            self.assertTrue(mutated)
        finally:
            for mapping in mappings:
                mapping.close()

    def test_frozen_directory_aba_is_rejected(self):
        reached_frozen = False
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                lease.freeze(self.lane, self.roots)
                reached_frozen = True
                root = self.lane / "v19-artifacts"
                root.rename(self.lane / "temporary-name")
                (self.lane / "temporary-name").rename(root)
                lease.validate_current()
        self.assertTrue(reached_frozen)

    def test_preexisting_output_is_not_reused_or_removed(self):
        destination = self.lane / "v19-artifacts"
        destination.mkdir(mode=0o700)
        (destination / "user-file").write_bytes(b"preserve")
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                lease.freeze(self.lane, self.roots)
        self.assertEqual((destination / "user-file").read_bytes(), b"preserve")

    def test_partial_write_failure_poisoned_and_descriptors_closed(self):
        before = set(os.listdir(f"/proc/{os.getpid()}/fd"))
        lease = self.lease()
        with self.assertRaises(owner.host.PreflightError):
            with lease:
                with mock.patch.object(owner.os, "write", side_effect=OSError("injected disk failure")):
                    with self.assertRaisesRegex(OSError, "injected disk failure"):
                        lease.freeze(self.lane, self.roots)
                with self.assertRaises(owner.host.PreflightError): lease.record()
                with self.assertRaises(owner.host.PreflightError): lease.freeze(self.lane, self.roots)
                # Suppression cannot turn a failed lease into successful exit.
        self.assertEqual(set(os.listdir(f"/proc/{os.getpid()}/fd")), before)
        self.assertTrue((self.lane / "v19-artifacts/candidate_executable").exists())
        with owner.CanonicalLock(_path=str(self.root / "canonical.lock")):
            pass

    def test_terminal_host_change_is_rejected(self):
        reached_frozen = False
        with self.assertRaises(FAILURES):
            with self.lease() as lease:
                lease.freeze(self.lane, self.roots)
                reached_frozen = True
                self.reader.files[fixtures.SCOPE + "/memory.swap.max"] = "1\n"
        self.assertTrue(reached_frozen)

    def test_early_execution_and_duplicate_freeze_are_rejected(self):
        with self.lease() as lease:
            with self.assertRaises(owner.host.PreflightError): lease.executable_descriptor("candidate")
            lease.freeze(self.lane, self.roots)
            with self.assertRaises(owner.host.PreflightError): lease.freeze(self.lane, self.roots)
            with self.assertRaises(owner.host.PreflightError): lease.executable_descriptor("foreign")


if __name__ == "__main__":
    unittest.main()
