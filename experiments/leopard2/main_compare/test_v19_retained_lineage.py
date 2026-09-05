#!/usr/bin/python3
"""Real-fd synthetic archive tests; leopard-79h.38.5.4.8.2.2.2.2."""
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import mmap
import os
from pathlib import Path
import stat
import subprocess
import tempfile
import tracemalloc
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_retained_lineage", HERE / "v19_retained_lineage.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


def encoded(value):
    return module.host.canonical_bytes(value) + b"\n"


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(65536), b""): result.update(part)
    return result.hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): path.chmod(0o600)
    path.write_bytes(data if type(data) is bytes else encoded(data))
    path.chmod(0o400)


def manifest(root, *, core=False):
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != root / "SHA256SUMS" and (not core or path.name != "SHA256SUMS"):
            rows.append((path.relative_to(root).as_posix(), digest(path)))
    return "".join(sha + "  ./" + name + "\n" for name, sha in sorted(rows)).encode()


class LineageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preregistration = subprocess.run([str(HERE / "run_authoritative_v17_gfni_main_compare.sh"),
            "--print-conditioned-v19-preregistration"], check=True, stdout=subprocess.PIPE).stdout
        cls.real_contract = module.host.load_preregistration(cls.preregistration)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-lineage-test-")
        self.addCleanup(temporary.cleanup)
        self.parent = Path(temporary.name)
        self.addCleanup(self.make_writable)
        self.archives = []
        for attempt, name, _sha in module.ARCHIVES:
            root = self.parent / name
            root.mkdir()
            common = {"status": "failed", "acquisition_generation": "passive-v2", "attempt": attempt,
                      "attempt_budget": 3, "source_commit": module.SOURCE_COMMIT, "source_tree": module.SOURCE_TREE,
                      "promotion_passed": False, "campaign_exit_status": 1, "failure_verified": True}
            priors = [{"attempt": number, "envelope": "/home/catid/leopard/.research/leopard-79h/" + prior_name,
                       "envelope_sha256sums_sha256": sha, "terminal": "FAILED.json",
                       "terminal_schema": "leopard2-v18-gfni-main-failed-envelope/v1"}
                      for number, prior_name, sha in self.archives]
            write(root / "core/attempt-lineage.json", {"schema": "leopard2-v18-gfni-main-attempt-lineage/v1",
                "acquisition_generation": "passive-v2", "attempt": attempt, "attempt_budget": 3,
                "source_commit": module.SOURCE_COMMIT, "source_tree": module.SOURCE_TREE, "prior_attempts": priors})
            common["attempt_lineage_sha256"] = digest(root / "core/attempt-lineage.json")
            write(root / "core/campaign/failure.json", {"schema": "leopard2-main-compare-failure/v18",
                "status": "failed", "valid": False, "error_type": "EvidenceError"})
            write(root / "core/manifest.json", {**common, "schema": "leopard2-v18-gfni-main-passive-failed-core-manifest/v1",
                "failure_verify_status": 0, "baseline_commit": self.real_contract["build_preflight"]["baseline_commit"],
                "failure_sha256": digest(root / "core/campaign/failure.json"), "canonical_lock": module.owners.LOCK_PATH,
                "cpu": 52, "sibling": 116})
            write(root / "core/nested/payload", b"payload")
            write(root / "core/nested/SHA256SUMS", b"historical nested manifest is outer-bound\n")
            write(root / "core/SHA256SUMS", manifest(root / "core", core=True))
            write(root / "FAILED.json", {**common, "schema": "leopard2-v18-gfni-main-failed-envelope/v1",
                "core_sha256sums_sha256": digest(root / "core/SHA256SUMS")})
            self.seal(root)
            self.archives.append((attempt, name, digest(root / "SHA256SUMS")))
        patcher = mock.patch.object(module, "ARCHIVES", self.archives)
        patcher.start()
        self.addCleanup(patcher.stop)
        original = module.host.load_preregistration
        def fixture_contract(data):
            result = original(data)
            result["attempt_contract"]["failure_lineage_sha256"] = hashlib.sha256(encoded(module.lineage_record())).hexdigest()
            return result
        patcher = mock.patch.object(module.host, "load_preregistration", new=fixture_contract)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.last = self.parent / self.archives[-1][1]

    def make_writable(self):
        for path in self.parent.rglob("*"):
            if not path.is_symlink(): path.chmod(0o700 if path.is_dir() else 0o600)

    def seal(self, root):
        root.chmod(0o700)
        write(root / "SHA256SUMS", b"")
        entries = []
        for path in [root, *sorted(root.rglob("*"))]:
            if path == root / "TREE-METADATA.json": continue
            value = path.stat()
            entries.append({"path": "." if path == root else path.relative_to(root).as_posix(),
                "type": "directory" if path.is_dir() else "file", "mode": format(stat.S_IMODE(value.st_mode) & ~0o222, "04o"),
                "nlink": value.st_nlink, "uid": value.st_uid, "gid": value.st_gid})
        write(root / "TREE-METADATA.json", {"entries": sorted(entries, key=lambda row: row["path"]),
            "schema": "leopard2-authoritative-tree-metadata/v1", "root": ".", "excluded_paths": ["TREE-METADATA.json"],
            "final_mode_policy": "observed mode with all write bits removed",
            "uid_gid_policy": {"uid": os.geteuid(), "gid": os.getegid(),
                               "rule": "every retained node has the invoking effective uid and gid"},
            "self_policy": {"uid": os.geteuid(), "gid": os.getegid(), "mode": "0400", "nlink": 1,
                            "type": "file", "sha256_binding": "exactly one ./TREE-METADATA.json checksum entry"}})
        write(root / "SHA256SUMS", manifest(root))
        for path in [root, *root.rglob("*")]: path.chmod(stat.S_IMODE(path.stat().st_mode) & ~0o222)

    def repin_last(self):
        root = self.last
        for path in root.rglob("*"):
            if path.is_dir(): path.chmod(0o700)
        root.chmod(0o700)
        write(root / "core/SHA256SUMS", manifest(root / "core", core=True))
        terminal = json.loads((root / "FAILED.json").read_bytes())
        terminal["core_sha256sums_sha256"] = digest(root / "core/SHA256SUMS")
        write(root / "FAILED.json", terminal)
        self.seal(root)
        attempt, name, _sha = self.archives[-1]
        self.archives[-1] = attempt, name, digest(root / "SHA256SUMS")

    @contextmanager
    def held(self):
        with module.PinnedV18Lineage(self.preregistration, self.parent) as owner: yield owner

    def test_all_archives_held_and_record_has_no_replay_or_runtime_authority(self):
        with self.held() as owner:
            record = owner.record()
            self.assertEqual(record["file_count"], 27)
            self.assertEqual(record["directory_count"], 12)
            self.assertTrue(record["physical_archives_verified"])
            for key in ("historical_failures_replayed", "runtime_closure_verified", "atomic_snapshot", "live_acquisition_armed"):
                self.assertIs(record[key], False)
            record["lineage"]["attempts"].clear()
            self.assertEqual(len(owner.record()["lineage"]["attempts"]), 3)
            self.assertTrue(all(not os.get_inheritable(entry[0]) for entry in owner._files.values()))
        with self.assertRaises(FAILURES): owner.record()
        with self.assertRaises(FAILURES): owner.__enter__()

    def test_original_contract_digest_matches_frozen_v19_authority(self):
        expected = [(1, "c8f825d-v18-passive-main-a1", "ce65c3a49ef1c1d89ba51ea03d0af4742d6790e6f2ea2662917d9ef9a9d945d7"),
                    (2, "c8f825d-v18-passive-main-a2", "a1bf0eda157c251f33f7260ebd76931d88054d460bd07a97bcba2811384b2c10"),
                    (3, "c8f825d-v18-passive-main-a3", "fe5b40cc98753cbd794ee019cb0e2643d0ccee0aca4c5fd7b2e0b27df8a86139")]
        with mock.patch.object(module, "ARCHIVES", expected):
            self.assertEqual(hashlib.sha256(encoded(module.lineage_record())).hexdigest(),
                             self.real_contract["attempt_contract"]["failure_lineage_sha256"])

    def test_parent_permissions_and_symlinks_are_rejected(self):
        self.parent.chmod(0o775)
        with self.assertRaises(FAILURES), self.held(): pass
        self.parent.chmod(0o700)
        root = self.parent / self.archives[0][1]
        target = self.parent / "moved"
        root.rename(target)
        root.symlink_to(target, target_is_directory=True)
        with self.assertRaises(FAILURES), self.held(): pass

    def test_wrong_file_bytes_are_rejected(self):
        write(self.last / "core/nested/payload", b"changed")
        with self.assertRaisesRegex(module.host.PreflightError, "inventory or hashes"), self.held(): pass

    def test_missing_and_extra_files_are_rejected(self):
        path = self.last / "core/nested/payload"
        path.parent.chmod(0o700)
        path.unlink()
        path.parent.chmod(0o500)
        with self.assertRaises(FAILURES), self.held(): pass
        path.parent.chmod(0o700)
        write(path, b"payload")
        path.parent.chmod(0o500)
        self.last.chmod(0o700)
        write(self.last / "extra", b"extra")
        self.last.chmod(0o500)
        with self.assertRaises(FAILURES), self.held(): pass

    def test_hardlinked_file_rejected(self):
        path = self.last / "core/nested/payload"
        os.link(path, self.parent / "alias")
        with self.assertRaisesRegex(module.host.PreflightError, "single-link"), self.held(): pass

    def test_missing_required_manifest_rejected_cleanly(self):
        self.last.chmod(0o700)
        (self.last / "SHA256SUMS").unlink()
        self.last.chmod(0o500)
        with self.assertRaisesRegex(module.host.PreflightError, "required evidence"), self.held(): pass

    def test_unlisted_empty_directory_rejected(self):
        self.last.chmod(0o700)
        extra = self.last / "extra-directory"
        extra.mkdir(mode=0o500)
        self.last.chmod(0o500)
        with self.assertRaisesRegex(module.host.PreflightError, "tree metadata"), self.held(): pass

    def test_symlink_leaf_rejected(self):
        path = self.last / "core/nested/payload"
        path.parent.chmod(0o700)
        path.unlink()
        path.symlink_to("/etc/hostname")
        path.parent.chmod(0o500)
        with self.assertRaises(FAILURES), self.held(): pass

    def test_writable_file_and_wrong_sealed_mode_rejected(self):
        path = self.last / "core/nested/payload"
        path.chmod(0o600)
        with self.assertRaises(FAILURES), self.held(): pass
        path.chmod(0o444)
        with self.assertRaisesRegex(module.host.PreflightError, "tree metadata"), self.held(): pass

    def test_write_restore_history_rejected(self):
        with self.assertRaises(FAILURES), self.held() as owner:
            path = self.last / "core/nested/payload"
            write(path, b"changed")
            write(path, b"payload")
            owner.validate_current()

    def test_persistent_mmap_drift_rejected_and_latched(self):
        path = self.last / "core/nested/payload"
        path.chmod(0o600)
        with path.open("r+b") as file, mmap.mmap(file.fileno(), 0) as view:
            # Fault the writable mapping before capture. Later stores can then
            # leave mtime/ctime and inotify unchanged: only a byte read catches it.
            view[0] = view[0]
            path.chmod(0o400)
            with self.assertRaisesRegex(module.host.PreflightError, "failed lineage owner"):
                with self.held() as owner:
                    before = module.provenance._stable_fields(path.stat())
                    view[0] = ord("X")
                    self.assertEqual(module.provenance._stable_fields(path.stat()), before)
                    owner._guard.verify()
                    with self.assertRaisesRegex(module.host.PreflightError, "archive bytes"):
                        owner.validate_current()
                    view[0] = ord("p")
                    with self.assertRaisesRegex(module.host.PreflightError, "not live"):
                        owner.record()

    def test_directory_rename_restore_rejected(self):
        with self.assertRaises(FAILURES), self.held() as owner:
            name = self.parent / self.archives[-1][1]
            moved = self.parent / "moved"
            name.rename(moved)
            moved.rename(name)
            owner.validate_current()

    def test_parent_permission_restore_rejected(self):
        with self.assertRaises(FAILURES), self.held() as owner:
            self.parent.chmod(0o777)
            self.parent.chmod(0o700)
            owner.validate_current()

    def test_lost_guard_rejected_and_latched(self):
        with self.assertRaises(FAILURES), self.held() as owner:
            owner._guard._close_without_verification()
            with self.assertRaises(FAILURES): owner.validate_current()
            self.assertEqual(owner._state, "failed")

    def test_inheritable_file_descriptor_rejected(self):
        with self.assertRaises(FAILURES), self.held() as owner:
            os.set_inheritable(next(iter(owner._files.values()))[0], True)
            owner.validate_current()

    def test_resigned_terminal_cannot_claim_promotion_or_success(self):
        terminal = json.loads((self.last / "FAILED.json").read_bytes())
        terminal["promotion_passed"] = True
        write(self.last / "FAILED.json", terminal)
        self.repin_last()
        with self.assertRaisesRegex(module.host.PreflightError, "failure terminal"), self.held(): pass

    def test_resigned_prior_attempt_cannot_change_label_or_digest(self):
        path = self.last / "core/attempt-lineage.json"
        value = json.loads(path.read_bytes())
        value["prior_attempts"][0]["envelope"] = str(self.parent / self.archives[0][1])
        write(path, value)
        new_digest = digest(path)
        for relative in ("FAILED.json", "core/manifest.json"):
            record = json.loads((self.last / relative).read_bytes())
            record["attempt_lineage_sha256"] = new_digest
            write(self.last / relative, record)
        self.repin_last()
        with self.assertRaisesRegex(module.host.PreflightError, "prior-attempt"), self.held(): pass

    def test_core_excludes_nested_checksums_but_outer_still_binds_them(self):
        with self.held() as owner: owner.validate_current()
        write(self.last / "core/nested/SHA256SUMS", b"changed")
        with self.assertRaisesRegex(module.host.PreflightError, "inventory or hashes"), self.held(): pass

    def test_manifest_parser_rejects_duplicates_unsorted_and_unsafe_names(self):
        digest_text = "a" * 64
        for data in (f"{digest_text}  ./x\n{digest_text}  ./x\n", f"{digest_text}  ./z\n{digest_text}  ./a\n",
                     f"{digest_text}  ./../escape\n", f"{digest_text}  ./SHA256SUMS\n", f"{digest_text}  ./x"):
            with self.subTest(data=data), self.assertRaises(FAILURES): module.checksums(data.encode())

    def test_metadata_duplicate_key_rejected(self):
        with self.assertRaises(FAILURES): module._json(b'{"schema":1,"schema":2}')

    def test_bounds_and_failed_constructor_release_descriptors(self):
        before = set(os.listdir("/proc/self/fd"))
        for key, value in (("MAX_FILES", 2), ("MAX_DIRECTORIES", 1), ("MAX_TOTAL_BYTES", 3),
                           ("MAX_FILE_BYTES", 2), ("MAX_JSON_BYTES", 8)):
            with self.subTest(key=key), mock.patch.object(module, key, value), self.assertRaises(FAILURES), self.held(): pass
            self.assertEqual(set(os.listdir("/proc/self/fd")), before)

    def test_large_archive_is_streamed_without_retaining_its_body(self):
        self.last.chmod(0o700)
        path = self.last / "large-archive.tar"
        with path.open("wb") as file: file.truncate(85 << 20)
        path.chmod(0o400)
        self.repin_last()
        maximum = 0
        pread = os.pread
        def bounded(descriptor, count, offset):
            nonlocal maximum
            if os.fstat(descriptor).st_size == 85 << 20:
                maximum = max(maximum, count)
            return pread(descriptor, count, offset)
        tracemalloc.start()
        try:
            with mock.patch.object(module.os, "pread", new=bounded), self.held() as owner:
                self.assertGreater(owner.record()["total_bytes"], 85 << 20)
                self.assertLess(tracemalloc.get_traced_memory()[1], 4 << 20)
        finally: tracemalloc.stop()
        self.assertEqual(maximum, 64 << 10)


if __name__ == "__main__":
    unittest.main()
