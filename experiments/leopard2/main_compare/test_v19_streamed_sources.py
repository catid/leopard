#!/usr/bin/python3
"""Offline lifetime/buffer tests for leopard-79h.38.5.4.8.2.2.2.1."""
import hashlib
import importlib.util
import mmap
import os
from pathlib import Path
import tempfile
import tracemalloc
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_streamed_sources", HERE / "v19_streamed_sources.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError)


class StreamedSourceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-stream-test-")
        self.addCleanup(temporary.cleanup)
        self.parent = Path(temporary.name)
        self.root = self.parent / "workspace"
        self.root.mkdir(mode=0o700)
        for role in module.SOURCE_ROOTS:
            source = self.root / role
            (source / ".git").mkdir(mode=0o700, parents=True)
            (source / ".git/HEAD").write_text("a" * 40 + "\n")
            (source / "source.cpp").write_bytes(b"source payload\n")
        self.source = self.root / "candidate-source/source.cpp"

    def owner(self):
        return module.StreamingSourceOwner(self.root)

    def test_valid_owner_binds_inventory_without_claiming_git_authenticity(self):
        with self.owner() as owner:
            record = owner.record()
            self.assertEqual(record["file_count"], 4)
            self.assertEqual(record["directory_count"], 4)
            self.assertTrue(record["staged_file_ownership_verified"])
            for key in ("git_identity_verified", "source_creation_verified", "runtime_closure_verified",
                        "atomic_snapshot", "live_acquisition_armed", "cache_eviction_requested"):
                self.assertIs(record[key], False)
            expected = next(row for row in record["files"] if row["path"].endswith("source.cpp"))
            self.assertEqual(expected["sha256"], hashlib.sha256(b"source payload\n").hexdigest())
            record["files"].clear()
            self.assertEqual(owner.record()["file_count"], 4)
            (self.root / "candidate-build").mkdir()
            (self.root / "candidate-build/output").write_bytes(b"unrelated build output")
            owner.validate_current()
        with self.assertRaises(module.host.PreflightError): owner.record()
        with self.assertRaises(module.host.PreflightError): owner.__enter__()

    def test_source_permissions_inside_private_boundary_are_preserved(self):
        self.source.chmod(0o664)
        self.source.parent.chmod(0o775)
        before = self.source.stat()
        with self.owner() as owner:
            owner.validate_current()
        self.assertEqual(module.provenance._stable_fields(before),
                         module.provenance._stable_fields(self.source.stat()))

    def test_workspace_must_be_private_and_not_a_symlink(self):
        self.root.chmod(0o775)
        with self.assertRaises(FAILURES):
            with self.owner(): pass
        self.root.chmod(0o700)
        link = self.parent / "alias"
        link.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(FAILURES):
            with module.StreamingSourceOwner(link): pass

    def test_linked_or_missing_git_metadata_is_rejected(self):
        directory = self.root / "candidate-source/.git"
        saved = self.parent / "saved-git"
        directory.rename(saved)
        for kind in ("missing", "gitfile", "symlink"):
            if kind == "gitfile": directory.write_text("gitdir: " + str(saved) + "\n")
            if kind == "symlink": directory.symlink_to(saved, target_is_directory=True)
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.owner(): pass
            if kind != "missing": directory.unlink()

    def test_unsafe_file_types_modes_and_links_are_rejected(self):
        other = self.parent / "foreign"
        other.write_bytes(b"source payload\n")
        for kind in ("hardlink", "symlink", "fifo", "world-writable", "setuid"):
            self.source.unlink()
            if kind == "hardlink": os.link(other, self.source)
            elif kind == "symlink": self.source.symlink_to(other)
            elif kind == "fifo": os.mkfifo(self.source)
            else:
                self.source.write_bytes(b"source payload\n")
                self.source.chmod(0o666 if kind == "world-writable" else 0o4644)
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.owner(): pass

    def test_inventory_and_byte_limits_precede_unbounded_reads(self):
        for key, value in (("MAX_FILES", 1), ("MAX_DIRECTORIES", 1), ("MAX_DEPTH", 1),
                           ("MAX_TOTAL_BYTES", 1), ("MAX_FILE_BYTES", 1), ("MAX_PATH_BYTES", 1)):
            with self.subTest(key=key), mock.patch.object(module, key, value), self.assertRaises(FAILURES):
                with self.owner(): pass

    def test_path_and_metadata_changes_are_rejected_while_held(self):
        for kind in ("rename-aba", "write", "metadata-aba", "extra", "delete"):
            entered = False
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    entered = True
                    if kind == "rename-aba":
                        self.source.rename(self.parent / "saved")
                        (self.parent / "saved").rename(self.source)
                    elif kind == "write": self.source.write_bytes(b"changed payload")
                    elif kind == "metadata-aba":
                        mode = self.source.stat().st_mode & 0o777
                        self.source.chmod(0o400)
                        self.source.chmod(mode)
                    elif kind == "extra": (self.source.parent / "extra").write_bytes(b"x")
                    else: self.source.unlink()
                    owner.validate_current()
            self.assertTrue(entered)
            extra = self.source.parent / "extra"
            if extra.exists(): extra.unlink()
            self.source.write_bytes(b"source payload\n")

    def test_workspace_and_subdirectory_aba_are_rejected(self):
        for path in (self.root, self.root / "candidate-source/.git"):
            entered = False
            with self.subTest(path=path), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    entered = True
                    saved = self.parent / "temporary-directory"
                    path.rename(saved)
                    saved.rename(path)
                    owner.validate_current()
            self.assertTrue(entered)

    def test_preexisting_writable_mapping_is_rehashed_and_failure_latches(self):
        with self.source.open("r+b") as stream:
            mapping = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_WRITE)
        mapping[0] = mapping[0]
        try:
            with self.assertRaisesRegex(module.host.PreflightError, "failed source owner"):
                with self.owner() as owner:
                    original = mapping[0]
                    mapping[0] ^= 1
                    with self.assertRaisesRegex(module.host.PreflightError, "current source bytes"):
                        owner.validate_current()
                    mapping[0] = original
                    with self.assertRaisesRegex(module.host.PreflightError, "not usable"):
                        owner.record()
        finally:
            mapping.close()

    def test_descriptor_inheritance_and_foreign_process_are_rejected(self):
        for kind in ("inherit", "pid"):
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    if kind == "inherit": os.set_inheritable(next(iter(owner._files.values()))[0], True)
                    else: owner._pid = -1
                    owner.validate_current()

    def test_eviction_uses_only_exact_held_file_descriptors(self):
        with self.owner() as owner:
            expected = {entry[0] for entry in owner._files.values()}
            with mock.patch.object(module.os, "fsync") as sync, \
                 mock.patch.object(module.os, "posix_fadvise") as advise:
                record = owner.record(evict_cache=True)
            self.assertTrue(record["cache_eviction_requested"])
            self.assertEqual({call.args[0] for call in sync.call_args_list}, expected)
            self.assertEqual({call.args for call in advise.call_args_list},
                             {(fd, 0, 0, os.POSIX_FADV_DONTNEED) for fd in expected})

    def test_real_cache_eviction_preserves_content_and_metadata(self):
        with self.owner() as owner:
            before = owner.record()
            owner.validate_current(evict_cache=True)
            self.assertEqual(owner.record(), before)

    def test_changed_git_metadata_cannot_reach_cache_eviction(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                (self.root / "candidate-source/.git/HEAD").write_text("b" * 40 + "\n")
                with mock.patch.object(module.os, "posix_fadvise") as advise:
                    try:
                        owner.validate_current(evict_cache=True)
                    finally:
                        advise.assert_not_called()

    def test_eviction_failure_closes_owner_and_never_uses_global_fallback(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                with mock.patch.object(module.os, "posix_fadvise", side_effect=OSError("unsupported")):
                    owner.validate_current(evict_cache=True)
        with self.assertRaises(module.host.PreflightError): owner.record()

    def test_failed_capture_and_exit_close_all_owned_descriptors(self):
        before = set(os.listdir("/proc/self/fd"))
        self.source.chmod(0o666)
        with self.assertRaises(FAILURES):
            with self.owner(): pass
        self.assertEqual(set(os.listdir("/proc/self/fd")), before)
        self.source.chmod(0o644)
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                self.source.write_bytes(b"changed")
                owner.validate_current()
        self.assertEqual(set(os.listdir("/proc/self/fd")), before)

    def test_large_file_uses_bounded_buffers_not_retained_content(self):
        large = self.root / "candidate-source/large.dat"
        with large.open("wb") as stream:
            stream.truncate(32 << 20)
        read_sizes = []
        pread = os.pread
        def bounded_read(descriptor, count, offset):
            read_sizes.append(count)
            return pread(descriptor, count, offset)
        tracemalloc.start()
        try:
            # A recording Mock would retain every call and inflate the measured
            # peak itself. Only retain the requested read sizes, not call objects.
            with mock.patch.object(module.os, "pread", new=bounded_read):
                with self.owner() as owner:
                    self.assertGreater(owner.record()["total_bytes"], 32 << 20)
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertLess(peak, 3 << 20)
        self.assertLessEqual(max(read_sizes), module.HASH_BLOCK_BYTES)


if __name__ == "__main__":
    unittest.main()
