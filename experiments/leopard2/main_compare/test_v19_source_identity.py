#!/usr/bin/python3
"""Bounded real-Git source authentication tests; leopard-79h.38.5.4.8.2.2.2.1."""
from contextlib import contextmanager
import copy
import hashlib
import importlib.util
import json
import mmap
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import tracemalloc
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_source_identity", HERE / "v19_source_identity.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError)


def git(root, *arguments):
    return subprocess.run(["/usr/bin/git", "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgSign=false",
        "-c", "user.name=Source Fixture", "-c", "user.email=fixture@example.invalid", "-C", str(root),
        *arguments], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=dict(module.provenance.GIT_ENVIRONMENT), timeout=15).stdout


class SourceIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-auth-template-")
        cls.addClassCleanup(temporary.cleanup)
        cls.template = Path(temporary.name)
        for role in module.REPOSITORIES:
            root = cls.template / role
            root.mkdir(parents=True)
            git(root, "init", "--initial-branch=master")
            (root / "source.cpp").write_bytes((role + " source\n").encode())
            (root / "foo").mkdir()
            (root / "foo/item").write_bytes(b"nested")
            (root / "foo.bar").write_bytes(b"file before directory in Git ordering")
            (root / "run.sh").write_bytes(b"#!/bin/sh\nexit 0\n")
            (root / "run.sh").chmod(0o755)
        for role in ("candidate-source/sse2neon", "candidate-source", "leopard1-source"):
            root = cls.template / role
            git(root, "add", "--all")
            if role == "leopard1-source":
                submodule = git(cls.template / "candidate-source/sse2neon", "rev-parse", "HEAD").decode().strip()
                (root / "sse2neon").mkdir()
                git(root, "update-index", "--add", "--cacheinfo", "160000," + submodule + ",sse2neon")
            git(root, "commit", "-m", "fixture")
            git(root, "checkout", "--detach")
        cls.identities = {role: git(cls.template / role, "rev-parse", "HEAD", "HEAD^{tree}").decode().splitlines()
                          for role in module.REPOSITORIES}

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-auth-test-")
        self.addCleanup(temporary.cleanup)
        self.parent = Path(temporary.name)
        self.root = self.parent / "workspace"
        shutil.copytree(self.template, self.root)
        self.root.chmod(0o700)
        self.canonical = str(self.parent / "historical")
        for role in module.REPOSITORIES:
            (self.root / role / ".git/config").write_bytes(module.clone_config(role, self.canonical))
        self.source = self.root / "candidate-source/source.cpp"
        self.manifest = self.current_manifest()
        self.pins = {"source_commit": self.identities["candidate-source"][0],
                     "source_tree": self.identities["candidate-source"][1],
                     "baseline_commit": self.identities["leopard1-source"][0],
                     "baseline_tree": self.identities["leopard1-source"][1]}
        # Only the retained-preflight input is synthetic. Source fds, Git object
        # stores, commit/index queries, sealed execution and mutations are real.
        self.retained = object.__new__(module.preflight.PinnedPreflight)
        self.preflight_live = True
        def validate():
            module.require(self.preflight_live, "fixture preflight is not held")
        self.retained.validate_current = validate
        self.retained.record = lambda: {"source_pins": copy.deepcopy(self.pins), "canonical_run_root": self.canonical}
        self.retained._bytes = lambda name: json.dumps({"tracked_source_manifest": self.manifest}).encode()
        patcher = mock.patch.object(module, "SUBMODULE_COMMIT", self.identities["candidate-source/sse2neon"][0])
        patcher.start()
        self.addCleanup(patcher.stop)

    def current_manifest(self):
        files = []
        root = self.root / "candidate-source"
        for path in sorted(root.rglob("*")):
            relative = path.relative_to(root).as_posix()
            if path.is_file() and ".git" not in relative.split("/"):
                data = path.read_bytes()
                files.append({"path": relative, "mode": path.stat().st_mode & 0o777,
                              "size": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        files.sort(key=lambda row: row["path"])
        return {"schema": "leopard2-tracked-source-tree/v1", "files": files,
                "total_bytes": sum(row["size"] for row in files)}

    @contextmanager
    def authenticated(self):
        with module.streamed.StreamingSourceOwner(self.root) as source:
            with module.PinnedSourceIdentity(source, self.retained) as identity:
                yield identity

    def test_real_commit_tree_index_and_sealed_git_identity(self):
        run = module.provenance._run
        queries = []
        def checked(command, label, **keywords):
            self.assertEqual(keywords["maximum_bytes"], 1 << 20)
            self.assertEqual(keywords["timeout"], 30)
            self.assertEqual(len(keywords["inherited_descriptors"]), 2)
            self.assertGreaterEqual(keywords["executable_descriptor"], 0)
            queries.append(tuple(command))
            return run(command, label, **keywords)
        with mock.patch.object(module.provenance, "_run", new=checked), self.authenticated() as identity:
            record = identity.record()
            self.assertTrue(record["source_content_verified"])
            self.assertTrue(record["commit_and_index_verified"])
            for key in ("source_creation_verified", "runtime_closure_verified", "atomic_snapshot", "live_acquisition_armed"):
                self.assertIs(record[key], False)
            self.assertEqual(record["git_executable"]["seals"] & 15, 15)
            for role in module.REPOSITORIES:
                self.assertEqual(record["repositories"][role]["commit"], self.identities[role][0])
                self.assertEqual(record["repositories"][role]["tree"], self.identities[role][1])
            record["repositories"].clear()
            self.assertEqual(len(identity.record()["repositories"]), 3)
            identity.validate_current(evict_cache=True)
        self.assertEqual(len(queries), 12)
        with self.assertRaises(FAILURES): identity.record()
        with self.assertRaises(FAILURES): identity.__enter__()

    def test_live_handles_are_required(self):
        with self.assertRaises(FAILURES): module.PinnedSourceIdentity({}, self.retained)
        with module.streamed.StreamingSourceOwner(self.root) as source:
            with self.assertRaises(FAILURES): module.PinnedSourceIdentity(source, {})

    def test_wrong_initial_candidate_content_or_mode_is_rejected(self):
        original = self.source.read_bytes()
        mode = self.source.stat().st_mode & 0o777
        for kind in ("bytes", "mode", "extra", "missing"):
            with self.subTest(kind=kind):
                extra = self.source.parent / "ignored-untracked"
                if kind == "bytes": self.source.write_bytes(b"incorrect initial content")
                elif kind == "mode": self.source.chmod(0o755)
                elif kind == "extra": extra.write_bytes(b"untracked")
                else: self.source.unlink()
                with self.assertRaisesRegex(module.host.PreflightError, "candidate worktree differs"):
                    with self.authenticated(): pass
                if extra.exists(): extra.unlink()
                self.source.write_bytes(original)
                self.source.chmod(mode)

    def test_wrong_initial_baseline_content_is_rejected(self):
        (self.root / "leopard1-source/source.cpp").write_bytes(b"wrong baseline")
        with self.assertRaisesRegex(module.host.PreflightError, "baseline worktree differs"):
            with self.authenticated(): pass

    def test_baseline_owner_execute_mode_matches_git_not_other_execute_bits(self):
        baseline = self.root / "leopard1-source"
        (baseline / "run.sh").chmod(0o654)
        self.assertNotEqual(git(baseline, "diff-files", "--raw"), b"")
        with self.assertRaisesRegex(module.host.PreflightError, "baseline worktree differs"):
            with self.authenticated(): pass

    def test_submodule_bytes_are_authenticated_even_with_a_coherent_manifest(self):
        (self.root / "candidate-source/sse2neon/source.cpp").write_bytes(b"wrong submodule")
        self.manifest = self.current_manifest()
        with self.assertRaisesRegex(module.host.PreflightError, "raw commit names a different content tree"):
            with self.authenticated(): pass

    def test_candidate_sha256_manifest_cannot_disagree_with_git_tree(self):
        self.manifest["files"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(module.host.PreflightError, "pinned SHA-256 manifest"):
            with self.authenticated(): pass

    def test_index_cannot_hide_assume_unchanged_or_skip_worktree(self):
        for flag in ("assume-unchanged", "skip-worktree"):
            with self.subTest(flag=flag):
                git(self.source.parent, "update-index", "--" + flag, "source.cpp")
                with self.assertRaisesRegex(module.host.PreflightError, "non-default flags"):
                    with self.authenticated(): pass
                git(self.source.parent, "update-index", "--no-" + flag, "source.cpp")

    def test_index_object_must_match_actual_source_bytes(self):
        wrong = git(self.source.parent, "rev-parse", "HEAD:foo.bar").decode().strip()
        git(self.source.parent, "update-index", "--cacheinfo", "100644," + wrong + ",source.cpp")
        with self.assertRaisesRegex(module.host.PreflightError, "index differs"):
            with self.authenticated(): pass

    def test_symbolic_or_wrong_head_is_rejected_before_git_queries(self):
        head = self.source.parent / ".git/HEAD"
        for contents in (b"ref: refs/heads/master\n", b"0" * 40 + b"\n"):
            head.write_bytes(contents)
            with self.assertRaisesRegex(module.host.PreflightError, "pinned detached commit"):
                with self.authenticated(): pass

    def test_unsafe_git_configuration_and_metadata_never_reaches_queries(self):
        config = self.source.parent / ".git/config"
        original = config.read_bytes()
        for relative, content in (("config", original + b"[include]\n\tpath=/tmp/external\n"),
                ("objects/info/alternates", b"/tmp/elsewhere\n"), ("info/grafts", b"redirect\n"),
                ("hooks/post-checkout", b"#!/bin/sh\nexit 0\n"),
                ("packed-refs", b"0" * 40 + b" refs/replace/" + b"1" * 40 + b"\n")):
            path = self.source.parent / ".git" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            with self.subTest(relative=relative), \
                 mock.patch.object(module.PinnedSourceIdentity, "_query", side_effect=AssertionError("query reached")):
                with self.assertRaises(module.host.PreflightError):
                    with self.authenticated(): pass
            if relative == "config": config.write_bytes(original)
            else: path.unlink()

    def test_additional_git_root_is_rejected(self):
        (self.source.parent / "unexpected/.git").mkdir(parents=True)
        with self.assertRaisesRegex(module.host.PreflightError, "Git root inventory"):
            with self.authenticated(): pass

    def test_gitfiles_outside_known_metadata_are_not_ignored(self):
        (self.source.parent / "unexpected").mkdir()
        (self.source.parent / "unexpected/.git").write_bytes(b"gitdir: /tmp/elsewhere\n")
        with self.assertRaisesRegex(module.host.PreflightError, "source tree path"):
            with self.authenticated(): pass

    def test_untracked_empty_directories_are_rejected(self):
        (self.source.parent / "unexpected").mkdir()
        with self.assertRaisesRegex(module.host.PreflightError, "untracked empty directories"):
            with self.authenticated(): pass

    def test_baseline_uninitialized_gitlink_is_required_and_must_stay_empty(self):
        empty = self.root / "leopard1-source/sse2neon"
        empty.rmdir()
        with self.assertRaisesRegex(module.host.PreflightError, "exact uninitialized directory"):
            with self.authenticated(): pass
        empty.mkdir()
        (empty / "unexpected").write_bytes(b"not the pinned empty state")
        with self.assertRaisesRegex(module.host.PreflightError, "exact uninitialized directory"):
            with self.authenticated(): pass

    def test_commit_transcript_is_verified_not_trusted(self):
        query = module.PinnedSourceIdentity._query
        def altered(identity, role, arguments):
            data = query(identity, role, arguments)
            return data + b"tampered" if arguments[0] == "cat-file" else data
        with mock.patch.object(module.PinnedSourceIdentity, "_query", new=altered):
            with self.assertRaisesRegex(module.host.PreflightError, "raw commit does not authenticate"):
                with self.authenticated(): pass

    def test_caught_mmap_drift_remains_latched_after_restore(self):
        with self.source.open("r+b") as stream:
            mapping = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_WRITE)
        mapping[0] = mapping[0]
        try:
            with self.assertRaisesRegex(module.host.PreflightError, "failed source identity owner"):
                with self.authenticated() as identity:
                    original = mapping[0]
                    mapping[0] ^= 1
                    with self.assertRaises(FAILURES): identity.validate_current()
                    mapping[0] = original
                    with self.assertRaises(FAILURES): identity.record()
        finally:
            mapping.close()

    def test_preflight_loss_fails_closed_and_releases_git_descriptors(self):
        before = set(os.listdir("/proc/self/fd"))
        with self.assertRaises(FAILURES):
            with self.authenticated() as identity:
                self.preflight_live = False
                identity.record()
        self.assertEqual(set(os.listdir("/proc/self/fd")), before)

    def test_failed_authentication_releases_git_descriptors(self):
        before = set(os.listdir("/proc/self/fd"))
        self.pins["baseline_tree"] = "0" * 40
        with self.assertRaises(FAILURES):
            with self.authenticated(): pass
        self.assertEqual(set(os.listdir("/proc/self/fd")), before)

    def test_git_tree_serializer_rejects_ambiguous_or_unbounded_entries(self):
        entry = {"path": "a", "mode": "100644", "object_id": "a" * 40}
        invalid = [[], [entry, entry], [entry, dict(entry, path="a/b")],
                   [dict(entry, path="a/b"), entry], [dict(entry, mode="120000")],
                   [dict(entry, object_id="Z" * 40)], [dict(entry, extra=1)]]
        invalid += [[dict(entry, path=path)] for path in ("../a", "a//b", "/a", ".git/HEAD", "a\\b", "a/./b")]
        for entries in invalid:
            with self.subTest(entries=entries), self.assertRaises(FAILURES): module.tree_id(entries)
        with mock.patch.object(module, "MAX_PATH_TOTAL_BYTES", 0), self.assertRaises(FAILURES):
            module.tree_id([entry])

    def test_blob_authentication_uses_bounded_buffers(self):
        large = self.source.parent / "large.bin"
        with large.open("wb") as stream: stream.truncate(32 << 20)
        pread = os.pread
        sizes = []
        def bounded(descriptor, count, offset):
            sizes.append(count)
            return pread(descriptor, count, offset)
        with module.streamed.StreamingSourceOwner(self.root) as owner:
            identity = module.PinnedSourceIdentity(owner, self.retained)
            tracemalloc.start()
            try:
                with mock.patch.object(module.os, "pread", new=bounded):
                    entries, _files = identity._worktree("candidate-source")
                _current, peak = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
            self.assertTrue(any(row["path"] == "large.bin" for row in entries))
            self.assertLess(peak, 3 << 20)
            self.assertLessEqual(max(sizes), module.streamed.HASH_BLOCK_BYTES)


if __name__ == "__main__":
    unittest.main()
