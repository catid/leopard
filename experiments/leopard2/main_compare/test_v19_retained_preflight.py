#!/usr/bin/python3
"""Offline retained-preflight tests for leopard-79h.38.5.4.8.2.2.2."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import mmap
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("v19_preflight_tested", HERE / "v19_retained_preflight.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)
RUN_ROOT = "/home/catid/leopard-v19-build-preflight.QSq8egne"


def encoded(value):
    return module.host.canonical_bytes(value) + b"\n"


class RetainedPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preregistration = subprocess.run(
            [str(HERE / "run_authoritative_v17_gfni_main_compare.sh"),
             "--print-conditioned-v19-preregistration"], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        cls.original = module.host.load_preregistration(cls.preregistration)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-preflight-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name) / "preflight"
        self.root.mkdir(mode=0o700)
        (self.root / "artifacts").mkdir(mode=0o700)
        self.contract = copy.deepcopy(self.original)
        preflight = self.contract["build_preflight"]
        self.source = {"candidate_commit": preflight["source_commit"], "candidate_tree": preflight["source_tree"],
                       "baseline_commit": preflight["baseline_commit"], "baseline_tree": preflight["baseline_tree"]}
        source_bytes = "".join(f"{key}={value}\n" for key, value in self.source.items()).encode("ascii")
        events = b"low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n"
        self.files = {"source-identities.txt": source_bytes, "run-root.txt": (RUN_ROOT + "\n").encode(),
            "terminal.env": b"exit_status=0\nstage=complete\n" + source_bytes +
                b"memory_max=536870912\nmemory_swap_max=0\nbenchmark_invoked=false\n",
            "build-stage.txt": b"complete\n", "provenance-validator-nofile65536.status": b"0\n",
            "memory.events": events, "memory.events.local": events,
            "memory.peak": b"391512064\n", "memory.current": b"123456\n", "empty-diagnostic.log": b""}
        artifacts = {}
        for role, (build_role, _filename, key, _mode) in module.owners.ARTIFACTS.items():
            data = (b"\x7fELF" if role.endswith("executable") else b"!<arch>\n") + role.encode()
            preflight[key] = hashlib.sha256(data).hexdigest()
            artifacts[role + "_sha256"] = preflight[key]
            name = build_role + ("" if role.endswith("executable") else "-archive")
            for generation in ("first", "replay"):
                self.files[f"artifacts/{name}-{generation}"] = data
        self.provenance = {"schema": "leopard2-production-build-closure/v2",
            "source_root": RUN_ROOT + "/candidate-source", "physical_source_root": RUN_ROOT + "/candidate-source",
            "build_root": RUN_ROOT + "/candidate-build", "executable_target": "bench_leopard2",
            "tracked_source_manifest": {"git": {"commit": self.source["candidate_commit"],
                "tree": self.source["candidate_tree"], "dirty": False,
                "status_sha256": hashlib.sha256(b"").hexdigest()}},
            "executable": {"path": RUN_ROOT + "/candidate-build/bench_leopard2", "sha256": preflight["candidate_binary_sha256"]},
            "archive": {"path": RUN_ROOT + "/candidate-build/libleopard.a", "sha256": preflight["candidate_archive_sha256"]}}
        self.result = {"schema": "leopard2-v19-gf16-build-preflight-result/v1", "status": "passed",
            "canonical_lock_held": True, "v19_attempt_consumed": False, "source": copy.deepcopy(self.source),
            "host": {"hostname": "ripper", "shared_host": True, "swap_rows": 0},
            "build": {"configuration": "Release", "generator": "Unix Makefiles", "compiler": "c++ 13.3.0",
                "parallelism": 1, "same_path_replay": True, "candidate_target": "bench_leopard2",
                "baseline_target": "leopard_main_benchmark",
                **{role + "_byte_identical": True for role in module.owners.ARTIFACTS}},
            "resource_envelope": {"benchmark_invoked": False, "memory_max_bytes": 536870912,
                "memory_swap_max_bytes": 0, "memory_peak_bytes": 391512064, "oom": 0, "oom_kill": 0,
                "memory_events_sha256": hashlib.sha256(events).hexdigest()},
            "provenance": {"schema": "leopard2-production-build-closure/v2", "validated": True, "validated_nofile": 65536},
            "artifacts": artifacts}
        self.regenerate()
        original_loader = module.host.load_preregistration
        def fixture_loader(data):
            original_loader(data)  # Still rejects modified or armed authority bytes.
            return copy.deepcopy(self.contract)
        patcher = mock.patch.object(module.host, "load_preregistration", side_effect=fixture_loader)
        patcher.start()
        self.addCleanup(patcher.stop)

    def regenerate(self):
        self.files["candidate-build-provenance.json"] = encoded(self.provenance)
        digest = hashlib.sha256(self.files["candidate-build-provenance.json"]).hexdigest()
        self.contract["build_preflight"]["provenance_record_sha256"] = digest
        self.result["provenance"]["candidate_record_sha256"] = digest
        self.files["RESULT.json"] = encoded(self.result)
        self.reseal()

    def reseal(self):
        manifest = b"".join((hashlib.sha256(data).hexdigest() + "  ./" + relative + "\n").encode()
                            for relative, data in sorted(self.files.items()))
        self.contract["build_preflight"]["outer_sha256sums_sha256"] = hashlib.sha256(manifest).hexdigest()
        for relative, data in dict(self.files, SHA256SUMS=manifest).items():
            path = self.root / relative
            if path.exists(): path.chmod(0o600)
            path.write_bytes(data)
            path.chmod(0o444)

    def capture(self):
        return module.PinnedPreflight(self.preregistration, _root=self.root)

    def test_valid_retained_tree_is_read_only_and_deep_copied(self):
        before = {str(path.relative_to(self.root)): path.stat().st_mtime_ns for path in self.root.rglob("*")}
        with mock.patch.object(module.os, "mkdir", side_effect=AssertionError("write")), \
             mock.patch.object(subprocess, "run", side_effect=AssertionError("execution")):
            with self.capture() as retained:
                result = retained.record()
                self.assertEqual(result["file_count"], len(self.files) + 1)
                self.assertEqual(result["canonical_run_root"], RUN_ROOT)
                self.assertTrue(result["retained_preflight_verified"])
                for key in ("new_build_verified", "live_acquisition_armed", "v18_archives_verified"):
                    self.assertFalse(result[key])
                result["source_pins"].clear()
                self.assertEqual(len(retained.record()["source_pins"]), 4)
        self.assertEqual(before, {str(path.relative_to(self.root)): path.stat().st_mtime_ns for path in self.root.rglob("*")})
        with self.assertRaises(module.host.PreflightError): retained.record()
        with self.assertRaises(module.host.PreflightError): retained.__enter__()

    def test_preregistration_and_outer_manifest_are_pinned(self):
        with self.assertRaises(module.host.PreflightError):
            module.PinnedPreflight(self.preregistration + b"\n", _root=self.root)
        path = self.root / "SHA256SUMS"
        path.chmod(0o600)
        path.write_bytes(path.read_bytes() + b"\n")
        path.chmod(0o444)
        with self.assertRaisesRegex(module.host.PreflightError, "hash differs"):
            with self.capture(): pass

    def test_checksum_parser_is_closed_world(self):
        prefix = "0" * 64 + "  ./"
        for text in (prefix + "../escape\n", prefix + "x//y\n", prefix + "SHA256SUMS\n",
                     prefix + "x\n" + prefix + "x\n", prefix + "z\n" + prefix + "a\n",
                     prefix + "x\n" + prefix + "x/y\n", prefix + "x", prefix + "x\\y\n",
                     prefix + "x\r\n", prefix + "x\v" + prefix + "y\n",
                     prefix + "x\r" + prefix + "y\n",
                     "0" * 64 + " */x\n"):
            with self.subTest(text=text), self.assertRaises(module.host.PreflightError):
                module.checksum_entries(text.encode())
        with self.assertRaises(module.host.PreflightError):
            module.checksum_entries(b"x" * 16385)

    def test_inventory_and_byte_bounds_are_enforced(self):
        prefix = "0" * 64 + "  ./"
        for text in ("".join(prefix + f"{index:03d}\n" for index in range(65)),
                     prefix + "/".join("a" for _ in range(9)) + "\n"):
            with self.subTest(text=text), self.assertRaises(module.host.PreflightError):
                module.checksum_entries(text.encode())
        with mock.patch.object(module, "MAX_TOTAL_BYTES", 1), \
             self.assertRaisesRegex(module.host.PreflightError, "total exceeds"):
            with self.capture(): pass
        with mock.patch.object(module, "MAX_JSON_BYTES", 1), \
             self.assertRaisesRegex(module.host.PreflightError, "JSON exceeds"):
            with self.capture(): pass

    def test_directory_permissions_and_replacement_are_rejected(self):
        directory = self.root / "artifacts"
        directory.chmod(0o755)
        with self.assertRaisesRegex(module.host.PreflightError, "owner-only"):
            with self.capture(): pass
        directory.chmod(0o700)
        entered = False
        with self.assertRaises(FAILURES):
            with self.capture() as retained:
                entered = True
                moved = self.root.parent / "moved-artifacts"
                directory.rename(moved)
                moved.rename(directory)
                retained.validate_current()
        self.assertTrue(entered)

    def test_missing_or_extra_tree_entries_fail(self):
        for kind in ("missing", "extra"):
            path = self.root / ("build-stage.txt" if kind == "missing" else "unexpected")
            if kind == "missing": path.unlink()
            else: path.write_bytes(b"unexpected")
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.capture(): pass
            if kind == "extra": path.unlink()
            self.reseal()

    def test_unsafe_metadata_and_aliases_fail(self):
        path = self.root / "build-stage.txt"
        for kind in ("writable", "symlink", "hardlink", "fifo", "oversized"):
            path.unlink()
            other = self.root.parent / "other"
            try:
                if kind in ("symlink", "hardlink"):
                    other.write_bytes(b"complete\n")
                    other.chmod(0o444)
                    if kind == "symlink": path.symlink_to(other)
                    else: os.link(other, path)
                elif kind == "fifo": os.mkfifo(path)
                else:
                    path.write_bytes(b"complete\n")
                    if kind == "oversized":
                        with path.open("r+b") as stream: stream.truncate(module.MAX_FILE_BYTES + 1)
                    path.chmod(0o600 if kind == "writable" else 0o444)
                with self.subTest(kind=kind), self.assertRaises(FAILURES):
                    with self.capture(): pass
            finally:
                path.unlink()
                if other.exists(): other.unlink()
                self.reseal()

    def test_coherently_resealed_result_cannot_relax_contract(self):
        original = copy.deepcopy(self.result)
        for field, key, replacement in (("build", "parallelism", True), ("build", "parallelism", 2),
                ("host", "hostname", "foureyes"), ("host", "swap_rows", False),
                ("source", "candidate_commit", "a" * 40),
                ("resource_envelope", "memory_swap_max_bytes", 1),
                ("resource_envelope", "benchmark_invoked", True),
                ("provenance", "validated_nofile", 1024),
                ("artifacts", "baseline_archive_sha256", "0" * 64)):
            self.result = copy.deepcopy(original)
            self.result[field][key] = replacement
            self.regenerate()
            with self.subTest(field=field, key=key, replacement=replacement), self.assertRaises(FAILURES):
                with self.capture(): pass

    def test_resealed_provenance_cannot_splice_paths_or_source(self):
        original = copy.deepcopy(self.provenance)
        for kind in ("build path", "source path", "source commit", "source tree", "dirty", "archive"):
            self.provenance = copy.deepcopy(original)
            if kind == "build path": self.provenance["build_root"] = RUN_ROOT + "/foreign-build"
            if kind == "source path": self.provenance["source_root"] = "/tmp/foreign-source"
            if kind == "source commit": self.provenance["tracked_source_manifest"]["git"]["commit"] = "0" * 40
            if kind == "source tree": self.provenance["tracked_source_manifest"]["git"]["tree"] = "0" * 40
            if kind == "dirty": self.provenance["tracked_source_manifest"]["git"]["dirty"] = True
            if kind == "archive": self.provenance["archive"]["sha256"] = "0" * 64
            self.regenerate()
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.capture(): pass

    def test_resealed_control_files_and_replay_bytes_fail(self):
        originals = dict(self.files)
        for path, data in (("terminal.env", self.files["terminal.env"].replace(b"exit_status=0", b"exit_status=1")),
                ("source-identities.txt", self.files["source-identities.txt"].replace(b"cf7a705", b"0000000")),
                ("provenance-validator-nofile65536.status", b"1\n"), ("build-stage.txt", b"failed\n"),
                ("memory.peak", b"536870913\n"), ("memory.current", b"536870913\n"),
                ("memory.events", b"low 0\nhigh 0\nmax 0\noom 0\noom_kill 1\n"),
                ("run-root.txt", b"/tmp/other\n"), ("artifacts/baseline-replay", b"different replay")):
            self.files = dict(originals, **{path: data})
            self.reseal()
            with self.subTest(path=path), self.assertRaises(FAILURES):
                with self.capture(): pass

    def test_tree_aba_and_content_change_are_rejected_while_held(self):
        for kind in ("rename", "write", "extra"):
            entered = False
            with self.subTest(kind=kind), self.assertRaises(FAILURES):
                with self.capture() as retained:
                    entered = True
                    path = self.root / "build-stage.txt"
                    if kind == "rename":
                        path.rename(self.root / "temporary")
                        (self.root / "temporary").rename(path)
                    elif kind == "write":
                        path.chmod(0o600)
                        path.write_bytes(b"complete\n")
                        path.chmod(0o444)
                    else: (self.root / "unexpected").write_bytes(b"x")
                    retained.validate_current()
            self.assertTrue(entered)
            extra = self.root / "unexpected"
            if extra.exists(): extra.unlink()
            self.reseal()

    def test_existing_writable_mapping_is_rehashed(self):
        path = self.root / "artifacts/candidate-first"
        path.chmod(0o600)
        with path.open("r+b") as stream:
            mapping = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_WRITE)
        mapping[0] = mapping[0]
        path.chmod(0o444)
        entered = False
        try:
            with self.assertRaises(FAILURES):
                with self.capture() as retained:
                    entered = True
                    mapping[4] ^= 1
                    retained.validate_current()
            self.assertTrue(entered)
        finally:
            mapping.close()

    def test_failed_capture_closes_descriptors(self):
        before = set(os.listdir(f"/proc/{os.getpid()}/fd"))
        self.result["status"] = "failed"
        self.regenerate()
        with self.assertRaises(FAILURES):
            with self.capture(): pass
        self.assertEqual(set(os.listdir(f"/proc/{os.getpid()}/fd")), before)

    def test_duplicate_json_keys_are_rejected_even_with_matching_manifest(self):
        self.files["RESULT.json"] = b'{"status":"passed","status":"failed"}\n'
        self.reseal()
        with self.assertRaisesRegex(module.host.PreflightError, "duplicate keys"):
            with self.capture(): pass


if __name__ == "__main__":
    unittest.main()
