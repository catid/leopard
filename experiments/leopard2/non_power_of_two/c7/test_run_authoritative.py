#!/usr/bin/env python3
"""Portable and mutation tests for the authoritative C7 runner."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


MODULE_PATH = Path(__file__).with_name("run_authoritative.py")
SPEC = importlib.util.spec_from_file_location("c7_run_authoritative", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def resign(value: dict) -> dict:
    payload = copy.deepcopy(value)
    payload.pop("digest", None)
    return runner.signed(payload)


class AuthoritativeRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(
            prefix="leopard2-c7-authoritative-test-")
        self.root = Path(self.temporary.name)
        self.manifest, self.raw, self.result = runner.synthetic_bundle(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def assert_raw_rejected(self, value: dict) -> None:
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(resign(value), self.root)

    def test_valid_fixture_and_portable_verify_mode(self) -> None:
        normalized = runner.validate_manifest(self.manifest, self.root)
        self.assertEqual(normalized["cell_count"], len(runner.EXPECTED_CELLS))
        manifest_path = self.root / "manifest.json"
        runner.write_json_exclusive(manifest_path, self.manifest)
        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify", "--manifest",
             str(manifest_path)], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["cells"], 12)

    def test_synthetic_raw_and_manifest_are_path_independent(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard2-c7-authoritative-peer-") as directory:
            peer_root = Path(directory)
            peer_manifest, peer_raw, peer_result = runner.synthetic_bundle(peer_root)
            self.assertEqual(runner.canonical_bytes(peer_manifest),
                             runner.canonical_bytes(self.manifest))
            self.assertEqual(runner.canonical_bytes(peer_raw),
                             runner.canonical_bytes(self.raw))
            self.assertEqual(runner.canonical_bytes(peer_result),
                             runner.canonical_bytes(self.result))
            self.assertEqual((peer_root / "raw.json").read_bytes(),
                             (self.root / "raw.json").read_bytes())

    def test_child_matrix_and_numeric_mutations_rejected(self) -> None:
        mutations = (
            lambda value: value.update(runtime_backend="scalar"),
            lambda value: value.update(affinity=[1]),
            lambda value: value.update(core_git_sha="f" * 40),
            lambda value: value["profile"].update(version=True),
            lambda value: value["correctness"].update(hot_path_allocations=False),
            lambda value: value["benchmarks"].reverse(),
            lambda value: value["benchmarks"][0].update(exact_field=True),
            lambda value: value["benchmarks"][0].update(
                padded_decode_scratch=False),
            lambda value: value["benchmarks"][0][
                "exact_encode_samples_us"].pop(),
            lambda value: value["benchmarks"][0]["exact_encode"].update(
                median_us=123.0),
        )
        inputs = runner.synthetic_inputs()
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(self.result)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                runner.validate_child_result(changed, 0, inputs)

    def test_request_child_and_input_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.raw)
        changed["request"]["cell_geometry"].reverse()
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["request"]["command"][2] = "1"
        changed["child"]["command"][2] = "1"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["returncode"] = False
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["environment"]["OMP_NUM_THREADS"] = "2"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["inputs_after"]["archive"]["sha256"] = "e" * 64
        changed["inputs_after"]["binding_sha256"] = runner.sha256_bytes(
            runner.canonical_bytes({key: value for key, value in
                                    changed["inputs_after"].items()
                                    if key != "binding_sha256"}))
        self.assert_raw_rejected(changed)

    def test_core_source_closure_mutations_rejected(self) -> None:
        for mutation in (
            lambda value: value["core_source_closure"].reverse(),
            lambda value: value["core_source_closure"][0].update(size=True),
            lambda value: value["core_source_closure"][0].update(
                path="LeopardFF8.cpp"),
        ):
            inputs = runner.synthetic_inputs()
            mutation(inputs)
            inputs["binding_sha256"] = runner.sha256_bytes(
                runner.canonical_bytes({key: value for key, value in inputs.items()
                                        if key != "binding_sha256"}))
            with self.assertRaises(runner.EvidenceError):
                runner.validate_input_snapshot(inputs)

    def test_host_isolation_and_reservation_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.raw)
        changed["host_after"]["timing_cpu"]["frequency_policy"][
            "scaling_governor"] = "powersave"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        before = changed["isolation"]["before"]
        after = copy.deepcopy(before)
        after["1"]["system"] += 1
        changed["isolation"] = runner.isolation_record(0, 1, before, after, 1, 2)
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["reservation"]["payload"]["status"] = "released"
        encoded = runner.canonical_bytes(changed["reservation"]["payload"])
        changed["reservation"]["bytes"] = len(encoded)
        changed["reservation"]["sha256"] = runner.sha256_bytes(encoded)
        self.assert_raw_rejected(changed)

    def test_retained_artifact_path_hash_and_bytes_mutations_rejected(self) -> None:
        result_path = self.root / "child/result.json"
        result_path.write_bytes(b"{}\n")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)

        result_path.write_bytes(runner.canonical_bytes(self.result) + b"\n")
        for path in ("../result.json", "child//result.json", "./child/result.json"):
            changed = copy.deepcopy(self.raw)
            changed["child"]["result"]["path"] = path
            with self.subTest(path=path):
                self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["stderr"]["size"] = True
        self.assert_raw_rejected(changed)

    def test_manifest_and_raw_coordinated_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.manifest)
        changed["validated_output"]["cell_count"] = 11
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(resign(changed), self.root)

        raw_path = self.root / "raw.json"
        pretty = (json.dumps(self.raw, indent=2, sort_keys=True) + "\n").encode()
        raw_path.write_bytes(pretty)
        changed = copy.deepcopy(self.manifest)
        changed["raw"] = {"path": "raw.json", "size": len(pretty),
                          "sha256": runner.sha256_bytes(pretty)}
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(resign(changed), self.root)

    def test_verify_rejects_noncanonical_manifest_bytes(self) -> None:
        manifest_path = self.root / "manifest.json"
        manifest_path.write_text(json.dumps(self.manifest, indent=2),
                                 encoding="utf-8")
        options = argparse.Namespace(manifest=manifest_path)
        with self.assertRaises(runner.EvidenceError):
            runner.verify_campaign(options)

    def test_reservation_lock_canonicality_and_current_bytes(self) -> None:
        payload = {"benchmark_cpu": 0, "nonce": "fixture-nonce",
                   "owner": "unit-test", "reserved_sibling": 1,
                   "schema": runner.RESERVATION_SCHEMA, "status": "held"}
        path = self.root / "reservation.json"
        path.write_bytes(runner.canonical_bytes(payload))
        guard = runner.Reservation(path, 0, 1)
        with guard:
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(path, 0, 1):
                    pass
            path.write_bytes(runner.canonical_bytes({**payload, "nonce": "changed-nonce"}))
            with self.assertRaises(runner.EvidenceError):
                guard.validate_current()

        path.write_bytes(runner.canonical_bytes(payload) + b"\n")
        with self.assertRaises(runner.EvidenceError):
            with runner.Reservation(path, 0, 1):
                pass
        path.write_bytes(runner.canonical_bytes(payload))
        with runner.Reservation(path, 0, 1):
            pass

    def test_git_identity_uses_exact_tree_and_rejects_untracked_files(self) -> None:
        repository = self.root / "repository"
        repository.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
        subprocess.run(["git", "config", "user.email", "fixture@example.com"],
                       cwd=repository, check=True)
        subprocess.run(["git", "config", "user.name", "Fixture"],
                       cwd=repository, check=True)
        (repository / "tracked.txt").write_text("fixture\n", encoding="utf-8")
        subprocess.run(["git", "add", "tracked.txt"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "fixture"],
                       cwd=repository, check=True)
        core_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        core_tree = subprocess.check_output(
            ["git", "rev-parse", "HEAD^{tree}"], cwd=repository, text=True).strip()
        (repository / "tooling.txt").write_text("tooling\n", encoding="utf-8")
        subprocess.run(["git", "add", "tooling.txt"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "tooling"],
                       cwd=repository, check=True)
        tooling_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        identity = runner.git_identity(repository, tooling_commit, core_commit)
        self.assertEqual(identity["core_tree"], core_tree)
        self.assertTrue(identity["core_is_ancestor"])
        with self.assertRaises(runner.EvidenceError):
            runner.git_identity(repository, tooling_commit, "f" * 40)
        (repository / "untracked.txt").write_text("untracked\n", encoding="utf-8")
        with self.assertRaises(runner.EvidenceError):
            runner.git_identity(repository, tooling_commit, core_commit)

    def test_preflight_failure_is_retained_without_running_child(self) -> None:
        executable = self.root / "executable"
        archive = self.root / "archive.a"
        executable.write_bytes(b"fixture")
        archive.write_bytes(b"!<arch>\n")
        executable.chmod(0o755)
        output = self.root / "failed-run"
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            source_root=self.root, cpu=0, sibling=1, timeout=1.0,
            expected_tooling_commit="a" * 40,
            expected_core_commit="a" * 40,
            reservation_file=self.root / "missing-reservation.json")
        with mock.patch.object(
                runner, "validate_topology",
                side_effect=runner.EvidenceError("synthetic topology failure")):
            with self.assertRaises(runner.EvidenceError):
                runner.run_campaign(options)
        failure_path = output / "failure.json"
        self.assertTrue(failure_path.is_file())
        failure = runner.strict_json(failure_path.read_bytes(), "failure")
        runner.validate_failure(failure)
        self.assertEqual(failure["schema"], runner.FAILURE_SCHEMA)
        self.assertIn("synthetic topology failure", failure["error"])
        self.assertIsNone(failure["child"])
        self.assertFalse((output / "manifest.json").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
