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

        changed = copy.deepcopy(self.raw)
        changed["inputs_after"]["executable"]["sha256"] = "f" * 64
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

    def test_v4_build_attestation_rejects_build_and_binary_mutations(self) -> None:
        inputs = runner.synthetic_inputs()
        original = runner.synthetic_build_manifest(inputs)

        def validate(value: dict) -> None:
            runner.derive_build_attestation(
                value, inputs["build_manifest"], inputs["git"]["tooling_commit"],
                inputs["git"]["core_commit"], inputs["source"],
                inputs["build_runner"], inputs["build_validator"],
                inputs["archive"], inputs["executable"])

        mutations = (
            lambda value: value.update(schema="leopard2-c7-build-run-manifest/v3"),
            lambda value: value["reproducibility"].update(
                comparison={"status": "not-run"}),
            lambda value: value["builds"][2]["compile_argv"].remove("-O2"),
            lambda value: value["builds"][2].update(sanitizer=True),
            lambda value: value["builds"][2]["executable"].update(
                sha256="f" * 64),
        )
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(original)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                validate(changed)

    def test_isolation_requires_observed_timing_and_duration_binding(self) -> None:
        pair = runner.synthetic_pair_lease(0, 1)
        inactive = runner.isolation_record(
            0, 1, pair, runner.synthetic_cpu_stat(0, user=0, idle=0),
            runner.synthetic_cpu_stat(0, user=0, idle=1),
            runner.synthetic_cpu_stat(1, user=0, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=1), 1, 2)
        with self.assertRaises(runner.EvidenceError):
            runner.validate_isolation(inactive, 0, 1)

        changed = copy.deepcopy(self.raw)
        changed["child"]["duration_ns"] += 1
        self.assert_raw_rejected(changed)

    def test_host_isolation_and_reservation_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.raw)
        changed["host_after"]["timing_cpu"]["frequency_policy"][
            "scaling_governor"] = "powersave"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        isolation = changed["isolation"]
        after_sibling = copy.deepcopy(isolation["after"]["sibling_cpu"])
        after_sibling["fields"]["system"] += 1
        after_sibling["total_jiffies"] += 1
        changed["isolation"] = runner.isolation_record(
            0, 1, isolation["pair_lease"], isolation["before"]["timing_cpu"],
            isolation["after"]["timing_cpu"],
            isolation["before"]["sibling_cpu"], after_sibling, 1, 2)
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

        build_path = self.root / "provenance/build-run-manifest-v4.json"
        build_original = build_path.read_bytes()
        build_path.write_bytes(b"{}\n")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)
        build_path.write_bytes(build_original)

    def test_live_build_validator_invocation_is_fail_closed(self) -> None:
        validator = self.root / runner.BUILD_VALIDATOR_RELATIVE
        validator.parent.mkdir(parents=True)
        validator.write_text("fixture\n", encoding="utf-8")
        build_manifest = self.root / "build-manifest-v4.json"
        build_manifest.write_text("{}\n", encoding="utf-8")
        successful = subprocess.CompletedProcess(
            ["validator"], 0, b"C7 evidence validation passed (live)\n", b"")
        with mock.patch.object(runner.subprocess, "run", return_value=successful) as call:
            runner.run_build_validator(self.root, build_manifest)
        argv = call.call_args.args[0]
        self.assertIn("--live", argv)
        self.assertIn("--require-checkout-head", argv)
        self.assertEqual(call.call_args.kwargs["env"], runner.VALIDATOR_ENVIRONMENT)

        changed = subprocess.CompletedProcess(["validator"], 0, b"portable only\n", b"")
        with mock.patch.object(runner.subprocess, "run", return_value=changed), \
             self.assertRaises(runner.EvidenceError):
            runner.run_build_validator(self.root, build_manifest)

    def test_verified_build_attestation_binds_exact_local_avx2_files(self) -> None:
        source = self.root / runner.SOURCE_RELATIVE
        build_runner = self.root / runner.BUILD_RUNNER_RELATIVE
        build_validator = self.root / runner.BUILD_VALIDATOR_RELATIVE
        archive = self.root / ".research/fixture/build/core-avx2/liblibleopard.a"
        executable = self.root / ".research/fixture/build/c7-avx2"
        for path, data in ((source, b"source\n"), (build_runner, b"runner\n"),
                           (build_validator, b"validator\n"),
                           (archive, b"!<arch>\n"), (executable, b"ELF fixture\n")):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        executable.chmod(0o755)
        identities = {
            "git": {"tooling_commit": "b" * 40, "core_commit": "e" * 40},
            "source": {**runner.file_identity(source, "source"),
                       "path": runner.SOURCE_RELATIVE.as_posix()},
            "build_runner": {**runner.file_identity(build_runner, "build-runner"),
                             "path": runner.BUILD_RUNNER_RELATIVE.as_posix()},
            "build_validator": {
                **runner.file_identity(build_validator, "build-validator"),
                "path": runner.BUILD_VALIDATOR_RELATIVE.as_posix()},
            "archive": runner.file_identity(archive, "archive"),
            "executable": runner.file_identity(executable, "executable"),
        }
        build_manifest = self.root / ".research/fixture/results/build-run-manifest.json"
        build_manifest.parent.mkdir(parents=True)
        build_manifest.write_bytes(runner.canonical_pretty_bytes(
            runner.synthetic_build_manifest(identities)))
        with mock.patch.object(runner, "run_build_validator"):
            identity, attestation = runner.verified_build_attestation(
                self.root, build_manifest, "b" * 40, "e" * 40,
                identities["source"], identities["build_runner"],
                identities["build_validator"], identities["archive"],
                identities["executable"], archive, executable)
        self.assertEqual(attestation["manifest"]["sha256"], identity["sha256"])
        self.assertEqual(attestation["avx2"]["optimization"], "-O2")
        self.assertEqual(identities["build_runner"]["mode"] & 0o111, 0)
        self.assertEqual(identities["build_validator"]["mode"] & 0o111, 0)

        with executable.open("ab") as stream:
            stream.write(b"trailing mutation")
        with mock.patch.object(runner, "run_build_validator"), \
             self.assertRaises(runner.EvidenceError):
            runner.verified_build_attestation(
                self.root, build_manifest, "b" * 40, "e" * 40,
                identities["source"], identities["build_runner"],
                identities["build_validator"], identities["archive"],
                runner.file_identity(executable, "executable"), archive, executable)

    def test_committed_python_tool_modes_are_nonexecutable_and_enforced(self) -> None:
        source_root = MODULE_PATH.parents[4]
        head = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            text=True).strip()
        for relative, kind in (
                (runner.BUILD_RUNNER_RELATIVE, "build-runner"),
                (runner.BUILD_VALIDATOR_RELATIVE, "build-validator")):
            entry = subprocess.check_output(
                ["git", "-C", str(source_root), "ls-tree", head, "--",
                 relative.as_posix()], text=True)
            self.assertTrue(entry.startswith("100644 blob "), entry)
            identity = runner.committed_file_identity(
                source_root, head, relative, kind)
            self.assertEqual(identity["mode"] & 0o111, 0)

        repository = self.root / "mode-repository"
        repository.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
        subprocess.run(["git", "config", "user.email", "fixture@example.com"],
                       cwd=repository, check=True)
        subprocess.run(["git", "config", "user.name", "Fixture"],
                       cwd=repository, check=True)
        tool = repository / "tool.py"
        tool.write_text("print('fixture')\n", encoding="utf-8")
        tool.chmod(0o644)
        subprocess.run(["git", "add", "tool.py"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "fixture"],
                       cwd=repository, check=True)
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        identity = runner.committed_file_identity(
            repository, commit, Path("tool.py"), "build-runner")
        self.assertEqual(identity["mode"] & 0o111, 0)
        tool.chmod(0o755)
        with self.assertRaises(runner.EvidenceError):
            runner.committed_file_identity(
                repository, commit, Path("tool.py"), "build-runner")

        inputs = runner.synthetic_inputs()
        self.assertEqual(inputs["build_runner"]["mode"] & 0o111, 0)
        self.assertEqual(inputs["build_validator"]["mode"] & 0o111, 0)
        runner.validate_input_snapshot(inputs)
        inputs["build_runner"]["mode"] = 0o755
        inputs["binding_sha256"] = runner.sha256_bytes(runner.canonical_bytes(
            {key: value for key, value in inputs.items()
             if key != "binding_sha256"}))
        with self.assertRaises(runner.EvidenceError):
            runner.validate_input_snapshot(inputs)

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
        path.chmod(0o600)
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

    def test_reservation_rejects_path_replacement(self) -> None:
        payload = {"benchmark_cpu": 0, "nonce": "fixture-nonce",
                   "owner": "unit-test", "reserved_sibling": 1,
                   "schema": runner.RESERVATION_SCHEMA, "status": "held"}
        encoded = runner.canonical_bytes(payload)
        path = self.root / "replacement-reservation.json"
        path.write_bytes(encoded)
        path.chmod(0o600)
        guard = runner.Reservation(path, 0, 1)
        with guard:
            path.rename(self.root / "old-reservation.json")
            path.write_bytes(encoded)
            path.chmod(0o600)
            with self.assertRaises(runner.EvidenceError):
                guard.validate_current()

    def test_pair_lease_serializes_pair_and_rejects_path_replacement(self) -> None:
        runtime_root = self.root / "runtime-root"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        first = runner.PairLease(7, 8, runtime_root=runtime_root)
        with first as identity:
            self.assertEqual(set(identity), {
                "device", "directory_device", "directory_inode", "inode", "lock",
                "path", "payload", "sha256"})
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(7, 8, runtime_root=runtime_root):
                    pass
            old = first.path.with_suffix(".old")
            first.path.rename(old)
            first.path.write_bytes(runner.canonical_bytes(
                runner.pair_lease_payload(7, 8)))
            first.path.chmod(0o600)
            with self.assertRaises(runner.EvidenceError):
                first.validate_current()

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
        build_manifest = self.root / "build-manifest.json"
        build_manifest.write_bytes(b"{}\n")
        executable.chmod(0o755)
        output = self.root / "failed-run"
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            source_root=self.root, cpu=0, sibling=1, timeout=1.0,
            expected_tooling_commit="a" * 40,
            expected_core_commit="a" * 40,
            reservation_file=self.root / "missing-reservation.json",
            build_manifest=build_manifest)
        with mock.patch.object(
                runner, "validate_topology",
                side_effect=runner.EvidenceError("synthetic topology failure")):
            with self.assertRaises(runner.EvidenceError):
                runner.run_campaign(options)
        failure_path = output / "failure.json"
        self.assertTrue(failure_path.is_file())
        failure = runner.strict_json(failure_path.read_bytes(), "failure")
        runner.validate_failure(failure, output, check_files=True)
        self.assertEqual(failure["schema"], runner.FAILURE_SCHEMA)
        self.assertEqual(failure["status"], "failed")
        self.assertEqual(failure["failure_kind"], "validation")
        self.assertIn("synthetic topology failure", failure["error"])
        self.assertIsNone(failure["child"])
        self.assertIsNotNone(failure["request"])
        self.assertFalse((output / "manifest.json").exists())

    def _run_child_failure(self, *, timeout: bool) -> tuple[dict, Path]:
        inputs = runner.synthetic_inputs()
        executable = self.root / ("timeout-executable" if timeout else "exit-executable")
        archive = self.root / ("timeout-archive.a" if timeout else "exit-archive.a")
        build_manifest = self.root / (
            "timeout-build-manifest.json" if timeout else "exit-build-manifest.json")
        reservation = self.root / (
            "timeout-reservation.json" if timeout else "exit-reservation.json")
        output = self.root / ("timeout-run" if timeout else "exit-run")
        executable.write_bytes(b"fixture")
        executable.chmod(0o755)
        archive.write_bytes(b"!<arch>\n")
        build_manifest.write_bytes(runner.canonical_pretty_bytes(
            runner.synthetic_build_manifest(inputs)))
        reservation_payload = {
            "benchmark_cpu": 0, "nonce": "failure-fixture", "owner": "unit-test",
            "reserved_sibling": 1, "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        reservation.write_bytes(runner.canonical_bytes(reservation_payload))
        reservation.chmod(0o600)
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            build_manifest=build_manifest, source_root=self.root, cpu=0, sibling=1,
            timeout=1.0, expected_tooling_commit="b" * 40,
            expected_core_commit="e" * 40, reservation_file=reservation)
        pair_guard = mock.MagicMock()
        pair_guard.__enter__.return_value = runner.synthetic_pair_lease(0, 1)
        pair_guard.__exit__.return_value = None
        snapshots = [
            runner.synthetic_cpu_stat(0, user=0, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=0),
            runner.synthetic_cpu_stat(0, user=1, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=1),
        ]

        if timeout:
            def child_effect(*unused_args: object, **unused_kwargs: object) -> object:
                raise subprocess.TimeoutExpired(
                    ["fixture"], 1.0, output=b"partial", stderr=b"timed out")
        else:
            def child_effect(*unused_args: object, **unused_kwargs: object) -> object:
                return subprocess.CompletedProcess(
                    ["fixture"], 7, b"", b"failed\n")
        with mock.patch.object(runner, "validate_topology", return_value=({0, 1, 2}, {2})), \
             mock.patch.object(runner, "host_identity",
                               return_value=runner.synthetic_host(0, 1)), \
             mock.patch.object(runner, "PairLease", return_value=pair_guard), \
             mock.patch.object(runner, "input_snapshot", return_value=inputs), \
             mock.patch.object(runner, "cpu_stat_snapshot", side_effect=snapshots), \
             mock.patch.object(runner.os, "sched_getaffinity", return_value={0, 1, 2}), \
             mock.patch.object(runner.os, "sched_setaffinity"), \
             mock.patch.object(runner.subprocess, "run", side_effect=child_effect):
            with self.assertRaises(runner.EvidenceError):
                runner.run_campaign(options)
        failure = runner.strict_json((output / "failure.json").read_bytes(), "failure")
        runner.validate_failure(failure, output, check_files=True)
        self.assertIsNotNone(failure["request"])
        self.assertIsNotNone(failure["inputs_before"])
        self.assertIsNotNone(failure["host_before"])
        self.assertIsNotNone(failure["reservation"])
        self.assertIsNotNone(failure["build_provenance"])
        self.assertIsNotNone(failure["isolation"])
        self.assertFalse((output / "manifest.json").exists())
        stdout = runner.read_artifact(output, failure["child"]["stdout"], "stdout")
        stderr = runner.read_artifact(output, failure["child"]["stderr"], "stderr")
        self.assertEqual(stdout, b"partial" if timeout else b"")
        self.assertEqual(stderr, b"timed out" if timeout else b"failed\n")
        return failure, output

    def test_nonzero_child_retains_failure_evidence(self) -> None:
        failure, _output = self._run_child_failure(timeout=False)
        self.assertEqual(failure["child"]["returncode"], 7)
        self.assertFalse(failure["child"]["timed_out"])
        self.assertEqual(failure["failure_kind"], "child-exit")
        self.assertEqual(failure["error"], "C7 child exited 7")

    def test_timeout_child_retains_failure_evidence(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        self.assertEqual(failure["child"]["returncode"], 124)
        self.assertTrue(failure["child"]["timed_out"])
        self.assertEqual(failure["failure_kind"], "timeout")
        self.assertEqual(failure["error"], "C7 child timed out")

        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
             str(output / "failure.json")], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["status"], "VERIFIED_FAILURE")
        self.assertEqual(report["failure_kind"], "timeout")

    def test_failure_coordinated_resigned_mutations_are_rejected(self) -> None:
        failure, output = self._run_child_failure(timeout=True)

        mutations = (
            lambda value: value.update(status="pass"),
            lambda value: value.update(failure_kind="child-exit"),
            lambda value: value.update(error="C7 child exited 124"),
            lambda value: value["child"].update(returncode=0),
            lambda value: value["child"].update(timed_out=False),
            lambda value: (
                value.update(failure_kind="validation"),
                value["child"].update(returncode=0, timed_out=False)),
            lambda value: value["child"]["stdout"].update(sha256="f" * 64),
            lambda value: value["child"]["stderr"].update(path="../outside"),
            lambda value: value["build_provenance"]["attestation"]["avx2"].update(
                optimization="-O0"),
            lambda value: value["reservation"]["payload"].update(status="released"),
            lambda value: value["pair_lease"].update(inode=99),
        )
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(failure)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                runner.validate_failure(
                    resign(changed), output, check_files=True)

    def test_malformed_nested_failure_is_caught_by_cli(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        malformed_values = []
        changed = copy.deepcopy(failure)
        changed["isolation"] = {"timing_cpu": 0, "sibling_cpu": 1}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["isolation"]["before"]["timing_cpu"] = None
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["child"] = {"timed_out": True}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["build_provenance"] = {"manifest": None}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["pair_lease"] = {"payload": None}
        malformed_values.append(changed)

        for index, changed in enumerate(malformed_values):
            path = output / f"malformed-{index}.json"
            path.write_bytes(runner.canonical_bytes(resign(changed)) + b"\n")
            completed = subprocess.run(
                [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
                 str(path)], check=False, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                     "PYTHONHASHSEED": "0"})
            with self.subTest(index=index):
                self.assertEqual(completed.returncode, 1)
                self.assertIn("C7 authoritative evidence error:", completed.stderr)
                self.assertNotIn("Traceback (most recent call last)", completed.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
