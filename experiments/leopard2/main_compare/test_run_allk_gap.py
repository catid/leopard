#!/usr/bin/env python3
"""Fail-closed source-identity tests for the all-K diagnostic runner."""

from __future__ import annotations

import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("run_allk_gap.py")
SPEC = importlib.util.spec_from_file_location("leopard2_run_allk_gap", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load run_allk_gap.py")
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def expect_rejected(test: unittest.TestCase, callback, text: str) -> None:
    with test.assertRaisesRegex(RuntimeError, text):
        callback()


class AllKIdentityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.main_commit = "a" * 40
        self.current_source = {
            "path": "/source", "head": "b" * 40, "tree": "c" * 40,
            "tracked_status": "clean",
        }
        self.main_snapshot = {"sha256": "d" * 64}
        self.current_snapshot = {"sha256": "e" * 64}

    def records(self):
        records = []
        for role in runner.ORDER:
            if role == "main":
                build = {"main_source_commit": self.main_commit}
            else:
                build = {
                    "source_commit": self.current_source["head"],
                    "source_tree": self.current_source["tree"],
                    "source_tracked_dirty": False,
                }
            result = {"build": build}
            if role == "leopard2":
                result.update({
                    "schema": "leopard2-benchmark-v3",
                    "parameters": {"report_decode_path": True},
                })
            snapshot = (self.main_snapshot if role == "main" else
                        self.current_snapshot)
            records.append({"role": role,
                            "executable_snapshot_sha256": snapshot["sha256"],
                            "returncode": 0,
                            "result": result})
        return records

    def attestation_record(self):
        return {
            "role": "leopard2", "returncode": 0,
            "executable_snapshot_sha256": self.current_snapshot["sha256"],
            "result": {
                "schema": "leopard2-benchmark-v5",
                "parameters": {"attest_source": True},
                "resolved": {
                    "profile": "legacy_high_v1", "field": "gf8",
                    "backend": "avx2",
                },
                "build": {
                    "source_commit": self.current_source["head"],
                    "source_tree": self.current_source["tree"],
                    "source_tracked_dirty": False,
                },
                "correctness": {"leopard2_round_trip": True},
            },
        }

    def correctness_records(self, cell):
        digest = {
            "algorithm": "fnv1a64",
            "original_data": "0123456789abcdef",
            "transmitted_parity": "123456789abcdef0",
            "recovered_originals": "23456789abcdef01",
        }
        records = self.records()
        for record in records:
            record["result"]["parameters"] = {
                "K": cell.k, "R": cell.r,
                "shard_bytes": cell.shard_bytes,
                "loss_count": cell.losses, "batch": 1,
                "reuse": cell.reuse, "iterations": cell.iterations,
                "warmup": cell.warmup, "thread_count": 1,
                "seed": cell.seed,
            }
            record["result"]["correctness"] = (
                {"round_trip": True} if record["role"] == "main" else
                {"leopard2_round_trip": True})
            record["result"]["workload_digests"] = copy.deepcopy(digest)
        return records

    def test_embedded_identities_accept_exact_values(self) -> None:
        runner.validate_invocation_identities(
            self.records(), self.main_commit, self.current_source,
            self.main_snapshot, self.current_snapshot)
        runner.validate_source_attestation(
            self.attestation_record(), self.current_source,
            self.current_snapshot)

    def test_pure_avx2_disassembly_gate(self) -> None:
        avx2 = """
  10: c5 fd ef c0          vpxor  ymm0,ymm0,ymm0
  14: c5 fe 7f 07          vmovdqu YMMWORD PTR [rdi],ymm0
"""
        identity = runner.inspect_isa_disassembly(avx2)
        self.assertEqual(identity["evex_prefixed_instruction_count"], 0)
        self.assertEqual(identity["ymm_operand_instruction_count"], 2)
        evex = avx2 + \
            "  20: 62 f1 7d 48 ef c0    vpxord zmm0,zmm0,zmm0\n"
        self.assertEqual(
            runner.inspect_isa_disassembly(evex)[
                "evex_prefixed_instruction_count"], 1)

    def test_gf8_only_matrix_is_all_k_and_stable(self) -> None:
        cells = runner.make_cells(gf8_only=True)
        self.assertEqual(len(cells), 2522)
        self.assertEqual({cell.k for cell in cells}, set(range(1, 256)))
        self.assertEqual({cell.family for cell in cells}, {"gf8-all-k"})
        self.assertEqual(runner.make_cells()[:len(cells)], cells)

    def test_current_command_forces_profile_field_and_avx2(self) -> None:
        gf8 = runner.Cell("one", "gf8-all-k", 3, 2, 4096, 1,
                          "low-R", "one-loss", 7, 5, 16, 1)
        command = runner.command(
            "leopard2", Path("/tmp/bench"), gf8, 2, False)
        self.assertEqual(command[command.index("--profile") + 1], "high")
        self.assertEqual(command[command.index("--field") + 1], "gf8")
        self.assertEqual(command[command.index("--backend") + 1], "avx2")
        self.assertIn("--report-decode-path", command)
        self.assertNotIn("--attest-source", command)
        self.assertIn("--skip-legacy", command)

    def test_path_classification_consumes_reported_route(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 65, 64, 4096, 4,
                           "max-GF8-R", "max-loss", 7, 5, 16, 1)
        result = {"resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "avx2", "parent_count": 256,
            "padded_side": 64, "selected_decode_path": "tiled",
            "selected_decode_rule": "workspace_tiled",
            "decode_required_work_slots": 132,
        }}
        selected = runner.classify_paths(cell, result)
        self.assertEqual(selected["decode_path"], "tiled")
        self.assertEqual(selected["decode_rule"], "workspace_tiled")
        self.assertEqual(selected["decode_required_work_slots"], 132)
        del result["resolved"]["selected_decode_path"]
        expect_rejected(
            self, lambda: runner.classify_paths(cell, result),
            "selected decode path")

    def test_embedded_identities_reject_every_mismatch(self) -> None:
        mutations = []
        wrong_main = self.records()
        wrong_main[0]["result"]["build"]["main_source_commit"] = "d" * 40
        mutations.append((wrong_main, "exact-main"))
        wrong_snapshot = self.records()
        wrong_snapshot[2]["executable_snapshot_sha256"] = "f" * 64
        mutations.append((wrong_snapshot, "snapshot mismatch"))
        wrong_role = self.records()
        wrong_role[0]["role"], wrong_role[1]["role"] = \
            wrong_role[1]["role"], wrong_role[0]["role"]
        mutations.append((wrong_role, "role order"))
        missing = self.records()[:-1]
        mutations.append((missing, "exactly four"))
        for records, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda records=records: runner.validate_invocation_identities(
                        records, self.main_commit, self.current_source,
                        self.main_snapshot, self.current_snapshot),
                    message)

        attestation_mutations = []
        wrong_commit = self.attestation_record()
        wrong_commit["result"]["build"]["source_commit"] = "d" * 40
        attestation_mutations.append((wrong_commit, "embedded commit"))
        wrong_tree = self.attestation_record()
        wrong_tree["result"]["build"]["source_tree"] = "d" * 40
        attestation_mutations.append((wrong_tree, "embedded tree"))
        dirty = self.attestation_record()
        dirty["result"]["build"]["source_tracked_dirty"] = True
        attestation_mutations.append((dirty, "tracked-dirty"))
        wrong_backend = self.attestation_record()
        wrong_backend["result"]["resolved"]["backend"] = "auto"
        attestation_mutations.append((wrong_backend, "codec identity"))
        for record, message in attestation_mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda record=record: runner.validate_source_attestation(
                        record, self.current_source, self.current_snapshot),
                    message)

    def test_git_identity_requires_requested_clean_top_level(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-git-") as directory:
            root = Path(directory)
            subprocess.run(["/usr/bin/git", "init", "-q", str(root)], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                            "user.name", "Leopard2 Test"], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                            "user.email", "test@example.invalid"], check=True)
            tracked = root / "tracked.txt"
            tracked.write_text("one\n", encoding="utf-8")
            subprocess.run(["/usr/bin/git", "-C", str(root), "add",
                            "tracked.txt"], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "commit", "-q",
                            "-m", "fixture"], check=True)
            commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
                text=True).strip()
            identity = runner.git_identity(root, commit)
            self.assertEqual(identity["head"], commit)
            self.assertEqual(identity["tracked_status"], "clean")
            expect_rejected(
                self, lambda: runner.git_identity(root, "0" * 40),
                "HEAD mismatch")
            child = root / "child"
            child.mkdir()
            expect_rejected(
                self, lambda: runner.git_identity(child, commit),
                "not the Git top level")
            child.rmdir()
            untracked = root / "untracked.txt"
            untracked.write_text("untracked\n", encoding="utf-8")
            expect_rejected(
                self, lambda: runner.git_identity(root, commit),
                "tracked or untracked modifications")
            untracked.unlink()
            tracked.write_text("two\n", encoding="utf-8")
            expect_rejected(
                self, lambda: runner.git_identity(root, commit),
                "tracked or untracked modifications")

    def test_git_identity_rejects_hidden_index_and_replace_state(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-git-flags-") as directory:
            root = Path(directory)
            subprocess.run(["/usr/bin/git", "init", "-q", str(root)], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                            "user.name", "Leopard2 Test"], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                            "user.email", "test@example.invalid"], check=True)
            tracked = root / "tracked.txt"
            tracked.write_text("one\n", encoding="utf-8")
            subprocess.run(["/usr/bin/git", "-C", str(root), "add",
                            "tracked.txt"], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "commit", "-q",
                            "-m", "one"], check=True)
            first_commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
                text=True).strip()

            subprocess.run(["/usr/bin/git", "-C", str(root), "update-index",
                            "--assume-unchanged", "tracked.txt"], check=True)
            tracked.write_text("hidden\n", encoding="utf-8")
            expect_rejected(
                self, lambda: runner.git_identity(root, first_commit),
                "non-default flag")
            subprocess.run(["/usr/bin/git", "-C", str(root), "update-index",
                            "--no-assume-unchanged", "tracked.txt"], check=True)
            tracked.write_text("one\n", encoding="utf-8")

            subprocess.run(["/usr/bin/git", "-C", str(root), "update-index",
                            "--skip-worktree", "tracked.txt"], check=True)
            expect_rejected(
                self, lambda: runner.git_identity(root, first_commit),
                "non-default flag")
            subprocess.run(["/usr/bin/git", "-C", str(root), "update-index",
                            "--no-skip-worktree", "tracked.txt"], check=True)

            tracked.write_text("two\n", encoding="utf-8")
            subprocess.run(["/usr/bin/git", "-C", str(root), "add",
                            "tracked.txt"], check=True)
            subprocess.run(["/usr/bin/git", "-C", str(root), "commit", "-q",
                            "-m", "two"], check=True)
            second_commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
                text=True).strip()
            subprocess.run(["/usr/bin/git", "-C", str(root), "replace",
                            first_commit, second_commit], check=True)
            branch = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "symbolic-ref", "HEAD"],
                text=True).strip()
            subprocess.run(["/usr/bin/git", "-C", str(root), "update-ref",
                            branch, first_commit], check=True)
            expect_rejected(
                self, lambda: runner.git_identity(root, first_commit),
                "tracked or untracked modifications")

    def test_correctness_requires_exact_workload_digests(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        exact = self.correctness_records(cell)
        runner.validate_correctness(cell, exact)

        mutations = []
        empty = self.correctness_records(cell)
        empty[0]["result"]["workload_digests"] = {}
        mutations.append((empty, "structure"))
        wrong_algorithm = self.correctness_records(cell)
        wrong_algorithm[1]["result"]["workload_digests"]["algorithm"] = "sha256"
        mutations.append((wrong_algorithm, "algorithm"))
        wrong_width = self.correctness_records(cell)
        wrong_width[2]["result"]["workload_digests"]["original_data"] = "0" * 15
        mutations.append((wrong_width, "original_data"))
        uppercase = self.correctness_records(cell)
        uppercase[3]["result"]["workload_digests"]["recovered_originals"] = \
            "ABCDEF0123456789"
        mutations.append((uppercase, "recovered_originals"))
        for records, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda records=records: runner.validate_correctness(
                        cell, records),
                    message)

    def test_manifest_rejects_stale_or_resigned_contract(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 1, 1, 64, 1,
                           "low-R", "one-loss", 1, 1, 1, 0)
        contract = {"source": self.current_source, "main": self.main_commit}
        digest = runner.canonical_digest(contract)
        manifest = {
            "schema": "leopard2-all-k-gap-manifest/v2",
            "run_contract": contract,
            "run_contract_sha256": digest,
            "cells": [runner.dataclasses.asdict(cell)],
            "completion": None,
        }
        runner.validate_manifest(manifest, contract, digest, [cell])

        stale = copy.deepcopy(manifest)
        stale["run_contract"]["main"] = "d" * 40
        stale["run_contract_sha256"] = runner.canonical_digest(
            stale["run_contract"])
        expect_rejected(
            self, lambda: runner.validate_manifest(stale, contract, digest, [cell]),
            "run contract")

        corrupted = copy.deepcopy(manifest)
        corrupted["run_contract_sha256"] = "e" * 64
        expect_rejected(
            self,
            lambda: runner.validate_manifest(corrupted, contract, digest, [cell]),
            "digest mismatch")

    def test_executable_identity_binds_bytes_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-exe-") as directory:
            executable = Path(directory) / "benchmark"
            executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            executable.chmod(0o755)
            first = runner.file_identity(executable, "fixture")
            self.assertEqual(len(first["sha256"]), 64)
            executable.write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
            second = runner.file_identity(executable, "fixture")
            self.assertNotEqual(first, second)
            expect_rejected(
                self, lambda: runner.snapshot_executable(executable, "fixture"),
                "not an ELF")
            executable.chmod(0o644)
            expect_rejected(
                self, lambda: runner.file_identity(executable, "fixture"),
                "not executable")

    def test_sealed_snapshot_executes_captured_bytes_after_path_swap(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-snapshot-") as directory:
            executable = Path(directory) / "benchmark"
            saved = Path(directory) / "saved"
            shutil.copy2("/bin/echo", executable)
            source, descriptor, snapshot = runner.snapshot_executable(
                executable, "fixture")
            try:
                executable.rename(saved)
                shutil.copy2("/bin/false", executable)
                command = [
                    "/usr/bin/taskset", "-c",
                    str(min(os.sched_getaffinity(0))), str(executable),
                    '{"value":"A"}',
                ]
                record = runner.run_one(
                    "main", command, 5.0, descriptor, snapshot)
                self.assertEqual(record["returncode"], 0)
                self.assertEqual(record["result"], {"value": "A"})
                self.assertEqual(record["command"], command)
                self.assertEqual(record["executable_snapshot_sha256"],
                                 snapshot["sha256"])
                with self.assertRaises(OSError):
                    os.write(descriptor, b"mutation")
                executable.unlink()
                saved.rename(executable)
                restored = runner.file_identity(executable, "fixture")
                for name in ("sha256", "device", "inode", "mode", "size",
                             "mtime_ns"):
                    self.assertEqual(restored[name], source[name])
                self.assertEqual(
                    runner.sealed_snapshot_identity(descriptor, "fixture"),
                    snapshot)
            finally:
                os.close(descriptor)

    def test_manifest_json_round_trip_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-json-") as directory:
            path = Path(directory) / "manifest.json"
            value = {"z": [3, 2, 1], "a": self.current_source}
            runner.write_json_atomic(path, value)
            self.assertEqual(json.loads(path.read_text(encoding="utf-8")), value)
            self.assertFalse(path.with_name(path.name + ".tmp").exists())


if __name__ == "__main__":
    unittest.main()
