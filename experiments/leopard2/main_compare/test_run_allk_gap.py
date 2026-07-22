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


class ProductionBuildFixture:
    """Minimal byte-valid CMake build closure; no compiler execution needed."""

    def __init__(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(
            prefix="leo2-build-closure-")
        self.root = Path(self.temporary.name)
        self.source = self.root / "source"
        self.build = self.root / "build"
        self.source.mkdir()
        self.build.mkdir()
        sources = set(runner.CORE_LIBRARY_SOURCES)
        sources.update({
            "LeopardFF8.cpp", "LeopardFF16.cpp",
            "Leopard2BackendSSSE3.cpp",
            "Leopard2BackendAVX2.cpp", "bench/leopard2/benchmark.cpp",
        })
        for relative in sorted(sources):
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("// " + relative + "\n", encoding="utf-8")
        subprocess.run(["/usr/bin/git", "init", "-q", str(self.source)],
                       check=True)
        subprocess.run(["/usr/bin/git", "-C", str(self.source), "config",
                        "user.name", "Leopard2 Test"], check=True)
        subprocess.run(["/usr/bin/git", "-C", str(self.source), "config",
                        "user.email", "test@example.invalid"], check=True)
        subprocess.run(["/usr/bin/git", "-C", str(self.source), "add", "."],
                       check=True)
        subprocess.run(["/usr/bin/git", "-C", str(self.source), "commit",
                        "-qm", "fixture"], check=True)
        self.executable = self.build / "bench_leopard2"
        self.entries = []
        self.archive_objects = []
        compiler = str(Path("/usr/bin/c++").resolve(strict=True))
        for relative in sorted(sources - {"bench/leopard2/benchmark.cpp"}):
            if relative == "Leopard2BackendSSSE3.cpp":
                output = Path("CMakeFiles/leopard2_backend_ssse3.dir") / \
                    (Path(relative).name + ".o")
                flags = ["-mssse3", "-mno-avx"]
            elif relative.startswith("Leopard2BackendAVX2"):
                output = Path("CMakeFiles/leopard2_backend_avx2.dir") / \
                    (Path(relative).name + ".o")
                flags = ["-mavx2", "-mno-avx512f"]
            else:
                output = Path("CMakeFiles/leopard.dir") / \
                    (Path(relative).name + ".o")
                flags = []
            self._add_entry(relative, output, compiler, flags)
            self.archive_objects.append(output)
        benchmark_output = Path("CMakeFiles/bench_leopard2.dir/bench/leopard2") / \
            "benchmark.cpp.o"
        self._add_entry("bench/leopard2/benchmark.cpp", benchmark_output,
                        compiler, [])
        self.benchmark_object = benchmark_output
        self._write_cache()
        self._write_commands()
        archive_recipe = self.build / "CMakeFiles/leopard.dir/link.txt"
        archive_recipe.parent.mkdir(parents=True, exist_ok=True)
        archive_recipe.write_text(
            "/usr/bin/ar qc libleopard.a " +
            " ".join(path.as_posix() for path in self.archive_objects) +
            "\n/usr/bin/ranlib libleopard.a\n", encoding="utf-8")
        executable_recipe = self.build / \
            "CMakeFiles/bench_leopard2.dir/link.txt"
        executable_recipe.parent.mkdir(parents=True, exist_ok=True)
        executable_recipe.write_text(
            compiler + " -O3 " + self.benchmark_object.as_posix() +
            " -o bench_leopard2 libleopard.a\n", encoding="utf-8")
        self._relink()

    def close(self) -> None:
        self.temporary.cleanup()

    def _add_entry(self, relative, output, compiler, flags) -> None:
        source = (self.source / relative).resolve()
        target = self.build / output
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((relative + " object\n").encode())
        arguments = [compiler, "-O3", *flags, "-c", str(source),
                     "-o", output.as_posix()]
        self.entries.append({
            "directory": str(self.build.resolve()),
            "arguments": arguments, "file": str(source),
            "output": output.as_posix(),
        })

    def _write_cache(self, **overrides) -> None:
        values = {
            "CMAKE_BUILD_TYPE:STRING": "Release",
            "CMAKE_EXPORT_COMPILE_COMMANDS:BOOL": "ON",
            "CMAKE_GENERATOR:INTERNAL": "Unix Makefiles",
            "CMAKE_HOME_DIRECTORY:INTERNAL": str(self.source.resolve()),
            "CMAKE_CXX_COMPILER:FILEPATH": "/usr/bin/c++",
            "CMAKE_CXX_FLAGS:STRING": "",
            "CMAKE_CXX_FLAGS_RELEASE:STRING": "-O3 -DNDEBUG",
            "CMAKE_AR:FILEPATH": "/usr/bin/ar",
            "ENABLE_OPENMP:BOOL": "ON",
            "LEO2_BUILD_BENCHMARKS:BOOL": "ON",
            "LEO2_BUILD_FUZZERS:BOOL": "OFF",
            "LEO2_ENABLE_CUDA:BOOL": "OFF",
            "LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT:BOOL": "OFF",
            "LEOPARD_ENABLE_GF8:BOOL": "ON",
            "LEOPARD_ENABLE_GF16:BOOL": "ON",
            "LEO2_BACKEND_VARIANT:STRING": "auto",
            "LEO2_BUILD_ALLK_DIAGNOSTIC:BOOL": "OFF",
            "LEO2_BUILD_TESTS:BOOL": "ON",
            "LEO2_FLAG_MAVX2:INTERNAL": "1",
            "LEO2_FLAG_MNO_AVX512F:INTERNAL": "1",
            "LEO2_FLAG_MAVX512F:UNINITIALIZED": "FALSE",
            "LEO2_FLAG_MAVX512BW:UNINITIALIZED": "FALSE",
            "LEO2_FLAG_MAVX512VL:UNINITIALIZED": "FALSE",
        }
        values.update(overrides)
        (self.build / "CMakeCache.txt").write_text(
            "".join(f"{name}={value}\n" for name, value in values.items()),
            encoding="utf-8")

    def _write_commands(self) -> None:
        (self.build / "compile_commands.json").write_text(
            json.dumps(self.entries), encoding="utf-8")

    def _relink(self) -> None:
        archive = self.build / "libleopard.a"
        if archive.exists():
            archive.unlink()
        subprocess.run(
            ["/usr/bin/ar", "qc", str(archive),
             *(str(self.build / path) for path in self.archive_objects)],
            check=True)
        subprocess.run(["/usr/bin/ranlib", str(archive)], check=True)
        shutil.copyfile("/bin/true", self.executable)
        self.executable.chmod(0o755)


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
            subprocess.run(["/usr/bin/git", "-C", str(root), "reset",
                            "--hard", "-q", first_commit], check=True)
            expect_rejected(
                self, lambda: runner.git_identity(root, first_commit),
                "replace refs")

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


class ProductionBuildClosureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = ProductionBuildFixture()

    def tearDown(self) -> None:
        self.fixture.close()

    def provenance(self):
        return runner.candidate_build_provenance(
            self.fixture.build, self.fixture.source,
            self.fixture.executable, "bench_leopard2")

    def test_benchmark_target_source_contract_is_exact(self) -> None:
        self.assertEqual(runner.BENCHMARK_SOURCE_BY_TARGET, {
            "bench_leopard2": "bench/leopard2/benchmark.cpp",
            "bench_leopard2_allk": "bench/leopard2/benchmark.cpp",
            "bench_leopard2_no_loss_setup":
                "bench/leopard2/benchmark_no_loss_setup.cpp",
        })

    def test_exact_release_pure_avx2_closure_is_accepted(self) -> None:
        result = self.provenance()
        self.assertEqual(result["schema"],
                         "leopard2-production-build-closure/v1")
        self.assertEqual(result["archive_members"], [
            path.name for path in self.fixture.archive_objects])
        profiles = {record["flag_profile"] for record in
                    result["source_object_compile_closure"]}
        self.assertEqual(profiles, {
            "portable-core", "ssse3-no-avx", "avx2-no-avx512"})

    def test_source_attestation_alone_can_accept_a_stale_library_object(self) \
            -> None:
        source = runner.git_identity(
            self.fixture.source,
            subprocess.check_output([
                "/usr/bin/git", "-C", str(self.fixture.source),
                "rev-parse", "HEAD"], text=True).strip())
        snapshot = {"sha256": "a" * 64}
        record = {
            "role": "leopard2", "returncode": 0,
            "executable_snapshot_sha256": snapshot["sha256"],
            "result": {
                "schema": "leopard2-benchmark-v5",
                "parameters": {"attest_source": True},
                "resolved": {"profile": "legacy_high_v1", "field": "gf8",
                             "backend": "avx2"},
                "build": {"source_commit": source["head"],
                          "source_tree": source["tree"],
                          "source_tracked_dirty": False},
                "correctness": {"leopard2_round_trip": True},
            },
        }
        runner.validate_source_attestation(record, source, snapshot)
        victim = self.fixture.source / "Leopard2BackendAVX2.cpp"
        object_path = self.fixture.build / next(
            path for path in self.fixture.archive_objects
            if path.name == "Leopard2BackendAVX2.cpp.o")
        future = object_path.stat().st_mtime_ns + 2_000_000_000
        os.utime(victim, ns=(future, future))
        self.assertEqual(runner.git_identity(
            self.fixture.source, source["head"])["tracked_status"], "clean")
        expect_rejected(
            self, self.provenance, "object predates source")

    def test_avx512_cache_or_compile_flags_are_rejected(self) -> None:
        self.fixture._write_cache(
            **{"LEO2_FLAG_MAVX512F:UNINITIALIZED": "TRUE"})
        expect_rejected(self, self.provenance, "did not disable")
        self.fixture._write_cache()
        entry = next(item for item in self.fixture.entries
                     if item["file"].endswith("Leopard2BackendAVX2.cpp"))
        entry["arguments"].append("-mavx512f")
        self.fixture._write_commands()
        expect_rejected(self, self.provenance, "AVX-512")

    def test_noncanonical_generator_fails_before_recipe_inference(self) -> None:
        self.fixture._write_cache(
            **{"CMAKE_GENERATOR:INTERNAL": "Ninja"})
        expect_rejected(
            self, self.provenance,
            "CMAKE_GENERATOR is 'Ninja', expected 'Unix Makefiles'")

    def test_mixed_compile_source_and_archive_member_are_rejected(self) -> None:
        entry = next(item for item in self.fixture.entries
                     if item["file"].endswith("Leopard2Plan.cpp"))
        wrong = str((self.fixture.source / "leopard2.cpp").resolve())
        entry["file"] = wrong
        entry["arguments"][entry["arguments"].index("-c") + 1] = wrong
        self.fixture._write_commands()
        expect_rejected(self, self.provenance, "source closure differs")
        entry["file"] = str(
            (self.fixture.source / "Leopard2Plan.cpp").resolve())
        entry["arguments"][entry["arguments"].index("-c") + 1] = entry["file"]
        self.fixture._write_commands()
        extra = self.fixture.build / "unexpected.o"
        extra.write_bytes(b"unexpected\n")
        subprocess.run([
            "/usr/bin/ar", "q", str(self.fixture.build / "libleopard.a"),
            str(extra)], check=True)
        shutil.copyfile("/bin/true", self.fixture.executable)
        self.fixture.executable.chmod(0o755)
        expect_rejected(self, self.provenance, "archive members differ")

    def test_archive_member_bytes_are_bound_to_external_objects(self) -> None:
        object_path = self.fixture.build / self.fixture.archive_objects[0]
        old_status = object_path.stat()
        old_bytes = object_path.read_bytes()
        replacement = bytes([old_bytes[0] ^ 1]) + old_bytes[1:]
        object_path.write_bytes(replacement)
        os.utime(object_path, ns=(old_status.st_atime_ns,
                                  old_status.st_mtime_ns))
        expect_rejected(
            self, self.provenance, "archive member bytes differ from object")

    def test_rebuild_comparison_rejects_stale_object_archive_and_executable(self) \
            -> None:
        exact = self.provenance()
        proof = runner.compare_reproducible_builds(
            exact, copy.deepcopy(exact))
        self.assertEqual(proof["archive_sha256"],
                         exact["archive"]["sha256"])

        mutations = []
        stale_object = copy.deepcopy(exact)
        stale_object["source_object_compile_closure"][0]["object"][
            "sha256"] = "0" * 64
        mutations.append((stale_object, "object bytes differ"))
        stale_member = copy.deepcopy(exact)
        stale_member["archive_member_identities"][0]["sha256"] = "1" * 64
        mutations.append((stale_member, "archive member bytes differ"))
        stale_archive = copy.deepcopy(exact)
        stale_archive["archive"]["sha256"] = "2" * 64
        mutations.append((stale_archive, "archive bytes differ"))
        stale_executable = copy.deepcopy(exact)
        stale_executable["executable"]["sha256"] = "3" * 64
        mutations.append((stale_executable, "executable bytes differ"))
        compare = runner.compare_reproducible_builds
        for stale, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self, lambda stale=stale: compare(stale, exact), message)


if __name__ == "__main__":
    unittest.main()
