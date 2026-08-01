#!/usr/bin/env python3
"""Fail-closed source-identity tests for the all-K diagnostic runner."""

from __future__ import annotations

import copy
import errno
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import signal
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock


MODULE_PATH = Path(__file__).with_name("run_allk_gap.py")
SPEC = importlib.util.spec_from_file_location("leopard2_run_allk_gap", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load run_allk_gap.py")
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)
provenance_module = sys.modules[
    runner.candidate_build_provenance.__module__]


def expect_rejected(test: unittest.TestCase, callback, text: str) -> None:
    with test.assertRaisesRegex(RuntimeError, text):
        callback()


def wait_for_path(
        test: unittest.TestCase, path: Path, label: str,
        timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while not path.is_file() and time.monotonic() < deadline:
        time.sleep(0.01)
    test.assertTrue(path.is_file(), f"{label} did not start")


def require_lock_blocked(
        test: unittest.TestCase, descriptor: int, label: str) -> None:
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return
    fcntl.flock(descriptor, fcntl.LOCK_UN)
    test.fail(f"{label} did not retain the campaign lock")


def wait_for_lock_release(
        test: unittest.TestCase, descriptor: int, label: str,
        timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return
        except BlockingIOError:
            time.sleep(0.02)
    test.fail(f"{label} leaked the inherited campaign lock")


def initialize_git_fixture(root: Path, content: str = "one\n") -> str:
    subprocess.run(["/usr/bin/git", "init", "-q", str(root)], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                    "user.name", "Leopard2 Test"], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                    "user.email", "test@example.invalid"], check=True)
    (root / "tracked.txt").write_text(content, encoding="utf-8")
    subprocess.run(["/usr/bin/git", "-C", str(root), "add", "."], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "commit", "-qm",
                    "fixture"], check=True)
    return subprocess.check_output(
        ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
        text=True).strip()


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
        c_test_source = "tests/leopard2/test_codec_options_abi.c"
        sources.update({
            "LeopardFF8.cpp", "LeopardFF16.cpp",
            "Leopard2BackendSSSE3.cpp",
            "Leopard2BackendAVX2.cpp", "bench/leopard2/benchmark.cpp",
            c_test_source,
        })
        self.library_source_names = (
            sources - {"bench/leopard2/benchmark.cpp", c_test_source})
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
        # compile_commands.json retains the exact configured driver spelling;
        # resolving these aliases would no longer match CMakeCache.txt.
        compiler = "/usr/bin/c++"
        c_compiler = "/usr/bin/cc"
        for relative in sorted(
                sources - {"bench/leopard2/benchmark.cpp", c_test_source}):
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
        c_test_output = Path(
            "CMakeFiles/leopard2_codec_options_abi_test.dir/tests/leopard2"
        ) / "test_codec_options_abi.c.o"
        self._add_entry(c_test_source, c_test_output, c_compiler, [])
        self.c_test_object = c_test_output
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
        libgomp = subprocess.check_output(
            [compiler, "-print-file-name=libgomp.so"],
            text=True).strip()
        libpthread = subprocess.check_output(
            [compiler, "-print-file-name=libpthread.a"],
            text=True).strip()
        executable_recipe.write_text(
            compiler + " -Wall -Wextra -fopenmp "
            "-O3 -DNDEBUG -O3 " + self.benchmark_object.as_posix() +
            " -o bench_leopard2 libleopard.a " +
            libgomp + " " + libpthread + "\n", encoding="utf-8")
        self._relink()

    def close(self) -> None:
        self.temporary.cleanup()

    def _add_entry(self, relative, output, compiler, flags) -> None:
        source = (self.source / relative).resolve()
        target = self.build / output
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((relative + " object\n").encode())
        source_name = Path(relative).name
        enhanced_backend = source_name.startswith(
            ("Leopard2BackendSSSE3", "Leopard2BackendAVX2",
             "Leopard2BackendGFNI"))
        definitions = [
            "-DNDEBUG",
            "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
            "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
        ]
        if source_name in {
                "leopard2.cpp", "Leopard2Backend.cpp",
                "Leopard2BackendAVX2.cpp"}:
            definitions.append(
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
        if relative in self.library_source_names and not enhanced_backend:
            definitions.extend((
                "-DLEO2_DISABLE_AVX2_CODEGEN=1",
                "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_HAVE_SSSE3_BACKEND=1",
                "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
                "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            ))
        elif source_name.startswith("Leopard2BackendAVX2"):
            definitions.append(
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1")
        elif relative == "bench/leopard2/benchmark.cpp":
            generated_header = (
                self.build /
                "generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            definitions.extend((
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="
                f'"{generated_header}"',
                "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="
                f'"{"0" * 64}"',
            ))
        options = [
            f"-I{self.source.resolve()}", "-Wall", "-Wextra",
            *(["-fopenmp"] if enhanced_backend else
              ["-fopenmp", "-fopenmp"]),
            "-O3", "-O3", "-std=gnu++11", *flags,
        ]
        if source_name.startswith(
                ("Leopard2BackendAVX2", "Leopard2BackendGFNI")):
            options.append("-falign-functions=64")
        arguments = [
            compiler, *definitions, *options, "-o", output.as_posix(),
            "-c", str(source)]
        self.entries.append({
            "directory": str(self.build.resolve()),
            "arguments": arguments, "file": str(source),
            "output": output.as_posix(),
        })

    def _write_cache(self, **overrides) -> None:
        values = {
            "CMAKE_BUILD_TYPE:STRING": "Release",
            "CMAKE_EXPORT_COMPILE_COMMANDS:UNINITIALIZED": "ON",
            "CMAKE_GENERATOR:INTERNAL": "Unix Makefiles",
            "CMAKE_HOME_DIRECTORY:INTERNAL": str(self.source.resolve()),
            "CMAKE_C_COMPILER:FILEPATH": "/usr/bin/cc",
            "CMAKE_CXX_COMPILER:FILEPATH": "/usr/bin/c++",
            "CMAKE_CXX_FLAGS:STRING": "",
            "CMAKE_CXX_FLAGS_RELEASE:STRING": "-O3 -DNDEBUG",
            "CMAKE_AR:FILEPATH": "/usr/bin/ar",
            "CMAKE_LINKER:FILEPATH": "/usr/bin/ld",
            "CMAKE_MAKE_PROGRAM:FILEPATH": "/usr/bin/make",
            "CMAKE_RANLIB:FILEPATH": "/usr/bin/ranlib",
            "ENABLE_OPENMP:BOOL": "ON",
            "LEO2_BUILD_BENCHMARKS:BOOL": "ON",
            "LEO2_BUILD_FUZZERS:BOOL": "OFF",
            "LEO2_ENABLE_CUDA:BOOL": "OFF",
            "LEO2_BENCHMARK_GIT_EXECUTABLE:FILEPATH": "/usr/bin/git",
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:INTERNAL":
                runner.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:INTERNAL":
                "0" * 64,
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN:BOOL": "OFF",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE:BOOL": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR:BOOL": "OFF",
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE:BOOL": "ON",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT:BOOL": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING:BOOL": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING:BOOL": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING:BOOL": "ON",
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT:BOOL": "ON",
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE:STRING": "0",
            "LEOPARD_ENABLE_GF8:BOOL": "ON",
            "LEOPARD_ENABLE_GF16:BOOL": "ON",
            "LEO2_BACKEND_VARIANT:STRING": "auto",
            "LEO2_BUILD_ALLK_DIAGNOSTIC:BOOL": "OFF",
            "LEO2_BUILD_TESTS:BOOL": "ON",
            "LEO2_FLAG_FALIGN_FUNCTIONS_64:INTERNAL": "1",
            "LEO2_FLAG_MAVX2:INTERNAL": "1",
            "LEO2_FLAG_MNO_AVX512F:INTERNAL": "1",
            "LEO2_FLAG_MAVX512F:UNINITIALIZED": "FALSE",
            "LEO2_FLAG_MAVX512BW:UNINITIALIZED": "FALSE",
            "LEO2_FLAG_MAVX512VL:UNINITIALIZED": "FALSE",
            "LEO2_LOCATOR_GIT_EXECUTABLE:FILEPATH": "/usr/bin/git",
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

    @staticmethod
    def run_contract(
            schema: str = runner.RUN_CONTRACT_SCHEMA) -> dict:
        expected = runner.ALL_K_EVIDENCE_CONTRACTS[schema]
        cache_keys = {
            runner.RUN_CONTRACT_SCHEMA_V4:
                runner.ALL_K_BUILD_CACHE_KEYS_V2,
            runner.RUN_CONTRACT_SCHEMA_V5:
                runner.ALL_K_BUILD_CACHE_KEYS_V3,
            runner.RUN_CONTRACT_SCHEMA:
                runner.ALL_K_BUILD_CACHE_KEYS,
        }[schema]
        cache = {key: "fixture" for key in cache_keys}
        cache.update({
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                expected["configuration"],
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "1" * 64,
            "LEO2_BUILD_ALLK_DIAGNOSTIC": "ON",
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
        })
        if schema == runner.RUN_CONTRACT_SCHEMA_V5:
            cache["LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = "OFF"
        elif schema == runner.RUN_CONTRACT_SCHEMA:
            cache.update({
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
                "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
                "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            })
        return {
            "schema": schema,
            "main_commit": "a" * 40,
            "current_commit": "b" * 40,
            "expected_main_sha256": "c" * 64,
            "current_source_initial": {"head": "b" * 40},
            "current_build_initial": {
                "schema": expected["closure"],
                "validated_cache": cache,
            },
            "current_reproducible_build": {
                "schema": expected["proof"],
                "immutable_replay": {
                    "schema": "leopard2-immutable-replay-tools/v2",
                    "recipe_transport": {
                        "schema": expected["replay_plan"],
                    },
                    "invocations": {
                        name: {"schema": expected["replay_invocation"]}
                        for name in ("configure", "build")
                    },
                },
            },
            "main_executable_initial": {"sha256": "d" * 64},
            "current_executable_initial": {"sha256": "e" * 64},
            "main_executable_snapshot": {"sha256": "d" * 64},
            "current_executable_snapshot": {"sha256": "e" * 64},
            "current_source_attestation_identity": {"sha256": "f" * 64},
            "main_isa": {"profile": "pure-avx2"},
            "current_isa": {"profile": "pure-avx2"},
            "benchmark_lock": {"path": "/tmp/fixture.lock"},
            "allowed_cpus": [0, 1, 2],
            "used_cpus": [0],
            "workers": 1,
            "order": list(runner.ORDER),
            "timeout_seconds": 30.0,
            "with_current_legacy": False,
            "matrix": {"cell_count": 1},
            "measurement_note": "fixture",
        }

    @staticmethod
    def bind_raw_output(record):
        payload = json.dumps(
            record["result"], sort_keys=True, separators=(",", ":"),
            allow_nan=False) + "\n"
        record.update({
            "command": record.get("command", ["fixture-benchmark"]),
            "duration_ns": record.get("duration_ns", 1),
            "stdout": payload,
            "stdout_sha256": runner.hashlib.sha256(
                payload.encode("utf-8")).hexdigest(),
            "stderr": "",
        })
        return record

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
            records.append(self.bind_raw_output({
                "role": role,
                "executable_snapshot_sha256": snapshot["sha256"],
                "returncode": 0,
                "result": result,
            }))
        return records

    def attestation_record(self):
        return self.bind_raw_output({
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
                "metrics": {
                    "codec_setup": {"median_us": 1.0},
                    "encode_execution": {
                        "median_us_per_batch_call": 2.0,
                    },
                },
            },
        })

    def correctness_records(self, cell):
        parent, padded = runner.parent_for(cell.k, cell.r)
        field = "gf8" if cell.family == "gf8-all-k" else "gf16"
        digest = {
            "algorithm": "fnv1a64",
            "original_data": "0123456789abcdef",
            "transmitted_parity": "123456789abcdef0",
            # Recovery fingerprints only the missing-original subset.  Keep
            # this distinct from the all-original digest so partial-loss test
            # fixtures exercise the real benchmark contract.
            "recovered_originals": "fedcba9876543210",
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
                "missing_original_indices": runner.expected_losses(cell),
            }
            if record["role"] == "leopard2":
                record["result"]["parameters"].update({
                    "requested_profile": "legacy_high_v1",
                    "requested_field":
                        "gf8" if cell.family == "gf8-all-k" else "gf16",
                    "requested_backend": "avx2",
                    "force_generic_decode": False,
                    "force_specialized_decode": False,
                    "force_tiled_decode": False,
                    "force_materialized_decode": False,
                    "skip_legacy": False,
                    "retain_samples": True,
                    "report_decode_path": True,
                })
            record["result"]["correctness"] = (
                {"round_trip": True} if record["role"] == "main" else
                {"leopard2_round_trip": True})
            record["result"]["resolved"] = {
                "profile": "legacy_high_v1", "field": field,
                "thread_count": 1, "parent_count": parent,
                "padded_side": padded,
            }
            record["result"]["workload_digests"] = copy.deepcopy(digest)
        return records

    def analyzed_cell(self, cell, main, current, cpu, contract_digest):
        records = self.correctness_records(cell)
        for record in records:
            executable = main if record["role"] == "main" else current
            record["command"] = runner.command(
                record["role"], executable, cell, cpu, False)
            if record["role"] == "main":
                record["result"]["metrics"] = {
                    "encode_execution": {
                        "median_us_per_batch_call": 12.0,
                    },
                    "decode_including_setup": {
                        "median_us_per_batch_call": 18.0,
                    },
                }
            else:
                record["result"]["parameters"]["report_decode_path"] = True
                record["result"]["parameters"]["skip_legacy"] = True
                record["result"]["resolved"] = {
                    "profile": "legacy_high_v1",
                    "field": "gf8",
                    "backend": "avx2",
                    "thread_count": 1,
                    "parent_count": 8,
                    "padded_side": 2,
                    "selected_decode_path": "direct",
                    "selected_decode_rule": "direct",
                    "decode_required_work_slots": 0,
                    "decode_aligned_prefix_bytes": cell.shard_bytes,
                    "decode_tail_bytes": 0,
                    "decode_rounded_bytes": cell.shard_bytes,
                    "decode_multi_item_batch": False,
                }
                record["result"]["metrics"] = {
                    "codec_setup": {"median_us": 1.0},
                    "encode_execution": {
                        "median_us_per_batch_call": 6.0,
                    },
                    "decode_plan_setup": {"median_us": 2.0},
                    "decode_execution": {
                        "median_us_per_batch_call": 3.0,
                    },
                }
                record["result"]["memory"] = {
                    "decode_scratch_bytes_per_stripe": 0,
                }
                record["result"]["legacy"] = {"available": False}
            self.bind_raw_output(record)
        return runner.analyze_cell(
            cell, records, cpu, contract_digest, self.main_commit,
            self.current_source, self.main_snapshot, self.current_snapshot)

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

    def test_cpu_lane_scheduler_serializes_unequal_jobs_per_cpu(self) -> None:
        guard = threading.Lock()
        active_by_cpu = {11: 0, 12: 0}
        maximum_by_cpu = {11: 0, 12: 0}
        total_active = 0
        maximum_total = 0
        first_job_barrier = threading.Barrier(2)

        def execute(item, cpu):
            nonlocal total_active, maximum_total
            if item in (0, 1):
                first_job_barrier.wait(timeout=2.0)
            with guard:
                active_by_cpu[cpu] += 1
                total_active += 1
                maximum_by_cpu[cpu] = max(
                    maximum_by_cpu[cpu], active_by_cpu[cpu])
                maximum_total = max(maximum_total, total_active)
            try:
                time.sleep(0.08 if item == 0 else 0.005)
                return {"item": item, "cpu": cpu}
            finally:
                with guard:
                    active_by_cpu[cpu] -= 1
                    total_active -= 1

        results = runner.run_serial_cpu_lanes(
            list(range(8)), [11, 12], execute)
        self.assertEqual(
            results,
            [{"item": item, "cpu": 11 + item % 2}
             for item in range(8)])
        self.assertEqual(maximum_by_cpu, {11: 1, 12: 1})
        self.assertGreaterEqual(maximum_total, 2)

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
            "backend": "avx2", "thread_count": 1, "parent_count": 256,
            "padded_side": 64, "selected_decode_path": "tiled",
            "selected_decode_rule": "workspace_tiled",
            "decode_required_work_slots": 132,
            "decode_aligned_prefix_bytes": 4096,
            "decode_tail_bytes": 0,
            "decode_rounded_bytes": 4096,
            "decode_multi_item_batch": False,
        }}
        selected = runner.classify_paths(cell, result)
        self.assertEqual(selected["decode_path"], "tiled")
        self.assertEqual(selected["decode_rule"], "workspace_tiled")
        self.assertEqual(selected["decode_required_work_slots"], 132)
        del result["resolved"]["selected_decode_path"]
        expect_rejected(
            self, lambda: runner.classify_paths(cell, result),
            "selected decode path")

        def classified_mutation(name, value, message):
            mutated = copy.deepcopy(result)
            mutated["resolved"]["selected_decode_path"] = "tiled"
            mutated["resolved"][name] = value
            expect_rejected(
                self, lambda: runner.classify_paths(cell, mutated), message)

        classified_mutation(
            "parent_count", 256.0, "parent geometry")
        classified_mutation(
            "decode_required_work_slots", 256, "decode workspace")
        classified_mutation(
            "decode_multi_item_batch", True, "batch geometry")
        classified_mutation(
            "decode_rounded_bytes", 4160, "byte geometry")
        classified_mutation(
            "selected_decode_rule", "measured_batch_tiled",
            "incoherent decode path/rule")

        translated_cell = runner.Cell(
            "translated", "gf8-all-k", 4, 4, 65, 2,
            "max-GF8-R", "max-loss", 9, 5, 16, 1)
        for path in ("materialized", "tiled"):
            translated = {"resolved": {
                "profile": "legacy_high_v1", "field": "gf8",
                "backend": "avx2", "thread_count": 1,
                "parent_count": 8, "padded_side": 4,
                "selected_decode_path": path,
                "selected_decode_rule": "translated_low",
                "decode_required_work_slots": 8,
                "decode_aligned_prefix_bytes": 64,
                "decode_tail_bytes": 1,
                "decode_rounded_bytes": 128,
                "decode_multi_item_batch": False,
            }}
            with self.subTest(translated_path=path):
                self.assertEqual(
                    runner.classify_paths(
                        translated_cell, translated)["decode_rule"],
                    "translated_low")

        nondual = copy.deepcopy(translated)
        nondual_cell = runner.Cell(
            "nondual", "gf8-all-k", 5, 4, 65, 2,
            "max-GF8-R", "max-loss", 9, 5, 16, 1)
        nondual["resolved"]["parent_count"] = 16
        expect_rejected(
            self, lambda: runner.classify_paths(nondual_cell, nondual),
            "incoherent decode path/rule")

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

    def test_invocation_records_bind_exact_raw_stdout(self) -> None:
        record = self.attestation_record()
        runner.validate_invocation_record(record, "fixture invocation")

        mutations = []
        missing_stdout = copy.deepcopy(record)
        del missing_stdout["stdout"]
        mutations.append((missing_stdout, "keys differ"))
        missing_digest = copy.deepcopy(record)
        del missing_digest["stdout_sha256"]
        mutations.append((missing_digest, "keys differ"))
        changed_result = copy.deepcopy(record)
        changed_result["result"]["build"]["source_commit"] = "d" * 40
        mutations.append((changed_result, "parsed result differs"))
        changed_stdout = copy.deepcopy(record)
        changed_stdout["stdout"] += " "
        mutations.append((changed_stdout, "SHA-256 differs"))
        for mutation, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda mutation=mutation: runner.validate_invocation_record(
                        mutation, "fixture invocation"),
                    message)

        for raw_value, retained_value in ((1.0, 1), (1, True)):
            numeric = copy.deepcopy(record)
            numeric["result"] = {"value": retained_value}
            numeric["stdout"] = json.dumps(
                {"value": raw_value}, sort_keys=True,
                separators=(",", ":")) + "\n"
            numeric["stdout_sha256"] = runner.hashlib.sha256(
                numeric["stdout"].encode("utf-8")).hexdigest()
            with self.subTest(
                    raw_value=raw_value, retained_value=retained_value):
                expect_rejected(
                    self,
                    lambda numeric=numeric:
                        runner.validate_invocation_record(
                            numeric, "numeric fixture invocation"),
                    "parsed result differs")

        failed = copy.deepcopy(record)
        failed["returncode"] = 1
        del failed["result"]
        failed["stdout"] = ""
        failed["stdout_sha256"] = runner.hashlib.sha256(b"").hexdigest()
        failed["stderr"] = "expected failure"
        runner.validate_invocation_record(failed, "failed fixture invocation")
        failed["unexpected"] = True
        expect_rejected(
            self,
            lambda: runner.validate_invocation_record(
                failed, "failed fixture invocation"),
            "keys differ")

        for payload, message in (
                (b'{"value":NaN}\\n', "non-finite"),
                (b'{"value":Infinity}\\n', "non-finite"),
                (b'{"value":1e9999}\\n', "non-finite"),
                (b'{"value":-1e9999}\\n', "non-finite"),
                (b'{"value":1,"value":1}\\n', "duplicate"),
                (b'\\xff', "strict finite JSON")):
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda payload=payload: runner.strict_json_bytes(
                        payload, "fixture JSON"),
                    message)
        deeply_nested = b"[" * 10000 + b"0" + b"]" * 10000
        expect_rejected(
            self,
            lambda: runner.strict_json_bytes(
                deeply_nested, "deeply nested fixture JSON"),
            "strict finite JSON")

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

    def test_direct_script_help_loads_shared_git_capture(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-direct-script-") as directory:
            completed = subprocess.run(
                [sys.executable, str(MODULE_PATH), "--help"],
                cwd=directory, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={**os.environ, "PYTHONWARNINGS":
                     "error::ResourceWarning"},
                timeout=30, check=False)
        self.assertEqual(
            completed.returncode, 0,
            completed.stderr.decode("utf-8", errors="replace"))

    def test_git_identity_accepts_linked_worktree_and_superproject(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-git-layouts-") as directory:
            root = Path(directory)
            repository = root / "repository"
            worktree = root / "worktree"
            submodule_source = root / "submodule-source"
            repository.mkdir()
            submodule_source.mkdir()
            commit = initialize_git_fixture(repository, "parent\n")
            subprocess.run(
                ["/usr/bin/git", "-C", str(repository), "worktree", "add",
                 "-q", "-b", "fixture-worktree", str(worktree)],
                check=True)
            worktree_identity = runner.git_identity(worktree, commit)
            self.assertEqual(
                worktree_identity["git_metadata"]["layout"],
                "linked-worktree")

            initialize_git_fixture(submodule_source, "submodule\n")
            subprocess.run(
                ["/usr/bin/git", "-c", "protocol.file.allow=always",
                 "-C", str(repository), "submodule", "add", "-q",
                 str(submodule_source), "submodule"], check=True)
            subprocess.run(
                ["/usr/bin/git", "-C", str(repository), "commit", "-qam",
                 "add submodule"], check=True)
            super_commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(repository), "rev-parse", "HEAD"],
                text=True).strip()
            super_identity = runner.git_identity(repository, super_commit)
            self.assertEqual(len(super_identity["submodules"]), 1)
            nested = super_identity["submodules"][0]
            self.assertEqual(nested["path"], "submodule")
            self.assertEqual(
                nested["object_id"], nested["identity"]["head"])

    def test_git_identity_rejects_deterministic_mixed_command_state(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-git-mixed-") as directory:
            root = Path(directory)
            first = initialize_git_fixture(root)
            (root / "tracked.txt").write_text("two\n", encoding="utf-8")
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "commit", "-qam",
                 "second"], check=True)
            alternate_tree = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse",
                 "HEAD^{tree}"], text=True).strip()
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "reset", "--hard", "-q",
                 first], check=True)
            invoke = runner.git_capture._invoke_git

            def mixed(*args, **kwargs):
                if tuple(args[3]) == (
                        "rev-parse", "--verify", "HEAD^{tree}"):
                    return (alternate_tree + "\n").encode("ascii")
                return invoke(*args, **kwargs)

            with mock.patch.object(
                    runner.git_capture, "_invoke_git", side_effect=mixed):
                expect_rejected(
                    self, lambda: runner.git_identity(root, first),
                    "commit object names a different root tree")

    def test_git_identity_rejects_actual_worktree_gitfile_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-gitfile-aba-") as directory:
            root = Path(directory)
            repository = root / "repository"
            worktree = root / "worktree"
            repository.mkdir()
            commit = initialize_git_fixture(repository)
            subprocess.run(
                ["/usr/bin/git", "-C", str(repository), "worktree", "add",
                 "-q", "-b", "aba-worktree", str(worktree)], check=True)
            gitfile = worktree / ".git"
            alternate = root / "alternate.gitfile"
            saved = root / "saved.gitfile"
            shutil.copy2(gitfile, alternate)
            invoke = runner.git_capture._invoke_git
            triggered = False

            def aba(*args, **kwargs):
                nonlocal triggered
                output = invoke(*args, **kwargs)
                if not triggered:
                    triggered = True
                    gitfile.rename(saved)
                    alternate.rename(gitfile)
                    gitfile.rename(alternate)
                    saved.rename(gitfile)
                return output

            with mock.patch.object(
                    runner.git_capture, "_invoke_git", side_effect=aba):
                expect_rejected(
                    self, lambda: runner.git_identity(worktree, commit),
                    "changed")
            self.assertTrue(triggered)

    def test_git_identity_propagates_campaign_lock_to_every_git_child(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-git-lock-") as directory:
            root = Path(directory)
            commit = initialize_git_fixture(root)
            lock = root.parent / "campaign.lock"
            descriptor = os.open(
                lock, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o600)
            observed = 0
            run = runner.git_capture._run

            def retained(*args, **kwargs):
                nonlocal observed
                self.assertIn(
                    descriptor, kwargs.get("inherited_descriptors", ()))
                observed += 1
                return run(*args, **kwargs)

            try:
                with mock.patch.object(
                        runner.git_capture, "_run", side_effect=retained):
                    runner.git_identity(
                        root, commit,
                        inherited_descriptors=(descriptor,))
            finally:
                os.close(descriptor)
            self.assertGreater(observed, 1)

    def test_git_identity_rejects_symlinked_and_hardlink_mutated_metadata(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-git-metadata-") as directory:
            parent = Path(directory)
            root = parent / "source"
            root.mkdir()
            commit = initialize_git_fixture(root)
            config = root / ".git/config"
            external_config = parent / "external.config"
            config.rename(external_config)
            config.symlink_to(external_config)
            expect_rejected(
                self, lambda: runner.git_identity(root, commit),
                "not a regular file")
            config.unlink()
            external_config.rename(config)

            head = root / ".git/HEAD"
            head_link = parent / "HEAD.hardlink"
            os.link(head, head_link)
            original = head.read_bytes()
            invoke = runner.git_capture._invoke_git
            triggered = False

            def hardlink_aba(*args, **kwargs):
                nonlocal triggered
                output = invoke(*args, **kwargs)
                if not triggered:
                    triggered = True
                    descriptor = os.open(
                        head_link, os.O_WRONLY | os.O_TRUNC | os.O_CLOEXEC)
                    try:
                        os.write(descriptor, b"x" * len(original))
                        os.fsync(descriptor)
                        os.lseek(descriptor, 0, os.SEEK_SET)
                        os.ftruncate(descriptor, 0)
                        os.write(descriptor, original)
                        os.fsync(descriptor)
                    finally:
                        os.close(descriptor)
                return output

            with mock.patch.object(
                    runner.git_capture, "_invoke_git",
                    side_effect=hardlink_aba):
                expect_rejected(
                    self, lambda: runner.git_identity(root, commit),
                    "changed")
            self.assertTrue(triggered)

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
        wrong_parameter_type = self.correctness_records(cell)
        wrong_parameter_type[0]["result"]["parameters"]["K"] = True
        mutations.append((wrong_parameter_type, "parameter K"))
        wrong_resolved_type = self.correctness_records(cell)
        wrong_resolved_type[0]["result"]["resolved"]["parent_count"] = 8.0
        mutations.append((wrong_resolved_type, "resolved geometry"))
        float_batch = self.correctness_records(cell)
        float_batch[1]["result"]["parameters"]["batch"] = 1.0
        mutations.append((float_batch, "parameter batch"))
        wrong_losses = self.correctness_records(cell)
        wrong_losses[2]["result"]["parameters"][
            "missing_original_indices"] = [2]
        mutations.append((wrong_losses, "missing-original indices"))
        coercible_losses = self.correctness_records(cell)
        coercible_losses[3]["result"]["parameters"][
            "missing_original_indices"] = [
                float(value) for value in runner.expected_losses(cell)]
        mutations.append((coercible_losses, "missing-original indices"))
        one_sided_bad_recovery = self.correctness_records(cell)
        one_sided_bad_recovery[2]["result"]["workload_digests"][
            "recovered_originals"] = "1111111111111111"
        mutations.append(
            (one_sided_bad_recovery,
             "Leopard1 and Leopard2 workload digests differ"))
        wrong_option_type = self.correctness_records(cell)
        wrong_option_type[1]["result"]["parameters"]["skip_legacy"] = 0
        mutations.append((wrong_option_type, "option skip_legacy"))
        for records, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda records=records: runner.validate_correctness(
                        cell, records),
                    message)

        full_loss = runner.Cell("full", "gf8-all-k", 3, 3, 64, 3,
                                "max-GF8-R", "max-loss", 8, 1, 1, 0)
        inconsistent_full_loss = self.correctness_records(full_loss)
        expect_rejected(
            self,
            lambda: runner.validate_correctness(
                full_loss, inconsistent_full_loss),
            "full-loss recovered-original digest")
        for record in inconsistent_full_loss:
            digests = record["result"]["workload_digests"]
            digests["recovered_originals"] = digests["original_data"]
        runner.validate_correctness(full_loss, inconsistent_full_loss)

    def test_manifest_rejects_stale_or_resigned_contract(self) -> None:
        proof_validator = mock.patch.object(
            runner, "validate_reproducible_build_proof")
        validated_proof = proof_validator.start()
        self.addCleanup(proof_validator.stop)
        cell = runner.Cell("one", "gf8-all-k", 1, 1, 64, 1,
                           "low-R", "one-loss", 1, 1, 1, 0)
        contract = self.run_contract()
        digest = runner.canonical_digest(contract)
        preflight = self.attestation_record()
        attestation_identity = runner.source_attestation_identity(
            preflight, self.current_source, self.current_snapshot)
        manifest = {
            "schema": runner.MANIFEST_SCHEMA,
            "run_contract": contract,
            "run_contract_sha256": digest,
            "cells": [runner.dataclasses.asdict(cell)],
            "current_source_attestation_preflights": [preflight],
            "completion": None,
        }
        runner.validate_manifest(
            manifest, contract, digest, [cell], self.current_source,
            self.current_snapshot, attestation_identity)

        stale = copy.deepcopy(manifest)
        stale["run_contract"]["main"] = "d" * 40
        stale["run_contract_sha256"] = runner.canonical_digest(
            stale["run_contract"])
        expect_rejected(
            self, lambda: runner.validate_manifest(
                stale, contract, digest, [cell], self.current_source,
                self.current_snapshot, attestation_identity),
            "run contract")

        corrupted = copy.deepcopy(manifest)
        corrupted["run_contract_sha256"] = "e" * 64
        expect_rejected(
            self,
            lambda: runner.validate_manifest(
                corrupted, contract, digest, [cell], self.current_source,
                self.current_snapshot, attestation_identity),
            "digest mismatch")

        historical = copy.deepcopy(manifest)
        historical["schema"] = runner.MANIFEST_SCHEMA_V4
        expect_rejected(
            self,
            lambda: runner.validate_manifest(
                historical, contract, digest, [cell], self.current_source,
                self.current_snapshot, attestation_identity),
            "schema tuple")

        coercible_cell = copy.deepcopy(manifest)
        coercible_cell["cells"][0]["k"] = True
        expect_rejected(
            self,
            lambda: runner.validate_manifest(
                coercible_cell, contract, digest, [cell],
                self.current_source, self.current_snapshot,
                attestation_identity),
            "cell matrix")
        self.assertGreater(validated_proof.call_count, 0)

    def test_all_k_outer_schema_binds_exact_replay_generation(self) -> None:
        proof_validator = mock.patch.object(
            runner, "validate_reproducible_build_proof")
        validated_proof = proof_validator.start()
        self.addCleanup(proof_validator.stop)
        current = self.run_contract()
        v5 = self.run_contract(runner.RUN_CONTRACT_SCHEMA_V5)
        v4 = self.run_contract(runner.RUN_CONTRACT_SCHEMA_V4)
        self.assertIs(
            runner.validate_run_contract_evidence(current), current)
        self.assertIs(
            runner.validate_run_contract_evidence(v5), v5)
        self.assertIs(
            runner.validate_run_contract_evidence(v4), v4)

        # Each body remains coherent under its own generation.  Relabeling only
        # the outer contract cannot upgrade or downgrade its nested closure.
        for body, schema in (
                (v4, runner.RUN_CONTRACT_SCHEMA_V5),
                (v4, runner.RUN_CONTRACT_SCHEMA),
                (v5, runner.RUN_CONTRACT_SCHEMA_V4),
                (v5, runner.RUN_CONTRACT_SCHEMA),
                (current, runner.RUN_CONTRACT_SCHEMA_V4),
                (current, runner.RUN_CONTRACT_SCHEMA_V5)):
            with self.subTest(body=body["schema"], label=schema):
                relabeled = copy.deepcopy(body)
                relabeled["schema"] = schema
                expect_rejected(
                    self,
                    lambda relabeled=relabeled:
                        runner.validate_run_contract_evidence(relabeled),
                    "schema tuple")

        extended = copy.deepcopy(v4)
        extended["current_build_initial"]["validated_cache"][
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = "OFF"
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(extended),
            "schema tuple")

        for selector in (
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"):
            with self.subTest(v5_extra=selector):
                extended_v5 = copy.deepcopy(v5)
                extended_v5["current_build_initial"]["validated_cache"][
                    selector] = "ON"
                expect_rejected(
                    self,
                    lambda extended_v5=extended_v5:
                        runner.validate_run_contract_evidence(extended_v5),
                    "schema tuple")

        for label, variable, value in (
            ("direct-source", "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN", "ON"),
            ("high-direct", "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE", "ON"),
            ("small-direct", "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE", "3"),
            ("general-one-loss",
             "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT", "OFF"),
            ("one-shot-equal-rounded",
             "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT", "OFF"),
            ("Cauchy-log-reuse",
             "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE", "OFF"),
            ("T8-partial", "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING", "OFF"),
            ("T8-two-block",
             "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING", "OFF"),
            ("T8-ragged", "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING", "OFF"),
            ("T8-disable", "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR", "ON"),
        ):
            with self.subTest(selector=label):
                changed = copy.deepcopy(current)
                changed["current_build_initial"]["validated_cache"][
                    variable] = value
                expect_rejected(
                    self,
                    lambda changed=changed:
                        runner.validate_run_contract_evidence(changed),
                    "selector tuple")

        missing_selector = copy.deepcopy(current)
        missing_selector["current_build_initial"]["validated_cache"].pop(
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE")
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(missing_selector),
            "schema tuple")

        extended_cache = copy.deepcopy(current)
        extended_cache["current_build_initial"]["validated_cache"][
            "LEO2_UNVERSIONED_SELECTOR"] = "OFF"
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(extended_cache),
            "schema tuple")

        wrong_plan = copy.deepcopy(current)
        wrong_plan["current_reproducible_build"]["immutable_replay"][
            "recipe_transport"]["schema"] = \
                runner.LEGACY_REPLAY_RECIPE_SCHEMA
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(wrong_plan),
            "schema tuple")

        wrong_invocation = copy.deepcopy(current)
        wrong_invocation["current_reproducible_build"]["immutable_replay"][
            "invocations"]["configure"]["schema"] = \
                runner.REPLAY_INVOCATION_SCHEMA_V1
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(wrong_invocation),
            "schema tuple")

        validated_proof.side_effect = RuntimeError("fixture proof drift")
        expect_rejected(
            self,
            lambda: runner.validate_run_contract_evidence(current),
            "reproducible-build proof is invalid")
        validated_proof.side_effect = None

        cell = runner.Cell(
            "one", "gf8-all-k", 1, 1, 64, 1,
            "low-R", "one-loss", 1, 1, 1, 0)
        preflight = self.attestation_record()
        attestation_identity = runner.source_attestation_identity(
            preflight, self.current_source, self.current_snapshot)
        v5_digest = runner.canonical_digest(v5)
        v5_manifest = {
            "schema": runner.MANIFEST_SCHEMA_V5,
            "run_contract": v5,
            "run_contract_sha256": v5_digest,
            "cells": [runner.dataclasses.asdict(cell)],
            "current_source_attestation_preflights": [preflight],
            "completion": None,
        }
        runner.validate_manifest(
            v5_manifest, v5, v5_digest, [cell],
            self.current_source, self.current_snapshot,
            attestation_identity)

        relabeled_manifest = copy.deepcopy(v5_manifest)
        relabeled_manifest["schema"] = runner.MANIFEST_SCHEMA
        expect_rejected(
            self,
            lambda: runner.validate_manifest(
                relabeled_manifest, v5, v5_digest, [cell],
                self.current_source, self.current_snapshot,
                attestation_identity),
            "schema tuple")
        self.assertGreaterEqual(validated_proof.call_count, 4)

    def test_attestation_contract_identity_excludes_only_volatile_timing(
            self) -> None:
        first = self.attestation_record()
        first_identity = runner.source_attestation_identity(
            first, self.current_source, self.current_snapshot)

        second = copy.deepcopy(first)
        second["duration_ns"] = 987654321
        second["result"]["metrics"]["codec_setup"]["median_us"] = 1234.5
        second["result"]["metrics"]["encode_execution"][
            "median_us_per_batch_call"] = 6789.0
        self.bind_raw_output(second)
        second_identity = runner.source_attestation_identity(
            second, self.current_source, self.current_snapshot)
        self.assertEqual(first_identity, second_identity)
        self.assertEqual(
            runner.canonical_digest({"attestation": first_identity}),
            runner.canonical_digest({"attestation": second_identity}))

        identity_change = copy.deepcopy(second)
        identity_change["result"]["build"]["compiler"] = "different compiler"
        self.bind_raw_output(identity_change)
        changed_identity = runner.source_attestation_identity(
            identity_change, self.current_source, self.current_snapshot)
        self.assertNotEqual(first_identity, changed_identity)

    def test_stored_cell_resumes_across_new_attestation_timings(self) -> None:
        first = self.attestation_record()
        second = copy.deepcopy(first)
        second["duration_ns"] = first["duration_ns"] + 1000
        second["result"]["metrics"]["codec_setup"]["median_us"] = 9.0
        self.bind_raw_output(second)
        first_identity = runner.source_attestation_identity(
            first, self.current_source, self.current_snapshot)
        second_identity = runner.source_attestation_identity(
            second, self.current_source, self.current_snapshot)
        self.assertEqual(first_identity, second_identity)
        digest = runner.canonical_digest({
            "schema": "fixture-contract/v1",
            "source_attestation_identity": second_identity,
        })
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-resume-") as directory:
            output = Path(directory)
            path = output / "cells" / "one.json"
            runner.write_json_atomic(
                path, self.analyzed_cell(cell, main, current, 2, digest))
            resumed = runner.run_cell(
                cell, 2, main, current, output, 1.0, False, digest,
                self.main_commit, self.current_source,
                0, 0, self.main_snapshot, 0, self.current_snapshot)
            self.assertEqual(resumed, runner.load_strict_json(
                path, "resumed fixture cell",
                runner.MAX_PERSISTED_CELL_BYTES))

    def test_stored_cell_requires_exact_recomputation(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        cpu = 2
        contract_digest = "f" * 64
        document = self.analyzed_cell(
            cell, main, current, cpu, contract_digest)
        runner.validate_cell_document(
            document, cell, cpu, contract_digest, main, current, False,
            self.main_commit, self.current_source,
            self.main_snapshot, self.current_snapshot)

        mutations = []
        fabricated_speedup = copy.deepcopy(document)
        fabricated_speedup["speedup_main_over_leopard2"]["encode"] = 999999.0
        mutations.append((fabricated_speedup, "exact recomputation"))
        omitted_derived_field = copy.deepcopy(document)
        del omitted_derived_field["selected"]
        mutations.append((omitted_derived_field, "keys differ"))
        unexpected_field = copy.deepcopy(document)
        unexpected_field["accepted_without_recomputation"] = True
        mutations.append((unexpected_field, "keys differ"))
        stale_derived_data = copy.deepcopy(document)
        stale_derived_data["invocations"][1]["result"]["metrics"][
            "decode_execution"]["median_us_per_batch_call"] = 3000.0
        mutations.append((stale_derived_data, "parsed result differs"))

        coherently_changed = copy.deepcopy(document)
        changed_invocation = coherently_changed["invocations"][1]
        changed_invocation["result"]["metrics"]["decode_execution"][
            "median_us_per_batch_call"] = 3000.0
        self.bind_raw_output(changed_invocation)
        recomputed_after_tamper = runner.analyze_cell(
            cell, coherently_changed["invocations"], cpu, contract_digest,
            self.main_commit, self.current_source,
            self.main_snapshot, self.current_snapshot)
        recomputed_after_tamper["invocations"][1]["stdout"] = \
            document["invocations"][1]["stdout"]
        recomputed_after_tamper["invocations"][1]["stdout_sha256"] = \
            document["invocations"][1]["stdout_sha256"]
        mutations.append((recomputed_after_tamper, "parsed result differs"))
        for mutation, message in mutations:
            with self.subTest(message=message):
                expect_rejected(
                    self,
                    lambda mutation=mutation: runner.validate_cell_document(
                        mutation, cell, cpu, contract_digest, main, current,
                        False, self.main_commit, self.current_source,
                        self.main_snapshot, self.current_snapshot),
                    message)

    def test_derived_speedup_overflow_is_rejected(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        document = self.analyzed_cell(
            cell, main, current, 2, "f" * 64)
        for record in document["invocations"]:
            if record["role"] == "main":
                record["result"]["metrics"]["encode_execution"][
                    "median_us_per_batch_call"] = 1.7976931348623157e308
            else:
                record["result"]["metrics"]["encode_execution"][
                    "median_us_per_batch_call"] = 5e-324
            self.bind_raw_output(record)
        expect_rejected(
            self,
            lambda: runner.analyze_cell(
                cell, document["invocations"], 2, "f" * 64,
                self.main_commit, self.current_source,
                self.main_snapshot, self.current_snapshot),
            "encode speedup is not finite")

    def test_metric_requires_exact_finite_positive_json_number(self) -> None:
        path = ("metrics", "encode_execution",
                "median_us_per_batch_call")
        for case_name, value, message in (
                ("boolean", True, "exact JSON number"),
                ("numeric string", "1.0", "exact JSON number"),
                ("zero", 0, "finite and positive"),
                ("negative", -1.0, "finite and positive"),
                ("NaN", float("nan"), "finite and positive"),
                ("infinity", float("inf"), "finite and positive"),
                ("huge integer", 10 ** 10000, "finite float range")):
            with self.subTest(case=case_name):
                record = {"result": {"metrics": {
                    "encode_execution": {
                        "median_us_per_batch_call": value,
                    },
                }}}
                expect_rejected(
                    self, lambda record=record: runner.metric(record, path),
                    message)

    def test_analyze_rejects_coercible_or_nonpositive_consumed_metrics(
            self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        base = self.analyzed_cell(cell, main, current, 2, "f" * 64)
        mutations = (
            (0, ("metrics", "encode_execution",
                 "median_us_per_batch_call"), True, "exact JSON number"),
            (1, ("metrics", "decode_plan_setup", "median_us"),
             "2.0", "exact JSON number"),
            (3, ("metrics", "decode_including_setup",
                 "median_us_per_batch_call"), 0, "finite and positive"),
            (2, ("metrics", "codec_setup", "median_us"),
             -1.0, "finite and positive"),
        )
        for record_index, path, value, message in mutations:
            with self.subTest(record_index=record_index, path=path):
                document = copy.deepcopy(base)
                record = document["invocations"][record_index]
                target = record["result"]
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                self.bind_raw_output(record)
                expect_rejected(
                    self,
                    lambda document=document: runner.analyze_cell(
                        cell, document["invocations"], 2, "f" * 64,
                        self.main_commit, self.current_source,
                        self.main_snapshot, self.current_snapshot),
                    message)

        legacy = copy.deepcopy(base)
        for record in legacy["invocations"]:
            if record["role"] == "leopard2":
                record["command"].remove("--skip-legacy")
                record["result"]["parameters"]["skip_legacy"] = False
                record["result"]["legacy"] = {
                    "available": True,
                    "encode_execution": {
                        "median_us_per_batch_call": True,
                    },
                    "decode_including_setup": {
                        "median_us_per_batch_call": 4.0,
                    },
                }
                self.bind_raw_output(record)
        expect_rejected(
            self,
            lambda: runner.analyze_cell(
                cell, legacy["invocations"], 2, "f" * 64,
                self.main_commit, self.current_source,
                self.main_snapshot, self.current_snapshot),
            "exact JSON number")

    def test_analyze_requires_exact_command_bound_legacy_availability(
            self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        base = self.analyzed_cell(cell, main, current, 2, "f" * 64)

        truthy_text = copy.deepcopy(base["invocations"])
        for record in truthy_text:
            if record["role"] == "leopard2":
                record["result"]["legacy"]["available"] = "false"
                self.bind_raw_output(record)
        expect_rejected(
            self,
            lambda: runner.analyze_cell(
                cell, truthy_text, 2, "f" * 64, self.main_commit,
                self.current_source, self.main_snapshot,
                self.current_snapshot),
            "not an exact boolean")

        enabled = copy.deepcopy(base["invocations"])
        for record in enabled:
            if record["role"] == "leopard2":
                record["command"].remove("--skip-legacy")
                record["result"]["parameters"]["skip_legacy"] = False
                record["result"]["legacy"] = {
                    "available": True,
                    "encode_execution": {
                        "median_us_per_batch_call": 7.0,
                    },
                    "decode_including_setup": {
                        "median_us_per_batch_call": 8.0,
                    },
                }
                self.bind_raw_output(record)
        attributed = runner.analyze_cell(
            cell, enabled, 2, "f" * 64, self.main_commit,
            self.current_source, self.main_snapshot, self.current_snapshot)
        self.assertIn("attribution_speedup", attributed)

        mixed = copy.deepcopy(enabled)
        mixed[1]["result"]["legacy"]["available"] = False
        self.bind_raw_output(mixed[1])
        expect_rejected(
            self,
            lambda: runner.analyze_cell(
                cell, mixed, 2, "f" * 64, self.main_commit,
                self.current_source, self.main_snapshot,
                self.current_snapshot),
            "contradicts its command")

        enabled_under_skip = copy.deepcopy(base["invocations"])
        enabled_under_skip[1]["result"]["legacy"]["available"] = True
        self.bind_raw_output(enabled_under_skip[1])
        expect_rejected(
            self,
            lambda: runner.analyze_cell(
                cell, enabled_under_skip, 2, "f" * 64,
                self.main_commit, self.current_source,
                self.main_snapshot, self.current_snapshot),
            "contradicts its command")

    def test_summary_rejects_coercions_and_partial_attribution(self) -> None:
        cell = runner.Cell("one", "gf8-all-k", 3, 2, 64, 1,
                           "low-R", "one-loss", 7, 1, 1, 0)
        main = Path("/tmp/exact-main")
        current = Path("/tmp/current-leopard2")
        ordinary = self.analyzed_cell(
            cell, main, current, 2, "f" * 64)
        metadata = {"run_contract": {"with_current_legacy": False}}
        summary = runner.summarize([ordinary], metadata)
        self.assertEqual(summary["valid_cell_count"], 1)

        coercible_valid = copy.deepcopy(ordinary)
        coercible_valid["valid"] = 1
        expect_rejected(
            self, lambda: runner.summarize([coercible_valid], metadata),
            "validity is not an exact boolean")

        coercible_speedup = copy.deepcopy(ordinary)
        coercible_speedup["speedup_main_over_leopard2"]["encode"] = True
        expect_rejected(
            self, lambda: runner.summarize([coercible_speedup], metadata),
            "not an exact JSON number")

        malformed_tags = copy.deepcopy(ordinary)
        malformed_tags["selected"]["gap_tags"] = ["valid", 1]
        expect_rejected(
            self, lambda: runner.summarize([malformed_tags], metadata),
            "gap tags are malformed")

        attributed = copy.deepcopy(ordinary)
        attributed["attribution_speedup"] = {
            "exact_main_over_current_tree_legacy": {
                "encode": 1.0, "decode_including_setup": 1.0,
            },
            "current_tree_legacy_over_leopard2": {
                "encode": 1.0, "decode_first_use": 1.0,
                "decode_at_reuse": 1.0,
                "decode_execution_only_optimistic": 1.0,
            },
        }
        attribution_metadata = {
            "run_contract": {"with_current_legacy": True},
        }
        expect_rejected(
            self,
            lambda: runner.summarize(
                [attributed, ordinary], attribution_metadata),
            "inconsistent legacy attribution")
        expect_rejected(
            self, lambda: runner.summarize([attributed], metadata),
            "contradicts the run contract")

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

    def test_run_closes_all_snapshots_after_injected_failure(self) -> None:
        original = runner._run_with_snapshot_owner
        captured: list[int] = []
        owners = []

        def fail_after_snapshots(options, lock_guard, snapshots):
            del options, lock_guard
            owners.append(snapshots)
            for executable, label in (
                    (Path("/bin/true"), "main fixture"),
                    (Path("/bin/false"), "current fixture")):
                unused_source, descriptor, unused_snapshot = \
                    snapshots.capture(executable, label)
                del unused_source, unused_snapshot
                captured.append(descriptor)
            raise RuntimeError("injected failure after snapshots")

        runner._run_with_snapshot_owner = fail_after_snapshots
        try:
            expect_rejected(
                self, lambda: runner.run(object(), object()),
                "injected failure after snapshots")
        finally:
            runner._run_with_snapshot_owner = original
        self.assertEqual(len(owners), 1)
        self.assertEqual(owners[0].descriptors, [])
        self.assertEqual(len(captured), 2)
        for descriptor in captured:
            with self.subTest(descriptor=descriptor):
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    def test_identity_bound_descriptor_handles_close_interrupt_and_recycle(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        real_close = os.close
        owner = runner.IdentityBoundDescriptor("close-before-effect fixture")
        descriptor = owner.open("/bin/true", os.O_RDONLY)
        interrupted = False

        def interrupt_once(value: int) -> None:
            nonlocal interrupted
            if value == descriptor and not interrupted:
                interrupted = True
                raise InjectedAbort("injected close-before-effect")
            real_close(value)

        with mock.patch.object(
                runner.os, "close", side_effect=interrupt_once), \
             self.assertRaisesRegex(
                 InjectedAbort, "close-before-effect"):
            owner.close()
        self.assertTrue(owner.closed)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

        owner = runner.IdentityBoundDescriptor("recycled descriptor fixture")
        descriptor = owner.open("/bin/true", os.O_RDONLY)
        replacement = -1

        def close_then_recycle(value: int) -> None:
            nonlocal replacement
            self.assertEqual(value, descriptor)
            real_close(value)
            replacement = os.open("/dev/null", os.O_RDONLY)
            self.assertEqual(replacement, descriptor)
            raise InjectedAbort("injected post-close recycle")

        with mock.patch.object(
                runner.os, "close", side_effect=close_then_recycle), \
             self.assertRaisesRegex(
                 InjectedAbort, "post-close recycle"):
            owner.close()
        try:
            self.assertTrue(owner.closed)
            os.fstat(replacement)
            owner.close()
            os.fstat(replacement)
        finally:
            real_close(replacement)

    def test_descriptor_adoption_fstat_interruption_is_self_cleaning(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        def descriptor_count() -> int:
            return len(os.listdir("/proc/self/fd"))

        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-adopt-fstat-") as directory:
            canonical = Path(directory) / "canonical.lock"
            original_canonical = runner.AUTHORITATIVE_LOCK
            runner.AUTHORITATIVE_LOCK = canonical
            guard = runner.CanonicalBenchmarkLock(canonical)
            real_fstat = os.fstat
            interrupted = False

            def interrupt_identity_once(descriptor: int):
                nonlocal interrupted
                metadata = real_fstat(descriptor)
                try:
                    target = os.readlink(f"/proc/self/fd/{descriptor}")
                except OSError:
                    target = ""
                if target == str(canonical) and not interrupted:
                    interrupted = True
                    raise InjectedAbort(
                        "injected fstat after canonical-lock open")
                return metadata

            before = descriptor_count()
            try:
                with mock.patch.object(
                        runner.os, "fstat",
                        side_effect=interrupt_identity_once), \
                     self.assertRaisesRegex(
                         InjectedAbort,
                         "fstat after canonical-lock open"):
                    guard.__enter__()
                self.assertTrue(interrupted)
                self.assertIsNone(guard.descriptor)
                self.assertEqual(descriptor_count(), before)
            finally:
                runner.AUTHORITATIVE_LOCK = original_canonical
                if guard.descriptor is not None:
                    guard.__exit__(None, None, None)

    def test_snapshot_owner_transfer_and_close_interruptions_do_not_leak(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        owner = runner.ExecutableSnapshotOwner()
        real_snapshot = runner.snapshot_executable
        captured: list[int] = []

        def interrupt_after_snapshot(*arguments, **keywords):
            result = real_snapshot(*arguments, **keywords)
            captured.append(result[1])
            raise InjectedAbort("injected snapshot transfer interruption")

        with mock.patch.object(
                runner, "snapshot_executable",
                side_effect=interrupt_after_snapshot), \
             self.assertRaisesRegex(
                 InjectedAbort, "snapshot transfer interruption"):
            owner.capture(Path("/bin/true"), "transfer fixture")
        self.assertEqual(owner.descriptors, [])
        self.assertEqual(len(captured), 1)
        with self.assertRaises(OSError):
            os.fstat(captured[0])

        unused_source, descriptor, unused_snapshot = owner.capture(
            Path("/bin/true"), "close fixture")
        del unused_source, unused_snapshot
        real_close = os.close
        interrupted = False

        def interrupt_close_once(value: int) -> None:
            nonlocal interrupted
            if value == descriptor and not interrupted:
                interrupted = True
                raise InjectedAbort("injected snapshot close-before-effect")
            real_close(value)

        with mock.patch.object(
                runner.os, "close", side_effect=interrupt_close_once), \
             self.assertRaisesRegex(
                 RuntimeError, "snapshot close-before-effect"):
            owner.close()
        self.assertEqual(owner.descriptors, [])
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_benchmark_lock_is_canonical_and_rejects_substitution(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-lock-") as directory:
            root = Path(directory)
            alternate = root / "alternate.lock"
            expect_rejected(
                self, lambda: runner.benchmark_lock(alternate).__enter__(),
                "require canonical lock")

            original_canonical = runner.AUTHORITATIVE_LOCK
            canonical = root / "canonical.lock"
            runner.AUTHORITATIVE_LOCK = canonical
            try:
                target = root / "symlink-target"
                target.write_bytes(b"")
                target.chmod(0o600)
                canonical.symlink_to(target.name)
                expect_rejected(
                    self,
                    lambda: runner.benchmark_lock(canonical).__enter__(),
                    "cannot open canonical benchmark lock")
                canonical.unlink()

                canonical.write_bytes(b"")
                canonical.chmod(0o644)
                expect_rejected(
                    self,
                    lambda: runner.benchmark_lock(canonical).__enter__(),
                    "unsafe metadata")
                canonical.chmod(0o600)

                alias = root / "hard-link"
                os.link(canonical, alias)
                expect_rejected(
                    self,
                    lambda: runner.benchmark_lock(canonical).__enter__(),
                    "unsafe metadata")
                alias.unlink()

                with runner.benchmark_lock(canonical) as guard:
                    identity = guard.validate_current()
                    self.assertEqual(identity["path"], str(canonical))
                    self.assertEqual(identity["uid"], os.getuid())
                    self.assertEqual(identity["mode"], 0o600)
                    self.assertEqual(identity["nlink"], 1)

                guard = runner.benchmark_lock(canonical)
                guard.__enter__()
                retained = root / "retained-lock"
                canonical.rename(retained)
                canonical.write_bytes(b"replacement")
                canonical.chmod(0o600)
                expect_rejected(
                    self, guard.validate_current, "replaced or changed")
                expect_rejected(
                    self, lambda: guard.__exit__(None, None, None),
                    "replaced or changed")
                self.assertIsNone(guard.descriptor)
            finally:
                runner.AUTHORITATIVE_LOCK = original_canonical

    def test_benchmark_lock_close_interruption_does_not_wedge_pair(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-lock-close-") as directory:
            canonical = Path(directory) / "canonical.lock"
            original_canonical = runner.AUTHORITATIVE_LOCK
            runner.AUTHORITATIVE_LOCK = canonical
            guard = runner.CanonicalBenchmarkLock(canonical)
            contender = -1
            try:
                guard.__enter__()
                descriptor = guard.descriptor
                self.assertIsNotNone(descriptor)
                real_close = os.close
                interrupted = False

                def interrupt_once(value: int) -> None:
                    nonlocal interrupted
                    if value == descriptor and not interrupted:
                        interrupted = True
                        raise InjectedAbort(
                            "injected campaign-lock close-before-effect")
                    real_close(value)

                with mock.patch.object(
                        runner.os, "close", side_effect=interrupt_once), \
                     self.assertRaisesRegex(
                         InjectedAbort, "campaign-lock close-before-effect"):
                    guard.__exit__(None, None, None)
                self.assertIsNone(guard.descriptor)
                self.assertIsNone(guard.identity)
                with self.assertRaises(OSError):
                    os.fstat(descriptor)
                contender = os.open(canonical, os.O_RDWR)
                fcntl.flock(
                    contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                runner.AUTHORITATIVE_LOCK = original_canonical
                if guard.descriptor is not None:
                    guard.__exit__(None, None, None)
                if contender >= 0:
                    os.close(contender)

    def test_benchmark_child_retains_lock_after_coordinator_sigkill(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-lock-inherit-") as directory:
            root = Path(directory)
            lock_path = root / "campaign.lock"
            marker = root / "benchmark-started"
            coordinator = os.fork()
            if coordinator == 0:
                try:
                    lock_descriptor = os.open(
                        lock_path, os.O_CREAT | os.O_RDWR, 0o600)
                    fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                    snapshot_descriptor = os.open(
                        "/bin/true", os.O_RDONLY)
                    child_code = (
                        "from pathlib import Path; import time; "
                        f"Path({str(marker)!r}).write_text('ready'); "
                        "time.sleep(1.5)"
                    )
                    runner.run_process_bounded(
                        [sys.executable, "-c", child_code], 5.0,
                        snapshot_descriptor, lock_descriptor)
                    os.close(snapshot_descriptor)
                    os.close(lock_descriptor)
                    os._exit(0)
                except BaseException:
                    os._exit(91)

            contender = -1
            try:
                deadline = time.monotonic() + 5.0
                while not marker.is_file() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertTrue(marker.is_file(),
                                "benchmark child did not start")
                contender = os.open(lock_path, os.O_RDWR)
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)

                os.kill(coordinator, signal.SIGKILL)
                unused_pid, status = os.waitpid(coordinator, 0)
                coordinator = -1
                self.assertTrue(os.WIFSIGNALED(status))
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)

                acquired = False
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline:
                    try:
                        fcntl.flock(
                            contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        acquired = True
                        break
                    except BlockingIOError:
                        time.sleep(0.02)
                self.assertTrue(
                    acquired,
                    "benchmark child leaked the inherited campaign lock")
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                if coordinator > 0:
                    try:
                        os.kill(coordinator, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    os.waitpid(coordinator, 0)
                if contender >= 0:
                    os.close(contender)

    def test_supervisor_identity_failure_cannot_launch_command(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-launch-gate-") as directory:
            root = Path(directory)
            marker = root / "must-not-start"
            lock_descriptor = os.open(
                root / "campaign.lock", os.O_CREAT | os.O_RDWR, 0o600)
            snapshot_descriptor = os.open("/bin/true", os.O_RDONLY)
            try:
                fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                with mock.patch.object(
                        runner, "_open_supervisor_proc_directory",
                        side_effect=FileNotFoundError(
                            errno.ENOENT, "gone")):
                    with self.assertRaisesRegex(
                            RuntimeError, "disappeared during launch"):
                        runner.run_process_bounded(
                            [
                                sys.executable, "-c",
                                "from pathlib import Path; "
                                f"Path({str(marker)!r}).write_text('unsafe')",
                            ],
                            5.0, snapshot_descriptor, lock_descriptor)
                self.assertFalse(
                    marker.exists(),
                    "unbound supervisor launched the benchmark command")

                with mock.patch.object(
                        runner, "_open_supervisor_proc_directory",
                        side_effect=FileNotFoundError(
                            errno.ENOENT, "gone")), \
                     mock.patch.object(
                         runner, "_collect_supervisor",
                         side_effect=RuntimeError(
                             "injected launch collection failure")):
                    with self.assertRaisesRegex(
                            RuntimeError,
                            "launch and collection both failed"):
                        runner.run_process_bounded(
                            [
                                sys.executable, "-c",
                                "from pathlib import Path; "
                                f"Path({str(marker)!r}).write_text('unsafe')",
                            ],
                            5.0, snapshot_descriptor, lock_descriptor)
                self.assertFalse(
                    marker.exists(),
                    "unbound supervisor with failed collection launched the "
                    "benchmark command")
            finally:
                os.close(snapshot_descriptor)
                os.close(lock_descriptor)

    def test_second_supervisor_pipe_failure_closes_first_pair(self) -> None:
        real_pipe2 = os.pipe2
        first_pair: list[int] = []

        def fail_second_pipe(flags: int) -> tuple[int, int]:
            if not first_pair:
                pair = real_pipe2(flags)
                first_pair.extend(pair)
                return pair
            raise OSError(errno.EMFILE, "injected descriptor exhaustion")

        with mock.patch.object(runner.os, "pipe2", side_effect=fail_second_pipe):
            with self.assertRaisesRegex(
                    OSError, "injected descriptor exhaustion"):
                runner.run_helper_bounded(["/bin/true"], timeout=1.0)
        self.assertEqual(len(first_pair), 2)
        for descriptor in first_pair:
            with self.subTest(descriptor=descriptor):
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    def test_unexpected_collection_failure_reaps_tree_and_releases_lock(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-collect-failure-") as directory:
            root = Path(directory)
            marker = root / "child.pid"
            lock_path = root / "campaign.lock"
            lock_descriptor = os.open(
                lock_path, os.O_CREAT | os.O_RDWR, 0o600)
            snapshot_descriptor = os.open("/bin/true", os.O_RDONLY)
            observed_identity: list[tuple[int, int, int, int]] = []

            def fail_after_child_started(
                    unused_supervisor, unused_control_read, unused_root,
                    unused_pidfds, unused_deadline, unused_stdout,
                    unused_stderr):
                deadline = time.monotonic() + 5.0
                while not marker.is_file() and time.monotonic() < deadline:
                    time.sleep(0.01)
                if not marker.is_file():
                    raise RuntimeError(
                        "injected collection failure before child marker")
                pid = int(marker.read_text(encoding="ascii"))
                record = runner._supervisor_proc_record(pid)
                if record is not None:
                    observed_identity.append(
                        (pid, record[1], record[3], record[4]))
                raise RuntimeError("injected collection failure")

            contender = -1
            try:
                fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                child_code = (
                    "import os,time; from pathlib import Path; "
                    f"marker=Path({str(marker)!r}); "
                    "temporary=marker.with_name("
                    "marker.name+'.tmp-'+str(os.getpid())); "
                    "temporary.write_text(str(os.getpid()),encoding='ascii'); "
                    "os.replace(temporary,marker); "
                    "time.sleep(60)"
                )
                with mock.patch.object(
                        runner, "_collect_supervisor",
                        side_effect=fail_after_child_started):
                    expect_rejected(
                        self,
                        lambda: runner.run_process_bounded(
                            [sys.executable, "-c", child_code], 30.0,
                            snapshot_descriptor, lock_descriptor),
                        "injected collection failure")
                os.close(snapshot_descriptor)
                snapshot_descriptor = -1
                os.close(lock_descriptor)
                lock_descriptor = -1

                contender = os.open(lock_path, os.O_RDWR)
                fcntl.flock(contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
                self.assertTrue(observed_identity)
                pid, starttime, proc_dev, proc_inode = observed_identity[0]
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline:
                    record = runner._supervisor_proc_record(pid)
                    if (record is None or
                            (record[1], record[3], record[4]) !=
                            (starttime, proc_dev, proc_inode)):
                        break
                    time.sleep(0.01)
                else:
                    self.fail(
                        "unexpected collection failure left the exact child "
                        "process alive")
            finally:
                if observed_identity:
                    identity = observed_identity[0]
                    record = runner._supervisor_proc_record(identity[0])
                    if (record is not None and
                            (record[1], record[3], record[4]) ==
                            (identity[1], identity[2], identity[3])):
                        retained = runner.RetainedPidfds()
                        try:
                            if retained.retain(identity[0], identity) is not None:
                                retained.signal(identity, signal.SIGKILL)
                        finally:
                            retained.close()
                if contender >= 0:
                    os.close(contender)
                if snapshot_descriptor >= 0:
                    os.close(snapshot_descriptor)
                if lock_descriptor >= 0:
                    os.close(lock_descriptor)

    def test_pidfd_retention_rejects_same_tick_proc_inode_reuse(self) -> None:
        proc_descriptor = os.open("/dev/null", os.O_RDONLY)
        pidfd_descriptor = os.open("/dev/null", os.O_RDONLY)
        records = iter((
            (1, 100, "S", 7, 11),
            (1, 100, "S", 7, 12),
        ))
        registry = runner.RetainedPidfds()
        try:
            with mock.patch.object(
                    runner, "_open_supervisor_proc_directory",
                    return_value=proc_descriptor), \
                 mock.patch.object(
                     runner, "_read_supervisor_proc_record",
                     side_effect=lambda unused_pid, unused_descriptor:
                     next(records)), \
                 mock.patch.object(
                     runner, "_open_pidfd", return_value=pidfd_descriptor):
                self.assertIsNone(registry.retain(
                    1234, (1234, 100, 7, 11)))
            self.assertEqual(registry.descriptors, {})
            with self.assertRaises(OSError):
                os.fstat(proc_descriptor)
            proc_descriptor = -1
            with self.assertRaises(OSError):
                os.fstat(pidfd_descriptor)
            pidfd_descriptor = -1
        finally:
            registry.close()
            if proc_descriptor >= 0:
                os.close(proc_descriptor)
            if pidfd_descriptor >= 0:
                os.close(pidfd_descriptor)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux retained pidfd ownership test")
    def test_pidfd_transfer_close_and_recycle_are_identity_bound(self) -> None:
        class InjectedAbort(BaseException):
            pass

        registry = runner.RetainedPidfds()

        class InsertThenInterrupt(dict):
            def __setitem__(self, key, value) -> None:
                super().__setitem__(key, value)
                raise InjectedAbort(
                    "injected pidfd transfer interruption")

        registry.descriptors = InsertThenInterrupt()
        with self.assertRaisesRegex(
                InjectedAbort, "pidfd transfer interruption"):
            registry.retain(os.getpid())
        self.assertEqual(len(registry.descriptors), 1)
        transferred = next(iter(registry.descriptors.values()))
        descriptor = transferred.descriptor
        self.assertIsNotNone(descriptor)
        os.fstat(descriptor)
        registry.close()
        self.assertEqual(registry.descriptors, {})
        with self.assertRaises(OSError):
            os.fstat(descriptor)

        registry = runner.RetainedPidfds()
        identity = registry.retain(os.getpid())
        self.assertIsNotNone(identity)
        retained = registry.descriptors[identity]
        descriptor = retained.descriptor
        self.assertIsNotNone(descriptor)
        real_close = os.close
        interrupted = False

        def interrupt_once(value: int) -> None:
            nonlocal interrupted
            if value == descriptor and not interrupted:
                interrupted = True
                raise InjectedAbort("injected pidfd close-before-effect")
            real_close(value)

        with mock.patch.object(
                runner.os, "close", side_effect=interrupt_once), \
             self.assertRaisesRegex(
                 RuntimeError, "pidfd close-before-effect"):
            registry.close()
        self.assertEqual(registry.descriptors, {})
        with self.assertRaises(OSError):
            os.fstat(descriptor)

        registry = runner.RetainedPidfds()
        identity = registry.retain(os.getpid())
        self.assertIsNotNone(identity)
        retained = registry.descriptors[identity]
        descriptor = retained.descriptor
        self.assertIsNotNone(descriptor)
        real_close(descriptor)
        replacements: list[int] = []
        while not replacements or replacements[-1] != descriptor:
            replacements.append(os.open("/dev/null", os.O_RDONLY))
            self.assertLessEqual(
                replacements[-1], descriptor,
                "descriptor recycling skipped the retained pidfd number")
        replacement = replacements[-1]
        try:
            with self.assertRaisesRegex(
                    RuntimeError, "recycled before close"):
                registry.close()
            self.assertEqual(registry.descriptors, {})
            os.fstat(replacement)
        finally:
            for replacement_descriptor in replacements:
                real_close(replacement_descriptor)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux retained pidfd ownership test")
    def test_pidfd_open_interruption_closes_held_proc_directory(self) -> None:
        class InjectedAbort(BaseException):
            pass

        before = len(os.listdir("/proc/self/fd"))
        registry = runner.RetainedPidfds()
        with mock.patch.object(
                runner, "_open_pidfd",
                side_effect=InjectedAbort(
                    "injected pidfd-open interruption")), \
             self.assertRaisesRegex(
                 InjectedAbort, "pidfd-open interruption"):
            registry.retain(os.getpid())
        self.assertEqual(registry.descriptors, {})
        self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_retained_pidfd_signal_never_reopens_numeric_pid(self) -> None:
        proc_descriptor = os.open("/dev/null", os.O_RDONLY)
        pidfd_descriptor = os.open("/dev/null", os.O_RDONLY)
        registry = runner.RetainedPidfds()
        signaled = []
        records = iter((
            (1, 100, "S", 7, 11),
            (1, 100, "S", 7, 11),
        ))
        try:
            with mock.patch.object(
                    runner, "_open_supervisor_proc_directory",
                    return_value=proc_descriptor), \
                 mock.patch.object(
                     runner, "_read_supervisor_proc_record",
                     side_effect=lambda unused_pid, unused_descriptor:
                     next(records)), \
                 mock.patch.object(
                     runner, "_open_pidfd", return_value=pidfd_descriptor):
                identity = registry.retain(
                    1234, (1234, 100, 7, 11))
            self.assertEqual(identity, (1234, 100, 7, 11))
            proc_descriptor = -1
            pidfd_descriptor = -1
            with mock.patch.object(
                    runner, "_open_supervisor_proc_directory",
                    side_effect=AssertionError(
                        "signal path reopened numeric PID")), \
                 mock.patch.object(
                     runner, "_send_pidfd_signal",
                     side_effect=lambda retained, number:
                     signaled.append((retained, number))):
                registry.signal(identity, signal.SIGKILL)
            self.assertEqual(len(signaled), 1)
            self.assertEqual(signaled[0][1], signal.SIGKILL)
        finally:
            registry.close()
            if proc_descriptor >= 0:
                os.close(proc_descriptor)
            if pidfd_descriptor >= 0:
                os.close(pidfd_descriptor)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux bounded helper supervision test")
    def test_helper_supervisor_drains_small_pipes_concurrently(self) -> None:
        real_popen = subprocess.Popen
        real_pipe2 = os.pipe2

        def small_popen(*arguments, **keywords):
            process = real_popen(*arguments, **keywords)
            self.assertIsNotNone(process.stdout)
            self.assertIsNotNone(process.stderr)
            fcntl.fcntl(process.stdout.fileno(), fcntl.F_SETPIPE_SZ, 4096)
            fcntl.fcntl(process.stderr.fileno(), fcntl.F_SETPIPE_SZ, 4096)
            return process

        def small_pipe2(flags):
            descriptors = real_pipe2(flags)
            fcntl.fcntl(descriptors[0], fcntl.F_SETPIPE_SZ, 4096)
            return descriptors

        count = 256 * 1024
        child = (
            "import os\n"
            f"first=b'a'*{count}\n"
            "offset=0\n"
            "while offset<len(first): "
            " offset+=os.write(1,first[offset:])\n"
            "second=b'b'*len(first)\n"
            "offset=0\n"
            "while offset<len(second): "
            " offset+=os.write(2,second[offset:])\n"
        )
        with mock.patch.object(
                runner.subprocess, "Popen", side_effect=small_popen), \
             mock.patch.object(runner.os, "pipe2", side_effect=small_pipe2):
            completed = runner.run_helper_bounded(
                [sys.executable, "-c", child], timeout=5.0,
                max_stdout=count, max_stderr=count)
        self.assertEqual(completed.returncode, 0)
        self.assertEqual(completed.stdout, b"a" * count)
        self.assertEqual(completed.stderr, b"b" * count)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux bounded helper supervision test")
    def test_helper_timeout_output_limit_and_retained_pipe(self) -> None:
        with self.assertRaises(subprocess.TimeoutExpired):
            runner.run_helper_bounded(
                [sys.executable, "-c", "import time;time.sleep(30)"],
                timeout=0.2, max_stdout=1024, max_stderr=1024)

        with self.assertRaisesRegex(
                RuntimeError, "output exceeded|output byte limit"):
            runner.run_helper_bounded(
                [sys.executable, "-c",
                 "import os;os.write(1,b'x'*8192)"],
                timeout=2.0, max_stdout=4096, max_stderr=1024)

        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-retained-pipe-") as directory:
            pid_path = Path(directory) / "descendant.pid"
            child = (
                "import os,pathlib,sys,time\n"
                "pid=os.fork()\n"
                "if pid==0:\n"
                " os.setsid()\n"
                " marker=pathlib.Path(sys.argv[1])\n"
                " temporary=marker.with_name("
                "marker.name+'.tmp-'+str(os.getpid()))\n"
                " temporary.write_text("
                "str(os.getpid()),encoding='ascii')\n"
                " os.replace(temporary,marker)\n"
                " time.sleep(30)\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+3\n"
                "while not os.path.exists(sys.argv[1]) and "
                "time.monotonic()<deadline: time.sleep(.01)\n"
                "sys.exit(0 if os.path.exists(sys.argv[1]) else 92)\n"
            )
            started = time.monotonic()
            completed = runner.run_helper_bounded(
                [sys.executable, "-c", child, str(pid_path)],
                timeout=4.0, max_stdout=1024, max_stderr=1024)
            elapsed = time.monotonic() - started
            descendant = int(pid_path.read_text(encoding="ascii"))
            self.assertEqual(completed.returncode, 0)
            self.assertLess(elapsed, 2.0)
            self.assertFalse(Path("/proc", str(descendant)).exists())

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux supervisor fallback containment test")
    def test_supervisor_fallback_kills_separate_session_descendant(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-supervisor-fallback-") as directory:
            root = Path(directory)
            lock_path = root / "campaign.lock"
            marker = root / "detached.pid"
            lock_descriptor = os.open(
                lock_path, os.O_CREAT | os.O_RDWR, 0o600)
            supervisor = None
            contender = -1
            detached_pid = -1
            try:
                fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                program = (
                    "import ctypes,os,pathlib,signal,sys,time\n"
                    "libc=ctypes.CDLL(None,use_errno=True)\n"
                    "assert libc.prctl(36,1,0,0,0)==0\n"
                    "lockfd=int(sys.argv[1])\n"
                    "pid=os.fork()\n"
                    "if pid==0:\n"
                    " os.setsid()\n"
                    " daemon=os.fork()\n"
                    " if daemon!=0: os._exit(0)\n"
                    " for fd in range(3,256):\n"
                    "  if fd!=lockfd:\n"
                    "   try: os.close(fd)\n"
                    "   except OSError: pass\n"
                    " os.fstat(lockfd)\n"
                    " marker=pathlib.Path(sys.argv[2])\n"
                    " temporary=marker.with_name("
                    "marker.name+'.tmp-'+str(os.getpid()))\n"
                    " temporary.write_text("
                    "str(os.getpid()),encoding='ascii')\n"
                    " os.replace(temporary,marker)\n"
                    " signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
                    " time.sleep(30)\n"
                    "deadline=time.monotonic()+5\n"
                    "while not os.path.exists(sys.argv[2]) and "
                    "time.monotonic()<deadline: time.sleep(.01)\n"
                    "time.sleep(30)\n"
                )
                supervisor = subprocess.Popen(
                    [
                        sys.executable, "-c", program,
                        str(lock_descriptor), str(marker),
                    ],
                    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, start_new_session=True,
                    pass_fds=(lock_descriptor,))
                wait_for_path(self, marker, "detached supervisor child")
                detached_pid = int(marker.read_text(encoding="ascii"))
                record = runner._supervisor_proc_record(supervisor.pid)
                self.assertIsNotNone(record)
                runner.terminate_supervisor_tree_bounded(
                    supervisor, (
                        supervisor.pid, record[1], record[3], record[4]))
                self.assertIsNotNone(supervisor.returncode)

                os.close(lock_descriptor)
                lock_descriptor = -1
                contender = os.open(lock_path, os.O_RDWR)
                fcntl.flock(
                    contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
                detached = runner._supervisor_proc_record(detached_pid)
                self.assertTrue(
                    detached is None or detached[2] in {"Z", "X", "x"})
            finally:
                if supervisor is not None and supervisor.poll() is None:
                    try:
                        os.killpg(supervisor.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    supervisor.wait(timeout=5)
                if lock_descriptor >= 0:
                    os.close(lock_descriptor)
                if contender >= 0:
                    os.close(contender)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_detached_grandchild_cannot_retain_benchmark_lock(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-detached-lock-") as directory:
            root = Path(directory)
            lock_path = root / "campaign.lock"
            ready = root / "detached-ready"
            escaped = root / "detached-escaped"
            lock_descriptor = os.open(
                lock_path, os.O_CREAT | os.O_RDWR, 0o600)
            snapshot_descriptor = os.open("/bin/true", os.O_RDONLY)
            contender = -1
            detached_pid = -1
            try:
                fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                child = (
                    "import os,pathlib,sys,time\n"
                    "lock_fd=int(sys.argv[1])\n"
                    "pid=os.fork()\n"
                    "if pid == 0:\n"
                    " os.setsid()\n"
                    " daemon=os.fork()\n"
                    " if daemon != 0: os._exit(0)\n"
                    " null=os.open('/dev/null', os.O_WRONLY)\n"
                    " os.dup2(null,1);os.dup2(null,2);os.close(null)\n"
                    " os.fstat(lock_fd)\n"
                    " ready=pathlib.Path(sys.argv[2])\n"
                    " temporary=ready.with_name("
                    "ready.name+'.tmp-'+str(os.getpid()))\n"
                    " temporary.write_text("
                    "str(os.getpid()),encoding='ascii')\n"
                    " os.replace(temporary,ready)\n"
                    " time.sleep(1.0)\n"
                    " open(sys.argv[3],'w').write('escaped')\n"
                    " os._exit(0)\n"
                    "deadline=time.monotonic()+5\n"
                    "while not os.path.exists(sys.argv[2]) and "
                    "time.monotonic()<deadline: time.sleep(.01)\n"
                    "sys.exit(0 if os.path.exists(sys.argv[2]) else 93)\n"
                )
                completed = runner.run_process_bounded(
                    [
                        sys.executable, "-c", child, str(lock_descriptor),
                        str(ready), str(escaped),
                    ],
                    3.0, snapshot_descriptor, lock_descriptor)
                self.assertEqual(completed.returncode, 0)
                detached_pid = int(ready.read_text(encoding="utf-8"))

                # Only an escaped descendant can keep this open-file-description
                # lock alive after the coordinator's descriptor is closed.
                os.close(lock_descriptor)
                lock_descriptor = -1
                contender = os.open(lock_path, os.O_RDWR)
                fcntl.flock(
                    contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
                time.sleep(1.05)
                self.assertFalse(escaped.exists())
                self.assertFalse(
                    Path("/proc", str(detached_pid)).exists())
            finally:
                if lock_descriptor >= 0:
                    os.close(lock_descriptor)
                os.close(snapshot_descriptor)
                if contender >= 0:
                    os.close(contender)

    def test_provenance_helper_retains_lock_after_coordinator_sigkill(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-helper-lock-inherit-") as directory:
            root = Path(directory)
            lock_path = root / "campaign.lock"
            marker = root / "helper-started"
            coordinator = os.fork()
            if coordinator == 0:
                try:
                    lock_descriptor = os.open(
                        lock_path, os.O_CREAT | os.O_RDWR, 0o600)
                    fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                    child_code = (
                        "from pathlib import Path; import os, time; "
                        f"os.fstat({lock_descriptor}); "
                        f"Path({str(marker)!r}).write_text('ready'); "
                        "time.sleep(1.5)"
                    )
                    provenance_module._run(
                        [sys.executable, "-c", child_code],
                        "lock-retaining provenance helper", timeout=5.0,
                        inherited_descriptors=(lock_descriptor,))
                    os.close(lock_descriptor)
                    os._exit(0)
                except BaseException:
                    os._exit(91)

            contender = -1
            try:
                wait_for_path(self, marker, "provenance helper")
                contender = os.open(lock_path, os.O_RDWR)
                require_lock_blocked(
                    self, contender, "live provenance helper")
                os.kill(coordinator, signal.SIGKILL)
                unused_pid, status = os.waitpid(coordinator, 0)
                coordinator = -1
                self.assertTrue(os.WIFSIGNALED(status))
                require_lock_blocked(
                    self, contender, "orphaned provenance helper")
                wait_for_lock_release(
                    self, contender, "provenance helper")
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                if coordinator > 0:
                    try:
                        os.kill(coordinator, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    os.waitpid(coordinator, 0)
                if contender >= 0:
                    os.close(contender)

    def test_nested_build_child_retains_lock_after_coordinator_sigkill(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-build-child-lock-") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            source.mkdir()
            lock_path = root / "campaign.lock"
            marker = root / "nested-build-child-started"
            lock_descriptor = os.open(
                lock_path, os.O_CREAT | os.O_RDWR, 0o600)
            child_script = source / "lock_child.py"
            child_script.write_text(
                "import os\n"
                "from pathlib import Path\n"
                "import sys\n"
                "import time\n"
                "os.fstat(int(sys.argv[1]))\n"
                "Path(sys.argv[2]).write_text('ready')\n"
                "time.sleep(1.5)\n",
                encoding="utf-8")
            (source / "CMakeLists.txt").write_text(
                "cmake_minimum_required(VERSION 3.16)\n"
                "project(lock_probe NONE)\n"
                "add_custom_target(lock_probe ALL\n"
                f"  COMMAND {sys.executable} {child_script} "
                f"{lock_descriptor} {marker}\n"
                "  VERBATIM)\n",
                encoding="utf-8")
            subprocess.run(
                ["/usr/bin/cmake", "-S", str(source), "-B", str(build),
                 "-G", "Unix Makefiles"],
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=True)

            coordinator = os.fork()
            if coordinator == 0:
                try:
                    fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                    provenance_module._run(
                        ["/usr/bin/cmake", "--build", str(build),
                         "--parallel", "1"],
                        "nested lock-retaining build", timeout=10.0,
                        inherited_descriptors=(lock_descriptor,))
                    os.close(lock_descriptor)
                    os._exit(0)
                except BaseException:
                    os._exit(91)
            os.close(lock_descriptor)

            contender = -1
            try:
                wait_for_path(self, marker, "nested build child")
                contender = os.open(lock_path, os.O_RDWR)
                require_lock_blocked(
                    self, contender, "live nested build child")
                os.kill(coordinator, signal.SIGKILL)
                unused_pid, status = os.waitpid(coordinator, 0)
                coordinator = -1
                self.assertTrue(os.WIFSIGNALED(status))
                require_lock_blocked(
                    self, contender, "orphaned nested build child")
                wait_for_lock_release(
                    self, contender, "nested build child")
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                if coordinator > 0:
                    try:
                        os.kill(coordinator, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    os.waitpid(coordinator, 0)
                if contender >= 0:
                    os.close(contender)

    def test_normal_helper_completion_releases_lock_without_leak(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-helper-lock-cleanup-") as directory:
            root = Path(directory)
            canonical = root / "canonical.lock"
            original_canonical = runner.AUTHORITATIVE_LOCK
            runner.AUTHORITATIVE_LOCK = canonical
            contender = -1
            try:
                with runner.benchmark_lock(canonical) as guard:
                    lock_descriptor = guard.descriptor
                    self.assertIsNotNone(lock_descriptor)
                    provenance_module._run(
                        [sys.executable, "-c",
                         f"import os; os.fstat({lock_descriptor})"],
                        "normally completing provenance helper",
                        inherited_descriptors=(lock_descriptor,))
                    contender = os.open(canonical, os.O_RDWR)
                    require_lock_blocked(
                        self, contender, "campaign owner")
                wait_for_lock_release(
                    self, contender, "normally completing helper", 1.0)
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                runner.AUTHORITATIVE_LOCK = original_canonical
                if contender >= 0:
                    os.close(contender)

    def test_timeout_and_error_descendants_are_reaped_before_return(
            self) -> None:
        for mode in ("timeout", "error"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                    prefix=f"leo2-helper-{mode}-descendant-") as directory:
                root = Path(directory)
                canonical = root / "canonical.lock"
                nested_marker = root / "nested-started"
                original_canonical = runner.AUTHORITATIVE_LOCK
                runner.AUTHORITATIVE_LOCK = canonical
                contender = -1
                try:
                    with runner.benchmark_lock(canonical) as guard:
                        lock_descriptor = guard.descriptor
                        self.assertIsNotNone(lock_descriptor)
                        nested_code = (
                            "from pathlib import Path; import os, time; "
                            f"os.fstat({lock_descriptor}); "
                            f"marker=Path({str(nested_marker)!r}); "
                            "temporary=marker.with_name("
                            "marker.name+'.tmp-'+str(os.getpid())); "
                            "temporary.write_text("
                            "str(os.getpid()),encoding='ascii'); "
                            "os.replace(temporary,marker); "
                            "time.sleep(30)"
                        )
                        nested_command = [
                            sys.executable, "-c", nested_code]
                        final_action = (
                            "time.sleep(10)"
                            if mode == "timeout" else "sys.exit(7)")
                        helper_code = (
                            "from pathlib import Path; "
                            "import os, subprocess, sys, time; "
                            "subprocess.Popen("
                            f"{nested_command!r}, stdin=subprocess.DEVNULL, "
                            "stdout=subprocess.DEVNULL, "
                            "stderr=subprocess.DEVNULL, "
                            f"pass_fds=({lock_descriptor},), "
                            "start_new_session=True); "
                            f"marker=Path({str(nested_marker)!r}); "
                            "deadline=time.monotonic()+2; "
                            "\nwhile not marker.is_file() and "
                            "time.monotonic()<deadline: time.sleep(.01)\n"
                            "if not marker.is_file(): sys.exit(93)\n"
                            f"{final_action}"
                        )
                        if mode == "timeout":
                            with self.assertRaisesRegex(
                                    provenance_module.BuildProvenanceError,
                                    "exceeded 0.250 seconds"):
                                provenance_module._run(
                                    [sys.executable, "-c", helper_code],
                                    "timed-out provenance helper",
                                    timeout=0.25,
                                    inherited_descriptors=(
                                        lock_descriptor,))
                        else:
                            with self.assertRaisesRegex(
                                    provenance_module.BuildProvenanceError,
                                    "rc=7"):
                                provenance_module._run(
                                    [sys.executable, "-c", helper_code],
                                    "failed provenance helper", timeout=5.0,
                                    inherited_descriptors=(
                                        lock_descriptor,))
                        wait_for_path(
                            self, nested_marker,
                            f"{mode} helper descendant")
                        nested_pid = int(nested_marker.read_text(
                            encoding="utf-8"))
                        self.assertFalse(
                            Path("/proc", str(nested_pid)).exists(),
                            f"{mode} helper descendant survived _run")
                    contender = os.open(canonical, os.O_RDWR)
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fcntl.flock(contender, fcntl.LOCK_UN)
                finally:
                    runner.AUTHORITATIVE_LOCK = original_canonical
                    if contender >= 0:
                        os.close(contender)

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
                    "main", command, 5.0, descriptor, snapshot, descriptor)
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

    def test_isa_gate_inspects_sealed_snapshot_after_evex_path_swap(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-isa-snapshot-") as directory:
            root = Path(directory)
            avx2_source = root / "avx2.s"
            evex_source = root / "evex.s"
            avx2_executable = root / "benchmark"
            evex_executable = root / "evex"
            common_tail = (
                "    mov $60, %eax\n"
                "    xor %edi, %edi\n"
                "    syscall\n"
            )
            avx2_source.write_text(
                ".text\n.global _start\n_start:\n"
                "    vpxor %ymm0, %ymm0, %ymm0\n" + common_tail,
                encoding="utf-8")
            evex_source.write_text(
                ".text\n.global _start\n_start:\n"
                "    vpxord %zmm0, %zmm0, %zmm0\n"
                "    vpxor %ymm0, %ymm0, %ymm0\n" + common_tail,
                encoding="utf-8")
            for source, executable in (
                    (avx2_source, avx2_executable),
                    (evex_source, evex_executable)):
                subprocess.run([
                    "/usr/bin/cc", "-nostdlib", "-no-pie",
                    "-Wl,-e,_start", "-o", str(executable), str(source),
                ], check=True, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)

            unused_source, avx2_descriptor, unused_snapshot = \
                runner.snapshot_executable(
                    avx2_executable, "AVX2 fixture")
            del unused_source, unused_snapshot
            try:
                shutil.copyfile(evex_executable, avx2_executable)
                avx2_executable.chmod(0o755)
                sealed = runner.inspect_pure_avx2(
                    avx2_descriptor, "sealed AVX2 fixture")
                self.assertEqual(
                    sealed["evex_prefixed_instruction_count"], 0)
                self.assertGreater(
                    sealed["ymm_operand_instruction_count"], 0)
            finally:
                os.close(avx2_descriptor)

            unused_source, evex_descriptor, unused_snapshot = \
                runner.snapshot_executable(
                    avx2_executable, "EVEX substitution fixture")
            del unused_source, unused_snapshot
            try:
                expect_rejected(
                    self,
                    lambda: runner.inspect_pure_avx2(
                        evex_descriptor, "EVEX substitution fixture"),
                    "EVEX-prefixed")
            finally:
                os.close(evex_descriptor)

    def test_manifest_json_round_trip_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leo2-allk-json-") as directory:
            path = Path(directory) / "manifest.json"
            value = {"z": [3, 2, 1], "a": self.current_source}
            runner.write_json_atomic(path, value)
            self.assertEqual(json.loads(path.read_text(encoding="utf-8")), value)
            self.assertFalse(path.with_name(path.name + ".tmp").exists())
            with self.assertRaises(ValueError):
                runner.write_json_atomic(path, {"invalid": float("inf")})

    def test_atomic_json_write_never_follows_symlink_victims(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-json-symlink-") as directory:
            root = Path(directory)
            victim = root / "victim"
            victim.write_bytes(b"do-not-touch")
            path = root / "manifest.json"

            predictable = path.with_name(path.name + ".tmp")
            predictable.symlink_to(victim.name)
            runner.write_json_atomic(path, {"safe": True})
            self.assertEqual(victim.read_bytes(), b"do-not-touch")
            self.assertTrue(predictable.is_symlink())
            self.assertEqual(
                json.loads(path.read_text(encoding="utf-8")),
                {"safe": True})

            redirected = root / "redirected.json"
            redirected.symlink_to(victim.name)
            expect_rejected(
                self,
                lambda: runner.write_json_atomic(
                    redirected, {"unsafe": True}),
                "unsafe evidence destination")
            self.assertEqual(victim.read_bytes(), b"do-not-touch")

    def test_atomic_output_directory_lifetime_binding_and_permissions(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-binding-") as directory:
            root = Path(directory)
            output = root / "output"
            previous_umask = os.umask(0o777)
            try:
                evidence = runner.AtomicEvidenceDirectory.open_or_create(
                    output)
            finally:
                os.umask(previous_umask)
            try:
                runner.write_json_atomic(
                    output / "cells/one.json", {"value": 1},
                    evidence_directory=evidence)
                self.assertEqual(
                    stat.S_IMODE(output.stat().st_mode), 0o700)
                self.assertEqual(
                    stat.S_IMODE((output / "cells").stat().st_mode), 0o700)
                self.assertEqual(
                    stat.S_IMODE(
                        (output / "cells/one.json").stat().st_mode),
                    0o600)

                retained = root / "retained"
                output.rename(retained)
                output.mkdir(mode=0o700)
                with self.assertRaisesRegex(
                        RuntimeError, "replaced or changed"):
                    runner.write_json_atomic(
                        output / "after.json", {"unsafe": True},
                        evidence_directory=evidence)
                self.assertFalse((output / "after.json").exists())
                self.assertFalse((retained / "after.json").exists())
            finally:
                evidence.close()

    def test_atomic_output_directory_validation_failures_close_child_fds(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        def descriptor_count() -> int:
            return len(os.listdir("/proc/self/fd"))

        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-fd-") as directory:
            root = Path(directory)
            absolute_child = root / "absolute-child"
            absolute_child.mkdir(mode=0o700)
            real_fstat = os.fstat

            def abort_on_absolute_child(descriptor: int):
                metadata = real_fstat(descriptor)
                try:
                    target = os.readlink(f"/proc/self/fd/{descriptor}")
                except OSError:
                    target = ""
                if target == str(absolute_child):
                    raise InjectedAbort("injected absolute-child validation")
                return metadata

            before = descriptor_count()
            with mock.patch.object(
                    runner.os, "fstat",
                    side_effect=abort_on_absolute_child), \
                 self.assertRaisesRegex(
                     InjectedAbort, "absolute-child validation"):
                runner.AtomicEvidenceDirectory._open_absolute(absolute_child)
            self.assertEqual(descriptor_count(), before)

            output = root / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            try:
                unsafe = output / "unsafe"
                unsafe.mkdir(mode=0o755)
                before = descriptor_count()
                with self.assertRaisesRegex(
                        RuntimeError, "evidence child directory is unsafe"):
                    evidence._open_parent("unsafe/result.json", create=False)
                self.assertEqual(descriptor_count(), before)
            finally:
                evidence.close()

            before = descriptor_count()
            with mock.patch.object(
                    runner.AtomicEvidenceDirectory, "validate_current",
                    side_effect=InjectedAbort(
                        "injected result validation")), \
                 self.assertRaisesRegex(
                     InjectedAbort, "result validation"):
                runner.AtomicEvidenceDirectory.open_or_create(
                    root / "interrupted-output")
            self.assertEqual(descriptor_count(), before)

    def test_atomic_directory_traversal_close_interruptions_do_not_leak(
            self) -> None:
        class InjectedAbort(BaseException):
            pass

        def descriptor_count() -> int:
            return len(os.listdir("/proc/self/fd"))

        real_close = os.close

        def close_interrupter(
                target: Path | str, message: str, *,
                matching_close: int = 1):
            target_text = str(target)
            matches = 0
            interrupted = False

            def interrupt(value: int) -> None:
                nonlocal matches, interrupted
                try:
                    current = os.readlink(f"/proc/self/fd/{value}")
                except OSError:
                    current = ""
                if current == target_text:
                    matches += 1
                    if matches == matching_close and not interrupted:
                        interrupted = True
                        raise InjectedAbort(message)
                real_close(value)

            return interrupt

        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-close-") as directory:
            root = Path(directory)
            absolute_child = root / "absolute-child"
            absolute_child.mkdir(mode=0o700)

            before = descriptor_count()
            with mock.patch.object(
                    runner.os, "close",
                    side_effect=close_interrupter(
                        "/", "injected absolute traversal close")), \
                 self.assertRaisesRegex(
                     InjectedAbort, "absolute traversal close"):
                runner.AtomicEvidenceDirectory._open_absolute(absolute_child)
            self.assertEqual(descriptor_count(), before)

            before = descriptor_count()
            with mock.patch.object(
                    runner.os, "close",
                    side_effect=close_interrupter(
                        "/", "injected create traversal close")), \
                 self.assertRaisesRegex(
                     InjectedAbort, "create traversal close"):
                runner.AtomicEvidenceDirectory.open_or_create(
                    root / "created-output")
            self.assertEqual(descriptor_count(), before)

            output = root / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            try:
                before = descriptor_count()
                with mock.patch.object(
                        runner.os, "close",
                        side_effect=close_interrupter(
                            output, "injected parent traversal close",
                            matching_close=2)), \
                     self.assertRaisesRegex(
                         InjectedAbort, "parent traversal close"):
                    evidence._open_parent(
                        "nested/result.json", create=True)
                self.assertEqual(descriptor_count(), before)
                evidence.validate_current()

                unsafe = output / "unsafe"
                unsafe.mkdir(mode=0o755)
                before = descriptor_count()
                with mock.patch.object(
                        runner.os, "close",
                        side_effect=close_interrupter(
                            unsafe, "injected child cleanup close")), \
                     self.assertRaisesRegex(
                         RuntimeError,
                         "evidence child directory is unsafe"):
                    evidence._open_parent(
                        "unsafe/result.json", create=False)
                self.assertEqual(descriptor_count(), before)
                evidence.validate_current()
            finally:
                evidence.close()

    def test_atomic_publication_detects_parent_replacement_during_replace(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-race-") as directory:
            root = Path(directory)
            output = root / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            retained = root / "retained"
            real_replace = os.replace

            def replace_after_substitution(
                    source, destination, *arguments, **keywords):
                output.rename(retained)
                output.mkdir(mode=0o700)
                return real_replace(
                    source, destination, *arguments, **keywords)

            try:
                with mock.patch.object(
                        runner.os, "replace",
                        side_effect=replace_after_substitution), \
                     self.assertRaisesRegex(
                         RuntimeError, "replaced"):
                    runner.write_json_atomic(
                        output / "result.json", {"value": 1},
                        evidence_directory=evidence)
                self.assertFalse((output / "result.json").exists())
                self.assertTrue((retained / "result.json").is_file())
            finally:
                evidence.close()

    def test_atomic_publication_cleanup_faults_preserve_diagnostics_and_fds(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-cleanup-") as directory:
            output = Path(directory) / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            try:
                real_open = os.open
                real_close = os.close
                real_fsync = os.fsync
                temporary_descriptor = None
                fsync_failed = False
                close_failed = False

                def tracked_open(path, flags, *arguments, **keywords):
                    nonlocal temporary_descriptor
                    descriptor = real_open(
                        path, flags, *arguments, **keywords)
                    if (isinstance(path, str) and
                            path.startswith(".close.json.tmp-")):
                        temporary_descriptor = descriptor
                    return descriptor

                def fail_primary_fsync(descriptor):
                    nonlocal fsync_failed
                    if (descriptor == temporary_descriptor and
                            not fsync_failed):
                        fsync_failed = True
                        raise OSError(
                            errno.ENOSPC, "injected temporary fsync failure")
                    return real_fsync(descriptor)

                def fail_close_after_release(descriptor):
                    nonlocal close_failed
                    real_close(descriptor)
                    if (descriptor == temporary_descriptor and
                            not close_failed):
                        close_failed = True
                        raise OSError(
                            errno.EIO, "injected temporary close failure")

                descriptors_before = len(os.listdir("/proc/self/fd"))
                with mock.patch.object(
                        runner.os, "open", side_effect=tracked_open), \
                     mock.patch.object(
                         runner.os, "fsync",
                         side_effect=fail_primary_fsync), \
                     mock.patch.object(
                         runner.os, "close",
                         side_effect=fail_close_after_release), \
                     self.assertRaisesRegex(
                         RuntimeError,
                         "injected temporary fsync failure.*cleanup also "
                         "failed.*injected temporary close failure"):
                    runner.write_json_atomic(
                        output / "close.json", {"value": 1},
                        evidence_directory=evidence)
                self.assertTrue(fsync_failed)
                self.assertTrue(close_failed)
                self.assertEqual(
                    len(os.listdir("/proc/self/fd")), descriptors_before)
                self.assertFalse((output / "close.json").exists())
                self.assertEqual(
                    list(output.glob(".close.json.tmp-*")), [])
                runner.write_json_atomic(
                    output / "close.json", {"value": 2},
                    evidence_directory=evidence)

                temporary_descriptor = None
                fsync_failed = False
                unlink_failed = False
                real_unlink = os.unlink

                def tracked_unlink(path, *arguments, **keywords):
                    nonlocal unlink_failed
                    if (isinstance(path, str) and
                            path.startswith(".unlink.json.tmp-") and
                            not unlink_failed):
                        unlink_failed = True
                        raise OSError(
                            errno.EIO, "injected temporary unlink failure")
                    return real_unlink(path, *arguments, **keywords)

                def tracked_unlink_open(
                        path, flags, *arguments, **keywords):
                    nonlocal temporary_descriptor
                    descriptor = real_open(
                        path, flags, *arguments, **keywords)
                    if (isinstance(path, str) and
                            path.startswith(".unlink.json.tmp-")):
                        temporary_descriptor = descriptor
                    return descriptor

                descriptors_before = len(os.listdir("/proc/self/fd"))
                with mock.patch.object(
                        runner.os, "open",
                        side_effect=tracked_unlink_open), \
                     mock.patch.object(
                         runner.os, "fsync",
                         side_effect=fail_primary_fsync), \
                     mock.patch.object(
                         runner.os, "unlink",
                         side_effect=tracked_unlink), \
                     self.assertRaisesRegex(
                         RuntimeError,
                         "injected temporary fsync failure.*cleanup also "
                         "failed.*injected temporary unlink failure"):
                    runner.write_json_atomic(
                        output / "unlink.json", {"value": 3},
                        evidence_directory=evidence)
                self.assertTrue(fsync_failed)
                self.assertTrue(unlink_failed)
                self.assertEqual(
                    len(os.listdir("/proc/self/fd")), descriptors_before)
                self.assertFalse((output / "unlink.json").exists())
                leftovers = list(output.glob(".unlink.json.tmp-*"))
                self.assertEqual(len(leftovers), 1)
                leftovers[0].unlink()
                runner.write_json_atomic(
                    output / "unlink.json", {"value": 4},
                    evidence_directory=evidence)
            finally:
                evidence.close()

    def test_held_evidence_reads_require_exact_permissions(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-held-read-mode-") as directory:
            output = Path(directory) / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            try:
                runner.write_json_atomic(
                    output / "manifest.json", {"value": 1},
                    evidence_directory=evidence)
                self.assertEqual(
                    runner.strict_json_bytes(
                        evidence.read_optional(
                            "manifest.json", 1024, "held manifest"),
                        "held manifest"),
                    {"value": 1})
                (output / "manifest.json").chmod(0o666)
                with self.assertRaisesRegex(
                        RuntimeError, "exact owner-only"):
                    evidence.read_optional(
                        "manifest.json", 1024, "permissive manifest")
                (output / "manifest.json").chmod(0o600)
                (output / "cells").mkdir(mode=0o700)
                (output / "cells").chmod(0o755)
                with self.assertRaisesRegex(
                        RuntimeError, "unsafe"):
                    evidence.read_optional(
                        "cells/missing.json", 1024,
                        "permissive child directory")
            finally:
                evidence.close()

    def test_held_evidence_read_cannot_follow_root_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-held-read-aba-") as directory:
            root = Path(directory)
            output = root / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            try:
                runner.write_json_atomic(
                    output / "manifest.json", {"source": "A"},
                    evidence_directory=evidence)
                retained = root / "retained-A"
                replacement = root / "replacement-B"

                def swap_root_aba() -> None:
                    output.rename(retained)
                    output.mkdir(mode=0o700)
                    (output / "manifest.json").write_text(
                        '{"source":"B"}\n', encoding="utf-8")
                    (output / "manifest.json").chmod(0o600)
                    output.rename(replacement)
                    retained.rename(output)

                payload = evidence.read_optional(
                    "manifest.json", 1024, "ABA-held manifest",
                    mutation_hook=swap_root_aba)
                self.assertEqual(
                    runner.strict_json_bytes(payload, "ABA-held manifest"),
                    {"source": "A"})
                self.assertEqual(
                    json.loads(
                        (replacement / "manifest.json").read_text(
                            encoding="utf-8")),
                    {"source": "B"})
            finally:
                evidence.close()

    def test_atomic_output_rejects_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-symlink-") as directory:
            root = Path(directory)
            victim = root / "victim"
            victim.mkdir(mode=0o700)
            redirected = root / "redirected"
            redirected.symlink_to(victim.name)
            with self.assertRaises(RuntimeError):
                runner.AtomicEvidenceDirectory.open_or_create(
                    redirected / "output")
            self.assertEqual(list(victim.iterdir()), [])

    def test_atomic_output_parallel_child_creation_is_race_free(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-output-parallel-") as directory:
            output = Path(directory) / "output"
            evidence = runner.AtomicEvidenceDirectory.open_or_create(output)
            barrier = threading.Barrier(16)
            failures = []

            def publish(index: int) -> None:
                try:
                    barrier.wait(timeout=2.0)
                    runner.write_json_atomic(
                        output / "cells" / f"{index}.json",
                        {"index": index}, evidence_directory=evidence)
                except BaseException as error:
                    failures.append(error)

            workers = [
                threading.Thread(target=publish, args=(index,))
                for index in range(16)
            ]
            try:
                for worker in workers:
                    worker.start()
                for worker in workers:
                    worker.join(timeout=5.0)
                self.assertFalse(
                    any(worker.is_alive() for worker in workers))
                self.assertEqual(failures, [])
                self.assertEqual(
                    len(list((output / "cells").glob("*.json"))), 16)
            finally:
                evidence.close()

    def test_1000_digit_json_number_has_bounded_evidence_error(self) -> None:
        parsed = runner.strict_json_bytes(
            b'{"value":' + b"9" * 1000 + b"}",
            "1000-digit fixture")
        self.assertIs(type(parsed["value"]), int)
        expect_rejected(
            self,
            lambda: runner.exact_positive_finite(
                parsed["value"], "1000-digit fixture"),
            "outside the finite float range")

    def test_persisted_json_loader_rejects_ambiguous_or_nonfinite_bytes(
            self) -> None:
        cases = (
            ("duplicate-schema",
             b'{"schema":"leopard2-all-k-gap-manifest/v4",'
             b'"schema":"leopard2-all-k-gap-manifest/v4"}\n',
             "duplicate"),
            ("duplicate-valid", b'{"valid":true,"valid":true}\n',
             "duplicate"),
            ("nan", b'{"valid":NaN}\n', "non-finite"),
            ("infinity", b'{"valid":Infinity}\n', "non-finite"),
            ("positive-overflow", b'{"valid":1e9999}\n', "non-finite"),
            ("negative-overflow", b'{"valid":-1e9999}\n', "non-finite"),
            ("invalid-utf8", b'{"valid":"\xff"}\n', "strict finite JSON"),
            ("trailing-text", b'{"valid":true} trailing\n',
             "strict finite JSON"),
            ("multiple-documents", b'{"valid":true}\n{"valid":true}\n',
             "strict finite JSON"),
        )
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-strict-json-") as directory:
            path = Path(directory) / "evidence.json"
            for name, payload, message in cases:
                with self.subTest(name=name):
                    path.write_bytes(payload)
                    path.chmod(0o600)
                    expect_rejected(
                        self,
                        lambda: runner.load_strict_json(
                            path, "persisted fixture", 4096),
                        message)

    def test_persisted_json_loader_rejects_oversized_sparse_file(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-large-json-") as directory:
            path = Path(directory) / "evidence.json"
            with path.open("wb") as output:
                output.truncate(4097)
            path.chmod(0o600)
            expect_rejected(
                self,
                lambda: runner.load_strict_json(
                    path, "persisted fixture", 4096),
                "exceeds 4096 bytes")

    def test_persisted_json_loader_requires_exact_owner_only_modes(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-json-modes-") as directory:
            root = Path(directory)
            path = root / "evidence.json"
            path.write_bytes(b'{"valid":true}\n')
            path.chmod(0o644)
            expect_rejected(
                self,
                lambda: runner.load_strict_json(
                    path, "permissive persisted fixture", 4096),
                "exact owner-only evidence")

            path.chmod(0o600)
            root.chmod(0o755)
            try:
                expect_rejected(
                    self,
                    lambda: runner.load_strict_json(
                        path, "permissive persisted fixture", 4096),
                    "directory is not exact owner-only evidence")
            finally:
                root.chmod(0o700)

    def test_persisted_json_loader_binds_parent_and_stable_file_snapshot(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-allk-json-race-") as directory:
            root = Path(directory)
            victim = root / "victim"
            victim.mkdir(mode=0o700)
            victim_file = victim / "evidence.json"
            victim_file.write_bytes(b'{"victim":true}\n')
            redirected = root / "redirected"
            redirected.symlink_to(victim.name)
            expect_rejected(
                self,
                lambda: runner.load_strict_json(
                    redirected / "evidence.json",
                    "redirected persisted fixture", 4096),
                "cannot open")
            self.assertEqual(
                victim_file.read_bytes(), b'{"victim":true}\n')

            path = root / "mutating.json"
            path.write_bytes(
                b'{"payload":"' + b"A" * (2 * 1024 * 1024) + b'"}\n')
            path.chmod(0o600)

            def mutate_during_read() -> None:
                descriptor = os.open(
                    path, os.O_WRONLY | getattr(os, "O_CLOEXEC", 0))
                try:
                    os.pwrite(descriptor, b"B", 32)
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)

            expect_rejected(
                self,
                lambda: runner.load_strict_json(
                    path, "mutating persisted fixture", 4 * 1024 * 1024,
                    mutation_hook=mutate_during_read),
                "changed while it was being read")


class ProductionBuildClosureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = ProductionBuildFixture()

    def tearDown(self) -> None:
        self.fixture.close()

    def provenance(self):
        return runner.candidate_build_provenance(
            self.fixture.build, self.fixture.source,
            self.fixture.executable, "bench_leopard2")

    def test_exact_release_pure_avx2_closure_is_accepted(self) -> None:
        result = self.provenance()
        self.assertEqual(result["schema"],
                         runner.PRODUCTION_BUILD_CLOSURE_SCHEMA)
        self.assertEqual(
            result["validated_cache"][
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"],
            runner.BENCHMARK_BUILD_CONFIGURATION_SCHEMA)
        self.assertEqual(
            result["validated_cache"][
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"], "ON")
        self.assertEqual(
            result["validated_cache"][
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT"], "ON")
        self.assertEqual(
            result["validated_cache"][
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"], "ON")
        self.assertEqual(result["archive_members"], [
            path.name for path in self.fixture.archive_objects])
        self.assertTrue(
            (self.fixture.build / self.fixture.c_test_object).is_file())
        self.assertNotIn(
            self.fixture.c_test_object, self.fixture.archive_objects)
        self.assertNotIn(
            self.fixture.c_test_object.name, result["archive_members"])
        profiles = {record["flag_profile"] for record in
                    result["source_object_compile_closure"]}
        self.assertEqual(profiles, {
            "portable-core", "ssse3-no-avx", "avx2-no-avx512"})

    def test_timed_snapshot_rejects_restored_symlink_swap(self) -> None:
        first = self.fixture.build / "candidate-a"
        second = self.fixture.build / "candidate-b"
        exact = self.provenance()
        proof = runner.compare_reproducible_builds(
            exact, copy.deepcopy(exact))

        with mock.patch.object(
                runner, "validate_reproducible_build_proof"):
            source, descriptor, snapshot = \
                runner.snapshot_reproduced_executable(
                    self.fixture.executable, exact, proof)
            try:
                self.assertEqual(source, exact["executable"])
                self.assertEqual(snapshot["sha256"],
                                 proof["executable_sha256"])
            finally:
                os.close(descriptor)

            shutil.copy2("/bin/false", second)
            self.fixture.executable.rename(first)
            self.fixture.executable.symlink_to(second.name)
            try:
                expect_rejected(
                    self,
                    lambda: runner.snapshot_reproduced_executable(
                        self.fixture.executable, exact, proof),
                    "snapshot source differs")
            finally:
                self.fixture.executable.unlink()
                first.rename(self.fixture.executable)

        # Restoring the original inode necessarily changes its ctime, but the
        # validated source/object/archive/link closure and executable bytes
        # reproduce exactly.  This demonstrates why end-of-campaign byte
        # recapture alone cannot catch the injected target swap.
        restored = self.provenance()
        self.assertEqual(
            runner.compare_reproducible_builds(exact, restored), proof)

    def test_timed_snapshot_rejects_incomplete_replay_proof(self) -> None:
        exact = self.provenance()
        legacy_core = runner.compare_reproducible_builds(
            exact, copy.deepcopy(exact))
        expect_rejected(
            self,
            lambda: runner.snapshot_reproduced_executable(
                self.fixture.executable, exact, legacy_core),
            "no valid reproducible executable proof")

    def test_c_compile_entry_requires_cached_c_driver(self) -> None:
        entry = next(item for item in self.fixture.entries
                     if item["file"].endswith(".c"))
        entry["arguments"][0] = str(
            Path("/usr/bin/c++").resolve(strict=True))
        self.fixture._write_commands()
        expect_rejected(
            self, self.provenance, "another command driver")

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
        AllKIdentityTests.bind_raw_output(record)
        runner.validate_source_attestation(record, source, snapshot)
        exact = self.provenance()
        object_path = self.fixture.build / next(
            path for path in self.fixture.archive_objects
            if path.name == "Leopard2BackendAVX2.cpp.o")
        old_bytes = object_path.read_bytes()
        object_path.write_bytes(
            bytes([old_bytes[0] ^ 1]) + old_bytes[1:])
        self.fixture._relink()
        self.assertEqual(runner.git_identity(
            self.fixture.source, source["head"])["tracked_status"], "clean")
        stale = self.provenance()
        self.assertNotEqual(
            exact["source_object_compile_closure"],
            stale["source_object_compile_closure"])
        expect_rejected(
            self,
            lambda: runner.compare_reproducible_builds(exact, stale),
            "object bytes differ")

    def test_avx512_cache_or_compile_flags_are_rejected(self) -> None:
        self.fixture._write_cache(
            **{"LEO2_FLAG_MAVX512F:UNINITIALIZED": "TRUE"})
        expect_rejected(self, self.provenance, "did not disable")
        self.fixture._write_cache()
        entry = next(item for item in self.fixture.entries
                     if item["file"].endswith("Leopard2BackendAVX2.cpp"))
        entry["arguments"].insert(
            entry["arguments"].index("-o"), "-mavx512f")
        self.fixture._write_commands()
        expect_rejected(
            self, self.provenance,
            "non-canonical or indirect compile option")

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
        general_one_loss = \
            "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1"
        entry["file"] = wrong
        entry["arguments"][entry["arguments"].index("-c") + 1] = wrong
        entry["arguments"].insert(
            entry["arguments"].index("-o"), general_one_loss)
        self.fixture._write_commands()
        expect_rejected(self, self.provenance, "source closure differs")
        entry["file"] = str(
            (self.fixture.source / "Leopard2Plan.cpp").resolve())
        entry["arguments"][entry["arguments"].index("-c") + 1] = entry["file"]
        entry["arguments"].remove(general_one_loss)
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

    def test_clean_replay_preserves_real_cmake_cache_types(self) -> None:
        exact = self.provenance()
        source = Path(__file__).resolve().parents[3]
        with tempfile.TemporaryDirectory(
                prefix="leo2-replay-cache-types-") as directory:
            build = Path(directory) / "build"
            configure = provenance_module._reproducible_configure_argv(
                source, build, exact["validated_cache"])
            completed = subprocess.run(
                configure, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=provenance_module.GIT_ENVIRONMENT,
                timeout=300, check=False)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode(errors="replace"))
            cache_bytes = (build / "CMakeCache.txt").read_bytes()
            parsed = provenance_module.parse_cmake_cache(cache_bytes)
            self.assertEqual(parsed["CMAKE_CXX_COMPILER"], "/usr/bin/c++")
            self.assertEqual(parsed["LEO2_FLAG_MAVX2"], "1")
            self.assertEqual(parsed["LEO2_FLAG_MNO_AVX512F"], "1")
            cache_text = cache_bytes.decode("utf-8", errors="strict")
            for expected in (
                    "CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\n",
                    "LEO2_FLAG_MAVX2:INTERNAL=1\n",
                    "LEO2_FLAG_MNO_AVX512F:INTERNAL=1\n"):
                with self.subTest(expected=expected):
                    self.assertIn(expected, cache_text)


if __name__ == "__main__":
    unittest.main()
