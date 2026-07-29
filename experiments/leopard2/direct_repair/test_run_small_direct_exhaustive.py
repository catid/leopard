#!/usr/bin/env python3
"""Adversarial contract tests for run_small_direct_exhaustive.py."""

from __future__ import annotations

import copy
import fcntl
import importlib.util
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


RUNNER_PATH = Path(__file__).with_name("run_small_direct_exhaustive.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_test_small_direct_exhaustive_runner", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import exhaustive direct-repair runner")
RUNNER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def write(path: Path, payload: bytes) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return RUNNER.file_identity(path)


def rich(identity: dict) -> dict:
    status = Path(identity["path"]).stat()
    return {
        **identity,
        "device": status.st_dev,
        "inode": status.st_ino,
        "mode": status.st_mode,
        "mtime_ns": status.st_mtime_ns,
        "ctime_ns": status.st_ctime_ns,
    }


def make_build_record(origin: Path, binary_identity: dict, mode: int) -> dict:
    source_root = origin / "source"
    source_files = {}
    for index, relative in enumerate(RUNNER.BOUND_SOURCE_FILES):
        source_files[relative] = write(
            source_root / relative,
            ("bound-%d-%s\n" % (index, relative)).encode())
    source = {
        "root": str(source_root.resolve()),
        "head": "1" * 40, "head_tree": "2" * 40,
        "status_short": [], "files": source_files,
    }
    source["digest"] = RUNNER.record_digest(source)

    definition = RUNNER.MODE_DEFINITIONS[mode]
    archive_source = source_files["leopard2.cpp"]
    unselected_source = source_files["leopard.cpp"]
    benchmark_source = source_files["bench/leopard2/benchmark.cpp"]
    archive_object = write(origin / "build/leopard2.cpp.o", b"archive object")
    unselected_object = write(
        origin / "build/leopard.cpp.o", b"unselected archive object")
    benchmark_object = write(
        origin / "build/benchmark.cpp.o", b"benchmark object")
    archive = write(origin / "build/libleopard.a", b"archive")
    benchmark = write(origin / "build/bench_leopard2", b"benchmark")
    production = {
        "schema": "leopard2-production-build-closure/v1",
        "archive": archive, "executable": benchmark,
        "source_object_compile_closure": [
            {
                "role": "archive", "source": archive_source,
                "object": archive_object,
                "compile_entry": {
                    "arguments": ["c++", "-O3", definition],
                },
                "flag_profile": "portable-core",
            },
            {
                "role": "archive", "source": unselected_source,
                "object": unselected_object,
                "compile_entry": {"arguments": ["c++", "-O3"]},
                "flag_profile": "portable-core",
            },
            {
                "role": "benchmark", "source": benchmark_source,
                "object": benchmark_object,
                "compile_entry": {"arguments": ["c++", "-O3"]},
                "flag_profile": "portable-core",
            },
        ],
        "build_root": str((origin / "build").resolve()),
        "compiler": RUNNER.file_identity(
            Path("/usr/bin/c++").resolve(strict=True)),
        "compiler_version_sha256": "3" * 64,
        "archiver": RUNNER.file_identity(
            Path("/usr/bin/ar").resolve(strict=True)),
        "archive_members": [
            Path(archive_object["path"]).name,
            Path(unselected_object["path"]).name,
        ],
        "archive_member_identities": [
            {
                "member": Path(archive_object["path"]).name,
                "size": archive_object["size"],
                "sha256": archive_object["sha256"],
            },
            {
                "member": Path(unselected_object["path"]).name,
                "size": unselected_object["size"],
                "sha256": unselected_object["sha256"],
            },
        ],
    }

    target_records = []
    for name in ("direct_oracle.cpp", "test_small_direct_exhaustive.cpp"):
        source_identity = source_files["tests/leopard2/" + name]
        object_identity = write(
            origin / "build" / (name + ".o"), name.encode())
        target_records.append({
            "source": rich(source_identity), "object": rich(object_identity),
            "compile_entry": {
                "directory": str((origin / "build").resolve()),
                "file": source_identity["path"],
                "output": object_identity["path"],
                "arguments": (
                    ["c++", "-O3", definition]
                    if name == "test_small_direct_exhaustive.cpp"
                    else ["c++", "-O3"]),
            },
            "flag_profile": "portable-core",
        })
    metadata = rich(write(origin / "build/CMakeCache.txt", b"cache"))
    commands = rich(write(origin / "build/compile_commands.json", b"[]"))
    recipe = rich(write(origin / "build/link.txt", b"link"))
    sidecar = write(origin / "build/configuration.txt", b"configuration")
    helper = source_files["tools/leopard2_build_provenance.py"]
    reader = source_files["tools/leopard2_direct_encode_crossover.py"]
    attestation_entries = {
        name: "" for name in RUNNER.BUILD_CONFIGURATION_VARIABLES
    }
    attestation_entries[RUNNER.MODE_CACHE_VARIABLE] = str(mode)
    for name in RUNNER.FORBIDDEN_EXPERIMENT_MACROS:
        attestation_entries[name] = "OFF"
    controlled = {
        name: "" for name in RUNNER.CONTROLLED_CONFIGURATION_KEYS
    }
    controlled.update({
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_AR": "/usr/bin/ar",
        "CMAKE_COMMAND": "/usr/bin/cmake",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_C_COMPILER": "/usr/bin/cc",
        "CMAKE_CXX_COMPILER": "/usr/bin/c++",
        "CMAKE_CXX_FLAGS": "",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
        "CMAKE_LINKER": "/usr/bin/ld",
        "CMAKE_MAKE_PROGRAM": "/usr/bin/gmake",
        "CMAKE_RANLIB": "/usr/bin/ranlib",
        "ENABLE_OPENMP": "ON",
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
        "LEO2_BACKEND_VARIANT": "avx2",
        "LEO2_BUILD_ALLK_DIAGNOSTIC": "OFF",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "ON",
        "LEO2_ENABLE_CUDA": "OFF",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
        RUNNER.MODE_CACHE_VARIABLE: str(mode),
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_FLAG_MAVX512F": "FALSE",
        "LEO2_FLAG_MAVX512BW": "FALSE",
        "LEO2_FLAG_MAVX512VL": "FALSE",
    })
    production["validated_cache"] = dict(controlled)
    result = {
        "schema": RUNNER.BUILD_SCHEMA, "mode": mode,
        "source": source, "production_closure": production,
        "controlled_configuration": controlled,
        "effective_configuration": {
            "binding_scope": "synthetic exhaustive-runner contract fixture",
            "build_root": str((origin / "build").resolve()),
            "cmake_cache": metadata["path"],
            "cmake_cache_sha256": metadata["sha256"],
            "effective_configuration_attestation": {
                "schema":
                    "leopard2-benchmark-build-configuration-attestation/v2",
                "entries": attestation_entries,
                "path": sidecar["path"],
                "sha256":
                    RUNNER.build_configuration_digest(attestation_entries),
            },
            "entries": {
                RUNNER.MODE_CACHE_VARIABLE: str(mode),
                **{
                    name: "OFF"
                    for name in RUNNER.FORBIDDEN_EXPERIMENT_MACROS
                },
            },
            "executable": {
                "mtime_ns": Path(benchmark["path"]).stat().st_mtime_ns,
                "path": benchmark["path"], "sha256": benchmark["sha256"],
                "size": benchmark["size"],
            },
            "extra_file_sha256": {},
        },
        "cache": metadata, "compile_commands": commands,
        "target_link_recipe": recipe,
        "target_link_arguments": [
            "/usr/bin/c++", "<BUILD>/direct_oracle.cpp.o",
            "<BUILD>/test_small_direct_exhaustive.cpp.o",
            "<BUILD>/libleopard.a",
        ],
        "target_compile_closure": target_records,
        "mode_independent_object_sha256": next(
            record["object"]["sha256"] for record in target_records
            if Path(record["source"]["path"]).name ==
                "direct_oracle.cpp"),
        "exhaustive_executable": rich(binary_identity),
        "configuration_sidecar": sidecar,
        "provenance_helper": helper,
        "configuration_reader": reader,
    }
    material = RUNNER.build_reproducible_material(result)
    proof_tools = {
        label: {
            "configured_path": controlled[key],
            "identity": RUNNER.file_identity(
                Path(controlled[key]).resolve(strict=True)),
        }
        for label, key in {
            "cmake": "CMAKE_COMMAND",
            "make": "CMAKE_MAKE_PROGRAM",
            "linker": "CMAKE_LINKER",
            "ranlib": "CMAKE_RANLIB",
            "archiver": "CMAKE_AR",
        }.items()
    }
    result["reproducible_rebuild"] = {
        "schema": "leopard2-small-direct-clean-rebuild/v1",
        "mode": mode,
        "tools": proof_tools,
        "configure_profile":
            RUNNER.controlled_configure_profile(controlled, mode),
        "configure_environment":
            RUNNER.controlled_configure_environment(controlled),
        "candidate_material": material,
        "rebuilt_material": material,
    }
    result["reproducible_rebuild"]["digest"] = RUNNER.record_digest(
        result["reproducible_rebuild"])
    result["digest"] = RUNNER.record_digest(result)
    RUNNER.validate_build_record(result)
    return result


def make_campaign(root: Path, workers: int = 30, mode: int = 1) -> Path:
    campaign = root / "campaign"
    origin = root / "origin"
    campaign.mkdir(parents=True, mode=0o700)
    origin.mkdir(parents=True)
    binary = origin / RUNNER.TARGET_NAME
    binary_identity = write(binary, b"synthetic exhaustive executable\n")
    os.chmod(binary, 0o555)
    binary_identity = RUNNER.file_identity(binary)
    build = make_build_record(origin, binary_identity, mode)
    artifacts = campaign / "artifacts"
    artifacts.mkdir(mode=0o700)
    frozen = RUNNER.freeze_file(
        binary_identity, artifacts / RUNNER.TARGET_NAME, 0o555, campaign)
    runner_identity = build["source"]["files"][
        "experiments/leopard2/direct_repair/"
        "run_small_direct_exhaustive.py"]
    taskset_identity = RUNNER.file_identity(RUNNER.CANONICAL_TASKSET)
    interpreter_identity = RUNNER.process_executable_identity()
    bundle = RUNNER.freeze_bundle(
        build, (runner_identity, taskset_identity, interpreter_identity),
        campaign)
    expected = [
        dict(value) for value in RUNNER.reconstruct_expected_shards(workers)
    ]
    request = {
        "schema": RUNNER.SCHEMA,
        "source_binary": binary_identity,
        "frozen_binary": frozen["frozen"],
        "build_provenance": build,
        "provenance_bundle": bundle,
        "runner": runner_identity,
        "taskset": taskset_identity,
        "interpreter": interpreter_identity,
        "workers": workers, "shard_count": workers,
        "expected_patterns": RUNNER.EXPECTED_PATTERNS,
        "timeout_seconds_per_shard": 10.0,
        "allowed_cpus": list(range(workers)),
        "assignment": "global_ordinal_mod_shard_count",
        "basis_seed": 0, "expected_shards": expected,
        "child_environment": RUNNER.CHILD_ENV,
        "execution_policy": RUNNER.EXECUTION_POLICY,
        "lock": {
            "schema": "leopard2-global-benchmark-lock/v1",
            "path": str(RUNNER.DEFAULT_LOCK), "device": 1, "inode": 2,
            "uid": os.getuid(), "mode": 0o600, "nlink": 1,
            "mechanism": "fcntl-flock-exclusive",
        },
        "resume_policy": "validated deterministic shard resume",
    }
    request["digest"] = RUNNER.record_digest(request)
    RUNNER.atomic_json(campaign / "request.json", request)

    shards = []
    for index, expected_shard in enumerate(expected):
        directory = campaign / ("shard-%04d" % index)
        directory.mkdir(mode=0o700)
        shard = {"schema": RUNNER.SHARD_SCHEMA, "mode": mode,
                 **expected_shard}
        RUNNER.atomic_json(directory / "stdout.json", shard)
        (directory / "stderr.txt").write_bytes(b"")
        command = [
            taskset_identity["path"], "-c", str(index),
            frozen["frozen"]["path"],
            "--shard-index", str(index), "--shard-count", str(workers),
        ]
        envelope = {
            "shard": shard,
            "stdout": RUNNER.file_identity(directory / "stdout.json"),
            "stderr": RUNNER.file_identity(directory / "stderr.txt"),
            "binary": frozen["frozen"], "cpu": index, "command": command,
            "execution": {
                "policy": RUNNER.EXECUTION_POLICY,
                "interpreter": interpreter_identity,
                "binary_snapshot": RUNNER.sealed_snapshot_identity(
                    frozen["frozen"], 0o500),
                "runner_snapshot": RUNNER.sealed_snapshot_identity(
                    runner_identity, 0o400),
                "taskset_snapshot": RUNNER.sealed_snapshot_identity(
                    taskset_identity, 0o500),
            },
        }
        RUNNER.atomic_json(directory / "result.json", envelope)
        shards.append(shard)
    manifest = {
        "schema": RUNNER.SCHEMA, "request": request, "mode": mode,
        "shard_count": workers,
        "verified_patterns": sum(
            value["assigned_patterns"] for value in shards),
        "recovered_shards": sum(
            value["recovered_shards"] for value in shards),
        "verified_basis_symbols": sum(
            value["verified_basis_symbols"] for value in shards),
        "shards": shards,
    }
    manifest["digest"] = RUNNER.record_digest(manifest)
    RUNNER.atomic_json(campaign / "manifest.json", manifest)
    return campaign / "manifest.json"


def executable_script(path: Path, source: str) -> dict:
    path.write_text("#!/usr/bin/python3\n" + source, encoding="utf-8")
    path.chmod(0o555)
    return RUNNER.file_identity(path)


def run_supervisor_fixture(
        root: Path, executable: Path, timeout: float,
        after_snapshot=None,
        close_parent_watch: bool = False,
        ready_path: Path | None = None) -> tuple[int, bytes, bytes]:
    """Run a fixture through the same sealed/cgroup supervisor as shards."""
    output = root / "campaign"
    output.mkdir(mode=0o700)
    shard = output / "shard-0000"
    shard.mkdir(mode=0o700)
    stdout_path = shard / "stdout.json"
    stderr_path = shard / "stderr.txt"
    binary_identity = RUNNER.file_identity(executable)
    runner_identity = RUNNER.file_identity(RUNNER_PATH)
    taskset_identity = RUNNER.file_identity(RUNNER.CANONICAL_TASKSET)
    interpreter_identity = RUNNER.process_executable_identity()
    binary_fd = -1
    runner_fd = -1
    taskset_fd = -1
    parent_watch_read = -1
    parent_watch_write = -1
    process = None
    campaign_scope = None
    scope = None
    try:
        binary_fd, binary_snapshot = RUNNER.sealed_file_snapshot(
            executable, binary_identity, "test-binary", 0o500)
        runner_fd, unused_runner_snapshot = RUNNER.sealed_file_snapshot(
            RUNNER_PATH, runner_identity, "test-runner", 0o400)
        taskset_fd, taskset_snapshot = RUNNER.sealed_file_snapshot(
            RUNNER.CANONICAL_TASKSET, taskset_identity,
            "test-taskset", 0o500)
        parent_watch_read, parent_watch_write = os.pipe2(
            os.O_CLOEXEC | os.O_NONBLOCK)
        if after_snapshot is not None:
            after_snapshot()
        if not RUNNER.process_executable_identity_matches(
                interpreter_identity):
            raise RuntimeError("fixture interpreter inode changed")
        try:
            parent = RUNNER.cgroup_v2_current()
            campaign_scope = parent / (
                "leopard2-exhaustive-%d-%x" %
                (os.getpid(), time.time_ns()))
            os.mkdir(campaign_scope)
            scope = campaign_scope / "shard-0000"
            os.mkdir(scope)
        except (OSError, RUNNER.EvidenceError) as error:
            raise unittest.SkipTest(
                "delegated cgroup-v2 fixture unavailable: %s" % error)
        command = [
            *RUNNER.supervisor_command_prefix(runner_fd),
            "--scope", str(scope),
            "--parent-cgroup", str(parent),
            "--coordinator-pid", str(os.getpid()),
            "--binary-fd", str(binary_fd),
            "--binary-size", str(binary_snapshot["size"]),
            "--binary-sha256", binary_snapshot["sha256"],
            "--taskset-fd", str(taskset_fd),
            "--taskset-size", str(taskset_snapshot["size"]),
            "--taskset-sha256", taskset_snapshot["sha256"],
            "--cpu", str(min(os.sched_getaffinity(0))),
            "--shard-index", "0", "--shard-count", "1",
            "--stdout", str(stdout_path),
            "--stderr", str(stderr_path),
            "--output-root", str(output),
            "--parent-watch-fd", str(parent_watch_read),
            "--timeout", str(timeout),
        ]
        process = RUNNER.launch_current_interpreter(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, start_new_session=True,
            pass_fds=(
                runner_fd, binary_fd, taskset_fd, parent_watch_read),
            env=RUNNER.CHILD_ENV)
        if ready_path is not None:
            ready_path.write_text(json.dumps({
                "campaign_scope": str(campaign_scope),
                "scope": str(scope),
                "supervisor_pid": process.pid,
            }), encoding="utf-8")
        if close_parent_watch:
            os.close(parent_watch_write)
            parent_watch_write = -1
        return_code = process.wait(
            timeout=timeout + RUNNER.SUPERVISOR_GRACE_SECONDS)
        if scope.exists():
            RUNNER.terminate_supervisor(process, scope)
        if not RUNNER.process_executable_identity_matches(
                interpreter_identity):
            raise RuntimeError("fixture interpreter inode changed")
        return (
            return_code,
            stdout_path.read_bytes() if stdout_path.exists() else b"",
            stderr_path.read_bytes() if stderr_path.exists() else b"",
        )
    finally:
        if process is not None and process.poll() is None and scope is not None:
            RUNNER.terminate_supervisor(process, scope)
        if scope is not None and scope.exists():
            RUNNER.signal_cgroup(scope)
            RUNNER.remove_empty_cgroup(scope)
        if campaign_scope is not None and campaign_scope.exists():
            campaign_scope.rmdir()
        for descriptor in (
                parent_watch_write, parent_watch_read,
                taskset_fd, runner_fd, binary_fd):
            if descriptor >= 0:
                os.close(descriptor)


class ExhaustiveRunnerTests(unittest.TestCase):
    def test_process_executable_identity_rejects_unrelated_path(self) -> None:
        with self.assertRaisesRegex(
                RUNNER.EvidenceError,
                "retained interpreter pathname differs"):
            RUNNER.process_executable_identity("/usr/bin/false")

    def test_process_executable_identity_record_is_canonical(self) -> None:
        identity = RUNNER.process_executable_identity()
        for field, value in (
                ("path", "/tmp/.."),
                ("device", -1),
                ("inode", -1),
                ("mode", 0),
                ("size", -1),
                ("sha256", "A" * 64)):
            malformed = dict(identity)
            malformed[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(
                    RUNNER.EvidenceError,
                    "retained interpreter identity is malformed"):
                RUNNER.process_executable_identity_matches(malformed)

    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory(
            prefix="leo2-exhaustive-evidence-")
        cls.root = Path(cls.temporary.name)
        cls.manifest = make_campaign(cls.root)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_thirty_shards_cover_exact_matrix_once(self) -> None:
        records = list(RUNNER.reconstruct_expected_shards(30))
        self.assertEqual(len(records), 30)
        self.assertEqual(
            sum(value["assigned_patterns"] for value in records),
            1_982_812)
        self.assertLessEqual(
            max(value["assigned_patterns"] for value in records) -
            min(value["assigned_patterns"] for value in records), 1)
        self.assertEqual(
            RUNNER.sha256_bytes(RUNNER.canonical_bytes(records)),
            "6bab1c0997e02000775bf9241b956d9c2689bf2b16c85911679c56ddf3c071f1")
        self.assertEqual(records[0]["digest_fnv1a64"],
                         "ebf0766317565412")
        self.assertEqual(records[-1]["ordinal_digest_fnv1a64"],
                         "bffe450894fab3b1")

    def test_exact_count_and_both_digests_are_required(self) -> None:
        expected = dict(RUNNER.reconstruct_expected_shards(30)[0])
        value = {"schema": RUNNER.SHARD_SCHEMA, "mode": 1, **expected}
        RUNNER.validate_shard(value, 0, 30, expected)
        for key in ("assigned_patterns", "recovered_shards",
                    "verified_basis_symbols"):
            forged = dict(value)
            forged[key] += 1
            with self.subTest(key=key), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "identity/counts/digests"):
                RUNNER.validate_shard(forged, 0, 30, expected)
        for key in ("digest_fnv1a64", "ordinal_digest_fnv1a64"):
            forged = dict(value)
            forged[key] = "0" * 16
            with self.subTest(key=key), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "identity/counts/digests"):
                RUNNER.validate_shard(forged, 0, 30, expected)

    def test_offline_manifest_reconstructs_all_artifacts(self) -> None:
        self.assertEqual(
            RUNNER.verify_manifest(self.manifest, True), 0)

    def test_offline_verifier_is_identical_under_python_optimized(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-O", str(RUNNER_PATH), "--verify",
             str(self.manifest), "--no-current-input-check"],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        self.assertEqual(
            completed.returncode, 0,
            completed.stderr.decode(errors="replace"))
        self.assertIn(b"verified 1982812 exact coefficient matrices",
                      completed.stdout)

    def test_execution_policy_is_exact_and_tamper_evident(self) -> None:
        request = copy.deepcopy(
            RUNNER.load_json(self.manifest)["request"])
        request["execution_policy"]["containment"] = "process-group-only"
        request["digest"] = RUNNER.record_digest(request)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "request policy"):
            RUNNER.validate_request(request, self.manifest.parent)

        result_path = self.manifest.parent / "shard-0000/result.json"
        saved = result_path.read_bytes()
        result = RUNNER.decode_json_bytes(saved, "fixture result")
        result["execution"]["binary_snapshot"]["sha256"] = "0" * 64
        result_path.write_bytes(
            RUNNER.canonical_bytes(result))
        try:
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "execution identity"):
                RUNNER.verify_manifest(self.manifest, True)
        finally:
            result_path.write_bytes(saved)

    def test_sealed_snapshot_executes_original_after_path_swap(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-sealed-exec-") as directory:
            root = Path(directory)
            executable = root / "fixture"
            replacement = root / "replacement"
            executable_script(
                executable, 'import os\nos.write(1, b"sealed-A\\n")\n')
            executable_script(
                replacement, 'import os\nos.write(1, b"swapped-B\\n")\n')

            def swap() -> None:
                os.replace(replacement, executable)

            return_code, stdout, stderr = run_supervisor_fixture(
                root, executable, 5.0, swap)
            self.assertEqual(return_code, 0, stderr.decode(errors="replace"))
            self.assertEqual(stdout, b"sealed-A\n")

    def test_supervisor_uses_current_interpreter_inode_after_path_swap(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-interpreter-swap-") as directory:
            root = Path(directory)
            copied_interpreter = root / "python-copy"
            replacement = root / "replacement"
            shutil.copy2(
                Path("/proc/self/exe").resolve(strict=True),
                copied_interpreter)
            shutil.copy2(Path("/usr/bin/false"), replacement)
            copied_interpreter.chmod(0o755)
            replacement.chmod(0o755)
            driver = root / "driver.py"
            driver.write_text("""
import importlib.util
import os
import pathlib
import sys

test_path = pathlib.Path({test_path!r})
spec = importlib.util.spec_from_file_location("copied_runner_test", test_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
root = pathlib.Path({root!r})
fixture = root / "fixture"
module.executable_script(
    fixture, 'import os\\nos.write(1, b"copied-inode\\\\n")\\n')
replacement = root / "replacement"
def swap_interpreter():
    os.replace(replacement, pathlib.Path(sys.executable))
return_code, stdout, stderr = module.run_supervisor_fixture(
    root, fixture, 5.0, swap_interpreter)
if return_code != 0 or stdout != b"copied-inode\\n":
    raise RuntimeError((return_code, stdout, stderr))
""".format(test_path=str(Path(__file__).resolve()), root=str(root)),
                              encoding="utf-8")
            completed = subprocess.run(
                [str(copied_interpreter), str(driver)],
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False,
                env={
                    **os.environ,
                    "PYTHONHOME": sys.base_prefix,
                    "PYTHONDONTWRITEBYTECODE": "1",
                },
                timeout=RUNNER.SUPERVISOR_GRACE_SECONDS + 10)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode(errors="replace"))

    def test_supervisor_startup_isolated_from_late_sitecustomize(
            self) -> None:
        discovery_argv0 = str(
            Path(sys.base_prefix) / "bin" /
            ("python%d.%d%s" % (
                sys.version_info.major, sys.version_info.minor,
                sys.abiflags)))
        self.assertEqual(
            RUNNER.supervisor_command_prefix(17),
            [discovery_argv0, "-I", "-S", "-B",
             "/proc/self/fd/17", "--internal-shard-supervisor"])
        self.assertEqual(
            RUNNER.EXECUTION_POLICY["supervisor_python_flags"],
            ["-I", "-S", "-B"])
        with mock.patch.object(
                RUNNER.subprocess, "Popen") as popen:
            RUNNER.launch_current_interpreter(
                [discovery_argv0, "-I", "-S", "-B", "-c", "pass"],
                stdin=subprocess.DEVNULL)
        self.assertEqual(
            popen.call_args.kwargs["executable"], "/proc/self/exe")
        with tempfile.TemporaryDirectory(
                prefix="leo2-isolated-supervisor-") as directory:
            root = Path(directory)
            marker = root / "sitecustomize-loaded"
            python_path = root / "python-path"
            user_base = root / "user-base"
            user_site = user_base / "lib" / (
                "python%d.%d" %
                (sys.version_info.major, sys.version_info.minor)
            ) / "site-packages"
            payload = (
                "from pathlib import Path\n"
                "Path(%r).write_text('loaded', encoding='ascii')\n" %
                str(marker)
            )

            def inject_after_snapshot() -> None:
                python_path.mkdir(parents=True)
                user_site.mkdir(parents=True)
                (python_path / "sitecustomize.py").write_text(
                    payload, encoding="utf-8")
                (user_site / "sitecustomize.py").write_text(
                    payload, encoding="utf-8")

            environment = {
                **RUNNER.CHILD_ENV,
                "PYTHONPATH": str(python_path),
                "PYTHONUSERBASE": str(user_base),
            }
            with mock.patch.object(RUNNER, "CHILD_ENV", environment):
                return_code, stdout, stderr = run_supervisor_fixture(
                    root, Path("/usr/bin/true"), 5.0,
                    inject_after_snapshot)
            self.assertEqual(
                return_code, 0, stderr.decode(errors="replace"))
            self.assertEqual(stdout, b"")
            self.assertFalse(
                marker.exists(),
                "isolated supervisor imported late sitecustomize")

    def test_guardian_exit_before_readiness_is_rejected(self) -> None:
        read_descriptor, write_descriptor = os.pipe2(
            os.O_CLOEXEC | os.O_NONBLOCK)
        process = subprocess.Popen(
            ["/usr/bin/false"], stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            os.close(write_descriptor)
            write_descriptor = -1
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError,
                    "did not establish lifecycle ownership"):
                RUNNER.wait_guardian_ready(
                    process, read_descriptor, 1.0)
        finally:
            process.wait(timeout=1.0)
            os.close(read_descriptor)
            if write_descriptor >= 0:
                os.close(write_descriptor)

    def test_setsided_double_fork_descendants_are_killed_and_reaped(
            self) -> None:
        for mode, expected_return, timeout in (
                ("success", 0, 5.0),
                ("nonzero", 125, 5.0),
                ("timeout", 124, 0.2),
                ("escaped-cgroup", 0, 5.0)):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                    prefix="leo2-descendant-%s-" % mode) as directory:
                root = Path(directory)
                executable = root / "fixture"
                pid_file = root / "descendant.pid"
                migration_file = root / "descendant.cgroup"
                marker = root / "escaped.marker"
                source = """
import os
import pathlib
import sys
import time

pid_file = {pid_file!r}
migration_file = {migration_file!r}
marker = {marker!r}
mode = {mode!r}
child = os.fork()
if child == 0:
    os.setsid()
    grandchild = os.fork()
    if grandchild != 0:
        os._exit(0)
    if mode == "escaped-cgroup":
        membership = next(
            line.split(":", 2)[2]
            for line in pathlib.Path("/proc/self/cgroup").read_text(
                encoding="ascii").splitlines()
            if line.startswith("0::"))
        scope = pathlib.Path("/sys/fs/cgroup") / membership.lstrip("/")
        with (scope.parent.parent / "cgroup.procs").open(
                "w", encoding="ascii") as stream:
            stream.write("0\\n")
        migrated = next(
            line.split(":", 2)[2]
            for line in pathlib.Path("/proc/self/cgroup").read_text(
                encoding="ascii").splitlines()
            if line.startswith("0::"))
        with open(migration_file, "w", encoding="ascii") as stream:
            stream.write(migrated)
            stream.flush()
            os.fsync(stream.fileno())
    os.close(1)
    os.close(2)
    with open(pid_file, "w", encoding="ascii") as stream:
        stream.write(str(os.getpid()))
        stream.flush()
        os.fsync(stream.fileno())
    time.sleep(0.5)
    with open(marker, "w", encoding="ascii") as stream:
        stream.write("escaped\\n")
    time.sleep(30)
    os._exit(0)
os.waitpid(child, 0)
deadline = time.monotonic() + 2
while not os.path.exists(pid_file):
    if time.monotonic() >= deadline:
        raise RuntimeError("double-fork descendant did not start")
    time.sleep(0.005)
os.write(1, b"leader-output\\n")
if mode == "nonzero":
    sys.exit(7)
if mode == "timeout":
    time.sleep(30)
""".format(
                    pid_file=str(pid_file),
                    migration_file=str(migration_file),
                    marker=str(marker), mode=mode)
                executable_script(executable, source)
                descendant_pid = None
                try:
                    return_code, unused_stdout, stderr = \
                        run_supervisor_fixture(root, executable, timeout)
                    self.assertEqual(
                        return_code, expected_return,
                        stderr.decode(errors="replace"))
                    descendant_pid = int(
                        pid_file.read_text(encoding="ascii"))
                    if mode == "escaped-cgroup":
                        migrated = migration_file.read_text(
                            encoding="ascii")
                        self.assertNotIn("leopard2-exhaustive-", migrated)
                    time.sleep(0.7)
                    self.assertFalse(marker.exists())
                    self.assertFalse(
                        Path("/proc/%d" % descendant_pid).exists(),
                        "double-fork descendant was not reaped")
                finally:
                    if (descendant_pid is None and pid_file.exists()):
                        descendant_pid = int(
                            pid_file.read_text(encoding="ascii"))
                    if (descendant_pid is not None and
                            Path("/proc/%d" % descendant_pid).exists()):
                        try:
                            os.kill(descendant_pid, 9)
                        except ProcessLookupError:
                            pass

    def test_nested_cgroup_tree_is_removed_on_normal_completion(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-nested-cgroup-normal-") as directory:
            root = Path(directory)
            executable = root / "fixture"
            scope_file = root / "nested.scope"
            executable_script(executable, """
import os
import pathlib

membership = next(
    line.split(":", 2)[2]
    for line in pathlib.Path("/proc/self/cgroup").read_text(
        encoding="ascii").splitlines()
    if line.startswith("0::"))
scope = pathlib.Path("/sys/fs/cgroup") / membership.lstrip("/")
nested = scope / "adversarial-nested"
deeper = nested / "deeper"
nested.mkdir()
deeper.mkdir()
with (deeper / "cgroup.procs").open("w", encoding="ascii") as stream:
    stream.write("0\\n")
with pathlib.Path({scope_file!r}).open(
        "w", encoding="ascii") as stream:
    stream.write(str(deeper))
    stream.flush()
    os.fsync(stream.fileno())
""".format(scope_file=str(scope_file)))
            return_code, unused_stdout, stderr = run_supervisor_fixture(
                root, executable, 5.0)
            self.assertEqual(
                return_code, 0, stderr.decode(errors="replace"))
            nested_scope = Path(scope_file.read_text(encoding="ascii"))
            self.assertFalse(
                nested_scope.exists(),
                "normal shard completion leaked a nested cgroup subtree")

    def test_output_is_bounded_while_child_is_running(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-output-bound-") as directory:
            root = Path(directory)
            executable = root / "fixture"
            executable_script(executable, """
import os
payload = b"x" * ({limit} + 65536)
offset = 0
while offset < len(payload):
    offset += os.write(1, payload[offset:])
""".format(limit=RUNNER.MAX_SHARD_JSON_BYTES))
            return_code, stdout, stderr = run_supervisor_fixture(
                root, executable, 5.0)
            self.assertEqual(return_code, 124)
            self.assertEqual(len(stdout), RUNNER.MAX_SHARD_JSON_BYTES)
            self.assertLessEqual(len(stderr), RUNNER.MAX_SHARD_JSON_BYTES)
            self.assertIn(b"exceeded its byte limit", stderr)
            shard = root / "campaign/shard-0000"
            self.assertFalse(any(".tmp-" in path.name
                                 for path in shard.iterdir()))

    def test_parent_exit_signal_tears_down_shard(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-parent-exit-") as directory:
            root = Path(directory)
            executable = root / "fixture"
            executable_script(
                executable,
                "import os\nimport time\n"
                "os.close(1)\nos.close(2)\ntime.sleep(30)\n")
            return_code, unused_stdout, stderr = run_supervisor_fixture(
                root, executable, 5.0, close_parent_watch=True)
            self.assertEqual(return_code, 124)
            self.assertIn(b"campaign parent exited", stderr)

    def test_parent_sigkill_leaves_no_campaign_cgroup(self) -> None:
        """A real parent death must be cleaned before the harness intervenes."""
        with tempfile.TemporaryDirectory(
                prefix="leo2-parent-sigkill-") as directory:
            root = Path(directory)
            executable = root / "fixture"
            pid_file = root / "fixture.pid"
            ready = root / "ready.json"
            executable_script(executable, """
import os
import time
with open({pid_file!r}, "w", encoding="ascii") as stream:
    stream.write(str(os.getpid()))
    stream.flush()
    os.fsync(stream.fileno())
time.sleep(30)
""".format(pid_file=str(pid_file)))
            driver = root / "driver.py"
            driver.write_text("""
import importlib.util
import pathlib
import sys

test_path = pathlib.Path({test_path!r})
spec = importlib.util.spec_from_file_location(
    "sigkill_exhaustive_test", test_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
module.run_supervisor_fixture(
    pathlib.Path({root!r}), pathlib.Path({executable!r}), 30.0,
    ready_path=pathlib.Path({ready!r}))
""".format(
                test_path=str(Path(__file__).resolve()),
                root=str(root), executable=str(executable),
                ready=str(ready)), encoding="utf-8")
            driver_process = subprocess.Popen(
                [sys.executable, "-I", "-S", "-B", str(driver)],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                env=RUNNER.CHILD_ENV)
            record = None
            fixture_pid = None
            try:
                deadline = time.monotonic() + 10.0
                while time.monotonic() < deadline:
                    if ready.exists() and pid_file.exists():
                        record = json.loads(
                            ready.read_text(encoding="utf-8"))
                        fixture_pid = int(
                            pid_file.read_text(encoding="ascii"))
                        break
                    if driver_process.poll() is not None:
                        break
                    time.sleep(0.01)
                if record is None:
                    driver_process.wait(timeout=1)
                    self.fail("SIGKILL fixture did not start")
                os.kill(driver_process.pid, signal.SIGKILL)
                driver_process.wait(timeout=5.0)
                campaign_scope = Path(record["campaign_scope"])
                supervisor_pid = int(record["supervisor_pid"])
                deadline = time.monotonic() + \
                    RUNNER.SUPERVISOR_GRACE_SECONDS
                while (campaign_scope.exists() or
                       Path("/proc/%d" % supervisor_pid).exists() or
                       Path("/proc/%d" % fixture_pid).exists()):
                    if time.monotonic() >= deadline:
                        break
                    time.sleep(0.02)

                # These assertions intentionally precede all emergency
                # harness cleanup: the last supervisor owns reclamation.
                self.assertFalse(
                    campaign_scope.exists(),
                    "parent SIGKILL leaked the campaign cgroup")
                self.assertFalse(
                    Path("/proc/%d" % supervisor_pid).exists(),
                    "parent SIGKILL leaked the shard supervisor")
                self.assertFalse(
                    Path("/proc/%d" % fixture_pid).exists(),
                    "parent SIGKILL leaked the shard process")
            finally:
                if driver_process.poll() is None:
                    driver_process.kill()
                    driver_process.wait(timeout=5.0)
                if record is not None:
                    for pid in (
                            int(record["supervisor_pid"]),
                            fixture_pid):
                        if pid is None:
                            continue
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    campaign_scope = Path(record["campaign_scope"])
                    if campaign_scope.exists():
                        for child in tuple(campaign_scope.iterdir()):
                            if child.is_dir():
                                RUNNER.signal_cgroup(child)
                                RUNNER.remove_empty_cgroup(child)
                        RUNNER.remove_empty_campaign_cgroup(
                            campaign_scope,
                            RUNNER.cgroup_v2_current())

    def test_parent_sigkill_during_guarded_launch_windows_cleans_scopes(
            self) -> None:
        """Fork-server guardian owns pre-launch and killed-supervisor windows."""
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("guardian launch-window test needs two CPUs")
        expected_two = [
            dict(value)
            for value in RUNNER.reconstruct_expected_shards(2)
        ]
        for launch_mode, label in (
                ("before-control", "before-first-supervisor"),
                ("kill-supervisor", "after-supervisor-kill")):
            with self.subTest(window=label), tempfile.TemporaryDirectory(
                    prefix="leo2-guardian-%s-" % label) as directory:
                root = Path(directory)
                output = root / "campaign"
                output.mkdir(mode=0o700)
                artifacts = output / "artifacts"
                artifacts.mkdir(mode=0o700)
                executable = artifacts / RUNNER.TARGET_NAME
                executable_script(executable, """
import os
import pathlib
import signal
import sys
import time
index = int(sys.argv[sys.argv.index("--shard-index") + 1])
root = pathlib.Path({root!r})
if index != 0:
    time.sleep(30)
    sys.exit(0)
supervisor_pid = os.getppid()
first = os.fork()
if first == 0:
    os.setsid()
    second = os.fork()
    if second != 0:
        os._exit(0)
    membership = next(
        line.split(":", 2)[2]
        for line in pathlib.Path("/proc/self/cgroup").read_text(
            encoding="ascii").splitlines()
        if line.startswith("0::"))
    scope = pathlib.Path("/sys/fs/cgroup") / membership.lstrip("/")
    nested = scope / "adversarial-nested"
    deeper = nested / "deeper"
    nested.mkdir()
    deeper.mkdir()
    with (scope.parent.parent / "cgroup.procs").open(
            "w", encoding="ascii") as stream:
        stream.write("0\\n")
    for descriptor in (1, 2):
        os.close(descriptor)
    with (root / "migrated.pid").open("w", encoding="ascii") as stream:
        stream.write(str(os.getpid()))
        stream.flush()
        os.fsync(stream.fileno())
    with (root / "supervisor.pid").open(
            "w", encoding="ascii") as stream:
        stream.write(str(supervisor_pid))
        stream.flush()
        os.fsync(stream.fileno())
    with (root / "migrated.ready").open(
            "w", encoding="ascii") as stream:
        stream.write("ready\\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.kill(supervisor_pid, signal.SIGKILL)
    time.sleep(1)
    (root / "escaped.marker").write_text(
        "escaped\\n", encoding="ascii")
    time.sleep(30)
    os._exit(0)
os.waitpid(first, 0)
time.sleep(30)
""".format(root=str(root)))
                expected = expected_two
                checkpoint = root / "checkpoint.json"
                driver = root / "driver.py"
                driver.write_text("""
import importlib.util
import json
import os
import pathlib
import sys
import time

runner_path = pathlib.Path({runner_path!r})
spec = importlib.util.spec_from_file_location(
    "guarded_launch_runner", runner_path)
runner = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = runner
spec.loader.exec_module(runner)
root = pathlib.Path({root!r})
output = root / "campaign"
executable = output / "artifacts" / runner.TARGET_NAME
expected = {expected!r}
request = {{
    "workers": 2,
    "shard_count": 2,
    "frozen_binary": runner.file_identity(executable),
    "runner": runner.file_identity(runner_path),
    "taskset": runner.file_identity(runner.CANONICAL_TASKSET),
    "allowed_cpus": sorted(os.sched_getaffinity(0))[:2],
    "timeout_seconds_per_shard": 30.0,
    "build_provenance": {{"mode": 1}},
}}
real_start = runner.start_campaign_guardian
def controlled_start(descriptor, campaign_scope):
    with pathlib.Path({checkpoint!r}).open(
            "w", encoding="utf-8") as stream:
        json.dump({{
            "campaign_scope": str(campaign_scope),
            "mode": {launch_mode!r},
        }}, stream)
        stream.flush()
        os.fsync(stream.fileno())
    if {launch_mode!r} == "before-control":
        time.sleep(30)
    real_start(descriptor, campaign_scope)
runner.start_campaign_guardian = controlled_start
try:
    runner.launch_shards(output, request, expected)
except BaseException:
    # Keep the coordinator alive until the harness explicitly SIGKILLs it.
    if {launch_mode!r} == "kill-supervisor":
        time.sleep(30)
    raise
""".format(
                    runner_path=str(RUNNER_PATH.resolve()),
                    root=str(root), launch_mode=launch_mode,
                    checkpoint=str(checkpoint),
                    expected=expected), encoding="utf-8")
                driver_process = subprocess.Popen(
                    [sys.executable, "-I", "-S", "-B", str(driver)],
                    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, env=RUNNER.CHILD_ENV)
                record = None
                migrated_pid = None
                supervisor_pid = None
                try:
                    deadline = time.monotonic() + 10.0
                    while time.monotonic() < deadline:
                        if checkpoint.exists():
                            record = json.loads(
                                checkpoint.read_text(encoding="utf-8"))
                            break
                        if driver_process.poll() is not None:
                            break
                        time.sleep(0.01)
                    if record is None:
                        driver_process.wait(timeout=1.0)
                        self.fail(
                            "guardian launch-window fixture did not start")
                    if launch_mode == "kill-supervisor":
                        migrated_ready = root / "migrated.ready"
                        fixture_deadline = time.monotonic() + 5.0
                        while (not migrated_ready.exists() and
                               time.monotonic() < fixture_deadline):
                            time.sleep(0.01)
                        self.assertTrue(
                            migrated_ready.exists(),
                            "partial-launch migrated shard never started")
                        migrated_pid = int(
                            (root / "migrated.pid").read_text(
                                encoding="ascii"))
                        supervisor_pid = int(
                            (root / "supervisor.pid").read_text(
                                encoding="ascii"))
                        while (Path("/proc/%d" % supervisor_pid).exists() and
                               time.monotonic() < fixture_deadline):
                            time.sleep(0.01)
                        self.assertFalse(
                            Path("/proc/%d" % supervisor_pid).exists(),
                            "migrated shard did not kill its supervisor "
                            "before coordinator death")
                    os.kill(driver_process.pid, signal.SIGKILL)
                    driver_process.wait(timeout=5.0)
                    campaign_scope = Path(record["campaign_scope"])
                    deadline = time.monotonic() + \
                        RUNNER.SUPERVISOR_GRACE_SECONDS
                    while (campaign_scope.exists() or
                           (migrated_pid is not None and
                            Path("/proc/%d" % migrated_pid).exists()) or
                           (supervisor_pid is not None and
                            Path("/proc/%d" % supervisor_pid).exists())):
                        if time.monotonic() >= deadline:
                            break
                        time.sleep(0.02)

                    # No test cleanup has run before this ownership assertion.
                    self.assertFalse(
                        campaign_scope.exists(),
                        "%s leaked the guarded campaign cgroup" % label)
                    if migrated_pid is not None:
                        self.assertFalse(
                            Path("/proc/%d" % supervisor_pid).exists(),
                            "guardian retained the killed supervisor")
                        self.assertFalse(
                            Path("/proc/%d" % migrated_pid).exists(),
                            "guardian did not reap its adopted "
                            "upward-migrated child")
                        self.assertFalse(
                            (root / "escaped.marker").exists(),
                            "upward-migrated child escaped guardian cleanup")
                finally:
                    if driver_process.poll() is None:
                        driver_process.kill()
                        driver_process.wait(timeout=5.0)
                    if record is not None:
                        campaign_scope = Path(record["campaign_scope"])
                        if campaign_scope.exists():
                            for child in tuple(campaign_scope.iterdir()):
                                if child.is_dir():
                                    RUNNER.signal_cgroup(child)
                                    RUNNER.remove_empty_cgroup(child)
                            RUNNER.remove_empty_campaign_cgroup(
                                campaign_scope,
                                RUNNER.cgroup_v2_current())
                    if (migrated_pid is not None and
                            Path("/proc/%d" % migrated_pid).exists()):
                        try:
                            os.kill(migrated_pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    def test_parent_sigkill_cleans_late_sibling_after_inventory(self) -> None:
        """A resumed tracked PID cannot strand a new campaign child."""
        with tempfile.TemporaryDirectory(
                prefix="leo2-guardian-late-sibling-") as directory:
            root = Path(directory)
            output = root / "campaign"
            output.mkdir(mode=0o700)
            artifacts = output / "artifacts"
            artifacts.mkdir(mode=0o700)
            executable = artifacts / RUNNER.TARGET_NAME
            executable_script(
                executable, "import time\ntime.sleep(30)\n")
            checkpoint = root / "checkpoint.json"
            driver = root / "driver.py"
            expected = [
                dict(value)
                for value in RUNNER.reconstruct_expected_shards(1)
            ]
            driver.write_text("""
import importlib.util
import json
import os
import pathlib
import sys
import time

runner_path = pathlib.Path({runner_path!r})
spec = importlib.util.spec_from_file_location(
    "late_sibling_runner", runner_path)
runner = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = runner
spec.loader.exec_module(runner)
root = pathlib.Path({root!r})
output = root / "campaign"
executable = output / "artifacts" / runner.TARGET_NAME
request = {{
    "workers": 1,
    "shard_count": 1,
    "frozen_binary": runner.file_identity(executable),
    "runner": runner.file_identity(runner_path),
    "taskset": runner.file_identity(runner.CANONICAL_TASKSET),
    "allowed_cpus": [min(os.sched_getaffinity(0))],
    "timeout_seconds_per_shard": 30.0,
    "build_provenance": {{"mode": 1}},
}}
real_start = runner.start_campaign_guardian
def controlled_start(descriptor, campaign_scope):
    children = runner.adopted_child_pids()
    if len(children) != 1:
        raise RuntimeError("guardian PID inventory changed")
    with pathlib.Path({checkpoint!r}).open(
            "w", encoding="utf-8") as stream:
        json.dump({{
            "campaign_scope": str(campaign_scope),
            "guardian_pid": children[0],
        }}, stream)
        stream.flush()
        os.fsync(stream.fileno())
    # The harness installs a stopped external PID before killing this
    # coordinator, so launch authorization must not be sent first.
    time.sleep(30)
    real_start(descriptor, campaign_scope)
runner.start_campaign_guardian = controlled_start
runner.launch_shards(output, request, {expected!r})
""".format(
                runner_path=str(RUNNER_PATH.resolve()),
                root=str(root), checkpoint=str(checkpoint),
                expected=expected), encoding="utf-8")
            driver_process = subprocess.Popen(
                [sys.executable, "-I", "-S", "-B", str(driver)],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, env=RUNNER.CHILD_ENV)
            adversary_process = None
            record = None
            adversary_pid = None
            try:
                deadline = time.monotonic() + 10.0
                while time.monotonic() < deadline:
                    if checkpoint.exists():
                        record = json.loads(
                            checkpoint.read_text(encoding="utf-8"))
                        break
                    if driver_process.poll() is not None:
                        break
                    time.sleep(0.01)
                if record is None:
                    driver_process.wait(timeout=1.0)
                    self.fail("late-sibling guardian fixture did not start")
                campaign_scope = Path(record["campaign_scope"])
                ready = root / "adversary.ready"
                pid_file = root / "adversary.pid"
                marker = root / "late-child.created"
                adversary = root / "adversary.py"
                adversary.write_text("""
import os
import pathlib
import signal

campaign_scope = pathlib.Path({campaign_scope!r})
ready = pathlib.Path({ready!r})
pid_file = pathlib.Path({pid_file!r})
marker = pathlib.Path({marker!r})
worker = os.fork()
if worker != 0:
    os.waitpid(worker, 0)
    raise SystemExit(0)
scope = campaign_scope / "shard-0000"
with (scope / "cgroup.procs").open(
        "w", encoding="ascii") as stream:
    stream.write("0\\n")
with pid_file.open("w", encoding="ascii") as stream:
    stream.write(str(os.getpid()))
    stream.flush()
    os.fsync(stream.fileno())
with ready.open("w", encoding="ascii") as stream:
    stream.write("stopping\\n")
    stream.flush()
    os.fsync(stream.fileno())
os.kill(os.getpid(), signal.SIGSTOP)
late = campaign_scope / "late-child"
late.mkdir()
with (late / "cgroup.procs").open(
        "w", encoding="ascii") as stream:
    stream.write("0\\n")
with marker.open("w", encoding="ascii") as stream:
    stream.write("created-after-sigcont\\n")
    stream.flush()
    os.fsync(stream.fileno())
os._exit(0)
""".format(
                    campaign_scope=str(campaign_scope),
                    ready=str(ready), pid_file=str(pid_file),
                    marker=str(marker)), encoding="utf-8")
                adversary_process = subprocess.Popen(
                    [sys.executable, "-I", "-S", "-B", str(adversary)],
                    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, env=RUNNER.CHILD_ENV)
                deadline = time.monotonic() + 5.0
                while (not ready.exists() and
                       time.monotonic() < deadline):
                    time.sleep(0.01)
                self.assertTrue(
                    ready.exists(), "late-sibling adversary did not start")
                adversary_pid = int(
                    pid_file.read_text(encoding="ascii"))
                state = None
                while time.monotonic() < deadline:
                    status_path = Path(
                        "/proc/%d/status" % adversary_pid)
                    if status_path.exists():
                        state = next(
                            line for line in
                            status_path.read_text(
                                encoding="ascii").splitlines()
                            if line.startswith("State:")).split()[1]
                        if state in ("T", "t"):
                            break
                    time.sleep(0.01)
                self.assertIn(state, ("T", "t"),
                              "late-sibling adversary did not stop")

                os.kill(driver_process.pid, signal.SIGKILL)
                driver_process.wait(timeout=5.0)
                adversary_process.wait(timeout=10.0)
                guardian_pid = int(record["guardian_pid"])
                deadline = time.monotonic() + \
                    RUNNER.SUPERVISOR_GRACE_SECONDS
                while (campaign_scope.exists() or
                       Path("/proc/%d" % guardian_pid).exists()):
                    if time.monotonic() >= deadline:
                        break
                    time.sleep(0.02)

                # All assertions precede emergency harness cleanup.  The
                # marker proves creation happened only after guardian SIGCONT.
                self.assertTrue(
                    marker.exists(),
                    "tracked PID never created its late campaign child")
                self.assertFalse(
                    campaign_scope.exists(),
                    "late campaign child escaped guardian cleanup")
                self.assertFalse(
                    Path("/proc/%d" % guardian_pid).exists(),
                    "late-sibling guardian did not exit")
            finally:
                if driver_process.poll() is None:
                    driver_process.kill()
                    driver_process.wait(timeout=5.0)
                if adversary_pid is None and (root / "adversary.pid").exists():
                    adversary_pid = int(
                        (root / "adversary.pid").read_text(
                            encoding="ascii"))
                if (adversary_pid is not None and
                        Path("/proc/%d" % adversary_pid).exists()):
                    try:
                        os.kill(adversary_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                if adversary_process is not None:
                    try:
                        adversary_process.wait(timeout=5.0)
                    except subprocess.TimeoutExpired:
                        adversary_process.kill()
                        adversary_process.wait(timeout=5.0)
                if record is not None:
                    campaign_scope = Path(record["campaign_scope"])
                    if campaign_scope.exists():
                        RUNNER.signal_cgroup(campaign_scope)
                        RUNNER.remove_empty_cgroup(campaign_scope)

    def test_parallel_launch_uses_contained_sealed_shards(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-parallel-launch-") as directory:
            root = Path(directory)
            output = root / "campaign"
            output.mkdir(mode=0o700)
            artifacts = output / "artifacts"
            artifacts.mkdir(mode=0o700)
            workers = min(2, len(os.sched_getaffinity(0)))
            expected = [
                dict(value)
                for value in RUNNER.reconstruct_expected_shards(workers)
            ]
            executable = artifacts / RUNNER.TARGET_NAME
            source = """
import json
import sys

expected = {expected!r}
index = int(sys.argv[sys.argv.index("--shard-index") + 1])
result = {{
    "schema": {schema!r},
    "mode": 1,
    **expected[index],
}}
print(json.dumps(result, sort_keys=True))
""".format(expected=expected, schema=RUNNER.SHARD_SCHEMA)
            executable_script(executable, source)
            request = {
                "workers": workers,
                "shard_count": workers,
                "frozen_binary": RUNNER.file_identity(executable),
                "runner": RUNNER.file_identity(RUNNER_PATH),
                "taskset": RUNNER.file_identity(
                    RUNNER.CANONICAL_TASKSET),
                "allowed_cpus":
                    sorted(os.sched_getaffinity(0))[:workers],
                "timeout_seconds_per_shard": 5.0,
                "build_provenance": {"mode": 1},
            }
            results = RUNNER.launch_shards(output, request, expected)
            self.assertEqual(results, [
                {"schema": RUNNER.SHARD_SCHEMA, "mode": 1, **value}
                for value in expected
            ])
            for index in range(workers):
                shard = output / ("shard-%04d" % index)
                self.assertEqual(
                    {path.name for path in shard.iterdir()},
                    {"stdout.json", "stderr.txt", "result.json"})

    def test_other_shard_failure_reaps_upward_migrated_descendant(
            self) -> None:
        """Forced supervisor death falls back to the campaign subreaper."""
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("two-shard containment regression needs two CPUs")
        with tempfile.TemporaryDirectory(
                prefix="leo2-two-shard-migration-") as directory:
            root = Path(directory)
            output = root / "campaign"
            output.mkdir(mode=0o700)
            artifacts = output / "artifacts"
            artifacts.mkdir(mode=0o700)
            ready = root / "migrated.ready"
            pid_file = root / "migrated.pid"
            marker = root / "escaped.marker"
            executable = artifacts / RUNNER.TARGET_NAME
            source = """
import os
import pathlib
import signal
import sys
import time

index = int(sys.argv[sys.argv.index("--shard-index") + 1])
ready = pathlib.Path({ready!r})
pid_file = pathlib.Path({pid_file!r})
marker = pathlib.Path({marker!r})
if index == 0:
    deadline = time.monotonic() + 5
    while not ready.exists():
        if time.monotonic() >= deadline:
            raise RuntimeError("migrated shard did not become ready")
        time.sleep(0.005)
    sys.exit(7)

supervisor_pid = os.getppid()
first = os.fork()
if first == 0:
    os.setsid()
    second = os.fork()
    if second != 0:
        os._exit(0)
    membership = next(
        line.split(":", 2)[2]
        for line in pathlib.Path("/proc/self/cgroup").read_text(
            encoding="ascii").splitlines()
        if line.startswith("0::"))
    scope = pathlib.Path("/sys/fs/cgroup") / membership.lstrip("/")
    with (scope.parent.parent / "cgroup.procs").open(
            "w", encoding="ascii") as stream:
        stream.write("0\\n")
    for descriptor in (1, 2):
        os.close(descriptor)
    with pid_file.open("w", encoding="ascii") as stream:
        stream.write(str(os.getpid()))
        stream.flush()
        os.fsync(stream.fileno())
    os.kill(supervisor_pid, signal.SIGSTOP)
    with ready.open("w", encoding="ascii") as stream:
        stream.write("ready\\n")
        stream.flush()
        os.fsync(stream.fileno())
    time.sleep(0.5)
    marker.write_text("escaped\\n", encoding="ascii")
    time.sleep(30)
    os._exit(0)
os.waitpid(first, 0)
time.sleep(30)
""".format(
                ready=str(ready), pid_file=str(pid_file),
                marker=str(marker))
            executable_script(executable, source)
            expected = [
                dict(value)
                for value in RUNNER.reconstruct_expected_shards(2)
            ]
            request = {
                "workers": 2,
                "shard_count": 2,
                "frozen_binary": RUNNER.file_identity(executable),
                "runner": RUNNER.file_identity(RUNNER_PATH),
                "taskset": RUNNER.file_identity(
                    RUNNER.CANONICAL_TASKSET),
                "allowed_cpus":
                    sorted(os.sched_getaffinity(0))[:2],
                "timeout_seconds_per_shard": 10.0,
                "build_provenance": {"mode": 1},
            }
            migrated_pid = None
            cgroup_parent = RUNNER.cgroup_v2_current()
            before_subreaper = RUNNER.child_subreaper_state()
            prefix = "leopard2-exhaustive-%d-" % os.getpid()
            before = {
                path.name for path in cgroup_parent.iterdir()
                if path.is_dir() and path.name.startswith(prefix)
            }
            try:
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError, "shard 0 failed"):
                    RUNNER.launch_shards(output, request, expected)
                migrated_pid = int(pid_file.read_text(encoding="ascii"))
                time.sleep(0.7)
                self.assertFalse(marker.exists())
                self.assertFalse(
                    Path("/proc/%d" % migrated_pid).exists(),
                    "campaign fallback did not reap migrated descendant")
                after = {
                    path.name for path in cgroup_parent.iterdir()
                    if path.is_dir() and path.name.startswith(prefix)
                }
                self.assertEqual(after, before)
                self.assertEqual(
                    RUNNER.child_subreaper_state(), before_subreaper)
            finally:
                if migrated_pid is None and pid_file.exists():
                    migrated_pid = int(
                        pid_file.read_text(encoding="ascii"))
                if (migrated_pid is not None and
                        Path("/proc/%d" % migrated_pid).exists()):
                    try:
                        os.kill(migrated_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

    def test_missing_corrupt_extra_and_swapped_artifacts_fail(self) -> None:
        campaign = self.manifest.parent
        result = campaign / "shard-0000/result.json"
        saved = result.read_bytes()
        result.unlink()
        try:
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.verify_manifest(self.manifest, True)
        finally:
            result.write_bytes(saved)

        frozen = campaign / "artifacts/provenance/file-0000"
        saved = frozen.read_bytes()
        os.chmod(frozen, 0o644)
        frozen.write_bytes(saved + b"x")
        try:
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.verify_manifest(self.manifest, True)
        finally:
            frozen.write_bytes(saved)
            os.chmod(frozen, 0o444)

        extra = campaign / "unexpected"
        extra.write_text("extra\n")
        try:
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.verify_manifest(self.manifest, True)
        finally:
            extra.unlink()

        first = campaign / "shard-0000/stdout.json"
        second = campaign / "shard-0001/stdout.json"
        left, right = first.read_bytes(), second.read_bytes()
        first.write_bytes(right)
        second.write_bytes(left)
        try:
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.verify_manifest(self.manifest, True)
        finally:
            first.write_bytes(left)
            second.write_bytes(right)

    def test_protected_lock_rejects_symlink_hardlink_and_replacement(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-exhaustive-lock-") as directory:
            root = Path(directory)
            lock_path = root / "lock"
            target = root / "target"
            target.write_text("")
            lock_path.symlink_to(target)
            with mock.patch.object(RUNNER, "DEFAULT_LOCK", lock_path):
                with self.assertRaises(OSError):
                    RUNNER.acquire_lock(lock_path, 0)
            lock_path.unlink()
            with mock.patch.object(RUNNER, "DEFAULT_LOCK", lock_path):
                stream = RUNNER.acquire_lock(lock_path, 0)
                try:
                    linked = root / "hardlink"
                    os.link(lock_path, linked)
                    with self.assertRaisesRegex(
                            RUNNER.EvidenceError, "unsafe metadata"):
                        RUNNER.lock_identity(stream, lock_path)
                    linked.unlink()
                    moved = root / "old"
                    lock_path.rename(moved)
                    lock_path.write_text("")
                    os.chmod(lock_path, 0o600)
                    with self.assertRaisesRegex(
                            RUNNER.EvidenceError, "replaced"):
                        RUNNER.lock_identity(stream, lock_path)
                finally:
                    stream.close()

    def test_close_only_lock_stays_held_by_surviving_guardian_copy(
            self) -> None:
        """Coordinator cleanup must not unlock a surviving shared descriptor."""
        with tempfile.TemporaryDirectory(
                prefix="leo2-exhaustive-lock-survivor-") as directory:
            lock_path = Path(directory) / "lock"
            survivor_descriptor = -1
            probe_descriptor = -1
            with mock.patch.object(RUNNER, "DEFAULT_LOCK", lock_path):
                stream = RUNNER.acquire_lock(lock_path, 0)
                try:
                    # os.dup models the descriptor inherited by the lifecycle
                    # guardian: it refers to the same open-file description.
                    survivor_descriptor = os.dup(stream.fileno())
                    probe_descriptor = os.open(
                        lock_path, os.O_RDWR | os.O_CLOEXEC)
                    RUNNER.close_campaign_lock(stream)
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(
                            probe_descriptor,
                            fcntl.LOCK_EX | fcntl.LOCK_NB)

                    os.close(survivor_descriptor)
                    survivor_descriptor = -1
                    fcntl.flock(
                        probe_descriptor,
                        fcntl.LOCK_EX | fcntl.LOCK_NB)
                finally:
                    if not stream.closed:
                        stream.close()
                    if survivor_descriptor >= 0:
                        os.close(survivor_descriptor)
                    if probe_descriptor >= 0:
                        os.close(probe_descriptor)

    def test_worker_count_rejected_before_exact_reconstruction(self) -> None:
        with mock.patch.object(
                RUNNER.os, "sched_getaffinity",
                return_value={2, 4, 6}):
            self.assertEqual(RUNNER.selected_worker_cpus(3), [2, 4, 6])
            for workers in (0, -1, 4, 10**9):
                with self.subTest(workers=workers), self.assertRaisesRegex(
                        RUNNER.EvidenceError, "at most 128"):
                    RUNNER.selected_worker_cpus(workers)
        for workers in (0, -1, 129, 10**9):
            with self.subTest(
                    reconstruction_workers=workers), self.assertRaisesRegex(
                        RUNNER.EvidenceError, "1..128"):
                RUNNER.reconstruct_expected_shards(workers)
        with mock.patch.object(
                RUNNER.os, "sched_getaffinity",
                return_value=set(range(256))):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "at most 128"):
                RUNNER.selected_worker_cpus(129)

    def test_nonfinite_and_invalid_timeouts_fail_before_deadlines(self) -> None:
        with mock.patch.object(
                RUNNER.os, "sched_getaffinity", return_value={0}):
            RUNNER.validate_run_limits(1, 1.0, 0.0)
            for timeout in (0.0, -1.0, float("nan"), float("inf")):
                with self.subTest(timeout=timeout), self.assertRaisesRegex(
                        RUNNER.EvidenceError, "timeouts are invalid"):
                    RUNNER.validate_run_limits(1, timeout, 0.0)
            for lock_timeout in (
                    -1.0, float("nan"), float("inf")):
                with self.subTest(
                        lock_timeout=lock_timeout), self.assertRaisesRegex(
                            RUNNER.EvidenceError, "timeouts are invalid"):
                    RUNNER.validate_run_limits(1, 1.0, lock_timeout)

    def test_build_record_rejects_wrong_mode_provenance(self) -> None:
        manifest = RUNNER.load_json(self.manifest)
        build = dict(manifest["request"]["build_provenance"])
        closure = dict(build["production_closure"])
        records = [dict(value)
                   for value in closure["source_object_compile_closure"]]
        records[0] = dict(records[0])
        records[0]["compile_entry"] = {
            "arguments": ["c++", "-O3",
                          RUNNER.MODE_DEFINITIONS[2]],
        }
        closure["source_object_compile_closure"] = records
        build["production_closure"] = closure
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "diagnostic mode definition"):
            RUNNER.validate_build_record(build)

        build = dict(manifest["request"]["build_provenance"])
        closure = dict(build["production_closure"])
        records = [dict(value)
                   for value in closure["source_object_compile_closure"]]
        records[1] = dict(records[1])
        records[1]["compile_entry"] = {
            "arguments": ["c++", "-O3",
                          RUNNER.MODE_DEFINITIONS[1]],
        }
        closure["source_object_compile_closure"] = records
        build["production_closure"] = closure
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "diagnostic mode definition"):
            RUNNER.validate_build_record(build)

        build = dict(manifest["request"]["build_provenance"])
        targets = [dict(value)
                   for value in build["target_compile_closure"]]
        oracle = next(index for index, record in enumerate(targets)
                      if Path(record["source"]["path"]).name ==
                      "direct_oracle.cpp")
        targets[oracle] = dict(targets[oracle])
        targets[oracle]["compile_entry"] = dict(
            targets[oracle]["compile_entry"])
        targets[oracle]["compile_entry"]["arguments"] = [
            "c++", "-O3", RUNNER.MODE_DEFINITIONS[1],
        ]
        build["target_compile_closure"] = targets
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "diagnostic mode definition"):
            RUNNER.validate_build_record(build)

        build = dict(manifest["request"]["build_provenance"])
        build["mode_independent_object_sha256"] = "0" * 64
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "object identity changed"):
            RUNNER.validate_build_record(build)

        for field, replacement in (
                ("head", "not-a-commit"),
                ("head_tree", "f" * 39)):
            build = dict(manifest["request"]["build_provenance"])
            source = dict(build["source"])
            source[field] = replacement
            source["digest"] = RUNNER.record_digest(source)
            build["source"] = source
            build["digest"] = RUNNER.record_digest(build)
            with self.subTest(field=field), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "source closure"):
                RUNNER.validate_build_record(build)

        build = dict(manifest["request"]["build_provenance"])
        source = dict(build["source"])
        files = dict(source["files"])
        runner_key = (
            "experiments/leopard2/direct_repair/"
            "run_small_direct_exhaustive.py")
        files[runner_key] = dict(files[runner_key])
        files[runner_key]["path"] = "/tmp/escaped-runner.py"
        source["files"] = files
        source["digest"] = RUNNER.record_digest(source)
        build["source"] = source
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "bound source identity"):
            RUNNER.validate_build_record(build)

    def test_build_record_rejects_schema_extensions_and_bad_identities(
            self) -> None:
        manifest = RUNNER.load_json(self.manifest)
        original = manifest["request"]["build_provenance"]

        build = copy.deepcopy(original)
        build["unbound_extension"] = "forged"
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "build provenance"):
            RUNNER.validate_build_record(build)

        build = copy.deepcopy(original)
        build["target_compile_closure"][0]["unbound_extension"] = "forged"
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "compile record"):
            RUNNER.validate_build_record(build)

        for mutation in ("missing", "extra"):
            build = copy.deepcopy(original)
            identity = build["target_compile_closure"][0]["object"]
            if mutation == "missing":
                identity.pop("ctime_ns")
            else:
                identity["unbound_extension"] = 1
            build["digest"] = RUNNER.record_digest(build)
            with self.subTest(mutation=mutation), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "object identity"):
                RUNNER.validate_build_record(build)

    def test_duplicate_json_object_keys_are_rejected(self) -> None:
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "duplicate key"):
            RUNNER.decode_json_bytes(
                b'{"schema":"first","schema":"second"}', "test fixture")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "invalid JSON"):
            RUNNER.decode_json_bytes(
                b"[" * 10000 + b"0" + b"]" * 10000,
                "excessively nested fixture")

    def test_compile_source_and_all_mode_spellings_are_fail_closed(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-compile-source-") as directory:
            root = Path(directory)
            declared = root / "declared.cpp"
            other = root / "other.cpp"
            declared.write_text("declared\n")
            other.write_text("other\n")
            entry = {"directory": str(root), "file": str(declared)}
            RUNNER.require_compile_source_match(
                entry, ["c++", "-c", str(declared)], declared)
            for arguments in (
                    ["c++", "-c", str(other)],
                    ["c++", "-c", str(declared), "-c", str(other)],
                    ["c++", str(declared)]):
                with self.subTest(arguments=arguments), \
                        self.assertRaises(RUNNER.EvidenceError):
                    RUNNER.require_compile_source_match(
                        entry, arguments, declared)

        name = RUNNER.MODE_CACHE_VARIABLE
        expected = RUNNER.MODE_DEFINITIONS[1]
        RUNNER.require_mode_macro(["c++", expected], 1, "fixture")
        for arguments in (
                ["c++", "-D%s=0" % name],
                ["c++", "-D", "%s=1" % name],
                ["c++", "-U%s" % name],
                ["c++", "-U", name],
                ["c++", "/D%s=1" % name],
                ["c++", expected, RUNNER.MODE_DEFINITIONS[2]]):
            with self.subTest(arguments=arguments), \
                    self.assertRaisesRegex(
                        RUNNER.EvidenceError, "alternate"):
                RUNNER.require_mode_macro(arguments, 1, "fixture")
        with self.assertRaisesRegex(RUNNER.EvidenceError, "alternate"):
            RUNNER.require_mode_macro(["c++", expected], None, "fixture")
        for experiment in RUNNER.FORBIDDEN_EXPERIMENT_MACROS:
            for arguments in (
                    ["c++", "-D%s=1" % experiment],
                    ["c++", "-D", "%s=1" % experiment],
                    ["c++", "-U%s" % experiment]):
                with self.subTest(
                        experiment=experiment, arguments=arguments), \
                        self.assertRaisesRegex(
                            RUNNER.EvidenceError,
                            "forbidden experiment macro"):
                    RUNNER.require_forbidden_experiment_macros_absent(
                        arguments, "fixture")

    def test_bool_mode_controlled_experiments_and_proof_tamper_fail(
            self) -> None:
        original = RUNNER.load_json(
            self.manifest)["request"]["build_provenance"]

        build = copy.deepcopy(original)
        build["mode"] = True
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "build provenance"):
            RUNNER.validate_build_record(build)

        for name in ("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
                     "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE"):
            build = copy.deepcopy(original)
            build["controlled_configuration"][name] = "ON"
            build["digest"] = RUNNER.record_digest(build)
            with self.subTest(name=name), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "controlled build configuration"):
                RUNNER.validate_build_record(build)

        build = copy.deepcopy(original)
        proof = build["reproducible_rebuild"]
        proof["rebuilt_material"]["archive"]["sha256"] = "0" * 64
        proof["digest"] = RUNNER.record_digest(proof)
        build["digest"] = RUNNER.record_digest(build)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "does not reproduce"):
            RUNNER.validate_build_record(build)

    def test_output_root_and_subdirectory_creation_are_fail_closed(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-output-root-") as directory:
            parent = Path(directory)
            source = parent / "source"
            source.mkdir(mode=0o700)
            output = RUNNER.prepare_output_root(
                parent / "campaign", source, 2)
            self.assertEqual(
                stat.S_IMODE(output.lstat().st_mode), 0o700)
            (output / "unknown").write_text("unknown\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "unknown top-level"):
                RUNNER.prepare_output_root(output, source, 2)
            (output / "unknown").unlink()

            outside = parent / "outside"
            outside.mkdir(mode=0o700)
            linked = output / "linked"
            linked.symlink_to(outside, target_is_directory=True)
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.ensure_directory(linked / "child", output)
            self.assertFalse((outside / "child").exists())
            linked.unlink()

            unsafe = parent / "unsafe"
            unsafe.mkdir(mode=0o755)
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.prepare_output_root(unsafe, source, 1)

    def test_bounded_canonical_files_and_finite_json_only(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-bounded-json-") as directory:
            root = Path(directory)
            oversized = root / "oversized.json"
            with oversized.open("wb") as stream:
                stream.truncate(RUNNER.MAX_SHARD_JSON_BYTES + 1)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "byte bound"):
                RUNNER.load_json(
                    oversized, RUNNER.MAX_SHARD_JSON_BYTES)

            target = root / "target"
            target.write_bytes(b"target")
            linked = root / "linked"
            linked.symlink_to(target)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "canonical regular file"):
                RUNNER.file_identity(linked)

            for constant in (
                    b"NaN", b"Infinity", b"-Infinity", b"1e999"):
                with self.subTest(constant=constant), \
                        self.assertRaisesRegex(
                            RUNNER.EvidenceError, "non-finite"):
                    RUNNER.decode_json_bytes(
                        b'{"value":' + constant + b"}", "fixture")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "invalid JSON"):
                RUNNER.decode_json_bytes(
                    b'{"value":' + b"9" * 5000 + b"}",
                    "oversized integer fixture")
            with self.assertRaises(ValueError):
                RUNNER.canonical_bytes({"value": float("nan")})
            with self.assertRaises(ValueError):
                RUNNER.atomic_json(
                    root / "nan.json", {"value": float("nan")})
            self.assertFalse((root / "nan.json").exists())

    def test_atomic_replace_failure_is_restart_safe(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-atomic-restart-") as directory:
            root = Path(directory)
            destination = root / "request.json"
            with mock.patch.object(
                    RUNNER.os, "replace",
                    side_effect=OSError("injected replace failure")):
                with self.assertRaisesRegex(
                        OSError, "injected replace failure"):
                    RUNNER.atomic_json(destination, {"generation": 1})
            self.assertFalse(destination.exists())
            self.assertFalse(tuple(root.glob("request.json.tmp-*")))
            RUNNER.atomic_json(destination, {"generation": 2})
            self.assertEqual(
                RUNNER.load_json(destination), {"generation": 2})

            source = root / "source"
            source.mkdir(mode=0o700)
            campaign = root / "campaign"
            campaign.mkdir(mode=0o700)
            stale = campaign / (
                "request.json.tmp-%d" % (os.getpid() + 1000))
            stale.write_bytes(b'{"interrupted":')
            stale.chmod(0o600)
            RUNNER.prepare_output_root(campaign, source, 1)
            self.assertFalse(stale.exists())

    def test_run_checked_reaps_detached_descendants_after_timeout(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-run-checked-descendant-") as directory:
            pid_path = Path(directory) / "grandchild.pid"
            program = (
                "import os,time\n"
                "pid=os.fork()\n"
                "if pid==0:\n"
                " os.setsid()\n"
                " open(%r,'w').write(str(os.getpid()))\n"
                " time.sleep(30)\n"
                " os._exit(0)\n"
                "time.sleep(30)\n" % str(pid_path))
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "timed out|failed to execute"):
                RUNNER.run_checked(
                    ["/usr/bin/python3", "-I", "-S", "-B", "-c", program],
                    "detached descendant fixture", timeout=0.5)
            self.assertTrue(pid_path.is_file())
            grandchild = int(pid_path.read_text(encoding="ascii"))
            with self.assertRaises(ProcessLookupError):
                os.kill(grandchild, 0)

    def test_file_snapshot_regular_to_fifo_swap_never_blocks(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-fifo-swap-") as directory:
            victim = Path(directory) / "victim"
            victim.write_bytes(b"regular")
            real_open = os.open
            observed_flags = []

            def swapping_open(path, flags, *arguments):
                if Path(path) == victim and not observed_flags:
                    observed_flags.append(flags)
                    victim.unlink()
                    os.mkfifo(victim, 0o600)
                return real_open(path, flags, *arguments)

            with mock.patch.object(RUNNER.os, "open", swapping_open):
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError, "not a regular file"):
                    RUNNER.file_snapshot(victim)
            self.assertTrue(
                observed_flags[0] & getattr(os, "O_NONBLOCK", 0),
                "file snapshot open omitted O_NONBLOCK")


if __name__ == "__main__":
    unittest.main()
