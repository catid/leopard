#!/usr/bin/env python3
"""Synthetic no-build tests for the exact-main per-role producer stage."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
import leopard2_exact_main_baseline as identity_contract  # noqa: E402
import leopard2_exact_main_baseline_acquire as acquire  # noqa: E402
import leopard2_exact_main_baseline_record as record_contract  # noqa: E402
import leopard2_exact_main_baseline_verifier as verifier  # noqa: E402


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


identity_fixtures = load_module(
    "leopard2_exact_main_build_stage_identity_fixtures",
    TOOLS / "leopard2_exact_main_baseline_test.py",
)
record_fixtures = load_module(
    "leopard2_exact_main_build_stage_record_fixtures",
    TOOLS / "leopard2_exact_main_baseline_record_test.py",
)
verifier_fixtures = load_module(
    "leopard2_exact_main_build_stage_verifier_fixtures",
    TOOLS / "leopard2_exact_main_baseline_verifier_test.py",
)
acquire_fixtures = load_module(
    "leopard2_exact_main_build_stage_acquire_fixtures",
    TOOLS / "leopard2_exact_main_baseline_acquire_test.py",
)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def byte_identity(path: str, content: bytes) -> dict:
    return {"relative_path": path, "size": len(content),
            "sha256": sha256(content)}


class ScriptedBuildEnvironment(acquire.HostEnvironment):
    """Materialize deterministic fixture outputs for exact frozen argv only."""

    def __init__(self, plan: acquire.LanePlan, *, contaminated_role: str | None = None,
                 failure_role: str | None = None,
                 dependency_path: str | None = None,
                 mutate_attested_executable: bool = False) -> None:
        super().__init__(anchor_path="/tmp",
                         canonical_lock_path="/tmp/leopard-build-stage.lock")
        self.plan = plan
        self.contaminated_role = contaminated_role
        self.failure_role = failure_role
        self.dependency_path = dependency_path
        self.mutate_attested_executable = mutate_attested_executable
        self.calls: list[list[str]] = []

    def now_utc(self) -> str:
        return "2026-08-30T08:00:00Z"

    def _role(self, build_root: str) -> str:
        for role in record_contract.BUILD_ROLES:
            if acquire.build_role_roots(self.plan, role)["build_root"] == \
                    build_root:
                return role
        raise AssertionError(f"unknown scripted build root {build_root}")

    def _configure(self, arguments: list[str]) -> acquire.CommandResult:
        build_root = arguments[arguments.index("-B") + 1]
        role = self._role(build_root)
        if self.failure_role == role:
            return acquire.CommandResult(7, b"", b"scripted configure failure\n")
        roots = acquire.build_role_roots(self.plan, role)
        Path(build_root, "CMakeCache.txt").write_bytes(
            verifier_fixtures.cache_bytes(roots))
        Path(build_root, "compile_commands.json").write_bytes(
            verifier_fixtures.compile_commands_bytes(roots))
        return acquire.CommandResult(0, b"configure complete\n", b"")

    def _build(self, arguments: list[str]) -> acquire.CommandResult:
        build_root = arguments[2]
        role = self._role(build_root)
        roots = acquire.build_role_roots(self.plan, role)
        rodata = b"portable selected-section payload"
        if self.contaminated_role == role:
            rodata += b"|" + roots["build_root"].encode("ascii")
        debug = b"canonical debug root" if role != "path_variant" else \
            b"path variant debug root differs"
        executable = identity_fixtures.synthetic_elf(
            rodata=rodata, debug=debug)
        executable_path = Path(build_root, "leopard_main_benchmark")
        executable_path.write_bytes(executable)
        executable_path.chmod(0o755)
        archive = (b"!<arch>\ncanonical archive\n" if role != "path_variant"
                   else b"!<arch>\npath variant archive\n")
        Path(build_root, "libleopard_main_exact.a").write_bytes(archive)
        cmake_files = Path(build_root, "CMakeFiles")
        cmake_files.mkdir(mode=0o755, exist_ok=True)
        Path(cmake_files, "fixture.make").write_bytes(b"fixture closure\n")
        return acquire.CommandResult(0, b"build complete\n", b"")

    def run(self, argv, *, cwd, env, timeout, maximum_bytes):
        arguments = list(argv)
        self.calls.append(arguments)
        if "-S" in arguments and "-B" in arguments:
            return self._configure(arguments)
        if len(arguments) >= 3 and arguments[1] == "--build":
            return self._build(arguments)
        if arguments[0] == record_fixtures.TOOL_PATHS["ldd"]:
            dependency_target = Path(
                self.plan.scratch_root, "fixture-libc.so.6")
            if not dependency_target.exists():
                dependency_target.write_bytes(b"fixture libc bytes\n")
            loader_target = Path(
                self.plan.scratch_root, "fixture-ld-linux.so.2")
            if not loader_target.exists():
                loader_target.write_bytes(b"fixture loader bytes\n")
            dependency_path = self.dependency_path or str(dependency_target)
            return acquire.CommandResult(
                0,
                b"linux-vdso.so.1 (0x00007fff00000000)\n"
                + f"libc.so.6 => {dependency_path} ".encode("ascii")
                + b"(0x0000700000000000)\n"
                + f"{loader_target} ".encode("ascii")
                + b"(0x0000700000100000)\n",
                b"",
            )
        if arguments[0] == record_fixtures.TOOL_PATHS["ctest"]:
            return acquire.CommandResult(
                0,
                b"Test project fixture\n"
                b"100% tests passed, 0 tests failed out of 1\n"
                b"Total Test time (real) = 9999.99 sec\n",
                b"",
            )
        if arguments[0].endswith("/leopard_main_benchmark"):
            if self.mutate_attested_executable:
                Path(arguments[0]).write_bytes(
                    b"mutated after benchmark attestation\n")
            return acquire.CommandResult(
                0, verifier_fixtures.benchmark_stdout(), b"")
        raise AssertionError(f"unexpected scripted command: {arguments!r}")


class ExactMainBuildStageTest(unittest.TestCase):
    def patched_baseline(self):
        stack = unittest.mock.patch.multiple(
            record_contract,
            BASELINE_COMMIT=verifier_fixtures._TEST_BASELINE_COMMIT,
            BASELINE_TREE=verifier_fixtures._TEST_BASELINE_TREE,
        )
        return stack

    def run_three_roles(
        self, plan: acquire.LanePlan, prepared: acquire.PreparedAcquisitionRoots,
        environment: ScriptedBuildEnvironment, toolchain: dict,
        controller_sha256: str,
    ) -> list[acquire.BuildStageResult]:
        results = []
        for index, role in enumerate(record_contract.BUILD_ROLES):
            if index == 1:
                prepared.reset_root("canonical_build_root")
            results.append(acquire.acquire_build_stage(
                environment, plan, prepared, role=role, toolchain=toolchain,
                controller_sha256=controller_sha256))
        return results

    def test_build_stage_records_are_accepted_and_normalize_paths(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-test.") as temporary, \
                self.patched_baseline():
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            toolchain = record_fixtures.toolchain_fixture()
            environment = ScriptedBuildEnvironment(plan)
            with acquire.prepare_acquisition_roots(plan) as prepared:
                results = self.run_three_roles(
                    plan, prepared, environment, toolchain, "a" * 64)
                for role, result in zip(record_contract.BUILD_ROLES, results):
                    self.assertEqual(record_contract.validate_exact_main_build(
                        result.build, role=role,
                        tools={item["role"]: item["resolved_path"]
                               for item in toolchain["tools"]}), result.build)
                    self.assertTrue(all(
                        root["occurrences"] == 0
                        for root in result.identity["path_string_census"][
                            "roots"]))
                self.assertEqual(
                    results[0].build["executable"]["sha256"],
                    results[1].build["executable"]["sha256"])
                self.assertNotEqual(
                    results[0].build["executable"]["sha256"],
                    results[2].build["executable"]["sha256"])
                self.assertEqual(len({result.identity["combined_sha256"]
                                      for result in results}), 1)

    def test_nonzero_selected_section_census_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-census.") as temporary, \
                self.patched_baseline():
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            environment = ScriptedBuildEnvironment(
                plan, contaminated_role="canonical_first")
            with acquire.prepare_acquisition_roots(plan) as prepared:
                with self.assertRaisesRegex(
                        acquire.BuildStageError, "retain an acquisition root"):
                    acquire.acquire_build_stage(
                        environment, plan, prepared,
                        role="canonical_first",
                        toolchain=record_fixtures.toolchain_fixture(),
                        controller_sha256="a" * 64)

    def test_runtime_dependency_symlink_preserves_ldd_path_and_hashes_target(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-dependency.") as temporary, \
                self.patched_baseline():
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            with acquire.prepare_acquisition_roots(plan) as prepared:
                dependency_target = Path(
                    plan.scratch_root, "fixture-soname-target.so.6")
                dependency_target.write_bytes(b"fixture soname target\n")
                dependency_link = Path(plan.scratch_root, "libc-soname.so.6")
                dependency_link.symlink_to(dependency_target)
                environment = ScriptedBuildEnvironment(
                    plan, dependency_path=str(dependency_link))
                result = acquire.acquire_build_stage(
                    environment, plan, prepared,
                    role="canonical_first",
                    toolchain=record_fixtures.toolchain_fixture(),
                    controller_sha256="a" * 64)
                dependency = next(
                    item for item in result.runtime["dependencies"]
                    if item["soname"] == "libc.so.6")
                self.assertEqual(dependency["path"], str(dependency_link))
                self.assertEqual(
                    dependency["sha256"], sha256(dependency_target.read_bytes()))

    def test_attested_executable_is_rebound_after_all_attestations(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-attestation.") as temporary, \
                self.patched_baseline():
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            environment = ScriptedBuildEnvironment(
                plan, mutate_attested_executable=True)
            with acquire.prepare_acquisition_roots(plan) as prepared:
                with self.assertRaisesRegex(
                        acquire.BuildStageError,
                        "attested executable changed during attestation"):
                    acquire.acquire_build_stage(
                        environment, plan, prepared,
                        role="canonical_first",
                        toolchain=record_fixtures.toolchain_fixture(),
                        controller_sha256="a" * 64)

    def test_ctest_summary_is_strict_and_never_interprets_duration(self) -> None:
        prefix = b"100% tests passed, 0 tests failed out of 1\n"
        self.assertEqual(acquire.parse_ctest_success_summary(
            prefix + b"Total Test time (real) = 0.01 sec\n"), (1, 0))
        self.assertEqual(acquire.parse_ctest_success_summary(
            prefix + b"Total Test time (real) = 999999 sec\n"), (1, 0))
        for changed in (
                b"0% tests passed, 1 tests failed out of 1\n",
                b"100% tests passed, 0 tests failed out of 2\n",
                b"Total Test time (real) = 0.01 sec\n"):
            with self.subTest(changed=changed):
                with self.assertRaises(acquire.AcquisitionError):
                    acquire.parse_ctest_success_summary(changed)

    def test_build_root_reset_and_artifact_staging_preserve_contracts(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-reset.") as temporary:
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            with acquire.prepare_acquisition_roots(plan) as prepared:
                build_root = Path(plan.canonical_build_root)
                nested = build_root / "nested"
                nested.mkdir()
                (nested / "file").write_bytes(b"content")
                (build_root / "link").symlink_to("nested/file")
                identity = prepared.identities["canonical_build_root"]
                with mock.patch.object(
                        acquire.shutil, "rmtree",
                        wraps=acquire.shutil.rmtree) as rmtree:
                    prepared.reset_root("canonical_build_root")
                self.assertEqual(rmtree.call_count, 1)
                self.assertEqual(
                    rmtree.call_args.args, (str(nested),))
                self.assertNotIn("dir_fd", rmtree.call_args.kwargs)
                self.assertEqual(os.listdir(build_root), [])
                self.assertEqual(
                    prepared.identities["canonical_build_root"], identity)
                with self.assertRaises(acquire.AcquisitionError):
                    prepared.reset_root("scratch_root")

                source = build_root / "executable"
                source.write_bytes(b"fixture executable")
                source.chmod(0o755)
                destination = Path(plan.scratch_root) / "staged-executable"
                result = acquire.stage_build_output(
                    str(source), str(destination), "fixture executable")
                status = destination.stat()
                self.assertEqual(stat.S_IMODE(status.st_mode), 0o600)
                self.assertEqual(status.st_nlink, 1)
                self.assertEqual(result["sha256"], sha256(source.read_bytes()))
                staged_content = destination.read_bytes()
                with self.assertRaisesRegex(
                        acquire.AcquisitionError, "already exists"):
                    acquire.stage_build_output(
                        str(source), str(destination), "fixture executable")
                self.assertEqual(destination.read_bytes(), staged_content)

    def test_build_stage_failure_prefixes_seal_and_verify(self) -> None:
        for stage_index in (1, 2, 3):
            with self.subTest(stage_index=stage_index), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-build-stage-failure.") as temporary:
                plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
                environment = ScriptedBuildEnvironment(plan)
                stage_names = record_contract.STAGE_NAMES[:stage_index + 1]
                stage_logs = {
                    name: (f"{name} {'failed' if index == stage_index else 'complete'}\n").
                    encode("ascii")
                    for index, name in enumerate(stage_names)
                }
                role = record_contract.BUILD_ROLES[stage_index - 1]
                command_error = acquire.CommandExecutionError(
                    "scripted role failure", 9)
                error = acquire.BuildStageError(
                    role, command_error, stage_logs[stage_names[-1]])
                retained_bytes = ({
                    "attestations/canonical_first/benchmark.stderr": b"",
                } if stage_index >= 2 else None)
                result = acquire.seal_stage_failure(
                    environment, plan, error, stage=stage_names[-1],
                    stage_logs=stage_logs,
                    retained_bytes=retained_bytes,
                    diagnostics={"diagnostics/stderr": b"failed\n"})
                self.assertEqual(result["terminal"], "FAILED.json")
                verdict = verifier.verify_sealed_lane(plan.lane_root)
                self.assertEqual(verdict["outcome"], "verified_failure")
                self.assertFalse(verdict["promoted"])
                process = subprocess.run(
                    [sys.executable, "-I", "-S", "-B",
                     str(TOOLS / "leopard2_exact_main_baseline_verifier.py"),
                     plan.lane_root],
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=False)
                self.assertEqual(process.returncode, 3, process.stderr)
                failure = json.loads(
                    Path(plan.lane_root, "FAILED.json").read_bytes())
                self.assertEqual(failure["stage"], stage_names[-1])
                if retained_bytes is not None:
                    self.assertEqual(
                        Path(plan.lane_root,
                             "attestations/canonical_first/benchmark.stderr").
                        read_bytes(), b"")
                self.assertFalse(any(
                    Path(plan.lane_root, f"logs/{index:02d}-{name}.log").exists()
                    for index, name in enumerate(record_contract.STAGE_NAMES)
                    if index > stage_index))

    def test_synthetic_authority_from_three_stages_verifies_offline(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-authority.") as temporary, \
                self.patched_baseline(), \
                mock.patch.object(
                    verifier, "BASELINE_COMMIT",
                    verifier_fixtures._TEST_BASELINE_COMMIT), \
                mock.patch.object(
                    verifier, "BASELINE_TREE",
                    verifier_fixtures._TEST_BASELINE_TREE):
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            roots_patch = mock.patch.object(
                verifier_fixtures, "roots",
                side_effect=lambda role: acquire.build_role_roots(plan, role))
            with roots_patch:
                fixture_files, fixture_record = \
                    verifier_fixtures.authority_files_and_record()
            toolchain = fixture_record["toolchain"]
            controller = fixture_record["attestation"]["test_controller"]
            environment = ScriptedBuildEnvironment(plan)
            with acquire.prepare_acquisition_roots(plan) as prepared:
                results = self.run_three_roles(
                    plan, prepared, environment, toolchain,
                    controller["sha256"])
                retained_bytes = {
                    path: content for path, content in fixture_files.items()
                    if path.startswith(("source/", "toolchain/", "controllers/"))
                }
                stage_logs = {
                    "source_acquisition":
                        fixture_files["logs/00-source_acquisition.log"],
                    "canonical_first_build": results[0].log,
                    "canonical_second_build": results[1].log,
                    "path_variant_build": results[2].log,
                    "independent_verification": fixture_files[
                        "logs/04-independent_verification.log"],
                    "seal": fixture_files["logs/05-seal.log"],
                }
                retained_paths: dict[str, str] = {}
                for result in results:
                    retained_bytes.update(result.retained_bytes)
                    retained_paths.update(result.retained_paths)
                for index, name in enumerate(record_contract.STAGE_NAMES):
                    retained_bytes[f"logs/{index:02d}-{name}.log"] = \
                        stage_logs[name]
                lane = {
                    "root": plan.lane_root,
                    "attempt": plan.attempt,
                    "attempt_budget": 3,
                    "record_relative_path": "baseline-authority.json",
                    "seal_protocol": record_contract.SEAL_PROTOCOL,
                    "stages": [{
                        "name": name,
                        "status": "complete",
                        "log": byte_identity(
                            f"logs/{index:02d}-{name}.log", stage_logs[name]),
                    } for index, name in enumerate(record_contract.STAGE_NAMES)],
                }
                identities = {result.role: result.identity for result in results}
                authority = record_contract.baseline_authority_record(
                    created_utc=environment.now_utc(),
                    host=fixture_record["host"], lane=lane,
                    source=fixture_record["source"],
                    adapter=fixture_record["adapter"], toolchain=toolchain,
                    builds={result.role: result.build for result in results},
                    runtime_closure={
                        "schema": record_contract.RUNTIME_CLOSURE_SCHEMA,
                        "normalization":
                            record_contract.CANONICAL_LDD_NORMALIZATION,
                        "records": [result.runtime for result in results],
                    },
                    attestation={
                        "schema": record_contract.ATTESTATION_SCHEMA,
                        "test_controller": controller,
                        "records": [result.attestation for result in results],
                    },
                    identity={
                        **identities,
                        "combined_sha256":
                            results[0].identity["combined_sha256"],
                        "canonical_raw_bytes_identical": True,
                        "path_variant_raw_bytes_differ": True,
                        "normalized_match": True,
                    },
                )
                with acquire.LaneWriter(plan.lane_root) as writer:
                    writer.seal_record(
                        authority, retained_bytes,
                        retained_paths=retained_paths)

                verdict = verifier.verify_sealed_lane(plan.lane_root)
                self.assertEqual(verdict["outcome"], "promoted_authority")
                self.assertTrue(verdict["promoted"])
                process = subprocess.run(
                    verifier_fixtures.verifier_cli_command(Path(plan.lane_root)),
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=False)
                self.assertEqual(process.returncode, 0, process.stderr)

    def test_stage_logs_have_no_timing_fields(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-build-stage-log.") as temporary, \
                self.patched_baseline():
            plan = acquire_fixtures.temporary_lane_plan(Path(temporary))
            with acquire.prepare_acquisition_roots(plan) as prepared:
                result = acquire.acquire_build_stage(
                    ScriptedBuildEnvironment(plan), plan, prepared,
                    role="canonical_first",
                    toolchain=record_fixtures.toolchain_fixture(),
                    controller_sha256="a" * 64)
                stage = json.loads(result.log)
                rendered = json.dumps(stage, sort_keys=True)
                self.assertNotIn("duration", rendered)
                self.assertNotIn("timestamp", rendered)
                for path in ("builds/canonical-first/configure.log",
                             "builds/canonical-first/build.log"):
                    transcript = json.loads(result.retained_bytes[path])
                    self.assertNotIn("duration", transcript)
                    self.assertNotIn("timestamp", transcript)


if __name__ == "__main__":
    unittest.main()
