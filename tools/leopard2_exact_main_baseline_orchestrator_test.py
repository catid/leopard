#!/usr/bin/env python3
"""Deterministic end-to-end tests for exact-main acquisition orchestration."""

from __future__ import annotations

import contextlib
import copy
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
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


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


verifier_fixtures = load_module(
    "leopard2_exact_main_orchestrator_verifier_fixtures",
    TOOLS / "leopard2_exact_main_baseline_verifier_test.py",
)
acquire_fixtures = load_module(
    "leopard2_exact_main_orchestrator_acquire_fixtures",
    TOOLS / "leopard2_exact_main_baseline_acquire_test.py",
)
build_fixtures = load_module(
    "leopard2_exact_main_orchestrator_build_fixtures",
    TOOLS / "leopard2_exact_main_baseline_build_stage_test.py",
)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


class OrchestratorEnvironment(build_fixtures.ScriptedBuildEnvironment):
    """Use synthetic builds but a real isolated verifier subprocess."""

    def __init__(
        self,
        plan: acquire.LanePlan,
        parent: Path,
        host: dict,
        *,
        failure_role: str | None = None,
        authority_verifier_status: int | None = None,
    ) -> None:
        super().__init__(plan, failure_role=failure_role)
        anchor = parent / "lease-anchor"
        anchor.mkdir(mode=0o700)
        self.anchor_path = str(anchor)
        self.canonical_lock_path = str(parent / "canonical.lock")
        self.fixture_host = copy.deepcopy(host)
        self.authority_verifier_status = authority_verifier_status
        self.verifier_lanes: list[str] = []

    def host_facts(self) -> dict:
        return copy.deepcopy(self.fixture_host)

    def independent_verifier_argv(
        self, *, interpreter: str, verifier: str, lane_root: str,
    ) -> list[str]:
        arguments = verifier_fixtures.verifier_cli_command(Path(lane_root))
        arguments[0] = interpreter
        return arguments

    def run(self, argv, *, cwd, env, timeout, maximum_bytes):
        arguments = list(argv)
        lane = arguments[-1] if arguments else ""
        if lane in (self.plan.lane_root, self.plan.failure_lane_root) and \
                "-c" in arguments:
            self.verifier_lanes.append(lane)
            if lane == self.plan.lane_root and \
                    Path(lane, "baseline-authority.json").is_file() and \
                    self.authority_verifier_status is not None:
                status = self.authority_verifier_status
                return acquire.CommandResult(
                    status, b"", b"scripted authority verifier failure\n")
            return acquire.HostEnvironment.run(
                self, arguments, cwd=cwd, env=env, timeout=timeout,
                maximum_bytes=maximum_bytes)
        return super().run(
            arguments, cwd=cwd, env=env, timeout=timeout,
            maximum_bytes=maximum_bytes)


class FixtureAcquisition:
    def __init__(self, parent: Path):
        self.parent = parent
        self.plan = acquire_fixtures.temporary_lane_plan(parent)
        roots_patch = mock.patch.object(
            verifier_fixtures, "roots",
            side_effect=lambda role: acquire.build_role_roots(self.plan, role))
        anchors_patch = mock.patch.multiple(
            record_contract,
            BASELINE_COMMIT=verifier_fixtures._TEST_BASELINE_COMMIT,
            BASELINE_TREE=verifier_fixtures._TEST_BASELINE_TREE,
        )
        with roots_patch, anchors_patch:
            self.files, self.fixture_record = \
                verifier_fixtures.authority_files_and_record()
        self.source = copy.deepcopy(self.fixture_record["source"])
        self.adapter = copy.deepcopy(self.fixture_record["adapter"])
        self.toolchain = copy.deepcopy(self.fixture_record["toolchain"])
        program = acquire.verification_program_identity(self.plan, None)
        python = program["files"][-1]
        tool = next(item for item in self.toolchain["tools"]
                    if item["role"] == "python")
        tool.update({
            "path": python["path"], "resolved_path": python["path"],
            "size": python["size"], "mode": python["mode"],
            "sha256": python["sha256"],
        })
        version = next(item for item in self.toolchain["versions"]
                       if item["role"] == "python")
        version["argv"] = [python["path"], "--version"]
        self.archive_paths = {
            self.source["baseline"]["archive"]["relative_path"],
            self.source["adapter_repository"]["archive"]["relative_path"],
        }
        self.retained_bytes = {
            path: content for path, content in self.files.items()
            if path.startswith(("source/", "toolchain/", "controllers/"))
            and path not in self.archive_paths
        }

    @contextlib.contextmanager
    def anchors(self):
        with mock.patch.multiple(
                record_contract,
                BASELINE_COMMIT=verifier_fixtures._TEST_BASELINE_COMMIT,
                BASELINE_TREE=verifier_fixtures._TEST_BASELINE_TREE):
            yield

    def source_stage(self, _environment, plan, prepared):
        self.assert_plan(plan, prepared)
        retained_paths: dict[str, str] = {}
        for relative in sorted(self.archive_paths):
            destination = Path(plan.scratch_root) / (
                "fixture-" + Path(relative).name)
            destination.write_bytes(self.files[relative])
            destination.chmod(0o600)
            retained_paths[relative] = str(destination)
        return acquire.SourceStageResult(
            source=copy.deepcopy(self.source),
            adapter=copy.deepcopy(self.adapter),
            toolchain=copy.deepcopy(self.toolchain),
            retained_bytes=copy.deepcopy(self.retained_bytes),
            retained_paths=retained_paths,
            log=self.files["logs/00-source_acquisition.log"],
        )

    def assert_plan(self, plan, prepared) -> None:
        if plan != self.plan or prepared.plan != self.plan:
            raise AssertionError("orchestrator changed its source-stage plan")


class ExactMainBaselineOrchestratorTest(unittest.TestCase):
    def run_fixture(
        self,
        fixture: FixtureAcquisition,
        environment: OrchestratorEnvironment,
        *,
        build_side_effect=None,
    ) -> acquire.AcquisitionOutcome:
        source_patch = mock.patch.object(
            acquire, "acquire_source_stage",
            side_effect=fixture.source_stage)
        build_patch = contextlib.nullcontext()
        if build_side_effect is not None:
            build_patch = mock.patch.object(
                acquire, "acquire_build_stage",
                side_effect=build_side_effect)
        with fixture.anchors(), source_patch, build_patch:
            return acquire.acquire_exact_main_baseline(
                environment, fixture.plan, mode="rehearsal")

    def test_success_seals_and_externally_verifies_without_timing(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-success.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            outcome = self.run_fixture(fixture, environment)
            self.assertEqual(outcome.outcome, "promoted_authority")
            self.assertEqual(outcome.verification_exit_status, 0)
            self.assertEqual(environment.verifier_lanes,
                             [fixture.plan.lane_root])
            self.assertTrue(Path(
                fixture.plan.lane_root, "baseline-authority.json").is_file())
            self.assertFalse(Path(fixture.plan.failure_lane_root).exists())
            for index, name in enumerate(record_contract.STAGE_NAMES):
                self.assertTrue(Path(
                    fixture.plan.lane_root,
                    f"logs/{index:02d}-{name}.log").is_file())
            for field in acquire._PREPARED_ROOT_FIELDS:
                self.assertFalse(Path(getattr(fixture.plan, field)).exists())
            report = acquire.acquisition_report(
                fixture.plan, outcome, mode="rehearsal")
            rendered = identity_contract.canonical_json_bytes(report)
            self.assertNotIn(b"duration", rendered)
            self.assertNotIn(b"elapsed", rendered)
            self.assertNotIn(b"throughput", rendered)
            for path in Path(fixture.plan.lane_root).rglob("*"):
                if path.is_file():
                    self.assertNotIn(b"rehearsal", path.read_bytes())

    def test_build_failures_seal_exact_prefix_and_verify(self) -> None:
        for role_index, role in enumerate(record_contract.BUILD_ROLES, 1):
            with self.subTest(role=role), tempfile.TemporaryDirectory(
                    prefix="leopard-orchestrator-build-failure.") as temporary:
                parent = Path(temporary)
                fixture = FixtureAcquisition(parent)
                environment = OrchestratorEnvironment(
                    fixture.plan, parent, fixture.fixture_record["host"])
                original = acquire.acquire_build_stage

                def fail_selected(*arguments, **keywords):
                    if keywords["role"] == role:
                        command_error = acquire.CommandExecutionError(
                            "scripted role failure", 7)
                        raise acquire.BuildStageError(
                            role, command_error,
                            f"{role} failed\n".encode("ascii"))
                    return original(*arguments, **keywords)

                outcome = self.run_fixture(
                    fixture, environment, build_side_effect=fail_selected)
                self.assertEqual(outcome.outcome, "acquisition_failure")
                self.assertEqual(outcome.verification_exit_status, 3)
                self.assertEqual(environment.verifier_lanes,
                                 [fixture.plan.lane_root])
                failure = record_contract.load_baseline_failure_record(
                    Path(fixture.plan.lane_root, "FAILED.json").read_bytes())
                self.assertEqual(failure["stage"], role + "_build")
                self.assertEqual(len(failure["lane"]["stages"]),
                                 role_index + 1)
                self.assertFalse(Path(fixture.plan.failure_lane_root).exists())

    def test_gate_failure_never_publishes_authority(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-gate-failure.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            original = acquire.acquire_build_stage
            canonical_first_identity: dict | None = None

            def changed_identity(*arguments, **keywords):
                nonlocal canonical_first_identity
                result = original(*arguments, **keywords)
                if result.role == "canonical_first":
                    canonical_first_identity = copy.deepcopy(result.identity)
                    return result
                if result.role == "canonical_second":
                    if canonical_first_identity is None:
                        raise AssertionError("canonical first role was skipped")
                    changed = copy.deepcopy(result.identity)
                    changed["record_sha256"] = "0" * 64
                    return acquire.BuildStageResult(
                        role=result.role, build=result.build,
                        runtime=result.runtime, attestation=result.attestation,
                        identity=changed,
                        retained_bytes=result.retained_bytes,
                        retained_paths=result.retained_paths, log=result.log)
                return result

            outcome = self.run_fixture(
                fixture, environment, build_side_effect=changed_identity)
            self.assertEqual(outcome.outcome, "acquisition_failure")
            self.assertFalse(Path(
                fixture.plan.lane_root, "baseline-authority.json").exists())
            failure = record_contract.load_baseline_failure_record(
                Path(fixture.plan.lane_root, "FAILED.json").read_bytes())
            self.assertEqual(failure["stage"], "path_variant_build")
            self.assertTrue(Path(
                fixture.plan.lane_root,
                "diagnostics/completed-path-variant-build.log").is_file())

    def test_verifier_failure_seals_bound_second_lane(self) -> None:
        for status in (1, 2, 3):
            with self.subTest(status=status), tempfile.TemporaryDirectory(
                    prefix="leopard-orchestrator-verify-failure.") as temporary:
                parent = Path(temporary)
                fixture = FixtureAcquisition(parent)
                environment = OrchestratorEnvironment(
                    fixture.plan, parent, fixture.fixture_record["host"],
                    authority_verifier_status=status)
                outcome = self.run_fixture(fixture, environment)
                self.assertEqual(outcome.outcome, "verification_failure")
                self.assertEqual(environment.verifier_lanes, [
                    fixture.plan.lane_root, fixture.plan.failure_lane_root])
                with fixture.anchors():
                    authority = record_contract.load_baseline_authority_record(
                        Path(fixture.plan.lane_root,
                             "baseline-authority.json").read_bytes())
                    failure = record_contract.load_baseline_failure_record(
                        Path(fixture.plan.failure_lane_root,
                             "FAILED.json").read_bytes())
                self.assertEqual(
                    failure["authority_record_sha256"],
                    authority["record_sha256"])
                self.assertEqual(failure["stage"],
                                 "independent_verification")
                self.assertEqual(outcome.verification_exit_status, 3)

    def test_source_failure_is_a_verified_stage_zero_terminal(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-source-failure.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            error = acquire.SourceStageError(
                acquire.AcquisitionError("scripted source failure"),
                b"source stage failed\n")
            with fixture.anchors(), mock.patch.object(
                    acquire, "acquire_source_stage", side_effect=error):
                outcome = acquire.acquire_exact_main_baseline(
                    environment, fixture.plan, mode="rehearsal")
            self.assertEqual(outcome.outcome, "acquisition_failure")
            self.assertEqual(outcome.verification_exit_status, 3)
            self.assertEqual(environment.verifier_lanes,
                             [fixture.plan.lane_root])
            failure = record_contract.load_baseline_failure_record(
                Path(fixture.plan.lane_root, "FAILED.json").read_bytes())
            self.assertEqual(failure["stage"], "source_acquisition")

    def test_authority_seal_failure_uses_second_lane_once(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-seal-failure.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            original = acquire.LaneWriter.seal_record

            def fail_primary(writer, record, retained_files, **keywords):
                if writer.root_path == fixture.plan.lane_root:
                    writer.publish_bytes("torn-primary", b"torn\n")
                    raise acquire.AcquisitionError(
                        "scripted primary seal failure")
                return original(writer, record, retained_files, **keywords)

            with mock.patch.object(
                    acquire.LaneWriter, "seal_record",
                    autospec=True, side_effect=fail_primary):
                outcome = self.run_fixture(fixture, environment)
            self.assertEqual(outcome.outcome, "acquisition_failure")
            self.assertEqual(outcome.lane_root,
                             fixture.plan.failure_lane_root)
            self.assertEqual(environment.verifier_lanes,
                             [fixture.plan.failure_lane_root])
            self.assertTrue(Path(
                fixture.plan.lane_root, "torn-primary").is_file())
            self.assertFalse(Path(
                fixture.plan.lane_root, "baseline-authority.json").exists())
            failure = record_contract.load_baseline_failure_record(
                Path(fixture.plan.failure_lane_root,
                     "FAILED.json").read_bytes())
            self.assertEqual(failure["stage"], "path_variant_build")

    def test_late_writer_error_recovers_only_verified_authority(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-late-seal.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            original = acquire.LaneWriter.seal_record

            def fail_after_complete(writer, record, retained_files, **keywords):
                result = original(
                    writer, record, retained_files, **keywords)
                if writer.root_path == fixture.plan.lane_root:
                    raise acquire.AcquisitionError(
                        "scripted post-seal self-check failure")
                return result

            with mock.patch.object(
                    acquire.LaneWriter, "seal_record", autospec=True,
                    side_effect=fail_after_complete):
                outcome = self.run_fixture(fixture, environment)
            self.assertEqual(outcome.outcome, "promoted_authority")
            self.assertEqual(environment.verifier_lanes,
                             [fixture.plan.lane_root])
            self.assertFalse(Path(fixture.plan.failure_lane_root).exists())
            self.assertEqual(
                outcome.seal["file_count"],
                outcome.verdict["seal"]["file_count"])

    def test_unclassifiable_late_seal_never_emits_contrary_failure(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-unclassified-seal.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            original_seal = acquire.LaneWriter.seal_record
            original_run = environment.run

            def fail_after_complete(writer, record, retained_files, **keywords):
                result = original_seal(
                    writer, record, retained_files, **keywords)
                if writer.root_path == fixture.plan.lane_root:
                    raise acquire.AcquisitionError(
                        "scripted post-seal self-check failure")
                return result

            def verifier_cannot_run(
                    argv, *, cwd, env, timeout, maximum_bytes):
                arguments = list(argv)
                if arguments[-1] == fixture.plan.lane_root and \
                        "-c" in arguments:
                    raise acquire.AcquisitionError(
                        "scripted verifier spawn failure")
                return original_run(
                    arguments, cwd=cwd, env=env, timeout=timeout,
                    maximum_bytes=maximum_bytes)

            with mock.patch.object(
                    acquire.LaneWriter, "seal_record", autospec=True,
                    side_effect=fail_after_complete), mock.patch.object(
                        environment, "run", side_effect=verifier_cannot_run), \
                    self.assertRaisesRegex(
                        acquire.AcquisitionError,
                        "could not be independently classified"):
                self.run_fixture(fixture, environment)
            self.assertTrue(Path(
                fixture.plan.lane_root, "baseline-authority.json").is_file())
            self.assertFalse(Path(fixture.plan.failure_lane_root).exists())
            completed = subprocess.run(
                verifier_fixtures.verifier_cli_command(
                    Path(fixture.plan.lane_root)),
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_plain_post_seal_exception_becomes_verification_failure(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-post-seal.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            original = acquire.run_independent_verification

            def fail_authority(*arguments, **keywords):
                if keywords.get("lane_root") is None:
                    raise acquire.AcquisitionError(
                        "scripted post-seal invariant failure")
                return original(*arguments, **keywords)

            with mock.patch.object(
                    acquire, "run_independent_verification",
                    side_effect=fail_authority):
                outcome = self.run_fixture(fixture, environment)
            self.assertEqual(outcome.outcome, "verification_failure")
            self.assertTrue(Path(
                fixture.plan.lane_root, "baseline-authority.json").is_file())
            failure = record_contract.load_baseline_failure_record(
                Path(fixture.plan.failure_lane_root,
                     "FAILED.json").read_bytes())
            self.assertEqual(failure["stage"],
                             "independent_verification")
            self.assertEqual(failure["authority_record_sha256"],
                             outcome.authority_record_sha256)

    def test_existing_lane_or_failure_root_is_never_reused(self) -> None:
        for field in ("lane_root", "failure_lane_root"):
            with self.subTest(field=field), tempfile.TemporaryDirectory(
                    prefix="leopard-orchestrator-no-reuse.") as temporary:
                parent = Path(temporary)
                fixture = FixtureAcquisition(parent)
                Path(getattr(fixture.plan, field)).mkdir(mode=0o700)
                environment = OrchestratorEnvironment(
                    fixture.plan, parent, fixture.fixture_record["host"])
                with fixture.anchors(), mock.patch.object(
                        acquire, "acquire_source_stage",
                        side_effect=fixture.source_stage), self.assertRaises(
                            acquire.AcquisitionError):
                    acquire.acquire_exact_main_baseline(
                        environment, fixture.plan, mode="rehearsal")
                self.assertEqual(environment.calls, [])
                self.assertEqual(environment.verifier_lanes, [])

    def test_cli_is_strict_and_report_is_canonical(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-orchestrator-cli.") as temporary:
            parent = Path(temporary)
            fixture = FixtureAcquisition(parent)
            environment = OrchestratorEnvironment(
                fixture.plan, parent, fixture.fixture_record["host"])
            outcome = self.run_fixture(fixture, environment)
            arguments = [
                "--mode", "rehearsal", "--attempt", "1",
                "--lane-root", fixture.plan.lane_root,
                "--failure-lane-root", fixture.plan.failure_lane_root,
                "--repository", fixture.plan.repository,
                "--verifier", fixture.plan.verifier,
                "--canonical-adapter-root",
                fixture.plan.canonical_adapter_root,
                "--canonical-baseline-root",
                fixture.plan.canonical_baseline_root,
                "--canonical-build-root", fixture.plan.canonical_build_root,
                "--variant-adapter-root", fixture.plan.variant_adapter_root,
                "--variant-baseline-root", fixture.plan.variant_baseline_root,
                "--variant-build-root", fixture.plan.variant_build_root,
                "--scratch-root", fixture.plan.scratch_root,
            ]

            class Stdout:
                def __init__(self):
                    self.buffer = io.BytesIO()

            stdout = Stdout()
            stderr = io.StringIO()
            with mock.patch.object(
                    acquire, "acquire_exact_main_baseline",
                    return_value=outcome), \
                    mock.patch.object(sys, "stdout", stdout), \
                    contextlib.redirect_stderr(stderr):
                status = acquire.main(arguments)
            self.assertEqual(status, 0, stderr.getvalue())
            report = identity_contract.strict_json_loads(
                stdout.buffer.getvalue(), "acquisition CLI report")
            self.assertEqual(
                identity_contract.canonical_json_bytes(report),
                stdout.buffer.getvalue())
            self.assertEqual(report["mode"], "rehearsal")
            for changed in (
                    arguments[:-2],
                    ["--unknown", "x"] + arguments[2:],
                    arguments + ["--mode", "rehearsal"]):
                with self.subTest(changed=changed[:2]):
                    with mock.patch.object(sys, "stderr", io.StringIO()):
                        self.assertEqual(acquire.main(changed), 2)


if __name__ == "__main__":
    unittest.main()
