#!/usr/bin/env python3
"""Adversarial tests for the closed K65 generation-3 live authority.

These tests launch no benchmark and acquire no qualification scan or timing.
The production pipeline is exercised with pure retained fixtures and a spawn
canary; the sequential executor uses an in-memory child result canary.
"""

from __future__ import annotations

import contextlib
import copy
import errno
import fcntl
import hashlib
import importlib.util
import inspect
import os
from pathlib import Path
import runpy
import shutil
import signal
import tempfile
import threading
import types
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
RUNNER_PATH = HERE / "run_k65r65_b64_packed_terminal_gen3_acquire.py"
BOOTSTRAP_PATH = HERE / "k65_gen3_exact_source_bootstrap.py"
PURE_TEST_PATH = HERE / "test_k65_gen3_execution_contract.py"
SEALED_AUTHORITY = Path(
    "/home/catid/leopard-exact-main-authority-6d4f690-a1/Aauthority")
AUTHORITY_AVAILABLE = all((
    (SEALED_AUTHORITY / "baseline-authority.json").is_file(),
    (SEALED_AUTHORITY / "SHA256SUMS").is_file(),
    (SEALED_AUTHORITY /
     "artifacts/path-variant/leopard_main_benchmark").is_file(),
))


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


bootstrap_namespace = runpy.run_path(
    str(BOOTSTRAP_PATH), run_name="k65_gen3_direct_source_bootstrap_under_test",
)
load_live_armer = bootstrap_namespace["load_live_armer"]
runner = load_live_armer("k65_gen3_live_armer_under_test")
pure_fixtures = load_module("k65_gen3_live_armer_pure_fixtures", PURE_TEST_PATH)
contract = runner.contract
execution = runner.execution


def private_lane(directory: str, name: str = "lane") -> Path:
    lane = (Path(directory) / name).resolve()
    lane.mkdir(mode=0o700)
    lane.chmod(0o700)
    return lane


def create_regular_slot_for_test(
    directory_fd: int, name: str, binding_factory,
) -> dict:
    descriptor = os.open(
        name, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
        os.O_NOFOLLOW, 0o600, dir_fd=directory_fd)
    try:
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        return binding_factory(os.fstat(descriptor), name)
    finally:
        os.close(descriptor)


def create_campaign_marker_slot_for_test(
    directory_fd: int, name: str,
) -> dict:
    return create_regular_slot_for_test(
        directory_fd, name, runner._campaign_marker_binding_record)


def create_prearmed_history_slot_for_test(
    directory_fd: int, name: str,
) -> dict:
    return create_regular_slot_for_test(
        directory_fd, name, runner._prearmed_history_marker_binding_record)


def create_prearmed_boundary_slot_for_test(
    directory_fd: int, name: str,
) -> dict:
    return create_regular_slot_for_test(
        directory_fd, name, runner._prearmed_boundary_marker_binding_record)


def launch_context() -> dict:
    return runner._validate_launch_context(
        runner.prereg.recommended_launch_context_record())


def synthetic_elf(size: int, marker: int = 0x5a) -> bytes:
    if size < 4:
        raise ValueError("synthetic ELF size is too small")
    return b"\x7fELF" + bytes((marker,)) * (size - 4)


_authority_cache = None


def retained_authority():
    global _authority_cache
    if _authority_cache is None:
        _authority_cache = runner.bind_exact_main_path_variant(SEALED_AUTHORITY)
    return _authority_cache


class FakeCanonicalLock:
    def __init__(self, unused_path: str):
        self.held = False
        self.validations = 0

    def __enter__(self):
        if self.held:
            raise ValueError("already held")
        self.held = True
        return self

    def validate_current(self) -> None:
        if not self.held:
            raise ValueError("not held")
        self.validations += 1

    def __exit__(self, unused_kind, unused_value, unused_traceback) -> None:
        self.held = False


class FakeRuntimeGuard:
    def __init__(self):
        self.closed = False
        self.validations = 0

    def validate_current(self, unused_label: str) -> None:
        if self.closed:
            raise ValueError("runtime guard is closed")
        self.validations += 1

    def close(self) -> None:
        self.closed = True


_lane_binding_cache: dict[tuple[str, int, int], dict] = {}


def lane_binding_for_test(lane: Path) -> dict:
    resolved = lane.resolve(strict=True)
    metadata = lane.stat(follow_symlinks=False)
    key = (str(resolved), metadata.st_dev, metadata.st_ino)
    binding = _lane_binding_cache.get(key)
    if binding is None:
        binding = runner.prereg.output_lane_binding_record(
            resolved, setup_invalid_budget=5,
            environment_rejected_budget=8, evidence_attempt_budget=3)
        _lane_binding_cache[key] = copy.deepcopy(binding)
    return copy.deepcopy(binding)


def live_fixture(
    expected_frozen_pair: dict | None = None, *, output_lane: Path | None = None,
) -> dict:
    candidate_data = synthetic_elf(32_768, 0x61)
    candidate_sha = hashlib.sha256(candidate_data).hexdigest()
    candidate_record = {"schema": "unit-candidate-authority/v1", "fixed": True}
    candidate_record_data = contract.canonical_json_bytes(candidate_record)
    ledger_sha = "7d" * 32
    controllers = runner.prereg.current_controller_bindings()
    pure = pure_fixtures.fixture(
        candidate_sha,
        candidate_authority_record_sha256=
            hashlib.sha256(candidate_record_data).hexdigest(),
        candidate_authority_ledger_sha256=ledger_sha,
        controller_bindings=controllers,
        output_lane_binding=(
            lane_binding_for_test(output_lane)
            if output_lane is not None else None),
        expected_frozen_pair=expected_frozen_pair)
    candidate_authority = {
        "root": Path("/candidate-authority"),
        "record": candidate_record,
        "record_data": candidate_record_data,
        "ledger": {"sha256": ledger_sha},
        "artifact_data": candidate_data,
        "runtime_closure": {
            "schema": "leopard2-k65-gen3-runtime-closure/v1",
            "dependencies": [], "launchers": [],
        },
        "source_authority": pure["source_authority"],
    }
    exact_main = {
        "root": Path("/exact-main-authority"),
        "record": {"schema": "unit-exact-main-authority/v1"},
        "record_data": b"{}",
        "verdict": {"verdict_sha256": "6" * 64},
        "artifact_data": b"unit-main-bytes",
    }
    qualification = {
        "host_instance": pure["host"],
        "qualification_evidence": pure["evidence"],
        "qualification_binding": pure["qualification"],
    }
    return {
        "pure": pure,
        "pure_fixture_args": {
            "candidate_sha256": candidate_sha,
            "candidate_authority_record_sha256":
                hashlib.sha256(candidate_record_data).hexdigest(),
            "candidate_authority_ledger_sha256": ledger_sha,
            "controller_bindings": controllers,
            "expected_frozen_pair": copy.deepcopy(expected_frozen_pair),
        },
        "candidate_authority": candidate_authority,
        "exact_main": exact_main,
        "qualification": qualification,
    }


def fake_capture_factory():
    counter = 0

    def capture(
        unused_data: bytes, *, expected_sha256: str, expected_size: int,
        label: str, authority_relative_path: str,
    ):
        nonlocal counter
        counter += 1
        descriptor = os.memfd_create(
            f"gen3-test-{counter}", os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING)
        os.write(descriptor, b"\x7fELF")
        os.fchmod(descriptor, 0o500)
        metadata = os.fstat(descriptor)
        snapshot = {
            "protocol": runner.SEALED_EXECUTABLE_PROTOCOL,
            "boot_id": "12345678-1234-1234-1234-123456789abc",
            "session_nonce": f"{counter:064x}",
            "creator_pid": os.getpid(),
            "creator_start_ticks": runner._process_start_ticks(),
            "device": metadata.st_dev, "inode": metadata.st_ino,
            "uid": metadata.st_uid, "gid": metadata.st_gid,
            "size": expected_size, "mode": 0o500,
            "raw_sha256": expected_sha256,
            "seals": runner.REQUIRED_SEALS, "elf": True,
        }
        return descriptor, {
            "source": {
                "authority_relative_path": authority_relative_path,
                "size": expected_size, "raw_sha256": expected_sha256,
            },
            "snapshot": snapshot,
        }

    return capture


@contextlib.contextmanager
def patched_live_acquisition(
    fixture: dict, *, enforce_output_lane_binding: bool = False,
):
    def current_monotonic_ns():
        bridge = fixture["pure"]["bridge"]
        return (bridge["bridge_finished_monotonic_ns"] +
                bridge["bridge_deadline_monotonic_ns"]) // 2
    original_controller_validator = \
        runner._validate_controller_bindings_current.__wrapped__
    (unused_armer_modules, unused_controller_modules,
     executable_references, class_member_references,
     unused_function_states) = \
        runner._frozen_controller_execution_graph()
    (live_executable_references, live_class_member_references,
     unused_live_function_states) = runner._frozen_live_execution_graph()
    public_callable_references = runner._frozen_live_public_surface()
    original_open_output_lane = \
        runner.prereg.open_output_lane_identity_for_preregistration
    original_validate_output_lane_descriptor = \
        runner.prereg.validate_output_lane_descriptor_identity_for_preregistration
    original_validate_preregistration = runner.prereg.validate_preregistration
    bound_test_lane: list[Path | None] = [None]

    def registration_bound_to_test_lane(registration, lane):
        adjusted = copy.deepcopy(registration)
        adjusted["output_lane"] = lane_binding_for_test(lane)
        return adjusted

    def open_test_output_lane(output_lane, registration, repo_root=runner.REPO_ROOT):
        bound_test_lane[0] = output_lane
        lane_binding = lane_binding_for_test(output_lane)
        arguments = fixture["pure_fixture_args"]
        retargeted = pure_fixtures.fixture(
            arguments["candidate_sha256"],
            candidate_authority_record_sha256=
                arguments["candidate_authority_record_sha256"],
            candidate_authority_ledger_sha256=
                arguments["candidate_authority_ledger_sha256"],
            controller_bindings=arguments["controller_bindings"],
            output_lane_binding=lane_binding,
            expected_frozen_pair=arguments["expected_frozen_pair"])
        fixture["pure"] = retargeted
        fixture["candidate_authority"]["source_authority"] = \
            retargeted["source_authority"]
        fixture["qualification"] = {
            "host_instance": retargeted["host"],
            "qualification_evidence": retargeted["evidence"],
            "qualification_binding": retargeted["qualification"],
        }
        registration["output_lane"] = copy.deepcopy(lane_binding)
        return original_open_output_lane(
            output_lane,
            registration_bound_to_test_lane(registration, output_lane),
            repo_root)

    def validate_test_preregistration(value, *args, **kwargs):
        validated = original_validate_preregistration(value, *args, **kwargs)
        if bound_test_lane[0] is None:
            return validated
        adjusted = copy.deepcopy(validated)
        adjusted["output_lane"] = lane_binding_for_test(bound_test_lane[0])
        return adjusted

    def validate_test_output_lane_descriptor(
        output_lane, lane_fd, registration, repo_root=runner.REPO_ROOT,
    ):
        return original_validate_output_lane_descriptor(
            output_lane, lane_fd,
            registration_bound_to_test_lane(registration, output_lane),
            repo_root)

    def validate_with_intentional_test_mocks(registration):
        missing = object()
        restored = []
        for owner, name, expected in (
                *public_callable_references, *executable_references,
                *class_member_references, *live_executable_references,
                *live_class_member_references):
            current = vars(owner).get(name, missing)
            if current is expected:
                continue
            restored.append((owner, name, current))
            setattr(owner, name, expected)
        try:
            # The plan-only suite deliberately loads some of the same source
            # files under ordinary aliases.  Production validation is wrapped
            # in this scope, which shadows those ambient aliases with the
            # retained exact-source modules.  The live-acquisition fixture
            # unwraps the public API so it can install child/host canaries, so
            # reproduce that production alias isolation around its raw
            # validator while the intentional callable mocks are restored.
            with runner._ExactImportScope():
                return original_controller_validator(registration)
        finally:
            for owner, name, current in reversed(restored):
                if current is missing:
                    delattr(owner, name)
                else:
                    setattr(owner, name, current)

    with contextlib.ExitStack() as stack:
        lock = FakeCanonicalLock(runner.AUTHORITATIVE_LOCK_PATH)
        runtime_guard = FakeRuntimeGuard()
        stack.enter_context(mock.patch.object(
            runner.baseline_acquire, "CanonicalFileLock",
            return_value=lock))
        stack.enter_context(mock.patch.object(
            runner, "bind_candidate_authority_lane",
            return_value=fixture["candidate_authority"]))
        stack.enter_context(mock.patch.object(
            runner, "bind_exact_main_path_variant",
            return_value=fixture["exact_main"]))
        stack.enter_context(mock.patch.object(
            runner, "capture_sealed_executable_bytes",
            side_effect=fake_capture_factory()))
        qualification = stack.enter_context(mock.patch.object(
            runner, "_acquire_track_a_qualification",
            side_effect=lambda *args, **kwargs: fixture["qualification"]))
        stack.enter_context(mock.patch.object(
            runner, "_capture_live_host_instance",
            side_effect=lambda *args, **kwargs: fixture["pure"]["host"]))
        stack.enter_context(mock.patch.object(
            runner, "_acquire_runtime_guard", return_value=runtime_guard))
        stack.enter_context(mock.patch.object(
            runner, "revalidate_sealed_descriptor",
            return_value={}))
        stack.enter_context(mock.patch.object(
            runner.time, "monotonic_ns",
            side_effect=lambda: current_monotonic_ns()))
        popen = stack.enter_context(mock.patch.object(
            runner.subprocess, "Popen",
            side_effect=AssertionError(
                "child process ran before durable ARMED")))
        stack.enter_context(mock.patch.object(
            runner, "_validate_controller_bindings_current",
            side_effect=validate_with_intentional_test_mocks))
        stack.enter_context(mock.patch.object(
            runner, "acquire_and_arm", runner.acquire_and_arm.__wrapped__))
        stack.enter_context(mock.patch.object(
            runner.ArmedCampaign, "run_all",
            runner.ArmedCampaign.run_all.__wrapped__))
        if not enforce_output_lane_binding:
            stack.enter_context(mock.patch.object(
                runner.prereg, "validate_preregistration",
                side_effect=validate_test_preregistration))
            stack.enter_context(mock.patch.object(
                runner.prereg,
                "open_output_lane_identity_for_preregistration",
                side_effect=open_test_output_lane))
            stack.enter_context(mock.patch.object(
                runner.prereg,
                "validate_output_lane_descriptor_identity_for_preregistration",
                side_effect=validate_test_output_lane_descriptor))
        yield {
            "lock": lock, "runtime_guard": runtime_guard,
            "qualification": qualification, "popen": popen,
        }


class FakeSealedTree:
    def __init__(self, files: dict[str, bytes]):
        self.payload = files
        self.files = tuple(sorted(files))
        self.reverified = False

    def __enter__(self):
        return self

    def __exit__(self, unused_kind, unused_value, unused_traceback):
        return None

    def read_file(self, name: str, maximum_bytes: int):
        value = self.payload[name]
        if len(value) > maximum_bytes:
            raise ValueError("oversized fake retained file")
        return value

    def reverify(self):
        self.reverified = True


def candidate_authority_fixture() -> tuple[dict, dict[str, bytes]]:
    candidate = synthetic_elf(16_384, 0x71)
    candidate_sha = hashlib.sha256(candidate).hexdigest()
    commit = "1" * 40
    tree = "2" * 40
    git_data = contract.canonical_json_bytes({})
    archive = b"unit source archive"
    provenance = {
        "tracked_source_manifest": {
            "git": {"commit": commit, "tree": tree, "dirty": False}},
        "executable": {"sha256": candidate_sha, "size": len(candidate)},
    }
    provenance_data = contract.canonical_json_bytes(provenance)
    core = {"schema": "unit-reproducible-core/v1"}
    core_data = contract.canonical_json_bytes(core)
    proof_data = contract.canonical_json_bytes({})
    header = b"#define UNIT_ATTESTATION 1\n"
    launchers = [{
        "path": path, "sha256": hashlib.sha256(path.encode()).hexdigest(),
        "size": 1, "mode": 0o755, "role": "launcher",
    } for path in ("/usr/bin/prlimit", "/usr/bin/taskset")]
    runtime = {
        "schema": "leopard2-k65-gen3-runtime-closure/v1",
        "dependencies": [], "launchers": launchers,
    }
    runtime_data = contract.canonical_json_bytes(runtime)
    payload = {
        "artifacts/bench_leopard2": candidate,
        "build/benchmark-source-attestation.h": header,
        "build/build-provenance.json": provenance_data,
        "build/reproducible-build-core.json": core_data,
        "build/reproducible-build-proof.json": proof_data,
        "runtime/runtime-closure.json": runtime_data,
        "source/candidate-source.tar": archive,
        "source/git-capture.json": git_data,
    }
    build_provenance_sha = hashlib.sha256(provenance_data).hexdigest()
    core_sha = hashlib.sha256(core_data).hexdigest()
    controllers = runner.prereg.current_controller_bindings()
    preliminary = pure_fixtures.preregistration(
        candidate_sha,
        candidate_size=len(candidate),
        candidate_build_provenance_sha256=build_provenance_sha,
        candidate_reproducible_build_core_sha256=core_sha,
        controller_bindings=controllers)
    record = {
        "schema": runner.CANDIDATE_AUTHORITY_SCHEMA,
        "generation": 3, "status": "promoted_authority",
        "source": {
            "commit": commit, "tree": tree,
            "git_capture_sha256": hashlib.sha256(git_data).hexdigest(),
            "source_archive_sha256": hashlib.sha256(archive).hexdigest(),
        },
        "candidate": {
            "relative_path": "artifacts/bench_leopard2",
            "sha256": candidate_sha, "size": len(candidate),
        },
        "build_closure": {
            "build_provenance_sha256": build_provenance_sha,
            "reproducible_build_core_sha256": core_sha,
            "reproducible_build_proof_sha256":
                hashlib.sha256(proof_data).hexdigest(),
            "attestation_header_sha256": hashlib.sha256(header).hexdigest(),
            "runtime_closure_sha256": hashlib.sha256(runtime_data).hexdigest(),
        },
        "controller_bindings_sha256": contract.canonical_sha256(
            preliminary["controller_bindings"]),
        "inventory": [{
            "relative_path": path,
            "sha256": hashlib.sha256(payload[path]).hexdigest(),
            "size": len(payload[path]),
        } for path in runner._CANDIDATE_PAYLOAD_PATHS],
    }
    record_data = contract.canonical_json_bytes(record)
    ledger_sha = "7d" * 32
    registration = pure_fixtures.preregistration(
        candidate_sha,
        candidate_size=len(candidate),
        candidate_build_provenance_sha256=build_provenance_sha,
        candidate_reproducible_build_core_sha256=core_sha,
        candidate_authority_record_sha256=
            hashlib.sha256(record_data).hexdigest(),
        candidate_authority_ledger_sha256=ledger_sha,
        controller_bindings=controllers)
    payload["candidate-authority.json"] = record_data
    payload["SHA256SUMS"] = b"unit ledger\n"
    payload["TREE-METADATA.json"] = b"{}"
    return registration, payload


class K65Generation3LiveArmerTests(unittest.TestCase):
    def assertRejected(self, function, pattern: str | None = None) -> None:
        with self.assertRaises(ValueError) as raised:
            function()
        if pattern is not None:
            self.assertIn(pattern, str(raised.exception))

    def test_production_api_has_no_authority_injection_surface(self) -> None:
        self.assertEqual(runner.__all__, (
            "ArmedCampaign", "ArmingError", "acquire_and_arm"))
        signature = inspect.signature(runner.acquire_and_arm)
        self.assertEqual(list(signature.parameters), [
            "lane", "preregistration_bytes", "candidate_authority_lane",
            "exact_main_authority_lane",
        ])
        self.assertTrue(all(parameter.kind is inspect.Parameter.KEYWORD_ONLY
                            for parameter in signature.parameters.values()))
        for forbidden in (
                "reader", "evidence", "frozen_pair", "candidate", "launcher",
                "fault_hook", "record_factory", "subprocess_kwargs"):
            self.assertNotIn(forbidden, signature.parameters)
        self.assertFalse(hasattr(runner, "arm_from_fresh_qualification"))
        self.assertFalse(hasattr(runner, "publish_armed_record"))
        self.assertFalse(hasattr(runner.ArmedCampaign, "materialize_child"))
        self.assertFalse(hasattr(runner.ArmedCampaign, "launch_child"))

    def test_exact_controller_graph_binds_actual_transitive_modules(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        bindings = {
            (runner.REPO_ROOT / item["path"]).resolve(): item["sha256"]
            for item in registration["controller_bindings"]
        }
        required = {
            (runner.REPO_ROOT / relative).resolve()
            for relative in runner._REQUIRED_CONTROLLER_RELATIVE_PATHS
        }
        self.assertEqual(set(bindings), required)
        for relative in (
                "experiments/leopard2/gf8_high_encode/"
                "run_k66r16_b64_tail_abba.py",
                "experiments/leopard2/gf8_high_encode/"
                "run_k9r5_b1024_terminal_abba.py",
                "experiments/leopard2/gf8_high_encode/"
                "run_k5r5_b64_terminal_abba.py",
                "experiments/leopard2/gf8_high_encode/run_t8_two_block_abba.py",
                "experiments/leopard2/main_compare/run_abba.py",
                "experiments/leopard2/decoder_dispatch/"
                "balanced_evidence_common.py"):
            self.assertIn((runner.REPO_ROOT / relative).resolve(), bindings)
        for path, digest in bindings.items():
            self.assertEqual(
                runner._LOADED_CONTROLLER_SOURCE_SHA256[path], digest)
        self.assertIs(runner.git_capture._build_provenance,
                      runner.build_provenance)
        before = len(runner._EXACT_CONTROLLER_LOADS)
        for unused_index in range(3):
            runner.prereg.validate_preregistration(
                registration, verify_files=True)
            runner._validate_controller_bindings_current(registration)
        self.assertEqual(len(runner._EXACT_CONTROLLER_LOADS), before)

    def test_controller_module_references_are_exact_identity_joins(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        ledger_ids = {
            id(module) for _, module, _, _ in runner._EXACT_CONTROLLER_LOADS
        }
        reviewed_edges = [
            (owner, name, expected)
            for owner, name, expected in
                runner._CONTROLLER_MODULE_REFERENCE_IDENTITIES
            if id(expected) in ledger_ids
        ]
        self.assertGreaterEqual(len(reviewed_edges), 20)
        for index, (owner, name, expected) in enumerate(reviewed_edges):
            rogue = types.ModuleType(f"rogue_controller_dependency_{index}")
            self.assertFalse(hasattr(rogue, "__file__"))
            setattr(owner, name, rogue)
            try:
                with self.assertRaisesRegex(
                        RuntimeError, "controller module reference changed"):
                    runner._validate_controller_bindings_current(registration)
            finally:
                setattr(owner, name, expected)
        rogue = types.ModuleType("rogue_live_armer_dependency")
        expected_result_contract = runner.result_contract
        runner.result_contract = rogue
        try:
            with self.assertRaisesRegex(
                    RuntimeError,
                    "live armer module reference changed: result_contract"):
                runner._validate_controller_bindings_current(registration)
        finally:
            runner.result_contract = expected_result_contract

    def test_exact_controller_ledgers_freeze_before_api_use(self):
        self.assertIsInstance(runner._EXACT_CONTROLLER_LOADS, tuple)
        self.assertIsInstance(
            runner._EXACT_CONTROLLER_MODULES, types.MappingProxyType)
        fake_path = (
            runner.HERE / "run_k65r65_b64_packed_terminal_abba.py").resolve()
        fake_digest = runner._LOADED_CONTROLLER_SOURCE_SHA256[fake_path]
        fake = types.ModuleType("forged_unexecuted_controller")
        fake.__file__ = str(fake_path)
        fake.__exact_source_sha256__ = fake_digest
        fake_entry = (
            "forged_unexecuted_controller", fake, fake_path, fake_digest)
        with self.assertRaises(AttributeError):
            runner._EXACT_CONTROLLER_LOADS.append(fake_entry)
        with self.assertRaises(TypeError):
            runner._EXACT_CONTROLLER_MODULES[fake_entry[0]] = (
                fake, fake_path, fake_digest)
        with self.assertRaisesRegex(RuntimeError, "ledger is frozen"):
            runner._ExactSourceLoader(
                fake_entry[0], fake_path).exec_module(fake)
        with self.assertRaisesRegex(RuntimeError, "ledger is frozen"):
            runner._load(fake_entry[0], fake_path)

        frozen_loads = runner._EXACT_CONTROLLER_LOADS
        frozen_modules = runner._EXACT_CONTROLLER_MODULES
        forged_loads = (*frozen_loads, fake_entry)
        forged_modules = types.MappingProxyType({
            **dict(frozen_modules),
            fake_entry[0]: (fake, fake_path, fake_digest),
        })
        runner._EXACT_CONTROLLER_LOADS = forged_loads
        runner._EXACT_CONTROLLER_MODULES = forged_modules
        try:
            with self.assertRaisesRegex(
                    RuntimeError, "ledger binding changed"):
                with runner._ExactImportScope():
                    self.fail("rebound controller ledger was accepted")
        finally:
            runner._EXACT_CONTROLLER_LOADS = frozen_loads
            runner._EXACT_CONTROLLER_MODULES = frozen_modules

        alias = fake_entry[0]
        old_alias = runner.sys.modules.get(alias, None)
        self.assertIsNone(old_alias)
        with self.assertRaisesRegex(
                RuntimeError, "new controller alias appeared in scope"):
            with runner._ExactImportScope():
                runner.sys.modules[alias] = fake
        self.assertNotIn(alias, runner.sys.modules)
        fixture = live_fixture()
        runner._validate_controller_bindings_current(
            fixture["pure"]["preregistration"])

        expected_callable = runner.prereg.validate_preregistration

        class ForgedCallable:
            def __call__(self, *unused_arguments, **unused_keywords):
                return fixture["pure"]["preregistration"]

        runner.prereg.validate_preregistration = ForgedCallable()
        try:
            self.assertRejected(
                lambda: runner._validate_controller_bindings_current(
                    fixture["pure"]["preregistration"]),
                "controller executable binding changed")
        finally:
            runner.prereg.validate_preregistration = expected_callable

    def test_exact_scope_rebuilds_hostile_module_registry_keys(self):
        class FlakyKey(str):
            fail_hash = False

            def __hash__(self):
                if self.fail_hash:
                    raise RuntimeError("hostile key hash failure")
                return super().__hash__()

        key = FlakyKey("mid_scope_hostile_controller_key")
        original_registry = runner.sys.modules
        original_snapshot = dict(runner.sys.modules)
        with self.assertRaisesRegex(RuntimeError, "non-string key"):
            with runner._ExactImportScope():
                runner.sys.modules[key] = types.ModuleType(str(key))
                key.fail_hash = True
        self.assertIs(runner.sys.modules, original_registry)
        self.assertEqual(len(runner.sys.modules), len(original_snapshot))
        self.assertTrue(all(type(alias) is str for alias in runner.sys.modules))
        for alias, expected in original_snapshot.items():
            self.assertIs(runner.sys.modules[alias], expected)

        acquired: list[bool] = []

        def probe_lock():
            locked = runner._EXACT_IMPORT_LOCK.acquire(timeout=1)
            acquired.append(locked)
            if locked:
                runner._EXACT_IMPORT_LOCK.release()

        worker = threading.Thread(target=probe_lock)
        worker.start()
        worker.join(timeout=2)
        self.assertFalse(worker.is_alive())
        self.assertEqual(acquired, [True])

    def test_exact_scope_shadows_and_restores_ambient_alias_pollution(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        alias = "t8_two_block_main_compare_support"
        expected = runner._EXACT_CONTROLLER_MODULES[alias][0]
        ambient_alias = "ambient_k65_v2_plan_comparison"
        ambient = types.ModuleType(ambient_alias)
        required_path = (
            runner.HERE / "run_k65r65_b64_packed_terminal_abba.py")
        path_calls = []

        class CustomPathLike:
            def __fspath__(self):
                path_calls.append("__fspath__")
                return str(required_path)

        ambient.__file__ = CustomPathLike()
        self.assertIsNone(runner._controller_module_path(ambient))
        self.assertEqual(path_calls, [])
        missing = object()
        original = runner.sys.modules.get(alias, missing)
        original_ambient = runner.sys.modules.get(ambient_alias, missing)
        foreign = types.ModuleType("ambient_foreign_t8_support")
        try:
            runner.sys.modules[alias] = foreign
            runner.sys.modules[ambient_alias] = ambient
            runner._validate_controller_bindings_current(registration)
            self.assertIs(runner.sys.modules[alias], foreign)
            self.assertIs(runner.sys.modules[ambient_alias], ambient)
            with self.assertRaisesRegex(
                    RuntimeError, "exact controller alias changed in scope"):
                with runner._ExactImportScope():
                    self.assertIs(runner.sys.modules[alias], expected)
                    self.assertNotIn(ambient_alias, runner.sys.modules)
                    runner.sys.modules[alias] = types.ModuleType(
                        "mid_scope_foreign_t8_support")
            self.assertIs(runner.sys.modules[alias], foreign)
            self.assertIs(runner.sys.modules[ambient_alias], ambient)

            del runner.sys.modules[alias]
            with runner._ExactImportScope():
                self.assertIs(runner.sys.modules[alias], expected)
                self.assertNotIn(ambient_alias, runner.sys.modules)
            self.assertNotIn(alias, runner.sys.modules)
            self.assertIs(runner.sys.modules[ambient_alias], ambient)
        finally:
            if original is missing:
                runner.sys.modules.pop(alias, None)
            else:
                runner.sys.modules[alias] = original
            if original_ambient is missing:
                runner.sys.modules.pop(ambient_alias, None)
            else:
                runner.sys.modules[ambient_alias] = original_ambient

    def test_exact_scope_rejects_new_alias_and_self_identity_spoof(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        new_alias = "mid_scope_new_required_controller"
        spoof_alias = "ambient_live_armer_identity_spoof"
        missing = object()
        old_new = runner.sys.modules.get(new_alias, missing)
        old_spoof = runner.sys.modules.get(spoof_alias, missing)
        required_path = str(
            runner.HERE / "run_k65r65_b64_packed_terminal_abba.py")

        class SelfIdentitySpoof:
            __file__ = str(RUNNER_PATH)
            __exact_source_sha256__ = runner.__exact_source_sha256__

            @property
            def __dict__(self):
                return vars(runner)

        spoof = SelfIdentitySpoof()
        try:
            runner.sys.modules[spoof_alias] = spoof
            with runner._ExactImportScope():
                self.assertNotIn(spoof_alias, runner.sys.modules)
            self.assertIs(runner.sys.modules[spoof_alias], spoof)
            runner._validate_controller_bindings_current(registration)
            self.assertIs(runner.sys.modules[spoof_alias], spoof)

            new_module = types.ModuleType(new_alias)
            new_module.__file__ = required_path
            with self.assertRaisesRegex(
                    RuntimeError, "new controller alias appeared in scope"):
                with runner._ExactImportScope():
                    runner.sys.modules[new_alias] = new_module
            self.assertNotIn(new_alias, runner.sys.modules)
            runner._validate_controller_bindings_current(registration)
        finally:
            if old_new is missing:
                runner.sys.modules.pop(new_alias, None)
            else:
                runner.sys.modules[new_alias] = old_new
            if old_spoof is missing:
                runner.sys.modules.pop(spoof_alias, None)
            else:
                runner.sys.modules[spoof_alias] = old_spoof

    def test_exact_scope_uses_identity_for_finder_and_loader_hook(self):
        class HostileEquality:
            def __eq__(self, unused_other):
                raise AssertionError("equality must not control import scope")

            def __call__(self, *unused_arguments, **unused_keywords):
                raise AssertionError("hostile loader hook was called")

        hostile = HostileEquality()
        previous_hook = runner.importlib.util.spec_from_file_location
        runner.importlib.util.spec_from_file_location = hostile
        try:
            with self.assertRaisesRegex(RuntimeError, "source loader"):
                with runner._ExactImportScope():
                    self.fail("hostile source loader was accepted")
        finally:
            runner.importlib.util.spec_from_file_location = previous_hook

        previous_module_hook = runner.importlib.util.module_from_spec
        runner.importlib.util.module_from_spec = hostile
        try:
            with self.assertRaisesRegex(RuntimeError, "source loader"):
                with runner._ExactImportScope():
                    self.fail("hostile module constructor was accepted")
        finally:
            runner.importlib.util.module_from_spec = previous_module_hook

        acquired: list[bool] = []

        def probe_lock():
            locked = runner._EXACT_IMPORT_LOCK.acquire(timeout=1)
            acquired.append(locked)
            if locked:
                runner._EXACT_IMPORT_LOCK.release()

        worker = threading.Thread(target=probe_lock)
        worker.start()
        worker.join(timeout=2)
        self.assertFalse(worker.is_alive())
        self.assertEqual(acquired, [True])

        original_meta_path = tuple(runner.sys.meta_path)
        first_equal = HostileEquality()
        second_equal = HostileEquality()
        runner.sys.meta_path.insert(0, first_equal)
        polluted_snapshot = tuple(runner.sys.meta_path)
        try:
            with self.assertRaisesRegex(
                    RuntimeError, "finder identity or order changed"):
                with runner._ExactImportScope():
                    self.assertIs(
                        runner.sys.meta_path[0],
                        runner._EXACT_REPOSITORY_FINDER)
                    self.assertTrue(any(
                        entry is first_equal for entry in runner.sys.meta_path))
                    runner.sys.meta_path.insert(0, second_equal)
            self.assertEqual(
                len(runner.sys.meta_path), len(polluted_snapshot))
            self.assertTrue(all(
                actual is expected
                for actual, expected in zip(
                    runner.sys.meta_path, polluted_snapshot)))
        finally:
            runner.sys.meta_path[:] = original_meta_path

    def test_exact_finder_discards_hostile_pathfinder_specification(self):
        required = (
            runner.HERE / "k65_gen3_execution_contract.py").resolve()

        class WrappedLoader:
            def __init__(self, delegate):
                self.delegate = delegate

            def exec_module(self, module):
                self.delegate.exec_module(module)
                module.validate_armed_record = lambda unused_value: None

        class HostileSpecification(runner._MODULE_SPEC_TYPE):
            def __setattr__(self, name, value):
                if (name == "loader" and
                        type(value) is runner._ExactSourceLoader):
                    value = WrappedLoader(value)
                super().__setattr__(name, value)

        hostile = HostileSpecification(
            "hostile_pathfinder_result", None, origin=str(required))
        with mock.patch.object(
                runner.importlib.machinery.PathFinder, "find_spec",
                return_value=hostile):
            specification = runner._EXACT_REPOSITORY_FINDER.find_spec(
                "discarded_hostile_pathfinder_spec")
        self.assertIs(type(specification), runner._MODULE_SPEC_TYPE)
        self.assertIs(type(specification.loader), runner._ExactSourceLoader)
        self.assertEqual(specification.origin, str(required))

    def test_exact_scope_exit_failure_revokes_transferred_campaign(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        alias = "t8_two_block_main_compare_support"
        descriptors: list[int] = []

        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                baseline = len(os.listdir("/proc/self/fd"))
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=
                        contract.canonical_json_bytes(registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            descriptors.extend(getattr(campaign, name) for name in (
                "_main_fd", "_candidate_fd", "_manifest_fd", "_bundle_fd",
                "_armed_fd", "_launch_marker_fd", "_journal_fd",
                "_campaign_marker_fd", "_attempt_fd", "_attempts_fd",
                "_lane_binding_fd", "_lane_lock_fd", "_lane_fd"))

            @runner._with_exact_source_imports
            def transfer_then_mutate_alias():
                runner.sys.modules[alias] = types.ModuleType(
                    "post_transfer_alias_mutation")
                return campaign

            with self.assertRaisesRegex(
                    RuntimeError, "exact controller alias changed in scope"):
                transfer_then_mutate_alias()
            self.assertTrue(campaign._closed)
            self.assertFalse(controls["lock"].held)
            self.assertTrue(controls["runtime_guard"].closed)
            self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
            self.assertEqual(len(descriptors), 13)
            for descriptor in descriptors:
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    def test_ordinary_self_import_is_rejected_and_seeded_aliases_are_purged(
        self,
    ):
        with self.assertRaisesRegex(RuntimeError, "exact-source bootstrap"):
            load_module("k65_gen3_forbidden_ordinary_import", RUNNER_PATH)
        script = f"""
import importlib.util
from pathlib import Path
import runpy
import sys

def ordinary(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

root = Path({str(runner.REPO_ROOT)!r})
pair_name = 'leopard2_pair_qualification_contract_for_k65_gen3_prereg'
old_pair = ordinary(pair_name, root / 'experiments/leopard2/main_compare/pair_qualification_contract.py')
old_pair.__file__ = Path(old_pair.__file__)
old_build = ordinary('leopard2_build_provenance', root / 'tools/leopard2_build_provenance.py')
old_build.__file__ = bytes(old_build.__file__, 'utf-8')
bootstrap = runpy.run_path(
    {str(BOOTSTRAP_PATH)!r}, run_name='seeded_alias_direct_bootstrap')
live = bootstrap['load_live_armer']('seeded_alias_live_armer')
assert sys.modules[pair_name] is not old_pair
assert sys.modules['leopard2_build_provenance'] is not old_build
assert live.git_capture._build_provenance is live.build_provenance
"""
        completed = runner.subprocess.run(
            ["/usr/bin/python3", "-c", script], cwd=str(HERE),
            stdout=runner.subprocess.PIPE, stderr=runner.subprocess.PIPE,
            check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())

    def test_exact_scope_never_invokes_ambient_module_descriptors(self):
        saved_acquire = runner.acquire_and_arm
        spoofed = lambda **unused_keywords: {"spoofed": True}
        getter_calls = []

        class PassiveModule(types.ModuleType):
            def __getattribute__(self, name):
                if name == "__file__":
                    getter_calls.append(name)
                    vars(runner)["acquire_and_arm"] = spoofed
                    return str(RUNNER_PATH)
                return super().__getattribute__(name)

        alias = "unrelated_passive_module_descriptor"
        passive = PassiveModule(alias)
        old_alias = runner.sys.modules.get(alias, None)
        try:
            runner.sys.modules[alias] = passive
            with runner._ExactImportScope():
                self.assertNotIn(alias, runner.sys.modules)
            self.assertEqual(getter_calls, [])
            self.assertIs(runner.acquire_and_arm, saved_acquire)

            runner.acquire_and_arm = spoofed
            with self.assertRaisesRegex(
                    RuntimeError, "live public callable identity changed"):
                with runner._ExactImportScope():
                    self.fail("changed public callable entered exact scope")
            self.assertIs(runner.acquire_and_arm, saved_acquire)
        finally:
            runner.acquire_and_arm = saved_acquire
            if old_alias is None:
                runner.sys.modules.pop(alias, None)
            else:
                runner.sys.modules[alias] = old_alias

    def test_bootstrap_rejects_compile_equivalent_stale_bytes_and_handoff_swap(
        self,
    ):
        source = BOOTSTRAP_PATH.read_bytes()
        old_comment = b"# Exact-byte self-digest trust anchor."
        new_comment = b"# exact-byte self-digest trust anchor."
        self.assertEqual(len(old_comment), len(new_comment))
        altered = source.replace(old_comment, new_comment, 1)
        self.assertNotEqual(altered, source)
        self.assertEqual(
            compile(source, str(BOOTSTRAP_PATH), "exec"),
            compile(altered, str(BOOTSTRAP_PATH), "exec"))
        namespace = load_live_armer.__globals__
        with mock.patch.dict(namespace, {
                "_read_exact_source_bytes": lambda unused_path: (
                    altered, hashlib.sha256(altered).hexdigest())}):
            with self.assertRaisesRegex(RuntimeError, "exact bytes changed"):
                namespace["_validate_direct_source_execution"]()

        stale_digest = "0" * 64
        with mock.patch.dict(namespace, {
                "_validate_direct_source_execution": lambda: stale_digest,
                "_DIRECT_SOURCE_BOOTSTRAP_SHA256": stale_digest}):
            with self.assertRaisesRegex(RuntimeError, "bootstrap changed"):
                load_live_armer("k65_gen3_bootstrap_handoff_swap")

    def test_bootstrap_rejects_passive_import_helper_pollution(self):
        for attribute in ("spec_from_file_location", "module_from_spec"):
            script = f"""
import importlib.util
from pathlib import Path
import runpy

original = getattr(importlib.util, {attribute!r})
def delegating_hook(*arguments, **keywords):
    return original(*arguments, **keywords)
setattr(importlib.util, {attribute!r}, delegating_hook)
bootstrap = runpy.run_path(
    {str(BOOTSTRAP_PATH)!r}, run_name='polluted_import_helper_bootstrap')
try:
    bootstrap['load_live_armer']('polluted_import_helper_live_armer')
except RuntimeError as error:
    assert 'polluted before bootstrap' in str(error), repr(error)
else:
    raise AssertionError('passive import helper pollution was accepted')
"""
            completed = runner.subprocess.run(
                ["/usr/bin/python3", "-c", script], cwd=str(HERE),
                stdout=runner.subprocess.PIPE, stderr=runner.subprocess.PIPE,
                check=False)
            self.assertEqual(
                completed.returncode, 0,
                f"{attribute}: {completed.stderr.decode()}")

    def test_bootstrap_rejects_passive_sys_module_subclass_alias(self):
        script = f"""
import runpy
import sys as real_sys
import types

class ProxyModule(types.ModuleType):
    def __getattribute__(self, name):
        if name == 'acquire_and_arm':
            return lambda **unused_keywords: {{'spoofed': True}}
        return super().__getattribute__(name)

proxy_sys = ProxyModule('sys')
proxy_sys.__dict__.update(real_sys.__dict__)
real_sys.modules['sys'] = proxy_sys
try:
    bootstrap = runpy.run_path(
        {str(BOOTSTRAP_PATH)!r}, run_name='polluted_sys_alias_bootstrap')
    bootstrap['load_live_armer']('polluted_sys_alias_live_armer')
except RuntimeError as error:
    assert 'module alias was replaced' in str(error), repr(error)
else:
    raise AssertionError('passive sys module subclass alias was accepted')
finally:
    real_sys.modules['sys'] = real_sys
"""
        completed = runner.subprocess.run(
            ["/usr/bin/python3", "-c", script], cwd=str(HERE),
            stdout=runner.subprocess.PIPE, stderr=runner.subprocess.PIPE,
            check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())

    def test_bootstrap_rejects_hostile_module_registry_key_without_hashing(self):
        script = f"""
from pathlib import Path
import runpy
import sys
import types

callbacks = []
class HostileKey(str):
    armed = False
    def __hash__(self):
        if self.armed:
            callbacks.append('hash')
        return super().__hash__()

key = HostileKey('passive_prebootstrap_controller_alias')
module = types.ModuleType(str(key))
module.__file__ = str(
    Path({str(runner.REPO_ROOT)!r}) /
    'experiments/leopard2/gf8_high_encode/k65_gen3_preregistration.py')
sys.modules[key] = module
key.armed = True
try:
    bootstrap = runpy.run_path(
        {str(BOOTSTRAP_PATH)!r}, run_name='hostile_key_bootstrap')
    bootstrap['load_live_armer']('hostile_key_live_armer')
except RuntimeError as error:
    assert 'non-string key' in str(error) or 'unsafe keys' in str(error), repr(error)
else:
    raise AssertionError('hostile module-registry key was accepted')
assert callbacks == [], callbacks
"""
        completed = runner.subprocess.run(
            ["/usr/bin/python3", "-c", script], cwd=str(HERE),
            stdout=runner.subprocess.PIPE, stderr=runner.subprocess.PIPE,
            check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())

    def test_controller_load_mismatch_rejects_before_lane_or_qualification(self):
        path = (runner.REPO_ROOT / "tools/leopard2_build_provenance.py").resolve()
        original = runner._LOADED_CONTROLLER_SOURCE_SHA256[path]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            fixture = live_fixture(output_lane=lane)
            registration = fixture["pure"]["preregistration"]
            before = sorted(path.name for path in lane.iterdir())
            try:
                runner._LOADED_CONTROLLER_SOURCE_SHA256[path] = "0" * 64
                with patched_live_acquisition(
                        fixture,
                        enforce_output_lane_binding=True) as controls:
                    self.assertRejected(lambda: runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=contract.canonical_json_bytes(
                            registration),
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path("/exact-main-authority")),
                        "loaded controller bytes differ")
                    self.assertEqual(controls["qualification"].call_count, 0)
            finally:
                runner._LOADED_CONTROLLER_SOURCE_SHA256[path] = original
            self.assertEqual(sorted(path.name for path in lane.iterdir()), before)

    def test_acquire_and_run_revalidation_never_uses_sourcefileloader(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        before = len(runner._EXACT_CONTROLLER_LOADS)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture), mock.patch.object(
                    runner.importlib.machinery.SourceFileLoader, "exec_module",
                    side_effect=AssertionError("ordinary loader used")):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                try:
                    campaign._revalidate_authority()
                    campaign._revalidate_authority()
                finally:
                    campaign.close()
        self.assertEqual(len(runner._EXACT_CONTROLLER_LOADS), before)

    def test_overlap_rejection_leaves_lane_inventory_untouched(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            sentinel = lane / "sentinel"
            sentinel.write_bytes(b"unchanged")
            before = sorted(path.name for path in lane.iterdir())
            with patched_live_acquisition(fixture) as controls, \
                    mock.patch.object(
                        runner.prereg,
                        "open_output_lane_identity_for_preregistration",
                        side_effect=runner.ArmingError(
                            "injected authority-lane overlap")):
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "overlap")
                self.assertEqual(controls["qualification"].call_count, 0)
            self.assertEqual(sorted(path.name for path in lane.iterdir()), before)
            self.assertEqual(sentinel.read_bytes(), b"unchanged")
            self.assertFalse((lane / ".arming.lock").exists())

    def test_preregistration_rejects_a_second_output_lane_before_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            ratified_lane = private_lane(directory, "ratified-lane")
            other_lane = private_lane(directory, "other-lane")
            fixture = live_fixture(output_lane=ratified_lane)
            preregistration_bytes = contract.canonical_json_bytes(
                fixture["pure"]["preregistration"])
            with patched_live_acquisition(
                    fixture, enforce_output_lane_binding=True):
                campaign = runner.acquire_and_arm(
                    lane=ratified_lane,
                    preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            campaign.close()

            with patched_live_acquisition(
                    fixture, enforce_output_lane_binding=True) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=other_lane,
                    preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "ratified canonical path")
                self.assertFalse(controls["lock"].held)
                self.assertEqual(controls["qualification"].call_count, 0)
            self.assertEqual(list(other_lane.iterdir()), [])
            self.assertEqual(other_lane.stat().st_mode & 0o777, 0o700)

    def test_recreated_ratified_path_cannot_reset_campaign_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory, "ratified-lane")
            fixture = live_fixture(output_lane=lane)
            registration = fixture["pure"]["preregistration"]
            preregistration_bytes = contract.canonical_json_bytes(
                registration)
            with patched_live_acquisition(
                    fixture, enforce_output_lane_binding=True):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            self.assertEqual(campaign.record["evidence_attempt"], 1)
            campaign.close()

            for child in sorted(
                    (path for path in lane.rglob("*") if path.is_dir()),
                    key=lambda path: len(path.parts), reverse=True):
                child.chmod(0o700)
            lane.chmod(0o700)
            shutil.rmtree(lane)
            replacement_inode = None
            for unused_attempt in range(4096):
                lane.mkdir(mode=0o700)
                lane.chmod(0o700)
                replacement_inode = lane.stat(follow_symlinks=False).st_ino
                if replacement_inode == \
                        registration["output_lane"]["inode"]:
                    break
                lane.rmdir()
            if replacement_inode != registration["output_lane"]["inode"]:
                self.skipTest(
                    "temporary filesystem did not recycle the lane inode")
            lane.chmod(0o700)
            with patched_live_acquisition(
                    fixture, enforce_output_lane_binding=True) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "ratified")
                self.assertFalse(controls["lock"].held)
                self.assertEqual(controls["qualification"].call_count, 0)
            self.assertEqual(list(lane.iterdir()), [])
            self.assertEqual(lane.stat().st_mode & 0o777, 0o700)

    def test_bound_authority_roots_never_overlap_the_output_lane(self):
        cases = (
            "lane-is-candidate", "lane-parent-of-candidate",
            "lane-child-of-candidate", "lane-is-main",
            "lane-parent-of-main", "lane-child-of-main",
            "candidate-parent-of-main",
        )

        def snapshot(root):
            observed = []
            for path in sorted(root.rglob("*")):
                relative = str(path.relative_to(root))
                mode = path.stat().st_mode & 0o777
                observed.append((
                    relative, mode,
                    path.read_bytes() if path.is_file() else None))
            return observed

        for case in cases:
            with self.subTest(case=case), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                root.chmod(0o700)
                candidate = root / "candidate"
                main = root / "main"
                lane = root / "lane"
                if case == "lane-is-candidate":
                    lane = candidate
                elif case == "lane-parent-of-candidate":
                    candidate = lane / "candidate"
                elif case == "lane-child-of-candidate":
                    lane = candidate / "lane"
                elif case == "lane-is-main":
                    lane = main
                elif case == "lane-parent-of-main":
                    main = lane / "main"
                elif case == "lane-child-of-main":
                    lane = main / "lane"
                elif case == "candidate-parent-of-main":
                    main = candidate / "main"
                for path in (candidate, main, lane):
                    path.mkdir(parents=True, exist_ok=True, mode=0o700)
                for path in (candidate, main, lane):
                    path.chmod(0o700)
                (candidate / "candidate-sentinel").write_bytes(b"candidate")
                (main / "main-sentinel").write_bytes(b"main")

                lane_was_empty = not any(lane.iterdir())
                fixture = live_fixture(
                    output_lane=lane if lane_was_empty else None)
                fixture["candidate_authority"]["root"] = candidate
                fixture["exact_main"]["root"] = main
                registration = fixture["pure"]["preregistration"]
                before = snapshot(root)
                with patched_live_acquisition(
                        fixture,
                        enforce_output_lane_binding=lane_was_empty) as controls, \
                        mock.patch.object(
                            runner, "_open_lane_and_lock") as open_lane:
                    self.assertRejected(lambda: runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=contract.canonical_json_bytes(
                            registration),
                        candidate_authority_lane=candidate,
                        exact_main_authority_lane=main))
                    open_lane.assert_not_called()
                    self.assertEqual(controls["qualification"].call_count, 0)
                self.assertEqual(snapshot(root), before)
                self.assertEqual(
                    list(root.rglob(".arming.lock")),
                    [lane / ".arming.lock"] if lane_was_empty else [])

    def test_nonempty_unbound_lane_rejects_before_lock_creation(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            sentinel = lane / "sentinel"
            sentinel.write_bytes(b"unchanged")
            before_mode = lane.stat().st_mode & 0o777
            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "empty owner-private")
                self.assertEqual(controls["qualification"].call_count, 0)
            self.assertEqual(sentinel.read_bytes(), b"unchanged")
            self.assertEqual(lane.stat().st_mode & 0o777, before_mode)
            self.assertEqual(sorted(path.name for path in lane.iterdir()),
                             ["sentinel"])

    def test_attempt_recovery_topology_rejects_before_any_mutation(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        cases = (
            ("out-of-budget-staging",
             (".attempt-999.staging-1-" + "0" * 32,), (),
             "unique next in-budget attempt"),
            ("multiple-staging",
             (".attempt-001.staging-1-" + "1" * 32,
              ".attempt-001.staging-2-" + "2" * 32), (),
             "multiple staging remnants"),
            ("wrong-next-staging",
             (".attempt-002.staging-1-" + "3" * 32,), (),
             "unique next in-budget attempt"),
            ("final-gap", (), ("attempt-002",),
             "gapless in-budget prefix"),
            ("final-over-budget", (),
             ("attempt-001", "attempt-002", "attempt-003", "attempt-004"),
             "gapless in-budget prefix"),
        )

        def snapshot(root):
            return [
                (str(path.relative_to(root)), path.lstat().st_mode & 0o777,
                 path.read_bytes() if path.is_file() else None)
                for path in sorted(root.rglob("*"))
            ]

        for label, staging_names, final_names, unused_diagnostic in cases:
            with self.subTest(case=label), \
                    tempfile.TemporaryDirectory() as directory:
                lane = private_lane(directory)
                lock = lane / ".arming.lock"
                lock.write_bytes(b"")
                lock.chmod(0o600)
                attempts = lane / runner.ATTEMPTS_DIRECTORY
                attempts.mkdir(mode=0o700)
                for name in staging_names:
                    (attempts / name).mkdir(mode=0o700)
                for name in final_names:
                    (attempts / name).mkdir(mode=0o500)
                before = snapshot(lane)
                with patched_live_acquisition(fixture) as controls:
                    self.assertRejected(lambda: runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority")))
                    self.assertEqual(controls["qualification"].call_count, 0)
                self.assertEqual(snapshot(lane), before)
                self.assertEqual(lane.stat().st_mode & 0o777, 0o700)

    def test_creation_umask_is_gated_before_lane_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            for masked in (0o700, 0o777):
                with self.subTest(mask=oct(masked)):
                    lane = private_lane(directory, f"masked-{masked:o}")
                    fixture = live_fixture(output_lane=lane)
                    preregistration_bytes = contract.canonical_json_bytes(
                        fixture["pure"]["preregistration"])
                    before = sorted(path.name for path in lane.iterdir())
                    previous = os.umask(masked)
                    try:
                        with patched_live_acquisition(
                                fixture,
                                enforce_output_lane_binding=True) as controls:
                            self.assertRejected(lambda: runner.acquire_and_arm(
                                lane=lane,
                                preregistration_bytes=preregistration_bytes,
                                candidate_authority_lane=Path(
                                    "/candidate-authority"),
                                exact_main_authority_lane=Path(
                                    "/exact-main-authority")),
                                "umask masks required owner permissions")
                            self.assertFalse(controls["lock"].held)
                            self.assertEqual(
                                controls["qualification"].call_count, 0)
                    finally:
                        os.umask(previous)
                    self.assertEqual(
                        sorted(path.name for path in lane.iterdir()), before)
                    self.assertEqual(lane.stat().st_mode & 0o777, 0o500)

            permitted = private_lane(directory, "permitted-077")
            permitted_fixture = live_fixture(output_lane=permitted)
            permitted_bytes = contract.canonical_json_bytes(
                permitted_fixture["pure"]["preregistration"])
            previous = os.umask(0o077)
            try:
                with patched_live_acquisition(
                        permitted_fixture,
                        enforce_output_lane_binding=True):
                    campaign = runner.acquire_and_arm(
                        lane=permitted,
                        preregistration_bytes=permitted_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
            finally:
                os.umask(previous)
            try:
                self.assertEqual(
                    (permitted / ".arming.lock").stat().st_mode & 0o777,
                    0o600)
                self.assertEqual(permitted.stat().st_mode & 0o777, 0o500)
            finally:
                campaign.close()

    def test_lane_identity_is_rechecked_after_preflight_and_at_initializer(self):
        with tempfile.TemporaryDirectory() as directory:
            failed_lane = private_lane(directory, "canonicalization-failure")
            baseline = len(os.listdir("/proc/self/fd"))
            retained_failed_fd = os.open(
                failed_lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            with mock.patch.object(
                    runner, "_canonical_lane",
                    side_effect=runner.ArmingError(
                        "injected canonical-lane failure")):
                self.assertRejected(lambda: runner._open_lane_and_lock(
                    failed_lane, retained_failed_fd, {},
                    evidence_attempt_limit=3, prearmed_attempt_limit=16),
                    "injected canonical-lane failure")
            self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)

            lane = private_lane(directory, "preflight-swap")
            lane_binding_for_test(lane)
            detached = Path(directory) / "preflight-detached"
            calls = 0
            real_preflight = runner._preflight_lane_inventory

            def swapping_preflight(*arguments, **keywords):
                nonlocal calls
                calls += 1
                result = real_preflight(*arguments, **keywords)
                if calls == 2:
                    lane.rename(detached)
                    lane.mkdir(mode=0o700)
                    replacement = lane / ".arming.lock"
                    replacement.write_bytes(b"")
                    replacement.chmod(0o600)
                return result

            baseline = len(os.listdir("/proc/self/fd"))
            retained_lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            with mock.patch.object(
                    runner, "_preflight_lane_inventory",
                    side_effect=swapping_preflight), mock.patch.object(
                        runner.prereg,
                        "validate_output_lane_descriptor_identity_for_preregistration",
                        return_value=lane):
                self.assertRejected(lambda: runner._open_lane_and_lock(
                    lane, retained_lane_fd, {}, evidence_attempt_limit=3,
                    prearmed_attempt_limit=16),
                    "changed while acquiring authority")
            self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
            self.assertEqual(detached.stat().st_mode & 0o777, 0o500)

            initial = private_lane(directory, "initializer-swap")
            lane_binding_for_test(initial)
            retained_lane_fd = os.open(
                initial, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            with mock.patch.object(
                    runner.prereg,
                    "validate_output_lane_descriptor_identity_for_preregistration",
                    return_value=initial):
                canonical, lane_fd, lock_fd = runner._open_lane_and_lock(
                    initial, retained_lane_fd, {}, evidence_attempt_limit=3,
                    prearmed_attempt_limit=16)
            initializer_detached = Path(directory) / "initializer-detached"
            initial.rename(initializer_detached)
            initial.mkdir(mode=0o700)
            try:
                self.assertRejected(lambda: runner._initialize_lane_locked(
                    canonical, lane_fd, lock_fd, {}, {}, {}, {}),
                    "changed while acquiring authority")
                self.assertEqual(
                    initializer_detached.stat().st_mode & 0o777, 0o500)
            finally:
                fcntl.flock(lock_fd, fcntl.LOCK_UN)
                os.close(lock_fd)
                os.close(lane_fd)

    def test_blocked_lane_lock_revalidates_path_before_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            lock_path = lane / ".arming.lock"
            lock_path.write_bytes(b"")
            lock_path.chmod(0o600)
            holder = os.open(lock_path, os.O_RDWR | os.O_CLOEXEC)
            real_flock = fcntl.flock
            real_flock(holder, fcntl.LOCK_EX)
            blocked = threading.Event()
            result: list[object] = []
            contender_lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)

            def observing_flock(descriptor, operation):
                if operation == fcntl.LOCK_EX:
                    blocked.set()
                return real_flock(descriptor, operation)

            def contender():
                try:
                    result.append(runner._open_lane_and_lock(
                        lane, contender_lane_fd, {}, evidence_attempt_limit=3,
                        prearmed_attempt_limit=16))
                except BaseException as error:
                    result.append(error)

            try:
                with mock.patch.object(
                        runner.fcntl, "flock", side_effect=observing_flock), \
                        mock.patch.object(
                            runner.prereg,
                            "validate_output_lane_descriptor_identity_for_preregistration",
                            return_value=lane):
                    worker = threading.Thread(target=contender)
                    worker.start()
                    self.assertTrue(blocked.wait(timeout=5))
                    orphan = Path(directory) / "orphaned-lane"
                    lane.rename(orphan)
                    lane.mkdir(mode=0o700)
                    (lane / "new-sentinel").write_bytes(b"new")
                    (lane / ".arming.lock").write_bytes(b"")
                    (lane / ".arming.lock").chmod(0o600)
                    real_flock(holder, fcntl.LOCK_UN)
                    worker.join(timeout=5)
                    self.assertFalse(worker.is_alive())
            finally:
                try:
                    real_flock(holder, fcntl.LOCK_UN)
                finally:
                    os.close(holder)
            self.assertEqual(len(result), 1)
            self.assertIsInstance(result[0], ValueError)
            self.assertIn("changed while acquiring", str(result[0]))
            self.assertEqual(
                sorted(path.name for path in orphan.iterdir()),
                [".arming.lock"])
            self.assertEqual(
                sorted(path.name for path in lane.iterdir()),
                [".arming.lock", "new-sentinel"])

    def test_campaign_marker_open_fsync_faults_reseal_and_close_every_fd(self):
        marker_name = "attempt-001-transcript.jsonl"
        with tempfile.TemporaryDirectory() as directory:
            lane = Path(directory)
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                for fail_at in (1, 2, 3):
                    binding = create_campaign_marker_slot_for_test(
                        lane_fd, marker_name)
                    target = (binding["device"], binding["inode"])
                    baseline = len(os.listdir("/proc/self/fd"))
                    calls = 0
                    real_fsync = os.fsync

                    def failing_fsync(descriptor):
                        nonlocal calls
                        metadata = os.fstat(descriptor)
                        if (metadata.st_dev, metadata.st_ino) == target:
                            calls += 1
                            if calls == fail_at:
                                raise OSError(errno.EIO, "injected fsync failure")
                        return real_fsync(descriptor)

                    with mock.patch.object(
                            runner.os, "fsync", side_effect=failing_fsync):
                        with self.assertRaisesRegex(OSError, "injected"):
                            runner._open_campaign_marker_for_append(
                                lane_fd, binding)
                    self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
                    self.assertEqual(
                        (lane / marker_name).stat().st_mode & 0o777, 0o400)
                    os.unlink(marker_name, dir_fd=lane_fd)
            finally:
                os.close(lane_fd)

    def test_acquisition_cleanup_survives_lane_unlock_failure(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            captured: list[int] = []
            real_open = runner._open_lane_and_lock
            real_flock = fcntl.flock

            def capture_open(
                path, lane_fd, registration, *,
                evidence_attempt_limit, prearmed_attempt_limit,
            ):
                value = real_open(
                    path, lane_fd, registration,
                    evidence_attempt_limit=evidence_attempt_limit,
                    prearmed_attempt_limit=prearmed_attempt_limit)
                captured.extend(value[1:])
                return value

            def failing_unlock(descriptor, operation):
                if operation == fcntl.LOCK_UN:
                    raise OSError(errno.EIO, "injected unlock failure")
                return real_flock(descriptor, operation)

            with patched_live_acquisition(fixture) as controls:
                baseline = len(os.listdir("/proc/self/fd"))
                with mock.patch.object(
                        runner, "_open_lane_and_lock",
                        side_effect=capture_open), \
                        mock.patch.object(
                            runner, "_initialize_lane_locked",
                            side_effect=RuntimeError(
                                "injected initializer failure")), \
                        mock.patch.object(
                            runner.fcntl, "flock",
                            side_effect=failing_unlock):
                    with self.assertRaisesRegex(RuntimeError, "initializer"):
                        runner.acquire_and_arm(
                            lane=lane,
                            preregistration_bytes=contract.canonical_json_bytes(
                                registration),
                            candidate_authority_lane=Path(
                                "/candidate-authority"),
                            exact_main_authority_lane=Path(
                                "/exact-main-authority"))
                self.assertFalse(controls["lock"].held)
                self.assertTrue(controls["runtime_guard"].closed)
                self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
            self.assertEqual(len(captured), 2)
            for descriptor in captured:
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    def test_close_identity_failure_still_releases_all_capabilities(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            canonical_path = Path(directory) / "canonical.lock"

            class DescriptorCanonicalLock:
                def __init__(self, unused_path):
                    self.descriptor = -1
                    self.identity = None

                def __enter__(self):
                    self.descriptor = os.open(
                        canonical_path,
                        os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o600)
                    fcntl.flock(self.descriptor, fcntl.LOCK_EX)
                    self.identity = os.fstat(self.descriptor)
                    return self

                def validate_current(self):
                    os.fstat(self.descriptor)

                def __exit__(self, unused_kind, unused_value, unused_traceback):
                    fcntl.flock(self.descriptor, fcntl.LOCK_UN)
                    os.close(self.descriptor)
                    self.descriptor = -1
                    self.identity = None

            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls, \
                    mock.patch.object(
                        runner.baseline_acquire, "CanonicalFileLock",
                        DescriptorCanonicalLock):
                baseline = len(os.listdir("/proc/self/fd"))
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                retained = [
                    campaign._main_fd, campaign._candidate_fd,
                    campaign._manifest_fd, campaign._bundle_fd,
                    campaign._armed_fd, campaign._launch_marker_fd,
                    campaign._journal_fd, campaign._campaign_marker_fd,
                    campaign._campaign_progress_fd,
                    campaign._attempt_fd, campaign._attempts_fd,
                    campaign._lane_binding_fd, campaign._lane_lock_fd,
                    campaign._lane_fd,
                ]
                with mock.patch.object(
                        runner, "_process_start_ticks",
                        side_effect=OSError(errno.EIO, "identity read failed")):
                    self.assertRejected(campaign.close, "cleanup failed")
                self.assertTrue(controls["runtime_guard"].closed)
                self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
                for descriptor in retained:
                    with self.assertRaises(OSError):
                        os.fstat(descriptor)
                contender = os.open(
                    canonical_path, os.O_RDWR | os.O_CLOEXEC)
                try:
                    fcntl.flock(contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fcntl.flock(contender, fcntl.LOCK_UN)
                finally:
                    os.close(contender)

    def test_transcript_recovery_covers_every_journal_entry_byte_boundary(self):
        armed = {
            "schema": "unit-test-armed-record",
            "evidence_attempt": 1,
            "authority_bundle_sha256": "1" * 64,
            "attempt_manifest_sha256": "2" * 64,
            "lane_binding_sha256": "3" * 64,
            "host_instance_sha256": "4" * 64,
            "selected_pair": {
                "benchmark_cpu": 8, "reserved_sibling": 72,
            },
        }
        journal_records = [
            ("intent-0000.json", hashlib.sha256(b"intent").hexdigest()),
            ("result-0000.json", hashlib.sha256(b"result").hexdigest()),
            ("complete.json", hashlib.sha256(b"complete").hexdigest()),
        ]
        expected = runner._expected_campaign_transcript_bytes(
            armed, journal_records)
        blocks = runner._expected_campaign_transcript_blocks(
            armed, journal_records)
        checkpoint_state = {
            "binding": {"markers": [{} for unused_block in blocks]},
            "markers": [
                {"mode": 0o400, "data": block} for block in blocks
            ],
            "frontier": len(blocks),
            "blocks": blocks,
            "payload_prefix": expected,
        }
        marker_name = "attempt-001-transcript.jsonl"
        with tempfile.TemporaryDirectory() as directory:
            lane = Path(directory)
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                for prefix_size in range(len(expected) + 1):
                    binding = create_campaign_marker_slot_for_test(
                        lane_fd, marker_name)
                    if prefix_size:
                        descriptor = runner._open_campaign_marker_for_append(
                            lane_fd, binding)
                        try:
                            self.assertEqual(
                                os.pwrite(
                                    descriptor, expected[:prefix_size], 0),
                                prefix_size)
                            os.fsync(descriptor)
                        finally:
                            os.close(descriptor)
                    if prefix_size % 2:
                        os.chmod(lane / marker_name, 0o600)
                    with mock.patch.object(
                            runner, "_validate_campaign_checkpoint_payloads",
                            return_value=checkpoint_state):
                        runner._reconcile_campaign_transcript_prefix(
                            lane_fd, binding, armed=armed,
                            journal_records=journal_records,
                            label=f"journal transcript prefix {prefix_size}",
                            progress_directory_fd=-1,
                            progress_directory_binding={})
                    self.assertEqual((lane / marker_name).read_bytes(), expected)
                    self.assertEqual(
                        (lane / marker_name).stat().st_mode & 0o777, 0o400)
                    os.unlink(marker_name, dir_fd=lane_fd)

                binding = create_campaign_marker_slot_for_test(
                    lane_fd, marker_name)
                descriptor = runner._open_campaign_marker_for_append(
                    lane_fd, binding)
                try:
                    os.pwrite(descriptor, b"!", 0)
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
                with mock.patch.object(
                        runner, "_validate_campaign_checkpoint_payloads",
                        return_value=checkpoint_state):
                    self.assertRejected(lambda:
                        runner._reconcile_campaign_transcript_prefix(
                            lane_fd, binding, armed=armed,
                            journal_records=journal_records,
                            label="non-prefix journal transcript",
                            progress_directory_fd=-1,
                            progress_directory_binding={}), "prefix")
                self.assertEqual((lane / marker_name).read_bytes(), b"!")
            finally:
                os.close(lane_fd)

    def test_launch_marker_recovery_covers_every_intent_byte_boundary(self):
        expected = contract.canonical_json_bytes({
            "schema": runner.JOURNAL_INTENT_SCHEMA,
            "child_index": 0,
            "durable_intent": "exact-prefix-recovery",
            "padding": "x" * 257,
        })
        with tempfile.TemporaryDirectory() as directory:
            attempt = Path(directory)
            attempt_fd = os.open(
                attempt, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                for prefix_size in range(len(expected) + 1):
                    descriptor = os.open(
                        runner.LAUNCH_CONSUMED_FILE,
                        os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
                        os.O_NOFOLLOW,
                        0o600, dir_fd=attempt_fd)
                    os.fsync(descriptor)
                    os.fchmod(descriptor, 0o400)
                    os.fsync(descriptor)
                    binding = runner._launch_marker_binding_record(
                        os.fstat(descriptor))
                    if prefix_size:
                        os.pwrite(descriptor, expected[:prefix_size], 0)
                        os.fsync(descriptor)
                    os.close(descriptor)
                    if prefix_size % 2:
                        os.chmod(
                            attempt / runner.LAUNCH_CONSUMED_FILE, 0o600)
                    runner._reconcile_launch_marker_prefix(
                        attempt_fd, binding, expected_data=expected,
                        label=f"launch marker prefix {prefix_size}")
                    marker = attempt / runner.LAUNCH_CONSUMED_FILE
                    self.assertEqual(marker.read_bytes(), expected)
                    self.assertEqual(marker.stat().st_mode & 0o777, 0o400)
                    marker.unlink()

                descriptor = os.open(
                    runner.LAUNCH_CONSUMED_FILE,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
                    os.O_NOFOLLOW,
                    0o600, dir_fd=attempt_fd)
                os.fsync(descriptor)
                os.fchmod(descriptor, 0o400)
                os.fsync(descriptor)
                binding = runner._launch_marker_binding_record(
                    os.fstat(descriptor))
                os.pwrite(descriptor, b"!", 0)
                os.fsync(descriptor)
                os.close(descriptor)
                self.assertRejected(lambda:
                    runner._reconcile_launch_marker_prefix(
                        attempt_fd, binding, expected_data=expected,
                        label="non-prefix launch marker"), "prefix")
            finally:
                os.close(attempt_fd)

    def test_atomic_journal_publication_faults_never_leave_partial_final(self):
        payloads = {
            "intent-0000.json": contract.canonical_json_bytes({
                "kind": "intent", "padding": "i" * 73}),
            "result-0000.json": contract.canonical_json_bytes({
                "kind": "result", "padding": "r" * 79}),
            "complete.json": contract.canonical_json_bytes({
                "kind": "complete", "padding": "c" * 83}),
        }
        with tempfile.TemporaryDirectory() as directory:
            journal = Path(directory)
            journal_fd = os.open(
                journal, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                real_write = os.write
                for name, data in payloads.items():
                    for prefix_size in range(len(data)):
                        calls = 0

                        def crashing_write(descriptor, content):
                            nonlocal calls
                            calls += 1
                            if calls == 1 and prefix_size:
                                return real_write(
                                    descriptor, content[:prefix_size])
                            raise OSError(errno.EIO, "injected journal write")

                        baseline = len(os.listdir("/proc/self/fd"))
                        with mock.patch.object(
                                runner.os, "write",
                                side_effect=crashing_write):
                            with self.assertRaisesRegex(OSError, "injected"):
                                runner._write_atomic_journal_file_at(
                                    journal_fd, name, data,
                                    f"faulted {name}")
                        self.assertEqual(os.listdir(journal), [])
                        self.assertEqual(
                            len(os.listdir("/proc/self/fd")), baseline)

                for fault in ("fchmod", "rename"):
                    name = "intent-0000.json"
                    data = payloads[name]
                    target = (runner.os if fault == "fchmod" else runner)
                    attribute = ("fchmod" if fault == "fchmod" else
                                 "_renameat2_noreplace")
                    with mock.patch.object(
                            target, attribute,
                            side_effect=OSError(
                                errno.EIO, f"injected {fault} failure")):
                        with self.assertRaisesRegex(OSError, "injected"):
                            runner._write_atomic_journal_file_at(
                                journal_fd, name, data,
                                f"faulted {fault}")
                    self.assertEqual(os.listdir(journal), [])

                real_fsync = os.fsync
                for fail_at in (1, 2, 3):
                    with tempfile.TemporaryDirectory() as fault_directory:
                        fault_journal = Path(fault_directory)
                        fault_fd = os.open(
                            fault_journal,
                            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
                        calls = 0

                        def failing_fsync(descriptor):
                            nonlocal calls
                            calls += 1
                            if calls == fail_at:
                                raise OSError(
                                    errno.EIO, "injected journal fsync")
                            return real_fsync(descriptor)

                        try:
                            with mock.patch.object(
                                    runner.os, "fsync",
                                    side_effect=failing_fsync):
                                with self.assertRaisesRegex(
                                        OSError, "injected"):
                                    runner._write_atomic_journal_file_at(
                                        fault_fd, "intent-0000.json",
                                        payloads["intent-0000.json"],
                                        "faulted fsync")
                            names = os.listdir(fault_journal)
                            if fail_at < 3:
                                self.assertEqual(names, [])
                            else:
                                self.assertEqual(names, ["intent-0000.json"])
                                final = fault_journal / "intent-0000.json"
                                self.assertEqual(
                                    final.read_bytes(),
                                    payloads["intent-0000.json"])
                                self.assertEqual(
                                    final.stat().st_mode & 0o777, 0o400)
                        finally:
                            os.close(fault_fd)
            finally:
                os.close(journal_fd)

    def test_durable_intent_prefixes_reconcile_before_consumed_rejection(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        for phase, partial_size in (("transcript", 17), ("marker", 19)):
            with self.subTest(phase=phase), \
                    tempfile.TemporaryDirectory() as directory:
                lane = private_lane(directory)
                with patched_live_acquisition(fixture):
                    campaign = runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
                target_fd = (
                    campaign._campaign_marker_fd
                    if phase == "transcript" else campaign._launch_marker_fd)
                target_calls = 0
                real_pwrite = os.pwrite

                def crash_after_durable_prefix(fd, data, offset):
                    nonlocal target_calls
                    if fd != target_fd:
                        return real_pwrite(fd, data, offset)
                    target_calls += 1
                    if target_calls == 1:
                        self.assertGreater(len(data), partial_size)
                        return real_pwrite(fd, data[:partial_size], offset)
                    raise OSError(
                        errno.EIO, f"injected {phase} prefix crash")

                try:
                    with patched_live_acquisition(fixture), mock.patch.object(
                            runner.ArmedCampaign, "_revalidate_authority",
                            return_value=None), \
                            mock.patch.object(
                                runner, "_run_child_process",
                                side_effect=AssertionError(
                                    "child ran before durable intent chain")) \
                                    as child_canary, \
                            mock.patch.object(
                                runner.subprocess, "Popen",
                                side_effect=AssertionError(
                                    "Popen ran before durable intent chain")) \
                                    as popen_canary, \
                            mock.patch.object(
                                runner.os, "pwrite",
                                side_effect=crash_after_durable_prefix):
                        with self.assertRaisesRegex(
                                OSError, f"injected {phase} prefix crash"):
                            campaign.run_all()
                    self.assertEqual(target_calls, 2)
                    self.assertEqual(child_canary.call_count, 0)
                    self.assertEqual(popen_canary.call_count, 0)
                finally:
                    record = copy.deepcopy(campaign.record)
                    campaign.close()

                intent_path = (
                    lane / runner.ATTEMPTS_DIRECTORY / "attempt-001" /
                    runner.JOURNAL_DIRECTORY / "intent-0000.json")
                intent_data = intent_path.read_bytes()
                records = [("intent-0000.json",
                            hashlib.sha256(intent_data).hexdigest())]
                expected_transcript = \
                    runner._expected_campaign_transcript_bytes(record, records)
                with patched_live_acquisition(fixture) as controls:
                    self.assertRejected(lambda: runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority")),
                        "permanently consumes")
                    self.assertEqual(controls["qualification"].call_count, 0)
                self.assertEqual(
                    (lane / "attempt-001-transcript.jsonl").read_bytes(),
                    expected_transcript)
                self.assertEqual(
                    (lane / runner.ATTEMPTS_DIRECTORY / "attempt-001" /
                     runner.LAUNCH_CONSUMED_FILE).read_bytes(), intent_data)

    def test_gen3_main_normalization_attributes_launched_raw_variant(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        logical_plan = runner.plan_contract.campaign_plan_record(
            preregistration=registration)
        child = next(
            child for child in logical_plan["child_plans"]
            if child["implementation"] == "main")
        cell = logical_plan["cells"][child["cell_index"]]
        canonical = {
            "validated": True,
            "exact_main_authority_attribution": {
                "pure_avx2": True,
                "main_source_commit":
                    runner.result_contract.EXACT_MAIN_SOURCE_COMMIT,
                "authority_record_sha256":
                    runner.prereg.EXACT_MAIN_AUTHORITY_RECORD_SHA256,
                "executable_sha256": runner.CANONICAL_RAW_SHA256,
            },
        }
        with mock.patch.object(
                runner.result_contract, "validate_result",
                side_effect=lambda *unused: copy.deepcopy(canonical)):
            normalized = runner._normalize_gen3_result(
                "main", {}, cell, registration)
            self.assertNotEqual(
                runner.PATH_VARIANT_RAW_SHA256, runner.CANONICAL_RAW_SHA256)
            self.assertEqual(
                normalized["exact_main_authority_attribution"]
                          ["executable_sha256"],
                runner.PATH_VARIANT_RAW_SHA256)
            forged = runner._journal_result_record(
                child_index=child["index"], intent_sha256="a" * 64,
                outcome="accepted", finished_monotonic_ns=2,
                elapsed_ns=1, returncode=0,
                stdout_sha256="b" * 64, stderr_sha256="c" * 64,
                payload={}, normalized=canonical, error=None)
            self.assertRejected(lambda: runner._validate_journal_result(
                forged, child=child, intent_sha256="a" * 64,
                registration=registration, logical_plan=logical_plan),
                "normalization")

    def test_bundle_descriptor_and_launch_context_fixed_points(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            try:
                artifacts = campaign._authority_bundle["artifact_bundle"]
                descriptors = copy.deepcopy(
                    campaign._authority_bundle["descriptor_binding"])
                runner._validate_descriptor_bundle(descriptors, artifacts)
                for field, value in (
                        ("seals", 0), ("size", 1), ("mode", 0o700),
                        ("creator_start_ticks", 0)):
                    malformed = copy.deepcopy(descriptors)
                    malformed["main"][field] = value
                    self.assertRejected(lambda malformed=malformed:
                        runner._validate_descriptor_bundle(
                            malformed, artifacts), "descriptor")
                extra = copy.deepcopy(descriptors)
                extra["candidate_control"]["unexpected"] = True
                self.assertRejected(lambda:
                    runner._validate_descriptor_bundle(extra, artifacts),
                    "descriptor")

                bundle = copy.deepcopy(campaign._authority_bundle)
                bundle["launch_context"]["rlimits"][1]["soft"] = 1
                bundle_data = contract.canonical_json_bytes(bundle)
                self.assertRejected(lambda:
                    runner._validate_authority_bundle_and_armed(
                        bundle, campaign.record, lane=campaign._lane,
                        registration=campaign.preregistration,
                        logical_plan=campaign.logical_plan,
                        candidate_authority=campaign._candidate_authority,
                        exact_main=campaign._exact_main,
                        lane_binding=campaign._lane_binding,
                        authority_bundle_data=bundle_data,
                        attempt_manifest_data=campaign._manifest_data),
                    "launch context differs")
            finally:
                campaign.close()

    def test_sealed_memfd_identity_rejects_same_byte_recreation_but_dup_works(self):
        content = synthetic_elf(16_384)
        expected = hashlib.sha256(content).hexdigest()
        first_fd, first = runner.capture_sealed_executable_bytes(
            content, expected_sha256=expected, expected_size=len(content),
            label="first", authority_relative_path="artifacts/first")
        second_fd, second = runner.capture_sealed_executable_bytes(
            content, expected_sha256=expected, expected_size=len(content),
            label="second", authority_relative_path="artifacts/second")
        duplicate = os.dup(first_fd)
        try:
            self.assertEqual(runner.revalidate_sealed_descriptor(
                duplicate, first["snapshot"], "duplicate"), first["snapshot"])
            self.assertNotEqual(
                (first["snapshot"]["device"], first["snapshot"]["inode"]),
                (second["snapshot"]["device"], second["snapshot"]["inode"]))
            self.assertNotEqual(contract.canonical_sha256(first["snapshot"]),
                                contract.canonical_sha256(second["snapshot"]))
            self.assertRejected(lambda: runner.revalidate_sealed_descriptor(
                second_fd, first["snapshot"], "separate same-byte memfd"),
                "changed")
            with self.assertRaises(OSError) as write_error:
                os.pwrite(first_fd, b"X", 4)
            self.assertEqual(write_error.exception.errno, errno.EPERM)
        finally:
            os.close(duplicate)
            os.close(second_fd)
            os.close(first_fd)

    @unittest.skipUnless(
        AUTHORITY_AVAILABLE, "retained exact-main authority is unavailable")
    def test_retained_exact_main_path_variant_replays_without_pooling(self) -> None:
        authority = retained_authority()
        self.assertEqual(authority["root"], SEALED_AUTHORITY)
        self.assertEqual(hashlib.sha256(authority["artifact_data"]).hexdigest(),
                         runner.PATH_VARIANT_RAW_SHA256)
        self.assertNotEqual(runner.PATH_VARIANT_RAW_SHA256,
                            runner.CANONICAL_RAW_SHA256)
        self.assertEqual(
            authority["record"]["identity"]["path_variant"]["combined_sha256"],
            runner.NORMALIZED_COMBINED_SHA256)
        with tempfile.TemporaryDirectory() as directory:
            linked = Path(directory) / "authority-link"
            linked.symlink_to(SEALED_AUTHORITY, target_is_directory=True)
            self.assertRejected(
                lambda: runner.bind_exact_main_path_variant(linked), "direct")

    def test_candidate_authority_replays_all_record_inventory_joins(self) -> None:
        registration, payload = candidate_authority_fixture()
        tree = FakeSealedTree(payload)
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(runner.authority_verifier,
                                  "read_sealed_tree", return_value=tree), \
                mock.patch.object(runner.authority_verifier,
                                  "verify_tree_metadata", return_value=None), \
                mock.patch.object(
                    runner.authority_verifier, "verify_sha256sums",
                    return_value={"sha256":
                                  registration["candidate_executable"]
                                              ["authority_ledger_sha256"]}), \
                mock.patch.object(
                    runner.git_capture, "validate_git_capture",
                    return_value={"tree":
                                  registration["candidate_source"]["tree"]}), \
                mock.patch.object(
                    runner, "_validate_candidate_source_archive",
                    return_value=None), \
                mock.patch.object(
                    runner.build_provenance,
                    "validate_reproducible_build_proof", return_value=None), \
                mock.patch.object(
                    runner.build_provenance, "compare_reproducible_builds",
                    return_value=contract.strict_json_loads(
                        payload["build/reproducible-build-core.json"])), \
                mock.patch.object(
                    runner.build_provenance,
                    "_canonical_replay_attestation_header_bytes",
                    return_value=payload[
                        "build/benchmark-source-attestation.h"]):
            lane = private_lane(directory, "candidate-authority")
            observed = runner.bind_candidate_authority_lane(
                lane, registration)
            self.assertEqual(observed["artifact_data"],
                             payload["artifacts/bench_leopard2"])
            self.assertEqual(observed["source_authority"]["commit"],
                             registration["candidate_source"]["commit"])
            self.assertTrue(tree.reverified)

            damaged = dict(payload)
            damaged["source/candidate-source.tar"] += b"damage"
            with mock.patch.object(
                    runner.authority_verifier, "read_sealed_tree",
                    return_value=FakeSealedTree(damaged)):
                self.assertRejected(lambda: runner.bind_candidate_authority_lane(
                    lane, registration), "payload differs")

    def test_candidate_source_archive_parser_rejects_arbitrary_bytes(self) -> None:
        captured = {"tracked_files": [{
            "path": "tracked.cpp", "kind": "regular", "git_mode": "100644",
            "object_id": "1" * 40,
        }]}
        self.assertRejected(
            lambda: runner._validate_candidate_source_archive(
                b"not a canonical tar archive", captured, "2" * 40),
            "archive")

    def test_runtime_guard_retains_files_and_detects_live_mutation(self) -> None:
        def runtime_entry(path: Path, role: str) -> dict:
            data = path.read_bytes()
            metadata = path.stat()
            return {
                "path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
                "size": len(data), "mode": metadata.st_mode & 0o7777,
                "role": role,
            }

        with tempfile.TemporaryDirectory() as directory:
            dependency = Path(directory).resolve() / "libunit.so"
            dependency.write_bytes(b"runtime-v1")
            candidate = {
                "schema": "leopard2-k65-gen3-runtime-closure/v1",
                "dependencies": [runtime_entry(dependency, "dependency")],
                "launchers": [
                    runtime_entry(Path("/usr/bin/prlimit"), "launcher"),
                    runtime_entry(Path("/usr/bin/taskset"), "launcher"),
                ],
            }
            duplicate = {
                **candidate,
                "launchers": [candidate["launchers"][0],
                              candidate["launchers"][0],
                              candidate["launchers"][1]],
            }
            self.assertRejected(
                lambda: runner._candidate_runtime_closure(duplicate),
                "metadata differs")
            oversized = {
                **candidate,
                "dependencies": [{
                    "path": f"/runtime/{index:03d}", "sha256": "0" * 64,
                    "size": 1, "mode": 0o444, "role": "dependency",
                } for index in range(runner.MAX_RUNTIME_CLOSURE_FILES - 1)],
            }
            self.assertRejected(
                lambda: runner._candidate_runtime_closure(oversized),
                "metadata differs")
            dep = runtime_entry(dependency, "dependency")
            exact = {"runtime_closure": {"records": [
                {"role": "canonical_first"},
                {"role": "canonical_second"},
                {
                    "role": "path_variant",
                    "executable_sha256": runner.PATH_VARIANT_RAW_SHA256,
                    "dependencies": [{
                        "kind": "file", "path": dep["path"],
                        "sha256": dep["sha256"], "size": dep["size"],
                    }],
                },
            ]}}
            runner._validate_runtime_closures_current(candidate, exact)
            guard = runner._acquire_runtime_guard(candidate, exact)
            try:
                guard.validate_current("unit stable runtime")
                dependency.write_bytes(b"runtime-v2")
                self.assertRejected(
                    lambda: guard.validate_current("unit changed runtime"),
                    "runtime")
            finally:
                guard.close()

    def test_atomic_arming_precedes_any_child_and_retains_locks(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                self.assertEqual(controls["popen"].call_count, 0)
                self.assertTrue(controls["lock"].held)
                self.assertGreater(controls["lock"].validations, 0)
                self.assertGreater(controls["runtime_guard"].validations, 0)
            try:
                attempt = lane / "attempts" / "attempt-001"
                self.assertEqual(campaign.attempt_path, attempt)
                self.assertEqual(set(path.name for path in attempt.iterdir()), {
                    runner.ARMED_FILE, runner.AUTHORITY_BUNDLE_FILE,
                    runner.ATTEMPT_MANIFEST_FILE, runner.JOURNAL_DIRECTORY,
                    runner.LAUNCH_CONSUMED_FILE,
                })
                self.assertEqual(attempt.stat().st_mode & 0o777, 0o500)
                self.assertEqual(lane.stat().st_mode & 0o777, 0o500)
                self.assertEqual(
                    (lane / runner.ATTEMPTS_DIRECTORY).stat().st_mode & 0o777,
                    0o500)
                self.assertEqual((attempt / runner.JOURNAL_DIRECTORY).stat().st_mode
                                 & 0o777, 0o700)
                for name in (runner.ARMED_FILE, runner.AUTHORITY_BUNDLE_FILE,
                             runner.ATTEMPT_MANIFEST_FILE):
                    self.assertEqual((attempt / name).stat().st_mode & 0o777,
                                     0o400)
                transcript_names = {
                    f"attempt-{index:03d}-transcript.jsonl"
                    for index in range(
                        1, registration["budgets"]["evidence_attempts"] + 1)
                }
                self.assertTrue(transcript_names <= {
                    path.name for path in lane.iterdir()})
                for index, name in enumerate(sorted(transcript_names), 1):
                    slot = lane / name
                    self.assertEqual(slot.stat().st_mode & 0o777, 0o400)
                    if index == 1:
                        self.assertGreater(slot.stat().st_size, 0)
                        transcript = runner._read_campaign_transcript(
                            campaign._lane_fd, campaign._campaign_marker_fd,
                            campaign._campaign_marker_binding,
                            "unit allocated transcript")
                        self.assertIsNotNone(transcript["allocation"])
                        self.assertEqual(transcript["journal_entries"], [])
                    else:
                        self.assertEqual(slot.stat().st_size, 0)
            finally:
                campaign.close()
            self.assertFalse(controls["lock"].held)
            self.assertTrue(controls["runtime_guard"].closed)

    def test_allocation_prefix_recovery_covers_every_byte_boundary(self):
        armed = {
            "schema": "unit-test-armed-record",
            "evidence_attempt": 1,
            "authority_bundle_sha256": "1" * 64,
            "attempt_manifest_sha256": "2" * 64,
            "lane_binding_sha256": "3" * 64,
            "host_instance_sha256": "4" * 64,
            "selected_pair": {
                "benchmark_cpu": 8, "reserved_sibling": 72,
            },
        }
        allocation_data = contract.canonical_json_bytes(
            runner._campaign_transcript_allocation_record(armed))
        marker_name = "attempt-001-transcript.jsonl"
        with tempfile.TemporaryDirectory() as directory:
            lane = Path(directory)
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                for prefix_size in range(len(allocation_data)):
                    binding = create_campaign_marker_slot_for_test(
                        lane_fd, marker_name)
                    if prefix_size:
                        marker_fd = runner._open_campaign_marker_for_append(
                            lane_fd, binding)
                        try:
                            self.assertEqual(os.pwrite(
                                marker_fd, allocation_data[:prefix_size], 0),
                                prefix_size)
                            os.fsync(marker_fd)
                        finally:
                            os.close(marker_fd)
                    if prefix_size % 2:
                        os.chmod(lane / marker_name, 0o600)
                    marker_fd = runner._open_campaign_marker_for_append(
                        lane_fd, binding,
                        associated_with_published_attempt=True)
                    try:
                        runner._append_campaign_transcript_allocation(
                            lane_fd, marker_fd, binding, armed)
                        transcript = runner._read_campaign_transcript(
                            lane_fd, marker_fd, binding,
                            f"allocation prefix {prefix_size}")
                    finally:
                        os.close(marker_fd)
                    self.assertEqual(
                        transcript["allocation"]["evidence_attempt"], 1)
                    self.assertEqual(transcript["journal_entries"], [])
                    os.unlink(marker_name, dir_fd=lane_fd)

                binding = create_campaign_marker_slot_for_test(
                    lane_fd, marker_name)
                marker_fd = runner._open_campaign_marker_for_append(
                    lane_fd, binding)
                try:
                    self.assertEqual(os.pwrite(marker_fd, b"!", 0), 1)
                    os.fsync(marker_fd)
                finally:
                    os.close(marker_fd)
                marker_fd = runner._open_campaign_marker_for_append(
                    lane_fd, binding,
                    associated_with_published_attempt=True)
                try:
                    self.assertRejected(lambda:
                        runner._append_campaign_transcript_allocation(
                            lane_fd, marker_fd, binding, armed),
                        "exact ARMED allocation prefix")
                finally:
                    os.close(marker_fd)
                self.assertEqual((lane / marker_name).read_bytes(), b"!")
            finally:
                os.close(lane_fd)

    def test_live_top_level_inventory_drift_fails_without_repair(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                junk = lane / ".lane.json.staging-attacker"
                lane.chmod(0o700)
                junk.write_bytes(b"must survive rejection")
                junk.chmod(0o400)
                lane.chmod(0o500)
                journal = campaign.attempt_path / runner.JOURNAL_DIRECTORY
                before = sorted(path.name for path in journal.iterdir())
                try:
                    self.assertRejected(campaign.run_all)
                    self.assertEqual(controls["popen"].call_count, 0)
                    self.assertEqual(
                        sorted(path.name for path in journal.iterdir()), before)
                    self.assertEqual(
                        junk.read_bytes(), b"must survive rejection")
                    self.assertEqual(junk.stat().st_mode & 0o777, 0o400)
                finally:
                    campaign.close()

    def test_live_older_checkpoint_rewrite_fails_before_result_publication(
            self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                progress = (
                    lane / campaign._campaign_progress_binding["name"])
                checkpoint = (
                    progress /
                    campaign._campaign_progress_binding["markers"][0]["name"])
                original = checkpoint.read_bytes()
                self.assertGreater(len(original), 1)

                def rewrite_allocation_checkpoint(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    self.assertGreater(descriptor, 0)
                    checkpoint.chmod(0o600)
                    checkpoint.write_bytes(b"!" + original[1:])
                    checkpoint.chmod(0o400)
                    return {
                        "outcome": "nonzero", "returncode": 7,
                        "stdout": b"", "stderr": b"rewritten",
                        "elapsed_ns": 1, "error": "rewritten checkpoint",
                    }

                try:
                    with mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=rewrite_allocation_checkpoint):
                        self.assertRejected(campaign.run_all, "checkpoint")
                    self.assertEqual(controls["popen"].call_count, 0)
                    self.assertFalse(
                        (campaign.attempt_path / runner.JOURNAL_DIRECTORY /
                         "result-0000.json").exists())
                    self.assertEqual(
                        checkpoint.read_bytes(), b"!" + original[1:])
                    self.assertEqual(checkpoint.stat().st_mode & 0o777, 0o400)
                finally:
                    campaign.close()

    def test_live_unlinked_prearmed_history_rewind_fails_without_repair(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as failed_controls:
                failed_controls["qualification"].side_effect = RuntimeError(
                    "injected environment rejection")
                with self.assertRaisesRegex(
                        RuntimeError, "injected environment rejection"):
                    runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
            history = lane / "prearmed-0001-history.jsonl"
            original = history.read_bytes()
            blocks = original.splitlines(keepends=True)
            self.assertGreaterEqual(len(blocks), 2)
            self.assertEqual(history.stat().st_mode & 0o777, 0o400)

            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                history.chmod(0o600)
                history.write_bytes(blocks[0])
                history.chmod(0o400)
                journal = campaign.attempt_path / runner.JOURNAL_DIRECTORY
                try:
                    self.assertRejected(
                        campaign.run_all, "live pre-ARMED history")
                    self.assertEqual(controls["popen"].call_count, 0)
                    self.assertEqual(list(journal.iterdir()), [])
                    self.assertEqual(history.read_bytes(), blocks[0])
                    self.assertEqual(history.stat().st_mode & 0o777, 0o400)
                finally:
                    campaign.close()

    def test_final_attempt_without_allocation_header_or_reseal_is_reconciled(
        self,
    ):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture), mock.patch.object(
                    runner, "_open_campaign_marker_for_append",
                    side_effect=runner.ArmingError(
                        "injected allocation-open failure")):
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "injected allocation-open failure")
            frozen_pair = copy.deepcopy(
                fixture["qualification"]["qualification_binding"]
                       ["selected_pair"])
            slot = lane / "attempt-001-transcript.jsonl"
            self.assertEqual(slot.stat().st_size, 0)
            slot.chmod(0o600)
            slot_fd = os.open(slot, os.O_RDONLY | os.O_CLOEXEC)
            try:
                os.fsync(slot_fd)
            finally:
                os.close(slot_fd)
            self.assertEqual(slot.stat().st_mode & 0o777, 0o600)
            self.assertEqual(lane.stat().st_mode & 0o777, 0o500)
            self.assertEqual(
                (lane / runner.ATTEMPTS_DIRECTORY).stat().st_mode & 0o777,
                0o500)

            retry_fixture = live_fixture(expected_frozen_pair=frozen_pair)
            with patched_live_acquisition(retry_fixture):
                retry = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            try:
                self.assertEqual(retry.record["evidence_attempt"], 2)
                first_fd = os.open(slot, os.O_RDONLY | os.O_CLOEXEC)
                try:
                    first_binding = retry._lane_binding["campaign_markers"][0]
                    transcript = runner._read_campaign_transcript(
                        retry._lane_fd, first_fd, first_binding,
                        "reconciled allocation transcript")
                finally:
                    os.close(first_fd)
                self.assertEqual(
                    transcript["allocation"]["evidence_attempt"], 1)
                self.assertEqual(transcript["journal_entries"], [])
            finally:
                retry.close()

    def test_partial_allocation_header_is_completed_before_retry(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        allocation_schema = \
            runner.CAMPAIGN_TRANSCRIPT_ALLOCATION_SCHEMA.encode("ascii")
        real_pwrite = os.pwrite
        target_fd = -1

        def crash_during_allocation(fd, data, offset):
            nonlocal target_fd
            descriptor_path = os.readlink(f"/proc/self/fd/{fd}")
            if (target_fd < 0 and offset == 0 and
                    descriptor_path.endswith(
                        "/attempt-001-transcript.jsonl") and
                    allocation_schema in data):
                target_fd = fd
                return real_pwrite(fd, data[:17], offset)
            if fd == target_fd and offset == 17:
                raise runner.ArmingError(
                    "injected crash during allocation header")
            return real_pwrite(fd, data, offset)

        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture), mock.patch.object(
                    runner.os, "pwrite", side_effect=crash_during_allocation):
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "injected crash during allocation header")
            attempts = lane / runner.ATTEMPTS_DIRECTORY
            staging_names = [
                path for path in attempts.iterdir()
                if runner.STAGING_NAME.fullmatch(path.name) is not None
            ]
            self.assertEqual(len(staging_names), 1)
            attempt = staging_names[0]
            slot = lane / "attempt-001-transcript.jsonl"
            armed = contract.strict_json_loads(
                (attempt / runner.ARMED_FILE).read_bytes(),
                "partially allocated unit ARMED record")
            allocation_data = contract.canonical_json_bytes(
                runner._campaign_transcript_allocation_record(armed))
            self.assertEqual(slot.read_bytes(), allocation_data[:17])
            self.assertEqual(slot.stat().st_mode & 0o777, 0o400)
            checkpoint_zero = list(
                (lane / "attempt-001-journal-checkpoints").glob(
                    "checkpoint-0000-*"))
            self.assertEqual(len(checkpoint_zero), 1)
            self.assertEqual(checkpoint_zero[0].read_bytes(), allocation_data)
            self.assertEqual(checkpoint_zero[0].stat().st_mode & 0o777, 0o400)

            frozen_pair = copy.deepcopy(
                fixture["qualification"]["qualification_binding"]
                       ["selected_pair"])
            retry_fixture = live_fixture(expected_frozen_pair=frozen_pair)
            with patched_live_acquisition(retry_fixture) as controls:
                retry = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            try:
                self.assertEqual(controls["qualification"].call_count, 1)
                self.assertEqual(retry.record["evidence_attempt"], 2)
                first_fd = os.open(slot, os.O_RDONLY | os.O_CLOEXEC)
                try:
                    first_binding = retry._lane_binding["campaign_markers"][0]
                    transcript = runner._read_campaign_transcript(
                        retry._lane_fd, first_fd, first_binding,
                        "completed partial allocation transcript")
                finally:
                    os.close(first_fd)
                self.assertEqual(
                    transcript["allocation"]["evidence_attempt"], 1)
                self.assertEqual(transcript["journal_entries"], [])
                self.assertFalse(staging_names[0].exists())
                self.assertTrue((attempts / "attempt-001").is_dir())
            finally:
                retry.close()

    def test_invalid_prior_chain_rejects_before_allocation_recovery_write(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                first = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            frozen_pair = copy.deepcopy(first.record["selected_pair"])
            first.close()

            real_armed_record = runner.execution.armed_record

            def wrong_prior_armed_record(*arguments, **keywords):
                changed = dict(keywords)
                if changed.get("prior_armed_chain_sha256") is not None:
                    changed["prior_armed_chain_sha256"] = "f" * 64
                return real_armed_record(*arguments, **changed)

            real_open_campaign_marker = \
                runner._open_campaign_marker_for_append

            def fail_attempt_two_allocation(lane_fd, binding, **keywords):
                if binding["name"] == "attempt-002-transcript.jsonl":
                    raise runner.ArmingError(
                        "injected attempt-2 allocation-open failure")
                return real_open_campaign_marker(
                    lane_fd, binding, **keywords)

            second_fixture = live_fixture(expected_frozen_pair=frozen_pair)
            with patched_live_acquisition(second_fixture), \
                    mock.patch.object(
                        runner.execution, "armed_record",
                        side_effect=wrong_prior_armed_record), \
                    mock.patch.object(
                        runner, "_open_campaign_marker_for_append",
                        side_effect=fail_attempt_two_allocation):
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "injected attempt-2 allocation-open failure")

            second_slot = lane / "attempt-002-transcript.jsonl"
            before = second_slot.read_bytes()
            before_mode = second_slot.stat().st_mode & 0o777
            self.assertEqual(before, b"")
            self.assertEqual(before_mode, 0o400)
            attempts = lane / runner.ATTEMPTS_DIRECTORY
            prepared = [
                path for path in attempts.iterdir()
                if (match := runner.STAGING_NAME.fullmatch(path.name))
                is not None and int(match.group(1)) == 2
            ]
            self.assertEqual(len(prepared), 1)
            with patched_live_acquisition(second_fixture) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "attempt 2 ARMED chain differs")
                self.assertEqual(controls["qualification"].call_count, 0)
            self.assertEqual(second_slot.read_bytes(), before)
            self.assertEqual(
                second_slot.stat().st_mode & 0o777, before_mode)
            self.assertTrue(prepared[0].is_dir())

    def test_mutated_committed_attempt_fails_before_new_qualification(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            campaign.close()
            bundle = (lane / "attempts" / "attempt-001" /
                      runner.AUTHORITY_BUNDLE_FILE)
            slot = lane / "attempt-001-transcript.jsonl"
            bundle.chmod(0o600)
            slot.chmod(0o600)
            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "immutable")
                self.assertEqual(controls["qualification"].call_count, 0)
                self.assertEqual(controls["popen"].call_count, 0)
            self.assertEqual(lane.stat().st_mode & 0o777, 0o500)
            self.assertEqual(
                (lane / runner.ATTEMPTS_DIRECTORY).stat().st_mode & 0o777,
                0o500)
            self.assertEqual(slot.stat().st_mode & 0o777, 0o600)

    def test_qualification_rejection_reseals_lane_and_attempts(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                controls["qualification"].side_effect = RuntimeError(
                    "simulated qualification rejection")
                with self.assertRaisesRegex(RuntimeError, "simulated"):
                    runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=contract.canonical_json_bytes(
                            registration),
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
            self.assertEqual(lane.stat().st_mode & 0o777, 0o500)
            self.assertEqual(
                (lane / runner.ATTEMPTS_DIRECTORY).stat().st_mode & 0o777,
                0o500)
            binding = contract.strict_json_loads(
                (lane / runner.LANE_BINDING_FILE).read_bytes(),
                "qualification-rejection lane binding")
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            try:
                bound_registration = fixture["pure"]["preregistration"]
                authority = runner._read_budget_authority_before_recovery(
                    lane_fd, registration=bound_registration,
                    logical_plan=fixture["pure"]["plan"],
                    lane_binding=binding, prior_attempts=[])
            finally:
                os.close(lane_fd)
            self.assertEqual(
                authority["ledger"]["environment_rejected_used"], 1)
            self.assertEqual(authority["ledger"]["setup_invalid_used"], 0)
            self.assertEqual(
                authority["histories"][0]["highest_state"], "QUALIFYING")

    def test_prearmed_setup_and_environment_budgets_are_durable(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)

        def acquire(lane):
            return runner.acquire_and_arm(
                lane=lane, preregistration_bytes=preregistration_bytes,
                candidate_authority_lane=Path("/candidate-authority"),
                exact_main_authority_lane=Path("/exact-main-authority"))

        with tempfile.TemporaryDirectory() as directory:
            environment_lane = private_lane(directory, "environment-budget")
            for index in range(1, 9):
                with self.subTest(environment_rejection=index), \
                        patched_live_acquisition(fixture) as controls:
                    controls["qualification"].side_effect = RuntimeError(
                        f"environment rejection {index}")
                    with self.assertRaisesRegex(
                            RuntimeError, "environment rejection"):
                        acquire(environment_lane)
            before = sorted(
                (path.name, path.stat().st_mode & 0o777, path.stat().st_size)
                for path in environment_lane.iterdir())
            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: acquire(environment_lane),
                                    "environment-rejected budget is exhausted")
                self.assertEqual(controls["qualification"].call_count, 0)
            after = sorted(
                (path.name, path.stat().st_mode & 0o777, path.stat().st_size)
                for path in environment_lane.iterdir())
            self.assertEqual(before, after)
            self.assertEqual(
                list((environment_lane / runner.ATTEMPTS_DIRECTORY).iterdir()),
                [])

            setup_lane = private_lane(directory, "setup-budget")
            for index in range(1, 6):
                with self.subTest(setup_rejection=index), \
                        patched_live_acquisition(fixture), \
                        mock.patch.object(
                            runner, "_validate_launch_context_current",
                            side_effect=runner.ArmingError(
                                f"setup rejection {index}")):
                    self.assertRejected(lambda: acquire(setup_lane),
                                        "setup rejection")
            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: acquire(setup_lane),
                                    "setup-invalid budget is exhausted")
                self.assertEqual(controls["qualification"].call_count, 0)

            for lane, setup_used, environment_used in (
                    (environment_lane, 0, 8), (setup_lane, 5, 0)):
                binding = contract.strict_json_loads(
                    (lane / runner.LANE_BINDING_FILE).read_bytes(),
                    "budget-exhaustion lane binding")
                arguments = fixture["pure_fixture_args"]
                lane_fixture = pure_fixtures.fixture(
                    arguments["candidate_sha256"],
                    candidate_authority_record_sha256=
                        arguments["candidate_authority_record_sha256"],
                    candidate_authority_ledger_sha256=
                        arguments["candidate_authority_ledger_sha256"],
                    controller_bindings=arguments["controller_bindings"],
                    output_lane_binding=lane_binding_for_test(lane),
                    expected_frozen_pair=arguments["expected_frozen_pair"])
                lane_fd = os.open(
                    lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
                try:
                    authority = runner._read_budget_authority_before_recovery(
                        lane_fd,
                        registration=lane_fixture["preregistration"],
                        logical_plan=lane_fixture["plan"],
                        lane_binding=binding, prior_attempts=[])
                finally:
                    os.close(lane_fd)
                self.assertEqual(
                    authority["ledger"]["setup_invalid_used"], setup_used)
                self.assertEqual(
                    authority["ledger"]["environment_rejected_used"],
                    environment_used)
                self.assertEqual(
                    authority["ledger"]["evidence_attempts_used"], 0)

    def test_first_prearmed_transition_crash_consumes_setup_at_every_byte(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        logical_plan = fixture["pure"]["plan"]
        prototype_lane_binding = {"schema": "unit-first-transition-lane/v1"}
        prototype_entry = runner._prearmed_state_entry_record(
            acquisition_index=1, sequence=1, from_state="INIT",
            to_state="PREREGISTERED", prior_entry_sha256=None,
            preregistration_sha256=contract.canonical_sha256(registration),
            plan_sha256=contract.canonical_sha256(logical_plan),
            lane_binding_sha256=contract.canonical_sha256(
                prototype_lane_binding))
        block_size = len(contract.canonical_json_bytes(prototype_entry))
        real_pwrite = os.pwrite
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for boundary in range(block_size):
                with self.subTest(boundary=boundary):
                    lane = root / f"boundary-{boundary:04d}"
                    lane.mkdir(mode=0o700)
                    lane_fd = os.open(
                        lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
                    binding = create_prearmed_history_slot_for_test(
                        lane_fd, "prearmed-0001-history.jsonl")
                    boundary_bindings = [
                        create_prearmed_boundary_slot_for_test(
                            lane_fd,
                            f"prearmed-0001-{state.lower()}-reached.marker")
                        for state in runner.PREARMED_BOUNDARY_STATES
                    ]
                    lane_binding = {
                        "schema": "unit-first-transition-lane/v1",
                        "prearmed_history_markers": [binding],
                        "prearmed_boundary_markers": boundary_bindings,
                    }
                    entry = runner._prearmed_state_entry_record(
                        acquisition_index=1, sequence=1,
                        from_state="INIT", to_state="PREREGISTERED",
                        prior_entry_sha256=None,
                        preregistration_sha256=
                            contract.canonical_sha256(registration),
                        plan_sha256=contract.canonical_sha256(logical_plan),
                        lane_binding_sha256=
                            contract.canonical_sha256(lane_binding))
                    block = contract.canonical_json_bytes(entry)
                    self.assertEqual(len(block), block_size)
                    runner._flip_prearmed_boundary_marker(
                        lane_fd, boundary_bindings[0], "PREREGISTERED")
                    marker_fd = runner._open_prearmed_history_for_allocation(
                        lane_fd, binding, registration=registration,
                        logical_plan=logical_plan,
                        lane_binding=lane_binding)

                    def crash_during_first_entry(
                        descriptor, data, offset, *, target=marker_fd,
                        prefix=boundary,
                    ):
                        if descriptor == target and offset == 0:
                            if prefix:
                                self.assertEqual(
                                    real_pwrite(
                                        descriptor, data[:prefix], offset),
                                    prefix)
                            raise OSError(errno.EIO, "injected first entry")
                        return real_pwrite(descriptor, data, offset)

                    try:
                        with mock.patch.object(
                                runner.os, "pwrite",
                                side_effect=crash_during_first_entry):
                            with self.assertRaisesRegex(
                                    OSError, "injected first entry"):
                                runner._append_prearmed_state(
                                    lane_fd, marker_fd, binding,
                                    "PREREGISTERED",
                                    registration=registration,
                                    logical_plan=logical_plan,
                                    lane_binding=lane_binding)
                        history = runner._read_prearmed_history(
                            lane_fd, marker_fd, binding,
                            registration=registration,
                            logical_plan=logical_plan,
                            lane_binding=lane_binding,
                            label="partially written first state")
                        self.assertTrue(history["used"])
                        self.assertEqual(history["mode"], 0o600)
                        self.assertEqual(history["data"], block[:boundary])
                        history["boundary"] = \
                            runner._read_prearmed_boundary_frontier(
                                lane_fd, lane_binding, 1)
                        history["effective_states"] = list(
                            history["boundary"]["states"])
                        history["inode_used"] = history["used"]
                        runner._recover_prearmed_histories(
                            lane_fd, {"observed": [history]},
                            registration=registration,
                            logical_plan=logical_plan,
                            lane_binding=lane_binding)
                        sealed = runner._read_prearmed_history(
                            lane_fd, marker_fd, binding,
                            registration=registration,
                            logical_plan=logical_plan,
                            lane_binding=lane_binding,
                            label="recovered first state")
                        self.assertTrue(sealed["used"])
                        self.assertEqual(sealed["mode"], 0o400)
                        self.assertEqual(
                            sealed["states"], ["INIT", "PREREGISTERED"])
                        self.assertEqual(sealed["data"], block)
                    finally:
                        os.close(marker_fd)
                        os.close(lane_fd)

    def test_zero_byte_first_transition_crash_survives_acquire_cleanup(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        schema = runner.PREARMED_STATE_ENTRY_SCHEMA.encode("ascii")
        real_pwrite = os.pwrite
        injected = False

        def crash_before_first_entry(descriptor, data, offset):
            nonlocal injected
            if not injected and offset == 0 and schema in data:
                injected = True
                raise OSError(errno.EIO, "injected zero-byte first entry")
            return real_pwrite(descriptor, data, offset)

        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture), mock.patch.object(
                    runner.os, "pwrite", side_effect=crash_before_first_entry):
                with self.assertRaisesRegex(
                        OSError, "injected zero-byte first entry"):
                    runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
            self.assertTrue(injected)
            first_history = lane / "prearmed-0001-history.jsonl"
            self.assertEqual(first_history.read_bytes(), b"")
            self.assertEqual(first_history.stat().st_mode & 0o777, 0o600)

            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            campaign.close()
            self.assertEqual(first_history.stat().st_mode & 0o777, 0o400)
            binding = contract.strict_json_loads(
                (lane / runner.LANE_BINDING_FILE).read_bytes(),
                "zero-byte crash lane binding")
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            history_fd = os.open(
                first_history, os.O_RDONLY | os.O_CLOEXEC)
            try:
                recovered = runner._read_prearmed_history(
                    lane_fd, history_fd,
                    binding["prearmed_history_markers"][0],
                    registration=fixture["pure"]["preregistration"],
                    logical_plan=fixture["pure"]["plan"],
                    lane_binding=binding,
                    label="zero-byte recovered first history")
            finally:
                os.close(history_fd)
                os.close(lane_fd)
            self.assertEqual(recovered["states"], ["INIT", "PREREGISTERED"])
            bundle = contract.strict_json_loads(
                (lane / runner.ATTEMPTS_DIRECTORY / "attempt-001" /
                 runner.AUTHORITY_BUNDLE_FILE).read_bytes(),
                "zero-byte recovery authority bundle")
            self.assertEqual(
                bundle["budget_commit"]["prospective_ledger"]
                      ["setup_invalid_used"],
                1)

    def test_prearmed_partial_must_be_exact_next_entry_prefix(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        logical_plan = fixture["pure"]["plan"]
        lane_binding = {"schema": "unit-partial-prefix-lane/v1"}
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            binding = create_prearmed_history_slot_for_test(
                lane_fd, "prearmed-0001-history.jsonl")
            marker_fd = runner._open_prearmed_history_for_allocation(
                lane_fd, binding, registration=registration,
                logical_plan=logical_plan, lane_binding=lane_binding)
            try:
                self.assertEqual(os.pwrite(marker_fd, b"!", 0), 1)
                os.fsync(marker_fd)
                with self.assertRaisesRegex(
                        ValueError, "non-canonical partial state entry"):
                    runner._read_prearmed_history(
                        lane_fd, marker_fd, binding,
                        registration=registration,
                        logical_plan=logical_plan,
                        lane_binding=lane_binding,
                        label="corrupt partial history")
            finally:
                os.close(marker_fd)
                os.close(lane_fd)

    def test_prearmed_partial_after_presampling_is_rejected(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        logical_plan = fixture["pure"]["plan"]
        lane_binding = {"schema": "unit-terminal-partial-lane/v1"}
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            lane_fd = os.open(
                lane, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
            binding = create_prearmed_history_slot_for_test(
                lane_fd, "prearmed-0001-history.jsonl")
            marker_fd = runner._open_prearmed_history_for_allocation(
                lane_fd, binding, registration=registration,
                logical_plan=logical_plan, lane_binding=lane_binding)
            try:
                history = None
                for to_state in runner.PREARMED_STATE_TRANSITIONS.values():
                    history = runner._append_prearmed_state(
                        lane_fd, marker_fd, binding, to_state,
                        registration=registration,
                        logical_plan=logical_plan,
                        lane_binding=lane_binding)
                self.assertIsNotNone(history)
                self.assertEqual(history["states"][-1], "PRESAMPLING")
                self.assertEqual(
                    os.pwrite(marker_fd, b"{", len(history["data"])), 1)
                os.fsync(marker_fd)
                with self.assertRaisesRegex(
                        ValueError, "extends terminal pre-ARMED state"):
                    runner._read_prearmed_history(
                        lane_fd, marker_fd, binding,
                        registration=registration,
                        logical_plan=logical_plan,
                        lane_binding=lane_binding,
                        label="terminal partial history")
            finally:
                os.close(marker_fd)
                os.close(lane_fd)

    def test_durable_launch_intent_permanently_consumes_campaign(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            child = campaign.logical_plan["child_plans"][0]
            role = child["implementation"]
            intent = runner._journal_intent_record(
                child=child, selected_pair=campaign.record["selected_pair"],
                artifact_handle_id=campaign._authority_bundle["artifact_bundle"]
                                                     ["roles"][role]["handle_id"],
                evidence_attempt=campaign.record["evidence_attempt"],
                prior_journal_sha256=None, started_monotonic_ns=1,
                launch_context=campaign._authority_bundle["launch_context"])
            campaign._write_journal(
                "intent-0000.json", intent, consume_launch=True)
            campaign.close()
            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "permanently consumes")
                self.assertEqual(controls["qualification"].call_count, 0)
                self.assertEqual(controls["popen"].call_count, 0)
            self.assertEqual(lane.stat().st_mode & 0o777, 0o500)
            self.assertEqual(
                (lane / runner.ATTEMPTS_DIRECTORY).stat().st_mode & 0o777,
                0o500)

    def test_replaced_journal_cannot_hide_a_durable_launch_intent(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                detached = Path(directory) / "detached-journal"

                def replace_journal(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    self.assertGreater(descriptor, 0)
                    self.assertGreater(timeout_seconds, 0)
                    self.assertEqual(
                        expected_launch_context,
                        campaign._authority_bundle["launch_context"])
                    campaign.attempt_path.chmod(0o700)
                    os.rename(
                        campaign.attempt_path / runner.JOURNAL_DIRECTORY,
                        detached)
                    (campaign.attempt_path / runner.JOURNAL_DIRECTORY).mkdir(
                        mode=0o700)
                    campaign.attempt_path.chmod(0o500)
                    return {
                        "outcome": "nonzero", "returncode": 7,
                        "stdout": b"", "stderr": b"injected",
                        "elapsed_ns": 1, "error": "injected child failure",
                    }

                try:
                    with mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=replace_journal):
                        self.assertRejected(campaign.run_all, "journal")
                finally:
                    campaign.close()
                self.assertTrue(any(detached.iterdir()))
                self.assertEqual(list(
                    (lane / "attempts" / "attempt-001" /
                     runner.JOURNAL_DIRECTORY).iterdir()), [])
                controls["qualification"].reset_mock()
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "manifest")
                self.assertEqual(controls["qualification"].call_count, 0)

    def test_deleted_first_intent_cannot_clear_monotone_launch_marker(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))

                def delete_intent_then_crash(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    self.assertGreater(descriptor, 0)
                    self.assertGreater(timeout_seconds, 0)
                    self.assertEqual(
                        expected_launch_context,
                        campaign._authority_bundle["launch_context"])
                    os.unlink("intent-0000.json", dir_fd=campaign._journal_fd)
                    os.fsync(campaign._journal_fd)
                    raise RuntimeError("simulated crash after child launch")

                try:
                    with mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=delete_intent_then_crash):
                        with self.assertRaisesRegex(RuntimeError, "simulated"):
                            campaign.run_all()
                finally:
                    campaign.close()
                attempt = lane / "attempts" / "attempt-001"
                self.assertGreater(
                    (attempt / runner.LAUNCH_CONSUMED_FILE).stat().st_size, 0)
                self.assertEqual(list(
                    (attempt / runner.JOURNAL_DIRECTORY).iterdir()), [])
                controls["qualification"].reset_mock()
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")))
                self.assertEqual(controls["qualification"].call_count, 0)

    def test_post_close_journal_suffix_rewrite_fails_checkpoint_replay(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            child = campaign.logical_plan["child_plans"][0]
            role = child["implementation"]
            intent = runner._journal_intent_record(
                child=child, selected_pair=campaign.record["selected_pair"],
                artifact_handle_id=campaign._authority_bundle["artifact_bundle"]
                                                     ["roles"][role]["handle_id"],
                evidence_attempt=campaign.record["evidence_attempt"],
                prior_journal_sha256=None, started_monotonic_ns=1,
                launch_context=campaign._authority_bundle["launch_context"])
            intent_data = campaign._write_journal(
                "intent-0000.json", intent, consume_launch=True)
            result = runner._journal_result_record(
                child_index=0,
                intent_sha256=hashlib.sha256(intent_data).hexdigest(),
                outcome="nonzero", finished_monotonic_ns=2,
                elapsed_ns=1, returncode=7,
                stdout_sha256=hashlib.sha256(b"").hexdigest(),
                stderr_sha256=hashlib.sha256(b"failed").hexdigest(),
                payload=None, normalized=None, error="unit failure")
            campaign._write_journal("result-0000.json", result)
            campaign.close()

            result_path = (
                lane / runner.ATTEMPTS_DIRECTORY / "attempt-001" /
                runner.JOURNAL_DIRECTORY / "result-0000.json")
            forged = copy.deepcopy(result)
            forged["finished_monotonic_ns"] += 1
            forged["elapsed_ns"] += 1
            result_path.unlink()
            result_path.write_bytes(contract.canonical_json_bytes(forged))
            result_path.chmod(0o400)

            with patched_live_acquisition(fixture) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "checkpoint")
                self.assertEqual(controls["qualification"].call_count, 0)

    def test_forked_close_cannot_unlock_parent_authority(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            canonical_path = str(Path(directory) / "canonical.lock")
            canonical_lock = runner.baseline_acquire.CanonicalFileLock(
                canonical_path)
            canonical_lock.__enter__()
            campaign._canonical_lock = canonical_lock
            runner._register_active_campaign(campaign)
            child = os.fork()
            if child == 0:
                try:
                    campaign.close()
                except runner.ArmingError:
                    os._exit(0)
                except BaseException:
                    os._exit(2)
                os._exit(1)
            waited, status = os.waitpid(child, 0)
            self.assertEqual(waited, child)
            self.assertTrue(os.WIFEXITED(status))
            self.assertEqual(os.WEXITSTATUS(status), 0)

            lane_contender = os.open(
                lane / ".arming.lock", os.O_RDWR | os.O_CLOEXEC)
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        lane_contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                os.close(lane_contender)
            canonical_contender = runner.baseline_acquire.CanonicalFileLock(
                canonical_path, blocking=False)
            with self.assertRaises(Exception):
                canonical_contender.__enter__()
            with patched_live_acquisition(fixture), mock.patch.object(
                    runner, "revalidate_sealed_descriptor",
                    return_value={}), mock.patch.object(
                    runner, "_capture_live_host_instance",
                    return_value=fixture["pure"]["host"]):
                campaign._revalidate_authority()
            campaign.close()

            lane_contender = os.open(
                lane / ".arming.lock", os.O_RDWR | os.O_CLOEXEC)
            try:
                fcntl.flock(
                    lane_contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(lane_contender, fcntl.LOCK_UN)
            finally:
                os.close(lane_contender)
            canonical_contender = runner.baseline_acquire.CanonicalFileLock(
                canonical_path, blocking=False)
            canonical_contender.__enter__()
            canonical_contender.__exit__(None, None, None)

    def test_uncontrolled_fork_revokes_worker_owned_campaign_before_child_code(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            fixture = live_fixture(output_lane=lane)
            registration = fixture["pure"]["preregistration"]
            canonical_path = str(Path(directory) / "canonical.lock")
            ready = threading.Event()
            release_owner = threading.Event()
            owner_finished = threading.Event()
            owner_errors: list[BaseException] = []
            descriptor_numbers: list[int] = []

            def own_campaign_without_exporting_object():
                try:
                    with patched_live_acquisition(
                            fixture,
                            enforce_output_lane_binding=True):
                        campaign = runner.acquire_and_arm(
                            lane=lane,
                            preregistration_bytes=contract.canonical_json_bytes(
                                registration),
                            candidate_authority_lane=Path(
                                "/candidate-authority"),
                            exact_main_authority_lane=Path(
                                "/exact-main-authority"))
                    canonical_lock = runner.baseline_acquire.CanonicalFileLock(
                        canonical_path)
                    canonical_lock.__enter__()
                    campaign._canonical_lock = canonical_lock
                    runtime_guard = runner._RetainedRuntimeGuard()
                    runtime_guard._open_inotify()
                    runtime_descriptor = os.open(
                        __file__, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
                    runtime_guard._files.append({
                        "descriptor": runtime_descriptor,
                    })
                    campaign._runtime_guard = runtime_guard
                    runner._register_active_campaign(campaign)
                    descriptor_numbers.extend(
                        getattr(campaign, name) for name in (
                            "_main_fd", "_candidate_fd", "_manifest_fd",
                            "_bundle_fd", "_armed_fd", "_launch_marker_fd",
                            "_journal_fd", "_campaign_marker_fd",
                            "_campaign_progress_fd",
                            "_attempt_fd", "_attempts_fd",
                            "_lane_binding_fd", "_lane_lock_fd", "_lane_fd"))
                    descriptor_numbers.extend((
                        runtime_descriptor, runtime_guard._inotify_fd,
                        canonical_lock.descriptor))
                    ready.set()
                    if not release_owner.wait(timeout=10):
                        raise AssertionError(
                            "timed out retaining worker campaign")
                    campaign.close()
                except BaseException as error:
                    owner_errors.append(error)
                    ready.set()
                finally:
                    owner_finished.set()

            owner = threading.Thread(
                target=own_campaign_without_exporting_object)
            owner.start()
            try:
                self.assertTrue(ready.wait(timeout=30))
                self.assertEqual(owner_errors, [])
                self.assertEqual(len(descriptor_numbers), 17)
                child = os.fork()
                if child == 0:
                    # No campaign method is called here.  The registered
                    # at-fork callback must already have revoked everything.
                    if runner._ACTIVE_CAMPAIGNS:
                        os._exit(2)
                    for descriptor in descriptor_numbers:
                        try:
                            os.fstat(descriptor)
                        except OSError as error:
                            if error.errno == errno.EBADF:
                                continue
                            os._exit(3)
                        os._exit(4)
                    os._exit(0)
                waited, status = os.waitpid(child, 0)
                self.assertEqual(waited, child)
                self.assertTrue(os.WIFEXITED(status))
                self.assertEqual(os.WEXITSTATUS(status), 0)
                for descriptor in descriptor_numbers:
                    os.fstat(descriptor)

                lane_contender = os.open(
                    lane / ".arming.lock", os.O_RDWR | os.O_CLOEXEC)
                try:
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(
                            lane_contender,
                            fcntl.LOCK_EX | fcntl.LOCK_NB)
                finally:
                    os.close(lane_contender)
                canonical_contender = \
                    runner.baseline_acquire.CanonicalFileLock(
                        canonical_path, blocking=False)
                with self.assertRaises(Exception):
                    canonical_contender.__enter__()
            finally:
                release_owner.set()
                owner.join(timeout=5)
            self.assertFalse(owner.is_alive())
            self.assertTrue(owner_finished.is_set())
            self.assertEqual(owner_errors, [])
            self.assertEqual(runner._ACTIVE_CAMPAIGNS, {})

    def test_exact_scope_handoff_registers_campaign_before_return(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            self.assertNotIn(id(campaign), runner._ACTIVE_CAMPAIGNS)

            def transfer_campaign():
                return campaign

            wrapped_transfer = runner._with_exact_source_imports(
                transfer_campaign)
            returned = wrapped_transfer()
            self.assertIs(returned, campaign)
            self.assertIs(
                runner._ACTIVE_CAMPAIGNS.get(id(campaign)), campaign)
            campaign.close()
            self.assertNotIn(id(campaign), runner._ACTIVE_CAMPAIGNS)

    def test_registered_campaign_popen_preserves_only_selected_exec_handoff(
            self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            runner._register_active_campaign(campaign)
            descriptor = campaign._candidate_fd
            try:
                result = runner._run_child_process(
                    [
                        "/usr/bin/python3", "-I", "-S", "-B", "-c",
                        f"import os; os.fstat({descriptor}); print('live')",
                    ],
                    descriptor=descriptor, timeout_seconds=5,
                    expected_launch_context=launch_context())
                self.assertEqual(result["outcome"], "accepted")
                self.assertEqual(result["stdout"], b"live\n")
                os.fstat(descriptor)
                self.assertIs(
                    runner._ACTIVE_CAMPAIGNS.get(id(campaign)), campaign)
                self.assertIsNone(runner._CONTROLLED_FORK_PERMIT[0])
            finally:
                campaign.close()

    def test_fork_child_cleanup_failure_is_fail_stop(self) -> None:
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))

            class FailingRuntimeGuard:
                def __init__(self):
                    self.fail = True

                def close(self):
                    if self.fail:
                        raise OSError(errno.EIO, "injected cleanup failure")

            guard = FailingRuntimeGuard()
            campaign._runtime_guard = guard
            runner._register_active_campaign(campaign)
            try:
                child = os.fork()
                if child == 0:
                    os._exit(99)
                waited, status = os.waitpid(child, 0)
                self.assertEqual(waited, child)
                self.assertTrue(os.WIFEXITED(status))
                self.assertEqual(os.WEXITSTATUS(status), 125)
                self.assertIs(
                    runner._ACTIVE_CAMPAIGNS.get(id(campaign)), campaign)
            finally:
                guard.fail = False
                campaign.close()

    def test_forked_methods_reject_before_inherited_python_locks(self) -> None:
        for method_name in ("close", "run_all"):
            with self.subTest(method=method_name), \
                    tempfile.TemporaryDirectory() as directory:
                lane = private_lane(directory)
                fixture = live_fixture()
                registration = fixture["pure"]["preregistration"]
                with patched_live_acquisition(fixture):
                    campaign = runner.acquire_and_arm(
                        lane=lane,
                        preregistration_bytes=contract.canonical_json_bytes(
                            registration),
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path(
                            "/exact-main-authority"))
                canonical_path = str(Path(directory) / "canonical.lock")
                canonical_lock = runner.baseline_acquire.CanonicalFileLock(
                    canonical_path)
                canonical_lock.__enter__()
                campaign._canonical_lock = canonical_lock
                runner._register_active_campaign(campaign)
                descriptors = [getattr(campaign, name) for name in (
                    "_main_fd", "_candidate_fd", "_manifest_fd",
                    "_bundle_fd", "_armed_fd", "_launch_marker_fd",
                    "_journal_fd", "_campaign_marker_fd", "_attempt_fd",
                    "_attempts_fd", "_lane_binding_fd", "_lane_lock_fd",
                    "_lane_fd")]
                locks_owned = threading.Event()
                release_locks = threading.Event()

                def hold_parent_locks():
                    runner._EXACT_IMPORT_LOCK.acquire()
                    campaign._operation_lock.acquire()
                    campaign._operation_owner = threading.get_ident()
                    locks_owned.set()
                    try:
                        if not release_locks.wait(timeout=10):
                            raise AssertionError(
                                "timed out holding parent campaign locks")
                    finally:
                        campaign._operation_owner = None
                        campaign._operation_lock.release()
                        runner._EXACT_IMPORT_LOCK.release()

                worker = threading.Thread(target=hold_parent_locks)
                worker.start()
                try:
                    self.assertTrue(locks_owned.wait(timeout=5))
                    child = os.fork()
                    if child == 0:
                        signal.signal(
                            signal.SIGALRM,
                            lambda unused_signal, unused_frame: os._exit(42))
                        signal.alarm(3)
                        try:
                            getattr(campaign, method_name)()
                        except runner.ArmingError:
                            for descriptor in descriptors:
                                try:
                                    os.fstat(descriptor)
                                except OSError:
                                    continue
                                os._exit(3)
                            os._exit(0)
                        except BaseException:
                            os._exit(2)
                        os._exit(1)
                    waited, status = os.waitpid(child, 0)
                    self.assertEqual(waited, child)
                    self.assertTrue(os.WIFEXITED(status))
                    self.assertEqual(os.WEXITSTATUS(status), 0)
                    self.assertTrue(worker.is_alive())
                    contender = runner.baseline_acquire.CanonicalFileLock(
                        canonical_path, blocking=False)
                    with self.assertRaises(Exception):
                        contender.__enter__()
                finally:
                    release_locks.set()
                    worker.join(timeout=5)
                    self.assertFalse(worker.is_alive())
                    campaign.close()

    def test_fork_during_parent_close_never_waits_on_inherited_lock(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            canonical_path = str(Path(directory) / "canonical.lock")
            canonical_lock = runner.baseline_acquire.CanonicalFileLock(
                canonical_path)
            canonical_lock.__enter__()
            close_reached_canonical = threading.Event()
            permit_parent_close = threading.Event()

            class BlockingCanonicalExit:
                def __init__(self, delegate):
                    self.delegate = delegate
                    self.descriptor = delegate.descriptor
                    self.identity = delegate.identity

                def validate_current(self):
                    return self.delegate.validate_current()

                def __exit__(self, kind, value, traceback):
                    close_reached_canonical.set()
                    if not permit_parent_close.wait(timeout=10):
                        raise AssertionError(
                            "timed out blocking parent canonical close")
                    result = self.delegate.__exit__(kind, value, traceback)
                    self.descriptor = self.delegate.descriptor
                    self.identity = self.delegate.identity
                    return result

            campaign._canonical_lock = BlockingCanonicalExit(canonical_lock)
            runner._register_active_campaign(campaign)
            parent_close_errors: list[BaseException] = []
            fork_errors: list[BaseException] = []
            child_statuses: list[int] = []
            fork_attempted = threading.Event()
            fork_finished = threading.Event()

            def close_in_parent():
                try:
                    campaign.close()
                except BaseException as error:
                    parent_close_errors.append(error)

            def fork_while_parent_is_closing():
                fork_attempted.set()
                try:
                    child = os.fork()
                    if child == 0:
                        signal.signal(
                            signal.SIGALRM,
                            lambda unused_signal, unused_frame: os._exit(42))
                        signal.alarm(3)
                        try:
                            campaign.close()
                        except runner.ArmingError:
                            descriptor = getattr(
                                campaign._canonical_lock, "descriptor", -1)
                            if descriptor >= 0:
                                try:
                                    os.fstat(descriptor)
                                except OSError:
                                    pass
                                else:
                                    os._exit(3)
                            os._exit(0)
                        except BaseException:
                            os._exit(2)
                        os._exit(1)
                    waited, status = os.waitpid(child, 0)
                    if waited != child:
                        raise AssertionError("waited for the wrong fork child")
                    child_statuses.append(status)
                except BaseException as error:
                    fork_errors.append(error)
                finally:
                    fork_finished.set()

            close_thread = threading.Thread(target=close_in_parent)
            close_thread.start()
            fork_thread: threading.Thread | None = None
            try:
                self.assertTrue(close_reached_canonical.wait(timeout=5))
                self.assertTrue(campaign._closed)
                self.assertFalse(campaign._cleanup_complete)
                fork_thread = threading.Thread(
                    target=fork_while_parent_is_closing)
                fork_thread.start()
                self.assertTrue(fork_attempted.wait(timeout=5))
                fork_thread.join(timeout=0.2)
                self.assertTrue(fork_thread.is_alive())
                self.assertFalse(fork_finished.is_set())
                contender = runner.baseline_acquire.CanonicalFileLock(
                    canonical_path, blocking=False)
                with self.assertRaises(Exception):
                    contender.__enter__()
            finally:
                permit_parent_close.set()
                close_thread.join(timeout=5)
                if fork_thread is not None:
                    fork_thread.join(timeout=5)
            self.assertFalse(close_thread.is_alive())
            self.assertIsNotNone(fork_thread)
            self.assertFalse(fork_thread.is_alive())
            self.assertEqual(parent_close_errors, [])
            self.assertEqual(fork_errors, [])
            self.assertEqual(len(child_statuses), 1)
            self.assertTrue(os.WIFEXITED(child_statuses[0]))
            self.assertEqual(os.WEXITSTATUS(child_statuses[0]), 0)
            self.assertTrue(campaign._cleanup_complete)

    def test_whole_attempt_rename_at_child_boundary_fails_closed(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                attempts = lane / runner.ATTEMPTS_DIRECTORY
                staging = attempts / (
                    f".attempt-001.staging-{os.getpid()}-" + "a" * 32)

                def rename_attempt(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    self.assertGreater(descriptor, 0)
                    attempts.chmod(0o700)
                    os.rename(campaign.attempt_path, staging)
                    attempts.chmod(0o500)
                    return {
                        "outcome": "nonzero", "returncode": 7,
                        "stdout": b"", "stderr": b"renamed",
                        "elapsed_ns": 1, "error": "renamed attempt",
                    }

                try:
                    with mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=rename_attempt):
                        self.assertRejected(campaign.run_all, "attempts")
                finally:
                    campaign.close()
                controls["qualification"].reset_mock()
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "stale unpublished journal is not empty")
                self.assertEqual(controls["qualification"].call_count, 0)
                self.assertTrue(staging.exists())
            self.assertEqual(attempts.stat().st_mode & 0o777, 0o500)

    def test_published_unrun_attempt_rename_is_recovered_without_reset(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            frozen_pair = copy.deepcopy(campaign.record["selected_pair"])
            campaign.close()
            attempts = lane / runner.ATTEMPTS_DIRECTORY
            staging = attempts / (
                f".attempt-001.staging-{os.getpid()}-" + "b" * 32)
            attempts.chmod(0o700)
            os.rename(attempts / "attempt-001", staging)
            attempts.chmod(0o500)
            retry_fixture = live_fixture(expected_frozen_pair=frozen_pair)
            with patched_live_acquisition(retry_fixture) as controls:
                retry = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            try:
                self.assertEqual(retry.record["evidence_attempt"], 2)
                self.assertEqual(controls["qualification"].call_count, 1)
            finally:
                retry.close()
            self.assertFalse(staging.exists())
            self.assertTrue((attempts / "attempt-001").is_dir())
            self.assertEqual(attempts.stat().st_mode & 0o777, 0o500)

    def test_attempt_entry_symlink_substitution_fails_closed(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            detached = lane / runner.ATTEMPTS_DIRECTORY / \
                "detached-attempt-001"
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                attempts = lane / runner.ATTEMPTS_DIRECTORY
                attempts.chmod(0o700)
                os.rename(campaign.attempt_path, detached)
                campaign.attempt_path.symlink_to(
                    detached, target_is_directory=True)
                attempts.chmod(0o500)
                try:
                    self.assertRejected(
                        campaign._revalidate_authority, "attempt inode")
                finally:
                    campaign.close()
                controls["qualification"].reset_mock()
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")))
                self.assertEqual(controls["qualification"].call_count, 0)
            attempts.chmod(0o700)
            campaign.attempt_path.unlink()
            os.rename(detached, campaign.attempt_path)
            attempts.chmod(0o500)

    def test_attempts_namespace_symlink_at_child_boundary_fails_closed(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            detached = Path(directory) / "detached-attempts"
            with patched_live_acquisition(fixture) as controls:
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
                attempts = lane / runner.ATTEMPTS_DIRECTORY

                def replace_attempts(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    self.assertGreater(descriptor, 0)
                    lane.chmod(0o700)
                    attempts.chmod(0o700)
                    os.rename(attempts, detached)
                    detached.chmod(0o500)
                    attempts.symlink_to(detached, target_is_directory=True)
                    lane.chmod(0o500)
                    return {
                        "outcome": "nonzero", "returncode": 7,
                        "stdout": b"", "stderr": b"replaced",
                        "elapsed_ns": 1, "error": "replaced namespace",
                    }

                try:
                    with mock.patch.object(
                        runner, "_run_child_process",
                            side_effect=replace_attempts):
                        self.assertRejected(campaign.run_all)
                finally:
                    campaign.close()
                controls["qualification"].reset_mock()
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")))
                self.assertEqual(controls["qualification"].call_count, 0)
            lane.chmod(0o700)
            attempts.unlink()
            detached.chmod(0o700)
            os.rename(detached, attempts)
            attempts.chmod(0o500)
            lane.chmod(0o500)

    def test_retry_host_drift_rejects_before_qualification_scan(self) -> None:
        first_fixture = live_fixture()
        registration = first_fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(first_fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            selected_pair = campaign.record["selected_pair"]
            campaign.close()
            retry_fixture = live_fixture(selected_pair)
            drifted = dict(retry_fixture["pure"]["host"])
            drifted["kernel_release"] += "-changed"
            with patched_live_acquisition(retry_fixture) as controls, \
                    mock.patch.object(
                        runner, "_capture_live_host_instance",
                        return_value=drifted):
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "frozen host")
                self.assertEqual(controls["qualification"].call_count, 0)

    def test_exhausted_evidence_budget_rejects_before_qualification_scan(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        preregistration_bytes = contract.canonical_json_bytes(registration)
        expected_pair = None
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            for unused_attempt in range(
                    registration["budgets"]["evidence_attempts"]):
                fixture = live_fixture(expected_pair)
                with patched_live_acquisition(fixture):
                    campaign = runner.acquire_and_arm(
                        lane=lane, preregistration_bytes=preregistration_bytes,
                        candidate_authority_lane=Path("/candidate-authority"),
                        exact_main_authority_lane=Path("/exact-main-authority"))
                expected_pair = campaign.record["selected_pair"]
                campaign.close()
            exhausted = live_fixture(expected_pair)
            with patched_live_acquisition(exhausted) as controls:
                self.assertRejected(lambda: runner.acquire_and_arm(
                    lane=lane, preregistration_bytes=preregistration_bytes,
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority")),
                    "budget is exhausted")
                self.assertEqual(controls["qualification"].call_count, 0)

    def test_safe_staging_cleanup_and_unsafe_orphan_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            attempts = Path(directory).resolve() / "attempts"
            attempts.mkdir(mode=0o700)
            attempts_fd = os.open(attempts, os.O_RDONLY | os.O_DIRECTORY)
            try:
                safe = attempts / (".attempt-001.staging-1-" + "0" * 32)
                safe.mkdir(mode=0o700)
                partial = safe / runner.AUTHORITY_BUNDLE_FILE
                partial.write_bytes(b"partial")
                partial.chmod(0o600)
                (safe / runner.JOURNAL_DIRECTORY).mkdir(mode=0o700)
                safe.chmod(0o500)
                runner._clean_staging_directories(attempts_fd)
                self.assertFalse(safe.exists())
                unsafe = attempts / (".attempt-001.staging-1-" + "1" * 32)
                unsafe.mkdir(mode=0o700)
                retained = unsafe / runner.ARMED_FILE
                retained.write_bytes(b"valid-shape")
                retained.chmod(0o400)
                (unsafe / "unexpected").write_bytes(b"x")
                self.assertRejected(
                    lambda: runner._clean_staging_directories(attempts_fd),
                    "unsafe entry")
                self.assertTrue(unsafe.exists())
                self.assertTrue(retained.exists())
            finally:
                os.close(attempts_fd)

    def test_run_all_is_exactly_once_and_sequential_without_launcher_knobs(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))
            commands = []
            journal_names = []
            journal_payloads = []

            def child_canary(
                command, *, descriptor, timeout_seconds,
                expected_launch_context,
            ):
                commands.append((list(command), descriptor, timeout_seconds))
                self.assertEqual(
                    expected_launch_context,
                    campaign._authority_bundle["launch_context"])
                return {
                    "outcome": "accepted", "returncode": 0,
                    "stdout": b"{}", "stderr": b"", "elapsed_ns": 1,
                    "error": None,
                }

            def journal_canary(name, value, **unused_options):
                journal_names.append(name)
                data = contract.canonical_json_bytes(value)
                journal_payloads.append(data)
                return data

            try:
                with patched_live_acquisition(fixture), mock.patch.object(
                        runner.ArmedCampaign,
                        "_revalidate_authority", return_value=None), \
                        mock.patch.object(
                            runner.ArmedCampaign, "_write_journal",
                            side_effect=journal_canary), \
                        mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=child_canary), \
                        mock.patch.object(
                            runner, "_normalize_gen3_result",
                            return_value={"validated": True}), \
                        mock.patch.object(
                            runner, "_process_start_ticks", return_value=1):
                    result = campaign.run_all()
                children = campaign.logical_plan["child_plans"]
                self.assertEqual(result["child_count"], 1650)
                self.assertEqual(len(commands), 1650)
                self.assertEqual(len(journal_names), 3301)
                self.assertEqual(journal_names[0], "intent-0000.json")
                self.assertEqual(journal_names[-1], "complete.json")
                self.assertEqual(
                    result["journal_chain_sha256"],
                    hashlib.sha256(journal_payloads[-1]).hexdigest())
                for index, (command, descriptor, timeout) in enumerate(commands):
                    child = children[index]
                    self.assertEqual(command, runner._actual_child_command(
                        child, campaign.record["selected_pair"], descriptor))
                    self.assertEqual(timeout,
                                     child["timeout_budget"]["timeout_seconds"])
                self.assertRejected(campaign.run_all, "exactly once")
            finally:
                campaign.close()

    def test_close_serializes_with_the_entire_run_operation(self):
        fixture = live_fixture()
        registration = fixture["pure"]["preregistration"]
        with tempfile.TemporaryDirectory() as directory:
            lane = private_lane(directory)
            with patched_live_acquisition(fixture):
                campaign = runner.acquire_and_arm(
                    lane=lane,
                    preregistration_bytes=contract.canonical_json_bytes(
                        registration),
                    candidate_authority_lane=Path("/candidate-authority"),
                    exact_main_authority_lane=Path("/exact-main-authority"))

                selected_calls = 0
                final_selection_reached = threading.Event()
                permit_launch = threading.Event()
                close_finished = threading.Event()
                child_observations: list[tuple[bool, bool]] = []
                run_errors: list[BaseException] = []
                close_errors: list[BaseException] = []
                original_selected = \
                    runner.ArmedCampaign._selected_executable_capability

                def blocking_selected(campaign_instance, role):
                    nonlocal selected_calls
                    result = original_selected(campaign_instance, role)
                    selected_calls += 1
                    if selected_calls == 2:
                        final_selection_reached.set()
                        if not permit_launch.wait(timeout=5):
                            raise AssertionError(
                                "timed out awaiting launch permission")
                    return result

                def rejected_child(
                    unused_command, *, descriptor, timeout_seconds,
                    expected_launch_context,
                ):
                    del timeout_seconds, expected_launch_context
                    child_observations.append((
                        close_finished.is_set(),
                        os.fstat(descriptor).st_ino > 0))
                    return {
                        "outcome": "output-rejected", "returncode": None,
                        "stdout": b"", "stderr": b"", "elapsed_ns": 1,
                        "error": "injected one-child stop",
                    }

                def journal_canary(unused_name, value, **unused_options):
                    return contract.canonical_json_bytes(value)

                def run_campaign():
                    try:
                        campaign.run_all()
                    except BaseException as error:
                        run_errors.append(error)

                def close_campaign():
                    try:
                        campaign.close()
                    except BaseException as error:
                        close_errors.append(error)
                    finally:
                        close_finished.set()

                with mock.patch.object(
                        runner.ArmedCampaign,
                        "_selected_executable_capability",
                        new=blocking_selected), mock.patch.object(
                            runner.ArmedCampaign,
                            "_revalidate_authority", return_value=None), \
                        mock.patch.object(
                            runner.ArmedCampaign, "_write_journal",
                            side_effect=journal_canary), mock.patch.object(
                            runner, "_run_child_process",
                            side_effect=rejected_child):
                    run_thread = threading.Thread(target=run_campaign)
                    run_thread.start()
                    self.assertTrue(final_selection_reached.wait(timeout=5))
                    close_thread = threading.Thread(target=close_campaign)
                    close_thread.start()
                    close_thread.join(timeout=0.2)
                    self.assertTrue(close_thread.is_alive())
                    self.assertFalse(close_finished.is_set())
                    permit_launch.set()
                    run_thread.join(timeout=5)
                    close_thread.join(timeout=5)
                    self.assertFalse(run_thread.is_alive())
                    self.assertFalse(close_thread.is_alive())

                self.assertEqual(close_errors, [])
                self.assertEqual(len(run_errors), 1)
                self.assertIn("injected one-child stop", str(run_errors[0]))
                self.assertEqual(child_observations, [(False, True)])
                self.assertTrue(close_finished.is_set())
                self.assertTrue(campaign._closed)

    def test_private_process_supervisor_fixes_environment_and_bounds_output(self):
        descriptor = os.open("/dev/null", os.O_RDONLY)
        real_popen = runner.subprocess.Popen
        real_kill = runner._kill_process_group
        calls = []
        kill_returncodes = []

        def recording_popen(*arguments, **keywords):
            calls.append((arguments, keywords))
            return real_popen(*arguments, **keywords)

        def recording_kill(process):
            kill_returncodes.append(process.returncode)
            return real_kill(process)

        try:
            with mock.patch.object(
                    runner.subprocess, "Popen", side_effect=recording_popen), \
                    mock.patch.object(
                        runner, "_kill_process_group",
                        side_effect=recording_kill):
                result = runner._run_child_process(
                    ["/usr/bin/printf", "{}"], descriptor=descriptor,
                    timeout_seconds=5,
                    expected_launch_context=launch_context())
            self.assertEqual(result["outcome"], "accepted")
            self.assertEqual(result["stdout"], b"{}")
            keywords = calls[0][1]
            self.assertEqual(keywords["env"], runner.CHILD_ENVIRONMENT)
            self.assertEqual(keywords["cwd"], "/")
            self.assertEqual(keywords["stdin"], runner.subprocess.DEVNULL)
            self.assertEqual(keywords["pass_fds"], (descriptor,))
            self.assertTrue(keywords["close_fds"])
            self.assertTrue(keywords["start_new_session"])
            self.assertEqual(keywords["preexec_fn"].__self__.__class__,
                             runner._ChildTreeScope)
            self.assertTrue(kill_returncodes)
            self.assertTrue(all(value is None for value in kill_returncodes))

            previous_mask = signal.pthread_sigmask(
                signal.SIG_UNBLOCK, {signal.SIGUSR1})
            try:
                child_mask = runner._run_child_process(
                    ["/usr/bin/python3", "-c",
                     "import json,signal; "
                    "m=signal.pthread_sigmask(signal.SIG_BLOCK,set()); "
                    "print(json.dumps({'blocked': signal.SIGUSR1 in m}))"],
                    descriptor=descriptor, timeout_seconds=5,
                    expected_launch_context=launch_context())
            finally:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
            self.assertEqual(child_mask["outcome"], "accepted")
            self.assertEqual(
                contract.strict_json_loads(child_mask["stdout"], "child mask"),
                {"blocked": False})

            previous_mask = signal.pthread_sigmask(
                signal.SIG_BLOCK, {signal.SIGUSR2})
            previous_disposition = signal.getsignal(signal.SIGUSR2)
            try:
                signal.signal(signal.SIGUSR2, signal.SIG_IGN)
                canonical_signals = runner._run_child_process(
                    ["/usr/bin/python3", "-c",
                     "import json,signal; "
                     "m=signal.pthread_sigmask(signal.SIG_BLOCK,set()); "
                     "print(json.dumps({"
                     "'blocked':signal.SIGUSR2 in m,"
                     "'default':signal.getsignal(signal.SIGUSR2) "
                     "is signal.SIG_DFL}))"],
                    descriptor=descriptor, timeout_seconds=5,
                    expected_launch_context=launch_context())
                parent_mask = signal.pthread_sigmask(
                    signal.SIG_BLOCK, set())
                self.assertIn(signal.SIGUSR2, parent_mask)
                self.assertIs(
                    signal.getsignal(signal.SIGUSR2), signal.SIG_IGN)
            finally:
                signal.signal(signal.SIGUSR2, previous_disposition)
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
            self.assertEqual(canonical_signals["outcome"], "accepted")
            self.assertEqual(contract.strict_json_loads(
                canonical_signals["stdout"], "canonical child signals"), {
                    "blocked": False, "default": True,
                })

            benchmark_cpu = min(os.sched_getaffinity(0))
            process_context = runner._run_child_process(
                ["/usr/bin/taskset", "-c", str(benchmark_cpu),
                 "/usr/bin/python3", "-c",
                 "import json,os,resource,signal; "
                 "u=os.umask(0o22); os.umask(u); "
                 "m=signal.pthread_sigmask(signal.SIG_BLOCK,set()); "
                 "print(json.dumps({"
                 "'affinity':sorted(os.sched_getaffinity(0)),"
                 "'as':list(resource.getrlimit(resource.RLIMIT_AS)),"
                 "'nice':os.getpriority(os.PRIO_PROCESS,0),"
                 "'policy':os.sched_getscheduler(0),"
                 "'priority':os.sched_getparam(0).sched_priority,"
                 "'signals_empty':not m,'umask':u}))"],
                descriptor=descriptor, timeout_seconds=5,
                expected_launch_context=launch_context())
            self.assertEqual(process_context["outcome"], "accepted")
            self.assertEqual(contract.strict_json_loads(
                process_context["stdout"], "canonical process context"), {
                    "affinity": [benchmark_cpu],
                    "as": [201_326_592, 201_326_592],
                    "nice": 0, "policy": os.SCHED_OTHER, "priority": 0,
                    "signals_empty": True, "umask": 0o022,
                })

            with mock.patch.object(runner, "MAX_CHILD_STDOUT_BYTES", 32):
                bounded = runner._run_child_process(
                    ["/usr/bin/python3", "-c",
                     "import os; os.write(1, b'x' * 4096)"],
                    descriptor=descriptor, timeout_seconds=5,
                    expected_launch_context=launch_context())
            self.assertEqual(bounded["outcome"], "output-limit")
            self.assertLessEqual(len(bounded["stdout"]), 32)

            descendant = runner._run_child_process(
                ["/usr/bin/python3", "-c",
                 "import os,time; p=os.fork(); "
                 "time.sleep(30) if p == 0 else os._exit(0)"],
                descriptor=descriptor, timeout_seconds=1,
                expected_launch_context=launch_context())
            self.assertEqual(descendant["outcome"], "output-rejected")
            self.assertIn("descendant", descendant["error"])

            escaped = runner._run_child_process(
                ["/usr/bin/python3", "-c",
                 "import os,time; p=os.fork(); "
                 "(os.setsid(),os.close(1),os.close(2),time.sleep(30)) "
                 "if p == 0 else os._exit(0)"],
                descriptor=descriptor, timeout_seconds=1,
                expected_launch_context=launch_context())
            self.assertEqual(escaped["outcome"], "output-rejected")
            self.assertIn("descendant", escaped["error"])
            self.assertEqual(runner._direct_child_pids(), [])
        finally:
            os.close(descriptor)

    def test_subreaper_scope_restores_state_on_every_failure_boundary(self):
        initial = runner._child_subreaper_state()
        with mock.patch.object(
                runner, "_direct_child_pids", side_effect=[[], [12345]]):
            self.assertRejected(
                lambda: runner._ChildTreeScope(
                    launch_context()).__enter__(), "raced")
        self.assertEqual(runner._child_subreaper_state(), initial)

        descriptor = os.open("/dev/null", os.O_RDONLY)
        try:
            with mock.patch.object(
                    runner.subprocess, "Popen",
                    side_effect=RuntimeError("injected Popen failure")):
                with self.assertRaisesRegex(RuntimeError, "injected"):
                    runner._run_child_process(
                        ["/usr/bin/printf", "{}"], descriptor=descriptor,
                        timeout_seconds=1,
                        expected_launch_context=launch_context())
            self.assertEqual(runner._child_subreaper_state(), initial)

            with mock.patch.object(
                    runner.selectors, "DefaultSelector",
                    side_effect=OSError(errno.EMFILE, "injected selector failure")):
                with self.assertRaisesRegex(OSError, "selector failure"):
                    runner._run_child_process(
                        ["/usr/bin/python3", "-c", "import time; time.sleep(30)"],
                        descriptor=descriptor, timeout_seconds=1,
                        expected_launch_context=launch_context())
            self.assertEqual(runner._child_subreaper_state(), initial)
            self.assertEqual(runner._direct_child_pids(), [])

            scope = runner._ChildTreeScope(launch_context()).__enter__()
            with mock.patch.object(
                    runner, "_direct_child_pids", return_value=[12345]):
                self.assertRejected(scope.close, "cleanup failed")
            self.assertEqual(runner._child_subreaper_state(), initial)
        finally:
            os.close(descriptor)

    def test_process_supervisor_rejects_nondefault_sigchld_disposition(self):
        descriptor = os.open("/dev/null", os.O_RDONLY)
        original = signal.getsignal(signal.SIGCHLD)
        try:
            for disposition in (signal.SIG_IGN, lambda unused_signum,
                                 unused_frame: None):
                signal.signal(signal.SIGCHLD, disposition)
                self.assertRejected(
                    lambda: runner._run_child_process(
                        ["/usr/bin/python3", "-c", "raise SystemExit(7)"],
                        descriptor=descriptor, timeout_seconds=1,
                        expected_launch_context=launch_context()),
                    "SIGCHLD")
                self.assertEqual(runner._direct_child_pids(), [])
        finally:
            signal.signal(signal.SIGCHLD, original)
            os.close(descriptor)

    def test_process_supervisor_rejects_midrun_sigchld_flip(self) -> None:
        descriptor = os.open("/dev/null", os.O_RDONLY)
        original_chld = signal.getsignal(signal.SIGCHLD)
        original_usr1 = signal.getsignal(signal.SIGUSR1)

        def flip_sigchld(unused_signum, unused_frame):
            signal.signal(signal.SIGCHLD, signal.SIG_IGN)

        try:
            signal.signal(signal.SIGUSR1, flip_sigchld)
            with self.assertRaisesRegex(ValueError, "signal authority"):
                runner._run_child_process(
                    ["/usr/bin/python3", "-c",
                     "import os,signal,sys; "
                    "os.kill(os.getppid(), signal.SIGUSR1); "
                    "os.write(1,b'{}'); sys.exit(7)"],
                    descriptor=descriptor, timeout_seconds=5,
                    expected_launch_context=launch_context())
            self.assertIs(signal.getsignal(signal.SIGCHLD), signal.SIG_DFL)
            self.assertEqual(runner._direct_child_pids(), [])
        finally:
            signal.signal(signal.SIGUSR1, original_usr1)
            signal.signal(signal.SIGCHLD, original_chld)
            os.close(descriptor)


if __name__ == "__main__":
    unittest.main()
