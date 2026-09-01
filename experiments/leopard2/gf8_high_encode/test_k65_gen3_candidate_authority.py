#!/usr/bin/env python3
"""Adversarial tests for the build-only K65 Gen3 authority producer.

The suite uses synthetic payloads and scripted child results.  It launches no
compiler, benchmark, timing runner, qualification scan, or campaign child.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k65_gen3_candidate_authority as authority  # noqa: E402
import k65_gen3_candidate_authority_acquire as acquire  # noqa: E402


contract = authority.contract


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def synthetic_elf(size: int = 4096) -> bytes:
    return b"\x7fELF" + b"C" * (size - 4)


def synthetic_discovery_rows() -> tuple[dict, ...]:
    return authority.normalize_runtime_discovery(
        b"linux-vdso.so.1 (0x00007fff00000000)\n"
        b"unit.so => /usr/bin/true (0x00007fff00001000)\n")


def runtime_closure() -> dict:
    rows = synthetic_discovery_rows()
    dependencies = [
        acquire._runtime_file_identity("/usr/bin/true", "dependency"),
        acquire._runtime_file_identity(
            authority.LDD_INTERPRETER_PATH,
            authority.LDD_INTERPRETER_ROLE),
        acquire._runtime_file_identity(
            authority.LDD_PATH,
            authority.LDD_DISCOVERY_ROLE_PREFIX +
                authority.runtime_discovery_sha256(rows)),
    ]
    dependencies.sort(key=lambda item: (item["role"], item["path"]))
    launchers = [
        acquire._runtime_file_identity(path, "launcher")
        for path in ("/usr/bin/prlimit", "/usr/bin/taskset")
    ]
    return authority.validate_runtime_closure({
        "schema": authority.RUNTIME_CLOSURE_SCHEMA,
        "dependencies": dependencies,
        "launchers": launchers,
    })


def runtime_discovery_verdict(closure: dict) -> dict:
    discovery = next(
        item for item in closure["dependencies"]
        if item["role"].startswith(authority.LDD_DISCOVERY_ROLE_PREFIX))
    interpreter = next(
        item for item in closure["dependencies"]
        if item["role"] == authority.LDD_INTERPRETER_ROLE)
    rows = synthetic_discovery_rows()
    return {
        "schema": authority.RUNTIME_DISCOVERY_SCHEMA,
        "transcript_sha256": authority.runtime_discovery_sha256(rows),
        "ldd_sha256": discovery["sha256"],
        "interpreter_sha256": interpreter["sha256"],
        "file_row_count": 1,
        "virtual_row_count": 1,
        "replay_count": 2,
    }


def synthetic_payloads() -> tuple[dict[str, bytes], dict, bytes]:
    commit = "1" * 40
    tree = "2" * 40
    artifact = synthetic_elf()
    git_data = contract.canonical_json_bytes({})
    archive = b"synthetic archive bytes"
    candidate_digest = digest(artifact)
    provenance = {
        "tracked_source_manifest": {
            "git": {"commit": commit, "tree": tree, "dirty": False},
        },
        "executable": {
            "sha256": candidate_digest, "size": len(artifact),
        },
    }
    provenance_data = contract.canonical_json_bytes(provenance)
    core = {"schema": "unit-reproducible-core/v1"}
    core_data = contract.canonical_json_bytes(core)
    proof_data = contract.canonical_json_bytes({})
    header = b"#define UNIT_K65_AUTHORITY 1\n"
    runtime_data = contract.canonical_json_bytes(runtime_closure())
    payload = {
        "artifacts/bench_leopard2": artifact,
        "build/benchmark-source-attestation.h": header,
        "build/build-provenance.json": provenance_data,
        "build/reproducible-build-core.json": core_data,
        "build/reproducible-build-proof.json": proof_data,
        "runtime/runtime-closure.json": runtime_data,
        "source/candidate-source.tar": archive,
        "source/git-capture.json": git_data,
    }
    identities = {
        path: authority.payload_identity(data) for path, data in payload.items()
    }
    record = authority.candidate_authority_record(
        source_commit=commit, source_tree=tree,
        controller_bindings_sha256="a" * 64,
        payload_identities=identities)
    return payload, record, core_data


def assignment_literal(path: Path, name: str):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(target, ast.Name) and target.id == name
                   for target in targets):
                return ast.literal_eval(node.value)
    raise AssertionError(f"assignment {name} is absent from {path}")


class CandidateAuthorityTests(unittest.TestCase):
    def assertAuthorityError(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(authority.CandidateAuthorityError):
            function(*args, **kwargs)

    def assertAcquireError(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(acquire.CandidateAcquisitionError):
            function(*args, **kwargs)

    def test_contract_constants_match_the_live_consumer(self) -> None:
        runner = HERE / "run_k65r65_b64_packed_terminal_gen3_acquire.py"
        self.assertEqual(
            assignment_literal(runner, "CANDIDATE_AUTHORITY_SCHEMA"),
            authority.CANDIDATE_AUTHORITY_SCHEMA)
        self.assertEqual(
            assignment_literal(runner, "_CANDIDATE_PAYLOAD_PATHS"),
            authority.PAYLOAD_PATHS)

    def test_unchanged_live_consumer_accepts_discovery_binding(self) -> None:
        closure = runtime_closure()
        encoded = contract.canonical_json_bytes(closure)
        bootstrap = HERE / "k65_gen3_exact_source_bootstrap.py"
        program = (
            "import json,runpy,sys;"
            "surface=runpy.run_path(sys.argv[1]);"
            "live=surface['load_live_armer']();"
            "value=json.loads(sys.stdin.buffer.read());"
            "validated=live._candidate_runtime_closure(value);"
            "sys.stdout.buffer.write(live.contract.canonical_json_bytes("
            "validated))"
        )
        result = subprocess.run(
            [sys.executable, "-I", "-S", "-B", "-c", program,
             str(bootstrap)], input=encoded, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, cwd=authority.REPO_ROOT,
            env={
                "HOME": "/nonexistent", "LANG": "C", "LC_ALL": "C",
                "PATH": "/usr/bin:/bin", "TZ": "UTC",
            }, timeout=30, check=False)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertEqual(result.stderr, b"")
        self.assertEqual(result.stdout, encoded)

    def test_record_constructor_closes_inventory_joins(self) -> None:
        payload, record, unused_core = synthetic_payloads()
        self.assertEqual(
            authority.validate_candidate_authority_record(record), record)
        self.assertEqual(
            [item["relative_path"] for item in record["inventory"]],
            list(authority.PAYLOAD_PATHS))
        self.assertEqual(
            record["candidate"]["size"],
            len(payload["artifacts/bench_leopard2"]))
        mutations = []
        wrong_size = copy.deepcopy(record)
        wrong_size["inventory"][0]["size"] += 1
        mutations.append(wrong_size)
        wrong_hash = copy.deepcopy(record)
        wrong_hash["inventory"][0]["sha256"] = "f" * 64
        mutations.append(wrong_hash)
        wrong_order = copy.deepcopy(record)
        wrong_order["inventory"][0], wrong_order["inventory"][1] = \
            wrong_order["inventory"][1], wrong_order["inventory"][0]
        mutations.append(wrong_order)
        extra = copy.deepcopy(record)
        extra["unexpected"] = True
        mutations.append(extra)
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.assertAuthorityError(
                    authority.validate_candidate_authority_record, mutation)

    def test_runtime_closure_rehashes_exact_launcher_targets(self) -> None:
        closure = runtime_closure()
        authority._verify_runtime_files(closure)
        self.assertEqual(
            [item["path"] for item in closure["launchers"]],
            ["/usr/bin/prlimit", "/usr/bin/taskset"])
        corrupted = copy.deepcopy(closure)
        corrupted["launchers"][0]["sha256"] = "0" * 64
        self.assertAuthorityError(authority._verify_runtime_files, corrupted)
        swapped = copy.deepcopy(closure)
        swapped["launchers"].reverse()
        self.assertAuthorityError(authority.validate_runtime_closure, swapped)
        missing = copy.deepcopy(closure)
        missing["dependencies"] = [
            item for item in missing["dependencies"]
            if not item["role"].startswith(
                authority.LDD_DISCOVERY_ROLE_PREFIX)]
        self.assertAuthorityError(
            authority.validate_runtime_closure, missing)
        malformed = copy.deepcopy(closure)
        discovery = next(
            item for item in malformed["dependencies"]
            if item["role"].startswith(authority.LDD_DISCOVERY_ROLE_PREFIX))
        discovery["role"] = authority.LDD_DISCOVERY_ROLE_PREFIX + "z" * 64
        self.assertAuthorityError(
            authority.validate_runtime_closure, malformed)

    def test_runtime_transcript_digest_binds_sonames_and_virtual_rows(self) \
            -> None:
        first = authority.normalize_runtime_discovery(
            b"linux-vdso.so.1 (0x1111)\n"
            b"libc.so.6 => /runtime/libc.so.6 (0x2222)\n")
        aslr_only = authority.normalize_runtime_discovery(
            b"linux-vdso.so.1 (0xaaaa)\n"
            b"libc.so.6 => /runtime/libc.so.6 (0xbbbb)\n")
        renamed = authority.normalize_runtime_discovery(
            b"linux-vdso.so.1 (0x1111)\n"
            b"other.so => /runtime/libc.so.6 (0x2222)\n")
        self.assertEqual(
            authority.runtime_discovery_sha256(first),
            authority.runtime_discovery_sha256(aslr_only))
        self.assertNotEqual(
            authority.runtime_discovery_sha256(first),
            authority.runtime_discovery_sha256(renamed))

    def test_runtime_discovery_is_repeated_and_filters_virtual_rows(self) -> None:
        raw = (
            b"\tlinux-vdso.so.1 (0x00007fff00000000)\n"
            b"\tlibc.so.6 => /runtime/libc.so.6 (0x00007fff00001000)\n"
            b"\t/runtime/ld-linux.so.2 (0x00007fff00002000)\n"
        )

        class Environment:
            def __init__(self):
                self.calls = []

            def run(self, argv, **kwargs):
                self.calls.append((argv, kwargs))
                return SimpleNamespace(
                    exit_status=0, stdout=raw, stderr=b"")

        environment = Environment()

        def identity(path, role):
            return {
                "path": path,
                "sha256": hashlib.sha256(path.encode("ascii")).hexdigest(),
                "size": len(path), "mode": 0o755, "role": role,
            }

        with mock.patch.object(
                acquire, "_runtime_file_identity", side_effect=identity):
            closure = acquire.capture_runtime_closure(
                environment, Path("/not-executed/bench_leopard2"),
                Path("/tmp"))
        self.assertEqual(len(environment.calls), 2)
        self.assertTrue(all(call[0][0] == "/usr/bin/ldd"
                            for call in environment.calls))
        self.assertEqual(
            [item["path"] for item in closure["dependencies"]
             if item["role"] == "dependency"],
            ["/runtime/ld-linux.so.2", "/runtime/libc.so.6"])
        discovery = next(
            item for item in closure["dependencies"]
            if item["role"].startswith(authority.LDD_DISCOVERY_ROLE_PREFIX))
        self.assertEqual(discovery["path"], authority.LDD_PATH)
        self.assertEqual(
            discovery["role"], authority.LDD_DISCOVERY_ROLE_PREFIX +
                authority.runtime_discovery_sha256(
                    authority.normalize_runtime_discovery(raw)))

    def test_generic_seal_replays_all_candidate_semantics(self) -> None:
        payload, record, core_data = synthetic_payloads()
        record_data = contract.canonical_json_bytes(record)
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-candidate-authority.") as temporary:
            lane = (Path(temporary) / "lane").resolve()
            with acquire.sealed_acquire.LaneWriter(str(lane)) as writer:
                seal = writer.seal_payload(
                    terminal=authority.TERMINAL_PATH,
                    terminal_content=record_data,
                    retained_files=payload)
            validated_git = {
                "tree": record["source"]["tree"],
                "detached": True,
                "tracked_files": [],
            }
            header = payload["build/benchmark-source-attestation.h"]
            core = contract.strict_json_loads(
                core_data, "unit reproducible core")
            runtime = contract.strict_json_loads(
                payload["runtime/runtime-closure.json"], "unit runtime")
            with mock.patch.object(
                    authority.git_capture, "validate_git_capture",
                    return_value=validated_git), mock.patch.object(
                    authority, "controller_bindings_sha256",
                    return_value="a" * 64), mock.patch.object(
                    authority, "_verify_source_archive"), mock.patch.object(
                    authority.build_provenance,
                    "validate_reproducible_build_proof"), mock.patch.object(
                    authority.build_provenance, "compare_reproducible_builds",
                    return_value=core), mock.patch.object(
                    authority.build_provenance,
                    "_canonical_replay_attestation_header_bytes",
                    return_value=header), mock.patch.object(
                    authority, "_verify_runtime_discovery",
                    return_value=runtime_discovery_verdict(runtime)):
                verdict = authority.verify_candidate_authority_lane(
                    lane, authority.REPO_ROOT,
                    expected_controller_bindings_sha256="a" * 64,
                    require_detached=True)
            self.assertEqual(
                verdict["record"]["sha256"], digest(record_data))
            self.assertEqual(
                verdict["seal"]["sha256sums_sha256"],
                seal["sha256sums_sha256"])
            self.assertEqual(
                authority.validate_candidate_authority_verdict(
                    verdict, expected_root=lane,
                    expected_record_sha256=digest(record_data),
                    expected_ledger_sha256=seal["sha256sums_sha256"],
                    expected_controller_bindings_sha256="a" * 64),
                verdict)

    def test_verifier_rejects_alias_lane_and_foreign_source_root(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-candidate-paths.") as temporary:
            root = Path(temporary).resolve()
            lane = root / "lane"
            lane.mkdir(mode=0o700)
            alias = root / "alias"
            alias.symlink_to(lane, target_is_directory=True)
            self.assertAuthorityError(
                authority.verify_candidate_authority_lane,
                alias, authority.REPO_ROOT)
            self.assertAuthorityError(
                authority.verify_candidate_authority_lane,
                lane, root)

    def test_fresh_receipt_rejects_noncanonical_or_unjoined_stdout(self) -> None:
        payload, record, unused_core = synthetic_payloads()
        record_data = contract.canonical_json_bytes(record)
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-candidate-verdict.") as temporary:
            lane = Path(temporary).resolve()
            verdict = {
                "schema": authority.VERDICT_SCHEMA,
                "outcome": "promoted_authority", "promoted": True,
                "root": str(lane),
                "record": {
                    "relative_path": authority.TERMINAL_PATH,
                    "schema": authority.CANDIDATE_AUTHORITY_SCHEMA,
                    "sha256": digest(record_data),
                },
                "seal": {
                    "sha256sums_sha256": "b" * 64,
                    "checksum_line_count": len(authority.SEALED_PATHS) - 1,
                    "metadata_entry_count": len(authority.SEALED_PATHS) + 4,
                },
                "source": {
                    "commit": record["source"]["commit"],
                    "tree": record["source"]["tree"], "detached": True,
                },
                "candidate": copy.deepcopy(record["candidate"]),
                "build_closure": copy.deepcopy(record["build_closure"]),
                "controller_bindings_sha256": "a" * 64,
                "inventory_file_count": len(authority.PAYLOAD_PATHS),
                "runtime_file_count": len(runtime_closure()["dependencies"]) +
                    len(runtime_closure()["launchers"]),
                "runtime_discovery":
                    runtime_discovery_verdict(runtime_closure()),
            }

            class Environment:
                def __init__(self, stdout):
                    self.stdout = stdout

                def run(self, unused_argv, **unused_kwargs):
                    return SimpleNamespace(
                        exit_status=0, stdout=self.stdout, stderr=b"")

            accepted = acquire._fresh_verdict(
                Environment(contract.canonical_json_bytes(verdict)),
                lane_root=lane, source_root=authority.REPO_ROOT,
                controller_sha256="a" * 64, record=record,
                record_sha256=digest(record_data), ledger_sha256="b" * 64)
            self.assertEqual(accepted, verdict)
            program = acquire.verification_program_identity()
            changed = copy.deepcopy(program)
            changed["files"][0]["sha256"] = "0" * 64
            with mock.patch.object(
                    acquire, "verification_program_identity",
                    side_effect=(program, changed)):
                self.assertAcquireError(
                    acquire._fresh_verdict,
                    Environment(contract.canonical_json_bytes(verdict)),
                    lane_root=lane, source_root=authority.REPO_ROOT,
                    controller_sha256="a" * 64, record=record,
                    record_sha256=digest(record_data),
                    ledger_sha256="b" * 64,
                    expected_program=program)
            corrupted = copy.deepcopy(verdict)
            corrupted["build_closure"]["build_provenance_sha256"] = "c" * 64
            self.assertAcquireError(
                acquire._fresh_verdict,
                Environment(contract.canonical_json_bytes(corrupted)),
                lane_root=lane, source_root=authority.REPO_ROOT,
                controller_sha256="a" * 64, record=record,
                record_sha256=digest(record_data), ledger_sha256="b" * 64)
            noncanonical = contract.canonical_json_bytes(verdict) + b" "
            self.assertAcquireError(
                acquire._fresh_verdict, Environment(noncanonical),
                lane_root=lane, source_root=authority.REPO_ROOT,
                controller_sha256="a" * 64, record=record,
                record_sha256=digest(record_data), ledger_sha256="b" * 64)

    def test_verifier_program_executes_bound_fds_without_repo_pyc(self) -> None:
        program = acquire.verification_program_identity()
        environment = acquire.sealed_acquire.HostEnvironment(
            canonical_lock_path=acquire.DEFAULT_CANONICAL_LOCK)
        with acquire.VerificationProgramGuard(program) as guard:
            interpreter = guard.descriptor("python")
            verifier = guard.descriptor("authority_verifier")
            result = environment.run(
                [f"/proc/self/fd/{interpreter}", "-v", "-I", "-S", "-B",
                 "-X", "pycache_prefix=/dev/null",
                 f"/proc/self/fd/{verifier}", "--help"],
                cwd=str(authority.REPO_ROOT),
                env=acquire.sealed_acquire.frozen_child_environment(),
                timeout=30, maximum_bytes=4 * 1024 * 1024,
                inherited_descriptors=(interpreter, verifier),
                executable_descriptor=interpreter)
            guard.validate_current()
        self.assertEqual(result.exit_status, 0)
        self.assertTrue(result.stdout.startswith(b"usage:"))
        self.assertNotIn(b"__pycache__", result.stderr)
        self.assertIn(b"leopard2_build_provenance.py", result.stderr)

    def test_entrypoints_ignore_self_deleting_repository_shadows(self) -> None:
        relative_files = (
            "experiments/leopard2/gf8_high_encode/"
                "k65_gen3_candidate_authority.py",
            "experiments/leopard2/gf8_high_encode/"
                "k65_gen3_candidate_authority_acquire.py",
            "experiments/leopard2/gf8_high_encode/"
                "k65_gen3_preregistration.py",
            "experiments/leopard2/main_compare/git_capture.py",
            "experiments/leopard2/main_compare/"
                "pair_qualification_contract.py",
            "experiments/leopard2/main_compare/"
                "pair_qualification_bridge_contract.py",
            "tools/leopard2_build_provenance.py",
            "tools/leopard2_exact_main_baseline.py",
            "tools/leopard2_exact_main_baseline_record.py",
            "tools/leopard2_exact_main_baseline_verifier.py",
            "tools/leopard2_exact_main_baseline_acquire.py",
        )
        cases = (
            ("k65_gen3_candidate_authority.py",
             "leopard2_exact_main_baseline_verifier.py"),
            ("k65_gen3_candidate_authority_acquire.py",
             "leopard2_exact_main_baseline_acquire.py"),
        )
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-exact-imports.") as temporary:
            copied_root = Path(temporary).resolve() / "repository"
            for relative in relative_files:
                source = authority.REPO_ROOT / relative
                destination = copied_root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
            entrypoint_root = copied_root / (
                "experiments/leopard2/gf8_high_encode")
            for index, (entrypoint_name, shadow_name) in enumerate(cases):
                marker = copied_root / f"shadow-executed-{index}"
                shadow = entrypoint_root / shadow_name
                expected = copied_root / "tools" / shadow_name
                shadow.write_text(
                    "from pathlib import Path\n"
                    "_origin = Path(__file__)\n"
                    f"Path({str(marker)!r}).write_text('executed')\n"
                    "_origin.unlink()\n"
                    f"__file__ = {str(expected)!r}\n",
                    encoding="utf-8")
                completed = subprocess.run(
                    [sys.executable, "-I", "-S", "-B", "-X",
                     "pycache_prefix=/dev/null",
                     str(entrypoint_root / entrypoint_name), "--help"],
                    cwd=copied_root, stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    timeout=30, check=False)
                self.assertEqual(
                    completed.returncode, 0, completed.stderr.decode(
                        "utf-8", errors="replace"))
                self.assertFalse(marker.exists())
                self.assertTrue(shadow.is_file())
                shadow.unlink()

    def test_verifier_program_guard_rejects_path_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-program-guard.") as temporary:
            root = Path(temporary).resolve()
            source = root / "verifier.py"
            alternate = root / "alternate.py"
            source.write_bytes(b"print('sealed')\n")
            source.chmod(0o500)
            program = {
                "schema": acquire.VERIFICATION_PROGRAM_SCHEMA,
                "files": [acquire._program_file_identity(
                    "authority_verifier", source)],
            }
            with acquire.VerificationProgramGuard(program) as guard:
                source.rename(alternate)
                alternate.rename(source)
                self.assertAcquireError(guard.validate_current)

    def test_directory_publication_is_atomic_and_no_replace(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-candidate-publish.") as temporary:
            parent = Path(temporary).resolve()
            parent.chmod(0o700)
            staging = parent / "staging"
            output = parent / "authority"
            staging.mkdir(mode=0o700)
            before = staging.stat()
            acquire._rename_directory_noreplace(staging, output)
            after = output.stat()
            self.assertEqual(
                (before.st_dev, before.st_ino),
                (after.st_dev, after.st_ino))
            collision = parent / "collision"
            collision.mkdir(mode=0o700)
            replacement = parent / "replacement"
            replacement.mkdir(mode=0o700)
            self.assertAcquireError(
                acquire._rename_directory_noreplace,
                replacement, collision)
            self.assertTrue(replacement.is_dir())
            self.assertTrue(collision.is_dir())

    def test_failed_postpublication_replay_restores_quarantine(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-candidate-rollback.") as temporary:
            parent = Path(temporary).resolve()
            parent.chmod(0o700)
            staging = parent / "staging"
            output = parent / "authority"
            staging.mkdir(mode=0o700)
            before = staging.stat()
            lock = SimpleNamespace(validate_current=lambda: None)
            with mock.patch.object(
                    acquire, "_fresh_verdict",
                    side_effect=acquire.CandidateAcquisitionError("rejected")):
                self.assertAcquireError(
                    acquire._publish_verified_receipt,
                    SimpleNamespace(), staging=staging, output_lane=output,
                    source_root=authority.REPO_ROOT,
                    controller_sha256="a" * 64, record={},
                    record_sha256="b" * 64, ledger_sha256="c" * 64,
                    program={}, runtime={}, lock=lock)
            self.assertFalse(output.exists())
            self.assertTrue(staging.is_dir())
            after = staging.stat()
            self.assertEqual(
                (before.st_dev, before.st_ino),
                (after.st_dev, after.st_ino))

            fsync_staging = parent / "fsync-staging"
            fsync_output = parent / "fsync-authority"
            fsync_staging.mkdir(mode=0o700)
            fsync_before = fsync_staging.stat()
            real_fsync = os.fsync
            calls = 0

            def fail_first_fsync(descriptor):
                nonlocal calls
                calls += 1
                if calls == 1:
                    raise OSError("injected post-rename fsync failure")
                return real_fsync(descriptor)

            with mock.patch.object(
                    acquire.os, "fsync", side_effect=fail_first_fsync):
                self.assertAcquireError(
                    acquire._publish_verified_receipt,
                    SimpleNamespace(), staging=fsync_staging,
                    output_lane=fsync_output,
                    source_root=authority.REPO_ROOT,
                    controller_sha256="a" * 64, record={},
                    record_sha256="b" * 64, ledger_sha256="c" * 64,
                    program={}, runtime={}, lock=lock)
            self.assertFalse(fsync_output.exists())
            self.assertTrue(fsync_staging.is_dir())
            fsync_after = fsync_staging.stat()
            self.assertEqual(
                (fsync_before.st_dev, fsync_before.st_ino),
                (fsync_after.st_dev, fsync_after.st_ino))

    def test_runtime_hashing_rejects_final_pathname_swap(self) -> None:
        def changed(metadata):
            names = (
                "st_dev", "st_ino", "st_mode", "st_uid", "st_gid",
                "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns",
            )
            values = {name: getattr(metadata, name) for name in names}
            values["st_ino"] += 1
            return SimpleNamespace(**values)

        target = Path("/usr/bin/taskset").resolve(strict=True)
        real_stat = os.stat
        calls = 0

        def swapped(path, *args, **kwargs):
            nonlocal calls
            result = real_stat(path, *args, **kwargs)
            try:
                is_target = Path(os.fspath(path)) == target
            except TypeError:
                is_target = False
            if is_target and kwargs.get("follow_symlinks") is False:
                calls += 1
                if calls == 2:
                    return changed(result)
            return result

        with mock.patch.object(acquire.os, "stat", side_effect=swapped):
            self.assertAcquireError(
                acquire._runtime_file_identity,
                "/usr/bin/taskset", "launcher")
        closure = runtime_closure()
        calls = 0
        with mock.patch.object(authority.os, "stat", side_effect=swapped):
            self.assertAuthorityError(authority._verify_runtime_files, closure)

    def test_producer_has_no_live_campaign_entrypoint(self) -> None:
        source = Path(acquire.__file__).read_text(encoding="utf-8")
        for forbidden in (
                "acquire_and_arm(", ".run_all(", "pair_qualification_acquire",
                "--timing", "benchmark_matrix"):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
