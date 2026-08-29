#!/usr/bin/python3

"""Frozen replay boundaries for the v18 windowed-observation generation."""

from __future__ import annotations

import importlib.util
import base64
import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MAIN_COMPARE = ROOT / "experiments/leopard2/main_compare"
SEALED_V17_FAILURE = (
    ROOT / ".research/leopard-79h/2cc900f-v17-passive-main-v1")


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load replay dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "v18_replay_fixtures", MAIN_COMPARE / "test_run_abba.py")
runner = fixtures.runner
auditor = load_module(
    "v18_replay_auditor",
    MAIN_COMPARE / "audit_v17_gfni_main_compare.py")


class V18ReplayCompatibilityTests(unittest.TestCase):
    def test_frozen_v17_auditor_self_test_stdout_is_preserved(self) -> None:
        completed = subprocess.run(
            ["/usr/bin/python3", "-I", "-S", "-B",
             str(MAIN_COMPARE / "audit_v17_gfni_main_compare.py"),
             "--self-test"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            completed.stdout,
            "v17 GFNI exact-main independent auditor self-test passed\n")

    def test_v1_through_v17_semantic_projections_are_frozen(self) -> None:
        expected = {
            **{generation:
               "2c617b922fab3dc1cba880507d660d1cace2430a58398ddecb3ae303b820e523"
               for generation in range(1, 10)},
            **{generation:
               "e8c816bbbc1460965ae0a7b9535ebc0ae34436a6704873ff33741f05bf30924c"
               for generation in (10, 11)},
            **{generation:
               "7b0b219de41910f1af717ff40a43fb66b592367510f56f0cc9366644896e03e4"
               for generation in range(12, 17)},
            17:
               "967b437ec057215c2f2adad40c048faad4aa353b0da45c10fe60e4d03fc5375d",
        }
        for generation, expected_digest in expected.items():
            with self.subTest(generation=generation):
                raw_schema = getattr(runner, f"RAW_SCHEMA_V{generation}")
                raw = fixtures.synthetic_raw(raw_schema=raw_schema)
                cell = runner.Cell(**raw["campaign"]["cells"][0])
                projection = {
                    "statistics": runner.statistics_policy(raw_schema),
                    "environment":
                        runner.child_environment_for_raw_schema(raw_schema),
                    "baseline_arguments": runner.benchmark_arguments(
                        "baseline", Path("/baseline"), cell,
                        raw["campaign"], raw_schema),
                    "candidate_arguments": runner.benchmark_arguments(
                        "candidate", Path("/candidate"), cell,
                        raw["campaign"], raw_schema),
                    "isolation": raw.get("isolation"),
                }
                self.assertEqual(
                    runner.sha256_bytes(runner.canonical_bytes(projection)),
                    expected_digest)

    def test_current_auditor_replays_a_frozen_v17_bundle(self) -> None:
        blob = b"production selector fixture\n"
        blob_id = hashlib.sha1(
            f"blob {len(blob)}\0".encode("ascii") + blob).hexdigest()
        tree_bytes = b"100644 leopard2.cpp\0" + bytes.fromhex(blob_id)
        tree_id = hashlib.sha1(
            f"tree {len(tree_bytes)}\0".encode("ascii") +
            tree_bytes).hexdigest()
        commit_bytes = (
            f"tree {tree_id}\n"
            "author Fixture <fixture@example.invalid> 0 +0000\n"
            "committer Fixture <fixture@example.invalid> 0 +0000\n\n"
            "fixture\n"
        ).encode("ascii")
        commit_id = hashlib.sha1(
            f"commit {len(commit_bytes)}\0".encode("ascii") +
            commit_bytes).hexdigest()
        original_source_fixture = fixtures.complete_source_fixture

        def source_fixture(role: str, raw_schema: str = runner.RAW_SCHEMA_V16
                           ) -> dict:
            if role == "baseline":
                return original_source_fixture(role, raw_schema)
            return fixtures.rich_git_source_fixture(
                fixtures.SPECIFICATION["candidate_source_root"],
                commit_id, commit_bytes,
                {tree_id: base64.b64encode(tree_bytes).decode("ascii")},
                detached=False)

        with tempfile.TemporaryDirectory(
                prefix="leopard-v18-v17-auditor-replay-") as temporary:
            with (
                mock.patch.multiple(
                    fixtures,
                    CANDIDATE_TREE=tree_id,
                    CANDIDATE_COMMIT_RAW=commit_bytes,
                    CANDIDATE_COMMIT=commit_id),
                mock.patch.dict(
                    fixtures.SPECIFICATION,
                    {"candidate_commit": commit_id}),
                mock.patch.object(
                    fixtures, "complete_source_fixture",
                    side_effect=source_fixture),
            ):
                raw = fixtures.synthetic_raw(
                    raw_schema=runner.RAW_SCHEMA_V17)

            specification = raw["input_specification"]
            runner_path = (
                specification["candidate_source_root"] +
                "/experiments/leopard2/main_compare/run_abba.py")
            specification["runner"] = runner_path
            raw["identities_initial"]["runner"]["path"] = runner_path
            raw["identities_initial"]["baseline_build"][
                "validated_cache"]["CMAKE_CXX_FLAGS_RELEASE"] = \
                "-g -O0 -O3"
            for role in ("baseline", "candidate"):
                build = raw["identities_initial"][f"{role}_build"]
                build["archiver"]["path"] = \
                    "/usr/bin/x86_64-linux-gnu-ar"
                build["ranlib"]["path"] = \
                    "/usr/bin/x86_64-linux-gnu-ranlib"
                build["archive_link_tool_invocations"]["archiver"][
                    "resolved_path"] = \
                    "/usr/bin/x86_64-linux-gnu-ar"
                build["archive_link_tool_invocations"]["ranlib"][
                    "resolved_path"] = \
                    "/usr/bin/x86_64-linux-gnu-ranlib"
                build["compiler"]["path"] = "/usr/bin/g++"
                build["compiler_invocation"]["resolved_path"] = "/usr/bin/g++"
            candidate_entries = raw["identities_initial"][
                "candidate_build"]["validated_compile_commands"][
                    "required_entries"]
            matching_entries = [
                entry for entry in candidate_entries
                if entry["file"].endswith(
                    "/Leopard2BackendAVX2T8K8B1024.cpp")]
            self.assertEqual(len(matching_entries), 1)
            arguments = matching_entries[0]["arguments"]
            arguments.insert(
                arguments.index("-falign-functions=64"),
                "-flive-range-shrinkage")

            raw["supervision"] = None
            for invocation in raw["invocations"]:
                invocation["duration_ns"] = 500_000_000
            raw["isolation"]["after"]["monotonic_ns"] = (
                raw["isolation"]["before"]["monotonic_ns"] +
                sum(item["duration_ns"] for item in raw["invocations"]))
            for cpu_name in ("benchmark_cpu", "reserved_sibling"):
                raw["host_initial"][cpu_name]["cpuinfo"]["flags"] = \
                    "avx2 gfni"
            raw["host_final"] = copy.deepcopy(raw["host_initial"])
            fixtures.synchronize_identity(raw)

            evidence_root = Path(temporary)
            cell_id = raw["campaign"]["cells"][0]["identifier"]
            for invocation in raw["invocations"]:
                prefix = (
                    f"invocations/{cell_id}/round-{invocation['round']}/"
                    f"slot-{invocation['slot']}-"
                    f"{invocation['implementation']}")
                stdout = json.dumps(invocation["result"]).encode("utf-8")
                for stream_name, stream_bytes in (
                        ("stdout", stdout), ("stderr", b"")):
                    relative = f"{prefix}.{stream_name}"
                    stream_path = evidence_root / relative
                    stream_path.parent.mkdir(parents=True, exist_ok=True)
                    stream_path.write_bytes(stream_bytes)
                    stream_path.chmod(0o600)
                    invocation[stream_name] = {
                        "path": relative,
                        "size": len(stream_bytes),
                        "sha256": runner.sha256_bytes(stream_bytes),
                    }

            raw = fixtures.resign(raw)
            raw_path = evidence_root / "raw.json"
            runner.write_json_exclusive(raw_path, raw)
            manifest = runner.signed({
                "schema": runner.MANIFEST_SCHEMA_V17,
                "created_utc": "2026-07-16T00:00:00Z",
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {
                    "path": "raw.json",
                    "size": raw_path.stat().st_size,
                    "sha256": runner.sha256_file(raw_path),
                    "payload_digest": raw["digest"],
                },
                "campaign": raw["campaign"],
                "host": raw["host_initial"],
                "isolation": raw["isolation"],
                "reservation": raw["reservation"],
                "identities": raw["identities_initial"],
                "analysis": raw["analysis"],
                "supervision": None,
                "executable_snapshots": raw["executable_snapshots"],
            })
            manifest_path = evidence_root / "manifest.json"
            runner.write_json_exclusive(manifest_path, manifest)
            result = auditor.replay(
                manifest_path, supervision_mode="absent")

        self.assertEqual(
            result["schema"],
            "leopard2-main-compare-v17-passive-independent-audit/v1")
        self.assertEqual(result["audit_mode"],
                         "passive-shared-host-observation")
        self.assertEqual(result["status"], "complete")
        self.assertIs(result["audit_passed"], True)
        self.assertEqual(result["replay"]["invocation_count"], 12)
        self.assertEqual(
            result["replay"]["candidate_route_attestation_count"], 6)
        result_without_auditor_identity = copy.deepcopy(result)
        result_without_auditor_identity.pop("auditor")
        self.assertEqual(
            auditor.sha256_bytes(auditor.canonical_bytes(
                result_without_auditor_identity)),
            "b34fcc7219a44a0795659967a61c5c60cfb9ce700d5ff069f37d019ea5801129")

    def test_decoder_promotion_remains_pinned_before_v17_and_v18(self) -> None:
        planner = fixtures.load_plan_runner()
        self.assertEqual(
            planner.EXACT_MANIFEST_SCHEMA,
            "leopard2-main-compare-manifest/v16")
        for generation in (17, 18):
            manifest = f"leopard2-main-compare-manifest/v{generation}"
            raw = f"leopard2-main-compare-raw/v{generation}"
            with self.subTest(generation=generation):
                self.assertNotIn(manifest, planner.EXACT_MANIFEST_SCHEMAS)
                self.assertNotIn(raw, planner.EXACT_RAW_SCHEMAS)
                self.assertNotIn(
                    (manifest, raw), planner.EXACT_SCHEMA_PAIRS)

    @unittest.skipUnless(
        (SEALED_V17_FAILURE / "FAILED.json").is_file(),
        "retained passive-v17 diagnostic envelope is unavailable")
    def test_current_runner_replays_retained_v17_failure_from_owner_only_copy(
            self) -> None:
        # The sealed source is intentionally 0500/0400.  Reconstruct the exact
        # campaign bytes under the producer's required owner-only 0700/0600
        # modes, then exercise the current v18-head runner's frozen v17 branch.
        with tempfile.TemporaryDirectory(
                prefix="leopard-v18-v17-replay-") as temporary:
            reconstructed = Path(temporary) / "campaign"
            shutil.copytree(
                SEALED_V17_FAILURE / "core/campaign", reconstructed)
            os.chmod(reconstructed, 0o700)
            for path in reconstructed.rglob("*"):
                os.chmod(path, 0o700 if path.is_dir() else 0o600)
            completed = subprocess.run(
                ["/usr/bin/python3", "-I", "-S", "-B",
                 str(MAIN_COMPARE / "run_abba.py"), "verify-failure",
                 "--failure", str(reconstructed / "failure.json")],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, text=True,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            completed.stdout.strip(),
            "failed exact-main campaign diagnostics and retained files "
            "verified; this is not valid performance evidence")

    @unittest.skipUnless(
        (SEALED_V17_FAILURE / "FAILED.json").is_file(),
        "retained passive-v17 diagnostic envelope is unavailable")
    def test_retained_v17_failure_verifies_through_owner_only_replay(
            self) -> None:
        # Sealed campaign nodes are deliberately 0500/0400.  The durable
        # wrapper verifier reconstructs their exact bytes as 0700/0600 before
        # invoking the frozen v17 producer, rather than weakening its mode gate.
        completed = subprocess.run(
            ["/usr/bin/bash", str(
                MAIN_COMPARE /
                "run_authoritative_v17_gfni_main_compare.sh"),
             "--verify", str(SEALED_V17_FAILURE)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("envelope verified", completed.stdout)


if __name__ == "__main__":
    unittest.main()
