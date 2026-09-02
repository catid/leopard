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


def synthetic_windowed_auditor_raw(raw_schema: str) -> dict:
    """Build a canonical v18/v19 fixture for the independent auditor."""
    if raw_schema not in (runner.RAW_SCHEMA_V18, runner.RAW_SCHEMA_V19):
        raise ValueError("windowed auditor fixture requires v18 or v19")
    blob = b"production selector fixture\n"
    blob_id = hashlib.sha1(
        f"blob {len(blob)}\0".encode("ascii") + blob,
        usedforsecurity=False).hexdigest()
    tree_bytes = b"100644 leopard2.cpp\0" + bytes.fromhex(blob_id)
    tree_id = hashlib.sha1(
        f"tree {len(tree_bytes)}\0".encode("ascii") +
        tree_bytes, usedforsecurity=False).hexdigest()
    commit_bytes = (
        f"tree {tree_id}\n"
        "author Fixture <fixture@example.invalid> 0 +0000\n"
        "committer Fixture <fixture@example.invalid> 0 +0000\n\n"
        "fixture\n"
    ).encode("ascii")
    commit_id = hashlib.sha1(
        f"commit {len(commit_bytes)}\0".encode("ascii") +
        commit_bytes, usedforsecurity=False).hexdigest()
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

    with (
            mock.patch.multiple(
                fixtures, CANDIDATE_TREE=tree_id,
                CANDIDATE_COMMIT_RAW=commit_bytes,
                CANDIDATE_COMMIT=commit_id),
            mock.patch.dict(
                fixtures.SPECIFICATION, {"candidate_commit": commit_id}),
            mock.patch.object(
                fixtures, "complete_source_fixture",
                side_effect=source_fixture),
    ):
        raw = fixtures.synthetic_raw(raw_schema=raw_schema)

    specification = raw["input_specification"]
    runner_path = (
        specification["candidate_source_root"] +
        "/experiments/leopard2/main_compare/run_abba.py")
    specification["runner"] = runner_path
    raw["identities_initial"]["runner"]["path"] = runner_path
    raw["identities_initial"]["baseline_build"]["validated_cache"][
        "CMAKE_CXX_FLAGS_RELEASE"] = "-g -O0 -O3"
    for role in ("baseline", "candidate"):
        build = raw["identities_initial"][f"{role}_build"]
        build["archiver"]["path"] = "/usr/bin/x86_64-linux-gnu-ar"
        build["ranlib"]["path"] = "/usr/bin/x86_64-linux-gnu-ranlib"
        build["archive_link_tool_invocations"]["archiver"][
            "resolved_path"] = "/usr/bin/x86_64-linux-gnu-ar"
        build["archive_link_tool_invocations"]["ranlib"][
            "resolved_path"] = "/usr/bin/x86_64-linux-gnu-ranlib"
        build["compiler"]["path"] = "/usr/bin/g++"
        build["compiler_invocation"]["resolved_path"] = "/usr/bin/g++"
    candidate_entries = raw["identities_initial"]["candidate_build"][
        "validated_compile_commands"]["required_entries"]
    matching_entries = [
        entry for entry in candidate_entries
        if entry["file"].endswith("/Leopard2BackendAVX2T8K8B1024.cpp")
    ]
    if len(matching_entries) != 1:
        raise RuntimeError("windowed fixture lacks its calibrated AVX2 entry")
    arguments = matching_entries[0]["arguments"]
    arguments.insert(
        arguments.index("-falign-functions=64"),
        "-flive-range-shrinkage")

    for name in ("benchmark_cpu", "reserved_sibling"):
        raw["host_initial"][name]["cpuinfo"]["flags"] = "avx2 gfni"
    raw["host_final"] = copy.deepcopy(raw["host_initial"])
    qualification = raw.get("pair_qualification")
    if raw_schema == runner.RAW_SCHEMA_V19:
        raw["pair_qualification"] = \
            runner.pair_v19.pair_qualified_v19_record(
                stage="complete", record_status="complete",
                terminal=runner.pair_v19.SUCCESS_TERMINAL,
                policy_value=qualification["policy"],
                expected_policy_sha256=qualification["policy_sha256"],
                host_identity_sha256=
                    runner.pair_v19.contract.canonical_sha256(
                        raw["host_initial"]),
                attempt_value=qualification["attempt"],
                acquisition_value=qualification["acquisition"],
                bridge_value=qualification["bridge"],
                first_window_before_value=
                    qualification["first_window_handoff"]
                        ["first_window_before"])
        qualification = raw["pair_qualification"]
    benchmark_cpu = raw["campaign"]["benchmark_cpu"]
    reserved_sibling = raw["campaign"]["reserved_sibling"]
    duration_ns = 500_000_000
    first_before = copy.deepcopy(
        qualification["first_window_handoff"]["first_window_before"]
        if qualification is not None else
        raw["invocations"][0]["cpu_window"]["before"])
    previous_after = None
    for index, invocation in enumerate(raw["invocations"]):
        if index == 0:
            before = copy.deepcopy(first_before)
        else:
            before = copy.deepcopy(previous_after)
            started = previous_after["read_finished_monotonic_ns"] + 10
            before.update({
                "monotonic_ns": started + 2,
                "read_started_monotonic_ns": started,
                "read_finished_monotonic_ns": started + 2,
            })
            for role in ("benchmark_cpu", "reserved_sibling"):
                before[role]["fields"]["idle"] += 1
                before[role]["total_jiffies"] += 1
        child_started = before["read_finished_monotonic_ns"] + 10
        after = copy.deepcopy(before)
        after["benchmark_cpu"]["fields"]["user"] += 1
        after["benchmark_cpu"]["fields"]["idle"] += 50
        after["benchmark_cpu"]["total_jiffies"] += 51
        after["reserved_sibling"]["fields"]["idle"] += 51
        after["reserved_sibling"]["total_jiffies"] += 51
        after_started = child_started + duration_ns
        after.update({
            "monotonic_ns": after_started,
            "read_started_monotonic_ns": after_started,
            "read_finished_monotonic_ns": after_started + 2,
        })
        invocation["duration_ns"] = duration_ns
        invocation["cpu_window"] = runner.cpu_window_record(
            benchmark_cpu, reserved_sibling, invocation["cell_id"],
            invocation["round"], invocation["slot"],
            invocation["implementation"], before, after,
            child_started, duration_ns)
        previous_after = after

    windows = [invocation["cpu_window"] for invocation in raw["invocations"]]
    if qualification is not None:
        outer_before = qualification["bridge"]["campaign_presample_before"]
        outer_before_ns = outer_before["read_finished_monotonic_ns"]
        outer_before_cpu = outer_before["cpus"][str(benchmark_cpu)]
        outer_before_sibling = outer_before["cpus"][str(reserved_sibling)]
    else:
        outer_before_ns = first_before["monotonic_ns"] - 10
        outer_before_cpu = first_before["benchmark_cpu"]
        outer_before_sibling = first_before["reserved_sibling"]
    raw["isolation"] = runner.isolation_record_v2(
        benchmark_cpu, reserved_sibling, raw["isolation"]["pair_lease"],
        outer_before_ns,
        previous_after["read_finished_monotonic_ns"] + 10,
        outer_before_cpu, previous_after["benchmark_cpu"],
        outer_before_sibling,
        previous_after["reserved_sibling"], windows)
    fixtures.synchronize_identity(raw)
    return raw


def materialize_windowed_auditor_bundle(
    evidence_root: Path, raw_value: dict, *,
    manifest_qualification: object = None,
) -> Path:
    raw = copy.deepcopy(raw_value)
    cell_id = raw["campaign"]["cells"][0]["identifier"]
    for invocation in raw["invocations"]:
        prefix = (
            f"invocations/{cell_id}/round-{invocation['round']}/"
            f"slot-{invocation['slot']}-{invocation['implementation']}")
        stdout = json.dumps(invocation["result"]).encode("utf-8")
        for stream_name, stream_bytes in (("stdout", stdout), ("stderr", b"")):
            relative = f"{prefix}.{stream_name}"
            stream_path = evidence_root / relative
            stream_path.parent.mkdir(parents=True, exist_ok=True)
            stream_path.write_bytes(stream_bytes)
            stream_path.chmod(0o600)
            invocation[stream_name] = {
                "path": relative, "size": len(stream_bytes),
                "sha256": runner.sha256_bytes(stream_bytes),
            }
    raw = fixtures.resign(raw)
    raw_path = evidence_root / "raw.json"
    runner.write_json_exclusive(raw_path, raw)
    raw_schema = raw["schema"]
    manifest_payload = {
        "schema": (
            runner.MANIFEST_SCHEMA_V19
            if raw_schema == runner.RAW_SCHEMA_V19
            else runner.MANIFEST_SCHEMA_V18),
        "created_utc": "2026-07-16T00:00:00Z",
        "valid": True, "validity_is_independent_of_speed": True,
        "raw": {
            "path": "raw.json", "size": raw_path.stat().st_size,
            "sha256": runner.sha256_file(raw_path),
            "payload_digest": raw["digest"],
        },
        "campaign": raw["campaign"], "host": raw["host_initial"],
        "isolation": raw["isolation"], "reservation": raw["reservation"],
        "identities": raw["identities_initial"],
        "analysis": raw["analysis"], "supervision": None,
        "executable_snapshots": raw["executable_snapshots"],
    }
    if raw_schema == runner.RAW_SCHEMA_V19:
        manifest_payload["pair_qualification"] = (
            raw["pair_qualification"] if manifest_qualification is None
            else manifest_qualification)
    manifest = runner.signed(manifest_payload)
    manifest_path = evidence_root / "manifest.json"
    runner.write_json_exclusive(manifest_path, manifest)
    return manifest_path


class V18ReplayCompatibilityTests(unittest.TestCase):
    def test_frozen_v17_auditor_self_test_stdout_is_preserved(self) -> None:
        for optimized in (False, True):
            with self.subTest(optimized=optimized):
                command = ["/usr/bin/python3"]
                if optimized:
                    command.append("-O")
                command.extend((
                    "-I", "-S", "-B",
                    str(MAIN_COMPARE / "audit_v17_gfni_main_compare.py"),
                    "--self-test"))
                completed = subprocess.run(
                    command, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=False, text=True)
                self.assertEqual(
                    completed.returncode, 0, completed.stderr)
                self.assertEqual(
                    completed.stdout,
                    "v17 GFNI exact-main independent auditor self-test "
                    "passed\n")

    def test_v19_auditor_replays_and_rejects_conditioning_splices(self) -> None:
        raw = synthetic_windowed_auditor_raw(runner.RAW_SCHEMA_V19)
        qualification = raw["pair_qualification"]
        expected_attempt = fixtures.expected_v19_attempt_one()
        host_hash = runner.pair_v19.contract.canonical_sha256(
            raw["host_initial"])
        validated, pair = auditor.validate_v19_qualification(
            qualification, expected_attempt,
            expected_host_identity_sha256=host_hash)
        self.assertEqual(pair, (1, 65))
        self.assertEqual(validated["terminal"], "NOT_PROMOTED")
        self.assertIs(validated["bridge"]["bridge_accepted"], True)
        self.assertIs(validated["first_window_handoff"]["accepted"], True)

        def changed(path: tuple[object, ...], value: object) -> dict:
            result = copy.deepcopy(qualification)
            target = result
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            return result

        mutations = {
            "policy": changed(
                ("policy", "candidate_primary_cpus"), list(range(1, 63))),
            "record selected pair": changed(
                ("selected_pair", "reserved_sibling"), 66),
            "bridge": changed(("bridge", "bridge_tail_sha256"), "0" * 64),
            "handoff": changed(
                ("first_window_handoff", "handoff_elapsed_ns"), 0),
            "lineage": changed(
                ("v18_failure_lineage", "attempts", 0,
                 "envelope_sha256sums_sha256"), "0" * 64),
            "claim ceiling": changed(
                ("shared_host_claim_ceiling", "promotion_eligible"), True),
            "status": changed(("record_status",), "failed"),
        }
        for label, mutation in mutations.items():
            with self.subTest(label=label):
                with self.assertRaises(auditor.AuditError):
                    auditor.validate_v19_qualification(
                        mutation, expected_attempt,
                        expected_host_identity_sha256=host_hash)
        wrong_attempt = auditor.v19_attempt_record(
            attempt=2, prior_attempt_failure_sha256="0" * 64,
            prior_attempt_acquisition_sha256="1" * 64,
            prior_attempt_selection_status="no-candidate-pair-qualified")
        with self.assertRaises(auditor.AuditError):
            auditor.validate_v19_qualification(
                qualification, wrong_attempt,
                expected_host_identity_sha256=host_hash)
        with self.assertRaises(auditor.AuditError):
            auditor.validate_v19_qualification(
                qualification, expected_attempt,
                expected_host_identity_sha256="0" * 64)
        with self.assertRaises(TypeError):
            auditor.validate_v19_qualification(
                qualification, expected_attempt)

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-auditor-replay-") as temporary:
            manifest_path = materialize_windowed_auditor_bundle(
                Path(temporary), raw)
            result = auditor.replay(
                manifest_path, supervision_mode="conditioned",
                expected_v19_attempt=expected_attempt)
        self.assertEqual(
            result["schema"],
            "leopard2-main-compare-v19-conditioned-passive-"
            "independent-audit/v1")
        self.assertEqual(result["acquisition_generation"], "passive-v3")
        self.assertEqual(result["contract"]["cpu"], 1)
        self.assertEqual(result["contract"]["sibling"], 65)
        self.assertEqual(
            result["pair_qualification"]["terminal"], "NOT_PROMOTED")
        self.assertEqual(
            result["pair_qualification"]["shared_host_claim_ceiling"], {
                "promotion_eligible": False,
                "host_exclusivity_proved": False,
                "whole_campaign_interval_observed": False,
                "causal_performance_claim_allowed": False,
            })
        self.assertIs(
            result["gate_results"][
                "pair_qualification_fixed_point_independently_replayed"],
            True)

        campaign_splice = copy.deepcopy(raw)
        campaign_splice["campaign"]["allowed_cpu_set_at_launch"].append(128)

        isolation_splice = copy.deepcopy(raw)
        splice_qualification = isolation_splice["pair_qualification"]
        splice_tail = splice_qualification["bridge"][
            "campaign_presample_before"]
        splice_windows = [
            invocation["cpu_window"]
            for invocation in isolation_splice["invocations"]
        ]
        splice_last = splice_windows[-1]["after"]
        isolation_splice["isolation"] = runner.isolation_record_v2(
            1, 65, isolation_splice["isolation"]["pair_lease"],
            splice_tail["read_finished_monotonic_ns"] + 1,
            isolation_splice["isolation"]["after"]["monotonic_ns"],
            splice_tail["cpus"]["1"], splice_last["benchmark_cpu"],
            splice_tail["cpus"]["65"], splice_last["reserved_sibling"],
            splice_windows)

        handoff_splice = copy.deepcopy(raw)
        splice_qualification = handoff_splice["pair_qualification"]
        different_first = copy.deepcopy(
            splice_qualification["first_window_handoff"][
                "first_window_before"])
        for name in ("monotonic_ns", "read_started_monotonic_ns",
                     "read_finished_monotonic_ns"):
            different_first[name] += 1
        handoff_splice["pair_qualification"] = \
            runner.pair_v19.pair_qualified_v19_record(
                stage="complete", record_status="complete",
                terminal=runner.pair_v19.SUCCESS_TERMINAL,
                policy_value=splice_qualification["policy"],
                expected_policy_sha256=splice_qualification["policy_sha256"],
                host_identity_sha256=
                    splice_qualification["host_identity_sha256"],
                attempt_value=splice_qualification["attempt"],
                acquisition_value=splice_qualification["acquisition"],
                bridge_value=splice_qualification["bridge"],
                first_window_before_value=different_first)

        manifest_splice = copy.deepcopy(qualification)
        manifest_splice["terminal"] = "campaign-rejected"
        replay_splices = (
            ("manifest/raw qualification", raw, manifest_splice),
            ("acquisition/campaign affinity", campaign_splice, None),
            ("bridge/isolation presample", isolation_splice, None),
            ("handoff/first retained window", handoff_splice, None),
        )
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-auditor-splices-") as temporary:
            root = Path(temporary)
            for index, (label, spliced_raw, manifest_override) in enumerate(
                    replay_splices):
                with self.subTest(replay_join=label):
                    evidence_root = root / str(index)
                    evidence_root.mkdir(mode=0o700)
                    manifest_path = materialize_windowed_auditor_bundle(
                        evidence_root, spliced_raw,
                        manifest_qualification=manifest_override)
                    with self.assertRaises(auditor.AuditError):
                        auditor.replay(
                            manifest_path, supervision_mode="conditioned",
                            expected_v19_attempt=expected_attempt)

    def test_current_auditor_preserves_frozen_v18_windowed_result(self) -> None:
        raw = synthetic_windowed_auditor_raw(runner.RAW_SCHEMA_V18)
        with tempfile.TemporaryDirectory(
                prefix="leopard-v18-auditor-replay-") as temporary:
            manifest_path = materialize_windowed_auditor_bundle(
                Path(temporary), raw)
            result = auditor.replay(
                manifest_path, supervision_mode="windowed")
        self.assertEqual(
            result["schema"],
            "leopard2-main-compare-v18-passive-independent-audit/v1")
        self.assertEqual(result["acquisition_generation"], "passive-v2")
        self.assertEqual(result["contract"]["cpu"], 52)
        self.assertEqual(result["contract"]["sibling"], 116)
        result_without_auditor_identity = copy.deepcopy(result)
        result_without_auditor_identity.pop("auditor")
        self.assertEqual(
            auditor.sha256_bytes(auditor.canonical_bytes(
                result_without_auditor_identity)),
            "e30bc1ef934a6ed018412fe73dd76a546b0acd3fdc132ff10ef808835435ed81")

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
            f"blob {len(blob)}\0".encode("ascii") + blob,
            usedforsecurity=False).hexdigest()
        tree_bytes = b"100644 leopard2.cpp\0" + bytes.fromhex(blob_id)
        tree_id = hashlib.sha1(
            f"tree {len(tree_bytes)}\0".encode("ascii") +
            tree_bytes, usedforsecurity=False).hexdigest()
        commit_bytes = (
            f"tree {tree_id}\n"
            "author Fixture <fixture@example.invalid> 0 +0000\n"
            "committer Fixture <fixture@example.invalid> 0 +0000\n\n"
            "fixture\n"
        ).encode("ascii")
        commit_id = hashlib.sha1(
            f"commit {len(commit_bytes)}\0".encode("ascii") +
            commit_bytes, usedforsecurity=False).hexdigest()
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
