#!/usr/bin/python3

"""Shared end-to-end fixtures for dormant conditioned v19 evidence."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MAIN_COMPARE = ROOT / "experiments/leopard2/main_compare"
RUNNER_PATH = MAIN_COMPARE / "run_abba.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load v19 end-to-end dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


replay_fixtures = load_module(
    "v19_end_to_end_replay_fixtures",
    MAIN_COMPARE / "test_v18_replay_compatibility.py",
)
passive_fixtures = load_module(
    "v19_end_to_end_passive_fixtures",
    MAIN_COMPARE / "test_passive_evidence.py",
)

runner = replay_fixtures.runner
auditor = replay_fixtures.auditor
census = passive_fixtures.census


def canonical_raw() -> dict:
    return replay_fixtures.synthetic_windowed_auditor_raw(
        runner.RAW_SCHEMA_V19)


def expected_attempt_one() -> dict:
    return replay_fixtures.fixtures.expected_v19_attempt_one()


def expected_attempt_two() -> dict:
    return census.v19_attempt_record(
        attempt=2,
        prior_attempt_failure_sha256="0" * 64,
        prior_attempt_acquisition_sha256="1" * 64,
        prior_attempt_selection_status="no-candidate-pair-qualified",
    )


def live_authority(raw: dict) -> dict:
    qualification = raw["pair_qualification"]
    return runner.prepare_v19_live_run_authority(
        campaign=raw["campaign"], host_identity_value=raw["host_initial"],
        expected_attempt=qualification["attempt"],
        acquisition_value=qualification["acquisition"],
        bridge_value=qualification["bridge"],
        expected_v18_failure_lineage_sha256=
            qualification["v18_failure_lineage_sha256"],
        expected_claim_ceiling_sha256=
            runner.pair_v19.contract.canonical_sha256(
                qualification["shared_host_claim_ceiling"]),
        resource_envelope_value=copy.deepcopy(
            runner.V19_LIVE_RESOURCE_ENVELOPE),
        controller_affinity_value=passive_fixtures.v19_controller(raw),
        environment_census_pre_value=
            passive_fixtures.v19_snapshot(raw, "pre"))


def handed_off_authority(raw: dict) -> dict:
    return runner.bind_v19_live_first_window(
        live_authority(raw), campaign=raw["campaign"],
        host_identity_value=raw["host_initial"],
        first_window_before_value=raw["invocations"][0]
            ["cpu_window"]["before"])


def complete_live_authority(raw: dict) -> dict:
    return runner.complete_v19_live_run_authority(
        handed_off_authority(raw), campaign=raw["campaign"],
        host_identity_value=raw["host_initial"])


def attempt_cli_arguments(attempt: int) -> tuple[str, ...]:
    if attempt == 1:
        return ("--v19-attempt", "1")
    if attempt == 2:
        return (
            "--v19-attempt", "2",
            "--v19-prior-failure-sha256", "0" * 64,
            "--v19-prior-acquisition-sha256", "1" * 64,
            "--v19-prior-selection-status",
            "no-candidate-pair-qualified",
        )
    raise ValueError("fixture attempt must be 1 or 2")


def run_runner_cli(*arguments: str) -> subprocess.CompletedProcess[str]:
    optimization = (
        ("-" + "O" * sys.flags.optimize,) if sys.flags.optimize else ())
    return subprocess.run(
        ["/usr/bin/python3", "-I", "-S", "-B", *optimization,
         str(RUNNER_PATH),
         *arguments],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, text=True, timeout=60,
    )


def consume_bundle(raw: dict, attempt: dict) -> tuple[dict, dict, dict, dict, dict]:
    with tempfile.TemporaryDirectory(
            prefix="leopard-v19-shared-consumers-") as temporary:
        manifest = replay_fixtures.materialize_windowed_auditor_bundle(
            Path(temporary), raw)
        retained_manifest, retained_raw, analysis, _manifest_snapshot = \
            runner.verified_campaign_bundle(
                manifest, no_current_input_check=True,
                expected_v19_attempt=attempt)
        audit_result = auditor.replay(
            manifest, supervision_mode="conditioned",
            expected_v19_attempt=attempt)
        census_result = compare_census(retained_raw, attempt)
        return (retained_manifest, retained_raw, analysis, audit_result,
                census_result)


def census_witnesses(raw: dict) -> tuple[dict, dict, dict]:
    return (
        passive_fixtures.v19_snapshot(raw, "pre"),
        passive_fixtures.v19_snapshot(raw, "post"),
        passive_fixtures.v19_controller(raw),
    )


def compare_census(
    raw: dict,
    attempt: dict,
    witnesses: tuple[dict, dict, dict] | None = None,
) -> dict:
    pre, post, controller = (
        census_witnesses(raw) if witnesses is None else
        tuple(copy.deepcopy(value) for value in witnesses)
    )
    return census.compare(
        pre, post, raw, controller,
        expected_v19_attempt=attempt)


def rebuild_isolation(raw: dict, *, before_ns: int | None = None) -> None:
    selected = raw["pair_qualification"]["selected_pair"]
    benchmark_cpu = selected["benchmark_cpu"]
    reserved_sibling = selected["reserved_sibling"]
    windows = [invocation["cpu_window"] for invocation in raw["invocations"]]
    current = raw["isolation"]
    raw["isolation"] = runner.isolation_record_v2(
        benchmark_cpu, reserved_sibling, current["pair_lease"],
        current["before"]["monotonic_ns"] if before_ns is None else before_ns,
        current["after"]["monotonic_ns"],
        current["before"]["benchmark_cpu"],
        current["after"]["benchmark_cpu"],
        current["before"]["reserved_sibling"],
        current["after"]["reserved_sibling"], windows)


def splice_first_window(raw: dict) -> None:
    invocation = raw["invocations"][0]
    old_window = invocation["cpu_window"]
    changed_before = copy.deepcopy(old_window["before"])
    for name in ("monotonic_ns", "read_started_monotonic_ns",
                 "read_finished_monotonic_ns"):
        changed_before[name] += 1
    selected = raw["pair_qualification"]["selected_pair"]
    invocation["cpu_window"] = runner.cpu_window_record(
        selected["benchmark_cpu"], selected["reserved_sibling"],
        invocation["cell_id"], invocation["round"], invocation["slot"],
        invocation["implementation"], changed_before, old_window["after"],
        old_window["child_started_monotonic_ns"], invocation["duration_ns"])
    rebuild_isolation(raw)


class V19EndToEndTests(unittest.TestCase):
    def assert_all_reject(
            self, raw: dict, attempt: dict, label: str,
            witnesses: tuple[dict, dict, dict]) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-shared-splice-") as temporary:
            root = Path(temporary)
            manifest = replay_fixtures.materialize_windowed_auditor_bundle(
                root, raw)
            retained_raw = runner.strict_json_loads(
                (root / "raw.json").read_bytes(), "shared retained raw")
            if type(retained_raw) is not dict:
                raise RuntimeError("shared retained raw is not an object")
            with self.subTest(consumer="runner", splice=label):
                with self.assertRaises(runner.EvidenceError):
                    runner.verified_campaign_bundle(
                        manifest, no_current_input_check=True,
                        expected_v19_attempt=attempt)
            with self.subTest(consumer="auditor", splice=label):
                with self.assertRaises(auditor.AuditError):
                    auditor.replay(
                        manifest, supervision_mode="conditioned",
                        expected_v19_attempt=attempt)
            with self.subTest(consumer="census", splice=label):
                with self.assertRaises(census.CensusError):
                    compare_census(retained_raw, attempt, witnesses)

    def test_one_canonical_raw_is_accepted_by_all_consumers(self) -> None:
        raw = canonical_raw()
        attempt = expected_attempt_one()
        (manifest, retained_raw, runner_analysis, audit_result,
         census_result) = consume_bundle(raw, attempt)

        self.assertEqual(
            audit_result["schema"],
            "leopard2-main-compare-v19-conditioned-passive-"
            "independent-audit/v1")
        self.assertEqual(census_result["schema"], census.POLICY_SCHEMA_V3)
        self.assertEqual(
            audit_result["acquisition_generation"], "passive-v3")
        self.assertEqual(
            census_result["acquisition_generation"], "passive-v3")
        self.assertEqual(
            audit_result["pair_qualification"],
            census_result["pair_qualification"])
        self.assertEqual(
            manifest["pair_qualification"], retained_raw["pair_qualification"])
        self.assertEqual(
            retained_raw["digest"], audit_result["raw"]["payload_digest"])
        self.assertEqual(runner_analysis, retained_raw["analysis"])
        self.assertEqual(
            runner_analysis, audit_result["replay"]["analysis"])
        self.assertEqual(
            (audit_result["contract"]["cpu"],
             audit_result["contract"]["sibling"]), (1, 65))
        self.assertEqual(
            census_result["pair_qualification"]["selected_pair"], {
                "benchmark_cpu": 1, "reserved_sibling": 65})
        self.assertTrue(audit_result["gate_results"]
                        ["pair_qualification_fixed_point_independently_replayed"])
        self.assertTrue(census_result["pair_qualification"]
                        ["fixed_point_replayed_independently"])
        self.assertFalse(audit_result["timing_performed"])
        self.assertFalse(audit_result["benchmark_executed"])
        for result in (audit_result, census_result):
            self.assertFalse(result["cpu_pair_exclusive"])
            self.assertFalse(result["causal_performance_claim_eligible"])
            self.assertFalse(result["promotion_eligible"])
            self.assertFalse(result["promotion_passed"])

    def test_runner_cli_binds_attempt_one_and_two_for_v19_replay(self) -> None:
        for attempt_number, attempt in (
                (1, expected_attempt_one()),
                (2, expected_attempt_two())):
            with self.subTest(kind="complete", attempt=attempt_number), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-v19-runner-cli-complete-") as temporary:
                raw = canonical_raw()
                raw["pair_qualification"]["attempt"] = copy.deepcopy(attempt)
                raw = replay_fixtures.fixtures.resign(raw)
                manifest = replay_fixtures.materialize_windowed_auditor_bundle(
                    Path(temporary), raw)
                completed = run_runner_cli(
                    "verify", "--manifest", str(manifest),
                    "--no-current-input-check",
                    *attempt_cli_arguments(attempt_number))
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(
                    completed.stdout,
                    "exact-main ABBA bundle structurally verified only; "
                    "current build/source closure was not revalidated\n")
                self.assertEqual(completed.stderr, "")

            with self.subTest(kind="failure", attempt=attempt_number), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-v19-runner-cli-failure-") as temporary:
                root = Path(temporary)
                failure = replay_fixtures.fixtures.\
                    synthetic_unselected_v19_failure(acquired=False)
                failure["pair_qualification"]["attempt"] = copy.deepcopy(
                    attempt)
                failure = replay_fixtures.fixtures.resign(failure)
                failure_path = root / "failure.json"
                runner.write_json_exclusive(failure_path, failure)
                completed = run_runner_cli(
                    "verify-failure", "--failure", str(failure_path),
                    *attempt_cli_arguments(attempt_number))
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(
                    completed.stdout,
                    "failed exact-main campaign diagnostics and retained "
                    "files verified; this is not valid performance evidence\n")
                self.assertEqual(completed.stderr, "")

    def test_runner_cli_rejects_missing_mismatched_and_legacy_authority(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-runner-cli-reject-") as temporary:
            root = Path(temporary)
            v19_manifest = \
                replay_fixtures.materialize_windowed_auditor_bundle(
                    root / "v19", canonical_raw())
            v19_manifest.parent.chmod(0o700)
            missing = run_runner_cli(
                "verify", "--manifest", str(v19_manifest),
                "--no-current-input-check")
            self.assertEqual(missing.returncode, 1)
            self.assertEqual(missing.stdout, "")
            self.assertEqual(
                missing.stderr,
                "main comparison evidence error: v19 attempt authority must "
                "be supplied exactly for a v19 manifest\n")

            mismatched = run_runner_cli(
                "verify", "--manifest", str(v19_manifest),
                "--no-current-input-check", *attempt_cli_arguments(2))
            self.assertEqual(mismatched.returncode, 1)
            self.assertEqual(mismatched.stdout, "")
            self.assertIn("v19 attempt differs", mismatched.stderr)

            v19_failure = replay_fixtures.fixtures.\
                synthetic_unselected_v19_failure(acquired=False)
            v19_failure_path = root / "v19-failure" / "failure.json"
            v19_failure_path.parent.mkdir(mode=0o700)
            runner.write_json_exclusive(v19_failure_path, v19_failure)
            missing_failure = run_runner_cli(
                "verify-failure", "--failure", str(v19_failure_path))
            self.assertEqual(missing_failure.returncode, 1)
            self.assertEqual(missing_failure.stdout, "")
            self.assertEqual(
                missing_failure.stderr,
                "main comparison evidence error: v19 attempt authority must "
                "be supplied exactly for a v19 failure\n")

            incomplete_attempt_failure = run_runner_cli(
                "verify-failure", "--failure", str(v19_failure_path),
                "--v19-attempt", "2")
            self.assertEqual(incomplete_attempt_failure.returncode, 1)
            self.assertEqual(incomplete_attempt_failure.stdout, "")
            self.assertEqual(
                incomplete_attempt_failure.stderr,
                "main comparison evidence error: v19 attempt authority is "
                "invalid: prior attempt failure hash is not exact lowercase "
                "SHA-256\n")

            mismatched_failure = run_runner_cli(
                "verify-failure", "--failure", str(v19_failure_path),
                *attempt_cli_arguments(2))
            self.assertEqual(mismatched_failure.returncode, 1)
            self.assertEqual(mismatched_failure.stdout, "")
            self.assertIn("v19 attempt differs", mismatched_failure.stderr)

            prior_without_attempt = run_runner_cli(
                "verify", "--manifest", str(v19_manifest),
                "--no-current-input-check",
                "--v19-prior-failure-sha256", "0" * 64)
            self.assertEqual(prior_without_attempt.returncode, 1)
            self.assertEqual(prior_without_attempt.stdout, "")
            self.assertEqual(
                prior_without_attempt.stderr,
                "main comparison evidence error: v19 prior-attempt flags "
                "require --v19-attempt\n")

            v18_raw = replay_fixtures.synthetic_windowed_auditor_raw(
                runner.RAW_SCHEMA_V18)
            v18_manifest = \
                replay_fixtures.materialize_windowed_auditor_bundle(
                    root / "v18", v18_raw)
            v18_manifest.parent.chmod(0o700)
            legacy_complete = run_runner_cli(
                "verify", "--manifest", str(v18_manifest),
                "--no-current-input-check", "--v19-attempt", "1")
            self.assertEqual(legacy_complete.returncode, 1)
            self.assertEqual(legacy_complete.stdout, "")
            self.assertEqual(
                legacy_complete.stderr,
                "main comparison evidence error: v19 attempt authority must "
                "be supplied exactly for a v19 manifest\n")

            v18_failure = replay_fixtures.fixtures.synthetic_failure(
                runner.RAW_SCHEMA_V18)
            v18_failure_path = root / "v18-failure" / "failure.json"
            v18_failure_path.parent.mkdir(mode=0o700)
            runner.write_json_exclusive(v18_failure_path, v18_failure)
            legacy_failure = run_runner_cli(
                "verify-failure", "--failure", str(v18_failure_path),
                "--v19-attempt", "1")
            self.assertEqual(legacy_failure.returncode, 1)
            self.assertEqual(legacy_failure.stdout, "")
            self.assertEqual(
                legacy_failure.stderr,
                "main comparison evidence error: v19 attempt authority must "
                "be supplied exactly for a v19 failure\n")

            malformed_pair_options = argparse.Namespace(
                v19_attempt=2,
                v19_prior_failure_sha256="0" * 64,
                v19_prior_acquisition_sha256="1" * 64,
                v19_prior_selection_status="selected-lowest-primary",
                v19_prior_selected_pair=object(),
            )
            with self.assertRaisesRegex(
                    runner.EvidenceError,
                    "must be BENCHMARK,SIBLING"):
                runner.v19_attempt_from_options(malformed_pair_options)

            with self.assertRaisesRegex(
                    runner.EvidenceError,
                    "supplied exactly for a v19 raw bundle"):
                runner.validate_raw(
                    replay_fixtures.fixtures.resign(v18_raw), None,
                    check_files=False,
                    check_current_inputs=False,
                    expected_v19_attempt=expected_attempt_one())

    def test_all_consumers_reject_common_authority_and_continuity_splices(
            self) -> None:
        pristine = canonical_raw()
        attempt = expected_attempt_one()
        witnesses = census_witnesses(pristine)

        self.assert_all_reject(
            pristine, expected_attempt_two(), "external authority argument",
            witnesses)

        def raw_attempt(raw: dict) -> None:
            raw["pair_qualification"]["attempt"] = expected_attempt_two()

        def policy(raw: dict) -> None:
            raw["pair_qualification"]["policy"] \
                ["candidate_primary_cpus"].pop()

        def pair(raw: dict) -> None:
            raw["pair_qualification"]["selected_pair"] \
                ["reserved_sibling"] = 66

        def campaign_pair(raw: dict) -> None:
            raw["campaign"]["reserved_sibling"] = 66

        def bridge(raw: dict) -> None:
            raw["pair_qualification"]["bridge"] \
                ["bridge_tail_sha256"] = "0" * 64

        def handoff(raw: dict) -> None:
            raw["pair_qualification"]["first_window_handoff"] \
                ["handoff_elapsed_ns"] = 0

        def lineage(raw: dict) -> None:
            raw["pair_qualification"]["v18_failure_lineage"] \
                ["attempts"][0]["envelope_sha256sums_sha256"] = "0" * 64

        def host(raw: dict) -> None:
            for name in ("host_initial", "host_final"):
                raw[name]["system"]["release"] = "changed-fixture"

        def host_continuity(raw: dict) -> None:
            raw["host_final"]["system"]["release"] = "changed-fixture"

        def claim(raw: dict) -> None:
            raw["pair_qualification"]["shared_host_claim_ceiling"] \
                ["promotion_eligible"] = True

        def status(raw: dict) -> None:
            raw["pair_qualification"]["record_status"] = "failed"

        def acquisition_campaign(raw: dict) -> None:
            raw["campaign"]["allowed_cpu_set_at_launch"].append(128)

        def bridge_isolation(raw: dict) -> None:
            rebuild_isolation(
                raw,
                before_ns=raw["isolation"]["before"]["monotonic_ns"] + 1)

        def pair_lease(raw: dict) -> None:
            lease = raw["isolation"]["pair_lease"]
            lease["payload"]["cpus"] = [1, 66]
            uid = lease["payload"]["uid"]
            lease["path"] = (
                f"/run/user/{uid}/leopard2-cpu-leases/"
                f"leopard2-cpu-pair-{uid}-1-66.lock")
            lease["sha256"] = runner.sha256_bytes(
                runner.canonical_bytes(lease["payload"]))

        mutations = (
            ("raw attempt authority", raw_attempt),
            ("qualification policy", policy),
            ("selected pair", pair),
            ("campaign selected pair", campaign_pair),
            ("qualification bridge", bridge),
            ("first-window handoff", handoff),
            ("v18 failure lineage", lineage),
            ("host identity", host),
            ("host endpoint continuity", host_continuity),
            ("claim ceiling", claim),
            ("qualification status", status),
            ("acquisition/campaign affinity", acquisition_campaign),
            ("selected pair/pair lease", pair_lease),
            ("bridge/isolation join", bridge_isolation),
            ("handoff/first-window join", splice_first_window),
        )
        for label, mutation in mutations:
            candidate = copy.deepcopy(pristine)
            mutation(candidate)
            self.assert_all_reject(
                candidate, attempt, label, witnesses)

    def test_dormant_live_producer_uses_one_dynamic_pair_end_to_end(
            self) -> None:
        raw = replay_fixtures.fixtures.synthetic_raw(
            raw_schema=runner.RAW_SCHEMA_V19, v19_primary=7)
        qualification = raw["pair_qualification"]
        self.assertEqual(qualification["selected_pair"], {
            "benchmark_cpu": 7, "reserved_sibling": 71,
        })

        forbidden = (
            "run_process_bounded", "validate_topology", "host_identity",
            "run_child",
        )
        patches = [
            mock.patch.object(
                runner, name,
                side_effect=AssertionError(
                    f"dormant v19 producer called {name}"))
            for name in forbidden
        ]
        with patches[0], patches[1], patches[2], patches[3]:
            authority = live_authority(raw)
            self.assertEqual(
                authority["acquisition_generation"], "passive-v3")
            self.assertEqual(
                authority["qualification"]["stage"], "bridged")
            handoff = runner.bind_v19_live_first_window(
                authority, campaign=raw["campaign"],
                host_identity_value=raw["host_initial"],
                first_window_before_value=raw["invocations"][0]
                    ["cpu_window"]["before"])
            launch = runner.require_v19_live_launch_authority(
                handoff, campaign=raw["campaign"],
                host_identity_value=raw["host_initial"])
            completed = runner.complete_v19_live_run_authority(
                launch, campaign=raw["campaign"],
                host_identity_value=raw["host_initial"])
            rebuilt_raw = runner.v19_live_raw_record(
                completed, created_utc=raw["created_utc"],
                campaign=raw["campaign"], host_initial=raw["host_initial"],
                isolation=raw["isolation"],
                reservation=raw["reservation"],
                input_specification=raw["input_specification"],
                identities_initial=raw["identities_initial"],
                executable_snapshots=raw["executable_snapshots"],
                invocations=raw["invocations"],
                identities_final=raw["identities_final"],
                host_final=raw["host_final"])
        self.assertEqual(rebuilt_raw, raw)
        self.assertEqual(rebuilt_raw["schema"], runner.RAW_SCHEMA_V19)
        self.assertNotEqual(
            rebuilt_raw["schema"], runner.RAW_SCHEMA)
        self.assertEqual(
            rebuilt_raw["isolation"]["before"]["monotonic_ns"],
            qualification["bridge"]["campaign_presample_before"]
                ["read_finished_monotonic_ns"])
        self.assertEqual(
            rebuilt_raw["invocations"][0]["cpu_window"]["before"],
            qualification["first_window_handoff"]["first_window_before"])

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-live-producer-") as temporary:
            manifest_path = \
                replay_fixtures.materialize_windowed_auditor_bundle(
                    Path(temporary), rebuilt_raw)
            stored_manifest = runner.strict_json_loads(
                manifest_path.read_bytes(), "stored v19 live manifest")
            stored_raw = runner.strict_json_loads(
                (manifest_path.parent / "raw.json").read_bytes(),
                "stored v19 live raw")
            self.assertIsInstance(stored_manifest, dict)
            self.assertIsInstance(stored_raw, dict)
            rebuilt_manifest = runner.v19_live_manifest_record(
                completed, stored_raw,
                created_utc=stored_manifest["created_utc"])
            self.assertEqual(rebuilt_manifest, stored_manifest)
            self.assertEqual(
                rebuilt_manifest["schema"], runner.MANIFEST_SCHEMA_V19)
            self.assertNotEqual(
                rebuilt_manifest["schema"], runner.MANIFEST_SCHEMA)

        failed_authority = runner.fail_v19_live_campaign_authority(
            handoff, campaign=raw["campaign"],
            host_identity_value=raw["host_initial"])
        before = raw["isolation"]["before"]
        after = raw["isolation"]["after"]
        failed_isolation = runner.isolation_record_v2(
            7, 71, raw["isolation"]["pair_lease"],
            before["monotonic_ns"], after["monotonic_ns"],
            before["benchmark_cpu"], after["benchmark_cpu"],
            before["reserved_sibling"], after["reserved_sibling"], [])
        failure = runner.v19_live_failure_record(
            failed_authority, created_utc=raw["created_utc"],
            error_type="EvidenceError", error="fixture failure",
            diagnostic_traceback="fixture traceback",
            campaign=raw["campaign"], host_initial=raw["host_initial"],
            reservation=raw["reservation"],
            pair_lease=raw["isolation"]["pair_lease"],
            isolation=failed_isolation,
            input_specification=raw["input_specification"],
            identities_initial=raw["identities_initial"],
            executable_snapshots=raw["executable_snapshots"],
            invocations=[], retained_files=[])
        self.assertEqual(
            failure,
            replay_fixtures.fixtures.synthetic_failure(
                runner.RAW_SCHEMA_V19, v19_primary=7))
        self.assertEqual(failure["schema"], runner.FAILURE_SCHEMA_V19)
        self.assertNotEqual(failure["schema"], runner.FAILURE_SCHEMA)
        self.assertEqual(
            failure["pair_qualification"]["selected_pair"], {
                "benchmark_cpu": 7, "reserved_sibling": 71,
            })
        self.assertEqual(
            (failure["pair_qualification"]["stage"],
             failure["pair_qualification"]["terminal"]),
            ("campaign", "campaign-rejected"))

    def test_dormant_live_prelaunch_and_handoff_fail_closed(self) -> None:
        raw = canonical_raw()
        qualification = raw["pair_qualification"]

        for label, mutate in (
            ("campaign pair", lambda values:
             values["campaign"].__setitem__("benchmark_cpu", 2)),
            ("boolean primary", lambda values:
             values["campaign"].__setitem__("benchmark_cpu", True)),
            ("float primary", lambda values:
             values["campaign"].__setitem__("benchmark_cpu", 1.0)),
            ("float sibling", lambda values:
             values["campaign"].__setitem__("reserved_sibling", 65.0)),
            ("launch affinity", lambda values:
             values["campaign"]["allowed_cpu_set_at_launch"].pop()),
            ("bridge", lambda values:
             values["bridge"].__setitem__("bridge_tail_sha256", "0" * 64)),
            ("attempt", lambda values:
             values["attempt"].__setitem__("attempt_budget", 3)),
            ("v18 lineage", lambda values:
             values.__setitem__("lineage", "0" * 64)),
            ("claim ceiling", lambda values:
             values.__setitem__("claim", "0" * 64)),
            ("resource envelope", lambda values:
             values["resource"].__setitem__("release_max_jobs", 2)),
            ("controller generation", lambda values:
             values["controller"].__setitem__(
                 "acquisition_generation", "passive-v2")),
            ("missing census wrapper", lambda values:
             values["controller"].__setitem__("wrapper_pid", 123456)),
            ("census wrapper affinity", lambda values:
             values["census"]["same_uid_processes"]["entries"][0]
                 .__setitem__("cpus_allowed", [0])),
            ("pre census", lambda values:
             values["census"]["reserved_cpus"].append(66)),
        ):
            values = {
                "campaign": copy.deepcopy(raw["campaign"]),
                "host": copy.deepcopy(raw["host_initial"]),
                "attempt": copy.deepcopy(qualification["attempt"]),
                "acquisition": copy.deepcopy(qualification["acquisition"]),
                "bridge": copy.deepcopy(qualification["bridge"]),
                "lineage": qualification["v18_failure_lineage_sha256"],
                "claim": runner.pair_v19.contract.canonical_sha256(
                    qualification["shared_host_claim_ceiling"]),
                "resource": copy.deepcopy(
                    runner.V19_LIVE_RESOURCE_ENVELOPE),
                "controller": passive_fixtures.v19_controller(raw),
                "census": passive_fixtures.v19_snapshot(raw, "pre"),
            }
            mutate(values)
            passive_fixtures.reseal_snapshot(values["census"])
            with self.subTest(splice=label), self.assertRaises(
                    runner.EvidenceError):
                runner.prepare_v19_live_run_authority(
                    campaign=values["campaign"],
                    host_identity_value=values["host"],
                    expected_attempt=values["attempt"],
                    acquisition_value=values["acquisition"],
                    bridge_value=values["bridge"],
                    expected_v18_failure_lineage_sha256=values["lineage"],
                    expected_claim_ceiling_sha256=values["claim"],
                    resource_envelope_value=values["resource"],
                    controller_affinity_value=values["controller"],
                    environment_census_pre_value=values["census"])

        authority = live_authority(raw)
        tail = qualification["bridge"]["campaign_presample_before"]
        late_before = copy.deepcopy(
            raw["invocations"][0]["cpu_window"]["before"])
        late_started = (
            tail["read_finished_monotonic_ns"] +
            runner.pair_v19.MAXIMUM_HANDOFF_ELAPSED_NS + 1)
        late_before.update({
            "read_started_monotonic_ns": late_started,
            "read_finished_monotonic_ns": late_started + 1_000_000,
            "monotonic_ns": late_started + 1_000_000,
        })
        rejected = runner.bind_v19_live_first_window(
            authority, campaign=raw["campaign"],
            host_identity_value=raw["host_initial"],
            first_window_before_value=late_before)
        self.assertEqual(
            rejected["qualification"]["terminal"],
            "first-window-handoff-late")
        launched = False
        with self.assertRaises(runner.EvidenceError):
            runner.require_v19_live_launch_authority(
                rejected, campaign=raw["campaign"],
                host_identity_value=raw["host_initial"])
            launched = True
        self.assertFalse(launched)

        mixed = copy.deepcopy(authority)
        mixed.pop("digest")
        mixed["acquisition_generation"] = "passive-v2"
        mixed = runner.signed(mixed)
        with self.assertRaises(runner.EvidenceError):
            runner.bind_v19_live_first_window(
                mixed, campaign=raw["campaign"],
                host_identity_value=raw["host_initial"],
                first_window_before_value=raw["invocations"][0]
                    ["cpu_window"]["before"])

        self.assertEqual(
            (runner.RAW_SCHEMA, runner.MANIFEST_SCHEMA,
             runner.FAILURE_SCHEMA),
            (runner.RAW_SCHEMA_V18, runner.MANIFEST_SCHEMA_V18,
             runner.FAILURE_SCHEMA_V18))

    def test_dormant_live_authority_rechecks_handoff_census(self) -> None:
        raw = canonical_raw()
        launch = handed_off_authority(raw)
        arguments = {
            "campaign": raw["campaign"],
            "host_identity_value": raw["host_initial"],
        }
        completed = runner.complete_v19_live_run_authority(launch, **arguments)
        failed = runner.fail_v19_live_campaign_authority(launch, **arguments)
        first = launch["qualification"]["first_window_handoff"]
        tail = launch["qualification"]["bridge"]["campaign_presample_before"]
        cpu = raw["campaign"]["benchmark_cpu"]

        for authority in (launch, completed, failed):
            for splice in ("late census", "before bridge", "after window"):
                changed = copy.deepcopy(authority)
                pre = changed["environment_census_pre"]
                if splice == "late census":
                    started = first["first_window_before_read_started_monotonic_ns"]
                    pre["scan_started_monotonic_ns"] = started + 1
                    pre["scan_finished_monotonic_ns"] = started + 2
                    pre["activity_boundary_monotonic_ns"] = started + 3
                else:
                    counter = pre["proc_stat"][str(cpu)]
                    if splice == "before bridge":
                        idle = tail["cpus"][str(cpu)]["fields"]["idle"] - 1
                    else:
                        idle = first["first_window_before"]["benchmark_cpu"] \
                            ["fields"]["idle"] + 1
                    counter["fields"]["idle"] = idle
                    counter["total_jiffies"] = sum(
                        counter["fields"][name]
                        for name in runner.CPU_STAT_FIELDS)
                passive_fixtures.reseal_snapshot(pre)
                changed["environment_census_pre_sha256"] = \
                    runner.sha256_bytes(runner.canonical_bytes(pre))
                changed.pop("digest")
                changed = runner.signed(changed)
                # The census and all hashes are individually valid; only
                # their relationship with the retained handoff is false.
                census.validate_snapshot(
                    pre, "pre", cpu_pair=(cpu, raw["campaign"]["reserved_sibling"]))
                with self.subTest(
                        stage=authority["qualification"]["stage"], splice=splice
                ), self.assertRaisesRegex(runner.EvidenceError, "pre-census"):
                    if authority is launch:
                        runner.require_v19_live_launch_authority(changed, **arguments)
                    else:
                        runner._validate_v19_live_run_authority(
                            changed, **arguments)

    def test_isolated_consumer_self_test_stdout_is_exact(self) -> None:
        cases = (
            (MAIN_COMPARE / "audit_v17_gfni_main_compare.py", "--self-test",
             "v17 GFNI exact-main independent auditor self-test passed\n"),
            (MAIN_COMPARE / "passive_environment_census.py", "self-test",
             "passive environment census self-test passed\n"),
        )
        for script, argument, expected_stdout in cases:
            for optimized in (False, True):
                command = [sys.executable]
                if optimized:
                    command.append("-O")
                command.extend(("-I", "-S", "-B", str(script), argument))
                completed = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=False, text=True)
                with self.subTest(script=script.name, optimized=optimized):
                    self.assertEqual(
                        completed.returncode, 0, completed.stderr)
                    self.assertEqual(completed.stdout, expected_stdout)


if __name__ == "__main__":
    unittest.main()
