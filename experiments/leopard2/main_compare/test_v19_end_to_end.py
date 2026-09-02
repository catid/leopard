#!/usr/bin/python3

"""Shared end-to-end fixtures for dormant conditioned v19 evidence."""

from __future__ import annotations

import copy
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MAIN_COMPARE = ROOT / "experiments/leopard2/main_compare"


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
