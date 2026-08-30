#!/usr/bin/env python3
"""Host-independent tests for the K65 generation-3 preregistration contract."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "k65_gen3_preregistration.py"
TEMPLATE_PATH = HERE / \
    "k65r65_b64_packed_terminal_gen3_preregistration.template.json"
PAIR_TEST_PATH = HERE.parent / "main_compare" / \
    "test_pair_qualification_contract.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


prereg = load_module("k65_gen3_preregistration_under_test", MODULE_PATH)
fixtures = load_module("k65_gen3_pair_fixtures", PAIR_TEST_PATH)


def retained_generation2_files_available() -> bool:
    lineage = prereg.closed_generation2_lineage_record(verify_files=False)
    paths = [prereg.REPO_ROOT / lineage["runner_path"]]
    for attempt in lineage["attempts"]:
        paths.extend((
            prereg.REPO_ROOT / attempt["terminal_path"],
            prereg.REPO_ROOT / attempt["failure_path"],
        ))
    paths.extend(
        prereg.REPO_ROOT / item["path"] for item in lineage["lineage_files"])
    return all(path.is_file() for path in paths)


RETAINED_GENERATION2_FILES_AVAILABLE = retained_generation2_files_available()


def final_preregistration(
    *, track_b_permitted: bool = True, scan_window_count: int = 2,
) -> dict:
    return prereg.preregistration_record(
        authority="unit-test-fixture-not-an-authorization",
        authorized_utc="2026-08-30T00:00:00Z",
        clock_ticks_per_second=100,
        candidate_primary_cpus=[7, 10, 30],
        excluded_pairs=[{
            "benchmark_cpu": 30,
            "reserved_sibling": 40,
            "reason": "fixture policy exclusion",
        }],
        track_b_permitted=track_b_permitted,
        setup_invalid_budget=5,
        environment_rejected_budget=8,
        evidence_attempt_budget=3,
        scan_window_count=scan_window_count,
        scan_nominal_window_ns=1_000_000_000,
        bridge_minimum_window_count=5,
        bridge_maximum_window_count=120,
        bridge_nominal_window_ns=1_000_000_000,
        maximum_handoff_elapsed_ns=120_000_000_000,
        freeze_point="armed",
        candidate_executable_mode="deferred-until-arming",
        candidate_executable_sha256=None,
        candidate_source_commit="1" * 40,
        candidate_source_tree="2" * 40,
        controller_bindings=[
            {"path": path, "sha256": f"{index + 1:064x}"}
            for index, path in enumerate(prereg.REQUIRED_CONTROLLER_PATHS)
        ],
        verify_files=False,
    )


def policy_scans(*, events: dict[int, str] | None = None,
                 frozen: dict | None = None) -> tuple[dict, dict]:
    policy_a = fixtures.policy_fixture("pair-and-domain")
    policy_b = fixtures.policy_fixture("pair-only")
    endpoints = fixtures.snapshots(events=events)
    return (
        fixtures.scan_fixture(
            policy=policy_a, endpoint_records=endpoints, frozen=frozen),
        fixtures.scan_fixture(
            policy=policy_b, endpoint_records=endpoints, frozen=frozen),
    )


class K65Generation3PreregistrationTests(unittest.TestCase):
    def assertRejected(self, function, pattern: str | None = None) -> None:
        with self.assertRaises((
            prereg.PreregistrationError,
            prereg.contract.QualificationError,
        )) as raised:
            function()
        if pattern is not None:
            self.assertIn(pattern, str(raised.exception))

    def test_committed_template_is_exact_and_never_executable(self) -> None:
        observed = prereg.load_preregistration_template(
            TEMPLATE_PATH.read_bytes(), verify_files=False)
        expected = prereg.preregistration_template_record(verify_files=False)
        self.assertEqual(observed, expected)
        self.assertFalse(observed["safe_to_execute"])
        self.assertTrue(all(
            decision["authorized_value"] is None
            for decision in observed["required_decisions"]
        ))
        self.assertRejected(
            lambda: prereg.validate_preregistration(observed, verify_files=False),
            "preregistration")
        self.assertEqual(
            prereg.contract.canonical_json_bytes(observed),
            TEMPLATE_PATH.read_bytes())

    @unittest.skipUnless(
        RETAINED_GENERATION2_FILES_AVAILABLE,
        "retained K65 generation-2 lanes are unavailable")
    def test_closed_generation2_lineage_replays_all_exact_files(self) -> None:
        lineage = prereg.closed_generation2_lineage_record(verify_files=True)
        self.assertTrue(lineage["generation_exhausted"])
        self.assertEqual(lineage["valid_ratio_count"], 0)
        self.assertEqual(
            [attempt["attempt"] for attempt in lineage["attempts"]],
            [1, 2, 3])
        self.assertEqual(
            [attempt["failure_class"] for attempt in lineage["attempts"]],
            ["environment_rejected", "setup_invalid",
             "contaminated_in_campaign"])
        self.assertEqual(len(lineage["lineage_files"]), 2)
        self.assertEqual(
            prereg.validate_closed_generation2_lineage(
                lineage, verify_files=True),
            lineage)

    def test_structural_template_validation_is_checkout_portable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            observed = prereg.load_preregistration_template(
                TEMPLATE_PATH.read_bytes(), root, verify_files=False)
            self.assertFalse(observed["safe_to_execute"])
            self.assertRejected(lambda: prereg.closed_generation2_lineage_record(
                root, verify_files=True), "not retained")

    def test_lineage_mutations_and_incomplete_attempts_reject(self) -> None:
        lineage = prereg.closed_generation2_lineage_record(verify_files=False)
        mutations = []
        missing = copy.deepcopy(lineage)
        missing["attempts"].pop()
        mutations.append(missing)
        relabeled = copy.deepcopy(lineage)
        relabeled["attempts"][1]["attempt"] = 4
        mutations.append(relabeled)
        promoted = copy.deepcopy(lineage)
        promoted["valid_ratio_count"] = 1
        mutations.append(promoted)
        changed_hash = copy.deepcopy(lineage)
        changed_hash["attempts"][0]["terminal_sha256"] = "0" * 64
        mutations.append(changed_hash)
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(lambda mutation=mutation: (
                    prereg.validate_closed_generation2_lineage(
                        mutation, verify_files=False)))

    def test_final_fixed_point_binds_policy_hashes_budgets_and_zero_claims(self) -> None:
        record = final_preregistration(track_b_permitted=True)
        self.assertEqual(
            prereg.validate_preregistration(record, verify_files=False), record)
        qualification = record["qualification"]
        self.assertEqual(
            qualification["policy_evaluation_order"],
            ["pair-and-domain", "pair-only"])
        self.assertEqual(
            qualification["policy_sha256s"],
            [prereg.contract.canonical_sha256(policy)
             for policy in qualification["policies"]])
        self.assertEqual(record["budgets"]["evidence_attempts"], 3)
        self.assertEqual(record["budgets"]["evidence_budget_commit_state"],
                         "ARMED")
        self.assertTrue(
            record["budgets"]["child_launch_requires_prior_budget_commit"])
        self.assertEqual(record["campaign"]["expected_child_process_count"], 1650)
        expected_claim_keys = {
            "promotion_eligible", "host_exclusivity_proved",
            "whole_campaign_interval_observed",
            "causal_performance_claim_allowed",
        }
        for track in ("track_a_claim_ceiling", "track_b_claim_ceiling"):
            self.assertEqual(set(record["reporting"][track]), expected_claim_keys)
            self.assertEqual(
                set(record["reporting"][track].values()), {False})
        self.assertEqual(record["reporting"]["track_a_ratio_field"],
                         "conditioned_nominal_ratio")
        self.assertEqual(record["reporting"]["track_b_ratio_field"],
                         "conditioned_nominal_ratio")
        decision_names = {
            item["name"]
            for item in prereg.preregistration_template_record(
                verify_files=False)["required_decisions"]
        }
        self.assertIn("clock_ticks_per_second", decision_names)
        self.assertIn("excluded_pairs", decision_names)
        self.assertEqual(
            decision_names,
            set(prereg._PREREGISTRATION_PARAMETER_DECISIONS.values()))

    def test_current_controller_binding_helper_covers_exact_reviewed_files(self) -> None:
        bindings = prereg.current_controller_bindings()
        self.assertEqual(
            [item["path"] for item in bindings],
            sorted(prereg.REQUIRED_CONTROLLER_PATHS))
        self.assertEqual(len({item["sha256"] for item in bindings}),
                         len(bindings))
        for item in bindings:
            self.assertRegex(item["sha256"], r"^[0-9a-f]{64}$")

    def test_policy_b_is_retained_but_not_consulted_when_a_selects(self) -> None:
        record = final_preregistration(track_b_permitted=True)
        scan_a, scan_b = policy_scans()
        selection = prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=None)
        self.assertEqual(selection["selected_track"], "pair-and-domain")
        self.assertEqual(selection["selected_pair"], {"benchmark_cpu": 7,
                                                       "reserved_sibling": 3})
        self.assertTrue(selection["policy_b_retained"])
        self.assertFalse(selection["policy_b_consulted"])
        self.assertEqual(set(selection["claim_ceiling"].values()), {False})

    def test_policy_b_is_consulted_only_after_a_has_no_candidate(self) -> None:
        record = final_preregistration(track_b_permitted=True)
        scan_a, scan_b = policy_scans(events={30: "user"})
        self.assertEqual(scan_a["selection_status"],
                         "no-candidate-pair-qualified")
        self.assertEqual(scan_b["selection_status"], "selected-lowest-primary")
        selection = prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=None)
        self.assertTrue(selection["policy_b_consulted"])
        self.assertEqual(selection["selected_track"], "pair-only")
        self.assertEqual(selection["selected_pair"], {"benchmark_cpu": 7,
                                                       "reserved_sibling": 3})

    def test_policy_a_and_b_must_share_the_exact_snapshot_chain(self) -> None:
        record = final_preregistration(track_b_permitted=True)
        scan_a, unused_scan_b = policy_scans()
        shifted = []
        for index in range(3):
            started = fixtures.BASE_NS + index * 250_100_000 + 10_000
            shifted.append(fixtures.snapshot(
                index, started_ns=started,
                finished_ns=started + fixtures.READ_NS))
        scan_b = fixtures.scan_fixture(
            policy=fixtures.policy_fixture("pair-only"),
            endpoint_records=shifted)
        self.assertRejected(lambda: prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=None), "share one retained snapshot chain")

    def test_scan_window_count_must_match_the_preregistration(self) -> None:
        record = final_preregistration(
            track_b_permitted=True, scan_window_count=3)
        scan_a, scan_b = policy_scans()
        self.assertEqual(scan_a["scan_window_count"], 2)
        self.assertRejected(lambda: prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=None), "scan window count")

    def test_frozen_pair_failure_never_falls_back_to_policy_b(self) -> None:
        record = final_preregistration(track_b_permitted=True)
        frozen = {"benchmark_cpu": 10, "reserved_sibling": 20}
        scan_a, scan_b = policy_scans(events={10: "system"}, frozen=frozen)
        selection = prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=frozen)
        self.assertEqual(selection["selection_status"],
                         "frozen-pair-did-not-requalify")
        self.assertFalse(selection["policy_b_consulted"])
        self.assertIsNone(selection["selected_pair"])

    def test_track_b_is_absent_when_not_authorized(self) -> None:
        record = final_preregistration(track_b_permitted=False)
        scan_a, scan_b = policy_scans(events={30: "user"})
        selection = prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=None,
            expected_frozen_pair=None)
        self.assertFalse(selection["policy_b_retained"])
        self.assertFalse(selection["policy_b_consulted"])
        self.assertIsNone(selection["selected_pair"])
        self.assertRejected(lambda: prereg.qualification_track_selection_record(
            record, policy_a_scan=scan_a, policy_b_scan=scan_b,
            expected_frozen_pair=None))

    def test_no_pooling_denylist_and_lane_overlap_are_fail_closed(self) -> None:
        record = final_preregistration(track_b_permitted=False)
        denylist = record["no_pooling"]["denylisted_artifacts"]
        self.assertEqual(len(denylist), 10)
        for item in denylist:
            with self.subTest(label=item["label"]):
                self.assertRejected(lambda item=item: (
                    prereg.reject_denylisted_evidence_hash(
                        item["sha256"], record)), "generation-2 artifact")
        prereg.reject_denylisted_evidence_hash("f" * 64, record)
        forged = copy.deepcopy(record)
        forged["no_pooling"]["denylisted_artifacts"] = []
        self.assertRejected(lambda: prereg.reject_denylisted_evidence_hash(
            "f" * 64, forged))

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake_repo = root / "repo"
            fake_repo.mkdir()
            empty = root / "empty"
            empty.mkdir()
            self.assertEqual(
                prereg.validate_output_lane_for_preregistration(
                    empty, record, fake_repo),
                empty.resolve())
            (empty / "occupied").write_text("x", encoding="ascii")
            self.assertRejected(lambda: (
                prereg.validate_output_lane_for_preregistration(
                    empty, record, fake_repo)),
                "not empty")
            target = root / "target"
            target.mkdir()
            link = root / "link"
            link.symlink_to(target, target_is_directory=True)
            self.assertRejected(lambda: (
                prereg.validate_output_lane_for_preregistration(
                    link, record, fake_repo)),
                "canonical directory")
            prior = fake_repo / \
                record["closed_generation_2"]["attempts"][0]["lane"]
            prior.mkdir(parents=True)
            self.assertRejected(lambda: (
                prereg.validate_output_lane_for_preregistration(
                    prior, record, fake_repo)), "overlaps")
            absent = root / "absent"
            self.assertRejected(lambda: (
                prereg.validate_output_lane_for_preregistration(
                    absent, record, fake_repo)), "cannot be resolved")
            missing_repo = root / "missing-repo"
            self.assertRejected(lambda: (
                prereg.validate_output_lane_for_preregistration(
                    target, record, missing_repo)), "cannot be resolved")

    def test_strict_json_and_extra_keys_reject(self) -> None:
        record = final_preregistration(track_b_permitted=False)
        encoded = prereg.contract.canonical_json_bytes(record)
        self.assertEqual(prereg.load_preregistration(
            encoded, verify_files=False), record)
        pretty = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
        self.assertRejected(lambda: prereg.load_preregistration(
            pretty, verify_files=False), "not canonical")
        for payload in (b'{}\n{}', b'{"x":1,"x":2}', b'{"x":NaN}'):
            with self.subTest(payload=payload):
                self.assertRejected(lambda payload=payload: (
                    prereg.load_preregistration(payload, verify_files=False)))
        extra = copy.deepcopy(record)
        extra["posthoc"] = False
        self.assertRejected(lambda: prereg.validate_preregistration(
            extra, verify_files=False))


if __name__ == "__main__":
    unittest.main()
