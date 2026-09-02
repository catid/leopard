#!/usr/bin/env python3
"""Host-free adversarial tests for the pair-qualified v19 contract."""

from __future__ import annotations

import ast
import copy
import importlib.util
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "pair_qualified_v19_contract.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "leopard2_pair_qualified_v19_acquire_fixtures",
    HERE / "test_pair_qualification_acquire.py")
live = load_module(
    "leopard2_pair_qualified_v19_bridge_fixtures",
    HERE / "pair_qualification_bridge_acquire.py")
v19 = load_module("leopard2_pair_qualified_v19_contract_test", MODULE_PATH)
contract = v19.contract


def policy_fixture() -> dict:
    return contract.qualification_policy_record(
        clock_ticks_per_second=100, candidate_primary_cpus=[2, 3],
        excluded_pairs=[], domain_mode="pair-only")


def mutate_proc_users(data: bytes, users: dict[int, int]) -> bytes:
    rows = []
    for row in data.splitlines():
        changed = row
        for cpu, user in users.items():
            prefix = f"cpu{cpu} ".encode("ascii")
            if row.startswith(prefix):
                fields = row.split()
                idle = int(fields[4])
                fields[1] = str(user).encode("ascii")
                fields[4] = str(idle - user).encode("ascii")
                changed = b" ".join(fields)
        rows.append(changed)
    return b"\n".join(rows) + b"\n"


def bump_proc_idle(data: bytes, amount: int) -> bytes:
    rows = []
    for row in data.splitlines():
        if any(row.startswith(f"cpu{cpu} ".encode("ascii"))
               for cpu in (2, 3, 6, 7)):
            fields = row.split()
            fields[4] = str(int(fields[4]) + amount).encode("ascii")
            row = b" ".join(fields)
        rows.append(row)
    return b"\n".join(rows) + b"\n"


def acquisition_fixture(
    *, frozen_pair=None, selected: bool = True,
    nonidle_cpus: tuple[int, ...] = (),
):
    indices = [28 * index for index in range(13)]
    records = [
        bump_proc_idle(fixtures.proc_stat(index), endpoint // 6)
        for endpoint, index in enumerate(indices)
    ]
    if not selected:
        nonidle_cpus = (2, 3)
    if nonidle_cpus:
        records[-1] = mutate_proc_users(
            records[-1], {cpu: 1 for cpu in nonidle_cpus})
    records.extend((
        bump_proc_idle(fixtures.proc_stat(340), 2),
        bump_proc_idle(fixtures.proc_stat(344), 2),
    ))
    reader = fixtures.FakeHostReader(proc_records=records)
    policy = policy_fixture()
    acquisition = live.acquire.acquire_pair_qualification(
        reader, policy_value=policy,
        window_count=v19.QUALIFICATION_WINDOW_COUNT,
        nominal_window_ns=v19.QUALIFICATION_NOMINAL_WINDOW_NS,
        frozen_pair_from_prior_attempt=frozen_pair)
    return acquisition, policy, reader


def accepted_bundle(*, frozen_pair=None):
    acquisition, policy, reader = acquisition_fixture(frozen_pair=frozen_pair)
    tail = acquisition["scan"]["windows"][-1]["after"]
    bridge_record = live.acquire_v19_pair_qualification_bridge(
        reader, acquisition, tail, expected_policy=policy,
        expected_policy_sha256=contract.canonical_sha256(policy),
        expected_frozen_pair=frozen_pair,
        expected_acquisition_window_count=v19.QUALIFICATION_WINDOW_COUNT,
        expected_acquisition_nominal_window_ns=
        v19.QUALIFICATION_NOMINAL_WINDOW_NS,
        expected_bridge_geometry=v19.bridge_geometry_record())
    return acquisition, bridge_record, policy


def first_window_before(
    bridge_record: dict, *, elapsed_ns: int = 1_000_000,
    nonidle_role: str | None = None,
) -> dict:
    selected = bridge_record["selected_pair"]
    tail = bridge_record["campaign_presample_before"]
    records = {}
    for role in ("benchmark_cpu", "reserved_sibling"):
        cpu = selected[role]
        record = copy.deepcopy(tail["cpus"][str(cpu)])
        if nonidle_role == role:
            record["fields"]["user"] += 1
        else:
            record["fields"]["idle"] += 1
        record["total_jiffies"] += 1
        records[role] = record
    started = tail["read_finished_monotonic_ns"] + elapsed_ns
    return {
        "benchmark_cpu": records["benchmark_cpu"],
        "monotonic_ns": started + 1_000_000,
        "read_finished_monotonic_ns": started + 1_000_000,
        "read_started_monotonic_ns": started,
        "reserved_sibling": records["reserved_sibling"],
    }


def attempt_one() -> dict:
    return v19.pair_qualified_attempt_record(attempt=1)


def complete_record():
    acquisition, bridge_record, policy = accepted_bundle()
    before = first_window_before(bridge_record)
    record = v19.pair_qualified_v19_record(
        stage="complete", record_status="complete",
        terminal=v19.SUCCESS_TERMINAL, policy_value=policy,
        expected_policy_sha256=contract.canonical_sha256(policy),
        host_identity_sha256="a" * 64, attempt_value=attempt_one(),
        acquisition_value=acquisition, bridge_value=bridge_record,
        first_window_before_value=before)
    return record, acquisition, bridge_record, policy, before


class PairQualifiedV19ContractTest(unittest.TestCase):
    def assertV19Failure(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(contract.QualificationError) as raised:
            function(*args, **kwargs)
        self.assertNotIsInstance(raised.exception, (KeyboardInterrupt, SystemExit))

    def validate(self, record, policy, attempt=None):
        return v19.validate_pair_qualified_v19_record(
            record, expected_policy=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            expected_host_identity_sha256="a" * 64,
            expected_attempt=attempt or attempt_one())

    def test_complete_record_binds_the_full_conditioned_chain(self) -> None:
        record, acquisition, bridge_record, policy, before = complete_record()
        self.assertEqual(self.validate(record, policy), record)
        self.assertEqual(record["qualification_geometry"], {
            "window_count": 12, "nominal_window_ns": 7_000_000_000,
        })
        self.assertEqual(record["selected_pair"], {
            "benchmark_cpu": 2, "reserved_sibling": 6,
        })
        self.assertEqual(record["acquisition_sha256"],
                         contract.canonical_sha256(acquisition))
        self.assertEqual(record["bridge_sha256"],
                         contract.canonical_sha256(bridge_record))
        self.assertEqual(
            record["first_window_handoff"]["first_window_before"], before)
        self.assertTrue(record["first_window_handoff"]["accepted"])
        self.assertEqual(record["terminal"], "NOT_PROMOTED")
        self.assertEqual(
            contract.canonical_sha256(record),
            "e322bdf85c2af666eefb1bff1c2a499cf5d101d944150a5c0bf07594b88d05cf")
        self.assertEqual(record["shared_host_claim_ceiling"], {
            "promotion_eligible": False,
            "host_exclusivity_proved": False,
            "whole_campaign_interval_observed": False,
            "causal_performance_claim_allowed": False,
        })
        self.assertEqual([
            entry["envelope_sha256sums_sha256"]
            for entry in record["v18_failure_lineage"]["attempts"]
        ], [
            "ce65c3a49ef1c1d89ba51ea03d0af4742d6790e6f2ea2662917d9ef9a9d945d7",
            "a1bf0eda157c251f33f7260ebd76931d88054d460bd07a97bcba2811384b2c10",
            "fe5b40cc98753cbd794ee019cb0e2643d0ccee0aca4c5fd7b2e0b27df8a86139",
        ])

    def test_attempt_two_freezes_or_freshly_selects_without_fallback(self) -> None:
        acquisition, unused_bridge, unused_policy = accepted_bundle()
        selected = acquisition["scan"]["selected"]
        frozen = v19.pair_qualified_attempt_record(
            attempt=2, prior_attempt_failure_sha256="b" * 64,
            prior_attempt_acquisition_sha256=
            contract.canonical_sha256(acquisition),
            prior_attempt_selection_status="selected-lowest-primary",
            prior_attempt_selected_pair=selected)
        self.assertFalse(frozen["fresh_selection_permitted"])
        self.assertEqual(frozen["frozen_pair_from_prior_attempt"], selected)
        frozen_acquisition, unused_policy, unused_reader = acquisition_fixture(
            frozen_pair=selected)
        self.assertEqual(frozen_acquisition["scan"]["selection_status"],
                         "selected-frozen-pair")
        no_fallback, unused_policy, unused_reader = acquisition_fixture(
            frozen_pair=selected, nonidle_cpus=(2,))
        self.assertEqual(no_fallback["scan"]["eligible_pair_count"], 1)
        self.assertEqual(no_fallback["scan"]["eligible_pairs"], [{
            "benchmark_cpu": 3, "reserved_sibling": 7,
        }])
        self.assertEqual(no_fallback["scan"]["selection_status"],
                         "frozen-pair-did-not-requalify")
        self.assertIsNone(no_fallback["scan"]["selected"])

        no_candidate, unused_policy, unused_reader = acquisition_fixture(
            selected=False)
        fresh = v19.pair_qualified_attempt_record(
            attempt=2, prior_attempt_failure_sha256="c" * 64,
            prior_attempt_acquisition_sha256=
            contract.canonical_sha256(no_candidate),
            prior_attempt_selection_status="no-candidate-pair-qualified")
        self.assertTrue(fresh["fresh_selection_permitted"])
        self.assertIsNone(fresh["frozen_pair_from_prior_attempt"])
        no_acquisition = v19.pair_qualified_attempt_record(
            attempt=2, prior_attempt_failure_sha256="d" * 64,
            prior_attempt_selection_status="not-acquired")
        self.assertTrue(no_acquisition["fresh_selection_permitted"])

        for keywords in (
                {"attempt": 3},
                {"attempt": 1, "prior_attempt_failure_sha256": "b" * 64},
                {
                    "attempt": 2, "prior_attempt_failure_sha256": "b" * 64,
                    "prior_attempt_selection_status": "selected-lowest-primary",
                    "prior_attempt_acquisition_sha256": "c" * 64,
                },
                {
                    "attempt": 2, "prior_attempt_failure_sha256": "b" * 64,
                    "prior_attempt_selection_status":
                    "no-candidate-pair-qualified",
                    "prior_attempt_acquisition_sha256": "c" * 64,
                    "prior_attempt_selected_pair": selected,
                }):
            with self.subTest(keywords=keywords):
                self.assertV19Failure(
                    v19.pair_qualified_attempt_record, **keywords)

    def test_attempt_two_records_and_dormant_states_are_fixed(self) -> None:
        first_acquisition, unused_first_bridge, unused_first_policy = \
            accepted_bundle()
        selected = first_acquisition["scan"]["selected"]
        attempt_two = v19.pair_qualified_attempt_record(
            attempt=2, prior_attempt_failure_sha256="b" * 64,
            prior_attempt_acquisition_sha256=
            contract.canonical_sha256(first_acquisition),
            prior_attempt_selection_status="selected-lowest-primary",
            prior_attempt_selected_pair=selected)
        acquisition, bridge_record, policy = accepted_bundle(
            frozen_pair=selected)
        self.assertEqual(acquisition["scan"]["selection_status"],
                         "selected-frozen-pair")
        common = {
            "policy_value": policy,
            "expected_policy_sha256": contract.canonical_sha256(policy),
            "host_identity_sha256": "a" * 64,
            "attempt_value": attempt_two,
            "acquisition_value": acquisition,
            "bridge_value": bridge_record,
        }
        complete = v19.pair_qualified_v19_record(
            stage="complete", record_status="complete",
            terminal=v19.SUCCESS_TERMINAL,
            first_window_before_value=first_window_before(bridge_record),
            **common)
        self.assertEqual(self.validate(complete, policy, attempt_two), complete)
        self.assertEqual(complete["selection_status"], "selected-frozen-pair")
        self.assertEqual(
            contract.canonical_sha256(complete),
            "9325366c9f197f03a48a8d96a6700605e000098cdac6aa563c3f5a1b83d634ed")
        self.assertV19Failure(self.validate, complete, policy, attempt_one())

        first_complete, unused_acquisition, unused_bridge, first_policy, \
            unused_before = complete_record()
        self.assertV19Failure(
            self.validate, first_complete, first_policy, attempt_two)

        in_progress = v19.pair_qualified_v19_record(
            stage="handoff", record_status="in-progress", terminal=None,
            first_window_before_value=first_window_before(bridge_record),
            **common)
        self.assertEqual(
            self.validate(in_progress, policy, attempt_two), in_progress)

        unavailable = v19.pair_qualified_v19_record(
            stage="bridged", record_status="failed",
            terminal="first-window-handoff-unavailable", **common)
        self.assertEqual(
            self.validate(unavailable, policy, attempt_two), unavailable)

        no_fallback, unused_policy, unused_reader = acquisition_fixture(
            frozen_pair=selected, nonidle_cpus=(2,))
        frozen_failure = v19.pair_qualified_v19_record(
            stage="acquired", record_status="failed",
            terminal="frozen-pair-did-not-requalify", policy_value=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            host_identity_sha256="a" * 64, attempt_value=attempt_two,
            acquisition_value=no_fallback)
        self.assertEqual(
            self.validate(frozen_failure, policy, attempt_two), frozen_failure)
        self.assertV19Failure(
            v19.pair_qualified_v19_record,
            stage="acquired", record_status="failed",
            terminal="no-candidate-pair-qualified", policy_value=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            host_identity_sha256="a" * 64, attempt_value=attempt_two,
            acquisition_value=no_fallback)

    def test_handoff_boundary_and_pair_nonidle_are_fixed(self) -> None:
        unused_acquisition, bridge_record, unused_policy = accepted_bundle()
        accepted = v19.first_window_handoff_record(
            bridge_record, bridge_record["selected_pair"],
            first_window_before(
                bridge_record, elapsed_ns=v19.MAXIMUM_HANDOFF_ELAPSED_NS))
        self.assertTrue(accepted["accepted"])
        self.assertIsNone(accepted["failure_terminal"])

        late = v19.first_window_handoff_record(
            bridge_record, bridge_record["selected_pair"],
            first_window_before(
                bridge_record,
                elapsed_ns=v19.MAXIMUM_HANDOFF_ELAPSED_NS + 1))
        self.assertFalse(late["accepted"])
        self.assertEqual(late["failure_terminal"],
                         "first-window-handoff-late")

        for role in ("benchmark_cpu", "reserved_sibling"):
            nonidle = v19.first_window_handoff_record(
                bridge_record, bridge_record["selected_pair"],
                first_window_before(bridge_record, nonidle_role=role))
            self.assertFalse(nonidle["accepted"])
            self.assertEqual(
                nonidle["selected_pair_nonidle_delta"][role], 1)
            self.assertEqual(
                nonidle["failure_terminal"],
                "first-window-handoff-selected-pair-nonidle")

    def test_valid_but_unaccepted_bridge_cannot_reach_a_bridged_stage(self) -> None:
        acquisition, unused_bridge, policy = accepted_bundle()
        tail = acquisition["scan"]["windows"][-1]["after"]

        def step(previous: dict, *, nonidle_cpu: int | None = None) -> dict:
            counters = {}
            for key, record in previous["cpus"].items():
                fields = dict(record["fields"])
                fields["idle"] += 100
                if int(key) == nonidle_cpu:
                    fields["user"] += 1
                counters[int(key)] = fields
            started = previous["read_finished_monotonic_ns"] + 1_000_000_000
            return contract.shared_snapshot_record(
                read_started_monotonic_ns=started,
                read_finished_monotonic_ns=started + 1_000_000,
                counters=counters)

        first = step(tail)
        second = step(first, nonidle_cpu=2)
        rejected = live.bridge.pair_qualification_bridge_record(
            acquisition, expected_policy=policy, expected_frozen_pair=None,
            expected_acquisition_window_count=
            v19.QUALIFICATION_WINDOW_COUNT,
            expected_acquisition_nominal_window_ns=
            v19.QUALIFICATION_NOMINAL_WINDOW_NS,
            expected_bridge_geometry=v19.bridge_geometry_record(),
            snapshots=[tail, first, second])
        self.assertFalse(rejected["bridge_accepted"])
        common = {
            "policy_value": policy,
            "expected_policy_sha256": contract.canonical_sha256(policy),
            "host_identity_sha256": "a" * 64,
            "attempt_value": attempt_one(),
            "acquisition_value": acquisition,
            "bridge_value": rejected,
        }
        self.assertV19Failure(
            v19.pair_qualified_v19_record, stage="bridged",
            record_status="in-progress", terminal=None, **common)
        self.assertV19Failure(
            v19.pair_qualified_v19_record, stage="bridged",
            record_status="failed",
            terminal="first-window-handoff-unavailable", **common)

    def test_failure_stage_prefixes_are_closed_and_canonical(self) -> None:
        acquisition, bridge_record, policy = accepted_bundle()
        common = {
            "policy_value": policy,
            "expected_policy_sha256": contract.canonical_sha256(policy),
            "host_identity_sha256": "a" * 64,
            "attempt_value": attempt_one(),
        }
        identity = v19.pair_qualified_v19_record(
            stage="identity", record_status="failed",
            terminal="identity-failed", **common)
        self.assertEqual(self.validate(identity, policy), identity)

        no_candidate, unused_policy, unused_reader = acquisition_fixture(
            selected=False)
        acquired = v19.pair_qualified_v19_record(
            stage="acquired", record_status="failed",
            terminal="no-candidate-pair-qualified",
            acquisition_value=no_candidate, **common)
        self.assertEqual(self.validate(acquired, policy), acquired)

        bridge_failed = v19.pair_qualified_v19_record(
            stage="acquired", record_status="failed",
            terminal="bridge-not-accepted", acquisition_value=acquisition,
            **common)
        self.assertEqual(self.validate(bridge_failed, policy), bridge_failed)

        late_before = first_window_before(
            bridge_record, elapsed_ns=v19.MAXIMUM_HANDOFF_ELAPSED_NS + 1)
        handoff_failed = v19.pair_qualified_v19_record(
            stage="handoff", record_status="failed",
            terminal="first-window-handoff-late",
            acquisition_value=acquisition, bridge_value=bridge_record,
            first_window_before_value=late_before, **common)
        self.assertEqual(self.validate(handoff_failed, policy), handoff_failed)

        campaign_failed = v19.pair_qualified_v19_record(
            stage="campaign", record_status="failed",
            terminal="campaign-rejected", acquisition_value=acquisition,
            bridge_value=bridge_record,
            first_window_before_value=first_window_before(bridge_record),
            **common)
        self.assertEqual(self.validate(campaign_failed, policy), campaign_failed)

        invalid = copy.deepcopy(campaign_failed)
        invalid["stage"] = "bridged"
        self.assertV19Failure(self.validate, invalid, policy)

    def test_every_bound_hash_claim_and_endpoint_is_recomputed(self) -> None:
        record, unused_acquisition, unused_bridge, policy, unused_before = \
            complete_record()
        mutations = []
        for key, value in (
                ("policy_sha256", "0" * 64),
                ("host_identity_sha256", "0" * 64),
                ("qualification_geometry_sha256", "0" * 64),
                ("acquisition_sha256", "0" * 64),
                ("bridge_sha256", "0" * 64),
                ("bridge_geometry_sha256", "0" * 64),
                ("v18_failure_lineage_sha256", "0" * 64),
                ("terminal", "PROMOTED"),
                ("host_mutation_performed", True),
                ("candidate_timing_performed", True)):
            mutation = copy.deepcopy(record)
            mutation[key] = value
            mutations.append(mutation)
        geometry = copy.deepcopy(record)
        geometry["qualification_geometry"]["window_count"] = 11
        geometry["qualification_geometry_sha256"] = contract.canonical_sha256(
            geometry["qualification_geometry"])
        mutations.append(geometry)
        claim = copy.deepcopy(record)
        claim["shared_host_claim_ceiling"]["promotion_eligible"] = True
        mutations.append(claim)
        lineage = copy.deepcopy(record)
        lineage["v18_failure_lineage"]["attempts"][2][
            "envelope_sha256sums_sha256"] = "0" * 64
        lineage["v18_failure_lineage_sha256"] = contract.canonical_sha256(
            lineage["v18_failure_lineage"])
        mutations.append(lineage)
        splice = copy.deepcopy(record)
        splice["first_window_handoff"]["first_window_before"][
            "benchmark_cpu"]["fields"]["idle"] += 1
        splice["first_window_handoff"]["first_window_before"][
            "benchmark_cpu"]["total_jiffies"] += 1
        mutations.append(splice)
        bridge_splice = copy.deepcopy(record)
        bridge_splice["bridge"]["bridge_accepted"] = False
        mutations.append(bridge_splice)
        for mutation in mutations:
            with self.subTest(keys={
                    key for key in record if
                    not contract.exact_json_equal(record[key], mutation[key])}):
                self.assertV19Failure(self.validate, mutation, policy)

    def test_strict_loader_and_external_expectations_fail_closed(self) -> None:
        record, unused_acquisition, unused_bridge, policy, unused_before = \
            complete_record()
        data = contract.canonical_json_bytes(record)
        loaded = v19.load_pair_qualified_v19_record(
            data, expected_policy=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            expected_host_identity_sha256="a" * 64,
            expected_attempt=attempt_one())
        self.assertEqual(loaded, record)
        for invalid in (b" " + data, data + data, b'{"x":NaN}\n'):
            with self.subTest(invalid=invalid[:24]):
                self.assertV19Failure(
                    v19.load_pair_qualified_v19_record, invalid,
                    expected_policy=policy,
                    expected_policy_sha256=contract.canonical_sha256(policy),
                    expected_host_identity_sha256="a" * 64,
                    expected_attempt=attempt_one())
        self.assertV19Failure(
            v19.validate_pair_qualified_v19_record, record,
            expected_policy=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            expected_host_identity_sha256="b" * 64,
            expected_attempt=attempt_one())

    def test_module_has_no_host_process_runner_or_writer_surface(self) -> None:
        tree = ast.parse(MODULE_PATH.read_text("utf-8"))
        imports = {
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        }
        from_imports = {
            node.module for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertEqual(imports, {"copy", "importlib.util", "sys"})
        self.assertEqual(from_imports, {"__future__", "pathlib", "typing"})
        text = MODULE_PATH.read_text("utf-8")
        for forbidden in (
                "run_abba", "SystemHostReader", "/proc/", "/sys/",
                "subprocess", "sched_setaffinity", "time.monotonic",
                "os.open", "os.write", "fork(", "execve("):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
