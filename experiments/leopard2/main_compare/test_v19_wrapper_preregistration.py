#!/usr/bin/python3

"""Deterministic tests for the dormant conditioned-v19 wrapper contract."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WRAPPER = ROOT / "experiments/leopard2/main_compare" / \
    "run_authoritative_v17_gfni_main_compare.sh"


class V19WrapperPreregistrationTests(unittest.TestCase):
    def test_preregistration_is_exact_and_dormant(self) -> None:
        self_test = subprocess.run(
            [str(WRAPPER), "--self-test-conditioned-v19-contract"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            text=True,
        )
        self.assertEqual(self_test.returncode, 0, self_test.stderr)
        self.assertEqual(
            self_test.stdout,
            "v19 conditioned wrapper preregistration self-test passed\n")
        self.assertEqual(self_test.stderr, "")

        printed = subprocess.run(
            [str(WRAPPER), "--print-conditioned-v19-preregistration"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        self.assertEqual(printed.returncode, 0, printed.stderr.decode())
        self.assertEqual(printed.stderr, b"")
        self.assertEqual(
            hashlib.sha256(printed.stdout).hexdigest(),
            "27c1b7d76a0ecdbe194d6e6b62c01e48b1c7d10fc8ef99ebad4d76238669f0c1")
        preregistration = json.loads(printed.stdout)
        self.assertEqual(
            preregistration["schema"],
            "leopard2-v19-conditioned-main-wrapper-preregistration/v1")
        self.assertIs(preregistration["live_acquisition_armed"], False)
        self.assertIs(preregistration["timing_performed"], False)
        self.assertEqual(
            preregistration["resource_envelope"]["memory_max_bytes"],
            536870912)
        self.assertEqual(
            preregistration["resource_envelope"]["memory_swap_max_bytes"],
            0)
        self.assertEqual(
            preregistration["resource_envelope"]
            ["release_authorized_max_jobs"], 2)
        self.assertEqual(
            preregistration["resource_envelope"]["release_max_jobs"], 1)
        self.assertEqual(
            preregistration["resource_envelope"]
            ["release_max_jobs_basis"],
            "retained-preflight-proven-cap")
        self.assertEqual(
            preregistration["build_preflight"]
            ["observed_memory_peak_bytes"], 391512064)
        self.assertEqual(
            preregistration["build_preflight"]["compiler"],
            "GNU c++ 13.3.0")
        self.assertEqual(
            preregistration["build_preflight"]["parallelism"], 1)
        self.assertEqual(
            preregistration["build_preflight"]["limit_nofile"], 65536)
        self.assertEqual(
            preregistration["build_preflight"]["execution_scope"],
            "systemd-user-service")
        self.assertIs(
            preregistration["build_preflight"]["zero_proc_swaps_rows"],
            True)
        self.assertIs(
            preregistration["build_preflight"]["builds_succeeded"], True)
        self.assertIs(
            preregistration["build_preflight"]
            ["provenance_validation_passed"], True)
        for name in (
                "benchmark_executable_executed",
                "qualification_scan_performed",
                "candidate_workload_executed",
                "v19_attempt_consumed",
                "foreign_process_mutation_performed",
                "timing_performed"):
            self.assertIs(preregistration["build_preflight"][name], False)
        self.assertEqual(
            preregistration["attempt_contract"]["attempt_budget"], 2)
        self.assertIs(preregistration["bridge"]["pair_lease_held"], False)
        self.assertIs(
            preregistration["bridge"]
            ["pair_lease_acquired_only_after_accepted_bridge"], True)
        for value in preregistration["claim_ceiling"].values():
            self.assertIs(value, False)

    def test_unarmed_v19_and_exhausted_v18_refuse_before_lane(self) -> None:
        wrapper_text = WRAPPER.read_text(encoding="utf-8")
        verify_dispatch = wrapper_text.index(
            'if [[ $# -eq 2 && $1 == --verify && $2 == /* ]]')
        v19_refusal = wrapper_text.index(
            "conditioned v19 authoritative acquisition is NOT ARMED")
        v18_refusal = wrapper_text.index(
            "fresh passive-v2 v18 acquisition is disabled")
        lane_creation = wrapper_text.index(
            '/usr/bin/mkdir -m 0700 "$envelope"')
        self.assertLess(verify_dispatch, v19_refusal)
        self.assertLess(v19_refusal, lane_creation)
        self.assertLess(v18_refusal, lane_creation)

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-dormant-") as temporary:
            root = Path(temporary)
            cases = (
                (("--conditioned-passive-v19", "--attempt", "1",
                  "--attempt-budget", "2"),
                 "conditioned v19 authoritative acquisition is NOT ARMED at "
                 "this preregistration commit; no scan, lane, or workload "
                 "was created"),
                (("--conditioned-passive-v19", "--attempt", "2",
                  "--attempt-budget", "2"),
                 "conditioned v19 authoritative acquisition is NOT ARMED at "
                 "this preregistration commit; no scan, lane, or workload "
                 "was created"),
                (("--conditioned-passive-v19", "--attempt", "3",
                  "--attempt-budget", "2"),
                 "conditioned v19 acquisition requires --attempt N (1..2) "
                 "and --attempt-budget 2"),
                (("--passive-shared-host", "--attempt", "1",
                  "--attempt-budget", "3"),
                 "fresh passive-v2 v18 acquisition is disabled because its "
                 "three-attempt budget is exhausted; retained v18 envelopes "
                 "remain verifiable"),
                (("--passive-shared-host", "--attempt", "2",
                  "--attempt-budget", "3"),
                 "fresh passive-v2 v18 acquisition is disabled because its "
                 "three-attempt budget is exhausted; retained v18 envelopes "
                 "remain verifiable"),
                (("--passive-shared-host", "--attempt", "3",
                  "--attempt-budget", "3"),
                 "fresh passive-v2 v18 acquisition is disabled because its "
                 "three-attempt budget is exhausted; retained v18 envelopes "
                 "remain verifiable"),
            )
            for index, (arguments, expected_stderr) in enumerate(cases):
                lane = root / f"lane-{index}"
                completed = subprocess.run(
                    [str(WRAPPER), *arguments, str(lane)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=False, text=True,
                )
                with self.subTest(arguments=arguments):
                    self.assertEqual(completed.returncode, 2)
                    self.assertEqual(completed.stdout, "")
                    self.assertEqual(
                        completed.stderr, expected_stderr + "\n")
                    self.assertFalse(lane.exists())
            bypass_lane = root / "mode-prefixed-verify"
            bypass = subprocess.run(
                [str(WRAPPER), "--conditioned-passive-v19", "--attempt", "1",
                 "--attempt-budget", "2", "--verify", str(bypass_lane)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, text=True,
            )
            self.assertEqual(bypass.returncode, 2)
            self.assertEqual(bypass.stdout, "")
            self.assertEqual(
                bypass.stderr,
                "--conditioned-passive-v19 cannot be combined with --verify; "
                "generation comes from the sealed envelope\n")
            self.assertFalse(bypass_lane.exists())
            passive_bypass = subprocess.run(
                [str(WRAPPER), "--passive-shared-host", "--attempt", "1",
                 "--attempt-budget", "3", "--verify", str(bypass_lane)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, text=True,
            )
            self.assertEqual(passive_bypass.returncode, 2)
            self.assertEqual(passive_bypass.stdout, "")
            self.assertEqual(
                passive_bypass.stderr,
                "--passive-shared-host cannot be combined with --verify; "
                "generation comes from the sealed envelope\n")
            self.assertFalse(bypass_lane.exists())
            bare_verify = subprocess.run(
                [str(WRAPPER), "--verify", str(root / "missing-envelope")],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, text=True,
            )
            self.assertEqual(bare_verify.returncode, 1)
            self.assertEqual(bare_verify.stdout, "")
            self.assertEqual(bare_verify.stderr, "")

            hidden_verify = (
                "--verify-v18-complete-core-semantics",
                str(root / "missing-envelope"),
                "/tmp/leopard-v18-complete-core-replay.fixture/"
                "NOT_PROMOTED.json",
                "/tmp/leopard-v18-complete-core-replay.fixture/"
                "postseal-audit.json",
            )
            mode_arguments = (
                ("--conditioned-passive-v19", "--attempt", "1",
                 "--attempt-budget", "2"),
                ("--passive-shared-host", "--attempt", "1",
                 "--attempt-budget", "3"),
            )
            guarded_dispatches = (
                hidden_verify,
                ("--print-conditioned-v19-preregistration",),
                ("--self-test-conditioned-v19-contract",),
            )
            for mode in mode_arguments:
                for dispatch in guarded_dispatches:
                    completed = subprocess.run(
                        [str(WRAPPER), *mode, *dispatch],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        check=False, text=True,
                    )
                    with self.subTest(mode=mode[0], dispatch=dispatch[0]):
                        self.assertEqual(completed.returncode, 1)
                        self.assertEqual(completed.stdout, "")
                        self.assertEqual(completed.stderr, "")


if __name__ == "__main__":
    unittest.main()
