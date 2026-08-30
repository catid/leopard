#!/usr/bin/env python3
"""Adversarial tests for the independent pair qualification verifier."""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
VERIFIER_PATH = HERE / "pair_qualification_verify.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "leopard2_pair_qualification_verify_fixtures",
    HERE / "test_pair_qualification_bridge_contract.py")
verifier = load_module(
    "leopard2_pair_qualification_independent_verifier_test", VERIFIER_PATH)
contract = verifier.contract


def bundle_fixture(*, nonidle_cpu: int | None = None):
    record, acquisition, policy, geometry = fixtures.bridge_fixture(
        nonidle_cpu=nonidle_cpu)
    keywords = {
        "expected_policy": policy,
        "expected_policy_sha256": contract.canonical_sha256(policy),
        "expected_frozen_pair": None,
        "expected_acquisition_window_count": 2,
        "expected_acquisition_nominal_window_ns": 250_000_000,
        "expected_bridge_geometry": geometry,
    }
    return (
        contract.canonical_json_bytes(acquisition),
        contract.canonical_json_bytes(record), policy, geometry, keywords,
    )


def python_command(*arguments: str) -> list[str]:
    command = [sys.executable, "-I", "-S", "-B"]
    if sys.flags.optimize:
        command.append("-O")
    command.extend(arguments)
    return command


class PairQualificationVerifierTest(unittest.TestCase):
    def assertVerificationFailure(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(contract.QualificationError):
            function(*args, **kwargs)

    def test_independent_verdict_has_a_frozen_canonical_projection(self) -> None:
        acquisition, bridge_data, unused_policy, unused_geometry, keywords = \
            bundle_fixture()
        verdict = verifier.require_accepted_pair_qualification_bundle(
            acquisition, bridge_data, **keywords)
        self.assertEqual(verdict["selected_pair"], {
            "benchmark_cpu": 2, "reserved_sibling": 6,
        })
        self.assertTrue(verdict["bridge_accepted"])
        self.assertTrue(verdict["fixed_point_verified"])
        self.assertFalse(verdict["verifier_host_access_performed"])
        self.assertFalse(verdict["candidate_timing_performed"])
        self.assertEqual(
            contract.canonical_sha256(verdict),
            "2a662b083dd2e464c7cb56bd2b49e5375462d216d2f46cca66db7d07831886bb")

    def test_valid_environment_rejection_is_verifiable_but_not_accepted(self) -> None:
        acquisition, bridge_data, unused_policy, unused_geometry, keywords = \
            bundle_fixture(nonidle_cpu=2)
        verdict = verifier.verify_pair_qualification_bundle(
            acquisition, bridge_data, **keywords)
        self.assertFalse(verdict["bridge_accepted"])
        self.assertVerificationFailure(
            verifier.require_accepted_pair_qualification_bundle,
            acquisition, bridge_data, **keywords)

    def test_producer_attestations_cannot_override_independent_replay(self) -> None:
        acquisition, bridge_data, unused_policy, unused_geometry, keywords = \
            bundle_fixture(nonidle_cpu=2)
        record = contract.strict_json_loads(bridge_data, "test bridge")
        record["bridge_accepted"] = True
        record["nonidle_guarded_cpus"] = []
        self.assertVerificationFailure(
            verifier.verify_pair_qualification_bundle,
            acquisition, contract.canonical_json_bytes(record), **keywords)

    def test_splices_relabels_and_malformed_acquisition_fail_closed(self) -> None:
        acquisition_data, bridge_data, policy, geometry, keywords = bundle_fixture()
        acquisition = contract.strict_json_loads(acquisition_data, "test acquisition")
        record = contract.strict_json_loads(bridge_data, "test bridge")
        mutations = []
        pair = copy.deepcopy(record)
        pair["selected_pair"] = {"benchmark_cpu": 3, "reserved_sibling": 7}
        mutations.append((acquisition_data, contract.canonical_json_bytes(pair), keywords))
        splice = copy.deepcopy(record)
        splice["windows"][0]["before"]["cpus"]["2"]["fields"]["idle"] += 1
        splice["windows"][0]["before"]["cpus"]["2"]["total_jiffies"] += 1
        mutations.append((acquisition_data, contract.canonical_json_bytes(splice), keywords))
        malformed = copy.deepcopy(acquisition)
        malformed["scan_sha256"] = "0" * 64
        mutations.append((contract.canonical_json_bytes(malformed), bridge_data, keywords))
        wrong_hash = dict(keywords)
        wrong_hash["expected_policy_sha256"] = "0" * 64
        mutations.append((acquisition_data, bridge_data, wrong_hash))
        changed_policy = contract.qualification_policy_record(
            clock_ticks_per_second=100, candidate_primary_cpus=[2, 3],
            excluded_pairs=[], domain_mode="pair-only")
        relabelled = dict(keywords)
        relabelled["expected_policy"] = changed_policy
        relabelled["expected_policy_sha256"] = contract.canonical_sha256(
            changed_policy)
        mutations.append((acquisition_data, bridge_data, relabelled))
        changed_geometry = dict(keywords)
        changed_geometry["expected_bridge_geometry"] = \
            fixtures.geometry_fixture(maximum=5)
        mutations.append((acquisition_data, bridge_data, changed_geometry))
        for mutated_acquisition, mutated_bridge, mutated_keywords in mutations:
            with self.subTest(
                    acquisition=hashlib.sha256(mutated_acquisition).hexdigest(),
                    bridge=hashlib.sha256(mutated_bridge).hexdigest()):
                self.assertVerificationFailure(
                    verifier.verify_pair_qualification_bundle,
                    mutated_acquisition, mutated_bridge, **mutated_keywords)

    def test_noncanonical_and_multivalue_wire_bytes_are_rejected(self) -> None:
        acquisition, bridge_data, unused_policy, unused_geometry, keywords = \
            bundle_fixture()
        for changed_acquisition, changed_bridge in (
                (b" " + acquisition, bridge_data),
                (acquisition, b" " + bridge_data),
                (acquisition + acquisition, bridge_data),
                (acquisition, bridge_data + bridge_data)):
            with self.subTest(
                    acquisition_size=len(changed_acquisition),
                    bridge_size=len(changed_bridge)):
                self.assertVerificationFailure(
                    verifier.verify_pair_qualification_bundle,
                    changed_acquisition, changed_bridge, **keywords)

    def test_verify_cli_reads_only_named_canonical_evidence(self) -> None:
        acquisition, bridge_data, policy, geometry, keywords = bundle_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            policy_path = root / "policy.json"
            acquisition_path = root / "acquisition.json"
            bridge_path = root / "bridge.json"
            policy_path.write_bytes(contract.canonical_json_bytes(policy))
            acquisition_path.write_bytes(acquisition)
            bridge_path.write_bytes(bridge_data)
            arguments = [
                str(VERIFIER_PATH), "verify",
                "--policy", str(policy_path),
                "--acquisition", str(acquisition_path),
                "--bridge", str(bridge_path),
                "--expected-policy-sha256", keywords["expected_policy_sha256"],
                "--expected-acquisition-window-count", "2",
                "--expected-acquisition-nominal-window-ns", "250000000",
                "--minimum-bridge-window-count",
                str(geometry["minimum_window_count"]),
                "--maximum-bridge-window-count",
                str(geometry["maximum_window_count"]),
                "--bridge-nominal-window-ns",
                str(geometry["nominal_window_ns"]),
                "--maximum-bridge-handoff-elapsed-ns",
                str(geometry["maximum_handoff_elapsed_ns"]),
            ]
            completed = subprocess.run(
                python_command(*arguments), check=False,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(completed.returncode, 0, completed.stderr.decode())
            verdict = contract.strict_json_loads(completed.stdout, "CLI verdict")
            self.assertEqual(
                contract.canonical_sha256(verdict),
                "2a662b083dd2e464c7cb56bd2b49e5375462d216d2f46cca66db7d07831886bb")

            one_role = subprocess.run(
                python_command(*arguments, "--frozen-benchmark-cpu", "2"),
                check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(one_role.returncode, 0)
            self.assertIn(b"requires both CPU roles", one_role.stderr)

    def test_cli_self_test_output_is_frozen(self) -> None:
        completed = subprocess.run(
            python_command(str(VERIFIER_PATH), "self-test"), check=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())
        self.assertEqual(
            completed.stdout,
            (verifier.SELF_TEST_OUTPUT + "\n").encode("utf-8"))

    def test_verifier_has_no_host_or_process_access_surface(self) -> None:
        tree = ast.parse(VERIFIER_PATH.read_text("utf-8"))
        imports = {
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        }
        from_imports = {
            node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
        }
        self.assertEqual(imports, {
            "argparse", "copy", "importlib.util", "sys",
        })
        self.assertEqual(from_imports, {"__future__", "pathlib", "typing"})
        text = VERIFIER_PATH.read_text("utf-8")
        for forbidden in (
                "/proc/", "/sys/", "sched_getaffinity", "sched_setaffinity",
                "subprocess", "time.monotonic", "sleep(", "fork(", "execve("):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
