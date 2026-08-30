#!/usr/bin/env python3
"""Freeze the K65 v2 evidence and pure qualification boundaries for v3 work."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
MAIN_COMPARE = HERE.parent / "main_compare"
RUNNER_PATH = HERE / "run_k65r65_b64_packed_terminal_abba.py"
CONTRACT_PATH = MAIN_COMPARE / "pair_qualification_contract.py"
CONTRACT_TEST_PATH = MAIN_COMPARE / "test_pair_qualification_contract.py"
ACQUISITION_PATH = MAIN_COMPARE / "pair_qualification_acquire.py"
BRIDGE_PATH = MAIN_COMPARE / "pair_qualification_bridge_contract.py"
BRIDGE_TEST_PATH = MAIN_COMPARE / "test_pair_qualification_bridge_contract.py"
VERIFIER_PATH = MAIN_COMPARE / "pair_qualification_verify.py"
VERIFIER_TEST_PATH = MAIN_COMPARE / "test_pair_qualification_verify.py"
PREREGISTRATION_PATH = HERE / "k65_gen3_preregistration.py"
PLAN_RUNNER_PATH = HERE / "run_k65r65_b64_packed_terminal_gen3_abba.py"
TEMPLATE_PATH = HERE / \
    "k65r65_b64_packed_terminal_gen3_preregistration.template.json"

FROZEN_SHA256 = {
    RUNNER_PATH: "2a5489ba7c1866135e5fc1577c3f4290e851bb41837d0e46e1118ac7699397ca",
    CONTRACT_PATH: "70a0d477ff947aaa301775ecee3370be42bb7b2e4d84a1ff372a411ee9eb8900",
    CONTRACT_TEST_PATH: "f670a4bcbaa3d6df784e81e163fa618d65780a135d234d8e5c27360e953199f1",
    ACQUISITION_PATH: "d26d926a9afabcedada5d96c3beb9d8b02c0751a169e00635a6080557f053f00",
    BRIDGE_PATH: "bbe5d306a33e58c67bf91b7a9e995dccce46e29969955f3447a8de5938dcc110",
    BRIDGE_TEST_PATH: "27bfaf5565ef51d711f1faac03a5e79e1c5d4263fe4c5dee27e8ed691d9bd59c",
    VERIFIER_PATH: "eb87be63bbc2865ae7d95c8de88e20c8606600d99be5efbd4d737461f3805d9f",
    VERIFIER_TEST_PATH: "780e837ce2d2f6ae631535697ce7a5fe7ec7fc2d436961b03eae4eb700b7f380",
}
GENERATION2_PROJECTION_SHA256 = \
    "50a56a31ff68fb1143347d1bef51701e16371a458e4b20a15e2c31999aae1ee6"


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "leopard2_k65_generation2_frozen_runner", RUNNER_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load frozen K65 generation-2 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


class K65Generation3FreezeBoundaryTest(unittest.TestCase):
    def test_v2_runner_and_pure_contract_bytes_are_frozen(self) -> None:
        for path, expected in FROZEN_SHA256.items():
            with self.subTest(path=str(path)):
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected)

    def test_generation2_policy_and_attempt_budget_are_frozen(self) -> None:
        runner = load_runner()
        self.assertEqual(runner.GENERATION, 2)
        self.assertEqual(runner.ATTEMPT_BUDGET, 3)
        self.assertEqual(runner.CONFIRMATORY_ROUNDS, 25)
        self.assertEqual(runner.FIXED_ITERATIONS, 31)
        self.assertEqual(runner.FIXED_WARMUP, 64)
        self.assertEqual(runner.FIXED_REUSE, 8192)
        self.assertEqual(runner.BASE.MAX_ISOLATION_ATTEMPTS, 3)
        self.assertEqual(runner.BASE.TARGET_CONTROL_FLOOR, 1.05)
        self.assertEqual(runner.BASE.TARGET_MAIN_FLOOR, 1.05)
        self.assertEqual(runner.BASE.NEIGHBOR_EQUIVALENCE_LOWER, 1.0 / 1.02)
        self.assertEqual(runner.BASE.NEIGHBOR_EQUIVALENCE_UPPER, 1.02)
        self.assertEqual(runner.BASE.RETAINED_MAIN_FLOOR, 0.98)
        self.assertEqual(runner.MATRIX_SHA256,
                         "ab9572c4101b2af5eda4b7cfab17e979239f698c4c1196e660cf6f5e3f4af27c")
        self.assertEqual(
            runner._canonical_sha256(runner.generation_projection(2)),
            GENERATION2_PROJECTION_SHA256)
        with self.assertRaises(runner.BASE.EvidenceError):
            runner.generation_projection(3)

    def test_new_acquisition_is_additive_and_does_not_import_the_v2_runner(self) -> None:
        additive_paths = (
            ACQUISITION_PATH,
            MAIN_COMPARE / "pair_qualification_bridge_contract.py",
            MAIN_COMPARE / "pair_qualification_verify.py",
        )
        for path in additive_paths:
            with self.subTest(path=str(path)):
                text = path.read_text("utf-8")
                self.assertNotIn("run_k65r65_b64_packed_terminal_abba", text)
                self.assertNotIn("generation_projection", text)
                self.assertNotIn("ATTEMPT_BUDGET", text)

    def test_gen3_template_and_plan_runner_are_non_executable(self) -> None:
        prereg = load_path_module(
            "leopard2_k65_gen3_preregistration_freeze", PREREGISTRATION_PATH)
        plan = load_path_module(
            "leopard2_k65_gen3_plan_runner_freeze", PLAN_RUNNER_PATH)
        template = prereg.load_preregistration_template(
            TEMPLATE_PATH.read_bytes(), verify_files=False)
        self.assertFalse(template["safe_to_execute"])
        self.assertEqual(template["ratification_status"],
                         "requires-explicit-user-or-external-authority")
        self.assertFalse(hasattr(plan, "BASE"))
        source = PLAN_RUNNER_PATH.read_text("utf-8")
        self.assertNotIn("subprocess", source)
        self.assertNotIn("sched_setaffinity", source)
        self.assertIn("--plan-only", source)


def load_path_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


if __name__ == "__main__":
    unittest.main()
