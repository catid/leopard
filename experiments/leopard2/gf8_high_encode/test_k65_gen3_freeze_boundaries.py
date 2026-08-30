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

FROZEN_SHA256 = {
    RUNNER_PATH: "2a5489ba7c1866135e5fc1577c3f4290e851bb41837d0e46e1118ac7699397ca",
    CONTRACT_PATH: "70a0d477ff947aaa301775ecee3370be42bb7b2e4d84a1ff372a411ee9eb8900",
    CONTRACT_TEST_PATH: "f670a4bcbaa3d6df784e81e163fa618d65780a135d234d8e5c27360e953199f1",
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
        acquisition = MAIN_COMPARE / "pair_qualification_acquire.py"
        text = acquisition.read_text("utf-8")
        self.assertNotIn("run_k65r65_b64_packed_terminal_abba", text)
        self.assertNotIn("generation_projection", text)
        self.assertNotIn("ATTEMPT_BUDGET", text)


if __name__ == "__main__":
    unittest.main()
