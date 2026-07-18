#!/usr/bin/env python3
"""Adversarial replay tests for the retained GF16 neighbor evidence."""

from __future__ import annotations

import copy
import gzip
import importlib.util
import json
import sys
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNNER_PATH = HERE / "run_neighbor_abba.py"
RAW_PATH = (HERE / "results" / "7921998-neighbor-amd9950x3d" /
            "raw.json.gz")


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "gf16_neighbor_runner_under_test", RUNNER_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {RUNNER_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()


class NeighborEvidenceValidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with gzip.open(RAW_PATH, "rt", encoding="utf-8") as source:
            cls.raw = json.load(source)

    def test_every_retained_result_validates(self) -> None:
        count = 0
        for cell_record in self.raw["cells"]:
            cell = cell_record["cell"]
            for round_record in cell_record["rounds"]:
                for invocation in round_record["invocations"]:
                    normalized = RUNNER.validate_result(
                        invocation["result"], cell)
                    self.assertEqual(
                        normalized["digests"], invocation["normalized"]["digests"])
                    self.assertEqual(
                        normalized["resolved"], invocation["normalized"]["resolved"])
                    self.assertEqual(
                        normalized["encode_us"], invocation["normalized"]["encode_us"])
                    self.assertEqual(
                        normalized["decode_us"], invocation["normalized"]["decode_us"])
                    count += 1
        self.assertEqual(count, 156)

    def retained_exemplar(self):
        cell_record = self.raw["cells"][0]
        invocation = cell_record["rounds"][0]["invocations"][0]
        return cell_record["cell"], copy.deepcopy(invocation["result"])

    def assert_rejected(self, mutator) -> None:
        cell, result = self.retained_exemplar()
        mutator(result)
        with self.assertRaises(RUNNER.EvidenceError):
            RUNNER.validate_result(result, cell)

    def test_wrong_requested_profile_is_rejected(self) -> None:
        self.assert_rejected(
            lambda result: result["parameters"].__setitem__(
                "requested_profile", "low_v1"))

    def test_wrong_selected_decode_path_is_rejected(self) -> None:
        self.assert_rejected(
            lambda result: result["resolved"].__setitem__(
                "selected_decode_path", "materialized"))

    def test_wrong_selected_decode_rule_is_rejected(self) -> None:
        self.assert_rejected(
            lambda result: result["resolved"].__setitem__(
                "selected_decode_rule", "workspace_materialized"))


if __name__ == "__main__":
    unittest.main()
