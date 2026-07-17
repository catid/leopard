#!/usr/bin/env python3
"""Focused fail-closed tests for the experiment-F checkpoint."""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "received_subset", HERE / "received_subset.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ReceivedSubsetTests(unittest.TestCase):
    def test_self_test_covers_dp_and_direct_mds(self) -> None:
        result = MODULE.self_test()
        self.assertEqual(result["field_checks"], 4367)
        self.assertEqual(result["dp_vs_bruteforce_patterns"], 546)
        self.assertEqual(result["representative_mds_subsets"], 140)

    def test_profile_coordinate_maps(self) -> None:
        high = MODULE.make_geometry("legacy_high_v1", 5, 3)
        self.assertEqual((high.side, high.parent, high.dimension), (4, 16, 12))
        self.assertEqual(high.original_coordinates, (4, 5, 6, 7, 8))
        self.assertEqual(high.recovery_coordinates, (0, 1, 2))
        self.assertEqual(high.source_coordinates, tuple(range(4, 16)))

        low = MODULE.make_geometry("low_v1", 5, 3)
        self.assertEqual((low.side, low.parent, low.dimension), (8, 16, 8))
        self.assertEqual(low.original_coordinates, (0, 1, 2, 3, 4))
        self.assertEqual(low.recovery_coordinates, (8, 9, 10))
        self.assertEqual(low.source_coordinates, tuple(range(8)))

    def test_lagrange_and_vandermonde_oracles_match_every_cell(self) -> None:
        entries = 0
        for profile, k, r in MODULE.valid_cells():
            geometry = MODULE.make_geometry(profile, k, r)
            self.assertEqual(MODULE.generator_rows(geometry),
                             MODULE.vandermonde_generator_rows(geometry))
            entries += (k + r) * k
        self.assertEqual(len(MODULE.valid_cells()), 170)
        self.assertEqual(entries, 9758)

    def test_exact_dp_matches_bruteforce_for_all_small_public_masks(self) -> None:
        comparisons = 0
        for profile, k, r in (
            ("legacy_high_v1", 2, 6),
            ("legacy_high_v1", 5, 3),
            ("low_v1", 2, 6),
            ("low_v1", 5, 3),
        ):
            geometry = MODULE.make_geometry(profile, k, r)
            public = geometry.public_coordinates
            for mask in range(1 << len(public)):
                available = tuple(
                    coordinate for index, coordinate in enumerate(public)
                    if mask & (1 << index))
                if len(available) < k:
                    continue
                self.assertEqual(
                    MODULE.select_exact_block_dp(geometry, available),
                    MODULE.brute_force_select(geometry, available))
                comparisons += 1
        self.assertEqual(comparisons, 680)

    def test_metric_ties_still_use_coordinate_objective(self) -> None:
        lower_coordinates = MODULE.selection_cost((0, 2), 4)
        higher_coordinates = MODULE.selection_cost((1, 2), 4)
        self.assertEqual(lower_coordinates.metric_key,
                         higher_coordinates.metric_key)
        self.assertLess(lower_coordinates.total_key,
                        higher_coordinates.total_key)

        geometry = MODULE.make_geometry("legacy_high_v1", 3, 3)
        available = (4, 5, 1, 2)
        greedy = MODULE.selection_cost(
            MODULE.select_block_greedy(geometry, available), geometry.side)
        exact = MODULE.selection_cost(
            MODULE.select_exact_block_dp(geometry, available), geometry.side)
        self.assertEqual(greedy.selected, (1, 4, 5))
        self.assertEqual(exact.selected, (1, 2, 4))
        self.assertEqual(greedy.metric_key, exact.metric_key)
        self.assertGreater(greedy.total_key, exact.total_key)

    def test_policies_are_deterministic_and_can_diverge(self) -> None:
        geometry = MODULE.make_geometry("legacy_high_v1", 5, 3)
        available = geometry.public_coordinates
        lowest = MODULE.select_lowest_parent(geometry, available)
        systematic = MODULE.select_prefer_systematic(geometry, available)
        self.assertEqual(lowest, (0, 1, 2, 4, 5))
        self.assertEqual(systematic, (4, 5, 6, 7, 8))
        self.assertNotEqual(lowest, systematic)
        for policy in MODULE.POLICIES:
            self.assertEqual(MODULE.select(policy, geometry, available),
                             MODULE.select(policy, geometry, available))

    def test_every_policy_rejects_insufficient_receives(self) -> None:
        geometry = MODULE.make_geometry("low_v1", 5, 3)
        available = geometry.public_coordinates[:4]
        for policy in MODULE.POLICIES:
            with self.subTest(policy=policy), self.assertRaises(ValueError):
                MODULE.select(policy, geometry, available)

    def test_exact_cost_never_loses_exhaustive_masks(self) -> None:
        checked = 0
        geometry = MODULE.make_geometry("low_v1", 4, 8)
        public = geometry.public_coordinates
        for erased_count in range(geometry.r + 1):
            for erased in itertools.combinations(range(len(public)), erased_count):
                erased_set = set(erased)
                available = tuple(coordinate for index, coordinate in enumerate(public)
                                  if index not in erased_set)
                exact = MODULE.selection_cost(
                    MODULE.select_exact_block_dp(geometry, available), geometry.side)
                for policy in MODULE.POLICIES:
                    candidate = MODULE.selection_cost(
                        MODULE.select(policy, geometry, available), geometry.side)
                    self.assertLessEqual(exact.metric_key, candidate.metric_key)
                checked += 1
        self.assertEqual(checked, 3797)

    def test_json_write_is_stable_and_atomic(self) -> None:
        value = {"z": [3, 2, 1], "a": {"field": "GF(2^4)"}}
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "nested" / "result.json"
            MODULE.write_json(path, value)
            first = path.read_bytes()
            MODULE.write_json(path, json.loads(first))
            self.assertEqual(path.read_bytes(), first)
            self.assertFalse(path.with_suffix(".json.tmp").exists())

    def test_retained_checkpoint_validates_and_rejects_mutation(self) -> None:
        retained = HERE / "results/checkpoint.json"
        validated = MODULE.validate_checkpoint(retained)
        self.assertEqual(validated["profile_cells"], 170)
        document = json.loads(retained.read_text(encoding="utf-8"))
        patterns = document["totals"]["availability_patterns"]
        for policy in MODULE.POLICIES:
            policy_total = document["policy_totals"][policy]
            comparison = policy_total["comparisons"]
            self.assertEqual(sum(comparison[key] for key in (
                "strictly_better_than_prefer_systematic",
                "equal_to_prefer_systematic",
                "strictly_worse_than_prefer_systematic",
            )), patterns)
            self.assertEqual(sum(comparison[key] for key in (
                "strictly_worse_than_exact", "equal_to_exact",
            )), patterns)
            self.assertEqual(sum(comparison[key] for key in (
                "metric_strictly_better_than_prefer_systematic",
                "metric_equal_to_prefer_systematic",
                "metric_strictly_worse_than_prefer_systematic",
            )), patterns)
            self.assertEqual(sum(comparison[key] for key in (
                "metric_strictly_worse_than_exact", "metric_equal_to_exact",
            )), patterns)
            self.assertEqual(
                comparison[MODULE.TIEBREAK_COMPARISON_KEY],
                comparison["metric_equal_to_exact"] -
                comparison["equal_to_exact"])
            self.assertEqual(
                policy_total["cost_sums"]["selected_subset_changes"],
                patterns - comparison["equal_to_prefer_systematic"])

        exact = document["policy_totals"]["exact_block_dp"]["comparisons"]
        self.assertEqual(exact["strictly_better_than_prefer_systematic"],
                         333421)
        self.assertEqual(exact["equal_to_prefer_systematic"], 425505)
        self.assertEqual(exact["metric_strictly_better_than_prefer_systematic"],
                         309047)
        self.assertEqual(exact["metric_equal_to_prefer_systematic"], 449879)
        self.assertEqual(
            document["policy_totals"]["exact_block_dp"]["cost_sums"][
                "selected_subset_changes"], 333421)

        greedy = document["policy_totals"]["block_aligned_greedy"][
            "comparisons"]
        self.assertEqual(greedy["equal_to_exact"], 667518)
        self.assertEqual(greedy["metric_equal_to_exact"], 702601)
        self.assertEqual(greedy[MODULE.TIEBREAK_COMPARISON_KEY], 35083)

        document["totals"]["mds_subsets"] += 1
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "mutated.json"
            MODULE.write_json(path, document)
            with self.assertRaises(ValueError):
                MODULE.validate_checkpoint(path)

    def test_validator_rejects_internally_remerged_tiebreak_error(self) -> None:
        retained = HERE / "results/checkpoint.json"
        document = json.loads(retained.read_text(encoding="utf-8"))
        cell = next(
            row for row in document["cells"]
            if row["comparison"]["block_aligned_greedy"][
                MODULE.TIEBREAK_COMPARISON_KEY] > 0)
        cell["comparison"]["block_aligned_greedy"][
            MODULE.TIEBREAK_COMPARISON_KEY] += 1
        malformed = MODULE.merge(document["cells"], document["source_sha256"])
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "remerged.json"
            MODULE.write_json(path, malformed)
            with self.assertRaisesRegex(ValueError,
                                        "metric-tie/tiebreak partition"):
                MODULE.validate_checkpoint(path)


if __name__ == "__main__":
    unittest.main()
