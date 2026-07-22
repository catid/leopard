#!/usr/bin/env python3
"""Contract tests for run_small_direct_exhaustive.py."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


RUNNER_PATH = Path(__file__).with_name("run_small_direct_exhaustive.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_test_small_direct_exhaustive_runner", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import exhaustive direct-repair runner")
RUNNER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def shard_record(index: int, count: int, mode: int = 1) -> dict:
    assigned = RUNNER.EXPECTED_PATTERNS // count + (
        1 if index < RUNNER.EXPECTED_PATTERNS % count else 0)
    return {
        "schema": RUNNER.SHARD_SCHEMA,
        "mode": mode,
        "shard_index": index,
        "shard_count": count,
        "total_patterns": RUNNER.EXPECTED_PATTERNS,
        "assigned_patterns": assigned,
        "recovered_shards": assigned * 5,
        "verified_basis_symbols": assigned * 25,
        "basis_seed": 0,
        "assignment": "global_ordinal_mod_shard_count",
        "digest_fnv1a64": "%016x" % (index + 1),
        "ordinal_digest_fnv1a64": "%016x" % (index + 101),
    }


class ExhaustiveRunnerTests(unittest.TestCase):
    def test_twenty_eight_shards_cover_exact_matrix_once(self) -> None:
        records = [
            RUNNER.validate_shard(shard_record(index, 28), index, 28)
            for index in range(28)
        ]
        self.assertEqual(
            sum(value["assigned_patterns"] for value in records),
            1_982_812)
        self.assertLessEqual(
            max(value["assigned_patterns"] for value in records) -
            min(value["assigned_patterns"] for value in records), 1)

    def test_wrong_count_or_mode_is_rejected(self) -> None:
        wrong_count = shard_record(0, 28)
        wrong_count["assigned_patterns"] -= 1
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "identity/counts"):
            RUNNER.validate_shard(wrong_count, 0, 28)

        wrong_mode = shard_record(0, 28, 0)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "identity/counts"):
            RUNNER.validate_shard(wrong_mode, 0, 28)

    def test_shard_identity_is_not_interchangeable(self) -> None:
        value = shard_record(0, 28)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "identity/counts"):
            RUNNER.validate_shard(value, 1, 28)

    def test_ordinal_digest_is_independently_reproducible(self) -> None:
        digests = RUNNER.expected_ordinal_digests(28)
        self.assertEqual(len(digests), 28)
        self.assertEqual(len(set(digests)), 28)
        value = shard_record(0, 28)
        value["ordinal_digest_fnv1a64"] = digests[0]
        RUNNER.validate_shard(value, 0, 28, digests[0])
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "identity/counts"):
            RUNNER.validate_shard(value, 0, 28, "0" * 16)


if __name__ == "__main__":
    unittest.main()
