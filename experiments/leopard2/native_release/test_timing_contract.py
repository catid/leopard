#!/usr/bin/env python3
"""Synthetic payload/collector and frozen-input mutations; no benchmark clocks."""
import copy
import hashlib
import os
from pathlib import Path
import tempfile
import unittest

from check_encode import require
from test_check_encode import ContractTests
from frozen_inputs import FrozenInputs, TIMING_ARTIFACTS
import timing_contract as contract
from run_timing import collect


def fixture_expected():
    fixture = ContractTests()
    result = []
    for cell in range(8):
        row = {"main": fixture.fixture(cell), "current": fixture.fixture(cell, False)}
        row["current"]["codec_commit"] = contract.CANDIDATE
        result.append(row)
    return result


def timed_record(expected, cell, implementation, mode="measure", ns=100000000):
    record = copy.deepcopy(expected[cell][implementation])
    record["schema"] = "leopard-native-release-encode-timing/v1"
    total = 6 + 9 * contract.GROUPS[cell]
    record["public_encode_calls"] = total
    record["grouped_timing"] = {
        "schema": "native-encode-groups/v1", "mode": mode,
        "clock_kind": "steady" if mode == "measure" else "synthetic" if mode == "synthetic" else "abort",
        "group_calls": contract.GROUPS[cell], "warmup_calls": 4, "total_public_calls": total,
        "clock_calls": 0 if mode == "exercise" else 18,
        "elapsed_ns": [ns] * 9 if mode == "measure" else
                      [31000000 + i * 100 for i in range(9)] if mode == "synthetic" else [0] * 9,
        "public_calls_at_clock": [6 + ((i + 1) // 2) * contract.GROUPS[cell] for i in range(18)]
                                  if mode == "synthetic" else [],
    }
    return record


def fixture_rows(expected):
    return [{"cell": cell, "round": round_id, "comparison": comparison, "slot": slot,
             "implementation": implementation, "sibling_delta": 0, "cpu_delta": 10,
             "record": timed_record(expected, cell, implementation,
                                    ns=100000000 if implementation == "main" else 50000000)}
            for cell, round_id, comparison, slot, implementation in contract.schedule()]


class TimingTests(unittest.TestCase):
    def setUp(self):
        self.expected = fixture_expected()

    def test_all_modes_and_cells(self):
        for cell in range(8):
            for implementation in ("main", "current"):
                for mode in ("measure", "synthetic", "exercise"):
                    record = timed_record(self.expected, cell, implementation, mode)
                    result = contract.validate_timing(record, cell, implementation,
                                                      self.expected[cell][implementation], mode)
                    self.assertEqual(len(result), 9)

    def test_timing_fields_bound(self):
        record = timed_record(self.expected, 0, "main")
        for key in record["grouped_timing"]:
            bad = copy.deepcopy(record)
            bad["grouped_timing"][key] = None
            with self.subTest(key=key), self.assertRaises(ValueError):
                contract.validate_timing(bad, 0, "main", self.expected[0]["main"])

    def test_minimum_window_and_call_counts(self):
        for value in (19999999, 0, -1, True, float("nan"), 9007199254740992):
            record = timed_record(self.expected, 0, "main")
            record["grouped_timing"]["elapsed_ns"][0] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                contract.validate_timing(record, 0, "main", self.expected[0]["main"])
        record = timed_record(self.expected, 0, "main", ns=20000000)
        contract.validate_timing(record, 0, "main", self.expected[0]["main"])
        record["public_encode_calls"] -= 1
        with self.assertRaises(ValueError):
            contract.validate_timing(record, 0, "main", self.expected[0]["main"])

    def test_synthetic_endpoint_and_wrong_mode(self):
        record = timed_record(self.expected, 5, "current", "synthetic")
        record["grouped_timing"]["public_calls_at_clock"][1] -= 1
        with self.assertRaises(ValueError):
            contract.validate_timing(record, 5, "current", self.expected[5]["current"], "synthetic")
        with self.assertRaises(ValueError):
            contract.validate_timing(record, 5, "current", self.expected[5]["current"], "unknown")

    def test_abba_success_and_aggregate_control_failure(self):
        rows = fixture_rows(self.expected)
        result = contract.analyze(rows, self.expected)
        self.assertEqual(result["decision"], "complete")
        self.assertTrue(all(cell["ratios"]["native_vs_current"] == 2 for cell in result["cells"]))
        self.assertTrue(all(cell["classification"] == "current_advantage" for cell in result["cells"]))
        rows[4]["record"]["grouped_timing"]["elapsed_ns"] = [150000000] * 9
        result = contract.analyze(rows, self.expected)
        self.assertEqual(result["decision"], "inconclusive_controls")
        self.assertTrue(all(cell["classification"] == "inconclusive_controls" for cell in result["cells"]))
        self.assertFalse(result["promotion"])

    def test_order_partial_and_contamination_rejected(self):
        rows = fixture_rows(self.expected)
        bad_rows = [rows[:-1], rows[1:] + rows[:1]]
        for key, value in (("sibling_delta", 1), ("cell", True), ("slot", 3)):
            bad = copy.deepcopy(rows)
            bad[0][key] = value
            bad_rows.append(bad)
        for bad in bad_rows:
            with self.assertRaises(ValueError):
                contract.analyze(bad, self.expected)

    def test_collector_complete_schedule(self):
        state = {"preflight": [], "invocations": [], "complete": False, "analysis": None}
        launches = []
        checkpoints = []
        def launch(implementation, cell, measured, label):
            launches.append((implementation, cell, measured, label))
            if label == "passive": return None
            value = timed_record(self.expected, cell, implementation,
                                 ns=100000000 if implementation == "main" else 50000000) if measured else self.expected[cell][implementation]
            return value, {"sibling_delta": 0, "cpu_delta": 10}
        collect(self.expected, launch, state,
                lambda value: checkpoints.append((len(value["preflight"]), len(value["invocations"]))))
        self.assertEqual(len(launches), 305)
        self.assertEqual(launches[16][3], "passive")
        self.assertEqual(len(state["preflight"]), 16)
        self.assertEqual(len(state["invocations"]), 288)
        self.assertTrue(state["complete"])
        self.assertEqual(checkpoints[-1], (16, 288))

    def test_collector_stops_without_analysis(self):
        state = {"preflight": [], "invocations": [], "complete": False, "analysis": None}
        def launch(implementation, cell, measured, label):
            if label == "passive": return None
            if measured:
                return timed_record(self.expected, cell, implementation), {"sibling_delta": 1, "cpu_delta": 10}
            return self.expected[cell][implementation], {"sibling_delta": 0, "cpu_delta": 1}
        with self.assertRaises(ValueError):
            collect(self.expected, launch, state, lambda value: None)
        self.assertEqual(len(state["invocations"]), 1)
        self.assertFalse(state["complete"])
        self.assertIsNone(state["analysis"])


class FrozenTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.lines = []
        for name in sorted(TIMING_ARTIFACTS):
            value = name.encode()
            (self.root / name).write_bytes(value)
            (self.root / name).chmod(0o444)
            self.lines.append(f"{hashlib.sha256(value).hexdigest()}  {name}")
        self.manifest(self.lines)

    def manifest(self, lines):
        path = self.root / "SHA256SUMS"
        if path.exists(): path.chmod(0o644)
        path.write_text("\n".join(lines) + "\n")
        path.chmod(0o444)

    def test_correct_and_unknown_execution(self):
        frozen = FrozenInputs(self.root, TIMING_ARTIFACTS)
        frozen.executable("main-steady")
        with self.assertRaises(ValueError): frozen.executable("unknown")

    def test_missing_replaced_duplicate_and_bad_hash(self):
        variants = (self.lines[:-1], self.lines + self.lines[:1],
                    [self.lines[0].replace("current-abort", "unrelated"), *self.lines[1:]],
                    ["not-a-sha  current-abort", *self.lines[1:]])
        for variant in variants:
            self.manifest(variant)
            with self.assertRaises(ValueError): FrozenInputs(self.root, TIMING_ARTIFACTS)

    def test_manifest_and_file_changes(self):
        frozen = FrozenInputs(self.root, TIMING_ARTIFACTS)
        self.manifest(list(reversed(self.lines)))
        with self.assertRaises(ValueError): frozen.verify()
        self.manifest(self.lines)
        frozen = FrozenInputs(self.root, TIMING_ARTIFACTS)
        target = self.root / "main-steady"
        target.chmod(0o644)
        target.write_bytes(b"changed")
        target.chmod(0o444)
        with self.assertRaises(ValueError): frozen.verify()


if __name__ == "__main__":
    unittest.main()
