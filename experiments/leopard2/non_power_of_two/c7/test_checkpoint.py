#!/usr/bin/env python3
"""Checkpoint and portability tests for C7 evidence manifest v3."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

import run_matrix
import validate_evidence
from test_evidence_portability import PortabilityTests  # noqa: F401


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE / "results"


class CheckpointTests(unittest.TestCase):
    def test_algebra_canonical_regeneration(self) -> None:
        retained = RESULTS / "algebra.json"
        self.assertTrue(retained.is_file())
        with tempfile.TemporaryDirectory(prefix="leopard2-c7-algebra-") as directory:
            regenerated = pathlib.Path(directory) / "algebra.json"
            subprocess.run([
                sys.executable, "-X", "dev", str(HERE / "algebra.py"),
                "--output", str(regenerated),
            ], cwd=ROOT, check=True, env={
                **dict(os.environ),
                "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
            })
            self.assertEqual(regenerated.read_bytes(), retained.read_bytes())

    def test_retained_manifest_is_current_v3_comparison_evidence(self) -> None:
        manifest_path = RESULTS / "build-run-manifest.json"
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(data.get("schema"),
                         "leopard2-c7-build-run-manifest/v3")
        self.assertEqual(
            data.get("reproducibility", {}).get("comparison", {}).get("status"),
            "pass")
        validate_evidence.validate_manifest(data)
        leaked = json.loads(json.dumps(data))
        compile_argv = leaked["builds"][0]["compile_argv"]
        index = next(i for i, value in enumerate(compile_argv)
                     if run_matrix.NORMALIZATION_TOKEN in value)
        compile_argv[index] = compile_argv[index].replace(
            run_matrix.NORMALIZATION_TOKEN, str(ROOT), 1)
        leaked["builds"][0]["argv_source_root_tokens"]["compile"] -= 1
        with self.assertRaises(ValueError):
            validate_evidence.validate_manifest(leaked)
        missing_token = json.loads(json.dumps(data))
        compile_argv = missing_token["builds"][0]["compile_argv"]
        index = next(i for i, value in enumerate(compile_argv)
                     if run_matrix.NORMALIZATION_TOKEN in value)
        compile_argv[index] = compile_argv[index].replace(
            run_matrix.NORMALIZATION_TOKEN, "LEO2_SOURCE_ROOT", 1)
        missing_token["builds"][0]["argv_source_root_tokens"]["compile"] -= 1
        with self.assertRaises(ValueError):
            validate_evidence.validate_manifest(missing_token)
        failed_comparison = json.loads(json.dumps(data))
        failed_comparison["reproducibility"]["comparison"]["status"] = "failed"
        with self.assertRaises(ValueError):
            validate_evidence.validate_manifest(failed_comparison)

    def test_current_attestation_constants_are_exact(self) -> None:
        self.assertEqual(
            sum(item["asan_lines"] for item in
                run_matrix.EXPECTED_ARCHIVE_MEMBER_COUNTS.values()), 329)
        self.assertEqual(
            sum(item["ubsan_lines"] for item in
                run_matrix.EXPECTED_ARCHIVE_MEMBER_COUNTS.values()), 87)
        self.assertEqual(
            run_matrix.EXPECTED_EXECUTABLE_SANITIZER_COUNTS,
            {"asan_lines": 320, "ubsan_lines": 54})
        self.assertEqual(
            validate_evidence.EXPECTED_CORRECTNESS[
                "decode_read_only_input_alias_calls"], 117)
        self.assertEqual(
            validate_evidence.EXPECTED_CORRECTNESS[
                "decode_read_only_input_alias_symbol_comparisons"], 6025)

    def test_optional_external_v3_manifest(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        path = pathlib.Path(requested)
        data = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(
            data.get("reproducibility", {}).get("comparison", {}).get("status"),
            "pass")
        validate_evidence.validate_manifest(data)


if __name__ == "__main__":
    unittest.main(verbosity=2)
