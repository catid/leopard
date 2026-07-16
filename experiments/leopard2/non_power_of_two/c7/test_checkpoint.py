#!/usr/bin/env python3
"""Checkpoint and portability tests for C7 evidence manifest v3."""

from __future__ import annotations

import hashlib
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

    def test_optional_authenticated_peer_rejects_exploits(self) -> None:
        current_path = os.environ.get("LEO2_C7_TEST_MANIFEST")
        peer_path_text = os.environ.get("LEO2_C7_TEST_PEER_MANIFEST")
        peer_root_text = os.environ.get("LEO2_C7_TEST_PEER_ROOT")
        if not current_path or not peer_path_text or not peer_root_text:
            self.skipTest("C7 peer evidence environment is not set")
        current = json.loads(pathlib.Path(current_path).read_text(encoding="utf-8"))
        peer_path = pathlib.Path(peer_path_text).resolve()
        peer_root = pathlib.Path(peer_root_text).resolve()
        peer = json.loads(peer_path.read_text(encoding="utf-8"))
        run_matrix.authenticate_peer(current, peer_path, peer_root)

        def rejected(candidate: dict) -> None:
            research = peer_root / ".research"
            research.mkdir(exist_ok=True)
            with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".json", prefix="c7-forged-peer-",
                    dir=research, encoding="utf-8", delete=False) as stream:
                json.dump(candidate, stream, indent=2, sort_keys=True)
                stream.write("\n")
                forged_path = pathlib.Path(stream.name)
            try:
                with self.assertRaises(RuntimeError):
                    run_matrix.authenticate_peer(current, forged_path, peer_root)
            finally:
                forged_path.unlink(missing_ok=True)

        redirected = json.loads(json.dumps(peer))
        for build in redirected["builds"]:
            for key in ("library", "executable"):
                build[key]["path"] = "/usr/bin/true"
        rejected(redirected)

        wrong_bytes = json.loads(json.dumps(peer))
        wrong_bytes["builds"][0]["library"]["bytes"] += 1
        rejected(wrong_bytes)

        for mutate in (
            lambda item: item.update(status="failed"),
            lambda item: item.update(runs=[]),
            lambda item: item["builds"][0].update(source_closure=[]),
            lambda item: item.update(padding=[{
                "path": "fake", "bytes": 0, "sha256": "0" * 64,
                "source_root_tokens": 0,
            }] * 64),
            lambda item: item.update(core_git_sha="0" * 40),
            lambda item: item["runner"].update(sha256="0" * 64),
        ):
            candidate = json.loads(json.dumps(peer))
            mutate(candidate)
            rejected(candidate)

    def test_optional_external_child_semantic_mutations(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_PEER_MANIFEST")
        root_text = os.environ.get("LEO2_C7_TEST_PEER_ROOT")
        if not requested or not root_text:
            self.skipTest("C7 peer evidence environment is not set")
        peer = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        peer_root = pathlib.Path(root_text).resolve()

        for run_index, key, value in (
            (0, "build", "auto"),
            (0, "kind", "authoritative"),
            (4, "environment", {
                "LC_ALL": "C", "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
                "ASAN_OPTIONS": "halt_on_error=1",
                "UBSAN_OPTIONS": "print_stacktrace=1",
            }),
        ):
            candidate = json.loads(json.dumps(peer))
            candidate["runs"][run_index][key] = value
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(
                    candidate, source_root=peer_root)

        def rehashed_child_rejected(run_index: int, mutation) -> None:
            candidate = json.loads(json.dumps(peer))
            record = candidate["runs"][run_index]["result"]
            path = peer_root / record["path"]
            original = path.read_bytes()
            child = json.loads(original)
            mutation(child)
            forged = (json.dumps(child, indent=2, sort_keys=True) + "\n").encode()
            try:
                path.write_bytes(forged)
                record["bytes"] = len(forged)
                record["sha256"] = hashlib.sha256(forged).hexdigest()
                with self.assertRaises(ValueError):
                    validate_evidence.validate_manifest(
                        candidate, source_root=peer_root)
            finally:
                path.write_bytes(original)

        for mutation in (
            lambda child: child.update(runtime_backend="scalar"),
            lambda child: child.update(profile={}),
            lambda child: child.update(timing_scope="benchmark"),
            lambda child: child.update(allocation_tracking="none"),
            lambda child: child.update(omp_num_threads="2"),
            lambda child: child.update(extra_key=True),
        ):
            rehashed_child_rejected(4, mutation)
        rehashed_child_rejected(
            5, lambda child: child["benchmarks"][0][
                "exact_encode_samples_us"].pop())
        rehashed_child_rejected(
            5, lambda child: child["benchmarks"][0]["exact_encode"].update(
                median_us=0))
        rehashed_child_rejected(
            5, lambda child: child["benchmarks"][0].update(extra_key=True))


if __name__ == "__main__":
    unittest.main(verbosity=2)
