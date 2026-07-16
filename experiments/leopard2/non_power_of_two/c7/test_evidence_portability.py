#!/usr/bin/env python3
"""Unit tests for C7 v3 evidence portability primitives."""

from __future__ import annotations

import copy
import hashlib
import pathlib
import subprocess
import tempfile
import unittest

import run_matrix
import validate_evidence


class PortabilityTests(unittest.TestCase):
    def test_normalizer_replaces_only_exact_root_prefix(self) -> None:
        root = str(run_matrix.ROOT)
        value = f"{root}/one {root} {root}=mapped {root}-sibling/two"
        normalized, count = run_matrix.normalize_source_root(value)
        self.assertEqual(count, 3)
        self.assertEqual(
            normalized,
            "${LEO2_SOURCE_ROOT}/one ${LEO2_SOURCE_ROOT} "
            "${LEO2_SOURCE_ROOT}=mapped "
            f"{root}-sibling/two")

    def test_affinity_defaults_and_explicit_validation(self) -> None:
        self.assertEqual(
            run_matrix.default_cpus({19, 2, 11, 5, 7, 13}),
            [2, 5, 7, 11, 13])
        run_matrix.validate_cpus([2, 5, 7, 11, 13], 19,
                                 {19, 2, 11, 5, 7, 13})
        for cpus, smoke, allowed in (
            ([2, 5, 7, 11], 2, {2, 5, 7, 11}),
            ([2, 5, 7, 11, 11], 2, {2, 5, 7, 11}),
            ([2, 5, 7, 11, 13], 99, {2, 5, 7, 11, 13}),
        ):
            with self.assertRaises(ValueError):
                run_matrix.validate_cpus(cpus, smoke, allowed)
        with self.assertRaises(ValueError):
            run_matrix.default_cpus({1, 3, 5, 7})

    def test_normalized_artifact_schema_tokens_and_path_leaks(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="c7-portability-", dir=run_matrix.ROOT) as directory:
            path = pathlib.Path(directory) / "log.txt"
            relative = path.relative_to(run_matrix.ROOT).as_posix()

            def record(contents: str, tokens: int) -> dict:
                path.write_text(contents, encoding="utf-8")
                data = path.read_bytes()
                return {
                    "path": relative, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "source_root_tokens": tokens,
                }

            valid = record("${LEO2_SOURCE_ROOT}/CMakeLists.txt\n", 1)
            validate_evidence.validate_normalized_text(
                valid, "valid", run_matrix.ROOT, require_token=True)
            missing = record("CMakeLists.txt\n", 0)
            with self.assertRaises(ValueError):
                validate_evidence.validate_normalized_text(
                    missing, "missing", run_matrix.ROOT, require_token=True)
            wrong_count = record("${LEO2_SOURCE_ROOT}/CMakeLists.txt\n", 0)
            with self.assertRaises(ValueError):
                validate_evidence.validate_normalized_text(
                    wrong_count, "count", run_matrix.ROOT, require_token=True)
            leaked = record(f"{run_matrix.ROOT}/CMakeLists.txt\n", 0)
            with self.assertRaises(ValueError):
                validate_evidence.validate_normalized_text(
                    leaked, "leak", run_matrix.ROOT, require_token=False)
            absolute = dict(valid, path=str(path))
            with self.assertRaises(ValueError):
                validate_evidence.validate_normalized_text(
                    absolute, "absolute", run_matrix.ROOT, require_token=True)

    def test_portable_program_record_does_not_require_or_execute_path(self) -> None:
        missing = {
            "path": "/definitely/not/installed/c7-tool",
            "sha256": "1" * 64,
            "version": "recorded tool 1.0\n",
        }
        self.assertEqual(
            validate_evidence.validate_program_record(missing, "missing"),
            pathlib.Path(missing["path"]))
        with self.assertRaises(ValueError):
            validate_evidence.validate_program_live(missing, "missing")

        true = pathlib.Path("/usr/bin/true")
        if not true.is_file():
            self.skipTest("/usr/bin/true is unavailable")
        version = subprocess.run(
            [str(true), "--version"], check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout
        live = {
            "path": str(true), "sha256": validate_evidence.sha256(true),
            "version": version,
        }
        self.assertEqual(
            validate_evidence.validate_program_live(live, "true"), true)
        forged = dict(live, sha256="0" * 64)
        with self.assertRaises(ValueError):
            validate_evidence.validate_program_live(forged, "true")

    def test_stale_unretained_output_is_ignored_only_portably(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="c7-stale-", dir=run_matrix.ROOT) as directory:
            path = pathlib.Path(directory) / "liblibleopard.a"
            path.write_bytes(b"stale")
            record = {
                "path": path.relative_to(run_matrix.ROOT).as_posix(),
                "bytes": 4,
                "sha256": hashlib.sha256(b"good").hexdigest(),
            }
            validate_evidence.validate_artifact(
                record, "portable stale output", run_matrix.ROOT,
                required=False, check_if_present=False)
            with self.assertRaises(ValueError):
                validate_evidence.validate_artifact(
                    record, "live stale output", run_matrix.ROOT,
                    required=True, check_if_present=True)

    def test_argv_token_and_checkout_path_mutations_rejected(self) -> None:
        argv = [
            "cc", "-I${LEO2_SOURCE_ROOT}",
            "${LEO2_SOURCE_ROOT}/LeopardFF8.cpp"]
        validate_evidence.validate_argv(
            argv, 2, "argv", require_token=True)
        with self.assertRaises(ValueError):
            validate_evidence.validate_argv(argv, 1, "argv", require_token=True)
        leaked = ["cc", f"-I{run_matrix.ROOT}",
                  f"{run_matrix.ROOT}/LeopardFF8.cpp"]
        with self.assertRaises(ValueError):
            validate_evidence.validate_argv(
                leaked, 0, "argv", require_token=False)

    def test_reproducibility_fingerprint_mutation_rejected(self) -> None:
        artifacts = {
            name: {
                "library_sha256": str(index) * 64,
                "executable_sha256": format(index + 8, "x") * 64,
            }
            for index, name in enumerate(run_matrix.BUILD_NAMES, start=1)
        }
        base = {
            "schema": "leopard2-c7-build-run-manifest/v3",
            "core_git_sha": "a" * 40,
            "source": {"path": "s", "sha256": "b" * 64, "bytes": 1},
            "runner": {"path": "r", "sha256": "c" * 64, "bytes": 1},
            "validator": {"path": "v", "sha256": "d" * 64, "bytes": 1},
            "reproducibility": {"fingerprints": artifacts},
        }
        run_matrix.require_reproducible_peer(base, copy.deepcopy(base))
        forged = copy.deepcopy(base)
        forged["reproducibility"]["fingerprints"]["scalar"][
            "library_sha256"] = "f" * 64
        with self.assertRaises(RuntimeError):
            run_matrix.require_reproducible_peer(base, forged)


if __name__ == "__main__":
    unittest.main(verbosity=2)
