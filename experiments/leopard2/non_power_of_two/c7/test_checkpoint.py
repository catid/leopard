#!/usr/bin/env python3
"""Checkpoint and portability tests for C7 evidence manifest v3."""

from __future__ import annotations

import hashlib
import gzip
import io
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import tarfile
import threading
import unittest
from unittest import mock

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
        candidate = json.loads(json.dumps(data))
        record = candidate["reproducibility"]["comparison"]["peer_attestation"]
        attestation_path = ROOT / record["path"]
        original = attestation_path.read_bytes()
        report = json.loads(original)
        report["checks"]["live_tools_and_outputs"] = "failed"
        forged = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
        try:
            attestation_path.write_bytes(forged)
            record["bytes"] = len(forged)
            record["sha256"] = hashlib.sha256(forged).hexdigest()
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(candidate)
        finally:
            attestation_path.write_bytes(original)

    def test_optional_retained_peer_bundle_rejects_exploits(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        manifest_path = pathlib.Path(requested).resolve()
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        comparison = data["reproducibility"]["comparison"]

        for key in ("peer_manifest", "peer_evidence_bundle", "peer_attestation"):
            candidate = json.loads(json.dumps(data))
            candidate["reproducibility"]["comparison"][key]["sha256"] = "0" * 64
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(candidate)

            record = comparison[key]
            path = ROOT / record["path"]
            original = path.read_bytes()
            path.unlink()
            try:
                with self.assertRaises(ValueError):
                    validate_evidence.validate_manifest(data)
            finally:
                path.write_bytes(original)

        def artifact_bytes_rejected(key: str, forged: bytes) -> None:
            candidate = json.loads(json.dumps(data))
            record = candidate["reproducibility"]["comparison"][key]
            path = ROOT / record["path"]
            original = path.read_bytes()
            try:
                path.write_bytes(forged)
                record["bytes"] = len(forged)
                record["sha256"] = hashlib.sha256(forged).hexdigest()
                if key == "peer_manifest":
                    candidate["reproducibility"]["comparison"][
                        "peer_manifest_sha256"] = record["sha256"]
                with self.assertRaises(ValueError):
                    validate_evidence.validate_manifest(candidate)
            finally:
                path.write_bytes(original)

        artifact_bytes_rejected("peer_manifest", b"{}\n")
        artifact_bytes_rejected("peer_evidence_bundle", b"not a gzip stream\n")
        artifact_bytes_rejected("peer_attestation", b"{}\n")
        # A main comparison manifest cannot be substituted as its own allegedly
        # independent peer, even when its new exact hash is recorded.
        artifact_bytes_rejected("peer_manifest", manifest_path.read_bytes())

        bundle_path = ROOT / comparison["peer_evidence_bundle"]["path"]
        valid_bundle = bundle_path.read_bytes()
        raw_tar = gzip.decompress(valid_bundle)
        artifact_bytes_rejected(
            "peer_evidence_bundle", gzip.compress(raw_tar, mtime=1))

        def unsafe_archive(
            infos: list[tuple[tarfile.TarInfo, bytes]], *,
            canonical_metadata: bool = True,
        ) -> bytes:
            output = io.BytesIO()
            with gzip.GzipFile(
                    filename="", mode="wb", fileobj=output, mtime=0
            ) as compressed:
                with tarfile.open(fileobj=compressed, mode="w|") as archive:
                    for info, payload in infos:
                        info.size = len(payload)
                        if canonical_metadata:
                            info.mode = 0o444
                            info.mtime = 0
                            info.uid = info.gid = 0
                            info.uname = info.gname = ""
                        archive.addfile(info, io.BytesIO(payload))
            return output.getvalue()

        traversal = tarfile.TarInfo("../escape")
        artifact_bytes_rejected(
            "peer_evidence_bundle", unsafe_archive([(traversal, b"x")]))
        for kind in (
                tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.CHRTYPE,
                tarfile.BLKTYPE, tarfile.FIFOTYPE):
            unsafe = tarfile.TarInfo("files/unsafe")
            unsafe.type = kind
            unsafe.linkname = "../escape"
            artifact_bytes_rejected(
                "peer_evidence_bundle", unsafe_archive([(unsafe, b"")]))
        duplicate_a = tarfile.TarInfo("index.json")
        duplicate_b = tarfile.TarInfo("index.json")
        artifact_bytes_rejected(
            "peer_evidence_bundle", unsafe_archive([
                (duplicate_a, b"{}\n"), (duplicate_b, b"{}\n")]))
        oversized = tarfile.TarInfo("files/oversized")
        oversized.size = run_matrix.PEER_BUNDLE_MAX_MEMBER_BYTES + 1
        oversized.mode = 0o444
        oversized.mtime = 0
        oversized.uid = oversized.gid = 0
        oversized.uname = oversized.gname = ""
        raw_oversized = oversized.tobuf(format=tarfile.USTAR_FORMAT) + b"\0" * 1024
        compressed_oversized = io.BytesIO()
        with gzip.GzipFile(
                filename="", mode="wb", fileobj=compressed_oversized, mtime=0
        ) as compressed:
            compressed.write(raw_oversized)
        artifact_bytes_rejected(
            "peer_evidence_bundle", compressed_oversized.getvalue())

        peer_path = ROOT / comparison["peer_manifest"]["path"]
        peer = json.loads(peer_path.read_text(encoding="utf-8"))
        expected = validate_evidence._peer_bundle_records(peer)
        files = validate_evidence._read_peer_bundle(valid_bundle, expected)
        index = {
            "schema": run_matrix.PEER_BUNDLE_SCHEMA,
            "files": [
                {
                    "path": path, "bytes": record["bytes"],
                    "sha256": record["sha256"],
                }
                for path, record in sorted(expected.items())
            ],
        }
        indexed_members = {
            "index.json": (
                json.dumps(index, indent=2, sort_keys=True) + "\n").encode(),
            **{f"files/{path}": contents for path, contents in files.items()},
        }
        artifact_bytes_rejected(
            "peer_evidence_bundle", unsafe_archive([
                (tarfile.TarInfo(name), contents)
                for name, contents in sorted(indexed_members.items())
            ], canonical_metadata=False))
        expected_payload = sum(item["bytes"] for item in expected.values())
        expansion_limit = min(
            run_matrix.PEER_BUNDLE_MAX_ARCHIVE_BYTES,
            expected_payload + (len(expected) + 8) * 2048 + 1024 * 1024)
        artifact_bytes_rejected(
            "peer_evidence_bundle",
            gzip.compress(b"\0" * (expansion_limit + 1), mtime=0))

    def test_optional_verified_peer_artifacts_are_single_snapshots(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        data = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        comparison = data["reproducibility"]["comparison"]
        real_read = validate_evidence._read_file_once
        for key, forged in (
            ("peer_manifest", b"{}\n"),
            ("peer_evidence_bundle", b"not a gzip stream\n"),
            ("peer_attestation", b"{}\n"),
        ):
            path = (ROOT / comparison[key]["path"]).resolve()
            original = path.read_bytes()
            swapped = threading.Event()

            def read_then_swap(
                candidate: pathlib.Path, maximum_bytes: int, *,
                target: pathlib.Path = path, replacement: bytes = forged,
            ) -> bytes:
                contents = real_read(candidate, maximum_bytes)
                if candidate.resolve() == target and not swapped.is_set():
                    worker = threading.Thread(
                        target=target.write_bytes, args=(replacement,))
                    worker.start()
                    worker.join()
                    swapped.set()
                return contents

            try:
                with mock.patch.object(
                        validate_evidence, "_read_file_once",
                        side_effect=read_then_swap):
                    validate_evidence.validate_manifest(data)
                self.assertTrue(swapped.is_set())
                self.assertEqual(path.read_bytes(), forged)
            finally:
                path.write_bytes(original)

    def test_optional_normalized_artifacts_are_verified_snapshots(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        data = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        record = data["builds"][0]["configure_stdout"]
        path = (ROOT / record["path"]).resolve()
        original = path.read_bytes()
        forged = original + b"\npost-read unauthenticated replacement\n"
        swapped = threading.Event()
        real_read = validate_evidence._read_file_once

        def read_then_swap(
            candidate: pathlib.Path, maximum_bytes: int,
        ) -> bytes:
            contents = real_read(candidate, maximum_bytes)
            if candidate.resolve() == path and not swapped.is_set():
                worker = threading.Thread(
                    target=path.write_bytes, args=(forged,))
                worker.start()
                worker.join()
                swapped.set()
            return contents

        try:
            with mock.patch.object(
                    validate_evidence, "_read_file_once",
                    side_effect=read_then_swap):
                validate_evidence.validate_manifest(data)
            self.assertTrue(swapped.is_set())
            self.assertEqual(path.read_bytes(), forged)
        finally:
            path.write_bytes(original)

    def test_optional_live_binary_scans_use_verified_snapshots(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        data = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        source_root = pathlib.Path(
            os.environ.get("LEO2_C7_TEST_PEER_ROOT", ROOT)).resolve()
        real_read = validate_evidence._read_file_once
        for output_key in ("library", "executable"):
            record = data["builds"][0][output_key]
            path = (ROOT / record["path"]).resolve()
            original = path.read_bytes()
            forged = b"binary replaced after authenticated read\n"
            swapped = False

            def read_then_swap(
                candidate: pathlib.Path, maximum_bytes: int,
            ) -> bytes:
                nonlocal swapped
                contents = real_read(candidate, maximum_bytes)
                if candidate.resolve() == path and not swapped:
                    path.write_bytes(forged)
                    swapped = True
                return contents

            try:
                with mock.patch.object(
                        validate_evidence, "_read_file_once",
                        side_effect=read_then_swap):
                    validate_evidence.validate_manifest(
                        data, source_root=source_root, evidence_root=ROOT,
                        live=True)
                self.assertTrue(swapped)
                self.assertEqual(path.read_bytes(), forged)
            finally:
                path.write_bytes(original)

    def test_optional_peer_binary_records_are_fully_cross_bound(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        main = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        comparison = main["reproducibility"]["comparison"]
        peer = json.loads(
            (ROOT / comparison["peer_manifest"]["path"]).read_text(
                encoding="utf-8"))
        bundle_bytes = (
            ROOT / comparison["peer_evidence_bundle"]["path"]).read_bytes()
        original_report = json.loads(
            (ROOT / comparison["peer_attestation"]["path"]).read_text(
                encoding="utf-8"))

        for build_index in range(len(run_matrix.BUILD_NAMES)):
            for output_key in ("library", "executable"):
                forged_peer = json.loads(json.dumps(peer))
                forged_peer["builds"][build_index][output_key]["bytes"] += 1
                peer_bytes = (
                    json.dumps(forged_peer, indent=2, sort_keys=True) +
                    "\n").encode()
                with tempfile.TemporaryDirectory(
                        prefix=".c7-binary-cross-bind-", dir=ROOT
                ) as directory:
                    root = pathlib.Path(directory)

                    def retain(name: str, contents: bytes) -> dict:
                        path = root / name
                        path.write_bytes(contents)
                        return {
                            "path": path.relative_to(ROOT).as_posix(),
                            "bytes": len(contents),
                            "sha256": hashlib.sha256(contents).hexdigest(),
                        }

                    peer_record = retain("peer-manifest.json", peer_bytes)
                    bundle_record = retain(
                        "peer-evidence-bundle.tar.gz", bundle_bytes)
                    report = json.loads(json.dumps(original_report))
                    report["peer_manifest_artifact"] = peer_record
                    report["peer_evidence_bundle"] = bundle_record
                    report_bytes = (
                        json.dumps(report, indent=2, sort_keys=True) +
                        "\n").encode()
                    attestation_record = retain(
                        "peer-reproducibility-attestation.json", report_bytes)
                    candidate = json.loads(json.dumps(main))
                    candidate_comparison = candidate[
                        "reproducibility"]["comparison"]
                    candidate_comparison["peer_manifest"] = peer_record
                    candidate_comparison["peer_manifest_sha256"] = (
                        peer_record["sha256"])
                    candidate_comparison["peer_evidence_bundle"] = bundle_record
                    candidate_comparison[
                        "peer_attestation"] = attestation_record
                    with self.assertRaises(ValueError):
                        validate_evidence.validate_manifest(candidate)

    def test_optional_retained_json_is_unambiguous_and_canonical(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_MANIFEST")
        if not requested:
            self.skipTest("LEO2_C7_TEST_MANIFEST is not set")
        main = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        comparison = main["reproducibility"]["comparison"]
        peer_bytes = (
            ROOT / comparison["peer_manifest"]["path"]).read_bytes()
        peer = json.loads(peer_bytes)
        bundle_bytes = (
            ROOT / comparison["peer_evidence_bundle"]["path"]).read_bytes()
        report_bytes = (
            ROOT / comparison["peer_attestation"]["path"]).read_bytes()
        report = json.loads(report_bytes)

        def forged_peer_rejected(forged_peer_bytes: bytes) -> None:
            with tempfile.TemporaryDirectory(
                    prefix=".c7-json-cross-bind-", dir=ROOT
            ) as directory:
                root = pathlib.Path(directory)

                def retain(name: str, contents: bytes) -> dict:
                    path = root / name
                    path.write_bytes(contents)
                    return {
                        "path": path.relative_to(ROOT).as_posix(),
                        "bytes": len(contents),
                        "sha256": hashlib.sha256(contents).hexdigest(),
                    }

                peer_record = retain("peer-manifest.json", forged_peer_bytes)
                bundle_record = retain(
                    "peer-evidence-bundle.tar.gz", bundle_bytes)
                forged_report = json.loads(json.dumps(report))
                forged_report["peer_manifest_artifact"] = peer_record
                forged_report["peer_evidence_bundle"] = bundle_record
                attestation_record = retain(
                    "peer-reproducibility-attestation.json",
                    (json.dumps(forged_report, indent=2, sort_keys=True) +
                     "\n").encode())
                candidate = json.loads(json.dumps(main))
                candidate_comparison = candidate[
                    "reproducibility"]["comparison"]
                candidate_comparison["peer_manifest"] = peer_record
                candidate_comparison[
                    "peer_manifest_sha256"] = peer_record["sha256"]
                candidate_comparison["peer_evidence_bundle"] = bundle_record
                candidate_comparison[
                    "peer_attestation"] = attestation_record
                with self.assertRaises(ValueError):
                    validate_evidence.validate_manifest(candidate)

        self.assertTrue(peer_bytes.startswith(b"{\n"))
        forged_peer_rejected(
            b'{\n  "status": "failed",\n' + peer_bytes[2:])
        forged_peer_rejected(
            json.dumps(peer, separators=(",", ":")).encode())

        attestation_path = ROOT / comparison["peer_attestation"]["path"]
        duplicate_attestation = (
            b'{\n  "status": "failed",\n' + report_bytes[2:])
        candidate = json.loads(json.dumps(main))
        candidate_record = candidate["reproducibility"]["comparison"][
            "peer_attestation"]
        try:
            attestation_path.write_bytes(duplicate_attestation)
            candidate_record["bytes"] = len(duplicate_attestation)
            candidate_record["sha256"] = hashlib.sha256(
                duplicate_attestation).hexdigest()
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(candidate)
        finally:
            attestation_path.write_bytes(report_bytes)

        result_record = main["runs"][0]["result"]
        result_path = ROOT / result_record["path"]
        result_bytes = result_path.read_bytes()
        self.assertTrue(result_bytes.startswith(b"{\n"))
        duplicate_result = b'{\n  "status":"failed",\n' + result_bytes[2:]
        candidate = json.loads(json.dumps(main))
        candidate_record = candidate["runs"][0]["result"]
        try:
            result_path.write_bytes(duplicate_result)
            candidate_record["bytes"] = len(duplicate_result)
            candidate_record["sha256"] = hashlib.sha256(
                duplicate_result).hexdigest()
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(candidate)
        finally:
            result_path.write_bytes(result_bytes)

    def test_optional_authenticated_peer_rejects_exploits(self) -> None:
        current_path = os.environ.get("LEO2_C7_TEST_MANIFEST")
        peer_path_text = os.environ.get("LEO2_C7_TEST_PEER_MANIFEST")
        peer_root_text = os.environ.get("LEO2_C7_TEST_PEER_ROOT")
        if not current_path or not peer_path_text or not peer_root_text:
            self.skipTest("C7 peer evidence environment is not set")
        current = json.loads(pathlib.Path(current_path).read_text(encoding="utf-8"))
        peer_path = pathlib.Path(peer_path_text).resolve()
        peer_root = pathlib.Path(peer_root_text).resolve()
        peer_bytes = peer_path.read_bytes()
        peer = json.loads(peer_bytes)

        def authenticated(candidate: dict) -> None:
            encoded = (json.dumps(candidate, indent=2, sort_keys=True) +
                       "\n").encode()
            with tempfile.TemporaryDirectory(
                    prefix=".c7-peer-auth-", dir=ROOT
            ) as directory:
                root = pathlib.Path(directory)
                manifest = root / "peer.json"
                manifest.write_bytes(encoded)
                evidence = root / "evidence"
                run_matrix.capture_peer_snapshot(candidate, peer_root, evidence)
                run_matrix.authenticate_peer(
                    current, encoded, manifest, peer_root, evidence)

        authenticated(peer)

        def rejected(candidate: dict) -> None:
            with self.assertRaises(RuntimeError):
                authenticated(candidate)

        redirected = json.loads(json.dumps(peer))
        for build in redirected["builds"]:
            for key in ("library", "executable"):
                build[key]["path"] = "/usr/bin/true"
        rejected(redirected)

        true_path = pathlib.Path("/usr/bin/true")
        if true_path.is_file():
            program_redirect = json.loads(json.dumps(peer))
            true_record = {
                "path": str(true_path),
                "sha256": validate_evidence.sha256(true_path),
                "version": subprocess.run(
                    [str(true_path), "--version"], check=True, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout,
            }
            program_redirect["taskset"] = true_record
            for run in program_redirect["runs"]:
                run["argv"][0] = str(true_path)
            rejected(program_redirect)

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

        # Once captured, peer-controlled manifest, log, and raw dependency
        # mutations cannot affect either portable or live authentication.
        with tempfile.TemporaryDirectory(
                prefix=".c7-peer-toctou-", dir=ROOT
        ) as directory:
            snapshot = pathlib.Path(directory)
            snapshot_manifest = snapshot / "peer.json"
            snapshot_manifest.write_bytes(peer_bytes)
            evidence = snapshot / "evidence"
            run_matrix.capture_peer_snapshot(peer, peer_root, evidence)
            log_path = peer_root / peer["builds"][0]["configure_stdout"]["path"]
            dependency_path = (
                peer_root /
                peer["builds"][0]["dependency_files"][0]["build_path"])
            originals = {
                peer_path: peer_path.read_bytes(),
                log_path: log_path.read_bytes(),
                dependency_path: dependency_path.read_bytes(),
            }
            try:
                peer_path.write_bytes(b"{}\n")
                log_path.write_bytes(b"peer log changed after capture\n")
                dependency_path.write_bytes(
                    b"peer dependency changed after capture\n")
                run_matrix.authenticate_peer(
                    current, peer_bytes, snapshot_manifest, peer_root, evidence)
            finally:
                for path, contents in originals.items():
                    path.write_bytes(contents)

    def test_optional_external_child_semantic_mutations(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_PEER_MANIFEST")
        root_text = os.environ.get("LEO2_C7_TEST_PEER_ROOT")
        if not requested or not root_text:
            self.skipTest("C7 peer evidence environment is not set")
        peer = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        peer_root = pathlib.Path(root_text).resolve()
        with self.assertRaises(ValueError):
            validate_evidence.validate_manifest(
                peer, source_root=peer_root / "experiments")

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

    def test_optional_external_closure_and_job_mutations(self) -> None:
        requested = os.environ.get("LEO2_C7_TEST_PEER_MANIFEST")
        root_text = os.environ.get("LEO2_C7_TEST_PEER_ROOT")
        if not requested or not root_text:
            self.skipTest("C7 peer evidence environment is not set")
        peer = json.loads(pathlib.Path(requested).read_text(encoding="utf-8"))
        peer_root = pathlib.Path(root_text).resolve()

        closure = peer["builds"][0]["source_closure"]
        self.assertEqual(
            [item["path"] for item in closure],
            list(run_matrix.EXPECTED_SOURCE_CLOSURE))
        for index in range(len(run_matrix.EXPECTED_SOURCE_CLOSURE)):
            missing = json.loads(json.dumps(peer))
            missing["builds"][0]["source_closure"].pop(index)
            with self.assertRaises(ValueError):
                validate_evidence.validate_source_closure_paths(
                    missing["builds"][0]["source_closure"])

            substituted = json.loads(json.dumps(peer))
            replacement = closure[(index + 1) % len(closure)]
            substituted["builds"][0]["source_closure"][index] = replacement
            with self.assertRaises(ValueError):
                validate_evidence.validate_source_closure_paths(
                    substituted["builds"][0]["source_closure"])

        for build_index, build in enumerate(peer["builds"]):
            retained_directory = (
                pathlib.PurePosixPath(build["configure_stdout"]["path"]).parent /
                f"{build['name']}-dependencies").as_posix()
            for bad_jobs in (0, 9, "4", True):
                candidate = json.loads(json.dumps(peer))
                candidate["builds"][build_index]["jobs"] = bad_jobs
                with self.assertRaises(ValueError):
                    validate_evidence.validate_manifest(
                        candidate, source_root=peer_root)
            dry_run = json.loads(json.dumps(peer))
            dry_run["builds"][build_index]["build_argv"][-1] = "--dry-run"
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(dry_run, source_root=peer_root)

            for dependency_index, dependency in enumerate(
                    build["dependency_files"]):
                relabeled = json.loads(json.dumps(build["dependency_files"]))
                relabeled[dependency_index]["build_path"] = (
                    f".research/leopard2/forged-dependencies/{build['name']}/"
                    f"{pathlib.PurePosixPath(dependency['build_path']).name}")
                with self.assertRaises(ValueError):
                    validate_evidence.dependency_source_closure(
                        relabeled, build["build_dir"], peer_root, peer_root,
                        retained_directory, live=False)

            swapped = json.loads(json.dumps(build["dependency_files"]))
            swapped[0]["retained"], swapped[1]["retained"] = (
                swapped[1]["retained"], swapped[0]["retained"])
            with self.assertRaises(ValueError):
                validate_evidence.dependency_source_closure(
                    swapped, build["build_dir"], peer_root, peer_root,
                    retained_directory, live=False)

            other_build = peer["builds"][
                (build_index + 1) % len(peer["builds"])]
            cross_build = json.loads(json.dumps(build["dependency_files"]))
            cross_build[0]["retained"] = json.loads(json.dumps(
                other_build["dependency_files"][0]["retained"]))
            with self.assertRaises(ValueError):
                validate_evidence.dependency_source_closure(
                    cross_build, build["build_dir"], peer_root, peer_root,
                    retained_directory, live=False)

        dependency = peer["builds"][0]["dependency_files"][0]
        retained_path = peer_root / dependency["retained"]["path"]
        original_retained = retained_path.read_bytes()
        forged_target = (
            b"CMakeFiles/forged.dir/forged.cpp.o:" +
            original_retained.split(b":", 1)[1])
        relabeled_target = json.loads(json.dumps(
            peer["builds"][0]["dependency_files"]))
        first_build = peer["builds"][0]
        first_retained_directory = (
            pathlib.PurePosixPath(
                first_build["configure_stdout"]["path"]).parent /
            f"{first_build['name']}-dependencies").as_posix()
        relabeled_record = relabeled_target[0]["retained"]
        try:
            retained_path.write_bytes(forged_target)
            relabeled_record["bytes"] = len(forged_target)
            relabeled_record["sha256"] = hashlib.sha256(
                forged_target).hexdigest()
            with self.assertRaises(ValueError):
                validate_evidence.dependency_source_closure(
                    relabeled_target, first_build["build_dir"], peer_root,
                    peer_root, first_retained_directory, live=False)
        finally:
            retained_path.write_bytes(original_retained)

        forged_retained = (
            f"object: {run_matrix.NORMALIZATION_TOKEN}/leopard.cpp\n".encode())
        candidate = json.loads(json.dumps(peer))
        candidate_record = candidate["builds"][0]["dependency_files"][0][
            "retained"]
        try:
            retained_path.write_bytes(forged_retained)
            candidate_record["bytes"] = len(forged_retained)
            candidate_record["sha256"] = hashlib.sha256(
                forged_retained).hexdigest()
            candidate_record["source_root_tokens"] = 1
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(
                    candidate, source_root=peer_root)
        finally:
            retained_path.write_bytes(original_retained)

        raw_path = peer_root / dependency["build_path"]
        original_raw = raw_path.read_bytes()
        try:
            raw_path.write_bytes(b"object: /tmp/substituted.cpp\n")
            with self.assertRaises(ValueError):
                validate_evidence.validate_manifest(
                    peer, source_root=peer_root, live=True)
        finally:
            raw_path.write_bytes(original_raw)


if __name__ == "__main__":
    unittest.main(verbosity=2)
