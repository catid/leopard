#!/usr/bin/env python3
"""Unit tests for C7 v5 evidence portability primitives."""

from __future__ import annotations

import copy
import hashlib
import os
import pathlib
import subprocess
import tempfile
import unittest
from unittest import mock

import run_matrix
import validate_evidence


class PortabilityTests(unittest.TestCase):
    def test_recursive_typed_equality_rejects_numeric_coercion(self) -> None:
        self.assertTrue(run_matrix.typed_equal(
            {"a": [0, 1, 1.0, {"flag": True}]},
            {"a": [0, 1, 1.0, {"flag": True}]}))
        for left, right in (
            (False, 0), (True, 1), (1, 1.0),
            ({"count": False}, {"count": 0}),
            ([{"field": True}], [{"field": 1}]),
            ({False: "value"}, {0: "value"}),
        ):
            self.assertFalse(run_matrix.typed_equal(left, right))

    def test_dual_git_identity_preflight_rejects_mutations(self) -> None:
        with tempfile.TemporaryDirectory(prefix="c7-dual-git-") as directory:
            root = pathlib.Path(directory)

            def git(*arguments: str) -> str:
                return subprocess.run(
                    ["git", *arguments], cwd=root, check=True, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                ).stdout.strip()

            git("init", "-q")
            git("config", "user.name", "C7 test")
            git("config", "user.email", "c7-test@example.invalid")
            (root / "core.txt").write_text("core-v1\n", encoding="utf-8")
            (root / "tool.txt").write_text("tool-v1\n", encoding="utf-8")
            git("add", "core.txt", "tool.txt")
            git("commit", "-q", "-m", "core")
            core_sha = git("rev-parse", "HEAD")
            (root / "tool.txt").write_text("tool-v2\n", encoding="utf-8")
            git("add", "tool.txt")
            git("commit", "-q", "-m", "tooling")
            tooling_sha = git("rev-parse", "HEAD")

            run_matrix.require_generation_identities(
                tooling_sha, core_sha, source_root=root,
                core_closure=("core.txt",), tooling_closure=("tool.txt",))
            validate_evidence.validate_repository(
                root, tooling_sha, core_sha, require_checkout_head=True)

            git("checkout", "-q", "--detach", core_sha)
            (root / "side.txt").write_text("side\n", encoding="utf-8")
            git("add", "side.txt")
            git("commit", "-q", "-m", "unrelated core candidate")
            nonancestor_sha = git("rev-parse", "HEAD")
            git("checkout", "-q", "--detach", tooling_sha)
            with self.assertRaisesRegex(RuntimeError, "not an ancestor"):
                run_matrix.require_generation_identities(
                    tooling_sha, nonancestor_sha, source_root=root,
                    core_closure=("core.txt",),
                    tooling_closure=("tool.txt",))
            with self.assertRaisesRegex(ValueError, "not an ancestor"):
                validate_evidence.validate_repository(
                    root, tooling_sha, nonancestor_sha,
                    require_checkout_head=False)

            (root / "untracked.txt").write_text("dirty\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "must be clean"):
                run_matrix.require_generation_identities(
                    tooling_sha, core_sha, source_root=root,
                    core_closure=("core.txt",),
                    tooling_closure=("tool.txt",))
            with self.assertRaisesRegex(ValueError, "not clean"):
                validate_evidence.validate_repository(
                    root, tooling_sha, core_sha, require_checkout_head=True)
            (root / "untracked.txt").unlink()

            for label, commit, filename, changed in (
                ("core", core_sha, "core.txt", "changed-core\n"),
                ("tooling", tooling_sha, "tool.txt", "changed-tool\n"),
            ):
                path = root / filename
                original = path.read_bytes()
                try:
                    path.write_text(changed, encoding="utf-8")
                    with self.assertRaisesRegex(
                            RuntimeError, f"{label} closure differs"):
                        run_matrix.require_committed_closure(
                            label, commit, (filename,), root)
                finally:
                    path.write_bytes(original)

    def test_historical_v4_source_replay_uses_recorded_git_bytes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="c7-v4-source-replay-") as directory:
            root = pathlib.Path(directory)

            def git(*arguments: str) -> str:
                return subprocess.run(
                    ["git", *arguments], cwd=root, check=True, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                ).stdout.strip()

            git("init", "-q")
            git("config", "user.name", "C7 test")
            git("config", "user.email", "c7-test@example.invalid")
            source = root / "leopard2.cpp"
            historical_bytes = b"historical-v4-core\n"
            source.write_bytes(historical_bytes)
            git("add", source.name)
            git("commit", "-q", "-m", "historical core")
            core_sha = git("rev-parse", "HEAD")
            record = {
                "bytes": len(historical_bytes), "path": source.name,
                "sha256": hashlib.sha256(historical_bytes).hexdigest(),
            }
            source.write_bytes(b"current-v5-checkout-bytes\n")
            git("add", source.name)
            git("commit", "-q", "-m", "current tooling checkout")

            validate_evidence.validate_committed_source_artifact(
                record, "historical v4 source", root, core_sha,
                current=False)
            with self.assertRaisesRegex(ValueError, "disagree with the manifest"):
                validate_evidence.validate_committed_source_artifact(
                    record, "current v5 source", root, core_sha,
                    current=True)

    def test_current_recipe_bytes_and_exact_semantics_are_bound(self) -> None:
        identity = validate_evidence.CANONICAL_CMAKE_IDENTITY
        build_dir = "build/core-avx2"
        records = {"ar": pathlib.Path("/usr/bin/ar"),
                   "ranlib": pathlib.Path("/usr/bin/ranlib")}
        objects = [
            f"CMakeFiles/{(
                'leopard2_backend_avx2.dir' if member ==
                'Leopard2BackendAVX2.cpp.o' else
                'leopard2_backend_ssse3.dir' if member ==
                'Leopard2BackendSSSE3.cpp.o' else
                identity['target_directory'])}/{member}"
            for member in validate_evidence.CURRENT_ARCHIVE_MEMBERS
        ]
        canonical = (
            f"/usr/bin/ar qc {identity['archive']} {' '.join(objects)}\n"
            f"/usr/bin/ranlib {identity['archive']}\n")

        def fixture(text: str) -> dict:
            data = text.encode("utf-8")
            digest = hashlib.sha256(data).hexdigest()
            return {
                "archive_link_recipe": {
                    "bytes": len(data),
                    "path": f"{build_dir}/CMakeFiles/leopard.dir/link.txt",
                    "sha256": digest,
                },
                "archive_link_recipe_content": {
                    "bytes": len(data), "encoding": "utf-8",
                    "sha256": digest, "text": text,
                },
            }

        with tempfile.TemporaryDirectory(prefix="c7-recipe-") as directory:
            root = pathlib.Path(directory)
            validate_evidence.validate_archive_link_recipe(
                fixture(canonical), "avx2", build_dir, identity, records,
                root, live=False)
            mutations = (
                canonical.replace("qc libleopard.a", "qc nested/libleopard.a", 1),
                canonical.replace("/usr/bin/ar", "/tmp/ar", 1),
                canonical.replace("LeopardCommon.cpp.o ",
                                  "LeopardCommon.cpp.o @objects.rsp ", 1),
                canonical.replace(
                    f"{objects[0]} {objects[1]}",
                    f"{objects[1]} {objects[0]}", 1),
                canonical.replace("/usr/bin/ranlib", "/tmp/ranlib", 1),
            )
            for text in mutations:
                with self.subTest(text=text[:48]), self.assertRaises(ValueError):
                    validate_evidence.validate_archive_link_recipe(
                        fixture(text), "avx2", build_dir, identity, records,
                        root, live=False)
            stale = fixture(canonical)
            stale["archive_link_recipe_content"]["text"] += "\n"
            with self.assertRaises(ValueError):
                validate_evidence.validate_archive_link_recipe(
                    stale, "avx2", build_dir, identity, records,
                    root, live=False)

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

    def test_compiler_launcher_rewrites_only_exact_checkout_paths(self) -> None:
        cwd = run_matrix.ROOT / ".research/leopard2/c7-launch/core-scalar"
        root = str(run_matrix.ROOT)
        relative_root = os.path.relpath(run_matrix.ROOT, cwd)
        argv = [
            "/usr/bin/c++", f"-I{root}", f"{root}/leopard.cpp",
            f"-ffile-prefix-map={root}=LEO2_SOURCE_ROOT",
            f"{root}-sibling/leopard.cpp",
        ]
        self.assertEqual(
            run_matrix.portable_compiler_argv(argv, cwd),
            [
                "/usr/bin/c++", f"-I{relative_root}",
                os.path.relpath(run_matrix.ROOT / "leopard.cpp", cwd),
                f"-ffile-prefix-map={root}=LEO2_SOURCE_ROOT",
                f"{root}-sibling/leopard.cpp",
            ])

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
            path = pathlib.Path(directory) / "libleopard.a"
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

    def test_all_artifact_paths_are_checkout_contained(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="c7-contained-", dir=run_matrix.ROOT) as directory:
            with tempfile.TemporaryDirectory(prefix="c7-outside-") as outside:
                directory_path = pathlib.Path(directory)
                outside_file = pathlib.Path(outside) / "artifact"
                outside_file.write_bytes(b"outside")
                record = {
                    "path": str(outside_file), "bytes": 7,
                    "sha256": hashlib.sha256(b"outside").hexdigest(),
                }
                with self.assertRaises(ValueError):
                    validate_evidence.validate_artifact(
                        record, "absolute", run_matrix.ROOT, required=True)
                record["path"] = "../escape"
                with self.assertRaises(ValueError):
                    validate_evidence.validate_artifact(
                        record, "parent", run_matrix.ROOT, required=False)
                for noncanonical in ("a//artifact", "C:/artifact", "./artifact"):
                    record["path"] = noncanonical
                    with self.assertRaises(ValueError):
                        validate_evidence.validate_artifact(
                            record, "noncanonical", run_matrix.ROOT,
                            required=False)
                link = directory_path / "escape"
                link.symlink_to(pathlib.Path(outside), target_is_directory=True)
                record["path"] = (
                    directory_path.relative_to(run_matrix.ROOT) /
                    "escape/artifact"
                ).as_posix()
                with self.assertRaises(ValueError):
                    validate_evidence.validate_artifact(
                        record, "symlink", run_matrix.ROOT, required=True)

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
        self.assertIsNone(validate_evidence.ABSOLUTE_PROJECT_PATH.search(
            "/usr/bin/python3;${LEO2_SOURCE_ROOT}/experiments/leopard2/"))
        self.assertIsNotNone(validate_evidence.ABSOLUTE_PROJECT_PATH.search(
            "/tmp/foreign-checkout/experiments/leopard2/c7/file.cpp"))
        for header in ("leopard.h", "leopard2.h"):
            self.assertIsNotNone(
                validate_evidence.ABSOLUTE_PROJECT_PATH.search(
                    f"/tmp/foreign-checkout/{header}"))

    def test_reproducibility_fingerprint_mutation_rejected(self) -> None:
        artifacts = {
            name: {
                "library_sha256": str(index) * 64,
                "executable_sha256": format(index + 8, "x") * 64,
            }
            for index, name in enumerate(run_matrix.BUILD_NAMES, start=1)
        }
        base = {
            "schema": run_matrix.MANIFEST_SCHEMA,
            "tooling_git_sha": "e" * 40,
            "core_git_sha": "a" * 40,
            "source": {"path": "s", "sha256": "b" * 64, "bytes": 1},
            "runner": {"path": "r", "sha256": "c" * 64, "bytes": 1},
            "validator": {"path": "v", "sha256": "d" * 64, "bytes": 1},
            "reproducibility": {"fingerprints": artifacts},
            "taskset": {"path": "/usr/bin/taskset"},
            "builds": [
                {
                    "name": name,
                    **{role: {"path": f"/usr/bin/{role}"}
                       for role in run_matrix.PROGRAM_ROLES},
                }
                for name in run_matrix.BUILD_NAMES
            ],
        }
        run_matrix.require_reproducible_peer(base, copy.deepcopy(base))
        forged = copy.deepcopy(base)
        forged["reproducibility"]["fingerprints"]["scalar"][
            "library_sha256"] = "f" * 64
        with self.assertRaises(RuntimeError):
            run_matrix.require_reproducible_peer(base, forged)
        forged = copy.deepcopy(base)
        forged["source"]["bytes"] = True
        with self.assertRaisesRegex(RuntimeError, "committed tooling"):
            run_matrix.require_reproducible_peer(base, forged)
        for key in ("tooling_git_sha", "core_git_sha"):
            forged = copy.deepcopy(base)
            forged[key] = "f" * 40
            with self.assertRaisesRegex(RuntimeError, key):
                run_matrix.require_reproducible_peer(base, forged)

    def test_comparison_attestation_is_machine_checkable(self) -> None:
        fingerprints = {
            name: {
                "library_sha256": format(index, "x") * 64,
                "executable_sha256": format(index + 8, "x") * 64,
            }
            for index, name in enumerate(run_matrix.BUILD_NAMES, start=1)
        }
        scan = {
            "normalized_text_records": 42,
            "archives": len(run_matrix.BUILD_NAMES),
            "executables": len(run_matrix.BUILD_NAMES),
        }
        comparison = {
            "status": "pass",
            "peer_manifest_sha256": "a" * 64,
            "fingerprints_sha256":
                run_matrix.canonical_json_sha256(fingerprints),
            "build_names": list(run_matrix.BUILD_NAMES),
            "checkout_roots_scanned": 2,
            "current_scan": scan,
            "peer_scan": scan,
            "peer_manifest": {
                "path": "results/peer-manifest.json",
                "bytes": 1, "sha256": "a" * 64,
            },
            "peer_evidence_bundle": {
                "path": "results/peer-evidence-bundle.tar.gz",
                "bytes": 1, "sha256": "c" * 64,
            },
            "peer_attestation": {
                "path": "results/peer-reproducibility-attestation.json",
                "bytes": 1, "sha256": "b" * 64,
            },
        }
        self.assertEqual(comparison["fingerprints_sha256"],
                         validate_evidence.canonical_json_sha256(fingerprints))
        validate_evidence.validate_comparison(comparison, fingerprints, 42)
        for key, value in (
            ("status", "not-pass"),
            ("peer_manifest_sha256", "0"),
            ("fingerprints_sha256", "f" * 64),
            ("checkout_roots_scanned", 1),
        ):
            forged = copy.deepcopy(comparison)
            forged[key] = value
            with self.assertRaises(ValueError):
                validate_evidence.validate_comparison(forged, fingerprints, 42)
        wrong_scan = copy.deepcopy(comparison)
        wrong_scan["peer_scan"]["normalized_text_records"] = 41
        with self.assertRaises(ValueError):
            validate_evidence.validate_comparison(
                wrong_scan, fingerprints, 42)
        boolean_one = copy.deepcopy(comparison)
        boolean_one["current_scan"]["normalized_text_records"] = True
        boolean_one["peer_scan"]["normalized_text_records"] = True
        with self.assertRaises(ValueError):
            validate_evidence.validate_comparison(
                boolean_one, fingerprints, 1)
        with self.assertRaises(ValueError):
            validate_evidence.validate_comparison(
                {"status": "not-run"}, fingerprints, False)

    def test_canonical_peer_bundle_round_trip(self) -> None:
        files = {
            "evidence/empty.txt": b"",
            "evidence/result.json": b'{"status":"pass"}\n',
        }
        records = {
            path: {
                "path": path, "bytes": len(contents),
                "sha256": hashlib.sha256(contents).hexdigest(),
                "source_root_tokens": 0,
            }
            for path, contents in files.items()
        }
        bundle = run_matrix.canonical_peer_bundle(files)
        self.assertEqual(
            bundle[:10], b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff")
        self.assertEqual(
            validate_evidence._read_peer_bundle(bundle, records), files)
        empty_member = (
            b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff\x03\x00"
            b"\x00\x00\x00\x00\x00\x00\x00\x00")
        for forged in (
                bundle + bundle, empty_member + bundle, bundle + empty_member,
                bundle + b"\0", bundle + b"\0" * 512):
            with self.assertRaises(ValueError):
                validate_evidence._read_peer_bundle(forged, records)

    def test_strict_json_rejects_ambiguous_objects_and_constants(self) -> None:
        self.assertEqual(
            validate_evidence.strict_json_loads('{"a":1,"b":2}', "test"),
            {"a": 1, "b": 2})
        for document in ('{"a":1,"a":1}', '{"nested":{"x":0,"x":0}}',
                         '{"value":NaN}', '{"value":Infinity}',
                         '{"value":1e9999}'):
            with self.assertRaises(ValueError):
                validate_evidence.strict_json_loads(document, "test")

    def test_root_scan_uses_authenticated_single_read_snapshots(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="c7-root-scan-", dir=run_matrix.ROOT
        ) as directory:
            root = pathlib.Path(directory)

            def retain(name: str, contents: bytes, *, normalized: bool) -> dict:
                path = root / name
                path.write_bytes(contents)
                record = {
                    "path": path.relative_to(run_matrix.ROOT).as_posix(),
                    "bytes": len(contents),
                    "sha256": hashlib.sha256(contents).hexdigest(),
                }
                if normalized:
                    record["source_root_tokens"] = 0
                return record

            text_record = retain("retained.txt", b"safe\n", normalized=True)
            builds = []
            for name in run_matrix.BUILD_NAMES:
                builds.append({
                    "name": name,
                    "library": retain(
                        f"{name}.a", f"archive-{name}\n".encode(),
                        normalized=False),
                    "executable": retain(
                        f"{name}.exe", f"executable-{name}\n".encode(),
                        normalized=False),
                })
            manifest = {"retained": text_record, "builds": builds}
            forbidden = root / "forbidden-checkout"
            forged = f"{forbidden}/leopard2.cpp\n".encode()
            real_read = run_matrix._read_file_once
            for record in (text_record, builds[0]["library"]):
                target = run_matrix.ROOT / record["path"]
                original = target.read_bytes()
                swapped = False

                def read_then_swap(
                    candidate: pathlib.Path, maximum_bytes: int,
                ) -> bytes:
                    nonlocal swapped
                    contents = real_read(candidate, maximum_bytes)
                    if candidate.resolve() == target.resolve() and not swapped:
                        target.write_bytes(forged)
                        swapped = True
                    return contents

                try:
                    with mock.patch.object(
                            run_matrix, "_read_file_once",
                            side_effect=read_then_swap):
                        run_matrix.require_no_root_bytes(
                            manifest, run_matrix.ROOT, (forbidden,))
                    self.assertTrue(swapped)
                    self.assertEqual(target.read_bytes(), forged)
                finally:
                    target.write_bytes(original)


if __name__ == "__main__":
    unittest.main(verbosity=2)
