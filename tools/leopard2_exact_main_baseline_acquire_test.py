#!/usr/bin/env python3
"""Deterministic tests for exact-main acquisition and lane sealing."""

from __future__ import annotations

import ast
import contextlib
import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
import leopard2_exact_main_baseline as identity_contract  # noqa: E402
import leopard2_exact_main_baseline_acquire as acquire  # noqa: E402
import leopard2_exact_main_baseline_record as record_contract  # noqa: E402
import leopard2_exact_main_baseline_verifier as verifier  # noqa: E402


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


verifier_fixtures = load_module(
    "leopard2_exact_main_baseline_acquire_verifier_fixtures",
    TOOLS / "leopard2_exact_main_baseline_verifier_test.py",
)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def resign_for_root(record: dict, root: Path) -> dict:
    value = copy.deepcopy(record)
    value["lane"]["root"] = str(root)
    value["record_sha256"] = identity_contract.canonical_sha256({
        key: copy.deepcopy(field)
        for key, field in value.items()
        if key != "record_sha256"
    })
    if value["schema"] == record_contract.AUTHORITY_SCHEMA:
        return record_contract.validate_baseline_authority_record(value)
    return record_contract.validate_baseline_failure_record(value)


@contextlib.contextmanager
def fixture_anchors():
    with contextlib.ExitStack() as stack:
        for module in (record_contract, verifier):
            stack.enter_context(mock.patch.object(
                module, "BASELINE_COMMIT",
                verifier_fixtures._TEST_BASELINE_COMMIT))
            stack.enter_context(mock.patch.object(
                module, "BASELINE_TREE",
                verifier_fixtures._TEST_BASELINE_TREE))
        yield


def fixture_for_root(root: Path, kind: str) -> tuple[dict[str, bytes], dict]:
    with fixture_anchors():
        if kind == "authority":
            files, record = verifier_fixtures.authority_files_and_record()
            terminal = "baseline-authority.json"
        else:
            files, record = verifier_fixtures.failure_files_and_record(
                verification_failure=kind == "verification_failure")
            terminal = "FAILED.json"
        relocated = resign_for_root(record, root)
    retained = copy.deepcopy(files)
    retained.pop(terminal)
    return retained, relocated


def lane_plan(root: str = "/tmp/leopard-exact-main-lane-a1") -> acquire.LanePlan:
    return acquire.LanePlan(
        lane_root=root,
        attempt=1,
        repository="/home/catid/leopard",
        verifier=(
            "/home/catid/leopard/tools/"
            "leopard2_exact_main_baseline_verifier.py"),
        canonical_adapter_root="/tmp/leopard-acquire-canonical-adapter",
        canonical_baseline_root="/tmp/leopard-acquire-canonical-baseline",
        canonical_build_root="/tmp/leopard-acquire-canonical-build",
        variant_adapter_root="/tmp/leopard-acquire-variant-adapter",
        variant_baseline_root="/tmp/leopard-acquire-variant-baseline",
        variant_build_root="/tmp/leopard-acquire-variant-build",
    )


class ExactMainBaselineAcquireTest(unittest.TestCase):
    def assertAcquireError(self, function, *arguments, **keywords) -> None:  # noqa: N802
        with self.assertRaises(acquire.AcquisitionError):
            function(*arguments, **keywords)

    def test_lane_plan_is_detached_and_roots_are_disjoint(self) -> None:
        plan = lane_plan()
        self.assertEqual(acquire.validate_lane_plan(plan), plan)
        for mutation in (
                {"attempt": 0}, {"attempt": 4},
                {"lane_root": "relative/lane"},
                {"variant_build_root": plan.canonical_build_root + "/child"},
                {"variant_build_root": "/var" + plan.canonical_build_root +
                 "-suffix"},
                {"canonical_adapter_root": plan.canonical_baseline_root},
                {"lane_root": plan.repository + "/lane"},
                {"canonical_build_root": plan.repository + "/build"},
                {"verifier": plan.canonical_build_root + "/verifier.py"},
                {"verifier": plan.repository}):
            with self.subTest(mutation=mutation):
                value = acquire.LanePlan(**{
                    **plan.__dict__, **mutation,
                })
                self.assertAcquireError(acquire.validate_lane_plan, value)
        self.assertAcquireError(acquire.validate_lane_plan, plan.__dict__)

    def test_canonical_ldd_and_build_closure_helpers(self) -> None:
        rows = [
            {"soname": "libc.so.6", "kind": "file",
             "path": "/usr/lib/libc.so.6"},
            {"soname": "linux-vdso.so.1", "kind": "virtual",
             "path": None},
        ]
        content = acquire.canonical_ldd_text(rows)
        self.assertEqual(
            list(record_contract.parse_canonical_ldd_output(content)), rows)
        self.assertAcquireError(acquire.canonical_ldd_text,
                                list(reversed(rows)))
        self.assertAcquireError(acquire.canonical_ldd_text, [rows[0], rows[0]])
        self.assertAcquireError(
            acquire.canonical_ldd_text,
            [rows[0]] * (record_contract.MAX_DEPENDENCIES + 1))
        oversized_soname = copy.deepcopy(rows[0])
        oversized_soname["soname"] = \
            "a" * (record_contract.MAX_TEXT_LENGTH + 1)
        self.assertAcquireError(
            acquire.canonical_ldd_text, [oversized_soname])
        files = [
            {"relative_path": "leopard_main_benchmark", "size": 5,
             "sha256": sha256(b"main\n")},
            {"relative_path": "libleopard_main_exact.a", "size": 8,
             "sha256": sha256(b"!<arch>\n")},
        ]
        files.sort(key=lambda item: item["relative_path"])
        closure = acquire.build_closure_document(
            "canonical_first", "/tmp/exact-main-build", files)
        self.assertEqual(closure, {
            "schema": acquire.BUILD_CLOSURE_SCHEMA,
            "role": "canonical_first",
            "build_root": "/tmp/exact-main-build",
            "files": files,
            "file_count": 2,
        })
        self.assertAcquireError(
            acquire.build_closure_document, "canonical_first",
            "/tmp/exact-main-build", list(reversed(files)))

    def test_authority_and_failure_lanes_verify_offline(self) -> None:
        for kind, expected_status in (
                ("authority", 0),
                ("acquisition_failure", 3),
                ("verification_failure", 3)):
            with self.subTest(kind=kind), fixture_anchors(), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-exact-main-acquire-test.") as temporary:
                root = Path(temporary) / "lane"
                retained, record = fixture_for_root(root, kind)
                with acquire.LaneWriter(str(root)) as writer:
                    seal = writer.seal_record(record, retained)
                self.assertEqual(seal["terminal"],
                                 "baseline-authority.json" if kind ==
                                 "authority" else "FAILED.json")
                self.assertEqual(oct(os.stat(root).st_mode & 0o7777), "0o500")
                self.assertFalse(any(
                    path.name.startswith(".leopard-stage-")
                    for path in root.rglob("*")))
                for path in root.rglob("*"):
                    mode = os.lstat(path).st_mode & 0o7777
                    self.assertEqual(
                        mode, acquire.LANE_DIRECTORY_MODE if path.is_dir()
                        else acquire.LANE_FILE_MODE)
                completed = subprocess.run(
                    verifier_fixtures.verifier_cli_command(root),
                    check=False, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, timeout=30)
                self.assertEqual(
                    completed.returncode, expected_status,
                    completed.stderr.decode("utf-8", "replace"))

    def test_lane_root_is_never_reused_and_publish_is_no_replace(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-collision.") as temporary:
            root = Path(temporary) / "lane"
            with acquire.LaneWriter(str(root)) as writer:
                collision = root / "collision.bin"
                collision.write_bytes(b"preexisting\n")
                self.assertAcquireError(
                    writer.publish_bytes, "collision.bin", b"replacement\n")
                self.assertEqual(collision.read_bytes(), b"preexisting\n")
                self.assertFalse(any(
                    path.name.startswith(".leopard-stage-")
                    for path in root.iterdir()))
            self.assertAcquireError(acquire.LaneWriter, str(root))

    def test_seal_rejects_inventory_and_byte_drift(self) -> None:
        for mutation in ("missing", "extra", "changed"):
            with self.subTest(mutation=mutation), fixture_anchors(), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-exact-main-acquire-inventory.") \
                    as temporary:
                root = Path(temporary) / "lane"
                retained, record = fixture_for_root(
                    root, "acquisition_failure")
                if mutation == "missing":
                    retained.pop(sorted(retained)[0])
                elif mutation == "extra":
                    retained["diagnostics/extra"] = b"extra\n"
                else:
                    retained[sorted(retained)[0]] += b"changed\n"
                with acquire.LaneWriter(str(root)) as writer:
                    self.assertAcquireError(
                        writer.seal_record, record, retained)

    def test_terminal_root_and_record_shape_are_bound(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-record.") as temporary:
            root = Path(temporary) / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            wrong_root = copy.deepcopy(record)
            wrong_root["lane"]["root"] = str(root) + "-other"
            wrong_root["record_sha256"] = identity_contract.canonical_sha256({
                key: copy.deepcopy(value)
                for key, value in wrong_root.items()
                if key != "record_sha256"
            })
            with acquire.LaneWriter(str(root)) as writer:
                self.assertAcquireError(
                    writer.seal_record, wrong_root, retained)

    def test_final_seal_rebinds_every_published_byte(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-substitution.") \
                as temporary:
            root = Path(temporary) / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            with acquire.LaneWriter(str(root)) as writer:
                apply_modes = writer._apply_final_modes

                def substitute_then_lock(files, directories):
                    terminal = root / "FAILED.json"
                    content = terminal.read_bytes()
                    terminal.write_bytes(b"X" + content[1:])
                    apply_modes(files, directories)

                with mock.patch.object(
                        writer, "_apply_final_modes",
                        side_effect=substitute_then_lock):
                    self.assertAcquireError(
                        writer.seal_record, record, retained)

    def test_production_boundary_is_independent_and_optimized_safe(self) -> None:
        source = (TOOLS / "leopard2_exact_main_baseline_acquire.py").read_text(
            encoding="utf-8")
        tree = ast.parse(source)
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(
                    alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.add((node.module or "").split(".")[0])
        self.assertNotIn("subprocess", imported)
        self.assertNotIn(
            "leopard2_exact_main_baseline_verifier", imported)
        contract_calls = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.Call) and
            isinstance(node.func, ast.Name) and
            node.func.id == "_load_local_contract"
        ]
        self.assertEqual(len(contract_calls), 2)
        loaded_contracts = set()
        for call in contract_calls:
            self.assertEqual(call.keywords, [])
            self.assertEqual(len(call.args), 2)
            self.assertTrue(all(
                isinstance(argument, ast.Constant) and
                type(argument.value) is str
                for argument in call.args))
            loaded_contracts.add(call.args[1].value)
        self.assertEqual(loaded_contracts, {
            "leopard2_exact_main_baseline.py",
            "leopard2_exact_main_baseline_record.py",
        })
        self.assertFalse(any(isinstance(node, ast.Assert)
                             for node in ast.walk(tree)))


if __name__ == "__main__":
    unittest.main()
