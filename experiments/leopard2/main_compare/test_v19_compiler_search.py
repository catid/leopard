#!/usr/bin/python3
"""Real absence/history/alias tests with a synthetic borrowed runtime."""
from contextlib import ExitStack
import importlib.util
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_search_dispatch", HERE / "v19_runtime_dispatch.py")
dispatch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dispatch)
module = dispatch.search_module
FAILURES = (module.builder.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


class SearchTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-search-test-")))
        self.system = self.root / "system"
        self.system.mkdir()
        self.missing = self.system / "specs"
        self.live = True
        test = self
        class Inventory:
            def validate_current(self): module.require(test.live, "fixture inventory lost")
        self.inventory = Inventory()

    def owner(self, paths=None):
        return module.CompilerSearch(self.inventory, _paths=(str(self.missing),) if paths is None else paths)

    def enter(self, paths=None):
        owner = self.owner(paths).__enter__()
        def cleanup():
            owner._guard._close_without_verification()
            owner._stack.close()
        self.stack.callback(cleanup)
        return owner

    def test_absence_record_is_detached_and_context_closes_descriptors(self):
        with self.owner() as owner:
            record = owner.record()
            self.assertEqual(record["missing"], [{"path": str(self.missing), "parent": str(self.system), "name": "specs"}])
            self.assertTrue(record["declared_absence_history_owned"])
            self.assertFalse(record["negative_search_closure_owned"])
            self.assertFalse(record["compiler_data_owned"])
            record["missing"][0]["name"] = "unrelated"
            self.assertEqual(owner.record()["missing"][0]["name"], "specs")
            descriptors = [fd for fd, _ in owner._directories.values()] + [owner._guard.descriptor]
        for fd in descriptors:
            with self.assertRaises(OSError): os.fstat(fd)
        with self.assertRaises(FAILURES): owner.__enter__()

    def test_persistent_file_creation_is_rejected(self):
        owner = self.enter()
        self.missing.write_bytes(b"*link:\n")
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_direct_absence_check_still_rejects_without_event_drain(self):
        owner = self.enter()
        self.missing.write_bytes(b"appeared")
        with mock.patch.object(owner._guard, "verify"), self.assertRaisesRegex(FAILURES, "input appeared"):
            owner.validate_current()

    def test_permission_error_is_not_treated_as_absence(self):
        original = os.stat
        def denied(path, *args, **kwargs):
            if path == "specs": raise PermissionError("fixture denied")
            return original(path, *args, **kwargs)
        before = len(os.listdir("/proc/self/fd"))
        with mock.patch.object(os, "stat", denied), self.assertRaises(PermissionError): self.owner().__enter__()
        self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_create_delete_history_latches_even_after_absence_is_restored(self):
        owner = self.enter()
        self.missing.write_bytes(b"*link:\n")
        self.missing.unlink()
        self.assertFalse(self.missing.exists())
        with self.assertRaises(FAILURES): owner.validate_current()
        owner._guard._close_without_verification()
        with self.assertRaises(FAILURES): owner.record()
        with self.assertRaises(FAILURES): owner.__exit__(None, None, None)

    def test_first_missing_component_is_held_before_cancelled_dotdot(self):
        (self.system / "present").write_bytes(b"exists")
        request = str(self.system) + "/missing/../present"
        owner = self.enter((request,))
        self.assertEqual(owner.record()["missing"][0]["name"], "missing")
        path = self.system / "missing"
        path.mkdir(); path.rmdir()
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_deep_missing_subtree_appearance_is_rejected(self):
        request = str(self.system / "absent/deep/specs")
        owner = self.enter((request,))
        path = self.system / "absent"
        path.mkdir(); path.rmdir()
        with self.assertRaises(FAILURES): owner.record()

    def test_relative_alias_dotdot_uses_resolved_directory(self):
        nested = self.system / "real/nested"
        nested.mkdir(parents=True)
        alias = self.system / "alias"
        alias.symlink_to("real/nested")
        request = str(alias) + "/../specs"
        owner = self.enter((request,))
        record = owner.record()
        self.assertEqual(record["missing"][0]["parent"], str(self.system / "real"))
        self.assertEqual(record["aliases"][str(alias)]["target"], "real/nested")
        (self.system / "specs").write_bytes(b"unrelated lexical path")
        owner.validate_current()

    def test_absolute_alias_is_retained_and_redirect_restore_is_rejected(self):
        alias = self.root / "alias"
        alias.symlink_to(self.system)
        owner = self.enter((str(alias / "specs"),))
        alias.unlink(); alias.symlink_to(self.root)
        alias.unlink(); alias.symlink_to(self.system)
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_parent_permission_restore_is_rejected(self):
        owner = self.enter()
        initial = self.system.stat().st_mode & 0o777
        self.system.chmod(0o700); self.system.chmod(initial)
        with self.assertRaises(FAILURES): owner.record()

    def test_parent_replacement_restore_is_rejected(self):
        owner = self.enter()
        moved = self.root / "old"
        self.system.rename(moved)
        self.system.mkdir(); self.system.rmdir()
        moved.rename(self.system)
        with self.assertRaises(FAILURES): owner.record()

    def test_unrelated_sibling_writes_are_not_a_failure(self):
        with self.owner() as owner:
            other = self.system / "unrelated"
            other.write_bytes(b"first"); other.write_bytes(b"second"); other.unlink()
            owner.validate_current()

    def test_existing_file_directory_and_nondirectory_traversal_rejected(self):
        self.missing.write_bytes(b"present")
        for request in (str(self.missing), str(self.system), str(self.missing / "child")):
            with self.subTest(request=request), self.assertRaises(FAILURES): self.owner((request,)).__enter__()

    def test_borrowed_inventory_loss_and_guard_loss_rejected(self):
        owner = self.enter()
        self.live = False
        with self.assertRaises(FAILURES): owner.record()
        self.live = True
        other = self.enter()
        other._guard._close_without_verification()
        with self.assertRaises(FAILURES): other.validate_current()

    def test_inheritable_directory_or_guard_rejected(self):
        for target in ("directory", "guard"):
            owner = self.enter()
            fd = owner._directories[self.system][0] if target == "directory" else owner._guard.descriptor
            os.set_inheritable(fd, True)
            with self.assertRaises(FAILURES): owner.validate_current()
            os.set_inheritable(fd, False)

    def test_changes_at_context_exit_rejected(self):
        owner = self.enter()
        self.missing.symlink_to("still-absent")
        with self.assertRaises(FAILURES): owner.__exit__(None, None, None)

    def test_invalid_paths_and_counts_rejected_before_capture(self):
        for paths in ([], (), ("relative",), ("/bad//path",), ("/bad\0path",), ("/",), (False,),
                      (str(self.missing),) * 2, tuple("/missing/" + str(index) for index in range(129)),
                      ("/" + "x" * 4096,), ("/" + "/".join("x" for _ in range(65)),)):
            with self.subTest(paths=paths), self.assertRaises(FAILURES): self.owner(paths)
        with self.assertRaises(FAILURES): module.CompilerSearch(self.inventory)

    def test_default_profile_uses_pinned_libexec_collect2_for_both_languages(self):
        tool = SimpleNamespace(path=Path("/usr/libexec/gcc/x86_64-linux-gnu/13/collect2"))
        self.inventory.phase = SimpleNamespace(language="c++", _tools={"collect2": tool})
        with mock.patch.object(module.runtime, "RuntimeInventory", type(self.inventory)):
            for language in ("c", "c++"):
                self.inventory.phase.language = language
                owner = module.CompilerSearch(self.inventory)
                self.assertEqual(len(owner.paths), 40)
                self.assertIn("/usr/lib/gcc/x86_64-linux-gnu/13/specs", owner.paths)
            tool.path = Path("/usr/lib/gcc/x86_64-linux-gnu/13/collect2")
            with self.assertRaises(FAILURES): module.CompilerSearch(self.inventory)

    def test_alias_loop_and_resource_bounds_fail_with_no_fd_leak(self):
        alias = self.system / "loop"
        alias.symlink_to("loop")
        before = len(os.listdir("/proc/self/fd"))
        with self.assertRaises(FAILURES): self.owner((str(alias / "specs"),)).__enter__()
        self.assertEqual(len(os.listdir("/proc/self/fd")), before)
        for constant in ("MAX_DIRECTORIES", "MAX_LINKS", "MAX_STEPS"):
            with mock.patch.object(module, constant, 0), self.assertRaises(FAILURES):
                self.owner((str(alias / "specs"),)).__enter__()
            self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_capture_race_after_parent_watch_is_rejected(self):
        add = module.provenance._InotifyMutationGuard._add_watch
        fired = False
        def race(guard, path, mask, names):
            nonlocal fired
            add(guard, path, mask, names)
            if path == self.system and names and b"specs" in names and not fired:
                fired = True
                self.missing.write_bytes(b"transient"); self.missing.unlink()
        with mock.patch.object(module.provenance._InotifyMutationGuard, "_add_watch", race), self.assertRaises(FAILURES):
            self.owner().__enter__()
        self.assertTrue(fired)


if __name__ == "__main__": unittest.main()
