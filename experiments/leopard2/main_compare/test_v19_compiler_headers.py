#!/usr/bin/python3
"""Real descriptor/header-view mutations, synthetic borrowed runtime only."""
from contextlib import ExitStack
import copy
import hashlib
import importlib.util
import mmap
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_header_dispatch", HERE / "v19_runtime_dispatch.py")
dispatch_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dispatch_module)
module = dispatch_module.header_module
FAILURES = (module.builder.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


class HeaderTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-header-test-")))
        self.parent = self.root / "new"
        self.parent.mkdir(mode=0o700)
        self.parentfd = os.open(self.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        self.stack.callback(os.close, self.parentfd)
        self.system = self.root / "sys"
        self.system.mkdir(mode=0o700)
        (self.system / "nested").mkdir()
        (self.system / "empty").mkdir()
        self.roots = tuple(str(self.system / path) for path in ("nested", "empty", ""))
        self.stack.enter_context(mock.patch.object(module, "CPP_INCLUDE_ROOTS", self.roots))
        self.first = self.system / "one.h"
        self.first.write_bytes(b"#define ONE 1\n")
        self.first.chmod(0o644)
        self.second = self.system / "nested/two.h"
        self.second.write_bytes(b'#include "../one.h"\n')
        self.second.chmod(0o644)
        self.pins = [self.pin(path) for path in (self.first, self.second)]
        self.live = True
        test = self
        class Phase:
            language = "c++"
            logical_driver = "/usr/bin/c++"
            parent, _parent_fd = test.parent, test.parentfd
            def validate_current(self): module.require(test.live, "fixture runtime lost")
        class Inventory:
            phase = Phase()
            def validate_current(self): self.phase.validate_current()
        self.inventory = Inventory()

    @staticmethod
    def pin(path):
        return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size": path.stat().st_size}

    def factory(self, path, **kwargs):
        return module.builder._StreamedTool(path, _trusted_owner=(os.getuid(), os.getgid()), **kwargs)

    def owner(self, pins=None):
        return module.CompilerHeaders(self.inventory, self.pins if pins is None else pins, self.roots, _file_factory=self.factory)

    def enter(self):
        owner = self.owner()
        owner.__enter__()
        # Failure tests must force-close without suppressing the tested error.
        self.stack.callback(self.cleanup_owner, owner)
        return owner

    @classmethod
    def cleanup_owner(cls, owner):
        # The failure has already been asserted. Close test guards without
        # re-reporting it, then restore fixture permissions for temp cleanup.
        owner._view_guard._close_without_verification()
        owner._source_guard._close_without_verification()
        owner._stack.close()
        cls.make_writable(owner.root)

    @staticmethod
    def make_writable(root):
        for directory, _dirs, _files in os.walk(root): os.chmod(directory, 0o700)

    def test_sealed_bytes_order_overlap_and_empty_root(self):
        owner = self.enter()
        record = owner.record()
        self.assertEqual(record["include_roots"], list(self.roots))
        self.assertEqual(record["header_bytes"], sum(row["size"] for row in self.pins))
        self.assertFalse(record["compiler_data_owned"])
        self.assertFalse(record["negative_search_closure_owned"])
        self.assertFalse(record["full_header_read_closure_owned"])
        self.assertEqual(len({id(file.guard) for file in owner._files.values()}), 1)
        for path, file in owner._files.items():
            self.assertEqual((owner.root / str(path).lstrip("/")).read_bytes(), path.read_bytes())
            with self.assertRaises(OSError): os.pwrite(file.executable_descriptor, b"x", 0)
        self.assertEqual(list((owner.root / self.roots[1].lstrip("/")).iterdir()), [])

    def test_arguments_preserve_logical_input_and_map_only_view_prefix(self):
        owner = self.enter()
        argv = ["/usr/bin/c++", "-g", "-O3", "-o", "/tmp/out.o", "-c", "/tmp/source.cpp"]
        original = list(argv)
        result = owner.arguments(argv)
        self.assertEqual(argv, original)
        self.assertEqual(result[-len(argv)+1:], argv[1:])
        self.assertIn(f"-ffile-prefix-map={owner.root}=", result)
        self.assertNotIn("-fno-canonical-system-headers", result)
        self.assertEqual(result.count("-isystem"), 3)
        self.assertEqual(result.count("-nostdinc"), 1)

    def test_recipe_overrides_latch_failure(self):
        for value in ("-isystem", "-include", "--sysroot=x", "-nostdinc++", "-Wp,-include,x", "-Xpreprocessor",
                      "-fno-canonical-system-headers", "-E"):
            owner = self.enter()
            with self.subTest(value=value), self.assertRaises(FAILURES):
                owner.arguments(["/usr/bin/c++", "-o", "out", "-c", "source.cpp", value])
            with self.assertRaises(FAILURES): owner.record()

    def test_noncompile_or_oversized_argv_rejected(self):
        for argv in (["/usr/bin/c++"], ["/usr/bin/cc", "-o", "out", "-c", "s.c"],
                     ["/usr/bin/c++", "-o", "out", "-c", "x\0.cpp"],
                     ["/usr/bin/c++", "-o", "out", "-c", "s.cpp", *(["-g"] * 480)]):
            owner = self.enter()
            with self.assertRaises(FAILURES): owner.arguments(argv)

    def test_invalid_pin_shapes_and_bounds(self):
        bad = [[], self.pins * 257, self.pins + [self.pins[0]]]
        for key, value in (("size", True), ("size", 0), ("size", module.MAX_FILE_BYTES + 1),
                           ("sha256", "bad"), ("path", "/outside.h"), ("path", str(self.system)+"/x/../one.h")):
            row = dict(self.pins[0]); row[key] = value; bad.append([row])
        for rows in bad:
            with self.assertRaises(FAILURES): self.owner(rows)
        with self.assertRaises(FAILURES): module.header_pins(self.pins, self.roots[::-1])
        with self.assertRaises(FAILURES): module.header_pins(self.pins, list(self.roots))
        with mock.patch.object(module, "MAX_TOTAL_BYTES", 1), self.assertRaises(FAILURES): self.owner()

    def test_changed_content_pin_rejects_entry(self):
        rows = copy.deepcopy(self.pins)
        rows[0]["sha256"] = "0" * 64
        with self.assertRaises(FAILURES): self.owner(rows).__enter__()

    def test_symlink_and_hardlink_aliases_rejected(self):
        alias = self.system / "alias.h"
        alias.symlink_to(self.first)
        with self.assertRaises(FAILURES): self.owner(self.pins + [self.pin(alias)]).__enter__()
        alias.unlink()
        os.link(self.first, alias)
        with self.assertRaises(FAILURES): self.owner(self.pins + [self.pin(alias)]).__enter__()

    def test_source_write_restore_is_latched(self):
        owner = self.enter()
        original = self.first.read_bytes()
        self.first.write_bytes(b"x" * len(original)); self.first.write_bytes(original)
        with self.assertRaises(FAILURES): owner.validate_current()
        with self.assertRaises(FAILURES): owner.record()

    def test_prefaulted_mmap_is_rehashed(self):
        with self.first.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[0] ^= 1
            mapping[0] ^= 1
            owner = self.enter()
            file = owner._files[self.first]
            with mock.patch.object(file, "_hash", wraps=file._hash) as hashed:
                mapping[0] ^= 1
                try:
                    with self.assertRaises(FAILURES): owner.validate_current()
                    self.assertTrue(hashed.called)
                finally: mapping[0] ^= 1

    def test_view_symlink_replacement_restore_is_latched(self):
        owner = self.enter()
        leaf = owner.root / str(self.first).lstrip("/")
        target = os.readlink(leaf)
        leaf.parent.chmod(0o700)
        leaf.unlink(); leaf.symlink_to("/proc/self/fd/0")
        leaf.unlink(); leaf.symlink_to(target)
        leaf.parent.chmod(0o500)
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_source_parent_permission_restore_is_latched(self):
        owner = self.enter()
        mode = self.system.stat().st_mode & 0o7777
        self.system.chmod(0o500); self.system.chmod(mode)
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_directory_bound_rejects_before_view_creation(self):
        with mock.patch.object(module, "MAX_DIRECTORIES", 1), self.assertRaises(FAILURES): self.owner().__enter__()
        self.assertEqual(list(self.parent.iterdir()), [])

    def test_missing_runtime_and_inherited_directory_fail(self):
        owner = self.enter()
        self.live = False
        with self.assertRaises(FAILURES): owner.validate_current()
        self.live = True
        owner = self.enter()
        fd = owner._view_dirs[Path(".")][0]
        os.set_inheritable(fd, True)
        try:
            with self.assertRaises(FAILURES): owner.descriptors()
        finally: os.set_inheritable(fd, False)

    def test_detached_records_and_exit_close_every_descriptor(self):
        owner = self.enter()
        record = owner.record()
        record["files"][0]["sha256"] = "corrupted"
        self.assertNotEqual(record, owner.record())
        descriptors = [*owner.descriptors(), *(row[0] for row in owner._source_dirs.values()),
                       *(row[0] for row in owner._view_dirs.values()), *(file.fd for file in owner._files.values())]
        owner.__exit__(None, None, None)
        for fd in descriptors:
            with self.assertRaises(OSError): os.fstat(fd)
        with self.assertRaises(FAILURES): owner.record()

    def test_file_close_does_not_close_borrowed_guard(self):
        guard = self.stack.enter_context(module.provenance._InotifyMutationGuard("shared header fixture"))
        first = self.factory(self.first, _guard=guard, permitted_modes=(0o644, 0o755))
        second = self.factory(self.second, _guard=guard, permitted_modes=(0o644, 0o755))
        self.stack.callback(second.close)
        first.close()
        second.validate_current()
        guard.verify()

    def test_dispatch_rejects_wrong_owner_or_inventory(self):
        owner = self.enter()
        with self.assertRaises(FAILURES): dispatch_module.RuntimeDispatch(self.inventory, headers=object(), _runner=lambda *a: b"")
        other = copy.copy(self.inventory)
        with self.assertRaises(FAILURES): dispatch_module.RuntimeDispatch(other, headers=owner, _runner=lambda *a: b"")


if __name__ == "__main__": unittest.main()
