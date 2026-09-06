#!/usr/bin/python3
"""Bounded link-data ownership tests; no codec execution."""
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
spec = importlib.util.spec_from_file_location("tested_link_dispatch", HERE / "v19_runtime_dispatch.py")
dispatch_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dispatch_module)
module = dispatch_module.link_module
FAILURES = (module.builder.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


class LinkerTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        # Keep disk writeback from changing the metadata of the rehash fixture.
        fixture_parent = "/dev/shm" if self._testMethodName == "test_preexisting_mmap_is_rehashed" else None
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-link-test-", dir=fixture_parent)))
        self.parent, self.gcc, self.system = (self.root / name for name in ("new", "gcc", "system"))
        for path in (self.parent, self.gcc, self.system): path.mkdir(mode=0o700)
        self.parentfd = os.open(self.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        self.stack.callback(os.close, self.parentfd)
        self.pins = []
        for original in module.LINK_INPUT_PATHS:
            path = (self.gcc if original.startswith(module.GCC_ROOT) else self.system) / Path(original).name
            path.write_bytes(b"fixture link input: " + path.name.encode() + b"\n")
            path.chmod(0o644)
            self.pins.append(self.pin(path))
        self.stack.enter_context(mock.patch.object(module, "LINK_INPUT_PATHS", tuple(row["path"] for row in self.pins)))
        self.stack.enter_context(mock.patch.object(module, "GCC_ROOT", str(self.gcc) + "/"))
        self.stack.enter_context(mock.patch.object(module, "SYSTEM_ROOT", str(self.system) + "/"))
        self.live = True
        test = self
        class Phase:
            language, logical_driver = "c++", "/usr/bin/c++"
            parent, _parent_fd = test.parent, test.parentfd
            def validate_current(self): module.require(test.live, "fixture runtime lost")
        class Inventory:
            phase = Phase()
            def validate_current(self): self.phase.validate_current()
        self.inventory = Inventory()

    @staticmethod
    def pin(path):
        return {"path": str(path), "size": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def factory(self, path, **kwargs):
        return module.builder._StreamedTool(path, _trusted_owner=(os.getuid(), os.getgid()), **kwargs)

    def owner(self, pins=None):
        return module.LinkerInputs(self.inventory, self.pins if pins is None else pins, _file_factory=self.factory)

    def enter(self):
        owner = self.owner().__enter__()
        def cleanup():
            owner._view_guard._close_without_verification()
            owner._source_guard._close_without_verification()
            owner._stack.close()
            for directory, _dirs, _files in os.walk(owner.root): os.chmod(directory, 0o700)
        self.stack.callback(cleanup)
        return owner

    def argv(self):
        return ["/usr/bin/c++", "-g", "-O0", "-O3", "adapter.o", "-o", "output", "library.a",
                str(self.gcc / "libgomp.so"), str(self.system / "libpthread.a")]

    def test_exact_aliases_seals_and_original_script_bytes(self):
        owner = self.enter()
        record = owner.record()
        self.assertEqual(len(record["files"]), 18)
        self.assertEqual(len(record["mappings"]), 53)
        self.assertTrue(record["original_script_bytes_preserved"])
        for name in ("compiler_data_owned", "full_link_read_closure_owned", "negative_search_closure_owned"):
            self.assertFalse(record[name])
        self.assertEqual(len({id(file.guard) for file in owner._files.values()}), 1)
        for relative, path in owner._aliases.items():
            self.assertEqual((owner.root / relative).read_bytes(), path.read_bytes())
        for file in owner._files.values():
            with self.assertRaises(OSError): os.pwrite(file.executable_descriptor, b"!", 0)
        for role in ("cc1", "cc1plus", "as", "collect2", "ld"):
            self.assertFalse((owner.root / "gcc-prefix" / role).exists())

    def test_arguments_preserve_logical_order_and_put_helpers_first(self):
        owner = self.enter()
        argv = self.argv()
        original = list(argv)
        effective = owner.arguments(argv, 123)
        self.assertEqual(argv, original)
        self.assertEqual(effective[:4], [argv[0], "-B/proc/self/fd/123/", "-B" + str(owner.root / "gcc-prefix") + "/",
                                        "--sysroot=" + str(owner.root)])
        self.assertEqual(effective[4:-2], argv[1:-2])
        self.assertEqual(effective[-2:], [str(owner.root / "gcc-prefix/libgomp.so"), str(owner.root / "gcc-prefix/libpthread.a")])

    def test_argument_overrides_and_compile_requests_fail_closed(self):
        for extra in ("-B/override", "@args", "-Wl,-T,bad", "-Xlinker", "-L/host", "-lbad", "--sysroot=/",
                      "-nostdlib", "-shared", "-c", "-E", "-S", "-specs=bad"):
            owner = self.enter()
            with self.subTest(extra=extra), self.assertRaises(FAILURES): owner.arguments(self.argv() + [extra], 123)
            with self.assertRaises(FAILURES): owner.record()

    def test_missing_or_duplicate_explicit_library_is_rejected(self):
        for argv in (self.argv()[:-1], self.argv() + [str(self.gcc / "libgomp.so")], self.argv() + ["-o", "second"]):
            owner = self.enter()
            with self.assertRaises(FAILURES): owner.arguments(argv, 123)

    def test_invalid_pins_counts_paths_and_byte_bounds(self):
        bad = [[], self.pins[:-1], self.pins + [self.pins[0]]]
        for key, value in (("size", True), ("size", 0), ("size", module.MAX_FILE_BYTES + 1), ("sha256", "bad"),
                           ("path", "/outside/input")):
            rows = copy.deepcopy(self.pins); rows[0][key] = value; bad.append(rows)
        rows = copy.deepcopy(self.pins); rows[-1] = rows[0]; bad.append(rows)
        for rows in bad:
            with self.assertRaises(FAILURES): self.owner(rows)
        with mock.patch.object(module, "MAX_TOTAL_BYTES", 1), self.assertRaises(FAILURES): self.owner()

    def test_changed_pin_and_aliased_file_rejected(self):
        path = Path(self.pins[0]["path"])
        path.write_bytes(b"changed")
        with self.assertRaises(FAILURES): self.owner().__enter__()
        path.unlink(); path.symlink_to(Path(self.pins[1]["path"]))
        with self.assertRaises(FAILURES): self.owner().__enter__()

    def test_source_write_restore_is_latched(self):
        owner = self.enter()
        path = Path(self.pins[0]["path"])
        original = path.read_bytes()
        path.write_bytes(b"!" * len(original)); path.write_bytes(original)
        with self.assertRaises(FAILURES): owner.validate_current()

    def test_preexisting_mmap_is_rehashed(self):
        path = Path(self.pins[0]["path"])
        with path.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[0:1] = mapping[0:1]
            owner = self.enter()
            file = owner._files[path]
            with mock.patch.object(file, "_hash", wraps=file._hash) as hashed:
                mapping[0] ^= 1
                try:
                    with self.assertRaises(FAILURES): owner.validate_current()
                    self.assertTrue(hashed.called)
                finally: mapping[0] ^= 1

    def test_view_permission_and_symlink_changes_are_rejected(self):
        owner = self.enter()
        directory = owner.root / "gcc-prefix"
        directory.chmod(0o700)
        path = directory / "libc.so"
        target = os.readlink(path)
        path.unlink(); path.symlink_to("/proc/self/fd/0")
        path.unlink(); path.symlink_to(target)
        directory.chmod(0o500)
        with self.assertRaises(FAILURES): owner.record()

    def test_runtime_loss_and_descriptor_inheritance_rejected(self):
        owner = self.enter()
        self.live = False
        with self.assertRaises(FAILURES): owner.record()
        self.live = True
        owner = self.enter()
        fd = owner.descriptors()[0]
        os.set_inheritable(fd, True)
        try:
            with self.assertRaises(FAILURES): owner.record()
        finally: os.set_inheritable(fd, False)

    def test_record_detachment_and_descriptor_closure(self):
        owner = self.enter()
        record = owner.record()
        record["mappings"].clear()
        self.assertNotEqual(record, owner.record())
        descriptors = [*owner.descriptors(), *(row[0] for row in owner._source_dirs.values()),
                       *(row[0] for row in owner._view_dirs.values())]
        owner.__exit__(None, None, None)
        for fd in descriptors:
            with self.assertRaises(OSError): os.fstat(fd)
        with self.assertRaises(FAILURES): owner.__enter__()

    def test_foreign_owner_or_runtime_rejected(self):
        owner = self.enter()
        with self.assertRaises(FAILURES): dispatch_module.RuntimeDispatch(self.inventory, link_inputs=object(), _runner=lambda *a: b"")
        with self.assertRaises(FAILURES): dispatch_module.RuntimeDispatch(copy.copy(self.inventory), link_inputs=owner, _runner=lambda *a: b"")


if __name__ == "__main__": unittest.main()
