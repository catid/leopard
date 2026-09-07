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
        # Keep disk writeback from changing the metadata of the rehash fixture.
        fixture_parent = "/dev/shm" if self._testMethodName in (
            "test_prefaulted_mmap_is_rehashed", "test_source_prefaulted_mmap_rehashes_quoted_input") else None
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-header-test-", dir=fixture_parent)))
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
        self.stack.enter_context(mock.patch.object(module, "C_PREDEFINITION_HEADER", str(self.first)))
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
        self.assertEqual(result.count("-include"), 1)
        self.assertEqual(result[result.index("-include") + 1], str(owner.root / str(self.first).lstrip("/")))
        self.assertEqual(owner.record()["schema"], "leopard2-v19-compiler-headers/v3")
        self.assertTrue(owner.record()["implicit_predefinition_sealed"])

    def test_recipe_overrides_latch_failure(self):
        for value in ("-isystem", "-include", "--sysroot=x", "-nostdinc++", "-Wp,-include,x", "-Xpreprocessor",
                      "-fno-canonical-system-headers", "-E"):
            owner = self.enter()
            with self.subTest(value=value), self.assertRaises(FAILURES):
                owner.arguments(["/usr/bin/c++", "-o", "out", "-c", "source.cpp", value])
            with self.assertRaises(FAILURES): owner.record()

    def c_phase(self):
        self.inventory.phase.language = "c"
        self.inventory.phase.logical_driver = "/usr/bin/cc"
        self.stack.enter_context(mock.patch.object(module, "C_INCLUDE_ROOTS", self.roots))
        self.stack.enter_context(mock.patch.object(module, "C_PREDEFINITION_HEADER", str(self.first)))

    def test_c_default_roots_and_compile_only_arguments(self):
        self.c_phase()
        owner = module.CompilerHeaders(self.inventory, self.pins, _file_factory=self.factory).__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        argv = ["/usr/bin/cc", "-c", "source.c", "-o", "out.o"]
        result = owner.arguments(argv)
        self.assertEqual(result[-4:], argv[1:])
        self.assertIn("-nostdinc", result)
        self.assertNotIn("-nostdinc++", result)
        self.assertEqual(result[result.index("-include") + 1], str(owner.root / str(self.first).lstrip("/")))
        self.assertEqual(owner.record()["language"], "c")
        self.assertEqual(owner.record()["include_roots"], list(self.roots))

    def test_c_combined_compile_link_preserves_arguments(self):
        self.c_phase()
        owner = self.enter()
        argv = ["/usr/bin/cc", "CMakeCCompilerId.c", "-o", "out"]
        original = list(argv)
        result = owner.arguments(argv, compile_only=False)
        self.assertEqual(argv, original)
        self.assertEqual(result[-3:], argv[1:])
        self.assertNotIn("-c", result)
        self.assertNotIn("-nostdinc++", result)
        self.assertEqual(result.count("-isystem"), len(self.roots))

    def test_c_stage_mismatch_or_nonboolean_mode_latches(self):
        self.c_phase()
        for mode, extra in ((True, []), (False, ["-c"]), (False, ["-E"]),
                            (False, ["-S"]), (0, []), (1, ["-c"]), (None, [])):
            owner = self.enter()
            with self.subTest(mode=mode, extra=extra), self.assertRaises(FAILURES):
                owner.arguments(["/usr/bin/cc", "s.c", "-o", "out", *extra], compile_only=mode)
            with self.assertRaises(FAILURES): owner.record()

    def test_c_missing_predefinition_pin_is_rejected(self):
        self.c_phase()
        with self.assertRaises(FAILURES): self.owner(self.pins[1:])

    def test_cpp_missing_predefinition_pin_is_rejected(self):
        with self.assertRaises(FAILURES): self.owner(self.pins[1:])

    def test_c_header_record_preserves_its_existing_schema(self):
        self.c_phase()
        owner = self.enter()
        self.assertEqual(owner.record()["schema"], "leopard2-v19-compiler-headers/v2")
        self.assertNotIn("implicit_predefinition_sealed", owner.record())

    def test_cpp_combined_compile_link_preserves_predefinition_and_cpp_roots(self):
        owner = self.enter()
        argv = ["/usr/bin/c++", "s.cpp", "-o", "out"]
        effective = owner.arguments(argv, compile_only=False)
        self.assertEqual(effective[-3:], argv[1:])
        self.assertIn("-nostdinc++", effective)
        self.assertNotIn("-c", effective)
        self.assertEqual(effective.count("-isystem"), len(self.roots))
        self.assertEqual(effective[effective.index("-include") + 1],
                         str(owner.root / module.C_PREDEFINITION_HEADER.lstrip("/")))

    def test_cpp_stage_mismatch_or_nonboolean_mode_latches(self):
        for mode, extra in ((True, []), (False, ["-c"]), (False, ["-E"]),
                            (False, ["-S"]), (0, []), (1, ["-c"]), (None, [])):
            owner = self.enter()
            with self.subTest(mode=mode, extra=extra), self.assertRaises(FAILURES):
                owner.arguments(["/usr/bin/c++", "s.cpp", "-o", "out", *extra], compile_only=mode)
            with self.assertRaises(FAILURES): owner.record()

    def test_c_rejects_cpp_roots_and_unknown_language(self):
        self.inventory.phase.language = "c"
        with self.assertRaises(FAILURES): self.owner()
        self.c_phase()
        with self.assertRaises(FAILURES): module.header_pins(self.pins, self.roots[::-1], "c")
        for language in ("fortran", None, True):
            with self.subTest(language=language), self.assertRaises(FAILURES):
                module.header_pins(self.pins, self.roots, language)

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

    def source_owner(self):
        self.c_phase()
        source = self.system / "unit.c"
        source.write_bytes(b'#include "one.h"\n#include "nested/two.h"\nint value = ONE;\n')
        source.chmod(0o644)
        owner = module.CompilerSourceInputs(self.inventory, self.pin(source), self.pins, _file_factory=self.factory).__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        return owner, source

    def test_source_unit_and_quoted_headers_are_sealed_together(self):
        owner, source = self.source_owner()
        record = owner.record()
        self.assertEqual(record["source"], str(source))
        self.assertEqual(record["input_bytes"], sum(row["size"] for row in record["files"]))
        self.assertEqual(len(record["files"]), 3)
        self.assertTrue(record["declared_source_inputs_sealed"])
        self.assertFalse(record["source_identity_owned"])
        self.assertFalse(record["full_source_read_closure_owned"])
        for path, file in owner._files.items():
            self.assertEqual((owner.root / str(path).lstrip("/")).read_bytes(), path.read_bytes())
            self.assertFalse(os.get_inheritable(file.executable_descriptor))
            with self.assertRaises(OSError): os.pwrite(file.executable_descriptor, b"!", 0)

    def test_source_argument_mapping_and_separate_link(self):
        owner, source = self.source_owner()
        for stage in (["-c"], []):
            argv = ["/usr/bin/cc", *stage, str(source), "-o", "out"]
            selected = owner.arguments(argv)
            self.assertIn(f"-ffile-prefix-map={owner.root}=", selected)
            self.assertNotIn(str(source), selected)
            self.assertIn(str(owner.root / str(source).lstrip("/")), selected)
            self.assertEqual(argv, ["/usr/bin/cc", *stage, str(source), "-o", "out"])
        link = ["/usr/bin/cc", "out.o", "-o", "out"]
        self.assertEqual(owner.arguments(link), link)
        self.assertIsNot(owner.arguments(link), link)

    def test_source_argument_mismatch_latches(self):
        for extra in (["another.c"], ["-E"], ["-S"], ["-c", "-c"], ["-o", "other"]):
            owner, source = self.source_owner()
            with self.subTest(extra=extra), self.assertRaises(FAILURES):
                owner.arguments(["/usr/bin/cc", str(source), "-o", "out", *extra])
            with self.assertRaises(FAILURES): owner.record()

    def test_source_absence_duplicate_and_wrong_driver_rejected(self):
        for case in ("absent", "duplicate", "driver"):
            owner, source = self.source_owner()
            argv = ["/usr/bin/cc", "-c", str(source), "-o", "out"]
            if case == "absent": argv[2] = "different.c"
            elif case == "duplicate": argv.append(str(source))
            else: argv[0] = "/usr/bin/c++"
            with self.subTest(case=case), self.assertRaises(FAILURES): owner.arguments(argv)

    def test_source_pin_scope_and_origin_rejected(self):
        self.c_phase()
        source = self.system / "unit.c"
        source.write_bytes(b"int value;\n"); source.chmod(0o644)
        pin = self.pin(source)
        for headers in ([pin], [dict(self.pins[0], path="/outside.h")], self.pins * 128):
            with self.assertRaises(FAILURES):
                module.CompilerSourceInputs(self.inventory, pin, headers, _file_factory=self.factory)
        for origin in (None, True, "untrusted"):
            with self.assertRaises(FAILURES):
                module.CompilerSourceInputs(self.inventory, pin, [], origin=origin, _file_factory=self.factory)
        with self.assertRaises(FAILURES):
            module.CompilerSourceInputs(self.inventory, self.pins[0], [], _file_factory=self.factory)
        with mock.patch.object(module, "MAX_TOTAL_BYTES", 1), self.assertRaises(FAILURES):
            module.CompilerSourceInputs(self.inventory, pin, [], _file_factory=self.factory)

    def test_source_generated_origin_uses_current_owner(self):
        self.c_phase()
        source = self.system / "unit.c"
        source.write_bytes(b"int value;\n"); source.chmod(0o644)
        with mock.patch.object(module.runtime, "RuntimeInventory", type(self.inventory)):
            owner = module.CompilerSourceInputs(self.inventory, self.pin(source), [], origin="generated").__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        self.assertEqual(owner.record()["origin"], "generated")
        value = os.fstat(owner._files[source].fd)
        self.assertEqual((value.st_uid, value.st_gid), (os.getuid(), os.getgid()))

    def test_source_cpp_mapping_and_language_mismatch(self):
        source = self.system / "unit.cpp"
        source.write_bytes(b'#include "one.h"\nint value = ONE;\n'); source.chmod(0o644)
        owner = module.CompilerSourceInputs(self.inventory, self.pin(source), self.pins, _file_factory=self.factory).__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        selected = owner.arguments(["/usr/bin/c++", "-c", str(source), "-o", "out.o"])
        self.assertEqual(owner.record()["language"], "c++")
        self.assertIn(str(owner.root / str(source).lstrip("/")), selected)
        self.inventory.phase.language = "c"
        with self.assertRaises(FAILURES):
            module.CompilerSourceInputs(self.inventory, self.pin(source), self.pins, _file_factory=self.factory)

    def test_source_system_origin_rejects_user_owned_inputs(self):
        if os.getuid() == 0 and os.getgid() == 0: self.skipTest("requires a non-root source owner")
        self.c_phase()
        source = self.system / "unit.c"
        source.write_bytes(b"int value;\n"); source.chmod(0o644)
        with mock.patch.object(module.runtime, "RuntimeInventory", type(self.inventory)), self.assertRaises(FAILURES):
            module.CompilerSourceInputs(self.inventory, self.pin(source), [], origin="system").__enter__()

    def test_source_and_quoted_header_write_restore_is_latched(self):
        for target in ("source", "quoted"):
            owner, source = self.source_owner()
            path = source if target == "source" else self.first
            original = path.read_bytes()
            path.write_bytes(b"!" * len(original)); path.write_bytes(original)
            with self.subTest(target=target), self.assertRaises(FAILURES): owner.record()

    def test_source_prefaulted_mmap_rehashes_quoted_input(self):
        with self.first.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[0:1] = mapping[0:1]
            owner, _ = self.source_owner()
            file = owner._files[self.first]
            with mock.patch.object(file, "_hash", wraps=file._hash) as hashed:
                mapping[0] ^= 1
                try:
                    with self.assertRaises(FAILURES): owner.record()
                    self.assertTrue(hashed.called)
                finally: mapping[0] ^= 1


if __name__ == "__main__": unittest.main()
