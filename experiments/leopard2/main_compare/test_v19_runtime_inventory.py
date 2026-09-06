#!/usr/bin/python3
"""Bounded synthetic ELF and real-fd runtime inventory tests; no compiler jobs."""
import hashlib
import importlib.util
import mmap
import os
from pathlib import Path
import stat
import struct
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_runtime_inventory", HERE / "v19_runtime_inventory.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.builder.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


def elf(needed=(), soname=None, interpreter=False):
    data = bytearray(2048)
    strings, entries = bytearray(b"\0"), []
    for name in needed:
        entries.append((1, len(strings)))
        strings += name.encode() + b"\0"
    if soname is not None:
        entries.append((14, len(strings)))
        strings += soname.encode() + b"\0"
    entries += [(5, 0x400400), (10, len(strings)), (0, 0)]
    phnum = 3 if interpreter else 2
    struct.pack_into("<16sHHIQQQIHHHHHH", data, 0, b"\x7fELF\x02\x01\x01" + bytes(9),
                     3, 62, 1, 0, 64, 0, 0, 64, 56, phnum, 0, 0, 0)
    struct.pack_into("<IIQQQQQQ", data, 64, 1, 5, 0, 0x400000, 0, len(data), len(data), 4096)
    struct.pack_into("<IIQQQQQQ", data, 120, 2, 4, 256, 0x400100, 0, len(entries)*16, len(entries)*16, 8)
    if interpreter:
        interp = module.INTERPRETER.encode() + b"\0"
        struct.pack_into("<IIQQQQQQ", data, 176, 3, 4, 768, 0x400300, 0, len(interp), len(interp), 1)
        data[768:768+len(interp)] = interp
    for index, entry in enumerate(entries): struct.pack_into("<qQ", data, 256 + index*16, *entry)
    data[1024:1024+len(strings)] = strings
    return bytes(data)


class ELFTests(unittest.TestCase):
    def parse(self, data):
        reads = []
        def read(offset, count): reads.append(count); return bytes(data[offset:offset+count])
        result = module.elf_startup(read, len(data))
        self.assertLessEqual(max(reads), 65536)
        return result

    def test_exact_transitive_metadata_and_interpreter(self):
        self.assertEqual(self.parse(elf(("liba.so.1", "libb.so.2"), "libroot.so", True)),
                         {"needed": ["liba.so.1", "libb.so.2"], "soname": "libroot.so", "interpreter": module.INTERPRETER})
        self.assertEqual(self.parse(elf()), {"needed": [], "soname": None, "interpreter": None})

    def test_header_and_table_bounds(self):
        for offset, fmt, value in ((4, "B", 1), (5, "B", 2), (18, "H", 3), (32, "Q", 1),
                                   (54, "H", 55), (56, "H", 0), (56, "H", 129), (32, "Q", 2030),
                                   (152, "Q", 65552), (152, "Q", 17)):
            with self.subTest(offset=offset, value=value):
                data = bytearray(elf())
                struct.pack_into("<"+fmt, data, offset, value)
                with self.assertRaises(FAILURES): self.parse(data)
        for size in (True, 0, 63, 65 << 20):
            with self.assertRaises(FAILURES): module.elf_startup(lambda *_: b"", size)

    def test_hidden_search_and_loading_tags_rejected(self):
        for tag in (15, 29, 0x6ffffefa, 0x6ffffefb, 0x6ffffefc, 0x7ffffffd, 0x7fffffff):
            data = bytearray(elf(("liba.so",)))
            struct.pack_into("<q", data, 256, tag)
            with self.subTest(tag=tag), self.assertRaises(FAILURES): self.parse(data)

    def test_string_metadata_and_termination(self):
        for names in (("../escape.so",), ("liba.so", "liba.so"), ("a"*256,), ("",)):
            with self.assertRaises(FAILURES): self.parse(elf(names))
        for offset, value in ((264, 4096), (296, 0), (280, 0x500000), (312, 2)):
            data = bytearray(elf(("liba.so",)))
            struct.pack_into("<Q", data, offset, value)
            with self.subTest(offset=offset), self.assertRaises(FAILURES): self.parse(data)
        with self.assertRaises(FAILURES): module.elf_startup(lambda *_: b"", 2048)

    def test_interpreter_must_match_supported_host_profile(self):
        data = bytearray(elf(interpreter=True))
        data[768] = ord("x")
        with self.assertRaises(FAILURES): self.parse(data)

    def test_gcc_sized_string_table_is_not_materialized(self):
        data = bytearray(elf(("liba.so",)))
        data.extend(bytes(2 << 20))
        struct.pack_into("<QQ", data, 96, len(data), len(data))
        struct.pack_into("<Q", data, 296, 1199001)
        reads = []
        def read(offset, count): reads.append(count); return bytes(data[offset:offset+count])
        self.assertEqual(module.elf_startup(read, len(data))["needed"], ["liba.so"])
        self.assertLessEqual(max(reads), 256)
        self.assertLess(sum(reads), 1024)


class InventoryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-runtime-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.libs = self.root / "libs"
        self.libs.mkdir(mode=0o700)
        self.parent = self.root / "new"
        self.parent.mkdir(mode=0o700)
        self.loader = self.libs / "ld-linux-x86-64.so.2"
        self.loader.write_bytes(elf(soname=self.loader.name))
        self.loader.chmod(0o755)
        self.leaf = self.libs / "libleaf.so.1.2"
        self.leaf.write_bytes(elf(soname="libleaf.so.1"))
        self.leaf.chmod(0o644)
        (self.libs / "libleaf.so.1").symlink_to(self.leaf.name)
        self.library = self.libs / "libtest.so.1"
        self.library.write_bytes(elf(("libleaf.so.1",), "libtest.so.1"))
        self.library.chmod(0o644)
        self.driver_path = self.root / "driver"
        self.driver_path.write_bytes(elf(("libtest.so.1",), interpreter=True))
        self.driver_path.chmod(0o755)
        self.plugin = self.root / "liblto_plugin.so"
        self.plugin.write_bytes(elf(("libtest.so.1",)))
        self.plugin.chmod(0o644)
        self.tools = []
        self.driver = self.factory(self.driver_path)
        self.addCleanup(self.driver.close)
        fd = os.open(self.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        self.addCleanup(os.close, fd)
        self.live = True
        test = self
        class Phase:
            _tools = {"driver": test.driver}
            _parent_fd, parent = fd, test.parent
            pins = {"collect2": {"path": str(test.root / "collect2")}}
            def validate_current(self):
                module.require(test.live, "fixture phase is no longer held")
                test.driver.validate_current()
        self.phase = Phase()

    def factory(self, path, **kwargs):
        tool = module.builder._StreamedTool(path, _trusted_owner=(os.geteuid(), os.getegid()), **kwargs)
        self.tools.append(tool)
        return tool

    def owner(self):
        return module.RuntimeInventory(self.phase, _tool_factory=self.factory, _library_root=self.libs,
                                       _loader=self.loader, _preload=self.root / "preload")

    def test_transitive_sealed_inventory_and_false_claims(self):
        with self.owner() as owner:
            record = owner.record()
            self.assertEqual(set(record["sonames"]), {"libtest.so.1", "libleaf.so.1"})
            self.assertEqual(len(record["files"]), 4)
            self.assertEqual(record["bytes"], 8192)
            for key, value in record.items():
                if type(value) is bool: self.assertIs(value, False, key)
            self.assertEqual(record["roots"]["gcc-link-plugin"], str(self.plugin))
            for name, key in owner._sonames.items():
                self.assertEqual(os.readlink(owner.root / name), f"/proc/self/fd/{owner._files[key].executable_descriptor}")
                self.assertEqual(owner._files[key].executable_record()["seals"], 15)
            record["files"].clear()
            self.assertEqual(len(owner.record()["files"]), 4)
        with self.assertRaises(FAILURES): owner.record()

    def test_listing_uses_sealed_loader_and_all_library_fds(self):
        with self.owner() as owner:
            output = (f"linux-vdso.so.1 (0x1000)\n"
                      f"libtest.so.1 => /proc/self/fd/{owner.prefix.descriptor}/libtest.so.1 (0x2000)\n"
                      f"libleaf.so.1 => /proc/self/fd/{owner.prefix.descriptor}/libleaf.so.1 (0x3000)\n"
                      f"{module.INTERPRETER} => {self.loader} (0x4000)\n").encode()
            def run(argv, label, **kwargs):
                self.assertEqual(argv[1:6], ["--list", "--inhibit-cache", "--glibc-hwcaps-mask", "", "--library-path"])
                self.assertEqual(kwargs["executable_descriptor"], owner._files[owner.loader_key].executable_descriptor)
                self.assertEqual(len(kwargs["inherited_descriptors"]), 6)
                return output
            with mock.patch.object(module.provenance, "_run", side_effect=run):
                self.assertEqual(owner.list_dependencies("driver"), output)
            self.assertEqual(owner.record()["loader_listings"]["driver"]["stdout_sha256"], hashlib.sha256(output).hexdigest())
            self.assertIs(owner.record()["loader_listings"]["driver"]["sealed_prefix_resolution_verified"], True)

    def test_missing_extra_disk_or_duplicate_loader_resolution_rejected(self):
        with self.owner() as owner:
            prefix = f"/proc/self/fd/{owner.prefix.descriptor}"
            good = ("linux-vdso.so.1 (0x1000)\n"
                    f"libtest.so.1 => {prefix}/libtest.so.1 (0x2000)\n"
                    f"libleaf.so.1 => {prefix}/libleaf.so.1 (0x3000)\n"
                    f"{module.INTERPRETER} => {self.loader} (0x4000)\n")
            for text in (good.replace(prefix, "/usr/lib"), good + good.splitlines()[1] + "\n",
                         good.replace(good.splitlines()[2] + "\n", ""), good.replace(str(self.loader), "/bad/loader"),
                         good + f"other.so => {prefix}/other.so (0x5000)\n", "libtest.so.1 => not found\n"):
                with self.assertRaises(FAILURES):
                    module.validate_listing(text.encode(), str(self.driver_path), owner._metadata, owner._sonames, prefix, str(self.loader))

    def test_bad_listing_latches_in_owner(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner, mock.patch.object(module.provenance, "_run", return_value=b"wrong listing\n"):
                with self.assertRaises(FAILURES): owner.list_dependencies("driver")
                with self.assertRaises(FAILURES): owner.record()

    def test_plugin_uses_direct_loader_display_and_cannot_overwrite_listing(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                prefix = f"/proc/self/fd/{owner.prefix.descriptor}"
                output = ("linux-vdso.so.1 (0x1000)\n"
                          f"libtest.so.1 => {prefix}/libtest.so.1 (0x2000)\n"
                          f"libleaf.so.1 => {prefix}/libleaf.so.1 (0x3000)\n"
                          f"{self.loader} (0x4000)\n").encode()
                with mock.patch.object(module.provenance, "_run", return_value=output) as run:
                    self.assertEqual(owner.list_dependencies("gcc-link-plugin"), output)
                    with self.assertRaises(FAILURES): owner.list_dependencies("gcc-link-plugin")
                    self.assertEqual(run.call_count, 1)

    def test_missing_library_wrong_soname_and_escape_rejected(self):
        self.leaf.unlink()
        with self.assertRaises(FAILURES):
            with self.owner(): self.fail("accepted missing library")
        self.leaf.write_bytes(elf(soname="wrong.so"))
        self.leaf.chmod(0o644)
        with self.assertRaises(FAILURES):
            with self.owner(): self.fail("accepted wrong SONAME")
        self.leaf.unlink()
        self.leaf.symlink_to(self.driver_path)
        with self.assertRaises(FAILURES):
            with self.owner(): self.fail("accepted escaped library")

    def test_alias_and_size_bounds(self):
        self.plugin.unlink()
        os.link(self.loader, self.plugin)
        with self.assertRaises(FAILURES):
            with self.owner(): self.fail("accepted aliased plugin")
        self.plugin.unlink()
        self.plugin.write_bytes(elf())
        self.plugin.chmod(0o644)
        for key in ("MAX_BYTES", "MAX_FILES"):
            with mock.patch.object(module, key, 1), self.assertRaises(FAILURES):
                with self.owner(): self.fail("accepted exceeded inventory bound")

    def test_absent_global_preload_history_is_retained(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                (self.root / "preload").touch()
                (self.root / "preload").unlink()
                with self.assertRaises(FAILURES): owner.validate_current()
                with self.assertRaises(FAILURES): owner.record()

    def test_dependency_symlink_and_prefix_mode_changes_latch(self):
        for target in ("library", "prefix", "descriptor"):
            with self.subTest(target=target), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    if target == "library":
                        path = self.libs / "libleaf.so.1"
                        path.unlink()
                        path.symlink_to(self.leaf.name)
                    elif target == "prefix":
                        owner.root.chmod(0o700)
                        owner.root.chmod(0o500)
                    else: os.set_inheritable(owner.prefix.descriptor, True)
                    owner.record()

    def test_prefaulted_library_mmap_drift(self):
        with self.leaf.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[2000] = mapping[2000]
            with self.assertRaises(FAILURES):
                with self.owner() as owner:
                    tool = owner._files[str(self.leaf)]
                    sealed = tool.executable_descriptor
                    fields = module.provenance._stable_fields(os.fstat(tool.fd))
                    mapping[2000] ^= 1
                    self.assertEqual(module.provenance._stable_fields(os.fstat(tool.fd)), fields)
                    tool.guard.verify()
                    self.assertEqual(os.pread(sealed, 1, 2000), b"\0")
                    owner.validate_current()

    def test_query_failure_and_borrowed_phase_loss(self):
        with self.assertRaises(RuntimeError):
            with self.owner() as owner, mock.patch.object(module.provenance, "_run", side_effect=RuntimeError("child failed")):
                owner.list_dependencies("driver")
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                self.live = False
                owner.record()

    def test_unknown_role_and_wrong_authority(self):
        with self.assertRaises(FAILURES): module.RuntimeInventory(self.phase)
        with self.assertRaises(FAILURES):
            with self.owner() as owner: owner.list_dependencies("unregistered")

    def test_original_launcher_mode_default_is_unchanged(self):
        with self.assertRaises(FAILURES): self.factory(self.library)
        for modes in ([], (0o777,), (True,), (0o755, 0o644)):
            with self.assertRaises(FAILURES): self.factory(self.driver_path, permitted_modes=modes)


if __name__ == "__main__": unittest.main()
