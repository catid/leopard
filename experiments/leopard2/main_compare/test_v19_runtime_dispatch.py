#!/usr/bin/python3
"""Static entry routing plus synthetic owner failures; no benchmark workload."""
from contextlib import ExitStack
import copy
import importlib.util
import mmap
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_runtime_dispatch", HERE / "v19_runtime_dispatch.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.builder.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


def static_fixture():
    data = bytearray(256)
    struct.pack_into("<16sHHIQQQIHHHHHH", data, 0, b"\x7fELF\x02\x01\x01" + bytes(9), 2, 62, 1,
                     0x4000c0, 64, 0, 0, 64, 56, 2, 0, 0, 0)
    struct.pack_into("<IIQQQQQQ", data, 64, 1, 5, 0, 0x400000, 0, 256, 256, 4096)
    struct.pack_into("<IIQQQQQQ", data, 120, 0x6474e551, 6, 0, 0, 0, 0, 0, 16)
    return bytes(data)


class SourceTests(unittest.TestCase):
    def setUp(self):
        self.bindings = {"loader_label": "/usr/lib/loader", "loader_descriptor": "/proc/self/fd/10",
                         "library_prefix": "/proc/self/fd/11", "plugin_original": "/usr/lib/plugin.so",
                         "plugin_descriptor": "/proc/self/fd/12",
                         "helpers": dict.fromkeys(("as", "cc1plus", "collect2", "ld"), "/proc/self/fd/13")}

    def test_bindings_render_both_languages_without_shell_syntax(self):
        template = module.TEMPLATE.read_bytes()
        for frontend in ("cc1", "cc1plus"):
            self.bindings["helpers"] = dict.fromkeys(("as", frontend, "collect2", "ld"), "/proc/self/fd/13")
            data = module.source_bytes(template, self.bindings)
            self.assertTrue(data.startswith(template))
            self.assertEqual(data.count(b".quad role_"), 4)
            self.assertIn(b'plugin_descriptor: .asciz "/proc/self/fd/12"', data)

    def test_unsafe_roles_descriptors_and_strings_rejected(self):
        for key, value in (("helpers", {}), ("loader_descriptor", "/proc/self/fd/0"),
                           ("plugin_descriptor", "/proc/self/fd/65536"), ("library_prefix", "/proc/self/fd/0011"),
                           ("plugin_original", '/bad/";evil'), ("loader_label", "/bad/../loader")):
            bindings = copy.deepcopy(self.bindings)
            bindings[key] = value
            with self.subTest(key=key), self.assertRaises(FAILURES): module.source_bytes(b"fixture\n", bindings)
        for fd in (True, -1, 2, 65536):
            with self.assertRaises(FAILURES): module.fd_path(fd)

    def test_static_elf_rejects_dynamic_writable_executable_or_bad_entry(self):
        module.validate_static_elf(static_fixture())
        for offset, fmt, value in ((16, "H", 3), (64, "I", 3), (64, "I", 2), (68, "I", 7),
                                   (124, "I", 7), (24, "Q", 0), (56, "H", 17), (96, "Q", 257)):
            data = bytearray(static_fixture())
            struct.pack_into("<"+fmt, data, offset, value)
            with self.subTest(offset=offset), self.assertRaises(FAILURES): module.validate_static_elf(bytes(data))

    def test_duplicate_executable_stack_declaration_is_rejected(self):
        data = bytearray(static_fixture())
        struct.pack_into("<H", data, 56, 3)
        struct.pack_into("<IIQQQQQQ", data, 176, 0x6474e551, 7, 0, 0, 0, 0, 0, 16)
        with self.assertRaises(FAILURES): module.validate_static_elf(bytes(data))


class RealTrampolineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.stack = ExitStack()
        cls.addClassCleanup(cls.stack.close)
        cls.root = Path(cls.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-entry-test-")))
        cls.loader = module.builder._StreamedTool(Path("/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2"))
        cls.stack.callback(cls.loader.close)
        cls.printf = module.builder._StreamedTool(Path("/usr/bin/printf"))
        cls.stack.callback(cls.printf.close)
        cls.printenv = module.builder._StreamedTool(Path("/usr/bin/printenv"))
        cls.stack.callback(cls.printenv.close)
        cls.libc = module.builder._StreamedTool(Path("/usr/lib/x86_64-linux-gnu/libc.so.6"), permitted_modes=(0o644, 0o755))
        cls.stack.callback(cls.libc.close)
        cls.libdir = cls.root / "libs"
        cls.libdir.mkdir(mode=0o700)
        (cls.libdir / "libc.so.6").symlink_to(module.fd_path(cls.libc.executable_descriptor))
        (cls.libdir / "ld-linux-x86-64.so.2").symlink_to(module.fd_path(cls.loader.executable_descriptor))
        cls.libfd = os.open(cls.libdir, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        cls.stack.callback(os.close, cls.libfd)
        cls.bindings = {"loader_label": str(cls.loader.path), "loader_descriptor": module.fd_path(cls.loader.executable_descriptor),
                        "library_prefix": module.fd_path(cls.libfd), "plugin_original": "/fixture/liblto_plugin.so",
                        "plugin_descriptor": module.fd_path(cls.libc.executable_descriptor),
                        "helpers": dict.fromkeys(("cc1plus", "as", "collect2", "ld"), module.fd_path(cls.printf.executable_descriptor))}
        cls.shim = cls.build_shim("fixture", cls.bindings)
        c_bindings = copy.deepcopy(cls.bindings)
        c_bindings["helpers"]["cc1"] = c_bindings["helpers"].pop("cc1plus")
        cls.c_shim = cls.build_shim("c-fixture", c_bindings)
        env_bindings = copy.deepcopy(cls.bindings)
        env_bindings["helpers"]["as"] = module.fd_path(cls.printenv.executable_descriptor)
        cls.env_shim = cls.build_shim("env-fixture", env_bindings)

    @classmethod
    def build_shim(cls, name, bindings):
        source, obj, binary = cls.root / (name+".S"), cls.root / (name+".o"), cls.root / name
        source.write_bytes(module.source_bytes(module.TEMPLATE.read_bytes(), bindings))
        # This small local assembler fixture tests the instruction stream;
        # pinned sealed bootstrap execution is checked by the native proof.
        subprocess.run(["/usr/bin/as", "--64", "-o", str(obj), str(source)], check=True, umask=0o022)
        subprocess.run(["/usr/bin/ld", "-static", "--build-id=none", "-z", "noexecstack", "-o", str(binary), str(obj)], check=True, umask=0o022)
        module.validate_static_elf(binary.read_bytes())
        shim = module.builder._StreamedTool(binary, _trusted_owner=(os.getuid(), os.getgid()))
        cls.stack.callback(shim.close)
        return shim

    def invoke(self, argv, shim=None):
        return module.provenance._run(argv, "static helper entry fixture", executable_descriptor=(shim or self.shim).executable_descriptor,
            inherited_descriptors=(self.libfd, self.loader.executable_descriptor, self.libc.executable_descriptor,
                                   self.printf.executable_descriptor, self.printenv.executable_descriptor),
            environment_overrides={"LEOPARD_TEST_MARKER": "literal marker"}, timeout=10)

    def test_all_roles_preserve_empty_unicode_whitespace_and_plugin_substitution(self):
        inputs = ["", "a b", "line\nbreak", "λ", "'quoted'", "/fixture/liblto_plugin.so", "/fixture/liblto_plugin.so-extra"]
        expected = [self.bindings["plugin_descriptor"] if value == self.bindings["plugin_original"] else value for value in inputs]
        for role in self.bindings["helpers"]:
            self.assertEqual(self.invoke(["/fixture/"+role, r"%s\0", *inputs]), b"".join(value.encode()+b"\0" for value in expected))

    def test_logical_argv0_is_preserved(self):
        output = self.invoke(["/fixture/as", "--help"])
        self.assertIn(b"Usage: /fixture/as FORMAT", output)

    def test_c_frontend_role_and_environment_are_preserved(self):
        self.assertEqual(self.invoke(["cc1", "%s", "c role"], self.c_shim), b"c role")
        with self.assertRaises(module.provenance.BuildProvenanceError):
            self.invoke(["cc1plus", "%s", "wrong language"], self.c_shim)
        self.assertEqual(self.invoke(["as", "LEOPARD_TEST_MARKER"], self.env_shim), b"literal marker\n")

    def test_invalid_role_long_name_and_argument_count_fail(self):
        for argv in (["/fixture/unregistered", "hello"], ["/" + "x"*4096 + "/as", "hello"],
                     ["/fixture/as", r"%s\0", *(["x"]*511)]):
            with self.assertRaises(module.provenance.BuildProvenanceError): self.invoke(argv)
        self.assertEqual(self.invoke(["/fixture/as", r"%s\0", *(["x"]*510)]), b"x\0"*510)


class OwnerTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-dispatch-test-")))
        self.parent = self.root / "new"
        self.parent.mkdir(mode=0o700)
        self.parentfd = os.open(self.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        self.stack.callback(os.close, self.parentfd)
        self.libfd = os.open(self.root, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        self.stack.callback(os.close, self.libfd)
        tools = {}
        for role in ("driver", "cc1plus", "as", "collect2", "ld", "loader", "plugin"):
            path = self.root / role
            path.write_bytes(b"fixture tool " + role.encode())
            path.chmod(0o755)
            tool = module.builder._StreamedTool(path, _trusted_owner=(os.getuid(), os.getgid()))
            self.stack.callback(tool.close)
            tools[role] = tool
        self.tools = tools
        self.live = True
        test = self
        class Phase:
            parent, _parent_fd = test.parent, test.parentfd
            _tools = {role: tools[role] for role in ("driver", "cc1plus", "as", "collect2", "ld")}
            pins = {role: {"path": str(tool.path)} for role, tool in _tools.items()}
            logical_driver = "/usr/bin/c++"
            def validate_current(self): module.require(test.live, "fixture phase lost")
        class Inventory:
            phase = Phase()
            loader_path, loader_key = tools["loader"].path, str(tools["loader"].path)
            _files = {str(tools[role].path): tools[role] for role in ("loader", "plugin")}
            _roots = {"gcc-link-plugin": str(tools["plugin"].path)}
            _lists = {}
            prefix = type("Prefix", (), {"descriptor": test.libfd})()
            def validate_current(self): self.phase.validate_current()
            def list_dependencies(self, role): self._lists[role] = {"fixture": True}
        self.inventory = Inventory()
        self.calls = []
        self.intercept = lambda *_: None

    def run_child(self, argv, label, **kwargs):
        self.calls.append((argv, kwargs))
        self.assertEqual(kwargs["executable_descriptor"], self.tools["loader"].executable_descriptor)
        self.assertEqual(argv[:8], [str(self.tools["loader"].path), "--inhibit-cache", "--glibc-hwcaps-mask", "",
                                   "--library-path", module.fd_path(self.libfd), "--argv0", argv[7]])
        role = label.rsplit(" ", 1)[-1]
        if role in ("as", "ld"):
            path = Path(argv[argv.index("-o")+1])
            path.write_bytes(b"fixture object" if role == "as" else static_fixture())
            path.chmod(0o644 if role == "as" else 0o755)
        self.intercept(role)
        return b""

    def owner(self): return module.RuntimeDispatch(self.inventory, _runner=self.run_child)

    def test_bootstrap_driver_route_false_claims_and_umask_restore(self):
        initial = os.umask(0o077)
        try:
            with self.owner() as owner:
                self.assertEqual(os.umask(0o077), 0o077)
                owner.run(["/usr/bin/c++", "-c", "source.cpp", "-o", "out.o"])
                record = owner.record()
                self.assertEqual([row["role"] for row in record["commands"]], ["as", "ld", "driver"])
                self.assertEqual(record["sealed_dispatcher"]["seals"], 15)
                self.assertEqual(len(record["artifacts"]), 4)
                for key, value in record.items():
                    if type(value) is bool: self.assertIs(value, False)
                self.assertEqual(record["commands"][-1]["logical_argv"], ["/usr/bin/c++", "-c", "source.cpp", "-o", "out.o"])
                self.assertIn(f"-B/proc/self/fd/{owner.prefix.descriptor}/", record["commands"][-1]["effective_argv"])
                self.assertEqual([row["umask_policy"] for row in record["commands"]], ["0022", "0022", "caller"])
                record["commands"].clear()
                self.assertEqual(len(owner.record()["commands"]), 3)
            with self.assertRaises(FAILURES): owner.record()
        finally: os.umask(initial)

    def test_bootstrap_failure_restores_umask_and_latches(self):
        def fail(role): raise RuntimeError("bootstrap failed")
        self.intercept = fail
        initial = os.umask(0o077)
        try:
            owner = self.owner()
            with self.assertRaises(RuntimeError): owner.__enter__()
            self.assertEqual(os.umask(0o077), 0o077)
            with self.assertRaises(FAILURES): owner.record()
        finally: os.umask(initial)

    def test_borrowed_owner_loss_during_job_fails(self):
        with self.assertRaises(FAILURES):
            with self.owner() as owner:
                def lose(role): self.live = False
                self.intercept = lose
                owner.run(["/usr/bin/c++"])

    def test_metadata_helper_prefix_and_source_drift(self):
        for target in ("root", "prefix", "source"):
            with self.subTest(target=target), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    if target == "source": (owner.root / "dispatch.S").write_bytes(b"changed")
                    else:
                        path = owner.root if target == "root" else owner.root / "helpers"
                        mode = path.stat().st_mode & 0o777
                        path.chmod(0o755)
                        path.chmod(mode)
                    owner.record()

    def test_bad_driver_arguments_do_not_launch(self):
        for argv in (None, ["/usr/bin/cc"], ["/usr/bin/c++", "@args"], ["/usr/bin/c++", "-B/tmp"], ["/usr/bin/c++"]*512):
            with self.subTest(argv=argv), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    count = len(self.calls)
                    with self.assertRaises(FAILURES): owner.run(argv)
                    self.assertEqual(len(self.calls), count)

    def test_wrong_authority_and_reuse_rejected(self):
        with self.assertRaises(FAILURES): module.RuntimeDispatch(self.inventory)
        with self.owner() as owner: pass
        with self.assertRaises(FAILURES): owner.__enter__()

    def test_source_fdopen_failure_closes_descriptor(self):
        real_fdopen, observed = os.fdopen, []
        def interrupt(fd, mode, *args, **kwargs):
            if mode == "wb":
                observed.append(fd)
                raise KeyboardInterrupt("source file-object construction interrupted")
            return real_fdopen(fd, mode, *args, **kwargs)
        with mock.patch.object(module.os, "fdopen", side_effect=interrupt), self.assertRaises(KeyboardInterrupt):
            self.owner().__enter__()
        self.assertEqual(len(observed), 1)
        with self.assertRaises(OSError): os.fstat(observed[0])

    def test_descriptor_inheritance_change_is_rejected(self):
        for descriptor in ("_root_fd", "_prefix_fd", "retained_prefix"):
            with self.subTest(descriptor=descriptor), self.assertRaises(FAILURES):
                with self.owner() as owner:
                    fd = owner.prefix.descriptor if descriptor == "retained_prefix" else getattr(owner, descriptor)
                    os.set_inheritable(fd, True)
                    owner.validate_current()

    def test_source_mmap_write_without_notification_is_rejected(self):
        with ExitStack() as maps:
            owner = self.owner()
            original_hold, observed = owner._hold, []
            def hold(path):
                if path.name == "dispatch.S":
                    stream = maps.enter_context(path.open("r+b"))
                    mapping = maps.enter_context(mmap.mmap(stream.fileno(), 0))
                    mapping[0:1] = mapping[0:1]  # Fault before snapshot acquisition.
                    observed.append(mapping)
                return original_hold(path)
            owner._hold = hold
            eventless = rejected_by_rehash = False
            with self.assertRaises(FAILURES):
                with owner:
                    source = owner.root / "dispatch.S"
                    fields = module.provenance._stable_fields(source.stat())
                    observed[0][0:1] = b"!"
                    self.assertEqual(module.provenance._stable_fields(source.stat()), fields)
                    owner._snapshots[str(source)]._path_guard.verify()
                    eventless = True
                    with self.assertRaisesRegex(module.builder.host.PreflightError, "current artifact bytes differ"):
                        owner.record()
                    rejected_by_rehash = True
            # Cleanup may also reject drift; it must not mask a broken fixture
            # or an assertion before the intended rehash check was reached.
            self.assertTrue(eventless)
            self.assertTrue(rejected_by_rehash)


if __name__ == "__main__": unittest.main()
