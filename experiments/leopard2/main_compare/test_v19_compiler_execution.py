#!/usr/bin/python3
"""Small real-fd fixtures; actual GCC equivalence requires the native probe."""
from contextlib import contextmanager
import copy
import fcntl
import hashlib
import importlib.util
import json
import mmap
import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_compiler_phase", HERE / "v19_compiler_execution.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


class CompilerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-compiler-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.parent = self.root / "new"
        self.parent.mkdir(mode=0o700)
        body = Path("/usr/bin/true").read_bytes()
        self.pins = {}
        for name in ("cc", "c++", "cc1", "cc1plus", "as", "collect2", "ld"):
            path = self.root / name
            path.write_bytes(body)
            path.chmod(0o755)
            self.pins[name] = {"path": str(path), "sha256": hashlib.sha256(body).hexdigest(),
                               "size": len(body), "uid": 0, "gid": 0, "mode": stat.S_IFREG | 0o755}
        self.pinned = {"compiler": self.pins["c++"], "c_compiler": self.pins["cc"],
                       "compiler_subtools": [{"language": language, "role": role, "identity": self.pins[role]}
                         for language, frontend in (("c", "cc1"), ("c++", "cc1plus"))
                         for role in (frontend, "as", "collect2", "ld")]}
        self.live = True
        self.tools = []
        test = self

        class Retained:
            def validate_current(self): module.require(test.live, "fixture preflight released")
            def _bytes(self, name):
                test.assertEqual(name, "candidate-build-provenance.json")
                return json.dumps(test.pinned).encode()

        self.retained = Retained()
        self.alias_bin, self.alias_alt = self.root / "aliases", self.root / "alternatives"
        self.alias_bin.mkdir()
        self.alias_alt.mkdir()
        for name in ("cc", "c++"):
            (self.alias_bin / name).symlink_to("../alternatives/" + name)
            (self.alias_alt / name).symlink_to("../" + name)
        retain_aliases = module._retain_driver_aliases
        def fixture_aliases(logical, tool):
            path = self.alias_bin / Path(logical).name if logical in ("/usr/bin/cc", "/usr/bin/c++") else logical
            return retain_aliases(path, tool)
        aliases_patch = mock.patch.object(module, "_retain_driver_aliases", side_effect=fixture_aliases)
        aliases_patch.start()
        self.addCleanup(aliases_patch.stop)

    def factory(self, path, **kwargs):
        tool = module.builder._StreamedTool(path, _trusted_owner=(os.geteuid(), os.getegid()), **kwargs)
        self.tools.append(tool)
        return tool

    def owner(self, language="c++"):
        return module.CompilerExecution(self.retained, self.parent, language, _tool_factory=self.factory)

    @contextmanager
    def opened(self, language="c++"):
        with self.owner(language) as owner:
            yield owner

    def test_exact_roles_for_both_languages_and_detached_inventory(self):
        for language, frontend in (("c", "cc1"), ("c++", "cc1plus")):
            logical, rows = module.phase_inventory(self.pinned, language)
            self.assertEqual(logical, "/usr/bin/cc" if language == "c" else "/usr/bin/c++")
            self.assertEqual(set(rows), {"driver", frontend, "as", "collect2", "ld"})
            rows["driver"]["size"] = 0
            self.assertGreater(self.pinned["compiler"]["size"], 0)

    def test_malformed_or_incomplete_inventory(self):
        for change in (
            lambda p: p["compiler_subtools"].pop(),
            lambda p: p["compiler_subtools"].__setitem__(0, p["compiler_subtools"][1]),
            lambda p: p["compiler_subtools"][0].update(language=[]),
            lambda p: p["compiler_subtools"][0].update(role={}),
            lambda p: p["compiler"].update(uid=False),
            lambda p: p["compiler"].update(size=True),
            lambda p: p["compiler"].update(sha256="0" * 63),
            lambda p: p["compiler"].update(path=self.pins["cc1plus"]["path"]),
            lambda p: p["compiler"].update(size=49 << 20),
            lambda p: p["compiler"].update(mode=0o755),
        ):
            with self.subTest(change=change):
                pinned = copy.deepcopy(self.pinned)
                change(pinned)
                with self.assertRaises(FAILURES): module.phase_inventory(pinned, "c++")
        for language in ([], None, "fortran"):
            with self.assertRaises(FAILURES): module.phase_inventory(self.pinned, language)

    def test_sealed_fd_route_exact_prefix_and_no_runtime_claim(self):
        for language in ("c", "c++"):
            with self.opened(language) as owner:
                argv = [owner.logical_driver, "-O3", "-c", "source.cpp", "-o", "out.o"]
                def run(effective, label, **kwargs):
                    self.assertEqual(effective, [argv[0], f"-B/proc/self/fd/{owner.prefix.descriptor}/", *argv[1:]])
                    self.assertEqual(kwargs["environment_overrides"], module.builder.ENVIRONMENT)
                    self.assertEqual(kwargs["executable_descriptor"], owner._tools["driver"].executable_descriptor)
                    self.assertEqual(len(set(kwargs["inherited_descriptors"])), 6)
                    for role, tool in owner._tools.items():
                        fd = tool.executable_descriptor
                        self.assertNotEqual(fd, tool.fd)
                        self.assertEqual(fcntl.fcntl(fd, getattr(fcntl, "F_GET_SEALS", 1034)), 15)
                        if role != "driver": self.assertEqual(os.readlink(owner.root / role), f"/proc/self/fd/{fd}")
                    return b"fixture\n"
                with mock.patch.object(module.provenance, "_run", side_effect=run):
                    self.assertEqual(owner.run(argv), b"fixture\n")
                record = owner.record()
                self.assertEqual(record["commands"][0]["logical_argv"], argv)
                for key in ("runtime_closure_verified", "compiler_data_owned", "fresh_build_recipe_integrated",
                            "live_acquisition_armed", "benchmark_executed", "compiler_subtool_execution_owned",
                            "atomic_snapshot"):
                    self.assertIs(record[key], False)
                record["commands"].clear()
                self.assertEqual(len(owner.record()["commands"]), 1)
            with self.assertRaises(FAILURES): owner.record()

    def test_real_fixture_driver_executes_but_does_not_claim_real_gcc(self):
        with self.opened() as owner:
            self.assertEqual(owner.run(["/usr/bin/c++"]), b"")
            self.assertEqual(owner.record()["commands"][0]["status"], "exit-zero")

    def test_override_rejection_latches_before_any_child(self):
        for option in ("@args", "-B/tmp/", "-specs=custom", "--specs=custom", "-wrapper", "-fplugin=x",
                       "-fuse-ld=gold", "-flto", "-xc", "bad\0arg"):
            with self.subTest(option=option), self.assertRaises(FAILURES):
                with self.opened() as owner, mock.patch.object(module.provenance, "_run") as run:
                    with self.assertRaises(FAILURES): owner.run(["/usr/bin/c++", option])
                    run.assert_not_called()
                    with self.assertRaises(FAILURES): owner.record()

    def test_driver_and_argv_type_rejected(self):
        for argv in ([], ("/usr/bin/c++",), ["/usr/bin/cc"], ["/usr/bin/c++", 3]):
            with self.subTest(argv=argv), self.assertRaises(FAILURES):
                with self.opened() as owner: owner.run(argv)

    def test_child_failure_latches_and_closes_descriptors(self):
        owner = self.owner()
        with self.assertRaises(RuntimeError):
            with owner, mock.patch.object(module.provenance, "_run", side_effect=RuntimeError("child failed")):
                owner.run(["/usr/bin/c++"])
        self.assertEqual(owner._commands[-1]["status"], "failed")
        for tool in self.tools:
            with self.assertRaises(OSError): os.fstat(tool.fd)

    def test_pinned_digest_drift_fails_and_closes(self):
        self.pinned["compiler"]["sha256"] = "0" * 64
        with self.assertRaises(FAILURES):
            with self.opened(): self.fail("accepted wrong compiler")
        for tool in self.tools:
            with self.assertRaises(OSError): os.fstat(tool.fd)

    def test_hardlinked_roles_rejected(self):
        alias = self.root / "alias"
        os.link(self.root / "c++", alias)
        self.pins["cc1plus"]["path"] = str(alias)
        with self.assertRaises(FAILURES):
            with self.opened(): self.fail("accepted helper inode alias")

    def test_preflight_loss_before_or_after_child(self):
        self.live = False
        with self.assertRaises(FAILURES):
            with self.opened(): self.fail("accepted released preflight")
        self.live = True
        with self.assertRaises(FAILURES):
            with self.opened() as owner:
                def release(*args, **kwargs): self.live = False; return b""
                with mock.patch.object(module.provenance, "_run", side_effect=release): owner.run(["/usr/bin/c++"])

    def test_prefix_parent_and_descriptor_mutations_latch(self):
        for target in ("parent", "prefix", "prefix-fd", "root-fd"):
            with self.subTest(target=target), self.assertRaises(FAILURES):
                with self.opened() as owner:
                    if target in ("parent", "prefix"):
                        path = self.parent if target == "parent" else owner.root
                        initial = stat.S_IMODE(path.stat().st_mode)
                        path.chmod(0o755)
                        path.chmod(initial)
                    else:
                        fd = owner.prefix.descriptor if target == "prefix-fd" else owner._root_fd
                        os.set_inheritable(fd, True)
                    with self.assertRaises(FAILURES): owner.validate_current()
                    with self.assertRaises(FAILURES): owner.record()

    def test_parent_replacement_cannot_redirect_creation(self):
        displaced = self.root / "old-parent"
        foreign = self.root / "foreign"
        foreign.mkdir(mode=0o700)
        mkdir = os.mkdir
        def replace(path, *args, **kwargs):
            self.parent.rename(displaced)
            self.parent.symlink_to(foreign, target_is_directory=True)
            return mkdir(path, *args, **kwargs)
        with self.assertRaises(FAILURES), mock.patch.object(module.os, "mkdir", side_effect=replace):
            with self.opened(): self.fail("accepted replaced parent")
        self.assertEqual(list(foreign.iterdir()), [])
        self.assertEqual(len(list(displaced.iterdir())), 1)

    def test_constructor_failure_and_interruption_release_tools(self):
        for failure in (RuntimeError("factory failed"), KeyboardInterrupt()):
            def factory(path, **kwargs):
                if len(self.tools) % 2: raise failure
                return self.factory(path, **kwargs)
            self.tools.clear()
            with self.assertRaises(type(failure)):
                with module.CompilerExecution(self.retained, self.parent, "c++", _tool_factory=factory):
                    self.fail("factory failure ignored")
            for tool in self.tools:
                with self.assertRaises(OSError): os.fstat(tool.fd)

    def test_process_guard_and_reuse(self):
        with self.assertRaises(FAILURES): module.CompilerExecution(self.retained, self.parent, "c++")
        owner = self.owner()
        owner._pid = -1
        with self.assertRaises(FAILURES): owner.__enter__()
        with self.assertRaises(FAILURES): owner.__enter__()
        with self.opened() as owner: pass
        with self.assertRaises(FAILURES): owner.__enter__()

    def test_prefaulted_mmap_source_drift_after_child_fails_closed(self):
        with (self.root / "cc1plus").open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[128] = mapping[128]
            with self.assertRaises(FAILURES):
                with self.opened() as owner:
                    tool = owner._tools["cc1plus"]
                    fields = module.provenance._stable_fields(os.fstat(tool.fd))
                    original = os.pread(tool.executable_descriptor, 1, 128)
                    def mutate(*args, **kwargs):
                        mapping[128] ^= 1
                        self.assertEqual(module.provenance._stable_fields(os.fstat(tool.fd)), fields)
                        tool.guard.verify()
                        self.assertEqual(os.pread(tool.executable_descriptor, 1, 128), original)
                        return b""
                    with mock.patch.object(module.provenance, "_run", side_effect=mutate):
                        owner.run(["/usr/bin/c++"])

    def test_helper_mapping_replacement_is_detected(self):
        with self.assertRaises(FAILURES):
            with self.opened() as owner:
                owner.root.chmod(0o700)
                link = owner.root / "as"
                link.unlink()
                link.symlink_to("/usr/bin/true")
                owner.root.chmod(0o500)
                owner.record()

    def test_closed_helper_descriptor_is_rejected(self):
        with self.assertRaises(FAILURES):
            with self.opened() as owner:
                owner._tools["as"].close()
                owner.record()

    def test_streamed_tool_bound_is_strict_and_enforced(self):
        for bound in (True, 0, -1, 1.5, 65 << 20, 1):
            with self.subTest(bound=bound), self.assertRaises(FAILURES):
                self.factory(self.pins["cc"]["path"], maximum_bytes=bound)
        tool = self.factory(self.pins["cc"]["path"], maximum_bytes=self.pins["cc"]["size"])
        try:
            self.assertGreaterEqual(tool.executable_descriptor, 0)
            self.assertEqual(tool.executable_record()["sha256"], self.pins["cc"]["sha256"])
        finally: tool.close()

    def test_intermediate_driver_alias_replacement_and_restore_is_rejected(self):
        logical, pins = module.phase_inventory(self.pinned, "c++")
        alias = self.alias_bin / "c++"
        unchanged = rejected = False
        with mock.patch.object(module, "phase_inventory", return_value=(str(alias), pins)):
            with self.assertRaises(FAILURES):
                with self.opened() as owner:
                    middle = self.alias_alt / "c++"
                    target = os.readlink(middle)
                    middle.unlink()
                    middle.symlink_to("../cc")
                    middle.unlink()
                    middle.symlink_to(target)
                    self.assertEqual(alias.resolve(), owner._tools["driver"].path)
                    owner._tools["driver"].validate_current()
                    unchanged = True
                    with self.assertRaises(FAILURES): owner.record()
                    rejected = True
        self.assertTrue(unchanged)
        self.assertTrue(rejected)

    def test_alias_target_mismatch_loop_and_directory_alias_are_rejected(self):
        alias = self.alias_bin / "c++"
        for target in ("../cc", "c++", "../redirect/c++"):
            with self.subTest(target=target):
                alias.unlink()
                alias.symlink_to(target)
                if target.startswith("../redirect"):
                    (self.root / "redirect").symlink_to(self.alias_alt, target_is_directory=True)
                with self.assertRaises(FAILURES):
                    with self.opened(): self.fail("unsafe compiler alias accepted")

    def test_alias_cancelled_directory_component_and_long_chain_are_rejected(self):
        alias = self.alias_bin / "c++"
        alias.unlink()
        # normpath would hide the uninspected directory component here.
        alias.symlink_to("../alternatives/../c++")
        self.assertEqual(alias.resolve(), Path(self.pins["c++"]["path"]))
        with self.assertRaises(FAILURES):
            with self.opened(): self.fail("cancelled directory component accepted")
        alias.unlink()
        alias.symlink_to("hop0")
        for index in range(32):
            (self.alias_bin / ("hop"+str(index))).symlink_to("hop"+str(index+1) if index < 31 else "../c++")
        with self.assertRaises(FAILURES):
            with self.opened(): self.fail("excessive alias chain accepted")

    def test_alias_drift_is_rejected_before_child_launch(self):
        rejected = False
        with self.assertRaises(FAILURES):
            with self.opened() as owner, mock.patch.object(module.provenance, "_run") as run:
                middle = self.alias_alt / "c++"
                middle.unlink()
                middle.symlink_to("../cc")
                with self.assertRaises(FAILURES): owner.run(["/usr/bin/c++"])
                run.assert_not_called()
                rejected = True
        self.assertTrue(rejected)

    def test_alias_parent_permission_aba_and_inheritable_fd_are_rejected(self):
        for mutation in ("mode", "fd"):
            reached = False
            with self.subTest(mutation=mutation), self.assertRaises(FAILURES):
                with self.opened() as owner:
                    if mutation == "mode":
                        initial = self.alias_alt.stat().st_mode & 0o777
                        self.alias_alt.chmod(0o700)
                        self.alias_alt.chmod(initial)
                    else: os.set_inheritable(owner.aliases._directories[self.alias_alt][0], True)
                    reached = True
                    owner.validate_current()
            self.assertTrue(reached)

    def test_alias_record_is_detached_and_parent_fds_close(self):
        with self.opened() as owner:
            record = owner.record()
            self.assertTrue(record["driver_aliases"]["file_alias_chain_retained"])
            self.assertEqual(len(record["driver_aliases"]["nodes"]), 3)
            record["driver_aliases"]["nodes"].clear()
            self.assertEqual(len(owner.record()["driver_aliases"]["nodes"]), 3)
            descriptors = [row[0] for row in owner.aliases._directories.values()]
        for descriptor in descriptors:
            with self.assertRaises(OSError): os.fstat(descriptor)


if __name__ == "__main__":
    unittest.main()
