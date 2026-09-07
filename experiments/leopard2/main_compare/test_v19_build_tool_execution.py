#!/usr/bin/python3
"""Bounded retained build-root/runtime tests; no external build jobs."""
from contextlib import ExitStack
import copy
import hashlib
import importlib.util
import json
import mmap
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("build_tool_runtime_fixtures", HERE / "test_v19_runtime_inventory.py")
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)
runtime, module = fixtures.module, fixtures.module.build_tools
FAILURES = fixtures.FAILURES


class BuildToolTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        parent = "/dev/shm" if "mmap" in self._testMethodName else None
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory(prefix="leopard-v19-build-tools-test-", dir=parent)))
        def permissions():
            for directory, _, _ in os.walk(self.root): os.chmod(directory, 0o700)
        self.stack.callback(permissions)
        self.parent, self.libs = self.root / "new", self.root / "libs"
        self.parent.mkdir(mode=0o700); self.libs.mkdir(mode=0o700)
        self.pinned, self.observed, roles, paths = {}, {}, {}, {}
        for role, (key, _) in module.ROLES.items():
            path = self.root / (role + "-tool")
            path.write_bytes(fixtures.elf(("libsample.so.1",), interpreter=True)); path.chmod(0o755)
            logical = self.root / role
            logical.symlink_to(path.name)
            row = {"path": str(path), "size": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            if key is None: self.observed[role], paths[role] = row, str(path)
            else: self.pinned[key] = dict(row, uid=0, gid=0, mode=0o100755)
            roles[role] = (key, str(logical))
        self.stack.enter_context(mock.patch.object(module, "ROLES", roles))
        self.stack.enter_context(mock.patch.object(module, "OBSERVED_PATHS", paths))
        self.loader = self.libs / "ld-linux-x86-64.so.2"
        self.loader.write_bytes(fixtures.elf(soname=self.loader.name)); self.loader.chmod(0o755)
        library = self.libs / "libsample.so.1"
        library.write_bytes(fixtures.elf(soname=library.name)); library.chmod(0o644)
        self.live, self.calls = True, []
        test = self
        class Retained:
            def validate_current(self): module.require(test.live, "fixture retained preflight lost")
            def _bytes(self, name):
                module.require(name == "candidate-build-provenance.json", "unexpected preflight input")
                return json.dumps(test.pinned).encode()
        self.retained = Retained()

    def factory(self, path, **kwargs):
        return module.builder._StreamedTool(path, _trusted_owner=(os.getuid(), os.getgid()), **kwargs)

    def cleanup_owner(self, owner):
        try: owner.__exit__(RuntimeError, RuntimeError("fixture teardown"), None)
        except FAILURES: pass

    def phase(self):
        owner = module.BuildToolExecution(self.retained, self.parent, self.observed, _tool_factory=self.factory).__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        return owner

    def inventory(self, phase=None):
        phase = self.phase() if phase is None else phase
        owner = runtime.RuntimeInventory(phase, _tool_factory=self.factory, _library_root=self.libs,
            _loader=self.loader, _preload=self.root / "preload").__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        return owner

    def run_child(self, argv, label, **kwargs):
        self.calls.append((argv, kwargs))
        return b"fixture output\n"

    def policy(self, phase):
        path = self.root / "gnutls-config"
        if not path.exists():
            path.write_bytes(b"[overrides]\ndisabled-version = tls1.0\n"); path.chmod(0o644)
        self.stack.enter_context(mock.patch.object(module, "PRIORITY_PATH", str(path)))
        pin = {"path": str(path), "size": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        owner = module.CMakePriorityConfiguration(phase, pin, _file_factory=self.factory).__enter__()
        self.stack.callback(self.cleanup_owner, owner)
        return owner

    def test_exact_roots_origin_aliases_and_seals(self):
        phase = self.phase()
        record = phase.record()
        self.assertEqual(set(record["tools"]), set(module.ROLES))
        self.assertEqual(set(record["original_preflight_roles"]), {"git", "make", "ar", "ranlib"})
        self.assertEqual(set(record["newly_observed_roles"]), {"cmake", "shell"})
        for role, tool in phase._tools.items():
            self.assertEqual(record["tools"][role]["seals"], 15)
            self.assertEqual(record["aliases"][role]["resolved_path"], str(tool.path))
            with self.assertRaises(OSError): os.pwrite(tool.executable_descriptor, b"!", 0)
        for flag in ("nested_tool_execution_owned", "build_tool_data_owned", "full_runtime_execution_owned",
                     "fresh_build_recipe_integrated", "live_acquisition_armed", "benchmark_executed"):
            self.assertFalse(record[flag])

    def test_inventory_rejects_missing_extra_mistyped_and_unpinned_roles(self):
        for observed in ({}, dict(self.observed, extra=self.observed["shell"])):
            with self.assertRaises(FAILURES): module.tool_inventory(self.pinned, observed)
        for key, value in (("size", True), ("size", 0), ("size", module.MAX_TOOL_BYTES + 1),
                           ("sha256", "bad"), ("path", "/unqualified/cmake")):
            observed = copy.deepcopy(self.observed); observed["cmake"][key] = value
            with self.subTest(key=key), self.assertRaises(FAILURES): module.tool_inventory(self.pinned, observed)
        for key in ("benchmark_git", "archiver", "ranlib", "make_program"):
            pinned = copy.deepcopy(self.pinned); del pinned[key]
            with self.assertRaises(FAILURES): module.tool_inventory(pinned, self.observed)
        with mock.patch.object(module, "MAX_PHASE_BYTES", 1), self.assertRaises(FAILURES):
            module.tool_inventory(self.pinned, self.observed)

    def test_wrong_hash_and_alias_endpoint_refuse_before_use(self):
        self.observed["cmake"]["sha256"] = "0" * 64
        with self.assertRaises(FAILURES): self.phase()
        self.observed["cmake"]["sha256"] = hashlib.sha256(Path(self.observed["cmake"]["path"]).read_bytes()).hexdigest()
        path = Path(module.ROLES["ar"][1]); path.unlink(); path.symlink_to(Path(self.observed["cmake"]["path"]).name)
        with self.assertRaises(FAILURES): self.phase()

    def test_alias_change_restore_is_latched(self):
        phase = self.phase()
        path = Path(phase.logical_paths["shell"]); target = os.readlink(path)
        path.unlink(); path.symlink_to("cmake-tool"); path.unlink(); path.symlink_to(target)
        with self.assertRaises(FAILURES): phase.record()

    def test_preexisting_mmap_is_rehashed(self):
        path = Path(self.observed["cmake"]["path"])
        with path.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[0:1] = mapping[0:1]
            phase = self.phase()
            with mock.patch.object(phase._tools["cmake"], "_hash", wraps=phase._tools["cmake"]._hash) as hashed:
                mapping[0] ^= 1
                try:
                    with self.assertRaises(FAILURES): phase.validate_current()
                    self.assertTrue(hashed.called)
                finally: mapping[0] ^= 1

    def test_parent_permission_restore_and_retained_loss_latch(self):
        phase = self.phase()
        self.parent.chmod(0o755); self.parent.chmod(0o700)
        with self.assertRaises(FAILURES): phase.record()
        phase = self.phase(); self.live = False
        with self.assertRaises(FAILURES): phase.record()

    def test_argument_role_response_and_type_rejections(self):
        for role, argv in (("other", []), ("ar", [module.ROLES["ar"][1], "@response"]),
                          ("git", [module.ROLES["make"][1]]), ("cmake", [module.ROLES["cmake"][1], True])):
            phase = self.phase()
            with self.subTest(role=role), self.assertRaises(FAILURES): phase.arguments(role, argv)
            with self.assertRaises(FAILURES): phase.record()

    def test_closure_retains_build_roots_without_gcc_plugin(self):
        inventory = self.inventory()
        record = inventory.record()
        self.assertEqual(record["schema"], "leopard2-v19-runtime-inventory/v2")
        self.assertEqual(record["root_profile"], "build-tools")
        self.assertEqual(record["commands"], [])
        self.assertNotIn("gcc-link-plugin", record["roots"])
        self.assertEqual(set(record["roots"]), set(module.ROLES))

    def test_root_job_routes_loader_and_inherits_only_sealed_files(self):
        inventory = self.inventory()
        for role in module.ROLES:
            argv = [inventory.phase.logical_paths[role], "--version"]
            self.assertEqual(inventory.run_tool(role, argv, _runner=self.run_child), b"fixture output\n")
            actual, kwargs = self.calls[-1]
            self.assertEqual(actual[:9], [str(self.loader), "--inhibit-cache", "--glibc-hwcaps-mask", "", "--library-path",
                f"/proc/self/fd/{inventory.prefix.descriptor}", "--argv0", argv[0],
                f"/proc/self/fd/{inventory.phase._tools[role].executable_descriptor}"])
            self.assertEqual(kwargs["executable_descriptor"], inventory._files[inventory.loader_key].executable_descriptor)
            for tool in [*inventory.phase._tools.values(), *inventory._files.values()]:
                self.assertIn(tool.executable_descriptor, kwargs["inherited_descriptors"])
                self.assertNotIn(tool.fd, kwargs["inherited_descriptors"])
            record = inventory.record()["commands"][-1]
            self.assertEqual(record["logical_argv"], argv)
            self.assertEqual(record["status"], "exit-zero")

    def test_unsealed_or_mistyped_input_descriptors_refuse_before_child(self):
        inventory = self.inventory()
        argv = [inventory.phase.logical_paths["ar"], "t", "archive.a"]
        for fds in ((True,), (1,), [5], (inventory.phase._tools["ar"].fd,)):
            owner = self.inventory()
            with self.subTest(fds=fds), self.assertRaises(FAILURES):
                owner.run_tool("ar", argv, input_descriptors=fds, _runner=self.run_child)
        self.assertEqual(self.calls, [])

    def test_input_inheritance_loss_during_job_latches(self):
        inventory = self.inventory()
        fd = inventory.phase._tools["git"].executable_descriptor
        def mutate(*args, **kwargs):
            os.set_inheritable(fd, True)
            return b""
        try:
            with self.assertRaises(FAILURES): inventory.run_tool("ar", [inventory.phase.logical_paths["ar"]],
                input_descriptors=(fd,), _runner=mutate)
            self.assertEqual(inventory._commands[-1]["status"], "failed")
        finally: os.set_inheritable(fd, False)

    def test_seal_query_works_when_python_omits_linux_constant(self):
        inventory = self.inventory()
        fd = inventory.phase._tools["git"].executable_descriptor
        with mock.patch.dict(module.builder.fcntl.__dict__):
            module.builder.fcntl.__dict__.pop("F_GET_SEALS", None)
            output = inventory.run_tool("ar", [inventory.phase.logical_paths["ar"]],
                input_descriptors=(fd,), _runner=self.run_child)
        self.assertEqual(output, b"fixture output\n")
        self.assertEqual(inventory.record()["commands"][-1]["status"], "exit-zero")

    def test_phase_loss_during_job_latches_and_records_failure(self):
        inventory = self.inventory()
        def mutate(*args, **kwargs):
            self.live = False
            return b""
        with self.assertRaises(FAILURES): inventory.run_tool("make", [inventory.phase.logical_paths["make"]], _runner=mutate)
        self.assertEqual(inventory._commands[-1]["status"], "failed")
        with self.assertRaises(FAILURES): inventory.record()

    def test_detached_records_and_phase_closure(self):
        phase = self.phase()
        record = phase.record(); record["tools"].clear()
        self.assertNotEqual(record, phase.record())
        fds = [phase._parent_fd, *[fd for tool in phase._tools.values() for fd in (tool.fd, tool.executable_descriptor)]]
        phase.__exit__(None, None, None)
        for fd in fds:
            with self.assertRaises(OSError): os.fstat(fd)
        with self.assertRaises(FAILURES): phase.__enter__()

    def test_policy_preserves_bytes_and_records_new_pin_authority(self):
        phase = self.phase(); policy = self.policy(phase)
        record = policy.record()
        self.assertEqual(os.pread(policy.descriptor(), record["size"], 0), Path(record["path"]).read_bytes())
        self.assertEqual(record["seals"], 15)
        self.assertEqual(record["source_mode"], 0o100644)
        self.assertTrue(record["declared_policy_bytes_sealed"])
        for name in ("original_preflight_pin", "nested_configuration_dependencies_owned", "full_runtime_execution_owned"):
            self.assertFalse(record[name])
        record["environment"].clear()
        self.assertNotEqual(record, policy.record())
        with self.assertRaises(OSError): os.pwrite(policy.descriptor(), b"!", 0)

    def test_policy_requires_exact_pin_and_build_phase(self):
        phase = self.phase(); policy = self.policy(phase)
        for name, value in (("path", "/another/config"), ("size", True), ("size", 0),
                            ("size", module.MAX_PRIORITY_BYTES + 1), ("sha256", "bad")):
            pin = dict(policy.pin); pin[name] = value
            with self.subTest(name=name), self.assertRaises(FAILURES): module.CMakePriorityConfiguration(phase, pin)
        with self.assertRaises(FAILURES): module.CMakePriorityConfiguration(object(), policy.pin)
        with self.assertRaises(FAILURES): module.CMakePriorityConfiguration(phase, dict(policy.pin, extra=1))

    def test_policy_wrong_hash_and_symlink_refuse(self):
        phase = self.phase(); policy = self.policy(phase)
        pin = dict(policy.pin, sha256="0" * 64)
        with self.assertRaises(FAILURES):
            module.CMakePriorityConfiguration(phase, pin, _file_factory=self.factory).__enter__()
        path = Path(policy.pin["path"]); data = path.read_bytes()
        target = self.root / "other-config"; target.write_bytes(data)
        path.unlink(); path.symlink_to(target.name)
        with self.assertRaises(FAILURES):
            module.CMakePriorityConfiguration(phase, policy.pin, _file_factory=self.factory).__enter__()

    def test_policy_mmap_mutation_is_rehashed_and_failure_latched(self):
        phase = self.phase()
        path = self.root / "gnutls-config"
        path.write_bytes(b"[overrides]\ndisabled-version = tls1.0\n"); path.chmod(0o644)
        with path.open("r+b") as stream, mmap.mmap(stream.fileno(), 0) as mapping:
            mapping[0:1] = mapping[0:1]
            policy = self.policy(phase)
            with mock.patch.object(policy._file, "_hash", wraps=policy._file._hash) as hashed:
                mapping[0] ^= 1
                try:
                    with self.assertRaises(FAILURES): policy.validate_current()
                    self.assertTrue(hashed.called)
                finally: mapping[0] ^= 1
            with self.assertRaises(FAILURES): policy.record()

    def test_policy_role_phase_and_lifetime_refuse_before_execution(self):
        inventory = self.inventory(); policy = self.policy(inventory.phase)
        other = self.inventory()
        cases = [(other, "cmake", policy), (inventory, "ar", policy), (self.inventory(), "cmake", object())]
        for owner, role, selected in cases:
            with self.subTest(role=role), self.assertRaises(FAILURES):
                owner.run_tool(role, [owner.phase.logical_paths[role]], configuration=selected, _runner=self.run_child)
        self.assertEqual(self.calls, [])

    def test_policy_routes_sealed_descriptor_only_and_keeps_default_environment(self):
        inventory = self.inventory(); policy = self.policy(inventory.phase)
        argv = [inventory.phase.logical_paths["cmake"], "--version"]
        baseline = dict(module.builder.ENVIRONMENT)
        inventory.run_tool("cmake", argv, configuration=policy, _runner=self.run_child)
        _, kwargs = self.calls[-1]
        self.assertEqual(kwargs["environment_overrides"], dict(baseline, **policy.record()["environment"]))
        self.assertEqual(module.builder.ENVIRONMENT, baseline)
        self.assertIn(policy.descriptor(), kwargs["inherited_descriptors"])
        self.assertNotIn(policy._file.fd, kwargs["inherited_descriptors"])
        self.assertEqual(inventory.record()["commands"][-1]["configuration"], policy.record())
        inventory.run_tool("cmake", argv, _runner=self.run_child)
        self.assertNotIn("configuration", inventory.record()["commands"][-1])
        self.assertEqual(self.calls[-1][1]["environment_overrides"], baseline)

    def test_policy_mutation_during_job_records_failure(self):
        inventory = self.inventory(); policy = self.policy(inventory.phase)
        def mutate(*args, **kwargs):
            path = Path(policy.pin["path"])
            path.chmod(0o755); path.chmod(0o644)
            return b""
        with self.assertRaises(FAILURES):
            inventory.run_tool("cmake", [inventory.phase.logical_paths["cmake"]], configuration=policy, _runner=mutate)
        self.assertEqual(inventory._commands[-1]["status"], "failed")
        with self.assertRaises(FAILURES): inventory.record()

    def test_policy_descriptor_and_phase_loss_latch(self):
        phase = self.phase(); policy = self.policy(phase)
        fd = policy.descriptor()
        os.set_inheritable(fd, True)
        try:
            with self.assertRaises(FAILURES): policy.validate_current()
        finally: os.set_inheritable(fd, False)
        with self.assertRaises(FAILURES): policy.record()
        policy = self.policy(phase); self.live = False
        with self.assertRaises(FAILURES): policy.record()

    def test_policy_close_releases_descriptors_and_refuses_job(self):
        inventory = self.inventory(); policy = self.policy(inventory.phase)
        fds = (policy._file.fd, policy.descriptor())
        policy.__exit__(None, None, None)
        for fd in fds:
            with self.assertRaises(OSError): os.fstat(fd)
        with self.assertRaises(FAILURES): policy.__enter__()
        with self.assertRaises(FAILURES):
            inventory.run_tool("cmake", [inventory.phase.logical_paths["cmake"]], configuration=policy, _runner=self.run_child)
        self.assertEqual(self.calls, [])

if __name__ == "__main__": unittest.main()
