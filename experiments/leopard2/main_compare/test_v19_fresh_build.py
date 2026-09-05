#!/usr/bin/python3
"""Bounded build-owner tests; real files/fds, synthetic host and compiler output."""
from contextlib import contextmanager
import copy
import hashlib
import importlib.util
import json
import mmap
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tested_fresh_build", HERE / "v19_fresh_build.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FAILURES = (module.host.PreflightError, module.provenance.BuildProvenanceError, OSError, ValueError)


def compile_bytes(entries):
    return json.dumps([{key: value for key, value in row.items() if key != "arguments"} |
                       {"command": shlex.join(row["arguments"])} for row in entries.values()]).encode()


class FreshBuildTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preregistration = subprocess.run([str(HERE / "run_authoritative_v17_gfni_main_compare.sh"),
            "--print-conditioned-v19-preregistration"], check=True, stdout=subprocess.PIPE).stdout

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="leopard-v19-builder-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.parent = self.root / "new"
        self.parent.mkdir(mode=0o700)
        self.canonical = self.root / "historical"
        old_build = self.canonical / "candidate-build"
        old_build.mkdir(parents=True)
        self.pins = {key: module.host.load_preregistration(self.preregistration)["build_preflight"][key]
                     for key in ("source_commit", "source_tree", "baseline_commit", "baseline_tree")}
        row = {"directory": str(old_build), "file": str(self.canonical / "candidate-source/source.cpp"),
               "output": "CMakeFiles/leopard.dir/source.cpp.o", "arguments": ["/usr/bin/c++", "-O3",
               "-o", "CMakeFiles/leopard.dir/source.cpp.o", "-c", str(self.canonical / "candidate-source/source.cpp")]}
        commands = compile_bytes({row["output"]: row})
        (old_build / "compile_commands.json").write_bytes(commands)
        (old_build / "compile_commands.json").chmod(0o600)
        self.pinned = {"compile_commands": {"sha256": hashlib.sha256(commands).hexdigest()},
                       "benchmark_git": {"sha256": hashlib.sha256(Path("/usr/bin/git").read_bytes()).hexdigest()},
                       "validated_cache": {"CMAKE_BUILD_TYPE": "Release", "CMAKE_CXX_FLAGS": ""},
                       "archive_link_commands": [["/usr/bin/ar", "qc", "libleopard.a", row["output"]]],
                       "executable_link_command": ["/usr/bin/c++", "-o", "bench_leopard2", "libleopard.a"]}
        self.events = []
        self.live = True
        self.owner = None
        self.intercept = lambda stage: None
        test = self

        class Lease:
            def __init__(self, data): test.assertEqual(data, test.preregistration)
            def __enter__(self): test.events.append("lock"); return self
            def validate_current(self): module.require(test.live, "fixture host or lock failed")
            def freeze(self, lane, roots):
                self.validate_current()
                test.assertEqual(len(test.owner._commands), 10)
                test.assertFalse(lane.is_relative_to(test.owner.workspace))
                test.assertEqual(set(roots), {"candidate", "baseline"})
                test.events.append("freeze")
                return self.record()
            def record(self): return {"fixture": True}
            def __exit__(self, *args): test.events.append("unlock")

        class Retained:
            def __init__(self, data): test.assertEqual(data, test.preregistration)
            def __enter__(self): test.events.append("preflight"); return self
            def validate_current(self): module.require(test.live, "fixture preflight failed")
            def record(self): return {"source_pins": copy.deepcopy(test.pins), "canonical_run_root": str(test.canonical)}
            def _bytes(self, name):
                test.assertEqual(name, "candidate-build-provenance.json")
                return json.dumps(test.pinned).encode()
            def __exit__(self, *args): test.events.append("release-preflight")

        class Authenticated:
            def __init__(self, source, retained): self.source = source
            def __enter__(self): test.events.append("authenticate"); return self
            def validate_current(self, **kwargs): self.source.validate_current(**kwargs)
            def record(self, **kwargs): self.validate_current(**kwargs); return {"fixture": True}
            def __exit__(self, *args): self.validate_current(); test.events.append("release-source")

        self.Lease, self.Retained, self.Authenticated = Lease, Retained, Authenticated

    def cache(self, build, values):
        lines = []
        for key, value in values.items():
            types = module.provenance.CMAKE_CACHE_REQUIRED_ENTRY_TYPES.get(key)
            kind = sorted(types)[0] if types else "STRING"
            if build.name == "baseline-build" and key == "CMAKE_EXPORT_COMPILE_COMMANDS": kind = "BOOL"
            lines.append(key + ":" + kind + "=" + value)
        (build / "CMakeCache.txt").write_text("\n".join(lines) + "\n")

    def child(self, argv, label, **keywords):
        owner = self.owner
        stage = owner._commands[-1]["name"]
        self.assertEqual(keywords["environment_overrides"], module.ENVIRONMENT)
        self.assertEqual(keywords["timeout"], 600)
        self.assertEqual(keywords["maximum_bytes"], 1 << 20)
        self.assertGreaterEqual(keywords["executable_descriptor"], 0)
        observed = os.umask(0o077)
        os.umask(observed)
        self.assertEqual(observed, owner._commands[-1]["umask"])
        if stage.endswith(":clone"):
            root = Path(argv[-1])
            root.mkdir(exist_ok=True)
            (root / ".git").mkdir()
            (root / ".git/config").write_bytes(b"fixture")
            (root / "source.cpp").write_bytes(b"source")
            if root.name == "candidate-source": (root / "sse2neon").mkdir()
        elif stage.endswith("-configure"):
            role = stage.split("-")[0]
            build = owner.workspace / (role + "-build")
            if role == "candidate":
                entries = owner._candidate_entries
                values = self.pinned["validated_cache"]
                links = {"leopard": self.pinned["archive_link_commands"],
                         "bench_leopard2": [self.pinned["executable_link_command"]]}
            else:
                entries = module.expected_baseline_entries(owner.workspace, owner.canonical, owner.pins)
                values = {"CMAKE_CXX_FLAGS": owner._profile["baseline_prefix_map"],
                          "CMAKE_CXX_FLAGS_RELEASE": "-g -O0 -O3", "CMAKE_BUILD_TYPE": "Release",
                          "CMAKE_GENERATOR": "Unix Makefiles", "CMAKE_CXX_COMPILER": "/usr/bin/c++",
                          "LEOPARD_MAIN_SOURCE_DIR": str(owner.workspace / "leopard1-source"),
                          "LEO_MAIN_PURE_AVX2": "OFF", "CMAKE_EXPORT_COMPILE_COMMANDS": "ON"}
                links = {"leopard_main_exact": [["/usr/bin/ar", "qc", "libleopard_main_exact.a", *list(entries)[:-1]],
                                                 ["/usr/bin/ranlib", "libleopard_main_exact.a"]],
                         "leopard_main_benchmark": [["/usr/bin/c++", owner._profile["baseline_prefix_map"],
                           "-g", "-O0", "-O3", list(entries)[-1], "-o", "leopard_main_benchmark",
                           "libleopard_main_exact.a", "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so",
                           "/usr/lib/x86_64-linux-gnu/libpthread.a"]]}
            (build / "compile_commands.json").write_bytes(compile_bytes(entries))
            self.cache(build, values)
            for target, commands in links.items():
                directory = build / ("CMakeFiles/" + target + ".dir")
                directory.mkdir(parents=True)
                (directory / "link.txt").write_text("\n".join(shlex.join(command) for command in commands) + "\n")
            for relative in ("Makefile", "CMakeFiles/Makefile2", "CMakeFiles/Makefile.cmake"):
                (build / relative).write_bytes(b"fixture make scheduler\n")
            for target in {"/".join(output.split("/")[:2]) for output in entries}:
                for name in ("build.make", "flags.make"):
                    (build / target / name).write_bytes(b"fixture compile recipe\n")
        self.events.append(stage)
        self.intercept(stage)
        return b"fixture output\n"

    @contextmanager
    def opened(self):
        with mock.patch.object(module.provenance, "_run", new=self.child), \
             mock.patch.object(module.identity, "PinnedSourceIdentity", self.Authenticated):
            owner = module.FreshBuildOwner(self.preregistration, self.parent,
                                           _lease_factory=self.Lease, _preflight_factory=self.Retained)
            self.owner = owner
            with owner:
                yield owner

    def test_complete_order_lifetimes_false_claims_and_umask_restore(self):
        initial_umask = os.umask(0o077)
        os.umask(initial_umask)
        with self.opened() as owner:
            self.assertEqual(self.events, ["lock", "preflight"])
            record = owner.build()
            self.assertEqual(len(record["commands"]), 10)
            self.assertEqual(self.events[-1], "freeze")
            for key in ("live_acquisition_armed", "benchmark_executed", "runtime_closure_verified",
                        "physical_v18_lineage_verified", "atomic_snapshot"):
                self.assertIs(record[key], False)
            self.assertEqual(len(record["metadata"]), 21)
            record["recipe"]["stages"].clear()
            self.assertEqual(len(owner.record()["recipe"]["stages"]), 10)
            with self.assertRaises(FAILURES): owner.build()
        self.assertEqual(self.events[-3:], ["release-source", "release-preflight", "unlock"])
        with self.assertRaises(FAILURES): owner.record()
        self.assertTrue(owner.root.exists())
        observed = os.umask(initial_umask)
        self.assertEqual(observed, initial_umask)

    def test_exact_recipe_rejects_maps_jobs_stages_environment_and_old_schema(self):
        workspace = self.parent / "work"
        expected = module.recipe(workspace, self.canonical, self.pins)
        mutations = [lambda p: p.update(schema="leopard2-v18-relocated-build-recipe/v1"),
                     lambda p: p["stages"][8]["argv"].pop(),
                     lambda p: p["stages"][8]["argv"].append(p["stages"][8]["argv"][-1]),
                     lambda p: p["stages"][6]["argv"].append(p["stages"][8]["argv"][-1]),
                     lambda p: p["stages"][8]["argv"].__setitem__(-1, "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=/wrong=/old"),
                     lambda p: p["stages"][7]["argv"].__setitem__(4, "2"),
                     lambda p: p["stages"].reverse(), lambda p: p["environment"].update(LD_PRELOAD="x"),
                     lambda p: p["stages"][0].update(umask=False), lambda p: p.update(live_acquisition_armed=True)]
        for mutation in mutations:
            value = copy.deepcopy(expected)
            mutation(value)
            with self.subTest(value=value), self.assertRaises(FAILURES):
                module.validate_recipe(value, workspace, self.canonical, self.pins)

    def test_ambiguous_mapping_roots_rejected(self):
        for path in (self.canonical, self.canonical / "nested", self.root, self.parent / "x=y", self.parent / "a b"):
            with self.subTest(path=path), self.assertRaises(FAILURES): module.recipe(path, self.canonical, self.pins)

    def test_duplicate_compile_entries_and_json_keys_rejected(self):
        for data in (b'[{"command":"x","command":"y"}]', b'[]', b'{}',
                     b'[{"directory":"x","command":"x","file":"x","output":"x"},'
                     b'{"directory":"x","command":"x","file":"x","output":"x"}]'):
            with self.subTest(data=data), self.assertRaises(FAILURES): module.compile_entries(data)

    def test_baseline_bool_cache_dialect_does_not_relax_candidate_parser(self):
        data = b"CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON\n"
        self.assertEqual(module.baseline_cache(data), {"CMAKE_EXPORT_COMPILE_COMMANDS": "ON"})
        with self.assertRaises(FAILURES): module.provenance.parse_cmake_cache(data)
        with self.assertRaises(FAILURES): module.baseline_cache(data.replace(b"BOOL", b"UNINITIALIZED"))
        with self.assertRaises(FAILURES): module.baseline_cache(data + data)

    def test_compile_extra_map_rejected_before_build(self):
        def change(stage):
            if stage == "candidate-configure":
                path = self.owner.workspace / "candidate-build/compile_commands.json"
                rows = json.loads(path.read_bytes())
                rows[0]["command"] += " -ffile-prefix-map=/wrong=/other"
                path.write_text(json.dumps(rows))
        self.intercept = change
        with self.assertRaisesRegex(module.host.PreflightError, "compile argv"), self.opened() as owner: owner.build()
        self.assertNotIn("candidate-build", self.events)
        self.assertNotIn("freeze", self.events)

    def test_source_mmap_change_during_build_rejected(self):
        def change(stage):
            if stage == "candidate-build":
                with (self.owner.workspace / "candidate-source/source.cpp").open("r+b") as file:
                    with mmap.mmap(file.fileno(), 0) as view: view[0] = ord("X")
        self.intercept = change
        with self.assertRaises(FAILURES), self.opened() as owner: owner.build()
        self.assertIn("candidate-build", self.events)
        self.assertNotIn("baseline-configure", self.events)

    def test_metadata_write_restore_during_build_rejected(self):
        def change(stage):
            if stage == "candidate-build":
                path = self.owner.workspace / "candidate-build/compile_commands.json"
                data = path.read_bytes()
                path.write_bytes(data + b" ")
                path.write_bytes(data)
        self.intercept = change
        with self.assertRaises(FAILURES), self.opened() as owner: owner.build()
        self.assertIn("candidate-build", self.events)
        self.assertNotIn("freeze", self.events)

    def test_makefile_mmap_drift_rejected(self):
        def change(stage):
            if stage == "candidate-build":
                path = self.owner.workspace / "candidate-build/CMakeFiles/leopard.dir/build.make"
                with path.open("r+b") as file:
                    with mmap.mmap(file.fileno(), 0) as view: view[0] = ord("X")
        self.intercept = change
        with self.assertRaises(FAILURES), self.opened() as owner: owner.build()
        self.assertIn("candidate-build", self.events)
        self.assertNotIn("freeze", self.events)

    def test_clone_destination_symlink_rejected_before_child(self):
        outside = self.root / "outside"
        outside.mkdir()
        def change(stage):
            if stage == "candidate-source:checkout":
                submodule = self.owner.workspace / "candidate-source/sse2neon"
                submodule.rmdir()
                submodule.symlink_to(outside, target_is_directory=True)
        self.intercept = change
        with self.assertRaisesRegex(module.host.PreflightError, "clone destination"), self.opened() as owner:
            owner.build()
        self.assertNotIn("candidate-source/sse2neon:clone", self.events)
        self.assertEqual(list(outside.iterdir()), [])

    def test_replaced_candidate_parent_cannot_redirect_submodule_clone(self):
        outside = self.root / "outside"
        (outside / "sse2neon").mkdir(parents=True)
        def change(stage):
            if stage == "candidate-source:checkout":
                candidate = self.owner.workspace / "candidate-source"
                candidate.rename(self.owner.workspace / "replaced-candidate")
                candidate.symlink_to(outside, target_is_directory=True)
        self.intercept = change
        with self.assertRaises(FAILURES), self.opened() as owner: owner.build()
        self.assertNotIn("candidate-source/sse2neon:clone", self.events)
        self.assertEqual(list((outside / "sse2neon").iterdir()), [])

    def test_lock_or_resource_loss_during_child_prevents_next_stage(self):
        def change(stage): self.live = False
        self.intercept = change
        with self.assertRaises(FAILURES), self.opened() as owner: owner.build()
        self.assertIn("leopard1-source:clone", self.events)
        self.assertNotIn("leopard1-source:checkout", self.events)

    def test_wrong_link_recipe_rejected_before_build(self):
        def change(stage):
            if stage == "candidate-configure":
                path = self.owner.workspace / "candidate-build/CMakeFiles/leopard.dir/link.txt"
                path.write_text(path.read_text() + "/usr/bin/true\n")
        self.intercept = change
        with self.assertRaisesRegex(module.host.PreflightError, "link argv"), self.opened() as owner: owner.build()
        self.assertNotIn("candidate-build", self.events)

    def test_effective_cache_mapping_mismatch_rejected_before_build(self):
        def change(stage):
            if stage == "baseline-configure":
                path = self.owner.workspace / "baseline-build/CMakeCache.txt"
                path.write_text(path.read_text().replace("-ffile-prefix-map=", "-fdebug-prefix-map="))
        self.intercept = change
        with self.assertRaisesRegex(module.host.PreflightError, "effective CMake cache"), self.opened() as owner: owner.build()
        self.assertNotIn("baseline-build", self.events)

    def test_child_failure_preserved_and_not_retried(self):
        def change(stage): raise module.host.PreflightError("injected child failure")
        self.intercept = change
        with self.assertRaisesRegex(module.host.PreflightError, "injected child failure"), self.opened() as owner:
            try: owner.build()
            except module.host.PreflightError:
                record = owner.failure_record()
                self.assertEqual(len(record["commands"]), 1)
                self.assertEqual(record["commands"][0]["status"], "failed")
                self.assertTrue(owner.root.exists())
                with self.assertRaises(FAILURES): owner.build()
                # Propagate the original failure through owner exit.
                raise

    def test_host_failure_before_creation(self):
        self.live = False
        with self.assertRaises(FAILURES), self.opened(): pass
        self.assertEqual(list(self.parent.iterdir()), [])
        self.assertNotIn("leopard1-source:clone", self.events)

    def test_historical_metadata_hash_mismatch_stops_before_clone(self):
        self.pinned["compile_commands"]["sha256"] = "0" * 64
        with self.assertRaises(FAILURES), self.opened(): pass
        self.assertNotIn("leopard1-source:clone", self.events)

    def test_symlink_parent_rejected(self):
        actual = self.root / "actual"
        self.parent.rename(actual)
        self.parent.symlink_to(actual, target_is_directory=True)
        with self.assertRaises(FAILURES), self.opened(): pass
        self.assertEqual(list(actual.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
