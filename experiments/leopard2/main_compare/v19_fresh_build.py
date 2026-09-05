#!/usr/bin/python3
"""Dormant owned fresh builds; leopard-79h.38.5.4.8.2.2.2.1.

This composes source authentication and artifact ownership, not acquisition.
No workload, qualification, arming, or historical-validator dispatch exists.
Compiler commands are checked against the pinned recipe. Runtime/toolchain
closure and physical v18 lineage remain separate required parent gates.
"""
from __future__ import annotations

from contextlib import ExitStack
import copy
import ctypes
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import stat
import tempfile
import threading

HERE = Path(__file__).resolve().parent
_dependency = HERE / "v19_source_identity.py"
if _dependency.resolve(strict=True) != _dependency:
    raise RuntimeError("fresh build dependency is not canonical")
_spec = importlib.util.spec_from_file_location("v19_build_identity", _dependency)
identity = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(identity)
streamed, preflight, owners = identity.streamed, identity.preflight, identity.owners
host, provenance, require = identity.host, identity.provenance, identity.require

SCHEMA = "leopard2-v19-relocated-build-recipe/v1"
MAX_METADATA_BYTES = 1 << 20
ENVIRONMENT = {"TZ": "UTC", "CMAKE_BUILD_PARALLEL_LEVEL": "1", "MAKEFLAGS": "",
               "PYTHONDONTWRITEBYTECODE": "1"}
BASELINE_FILES = ("leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp")


def canonical_root(value):
    host.canonical_path(str(value), nonroot=True)
    require(re.fullmatch(r"/[A-Za-z0-9_./-]+", str(value)) is not None,
            "build path cannot be represented by the exact relocation recipe")
    return Path(value)


def recipe(workspace: Path, canonical: Path, pins: dict) -> dict:
    """The sole supported relocation: one baseline-only debug/file prefix map."""
    workspace, canonical = canonical_root(workspace), canonical_root(canonical)
    require(workspace != canonical and not workspace.is_relative_to(canonical) and
            not canonical.is_relative_to(workspace), "build workspace overlaps historical sources")
    candidate, baseline = workspace / "candidate-source", workspace / "leopard1-source"
    candidate_build, baseline_build = workspace / "candidate-build", workspace / "baseline-build"
    git = ["/usr/bin/git", "-c", "core.hooksPath=/dev/null"]
    stages = []
    for role, commit in (("leopard1-source", pins["baseline_commit"]),
                         ("candidate-source", pins["source_commit"]),
                         ("candidate-source/sse2neon", identity.SUBMODULE_COMMIT)):
        require(type(commit) is str and identity.HEX40.fullmatch(commit), "invalid pinned source commit")
        stages.extend([
            {"name": role + ":clone", "umask": 0o002,
             "argv": git + ["clone", "--no-local", "--no-checkout", "--single-branch", "--depth", "1",
                            str(canonical / role), str(workspace / role)]},
            {"name": role + ":checkout", "umask": 0o002,
             "argv": git + ["-C", str(workspace / role), "checkout", "--detach", commit]},
        ])
    mapping = "-ffile-prefix-map=" + str(workspace) + "=" + str(canonical)
    stages.extend([
        {"name": "candidate-configure", "umask": 0o022,
         "argv": ["/usr/bin/cmake", "-S", str(candidate), "-B", str(candidate_build),
                  "-G", "Unix Makefiles", "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
                  "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=ON", "-DLEO2_BUILD_FUZZERS=OFF",
                  "-DLEO2_ENABLE_CUDA=OFF", "-DLEO2_BACKEND_VARIANT=auto", "-DLEOPARD_ENABLE_GF8=ON",
                  "-DLEOPARD_ENABLE_GF16=ON", "-DENABLE_OPENMP=ON", "-DLEO2_FLAG_MAVX512F=FALSE",
                  "-DLEO2_FLAG_MAVX512BW=FALSE", "-DLEO2_FLAG_MAVX512VL=FALSE",
                  "-DLEO2_FLAG_MPREFER_VECTOR_WIDTH_256=FALSE", "-DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git"]},
        {"name": "candidate-build", "umask": 0o022,
         "argv": ["/usr/bin/cmake", "--build", str(candidate_build), "--parallel", "1",
                  "--target", "bench_leopard2"]},
        {"name": "baseline-configure", "umask": 0o022,
         "argv": ["/usr/bin/cmake", "-S", str(candidate / "experiments/leopard2/main_compare"),
                  "-B", str(baseline_build), "-G", "Unix Makefiles", "-DCMAKE_BUILD_TYPE=Release",
                  "-DLEOPARD_MAIN_SOURCE_DIR=" + str(baseline), "-DLEO_MAIN_PURE_AVX2=OFF",
                  "-DCMAKE_CXX_FLAGS=" + mapping]},
        {"name": "baseline-build", "umask": 0o022,
         "argv": ["/usr/bin/cmake", "--build", str(baseline_build), "--parallel", "1",
                  "--target", "leopard_main_benchmark"]},
    ])
    return {"schema": SCHEMA, "workspace": str(workspace), "historical_root": str(canonical),
            "source_pins": copy.deepcopy(pins), "stages": stages, "environment": dict(ENVIRONMENT),
            "baseline_prefix_map": mapping, "live_acquisition_armed": False}


def validate_recipe(value: dict, workspace: Path, canonical: Path, pins: dict):
    require(type(value) is dict and host.canonical_bytes(value) ==
            host.canonical_bytes(recipe(workspace, canonical, pins)),
            "v19 relocated recipe differs from the exact stage/argv profile")


def compile_entries(data: bytes) -> dict:
    require(type(data) is bytes and 0 < len(data) <= MAX_METADATA_BYTES, "compile metadata exceeds bound")
    rows = json.loads(data, object_pairs_hook=preflight._unique_object)
    require(type(rows) is list and 0 < len(rows) <= 64, "compile inventory is empty or oversized")
    result = {}
    for row in rows:
        require(type(row) is dict and set(row) == {"directory", "command", "file", "output"} and
                all(type(value) is str for value in row.values()), "unexpected compile metadata representation")
        for value in row.values():
            provenance._require_safe_unicode(value, "compile metadata")
        output = row["output"]
        require(output not in result, "duplicate compile output")
        result[output] = {"directory": row["directory"], "file": row["file"],
                          "output": output, "arguments": shlex.split(row["command"])}
    return result


def baseline_cache(data: bytes) -> dict:
    """Parse this adapter's BOOL export key without relaxing candidate/old profiles."""
    require(type(data) is bytes and 0 < len(data) <= MAX_METADATA_BYTES, "baseline cache exceeds bound")
    text = data.decode("utf-8", errors="strict")
    require("\0" not in text and "\r" not in text, "baseline cache contains a forbidden delimiter")
    result = {}
    for line in text.split("\n"):
        provenance._require_safe_unicode(line, "baseline CMake cache")
        if not line or line.startswith(("#", "//")):
            continue
        require("=" in line, "baseline cache record is unframed")
        typed, value = line.split("=", 1)
        require(typed.count(":") == 1, "baseline cache typed key is malformed")
        key, kind = typed.split(":", 1)
        require(re.fullmatch(r"[A-Za-z0-9_.+-]+", key) is not None and key not in result and
                kind in provenance.CMAKE_CACHE_ENTRY_TYPES, "baseline cache key is invalid or duplicated")
        allowed = ({"BOOL"} if key == "CMAKE_EXPORT_COMPILE_COMMANDS" else
                   provenance.CMAKE_CACHE_REQUIRED_ENTRY_TYPES.get(key))
        require(allowed is None or kind in allowed, "baseline cache key has an unexpected type: " + key)
        result[key] = value
    return result


def expected_baseline_entries(workspace: Path, canonical: Path, pins: dict) -> dict:
    profile = recipe(workspace, canonical, pins)
    baseline = workspace / "leopard1-source"
    result = {}
    for filename in (*BASELINE_FILES, "legacy_main_benchmark.cpp"):
        adapter = filename == "legacy_main_benchmark.cpp"
        source = (workspace / "candidate-source/experiments/leopard2/main_compare" if adapter else baseline) / filename
        output = ("CMakeFiles/leopard_main_benchmark.dir/legacy_main_benchmark.cpp.o" if adapter else
                  "CMakeFiles/leopard_main_exact.dir/" + str(baseline).lstrip("/") + "/" + filename + ".o")
        defines = ([f'-DLEOPARD_MAIN_SOURCE_COMMIT="{pins["baseline_commit"]}"',
                    "-DLEO_MAIN_PURE_AVX2_PROFILE=0"] if adapter else [])
        argv = ["/usr/bin/c++", *defines, "-I" + str(baseline), profile["baseline_prefix_map"],
                "-g", "-O0", "-O3", "-std=gnu++11", "-march=native", "-Wall", "-Wextra", "-fopenmp",
                "-o", output, "-c", str(source)]
        result[output] = {"directory": str(workspace / "baseline-build"), "file": str(source),
                          "output": output, "arguments": argv}
    return result


def relocated(value, canonical: Path, workspace: Path):
    if type(value) is str:
        return value.replace(str(canonical) + "/", str(workspace) + "/")
    if type(value) is list:
        return [relocated(item, canonical, workspace) for item in value]
    if type(value) is dict:
        return {key: relocated(item, canonical, workspace) for key, item in value.items()}
    return value


class _StreamedTool:
    """Guard a tool fd without caching another 12 MiB CMake image in Python.

    This is before/after byte observation, not immutable execution or complete
    loader/subtool ownership. The final orchestrator must establish that gate.
    """
    def __init__(self, path):
        self.path = Path(path)
        self.stack = ExitStack()
        try:
            require(self.path.resolve(strict=True) == self.path, "build launcher is not canonical")
            self.guard = self.stack.enter_context(provenance._InotifyMutationGuard("v19 build tool"))
            self.guard.add_file_path(self.path)
            self.fd = os.open(self.path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_NONBLOCK)
            self.stack.callback(os.close, self.fd)
            value = os.fstat(self.fd)
            require(stat.S_ISREG(value.st_mode) and value.st_uid == 0 and value.st_gid == 0 and
                    value.st_nlink >= 1 and stat.S_IMODE(value.st_mode) == 0o755 and
                    0 < value.st_size <= 16 << 20, "unsafe build launcher")
            # Distribution installations may hardlink packaged launchers.
            # This is not a task-owned source file:
            # retain the actual link count, guard the inode, and rehash bytes.
            self.fields = provenance._stable_fields(value)
            self.sha256 = self._hash(value.st_size)
            self.validate_current()
        except BaseException:
            self.stack.close()
            raise

    def _hash(self, size):
        digest = hashlib.sha256()
        for offset in range(0, size, streamed.HASH_BLOCK_BYTES):
            part = os.pread(self.fd, min(streamed.HASH_BLOCK_BYTES, size - offset), offset)
            require(len(part) == min(streamed.HASH_BLOCK_BYTES, size - offset), "build launcher truncated")
            digest.update(part)
        return digest.hexdigest()

    def validate_current(self):
        self.guard.verify()
        value = os.fstat(self.fd)
        require(provenance._stable_fields(value) == self.fields == provenance._stable_fields(self.path.lstat())
                and not os.get_inheritable(self.fd) and self._hash(value.st_size) == self.sha256,
                "build launcher changed")
        self.guard.verify()

    def close(self):
        self.stack.close()


class FreshBuildOwner:
    """Create, build, and freeze once while retaining all borrowed lifetimes.

    The private parent must already exist; all new outputs are exclusive and
    retained on failure. Neither the parent nor historical sources are removed.
    The caller must keep this context open through later handoff and sealing.
    Test injections are private; there is no CLI or armed wrapper integration.
    """
    def __init__(self, preregistration_bytes: bytes, parent: Path, *, _lease_factory=None,
                 _preflight_factory=None):
        self.preregistration_bytes = preregistration_bytes
        self.contract = host.load_preregistration(preregistration_bytes)["build_preflight"]
        self.parent = canonical_root(parent)
        self._lease_factory = owners.BuildArtifactLease if _lease_factory is None else _lease_factory
        self._preflight_factory = preflight.PinnedPreflight if _preflight_factory is None else _preflight_factory
        self._stack = ExitStack()
        self._state = "new"
        self._pid = os.getpid()
        self._directories = {}
        self._metadata = {}
        self._tools = {}
        self._commands = []
        self.authenticated = None
        self.root = None
        self._profile = None

    def _directory(self, path, *, private=True):
        descriptor = host.LinuxReader.open_directory(str(path))
        self._stack.callback(os.close, descriptor)
        value = os.fstat(descriptor)
        require((stat.S_IMODE(value.st_mode) == 0o700 if private else not value.st_mode & 0o7002)
                and value.st_uid == os.getuid() and
                value.st_gid == os.getgid(), "fresh build directory is not private")
        # Drain under the old watched-name set before adding this newly created
        # child. Otherwise its queued creation is retrospectively classified as
        # a mutation of an already-held directory. Existing names stay guarded.
        self._guard.verify()
        self._guard.add_directory_path(path)
        self._directories[path] = descriptor, streamed._directory_identity(value)

    def __enter__(self):
        require(self._state == "new", "fresh build owner cannot be reused")
        self._state = "entering"
        try:
            self.lease = self._stack.enter_context(self._lease_factory(self.preregistration_bytes))
            self.retained = self._stack.enter_context(self._preflight_factory(self.preregistration_bytes))
            self.lease.validate_current()
            self.retained.validate_current()
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 build directories"))
            self._directory(self.parent)
            historical = self.retained.record()
            self.canonical = canonical_root(historical["canonical_run_root"])
            self.pins = historical["source_pins"]
            self.root = Path(tempfile.mkdtemp(prefix="v19-build-", dir=self.parent))
            self._directory(self.root)
            self.workspace, self.lane = self.root / "work", self.root / "artifact-lane"
            for path in (self.workspace, self.lane, self.workspace / "candidate-build", self.workspace / "baseline-build"):
                path.mkdir(mode=0o700)
                self._directory(path)
            self._profile = recipe(self.workspace, self.canonical, self.pins)
            self._pinned = preflight._json(self.retained._bytes("candidate-build-provenance.json"))
            original = self._retain_metadata(self.canonical / "candidate-build/compile_commands.json")
            require(original.identity["sha256"] == self._pinned["compile_commands"]["sha256"],
                    "historical compile metadata differs from pinned provenance")
            self._candidate_entries = relocated(compile_entries(original.content), self.canonical, self.workspace)
            for path in ("/usr/bin/git", "/usr/bin/cmake"):
                tool = _StreamedTool(path)
                self._stack.callback(tool.close)
                self._tools[path] = tool
            require(self._tools["/usr/bin/git"].sha256 == self._pinned["benchmark_git"]["sha256"],
                    "fresh staging Git differs from the pinned Git executable")
            self._state = "prepared"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def _retain_metadata(self, path):
        require(path not in self._metadata and len(self._metadata) < 128, "duplicate or excessive build metadata")
        snapshot = self._stack.enter_context(provenance._RetainedFileSnapshot(
            path, "v19 build metadata", maximum_bytes=MAX_METADATA_BYTES))
        require(snapshot.resolved == path and snapshot.identity["links"] == 1 and
                snapshot.identity["uid"] == os.getuid() and snapshot.identity["gid"] == os.getgid() and
                not snapshot.identity["mode"] & 0o7022, "unsafe build metadata")
        require(sum(row.identity["size"] for row in self._metadata.values()) +
                snapshot.identity["size"] <= 4 << 20, "total build metadata exceeds bound")
        self._metadata[path] = snapshot
        return snapshot

    def validate_current(self):
        require(self._state in ("prepared", "building", "frozen") and os.getpid() == self._pid,
                "fresh build owner is not live")
        try:
            require(threading.active_count() == 1, "fresh build requires a single-threaded owner")
            self.lease.validate_current()
            self.retained.validate_current()
            self._guard.verify()
            for path, (descriptor, fields) in self._directories.items():
                require(streamed._directory_identity(os.fstat(descriptor)) == fields ==
                        streamed._directory_identity(path.lstat()) and not os.get_inheritable(descriptor),
                        "fresh build directory identity changed")
            validate_recipe(self._profile, self.workspace, self.canonical, self.pins)
            for tool in self._tools.values():
                tool.validate_current()
            for snapshot in self._metadata.values():
                owners.verify_current_bytes(snapshot, snapshot.identity["sha256"])
            if self.authenticated is not None:
                self.authenticated.validate_current(evict_cache=True)
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def _run_stage(self, stage):
        self.validate_current()
        expected = self._profile["stages"][len(self._commands)]
        require(stage == expected, "fresh build command is duplicated or out of order")
        record = copy.deepcopy(stage)
        record["status"] = "running"
        self._commands.append(record)
        if stage["name"].endswith("-build"):
            libc = ctypes.CDLL(None)
            trim = libc.malloc_trim
            trim.argtypes, trim.restype = [ctypes.c_size_t], ctypes.c_int
            record["controller_heap_preparation"] = {"collected": gc.collect(), "malloc_trim_result": trim(0)}
        saved_umask = os.umask(stage["umask"])
        try:
            output = provenance._run(stage["argv"], "v19 " + stage["name"], maximum_bytes=MAX_METADATA_BYTES,
                timeout=600, environment_overrides=ENVIRONMENT,
                executable_descriptor=self._tools[stage["argv"][0]].fd)
            record.update(status="exit-zero", stdout=output.decode("utf-8", errors="strict"),
                          stdout_sha256=hashlib.sha256(output).hexdigest())
        except BaseException as error:
            record.update(status="failed", failure=type(error).__name__ + ": " + str(error))
            raise
        finally:
            os.umask(saved_umask)
        self.validate_current()

    def _validate_build_metadata(self, role):
        build = self.workspace / (role + "-build")
        entries = compile_entries(self._retain_metadata(build / "compile_commands.json").content)
        expected = (self._candidate_entries if role == "candidate" else
                    expected_baseline_entries(self.workspace, self.canonical, self.pins))
        require(entries == expected, role + " compile argv differs from the exact v19 recipe")
        content = self._retain_metadata(build / "CMakeCache.txt").content
        cache = provenance.parse_cmake_cache(content) if role == "candidate" else baseline_cache(content)
        if role == "candidate":
            required = relocated(self._pinned["validated_cache"], self.canonical, self.workspace)
            links = {"leopard": self._pinned["archive_link_commands"],
                     "bench_leopard2": [self._pinned["executable_link_command"]]}
        else:
            required = {"CMAKE_CXX_FLAGS": self._profile["baseline_prefix_map"],
                        "CMAKE_CXX_FLAGS_RELEASE": "-g -O0 -O3", "CMAKE_BUILD_TYPE": "Release",
                        "CMAKE_GENERATOR": "Unix Makefiles", "CMAKE_CXX_COMPILER": "/usr/bin/c++",
                        "LEOPARD_MAIN_SOURCE_DIR": str(self.workspace / "leopard1-source"),
                        "LEO_MAIN_PURE_AVX2": "OFF", "CMAKE_EXPORT_COMPILE_COMMANDS": "ON"}
            objects = list(expected)[:-1]
            libraries = ["/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so", "/usr/lib/x86_64-linux-gnu/libpthread.a"]
            links = {"leopard_main_exact": [["/usr/bin/ar", "qc", "libleopard_main_exact.a", *objects],
                                             ["/usr/bin/ranlib", "libleopard_main_exact.a"]],
                     "leopard_main_benchmark": [["/usr/bin/c++", self._profile["baseline_prefix_map"],
                         "-g", "-O0", "-O3", list(expected)[-1], "-o", "leopard_main_benchmark",
                         "libleopard_main_exact.a", *libraries]]}
        require(all(cache.get(key) == value for key, value in required.items()),
                role + " effective CMake cache differs from the pinned recipe")
        for target, commands in links.items():
            content = self._retain_metadata(build / ("CMakeFiles/" + target + ".dir/link.txt")).content
            actual = [shlex.split(line) for line in content.decode("utf-8").splitlines() if line.strip()]
            require(actual == relocated(commands, self.canonical, self.workspace),
                    role + " link argv differs from the exact v19 recipe")
        # Hold the generated scheduling/flag files too, so a later make-file
        # edit cannot silently diverge from the compile database we checked.
        # This is not a claim that arbitrary compiler/runtime behavior is owned.
        for relative in ("Makefile", "CMakeFiles/Makefile2", "CMakeFiles/Makefile.cmake"):
            self._retain_metadata(build / relative)
        for target_directory in sorted({"/".join(output.split("/")[:2]) for output in entries}):
            require(re.fullmatch(r"CMakeFiles/[A-Za-z0-9_]+\.dir", target_directory) is not None,
                    "compile output has an unexpected target directory")
            for name in ("build.make", "flags.make"):
                self._retain_metadata(build / target_directory / name)

    def build(self):
        require(self._state == "prepared", "fresh build is duplicate or out of order")
        self._state = "building"
        try:
            for stage in self._profile["stages"][:6]:
                destination = self.workspace / stage["name"].split(":")[0]
                if stage["name"].endswith(":clone"):
                    require(not destination.is_symlink() and
                            (not destination.exists() or (destination.is_dir() and not any(destination.iterdir()))),
                            "fresh clone destination is not absent or empty")
                self._run_stage(stage)
                if stage["name"].endswith(":clone"):
                    # Retain each created source directory before checkout or
                    # any later nested clone can use it as a writable parent.
                    # Full file-content ownership follows detached checkout.
                    self._directory(destination, private=False)
            source = self._stack.enter_context(streamed.StreamingSourceOwner(self.workspace))
            self.authenticated = self._stack.enter_context(identity.PinnedSourceIdentity(source, self.retained))
            for stage in self._profile["stages"][6:]:
                self._run_stage(stage)
                if stage["name"].endswith("-configure"):
                    self._validate_build_metadata(stage["name"].split("-")[0])
            self.validate_current()
            self.lease.freeze(self.lane, {role: self.workspace / (role + "-build")
                                         for role in ("candidate", "baseline")})
            self._state = "frozen"
            return self.record()
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        require(self._state == "frozen", "fresh build has not completed")
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-fresh-build-owner/v1", "status": "frozen",
            "recipe": self._profile, "commands": self._commands,
            "source_identity": self.authenticated.record(evict_cache=True), "artifacts": self.lease.record(),
            "metadata": {str(path): snapshot.identity for path, snapshot in self._metadata.items()},
            "launchers": {path: tool.sha256 for path, tool in self._tools.items()},
            "fresh_staging_completed": True, "generated_compile_and_link_argv_verified": True,
            "all_four_pinned_outputs_verified": True, "live_acquisition_armed": False,
            "benchmark_executed": False, "runtime_closure_verified": False,
            "physical_v18_lineage_verified": False, "atomic_snapshot": False})

    def failure_record(self):
        require(self._state == "failed", "fresh build has not failed")
        return copy.deepcopy({"schema": "leopard2-v19-fresh-build-failure/v1", "status": "failed",
            "root": str(self.root) if self.root else None, "commands": self._commands,
            "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state in ("prepared", "building", "frozen"):
                self.validate_current()
            elif value is None:
                raise host.PreflightError("failed fresh build cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
