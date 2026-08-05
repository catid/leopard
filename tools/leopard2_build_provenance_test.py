#!/usr/bin/env python3
"""Regression tests for the shared Leopard2 build-provenance helper."""

from __future__ import annotations

import copy
from contextlib import ExitStack
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import shlex
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


SOURCE_ROOT = Path(sys.argv.pop(1)).resolve(strict=True)
sys.path.insert(0, str(SOURCE_ROOT / "tools"))

import leopard2_build_provenance as provenance  # noqa: E402


def tool_path(name: str) -> Path:
    value = shutil.which(name, path="/usr/bin:/bin")
    if value is None:
        raise unittest.SkipTest(f"{name} is unavailable")
    return Path(value)


def distinct_c_driver() -> Path:
    default = Path(shutil.which("cc", path="/usr/bin:/bin") or
                   "/usr/bin/cc").resolve(strict=True)
    for name in ("clang-18", "clang", "gcc-12", "gcc-14", "gcc-11"):
        candidate = shutil.which(name, path="/usr/bin:/bin")
        if candidate is None:
            continue
        resolved = Path(candidate).resolve(strict=True)
        if resolved != default:
            return resolved
    raise unittest.SkipTest(
        "no C compiler distinct from the default /usr/bin/cc is available")


def compiler_runtime_path(name: str) -> Path:
    compiler = tool_path("c++")
    completed = subprocess.run(
        [str(compiler), f"-print-file-name={name}"],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=30, check=False)
    if completed.returncode != 0:
        raise unittest.SkipTest(
            f"cannot query {name}: " +
            completed.stderr.decode("utf-8", errors="replace"))
    value = completed.stdout.decode("utf-8", errors="strict").strip()
    candidate = Path(value)
    if not candidate.is_absolute() or candidate.name != name or \
            not candidate.exists():
        raise unittest.SkipTest(f"compiler did not resolve {name}")
    return candidate


def write_detached_helper(path: Path) -> None:
    """Write a helper that double-forks, calls setsid, and retains all pass_fds."""
    path.write_text(
        "import os,sys,time\n"
        "mode=sys.argv[2] if len(sys.argv)>2 else 'success'\n"
        "retain_pipes=mode.startswith('retained-pipes')\n"
        "ready_read,ready_write=os.pipe()\n"
        "child=os.fork()\n"
        "if child==0:\n"
        " os.close(ready_read)\n"
        " os.setsid()\n"
        " grandchild=os.fork()\n"
        " if grandchild:\n"
        "  os._exit(0)\n"
        " if not retain_pipes:\n"
        "  null=os.open('/dev/null',os.O_RDWR)\n"
        "  os.dup2(null,1); os.dup2(null,2); os.close(null)\n"
        " with open(sys.argv[1],'w',encoding='ascii') as stream:\n"
        "  stream.write(str(os.getpid())); stream.flush(); os.fsync(stream.fileno())\n"
        " if retain_pipes:\n"
        "  os.write(1,b'retained-pipe-output\\n')\n"
        "  os.write(2,b'retained-pipe-error\\n')\n"
        " os.write(ready_write,b'1'); os.close(ready_write)\n"
        " time.sleep(60); os._exit(0)\n"
        "os.close(ready_write)\n"
        "if os.read(ready_read,1)!=b'1': raise RuntimeError('no daemon')\n"
        "os.close(ready_read); os.waitpid(child,0)\n"
        "if mode=='timeout': time.sleep(60)\n"
        "if mode=='retained-pipes-failure': raise SystemExit(7)\n"
        "if mode=='failure': raise SystemExit(7)\n",
        encoding="utf-8")


def process_gone(pid: int, timeout: float = 2.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not Path("/proc", str(pid)).exists():
            return True
        time.sleep(0.01)
    return not Path("/proc", str(pid)).exists()


def canonical_replay_fixture(
    source: Path, build: Path, *,
    executable_target: str = "bench_leopard2",
) -> tuple[dict[str, object], dict[str, str]]:
    """Return minimal independently bound inputs for a canonical replay."""
    source = source.resolve(strict=True)
    build = build.resolve(strict=True)
    script = source / provenance._CANONICAL_REPLAY_ATTESTATION_SOURCE
    script.parent.mkdir(parents=True, exist_ok=True)
    if not script.exists():
        script.write_bytes(b"# canonical replay test attestation\n")
    library_source = source / "leopard2.cpp"
    benchmark_source = source / "bench/leopard2/benchmark.cpp"
    benchmark_source.parent.mkdir(parents=True, exist_ok=True)
    if not library_source.exists():
        library_source.write_bytes(b"// canonical replay library fixture\n")
    if not benchmark_source.exists():
        benchmark_source.write_bytes(
            b"// canonical replay benchmark fixture\n")

    library_object = "CMakeFiles/leopard.dir/leopard2.cpp.o"
    benchmark_object = (
        f"CMakeFiles/{executable_target}.dir/"
        "bench/leopard2/benchmark.cpp.o")
    (build / Path(library_object).parent).mkdir(
        parents=True, exist_ok=True)
    (build / Path(benchmark_object).parent).mkdir(
        parents=True, exist_ok=True)
    (build / "generated/leopard2-benchmark-attestation").mkdir(
        parents=True, exist_ok=True)

    compiler = str(tool_path("c++"))
    archiver = str(tool_path("ar"))
    ranlib = str(tool_path("ranlib"))
    common_options = [
        "-Wall", "-Wextra", "-fopenmp", "-fopenmp", "-O3",
        "-std=gnu++11", "-I${SOURCE_ROOT}",
    ]

    def compile_record(
        role: str, output: str, relative_source: str,
    ) -> dict[str, object]:
        definitions = ["-DNDEBUG"]
        if role == "archive":
            definitions.extend((
                "-DLEO2_DISABLE_AVX2_CODEGEN=1",
                "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            ))
        else:
            header = (
                "${BUILD_ROOT}/generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            definitions.extend((
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
                "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=" +
                '"' + ("a" * 64) + '"',
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER=" +
                '"' + header + '"',
            ))
        return {
            "role": role,
            "flag_profile": "portable-core",
            "compile_entry": {
                "normalized_arguments": [
                    compiler, *common_options, *definitions, "-O3",
                    "-o", output, "-c",
                    "${SOURCE_ROOT}/" + relative_source,
                ],
            },
        }

    def source_record(path: Path) -> dict[str, object]:
        content = path.read_bytes()
        return {
            "path": path.relative_to(source).as_posix(),
            "sha256": hashlib.sha256(content).hexdigest(),
            "size": len(content),
            "mode": 0o644,
        }

    candidate: dict[str, object] = {
        "executable_target": executable_target,
        "validated_cache": {
            "CMAKE_CXX_COMPILER": compiler,
            "CMAKE_AR": archiver,
            "CMAKE_RANLIB": ranlib,
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "a" * 64,
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
            "LEO2_FLAG_FALIGN_FUNCTIONS_64": "FALSE",
        },
        "tracked_source_manifest": {
            "git": {
                "commit": "1" * 40,
                "tree": "2" * 40,
                "dirty": False,
            },
            "files": sorted((
                source_record(script),
                source_record(library_source),
                source_record(benchmark_source),
            ), key=lambda record: str(record["path"])),
        },
        "source_object_compile_closure": [
            compile_record(
                "archive", library_object, "leopard2.cpp"),
            compile_record(
                "benchmark", benchmark_object,
                "bench/leopard2/benchmark.cpp"),
        ],
        "archive_link_commands": [
            [archiver, "qc", "libleopard.a", library_object],
            [ranlib, "libleopard.a"],
        ],
        "executable_link_command": [
            compiler, "-Wall", "-Wextra", "-fopenmp",
            "-O3", "-DNDEBUG", "-O3", benchmark_object,
            "-o", executable_target, "libleopard.a",
            str(compiler_runtime_path("libgomp.so")),
            str(compiler_runtime_path("libpthread.a")),
        ],
    }
    transports = {
        "cxx": "/proc/self/fd/101",
        "archiver": "/proc/self/fd/102",
        "ranlib": "/proc/self/fd/103",
        "cmake": "/proc/self/fd/104",
        "git": "/proc/self/fd/105",
    }
    return candidate, transports


class StrictMetadataTests(unittest.TestCase):
    def test_strict_json_accepts_canonical_compile_commands_shape(self) -> None:
        value = provenance._strict_json_loads(
            b'[{"directory":"/tmp/build","file":"a.cpp",'
            b'"arguments":["c++","-c","a.cpp"]}]',
            "compile commands")
        self.assertEqual(value[0]["file"], "a.cpp")

    def test_strict_json_rejects_duplicate_keys(self) -> None:
        for encoded in (
            b'{"file":"a.cpp","file":"b.cpp"}',
            b'[{"outer":{"file":"a.cpp","file":"b.cpp"}}]',
        ):
            with self.subTest(encoded=encoded):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError, "duplicate key"):
                    provenance._strict_json_loads(
                        encoded, "compile commands")

    def test_strict_json_rejects_non_finite_numbers(self) -> None:
        for encoded in (
            b'NaN', b'Infinity', b'-Infinity', b'1e999999',
        ):
            with self.subTest(encoded=encoded):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-(?:finite|standard)"):
                    provenance._strict_json_loads(
                        encoded, "compile commands")

    def test_strict_json_rejects_depth_and_malformed_utf8(self) -> None:
        too_deep = (
            b"[" * (provenance.MAX_METADATA_JSON_DEPTH + 2) +
            b"0" +
            b"]" * (provenance.MAX_METADATA_JSON_DEPTH + 2)
        )
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "(?:nesting depth|invalid)"):
            provenance._strict_json_loads(too_deep, "compile commands")
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "invalid"):
            provenance._strict_json_loads(
                b'["valid", "' + bytes((0xff,)) + b'"]',
                "compile commands")

    def test_strict_json_rejects_controls_and_surrogates(self) -> None:
        for encoded in (
            b'["embedded\\u0000nul"]',
            b'["embedded\\nnewline"]',
            b'["unpaired\\ud800surrogate"]',
            b'{"unsafe\\u0001key":"value"}',
            b'["delete\\u007fcontrol"]',
            b'["c1\\u0085control"]',
            b'["bidi\\u202eoverride"]',
            b'["zero\\u200bwidth"]',
            b'["unicode\\u2028line"]',
        ):
            with self.subTest(encoded=encoded):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "unsafe control, format, surrogate"):
                    provenance._strict_json_loads(
                        encoded, "compile commands")

    def test_cmake_cache_rejects_controls_and_format_characters(self) -> None:
        for value in ("\t", "\u2028", "\u2029", "\u202e", "\u200b"):
            with self.subTest(value=repr(value)):
                encoded = (
                    "UNUSED:STRING=left" + value + "right\n").encode(
                        "utf-8")
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError, "unsafe"):
                    provenance.parse_cmake_cache(encoded)
        self.assertEqual(
            provenance.parse_cmake_cache(
                b"CANONICAL:STRING=safe-value\n"),
            {"CANONICAL": "safe-value"})


class StableFileSnapshotTests(unittest.TestCase):
    def test_memfd_falls_back_when_python_omits_wrapper(self) -> None:
        native = getattr(provenance.os, "memfd_create", None)
        descriptor = -1
        try:
            if native is not None:
                delattr(provenance.os, "memfd_create")
            descriptor = provenance._linux_memfd_create(
                "leopard2-python-wrapper-fallback",
                getattr(provenance.os, "MFD_CLOEXEC", 0x0001) |
                getattr(provenance.os, "MFD_ALLOW_SEALING", 0x0002))
            self.assertGreaterEqual(descriptor, 0)
            os.fstat(descriptor)
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            if native is not None:
                provenance.os.memfd_create = native

    def test_replay_artifact_sink_is_frozen_against_same_inode_mutation(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-output-sink-") as directory, \
                provenance._ReplayArtifactSink(
                    Path(directory) / "artifact", "test artifact") as sink:
            os.write(sink.descriptor, b"trusted replay output")
            sink.seal()
            required_seals = (
                getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
                getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
                getattr(fcntl, "F_SEAL_GROW", 0x0004) |
                getattr(fcntl, "F_SEAL_WRITE", 0x0008))
            self.assertEqual(
                fcntl.fcntl(
                    sink.descriptor, getattr(fcntl, "F_GET_SEALS", 1034)) &
                required_seals,
                required_seals)
            with self.assertRaises(OSError):
                os.pwrite(sink.descriptor, b"untrusted", 0)
            self.assertEqual(sink.content, b"trusted replay output")
            sink.verify()

    def test_private_tmpfs_descriptor_survives_lexical_root_swap(
            self) -> None:
        """A post-mount pathname replacement cannot redirect stage writes."""
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-private-root-") as directory:
            workspace = Path(directory)
            stage = workspace / "stage"
            saved = workspace / "saved-stage"
            victim = workspace / "victim"
            stage.mkdir()
            victim.write_bytes(b"victim-is-unchanged")
            stage_descriptor = os.open(
                stage,
                getattr(os, "O_PATH", os.O_RDONLY) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_CLOEXEC", 0))
            ready_read, ready_write = os.pipe()
            continue_read, continue_write = os.pipe()
            error_read, error_write = os.pipe()
            child = os.fork()
            if child == 0:
                try:
                    os.close(ready_read)
                    os.close(continue_write)
                    os.close(error_read)
                    provenance._mount_private_replay_tmpfs(
                        stage, stage_descriptor, ("out",),
                        os.getuid(), os.getgid())
                    os.write(ready_write, b"1")
                    if os.read(continue_read, 1) != b"1":
                        raise RuntimeError("parent did not release child")
                    output = os.open(
                        f"/proc/self/fd/{stage_descriptor}/out/object",
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
                    try:
                        os.write(output, b"private-output")
                    finally:
                        os.close(output)
                    os._exit(0)
                except BaseException as error:
                    try:
                        os.write(
                            error_write,
                            (type(error).__name__ + ": " + str(error)).
                            encode("utf-8", errors="replace"))
                    finally:
                        os._exit(1)

            os.close(ready_write)
            os.close(continue_read)
            os.close(error_write)
            try:
                ready = os.read(ready_read, 1)
                if ready == b"1":
                    stage.rename(saved)
                    (stage / "out").mkdir(parents=True)
                    os.link(victim, stage / "out/object")
                    os.write(continue_write, b"1")
                else:
                    try:
                        os.write(continue_write, b"0")
                    except BrokenPipeError:
                        pass
                _pid, wait_status = os.waitpid(child, 0)
                error = os.read(error_read, 65536).decode(
                    "utf-8", errors="replace")
                self.assertEqual(
                    ready, b"1",
                    "private tmpfs child failed before synchronization: " +
                    error)
                self.assertTrue(
                    os.WIFEXITED(wait_status) and
                    os.WEXITSTATUS(wait_status) == 0,
                    "private tmpfs child failed: " + error)
                self.assertEqual(victim.read_bytes(), b"victim-is-unchanged")
                self.assertEqual(
                    (stage / "out/object").read_bytes(),
                    b"victim-is-unchanged")
                self.assertFalse((saved / "out/object").exists())
            finally:
                os.close(ready_read)
                os.close(continue_write)
                os.close(error_read)
                os.close(stage_descriptor)

    def test_landlocked_cat_populates_only_inherited_memfd_sink(self) -> None:
        """Private-stage aliases cannot redirect the retained sink write."""
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-landlock-sink-") as directory:
            workspace = Path(directory)
            stage = workspace / "stage"
            victim = workspace / "victim"
            stage.mkdir()
            victim.write_bytes(b"victim-is-unchanged")
            stage_descriptor = os.open(
                stage,
                getattr(os, "O_PATH", os.O_RDONLY) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_CLOEXEC", 0))
            try:
                with provenance._ReplayArtifactSink(
                        workspace / "artifact", "Landlock artifact") as sink, \
                        provenance._RetainedFileSnapshot(
                            tool_path("bash"), "Landlock test Bash") as bash, \
                        provenance._RetainedFileSnapshot(
                            tool_path("cat"), "Landlock test cat") as cat:
                    command = (
                        f"cd /proc/self/fd/{stage_descriptor}; "
                        f"if ln {shlex.quote(str(victim))} out/object; then "
                        f"printf corrupt > out/object; "
                        f"else printf sealed-output > out/object; fi; "
                        f"/proc/self/fd/{cat.executable_descriptor} "
                        f"-- out/object 1>&{sink.descriptor}")
                    provenance._run(
                        [str(tool_path("bash")), "-c", command],
                        "Landlocked inherited artifact sink",
                        inherited_descriptors=(
                            bash.executable_descriptor,
                            cat.executable_descriptor,
                            sink.descriptor,
                            stage_descriptor),
                        executable_descriptor=bash.executable_descriptor,
                        write_sandbox_root=workspace,
                        write_sandbox_descriptors=(sink.descriptor,),
                        private_tmpfs_root=stage,
                        private_tmpfs_directories=("out",),
                        private_tmpfs_descriptor=stage_descriptor,
                        timeout=10)
                    sink.seal()
                    self.assertEqual(sink.content, b"sealed-output")
                    self.assertEqual(
                        victim.read_bytes(), b"victim-is-unchanged")
                    self.assertFalse((stage / "out/object").exists())
            finally:
                os.close(stage_descriptor)

    def test_owned_descriptor_callers_close_post_open_interruptions(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-owned-fd-") as directory:
            root = Path(directory)

            def exercise(label: str, callback) -> None:
                captured: list[int] = []
                armed = False
                real_open = os.open

                def tracked_open(*args, **kwargs) -> int:
                    nonlocal armed
                    descriptor = real_open(*args, **kwargs)
                    captured.append(descriptor)
                    armed = True
                    return descriptor

                owner_code = provenance._OwnedDescriptor.open.__code__

                def interrupt_after_open(frame, event: str, argument):
                    nonlocal armed
                    del argument
                    if (armed and event == "line" and
                            frame.f_code is owner_code):
                        armed = False
                        raise KeyboardInterrupt(
                            f"injected {label} post-open interruption")
                    return interrupt_after_open

                previous_trace = sys.gettrace()
                try:
                    with mock.patch.object(
                            provenance.os, "open",
                            side_effect=tracked_open):
                        sys.settrace(interrupt_after_open)
                        with self.assertRaisesRegex(
                                KeyboardInterrupt, "post-open interruption"):
                            callback()
                finally:
                    sys.settrace(previous_trace)
                self.assertEqual(len(captured), 1)
                with self.assertRaises(OSError):
                    os.fstat(captured[0])

            open_root = root / "open-root"
            open_root.mkdir()
            exercise(
                "open tree root",
                lambda: provenance._OpenDirectoryTree(
                    open_root, provenance.ExitStack(), "open tree"))

            create_root = root / "created-root"
            exercise(
                "create tree root",
                lambda: provenance._CreateDirectoryTree(
                    create_root, provenance.ExitStack(), "create tree"))

            executable = root / "generated-executable"
            exercise(
                "private executable",
                lambda: provenance._write_private_executable(
                    executable, b"#!/bin/sh\nexit 0\n"))

            replace_target = root / "replace-target"
            replace_target.write_bytes(b"old")
            exercise(
                "private replacement",
                lambda: provenance._replace_private_file(
                    replace_target, b"new", 0o600))

            with provenance.ExitStack() as stack:
                tree_root = root / "write-tree"
                tree = provenance._CreateDirectoryTree(
                    tree_root, stack, "write tree")
                exercise(
                    "tree file",
                    lambda: tree.write_exclusive(
                        "generated/file", b"content", 0o600))

    def test_owned_descriptor_does_not_close_recycled_fd(self) -> None:
        target = tool_path("bash")
        owner = provenance._OwnedDescriptor()
        descriptor = owner.open(target, os.O_RDONLY)
        replacement = -1
        real_close = os.close

        def close_then_interrupt(value: int) -> None:
            nonlocal replacement
            self.assertEqual(value, descriptor)
            real_close(value)
            replacement = os.open("/dev/null", os.O_RDONLY)
            self.assertEqual(replacement, descriptor)
            raise KeyboardInterrupt("injected post-close interruption")

        with mock.patch.object(
                provenance.os, "close", side_effect=close_then_interrupt):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "post-close interruption"):
                owner.close()
        try:
            os.fstat(replacement)
            self.assertEqual(owner.descriptor, -1)
            owner.close()
            os.fstat(replacement)
        finally:
            real_close(replacement)

    def test_owned_descriptor_guard_closes_fd_when_close_did_not_run(
            self) -> None:
        target = tool_path("bash")
        owner = provenance._OwnedDescriptor()
        descriptor = owner.open(target, os.O_RDONLY)

        def interrupt_before_close(value: int) -> None:
            self.assertEqual(value, descriptor)
            raise KeyboardInterrupt("injected pre-close interruption")

        with mock.patch.object(
                provenance.os, "close",
                side_effect=interrupt_before_close):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "pre-close interruption"):
                owner.close()
        self.assertEqual(owner.descriptor, -1)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_owned_descriptor_same_inode_recycle_is_never_retry_closed(
            self) -> None:
        for failure_type in (OSError, KeyboardInterrupt, SystemExit):
            with self.subTest(failure=failure_type.__name__), \
                    tempfile.TemporaryDirectory(
                        prefix="leo2-owned-same-inode-aba-") as directory:
                path = Path(directory) / "same-inode"
                path.write_bytes(b"same inode, distinct descriptions")
                owner = provenance._OwnedDescriptor()
                descriptor = owner.open(path, os.O_RDONLY)
                replacement = -1
                guards: list[int] = []
                real_close = os.close
                real_dup = os.dup
                failure = failure_type(
                    f"injected same-inode {failure_type.__name__}")

                def tracked_dup(value: int) -> int:
                    self.assertEqual(value, descriptor)
                    guard = real_dup(value)
                    guards.append(guard)
                    return guard

                def close_reopen_same_inode(value: int) -> None:
                    nonlocal replacement
                    self.assertEqual(value, descriptor)
                    real_close(value)
                    replacement = os.open(path, os.O_RDONLY)
                    self.assertEqual(replacement, descriptor)
                    raise failure

                try:
                    with mock.patch.object(
                            provenance.os, "dup",
                            side_effect=tracked_dup), mock.patch.object(
                                provenance.os, "close",
                                side_effect=close_reopen_same_inode):
                        with self.assertRaises(failure_type) as raised:
                            owner.close()
                    self.assertIs(raised.exception, failure)
                    self.assertEqual(owner.descriptor, -1)
                    os.fstat(replacement)
                    owner.close()
                    os.fstat(replacement)
                    self.assertEqual(len(guards), 1)
                    with self.assertRaises(OSError):
                        os.fstat(guards[0])
                finally:
                    if replacement >= 0:
                        real_close(replacement)

    def test_owned_descriptor_clears_fd_when_close_preceded_interruption(
            self) -> None:
        target = tool_path("bash")
        owner = provenance._OwnedDescriptor()
        descriptor = owner.open(target, os.O_RDONLY)
        real_close = os.close

        def close_then_interrupt(value: int) -> None:
            self.assertEqual(value, descriptor)
            real_close(value)
            raise KeyboardInterrupt("injected post-close interruption")

        with mock.patch.object(
                provenance.os, "close", side_effect=close_then_interrupt):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "post-close interruption"):
                owner.close()
        self.assertEqual(owner.descriptor, -1)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_interrupted_inotify_initialization_closes_descriptor(self) -> None:
        captured: list[int] = []

        class InitFunction:
            argtypes = None
            restype = None

            def __call__(self, unused_flags) -> int:
                descriptor = os.open("/dev/null", os.O_RDONLY)
                captured.append(descriptor)
                return descriptor

        class Library:
            inotify_init1 = InitFunction()

        real_cdll = provenance.ctypes.CDLL
        served_inotify_library = False

        def load_library(*args, **kwargs):
            nonlocal served_inotify_library
            if not served_inotify_library:
                served_inotify_library = True
                return Library()
            # Constructor cleanup obtains libc again for its one-shot close;
            # do not let the inotify-only fixture replace that independent
            # primitive with an object that has no close symbol.
            return real_cdll(*args, **kwargs)

        # Raising from fdopen is the first operation after inotify_init1 has
        # returned its owned descriptor.  This deterministic boundary tests
        # the same constructor cleanup without a line tracer that can re-enter
        # the exception handler on newer CPython versions.
        with mock.patch.object(
                provenance.ctypes, "CDLL", side_effect=load_library), \
                mock.patch.object(
                    provenance.os, "fdopen",
                    side_effect=KeyboardInterrupt(
                        "injected post-inotify-init interruption")):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "post-inotify-init"):
                provenance._InotifyMutationGuard(
                    "interrupted inotify")
        self.assertEqual(len(captured), 1)
        with self.assertRaises(OSError):
            os.fstat(captured[0])

    def test_snapshot_accepts_a_stable_symlink_alias(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-snapshot-") as directory:
            root = Path(directory)
            target = root / "compiler"
            target.write_bytes(b"stable compiler")
            alias = root / "c++"
            alias.symlink_to(target.name)
            identity, content = provenance.file_snapshot(alias, "compiler")
            self.assertEqual(content, b"stable compiler")
            self.assertEqual(identity["path"], str(target.resolve(strict=True)))
            status = target.stat()
            self.assertEqual(identity["device"], status.st_dev)
            self.assertEqual(identity["inode"], status.st_ino)
            self.assertEqual(identity["uid"], status.st_uid)
            self.assertEqual(identity["gid"], status.st_gid)
            self.assertEqual(identity["links"], status.st_nlink)

    def test_snapshot_rejects_final_component_symlink_swap(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-snapshot-") as directory:
            root = Path(directory)
            target = root / "artifact"
            target.write_bytes(b"expected")
            replacement = root / "replacement"
            replacement.write_bytes(b"substituted")
            saved = root / "saved"
            real_open = os.open
            swapped = False

            def swapping_open(
                path: os.PathLike[str] | str, flags: int, *args: object,
                **kwargs: object,
            ) -> int:
                nonlocal swapped
                if not swapped and Path(path) == target:
                    swapped = True
                    target.rename(saved)
                    target.symlink_to(replacement.name)
                return real_open(path, flags, *args, **kwargs)

            with mock.patch.object(
                    provenance.os, "open", side_effect=swapping_open):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "cannot open stable artifact"):
                    provenance.file_snapshot(target, "artifact")

    def test_snapshot_rejects_parent_directory_swap(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-snapshot-") as directory:
            root = Path(directory)
            active = root / "active"
            alternate = root / "alternate"
            active.mkdir()
            alternate.mkdir()
            target = active / "artifact"
            target.write_bytes(b"expected")
            (alternate / "artifact").write_bytes(b"substituted")
            saved = root / "saved"
            real_open = os.open
            swapped = False

            def swapping_open(
                path: os.PathLike[str] | str, flags: int, *args: object,
                **kwargs: object,
            ) -> int:
                nonlocal swapped
                if not swapped and Path(path) == target:
                    swapped = True
                    active.rename(saved)
                    alternate.rename(active)
                return real_open(path, flags, *args, **kwargs)

            with mock.patch.object(
                    provenance.os, "open", side_effect=swapping_open):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "pathname changed"):
                    provenance.file_snapshot(target, "artifact")

    def test_retained_descriptor_defeats_executable_path_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-aba-") as directory:
            root = Path(directory)
            executable = root / "tool"
            saved = root / "saved-tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "retained executable")
            try:
                executable.rename(saved)
                shutil.copyfile(tool_path("false"), executable)
                executable.chmod(0o700)
                self.assertEqual(
                    provenance._run(
                        [str(executable)], "descriptor-bound executable",
                        executable_descriptor=snapshot.executable_descriptor,
                        timeout=3),
                    b"")
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:pathname|changed while retained)"):
                    snapshot.close()
            finally:
                snapshot._close_without_verification()

    def test_sealed_executable_defeats_same_inode_byte_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-byte-aba-") as directory:
            executable = Path(directory) / "tool"
            original = tool_path("true").read_bytes()
            alternate = tool_path("false").read_bytes()
            executable.write_bytes(original)
            executable.chmod(0o700)
            inode = executable.stat().st_ino
            snapshot = provenance._RetainedFileSnapshot(
                executable, "same-inode executable")
            try:
                sealed = snapshot.executable_descriptor
                self.assertEqual(
                    snapshot.executable_identity["sha256"],
                    hashlib.sha256(original).hexdigest())
                self.assertEqual(
                    snapshot.executable_identity["source_sha256"],
                    snapshot.executable_identity["sha256"])
                required_seals = (
                    getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
                    getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
                    getattr(fcntl, "F_SEAL_GROW", 0x0004) |
                    getattr(fcntl, "F_SEAL_WRITE", 0x0008))
                self.assertEqual(
                    snapshot.executable_identity["seals"] & required_seals,
                    required_seals)
                with executable.open("r+b") as stream:
                    stream.write(alternate)
                    stream.truncate()
                    stream.flush()
                    os.fsync(stream.fileno())
                self.assertEqual(executable.stat().st_ino, inode)
                self.assertEqual(
                    provenance._run(
                        [str(executable)], "sealed same-inode executable",
                        executable_descriptor=sealed, timeout=3),
                    b"")
                with executable.open("r+b") as stream:
                    stream.write(original)
                    stream.truncate()
                    stream.flush()
                    os.fsync(stream.fileno())
                self.assertEqual(executable.stat().st_ino, inode)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:pathname changed|changed while retained)"):
                    snapshot.close()
            finally:
                snapshot._close_without_verification()

    def test_interrupted_executable_snapshot_does_not_own_closed_fd(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-interrupt-") as directory:
            executable = Path(directory) / "tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "interrupted executable")
            interrupted_descriptors: list[int] = []

            def interrupt_read(
                    descriptor: int, unused_size: int, unused_label: str,
                    unused_maximum: int) -> bytes:
                interrupted_descriptors.append(descriptor)
                raise KeyboardInterrupt("injected immutable read interruption")

            replacement = -1
            try:
                with mock.patch.object(
                        provenance, "_read_bounded_descriptor",
                        side_effect=interrupt_read):
                    with self.assertRaisesRegex(
                            KeyboardInterrupt, "injected immutable"):
                        snapshot.executable_descriptor
                self.assertEqual(snapshot._executable_descriptor, -1)
                self.assertEqual(snapshot.executable_identity, {})
                self.assertEqual(len(interrupted_descriptors), 1)
                with self.assertRaises(OSError):
                    os.fstat(interrupted_descriptors[0])

                # Linux immediately reuses the lowest available descriptor.
                # Snapshot cleanup must not treat that recycled number as its
                # own immutable executable and close an unrelated file.
                replacement = os.open("/dev/null", os.O_RDONLY)
                self.assertEqual(replacement, interrupted_descriptors[0])
                snapshot.close()
                os.fstat(replacement)
            finally:
                snapshot._close_without_verification()
                if replacement >= 0:
                    try:
                        os.close(replacement)
                    except OSError:
                        pass

    def test_interrupted_fdopen_closes_unowned_memfd(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-fdopen-") as directory:
            executable = Path(directory) / "tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "fdopen-interrupted executable")
            captured: list[int] = []

            def interrupt_fdopen(
                    descriptor: int, unused_mode: str, *,
                    buffering: int) -> object:
                self.assertEqual(buffering, 0)
                captured.append(descriptor)
                raise KeyboardInterrupt("injected fdopen interruption")

            try:
                with mock.patch.object(
                        provenance.os, "fdopen",
                        side_effect=interrupt_fdopen):
                    with self.assertRaisesRegex(
                            KeyboardInterrupt, "injected fdopen"):
                        snapshot.executable_descriptor
                self.assertEqual(len(captured), 1)
                self.assertEqual(snapshot._executable_descriptor, -1)
                self.assertEqual(snapshot.executable_identity, {})
                with self.assertRaises(OSError):
                    os.fstat(captured[0])
                snapshot.close()
            finally:
                snapshot._close_without_verification()

    def test_interrupted_after_memfd_create_closes_raw_descriptor(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-memfd-") as directory:
            executable = Path(directory) / "tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "memfd-interrupted executable")
            real_memfd_create = provenance._linux_memfd_create
            captured: list[int] = []
            armed = False
            property_code = type(snapshot).executable_descriptor.fget.__code__

            def recording_memfd_create(
                    name: str, flags: int = 0) -> int:
                nonlocal armed
                descriptor = real_memfd_create(name, flags)
                captured.append(descriptor)
                armed = True
                return descriptor

            def interrupt_after_return(frame, event: str, argument):
                nonlocal armed
                del argument
                if armed and event == "line" and frame.f_code is property_code:
                    armed = False
                    sys.settrace(None)
                    raise KeyboardInterrupt(
                        "injected post-memfd-create interruption")
                return interrupt_after_return

            try:
                with mock.patch.object(
                        provenance, "_linux_memfd_create",
                        side_effect=recording_memfd_create):
                    previous_trace = sys.gettrace()
                    sys.settrace(interrupt_after_return)
                    try:
                        with self.assertRaisesRegex(
                                KeyboardInterrupt, "post-memfd-create"):
                            snapshot.executable_descriptor
                    finally:
                        sys.settrace(previous_trace)
                self.assertEqual(len(captured), 1)
                self.assertEqual(snapshot._executable_descriptor, -1)
                self.assertEqual(snapshot.executable_identity, {})
                with self.assertRaises(OSError):
                    os.fstat(captured[0])
                snapshot.close()
            finally:
                sys.settrace(None)
                snapshot._close_without_verification()

    def test_interrupted_close_retries_before_source_close(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-close-source-") as directory:
            executable = Path(directory) / "tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "source-close-interrupted executable")
            source_descriptor = snapshot.descriptor
            sealed_descriptor = snapshot.executable_descriptor
            close_code = type(snapshot).close.__code__
            injected = False

            def interrupt_before_source_close(frame, event: str, argument):
                nonlocal injected
                del argument
                source = frame.f_locals.get("source")
                if not injected and event == "line" and \
                        frame.f_code is close_code and \
                        source is snapshot._source_file and \
                        source is not None and not source.closed:
                    injected = True
                    raise KeyboardInterrupt(
                        "injected before retained source close")
                return interrupt_before_source_close

            try:
                previous_trace = sys.gettrace()
                sys.settrace(interrupt_before_source_close)
                try:
                    with self.assertRaisesRegex(
                            KeyboardInterrupt, "before retained source"):
                        snapshot.close()
                finally:
                    sys.settrace(previous_trace)
                self.assertTrue(injected)
                os.fstat(source_descriptor)
                os.fstat(sealed_descriptor)
                snapshot.close()
                with self.assertRaises(OSError):
                    os.fstat(source_descriptor)
                with self.assertRaises(OSError):
                    os.fstat(sealed_descriptor)
            finally:
                sys.settrace(None)
                snapshot._close_without_verification()

    def test_interrupted_close_retries_after_source_close(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-close-exec-") as directory:
            executable = Path(directory) / "tool"
            shutil.copyfile(tool_path("true"), executable)
            executable.chmod(0o700)
            snapshot = provenance._RetainedFileSnapshot(
                executable, "exec-close-interrupted executable")
            source_descriptor = snapshot.descriptor
            sealed_descriptor = snapshot.executable_descriptor
            close_code = type(snapshot).close.__code__
            injected = False

            def interrupt_after_source_close(frame, event: str, argument):
                nonlocal injected
                del argument
                source = snapshot._source_file
                if not injected and event == "line" and \
                        frame.f_code is close_code and \
                        source is not None and source.closed and \
                        snapshot._executable_descriptor >= 0:
                    injected = True
                    raise KeyboardInterrupt(
                        "injected after retained source close")
                return interrupt_after_source_close

            try:
                previous_trace = sys.gettrace()
                sys.settrace(interrupt_after_source_close)
                try:
                    with self.assertRaisesRegex(
                            KeyboardInterrupt, "after retained source"):
                        snapshot.close()
                finally:
                    sys.settrace(previous_trace)
                self.assertTrue(injected)
                with self.assertRaises(OSError):
                    os.fstat(source_descriptor)
                os.fstat(sealed_descriptor)
                snapshot.close()
                with self.assertRaises(OSError):
                    os.fstat(sealed_descriptor)
            finally:
                sys.settrace(None)
                snapshot._close_without_verification()

    def test_sealed_wrapper_preserves_all_colliding_inherited_fds(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-exec-fds-") as directory:
            root = Path(directory)
            marker_path = root / "marker"
            marker_path.write_bytes(b"collision-safe\n")
            fillers = [os.open("/dev/null", os.O_RDONLY) for _ in range(48)]
            marker_descriptor = os.open(marker_path, os.O_RDONLY)
            try:
                with provenance._RetainedFileSnapshot(
                        tool_path("bash"), "sealed wrapper Bash") as bash:
                    wrapper_path = root / "wrapper"
                    command = (
                        'IFS= read -r value < "/proc/self/fd/$1"; '
                        'printf %s "$value"')
                    provenance._write_replay_exec_wrapper(
                        wrapper_path, bash.executable_descriptor,
                        str(tool_path("bash")), bash.executable_descriptor,
                        ("-c", command, "leo2-wrapper"))
                    with provenance._RetainedFileSnapshot(
                            wrapper_path, "sealed descriptor wrapper"
                            ) as wrapper:
                        output = provenance._run(
                            [str(wrapper_path), str(marker_descriptor)],
                            "collision-safe inherited descriptors",
                            inherited_descriptors=(
                                marker_descriptor,
                                bash.executable_descriptor,
                                wrapper.executable_descriptor,
                                marker_descriptor,
                            ),
                            executable_descriptor=
                            wrapper.executable_descriptor,
                            timeout=3)
                        self.assertEqual(output, b"collision-safe")
            finally:
                os.close(marker_descriptor)
                for descriptor in fillers:
                    os.close(descriptor)

    def test_generated_executable_rejects_pre_retention_replacement(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-generated-") as directory:
            executable = Path(directory) / "generated"
            expected = b"#!/bin/sh\nexit 0\n"
            substituted = b"#!/bin/sh\nexit 7\n"
            provenance._write_private_executable(executable, expected)
            real_init = provenance._RetainedFileSnapshot.__init__
            replaced = False

            def replacing_init(
                    instance, path: Path | str, label: str, *,
                    maximum_bytes: int = provenance.MAX_FILE_BYTES) -> None:
                nonlocal replaced
                if not replaced:
                    Path(path).write_bytes(substituted)
                    Path(path).chmod(0o700)
                    replaced = True
                real_init(
                    instance, path, label, maximum_bytes=maximum_bytes)

            with mock.patch.object(
                    provenance._RetainedFileSnapshot, "__init__",
                    new=replacing_init):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "changed between generation and retention"):
                    provenance._retain_exact_generated_executable(
                        executable, expected, "generated replay executable")
            self.assertTrue(replaced)
            self.assertEqual(executable.read_bytes(), substituted)

    def test_compiler_prefix_rejects_pre_guard_symlink_replacement(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-prefix-") as directory:
            prefix = Path(directory) / "prefix"
            prefix.mkdir()
            role = prefix / "ld"
            expected = "/proc/self/fd/123"
            substituted = "/proc/self/fd/124"
            role.symlink_to(expected)
            real_init = provenance._RetainedDirectoryTree.__init__
            replaced = False

            def replacing_init(instance, path: Path | str, label: str) -> None:
                nonlocal replaced
                if not replaced:
                    role.unlink()
                    role.symlink_to(substituted)
                    replaced = True
                real_init(instance, path, label)

            with mock.patch.object(
                    provenance._RetainedDirectoryTree, "__init__",
                    new=replacing_init):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "mappings changed before retention"):
                    provenance._retain_exact_symlink_directory(
                        prefix, {"ld": expected}, "compiler prefix")
            self.assertTrue(replaced)
            self.assertEqual(os.readlink(role), substituted)

    def test_parent_directory_aba_is_retained_as_an_event(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-parent-aba-") as directory:
            root = Path(directory)
            active = root / "active"
            alternate = root / "alternate"
            saved = root / "saved"
            active.mkdir()
            alternate.mkdir()
            for target, marker in (
                    (active / "tool", "A"),
                    (alternate / "tool", "B")):
                target.write_text(
                    "#!/bin/sh\nprintf %s " + marker + "\n",
                    encoding="ascii")
                target.chmod(0o700)
            path = active / "tool"
            snapshot = provenance._RetainedFileSnapshot(
                path, "parent ABA tool")
            try:
                active.rename(saved)
                alternate.rename(active)
                self.assertEqual(
                    provenance._run(
                        [str(path)], "pathname tool during ABA", timeout=3),
                    b"B")
                active.rename(alternate)
                saved.rename(active)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "pathname changed while retained"):
                    snapshot.close()
            finally:
                snapshot._close_without_verification()

    def test_tracked_source_copy_is_complete_and_mutation_guarded(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-source-copy-") as directory:
            destination = Path(directory) / "source"
            snapshot = provenance._RetainedPrivateSourceTree(
                SOURCE_ROOT, destination)
            try:
                records = snapshot.manifest["files"]
                self.assertGreater(len(records), 100)
                relative = next(
                    record["path"] for record in records
                    if record["path"] == "leopard2.h")
                target = destination / relative
                original = target.read_bytes()
                target.write_bytes(original + b"\n")
                target.write_bytes(original)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:pathname changed|snapshot changed)"):
                    snapshot.close()
            finally:
                snapshot.guard._close_without_verification()

    def test_tracked_source_rejects_git_directory_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-gitdir-aba-") as directory:
            root = Path(directory) / "repo"
            root.mkdir()
            subprocess.run(
                [str(tool_path("git")), "init", "-q"], cwd=root, check=True)
            subprocess.run(
                [str(tool_path("git")), "config", "user.name", "Audit"],
                cwd=root, check=True)
            subprocess.run(
                [str(tool_path("git")), "config", "user.email",
                 "audit@example.invalid"], cwd=root, check=True)
            (root / "tracked.txt").write_text("A\n", encoding="ascii")
            subprocess.run(
                [str(tool_path("git")), "add", "tracked.txt"],
                cwd=root, check=True)
            subprocess.run(
                [str(tool_path("git")), "commit", "-q", "-m", "A"],
                cwd=root, check=True)
            alternate = root / ".git-alternate"
            saved = root / ".git-saved"
            shutil.copytree(root / ".git", alternate)

            original = provenance._git_source_state
            calls = 0

            def git_state_with_aba(*args: object, **kwargs: object):
                nonlocal calls
                calls += 1
                if calls == 1:
                    result = original(*args, **kwargs)
                    (root / ".git").rename(saved)
                    alternate.rename(root / ".git")
                    return result
                (root / ".git").rename(alternate)
                saved.rename(root / ".git")
                return original(*args, **kwargs)

            with mock.patch.object(
                    provenance, "_git_source_state",
                    side_effect=git_state_with_aba):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:pathname changed|directory pathname changed)"):
                    provenance._capture_tracked_source_tree(root)
            self.assertEqual(calls, 2)
            self.assertTrue((root / ".git").is_dir())

    def test_worktree_gitfile_is_retained_and_rejects_byte_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-worktree-gitfile-") as directory:
            main = Path(directory) / "main"
            worktree = Path(directory) / "worktree"
            main.mkdir()
            subprocess.run(
                [str(tool_path("git")), "init", "-q"], cwd=main, check=True)
            subprocess.run(
                [str(tool_path("git")), "config", "user.name", "Audit"],
                cwd=main, check=True)
            subprocess.run(
                [str(tool_path("git")), "config", "user.email",
                 "audit@example.invalid"], cwd=main, check=True)
            (main / "tracked.txt").write_text("A\n", encoding="ascii")
            subprocess.run(
                [str(tool_path("git")), "add", "tracked.txt"],
                cwd=main, check=True)
            subprocess.run(
                [str(tool_path("git")), "commit", "-q", "-m", "A"],
                cwd=main, check=True)
            subprocess.run(
                [str(tool_path("git")), "worktree", "add", "-q", "-b",
                 "audit-worktree", str(worktree)],
                cwd=main, check=True)
            manifest = provenance._capture_tracked_source_tree(worktree)
            self.assertEqual(
                [record["path"] for record in manifest["files"]],
                ["tracked.txt"])
            gitfile = worktree / ".git"
            original_content = gitfile.read_bytes()
            original = provenance._git_source_state
            calls = 0

            def git_state_with_byte_aba(*args: object, **kwargs: object):
                nonlocal calls
                calls += 1
                if calls == 1:
                    result = original(*args, **kwargs)
                    gitfile.write_bytes(original_content + b" ")
                    return result
                gitfile.write_bytes(original_content)
                return original(*args, **kwargs)

            with mock.patch.object(
                    provenance, "_git_source_state",
                    side_effect=git_state_with_byte_aba):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:pathname changed|changed while retained)"):
                    provenance._capture_tracked_source_tree(worktree)
            self.assertEqual(calls, 2)
            self.assertEqual(gitfile.read_bytes(), original_content)

    def test_source_tree_nested_directory_aba_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-tree-aba-") as directory:
            root = Path(directory) / "source"
            active = root / "nested"
            alternate = root / "alternate"
            saved = root / "saved"
            active.mkdir(parents=True)
            alternate.mkdir()
            (active / "input.h").write_bytes(b"A")
            (alternate / "input.h").write_bytes(b"B")
            guard = provenance._InotifyMutationGuard("source tree ABA")
            try:
                guard.add_tree(root)
                active.rename(saved)
                alternate.rename(active)
                self.assertEqual((active / "input.h").read_bytes(), b"B")
                active.rename(alternate)
                saved.rename(active)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "pathname changed while retained"):
                    guard.close()
            finally:
                guard._close_without_verification()

    def test_descriptor_exec_preserves_lexical_argv_zero(self) -> None:
        shell = tool_path("sh")
        lexical_name = "/synthetic/original/tool-name"
        with provenance._RetainedFileSnapshot(
                shell, "retained shell") as snapshot:
            output = provenance._run(
                [lexical_name, "-c", "printf %s \"$0\""],
                "descriptor argv zero", timeout=3,
                executable_descriptor=snapshot.executable_descriptor)
        self.assertEqual(output, lexical_name.encode("ascii"))


class CanonicalReplayPlanTests(unittest.TestCase):
    def _render_plan(
        self, root: Path,
    ) -> tuple[
        Path, Path, dict[str, object], dict[str, str], bytes, dict[str, int],
    ]:
        source = root / "source tree"
        build = root / "build tree"
        source.mkdir()
        (build / "CMakeFiles").mkdir(parents=True)
        candidate, transports = canonical_replay_fixture(source, build)
        content, counts = provenance._canonical_replay_makefile_bytes(
            candidate, source.resolve(), build.resolve(), transports)
        return (
            source.resolve(), build.resolve(), candidate, transports,
            content, counts)

    def _run_while_ignoring_generated_input(
        self, relative: str, payload: bytes,
    ) -> tuple[bytes, subprocess.CompletedProcess[bytes]]:
        with tempfile.TemporaryDirectory(
                prefix="leo2 canonical replay exploits ") as directory:
            root = Path(directory)
            (source, build, candidate, _transports,
             _content, _counts) = self._render_plan(root)
            helper = root / "unretained-helper"
            marker = root / "unexpected-marker"
            helper.write_text(
                "#!/bin/sh\n"
                f"printf executed > {shlex.quote(str(marker))}\n",
                encoding="utf-8")
            helper.chmod(0o700)
            generated = build / relative
            generated.parent.mkdir(parents=True, exist_ok=True)
            generated.write_bytes(
                payload.replace(b"@HELPER@", os.fsencode(helper)))

            true_descriptor = os.open(
                tool_path("true"), os.O_RDONLY |
                getattr(os, "O_CLOEXEC", 0))
            try:
                transports = {
                    role: f"/proc/self/fd/{true_descriptor}"
                    for role in ("cxx", "archiver", "ranlib", "cmake", "git")
                }
                content, _counts = \
                    provenance._canonical_replay_makefile_bytes(
                        candidate, source, build, transports)
                plan_path = (
                    build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE)
                provenance._write_private_executable(plan_path, content)
                with provenance._RetainedFileSnapshot(
                        plan_path, "canonical replay exploit test plan"
                        ) as plan:
                    completed = subprocess.run(
                        [
                            str(tool_path("make")), "-C", str(build),
                            "-f",
                            f"/proc/self/fd/{plan.executable_descriptor}",
                            "--",
                            "CMakeFiles/bench_leopard2.dir/replay",
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        env=dict(provenance.GIT_ENVIRONMENT),
                        pass_fds=(
                            true_descriptor,
                            plan.executable_descriptor,
                        ),
                        timeout=30, check=False)
                    plan.verify()
            finally:
                os.close(true_descriptor)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))
            self.assertFalse(
                marker.exists(),
                "an unbound generated CMake/Make input executed")
            self.assertNotIn(os.fsencode(helper), content)
            return content, completed

    def test_plan_is_candidate_derived_and_handles_space_paths(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2 canonical replay paths ") as directory:
            (source, _build, _candidate, _transports,
             content, counts) = self._render_plan(Path(directory))
            text = content.decode("utf-8", errors="strict")
            self.assertIn(shlex.quote(str(source / "leopard2.cpp")), text)
            self.assertNotIn("Makefile2", text)
            self.assertNotIn("DependInfo.cmake", text)
            self.assertNotIn("cmake_depends", text)
            self.assertNotIn("cmake_link_script", text)
            self.assertNotIn("link.txt", text)
            self.assertEqual(counts["recursive_make_count"], 0)
            self.assertEqual(counts["object_count"], 2)
            for line in text.splitlines():
                if line.startswith("CMakeFiles/") and ":" in line:
                    self.assertNotIn(str(source), line)

    def test_plan_ignores_selected_target_added_prerequisite(self) -> None:
        self._run_while_ignoring_generated_input(
            "CMakeFiles/Makefile2",
            b"CMakeFiles/bench_leopard2.dir/all: injected-helper\n"
            b"injected-helper:\n"
            b"\t@HELPER@\n")

    def test_plan_ignores_cmake_env_helper_recipe(self) -> None:
        self._run_while_ignoring_generated_input(
            "CMakeFiles/bench_leopard2.dir/build.make",
            b"CMAKE_COMMAND = /usr/bin/cmake\n"
            b"all:\n"
            b"\t$(CMAKE_COMMAND) -E env @HELPER@\n")

    def test_plan_ignores_dependinfo_execute_process(self) -> None:
        self._run_while_ignoring_generated_input(
            "CMakeFiles/leopard.dir/DependInfo.cmake",
            b"execute_process(COMMAND \"@HELPER@\")\n")

    def test_plan_ignores_exported_ld_preload(self) -> None:
        content, completed = self._run_while_ignoring_generated_input(
            "CMakeFiles/leopard.dir/build.make",
            b"export LD_PRELOAD = /tmp/unretained-leo2.so\n"
            b"all:\n"
            b"\t/bin/true\n")
        self.assertNotIn(b"LD_PRELOAD", content)
        self.assertNotIn(b"LD_PRELOAD", completed.stderr)

    def test_plan_rejects_response_file_indirection(self) -> None:
        def compile_response(candidate):
            candidate["source_object_compile_closure"][0][
                "compile_entry"]["normalized_arguments"].insert(
                    1, "@/tmp/unretained-compile-response")

        def archive_response(candidate):
            candidate["archive_link_commands"][0].append(
                "@/tmp/unretained-archive-response")

        def link_response(candidate):
            candidate["executable_link_command"].insert(
                1, "-Wl,@/tmp/unretained-link-response")

        for name, mutate in (
                ("compile", compile_response),
                ("archive", archive_response),
                ("link", link_response)):
            with self.subTest(location=name):
                with tempfile.TemporaryDirectory(
                        prefix="leo2-canonical-response-") as directory:
                    root = Path(directory)
                    source = root / "source"
                    build = root / "build"
                    source.mkdir()
                    build.mkdir()
                    candidate, transports = canonical_replay_fixture(
                        source, build)
                    mutate(candidate)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "response-file indirection"):
                        provenance._canonical_replay_makefile_bytes(
                            candidate, source.resolve(), build.resolve(),
                            transports)

    def test_plan_rejects_indirect_compile_and_link_options(self) -> None:
        mutations = (
            (
                "compile-wrapper",
                lambda candidate:
                    candidate["source_object_compile_closure"][0][
                        "compile_entry"]["normalized_arguments"].insert(
                            1, "-wrapper"),
            ),
            (
                "compile-specs",
                lambda candidate:
                    candidate["source_object_compile_closure"][0][
                        "compile_entry"]["normalized_arguments"].insert(
                            1, "-specs=/tmp/unretained.specs"),
            ),
            (
                "compile-plugin",
                lambda candidate:
                    candidate["source_object_compile_closure"][0][
                        "compile_entry"]["normalized_arguments"].insert(
                            1, "-fplugin=/tmp/unretained.so"),
            ),
            (
                "compile-prefix",
                lambda candidate:
                    candidate["source_object_compile_closure"][0][
                        "compile_entry"]["normalized_arguments"].insert(
                            1, "-B/tmp/unretained-prefix"),
            ),
            (
                "link-wrapper",
                lambda candidate:
                    candidate["executable_link_command"].insert(
                        1, "-wrapper"),
            ),
            (
                "link-specs",
                lambda candidate:
                    candidate["executable_link_command"].insert(
                        1, "-specs=/tmp/unretained.specs"),
            ),
            (
                "link-plugin",
                lambda candidate:
                    candidate["executable_link_command"].insert(
                        1, "-fplugin=/tmp/unretained.so"),
            ),
            (
                "link-prefix",
                lambda candidate:
                    candidate["executable_link_command"].insert(
                        1, "-B/tmp/unretained-prefix"),
            ),
        )
        for label, mutate in mutations:
            with self.subTest(location=label):
                with tempfile.TemporaryDirectory(
                        prefix="leo2-canonical-indirect-") as directory:
                    root = Path(directory)
                    source = root / "source"
                    build = root / "build"
                    source.mkdir()
                    build.mkdir()
                    candidate, transports = canonical_replay_fixture(
                        source, build)
                    mutate(candidate)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "non-canonical|canonical CMake"):
                        provenance._canonical_replay_makefile_bytes(
                            candidate, source.resolve(), build.resolve(),
                            transports)

    def test_plan_accepts_at_signs_inside_ordinary_paths(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2@canonical@") as directory:
            root = Path(directory)
            source = root / "source@tree"
            build = root / "build@tree"
            source.mkdir()
            build.mkdir()
            candidate, transports = canonical_replay_fixture(source, build)
            content, counts = provenance._canonical_replay_makefile_bytes(
                candidate, source.resolve(), build.resolve(), transports)
            self.assertTrue(content.endswith(b"\n"))
            self.assertEqual(counts["object_count"], 2)

    def test_plan_rejects_source_parent_traversal(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-source-traversal-") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            source.mkdir()
            build.mkdir()
            candidate, transports = canonical_replay_fixture(source, build)
            candidate["source_object_compile_closure"][0][
                "compile_entry"]["normalized_arguments"][-1] = (
                    "${SOURCE_ROOT}/../evil.cpp")
            candidate["tracked_source_manifest"]["files"].append({
                "path": "../evil.cpp",
                "sha256": "0" * 64,
                "size": 0,
                "mode": 0o644,
            })
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError, "escapes|unsafe"):
                provenance._canonical_replay_makefile_bytes(
                    candidate, source.resolve(), build.resolve(), transports)

    def test_plan_rejects_symlinked_build_operands(self) -> None:
        for location in ("parent", "final", "hardlink"):
            with self.subTest(location=location):
                with tempfile.TemporaryDirectory(
                        prefix="leo2-canonical-build-symlink-") as directory:
                    root = Path(directory)
                    source = root / "source"
                    build = root / "build"
                    outside = root / "outside"
                    source.mkdir()
                    build.mkdir()
                    outside.mkdir()
                    candidate, transports = canonical_replay_fixture(
                        source, build)
                    object_path = (
                        build / "CMakeFiles/leopard.dir/leopard2.cpp.o")
                    if location == "parent":
                        object_path.parent.rmdir()
                        object_path.parent.symlink_to(
                            outside, target_is_directory=True)
                    elif location == "final":
                        target = outside / "object"
                        target.write_bytes(b"outside")
                        object_path.symlink_to(target)
                    else:
                        target = outside / "object"
                        target.write_bytes(b"outside")
                        os.link(target, object_path)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "unstable or unsafe|symlink|hard link|non-regular"):
                        with ExitStack() as stack:
                            provenance._retain_canonical_replay_output_topology(
                                candidate, source.resolve(), build.resolve(),
                                stack)

    def test_live_output_topology_rejects_parent_aba_swap(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-output-aba-") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            outside = root / "outside"
            source.mkdir()
            build.mkdir()
            outside.mkdir()
            candidate, transports = canonical_replay_fixture(source, build)
            plan_content, _counts = \
                provenance._canonical_replay_makefile_bytes(
                    candidate, source.resolve(), build.resolve(), transports)
            provenance._write_private_executable(
                build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE,
                plan_content)
            parent = build / "CMakeFiles/leopard.dir"
            saved = build / "CMakeFiles/leopard.dir.saved"
            with ExitStack() as stack:
                retained = \
                    provenance._retain_canonical_replay_output_topology(
                        candidate, source.resolve(), build.resolve(), stack)
                parent.rename(saved)
                parent.symlink_to(outside, target_is_directory=True)
                parent.unlink()
                saved.rename(parent)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "pathname changed"):
                    provenance._verify_canonical_replay_output_topology(
                        retained, require_outputs=False)
                retained["guard"]._close_without_verification()
                retained["entry_guard"]._close_without_verification()

    def test_descriptor_materialization_cannot_modify_transient_victims(
            self) -> None:
        """Every mutable output role is populated through its retained FD."""
        roles = ("object", "archive", "executable", "header", "lock")
        for role in roles:
            with self.subTest(role=role), tempfile.TemporaryDirectory(
                    prefix=f"leo2-canonical-{role}-aba-") as directory:
                root = Path(directory)
                source = root / "source"
                build = root / "build"
                source.mkdir()
                build.mkdir()
                candidate, transports = canonical_replay_fixture(
                    source, build)
                plan_content, _counts = \
                    provenance._canonical_replay_makefile_bytes(
                        candidate, source.resolve(), build.resolve(),
                        transports)
                provenance._write_private_executable(
                    build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE,
                    plan_content)
                victim = root / "victim"
                victim.write_bytes(b"external-victim")
                with ExitStack() as stack:
                    retained = \
                        provenance._retain_canonical_replay_output_topology(
                            candidate, source.resolve(), build.resolve(),
                            stack)

                    def output_role(path: Path) -> str | None:
                        if path.suffix == ".o":
                            return "object"
                        if path.name == "libleopard.a":
                            return "archive"
                        if path.name == candidate["executable_target"]:
                            return "executable"
                        if path.name.endswith(".h.lock"):
                            return "lock"
                        if path.name.endswith(".h"):
                            return "header"
                        return None

                    record = next(
                        item for item in retained["outputs"]
                        if output_role(item["path"]) == role)
                    output = record["path"]
                    saved = output.with_name(output.name + ".reserved")
                    output.rename(saved)
                    os.link(victim, output)
                    descriptor = record["file"].fileno()
                    os.lseek(descriptor, 0, os.SEEK_SET)
                    os.write(descriptor, b"descriptor-only-output")
                    os.ftruncate(descriptor, len(b"descriptor-only-output"))
                    self.assertEqual(
                        victim.read_bytes(), b"external-victim",
                        f"{role} materialization followed a transient link")
                    output.unlink()
                    saved.rename(output)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "pathname changed"):
                        provenance._verify_canonical_replay_output_topology(
                            retained, require_outputs=True)
                    retained["guard"]._close_without_verification()
                    retained["entry_guard"]._close_without_verification()

    def test_sealed_artifacts_bind_every_private_output_identity(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-artifact-binding-") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            source.mkdir()
            build.mkdir()
            candidate, _transports = canonical_replay_fixture(source, build)
            object_digests = ("3" * 64, "4" * 64)
            for item, digest in zip(
                    candidate["source_object_compile_closure"],
                    object_digests):
                item["object"] = {"sha256": digest, "size": 101}
            candidate["archive"] = {"sha256": "5" * 64, "size": 202}
            candidate["executable"] = {"sha256": "6" * 64, "size": 303}
            header_relative = (
                "generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            header = provenance._canonical_replay_attestation_header_bytes(
                candidate)
            identities = {
                "CMakeFiles/leopard.dir/leopard2.cpp.o":
                    (object_digests[0], 101, 0o600),
                (
                    "CMakeFiles/bench_leopard2.dir/"
                    "bench/leopard2/benchmark.cpp.o"
                ): (object_digests[1], 101, 0o600),
                "libleopard.a": ("5" * 64, 202, 0o600),
                "bench_leopard2": ("6" * 64, 303, 0o700),
                header_relative: (
                    hashlib.sha256(header).hexdigest(), len(header), 0o600),
                header_relative + ".lock": (
                    hashlib.sha256(b"").hexdigest(), 0, 0o600),
            }
            artifacts = {
                path: {
                    "path": path,
                    "sha256": digest,
                    "size": size,
                    "mode": stat.S_IFREG | mode,
                }
                for path, (digest, size, mode) in identities.items()
            }
            provenance._validate_sealed_replay_artifact_bindings(
                candidate, artifacts, source.resolve(), build.resolve(),
                label="test")
            mutations = (
                lambda value: value[header_relative].update(
                    sha256="0" * 64),
                lambda value: value[
                    "CMakeFiles/leopard.dir/leopard2.cpp.o"].update(size=102),
                lambda value: value["bench_leopard2"].update(
                    mode=stat.S_IFREG | 0o600),
                lambda value: value.update({
                    "extra": {
                        "path": "extra", "sha256": "0" * 64,
                        "size": 0, "mode": stat.S_IFREG | 0o600,
                    },
                }),
            )
            for mutate in mutations:
                changed = copy.deepcopy(artifacts)
                mutate(changed)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "sealed replay artifact"):
                    provenance._validate_sealed_replay_artifact_bindings(
                        candidate, changed, source.resolve(), build.resolve(),
                        label="test")

    def test_plan_force_rebuilds_precreated_future_artifacts(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2 canonical force ") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            source.mkdir()
            build.mkdir()
            candidate, _transports = canonical_replay_fixture(source, build)
            log = root / "invocations.log"
            helper = root / "transport"
            helper.write_text(
                "#!/bin/sh\n"
                "set -eu\n"
                "printf '%s\\n' \"$*\" >> \"$LEO2_TEST_LOG\"\n"
                "if [ \"${1-}\" = -E ] && [ \"${2-}\" = rm ]; then\n"
                "  shift 3; rm -f -- \"$@\"; exit 0\n"
                "fi\n"
                "if [ \"${1-}\" = qc ]; then\n"
                "  : > \"$2\"; exit 0\n"
                "fi\n"
                "output=\n"
                "previous=\n"
                "for argument in \"$@\"; do\n"
                "  if [ \"$previous\" = -o ]; then output=$argument; break; fi\n"
                "  previous=$argument\n"
                "done\n"
                "if [ -n \"$output\" ]; then\n"
                "  mkdir -p -- \"$(dirname -- \"$output\")\"\n"
                "  : > \"$output\"\n"
                "fi\n",
                encoding="utf-8")
            helper.chmod(0o700)
            descriptor = os.open(
                helper, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
            try:
                transports = {
                    role: f"/proc/self/fd/{descriptor}"
                    for role in ("cxx", "archiver", "ranlib", "cmake", "git")
                }
                content, _counts = \
                    provenance._canonical_replay_makefile_bytes(
                        candidate, source.resolve(), build.resolve(),
                        transports)
                text = content.decode("utf-8", errors="strict")
                self.assertIn(
                    f".PHONY: CMakeFiles/bench_leopard2.dir/replay "
                    f"{provenance._CANONICAL_REPLAY_ATTESTATION_TARGET} "
                    f"{provenance._CANONICAL_REPLAY_FORCE_TARGET}",
                    text)
                for target in (
                        "CMakeFiles/leopard.dir/leopard2.cpp.o",
                        "CMakeFiles/bench_leopard2.dir/"
                        "bench/leopard2/benchmark.cpp.o",
                        "libleopard.a", "bench_leopard2"):
                    rule = next(
                        line for line in text.splitlines()
                        if line.startswith(target + ":"))
                    self.assertIn(
                        provenance._CANONICAL_REPLAY_FORCE_TARGET, rule)

                future = time.time() + 3600
                for artifact in (
                        build / "CMakeFiles/leopard.dir/leopard2.cpp.o",
                        build / "CMakeFiles/bench_leopard2.dir/"
                        "bench/leopard2/benchmark.cpp.o",
                        build / "libleopard.a",
                        build / "bench_leopard2"):
                    artifact.parent.mkdir(parents=True, exist_ok=True)
                    artifact.write_bytes(b"precreated")
                    os.utime(artifact, (future, future))
                plan = build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE
                provenance._write_private_executable(plan, content)
                completed = subprocess.run(
                    [
                        str(tool_path("make")), "-C", str(build),
                        "-f", str(plan), "-j1", "--",
                        "CMakeFiles/bench_leopard2.dir/replay",
                    ],
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env={
                        **provenance.GIT_ENVIRONMENT,
                        "LEO2_TEST_LOG": str(log),
                    },
                    pass_fds=(descriptor,), timeout=30, check=False)
            finally:
                os.close(descriptor)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))
            invocations = log.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                sum(" -c " in (" " + line + " ") for line in invocations), 2)
            self.assertEqual(
                sum(line.startswith("qc libleopard.a ") for line in invocations),
                1)
            self.assertEqual(
                sum(line == "libleopard.a" for line in invocations), 1)
            self.assertEqual(
                sum(line.startswith("-E rm -f libleopard.a")
                    for line in invocations), 0)
            self.assertEqual(
                sum(" -o bench_leopard2 " in (" " + line + " ")
                    for line in invocations), 1)

    def test_private_plan_parent_symlink_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-parent-symlink-") as directory:
            root = Path(directory)
            build = root / "build"
            outside = root / "outside"
            build.mkdir()
            outside.mkdir()
            (build / "CMakeFiles").symlink_to(
                outside, target_is_directory=True)
            destination = (
                build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "parent follows a symlink"):
                provenance._write_private_executable(
                    destination, b"# sealed plan\n")
            self.assertFalse(
                (outside / destination.name).exists(),
                "symlinked plan parent received replay bytes")

    def test_manifest_inventory_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-inventory-") as directory:
            (source, build, candidate, transports,
             content, counts) = self._render_plan(Path(directory))
            plan_path = build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE
            provenance._write_private_executable(plan_path, content)
            with provenance._RetainedFileSnapshot(
                    plan_path, "canonical replay inventory test plan"
                    ) as snapshot:
                manifest = provenance._canonical_replay_plan_manifest(
                    candidate, plan_path, content, snapshot, counts)
                descriptor = \
                    provenance._validate_canonical_replay_plan_manifest(
                        manifest, candidate, source, build, transports,
                        label="test")
                self.assertEqual(
                    descriptor, snapshot.executable_descriptor)

                def coherent_inventory(
                    mutation,
                ) -> dict[str, object]:
                    changed = copy.deepcopy(manifest)
                    mutation(changed["retained_inputs"])
                    changed["file_count"] = len(changed["retained_inputs"])
                    changed["total_bytes"] = sum(
                        record["size"]
                        for record in changed["retained_inputs"])
                    changed["inventory_sha256"] = \
                        provenance._canonical_replay_inventory_sha256(
                            changed["retained_inputs"])
                    return changed

                mutations: dict[str, dict[str, object]] = {
                    "omission": coherent_inventory(
                        lambda records: records.pop()),
                    "addition": coherent_inventory(
                        lambda records: records.append({
                            "role": "generated-plan",
                            "root": "build",
                            "path": "CMakeFiles/forged.make",
                            "sha256": "0" * 64,
                            "size": 0,
                            "mode": 0o700,
                        })),
                    "reorder": coherent_inventory(
                        lambda records: records.reverse()),
                    "file-count": copy.deepcopy(manifest),
                    "recursive-make-count": copy.deepcopy(manifest),
                    "command-count": copy.deepcopy(manifest),
                }
                mutations["file-count"]["file_count"] += 1
                mutations["recursive-make-count"]["counts"][
                    "recursive_make_count"] = 1
                mutations["command-count"]["counts"]["command_count"] += 1
                validator = (
                    provenance._validate_canonical_replay_plan_manifest)
                for name, changed in mutations.items():
                    with self.subTest(mutation=name):
                        with self.assertRaisesRegex(
                                provenance.BuildProvenanceError,
                                "retained input inventory differs"):
                            validator(
                                changed, candidate, source, build,
                                transports, label="test")

    def test_manifest_validation_is_offline_after_roots_are_deleted(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-offline-") as directory:
            root = Path(directory)
            (source_root, build_root, candidate, transports,
             content, counts) = self._render_plan(root)
            plan = build_root / provenance._CANONICAL_REPLAY_PLAN_RELATIVE
            provenance._write_private_executable(plan, content)
            with provenance._RetainedFileSnapshot(
                    plan, "offline canonical plan") as snapshot:
                manifest = provenance._canonical_replay_plan_manifest(
                    candidate, plan, content, snapshot, counts)
                expected_descriptor = snapshot.executable_descriptor
        self.assertFalse(source_root.exists())
        self.assertFalse(build_root.exists())
        self.assertEqual(
            provenance._validate_canonical_replay_plan_manifest(
                manifest, candidate, source_root, build_root, transports,
                label="offline test"),
            expected_descriptor)

    def test_manifest_validation_rejects_tampering_under_python_o(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-canonical-optimized-python-") as directory:
            root = Path(directory)
            (source, build, candidate, transports,
             content, counts) = self._render_plan(root)
            plan_path = build / provenance._CANONICAL_REPLAY_PLAN_RELATIVE
            provenance._write_private_executable(plan_path, content)
            with provenance._RetainedFileSnapshot(
                    plan_path, "optimized canonical replay test plan"
                    ) as snapshot:
                manifest = provenance._canonical_replay_plan_manifest(
                    candidate, plan_path, content, snapshot, counts)
            changed = copy.deepcopy(manifest)
            changed["counts"]["recursive_make_count"] = 1
            payload = root / "payload.json"
            payload.write_text(json.dumps({
                "candidate": candidate,
                "manifest": manifest,
                "changed": changed,
                "source": str(source),
                "build": str(build),
                "transports": transports,
            }), encoding="utf-8")
            program = (
                "import json, pathlib, sys\n"
                "sys.path.insert(0, sys.argv[1])\n"
                "import leopard2_build_provenance as p\n"
                "value=json.loads(pathlib.Path(sys.argv[2]).read_text("
                "encoding='utf-8'))\n"
                "arguments=(value['candidate'],"
                "pathlib.Path(value['source']),"
                "pathlib.Path(value['build']),value['transports'])\n"
                "try:\n"
                " p._validate_canonical_replay_plan_manifest("
                "value['manifest'],*arguments,label='python-o-valid')\n"
                "except p.BuildProvenanceError as error:\n"
                " print(error,file=sys.stderr); raise SystemExit(3)\n"
                "try:\n"
                " p._validate_canonical_replay_plan_manifest("
                "value['changed'],*arguments,label='python-o-tampered')\n"
                "except p.BuildProvenanceError:\n"
                " raise SystemExit(0)\n"
                "raise SystemExit(4)\n")
            environment = dict(os.environ)
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            completed = subprocess.run(
                [sys.executable, "-O", "-c", program,
                 str(SOURCE_ROOT / "tools"), str(payload)],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=environment, timeout=30, check=False)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))


class ReproducibleCompilerReplayTests(unittest.TestCase):
    def test_recursive_make_uses_sealed_image_after_same_inode_mutation(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-recursive-make-") as directory:
            root = Path(directory)
            mutable_make = root / "make"
            original = tool_path("make").read_bytes()
            alternate = tool_path("false").read_bytes()
            mutable_make.write_bytes(original)
            mutable_make.chmod(0o700)
            inode = mutable_make.stat().st_ino
            (root / "Makefile").write_text(
                "all:\n"
                "\t$(MAKE) --no-print-directory -f child.mk child\n",
                encoding="ascii")
            (root / "child.mk").write_text(
                "child:\n"
                "\t$(MAKE) --no-print-directory -f grandchild.mk "
                "grandchild\n",
                encoding="ascii")
            (root / "grandchild.mk").write_text(
                "grandchild:\n"
                "\t@printf 'recursive-sealed\\n'\n",
                encoding="ascii")

            make_snapshot = provenance._RetainedFileSnapshot(
                mutable_make, "recursive Make")
            try:
                with provenance._RetainedFileSnapshot(
                        tool_path("bash"), "recursive Make Bash") as bash, \
                        provenance._RetainedFileSnapshot(
                            tool_path("sh"),
                            "recursive Make command shell") as shell:
                    wrapper_path = root / "sealed-make"
                    provenance._write_replay_exec_wrapper(
                        wrapper_path, bash.executable_descriptor,
                        str(mutable_make),
                        make_snapshot.executable_descriptor)
                    with provenance._RetainedFileSnapshot(
                            wrapper_path, "recursive Make wrapper"
                            ) as wrapper:
                        with mutable_make.open("r+b") as stream:
                            stream.write(alternate)
                            stream.truncate()
                            stream.flush()
                            os.fsync(stream.fileno())
                        self.assertEqual(mutable_make.stat().st_ino, inode)
                        output = provenance._run(
                            [str(mutable_make), "-C", str(root),
                             "--no-print-directory",
                             f"SHELL=/proc/self/fd/"
                             f"{shell.executable_descriptor}",
                             f"MAKE=/proc/self/fd/"
                             f"{wrapper.executable_descriptor}",
                             "--", "all"],
                            "sealed recursive Make",
                            inherited_descriptors=(
                                bash.executable_descriptor,
                                shell.executable_descriptor,
                                make_snapshot.executable_descriptor,
                                wrapper.executable_descriptor,
                            ),
                            executable_descriptor=
                            make_snapshot.executable_descriptor,
                            timeout=10)
                        expected_command = (
                            f"/proc/self/fd/"
                            f"{wrapper.executable_descriptor} "
                            "--no-print-directory -f child.mk child\n"
                            f"/proc/self/fd/"
                            f"{wrapper.executable_descriptor} "
                            "--no-print-directory -f grandchild.mk "
                            "grandchild\n"
                            "recursive-sealed\n").encode("ascii")
                        self.assertEqual(output, expected_command)
            finally:
                # The mutation is intentional and must not obscure whether
                # both outer and recursive Make used the sealed snapshot.
                make_snapshot._close_without_verification()

    def test_clean_configure_preserves_distinct_c_driver(self) -> None:
        c_compiler = distinct_c_driver()
        cxx_compiler = Path(
            shutil.which("c++", path="/usr/bin:/bin") or
            "/usr/bin/c++").resolve(strict=True)
        self.assertNotEqual(c_compiler, cxx_compiler)
        cache = {
            "CMAKE_AR": str(tool_path("ar")),
            "CMAKE_C_COMPILER": str(c_compiler),
            "CMAKE_CXX_COMPILER": str(cxx_compiler),
            "CMAKE_CXX_FLAGS": "",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "CMAKE_LINKER": str(tool_path("ld")),
            "CMAKE_MAKE_PROGRAM": str(tool_path("make")),
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BENCHMARK_GIT_EXECUTABLE": str(tool_path("git")),
            "LEO2_BUILD_ALLK_DIAGNOSTIC": "OFF",
            "LEO2_BUILD_TESTS": "ON",
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
            "LEOPARD_ENABLE_GF16": "ON",
            "LEO2_FLAG_FALIGN_FUNCTIONS_64": "1",
            "LEO2_FLAG_MAVX2": "1",
            "LEO2_FLAG_MNO_AVX512F": "1",
            "LEO2_LOCATOR_GIT_EXECUTABLE": str(tool_path("git")),
            "CMAKE_RANLIB": str(tool_path("ranlib")),
        }
        missing_c_compiler = dict(cache)
        del missing_c_compiler["CMAKE_C_COMPILER"]
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "lacks a C compiler"):
            provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path("/tmp/unused-build"),
                missing_c_compiler)
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-cc-") as directory:
            build = Path(directory) / "build"
            configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, build, cache)
            self.assertNotIn(f"CC={c_compiler}", configure)
            self.assertNotIn(f"CXX={cxx_compiler}", configure)
            self.assertIn(
                "-DLEO2_BENCHMARK_GIT_EXECUTABLE:FILEPATH=" +
                str(tool_path("git")), configure)
            self.assertIn(
                "-DLEO2_LOCATOR_GIT_EXECUTABLE:FILEPATH=" +
                str(tool_path("git")), configure)
            self.assertIn(
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK:BOOL=OFF",
                configure)
            v8_cache = dict(cache)
            for selector in (
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
                v8_cache.pop(selector)
            v8_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8
            v8_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v8-build", v8_cache)
            self.assertIn(
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT:BOOL=ON",
                v8_configure)
            self.assertFalse(any(
                any(argument.startswith(f"-D{selector}:")
                    for selector in (
                        "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                        "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"))
                for argument in v8_configure))
            v7_cache = dict(cache)
            for selector in (
                    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
                v7_cache.pop(selector)
            v7_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7
            v7_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v7-build", v7_cache)
            self.assertFalse(any(
                argument.startswith(
                    "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT:")
                for argument in v7_configure))
            self.assertIn(
                "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED:BOOL=ON",
                configure)
            self.assertIn(
                "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED:BOOL=OFF",
                configure)
            v6_only = {
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
            }
            for selector in v6_only:
                self.assertTrue(any(
                    argument.startswith(f"-D{selector}:BOOL=")
                    for argument in configure))
            v5_cache = dict(cache)
            v5_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5
            v5_cache["LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED"] = "OFF"
            for selector in (
                    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
                v5_cache.pop(selector)
            for selector in v6_only:
                v5_cache.pop(selector)
            v5_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v5-build", v5_cache)
            self.assertFalse(any(
                any(argument.startswith(f"-D{selector}:")
                    for selector in v6_only)
                for argument in v5_configure))
            v5_only = {
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED",
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
            }
            v4_cache = dict(v5_cache)
            v4_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4
            for selector in v5_only:
                v4_cache.pop(selector)
            v4_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v4-build", v4_cache)
            self.assertFalse(any(
                any(argument.startswith(f"-D{selector}:")
                    for selector in v5_only)
                for argument in v4_configure))
            current_only = {
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE",
                "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING",
                "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING",
                "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
            }
            v3_cache = dict(v4_cache)
            v3_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3
            v3_cache["LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = "OFF"
            for selector in current_only:
                v3_cache.pop(selector)
            v3_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v3-build", v3_cache)
            self.assertIn(
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT:BOOL=OFF",
                v3_configure)
            self.assertFalse(any(
                any(argument.startswith(f"-D{selector}:")
                    for selector in current_only)
                for argument in v3_configure))

            v2_cache = dict(v3_cache)
            v2_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2
            v2_cache.pop("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT")
            v2_configure = provenance._reproducible_configure_argv(
                SOURCE_ROOT, Path(directory) / "v2-build", v2_cache)
            self.assertFalse(any(
                argument.startswith(
                    "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT:")
                for argument in v2_configure))
            configure_environment = dict(provenance.GIT_ENVIRONMENT)
            configure_environment.update({
                "CC": str(c_compiler), "CXX": str(cxx_compiler)})
            completed = subprocess.run(
                configure, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=configure_environment,
                timeout=300, check=False)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))

            cache_bytes = (build / "CMakeCache.txt").read_bytes()
            parsed = provenance.parse_cmake_cache(cache_bytes)
            self.assertEqual(parsed["CMAKE_AR"], cache["CMAKE_AR"])
            self.assertEqual(parsed["CMAKE_C_COMPILER"], str(c_compiler))
            self.assertEqual(parsed["CMAKE_CXX_COMPILER"], str(cxx_compiler))
            self.assertEqual(parsed["CMAKE_RANLIB"], cache["CMAKE_RANLIB"])
            cache_text = cache_bytes.decode("utf-8", errors="strict")
            self.assertIn(
                f"CMAKE_C_COMPILER:FILEPATH={c_compiler}\n", cache_text)
            self.assertIn(
                f"CMAKE_CXX_COMPILER:FILEPATH={cxx_compiler}\n", cache_text)

            commands = json.loads(
                (build / "compile_commands.json").read_text(
                    encoding="utf-8"))
            c_source = (
                SOURCE_ROOT / "tests/leopard2/test_codec_options_abi.c"
            ).resolve(strict=True)
            c_entries = [
                entry for entry in commands
                if Path(entry.get("file", "")).resolve(strict=False) ==
                c_source
            ]
            self.assertEqual(len(c_entries), 1)
            c_tokens = provenance._compile_tokens(c_entries[0])
            self.assertEqual(
                Path(c_tokens[0]).resolve(strict=True), c_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "another command driver"):
                provenance._require_compile_driver(
                    [str(cxx_compiler), "-c", str(c_source)],
                    c_source, c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "another command driver"):
                provenance._require_compile_driver(
                    ["/usr/bin/env", str(c_compiler), "-c", str(c_source)],
                    c_source, c_compiler, cxx_compiler)

    def test_replay_git_transport_preserves_configured_identity(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-git-transport-") as directory:
            build = Path(directory) / "build"
            recipe = (
                build /
                "CMakeFiles/"
                "leopard2_benchmark_source_attestation_refresh.dir/"
                "build.make")
            recipe.parent.mkdir(parents=True)
            identity_git = tool_path("git")
            transport_git = Path(directory) / "replay-git"
            transport_git.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            transport_git.chmod(0o700)
            configured_hash = (
                b"LEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git\n")
            original = (
                b"command -DLEO2_GIT_EXECUTABLE=" +
                str(identity_git).encode("utf-8") + b" rest\n")
            recipe.write_bytes(original)
            sidecar = build / "configuration.txt"
            sidecar.write_bytes(configured_hash)
            provenance._retarget_replay_git_transport(
                build, identity_git, transport_git)
            self.assertEqual(sidecar.read_bytes(), configured_hash)
            rewritten = recipe.read_bytes()
            self.assertNotIn(
                ("-DLEO2_GIT_EXECUTABLE=" +
                 str(identity_git)).encode("utf-8"), rewritten)
            self.assertIn(
                ("-DLEO2_GIT_EXECUTABLE=" +
                 str(transport_git)).encode("utf-8"), rewritten)

    def test_make_variable_shell_injection_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-make-value-injection-") as directory:
            build = Path(directory)
            helper = build / "unretained-helper"
            marker = build / "unexpected-marker"
            helper.write_text(
                "#!/bin/sh\n"
                f"printf executed > {shlex.quote(str(marker))}\n",
                encoding="utf-8")
            helper.chmod(0o700)
            cmake = tool_path("cmake")
            makefile = build / "Makefile"
            payload = (
                f"CMAKE_COMMAND = {cmake}\n"
                f"ARGS = ; {helper}\n"
                "all:\n"
                "\t$(CMAKE_COMMAND) -E echo retained $(ARGS)\n"
            ).encode("utf-8")
            makefile.write_bytes(payload)
            completed = subprocess.run(
                [str(tool_path("make")), "-f", str(makefile), "all"],
                cwd=build, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=10, check=False)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))
            self.assertTrue(
                marker.is_file(),
                "live GNU Make fixture did not demonstrate value injection")
            marker.unlink()
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "external-only variable"):
                provenance._validate_replay_make_include_closure(
                    makefile, build, payload, {"Makefile"}, str(cmake))
            self.assertFalse(marker.exists())

    def test_replay_compile_variables_are_bound_to_audited_argv(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-compile-variable-") as directory:
            root = Path(directory)
            build = root / "build"
            source = root / "source"
            flags = build / "CMakeFiles/target.dir/flags.make"
            build_make = build / "CMakeFiles/target.dir/build.make"
            cache = build / "CMakeCache.txt"
            flags.parent.mkdir(parents=True)
            source.mkdir()
            source_file = source / "input.cpp"
            source_file.write_text("int value;\n", encoding="ascii")
            compiler = str(tool_path("c++"))
            output = "CMakeFiles/target.dir/input.cpp.o"
            arguments = [
                compiler, "-DNDEBUG", f"-I{source}", "-Wall",
                "-o", output, "-c", str(source_file),
            ]
            commands = [{
                "directory": str(build),
                "file": str(source_file),
                "output": output,
                "arguments": arguments,
            }]
            expected_closure = [{
                "role": "archive",
                "compile_entry": {
                    "normalized_arguments":
                        provenance._normalize_build_argv(
                            arguments, build, source),
                },
            }]
            canonical_flags = (
                b"CXX_DEFINES = -DNDEBUG\n"
                + f"CXX_INCLUDES = -I{source}\n".encode("utf-8")
                + b"CXX_FLAGS = -Wall\n")
            canonical_controls = (
                f"CMAKE_COMMAND = {tool_path('cmake')}\n"
                f"CMAKE_SOURCE_DIR = {source}\n"
                f"CMAKE_BINARY_DIR = {build}\n").encode("utf-8")
            canonical_recipe = (
                f"{output}: CMakeFiles/target.dir/flags.make\n"
                f"{output}: {source_file}\n"
                f"{output}: CMakeFiles/target.dir/compiler_depend.ts\n"
                f"\t@$(CMAKE_COMMAND) -E cmake_echo_color "
                f"\"--switch=$(COLOR)\" --green "
                f"--progress-dir={build}/CMakeFiles "
                f"--progress-num=$(CMAKE_PROGRESS_1) "
                f"\"Building CXX object {output}\"\n"
                f"\t{compiler} $(CXX_DEFINES) $(CXX_INCLUDES) "
                f"$(CXX_FLAGS) -MD -MT {output} -MF {output}.d "
                f"-o {output} -c {source_file}\n").encode("utf-8")
            canonical_build_make = canonical_controls + canonical_recipe
            canonical_cache = (
                f"CMAKE_HOME_DIRECTORY:INTERNAL={source}\n").encode(
                    "utf-8")

            def validate(
                *,
                flags_content: bytes = canonical_flags,
                command_records: object = commands,
                build_make_content: bytes = canonical_build_make,
            ) -> None:
                recipe_contents = {
                    flags: flags_content,
                    build_make: build_make_content,
                    cache: canonical_cache,
                }
                provenance._validate_replay_compile_variable_closure(
                    build, source,
                    json.dumps(
                        command_records, sort_keys=True,
                        separators=(",", ":")).encode("ascii"),
                    recipe_contents, expected_closure,
                    str(tool_path("cmake")))

            validate()
            for injected in (
                    b"CXX_FLAGS = -Wall -B /tmp/unretained-prefix\n",
                    b"CXX_FLAGS = -Wall -fplugin=/tmp/unretained.so\n",
                    b"CXX_FLAGS = -Wall -include /tmp/unretained.h\n",
                    b"CXX_FLAGS = -Wall @/tmp/unretained-response\n",
                    b"CXX_FLAGS = -Wall -I~/unretained\n",
                    b"CXX_FLAGS = -Wall -DVALUE={unretained}\n",
                    b"CXX_FLAGS = -Wall # ignored remainder\n",
                    b"CXX_FLAGS = -Wextra\n"):
                with self.subTest(injected=injected):
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "(?:flags.make variables differ|"
                            "response-file indirection|"
                            "shell code or indirection)"):
                        validate(
                            flags_content=canonical_flags.replace(
                                b"CXX_FLAGS = -Wall\n", injected))

            mismatched_commands = json.loads(json.dumps(commands))
            mismatched_commands[0]["arguments"].insert(3, "-Wextra")
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "generated compile argv differs"):
                validate(command_records=mismatched_commands)
            recipe_injections = (
                b"-B /tmp/unretained-prefix",
                b"-fplugin=/tmp/unretained.so",
                b"-include /tmp/unretained.h",
                b"@/tmp/unretained-response",
            )
            for injected in recipe_injections:
                with self.subTest(recipe_injected=injected):
                    mutated = canonical_build_make.replace(
                        b" -MD -MT ", b" " + injected + b" -MD -MT ")
                    if injected.startswith(b"-include "):
                        # The generic Make syntax boundary deliberately
                        # permits ordinary compiler operands.  This reproduces
                        # why the semantic object-recipe binding is required.
                        provenance._validate_replay_make_include_closure(
                            build_make, build, mutated,
                            {"CMakeFiles/target.dir/build.make"},
                            str(tool_path("cmake")), str(source))
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "production compile recipe differs"):
                        validate(build_make_content=mutated)
            extra_recipe = (
                canonical_build_make +
                f"\t{compiler} -fplugin=/tmp/unretained.so\n".encode(
                    "utf-8"))
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "non-canonical recipe count"):
                validate(build_make_content=extra_recipe)
            inline_recipe = canonical_build_make.replace(
                f"{output}: CMakeFiles/target.dir/flags.make\n".encode(
                    "utf-8"),
                (f"{output}: CMakeFiles/target.dir/flags.make ; "
                 f"{compiler} -fplugin=/tmp/unretained.so\n").encode(
                     "utf-8"))
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "dependency rules differ"):
                validate(build_make_content=inline_recipe)

            injected_header = root / "unretained.h"
            injected_header.write_text(
                "#define LEO2_RECIPE_INJECTION_OBSERVED 1\n",
                encoding="ascii")
            injected_source = root / "injected.cpp"
            injected_source.write_text(
                "#ifndef LEO2_RECIPE_INJECTION_OBSERVED\n"
                "#error direct recipe option was not consumed\n"
                "#endif\n"
                "int injected_value;\n",
                encoding="ascii")
            injected_object = build / "injected.o"
            live_make = root / "live-injected.make"
            live_make.write_text(
                "all:\n"
                f"\t{compiler} -include {injected_header} "
                f"-o {injected_object} -c {injected_source}\n",
                encoding="utf-8")
            completed = subprocess.run(
                [str(tool_path("make")), "-f", str(live_make), "all"],
                cwd=build, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30, check=False)
            self.assertEqual(
                completed.returncode, 0,
                completed.stderr.decode("utf-8", errors="replace"))
            self.assertTrue(
                injected_object.is_file(),
                "live Make fixture did not execute the direct recipe option")
            for original, replacement in (
                    (f"CMAKE_SOURCE_DIR = {source}\n".encode("utf-8"),
                     b"CMAKE_SOURCE_DIR = /tmp/unretained-source\n"),
                    (f"CMAKE_BINARY_DIR = {build}\n".encode("utf-8"),
                     b"CMAKE_BINARY_DIR = /tmp/unretained-build\n"),
                    (f"CMAKE_COMMAND = {tool_path('cmake')}\n".encode(
                        "utf-8"),
                     b"CMAKE_COMMAND = /tmp/unretained-cmake\n")):
                with self.subTest(replacement=replacement):
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "control variable differs"):
                        validate(build_make_content=
                                 canonical_build_make.replace(
                                     original, replacement))

            flags.write_bytes(canonical_flags)
            build_make.write_bytes(canonical_build_make)
            cache.write_bytes(canonical_cache)
            (build / "compile_commands.json").write_text(
                json.dumps(commands, sort_keys=True, separators=(",", ":")),
                encoding="ascii")
            required_inputs = {
                "Makefile": canonical_controls,
                "CMakeFiles/Makefile2": b"# retained target graph\n",
                "CMakeFiles/Makefile.cmake":
                    b"# retained CMake inputs\n",
                "CMakeFiles/cmake.check_cache":
                    b"# retained cache check\n",
            }
            for relative, content in required_inputs.items():
                path = build / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
            replacements = {
                "cmake": (
                    str(tool_path("cmake")), "/proc/self/fd/101"),
            }
            with provenance._TemporaryReplayRecipeTransport(
                    build, replacements, set(),
                    expected_compile_closure=expected_closure,
                    expected_source=source) as transport:
                semantic = transport.manifest[
                    "semantic_compile_closure"]
                self.assertEqual(
                    semantic["schema"],
                    "leopard2-replay-compile-variable-closure/v2")
                self.assertEqual(semantic["object_count"], 1)
                self.assertEqual(semantic["target_count"], 1)
                self.assertEqual(
                    semantic["candidate_compile_closure_sha256"],
                    provenance._replay_candidate_compile_closure_sha256(
                        expected_closure))
                self.assertEqual(
                    semantic["candidate_target_variables_sha256"],
                    provenance._replay_candidate_target_variables_sha256(
                        expected_closure))
                self.assertEqual(
                    transport.manifest["semantic_inputs"][0]["path"],
                    "compile_commands.json")
                coherently_replaced = copy.deepcopy(semantic)
                coherently_replaced["compile_commands_sha256"] = "a" * 64
                coherently_replaced[
                    "candidate_compile_closure_sha256"] = "b" * 64
                coherently_replaced[
                    "candidate_target_variables_sha256"] = "c" * 64
                replaced_input = copy.deepcopy(
                    transport.manifest["semantic_inputs"][0])
                replaced_input["sha256"] = "a" * 64
                self.assertEqual(
                    coherently_replaced["compile_commands_sha256"],
                    replaced_input["sha256"])
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "candidate compile closure digest differs"):
                    binding = provenance.\
                        _require_replay_candidate_compile_closure_binding
                    binding(
                        coherently_replaced, expected_closure,
                        label="coherently replaced proof")

    def test_replay_recipe_transport_is_exact_and_temporary(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-recipe-transport-") as directory:
            build = Path(directory) / "build"
            link = build / "CMakeFiles/target.dir/link.txt"
            link.parent.mkdir(parents=True)
            progress = build / "CMakeFiles/target.dir/progress.make"
            build_make = build / "CMakeFiles/target.dir/build.make"
            depend_make = build / "CMakeFiles/target.dir/depend.make"
            compiler_depend_make = (
                build / "CMakeFiles/target.dir/compiler_depend.make")
            makefile_cmake = build / "CMakeFiles/Makefile.cmake"
            check_cache = build / "CMakeFiles/cmake.check_cache"
            cache = build / "CMakeCache.txt"
            makefile2 = build / "CMakeFiles/Makefile2"
            makefile = build / "Makefile"
            original_makefile = (
                b"# Command-line flag to silence nested $(MAKE).\n"
                b"/usr/bin/cmake --build .\n"
                b"/usr/bin/c++ -c input.cpp\n"
                b"/usr/bin/c++-unrelated -c other.cpp\n"
                b"\t$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 target\n")
            original_makefile2 = b"# retained recursive Make input\n"
            original_link = (
                b"/usr/bin/ar qc libleopard.a object.o\n"
                b"/usr/bin/ranlib libleopard.a\n"
                b"/usr/bin/git status\n")
            original_progress = b"CMAKE_PROGRESS_1 = 1\n"
            original_build_make = (
                b"include CMakeFiles/target.dir/depend.make\n"
                b"include CMakeFiles/target.dir/compiler_depend.make\n")
            original_makefile_cmake = b"set(CMAKE_DEPENDS_GENERATOR \"Unix Makefiles\")\n"
            makefile.write_bytes(original_makefile)
            makefile2.write_bytes(original_makefile2)
            link.write_bytes(original_link)
            progress.write_bytes(original_progress)
            build_make.write_bytes(original_build_make)
            depend_make.write_bytes(
                b"$(shell touch injected-dependency-marker)\n")
            compiler_depend_make.write_bytes(
                b"$(shell touch injected-compiler-dependency-marker)\n")
            makefile_cmake.write_bytes(original_makefile_cmake)
            check_cache.write_bytes(b"# retained cache check\n")
            cache.write_bytes(b"# retained CMake cache\n")
            replacements = {
                "cmake": ("/usr/bin/cmake", "/proc/self/fd/101"),
                "cxx": ("/usr/bin/c++", "/proc/self/fd/102"),
                "ar": ("/usr/bin/ar", "/proc/self/fd/103"),
                "ranlib": ("/usr/bin/ranlib", "/proc/self/fd/104"),
                "git": ("/usr/bin/git", "/proc/self/fd/105"),
                "make": ("/usr/bin/make", "/proc/self/fd/106"),
                "dependency-includes": (
                    "CMakeFiles/**/{compiler_,}depend.make",
                    "/proc/self/fd/107"),
            }
            with provenance._TemporaryReplayRecipeTransport(
                    build, replacements, set(replacements)) as transport:
                rewritten_makefile = makefile.read_bytes()
                rewritten_link = link.read_bytes()
                rewritten_build_make = build_make.read_bytes()
                self.assertIn(b"/proc/self/fd/101", rewritten_makefile)
                self.assertIn(b"/proc/self/fd/102", rewritten_makefile)
                self.assertIn(
                    b"/usr/bin/c++-unrelated", rewritten_makefile)
                self.assertIn(b"/proc/self/fd/103", rewritten_link)
                self.assertIn(b"/proc/self/fd/104", rewritten_link)
                self.assertIn(b"/proc/self/fd/105", rewritten_link)
                self.assertIn(
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 target",
                    rewritten_makefile)
                self.assertEqual(
                    rewritten_build_make,
                    b"include /proc/self/fd/107\n"
                    b"include /proc/self/fd/107\n")
                retained = {
                    record["path"]: record["rewritten"]
                    for record in transport.manifest["retained_recipes"]
                }
                self.assertEqual(
                    retained["CMakeFiles/target.dir/progress.make"], False)
                self.assertEqual(
                    retained["CMakeFiles/Makefile.cmake"], False)
                self.assertNotIn(
                    "CMakeFiles/target.dir/depend.make", retained)
                self.assertNotIn(
                    "CMakeFiles/target.dir/compiler_depend.make", retained)
                self.assertTrue(
                    all(record["replacement_counts"]
                        for record in transport.manifest["rewrites"]))
            self.assertEqual(makefile.read_bytes(), original_makefile)
            self.assertEqual(makefile2.read_bytes(), original_makefile2)
            self.assertEqual(link.read_bytes(), original_link)
            self.assertEqual(progress.read_bytes(), original_progress)
            self.assertEqual(build_make.read_bytes(), original_build_make)

            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "recipe changed while retained"):
                with provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements)):
                    makefile.write_bytes(
                        makefile.read_bytes() +
                        b"\t/usr/bin/true # injected after retention\n")
            self.assertEqual(makefile.read_bytes(), original_makefile)
            self.assertEqual(link.read_bytes(), original_link)

            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "recipe changed while retained"):
                with provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements)):
                    progress.write_bytes(
                        original_progress +
                        b"$(shell touch unexpected-marker)\n")
            progress.write_bytes(original_progress)

            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "recipe changed while retained"):
                with provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements)):
                    makefile_cmake.write_bytes(
                        original_makefile_cmake +
                        b"execute_process(COMMAND false)\n")
            makefile_cmake.write_bytes(original_makefile_cmake)

            bypass_build = (
                build / "CMakeFiles/bypass.dir/build.make")
            bypass_dependency = (
                build / "CMakeFiles/bypass.dir/depend.make")
            bypass_build.parent.mkdir()
            bypass_dependency.write_bytes(
                b"$(shell touch dependency-bypass-marker)\n")
            (bypass_build.parent / "extra.mk").write_bytes(
                b"$(shell touch dependency-bypass-marker)\n")
            for bypass in (
                    b"dependency_file := "
                    b"CMakeFiles/bypass.dir/depend.make\n"
                    b"include $(dependency_file)\n",
                    b"include CMakeFiles/bypass.dir/depend.mak?\n",
                    b"include CMakeFiles/bypass.dir/depen$()d.make\n",
                    b"include CMakeFiles/bypass.dir/extra.mk\n"):
                with self.subTest(bypass=bypass):
                    bypass_build.write_bytes(bypass)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "(?:dynamic include|non-canonical include|"
                            "(?:computed|non-canonical) variable reference)"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            makefile2.write_bytes(
                b"\t$(MAKE) $(MAKESILENT) -f "
                b"CMakeFiles/bypass.dir/extra.mk target\n")
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "recursive Make selects an unretained recipe"):
                with provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements)):
                    pass
            makefile2.write_bytes(original_makefile2)

            for noncanonical_make in (
                    b"\t/usr/bin/make target\n",
                    b"\tmake -f CMakeFiles/bypass.dir/extra.mk target\n",
                    b"\t/usr/bin/gmake "
                    b"-CCMakeFiles/bypass.dir target\n",
                    b"\t/usr/bin/env make "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"\t\"/usr/bin/make\" target\n",
                    b"\tcd . && make "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"\t/bin/sh -c \"make "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\"\n",
                    b"\teval make "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"\tprintf x | xargs make "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"\tX=make; $$X "
                    b"-fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"M = /usr/bin/make\n"
                    b"\t$(M) -f CMakeFiles/bypass.dir/extra.mk target\n",
                    b"\t$(UNTRUSTED_MAKE) target\n",
                    b"\t$(UNTRUSTED_MAKE)"
                    b" -fCMakeFiles/bypass.dir/extra.mk target\n",
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 "
                    b"-fCMakeFiles/bypass.dir/extra.mk\n"):
                with self.subTest(noncanonical_make=noncanonical_make):
                    makefile2.write_bytes(noncanonical_make)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "non-canonical recursive Make command"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            makefile2.write_bytes(original_makefile2)

            for value_injection in (
                    b"ARGS = ; /tmp/unretained-helper\n"
                    b"\t$(CMAKE_COMMAND) -E echo ok $(ARGS)\n",
                    b"COLOR = && /tmp/unretained-helper\n"
                    b"\t$(CMAKE_COMMAND) -E echo $(COLOR)\n",
                    b"VERBOSE = | /tmp/unretained-helper\n"
                    b"\t$(CMAKE_COMMAND) -E echo $(VERBOSE)\n",
                    b"CXX_FLAGS += ; /tmp/unretained-helper\n"
                    b"\t/usr/bin/c++ $(CXX_FLAGS) -c source.cpp\n",
                    b"target: CXX_FLAGS = ; /tmp/unretained-helper\n"
                    b"\t/usr/bin/c++ $(CXX_FLAGS) -c source.cpp\n",
                    b"CXX_FLAGS = @/tmp/unretained-response\n"
                    b"\t/usr/bin/c++ $(CXX_FLAGS) -c source.cpp\n",
                    b"CMAKE_PROGRESS_1 = 1; /tmp/unretained-helper\n"
                    b"\t$(CMAKE_COMMAND) -E echo $(CMAKE_PROGRESS_1)\n",
                    b"CMAKE_SOURCE_DIR = /tmp/source; /tmp/helper\n"
                    b"\t$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR)\n"):
                with self.subTest(value_injection=value_injection):
                    bypass_build.write_bytes(value_injection)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "(?:value-bearing variable assignment|"
                            "external-only variable|shell code or indirection|"
                            "response-file indirection|progress value|"
                            "source directory)"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            bypass_build.write_bytes(b"# canonical fixture\n")

            for control in (
                    b"override SHELL := /tmp/unretained-shell\n",
                    b"export MAKEFILES := CMakeFiles/bypass.dir/extra.mk\n",
                    b"MAKEFLAGS += --eval=include\\ extra.mk\n",
                    b"target: MAKE := /tmp/unretained-make\n",
                    b"CMAKE_COMMAND = /usr/bin/ma${EMPTY}ke\n"
                    b"\t$(CMAKE_COMMAND) target\n",
                    b"MAKESILENT = -f CMakeFiles/bypass.dir/extra.mk\n"):
                with self.subTest(control=control):
                    bypass_build.write_bytes(control)
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "execution-control variable"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            for computed in (
                    b"I = include CMakeFiles/bypass.dir/extra.mk\n"
                    b"$(call eval,$(I))\n",
                    b"$(call shell,true)\n",
                    b"M := MA$()KE\n"
                    b"\t$($(M)) -f CMakeFiles/bypass.dir/extra.mk target\n",
                    b"M := /usr/bin/ma$()ke\n"
                    b"\t$(M) target\n",
                    b"M := MAKEFILE$()S\n"
                    b"$(M) = CMakeFiles/bypass.dir/extra.mk\n"
                    b"export $(M)\n",
                    b"define MAKESILEN$()T\n"
                    b"-f CMakeFiles/bypass.dir/extra.mk\n"
                    b"endef\n",
                    b"-load CMakeFiles/bypass.dir/plugin.so\n",
                    b"$(guile (display \"unexpected\"))\n",
                    b"\t$(value MAKE) "
                    b"-f CMakeFiles/bypass.dir/extra.mk target\n",
                    b"A = g\n"
                    b"M = $Amake\n"
                    b"all:\n"
                    b"\t$M -f CMakeFiles/bypass.dir/extra.mk target\n",
                    b"all:\n\t$$M target\n",
                    b"all:\n\t$@ target\n",
                    b"all:\n\t$ target\n",
                    b"all:\n\t$(MAKE target\n"):
                with self.subTest(computed=computed):
                    bypass_build.write_bytes(computed)
                    with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:dynamic (?:code or shell function|definition "
                        "or load directive)|(?:computed|non-canonical) "
                        "(?:variable reference|assignment|directive))"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            for malformed_lines in (
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 target\r\n",
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 tar\0get\n",
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 target"):
                makefile2.write_bytes(malformed_lines)
                with self.subTest(malformed_lines=malformed_lines):
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "non-canonical line framing"):
                        with provenance._TemporaryReplayRecipeTransport(
                                build, replacements, set(replacements)):
                            pass
            makefile2.write_bytes(original_makefile2)
            self.assertFalse(
                (build / "dependency-bypass-marker").exists())
            for derived in (
                    "CMakeFiles/target.dir/depend.make",
                    "CMakeFiles/target.dir/compiler_depend.make"):
                self.assertFalse(
                    provenance._is_replay_recipe_relative_path(derived))
            cache.unlink()
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "inventory omits a required CMake input"):
                provenance._TemporaryReplayRecipeTransport(
                    build, replacements, set(replacements))

    def test_replay_recipe_transport_uses_bounded_guard_instances(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-recipe-guard-bound-") as directory:
            build = Path(directory) / "build"
            cmake_files = build / "CMakeFiles"
            cmake_files.mkdir(parents=True)
            required = {
                "CMakeCache.txt": b"CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\n",
                "Makefile": (
                    b"\t/usr/bin/cmake -E true\n"
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 nested\n"),
                "CMakeFiles/Makefile2": b"\t/usr/bin/c++ -c source.cpp\n",
                "CMakeFiles/Makefile.cmake":
                    b"set(CMAKE_MAKE_PROGRAM \"/usr/bin/make\")\n",
                "CMakeFiles/cmake.check_cache": b"# generated\n",
            }
            for relative, content in required.items():
                path = build / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
            link = cmake_files / "target.dir/link.txt"
            link.parent.mkdir(parents=True)
            link.write_bytes(
                b"/usr/bin/ar qc libleopard.a object.o\n"
                b"/usr/bin/ranlib libleopard.a\n"
                b"/usr/bin/git status\n")
            for index in range(140):
                path = cmake_files / f"bulk-{index}.dir/build.make"
                path.parent.mkdir()
                path.write_bytes(b"# retained generated recipe\n")

            replacements = {
                "cmake": ("/usr/bin/cmake", "/proc/self/fd/101"),
                "cxx": ("/usr/bin/c++", "/proc/self/fd/102"),
                "ar": ("/usr/bin/ar", "/proc/self/fd/103"),
                "ranlib": ("/usr/bin/ranlib", "/proc/self/fd/104"),
                "git": ("/usr/bin/git", "/proc/self/fd/105"),
                "make": ("/usr/bin/make", "/proc/self/fd/106"),
            }
            active: set[int] = set()
            maximum_active = 0
            real_init = provenance._InotifyMutationGuard.__init__
            real_close = provenance._InotifyMutationGuard.close

            def counted_init(instance, label: str) -> None:
                nonlocal maximum_active
                real_init(instance, label)
                active.add(id(instance))
                maximum_active = max(maximum_active, len(active))

            def counted_close(instance) -> None:
                try:
                    real_close(instance)
                finally:
                    if instance.descriptor < 0:
                        active.discard(id(instance))

            with mock.patch.object(
                    provenance._InotifyMutationGuard, "__init__",
                    new=counted_init), mock.patch.object(
                        provenance._InotifyMutationGuard, "close",
                        new=counted_close):
                with provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements)) as transport:
                    self.assertEqual(
                        len(transport.manifest["retained_recipes"]),
                        146)
            self.assertFalse(active)
            self.assertLessEqual(maximum_active, 4)

            with mock.patch.object(
                    provenance, "MAX_REPLAY_RECIPE_FILES", 5):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "inventory exceeds its file bound"):
                    provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements))
            with mock.patch.object(
                    provenance, "MAX_REPLAY_RECIPE_TOTAL_BYTES", 1):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "aggregate byte bound"):
                    provenance._TemporaryReplayRecipeTransport(
                        build, replacements, set(replacements))
            original_total = sum(
                path.stat().st_size
                for path in provenance._replay_recipe_candidates(build))
            expanding = {
                "cmake": ("/usr/bin/cmake", "/proc/self/fd/101"),
            }
            with mock.patch.object(
                    provenance, "MAX_REPLAY_RECIPE_TOTAL_BYTES",
                    original_total):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "transported private replay recipes exceed their "
                        "aggregate byte bound"):
                    provenance._TemporaryReplayRecipeTransport(
                        build, expanding, set(expanding))

    def test_replay_recipe_transport_close_recovers_from_interruption(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-recipe-close-interrupt-") as directory:
            build = Path(directory) / "build"
            cmake_files = build / "CMakeFiles"
            cmake_files.mkdir(parents=True)
            inputs = {
                "CMakeCache.txt": b"CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\n",
                "Makefile": (
                    b"\t/usr/bin/cmake -E true\n"
                    b"\t$(MAKE) $(MAKESILENT) -f "
                    b"CMakeFiles/Makefile2 target\n"),
                "CMakeFiles/Makefile2":
                    b"\t/usr/bin/c++ -c source.cpp\n",
                "CMakeFiles/Makefile.cmake":
                    b"set(CMAKE_MAKE_PROGRAM \"/usr/bin/make\")\n",
                "CMakeFiles/cmake.check_cache": b"# generated\n",
                "CMakeFiles/target.dir/link.txt": (
                    b"/usr/bin/ar qc libleopard.a object.o\n"
                    b"/usr/bin/ranlib libleopard.a\n"
                    b"/usr/bin/git status\n"),
            }
            for relative, content in inputs.items():
                path = build / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
            replacements = {
                "cmake": ("/usr/bin/cmake", "/proc/self/fd/101"),
                "cxx": ("/usr/bin/c++", "/proc/self/fd/102"),
                "ar": ("/usr/bin/ar", "/proc/self/fd/103"),
                "ranlib": ("/usr/bin/ranlib", "/proc/self/fd/104"),
                "git": ("/usr/bin/git", "/proc/self/fd/105"),
                "make": ("/usr/bin/make", "/proc/self/fd/106"),
            }
            transport = provenance._TemporaryReplayRecipeTransport(
                build, replacements, set(replacements))
            rewritten = (build / "Makefile").read_bytes()
            self.assertNotEqual(rewritten, inputs["Makefile"])
            cleanup_code = transport._cleanup_owned_state.__code__
            armed = True

            def interrupt_before_cleanup(frame, event: str, argument):
                nonlocal armed
                del argument
                if (armed and event == "line" and
                        frame.f_code is cleanup_code):
                    armed = False
                    raise KeyboardInterrupt(
                        "injected replay cleanup interruption")
                return interrupt_before_cleanup

            previous_trace = sys.gettrace()
            try:
                sys.settrace(interrupt_before_cleanup)
                with self.assertRaisesRegex(
                        KeyboardInterrupt, "cleanup interruption"):
                    transport.close()
            finally:
                sys.settrace(previous_trace)
            self.assertFalse(transport.originals)
            self.assertFalse(transport._guards)
            for relative, content in inputs.items():
                self.assertEqual((build / relative).read_bytes(), content)

            transport = provenance._TemporaryReplayRecipeTransport(
                build, replacements, set(replacements))

            class InterruptingGuards(list):
                armed = True

                def __iter__(self):
                    if self.armed:
                        self.armed = False
                        raise KeyboardInterrupt(
                            "injected replay verification interruption")
                    return super().__iter__()

            transport._guards = InterruptingGuards(transport._guards)
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "verification interruption"):
                transport.close()
            self.assertFalse(transport.originals)
            self.assertFalse(transport._guards)
            for relative, content in inputs.items():
                self.assertEqual((build / relative).read_bytes(), content)

            transport = provenance._TemporaryReplayRecipeTransport(
                build, replacements, set(replacements))
            real_replace = provenance._replace_private_file
            interrupted = False

            def interrupt_one_restore(path, content: bytes, mode: int) -> None:
                nonlocal interrupted
                if not interrupted:
                    interrupted = True
                    raise KeyboardInterrupt(
                        "injected caught cleanup interruption")
                real_replace(path, content, mode)

            with mock.patch.object(
                    provenance, "_replace_private_file",
                    side_effect=interrupt_one_restore):
                with self.assertRaisesRegex(
                        KeyboardInterrupt,
                        "caught cleanup interruption"):
                    transport.close()
            self.assertTrue(interrupted)
            self.assertFalse(transport.originals)
            self.assertFalse(transport._guards)
            for relative, content in inputs.items():
                self.assertEqual((build / relative).read_bytes(), content)

    def test_end_to_end_replay_rejects_tool_parent_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-replay-aba-") as directory:
            root = Path(directory)
            active = root / "active"
            alternate = root / "alternate"
            saved = root / "saved"
            active.mkdir()
            alternate.mkdir()
            for target, marker in (
                    (active / "tool", "A"),
                    (alternate / "tool", "B")):
                target.write_text(
                    "#!/bin/sh\nprintf %s " + marker + "\n",
                    encoding="ascii")
                target.chmod(0o700)
            tool = active / "tool"
            identity = provenance.file_identity(tool, "replay tool")
            manifest = {
                "schema": "leopard2-tracked-source-tree/v1",
                "total_bytes": 0,
                "files": [],
                "git": {
                    "commit": "1" * 40,
                    "tree": "2" * 40,
                    "dirty": False,
                    "status_sha256": hashlib.sha256(b"").hexdigest(),
                },
                "git_tool": identity,
            }
            candidate = {
                "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
                "source_root": str(SOURCE_ROOT),
                "executable_target": "bench_leopard2",
                "tracked_source_manifest": manifest,
                "validated_cache": {
                    **{name: str(tool) for name in (
                        "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER",
                        "CMAKE_AR", "CMAKE_RANLIB",
                        "CMAKE_MAKE_PROGRAM", "CMAKE_LINKER",
                        "LEO2_BENCHMARK_GIT_EXECUTABLE")},
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                        provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                        "a" * 64,
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
                    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
                    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
                    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
                    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
                    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
                    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
                },
                "c_compiler": identity,
                "compiler": identity,
                "archiver": identity,
                "ranlib": identity,
                "make_program": identity,
                "linker": identity,
                "benchmark_git": identity,
                "compiler_subtools": [{
                    "language": "c++", "role": "cc1plus",
                    "identity": identity,
                }],
            }
            observed: list[bytes] = []

            class FakeSourceTree:
                def __init__(
                    self, unused_source: Path, destination: Path,
                    **unused_kwargs: object,
                ) -> None:
                    destination.mkdir()
                    self.manifest = manifest

                def __enter__(self):
                    return self

                def verify(self) -> None:
                    pass

                def __exit__(self, *unused_args: object) -> None:
                    pass

            calls = 0

            def swapped_run(*unused_args: object, **unused_kwargs: object):
                nonlocal calls
                calls += 1
                if calls == 1:
                    active.rename(saved)
                    alternate.rename(active)
                    completed = subprocess.run(
                        [str(tool)], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, check=True)
                    observed.append(completed.stdout)
                    active.rename(alternate)
                    saved.rename(active)
                return b""

            with mock.patch.object(
                    provenance, "_RetainedPrivateSourceTree",
                    FakeSourceTree), \
                    mock.patch.object(
                        provenance, "_reproducible_configure_argv",
                        return_value=[str(tool)]), \
                    mock.patch.object(
                        provenance, "_run", side_effect=swapped_run), \
                    mock.patch.object(
                        provenance, "candidate_build_provenance",
                        return_value=candidate), \
                    mock.patch.object(
                        provenance, "compare_reproducible_builds",
                        return_value={"proof": True}):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "pathname changed while retained"):
                    provenance.verify_reproducible_candidate_build(
                        candidate, jobs=1)
            self.assertEqual(observed, [b"B"])


class ExactCommandValidationTests(unittest.TestCase):
    def test_current_production_selector_cache_is_exact(self) -> None:
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "ENABLE_OPENMP": "ON",
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "a" * 64,
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            "LEOPARD_ENABLE_GF8": "ON",
            "LEOPARD_ENABLE_GF16": "ON",
        }
        validated = provenance._validate_candidate_required_cache(cache)
        v6_cache = dict(cache)
        for selector in (
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v6_cache.pop(selector)
        v6_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6
        v6_validated = provenance._validate_candidate_required_cache(
            v6_cache,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6)
        self.assertEqual(
            v6_validated["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6)
        v7_cache = dict(cache)
        for selector in (
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v7_cache.pop(selector)
        v7_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7
        v7_validated = provenance._validate_candidate_required_cache(
            v7_cache,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7)
        self.assertEqual(
            v7_validated["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7)
        v7_with_current_selector = dict(v7_cache)
        v7_with_current_selector[
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "ON"
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "historical candidate CMake cache"):
            provenance._validate_candidate_required_cache(
                v7_with_current_selector,
                expected_configuration_schema=
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7)
        v8_cache = dict(cache)
        for selector in (
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v8_cache.pop(selector)
        v8_cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8
        v8_validated = provenance._validate_candidate_required_cache(
            v8_cache,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8)
        self.assertEqual(
            v8_validated["LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"], "ON")
        v8_with_current_selector = dict(v8_cache)
        v8_with_current_selector[
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS"] = "ON"
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "historical candidate CMake cache"):
            provenance._validate_candidate_required_cache(
                v8_with_current_selector,
                expected_configuration_schema=
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8)
        selectors = {
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
        }
        for selector, expected in selectors.items():
            self.assertEqual(validated[selector], expected)
            for value in (
                    "OFF" if expected == "ON" else "ON",
                    "TRUE", "", None, False):
                with self.subTest(selector=selector, value=value):
                    changed = dict(cache)
                    if value is None:
                        del changed[selector]
                    else:
                        changed[selector] = value
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            selector.removeprefix("LEO2_")):
                        provenance._validate_candidate_required_cache(changed)
            with self.subTest(selector_type=selector), \
                    self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "expected one of \\['BOOL'\\]"):
                provenance.parse_cmake_cache(
                    f"{selector}:STRING={expected}\n".encode("ascii"))

    def test_v5_replay_capture_uses_its_frozen_cache_contract(self) -> None:
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "ENABLE_OPENMP": "ON",
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "a" * 64,
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            "LEOPARD_ENABLE_GF8": "ON",
            "LEOPARD_ENABLE_GF16": "ON",
        }
        validated = provenance._validate_candidate_required_cache(
            cache,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5)
        self.assertEqual(
            validated["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5)
        for selector in (
                "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL"):
            self.assertNotIn(selector, validated)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "EFFECTIVE_CONFIGURATION_SCHEMA"):
            provenance._validate_candidate_required_cache(cache)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "configuration schema is unsupported"):
            provenance._validate_candidate_required_cache(
                cache,
                expected_configuration_schema=
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4)

    def test_active_replay_generation_rejects_pre_v5_contracts(self) -> None:
        digest = "a" * 64
        for schema, selectors in (
                (
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
                    {"LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF"},
                ),
                (
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
                    {
                        "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                        "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT":
                            "ON",
                        "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                    },
                )):
            with self.subTest(schema=schema):
                candidate = {
                    "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
                    "source_root": "/does/not/need/to/exist",
                    "validated_cache": {
                        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                            schema,
                        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                            digest,
                        **selectors,
                    },
                }
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "supports only v5 through v8 and current"):
                    provenance.verify_reproducible_candidate_build(
                        candidate, jobs=1)

    def test_verify_replay_capture_forwards_v5_configuration_schema(
            self) -> None:
        candidate = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "a" * 64,
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            },
        }
        replay_contract = provenance._reproducible_replay_contract(
            candidate)
        rebuilt = {"captured": True}
        build = Path("/tmp/replayed-build")
        source = Path("/tmp/replayed-source")
        manifest = {"manifest": True}
        sealed = {}
        with mock.patch.object(
                provenance, "candidate_build_provenance",
                return_value=rebuilt) as capture:
            self.assertIs(
                provenance._capture_replayed_candidate_provenance(
                    build, source, "bench_leopard2",
                    replay_contract=replay_contract,
                    inherited_descriptors=(17,),
                    tracked_source_manifest=manifest,
                    logical_source_root=SOURCE_ROOT,
                    sealed_artifacts=sealed),
                rebuilt)
        capture.assert_called_once_with(
            build, source, build / "bench_leopard2", "bench_leopard2",
            inherited_descriptors=(17,),
            tracked_source_manifest=manifest,
            logical_source_root=SOURCE_ROOT,
            _expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
            sealed_artifacts=sealed)

    def test_replay_contract_cannot_downgrade_current_closures(self) -> None:
        digest = "a" * 64
        current = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
            },
        }
        v5 = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            },
        }
        v6 = copy.deepcopy(current)
        for selector in (
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v6["validated_cache"].pop(selector)
        v6["validated_cache"][
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6
        v7 = copy.deepcopy(current)
        for selector in (
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v7["validated_cache"].pop(selector)
        v7["validated_cache"][
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7
        v8 = copy.deepcopy(current)
        for selector in (
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v8["validated_cache"].pop(selector)
        v8["validated_cache"][
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8
        v4 = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
            },
        }
        v3 = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
            },
        }
        historical = {
            "schema": provenance.PRODUCTION_BUILD_CLOSURE_SCHEMA_V1,
            "validated_cache": {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
            },
        }
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                current, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "invocation_schema"],
            provenance.REPLAY_INVOCATION_SCHEMA)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v6, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v7, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v8, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v5, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v4, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                v3, provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA)[
                    "configuration_schema"],
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3)
        self.assertEqual(
            provenance._require_reproducible_replay_artifact_contract(
                historical,
                provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2,
                provenance.LEGACY_REPLAY_RECIPE_SCHEMA)[
                    "invocation_schema"],
            provenance.REPLAY_INVOCATION_SCHEMA_V1)
        for candidate, proof, recipe in (
            (
                current,
                provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2,
                provenance.LEGACY_REPLAY_RECIPE_SCHEMA,
            ),
            (
                historical,
                provenance.REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                provenance.CANONICAL_REPLAY_RECIPE_SCHEMA,
            ),
        ):
            with self.subTest(schema=candidate["schema"]):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "incompatible with its closure contract"):
                    provenance._require_reproducible_replay_artifact_contract(
                        candidate, proof, recipe)

        for selector in (
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE",
                "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            with self.subTest(current_missing=selector):
                missing_selector = copy.deepcopy(current)
                del missing_selector["validated_cache"][selector]
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "current reproducible-build closure configuration"):
                    provenance._reproducible_replay_contract(
                        missing_selector)

        v5_with_v6_selector = copy.deepcopy(v5)
        v5_with_v6_selector["validated_cache"][
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED"] = "ON"
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "v5 reproducible-build closure configuration"):
            provenance._reproducible_replay_contract(v5_with_v6_selector)

        v7_with_current_selector = copy.deepcopy(v7)
        v7_with_current_selector["validated_cache"][
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "ON"
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "v6/v7/v8/current reproducible-build closure"):
            provenance._reproducible_replay_contract(
                v7_with_current_selector)

        for selector, value in (
                ("LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS", "ON"),
                ("LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK", "OFF")):
            with self.subTest(v8_forward_selector=selector):
                v8_with_current_selector = copy.deepcopy(v8)
                v8_with_current_selector["validated_cache"][selector] = value
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "v6/v7/v8/current reproducible-build closure"):
                    provenance._reproducible_replay_contract(
                        v8_with_current_selector)

        for selector in (
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED"):
            with self.subTest(v4_extra=selector):
                v4_with_selector = copy.deepcopy(v4)
                v4_with_selector["validated_cache"][selector] = "OFF"
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "v4 reproducible-build closure configuration"):
                    provenance._reproducible_replay_contract(v4_with_selector)

        for selector in (
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"):
            with self.subTest(v3_extra=selector):
                v3_with_selector = copy.deepcopy(v3)
                v3_with_selector["validated_cache"][selector] = "ON"
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "v3 reproducible-build closure configuration"):
                    provenance._reproducible_replay_contract(
                        v3_with_selector)

        for selector in (
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"):
            with self.subTest(historical_extra=selector):
                historical_with_selector = copy.deepcopy(historical)
                historical_with_selector["validated_cache"][selector] = "OFF"
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "historical reproducible-build closure configuration"):
                    provenance._reproducible_replay_contract(
                        historical_with_selector)

    def test_benchmark_compile_digest_is_bound_to_retained_cache(self) -> None:
        source = (
            SOURCE_ROOT / "bench/leopard2/benchmark.cpp").resolve(strict=True)
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-config-digest-") as directory:
            build = Path(directory).resolve(strict=True)
            digest = "a" * 64
            header = (
                build / "generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            tokens = [
                str(tool_path("c++")),
                "-Wall", "-Wextra", "-fopenmp", "-fopenmp",
                "-O3", "-DNDEBUG", "-O3", "-std=gnu++11",
                f"-I{SOURCE_ROOT}",
                "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
                "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                '-DLEO2_BENCHMARK_BUILD_TYPE=\"Release\"',
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER=" +
                f'\"{header}\"',
                "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=" +
                f'\"{digest}\"',
                "-o", "benchmark.cpp.o", "-c", str(source),
            ]
            self.assertEqual(
                provenance._validate_compile_flags(
                    tokens, source, source_root=SOURCE_ROOT,
                    build_root=build,
                    cache={
                        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                            digest,
                    },
                    benchmark_source=True,
                    lexical_build_output=True),
                "portable-core")
            changed_cache = {
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "b" * 64,
            }
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "configuration digest is not bound"):
                provenance._validate_compile_flags(
                    tokens, source, source_root=SOURCE_ROOT,
                    build_root=build, cache=changed_cache,
                    benchmark_source=True,
                    lexical_build_output=True)

    def test_build_root_normalization_requires_path_boundaries(self) -> None:
        build = Path("/tmp/leopard2-build-A")
        self.assertEqual(
            provenance._normalize_build_token(
                f"-I{build}/generated", build),
            "-I${BUILD_ROOT}/generated")
        self.assertEqual(
            provenance._normalize_build_token(
                f"-DROOT={build}", build),
            "-DROOT=${BUILD_ROOT}")
        for value in (
            f"-DVALUE={build}-suffix",
            f"-DVALUE=prefix{build}/suffix",
            f"{build}suffix",
        ):
            with self.subTest(value=value):
                self.assertEqual(
                    provenance._normalize_build_token(value, build),
                    value)

    def test_compile_argv_is_direct_and_canonical(self) -> None:
        compiler = tool_path("c++")
        source = (SOURCE_ROOT / "leopard2.cpp").resolve(strict=True)
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-argv-") as directory:
            build = Path(directory).resolve(strict=True)
            output = "CMakeFiles/leopard.dir/leopard2.cpp.o"
            entry = {
                "directory": str(build),
                "file": str(source),
                "output": output,
            }
            arguments = [
                str(compiler), f"-I{build}/generated", "-O3",
                "-o", output, "-c", str(source),
            ]
            provenance._require_compile_driver(
                arguments, source, tool_path("cc"), compiler)
            normalized = provenance._canonical_compile_argv(
                entry, arguments, source, build)
            self.assertEqual(
                normalized[1], "-I${BUILD_ROOT}/generated")

            string_entry = {
                **entry, "command": shlex.join(arguments),
            }
            string_tokens = provenance._compile_tokens(string_entry)
            self.assertEqual(string_tokens, arguments)
            string_identity = provenance._compile_recipe_identity(
                string_entry, build)
            self.assertEqual(string_identity["representation"], "command")
            self.assertIn(
                "${BUILD_ROOT}/generated",
                string_identity["normalized_command"])
            raw_adversaries = (
                string_entry["command"].replace(
                    " -O3 ", "\n-O3 ", 1),
                string_entry["command"].replace(
                    " -O3 ", "\r-O3 ", 1),
                string_entry["command"].replace(
                    " -O3 ", "\0-O3 ", 1),
                string_entry["command"].replace(
                    " -O3 ", " -DVALUE=ok&&/usr/bin/true -O3 ", 1),
                string_entry["command"].replace(
                    " -O3 ", " '-DVALUE=ok;still' -O3 ", 1),
                string_entry["command"].replace(
                    " -O3 ", r" -DVALUE=\$HOME -O3 ", 1),
            )
            for raw in raw_adversaries:
                with self.subTest(raw=raw):
                    rejected = {**entry, "command": raw}
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "shell control"):
                        provenance._compile_tokens(rejected)

            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "non-canonical source operand"):
                provenance._canonical_compile_argv(
                    entry, arguments + ["-DUNUSED=1"], source, build)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "shell control"):
                provenance._canonical_compile_argv(
                    entry,
                    [str(compiler), "&&", "/usr/bin/true"] + arguments[1:],
                    source, build)
            for fused in (
                    "-DVALUE=ok&&/usr/bin/true",
                    "-DVALUE=$(touch /tmp/forbidden)",
                    "-DVALUE=`/usr/bin/true`"):
                with self.subTest(fused=fused):
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "shell control"):
                        provenance._canonical_compile_argv(
                            entry, arguments[:-2] + [fused] + arguments[-2:],
                            source, build)

    def test_compile_driver_rejects_same_binary_under_another_spelling(
            self) -> None:
        compiler = tool_path("c++")
        source = (SOURCE_ROOT / "leopard2.cpp").resolve(strict=True)
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-driver-") as directory:
            alias = Path(directory) / "configured-cxx"
            alias.symlink_to(compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "another command driver"):
                provenance._require_compile_driver(
                    [str(compiler), "-c", str(source)],
                    source, tool_path("cc"), alias)

    def test_target_flag_profiles_are_exact_and_order_independent(self) -> None:
        canonical = {
            "leopard2.cpp": ((), "portable-core"),
            "Leopard2BackendSSSE3.cpp": (
                ("-mssse3", "-mno-avx"), "ssse3-no-avx"),
            "Leopard2BackendAVX2.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2BackendAVX2T2K4.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2BackendAVX2T8K8B1024.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2BackendAVX2T16B64.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2BackendAVX2T32B256.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2LowP32B64AVX2.cpp": (
                ("-mavx2", "-mno-avx512f"), "avx2-no-avx512"),
            "Leopard2BackendGFNI.cpp": (
                ("-mavx2", "-mgfni", "-mno-avx512f"),
                "avx2-gfni-no-avx512"),
        }

        def compile_argv(
            source_name: str, target_flags: tuple[str, ...],
            extra: tuple[str, ...] = (),
        ) -> list[str]:
            enhanced = source_name.startswith((
                "Leopard2BackendSSSE3", "Leopard2BackendAVX2",
                "Leopard2BackendGFNI")) or \
                source_name == "Leopard2LowP32B64AVX2.cpp"
            return [
                str(tool_path("c++")),
                "-Wall", "-Wextra",
                *(("-fopenmp",) if enhanced else
                  ("-fopenmp", "-fopenmp")),
                "-O3", "-DNDEBUG", "-O3", "-std=gnu++11",
                *target_flags, *extra,
                "-o", f"{source_name}.o", "-c", source_name,
            ]

        for source_name, (target_flags, profile) in canonical.items():
            with self.subTest(source=source_name):
                tokens = compile_argv(
                    source_name, tuple(reversed(target_flags)))
                self.assertEqual(
                    provenance._validate_compile_flags(
                        tokens, Path(source_name)),
                    profile)

        ff8_source = (SOURCE_ROOT / "LeopardFF8.cpp").resolve(strict=True)
        backend_sources = {
            ff8_source,
            (SOURCE_ROOT / "Leopard2BackendSSSE3.cpp").resolve(strict=True),
            (SOURCE_ROOT / "Leopard2BackendAVX2.cpp").resolve(strict=True),
            (SOURCE_ROOT / "Leopard2BackendGFNI.cpp").resolve(strict=True),
        }
        ff8_tokens = [
            str(tool_path("c++")),
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_GFNI_BACKEND=1",
            "-DLEO2_HAVE_SSSE3_BACKEND=1",
            "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
            f"-I{SOURCE_ROOT}",
            "-Wall", "-Wextra", "-fopenmp",
            "-O3", "-DNDEBUG", "-O3", "-std=gnu++11", "-fopenmp",
            "-o", "LeopardFF8.cpp.o", "-c", str(ff8_source),
        ]
        self.assertEqual(
            provenance._validate_compile_flags(
                ff8_tokens, ff8_source, source_root=SOURCE_ROOT,
                cache={
                    "LEO2_BACKEND_VARIANT": "auto",
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                        provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
                },
                library_sources=backend_sources),
            "portable-core")
        experimental_tokens = list(ff8_tokens)
        experimental_tokens.insert(
            1, "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "non-canonical compile definitions"):
            provenance._validate_compile_flags(
                experimental_tokens, ff8_source,
                source_root=SOURCE_ROOT,
                cache={
                    "LEO2_BACKEND_VARIANT": "auto",
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                        provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
                    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
                    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
                },
                library_sources=backend_sources)

        adversaries = {
            "portable umbrella": (
                "leopard2.cpp", ("-march=x86-64-v4",)),
            "portable feature": (
                "leopard2.cpp", ("-msse4.2",)),
            "AVX2 umbrella": (
                "Leopard2BackendAVX2.cpp",
                ("-mavx2", "-mno-avx512f", "-march=skylake-avx512")),
            "AVX2 contradiction": (
                "Leopard2BackendAVX2.cpp",
                ("-mavx2", "-mno-avx512f", "-mno-avx2")),
            "AVX2 duplicate": (
                "Leopard2BackendAVX2.cpp",
                ("-mavx2", "-mavx2", "-mno-avx512f")),
            "GFNI contradiction": (
                "Leopard2BackendGFNI.cpp",
                ("-mavx2", "-mgfni", "-mno-avx512f", "-mno-gfni")),
            "SSSE3 contradiction": (
                "Leopard2BackendSSSE3.cpp",
                ("-mssse3", "-mno-avx", "-mno-ssse3")),
            "Clang feature escape": (
                "Leopard2BackendAVX2.cpp",
                ("-mavx2", "-mno-avx512f", "-Xclang",
                 "-target-feature", "-Xclang", "+avx512f")),
            "GCC specs escape": (
                "leopard2.cpp", ("-specs=/tmp/evil.specs",)),
            "GCC attached B escape": (
                "leopard2.cpp", ("-B/tmp/evil",)),
            "GCC split B escape": (
                "leopard2.cpp", ("-B", "/tmp/evil")),
            "compiler plugin escape": (
                "leopard2.cpp", ("-fplugin=/tmp/evil.so",)),
            "sysroot escape": (
                "leopard2.cpp", ("--sysroot=/tmp/evil",)),
            "driver wrapper escape": (
                "leopard2.cpp", ("-wrapper", "/tmp/evil")),
        }
        for name, (source_name, target_flags) in adversaries.items():
            with self.subTest(name=name):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "(?:non-canonical|indirect compile option)"):
                    provenance._validate_compile_flags(
                        compile_argv(source_name, (), target_flags),
                        Path(source_name))

        response = compile_argv(
            "leopard2.cpp", (), ("@/tmp/evil.rsp",))
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "response-file"):
            provenance._validate_compile_flags(
                response, Path("leopard2.cpp"))

    def test_k8_live_range_flag_requires_unambiguous_gnu_driver(self) -> None:
        source = Path("Leopard2BackendAVX2T8K8B1024.cpp")
        current_cache = {
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
        }

        def compile_argv(*extra: str) -> list[str]:
            return [
                str(tool_path("c++")),
                "-Wall", "-Wextra", "-fopenmp",
                "-O3", "-DNDEBUG", "-O3", "-std=gnu++11",
                "-mavx2", "-mno-avx512f", *extra,
                "-o", f"{source.name}.o", "-c", str(source),
            ]

        gnu_drivers = (
            "/usr/bin/g++",
            "/usr/bin/g++-13",
            "/usr/bin/x86_64-linux-gnu-g++-13",
            "/opt/cross/bin/aarch64-linux-gnu-g++",
        )
        for compiler_path in gnu_drivers:
            with self.subTest(gnu=compiler_path):
                self.assertTrue(
                    provenance._resolved_compiler_is_gnu(compiler_path))
                self.assertEqual(
                    provenance._validate_compile_flags(
                        compile_argv("-flive-range-shrinkage"), source,
                        cache=current_cache, compiler_path=compiler_path),
                    "avx2-no-avx512")
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-canonical or indirect compile option"):
                    provenance._validate_compile_flags(
                        compile_argv(), source, cache=current_cache,
                        compiler_path=compiler_path)

        ambiguous_drivers = (
            "/usr/bin/compiler",
            "/usr/bin/clang++-18",
            "/usr/bin/ccache",
            "/usr/lib/ccache-g++",
            "/usr/bin/vendor-g++",
            "g++",
        )
        for compiler_path in ambiguous_drivers:
            with self.subTest(ambiguous=compiler_path):
                self.assertFalse(
                    provenance._resolved_compiler_is_gnu(compiler_path))
                self.assertEqual(
                    provenance._validate_compile_flags(
                        compile_argv(), source, cache=current_cache,
                        compiler_path=compiler_path),
                    "avx2-no-avx512")
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-canonical or indirect compile option"):
                    provenance._validate_compile_flags(
                        compile_argv("-flive-range-shrinkage"), source,
                        cache=current_cache, compiler_path=compiler_path)

    def test_isolated_avx2_member_definitions_are_exact(self) -> None:
        cache = {
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
            "LEO2_FLAG_FALIGN_FUNCTIONS_64": "1",
            "LEOPARD_ENABLE_GF16": "ON",
        }
        common_definitions = {
            "-DNDEBUG",
            "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
            "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
        }
        member_definitions = {
            "Leopard2BackendAVX2T2K4.cpp": {
                "-DLEO2_HAVE_AVX2_BACKEND=1",
            },
            "Leopard2BackendAVX2T8K8B1024.cpp": {
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
            },
            "Leopard2BackendAVX2T16B64.cpp": {
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            },
            "Leopard2BackendAVX2T32B256.cpp": {
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
                "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=0",
            },
            "Leopard2LowP32B64AVX2.cpp": {
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            },
        }
        member_sources = {
            name: (SOURCE_ROOT / name).resolve(strict=True)
            for name in member_definitions
        }
        library_sources = set(member_sources.values())

        def compile_argv(
            source: Path, definitions: set[str], *,
            target_flags: tuple[str, ...] = (
                "-mavx2", "-mno-avx512f", "-falign-functions=64"),
        ) -> list[str]:
            enhanced = source.name.startswith((
                "Leopard2BackendSSSE3", "Leopard2BackendAVX2",
                "Leopard2BackendGFNI")) or \
                source.name == "Leopard2LowP32B64AVX2.cpp"
            return [
                str(tool_path("c++")),
                *sorted(definitions),
                f"-I{SOURCE_ROOT}",
                "-Wall", "-Wextra",
                *(("-fopenmp",) if enhanced else
                  ("-fopenmp", "-fopenmp")),
                "-O3", "-O3", "-std=gnu++11",
                *target_flags,
                "-o", f"{source.name}.o", "-c", str(source),
            ]

        canonical_tokens: dict[str, list[str]] = {}
        for name, specific_definitions in member_definitions.items():
            source = member_sources[name]
            definitions = common_definitions | specific_definitions
            tokens = compile_argv(source, definitions)
            canonical_tokens[name] = tokens
            with self.subTest(member=name):
                self.assertEqual(
                    provenance._validate_compile_flags(
                        tokens, source, source_root=SOURCE_ROOT,
                        cache=cache, library_sources=library_sources),
                    "avx2-no-avx512")
            for definition in specific_definitions:
                missing = [token for token in tokens if token != definition]
                with self.subTest(member=name, missing=definition), \
                        self.assertRaisesRegex(
                            provenance.BuildProvenanceError,
                            "non-canonical compile definitions"):
                    provenance._validate_compile_flags(
                        missing, source, source_root=SOURCE_ROOT,
                        cache=cache, library_sources=library_sources)
            mature_selector = list(tokens)
            mature_selector.insert(
                1, "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1")
            with self.subTest(member=name, mature_selector=True), \
                    self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-canonical compile definitions"):
                provenance._validate_compile_flags(
                    mature_selector, source, source_root=SOURCE_ROOT,
                    cache=cache, library_sources=library_sources)
            missing_alignment = [
                token for token in tokens
                if token != "-falign-functions=64"
            ]
            with self.subTest(member=name, missing_alignment=True), \
                    self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-canonical or indirect compile option"):
                provenance._validate_compile_flags(
                    missing_alignment, source, source_root=SOURCE_ROOT,
                    cache=cache, library_sources=library_sources)

        # The retained replay cache predates the T32 two-block selector.  Its
        # archive member proves the default-on path while generated T32 is
        # retained as explicitly off.
        replay_cache = {
            key: value for key, value in cache.items()
            if key not in {
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
            }
        }
        t32_name = "Leopard2BackendAVX2T32B256.cpp"
        self.assertEqual(
            provenance._validate_compile_flags(
                canonical_tokens[t32_name], member_sources[t32_name],
                source_root=SOURCE_ROOT, cache=replay_cache,
                library_sources=library_sources),
            "avx2-no-avx512")
        p32_name = "Leopard2LowP32B64AVX2.cpp"
        self.assertEqual(
            provenance._validate_compile_flags(
                canonical_tokens[p32_name], member_sources[p32_name],
                source_root=SOURCE_ROOT, cache=replay_cache,
                library_sources=library_sources),
            "avx2-no-avx512")

        generated_cache = dict(cache)
        generated_cache.update({
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "OFF",
        })
        generated_definitions = common_definitions | {
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1",
            "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=1",
        }
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(
                    member_sources[t32_name], generated_definitions),
                member_sources[t32_name], source_root=SOURCE_ROOT,
                cache=generated_cache, library_sources=library_sources),
            "avx2-no-avx512")

        gf16_off_cache = dict(cache)
        gf16_off_cache["LEOPARD_ENABLE_GF16"] = "OFF"
        gf16_off_definitions = common_definitions | {
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DNO_LEO_HAS_FF16=1",
        }
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(
                    member_sources[p32_name], gf16_off_definitions),
                member_sources[p32_name], source_root=SOURCE_ROOT,
                cache=gf16_off_cache, library_sources=library_sources),
            "avx2-no-avx512")
        unexpected_field_disable = compile_argv(
            member_sources[p32_name], gf16_off_definitions)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "non-canonical compile definitions"):
            provenance._validate_compile_flags(
                unexpected_field_disable, member_sources[p32_name],
                source_root=SOURCE_ROOT, cache=cache,
                library_sources=library_sources)

        disabled_object_caches = {
            "Leopard2BackendAVX2T16B64.cpp": {
                "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "OFF",
            },
            "Leopard2BackendAVX2T32B256.cpp": {
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
                "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "OFF",
            },
            "Leopard2LowP32B64AVX2.cpp": {
                "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "OFF",
            },
        }
        for name, overrides in disabled_object_caches.items():
            disabled_cache = dict(cache)
            disabled_cache.update(overrides)
            with self.subTest(disabled_object=name), \
                    self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "contradicts its disabled selector"):
                provenance._validate_compile_flags(
                    canonical_tokens[name], member_sources[name],
                    source_root=SOURCE_ROOT, cache=disabled_cache,
                    library_sources=library_sources)

        lookalike = SOURCE_ROOT / "Leopard2BackendAVX2Lookalike.cpp"
        lookalike_tokens = compile_argv(
            lookalike, common_definitions | {
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1"})
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "not an allowlisted production AVX2"):
            provenance._validate_compile_flags(
                lookalike_tokens, lookalike, source_root=SOURCE_ROOT,
                cache=cache, library_sources=library_sources | {lookalike})

        mature_definitions = {
            "Leopard2BackendAVX2.cpp": common_definitions | {
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1",
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            },
            "Leopard2BackendAVX2Xor.cpp": common_definitions | {
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            },
        }
        for name, definitions in mature_definitions.items():
            source = (SOURCE_ROOT / name).resolve(strict=True)
            with self.subTest(mature=name):
                self.assertEqual(
                    provenance._validate_compile_flags(
                        compile_argv(source, definitions), source,
                        source_root=SOURCE_ROOT, cache=cache,
                        library_sources=library_sources | {source}),
                    "avx2-no-avx512")
            low_p32_selector = compile_argv(
                source, definitions | {
                    "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1"})
            with self.subTest(mature=name, low_p32_selector=True), \
                    self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "non-canonical compile definitions"):
                provenance._validate_compile_flags(
                    low_p32_selector, source, source_root=SOURCE_ROOT,
                    cache=cache,
                    library_sources=library_sources | {source})

        router = (SOURCE_ROOT / "leopard2.cpp").resolve(strict=True)
        mature = (SOURCE_ROOT / "Leopard2BackendAVX2.cpp").resolve(
            strict=True)
        routed_library_sources = library_sources | {router, mature}
        router_definitions = common_definitions | {
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1",
            "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
            "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
        }
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(router, router_definitions, target_flags=()),
                router, source_root=SOURCE_ROOT, cache=replay_cache,
                library_sources=routed_library_sources),
            "portable-core")
        dual_off_cache = dict(replay_cache)
        dual_off_cache["LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "OFF"
        dual_off_definitions = (
            router_definitions - {
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1"}) | {
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=0"}
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(
                    router, dual_off_definitions, target_flags=()),
                router, source_root=SOURCE_ROOT, cache=dual_off_cache,
                library_sources=routed_library_sources),
            "portable-core")
        noncanonical_dual_cache = dict(replay_cache)
        noncanonical_dual_cache[
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "TRUE"
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "exactly ON or OFF"):
            provenance._validate_compile_flags(
                compile_argv(router, router_definitions, target_flags=()),
                router, source_root=SOURCE_ROOT,
                cache=noncanonical_dual_cache,
                library_sources=routed_library_sources)
        for selector, definition, changed_value, changed_definition in (
                (
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                    "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
                    "OFF",
                    "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=0",
                ),
                (
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK",
                    "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0",
                    "ON",
                    "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=1",
                )):
            with self.subTest(router_selector=selector):
                changed_cache = dict(replay_cache)
                changed_cache[selector] = changed_value
                changed_definitions = (
                    router_definitions - {definition}) | {
                        changed_definition}
                self.assertEqual(
                    provenance._validate_compile_flags(
                        compile_argv(
                            router, changed_definitions, target_flags=()),
                        router, source_root=SOURCE_ROOT,
                        cache=changed_cache,
                        library_sources=routed_library_sources),
                    "portable-core")
                noncanonical_cache = dict(replay_cache)
                noncanonical_cache[selector] = "TRUE"
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "exactly ON or OFF"):
                    provenance._validate_compile_flags(
                        compile_argv(
                            router, router_definitions, target_flags=()),
                        router, source_root=SOURCE_ROOT,
                        cache=noncanonical_cache,
                        library_sources=routed_library_sources)
        v8_router_cache = dict(replay_cache)
        for selector in (
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v8_router_cache.pop(selector)
        v8_router_cache[
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(
                    router,
                    router_definitions - {
                        "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
                        "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0",
                    },
                    target_flags=()),
                router, source_root=SOURCE_ROOT, cache=v8_router_cache,
                library_sources=routed_library_sources),
            "portable-core")
        v7_router_cache = dict(replay_cache)
        for selector in (
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            v7_router_cache.pop(selector)
        v7_router_cache[
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = \
            provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(
                    router,
                    router_definitions - {
                        "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1",
                        "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
                        "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0",
                    },
                    target_flags=()),
                router, source_root=SOURCE_ROOT, cache=v7_router_cache,
                library_sources=routed_library_sources),
            "portable-core")
        missing_low_p32 = compile_argv(
            router,
            router_definitions - {
                "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1"},
            target_flags=())
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "non-canonical compile definitions"):
            provenance._validate_compile_flags(
                missing_low_p32, router, source_root=SOURCE_ROOT,
                cache=replay_cache,
                library_sources=routed_library_sources)

        ff8 = (SOURCE_ROOT / "LeopardFF8.cpp").resolve(strict=True)
        ff8_library_sources = routed_library_sources | {ff8}
        ff8_definitions = common_definitions | {
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
        }
        self.assertEqual(
            provenance._validate_compile_flags(
                compile_argv(ff8, ff8_definitions, target_flags=()),
                ff8, source_root=SOURCE_ROOT, cache=replay_cache,
                library_sources=ff8_library_sources),
            "portable-core")

    def test_low_p32_source_discovery_tracks_the_selector(self) -> None:
        p32_name = "Leopard2LowP32B64AVX2.cpp"
        t16_name = "Leopard2BackendAVX2T16B64.cpp"
        t32_name = "Leopard2BackendAVX2T32B256.cpp"
        tracked_names = set(provenance.CORE_LIBRARY_SOURCES) | {
            "LeopardFF8.cpp",
            "Leopard2BackendSSSE3.cpp",
            "Leopard2BackendAVX2.cpp",
            "Leopard2BackendAVX2Xor.cpp",
            "Leopard2BackendAVX2T2K4.cpp",
            "Leopard2BackendAVX2T8K8B1024.cpp",
            t16_name,
            t32_name,
            p32_name,
        }
        tracked = {
            (SOURCE_ROOT / name).resolve(strict=True)
            for name in tracked_names
        }
        lookalike = (
            SOURCE_ROOT / "Leopard2BackendAVX2Lookalike.cpp").resolve(
                strict=False)
        tracked.add(lookalike)
        cache = {
            "LEOPARD_ENABLE_GF16": "OFF",
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
        }
        expected = provenance._expected_library_sources(
            SOURCE_ROOT, cache, tracked)
        t2_path = (
            SOURCE_ROOT / "Leopard2BackendAVX2T2K4.cpp").resolve(
                strict=True)
        t8_k8_b1024_path = (
            SOURCE_ROOT / "Leopard2BackendAVX2T8K8B1024.cpp").resolve(
                strict=True)
        self.assertIn(t2_path, expected)
        self.assertIn(t8_k8_b1024_path, expected)
        self.assertIn((SOURCE_ROOT / p32_name).resolve(strict=True), expected)
        self.assertIn((SOURCE_ROOT / t16_name).resolve(strict=True), expected)
        self.assertIn((SOURCE_ROOT / t32_name).resolve(strict=True), expected)
        self.assertNotIn(lookalike, expected)

        historical_expected = provenance._expected_library_sources(
            SOURCE_ROOT, cache, tracked,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5)
        self.assertNotIn(t2_path, historical_expected)
        self.assertNotIn(t8_k8_b1024_path, historical_expected)
        self.assertIn(
            (SOURCE_ROOT / t16_name).resolve(strict=True),
            historical_expected)
        self.assertIn(
            (SOURCE_ROOT / t32_name).resolve(strict=True),
            historical_expected)

        v6_expected = provenance._expected_library_sources(
            SOURCE_ROOT, cache, tracked,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6)
        self.assertIn(t2_path, v6_expected)
        self.assertNotIn(t8_k8_b1024_path, v6_expected)

        v7_expected = provenance._expected_library_sources(
            SOURCE_ROOT, cache, tracked,
            expected_configuration_schema=
                provenance.BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7)
        self.assertIn(t2_path, v7_expected)
        self.assertIn(t8_k8_b1024_path, v7_expected)

        disabled_cache = dict(cache)
        disabled_cache["LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED"] = "OFF"
        disabled_cache["LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK"] = "OFF"
        disabled_cache["LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL"] = "OFF"
        expected_without_p32 = provenance._expected_library_sources(
            SOURCE_ROOT, disabled_cache, tracked)
        self.assertNotIn(
            (SOURCE_ROOT / p32_name).resolve(strict=True),
            expected_without_p32)
        self.assertNotIn(
            (SOURCE_ROOT / t16_name).resolve(strict=True),
            expected_without_p32)
        self.assertNotIn(
            (SOURCE_ROOT / t32_name).resolve(strict=True),
            expected_without_p32)

    def test_archive_and_executable_recipe_shapes_are_exact(self) -> None:
        archiver = tool_path("ar")
        ranlib = tool_path("ranlib")
        compiler = tool_path("c++")
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-recipe-") as directory:
            build = Path(directory).resolve(strict=True)
            archive = (build / "libleopard.a").resolve(strict=False)
            executable = (build / "bench_leopard2").resolve(strict=False)
            object_name = "CMakeFiles/leopard.dir/leopard2.cpp.o"
            benchmark_object = (
                "CMakeFiles/bench_leopard2.dir/"
                "bench/leopard2/benchmark.cpp.o")
            libgomp = compiler_runtime_path("libgomp.so")
            libpthread = compiler_runtime_path("libpthread.a")
            archive_recipe = (
                f"{archiver} qc libleopard.a {object_name}\n"
                f"{ranlib} libleopard.a\n").encode()
            archive_objects, archive_commands = \
                provenance._archive_recipe_semantics(
                    archive_recipe, build, archive, archiver, ranlib)
            self.assertEqual(
                archive_objects,
                [(build / object_name).resolve(strict=False)])
            self.assertEqual(len(archive_commands), 2)
            provenance._require_exact_compile_object_operand(
                {"output": object_name}, archive_commands[0][3],
                build, "archive object")
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "does not exactly name its compile output"):
                provenance._require_exact_compile_object_operand(
                    {"output": object_name},
                    "./" + archive_commands[0][3],
                    build, "archive object")

            executable_recipe = (
                f"{compiler} -Wall -Wextra -fopenmp -O3 -DNDEBUG -O3 "
                f"{benchmark_object} -o bench_leopard2 libleopard.a "
                f"{libgomp} {libpthread}\n").encode()
            executable_objects, executable_command = \
                provenance._executable_recipe_semantics(
                    executable_recipe, build, executable, archive, compiler)
            self.assertEqual(
                executable_objects,
                [(build / benchmark_object).resolve(strict=False)])
            self.assertEqual(executable_command[-2:], [
                str(libgomp), str(libpthread)])
            allk_executable = (
                build / "bench_leopard2_allk").resolve(strict=False)
            allk_object = (
                "CMakeFiles/bench_leopard2_allk.dir/"
                "bench/leopard2/benchmark.cpp.o")
            allk_recipe = (
                f"{compiler} -Wall -Wextra -fopenmp -O3 -DNDEBUG -O3 "
                f"{allk_object} -o bench_leopard2_allk libleopard.a "
                f"{libgomp} {libpthread}\n").encode()
            allk_objects, _ = provenance._executable_recipe_semantics(
                allk_recipe, build, allk_executable, archive, compiler)
            self.assertEqual(
                allk_objects,
                [(build / allk_object).resolve(strict=False)])

            for invalid in (
                    b"/usr/bin/true\n" + archive_recipe,
                    archive_recipe.replace(
                        b" qc libleopard.a", b" qcD libleopard.a"),
                    archive_recipe.replace(
                        str(archiver).encode(),
                        b"/usr/bin/env " + str(archiver).encode(), 1),
                    archive_recipe.replace(
                        object_name.encode(),
                        object_name.encode() + b" --plugin=/tmp/plugin"),
                    archive_recipe.replace(
                        object_name.encode(), b"--plugin=plugin.o"),
                    archive_recipe.replace(
                        object_name.encode(), b"--plugin=/tmp/plugin.o"),
                    archive_recipe.replace(
                        object_name.encode(), b"-plugin=plugin.o")):
                with self.subTest(invalid=invalid):
                    with self.assertRaises(
                            provenance.BuildProvenanceError):
                        provenance._archive_recipe_semantics(
                            invalid, build, archive, archiver, ranlib)

            for invalid in (
                    b"/usr/bin/true\n" + executable_recipe,
                    executable_recipe.replace(
                        str(compiler).encode(),
                        b"/usr/bin/env " + str(compiler).encode(), 1),
                    executable_recipe.replace(
                        b" -Wall ", b" -Wall&&/usr/bin/true ", 1),
                    executable_recipe.replace(
                        b" -Wall ", b" -Wl,-rpath,$ORIGIN -Wall ", 1),
                    executable_recipe.replace(
                        b" libleopard.a ", b" libleopard.a -levil ", 1),
                    executable_recipe.replace(
                        b" libleopard.a ", b" libleopard.a -Wl,--evil ", 1),
                    executable_recipe.replace(
                        b" libleopard.a ",
                        b" libleopard.a -Wl,-rpath,/tmp ", 1),
                    executable_recipe.replace(
                        b" libleopard.a ", b" libleopard.a evil.a ", 1),
                    executable_recipe.replace(
                        b" libleopard.a ", b" libleopard.a /tmp/evil.a ", 1)):
                with self.subTest(invalid=invalid):
                    with self.assertRaises(
                            provenance.BuildProvenanceError):
                        provenance._executable_recipe_semantics(
                            invalid, build, executable, archive, compiler)


class PidfdWrapperTests(unittest.TestCase):
    def test_mapping_descriptor_cleanup_reconciles_interruptions(self) -> None:
        descriptor = os.open("/dev/null", os.O_RDONLY)
        mapping = {"descriptor": descriptor}
        replacement = -1
        real_close = os.close

        def close_then_interrupt(value: int) -> None:
            nonlocal replacement
            self.assertEqual(value, descriptor)
            real_close(value)
            replacement = os.open("/dev/zero", os.O_RDONLY)
            self.assertEqual(replacement, descriptor)
            raise KeyboardInterrupt("injected mapped-close interruption")

        with mock.patch.object(
                provenance.os, "close", side_effect=close_then_interrupt):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "mapped-close interruption"):
                provenance._close_retained_mapping_descriptor(
                    mapping, "descriptor")
        try:
            self.assertFalse(mapping)
            os.fstat(replacement)
        finally:
            real_close(replacement)

        descriptor = os.open("/dev/null", os.O_RDONLY)
        mapping = {"descriptor": descriptor}
        with mock.patch.object(
                provenance.os, "close",
                side_effect=KeyboardInterrupt(
                    "injected pre-close interruption")):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "pre-close interruption"):
                provenance._close_retained_mapping_descriptor(
                    mapping, "descriptor")
        self.assertFalse(mapping)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_mapping_descriptor_same_inode_recycle_is_not_closed(
            self) -> None:
        for failure_type in (OSError, KeyboardInterrupt, SystemExit):
            with self.subTest(failure=failure_type.__name__), \
                    tempfile.TemporaryDirectory(
                        prefix="leo2-mapping-same-inode-aba-") as directory:
                path = Path(directory) / "same-inode"
                path.write_bytes(b"mapping descriptor")
                descriptor = os.open(path, os.O_RDONLY)
                mapping = {"descriptor": descriptor}
                replacement = -1
                guards: list[int] = []
                real_close = os.close
                real_dup = os.dup
                failure = failure_type(
                    f"injected mapping {failure_type.__name__}")

                def tracked_dup(value: int) -> int:
                    self.assertEqual(value, descriptor)
                    guard = real_dup(value)
                    guards.append(guard)
                    return guard

                def close_reopen_same_inode(value: int) -> None:
                    nonlocal replacement
                    self.assertEqual(value, descriptor)
                    real_close(value)
                    replacement = os.open(path, os.O_RDONLY)
                    self.assertEqual(replacement, descriptor)
                    raise failure

                try:
                    with mock.patch.object(
                            provenance.os, "dup",
                            side_effect=tracked_dup), mock.patch.object(
                                provenance.os, "close",
                                side_effect=close_reopen_same_inode):
                        with self.assertRaises(failure_type) as raised:
                            provenance._close_retained_mapping_descriptor(
                                mapping, "descriptor")
                    self.assertIs(raised.exception, failure)
                    self.assertFalse(mapping)
                    os.fstat(replacement)
                    self.assertEqual(len(guards), 1)
                    with self.assertRaises(OSError):
                        os.fstat(guards[0])
                finally:
                    if replacement >= 0:
                        real_close(replacement)

    def test_proc_record_binds_exact_task_directory_inode(self) -> None:
        record = provenance._proc_process_record(os.getpid())
        self.assertIsNotNone(record)
        status = os.stat(
            f"/proc/{os.getpid()}", follow_symlinks=False)
        self.assertEqual(record[5:], (status.st_dev, status.st_ino))

    def test_python_pidfd_wrappers_are_preferred(self) -> None:
        with mock.patch.object(
                provenance.os, "pidfd_open", create=True,
                return_value=73) as pidfd_open, \
                mock.patch.object(
                    provenance.ctypes, "CDLL",
                    side_effect=AssertionError("libc fallback used")):
            self.assertEqual(provenance._linux_pidfd_open(123), 73)
        pidfd_open.assert_called_once_with(123, 0)

        with mock.patch.object(
                provenance.signal, "pidfd_send_signal", create=True
                ) as pidfd_signal, \
                mock.patch.object(
                    provenance.ctypes, "CDLL",
                    side_effect=AssertionError("libc fallback used")):
            provenance._linux_pidfd_signal(73, provenance.signal.SIGKILL)
        pidfd_signal.assert_called_once_with(
            73, provenance.signal.SIGKILL, None, 0)

    def test_pidfd_libc_fallback_is_retained(self) -> None:
        library = mock.Mock()
        library.pidfd_open.return_value = 81
        library.pidfd_send_signal.return_value = 0
        with mock.patch.object(
                provenance.os, "pidfd_open", None, create=True), \
                mock.patch.object(
                    provenance.signal, "pidfd_send_signal", None,
                    create=True), \
                mock.patch.object(
                    provenance.ctypes, "CDLL", return_value=library):
            self.assertEqual(provenance._linux_pidfd_open(456), 81)
            provenance._linux_pidfd_signal(
                81, provenance.signal.SIGKILL)
        library.pidfd_open.assert_called_once()
        library.pidfd_send_signal.assert_called_once()

    def test_pidfd_wrapper_and_fallback_errors_fail_closed(self) -> None:
        vanished = ProcessLookupError(
            errno.ESRCH, os.strerror(errno.ESRCH))
        denied = PermissionError(
            errno.EPERM, os.strerror(errno.EPERM))
        with mock.patch.object(
                provenance.os, "pidfd_open", create=True,
                side_effect=vanished):
            self.assertIsNone(provenance._linux_pidfd_open(789))
        with mock.patch.object(
                provenance.os, "pidfd_open", create=True,
                side_effect=denied):
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "cannot open pidfd"):
                provenance._linux_pidfd_open(789)

        library = mock.Mock()
        library.pidfd_send_signal.return_value = -1
        with mock.patch.object(
                provenance.signal, "pidfd_send_signal", None,
                create=True), \
                mock.patch.object(
                    provenance.ctypes, "CDLL", return_value=library), \
                mock.patch.object(
                    provenance.ctypes, "get_errno",
                    return_value=errno.EINVAL):
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "cannot signal"):
                provenance._linux_pidfd_signal(
                    81, provenance.signal.SIGKILL)
        with mock.patch.object(
                provenance.signal, "pidfd_send_signal", create=True,
                side_effect=vanished):
            provenance._linux_pidfd_signal(
                81, provenance.signal.SIGKILL)
        with mock.patch.object(
                provenance.signal, "pidfd_send_signal", create=True,
                side_effect=denied):
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "cannot signal"):
                provenance._linux_pidfd_signal(
                    81, provenance.signal.SIGKILL)

        library = mock.Mock()
        library.pidfd_open.return_value = -1
        with mock.patch.object(
                provenance.os, "pidfd_open", None, create=True), \
                mock.patch.object(
                    provenance.ctypes, "CDLL", return_value=library), \
                mock.patch.object(
                    provenance.ctypes, "get_errno",
                    return_value=errno.ENOSYS):
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "cannot open pidfd"):
                provenance._linux_pidfd_open(789)


class OwnerExceptionPrecedenceTests(unittest.TestCase):
    """Every retained-owner family obeys one terminal-exception policy."""

    _SIMPLE_OWNER_EXITS = (
        ("descriptor", provenance._OwnedDescriptor.__exit__),
        ("inotify", provenance._InotifyMutationGuard.__exit__),
        ("artifact", provenance._ReplayArtifactSink.__exit__),
        ("file", provenance._RetainedFileSnapshot.__exit__),
        ("directory", provenance._RetainedDirectoryTree.__exit__),
        ("git", provenance._RetainedGitMetadata.__exit__),
        ("source", provenance._RetainedPrivateSourceTree.__exit__),
        ("proc", provenance._ProcProcessSnapshot.__exit__),
        ("recipe", provenance._TemporaryReplayRecipeTransport.__exit__),
    )

    @staticmethod
    def _simple_owner(cleanup: BaseException) -> object:
        class Owner:
            label = "injected owner"

            def close(self) -> None:
                raise cleanup

        return Owner()

    @staticmethod
    def _descriptor_cleanup_failure(cleanup: BaseException):
        def fail(
                descriptor: int, context: str,
                consumed: list[bool]) -> None:
            del descriptor, context
            # Model the current guarded-close contract: ownership has been
            # consumed even when the close operation reports a failure.
            consumed[0] = True
            raise cleanup

        return fail

    @classmethod
    def _patch_descriptor_cleanup_failure(cls, cleanup: BaseException):
        return mock.patch.object(
            provenance, "_close_descriptor_with_ofd_guard",
            side_effect=cls._descriptor_cleanup_failure(cleanup))

    def test_owner_precedence_is_idempotent_for_same_exception(self) -> None:
        failure = KeyboardInterrupt("same owner failure")
        failure.add_note("retained cleanup detail")
        notes_before = tuple(failure.__notes__)
        selected = provenance._owner_exception_precedence(
            failure, failure, "same owner exception")
        self.assertIs(selected, failure)
        self.assertEqual(tuple(failure.__notes__), notes_before)

    def test_every_simple_owner_retains_earlier_terminal(self) -> None:
        for name, owner_exit in self._SIMPLE_OWNER_EXITS:
            with self.subTest(owner=name, primary="KeyboardInterrupt"):
                primary = KeyboardInterrupt(f"{name} primary interrupt")
                cleanup = OSError(f"{name} ordinary cleanup")
                owner = self._simple_owner(cleanup)
                with self.assertRaises(KeyboardInterrupt) as raised:
                    owner_exit(
                        owner, type(primary), primary,
                        primary.__traceback__)
                self.assertIs(raised.exception, primary)
                self.assertTrue(any(
                    "later cleanup failure" in note and
                    "ordinary cleanup" in note
                    for note in getattr(primary, "__notes__", ())))

            with self.subTest(owner=name, primary="SystemExit"):
                primary = SystemExit(f"{name} primary exit")
                cleanup = KeyboardInterrupt(f"{name} later interrupt")
                owner = self._simple_owner(cleanup)
                with self.assertRaises(SystemExit) as raised:
                    owner_exit(
                        owner, type(primary), primary,
                        primary.__traceback__)
                self.assertIs(raised.exception, primary)
                self.assertTrue(any(
                    "later cleanup failure" in note and
                    "later interrupt" in note
                    for note in getattr(primary, "__notes__", ())))

    def test_every_simple_owner_promotes_later_terminal(self) -> None:
        for terminal_type in (KeyboardInterrupt, SystemExit):
            for name, owner_exit in self._SIMPLE_OWNER_EXITS:
                with self.subTest(
                        owner=name, cleanup=terminal_type.__name__):
                    primary = ValueError(f"{name} ordinary primary")
                    cleanup = terminal_type(
                        f"{name} terminal cleanup")
                    owner = self._simple_owner(cleanup)
                    with self.assertRaises(terminal_type) as raised:
                        owner_exit(
                            owner, type(primary), primary,
                            primary.__traceback__)
                    self.assertIs(raised.exception, cleanup)
                    self.assertTrue(any(
                        "earlier failure" in note and
                        "ordinary primary" in note
                        for note in getattr(cleanup, "__notes__", ())))

    def test_every_simple_owner_keeps_ordinary_failures_fail_closed(
            self) -> None:
        for name, owner_exit in self._SIMPLE_OWNER_EXITS:
            with self.subTest(owner=name):
                primary = ValueError(f"{name} ordinary primary")
                cleanup = OSError(f"{name} ordinary cleanup")
                owner = self._simple_owner(cleanup)
                with self.assertRaises(
                        provenance.BuildProvenanceError) as raised:
                    owner_exit(
                        owner, type(primary), primary,
                        primary.__traceback__)
                self.assertIn("ordinary primary", str(raised.exception))
                self.assertIn("ordinary cleanup", str(raised.exception))
                self.assertIs(raised.exception.__cause__, cleanup)

    def test_every_simple_owner_propagates_cleanup_without_primary(
            self) -> None:
        for cleanup_type in (OSError, KeyboardInterrupt, SystemExit):
            for name, owner_exit in self._SIMPLE_OWNER_EXITS:
                with self.subTest(
                        owner=name, cleanup=cleanup_type.__name__):
                    cleanup = cleanup_type(
                        f"{name} unaccompanied cleanup")
                    owner = self._simple_owner(cleanup)
                    with self.assertRaises(cleanup_type) as raised:
                        owner_exit(owner, None, None, None)
                    self.assertIs(raised.exception, cleanup)

    @staticmethod
    def _construct_descriptor(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        owner = provenance._OwnedDescriptor()
        with mock.patch.object(
                provenance.os, "open", return_value=91), \
                mock.patch.object(
                    provenance.os, "fstat", side_effect=primary), \
                OwnerExceptionPrecedenceTests.\
                    _patch_descriptor_cleanup_failure(cleanup):
            owner.open("/injected/descriptor", os.O_RDONLY)

    @staticmethod
    def _construct_inotify(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        library = mock.Mock()
        library.inotify_init1.return_value = 91
        with mock.patch.object(
                provenance.ctypes, "CDLL", return_value=library), \
                mock.patch.object(
                    provenance.os, "fdopen", side_effect=primary), \
                OwnerExceptionPrecedenceTests.\
                    _patch_descriptor_cleanup_failure(cleanup):
            provenance._InotifyMutationGuard("injected inotify")

    @staticmethod
    def _construct_artifact(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        with mock.patch.object(
                provenance, "_linux_memfd_create", return_value=91), \
                mock.patch.object(
                    provenance.os, "fdopen", side_effect=primary), \
                OwnerExceptionPrecedenceTests.\
                    _patch_descriptor_cleanup_failure(cleanup):
            provenance._ReplayArtifactSink(
                "/injected/artifact", "injected artifact")

    @staticmethod
    def _construct_file(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        class Guard:
            descriptor = 91

            @staticmethod
            def _absolute_lexical(path: Path) -> Path:
                return path

            @staticmethod
            def add_file_path(path: Path) -> None:
                del path
                raise primary

            @staticmethod
            def _close_without_verification() -> None:
                raise cleanup

        with tempfile.TemporaryDirectory(
                prefix="leo2-owner-file-constructor-") as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"input")
            with mock.patch.object(
                    provenance, "_InotifyMutationGuard",
                    return_value=Guard()):
                provenance._RetainedFileSnapshot(path, "injected file")

    @staticmethod
    def _construct_directory(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        class Guard:
            descriptor = 91

            @staticmethod
            def _absolute_lexical(path: Path) -> Path:
                return path

            @staticmethod
            def add_directory_path(path: Path) -> None:
                del path
                raise primary

            @staticmethod
            def _close_without_verification() -> None:
                raise cleanup

        with tempfile.TemporaryDirectory(
                prefix="leo2-owner-directory-constructor-") as directory:
            with mock.patch.object(
                    provenance, "_InotifyMutationGuard",
                    return_value=Guard()):
                provenance._RetainedDirectoryTree(
                    Path(directory), "injected directory")

    @staticmethod
    def _construct_git(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        class RetainedGitDirectory:
            descriptor = 91

            def __init__(self, path: Path, label: str) -> None:
                del label
                self.resolved = path.resolve(strict=True)

            def __enter__(self) -> "RetainedGitDirectory":
                return self

            def __exit__(
                self, exc_type: object, exc: object, tb: object,
            ) -> None:
                del exc_type, exc, tb
                raise cleanup

            @staticmethod
            def verify() -> None:
                raise primary

        with tempfile.TemporaryDirectory(
                prefix="leo2-owner-git-constructor-") as directory:
            source = Path(directory)
            (source / ".git").mkdir()
            with mock.patch.object(
                    provenance, "_RetainedDirectoryTree",
                    RetainedGitDirectory):
                provenance._RetainedGitMetadata(source)

    @staticmethod
    def _construct_source(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        class Guard:
            descriptor = 91

            @staticmethod
            def add_tree(path: Path) -> None:
                del path
                raise primary

            @staticmethod
            def _close_without_verification() -> None:
                raise cleanup

        with mock.patch.object(
                provenance, "_capture_tracked_source_tree",
                return_value={"files": [], "total_bytes": 0}), \
                mock.patch.object(
                    provenance, "_InotifyMutationGuard",
                    return_value=Guard()):
            provenance._RetainedPrivateSourceTree(
                Path("/injected/source"), Path("/injected/destination"))

    @staticmethod
    def _construct_proc(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        with mock.patch.object(
                provenance.os, "listdir", return_value=[]), \
                mock.patch.object(
                    provenance.os, "getpid", side_effect=primary), \
                mock.patch.object(
                    provenance._ProcProcessSnapshot, "close",
                    side_effect=cleanup):
            provenance._ProcProcessSnapshot()

    @staticmethod
    def _construct_recipe(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-owner-recipe-constructor-") as directory:
            with mock.patch.object(
                    provenance, "_replay_recipe_candidates",
                    side_effect=primary), mock.patch.object(
                        provenance._TemporaryReplayRecipeTransport, "close",
                        side_effect=cleanup):
                provenance._TemporaryReplayRecipeTransport(
                    Path(directory), {}, set())

    @staticmethod
    def _construct_containment(
        primary: BaseException, cleanup: BaseException,
    ) -> None:
        with mock.patch.object(
                provenance.os, "listdir", return_value=["1"]), \
                mock.patch.object(
                    provenance, "_validate_pidfd_support"), \
                mock.patch.object(
                    provenance, "_get_child_subreaper", return_value=0), \
                mock.patch.object(
                    provenance, "_set_child_subreaper",
                    side_effect=(primary, cleanup)):
            provenance._LinuxDescendantContainment().__enter__()

    _CONSTRUCTORS = (
        ("descriptor", _construct_descriptor),
        ("inotify", _construct_inotify),
        ("artifact", _construct_artifact),
        ("file", _construct_file),
        ("directory", _construct_directory),
        ("git", _construct_git),
        ("source", _construct_source),
        ("proc", _construct_proc),
        ("recipe", _construct_recipe),
        ("containment", _construct_containment),
    )

    def test_every_constructor_retains_earlier_terminal(self) -> None:
        for name, construct in self._CONSTRUCTORS:
            with self.subTest(owner=name):
                primary = KeyboardInterrupt(
                    f"{name} constructor interrupt")
                cleanup = OSError(f"{name} constructor cleanup")
                with self.assertRaises(KeyboardInterrupt) as raised:
                    construct(primary, cleanup)
                self.assertIs(raised.exception, primary)
                self.assertTrue(any(
                    "later cleanup failure" in note and
                    "constructor cleanup" in note
                    for note in getattr(primary, "__notes__", ())))

    def test_every_constructor_promotes_cleanup_terminal(self) -> None:
        for name, construct in self._CONSTRUCTORS:
            with self.subTest(owner=name):
                primary = ValueError(
                    f"{name} ordinary constructor failure")
                cleanup = SystemExit(
                    f"{name} constructor terminal cleanup")
                with self.assertRaises(SystemExit) as raised:
                    construct(primary, cleanup)
                self.assertIs(raised.exception, cleanup)
                self.assertTrue(any(
                    "earlier failure" in note and
                    "ordinary constructor failure" in note
                    for note in getattr(cleanup, "__notes__", ())))

    def test_every_constructor_retains_first_of_two_terminals(self) -> None:
        for name, construct in self._CONSTRUCTORS:
            with self.subTest(owner=name):
                primary = SystemExit(
                    f"{name} constructor primary exit")
                cleanup = KeyboardInterrupt(
                    f"{name} constructor later interrupt")
                with self.assertRaises(SystemExit) as raised:
                    construct(primary, cleanup)
                self.assertIs(raised.exception, primary)
                self.assertTrue(any(
                    "later cleanup failure" in note and
                    "later interrupt" in note
                    for note in getattr(primary, "__notes__", ())))

    def test_every_constructor_ordinary_failures_remain_fail_closed(
            self) -> None:
        for name, construct in self._CONSTRUCTORS:
            with self.subTest(owner=name):
                primary = ValueError(
                    f"{name} ordinary constructor primary")
                cleanup = OSError(
                    f"{name} ordinary constructor cleanup")
                with self.assertRaises(
                        provenance.BuildProvenanceError) as raised:
                    construct(primary, cleanup)
                self.assertIn(
                    "ordinary constructor primary",
                    str(raised.exception))
                self.assertIn(
                    "ordinary constructor cleanup",
                    str(raised.exception))

    @staticmethod
    def _containment() -> object:
        containment = object.__new__(
            provenance._LinuxDescendantContainment)
        containment.active = True
        containment.process = None
        containment.proven_empty = False
        containment.previous_subreaper = 0
        containment.leader = None
        containment.known = set()
        containment.pidfds = {}
        containment.procfds = {}
        return containment

    def test_containment_retains_first_terminal_across_every_phase(
            self) -> None:
        primary = KeyboardInterrupt("containment primary interrupt")
        cleanup = ValueError("containment ordinary cleanup")
        restore = SystemExit("containment later exit")
        descriptor = OSError("containment descriptor cleanup")
        containment = self._containment()
        with mock.patch.object(
                provenance._LinuxDescendantContainment,
                "terminate_unattached_and_reap",
                side_effect=cleanup), mock.patch.object(
                    provenance, "_set_child_subreaper",
                    side_effect=restore), mock.patch.object(
                        provenance._LinuxDescendantContainment,
                        "_close_pidfds", side_effect=descriptor):
            with self.assertRaises(KeyboardInterrupt) as raised:
                provenance._LinuxDescendantContainment.__exit__(
                    containment, type(primary), primary,
                    primary.__traceback__)
        self.assertIs(raised.exception, primary)
        notes = getattr(primary, "__notes__", ())
        self.assertTrue(any("ordinary cleanup" in note for note in notes))
        self.assertTrue(any("later exit" in note for note in notes))
        self.assertTrue(any("descriptor cleanup" in note for note in notes))

    def test_containment_promotes_first_cleanup_terminal(self) -> None:
        primary = ValueError("containment ordinary primary")
        cleanup = SystemExit("containment terminal cleanup")
        restore = KeyboardInterrupt("containment later interrupt")
        descriptor = OSError("containment descriptor cleanup")
        containment = self._containment()
        with mock.patch.object(
                provenance._LinuxDescendantContainment,
                "terminate_unattached_and_reap",
                side_effect=cleanup), mock.patch.object(
                    provenance, "_set_child_subreaper",
                    side_effect=restore), mock.patch.object(
                        provenance._LinuxDescendantContainment,
                        "_close_pidfds", side_effect=descriptor):
            with self.assertRaises(SystemExit) as raised:
                provenance._LinuxDescendantContainment.__exit__(
                    containment, type(primary), primary,
                    primary.__traceback__)
        self.assertIs(raised.exception, cleanup)
        notes = getattr(cleanup, "__notes__", ())
        self.assertTrue(any("ordinary primary" in note for note in notes))
        self.assertTrue(any("later interrupt" in note for note in notes))
        self.assertTrue(any("descriptor cleanup" in note for note in notes))

    def test_containment_ordinary_phases_remain_fail_closed(self) -> None:
        primary = ValueError("containment ordinary primary")
        cleanup = RuntimeError("containment ordinary cleanup")
        restore = OSError("containment ordinary restore")
        descriptor = LookupError("containment ordinary descriptor")
        containment = self._containment()
        with mock.patch.object(
                provenance._LinuxDescendantContainment,
                "terminate_unattached_and_reap",
                side_effect=cleanup), mock.patch.object(
                    provenance, "_set_child_subreaper",
                    side_effect=restore), mock.patch.object(
                        provenance._LinuxDescendantContainment,
                        "_close_pidfds", side_effect=descriptor):
            with self.assertRaises(
                    provenance.BuildProvenanceError) as raised:
                provenance._LinuxDescendantContainment.__exit__(
                    containment, type(primary), primary,
                    primary.__traceback__)
        message = str(raised.exception)
        for expected in (
                "ordinary primary", "ordinary cleanup",
                "ordinary restore", "ordinary descriptor"):
            self.assertIn(expected, message)

    def test_containment_cleanup_without_primary_is_fail_closed(self) -> None:
        for cleanup_type, expected_type in (
                (OSError, provenance.BuildProvenanceError),
                (KeyboardInterrupt, KeyboardInterrupt),
                (SystemExit, SystemExit)):
            with self.subTest(cleanup=cleanup_type.__name__):
                cleanup = cleanup_type(
                    "containment unaccompanied cleanup")
                containment = self._containment()
                with mock.patch.object(
                        provenance._LinuxDescendantContainment,
                        "terminate_unattached_and_reap",
                        side_effect=cleanup), mock.patch.object(
                            provenance, "_set_child_subreaper"), \
                        mock.patch.object(
                            provenance._LinuxDescendantContainment,
                            "_close_pidfds"):
                    with self.assertRaises(expected_type) as raised:
                        provenance._LinuxDescendantContainment.__exit__(
                            containment, None, None, None)
                if cleanup_type is not OSError:
                    self.assertIs(raised.exception, cleanup)


class BoundedDescendantContainmentTests(unittest.TestCase):
    def test_restores_both_prior_subreaper_states(self) -> None:
        original = provenance._get_child_subreaper()
        try:
            for prior in (0, 1):
                with self.subTest(prior=prior):
                    provenance._set_child_subreaper(prior)
                    self.assertEqual(
                        provenance._run(
                            [str(tool_path("true"))],
                            f"subreaper restoration from {prior}",
                            timeout=3),
                        b"")
                    self.assertEqual(
                        provenance._get_child_subreaper(), prior)
        finally:
            provenance._set_child_subreaper(original)

    def test_direct_argv_allows_multiline_data_on_every_exit_path(self) -> None:
        previous = provenance._get_child_subreaper()
        output = provenance._run(
            [sys.executable, "-c",
             "import sys\nsys.stdout.write('multiline-ok\\n')\n",
             ""],
            "multiline success", timeout=3)
        self.assertEqual(output, b"multiline-ok\n")
        self.assertEqual(provenance._get_child_subreaper(), previous)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "failed with rc=9"):
            provenance._run(
                [sys.executable, "-c",
                 "import sys\nsys.stderr.write('multiline-error\\n')\n"
                 "raise SystemExit(9)\n"],
                "multiline failure", timeout=3)
        self.assertEqual(provenance._get_child_subreaper(), previous)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "exceeded 0.200 seconds"):
            provenance._run(
                [sys.executable, "-c",
                 "import time\ntime.sleep(60)\n"],
                "multiline timeout", timeout=0.2)
        self.assertEqual(provenance._get_child_subreaper(), previous)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "argv is invalid"):
            provenance._run(
                [sys.executable, "bad\0argument"],
                "NUL argv", timeout=3)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError, "argv is invalid"):
            provenance._run(
                sys.executable, "string is not argv", timeout=3)

    def test_success_kills_detached_fd_holder_before_return(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-detached-success-") as directory:
            root = Path(directory)
            helper = root / "detach.py"
            pid_path = root / "daemon.pid"
            lock_path = root / "campaign.lock"
            write_detached_helper(helper)
            lock_descriptor = os.open(
                lock_path, os.O_RDWR | os.O_CREAT, 0o600)
            competitor = -1
            try:
                fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
                provenance._run(
                    [sys.executable, str(helper), str(pid_path)],
                    "detached success helper", timeout=3,
                    inherited_descriptors=(lock_descriptor,))
                daemon_pid = int(pid_path.read_text(encoding="ascii"))
                self.assertTrue(process_gone(daemon_pid))

                # Close the coordinator's copy.  If containment returned while
                # a detached descendant still held the inherited open-file
                # description, this independent lock attempt would fail.
                os.close(lock_descriptor)
                lock_descriptor = -1
                competitor = os.open(lock_path, os.O_RDWR)
                fcntl.flock(
                    competitor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                if competitor >= 0:
                    os.close(competitor)
                if lock_descriptor >= 0:
                    os.close(lock_descriptor)

    def test_success_and_failure_drain_retained_pipes_after_leader_exit(
            self) -> None:
        for mode, pattern in (
                ("retained-pipes", None),
                ("retained-pipes-failure", "retained-pipe-error")):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                    prefix="leo2-provenance-retained-pipes-") as directory:
                root = Path(directory)
                helper = root / "detach.py"
                pid_path = root / "daemon.pid"
                write_detached_helper(helper)
                if pattern is None:
                    output = provenance._run(
                        [sys.executable, str(helper), str(pid_path), mode],
                        "retained-pipe success", timeout=3)
                    self.assertEqual(output, b"retained-pipe-output\n")
                else:
                    with self.assertRaisesRegex(
                            provenance.BuildProvenanceError, pattern):
                        provenance._run(
                            [sys.executable, str(helper), str(pid_path), mode],
                            "retained-pipe failure", timeout=3)
                daemon_pid = int(pid_path.read_text(encoding="ascii"))
                self.assertTrue(process_gone(daemon_pid))

    def test_run_does_not_reap_leader_before_tree_cleanup(self) -> None:
        observed_returncodes: list[int | None] = []
        original = (
            provenance._LinuxDescendantContainment.terminate_and_reap)

        def checking_terminate_and_reap(
            containment: provenance._LinuxDescendantContainment,
        ) -> None:
            self.assertIsNotNone(containment.process)
            observed_returncodes.append(containment.process.returncode)
            original(containment)

        with mock.patch.object(
                provenance._LinuxDescendantContainment,
                "terminate_and_reap", new=checking_terminate_and_reap):
            provenance._run(
                [str(tool_path("true"))],
                "deferred leader reap", timeout=3)
        self.assertEqual(observed_returncodes, [None])

    def test_failure_and_timeout_kill_setsid_double_fork_descendants(
            self) -> None:
        for mode, timeout, pattern in (
                ("failure", 3.0, "failed with rc=7"),
                ("timeout", 0.2, "exceeded 0.200 seconds")):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                    prefix=f"leo2-provenance-detached-{mode}-") as directory:
                root = Path(directory)
                helper = root / "detach.py"
                pid_path = root / "daemon.pid"
                write_detached_helper(helper)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError, pattern):
                    provenance._run(
                        [sys.executable, str(helper), str(pid_path), mode],
                        f"detached {mode} helper", timeout=timeout)
                daemon_pid = int(pid_path.read_text(encoding="ascii"))
                self.assertTrue(process_gone(daemon_pid))

    def test_popen_exception_reaps_unattached_double_fork_and_restores(
            self) -> None:
        previous = provenance._get_child_subreaper()
        with tempfile.TemporaryDirectory(
                prefix="leo2-provenance-unattached-") as directory:
            pid_path = Path(directory) / "daemon.pid"

            def spawn_then_fail(*_args, **_kwargs):
                ready_read, ready_write = os.pipe()
                child = os.fork()
                if child == 0:
                    try:
                        os.close(ready_read)
                        os.setsid()
                        grandchild = os.fork()
                        if grandchild:
                            os._exit(0)
                        pid_path.write_text(
                            str(os.getpid()), encoding="ascii")
                        os.write(ready_write, b"1")
                        os.close(ready_write)
                        time.sleep(60)
                    finally:
                        os._exit(0)
                os.close(ready_write)
                ready = os.read(ready_read, 1)
                os.close(ready_read)
                os.waitpid(child, 0)
                if ready != b"1":
                    raise RuntimeError("detached helper was not ready")
                raise OSError("synthetic Popen failure after fork")

            with mock.patch.object(
                    provenance.subprocess, "Popen",
                    side_effect=spawn_then_fail):
                with self.assertRaisesRegex(
                        OSError, "synthetic Popen failure"):
                    provenance._run(
                        ["/unused/direct-command"],
                        "unattached Popen failure", timeout=3)
            daemon_pid = int(pid_path.read_text(encoding="ascii"))
            self.assertTrue(process_gone(daemon_pid))
            self.assertEqual(
                provenance._get_child_subreaper(), previous)

    def test_attach_exceptions_kill_leader_and_restore_subreaper(self) -> None:
        previous = provenance._get_child_subreaper()
        original_attach = provenance._LinuxDescendantContainment.attach
        for after_attach in (False, True):
            captured: list[int] = []

            def fail_attach(
                containment: provenance._LinuxDescendantContainment,
                process: subprocess.Popen[bytes],
            ) -> None:
                captured.append(process.pid)
                if after_attach:
                    original_attach(containment, process)
                raise RuntimeError("synthetic attach failure")

            with self.subTest(after_attach=after_attach), \
                    mock.patch.object(
                        provenance._LinuxDescendantContainment, "attach",
                        new=fail_attach):
                with self.assertRaisesRegex(
                        RuntimeError, "synthetic attach failure"):
                    provenance._run(
                        [str(tool_path("sleep")), "60"],
                        "attach failure", timeout=3)
            self.assertEqual(len(captured), 1)
            self.assertTrue(process_gone(captured[0]))
            self.assertEqual(
                provenance._get_child_subreaper(), previous)

    def test_selector_setup_exception_cleans_attached_process(self) -> None:
        previous = provenance._get_child_subreaper()
        captured: list[int] = []
        original_popen = subprocess.Popen

        def capturing_popen(*args: object, **kwargs: object):
            process = original_popen(*args, **kwargs)
            captured.append(process.pid)
            return process

        with mock.patch.object(
                provenance.subprocess, "Popen",
                side_effect=capturing_popen), \
                mock.patch.object(
                    provenance.selectors.DefaultSelector, "register",
                    side_effect=OSError("synthetic selector failure")):
            with self.assertRaisesRegex(
                    OSError, "synthetic selector failure"):
                provenance._run(
                    [str(tool_path("sleep")), "60"],
                    "selector setup failure", timeout=3)
        self.assertEqual(len(captured), 1)
        self.assertTrue(process_gone(captured[0]))
        self.assertEqual(provenance._get_child_subreaper(), previous)

    def test_pidfd_retention_rejects_same_tick_reuse_after_open(self) -> None:
        containment = provenance._LinuxDescendantContainment()
        pid = 424242
        starttime = 777
        identity = (pid, starttime, 7, 700)
        observed = (os.getpid(), 100, 100, starttime, "R", 7, 700)
        # Numeric starttime and proc inode deliberately collide.  The retained
        # old proc-directory FD nevertheless becomes unreadable after reuse.
        replacement = observed
        with mock.patch.object(
                provenance, "_read_proc_process_record_descriptor",
                side_effect=(observed, None)), \
                mock.patch.object(
                    provenance, "_open_proc_process_record",
                    return_value=(replacement, 82)), \
                mock.patch.object(provenance.os, "dup", return_value=81), \
                mock.patch.object(
                    provenance, "_linux_pidfd_open", return_value=91), \
                mock.patch.object(provenance.os, "close") as close:
            self.assertFalse(
                containment._retain_pidfd(
                    identity, observed, observed_descriptor=80))
        self.assertEqual(
            {call.args[0] for call in close.call_args_list},
            {81, 82, 91})
        self.assertNotIn(identity, containment.pidfds)
        self.assertNotIn(identity, containment.procfds)

        containment.pidfds[identity] = 92
        containment.procfds[identity] = 93

        def consume_fake_descriptor(
                descriptor: int, context: str,
                consumed: list[bool]) -> None:
            self.assertIn(descriptor, (92, 93))
            self.assertEqual(context, "retained mapping descriptor")
            consumed[0] = True

        with mock.patch.object(
                provenance, "_linux_pidfd_signal") as pidfd_signal, \
                mock.patch.object(
                    provenance, "_close_descriptor_with_ofd_guard",
                    side_effect=consume_fake_descriptor):
            containment._signal_retained({identity}, provenance.signal.SIGKILL)
            containment._close_pidfds()
        pidfd_signal.assert_called_once_with(92, provenance.signal.SIGKILL)
        self.assertFalse(containment.pidfds)

    def test_pidfd_retention_closes_handle_when_postcheck_raises(self) -> None:
        containment = provenance._LinuxDescendantContainment()
        pid = 424243
        starttime = 778
        identity = (pid, starttime, 8, 800)
        observed = (os.getpid(), 100, 100, starttime, "R", 8, 800)
        with mock.patch.object(
                provenance, "_read_proc_process_record_descriptor",
                side_effect=(
                    observed,
                    provenance.BuildProvenanceError(
                        "synthetic post-open procfs failure"),
                )), \
                mock.patch.object(provenance.os, "dup", return_value=83), \
                mock.patch.object(
                    provenance, "_linux_pidfd_open", return_value=93), \
                mock.patch.object(provenance.os, "close") as close:
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "synthetic post-open"):
                containment._retain_pidfd(
                    identity, observed, observed_descriptor=80)
        self.assertEqual(
            {call.args[0] for call in close.call_args_list}, {83, 93})
        self.assertNotIn(identity, containment.pidfds)
        self.assertNotIn(identity, containment.procfds)

    def test_proc_directory_remains_held_through_pidfd_open_and_cleanup(
            self) -> None:
        containment = provenance._LinuxDescendantContainment()
        pid = 424246
        starttime = 781
        observed = (
            containment.runner_pid, 102, 102, starttime, "R", 11, 1100)
        identity = (pid, starttime, 11, 1100)
        closed: list[int] = []

        def pidfd_open(unused_pid: int) -> int:
            self.assertNotIn(86, closed)
            return 96

        with mock.patch.object(
                provenance, "_read_proc_process_record_descriptor",
                side_effect=(observed, observed)), \
                mock.patch.object(
                    provenance, "_open_proc_process_record",
                    return_value=(observed, 87)), \
                mock.patch.object(provenance.os, "dup", return_value=86), \
                mock.patch.object(
                    provenance, "_linux_pidfd_open",
                    side_effect=pidfd_open), \
                mock.patch.object(
                    provenance.os, "close",
                    side_effect=lambda descriptor: closed.append(descriptor)):
            self.assertTrue(containment._retain_pidfd(
                identity, observed, observed_descriptor=80))
            self.assertEqual(containment.pidfds[identity], 96)
            self.assertEqual(containment.procfds[identity], 86)
            self.assertNotIn(86, closed)
            containment._close_pidfds()
        self.assertEqual(set(closed), {86, 87, 96})

    def test_unretained_observation_cannot_rebind_by_pid_number(self) -> None:
        containment = provenance._LinuxDescendantContainment()
        pid = 424244
        starttime = 779
        direct = (
            containment.runner_pid, 100, 100, starttime, "R", 9, 900)
        with mock.patch.object(
                containment, "_retain_pidfd", return_value=False):
            first = containment._discover({pid: direct})
            self.assertEqual(first, {(pid, starttime, 9, 900)})
            self.assertNotIn((pid, starttime, 9, 900), containment.known)

            unrelated = (
                containment.runner_pid + 1000, 100, 100, starttime, "R",
                9, 901)
            second = containment._discover({pid: unrelated})
            self.assertEqual(second, set())
            self.assertNotIn((pid, starttime, 9, 900), containment.known)

    def test_same_starttime_different_proc_inode_cannot_rebind(self) -> None:
        containment = provenance._LinuxDescendantContainment()
        pid = 424245
        starttime = 780
        original = (
            containment.runner_pid, 101, 101, starttime, "R", 10, 1000)
        replacement = (
            containment.runner_pid, 101, 101, starttime, "R", 10, 1001)
        identity = (pid, starttime, 10, 1000)
        with mock.patch.object(
                provenance, "_read_proc_process_record_descriptor",
                side_effect=(original, original)), \
                mock.patch.object(
                    provenance, "_open_proc_process_record",
                    return_value=(replacement, 84)), \
                mock.patch.object(provenance.os, "dup", return_value=85), \
                mock.patch.object(
                    provenance, "_linux_pidfd_open", return_value=94), \
                mock.patch.object(provenance.os, "close") as close:
            self.assertFalse(
                containment._retain_pidfd(
                    identity, original, observed_descriptor=80))
        self.assertEqual(
            {call.args[0] for call in close.call_args_list},
            {84, 85, 94})
        self.assertNotIn(identity, containment.pidfds)
        self.assertNotIn(identity, containment.procfds)


def reproducible_record() -> dict[str, object]:
    source = (SOURCE_ROOT / "leopard2.cpp").resolve(strict=True)
    build = Path("/tmp/leopard2-provenance-candidate")
    arguments = [
        "/usr/bin/c++", "-Wall", "-O3",
        "-o", "CMakeFiles/leopard.dir/leopard2.cpp.o",
        "-c", str(source),
    ]

    def identity(label: str) -> dict[str, object]:
        return {"sha256": label * 64, "size": len(label)}

    return {
        "schema": "leopard2-production-build-closure/v1",
        "build_root": str(build),
        "physical_source_root": str(SOURCE_ROOT),
        "source_root": str(SOURCE_ROOT),
        "executable_target": "bench_leopard2",
        "validated_cache": {"profile": "canonical"},
        "tracked_source_manifest": {
            "schema": "leopard2-tracked-source-tree/v1",
            "total_bytes": 123,
            "files": [{
                "path": "leopard2.cpp", "sha256": "s" * 64,
                "size": 123, "mode": 0o644,
            }],
            "git": {
                "commit": "1" * 40, "tree": "2" * 40, "dirty": True,
                "status_sha256": "3" * 64,
            },
            "git_tool": {
                "path": "/usr/bin/git", "sha256": "4" * 64,
                "size": 1, "device": 1, "inode": 1,
            },
        },
        "source_object_compile_closure": [{
            "role": "archive",
            "source": {
                "path": str(source), "sha256": "s" * 64, "size": 123,
            },
            "object": {
                "path": "${BUILD_ROOT}/leopard2.cpp.o",
                "sha256": "o" * 64, "size": 456,
            },
            "compile_entry": {
                "arguments": list(arguments),
                "normalized_arguments": provenance._normalize_build_argv(
                    arguments, build, SOURCE_ROOT),
                "representation": "arguments",
            },
            "flag_profile": "portable-core",
        }],
        "archive_link_recipe": identity("a"),
        "executable_link_recipe": identity("e"),
        "archive_link_commands": [
            ["/usr/bin/ar", "qc", "libleopard.a",
             "CMakeFiles/leopard.dir/leopard2.cpp.o"],
            ["/usr/bin/ranlib", "libleopard.a"],
        ],
        "executable_link_command": [
            "/usr/bin/c++", "-O3",
            "CMakeFiles/bench_leopard2.dir/benchmark.cpp.o",
            "-o", "bench_leopard2", "libleopard.a",
        ],
        "archive_member_identities": [{
            "member": "leopard2.cpp.o", "sha256": "o" * 64, "size": 456,
        }],
        "archive": identity("r"),
        "executable": identity("x"),
        "compiler": identity("c"),
        "compiler_version_sha256": "v" * 64,
        "c_compiler": identity("d"),
        "make_program": identity("m"),
        "linker": identity("l"),
        "compiler_subtools": [{
            "language": "c++", "role": "cc1plus",
            "identity": identity("u"),
        }],
        "archiver": identity("i"),
        "ranlib": identity("n"),
        "benchmark_git": identity("g"),
    }


class ExactRecipeComparisonTests(unittest.TestCase):
    def test_compile_argv_changes_fail_even_when_object_bytes_match(
            self) -> None:
        candidate = reproducible_record()
        mutations = {
            "unused definition": lambda argv: argv.insert(1, "-DUNUSED=1"),
            "prepended command": lambda argv: argv.insert(0, "/usr/bin/true"),
            "reordered flags": lambda argv: argv.__setitem__(
                slice(1, 3), reversed(argv[1:3])),
            "substituted driver": lambda argv: argv.__setitem__(
                0, "/usr/bin/clang++"),
        }
        for name, mutate in mutations.items():
            rebuilt = copy.deepcopy(candidate)
            entry = rebuilt["source_object_compile_closure"][0][
                "compile_entry"]
            mutate(entry["arguments"])
            mutate(entry["normalized_arguments"])
            with self.subTest(name=name):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "compile argv differs"):
                    provenance.compare_reproducible_builds(
                        candidate, rebuilt)

    def test_compile_argv_normalization_cannot_be_forged(self) -> None:
        candidate = reproducible_record()
        rebuilt = copy.deepcopy(candidate)
        rebuilt["source_object_compile_closure"][0]["compile_entry"][
            "arguments"].insert(1, "-DUNUSED=1")
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "compile argv normalization is invalid"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

    def test_build_root_substring_does_not_hide_argv_difference(self) -> None:
        candidate = reproducible_record()
        rebuilt = copy.deepcopy(candidate)
        candidate_root = Path(str(candidate["build_root"]))
        rebuilt_root = Path("/tmp/leopard2-provenance-rebuilt")
        rebuilt["build_root"] = str(rebuilt_root)
        candidate_entry = candidate[
            "source_object_compile_closure"][0]["compile_entry"]
        rebuilt_entry = rebuilt[
            "source_object_compile_closure"][0]["compile_entry"]
        candidate_argument = f"-DUNUSED={candidate_root}-suffix"
        rebuilt_argument = f"-DUNUSED={rebuilt_root}-suffix"
        candidate_entry["arguments"].insert(1, candidate_argument)
        rebuilt_entry["arguments"].insert(1, rebuilt_argument)
        candidate_entry["normalized_arguments"] = \
            provenance._normalize_build_argv(
                candidate_entry["arguments"], candidate_root,
                Path(candidate["physical_source_root"]))
        rebuilt_entry["normalized_arguments"] = \
            provenance._normalize_build_argv(
                rebuilt_entry["arguments"], rebuilt_root,
                Path(rebuilt["physical_source_root"]))
        self.assertNotEqual(
            candidate_entry["normalized_arguments"],
            rebuilt_entry["normalized_arguments"])
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "compile argv differs"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

    def test_compile_recipe_canonicalizes_safe_spelling_only(
            self) -> None:
        candidate = reproducible_record()
        entry = candidate["source_object_compile_closure"][0][
            "compile_entry"]
        entry["representation"] = "command"
        entry["command"] = shlex.join(entry["arguments"])
        entry["normalized_command"] = shlex.join(
            entry["normalized_arguments"])

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["command"] = rebuilt_entry["command"].replace(
            " -Wall ", "  -Wall ", 1)
        provenance.compare_reproducible_builds(candidate, rebuilt)

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["command"] = rebuilt_entry["command"].replace(
            " -Wall ", " -Wextra ", 1)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "compile recipe normalization is invalid"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["command"] = rebuilt_entry["command"].replace(
            " -Wall ", " '-Wall -Wextra' ", 1)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "compile recipe normalization is invalid"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["command"] = rebuilt_entry["command"].replace(
            " -Wall ", " -Wall; /usr/bin/true ", 1)
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "shell control"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["representation"] = "arguments"
        del rebuilt_entry["command"]
        del rebuilt_entry["normalized_command"]
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "representation differs"):
            provenance.compare_reproducible_builds(candidate, rebuilt)

    def test_archive_recipe_mutations_fail(self) -> None:
        candidate = reproducible_record()
        mutations = {
            "prepended command": lambda commands: commands.insert(
                0, ["/usr/bin/true"]),
            "reordered objects": lambda commands: commands[0].extend(
                commands[0][3:4]),
            "substituted driver": lambda commands: commands[0].__setitem__(
                0, "/usr/bin/llvm-ar"),
        }
        for name, mutate in mutations.items():
            rebuilt = copy.deepcopy(candidate)
            mutate(rebuilt["archive_link_commands"])
            with self.subTest(name=name):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "archive link recipe"):
                    provenance.compare_reproducible_builds(
                        candidate, rebuilt)

    def test_executable_recipe_mutations_fail(self) -> None:
        candidate = reproducible_record()
        mutations = {
            "appended linker flag": lambda command: command.append(
                "-Wl,--as-needed"),
            "wrapper": lambda command: command.insert(0, "/usr/bin/env"),
            "reordered inputs": lambda command: command.__setitem__(
                slice(1, 3), reversed(command[1:3])),
            "substituted driver": lambda command: command.__setitem__(
                0, "/usr/bin/clang++"),
        }
        for name, mutate in mutations.items():
            rebuilt = copy.deepcopy(candidate)
            mutate(rebuilt["executable_link_command"])
            with self.subTest(name=name):
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError,
                        "executable link recipe semantics differ"):
                    provenance.compare_reproducible_builds(
                        candidate, rebuilt)

    def test_raw_recipe_identity_is_compared(self) -> None:
        candidate = reproducible_record()
        rebuilt = copy.deepcopy(candidate)
        rebuilt["archive_link_recipe"]["sha256"] = "z" * 64
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "archive link recipe identity differs"):
            provenance.compare_reproducible_builds(candidate, rebuilt)


if __name__ == "__main__":
    unittest.main()
