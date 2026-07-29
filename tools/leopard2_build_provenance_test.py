#!/usr/bin/env python3
"""Regression tests for the shared Leopard2 build-provenance helper."""

from __future__ import annotations

import copy
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import shlex
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

    def test_owned_descriptor_retains_same_fd_when_close_did_not_run(
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
        self.assertEqual(owner.descriptor, descriptor)
        os.fstat(descriptor)
        owner.close()
        self.assertEqual(owner.descriptor, -1)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

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
        armed = False

        class InitFunction:
            argtypes = None
            restype = None

            def __call__(self, unused_flags) -> int:
                nonlocal armed
                descriptor = os.open("/dev/null", os.O_RDONLY)
                captured.append(descriptor)
                armed = True
                return descriptor

        class Library:
            inotify_init1 = InitFunction()

        init_code = provenance._InotifyMutationGuard.__init__.__code__

        def interrupt_after_init(frame, event: str, argument):
            del argument
            if armed and event == "line" and frame.f_code is init_code:
                raise KeyboardInterrupt(
                    "injected post-inotify-init interruption")
            return interrupt_after_init

        previous_trace = sys.gettrace()
        try:
            with mock.patch.object(
                    provenance.ctypes, "CDLL", return_value=Library()):
                sys.settrace(interrupt_after_init)
                with self.assertRaisesRegex(
                        KeyboardInterrupt, "post-inotify-init"):
                    provenance._InotifyMutationGuard(
                        "interrupted inotify")
        finally:
            sys.settrace(previous_trace)
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
                del argument
                if armed and event == "line" and frame.f_code is property_code:
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
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
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
                        provenance.BuildProvenanceError,
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
                "source_root": str(SOURCE_ROOT),
                "executable_target": "bench_leopard2",
                "tracked_source_manifest": manifest,
                "validated_cache": {
                    name: str(tool) for name in (
                        "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER",
                        "CMAKE_AR", "CMAKE_RANLIB",
                        "CMAKE_MAKE_PROGRAM", "CMAKE_LINKER",
                        "LEO2_BENCHMARK_GIT_EXECUTABLE")
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
                "Leopard2BackendGFNI"))
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
                    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
                },
                library_sources=backend_sources),
            "portable-core")

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
        self.assertEqual(mapping, {"descriptor": descriptor})
        provenance._close_retained_mapping_descriptor(
            mapping, "descriptor")
        self.assertFalse(mapping)

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
        with mock.patch.object(
                provenance, "_linux_pidfd_signal") as pidfd_signal, \
                mock.patch.object(provenance.os, "close"):
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

    def test_raw_compile_recipe_changes_and_representation_changes_fail(
            self) -> None:
        candidate = reproducible_record()
        entry = candidate["source_object_compile_closure"][0][
            "compile_entry"]
        entry["representation"] = "command"
        entry["command"] = shlex.join(entry["arguments"])
        entry["normalized_command"] = provenance._normalize_root_token(
            provenance._normalize_build_token(
                entry["command"], Path(candidate["build_root"])),
            Path(candidate["physical_source_root"]), "${SOURCE_ROOT}")

        rebuilt = copy.deepcopy(candidate)
        rebuilt_entry = rebuilt["source_object_compile_closure"][0][
            "compile_entry"]
        rebuilt_entry["command"] = rebuilt_entry["command"].replace(
            " -Wall ", "  -Wall ", 1)
        rebuilt_entry["normalized_command"] = provenance._normalize_root_token(
            provenance._normalize_build_token(
                rebuilt_entry["command"], Path(rebuilt["build_root"])),
            Path(rebuilt["physical_source_root"]), "${SOURCE_ROOT}")
        with self.assertRaisesRegex(
                provenance.BuildProvenanceError,
                "raw compile recipe differs"):
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
