#!/usr/bin/env python3
"""Reversibly keep same-user work off an authoritative benchmark CPU pair.

This Linux-only helper is deliberately outside the codec.  It cannot exclude
kernel work or another UID, and it never replaces the benchmark runner's
independent ``/proc/stat`` zero-sibling-jiffy acceptance gate.

``run`` is a durable affinity transaction.  It records exact original masks,
removes the reserved CPUs from same-UID threads, launches one descriptor-pinned
command through a fork/pipe gate, monitors safe newcomers, terminates every
observed surviving child, and makes a bounded best-effort exact-mask restore.
The child cannot execute until its PID/start-time/session identity is durable.
Because Linux offers no unprivileged persistent containment here, a run that
actually changes another mask is never accepted as authoritative evidence: an
in-flight clone can publish a copied restricted mask after any finite scan.

``restore`` is the crash-recovery operation.  It first proves that the journal
belongs to this boot and PID namespace, then cleans known child identities and
retries pending or failed affinity restores.  If a released child or affinity
mutation existed, recovery remains explicitly failed: after the original
subreaper crashes, an unobserved cross-session orphan cannot be rediscovered
reliably without persistent cgroup-style containment.

The remaining Linux limitation is explicit: sched_getaffinity(2) and
sched_setaffinity(2) address a numeric TID.  Identity checks immediately before
and after each syscall detect ordinary reuse, but the kernel offers no pidfd-like
handle for a non-leader TID affinity operation.  Any observed post-syscall reuse
marks the transaction uncertain and permanently unsuitable as benchmark
evidence.
"""

from __future__ import print_function

import argparse
import ctypes
import datetime
import errno
import fcntl
import hashlib
import json
import os
import re
import select
import secrets
import shutil
import signal
import stat
import sys
import tempfile
import time
from dataclasses import dataclass, replace
from pathlib import Path


REPORT_SCHEMA = "leopard2-affinity-supervisor/v6"
ACCEPTANCE_SCHEMA = "leopard2-affinity-acceptance/v1"
BINDING_SCHEMA = "leopard2-affinity-main-binding/v2"
MAIN_MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v5"
MAIN_RAW_SCHEMA = "leopard2-main-compare-raw/v5"
MAIN_SUPERVISION_SCHEMA = "leopard2-main-supervision/v1"
DEFAULT_POLL_MS = 25
MAX_STABILIZE_PASSES = 64
LAUNCH_READY_TIMEOUT_SECONDS = 10.0
SIGNAL_GRACE_SECONDS = 2.0
KILL_GRACE_SECONDS = 2.0
MAX_IDENTITY_BYTES = 64 * 1024 * 1024
MUTATED_RESTORE_UNPROVABLE = (
    "authoritative completion is unavailable after affinity mutation: finite "
    "procfs scans cannot exclude an in-flight clone publishing with the "
    "restricted mask; persistent task containment is required")
CRASH_SCOPE_UNPROVABLE = (
    "crash recovery cannot prove released-child scope empty: an unjournaled "
    "descendant may have changed sessions before the subreaper exited; "
    "persistent task containment is required")
HEX256 = re.compile(r"[0-9a-f]{64}")
BOOT_ID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
GATED_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)

REPORT_KEYS = {
    "accepted", "boot_id", "child", "child_processes", "child_returncode", "command",
    "command_identity", "command_sha256", "error", "events", "execution",
    "finished_utc", "global_lock", "launch_cpus", "pid_namespace", "poll_ms",
    "process_provenance", "received_signal", "records", "reserved_cpus",
    "schema", "started_utc", "state", "subreaper_pid", "uid", "uncertainty",
}
RECORD_KEYS = {
    "move_count", "move_status", "original_cpus", "pid", "ppid",
    "process_device", "process_inode", "process_starttime_ticks", "provenance",
    "restricted_cpus", "restore_status", "session", "starttime_ticks",
    "task_device", "task_inode", "tid", "uid",
}
PROVENANCE_KEYS = {
    "original_cpus", "pid", "ppid", "process_device", "process_inode",
    "process_starttime_ticks", "source",
}
EVENT_KEYS = {"event", "pid", "tid", "value"}
CHILD_KEYS = {
    "pid", "process_device", "process_inode", "released", "session",
    "starttime_ticks",
}
CHILD_PROCESS_KEYS = {
    "pid", "ppid", "process_device", "process_inode", "session",
    "starttime_ticks",
}
EXECUTION_KEYS = {
    "child_finished_monotonic_ns", "child_ready_monotonic_ns", "nonce",
}
EXECUTION_NONCE_ENV = "LEO2_AFFINITY_EXECUTION_NONCE"
PID_NAMESPACE_KEYS = {"device", "inode"}
GLOBAL_LOCK_KEYS = {
    "device", "directory_device", "directory_inode", "inode", "lock", "mode",
    "path", "uid",
}
ACCEPTANCE_KEYS = {
    "committed_utc", "digest", "report_path", "report_sha256", "report_size",
    "schema",
}
COMMAND_IDENTITY_KEYS = {"cwd", "executable", "script"}
PATH_IDENTITY_KEYS = {"device", "inode", "mode", "path"}
FILE_IDENTITY_KEYS = {
    "argument", "device", "inode", "mode", "path", "sha256", "size",
}

REPORT_STATES = {
    "created", "prepared", "moving", "isolated", "launch_gated", "running",
    "restoring", "complete", "failed", "recovering", "recovered",
    "recovery_failed",
}
MOVE_STATUSES = {
    "planned", "moved", "inherited_restricted", "gone_before_move",
    "reused_before_move", "uncertain_after_move", "move_failed",
}
RESTORE_STATUSES = {
    "pending", "restored", "restored_after_recovery", "gone", "failed",
    "reused", "uncertain",
}
PROVENANCE_SOURCES = {"baseline", "parent_process", "observed"}
RECORD_PROVENANCE = {"baseline", "inherited", "observed"}
EVENT_NAMES = {
    "child_descendants", "child_session_cleanup", "gone_before_capture",
    "process_terminated", "reused_before_capture", "signal",
    "tid_reuse_uncertainty", "unsafe_newcomer", "write_failure",
}


class IsolationError(RuntimeError):
    pass


class SelfTestError(RuntimeError):
    pass


class ThreadGone(Exception):
    pass


class ThreadReused(Exception):
    pass


class ThreadMutationUncertain(IsolationError):
    """A TID identity changed after a numeric-TID affinity syscall."""


class UnsafeNewcomer(IsolationError):
    def __init__(self, identity, reason):
        super().__init__(reason)
        self.identity = identity


def require(condition, message):
    if not condition:
        raise IsolationError(message)


def check(condition, message):
    if not condition:
        raise SelfTestError(message)


def expect_exception(exception_type, function, label):
    try:
        function()
    except exception_type:
        return
    except BaseException as error:
        raise SelfTestError(
            "{} raised {}, expected {}".format(
                label, type(error).__name__, exception_type.__name__)) from error
    raise SelfTestError("{} did not raise {}".format(label, exception_type.__name__))


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="microseconds").replace("+00:00", "Z")


def validate_utc(value, label, allow_none=False):
    if value is None and allow_none:
        return None
    require(isinstance(value, str) and value.endswith("Z"),
            "{} is not canonical UTC".format(label))
    try:
        parsed = datetime.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise IsolationError("{} is not canonical UTC".format(label)) from error
    require(parsed.tzinfo is not None and parsed.utcoffset() == datetime.timedelta(0),
            "{} is not UTC".format(label))
    require(parsed.isoformat(timespec="microseconds").replace("+00:00", "Z") == value,
            "{} is not canonical UTC".format(label))
    return value


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("utf-8")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_value(value):
    return sha256_bytes(canonical_bytes(value))


def exact_int(value, label, minimum=0, maximum=1048575):
    if type(value) is not int or value < minimum or value > maximum:
        raise IsolationError("{} is not a bounded integer".format(label))
    return value


def exact_bool(value, label):
    require(type(value) is bool, "{} is not a boolean".format(label))
    return value


def exact_nullable_string(value, label):
    require(value is None or isinstance(value, str),
            "{} is not null or a string".format(label))
    return value


def require_exact_keys(value, keys, label):
    require(isinstance(value, dict) and set(value) == set(keys),
            "{} has unexpected or missing fields".format(label))
    return value


def validate_cpu_list(values, label, allow_empty=False):
    require(isinstance(values, list), "{} is not a list".format(label))
    result = [exact_int(value, "{} CPU".format(label)) for value in values]
    require(result == sorted(set(result)),
            "{} is not sorted and unique".format(label))
    require(allow_empty or bool(result), "{} is empty".format(label))
    return result


def parse_cpu_list(text):
    result = set()
    require(isinstance(text, str) and bool(text), "CPU list is empty")
    for part in text.split(","):
        fields = part.split("-", 1)
        try:
            first = int(fields[0], 10)
            last = int(fields[1], 10) if len(fields) == 2 else first
        except ValueError as error:
            raise IsolationError("invalid CPU list: {}".format(text)) from error
        require(0 <= first <= last <= 1048575,
                "invalid CPU range: {}".format(part))
        result.update(range(first, last + 1))
    return sorted(result)


def format_cpu_list(values):
    values = validate_cpu_list(sorted(values), "formatted CPU list")
    ranges = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(first) if first == previous else
                      "{}-{}".format(first, previous))
        first = previous = value
    ranges.append(str(first) if first == previous else
                  "{}-{}".format(first, previous))
    return ",".join(ranges)


def open_bounded_file_identity(argument):
    """Open and hash one exact regular-file inode for later descriptor exec."""
    require(isinstance(argument, str) and argument and "\x00" not in argument,
            "file argument is invalid")
    path = Path(argument).resolve(strict=True)
    descriptor = os.open(
        str(path), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NONBLOCK", 0))
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode) and 0 <= before.st_size <= MAX_IDENTITY_BYTES,
                "identity file is not a bounded regular file: {}".format(path))
        digest = hashlib.sha256()
        retained = 0
        while True:
            block = os.read(descriptor, min(1024 * 1024,
                                             MAX_IDENTITY_BYTES + 1 - retained))
            if not block:
                break
            retained += len(block)
            require(retained <= MAX_IDENTITY_BYTES,
                    "identity file is too large: {}".format(path))
            digest.update(block)
        after = os.fstat(descriptor)
        current = os.stat(path)
        require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
                 before.st_ctime_ns) ==
                (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                 after.st_ctime_ns) and
                (before.st_dev, before.st_ino) == (current.st_dev, current.st_ino),
                "identity file changed while hashing: {}".format(path))
        identity = {
            "argument": argument,
            "path": str(path),
            "device": before.st_dev,
            "inode": before.st_ino,
            "mode": before.st_mode & 0o7777,
            "size": before.st_size,
            "sha256": digest.hexdigest(),
        }
        return descriptor, identity
    except BaseException:
        os.close(descriptor)
        raise


def bounded_file_identity(argument):
    descriptor, identity = open_bounded_file_identity(argument)
    os.close(descriptor)
    return identity


class PinnedCommand:
    """Descriptor-bound command/cwd snapshot used by the gated child.

    The report keeps the user's original argv and path identities.  Execution
    uses ``/proc/self/fd`` names for the already-hashed executable and optional
    Python script, so replacing either pathname after preparation cannot select
    different code.  The child also fchdir()s to the captured directory inode.
    """

    def __init__(self, command):
        require(isinstance(command, list) and command and
                all(isinstance(item, str) and item and "\x00" not in item
                    for item in command), "child command is malformed")
        self.command = list(command)
        self.descriptors = []
        self.cwd_descriptor = None
        self.executable_descriptor = None
        self.script_descriptor = None
        try:
            cwd = Path.cwd().resolve(strict=True)
            cwd_flags = (os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                         getattr(os, "O_CLOEXEC", 0))
            if hasattr(os, "O_NOFOLLOW"):
                cwd_flags |= os.O_NOFOLLOW
            self.cwd_descriptor = os.open(str(cwd), cwd_flags)
            self.descriptors.append(self.cwd_descriptor)
            cwd_metadata = os.fstat(self.cwd_descriptor)
            current_cwd = cwd.stat()
            require(stat.S_ISDIR(cwd_metadata.st_mode) and
                    (cwd_metadata.st_dev, cwd_metadata.st_ino) ==
                    (current_cwd.st_dev, current_cwd.st_ino),
                    "working directory changed while it was captured")
            cwd_identity = {
                "path": str(cwd), "device": cwd_metadata.st_dev,
                "inode": cwd_metadata.st_ino,
                "mode": cwd_metadata.st_mode & 0o7777,
            }

            executable_argument = self.command[0]
            executable_path = executable_argument
            if os.sep not in executable_argument:
                executable_path = shutil.which(executable_argument)
                require(executable_path is not None,
                        "child executable is not present on PATH")
            descriptor, executable = open_bounded_file_identity(executable_path)
            self.executable_descriptor = descriptor
            self.descriptors.append(descriptor)
            executable["argument"] = executable_argument
            require(executable["mode"] & 0o111,
                    "child executable has no execute permission bit")
            os.set_inheritable(descriptor, True)

            script = None
            if len(self.command) >= 2 and self.command[1].endswith(".py"):
                descriptor, script = open_bounded_file_identity(self.command[1])
                self.script_descriptor = descriptor
                self.descriptors.append(descriptor)
                os.set_inheritable(descriptor, True)
            self.identity = {
                "cwd": cwd_identity,
                "executable": executable,
                "script": script,
            }
        except BaseException:
            self.close()
            raise

    @staticmethod
    def _proc_fd(descriptor):
        return "/proc/self/fd/{}".format(descriptor)

    def executable_path(self):
        require(self.executable_descriptor is not None,
                "pinned executable is closed")
        return self._proc_fd(self.executable_descriptor)

    def execution_argv(self):
        result = list(self.command)
        if self.script_descriptor is not None:
            result[1] = self._proc_fd(self.script_descriptor)
        return result

    def close(self):
        for descriptor in reversed(getattr(self, "descriptors", [])):
            try:
                os.close(descriptor)
            except OSError:
                pass
        self.descriptors = []
        self.cwd_descriptor = None
        self.executable_descriptor = None
        self.script_descriptor = None

    def __del__(self):  # pragma: no cover - deterministic paths call close()
        self.close()
def runtime_identity(proc_root=Path("/proc")):
    proc_root = Path(proc_root)
    try:
        boot_id = (proc_root / "sys/kernel/random/boot_id").read_text(
            encoding="ascii").strip()
        namespace = (proc_root / "self/ns/pid").stat()
    except OSError as error:
        raise IsolationError("cannot identify boot and PID namespace") from error
    require(BOOT_ID.fullmatch(boot_id) is not None, "kernel boot ID is malformed")
    return {
        "boot_id": boot_id,
        "pid_namespace": {"device": namespace.st_dev, "inode": namespace.st_ino},
    }


def fsync_directory(path):
    descriptor = os.open(str(path), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


class GlobalSupervisorLock:
    """Serialize every same-UID run and recovery transaction."""

    def __init__(self, uid=None, runtime_directory=None):
        self.uid = os.getuid() if uid is None else exact_int(
            uid, "global lock UID", 0, 2**31 - 1)
        self.directory = Path(
            "/run/user/{}".format(self.uid) if runtime_directory is None
            else runtime_directory)
        self.path = self.directory / "leopard2-affinity-supervisor.lock"
        self.descriptor = None
        self.identity = None

    def __enter__(self):
        directory = self.directory.resolve(strict=True)
        metadata = directory.stat()
        require(stat.S_ISDIR(metadata.st_mode) and metadata.st_uid == self.uid,
                "global affinity-lock directory has unsafe ownership/type")
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0)
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(str(self.path), flags, 0o600)
        try:
            os.fchmod(descriptor, 0o600)
            file_metadata = os.fstat(descriptor)
            require(stat.S_ISREG(file_metadata.st_mode) and
                    file_metadata.st_uid == self.uid and file_metadata.st_nlink == 1,
                    "global affinity lock has unsafe ownership/type/link count")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise IsolationError(
                    "another same-UID affinity supervisor is active") from error
            current = os.stat(self.path)
            require((current.st_dev, current.st_ino) ==
                    (file_metadata.st_dev, file_metadata.st_ino),
                    "global affinity lock path changed while acquiring")
            os.fsync(descriptor)
            fsync_directory(directory)
            self.descriptor = descriptor
            self.identity = {
                "path": str(self.path.resolve()), "uid": self.uid,
                "device": file_metadata.st_dev, "inode": file_metadata.st_ino,
                "mode": file_metadata.st_mode & 0o7777,
                "directory_device": metadata.st_dev,
                "directory_inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_same_uid",
            }
            return self
        except BaseException:
            os.close(descriptor)
            raise

    def validate_current(self):
        require(self.descriptor is not None and self.identity is not None,
                "global affinity lock is not held")
        metadata = os.fstat(self.descriptor)
        current = os.stat(self.path)
        require((metadata.st_dev, metadata.st_ino) ==
                (self.identity["device"], self.identity["inode"]) ==
                (current.st_dev, current.st_ino),
                "global affinity lock identity changed")
        try:
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise IsolationError("global affinity lock was lost") from error
        return dict(self.identity)

    def __exit__(self, _exc_type, _exc, _traceback):
        descriptor = self.descriptor
        self.descriptor = None
        if descriptor is not None:
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)


def atomic_json(path, value, exclusive=False):
    """Write canonical JSON and durably bind the directory entry.

    The parent directory must already exist.  A directory-fsync failure is a
    transaction failure; it is never swallowed before affinity mutation.
    """
    path = Path(path)
    require(path.parent.is_dir(),
            "report parent directory does not exist: {}".format(path.parent))
    if exclusive and path.exists():
        raise IsolationError("report path already exists: {}".format(path))
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    installed = False
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(canonical_bytes(value))
            output.flush()
            os.fsync(output.fileno())
        if exclusive:
            try:
                os.link(temporary, str(path))
            except FileExistsError as error:
                raise IsolationError(
                    "report path already exists: {}".format(path)) from error
            os.unlink(temporary)
        else:
            os.replace(temporary, str(path))
        installed = True
        fsync_directory(path.parent)
    except BaseException as error:
        if installed and exclusive:
            # An exclusive evidence file is not usable unless its directory
            # entry was durably committed.  Remove a possibly-installed file
            # before returning the durability failure.
            try:
                os.unlink(path)
                fsync_directory(path.parent)
            except FileNotFoundError:
                pass
            except BaseException as cleanup_error:
                raise IsolationError(
                    "exclusive write failed and installed-file cleanup was "
                    "not durable: {}; {}".format(error, cleanup_error)) from error
        raise
    finally:
        if not installed:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def acceptance_path(report_path):
    report_path = Path(report_path)
    return report_path.with_name(report_path.name + ".accepted.json")


def ambiguity_path(report_path):
    report_path = Path(report_path)
    return report_path.with_name(report_path.name + ".ambiguous")


def durable_unlink(path):
    path = Path(path)
    try:
        os.unlink(path)
    except FileNotFoundError:
        return False
    fsync_directory(path.parent)
    return True


def stable_file_snapshot(path, label, maximum=MAX_IDENTITY_BYTES,
                         after_read=None):
    """Read one inode-bound regular-file snapshot without validate/reopen races."""
    path = Path(os.path.abspath(os.fspath(path)))
    try:
        path_before = os.lstat(path)
    except OSError as error:
        raise IsolationError("cannot inspect {}: {}".format(label, error)) from error
    require(stat.S_ISREG(path_before.st_mode) and path_before.st_nlink == 1,
            "{} path is not one singly-linked regular file".format(label))
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NONBLOCK", 0))
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = None
    try:
        descriptor = os.open(path, flags)
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                (before.st_dev, before.st_ino) ==
                (path_before.st_dev, path_before.st_ino),
                "{} is not one singly-linked regular file".format(label))
        blocks = []
        size = 0
        while True:
            block = os.read(descriptor, min(65536, maximum + 1 - size))
            if not block:
                break
            blocks.append(block)
            size += len(block)
            require(size <= maximum, "{} exceeds its size bound".format(label))
        if after_read is not None:
            after_read(path)
        after = os.fstat(descriptor)
        current = os.lstat(path)
        identity = (before.st_dev, before.st_ino, before.st_size,
                    before.st_mtime_ns, before.st_ctime_ns)
        require(identity == (after.st_dev, after.st_ino, after.st_size,
                             after.st_mtime_ns, after.st_ctime_ns) and
                (current.st_dev, current.st_ino) ==
                (before.st_dev, before.st_ino) and
                size == before.st_size,
                "{} changed while it was read".format(label))
        data = b"".join(blocks)
        return path, data, {
            "device": before.st_dev, "inode": before.st_ino,
            "size": size, "sha256": sha256_bytes(data),
        }
    except OSError as error:
        raise IsolationError("cannot snapshot {}: {}".format(label, error)) from error
    finally:
        if descriptor is not None:
            os.close(descriptor)


def stable_json_snapshot(path, label, maximum=MAX_IDENTITY_BYTES):
    path, data, identity = stable_file_snapshot(path, label, maximum)
    try:
        return path, json.loads(data.decode("utf-8")), data, identity
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise IsolationError("{} is not canonical JSON".format(label)) from error


def validate_acceptance(value, report_path, report_bytes):
    value = verify_signed_json(value, "affinity acceptance seal")
    require_exact_keys(value, ACCEPTANCE_KEYS, "affinity acceptance seal")
    require(value["schema"] == ACCEPTANCE_SCHEMA,
            "affinity acceptance seal has wrong schema")
    validate_utc(value["committed_utc"], "affinity acceptance commit")
    exact_int(value["report_size"], "acceptance report size", 0,
              MAX_IDENTITY_BYTES)
    require(isinstance(value["report_sha256"], str) and
            HEX256.fullmatch(value["report_sha256"]) is not None,
            "acceptance report hash is invalid")
    canonical_report = str(Path(report_path).resolve(strict=True))
    require(value["report_path"] == canonical_report and
            value["report_size"] == len(report_bytes) and
            value["report_sha256"] == sha256_bytes(report_bytes),
            "affinity acceptance seal does not bind the report snapshot")
    return value


def accepted_report_snapshot(path):
    """Return one sealed report snapshot or fail closed on ambiguity."""
    path = Path(os.path.abspath(os.fspath(path)))
    require(not ambiguity_path(path).exists(),
            "affinity report has a durable ambiguity marker")
    try:
        path, report, report_bytes, _report_identity = stable_json_snapshot(
            path, "affinity report")
        _seal_path, seal, seal_bytes, _seal_identity = stable_json_snapshot(
            acceptance_path(path), "affinity acceptance seal")
    except IsolationError as error:
        raise IsolationError("accepted affinity evidence is incomplete") from error
    require(len(report_bytes) <= MAX_IDENTITY_BYTES and
            len(seal_bytes) <= MAX_IDENTITY_BYTES,
            "accepted affinity evidence exceeds its size bound")
    try:
        report = validate_report(report)
        seal = validate_acceptance(seal, path, report_bytes)
    except IsolationError:
        raise
    require(report["accepted"] is True,
            "affinity report body is not accepted")
    return report, report_bytes, seal, seal_bytes


def invalidate_acceptance_path(report_path, replace_function=os.replace,
                               fsync_function=fsync_directory,
                               unlink_function=durable_unlink):
    """Atomically exchange the accepting name for a rejecting marker."""
    report_path = Path(report_path)
    seal = acceptance_path(report_path)
    marker = ambiguity_path(report_path)
    if not seal.exists():
        return False
    try:
        replace_function(seal, marker)
        fsync_function(report_path.parent)
    except BaseException as error:
        try:
            unlink_function(seal)
        except BaseException as unlink_error:
            raise IsolationError(
                "acceptance rename and unlink both failed: {}; {}".format(
                    error, unlink_error)) from error
        raise IsolationError(
            "acceptance invalidation was visible but not durably confirmed: "
            "{}".format(error)) from error
    return True


def parse_proc_stat(text, label):
    right = text.rfind(")")
    require(right >= 0 and right + 2 <= len(text),
            "{} has malformed stat data".format(label))
    fields = text[right + 2:].split()
    require(len(fields) > 19, "{} has truncated stat data".format(label))
    try:
        # fields[0] is proc field 3; ppid/session/starttime are 4/6/22.
        return int(fields[1], 10), int(fields[3], 10), int(fields[19], 10)
    except ValueError as error:
        raise IsolationError("{} has nonnumeric stat data".format(label)) from error


@dataclass(frozen=True)
class ThreadIdentity:
    pid: int
    tid: int
    starttime_ticks: int
    uid: int
    session: int
    ppid: int
    process_starttime_ticks: int
    task_device: int
    task_inode: int
    process_device: int
    process_inode: int

    def key(self):
        return (self.pid, self.tid, self.starttime_ticks,
                self.task_device, self.task_inode)

    def process_key(self):
        return (self.pid, self.process_starttime_ticks,
                self.process_device, self.process_inode)


class LinuxThreadBackend:
    def __init__(self, proc_root=Path("/proc"), uid=None):
        self.proc_root = Path(proc_root)
        self.uid = os.getuid() if uid is None else exact_int(
            uid, "UID", 0, 2**31 - 1)

    @staticmethod
    def validate_capabilities():
        require(hasattr(os, "pidfd_open") and hasattr(signal, "pidfd_send_signal"),
                "safe child/newcomer cleanup requires pidfd_open and "
                "pidfd_send_signal")
        require(hasattr(signal, "pthread_sigmask") and hasattr(signal, "sigpending"),
                "signal-race closure requires pthread_sigmask and sigpending")

    @staticmethod
    def enable_subreaper():
        # Linux prctl(2): orphaned benchmark grandchildren are adopted by this
        # supervisor instead of disappearing from the ancestry proof.
        pr_set_child_subreaper = 36
        pr_get_child_subreaper = 37
        libc = ctypes.CDLL(None, use_errno=True)
        if libc.prctl(pr_set_child_subreaper, 1, 0, 0, 0) != 0:
            error = ctypes.get_errno()
            raise IsolationError("cannot enable child subreaper: {}".format(
                os.strerror(error)))
        enabled = ctypes.c_int(0)
        if libc.prctl(pr_get_child_subreaper, ctypes.byref(enabled), 0, 0, 0) != 0 \
                or enabled.value != 1:
            raise IsolationError("kernel did not retain child-subreaper state")
        return os.getpid()

    def _read_stat(self, path, label):
        try:
            return parse_proc_stat(path.read_text(encoding="ascii"), label)
        except (FileNotFoundError, ProcessLookupError) as error:
            raise ThreadGone() from error
        except PermissionError as error:
            raise IsolationError("cannot inspect {}".format(label)) from error

    @staticmethod
    def _open_proc_directory(path, label):
        flags = (getattr(os, "O_PATH", os.O_RDONLY) |
                 getattr(os, "O_DIRECTORY", 0) |
                 getattr(os, "O_CLOEXEC", 0))
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            descriptor = os.open(str(path), flags)
            return descriptor, os.fstat(descriptor)
        except (FileNotFoundError, ProcessLookupError) as error:
            raise ThreadGone() from error
        except PermissionError as error:
            raise IsolationError("cannot inspect {}".format(label)) from error

    def _identity(self, pid, tid):
        task = self.proc_root / str(pid) / "task" / str(tid)
        process = self.proc_root / str(pid)
        process_fd = None
        task_fd = None
        try:
            process_fd, process_metadata = self._open_proc_directory(
                process, "process {}".format(pid))
            task_fd, task_metadata = self._open_proc_directory(
                task, "thread {}/{}".format(pid, tid))
            owner = task_metadata.st_uid
            if owner != self.uid:
                return None
            require(process_metadata.st_uid == owner,
                    "thread/process ownership differs for {}/{}".format(pid, tid))
            ppid, session_id, starttime = self._read_stat(
                task / "stat", "thread {}/{}".format(pid, tid))
            process_ppid, process_session, process_starttime = self._read_stat(
                process / "stat", "process {}".format(pid))
            require(ppid == process_ppid and session_id == process_session,
                    "thread/process identity changed while inspecting {}/{}".format(
                        pid, tid))
            current_task = task.stat()
            current_process = process.stat()
            require((current_task.st_dev, current_task.st_ino) ==
                    (task_metadata.st_dev, task_metadata.st_ino) and
                    (current_process.st_dev, current_process.st_ino) ==
                    (process_metadata.st_dev, process_metadata.st_ino),
                    "procfs identity changed while inspecting {}/{}".format(pid, tid))
            return ThreadIdentity(
                pid, tid, starttime, owner, session_id, ppid, process_starttime,
                task_metadata.st_dev, task_metadata.st_ino,
                process_metadata.st_dev, process_metadata.st_ino)
        except (FileNotFoundError, ProcessLookupError) as error:
            raise ThreadGone() from error
        finally:
            if task_fd is not None:
                os.close(task_fd)
            if process_fd is not None:
                os.close(process_fd)

    def _open_identity_guard(self, identity):
        process_path = self.proc_root / str(identity.pid)
        task_path = process_path / "task" / str(identity.tid)
        process_fd = None
        task_fd = None
        try:
            process_fd, process_metadata = self._open_proc_directory(
                process_path, "process {}".format(identity.pid))
            task_fd, task_metadata = self._open_proc_directory(
                task_path, "thread {}/{}".format(identity.pid, identity.tid))
            if ((process_metadata.st_dev, process_metadata.st_ino) !=
                    (identity.process_device, identity.process_inode) or
                    (task_metadata.st_dev, task_metadata.st_ino) !=
                    (identity.task_device, identity.task_inode)):
                raise ThreadReused()
            self.current_identity(identity)
            return process_fd, task_fd
        except BaseException:
            if task_fd is not None:
                os.close(task_fd)
            if process_fd is not None:
                os.close(process_fd)
            raise

    def list_threads(self, excluded_session=None):
        result = []
        try:
            processes = sorted(
                (entry for entry in self.proc_root.iterdir() if entry.name.isdigit()),
                key=lambda entry: int(entry.name))
        except OSError as error:
            raise IsolationError("cannot enumerate procfs") from error
        for process in processes:
            try:
                if process.stat().st_uid != self.uid:
                    continue
                tasks = sorted(
                    (entry for entry in (process / "task").iterdir()
                     if entry.name.isdigit()), key=lambda entry: int(entry.name))
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError as error:
                raise IsolationError(
                    "cannot enumerate same-UID process {}".format(process.name)) from error
            for task in tasks:
                try:
                    identity = self._identity(int(process.name), int(task.name))
                except (ThreadGone, ThreadReused):
                    continue
                if identity is not None and identity.session != excluded_session:
                    result.append(identity)
        return result

    def process_identity(self, pid):
        return self._identity(pid, pid)

    def current_identity(self, identity):
        current = self._identity(identity.pid, identity.tid)
        if current is None:
            raise ThreadReused()
        if (current.starttime_ticks != identity.starttime_ticks or
                current.process_starttime_ticks != identity.process_starttime_ticks or
                current.task_device != identity.task_device or
                current.task_inode != identity.task_inode or
                current.process_device != identity.process_device or
                current.process_inode != identity.process_inode):
            raise ThreadReused()
        return current

    def get_affinity(self, identity):
        process_fd, task_fd = self._open_identity_guard(identity)
        try:
            try:
                affinity = set(os.sched_getaffinity(identity.tid))
            except ProcessLookupError as error:
                raise ThreadGone() from error
            except PermissionError as error:
                raise IsolationError(
                    "cannot read affinity for same-UID thread {}/{}".format(
                        identity.pid, identity.tid)) from error
            try:
                self.current_identity(identity)
            except ThreadReused as error:
                raise ThreadMutationUncertain(
                    "TID identity changed across sched_getaffinity for {}/{}".format(
                        identity.pid, identity.tid)) from error
            return affinity
        finally:
            os.close(task_fd)
            os.close(process_fd)

    def set_affinity(self, identity, cpus):
        process_fd, task_fd = self._open_identity_guard(identity)
        try:
            try:
                os.sched_setaffinity(identity.tid, set(cpus))
            except ProcessLookupError as error:
                raise ThreadGone() from error
            except (PermissionError, OSError) as error:
                if isinstance(error, OSError) and error.errno == errno.ESRCH:
                    raise ThreadGone() from error
                raise IsolationError(
                    "cannot set affinity for same-UID thread {}/{}: {}".format(
                        identity.pid, identity.tid, error)) from error
            try:
                self.current_identity(identity)
            except ThreadReused as error:
                raise ThreadMutationUncertain(
                    "TID identity changed after sched_setaffinity for {}/{}".format(
                        identity.pid, identity.tid)) from error
        finally:
            os.close(task_fd)
            os.close(process_fd)
        actual = self.get_affinity(identity)
        if actual != set(cpus):
            raise IsolationError(
                "kernel did not retain exact affinity for thread {}/{}".format(
                    identity.pid, identity.tid))

    def session_processes(self, session_id):
        by_process = {}
        for identity in self.list_threads(excluded_session=None):
            if identity.session == session_id:
                by_process.setdefault(identity.process_key(), identity)
        result = []
        for identity in by_process.values():
            try:
                leader = self.process_identity(identity.pid)
            except ThreadGone:
                continue
            if leader is not None and leader.session == session_id:
                result.append(leader)
        return sorted(result, key=lambda item: (item.pid, item.starttime_ticks))

    def child_scope_processes(self, child_record, retained=(), subreaper_pid=None):
        identities = self.list_threads(excluded_session=None)
        keys = child_scope_process_keys(
            identities, child_record, retained=retained,
            subreaper_pid=subreaper_pid)
        by_process = {}
        for identity in identities:
            if identity.process_key() in keys:
                by_process.setdefault(identity.process_key(), identity)
        result = []
        for identity in by_process.values():
            try:
                leader = self.process_identity(identity.pid)
            except ThreadGone:
                continue
            if leader is not None and leader.process_key() == identity.process_key():
                result.append(leader)
        return sorted(result, key=lambda item: (item.pid, item.starttime_ticks))

    def signal_process(self, identity, signum):
        current = self.process_identity(identity.pid)
        if current is None or current.process_key() != identity.process_key():
            raise ThreadReused()
        self.validate_capabilities()
        pidfd = None
        try:
            pidfd = os.pidfd_open(identity.pid, 0)
            current = self.process_identity(identity.pid)
            if current.process_key() != identity.process_key():
                raise ThreadReused()
            signal.pidfd_send_signal(pidfd, signum)
        except ProcessLookupError as error:
            raise ThreadGone() from error
        finally:
            if pidfd is not None:
                os.close(pidfd)

    def reap_adopted(self, identities, excluded_pid=None):
        """Reap exact retained orphans without stealing the leader's status."""
        reaped = 0
        for identity in identities:
            if identity.pid == excluded_pid:
                continue
            try:
                current = self.process_identity(identity.pid)
            except ThreadGone:
                continue
            if current is None or current.process_key() != identity.process_key():
                continue
            try:
                waited, _status = os.waitpid(identity.pid, os.WNOHANG)
            except (ChildProcessError, ProcessLookupError):
                continue
            if waited == identity.pid:
                reaped += 1
        return reaped


def validate_path_identity(value, label):
    require_exact_keys(value, PATH_IDENTITY_KEYS, label)
    require(isinstance(value["path"], str) and value["path"],
            "{} path is invalid".format(label))
    for name in ("device", "inode", "mode"):
        exact_int(value[name], "{} {}".format(label, name), 0, 2**63 - 1)


def validate_file_identity(value, label):
    require_exact_keys(value, FILE_IDENTITY_KEYS, label)
    require(all(isinstance(value[name], str) and value[name]
                for name in ("argument", "path")),
            "{} paths are invalid".format(label))
    for name in ("device", "inode", "mode", "size"):
        exact_int(value[name], "{} {}".format(label, name), 0, 2**63 - 1)
    require(isinstance(value["sha256"], str) and
            HEX256.fullmatch(value["sha256"]) is not None,
            "{} hash is invalid".format(label))


def validate_global_lock(value, uid, label="global affinity lock"):
    require_exact_keys(value, GLOBAL_LOCK_KEYS, label)
    exact_int(value["uid"], "{} UID".format(label), 0, 2**31 - 1)
    require(value["uid"] == uid, "{} UID differs from report UID".format(label))
    require(isinstance(value["path"], str) and value["path"].startswith("/"),
            "{} path is not absolute".format(label))
    for name in ("device", "inode", "mode", "directory_device",
                 "directory_inode"):
        exact_int(value[name], "{} {}".format(label, name), 0, 2**63 - 1)
    require(value["lock"] == "exclusive_nonblocking_same_uid",
            "{} protocol is unknown".format(label))
    require(value["mode"] & 0o777 == 0o600,
            "{} mode is not owner-only".format(label))
    return value


def validate_record(record, label="affinity record"):
    require_exact_keys(record, RECORD_KEYS, label)
    for name in ("pid", "tid"):
        exact_int(record[name], "{} {}".format(label, name), 1, 2**31 - 1)
    for name in ("starttime_ticks", "process_starttime_ticks", "task_device",
                 "task_inode", "process_device", "process_inode"):
        exact_int(record[name], "{} {}".format(label, name), 0, 2**63 - 1)
    exact_int(record["uid"], "{} UID".format(label), 0, 2**31 - 1)
    exact_int(record["session"], "{} session".format(label), 0, 2**31 - 1)
    exact_int(record["ppid"], "{} parent PID".format(label), 0, 2**31 - 1)
    exact_int(record["move_count"], "{} move count".format(label), 0, 2**31 - 1)
    validate_cpu_list(record["original_cpus"], "{} original CPUs".format(label))
    validate_cpu_list(record["restricted_cpus"], "{} restricted CPUs".format(label),
                      allow_empty=False)
    require(record["move_status"] in MOVE_STATUSES,
            "{} move status is unknown".format(label))
    require(record["restore_status"] in RESTORE_STATUSES,
            "{} restore status is unknown".format(label))
    require(record["provenance"] in RECORD_PROVENANCE,
            "{} provenance is unknown".format(label))
    return record


def validate_provenance(value, label="process provenance"):
    require_exact_keys(value, PROVENANCE_KEYS, label)
    exact_int(value["pid"], "{} PID".format(label), 1, 2**31 - 1)
    exact_int(value["ppid"], "{} parent PID".format(label), 0, 2**31 - 1)
    for name in ("process_starttime_ticks", "process_device", "process_inode"):
        exact_int(value[name], "{} {}".format(label, name), 0, 2**63 - 1)
    validate_cpu_list(value["original_cpus"], "{} CPUs".format(label))
    require(value["source"] in PROVENANCE_SOURCES,
            "{} source is unknown".format(label))
    return value


def validate_child_process(value, label="child process"):
    require_exact_keys(value, CHILD_PROCESS_KEYS, label)
    for name in ("pid", "session"):
        exact_int(value[name], "{} {}".format(label, name), 1, 2**31 - 1)
    exact_int(value["ppid"], "{} parent PID".format(label), 0, 2**31 - 1)
    for name in ("starttime_ticks", "process_device", "process_inode"):
        exact_int(value[name], "{} {}".format(label, name), 0, 2**63 - 1)
    return value


def validate_execution(value):
    require_exact_keys(value, EXECUTION_KEYS, "supervised execution")
    require(isinstance(value["nonce"], str) and
            HEX256.fullmatch(value["nonce"]) is not None,
            "supervised execution nonce is invalid")
    ready = value["child_ready_monotonic_ns"]
    finished = value["child_finished_monotonic_ns"]
    require(ready is None or (type(ready) is int and 0 <= ready <= 2**63 - 1),
            "child-ready monotonic timestamp is invalid")
    require(finished is None or
            (type(finished) is int and 0 <= finished <= 2**63 - 1),
            "child-finished monotonic timestamp is invalid")
    require(finished is None or (ready is not None and finished >= ready),
            "child-finished timestamp predates child-ready")
    return value


def validate_event(value, label="event"):
    require_exact_keys(value, EVENT_KEYS, label)
    require(value["event"] in EVENT_NAMES, "{} name is unknown".format(label))
    for name in ("pid", "tid"):
        require(value[name] is None or
                (type(value[name]) is int and 0 <= value[name] <= 2**31 - 1),
                "{} {} is invalid".format(label, name))
    require(value["value"] is None or
            (isinstance(value["value"], str) and value["value"]) or
            (type(value["value"]) is int and 0 <= value["value"] <= 2**31 - 1),
            "{} value is invalid".format(label))
    return value


def uncontained_mutation_records(records, subreaper_pid):
    """Return mutations whose creator is not the supervisor's sole main thread."""
    records = list(records)
    self_records = [
        item for item in records
        if subreaper_pid is not None and item["pid"] == subreaper_pid and
        item["tid"] == subreaper_pid
    ]
    if len(self_records) > 1:
        return records
    return [item for item in records if item not in self_records]


def validate_report(report, require_accepted=False):
    require_exact_keys(report, REPORT_KEYS, "affinity report")
    require(report["schema"] == REPORT_SCHEMA, "affinity report has wrong schema")
    require(report["state"] in REPORT_STATES, "affinity report state is unknown")
    exact_int(report["uid"], "report UID", 0, 2**31 - 1)
    validate_global_lock(report["global_lock"], report["uid"])
    execution = validate_execution(report["execution"])
    reserved = validate_cpu_list(report["reserved_cpus"], "report reserved CPUs")
    launch = validate_cpu_list(report["launch_cpus"], "report launch CPUs")
    require(set(reserved).issubset(launch) and set(launch).difference(reserved),
            "report CPU sets are inconsistent")
    exact_int(report["poll_ms"], "report poll milliseconds", 1, 10000)
    validate_utc(report["started_utc"], "report start")
    validate_utc(report["finished_utc"], "report finish", allow_none=True)
    if report["finished_utc"] is not None:
        require(report["finished_utc"] >= report["started_utc"],
                "report finish predates its start")
    require(isinstance(report["boot_id"], str) and
            BOOT_ID.fullmatch(report["boot_id"]) is not None,
            "report boot ID is invalid")
    namespace = require_exact_keys(
        report["pid_namespace"], PID_NAMESPACE_KEYS, "report PID namespace")
    for name in PID_NAMESPACE_KEYS:
        exact_int(namespace[name], "PID namespace {}".format(name), 0, 2**63 - 1)
    require(isinstance(report["command"], list) and
            all(isinstance(item, str) and item and "\x00" not in item
                for item in report["command"]), "report command is malformed")
    if report["command"]:
        require(isinstance(report["command_sha256"], str) and
                HEX256.fullmatch(report["command_sha256"]) is not None and
                report["command_sha256"] == sha256_value(report["command"]),
                "report command hash is invalid")
        identity = require_exact_keys(
            report["command_identity"], COMMAND_IDENTITY_KEYS,
            "report command identity")
        validate_path_identity(identity["cwd"], "command working directory")
        validate_file_identity(identity["executable"], "command executable")
        require(identity["executable"]["argument"] == report["command"][0],
                "command executable identity has another argument")
        if identity["script"] is not None:
            validate_file_identity(identity["script"], "command script")
            require(len(report["command"]) >= 2 and
                    identity["script"]["argument"] == report["command"][1],
                    "command script identity has another argument")
    else:
        require(report["command_sha256"] is None and
                report["command_identity"] is None,
                "empty report command has an identity")
    child = report["child"]
    if child is not None:
        require_exact_keys(child, CHILD_KEYS, "report child")
        for name in ("pid", "session"):
            exact_int(child[name], "child {}".format(name), 1, 2**31 - 1)
        for name in ("starttime_ticks", "process_device", "process_inode"):
            exact_int(child[name], "child {}".format(name), 0, 2**63 - 1)
        exact_bool(child["released"], "child release flag")
        require(child["pid"] == child["session"],
                "child is not its recorded session leader")
    require(isinstance(report["child_processes"], list),
            "report child-process identities are malformed")
    child_processes = [
        validate_child_process(item, "child process {}".format(index))
        for index, item in enumerate(report["child_processes"])
    ]
    child_process_keys = [child_process_key(item) for item in child_processes]
    require(len(child_process_keys) == len(set(child_process_keys)),
            "report child-process identities are duplicated")
    if child is None:
        require(not child_processes,
                "report has child-process identities without a child")
    else:
        require((child["pid"], child["starttime_ticks"],
                 child["process_device"], child["process_inode"]) in
                set(child_process_keys),
                "report child root is absent from retained process identities")
    require(report["subreaper_pid"] is None or
            (type(report["subreaper_pid"]) is int and
             1 <= report["subreaper_pid"] <= 2**31 - 1),
            "report subreaper PID is invalid")
    require(report["child_returncode"] is None or
            (type(report["child_returncode"]) is int and
             -(2**31) <= report["child_returncode"] <= 2**31 - 1),
            "child return code is invalid")
    require(report["received_signal"] is None or
            (type(report["received_signal"]) is int and
             1 <= report["received_signal"] <= 255),
            "received signal is invalid")
    require(isinstance(report["records"], list), "report records are malformed")
    records = [validate_record(item, "affinity record {}".format(index))
               for index, item in enumerate(report["records"])]
    for item in records:
        require(item["uid"] == report["uid"],
                "affinity record UID differs from report UID")
        require(set(item["restricted_cpus"]) ==
                set(item["original_cpus"]).difference(reserved),
                "affinity record restricted mask is not exact subtraction")
        require(set(item["original_cpus"]).intersection(reserved),
                "affinity record original mask never used a reserved CPU")
        if item["move_status"] == "moved":
            require(item["move_count"] >= 1,
                    "moved affinity record has no successful move")
    record_keys = [(item["pid"], item["tid"], item["starttime_ticks"],
                    item["task_device"], item["task_inode"])
                   for item in records]
    require(len(record_keys) == len(set(record_keys)), "report records are duplicated")
    require(isinstance(report["process_provenance"], list),
            "report process provenance is malformed")
    provenance = [validate_provenance(item, "process provenance {}".format(index))
                  for index, item in enumerate(report["process_provenance"])]
    process_keys = [(item["pid"], item["process_starttime_ticks"],
                     item["process_device"], item["process_inode"])
                    for item in provenance]
    require(len(process_keys) == len(set(process_keys)),
            "report process provenance is duplicated")
    provenance_by_key = {
        (item["pid"], item["process_starttime_ticks"],
         item["process_device"], item["process_inode"]): item
        for item in provenance
    }
    for item in records:
        if item["provenance"] in {"baseline", "inherited"}:
            process = provenance_by_key.get(
                (item["pid"], item["process_starttime_ticks"],
                 item["process_device"], item["process_inode"]))
            require(process is not None and
                    process["original_cpus"] == item["original_cpus"],
                    "record original mask differs from proven process mask")
    require(isinstance(report["events"], list), "report events are malformed")
    for index, event in enumerate(report["events"]):
        validate_event(event, "event {}".format(index))
    signal_events = [item for item in report["events"] if item["event"] == "signal"]
    if report["received_signal"] is None:
        require(not signal_events, "report has a signal event without a signal")
    else:
        require(len(signal_events) == 1 and
                signal_events[0]["value"] == report["received_signal"],
                "report signal event differs from received signal")
    exact_bool(report["uncertainty"], "report uncertainty")
    exact_bool(report["accepted"], "report acceptance")
    exact_nullable_string(report["error"], "report error")
    terminal_states = {"complete", "failed", "recovered", "recovery_failed"}
    if report["state"] in terminal_states:
        require(report["finished_utc"] is not None,
                "terminal report has no finish timestamp")
    if report["state"] == "complete":
        require(report["accepted"] is True,
                "complete report is not accepted")
    if report["state"] in {"failed", "recovery_failed"}:
        require(isinstance(report["error"], str) and bool(report["error"]),
                "failed report has no error")
    if report["state"] in {"recovered", "recovery_failed"}:
        require(report["accepted"] is False,
                "recovery report cannot be accepted evidence")
    if any(item["restore_status"] in {"reused", "uncertain"} for item in records):
        require(report["uncertainty"] is True,
                "TID-reuse record does not mark report uncertainty")
    if report["accepted"]:
        require(report["state"] == "complete" and report["finished_utc"] is not None and
                report["error"] is None and not report["uncertainty"] and
                report["subreaper_pid"] is not None and
                not uncontained_mutation_records(
                    records, report["subreaper_pid"]) and
                report["child_returncode"] == 0 and
                report["received_signal"] is None and child is not None and
                child["released"] is True and
                execution["child_ready_monotonic_ns"] is not None and
                execution["child_finished_monotonic_ns"] is not None and
                all(item["restore_status"] in {
                    "restored", "restored_after_recovery", "gone"}
                    for item in records),
                "accepted report is not a complete clean transaction")
    if require_accepted:
        raise IsolationError(
            "accepted evidence is path-bound; use load_report(..., "
            "require_accepted=True)")
    return report


class ForkedChild:
    def __init__(self, pid, identity):
        self.pid = pid
        self.identity = identity
        self.returncode = None

    @staticmethod
    def _decode_status(status):
        if os.WIFEXITED(status):
            return os.WEXITSTATUS(status)
        if os.WIFSIGNALED(status):
            return -os.WTERMSIG(status)
        return 1

    def poll(self):
        if self.returncode is not None:
            return self.returncode
        try:
            waited, status = os.waitpid(self.pid, os.WNOHANG)
        except ChildProcessError:
            return self.returncode
        if waited == self.pid:
            self.returncode = self._decode_status(status)
        return self.returncode

    def wait(self, timeout):
        deadline = time.monotonic() + timeout
        while self.poll() is None:
            if time.monotonic() >= deadline:
                raise TimeoutError("child did not exit before timeout")
            time.sleep(0.01)
        return self.returncode


def execute_pinned_command(pinned, launch_cpus, execution_nonce,
                           fchdir_function=os.fchdir,
                           set_affinity_function=os.sched_setaffinity,
                           get_affinity_function=os.sched_getaffinity,
                           execve_function=os.execve, environment=None):
    """Enter the captured cwd/affinity and exec only captured code in a child."""
    require(isinstance(pinned, PinnedCommand),
            "gated execution requires a pinned command")
    launch_cpus = set(validate_cpu_list(sorted(launch_cpus), "launch CPUs"))
    require(isinstance(execution_nonce, str) and
            HEX256.fullmatch(execution_nonce) is not None,
            "execution nonce is invalid")
    require(pinned.cwd_descriptor is not None,
            "pinned working directory is closed")
    fchdir_function(pinned.cwd_descriptor)
    set_affinity_function(0, launch_cpus)
    require(set(get_affinity_function(0)) == launch_cpus,
            "gated child did not retain its exact launch affinity")
    child_environment = dict(os.environ if environment is None else environment)
    child_environment[EXECUTION_NONCE_ENV] = execution_nonce
    execve_function(
        pinned.executable_path(), pinned.execution_argv(), child_environment)
    raise IsolationError("execve unexpectedly returned")


def reset_gated_child_signals(old_mask, signal_function=signal.signal,
                              mask_function=signal.pthread_sigmask):
    """Install terminating child dispositions before unblocking launch signals."""
    for signum in GATED_SIGNALS:
        signal_function(signum, signal.SIG_DFL)
    mask_function(signal.SIG_SETMASK, old_mask)


class GatedLauncher:
    """Fork a new session whose command is blocked on a parent-owned pipe."""

    def __init__(self, ready_timeout=LAUNCH_READY_TIMEOUT_SECONDS):
        self.ready_timeout = float(ready_timeout)

    @staticmethod
    def _close(descriptor):
        try:
            os.close(descriptor)
        except OSError:
            pass

    def launch(self, command, launch_cpus, backend, ready_callback,
               release_callback, execution_nonce):
        require(isinstance(command, PinnedCommand),
                "gated launcher received an unpinned command")
        gate_read, gate_write = os.pipe2(os.O_CLOEXEC)
        ready_read, ready_write = os.pipe2(os.O_CLOEXEC)
        blocked = set(GATED_SIGNALS)
        old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked)
        if set(old_mask).intersection(blocked):
            signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
            self._close(gate_read)
            self._close(gate_write)
            self._close(ready_read)
            self._close(ready_write)
            raise IsolationError(
                "authoritative child launch requires INT/TERM/HUP unblocked")
        try:
            pid = os.fork()
        except BaseException:
            signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
            self._close(gate_read)
            self._close(gate_write)
            self._close(ready_read)
            self._close(ready_write)
            raise
        if pid == 0:  # pragma: no cover - behavior is asserted by the parent
            try:
                self._close(gate_write)
                self._close(ready_read)
                # Fork inherits the supervisor's Python handlers.  Reset them
                # while blocked, then restore the original unblocked mask.  A
                # signal sent to the gated child therefore terminates it; it
                # cannot run a copied Python handler and restart gate_read.
                reset_gated_child_signals(old_mask)
                os.setsid()
                os.write(ready_write, b"R")
                self._close(ready_write)
                token = os.read(gate_read, 1)
                self._close(gate_read)
                if token != b"G":
                    os._exit(125)
                execute_pinned_command(command, launch_cpus, execution_nonce)
            except BaseException:
                os._exit(126)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        self._close(gate_read)
        self._close(ready_write)
        child = None
        try:
            poller = select.poll()
            poller.register(ready_read, select.POLLIN | select.POLLHUP | select.POLLERR)
            events = poller.poll(max(1, int(self.ready_timeout * 1000)))
            require(events and os.read(ready_read, 1) == b"R",
                    "gated child did not establish its session")
            identity = backend.process_identity(pid)
            require(identity is not None and identity.pid == pid and
                    identity.tid == pid and identity.session == pid,
                    "gated child did not become its recorded session leader")
            child = ForkedChild(pid, identity)
            # This callback must durably record PID/starttime/session.  If it
            # raises, closing gate_write gives the child EOF and it cannot exec.
            ready_callback(identity)
            # The transaction owns the last signal check and the one-byte
            # write.  Its signal handler closes this exact descriptor if a
            # handled signal wins the race before write(2) commits the byte.
            release_callback(gate_write)
            self._close(gate_write)
            gate_write = -1
            return child
        except BaseException:
            self._close(gate_write)
            if child is None:
                try:
                    identity = backend.process_identity(pid)
                    child = ForkedChild(pid, identity)
                except BaseException:
                    child = ForkedChild(pid, None)
            try:
                child.wait(self.ready_timeout)
            except TimeoutError:
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                try:
                    child.wait(self.ready_timeout)
                except TimeoutError:
                    pass
            raise
        finally:
            self._close(ready_read)
            self._close(gate_write)


def thread_from_record(record):
    return ThreadIdentity(
        record["pid"], record["tid"], record["starttime_ticks"], record["uid"],
        record["session"], record["ppid"], record["process_starttime_ticks"],
        record["task_device"], record["task_inode"],
        record["process_device"], record["process_inode"])


def child_record_from_identity(identity, released):
    return {
        "pid": identity.pid,
        "starttime_ticks": identity.process_starttime_ticks,
        "process_device": identity.process_device,
        "process_inode": identity.process_inode,
        "session": identity.session,
        "released": bool(released),
    }


def child_process_record(identity):
    return {
        "pid": identity.pid,
        "starttime_ticks": identity.process_starttime_ticks,
        "process_device": identity.process_device,
        "process_inode": identity.process_inode,
        "session": identity.session,
        "ppid": identity.ppid,
    }


def child_process_key(record):
    return (record["pid"], record["starttime_ticks"],
            record["process_device"], record["process_inode"])


def child_scope_process_keys(identities, child_record, retained=(),
                             subreaper_pid=None):
    """Return the recorded session plus descendants that created new sessions."""
    if child_record is None:
        return set()
    by_pid = {}
    for identity in identities:
        by_pid.setdefault(identity.pid, identity)
    root = by_pid.get(child_record["pid"])
    if root is not None:
        require(root.process_starttime_ticks == child_record["starttime_ticks"] and
                root.process_device == child_record["process_device"] and
                root.process_inode == child_record["process_inode"] and
                root.session == child_record["session"],
                "recorded child PID was reused while identifying descendants")
    retained_keys = {child_process_key(item) for item in retained}
    selected_pids = {
        identity.pid for identity in by_pid.values()
        if (identity.session == child_record["session"] or
            identity.process_key() in retained_keys or
            (subreaper_pid is not None and identity.ppid == subreaper_pid and
             identity.process_starttime_ticks >= child_record["starttime_ticks"]))
    }
    # Retain the recorded root PID as an ancestry anchor even after its leader
    # has been reaped; a still-running direct child continues to expose ppid.
    selected_pids.add(child_record["pid"])
    changed = True
    while changed:
        changed = False
        for identity in by_pid.values():
            if identity.pid in selected_pids:
                continue
            if identity.ppid in selected_pids:
                require(identity.process_starttime_ticks >= child_record["starttime_ticks"],
                        "descendant predates recorded child identity")
                selected_pids.add(identity.pid)
                changed = True
    return {
        identity.process_key() for identity in by_pid.values()
        if identity.pid in selected_pids
    }


def wait_for_empty_child_scope(backend, child_record, deadline,
                               sleep_function=time.sleep, reap_callback=None,
                               retained=(), subreaper_pid=None):
    empty_scans = 0
    while time.monotonic() < deadline:
        if reap_callback is not None:
            reap_callback()
        members = backend.child_scope_processes(
            child_record, retained=retained, subreaper_pid=subreaper_pid)
        if hasattr(backend, "reap_adopted") and backend.reap_adopted(
                members, excluded_pid=child_record["pid"]):
            members = backend.child_scope_processes(
                child_record, retained=retained,
                subreaper_pid=subreaper_pid)
        if not members:
            empty_scans += 1
            if empty_scans >= 2:
                return []
        else:
            empty_scans = 0
        sleep_function(min(0.02, max(0.0, deadline - time.monotonic())))
    return backend.child_scope_processes(
        child_record, retained=retained, subreaper_pid=subreaper_pid)


def terminate_session(backend, child_record, initial_signal=signal.SIGTERM,
                      grace_seconds=SIGNAL_GRACE_SECONDS,
                      kill_grace_seconds=KILL_GRACE_SECONDS,
                      sleep_function=time.sleep, reap_callback=None,
                      retained=(), subreaper_pid=None):
    """Boundedly stop the child session and traceable cross-session descendants."""
    if child_record is None:
        return [], False
    session_id = child_record["session"]
    try:
        live_leader = backend.process_identity(child_record["pid"])
    except ThreadGone:
        live_leader = None
    if live_leader is not None:
        require(live_leader.process_starttime_ticks == child_record["starttime_ticks"] and
                live_leader.process_device == child_record["process_device"] and
                live_leader.process_inode == child_record["process_inode"] and
                live_leader.session == session_id,
                "recorded child PID was reused or changed session before cleanup")
    members = backend.child_scope_processes(
        child_record, retained=retained, subreaper_pid=subreaper_pid)
    had_members = bool(members)
    members = sorted(members, key=lambda item: item.pid == child_record["pid"])
    for identity in members:
        try:
            backend.signal_process(identity, initial_signal)
        except ThreadGone:
            pass
        except ThreadReused as error:
            raise IsolationError(
                "process identity changed while signaling child session") from error
    remaining = wait_for_empty_child_scope(
        backend, child_record, time.monotonic() + grace_seconds, sleep_function,
        reap_callback, retained, subreaper_pid)
    remaining = sorted(
        remaining, key=lambda item: item.pid == child_record["pid"])
    for identity in remaining:
        try:
            backend.signal_process(identity, signal.SIGKILL)
        except ThreadGone:
            pass
        except ThreadReused as error:
            raise IsolationError(
                "process identity changed while killing child session") from error
    remaining = wait_for_empty_child_scope(
        backend, child_record, time.monotonic() + kill_grace_seconds, sleep_function,
        reap_callback, retained, subreaper_pid)
    if remaining:
        return ["child session {} retained {} process(es) after SIGKILL".format(
            session_id, len(remaining))], had_members
    return [], had_members


def terminate_process(backend, identity, grace_seconds=SIGNAL_GRACE_SECONDS,
                      kill_grace_seconds=KILL_GRACE_SECONDS,
                      sleep_function=time.sleep):
    """Fail closed for an unprovable newcomer without leaving its inherited mask."""
    for signum, timeout in ((signal.SIGTERM, grace_seconds),
                            (signal.SIGKILL, kill_grace_seconds)):
        try:
            backend.signal_process(identity, signum)
        except ThreadGone:
            return []
        except ThreadReused as error:
            return ["unsafe newcomer PID identity was reused during cleanup: {}".format(error)]
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                current = backend.process_identity(identity.pid)
            except ThreadGone:
                return []
            if current is None:
                return []
            if current.process_key() != identity.process_key():
                return ["unsafe newcomer PID was reused during cleanup"]
            sleep_function(min(0.02, max(0.0, deadline - time.monotonic())))
    return ["unsafe newcomer process {} survived SIGKILL".format(identity.pid)]


class AffinityTransaction:
    def __init__(self, backend, reserved_cpus, launch_cpus, report_path,
                 poll_ms=DEFAULT_POLL_MS, writer=atomic_json,
                 host=None, sleep_function=time.sleep, launcher=None,
                 global_lock=None, seal_writer=atomic_json,
                 nonce_factory=lambda: secrets.token_hex(32)):
        self.backend = backend
        if hasattr(self.backend, "validate_capabilities"):
            self.backend.validate_capabilities()
        self.reserved = set(validate_cpu_list(
            sorted(reserved_cpus), "reserved CPUs"))
        self.launch = set(validate_cpu_list(sorted(launch_cpus), "launch CPUs"))
        require(self.reserved.issubset(self.launch),
                "reserved CPUs are outside the launch affinity")
        require(bool(self.launch.difference(self.reserved)),
                "removing reserved CPUs leaves no housekeeping CPU")
        self.report_path = Path(report_path)
        self.poll_ms = exact_int(poll_ms, "poll milliseconds", 1, 10000)
        self.writer = writer
        self.seal_writer = seal_writer
        self.sleep = sleep_function
        self.launcher = GatedLauncher() if launcher is None else launcher
        require(global_lock is not None,
                "affinity transaction requires the same-UID global lock")
        self.global_lock = global_lock
        lock_identity = global_lock.validate_current() if hasattr(
            global_lock, "validate_current") else dict(global_lock)
        validate_global_lock(lock_identity, backend.uid)
        nonce = nonce_factory()
        require(isinstance(nonce, str) and HEX256.fullmatch(nonce) is not None,
                "execution nonce source returned an invalid value")
        identity = runtime_identity() if host is None else host
        require_exact_keys(identity, {"boot_id", "pid_namespace"}, "runtime identity")
        self.records = {}
        self.record_order = []
        self.process_provenance = {}
        self.provenance_order = []
        self.seen_keys = set()
        self.events = []
        self.child = None
        self.pinned_command = None
        self.child_processes = {}
        self.child_process_order = []
        self.subreaper_pid = None
        self.received_signal = None
        self.signal_boundary_closed = False
        self.gate_descriptor = None
        self.initialized = False
        self.write_count = 0
        require(not acceptance_path(self.report_path).exists() and
                not ambiguity_path(self.report_path).exists(),
                "report acceptance/ambiguity sidecar already exists")
        self.report = {
            "schema": REPORT_SCHEMA,
            "state": "created",
            "uid": backend.uid,
            "reserved_cpus": sorted(self.reserved),
            "launch_cpus": sorted(self.launch),
            "poll_ms": self.poll_ms,
            "started_utc": utc_now(),
            "finished_utc": None,
            "boot_id": identity["boot_id"],
            "pid_namespace": dict(identity["pid_namespace"]),
            "global_lock": lock_identity,
            "command": [],
            "command_sha256": None,
            "command_identity": None,
            "child": None,
            "child_processes": [],
            "child_returncode": None,
            "subreaper_pid": None,
            "received_signal": None,
            "records": [],
            "process_provenance": [],
            "events": self.events,
            "execution": {
                "nonce": nonce,
                "child_ready_monotonic_ns": None,
                "child_finished_monotonic_ns": None,
            },
            "uncertainty": False,
            "accepted": False,
            "error": None,
        }
        self._write(exclusive=True)

    @classmethod
    def recovery_view(cls, report, report_path, backend, writer, global_lock,
                      sleep_function):
        """Rehydrate journal state without creating or replacing the report."""
        value = cls.__new__(cls)
        value.backend = backend
        value.reserved = set(report["reserved_cpus"])
        value.launch = set(report["launch_cpus"])
        value.report_path = Path(report_path)
        value.poll_ms = report["poll_ms"]
        value.writer = writer
        value.seal_writer = atomic_json
        value.sleep = sleep_function
        value.launcher = None
        value.global_lock = global_lock
        value.records = {}
        value.record_order = []
        for record in report["records"]:
            key = thread_from_record(record).key()
            value.records[key] = record
            value.record_order.append(key)
        value.process_provenance = {}
        value.provenance_order = []
        for provenance in report["process_provenance"]:
            key = (provenance["pid"], provenance["process_starttime_ticks"],
                   provenance["process_device"], provenance["process_inode"])
            value.process_provenance[key] = provenance
            value.provenance_order.append(key)
        value.seen_keys = set(value.record_order)
        value.events = report["events"]
        value.child = None
        value.pinned_command = None
        value.child_processes = {
            child_process_key(item): item for item in report["child_processes"]}
        value.child_process_order = [
            child_process_key(item) for item in report["child_processes"]]
        # Retain the original supervisor PID only to classify its single main
        # thread record.  Recovery never claims to be that process's subreaper.
        value.subreaper_pid = report["subreaper_pid"]
        value.received_signal = report["received_signal"]
        value.signal_boundary_closed = True
        value.gate_descriptor = None
        value.initialized = True
        value.write_count = 0
        value.report = report
        return value

    def _invalidate_acceptance(self):
        """Make a previously sealed report unusable before changing its body."""
        invalidate_acceptance_path(self.report_path)

    def _seal_acceptance(self):
        require(self.report["accepted"] is True and
                self.report["state"] == "complete",
                "cannot seal an unaccepted affinity report")
        require(not ambiguity_path(self.report_path).exists(),
                "cannot seal an affinity report with ambiguous durability")
        report_bytes = self.report_path.read_bytes()
        require(report_bytes == canonical_bytes(self.report),
                "retained report differs from the accepted in-memory snapshot")
        payload = {
            "schema": ACCEPTANCE_SCHEMA,
            "committed_utc": utc_now(),
            "report_path": str(self.report_path.resolve(strict=True)),
            "report_size": len(report_bytes),
            "report_sha256": sha256_bytes(report_bytes),
        }
        payload["digest"] = sha256_value(payload)
        seal = acceptance_path(self.report_path)
        try:
            self.seal_writer(seal, payload, exclusive=True)
            validate_acceptance(payload, self.report_path, report_bytes)
        except BaseException as commit_error:
            # A writer may report a post-install fsync error.  Never let that
            # possibly-installed seal disagree with the returned result.  If
            # durable rejection succeeds, propagate the commit error.  If
            # rejection itself is ambiguous, read back the exact path-bound
            # pair: a fully valid pair defines a committed success; every other
            # state remains rejected and raises.  Thus no caller-visible failure
            # can coexist with evidence that verify-report accepts.
            try:
                self._invalidate_acceptance()
            except BaseException as invalidation_error:
                try:
                    retained, retained_bytes, _seal, _seal_bytes = \
                        accepted_report_snapshot(self.report_path)
                    require(retained_bytes == report_bytes and
                            retained == self.report,
                            "read-back acceptance differs from committing snapshot")
                    return
                except BaseException as readback_error:
                    raise IsolationError(
                        "acceptance commit, invalidation, and read-back failed: "
                        "{}; {}; {}".format(
                            commit_error, invalidation_error, readback_error)) \
                        from commit_error
            raise

    def _write(self, exclusive=False):
        if hasattr(self.global_lock, "validate_current"):
            current_lock = self.global_lock.validate_current()
            require(current_lock == self.report["global_lock"],
                    "global affinity lock changed during transaction")
        self.report["records"] = [self.records[key] for key in self.record_order]
        self.report["process_provenance"] = [
            self.process_provenance[key] for key in self.provenance_order]
        self.report["child_processes"] = [
            self.child_processes[key] for key in self.child_process_order]
        validate_report(self.report)
        self.writer(self.report_path, self.report, exclusive=exclusive)
        self.write_count += 1

    def _event(self, name, identity=None, value=None):
        event = {
            "event": name,
            "pid": None if identity is None else identity.pid,
            "tid": None if identity is None else identity.tid,
            "value": value,
        }
        validate_event(event)
        self.events.append(event)

    def prepare_command(self, command):
        require(not self.report["command"], "child command was already prepared")
        command = list(command)
        pinned = PinnedCommand(command)
        try:
            self.report["command"] = command
            self.report["command_sha256"] = sha256_value(command)
            self.report["command_identity"] = json.loads(json.dumps(
                pinned.identity))
            self.report["state"] = "prepared"
            self._write()
            self.pinned_command = pinned
        except BaseException:
            pinned.close()
            raise

    @staticmethod
    def _provenance(identity, affinity, source):
        return {
            "pid": identity.pid,
            "process_starttime_ticks": identity.process_starttime_ticks,
            "process_device": identity.process_device,
            "process_inode": identity.process_inode,
            "ppid": identity.ppid,
            "original_cpus": sorted(affinity),
            "source": source,
        }

    @staticmethod
    def _record(identity, original, restricted, provenance, move_status="planned"):
        return {
            "pid": identity.pid,
            "tid": identity.tid,
            "starttime_ticks": identity.starttime_ticks,
            "uid": identity.uid,
            "session": identity.session,
            "ppid": identity.ppid,
            "process_starttime_ticks": identity.process_starttime_ticks,
            "task_device": identity.task_device,
            "task_inode": identity.task_inode,
            "process_device": identity.process_device,
            "process_inode": identity.process_inode,
            "original_cpus": sorted(original),
            "restricted_cpus": sorted(restricted),
            "provenance": provenance,
            "move_count": 0,
            "move_status": move_status,
            "restore_status": "pending",
        }

    def request_signal(self, signum):
        if self.signal_boundary_closed:
            return False
        if isinstance(signum, signal.Signals):
            signum = int(signum)
        signum = exact_int(signum, "signal", 1, 255)
        if self.received_signal is None:
            self.received_signal = signum
            self.report["received_signal"] = signum
            self._event("signal", value=signum)
        if self.gate_descriptor is not None:
            descriptor = self.gate_descriptor
            self.gate_descriptor = None
            try:
                os.close(descriptor)
            except OSError:
                pass
        return True

    def _snapshot(self, excluded_child=None):
        result = []
        identities = self.backend.list_threads(excluded_session=None)
        excluded_processes = child_scope_process_keys(
            identities, excluded_child,
            retained=[self.child_processes[key]
                      for key in self.child_process_order],
            subreaper_pid=self.subreaper_pid)
        for identity in identities:
            if identity.process_key() in excluded_processes:
                continue
            try:
                affinity = self.backend.get_affinity(identity)
            except ThreadGone:
                self._event("gone_before_capture", identity)
                continue
            except ThreadReused:
                self._event("reused_before_capture", identity)
                continue
            except ThreadMutationUncertain:
                self.report["uncertainty"] = True
                self._event("tid_reuse_uncertainty", identity, "get_affinity")
                raise
            result.append((identity, affinity))
        return result

    def _add_process_provenance(self, identity, affinity, source):
        key = identity.process_key()
        existing = self.process_provenance.get(key)
        if existing is not None:
            require(set(existing["original_cpus"]) == set(affinity),
                    "process provenance became nonuniform")
            return False
        value = self._provenance(identity, affinity, source)
        self.process_provenance[key] = value
        self.provenance_order.append(key)
        return True

    def initialize_provenance(self):
        require(not self.initialized, "process provenance was already initialized")
        snapshot = self._snapshot()
        by_process = {}
        for identity, affinity in snapshot:
            by_process.setdefault(identity.process_key(), []).append((identity, affinity))
        require(by_process, "same-UID process snapshot is empty")
        for key in sorted(by_process):
            members = by_process[key]
            unique = {tuple(sorted(affinity)) for _, affinity in members}
            require(len(unique) == 1,
                    "same-UID process {} has nonuniform thread affinities; "
                    "refusing mutation because newcomer creators are ambiguous".format(key[0]))
            original = set(next(iter(unique)))
            self._add_process_provenance(members[0][0], original, "baseline")
            self.seen_keys.update(identity.key() for identity, _ in members)
        self.initialized = True
        self._apply_snapshot(snapshot, baseline=True)

    def _parent_original(self, identity):
        existing = self.process_provenance.get(identity.process_key())
        if existing is not None:
            return set(existing["original_cpus"]), existing["source"], False
        try:
            parent = self.backend.process_identity(identity.ppid)
        except (ThreadGone, ThreadReused):
            parent = None
        if parent is not None:
            provenance = self.process_provenance.get(parent.process_key())
            if provenance is not None:
                original = set(provenance["original_cpus"])
                return original, "parent_process", True
        return None, None, False

    def _apply_snapshot(self, snapshot, baseline=False):
        pending_moves = []
        inherited_records = []
        events_before = len(self.events)

        # Preflight every newcomer before changing any mask in this scan.
        for identity, affinity in snapshot:
            key = identity.key()
            record = self.records.get(key)
            if record is not None:
                if affinity.intersection(self.reserved):
                    restricted = set(record["restricted_cpus"])
                    require(restricted,
                            "recorded thread {}/{} has no housekeeping CPU".format(
                                identity.pid, identity.tid))
                    pending_moves.append((identity, affinity, restricted, record, False))
                continue

            if baseline:
                original = set(affinity)
                provenance = "baseline"
            elif key in self.seen_keys:
                # A known thread deliberately changed its mask while the child
                # ran.  The observed pre-mutation mask is exact.
                original = set(affinity)
                provenance = "observed"
            else:
                original, source, is_new_process = self._parent_original(identity)
                if original is None:
                    if affinity.intersection(self.reserved):
                        # We observe the exact unrestricted mask before touching it.
                        original = set(affinity)
                        source = "observed"
                        is_new_process = True
                    else:
                        raise UnsafeNewcomer(
                            identity,
                            "new process {}/{} has an already restricted mask but no "
                            "safe uniform creator provenance".format(
                                identity.pid, identity.tid))
                expected = original.difference(self.reserved)
                if (affinity != original and affinity != expected):
                    raise UnsafeNewcomer(
                        identity,
                        "new thread {}/{} does not match either its proven original "
                        "or inherited restricted affinity".format(
                            identity.pid, identity.tid))
                provenance = "inherited" if affinity == expected and \
                    original.intersection(self.reserved) else "observed"
                if is_new_process:
                    self._add_process_provenance(identity, original, source)

            self.seen_keys.add(key)
            if not original.intersection(self.reserved):
                continue
            restricted = original.difference(self.reserved)
            require(restricted,
                    "thread {}/{} has no CPU outside the reserved set".format(
                        identity.pid, identity.tid))
            if affinity == restricted:
                record = self._record(
                    identity, original, restricted, provenance,
                    move_status="inherited_restricted")
                self.records[key] = record
                self.record_order.append(key)
                inherited_records.append(record)
            else:
                record = self._record(identity, original, restricted, provenance)
                self.records[key] = record
                self.record_order.append(key)
                pending_moves.append((identity, original, restricted, record, True))

        if not pending_moves and not inherited_records and \
                len(self.events) == events_before:
            # Provenance for a process that needs no affinity change is useful
            # in memory for its future children, but it is not crash-recovery
            # state yet.  Defer that audit-only update until the next necessary
            # write instead of fsyncing merely because an unrelated process was
            # created during a timed benchmark.
            return 0, 0

        # Every exact original is durable before the first syscall in this scan.
        self.report["state"] = "moving"
        self._write()
        for identity, _original, restricted, record, _new_record in pending_moves:
            record["restricted_cpus"] = sorted(restricted)
            record["move_status"] = "planned"
            record["restore_status"] = "pending"
            try:
                self.backend.set_affinity(identity, restricted)
            except ThreadGone:
                record["move_status"] = "gone_before_move"
                record["restore_status"] = "gone"
            except ThreadReused:
                record["move_status"] = "reused_before_move"
                record["restore_status"] = "reused"
                self.report["uncertainty"] = True
                self._event("tid_reuse_uncertainty", identity, "before_move")
                raise IsolationError("TID was reused before affinity mutation")
            except ThreadMutationUncertain:
                record["move_status"] = "uncertain_after_move"
                record["restore_status"] = "uncertain"
                self.report["uncertainty"] = True
                self._event("tid_reuse_uncertainty", identity, "after_move")
                raise
            except IsolationError:
                record["move_status"] = "move_failed"
                # A failed syscall may have partially taken effect.  Keep it
                # pending so normal cleanup and recovery both retry the original.
                record["restore_status"] = "pending"
                raise
            else:
                record["move_count"] += 1
                record["move_status"] = "moved"
        self._write()
        return len(pending_moves), len(inherited_records)

    def restrict_once(self, excluded_child=None):
        require(self.initialized, "process provenance is not initialized")
        return self._apply_snapshot(self._snapshot(excluded_child), baseline=False)

    def stabilize(self):
        require(bool(self.report["command"]), "child command is not prepared")
        if not self.initialized:
            self.initialize_provenance()
        stable = 0
        for _ in range(MAX_STABILIZE_PASSES):
            if self.received_signal is not None:
                raise IsolationError("received signal before child launch")
            moved, inherited = self.restrict_once()
            if moved == 0 and inherited == 0:
                stable += 1
                if stable == 2:
                    self.report["state"] = "isolated"
                    self._write()
                    return
            else:
                stable = 0
            self.sleep(self.poll_ms / 1000.0)
        raise IsolationError("same-UID thread set did not stabilize")

    def _child_ready(self, identity):
        require(self.received_signal is None,
                "signal arrived before gated child release")
        child = child_record_from_identity(identity, False)
        self.report["child"] = child
        self.report["execution"]["child_ready_monotonic_ns"] = time.monotonic_ns()
        process = child_process_record(identity)
        key = child_process_key(process)
        self.child_processes[key] = process
        self.child_process_order.append(key)
        self.report["state"] = "launch_gated"
        # This durable write closes the SIGKILL launch window.  Until it
        # succeeds, the only path to child execution remains blocked by a pipe
        # whose writer dies with this supervisor.
        self._write()

    def _observe_child_processes(self):
        if self.report["child"] is None:
            return 0
        retained = [self.child_processes[key] for key in self.child_process_order]
        members = self.backend.child_scope_processes(
            self.report["child"], retained=retained,
            subreaper_pid=self.subreaper_pid)
        added = 0
        for identity in members:
            process = child_process_record(identity)
            key = child_process_key(process)
            if key in self.child_processes:
                continue
            self.child_processes[key] = process
            self.child_process_order.append(key)
            added += 1
        if added:
            self._event("child_descendants", value=added)
            self._write()
        return added

    def _mark_child_released(self):
        require(self.report["child"] is not None,
                "cannot release an unrecorded child")
        self.report["child"]["released"] = True
        self.report["state"] = "running"
        self._write()

    def _authorize_child_release(self, descriptor):
        require(self.received_signal is None,
                "signal arrived before gated child release")
        self.gate_descriptor = descriptor
        try:
            # Recheck after publishing the descriptor to request_signal().  A
            # Python signal handler that runs after this check but before the
            # syscall closes the FD; a signal delivered during write(2) is
            # linearized either before the byte (EBADF/EINTR) or after release.
            require(self.received_signal is None,
                    "signal arrived before gated child release")
            if descriptor is not None:
                written = os.write(descriptor, b"G")
                require(written == 1, "gated child release write was incomplete")
        finally:
            if self.gate_descriptor == descriptor:
                self.gate_descriptor = None

    def _cleanup_unsafe_newcomer(self, error):
        self._event("unsafe_newcomer", error.identity, str(error))
        failures = terminate_process(
            self.backend, error.identity, sleep_function=self.sleep)
        if failures:
            self.report["uncertainty"] = True
            return failures
        self._event("process_terminated", error.identity, "unsafe_newcomer")
        return []

    def restore(self, restored_status="restored"):
        failures = []
        for key in reversed(self.record_order):
            record = self.records[key]
            if record["restore_status"] not in {"pending", "failed"}:
                continue
            identity = thread_from_record(record)
            try:
                self.backend.current_identity(identity)
            except ThreadGone:
                record["restore_status"] = "gone"
                continue
            except ThreadReused:
                record["restore_status"] = "reused"
                self.report["uncertainty"] = True
                failures.append(
                    "TID {}/{} was reused before restore".format(
                        identity.pid, identity.tid))
                continue
            try:
                self.backend.set_affinity(identity, set(record["original_cpus"]))
                record["restore_status"] = restored_status
            except ThreadGone:
                record["restore_status"] = "gone"
            except ThreadReused:
                record["restore_status"] = "reused"
                self.report["uncertainty"] = True
                failures.append(
                    "TID {}/{} was reused during restore".format(
                        identity.pid, identity.tid))
            except ThreadMutationUncertain as error:
                record["restore_status"] = "uncertain"
                self.report["uncertainty"] = True
                self._event("tid_reuse_uncertainty", identity, "after_restore")
                failures.append(str(error))
            except IsolationError as error:
                record["restore_status"] = "failed"
                failures.append(str(error))
        return failures

    def _plan_late_restores(self):
        """Journal restricted late arrivals before widening their exact masks."""
        additions = 0
        snapshot = self._snapshot()
        for identity, affinity in snapshot:
            key = identity.key()
            record = self.records.get(key)
            if record is not None:
                original = set(record["original_cpus"])
                restricted = set(record["restricted_cpus"])
                if affinity == original:
                    continue
                if affinity == restricted:
                    record["restore_status"] = "pending"
                    additions += 1
                    continue
                self.report["uncertainty"] = True
                record["restore_status"] = "uncertain"
                raise IsolationError(
                    "journaled thread {}/{} changed to an unprovable mask during "
                    "restoration".format(identity.pid, identity.tid))

            original, source, is_new_process = self._parent_original(identity)
            if original is None:
                if affinity.intersection(self.reserved):
                    original = set(affinity)
                    source = "observed"
                    is_new_process = True
                else:
                    raise UnsafeNewcomer(
                        identity,
                        "late process {}/{} retained a restricted mask without "
                        "exact creator provenance".format(identity.pid, identity.tid))
            expected = original.difference(self.reserved)
            if affinity == original:
                if is_new_process:
                    self._add_process_provenance(identity, original, source)
                self.seen_keys.add(key)
                continue
            if (original.intersection(self.reserved) and affinity == expected and
                    expected):
                if is_new_process:
                    self._add_process_provenance(identity, original, source)
                record = self._record(
                    identity, original, expected, "inherited",
                    move_status="inherited_restricted")
                self.records[key] = record
                self.record_order.append(key)
                self.seen_keys.add(key)
                additions += 1
                continue
            raise UnsafeNewcomer(
                identity,
                "late thread {}/{} has neither its proven original nor inherited "
                "restricted affinity".format(identity.pid, identity.tid))
        return additions

    def restore_reconciled(self, restored_status="restored"):
        """Restore and perform bounded best-effort late-arrival reconciliation.

        For an outside creator whose mask was changed, no finite scan count is
        a proof of complete restoration: clone may have captured a restricted
        mask before widening and publish after the final scan.  The supervisor's
        own single main thread is the one bounded exception: this implementation
        controls its fork point and creates no other thread.  Keep scanning for
        the full bound after any outside mutation, then return an explicit
        containment failure.
        """
        failures = self.restore(restored_status=restored_status)
        stable = 0
        for _ in range(MAX_STABILIZE_PASSES):
            try:
                additions = self._plan_late_restores()
            except UnsafeNewcomer as error:
                failures.extend(self._cleanup_unsafe_newcomer(error))
                if not failures:
                    failures.append(str(error))
                stable = 0
                continue
            if additions:
                # Exact originals and identities are durable before widening.
                self.report["state"] = "restoring" if \
                    restored_status == "restored" else "recovering"
                self._write()
                failures.extend(self.restore(restored_status=restored_status))
                stable = 0
            else:
                stable += 1
                if stable == 2 and not uncontained_mutation_records(
                        self.records.values(), self.subreaper_pid):
                    return failures
            self.sleep(self.poll_ms / 1000.0)
        if uncontained_mutation_records(
                self.records.values(), self.subreaper_pid):
            self.report["uncertainty"] = True
            if MUTATED_RESTORE_UNPROVABLE not in failures:
                failures.append(MUTATED_RESTORE_UNPROVABLE)
        else:
            failures.append(
                "same-UID thread set did not stabilize during restoration")
        return failures

    def _finish_child_wait(self):
        if self.child is None:
            return None
        result = self.child.poll()
        if result is None:
            try:
                result = self.child.wait(KILL_GRACE_SECONDS)
            except TimeoutError:
                result = self.child.poll()
        return result

    def _capture_blocked_signals(self):
        if not hasattr(signal, "sigpending"):
            return
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            if signum in signal.sigpending():
                self.request_signal(signum)

    def run(self, command):
        if not command:
            raise IsolationError("no child command was provided")
        child_returncode = None
        run_error = None
        cleanup_failures = []
        restore_failures = []
        write_failures = []
        descendants_observed = False
        try:
            if hasattr(self.backend, "enable_subreaper"):
                self.subreaper_pid = self.backend.enable_subreaper()
                self.report["subreaper_pid"] = self.subreaper_pid
            self.prepare_command(list(command))
            self.stabilize()
            require(self.pinned_command is not None,
                    "prepared command descriptors are unavailable")
            self.child = self.launcher.launch(
                self.pinned_command, self.launch, self.backend,
                self._child_ready, self._authorize_child_release,
                self.report["execution"]["nonce"])
            self._mark_child_released()
            self._observe_child_processes()
            while self.child.poll() is None:
                if self.received_signal is not None:
                    break
                self._observe_child_processes()
                self.restrict_once(excluded_child=self.report["child"])
                self.sleep(self.poll_ms / 1000.0)
            child_returncode = self.child.poll()
            if self.received_signal is None:
                self._observe_child_processes()
                self.restrict_once(excluded_child=self.report["child"])
                members = self.backend.child_scope_processes(
                    self.report["child"],
                    retained=[self.child_processes[key]
                              for key in self.child_process_order],
                    subreaper_pid=self.subreaper_pid)
                if child_returncode is not None and members:
                    descendants_observed = True
                    self._event("child_descendants", value=len(members))
                    run_error = IsolationError(
                        "child leader exited while its session retained descendants")
        except UnsafeNewcomer as error:
            run_error = error
            cleanup_failures.extend(self._cleanup_unsafe_newcomer(error))
        except BaseException as error:
            run_error = error
        finally:
            old_signal_mask = None
            if hasattr(signal, "pthread_sigmask"):
                old_signal_mask = signal.pthread_sigmask(
                    signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM, signal.SIGHUP})
                self._capture_blocked_signals()
            # No affinity is restored while a process in the benchmark session
            # remains alive.  Recovery uses the same ordering.
            if self.report["child"] is not None:
                initial_signal = self.received_signal or signal.SIGTERM
                try:
                    failures, had_members = terminate_session(
                        self.backend, self.report["child"], initial_signal,
                        sleep_function=self.sleep,
                        reap_callback=self.child.poll if self.child is not None else None,
                        retained=[self.child_processes[key]
                                  for key in self.child_process_order],
                        subreaper_pid=self.subreaper_pid)
                    cleanup_failures.extend(failures)
                    if had_members:
                        self._event("child_session_cleanup", value=initial_signal)
                except BaseException as error:
                    cleanup_failures.append(str(error))
                child_returncode = self._finish_child_wait()
            self.report["child_returncode"] = child_returncode
            if self.report["execution"]["child_ready_monotonic_ns"] is not None:
                self.report["execution"]["child_finished_monotonic_ns"] = \
                    time.monotonic_ns()
            self.report["state"] = "restoring"
            try:
                self._write()
            except BaseException as error:
                # A pre-restore journal write failure is evidence failure, but
                # it must never suppress the independently attempted restore.
                write_failures.append(str(error))
            if not cleanup_failures:
                try:
                    restore_failures = self.restore_reconciled()
                except BaseException as error:
                    restore_failures.append(str(error))
            else:
                restore_failures.append(
                    "affinity restore withheld until child-session cleanup succeeds")
            errors = []
            if run_error is not None:
                errors.append(str(run_error))
            errors.extend(cleanup_failures)
            errors.extend(restore_failures)
            errors.extend(write_failures)
            if child_returncode not in (None, 0):
                errors.append("child returned {}".format(child_returncode))
            if self.received_signal is not None:
                errors.append("received signal {}".format(self.received_signal))
            if descendants_observed:
                errors.append("child-session descendants required forced cleanup")
            self._capture_blocked_signals()
            if self.received_signal is not None and not any(
                    item.startswith("received signal ") for item in errors):
                errors.append("received signal {}".format(self.received_signal))
            self.report["error"] = "; ".join(errors) if errors else None
            self.report["accepted"] = (
                not errors and not self.report["uncertainty"] and
                child_returncode == 0 and self.received_signal is None)
            self.report["state"] = "complete" if self.report["accepted"] else "failed"
            self.report["finished_utc"] = utc_now()
            try:
                self._write()
            except BaseException as error:
                write_failures.append(str(error))
                self.report["accepted"] = False
                self.report["state"] = "failed"
            # A signal arriving during the durable final write is pending while
            # blocked.  Observe it and rewrite a failed report before defining
            # the transaction as complete.  Signals delivered after unblocking
            # occur after this durability boundary, not during the benchmark
            # transaction.
            previous_signal = self.received_signal
            self._capture_blocked_signals()
            if self.received_signal is not None and previous_signal is None:
                self.report["accepted"] = False
                self.report["state"] = "failed"
                self.report["finished_utc"] = utc_now()
                signal_error = "received signal {}".format(self.received_signal)
                self.report["error"] = signal_error if self.report["error"] is None else \
                    self.report["error"] + "; " + signal_error
                try:
                    self._write()
                except BaseException as error:
                    write_failures.append(str(error))
            # This is the signal-observation boundary for the supervised
            # transaction.  The benchmark child is gone, best-effort restoration
            # and the final pending-signal rewrite are complete, and signals stay
            # blocked until after the acceptance commit.  A later signal belongs
            # to the caller, not to the already-finished benchmark evidence.
            self.signal_boundary_closed = True
            if self.report["accepted"]:
                try:
                    self._seal_acceptance()
                except BaseException as error:
                    write_failures.append(str(error))
                    self.report["accepted"] = False
                    self.report["state"] = "failed"
                    self.report["finished_utc"] = utc_now()
                    self.report["error"] = "acceptance seal failed: {}".format(error)
                    try:
                        self._write()
                    except BaseException as rewrite_error:
                        write_failures.append(str(rewrite_error))
            if old_signal_mask is not None:
                signal.pthread_sigmask(signal.SIG_SETMASK, old_signal_mask)
            if self.pinned_command is not None:
                self.pinned_command.close()
                self.pinned_command = None
        if write_failures or run_error is not None or cleanup_failures or restore_failures:
            return 1
        if self.received_signal is not None:
            return 128 + self.received_signal
        return child_returncode if child_returncode is not None else 1


def load_json(path, label):
    path = Path(path)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise IsolationError("cannot read {}: {}".format(label, error)) from error


def load_report(path, require_accepted=False):
    if require_accepted:
        return accepted_report_snapshot(path)[0]
    return validate_report(load_json(path, "affinity report"))


def report_host_matches(report, host=None):
    current = runtime_identity() if host is None else host
    require(current == {
        "boot_id": report["boot_id"],
        "pid_namespace": report["pid_namespace"],
    }, "refusing affinity recovery across a boot or PID-namespace boundary")


def recover_report(path, backend=None, writer=atomic_json, host=None,
                   sleep_function=time.sleep, global_lock=None):
    path = Path(path)
    report = load_report(path)
    uid = exact_int(report["uid"], "report UID", 0, 2**31 - 1)
    require(uid == os.getuid() or backend is not None,
            "refusing to restore another user's report")
    report_host_matches(report, host)
    backend = LinuxThreadBackend(uid=uid) if backend is None else backend
    require(backend.uid == uid, "recovery backend UID differs from report UID")
    require(global_lock is not None,
            "affinity recovery requires the same-UID global lock")
    current_lock = global_lock.validate_current() if hasattr(
        global_lock, "validate_current") else dict(global_lock)
    validate_global_lock(current_lock, uid)
    require(current_lock == report["global_lock"],
            "recovery global lock differs from the journaled lock identity")
    failures = []
    write_failures = []
    try:
        invalidate_acceptance_path(path)
    except BaseException as error:
        write_failures.append("acceptance invalidation failed: {}".format(error))
    report["state"] = "recovering"
    report["accepted"] = False
    report["finished_utc"] = None
    report["error"] = None
    try:
        writer(path, report)
    except BaseException as error:
        write_failures.append(str(error))

    # A released or gated child may still exist after supervisor SIGKILL.
    # Verify and empty that exact session before widening any restored mask.
    cleanup_failures = []
    try:
        cleanup_failures, _ = terminate_session(
            backend, report["child"], signal.SIGTERM,
            sleep_function=sleep_function,
            retained=report["child_processes"])
    except BaseException as error:
        cleanup_failures.append(str(error))
    failures.extend(cleanup_failures)

    if not cleanup_failures:
        recovery = AffinityTransaction.recovery_view(
            report, path, backend, writer, global_lock, sleep_function)
        failures.extend(recovery.restore_reconciled(
            restored_status="restored_after_recovery"))
    else:
        failures.append(
            "affinity restore withheld until child-session cleanup succeeds")

    # A replacement supervisor is not the original subreaper.  Once the gate
    # was released, a double-forked descendant may have changed sessions and
    # been reparented elsewhere before its identity reached the journal.  Known
    # identities are still cleaned and masks are restored best-effort, but a
    # successful recovery claim would be unsound without persistent containment.
    if report["child"] is not None and report["child"]["released"]:
        report["uncertainty"] = True
        failures.append(CRASH_SCOPE_UNPROVABLE)
    if uncontained_mutation_records(
            report["records"], report["subreaper_pid"]):
        report["uncertainty"] = True
        if MUTATED_RESTORE_UNPROVABLE not in failures:
            failures.append(MUTATED_RESTORE_UNPROVABLE)

    unresolved = [
        "{}/{}:{}".format(item["pid"], item["tid"], item["restore_status"])
        for item in report["records"]
        if item["restore_status"] in {"pending", "failed", "reused", "uncertain"}
    ]
    if unresolved:
        failures.append("unresolved affinity records: {}".format(",".join(unresolved)))
    if report["uncertainty"] and not unresolved and not failures:
        failures.append("journal retains an unrecoverable uncertainty")

    failures.extend(write_failures)
    report["state"] = "recovered" if not failures else "recovery_failed"
    report["accepted"] = False
    report["finished_utc"] = utc_now()
    report["error"] = None if not failures else "; ".join(failures)
    try:
        validate_report(report)
        writer(path, report)
    except BaseException as error:
        failures.append(str(error))
    if failures:
        raise IsolationError("; ".join(failures))
    print("Leopard2 affinity report restored; benchmark evidence remains unaccepted")
    return 0


def verify_signed_json(value, label):
    require(isinstance(value, dict), "{} is not an object".format(label))
    payload = dict(value)
    digest = payload.pop("digest", None)
    require(isinstance(digest, str) and HEX256.fullmatch(digest) is not None and
            digest == sha256_value(payload), "{} digest is invalid".format(label))
    return value


def safe_relative(root, relative, label):
    require(isinstance(relative, str) and relative and not os.path.isabs(relative),
            "{} path is not relative".format(label))
    root = Path(root).resolve()
    result = (root / relative).resolve()
    try:
        result.relative_to(root)
    except ValueError as error:
        raise IsolationError("{} path escapes its evidence root".format(label)) from error
    return result


def parse_long_options(arguments):
    """Parse the exact value-taking option shape used by main_compare/run_abba.py."""
    allowed = {
        "--baseline", "--candidate", "--baseline-archive", "--candidate-archive",
        "--baseline-build-dir", "--candidate-build-dir", "--baseline-source-root",
        "--candidate-source-root", "--candidate-commit", "--candidate-mode",
        "--reservation-file", "--output", "--cpu", "--reserved-sibling",
        "--taskset", "--ldd", "--preset", "--cell", "--reuse", "--iterations",
        "--warmup", "--timeout",
    }
    values = {}
    index = 0
    while index < len(arguments):
        option = arguments[index]
        require(option in allowed and index + 1 < len(arguments),
                "supervised main-comparison command has an unknown/missing option")
        value = arguments[index + 1]
        require(not value.startswith("--"),
                "supervised main-comparison option has no value")
        if option == "--cell":
            values.setdefault(option, []).append(value)
        else:
            require(option not in values,
                    "supervised main-comparison option is duplicated: {}".format(option))
            values[option] = value
        index += 2
    required = allowed.difference({"--taskset", "--ldd", "--preset", "--cell",
                                   "--candidate-mode"})
    require(required.issubset(values),
            "supervised main-comparison command omits required options")
    require(not ("--preset" in values and "--cell" in values),
            "supervised command mixes preset and custom cells")
    return values


def resolve_command_path(value, cwd):
    path = Path(value)
    if not path.is_absolute():
        path = Path(cwd) / path
    return str(path.resolve())


def validate_main_manifest_binding(report, manifest_path):
    validate_report(report)
    require(report["accepted"] is True,
            "main-comparison binding requires an accepted report snapshot")
    manifest_path, manifest, manifest_bytes, _manifest_identity = \
        stable_json_snapshot(manifest_path, "main-comparison manifest")
    manifest = verify_signed_json(manifest, "main-comparison manifest")
    require(manifest.get("schema") == MAIN_MANIFEST_SCHEMA and
            manifest.get("valid") is True,
            "main-comparison manifest is not valid current-schema evidence")
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict) and
            set(raw_info) == {"path", "payload_digest", "sha256", "size"},
            "main-comparison manifest raw identity is malformed")
    raw_path = safe_relative(manifest_path.parent, raw_info["path"], "raw bundle")
    raw_path, raw, raw_bytes, _raw_identity = stable_json_snapshot(
        raw_path, "main-comparison raw bundle")
    require(type(raw_info["size"]) is int and raw_info["size"] == len(raw_bytes) and
            raw_info["sha256"] == sha256_bytes(raw_bytes),
            "main-comparison raw file identity differs from its manifest")
    raw = verify_signed_json(raw, "main-comparison raw bundle")
    require(raw.get("schema") == MAIN_RAW_SCHEMA,
            "main-comparison raw bundle is not current-schema evidence")
    require(raw_info["payload_digest"] == raw.get("digest"),
            "main-comparison payload digest differs from its manifest")

    command = report["command"]
    identity = report["command_identity"]
    require(identity["script"] is not None and len(command) >= 3 and
            resolve_command_path(command[1], identity["cwd"]["path"]) ==
            identity["script"]["path"] and command[2] == "run",
            "accepted supervisor command is not a direct Python main-comparison run")
    specification = raw.get("input_specification")
    identities = raw.get("identities_initial")
    campaign = raw.get("campaign")
    reservation = raw.get("reservation")
    require(isinstance(specification, dict) and isinstance(identities, dict) and
            isinstance(campaign, dict) and isinstance(reservation, dict),
            "main-comparison binding inputs are incomplete")
    runner_identity = identities.get("runner")
    require(isinstance(runner_identity, dict) and
            specification.get("runner") == identity["script"]["path"] and
            runner_identity.get("path") == identity["script"]["path"] and
            runner_identity.get("sha256") == identity["script"]["sha256"],
            "supervised runner identity differs from the retained campaign")
    options = parse_long_options(command[3:])
    cwd = identity["cwd"]["path"]
    path_options = {
        "--baseline": "baseline_executable",
        "--candidate": "candidate_executable",
        "--baseline-archive": "baseline_archive",
        "--candidate-archive": "candidate_archive",
        "--baseline-build-dir": "baseline_build_dir",
        "--candidate-build-dir": "candidate_build_dir",
        "--baseline-source-root": "baseline_source_root",
        "--candidate-source-root": "candidate_source_root",
    }
    for option, field in path_options.items():
        require(resolve_command_path(options[option], cwd) == specification.get(field),
                "supervised {} differs from retained campaign".format(option))
    for option, field, default in (
        ("--taskset", "taskset", "/usr/bin/taskset"),
        ("--ldd", "ldd", "/usr/bin/ldd"),
    ):
        require(resolve_command_path(options.get(option, default), cwd) ==
                specification.get(field),
                "supervised {} differs from retained campaign".format(option))
    require(options["--candidate-commit"] == specification.get("candidate_commit"),
            "supervised candidate commit differs from retained campaign")
    require(resolve_command_path(options["--reservation-file"], cwd) ==
            reservation.get("path"),
            "supervised reservation differs from retained campaign")
    require(resolve_command_path(options["--output"], cwd) ==
            str(manifest_path.parent),
            "supervised output differs from retained manifest location")
    require(int(options["--cpu"], 10) == campaign.get("benchmark_cpu") and
            int(options["--reserved-sibling"], 10) == campaign.get("reserved_sibling") and
            int(options["--reuse"], 10) == campaign.get("reuse") and
            int(options["--iterations"], 10) == campaign.get("iterations") and
            int(options["--warmup"], 10) == campaign.get("warmup") and
            float(options["--timeout"]) == campaign.get("timeout_seconds"),
            "supervised numeric campaign parameters differ from retained evidence")
    require(options.get("--candidate-mode", "auto") ==
            campaign.get("candidate_mode", "auto"),
            "supervised candidate mode differs from retained evidence")
    require(report["reserved_cpus"] == sorted({campaign["benchmark_cpu"],
                                                campaign["reserved_sibling"]}),
            "supervisor reserved set differs from the benchmark pair")
    supervision = raw.get("supervision")
    supervision_keys = {
        "campaign_sha256", "execution_nonce",
        "isolation_after_monotonic_ns", "isolation_before_monotonic_ns",
        "launch_cpus", "reservation_nonce", "reservation_sha256",
        "reserved_cpus", "runner_finished_monotonic_ns", "runner_pid",
        "runner_started_monotonic_ns", "schema",
    }
    require(isinstance(supervision, dict) and
            set(supervision) == supervision_keys and
            supervision.get("schema") == MAIN_SUPERVISION_SCHEMA,
            "main-comparison supervision handshake is missing or malformed")
    execution = report["execution"]
    require(supervision.get("execution_nonce") == execution["nonce"],
            "supervisor and main-comparison execution nonces differ")
    for name in ("runner_pid", "runner_started_monotonic_ns",
                 "runner_finished_monotonic_ns",
                 "isolation_before_monotonic_ns",
                 "isolation_after_monotonic_ns"):
        exact_int(supervision.get(name), "supervision {}".format(name), 0,
                  2**63 - 1)
    require(supervision["runner_pid"] >= 1 and
            supervision["runner_pid"] == report["child"]["pid"] and
            execution["child_ready_monotonic_ns"] <=
            supervision["runner_started_monotonic_ns"] <=
            supervision["isolation_before_monotonic_ns"] <=
            supervision["isolation_after_monotonic_ns"] <=
            supervision["runner_finished_monotonic_ns"] <=
            execution["child_finished_monotonic_ns"],
            "supervisor interval does not enclose the retained runner campaign")
    require(supervision.get("launch_cpus") == report["launch_cpus"] ==
            campaign.get("allowed_cpu_set_at_launch"),
            "supervisor launch set differs from the retained runner launch set")
    require(supervision.get("reserved_cpus") == report["reserved_cpus"],
            "supervision handshake uses another reserved CPU set")
    require(supervision.get("campaign_sha256") ==
            sha256_bytes(canonical_bytes(campaign)),
            "supervision handshake does not bind the campaign payload")
    reservation_payload = reservation.get("payload")
    require(isinstance(reservation_payload, dict) and
            supervision.get("reservation_sha256") == reservation.get("sha256") ==
            sha256_bytes(canonical_bytes(reservation_payload)) and
            supervision.get("reservation_nonce") == reservation_payload.get("nonce"),
            "supervision handshake does not bind the held reservation payload")
    reservation_path, reservation_bytes, _reservation_identity = \
        stable_file_snapshot(reservation.get("path"), "CPU reservation")
    require(str(reservation_path) == reservation.get("path") and
            reservation_bytes == canonical_bytes(reservation_payload) and
            sha256_bytes(reservation_bytes) == reservation.get("sha256"),
            "retained reservation file differs from the supervised payload")
    isolation = raw.get("isolation")
    require(isinstance(isolation, dict) and
            isinstance(isolation.get("before"), dict) and
            isinstance(isolation.get("after"), dict) and
            supervision["isolation_before_monotonic_ns"] ==
            isolation["before"].get("monotonic_ns") and
            supervision["isolation_after_monotonic_ns"] ==
            isolation["after"].get("monotonic_ns"),
            "supervision handshake differs from scheduler isolation interval")
    require(manifest.get("supervision") == supervision,
            "manifest supervision handshake differs from the raw bundle")
    return {
        "manifest_bytes": manifest_bytes,
        "manifest": manifest,
        "raw_bytes": raw_bytes,
        "raw": raw,
    }


def binding_payload(report_path, manifest_path):
    report_path = Path(report_path).resolve(strict=True)
    report, report_bytes, _seal, _seal_bytes = accepted_report_snapshot(report_path)
    inputs = validate_main_manifest_binding(report, manifest_path)
    manifest_path = Path(manifest_path).resolve(strict=True)
    return {
        "schema": BINDING_SCHEMA,
        "created_utc": utc_now(),
        "report": {
            "path": str(report_path), "size": len(report_bytes),
            "sha256": sha256_bytes(report_bytes), "schema": REPORT_SCHEMA,
            "command_sha256": report["command_sha256"],
        },
        "manifest": {
            "path": str(manifest_path), "size": len(inputs["manifest_bytes"]),
            "sha256": sha256_bytes(inputs["manifest_bytes"]),
            "schema": inputs["manifest"]["schema"],
            "payload_digest": inputs["manifest"]["digest"],
        },
    }


def validate_binding(value, binding_path=None, expected_manifest_path=None,
                     expected_manifest_sha256=None):
    value = verify_signed_json(value, "affinity/main binding")
    require(set(value) == {"schema", "created_utc", "report", "manifest", "digest"} and
            value["schema"] == BINDING_SCHEMA,
            "affinity/main binding has an unknown schema or fields")
    validate_utc(value["created_utc"], "binding creation time")
    report_ref = require_exact_keys(
        value["report"], {"path", "size", "sha256", "schema", "command_sha256"},
        "binding report reference")
    manifest_ref = require_exact_keys(
        value["manifest"], {"path", "size", "sha256", "schema", "payload_digest"},
        "binding manifest reference")
    require(report_ref["schema"] == REPORT_SCHEMA and
            isinstance(report_ref["command_sha256"], str) and
            HEX256.fullmatch(report_ref["command_sha256"]) is not None,
            "binding report schema/command identity is invalid")
    require(isinstance(manifest_ref["schema"], str) and manifest_ref["schema"] and
            isinstance(manifest_ref["payload_digest"], str) and
            HEX256.fullmatch(manifest_ref["payload_digest"]) is not None,
            "binding manifest schema/payload identity is invalid")
    for label, reference in (("report", report_ref), ("manifest", manifest_ref)):
        require(isinstance(reference["path"], str) and reference["path"] and
                type(reference["size"]) is int and reference["size"] >= 0 and
                isinstance(reference["sha256"], str) and
                HEX256.fullmatch(reference["sha256"]) is not None,
                "binding {} file identity is invalid".format(label))
    report_path = Path(report_ref["path"]).resolve(strict=True)
    manifest_path = Path(manifest_ref["path"]).resolve(strict=True)
    if expected_manifest_path is not None:
        require(manifest_path == Path(expected_manifest_path).resolve(strict=True),
                "affinity binding refers to another manifest path")
    if expected_manifest_sha256 is not None:
        require(isinstance(expected_manifest_sha256, str) and
                HEX256.fullmatch(expected_manifest_sha256) is not None and
                manifest_ref["sha256"] == expected_manifest_sha256,
                "affinity binding refers to another manifest snapshot")
    report, report_bytes, _seal, _seal_bytes = accepted_report_snapshot(report_path)
    inputs = validate_main_manifest_binding(report, manifest_path)
    manifest_bytes = inputs["manifest_bytes"]
    require(len(report_bytes) == report_ref["size"] and
            sha256_bytes(report_bytes) == report_ref["sha256"] and
            len(manifest_bytes) == manifest_ref["size"] and
            sha256_bytes(manifest_bytes) == manifest_ref["sha256"],
            "binding file identity changed")
    require(report["command_sha256"] == report_ref["command_sha256"],
            "binding command identity differs from accepted report")
    require(inputs["manifest"]["schema"] == manifest_ref["schema"] and
            inputs["manifest"]["digest"] == manifest_ref["payload_digest"],
            "binding manifest semantics changed")
    if binding_path is not None:
        require(Path(binding_path).resolve(strict=True).is_file(),
                "binding path is not a file")
    return value


def create_binding(report_path, manifest_path, output_path):
    output_path = Path(output_path)
    payload = binding_payload(report_path, manifest_path)
    payload["digest"] = sha256_value(payload)
    validate_binding_structure_only(payload)
    atomic_json(output_path, payload, exclusive=True)
    _path, retained, _bytes, _identity = stable_json_snapshot(
        output_path, "affinity/main binding")
    validate_binding(retained, output_path)
    print(output_path)
    return 0


def validate_binding_structure_only(value):
    # Used before installation so validation does not require output_path itself.
    verify_signed_json(value, "affinity/main binding")
    require(set(value) == {"schema", "created_utc", "report", "manifest", "digest"} and
            value["schema"] == BINDING_SCHEMA,
            "affinity/main binding has invalid fields")
    validate_utc(value["created_utc"], "binding creation time")
    return value


FAKE_HOST = {
    "boot_id": "00000000-0000-0000-0000-000000000001",
    "pid_namespace": {"device": 1, "inode": 2},
}

FAKE_LOCK_IDENTITY = {
    "path": "/run/user/{}/leopard2-affinity-supervisor.lock".format(os.getuid()),
    "uid": os.getuid(), "device": 10, "inode": 11, "mode": 0o600,
    "directory_device": 12, "directory_inode": 13,
    "lock": "exclusive_nonblocking_same_uid",
}


class FakeGlobalLock:
    """Deterministic same-UID lock model; never takes an operating-system lock."""

    held = False

    def __init__(self, identity=None):
        self.identity = dict(FAKE_LOCK_IDENTITY if identity is None else identity)
        self.active = False

    def __enter__(self):
        if FakeGlobalLock.held:
            raise IsolationError("another same-UID affinity supervisor is active")
        FakeGlobalLock.held = True
        self.active = True
        return self

    def validate_current(self):
        require(self.active, "fake global affinity lock is not held")
        return dict(self.identity)

    def __exit__(self, _exc_type, _exc, _traceback):
        self.active = False
        FakeGlobalLock.held = False


class FakeBackend:
    """Deterministic backend; it never calls a real affinity or signal syscall."""

    def __init__(self):
        self.uid = os.getuid()
        self.threads = {}
        self.fail_set_once = set()
        self.reuse_after_set = set()
        self.stubborn = set()
        self.immortal = set()
        self.set_log = []
        self.signal_log = []
        self.operation_log = []
        self.on_set = None
        self.on_list = None

    @staticmethod
    def enable_subreaper():
        return 500

    def add(self, pid, tid, start, session, cpus, ppid=1, process_start=None,
            task_device=101, task_inode=None, process_device=102,
            process_inode=None):
        if process_start is None:
            process_start = start if tid == pid else self._process_start(pid, start)
        if task_inode is None:
            task_inode = pid * 1000000000 + tid * 100000 + start
        if process_inode is None:
            existing = [value[0] for key, value in self.threads.items()
                        if key[0] == pid]
            process_inode = existing[0].process_inode if existing else \
                pid * 1000000000 + process_start
        identity = ThreadIdentity(
            pid, tid, start, self.uid, session, ppid, process_start,
            task_device, task_inode, process_device, process_inode)
        self.threads[(pid, tid)] = [identity, set(cpus)]
        return identity

    def _process_start(self, pid, fallback):
        members = [value[0] for key, value in self.threads.items() if key[0] == pid]
        return members[0].process_starttime_ticks if members else fallback

    def remove_process(self, pid):
        for key in [key for key in self.threads if key[0] == pid]:
            del self.threads[key]

    def list_threads(self, excluded_session=None):
        result = [value[0] for _, value in sorted(self.threads.items())
                  if value[0].session != excluded_session]
        if self.on_list is not None:
            self.on_list(list(result))
        return result

    def process_identity(self, pid):
        value = self.threads.get((pid, pid))
        if value is None:
            raise ThreadGone()
        return value[0]

    def current_identity(self, identity):
        value = self.threads.get((identity.pid, identity.tid))
        if value is None:
            raise ThreadGone()
        if value[0].key() != identity.key() or \
                value[0].process_key() != identity.process_key():
            raise ThreadReused()
        return value[0]

    def get_affinity(self, identity):
        self.current_identity(identity)
        return set(self.threads[(identity.pid, identity.tid)][1])

    def set_affinity(self, identity, cpus):
        self.current_identity(identity)
        key = (identity.pid, identity.tid)
        if key in self.fail_set_once:
            self.fail_set_once.remove(key)
            raise IsolationError("injected set failure")
        self.threads[key][1] = set(cpus)
        self.set_log.append((identity.key(), tuple(sorted(cpus))))
        self.operation_log.append(("set", identity.pid, tuple(sorted(cpus))))
        if self.on_set is not None:
            self.on_set(identity, set(cpus))
        if key in self.reuse_after_set:
            self.reuse_after_set.remove(key)
            old = self.threads[key][0]
            self.threads[key][0] = ThreadIdentity(
                old.pid, old.tid, old.starttime_ticks + 100000, old.uid,
                old.session, old.ppid, old.process_starttime_ticks + 100000,
                old.task_device, old.task_inode + 100000,
                old.process_device, old.process_inode + 100000)
            raise ThreadMutationUncertain("injected post-syscall TID reuse")

    def session_processes(self, session_id):
        result = []
        seen = set()
        for identity, _affinity in self.threads.values():
            key = identity.process_key()
            if identity.session == session_id and key not in seen:
                leader = self.threads.get((identity.pid, identity.pid))
                if leader is not None:
                    result.append(leader[0])
                    seen.add(key)
        return sorted(result, key=lambda item: item.pid)

    def child_scope_processes(self, child_record, retained=(), subreaper_pid=None):
        identities = self.list_threads()
        keys = child_scope_process_keys(
            identities, child_record, retained=retained,
            subreaper_pid=subreaper_pid)
        result = []
        seen = set()
        for identity in identities:
            if identity.process_key() in keys and identity.process_key() not in seen:
                leader = self.threads.get((identity.pid, identity.pid))
                if leader is not None:
                    result.append(leader[0])
                    seen.add(identity.process_key())
        return sorted(result, key=lambda item: item.pid)

    def signal_process(self, identity, signum):
        current = self.process_identity(identity.pid)
        if current.process_key() != identity.process_key():
            raise ThreadReused()
        self.signal_log.append((identity.pid, int(signum)))
        self.operation_log.append(("signal", identity.pid, int(signum)))
        if identity.pid in self.immortal:
            return
        if signum == signal.SIGKILL or identity.pid not in self.stubborn:
            self.remove_process(identity.pid)


class CountingWriter:
    def __init__(self, fail_state_once=None):
        self.calls = 0
        self.fail_state_once = fail_state_once
        self.failed = False

    def __call__(self, path, value, exclusive=False):
        self.calls += 1
        if (not self.failed and self.fail_state_once is not None and
                value.get("state") == self.fail_state_once):
            self.failed = True
            raise OSError("injected journal write failure")
        atomic_json(path, value, exclusive=exclusive)


class PostInstallFailWriter:
    def __call__(self, path, value, exclusive=False):
        atomic_json(path, value, exclusive=exclusive)
        raise OSError("injected post-install durability failure")


class FakeChild:
    def __init__(self, backend, identity, descendant=False, stay_alive=False,
                 reparent_on_exit=False):
        self.backend = backend
        self.identity = identity
        self.pid = identity.pid
        self.returncode = None
        self.descendant = descendant
        self.stay_alive = stay_alive
        self.reparent_on_exit = reparent_on_exit
        self.poll_count = 0

    def poll(self):
        if self.returncode is not None:
            return self.returncode
        self.poll_count += 1
        try:
            self.backend.process_identity(self.pid)
        except ThreadGone:
            self.returncode = -signal.SIGTERM if self.stay_alive else 0
            return self.returncode
        if not self.stay_alive:
            # The leader exits.  A deliberately retained descendant remains in
            # the same session to exercise bounded cleanup and evidence failure.
            if self.descendant:
                del self.backend.threads[(self.pid, self.pid)]
                if self.reparent_on_exit and (901, 901) in self.backend.threads:
                    adopted = self.backend.threads[(901, 901)][0]
                    self.backend.threads[(901, 901)][0] = replace(
                        adopted, ppid=500, session=901)
            else:
                self.backend.remove_process(self.pid)
            self.returncode = 0
        return self.returncode

    def wait(self, timeout):
        deadline = time.monotonic() + timeout
        while self.poll() is None and time.monotonic() < deadline:
            time.sleep(0.001)
        if self.returncode is None:
            raise TimeoutError("fake child stayed alive")
        return self.returncode


class FakeLauncher:
    def __init__(self, backend, descendant=False, stay_alive=False,
                 fail_before_release=False, cross_session=False):
        self.backend = backend
        self.descendant = descendant
        self.stay_alive = stay_alive
        self.fail_before_release = fail_before_release
        self.cross_session = cross_session
        self.reparent_on_exit = False
        self.launch_count = 0
        self.released = False
        self.before_release = None

    def launch(self, _command, launch_cpus, _backend, ready_callback,
               release_callback, execution_nonce):
        require(isinstance(execution_nonce, str) and
                HEX256.fullmatch(execution_nonce) is not None,
                "fake launcher received an invalid execution nonce")
        self.launch_count += 1
        leader = self.backend.add(900, 900, 9000, 900, set(launch_cpus), ppid=1)
        if self.descendant:
            self.backend.add(901, 901, 9010,
                             901 if self.cross_session else 900, set(launch_cpus),
                             ppid=900, process_start=9010)
        try:
            ready_callback(leader)
            if self.before_release is not None:
                self.before_release()
            if self.fail_before_release:
                raise IsolationError("injected gate abort")
            release_callback(None)
            self.released = True
            return FakeChild(
                self.backend, leader, self.descendant, self.stay_alive,
                self.reparent_on_exit)
        except BaseException:
            self.backend.remove_process(900)
            self.backend.remove_process(901)
            raise


def make_fake_transaction(directory, backend, writer=atomic_json, launcher=None,
                          sleep_function=lambda _seconds: None, global_lock=None,
                          seal_writer=atomic_json):
    if global_lock is None:
        global_lock = FakeGlobalLock()
        # Each fake transaction models a previously acquired coordinator lock;
        # the overlap regression below exercises acquisition itself.
        global_lock.active = True
    return AffinityTransaction(
        backend, (1, 3), (0, 1, 2, 3), Path(directory) / "report.json",
        poll_ms=1, writer=writer, host=FAKE_HOST,
        sleep_function=sleep_function, launcher=launcher,
        global_lock=global_lock, seal_writer=seal_writer,
        nonce_factory=lambda: "ab" * 32)


def fake_command():
    return [sys.executable, str(Path(__file__).resolve()), "self-test-child"]


def test_global_lock_serialization():
    FakeGlobalLock.held = False
    first = FakeGlobalLock()
    second = FakeGlobalLock()
    with first:
        expect_exception(IsolationError, second.__enter__,
                         "overlapping same-UID supervisor")
        check(first.validate_current() == FAKE_LOCK_IDENTITY,
              "held global lock identity changed")
    with second:
        check(second.validate_current() == FAKE_LOCK_IDENTITY,
              "global lock was not reusable after release")
    check(not FakeGlobalLock.held,
          "fake global lock remained held after context exit")


def test_inode_bound_evidence_snapshot():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        target = root / "evidence.json"
        replacement = root / "replacement.json"
        atomic_json(target, {"value": "original"}, exclusive=True)
        atomic_json(replacement, {"value": "replacement"}, exclusive=True)

        def replace_after_read(_path):
            os.replace(replacement, target)

        expect_exception(
            IsolationError,
            lambda: stable_file_snapshot(
                target, "racing evidence", after_read=replace_after_read),
            "validate-then-reopen replacement")
        fifo = root / "evidence.fifo"
        os.mkfifo(fifo, 0o600)
        expect_exception(
            IsolationError,
            lambda: stable_file_snapshot(fifo, "FIFO evidence"),
            "non-regular evidence path")


def test_exact_restore_and_newcomer():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        first = backend.add(10, 10, 100, 5, {0, 1, 2, 3}, ppid=1)
        second = backend.add(10, 11, 101, 5, {0, 1, 2, 3},
                             ppid=1, process_start=100)
        backend.add(20, 20, 200, 6, {0, 2}, ppid=1)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.stabilize()
        check(backend.get_affinity(first) == {0, 2} and
              backend.get_affinity(second) == {0, 2},
              "baseline masks were not restricted")
        newcomer = backend.add(10, 12, 102, 5, {0, 2},
                               ppid=1, process_start=100)
        moved, inherited = transaction.restrict_once()
        check(moved == 0 and inherited == 1,
              "safe inherited newcomer was not journaled without a syscall")
        before = transaction.write_count
        check(transaction.restrict_once() == (0, 0) and
              transaction.write_count == before,
              "unchanged poll wrote or changed journal state")
        backend.add(21, 21, 210, 7, {0, 2}, ppid=20)
        check(transaction.restrict_once() == (0, 0) and
              transaction.write_count == before,
              "audit-only newcomer provenance performed a benchmark-time fsync")
        check(not transaction.restore(), "exact restore reported a failure")
        check(backend.get_affinity(first) == {0, 1, 2, 3} and
              backend.get_affinity(second) == {0, 1, 2, 3} and
              backend.get_affinity(newcomer) == {0, 1, 2, 3},
              "exact inherited masks were not restored")


def test_late_inherited_restore_reconciliation():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(10, 10, 100, 10, {0, 1, 2, 3}, ppid=1)
        backend.add(10, 11, 101, 10, {0, 1, 2, 3}, ppid=1,
                    process_start=100)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.stabilize()
        spawned = {"done": False}

        def spawn_from_still_restricted_creator(identity, cpus):
            if (not spawned["done"] and identity.tid == 11 and
                    cpus == {0, 1, 2, 3}):
                spawned["done"] = True
                backend.add(10, 12, 102, 10, {0, 2}, ppid=1,
                            process_start=100)

        backend.on_set = spawn_from_still_restricted_creator
        failures = transaction.restore_reconciled()
        check(failures == [MUTATED_RESTORE_UNPROVABLE] and
              transaction.report["uncertainty"] is True,
              "finite restoration scan was certified without containment")
        late = backend.threads[(10, 12)][0]
        check(backend.get_affinity(late) == {0, 1, 2, 3} and
              transaction.records[late.key()]["restore_status"] == "restored",
              "late inherited thread was stranded on the restricted mask")
        retained = load_report(transaction.report_path)
        check(any(item["tid"] == 12 for item in retained["records"]),
              "late inherited original was not durable before widening")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(20, 20, 200, 20, {0, 1, 2, 3}, ppid=1)
        backend.add(20, 21, 201, 20, {0, 1, 2, 3}, ppid=1,
                    process_start=200)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        spawned = {"done": False}

        def spawn_during_recovery(identity, cpus):
            if (not spawned["done"] and identity.tid == 21 and
                    cpus == {0, 1, 2, 3}):
                spawned["done"] = True
                backend.add(20, 22, 202, 20, {0, 2}, ppid=1,
                            process_start=200)

        backend.on_set = spawn_during_recovery
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, host=FAKE_HOST,
                sleep_function=lambda _seconds: None,
                global_lock=transaction.global_lock),
            "uncontained recovery completion")
        late = backend.threads[(20, 22)][0]
        check(backend.get_affinity(late) == {0, 1, 2, 3},
              "recovery stranded a late inherited thread")

    # Reproduce publication after the historical second/final empty scan.  The
    # new bounded cleanup keeps scanning and repairs this concrete late thread,
    # but still rejects the unprovable global completion claim because another
    # clone could publish after any finite last scan.
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(30, 30, 300, 30, {0, 1, 2, 3}, ppid=1)
        backend.add(30, 31, 301, 30, {0, 1, 2, 3}, ppid=1,
                    process_start=300)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.stabilize()
        scans = {"count": 0, "published": False}

        def publish_after_second_snapshot(_snapshot):
            scans["count"] += 1
            if scans["count"] == 2:
                scans["published"] = True
                backend.add(30, 32, 302, 30, {0, 2}, ppid=1,
                            process_start=300)

        backend.on_list = publish_after_second_snapshot
        failures = transaction.restore_reconciled()
        late = backend.threads[(30, 32)][0]
        check(scans["published"] and scans["count"] == MAX_STABILIZE_PASSES and
              failures == [MUTATED_RESTORE_UNPROVABLE] and
              backend.get_affinity(late) == {0, 1, 2, 3} and
              late.key() in transaction.records,
              "post-final-scan thread publication was accepted or left unrepaired")


def test_nonuniform_and_unsafe_newcomer():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2}, ppid=0)
        backend.add(1, 2, 2, 1, {0, 2}, ppid=0, process_start=1)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        expect_exception(IsolationError, transaction.initialize_provenance,
                         "nonuniform baseline process")
        check(backend.threads[(1, 1)][1] == {0, 1, 2},
              "nonuniform preflight mutated a mask")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        backend.add(2, 2, 2, 2, {0, 2}, ppid=1)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.stabilize()
        # Remove the only known parent and disable global uniform fallback by
        # retaining two original process masks in the initial provenance.
        backend.remove_process(1)
        backend.add(3, 3, 3, 3, {0, 2}, ppid=999)
        try:
            transaction.restrict_once()
            raise SelfTestError("unsafe newcomer was accepted")
        except UnsafeNewcomer as error:
            failures = transaction._cleanup_unsafe_newcomer(error)
            check(not failures and (3, 3) not in backend.threads,
                  "unsafe inherited-mask newcomer survived fail-closed cleanup")


def test_write_failure_still_restores():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        writer = CountingWriter(fail_state_once="restoring")
        launcher = FakeLauncher(backend)
        transaction = make_fake_transaction(
            directory, backend, writer=writer, launcher=launcher)
        result = transaction.run(fake_command())
        check(result == 1 and writer.failed,
              "pre-restore write failure was not reported")
        check(backend.get_affinity(identity) == {0, 1, 2, 3},
              "pre-restore write failure suppressed mask restoration")


def test_acceptance_is_path_sealed_and_fail_closed():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        supervisor = backend.add(500, 500, 5000, 500,
                                 {0, 1, 2, 3}, ppid=1)
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        transaction = make_fake_transaction(
            directory, backend, launcher=FakeLauncher(backend))
        check(transaction.run(fake_command()) == 0 and
              backend.get_affinity(supervisor) == {0, 1, 2, 3} and
              transaction.report["accepted"] is True and
              [(item["pid"], item["tid"]) for item in
               transaction.report["records"]] == [(500, 500)] and
              load_report(transaction.report_path,
                          require_accepted=True)["accepted"] is True,
              "single-threaded supervisor self-mutation was not bounded")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(
            directory, backend, launcher=FakeLauncher(backend))
        check(transaction.run(fake_command()) == 1 and
              backend.get_affinity(identity) == {0, 1, 2, 3} and
              transaction.report["accepted"] is False and
              transaction.report["uncertainty"] is True and
              MUTATED_RESTORE_UNPROVABLE in transaction.report["error"] and
              not acceptance_path(transaction.report_path).exists(),
              "uncontained affinity mutation produced authoritative evidence")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        transaction = make_fake_transaction(
            directory, backend, launcher=FakeLauncher(backend),
            seal_writer=PostInstallFailWriter())
        check(transaction.run(fake_command()) == 1,
              "post-install seal failure did not fail the transaction")
        expect_exception(
            IsolationError,
            lambda: load_report(transaction.report_path, require_accepted=True),
            "post-install seal failure acceptance")
        check(transaction.report["accepted"] is False,
              "post-install seal failure left accepted in-memory state")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        transaction = make_fake_transaction(
            directory, backend, launcher=FakeLauncher(backend))
        check(transaction.run(fake_command()) == 0 and
              load_report(transaction.report_path, require_accepted=True)["accepted"],
              "valid transaction did not create a usable acceptance seal")

        def fail_fsync(_path):
            raise OSError("injected late invalidation fsync failure")

        expect_exception(
            IsolationError,
            lambda: invalidate_acceptance_path(
                transaction.report_path, fsync_function=fail_fsync),
            "late invalidation durability failure")
        expect_exception(
            IsolationError,
            lambda: load_report(transaction.report_path, require_accepted=True),
            "late invalidation ambiguity")
        check(ambiguity_path(transaction.report_path).exists() and
              not acceptance_path(transaction.report_path).exists(),
              "late invalidation did not atomically prefer rejection")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        body_writer = CountingWriter(fail_state_once="failed")
        transaction = make_fake_transaction(
            directory, backend, writer=body_writer,
            launcher=FakeLauncher(backend),
            seal_writer=PostInstallFailWriter())

        def fail_revoke():
            raise OSError("injected acceptance revocation failure")

        transaction._invalidate_acceptance = fail_revoke
        check(transaction.run(fake_command()) == 0 and
              transaction.report["accepted"] is True and
              body_writer.failed is False and
              load_report(transaction.report_path,
                          require_accepted=True)["accepted"] is True,
              "valid installed seal survived a failed return/body rewrite")


def test_gate_signal_and_descendants():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend)
        writer = CountingWriter(fail_state_once="launch_gated")
        transaction = make_fake_transaction(
            directory, backend, writer=writer, launcher=launcher)
        check(transaction.run(fake_command()) == 1 and not launcher.released,
              "failed durable child-identity write released gated child code")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend, fail_before_release=True)
        transaction = make_fake_transaction(directory, backend, launcher=launcher)
        check(transaction.run(fake_command()) == 1 and not launcher.released,
              "parent-side gate abort released child code")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend)
        transaction = make_fake_transaction(directory, backend, launcher=launcher)
        launcher.before_release = lambda: transaction.request_signal(signal.SIGTERM)
        check(transaction.run(fake_command()) == 1 and not launcher.released and
              transaction.report["child"]["released"] is False,
              "signal after durable-ready escaped the child release gate")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend, descendant=True)
        transaction = make_fake_transaction(directory, backend, launcher=launcher)
        check(transaction.run(fake_command()) == 1,
              "child-session descendant did not invalidate evidence")
        check(not backend.session_processes(900) and
              any(item[1] == signal.SIGTERM for item in backend.signal_log),
              "child-session descendant was not boundedly cleaned")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend, descendant=True, cross_session=True)
        launcher.reparent_on_exit = True
        transaction = make_fake_transaction(directory, backend, launcher=launcher)
        check(transaction.run(fake_command()) == 1 and
              (901, 901) not in backend.threads and
              any(item[0] == 901 for item in backend.signal_log),
              "retained cross-session orphan escaped cleanup after reparenting")
        retained = load_report(transaction.report_path)
        check(any(item["pid"] == 901 for item in retained["child_processes"]),
              "cross-session descendant identity was not durably retained")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        launcher = FakeLauncher(backend)
        transaction = make_fake_transaction(directory, backend, launcher=launcher)
        transaction.request_signal(signal.SIGTERM)
        check(transaction.run(fake_command()) == 1 and launcher.launch_count == 0,
              "prelaunch signal raced through the child gate")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        launcher = FakeLauncher(backend, stay_alive=True)
        holder = {}

        def signal_after_release(_seconds):
            if launcher.released and holder["transaction"].received_signal is None:
                holder["transaction"].request_signal(signal.SIGHUP)

        transaction = make_fake_transaction(
            directory, backend, launcher=launcher,
            sleep_function=signal_after_release)
        holder["transaction"] = transaction
        check(transaction.run(fake_command()) == 128 + signal.SIGHUP,
              "running signal did not produce signal exit status")
        retained = load_report(transaction.report_path)
        check(retained["accepted"] is False and
              retained["received_signal"] == signal.SIGHUP and
              not backend.session_processes(900),
              "running signal race was accepted or left child session alive")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        launcher = FakeLauncher(
            backend, descendant=True, stay_alive=True, cross_session=True)
        holder = {}

        def inspect_cross_session_then_signal(_seconds):
            if launcher.released and holder["transaction"].received_signal is None:
                check(backend.threads[(901, 901)][1] == {0, 1, 2, 3},
                      "legitimate cross-session child was affinity-restricted")
                holder["transaction"].request_signal(signal.SIGTERM)

        transaction = make_fake_transaction(
            directory, backend, launcher=launcher,
            sleep_function=inspect_cross_session_then_signal)
        holder["transaction"] = transaction
        check(transaction.run(fake_command()) == 128 + signal.SIGTERM and
              not backend.child_scope_processes(transaction.report["child"]),
              "cross-session child tree was not excluded and boundedly cleaned")

    backend = FakeBackend()
    child = backend.add(700, 700, 7000, 700, {0, 1}, ppid=1)
    backend.stubborn.add(700)
    failures, had_members = terminate_session(
        backend, child_record_from_identity(child, True),
        grace_seconds=0.001, kill_grace_seconds=0.001,
        sleep_function=lambda _seconds: None)
    check(not failures and had_members and
          backend.signal_log == [(700, signal.SIGTERM), (700, signal.SIGKILL)],
          "bounded signal grace did not escalate stubborn child to SIGKILL")


def test_descriptor_pinned_exec_and_child_signal_defaults():
    operations = []
    handlers = {signum: "copied-python-handler" for signum in GATED_SIGNALS}
    outcome = {"terminated": False, "executed": False}

    def fake_signal(signum, disposition):
        handlers[signum] = disposition
        operations.append(("handler", int(signum), disposition))

    def fake_mask(how, mask):
        operations.append(("mask", how, frozenset(mask)))
        # Model a SIGTERM delivered while the child was blocked in gate_read.
        # Unblocking with a copied Python handler would return to the read and
        # later execute; SIG_DFL terminates before any release is possible.
        if handlers[signal.SIGTERM] == signal.SIG_DFL:
            outcome["terminated"] = True
        else:
            outcome["executed"] = True

    reset_gated_child_signals(
        set(), signal_function=fake_signal, mask_function=fake_mask)
    check(outcome == {"terminated": True, "executed": False} and
          [item[0] for item in operations] ==
          ["handler", "handler", "handler", "mask"],
          "gated child unblocked before replacing copied signal handlers")

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        executable = root / "command"
        script = root / "runner.py"
        executable.write_bytes(b"#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)
        script.write_bytes(b"print('original')\n")
        pinned = PinnedCommand([str(executable), str(script), "argument"])
        try:
            original_identity = json.loads(json.dumps(pinned.identity))
            original_command_hash = sha256_value(pinned.command)
            saved_executable = root / "command.original"
            saved_script = root / "runner.original.py"
            os.replace(executable, saved_executable)
            os.replace(script, saved_script)
            executable.write_bytes(b"#!/bin/sh\nexit 99\n")
            executable.chmod(0o755)
            script.write_bytes(b"print('replacement')\n")
            captured = {}

            class CapturedExec(Exception):
                pass

            def fake_fchdir(descriptor):
                metadata = os.fstat(descriptor)
                captured["cwd"] = (metadata.st_dev, metadata.st_ino)

            def fake_set_affinity(pid, cpus):
                captured["set_affinity"] = (pid, set(cpus))

            def fake_get_affinity(pid):
                check(pid == 0, "pinned exec queried another PID")
                return {0, 1, 2, 3}

            def fake_execve(path, argv, environment):
                captured["exec_path"] = path
                captured["argv"] = list(argv)
                captured["executable_bytes"] = Path(path).read_bytes()
                captured["script_bytes"] = Path(argv[1]).read_bytes()
                captured["nonce"] = environment.get(EXECUTION_NONCE_ENV)
                raise CapturedExec()

            expect_exception(
                CapturedExec,
                lambda: execute_pinned_command(
                    pinned, {0, 1, 2, 3}, "cd" * 32,
                    fchdir_function=fake_fchdir,
                    set_affinity_function=fake_set_affinity,
                    get_affinity_function=fake_get_affinity,
                    execve_function=fake_execve, environment={}),
                "descriptor-pinned execution capture")

            # Model replacement code restoring the original pathname before
            # final evidence inspection.  Path-only pre/post hashing would miss
            # this attack; the captured descriptors never selected replacement.
            executable.unlink()
            script.unlink()
            os.replace(saved_executable, executable)
            os.replace(saved_script, script)
            final_executable = bounded_file_identity(str(executable))
            final_script = bounded_file_identity(str(script))
            check(captured["executable_bytes"] == b"#!/bin/sh\nexit 0\n" and
                  captured["script_bytes"] == b"print('original')\n" and
                  captured["argv"] == [str(executable),
                                        pinned._proc_fd(
                                            pinned.script_descriptor),
                                        "argument"] and
                  captured["exec_path"] == pinned._proc_fd(
                      pinned.executable_descriptor) and
                  captured["cwd"] == (
                      original_identity["cwd"]["device"],
                      original_identity["cwd"]["inode"]) and
                  captured["set_affinity"] == (0, {0, 1, 2, 3}) and
                  captured["nonce"] == "cd" * 32 and
                  original_command_hash == sha256_value(pinned.command) and
                  (final_executable["device"], final_executable["inode"],
                   final_executable["sha256"]) ==
                  (original_identity["executable"]["device"],
                   original_identity["executable"]["inode"],
                   original_identity["executable"]["sha256"]) and
                  (final_script["device"], final_script["inode"],
                   final_script["sha256"]) ==
                  (original_identity["script"]["device"],
                   original_identity["script"]["inode"],
                   original_identity["script"]["sha256"]),
                  "pathname replacement selected unpinned command code")
        finally:
            pinned.close()


def test_tid_reuse_and_recovery():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        backend.reuse_after_set.add((1, 1))
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        expect_exception(ThreadMutationUncertain,
                         transaction.initialize_provenance,
                         "post-syscall TID reuse")
        check(transaction.report["uncertainty"] is True,
              "post-syscall TID reuse did not mark uncertainty")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        # Model PID/TID reuse inside the same clock tick.  Tick-only identity
        # would accept this replacement; procfs directory inodes must not.
        backend.threads[(1, 1)][0] = replace(
            identity, task_inode=identity.task_inode + 1,
            process_inode=identity.process_inode + 1)
        failures = transaction.restore()
        record = transaction.records[identity.key()]
        check(failures and record["restore_status"] == "reused" and
              transaction.report["uncertainty"] is True,
              "same-tick procfs identity reuse was not rejected")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        check(backend.get_affinity(identity) == {0, 2},
              "recovery fixture was not restricted")
        transaction.records[identity.key()]["restore_status"] = "failed"
        transaction._write()
        wrong_host = {
            "boot_id": "00000000-0000-0000-0000-000000000002",
            "pid_namespace": dict(FAKE_HOST["pid_namespace"]),
        }
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, host=wrong_host,
                global_lock=transaction.global_lock),
            "cross-boot recovery")
        check(backend.get_affinity(identity) == {0, 2},
              "cross-boot recovery mutated a journaled mask")
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, host=FAKE_HOST,
                sleep_function=lambda _seconds: None,
                global_lock=transaction.global_lock),
            "uncontained recovery completion")
        check(backend.get_affinity(identity) == {0, 1, 2, 3},
              "best-effort recovery did not restore the exact original mask")

        malformed = load_json(transaction.report_path, "test report")
        malformed["records"][0]["restore_status"] = "mystery"
        atomic_json(transaction.report_path, malformed)
        expect_exception(IsolationError,
                         lambda: recover_report(
                             transaction.report_path, backend=backend, host=FAKE_HOST,
                             global_lock=transaction.global_lock),
                         "unknown recovery status")

    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        failing_writer = CountingWriter(fail_state_once="recovering")
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, writer=failing_writer,
                host=FAKE_HOST, sleep_function=lambda _seconds: None,
                global_lock=transaction.global_lock),
            "pre-recovery journal write failure")
        check(backend.get_affinity(identity) == {0, 1, 2, 3},
              "pre-recovery write failure suppressed exact restore")


def test_recovery_kills_session_before_restore():
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        child = backend.add(900, 900, 9000, 900, {0, 1, 2, 3}, ppid=1)
        descendant = backend.add(
            901, 901, 9010, 901, {0, 1, 2, 3}, ppid=900)
        transaction.report["child"] = child_record_from_identity(child, True)
        for member in (child, descendant):
            process = child_process_record(member)
            key = child_process_key(process)
            transaction.child_processes[key] = process
            transaction.child_process_order.append(key)
        transaction.report["state"] = "running"
        transaction._write()
        backend.remove_process(900)
        backend.threads[(901, 901)][0] = replace(
            descendant, ppid=1, session=901)
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, host=FAKE_HOST,
                sleep_function=lambda _seconds: None,
                global_lock=transaction.global_lock),
            "released-child crash recovery completion")
        check((901, 901) not in backend.threads and
              backend.get_affinity(identity) == {0, 1, 2, 3},
              "recovery restored masks before emptying retained descendants")
        signal_index = next(index for index, item in enumerate(backend.operation_log)
                            if item[0] == "signal")
        restore_index = next(
            index for index, item in enumerate(backend.operation_log)
            if item[0] == "set" and item[2] == (0, 1, 2, 3))
        check(signal_index < restore_index,
              "recovery ordering did not include child cleanup")

    # Reproduce a crash before the original subreaper observed a double-forked
    # child.  Once it has changed session and reparented outside the replacement
    # supervisor, PID/session/ancestry scans cannot rediscover it.  Recovery may
    # repair known masks, but must return failure instead of claiming the scope
    # was emptied.
    with tempfile.TemporaryDirectory() as directory:
        backend = FakeBackend()
        identity = backend.add(1, 1, 1, 1, {0, 1, 2, 3}, ppid=0)
        transaction = make_fake_transaction(directory, backend)
        transaction.prepare_command(fake_command())
        transaction.initialize_provenance()
        child = backend.add(900, 900, 9000, 900, {0, 1, 2, 3}, ppid=1)
        transaction.report["child"] = child_record_from_identity(child, True)
        process = child_process_record(child)
        key = child_process_key(process)
        transaction.child_processes[key] = process
        transaction.child_process_order.append(key)
        transaction.report["state"] = "running"
        transaction._write()
        backend.remove_process(900)
        backend.add(901, 901, 9010, 901, {0, 1, 2, 3}, ppid=1)
        expect_exception(
            IsolationError,
            lambda: recover_report(
                transaction.report_path, backend=backend, host=FAKE_HOST,
                sleep_function=lambda _seconds: None,
                global_lock=transaction.global_lock),
            "unobserved cross-session orphan recovery")
        retained = load_report(transaction.report_path)
        check((901, 901) in backend.threads and not backend.signal_log and
              backend.get_affinity(identity) == {0, 1, 2, 3} and
              retained["state"] == "recovery_failed" and
              retained["uncertainty"] is True and
              CRASH_SCOPE_UNPROVABLE in retained["error"],
              "crash-unobserved orphan was hidden by successful recovery")


def test_binding():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        evidence = root / "evidence"
        evidence.mkdir()
        backend = FakeBackend()
        # Authoritative v6 evidence is accepted only when the surrounding work
        # was already outside the reserved pair and no mask needed mutation.
        backend.add(1, 1, 1, 1, {0, 2}, ppid=0)
        report_path = root / "report.json"
        reservation_path = root / "reservation.json"
        reservation_payload = {"nonce": "reservation-nonce"}
        atomic_json(reservation_path, reservation_payload, exclusive=True)
        launcher = FakeLauncher(backend)
        global_lock = FakeGlobalLock()
        global_lock.active = True
        transaction = AffinityTransaction(
            backend, (1, 3), (0, 1, 2, 3), report_path, poll_ms=1,
            host=FAKE_HOST, sleep_function=lambda _seconds: None,
            launcher=launcher, global_lock=global_lock)
        artifact = str(Path(__file__).resolve())
        command = [
            sys.executable, artifact, "run",
            "--baseline", artifact, "--candidate", artifact,
            "--baseline-archive", artifact, "--candidate-archive", artifact,
            "--baseline-build-dir", str(root),
            "--candidate-build-dir", str(root),
            "--baseline-source-root", str(root),
            "--candidate-source-root", str(root),
            "--candidate-commit", "0" * 40,
            "--candidate-mode", "auto",
            "--reservation-file", str(root / "reservation.json"),
            "--output", str(evidence), "--cpu", "1",
            "--reserved-sibling", "3", "--preset", "smoke",
            "--reuse", "8", "--iterations", "9", "--warmup", "2",
            "--timeout", "120.0",
        ]
        check(transaction.run(command) == 0 and transaction.report["accepted"] is True,
              "binding fixture transaction was not accepted")
        script_identity = transaction.report["command_identity"]["script"]
        raw_payload = {
            "schema": MAIN_RAW_SCHEMA,
            "campaign": {
                "benchmark_cpu": 1, "reserved_sibling": 3, "reuse": 8,
                "iterations": 9, "warmup": 2, "timeout_seconds": 120.0,
                "candidate_mode": "auto",
                "allowed_cpu_set_at_launch": [0, 1, 2, 3],
            },
            "input_specification": {
                "runner": script_identity["path"],
                "baseline_executable": artifact,
                "candidate_executable": artifact,
                "baseline_archive": artifact,
                "candidate_archive": artifact,
                "baseline_build_dir": str(root),
                "candidate_build_dir": str(root),
                "baseline_source_root": str(root),
                "candidate_source_root": str(root),
                "candidate_commit": "0" * 40,
                "taskset": str(Path("/usr/bin/taskset").resolve()),
                "ldd": str(Path("/usr/bin/ldd").resolve()),
            },
            "identities_initial": {
                "runner": {"path": script_identity["path"],
                           "sha256": script_identity["sha256"]},
            },
            "reservation": {
                "path": str(reservation_path),
                "payload": reservation_payload,
            },
        }
        raw_payload["reservation"]["sha256"] = sha256_bytes(canonical_bytes(
            raw_payload["reservation"]["payload"]))
        timestamp = transaction.report["execution"]["child_ready_monotonic_ns"]
        raw_payload["isolation"] = {
            "before": {"monotonic_ns": timestamp},
            "after": {"monotonic_ns": timestamp},
        }
        raw_payload["supervision"] = {
            "schema": MAIN_SUPERVISION_SCHEMA,
            "execution_nonce": transaction.report["execution"]["nonce"],
            "runner_pid": 900,
            "runner_started_monotonic_ns": timestamp,
            "runner_finished_monotonic_ns": timestamp,
            "launch_cpus": [0, 1, 2, 3],
            "reserved_cpus": [1, 3],
            "campaign_sha256": sha256_bytes(canonical_bytes(
                raw_payload["campaign"])),
            "reservation_sha256": raw_payload["reservation"]["sha256"],
            "reservation_nonce": "reservation-nonce",
            "isolation_before_monotonic_ns": timestamp,
            "isolation_after_monotonic_ns": timestamp,
        }
        raw_payload["digest"] = sha256_value(raw_payload)
        raw_path = evidence / "raw.json"
        atomic_json(raw_path, raw_payload, exclusive=True)
        raw_bytes = raw_path.read_bytes()
        manifest_payload = {
            "schema": MAIN_MANIFEST_SCHEMA, "valid": True,
            "supervision": raw_payload["supervision"],
            "raw": {
                "path": "raw.json", "size": len(raw_bytes),
                "sha256": sha256_bytes(raw_bytes),
                "payload_digest": raw_payload["digest"],
            },
        }
        manifest_payload["digest"] = sha256_value(manifest_payload)
        manifest_path = evidence / "manifest.json"
        atomic_json(manifest_path, manifest_payload, exclusive=True)
        binding_path = evidence / "affinity-binding.json"
        create_binding(report_path, manifest_path, binding_path)
        binding = load_json(binding_path, "test binding")
        validate_binding(
            binding, binding_path, manifest_path,
            sha256_bytes(manifest_path.read_bytes()))
        expect_exception(
            IsolationError,
            lambda: validate_binding(
                binding, binding_path, manifest_path, "f" * 64),
            "binding against another manifest snapshot")

        def install_bundle(raw_value):
            raw_value = json.loads(json.dumps(raw_value))
            raw_value.pop("digest", None)
            raw_value["digest"] = sha256_value(raw_value)
            atomic_json(raw_path, raw_value)
            raw_bytes_value = raw_path.read_bytes()
            manifest_value = json.loads(json.dumps(manifest_payload))
            manifest_value["supervision"] = raw_value["supervision"]
            manifest_value["raw"] = {
                "path": "raw.json", "size": len(raw_bytes_value),
                "sha256": sha256_bytes(raw_bytes_value),
                "payload_digest": raw_value["digest"],
            }
            manifest_value.pop("digest", None)
            manifest_value["digest"] = sha256_value(manifest_value)
            atomic_json(manifest_path, manifest_value)

        mutations = []
        changed_nonce = json.loads(json.dumps(raw_payload))
        changed_nonce["supervision"]["execution_nonce"] = "cd" * 32
        mutations.append(("execution nonce", changed_nonce))
        changed_interval = json.loads(json.dumps(raw_payload))
        changed_interval["supervision"]["runner_started_monotonic_ns"] = timestamp - 1
        mutations.append(("execution interval", changed_interval))
        changed_pid = json.loads(json.dumps(raw_payload))
        changed_pid["supervision"]["runner_pid"] = 901
        mutations.append(("runner PID", changed_pid))
        changed_launch = json.loads(json.dumps(raw_payload))
        changed_launch["supervision"]["launch_cpus"] = [0, 1, 2]
        mutations.append(("launch CPU set", changed_launch))
        changed_campaign = json.loads(json.dumps(raw_payload))
        changed_campaign["campaign"]["reuse"] = 9
        changed_campaign["supervision"]["campaign_sha256"] = sha256_bytes(
            canonical_bytes(changed_campaign["campaign"]))
        mutations.append(("campaign payload", changed_campaign))
        changed_reservation = json.loads(json.dumps(raw_payload))
        changed_reservation["reservation"]["payload"]["nonce"] = "other-reservation"
        changed_reservation["reservation"]["sha256"] = sha256_bytes(canonical_bytes(
            changed_reservation["reservation"]["payload"]))
        changed_reservation["supervision"]["reservation_nonce"] = "other-reservation"
        changed_reservation["supervision"]["reservation_sha256"] = \
            changed_reservation["reservation"]["sha256"]
        mutations.append(("reservation payload", changed_reservation))
        for label, mutation in mutations:
            install_bundle(mutation)
            expect_exception(
                IsolationError,
                lambda path=manifest_path: validate_main_manifest_binding(
                    transaction.report, path),
                "changed {} binding".format(label))
        install_bundle(raw_payload)

        changed = load_json(report_path, "test report")
        changed["command_sha256"] = "f" * 64
        atomic_json(report_path, changed)
        expect_exception(IsolationError,
                         lambda: validate_binding(
                             load_json(binding_path, "test binding"), binding_path),
                         "binding command tamper")


def self_test():
    tests = (
        test_global_lock_serialization,
        test_inode_bound_evidence_snapshot,
        test_exact_restore_and_newcomer,
        test_late_inherited_restore_reconciliation,
        test_nonuniform_and_unsafe_newcomer,
        test_write_failure_still_restores,
        test_acceptance_is_path_sealed_and_fail_closed,
        test_gate_signal_and_descendants,
        test_descriptor_pinned_exec_and_child_signal_defaults,
        test_tid_reuse_and_recovery,
        test_recovery_kills_session_before_restore,
        test_binding,
    )
    for test in tests:
        test()
    check(parse_cpu_list("0-2,4,6-7") == [0, 1, 2, 4, 6, 7],
          "CPU list parser changed")
    check(format_cpu_list([0, 1, 2, 4, 6, 7]) == "0-2,4,6-7",
          "CPU list formatter changed")
    print("Leopard2 affinity supervisor self-test passed ({} cases)".format(
        len(tests)))
    return 0


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="operation", required=True)
    run = commands.add_parser("run", help="isolate a pair and run one command")
    run.add_argument("--report", required=True, type=Path)
    run.add_argument("--reserved-cpus", required=True,
                     help="logical CPU list, normally one SMT pair")
    run.add_argument("--poll-ms", type=int, default=DEFAULT_POLL_MS)
    run.add_argument("child", nargs=argparse.REMAINDER)
    restore = commands.add_parser(
        "restore", help="restore a journal left by an unclean exit")
    restore.add_argument("--report", required=True, type=Path)
    verify = commands.add_parser(
        "verify-report", help="validate an accepted affinity report")
    verify.add_argument("--report", required=True, type=Path)
    bind = commands.add_parser(
        "bind", help="bind an accepted report to a main-comparison manifest")
    bind.add_argument("--report", required=True, type=Path)
    bind.add_argument("--manifest", required=True, type=Path)
    bind.add_argument("--output", required=True, type=Path)
    verify_binding_parser = commands.add_parser(
        "verify-binding", help="verify an affinity/main evidence binding")
    verify_binding_parser.add_argument("--binding", required=True, type=Path)
    verify_binding_parser.add_argument("--manifest", type=Path)
    verify_binding_parser.add_argument("--manifest-sha256")
    commands.add_parser(
        "self-test", help="run deterministic fake tests without changing affinity")
    return result


def main(arguments=None):
    options = parser().parse_args(arguments)
    try:
        if options.operation == "self-test":
            return self_test()
        if options.operation == "restore":
            with GlobalSupervisorLock() as global_lock:
                return recover_report(options.report, global_lock=global_lock)
        if options.operation == "verify-report":
            load_report(options.report, require_accepted=True)
            print("Leopard2 affinity supervisor report verified")
            return 0
        if options.operation == "bind":
            return create_binding(options.report, options.manifest, options.output)
        if options.operation == "verify-binding":
            require((options.manifest is None) ==
                    (options.manifest_sha256 is None),
                    "manifest path and hash must be supplied together")
            binding_path, binding, _bytes, _identity = stable_json_snapshot(
                options.binding, "affinity/main binding")
            validate_binding(
                binding, binding_path, options.manifest,
                options.manifest_sha256)
            print("Leopard2 affinity/main evidence binding verified")
            return 0
        if sys.platform != "linux" or not hasattr(os, "sched_getaffinity"):
            raise IsolationError("the affinity supervisor requires Linux sched affinity")
        child = list(options.child)
        if child and child[0] == "--":
            child = child[1:]
        require(bool(child), "run requires a child command after --")
        with GlobalSupervisorLock() as global_lock:
            launch = sorted(os.sched_getaffinity(0))
            transaction = AffinityTransaction(
                LinuxThreadBackend(), parse_cpu_list(options.reserved_cpus), launch,
                options.report, options.poll_ms, global_lock=global_lock)
            previous_handlers = {}
            for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
                previous_handlers[signum] = signal.getsignal(signum)

                def supervised_signal(received, frame, tx=transaction,
                                      previous=previous_handlers[signum]):
                    if tx.request_signal(received):
                        return
                    # The transaction's explicit signal boundary has closed.
                    # Preserve the caller's disposition instead of swallowing a
                    # post-transaction signal in the temporary handler.
                    if previous == signal.SIG_IGN:
                        return
                    if previous == signal.SIG_DFL:
                        signal.signal(received, signal.SIG_DFL)
                        os.kill(os.getpid(), received)
                        return
                    previous(received, frame)

                signal.signal(signum, supervised_signal)
            result_code = 1
            late_mask = None
            try:
                result_code = transaction.run(child)
            finally:
                # Block while restoring the caller's dispositions.  Signals after
                # run()'s explicit boundary are delivered under those original
                # dispositions, never folded back into committed evidence.
                late_mask = signal.pthread_sigmask(
                    signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM, signal.SIGHUP})
                try:
                    for signum, handler in previous_handlers.items():
                        signal.signal(signum, handler)
                finally:
                    signal.pthread_sigmask(signal.SIG_SETMASK, late_mask)
            return result_code
    except (IsolationError, OSError, ValueError) as error:
        print("affinity supervisor error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
