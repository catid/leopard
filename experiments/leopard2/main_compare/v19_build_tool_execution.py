#!/usr/bin/python3
"""Retain non-GCC build roots; leopard-79h.38.5.4.8.2.2.2.3.

The caller owns host/resources, source/recipe inputs and later consumers.
This component does not route nested tools, scripts, plugins or build data.
"""
from __future__ import annotations
from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import stat
import threading

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_compiler_execution.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("build-tool execution dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_build_tool_compiler", dependency)
compiler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compiler)
builder, provenance, require = compiler.builder, compiler.provenance, compiler.require
ROLES = {"git": ("benchmark_git", "/usr/bin/git"), "cmake": (None, "/usr/bin/cmake"),
         "make": ("make_program", "/usr/bin/make"), "shell": (None, "/usr/bin/sh"),
         "ar": ("archiver", "/usr/bin/ar"), "ranlib": ("ranlib", "/usr/bin/ranlib")}
OBSERVED_PATHS = {"cmake": "/usr/bin/cmake", "shell": "/usr/bin/dash"}
MAX_TOOL_BYTES, MAX_PHASE_BYTES = 16 << 20, 24 << 20


def tool_inventory(pinned, observed):
    require(type(pinned) is dict and type(observed) is dict and set(observed) == set(OBSERVED_PATHS),
            "build tool supplemental inventory differs")
    result, paths, total = {}, set(), 0
    for role, (key, logical) in ROLES.items():
        row = pinned.get(key) if key is not None else observed[role]
        require(type(row) is dict and type(row.get("path")) is str and
                type(row.get("size")) is int and 0 < row["size"] <= MAX_TOOL_BYTES and
                type(row.get("sha256")) is str and builder.re.fullmatch(r"[0-9a-f]{64}", row["sha256"]),
                "build tool pin is invalid")
        builder.host.canonical_path(row["path"], nonroot=True)
        if key is not None:
            require(type(row.get("uid")) is int and type(row.get("gid")) is int and row["uid"] == row["gid"] == 0 and
                    type(row.get("mode")) is int and row["mode"] == (stat.S_IFREG | 0o755), "pinned build tool identity differs")
        else:
            require(set(row) == {"path", "size", "sha256"} and row["path"] == OBSERVED_PATHS[role],
                    "observed build tool role/path differs")
        require(row["path"] not in paths, "build tool paths alias")
        paths.add(row["path"])
        total += row["size"]
        result[role] = copy.deepcopy(row)
    require(total <= MAX_PHASE_BYTES, "build tool phase exceeds byte bound")
    return result


class BuildToolExecution:
    """Six retained roots, with original pin authority explicitly distinguished.

    CMake and dash are newly observed caller pins. Git, Make, ar and ranlib
    retain their original preflight authority. /usr/bin/sh is qualified here;
    Make's /bin/sh directory-alias route is a separate integration obligation.
    """
    def __init__(self, retained_preflight, parent, observed_pins, *, _tool_factory=None):
        require(type(retained_preflight) is builder.preflight.PinnedPreflight or _tool_factory is not None,
                "build tools require a live pinned preflight owner")
        self.retained = retained_preflight
        self.parent = builder.canonical_root(parent)
        self.observed = copy.deepcopy(observed_pins)
        self._factory = builder._StreamedTool if _tool_factory is None else _tool_factory
        self._stack, self._state, self._pid = ExitStack(), "new", os.getpid()
        self._tools, self._aliases = {}, {}
        self.logical_paths = {role: logical for role, (_, logical) in ROLES.items()}

    def __enter__(self):
        require(self._state == "new", "build tool phase cannot be reused")
        self._state = "entering"
        try:
            require(os.getpid() == self._pid and threading.active_count() == 1, "build tool owner process changed")
            self.retained.validate_current()
            pinned = builder.preflight._json(self.retained._bytes("candidate-build-provenance.json"))
            self.pins = tool_inventory(pinned, self.observed)
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 build tool parent"))
            self._guard.add_directory_path(self.parent)
            self._guard._add_watch(self.parent, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, {b""})
            self._parent_fd = builder.host.LinuxReader.open_directory(str(self.parent))
            self._stack.callback(os.close, self._parent_fd)
            value = os.fstat(self._parent_fd)
            require(value.st_uid == os.geteuid() and value.st_gid == os.getegid() and stat.S_IMODE(value.st_mode) == 0o700,
                    "build tool parent is not private")
            self._parent_fields = builder.streamed._directory_identity(value)
            for role, pin in self.pins.items():
                tool = self._factory(pin["path"], maximum_bytes=MAX_TOOL_BYTES)
                self._stack.callback(tool.close)
                require(tool.sha256 == pin["sha256"] and os.fstat(tool.fd).st_size == pin["size"],
                        "build tool differs from retained pin")
                self._tools[role] = tool
                self._aliases[role] = self._stack.enter_context(compiler._retain_driver_aliases(self.logical_paths[role], tool))
                tool.executable_descriptor
            inodes = {(os.fstat(tool.fd).st_dev, os.fstat(tool.fd).st_ino) for tool in self._tools.values()}
            require(len(inodes) == len(self._tools), "build tool inodes alias")
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def validate_current(self):
        require(self._state == "held", "build tool phase is not live")
        try:
            require(os.getpid() == self._pid and threading.active_count() == 1, "build tool owner process changed")
            self.retained.validate_current()
            self._guard.verify()
            require(not os.get_inheritable(self._parent_fd) and self.parent.resolve(strict=True) == self.parent and
                    builder.streamed._directory_identity(os.fstat(self._parent_fd)) == self._parent_fields ==
                    builder.streamed._directory_identity(self.parent.lstat()), "build tool parent changed")
            for aliases in self._aliases.values(): aliases.validate_current()
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def arguments(self, role, argv):
        self.validate_current()
        try:
            require(type(role) is str and role in self._tools and type(argv) is list and 0 < len(argv) <= 512 and
                    all(type(value) is str and "\0" not in value for value in argv) and argv[0] == self.logical_paths[role],
                    "build tool role or argument list differs")
            require(not any(value.startswith("@") for value in argv[1:]), "build tool response files are not qualified")
            return list(argv)
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-build-tool-execution/v1", "logical_paths": self.logical_paths,
            "tools": {role: {"path": self.pins[role]["path"], **tool.executable_record()} for role, tool in self._tools.items()},
            "aliases": {role: aliases.record() for role, aliases in self._aliases.items()},
            "original_preflight_roles": [role for role, (key, _) in ROLES.items() if key is not None],
            "newly_observed_roles": list(OBSERVED_PATHS), "declared_build_roots_sealed": True,
            "nested_tool_execution_owned": False, "build_tool_data_owned": False,
            "full_runtime_execution_owned": False, "fresh_build_recipe_integrated": False,
            "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise builder.host.PreflightError("failed build tools cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
