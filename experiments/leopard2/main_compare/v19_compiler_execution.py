#!/usr/bin/python3
"""Phase-owned GCC execution; leopard-79h.38.5.4.8.2.2.2.3.

Explicit component API only: no acquisition or fresh-build dispatch. The caller
must own host/resources, sources and its validated recipe. Shared libraries,
the ELF interpreter, compiler data and other build tools are not owned here.
"""
from __future__ import annotations

from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import secrets
import stat
import threading

HERE = Path(__file__).resolve().parent
_path = HERE / "v19_fresh_build.py"
if _path.resolve(strict=True) != _path:
    raise RuntimeError("compiler execution dependency is not canonical")
_spec = importlib.util.spec_from_file_location("v19_compiler_builder", _path)
builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(builder)
host, provenance, require = builder.host, builder.provenance, builder.require
MAX_PHASE_BYTES = 48 << 20
LANGUAGES = {"c": ("c_compiler", "/usr/bin/cc", "cc1"),
             "c++": ("compiler", "/usr/bin/c++", "cc1plus")}


def phase_inventory(pinned, language):
    require(type(language) is str and language in LANGUAGES, "unsupported GCC phase language")
    driver_key, logical, frontend = LANGUAGES[language]
    require(type(pinned) is dict and type(pinned.get("compiler_subtools")) is list and
            len(pinned["compiler_subtools"]) == 8, "pinned GCC subtool inventory differs")
    roles = {}
    for row in pinned["compiler_subtools"]:
        require(type(row) is dict and set(row) == {"language", "role", "identity"} and
                type(row["language"]) is str and row["language"] in LANGUAGES and
                type(row["role"]) is str and
                row["role"] in (LANGUAGES[row["language"]][2], "as", "collect2", "ld"),
                "pinned GCC subtool role differs")
        key = row["language"], row["role"]
        require(key not in roles, "duplicate pinned GCC subtool role")
        roles[key] = row["identity"]
    require(set(roles) == {(lang, role) for lang in LANGUAGES
                          for role in (LANGUAGES[lang][2], "as", "collect2", "ld")},
            "pinned GCC subtool coverage differs")
    result = {"driver": pinned.get(driver_key), **{role: roles[(language, role)]
              for role in (frontend, "as", "collect2", "ld")}}
    paths, total = set(), 0
    for role, row in result.items():
        require(type(row) is dict and type(row.get("path")) is str and
                type(row.get("size")) is int and 0 < row["size"] <= 64 << 20 and
                type(row.get("sha256")) is str and builder.re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) and
                type(row.get("uid")) is int and type(row.get("gid")) is int and
                row["uid"] == row["gid"] == 0 and
                type(row.get("mode")) is int and row["mode"] == (stat.S_IFREG | 0o755),
                "pinned compiler tool identity differs")
        host.canonical_path(row["path"], nonroot=True)
        require(row["path"] not in paths, "compiler phase tools alias one pathname")
        paths.add(row["path"])
        total += row["size"]
    require(total <= MAX_PHASE_BYTES, "compiler phase exceeds memory bound")
    return logical, copy.deepcopy(result)


def driver_arguments(argv, logical_driver, prefix_descriptor):
    require(type(argv) is list and argv and argv[0] == logical_driver and
            all(type(item) is str and "\0" not in item for item in argv), "compiler phase argv differs")
    require(type(prefix_descriptor) is int and prefix_descriptor >= 3, "compiler helper prefix descriptor differs")
    # The caller still owns the complete recipe and compiler inputs.
    require(not any(item.startswith(("@", "-B", "-specs", "--specs", "-wrapper", "-fplugin",
                                     "-fuse-ld", "-flto", "-x")) for item in argv[1:]),
            "compiler phase argv overrides owned helper selection")
    return [argv[0], f"-B/proc/self/fd/{prefix_descriptor}/", *argv[1:]]


class CompilerExecution:
    """Own exactly one GCC language phase and preserve failed scratch roots.

    The retained preflight supplies authority, never a caller-provided freeform
    tool hash list. A new helper directory is held through child execution;
    compiler argv gains only the explicit -B prefix. This component does not
    assert those injected arguments preserve a particular output: the caller
    must independently validate that before using a changed build recipe.
    """
    def __init__(self, retained_preflight, parent: Path, language: str, *, _tool_factory=None):
        require(type(retained_preflight) is builder.preflight.PinnedPreflight or _tool_factory is not None,
                "compiler execution requires a live pinned preflight owner")
        self.retained = retained_preflight
        self.parent = builder.canonical_root(parent)
        self.language = language
        self._factory = builder._StreamedTool if _tool_factory is None else _tool_factory
        self._stack = ExitStack()
        self._state = "new"
        self._pid = os.getpid()
        self._tools = {}
        self._commands = []
        self.root = None

    def __enter__(self):
        require(self._state == "new", "compiler phase cannot be reused")
        self._state = "entering"
        try:
            require(os.getpid() == self._pid and threading.active_count() == 1, "compiler owner process changed")
            self.retained.validate_current()
            pinned = builder.preflight._json(self.retained._bytes("candidate-build-provenance.json"))
            self.logical_driver, self.pins = phase_inventory(pinned, self.language)
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 compiler phase parent"))
            self._guard.add_directory_path(self.parent)
            self._guard._add_watch(self.parent, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, {b""})
            self._parent_fd = host.LinuxReader.open_directory(str(self.parent))
            self._stack.callback(os.close, self._parent_fd)
            value = os.fstat(self._parent_fd)
            require(value.st_uid == os.geteuid() and value.st_gid == os.getegid() and
                    stat.S_IMODE(value.st_mode) == 0o700, "compiler phase parent is not private")
            self._parent_fields = builder.streamed._directory_identity(value)
            self._validate_parent()
            # All creation is relative to retained descriptors. A renamed
            # parent cannot redirect our writes through a replacement symlink.
            name = "v19-gcc-" + secrets.token_hex(16)
            os.mkdir(name, mode=0o700, dir_fd=self._parent_fd)
            self.root = self.parent / name
            self._root_fd = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
                                    dir_fd=self._parent_fd)
            self._stack.callback(os.close, self._root_fd)
            self._validate_parent()
            value = os.fstat(self._root_fd)
            require(value.st_uid == os.geteuid() and value.st_gid == os.getegid() and
                    stat.S_IMODE(value.st_mode) == 0o700, "compiler prefix is not private")
            for role, pin in self.pins.items():
                tool = self._factory(pin["path"], maximum_bytes=64 << 20)
                self._stack.callback(tool.close)
                require(tool.sha256 == pin["sha256"] and os.fstat(tool.fd).st_size == pin["size"],
                        "compiler tool differs from pinned preflight")
                self._tools[role] = tool
            inodes = [(os.fstat(tool.fd).st_dev, os.fstat(tool.fd).st_ino) for tool in self._tools.values()]
            require(len(set(inodes)) == len(inodes), "compiler phase tool inodes alias")
            mappings = {}
            for role, tool in self._tools.items():
                descriptor = tool.executable_descriptor
                if role != "driver":
                    mappings[role] = f"/proc/self/fd/{descriptor}"
                    os.symlink(mappings[role], role, dir_fd=self._root_fd)
            os.fchmod(self._root_fd, 0o500)
            self._guard.verify()
            self._guard.add_directory_path(self.root)
            self._guard._add_watch(self.root, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, None)
            self.prefix = self._stack.enter_context(provenance._retain_exact_symlink_directory(
                self.root, mappings, "v19 GCC immutable helper prefix"))
            require(provenance._stable_fields(os.fstat(self.prefix.descriptor)) ==
                    provenance._stable_fields(os.fstat(self._root_fd)), "compiler prefix was replaced")
            self._mappings = mappings
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def _validate_parent(self):
        self._guard.verify()
        require(not os.get_inheritable(self._parent_fd) and
                builder.streamed._directory_identity(os.fstat(self._parent_fd)) == self._parent_fields ==
                builder.streamed._directory_identity(self.parent.lstat()), "compiler phase parent changed")
        self._guard.verify()

    def validate_current(self):
        require(self._state == "held", "compiler phase is not live")
        try:
            require(os.getpid() == self._pid and threading.active_count() == 1, "compiler owner process changed")
            self.retained.validate_current()
            self._validate_parent()
            self.prefix.verify()
            require(not os.get_inheritable(self.prefix.descriptor) and not os.get_inheritable(self._root_fd) and
                    provenance._stable_fields(os.fstat(self._root_fd)) ==
                    provenance._stable_fields(os.fstat(self.prefix.descriptor)), "compiler prefix descriptor changed")
            for tool in self._tools.values(): tool.validate_current()
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def run(self, argv, *, maximum_bytes=1 << 20, timeout=600):
        self.validate_current()
        record = None
        try:
            effective = driver_arguments(argv, self.logical_driver, self.prefix.descriptor)
            record = {"logical_argv": list(argv), "effective_argv": effective, "status": "running"}
            self._commands.append(record)
            output = provenance._run(effective, "v19 sealed GCC phase", maximum_bytes=maximum_bytes, timeout=timeout,
                inherited_descriptors=(self.prefix.descriptor, *(tool.executable_descriptor for tool in self._tools.values())),
                executable_descriptor=self._tools["driver"].executable_descriptor,
                environment_overrides=builder.ENVIRONMENT)
            self.validate_current()
            record.update(status="exit-zero", stdout_sha256=builder.hashlib.sha256(output).hexdigest())
            return output
        except BaseException as error:
            if record is not None:
                record.update(status="failed", failure=type(error).__name__ + ": " + str(error))
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-compiler-execution/v1", "language": self.language,
            "prefix": str(self.root), "mappings": self._mappings, "commands": self._commands,
            "tools": {role: {"path": self.pins[role]["path"], **tool.executable_record()}
                      for role, tool in self._tools.items()},
            "sealed_driver_and_retained_helper_prefix": True, "runtime_closure_verified": False,
            "compiler_subtool_execution_owned": False, "atomic_snapshot": False,
            "compiler_data_owned": False, "fresh_build_recipe_integrated": False,
            "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise host.PreflightError("failed compiler phase cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
