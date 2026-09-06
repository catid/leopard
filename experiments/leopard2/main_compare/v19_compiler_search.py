#!/usr/bin/python3
"""Retain declared GCC search inputs; leopard-79h.38.5.4.8.2.2.2.3.

This is a bounded search-input owner, not a sandbox or complete search trace.
It detects changes before accepting a child result; it cannot prevent a child
from observing a transient file. Callers must discard failed jobs and outputs.
"""
from collections import deque
from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import stat

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_linker_inputs.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("compiler search dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_search_link_inputs", dependency)
link_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(link_module)
runtime = link_module.header_module.runtime
builder, provenance, require = runtime.builder, runtime.provenance, runtime.require
GCC = "/usr/lib/gcc/x86_64-linux-gnu/13"
CROSS = GCC + "/../../../../x86_64-linux-gnu"
# Original syscall spellings matter: missing/x/../y is not equivalent to y.
GCC13_ABSENT_PATHS = tuple(sorted(
    [root + "/" + name for root in (GCC, str(Path(GCC).parent),
         "/usr/libexec/gcc/x86_64-linux-gnu/13", "/usr/libexec/gcc/x86_64-linux-gnu")
     for name in ("collect-ld", "real-ld", "gnm", "gstrip", "nm", "strip")]
    + [GCC + "/specs", str(Path(GCC).parent) + "/specs", "/usr/bin/gnm", "/usr/bin/gstrip",
       "/lib/x86_64-linux-gnu/13/.", "/usr/lib/x86_64-linux-gnu/13/.",
       GCC + "/../../../x86_64-linux-gnu/13/."]
    + [CROSS + "/" + name for name in ("bin/.", "bin/x86_64-linux-gnu/.",
       "bin/x86_64-linux-gnu/13/.", "lib/.", "lib/../lib/.", "lib/x86_64-linux-gnu/.",
       "lib/x86_64-linux-gnu/13/.", "lib/specs", "lib/x86_64-linux-gnu/13/specs")]))
MAX_PATHS, MAX_DIRECTORIES, MAX_LINKS, MAX_STEPS = 128, 256, 64, 128
PRESENT_PATHS = ("/usr/bin/nm", "/usr/bin/strip", "/usr/libexec/gcc/x86_64-linux-gnu/13/lto-wrapper")
MAX_PRESENT_FILES, MAX_PRESENT_FILE_BYTES, MAX_PRESENT_BYTES = 8, 2 << 20, 4 << 20


def components(value, *, absolute):
    require(type(value) is str and 0 < len(os.fsencode(value)) < 4096 and "\0" not in value,
            "compiler search path exceeds profile")
    require(not absolute or value.startswith("/"), "compiler search path is not absolute")
    parts = value.split("/")[1:] if value.startswith("/") else value.split("/")
    require(0 < len(parts) <= 64 and all(parts), "compiler search path components differ")
    return parts


class CompilerSearch:
    """Borrow the dispatcher's runtime; hold absent and present search inputs.

    Production uses the fixed observed GCC13 path set. The private paths seam
    permits small real-filesystem mutation tests, not arbitrary build recipes.
    """
    def __init__(self, inventory, present_pins, *, _paths=None, _file_factory=None):
        require(type(inventory) is runtime.RuntimeInventory or _paths is not None,
                "compiler search requires a live runtime inventory")
        self.inventory = inventory
        self.paths = GCC13_ABSENT_PATHS if _paths is None else _paths
        require(type(self.paths) is tuple and 0 < len(self.paths) <= MAX_PATHS and
                all(type(path) is str for path in self.paths) and len(set(self.paths)) == len(self.paths),
                "compiler search path inventory differs")
        for path in self.paths: components(path, absolute=True)
        require(type(present_pins) is list and len(present_pins) <= MAX_PRESENT_FILES,
                "present compiler search inventory exceeds bound")
        names, targets, total = set(), set(), 0
        for pin in present_pins:
            require(type(pin) is dict and set(pin) == {"path", "resolved_path", "size", "sha256"} and
                    type(pin["path"]) is str and type(pin["resolved_path"]) is str and
                    type(pin["size"]) is int and 0 < pin["size"] <= MAX_PRESENT_FILE_BYTES and
                    type(pin["sha256"]) is str and builder.re.fullmatch(r"[0-9a-f]{64}", pin["sha256"]),
                    "present compiler search pin differs")
            for key in ("path", "resolved_path"): builder.host.canonical_path(pin[key], nonroot=True)
            require(pin["path"] not in names and pin["resolved_path"] not in targets,
                    "present compiler searches duplicate a name or target")
            names.add(pin["path"]); targets.add(pin["resolved_path"])
            total += pin["size"]
        require(total <= MAX_PRESENT_BYTES and not (names | targets).intersection(self.paths),
                "present compiler search bytes or absence overlap differ")
        self._present_pins = copy.deepcopy(present_pins)
        self._factory = builder._StreamedTool if _file_factory is None else _file_factory
        if _paths is None:
            phase = inventory.phase
            require(phase.language in ("c", "c++") and
                    str(phase._tools["collect2"].path) == "/usr/libexec/gcc/x86_64-linux-gnu/13/collect2",
                    "compiler search requires the qualified GCC13 profile")
            require(names == set(PRESENT_PATHS), "qualified present compiler search coverage differs")
        self._stack, self._state = ExitStack(), "new"
        self._directories, self._links, self._missing = {}, {}, []
        self._files, self._file_aliases = {}, {}

    def _directory(self, path):
        if path not in self._directories:
            require(len(self._directories) < MAX_DIRECTORIES, "compiler search directory bound exceeded")
            fd = builder.host.LinuxReader.open_directory(str(path))
            self._stack.callback(os.close, fd)
            self._guard._add_watch(path, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, {b""})
            self._directories[path] = fd, builder.streamed._directory_identity(os.fstat(fd))
        return self._directories[path][0]

    def _capture(self, request):
        pending, parent = deque(components(request, absolute=True)), Path("/")
        for _ in range(MAX_STEPS):
            fd = self._directory(parent)
            if not pending: raise builder.host.PreflightError("compiler search input already exists: " + request)
            name = pending.popleft()
            if name == ".": continue
            if name == "..":
                parent = parent.parent
                continue
            self._guard._add_watch(parent, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR |
                                   provenance.IN_ATTRIB, {b"", os.fsencode(name)})
            try:
                value = os.stat(name, dir_fd=fd, follow_symlinks=False)
            except FileNotFoundError:
                self._missing.append({"path": request, "parent": str(parent), "name": name})
                return
            path = parent / name
            if stat.S_ISLNK(value.st_mode):
                require(len(self._links) < MAX_LINKS or path in self._links, "compiler search alias bound exceeded")
                self._guard._add_watch(path, self._guard._FILE_MASK, None)
                target = os.readlink(name, dir_fd=fd)
                require(provenance._stable_fields(os.stat(name, dir_fd=fd, follow_symlinks=False)) ==
                        provenance._stable_fields(value), "compiler search alias changed while captured")
                row = provenance._stable_fields(value), target
                require(path not in self._links or self._links[path] == row, "compiler search alias changed")
                self._links[path] = row
                pending.extendleft(reversed(components(target, absolute=False)))
                require(len(pending) <= MAX_STEPS, "compiler search alias expansion exceeds bound")
                if target.startswith("/"): parent = Path("/")
            else:
                require(stat.S_ISDIR(value.st_mode), "compiler search input exists or traversal is not a directory")
                parent = path
        raise builder.host.PreflightError("compiler search traversal exceeds bound")

    def __enter__(self):
        require(self._state == "new", "compiler search cannot be reused")
        self._state = "entering"
        try:
            self.inventory.validate_current()
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 compiler search inputs"))
            for path in self.paths: self._capture(path)
            inodes = set()
            for pin in self._present_pins:
                tool = self._factory(Path(pin["resolved_path"]), maximum_bytes=MAX_PRESENT_FILE_BYTES, _guard=self._guard)
                self._stack.callback(tool.close)
                value = os.fstat(tool.fd)
                require(tool.sha256 == pin["sha256"] and value.st_size == pin["size"],
                        "present compiler search bytes differ from pin")
                require((value.st_dev, value.st_ino) not in inodes, "present compiler search files alias one inode")
                inodes.add((value.st_dev, value.st_ino))
                aliases = self._stack.enter_context(runtime.compiler._retain_driver_aliases(pin["path"], tool))
                self._files[pin["path"]] = tool
                self._file_aliases[pin["path"]] = aliases
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    @staticmethod
    def _require_absent(path, **kwargs):
        try: os.stat(path, **kwargs)
        except FileNotFoundError: return
        raise builder.host.PreflightError("compiler search input appeared: " + str(path))

    def validate_current(self):
        require(self._state == "held", "compiler searches are not held")
        try:
            self.inventory.validate_current()
            self._guard.verify()
            require(not os.get_inheritable(self._guard.descriptor), "compiler search guard inheritance changed")
            for path, (fd, identity) in self._directories.items():
                require(not os.get_inheritable(fd) and path.resolve(strict=True) == path and
                        builder.streamed._directory_identity(os.fstat(fd)) == identity ==
                        builder.streamed._directory_identity(path.lstat()), "compiler search parent changed")
            for path, (fields, target) in self._links.items():
                require(provenance._stable_fields(path.lstat()) == fields and os.readlink(path) == target,
                        "compiler search alias changed")
            for row in self._missing:
                self._require_absent(row["name"], dir_fd=self._directories[Path(row["parent"])][0], follow_symlinks=False)
                self._require_absent(row["path"])
            for aliases in self._file_aliases.values(): aliases.validate_current()
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-compiler-search/v2", "missing": self._missing,
            "directories": {str(path): list(identity) for path, (_, identity) in self._directories.items()},
            "aliases": {str(path): {"stable_fields": list(fields), "target": target}
                        for path, (fields, target) in self._links.items()},
            "present_files": {path: {"resolved_path": str(tool.path), "sha256": tool.sha256,
                "size": os.fstat(tool.fd).st_size, "stable_fields": list(tool.fields),
                "aliases": self._file_aliases[path].record()} for path, tool in self._files.items()},
            "present_file_bytes": sum(pin["size"] for pin in self._present_pins),
            "declared_present_searches_retained": bool(self._present_pins),
            "present_files_pinned_by_original_preflight": False, "queried_file_execution_owned": False,
            "declared_absence_history_owned": True, "negative_search_closure_owned": False,
            "compiler_data_owned": False, "full_runtime_execution_owned": False,
            "fresh_build_recipe_integrated": False, "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise builder.host.PreflightError("failed compiler searches cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
