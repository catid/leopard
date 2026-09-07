#!/usr/bin/python3
"""Retained GCC/build-tool startup dependencies; leopard-79h.38.5.4.8.2.2.2.3.

Linux ELF64/x86-64 system-library profile only. This API inventories and seals
ELF startup inputs, lists resolutions and runs explicit build-root jobs. It does
can borrow an explicit CMake priority-policy owner. Other configuration, dlopen
plugins, compiler dispatch and full build integration remain separate obligations.
"""
from __future__ import annotations

from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import re
import secrets
import struct

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_build_tool_execution.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("runtime inventory dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_runtime_build_tools", dependency)
build_tools = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_tools)
compiler = build_tools.compiler
builder, provenance, require = compiler.builder, compiler.provenance, compiler.require
LIBRARY_ROOT = Path("/usr/lib/x86_64-linux-gnu")
LOADER = LIBRARY_ROOT / "ld-linux-x86-64.so.2"
INTERPRETER = "/lib64/ld-linux-x86-64.so.2"
MAX_FILES, MAX_BYTES = 64, 64 << 20
NAME = re.compile(r"[A-Za-z0-9_+.-]{1,255}")


def elf_startup(read, size):
    """Read bounded program/dynamic metadata, never an entire executable body.

    Constants/layouts follow the installed elf.h ABI. Unsupported search-path,
    filter, configuration and audit tags fail closed instead of being ignored.
    No section table, external inspector or executable invocation is required.
    """
    require(type(size) is int and 64 <= size <= 64 << 20, "ELF size is outside the profile")
    def get(offset, count):
        require(0 <= offset <= size and 0 <= count <= 65536 and count <= size - offset, "ELF read outside bounds")
        data = read(offset, count)
        require(type(data) is bytes and len(data) == count, "ELF read truncated")
        return data
    h = struct.unpack("<16sHHIQQQIHHHHHH", get(0, 64))
    require(h[0][:7] == b"\x7fELF\x02\x01\x01" and h[1] in (2, 3) and h[2] == 62 and h[3] == 1 and
            h[8] == 64 and h[9] == 56 and 0 < h[10] <= 128 and h[5] >= 64, "unsupported ELF header")
    table = get(h[5], h[9] * h[10])
    segments = [struct.unpack_from("<IIQQQQQQ", table, offset) for offset in range(0, len(table), 56)]
    loads, dynamic, interpreter = [], [], []
    for kind, flags, offset, address, physical, filesz, memsz, alignment in segments:
        if kind not in (1, 2, 3): continue
        require(filesz <= memsz and offset <= size and filesz <= size - offset, "ELF segment exceeds file")
        if kind == 1: loads.append((offset, address, filesz))
        elif kind == 2: dynamic.append((offset, filesz))
        else: interpreter.append((offset, filesz))
    require(loads and len(dynamic) == 1 and len(interpreter) <= 1, "ELF runtime segment coverage differs")
    interp = None
    if interpreter:
        offset, count = interpreter[0]
        require(1 < count <= 4096, "ELF interpreter exceeds bound")
        data = get(offset, count)
        require(data.endswith(b"\0") and b"\0" not in data[:-1], "ELF interpreter is not framed")
        interp = data[:-1].decode("ascii", errors="strict")
        require(interp == INTERPRETER, "ELF interpreter differs from supported profile")
    offset, count = dynamic[0]
    require(0 < count <= 65536 and count % 16 == 0, "ELF dynamic table exceeds bound")
    data = get(offset, count)
    needed, strings, ended = [], {}, False
    for offset in range(0, count, 16):
        tag, value = struct.unpack_from("<qQ", data, offset)
        if ended:
            require(tag == value == 0, "ELF dynamic data follows terminator")
            continue
        if tag == 0:
            require(value == 0, "ELF dynamic terminator value differs")
            ended = True
        elif tag == 1:
            needed.append(value)
            require(len(needed) <= 64, "ELF dependency count exceeds bound")
        elif tag in (5, 10, 14):
            require(tag not in strings, "duplicate ELF string metadata")
            strings[tag] = value
        else:
            require(tag not in (15, 29, 0x6ffffefa, 0x6ffffefb, 0x6ffffefc, 0x7ffffffd, 0x7fffffff),
                    "unsupported ELF search/configuration/audit/filter tag")
    # GCC frontends have >1 MiB dynamic string tables. Only the at-most-64
    # dependency names are read (256 bytes each), not that entire table.
    require(ended and 5 in strings and 10 in strings and 0 < strings[10] <= size,
            "ELF string table is absent or exceeds bound")
    ranges = [(offset + strings[5] - address, strings[10]) for offset, address, count in loads
              if address <= strings[5] and strings[5] - address <= count and strings[10] <= count - (strings[5] - address)]
    require(len(ranges) == 1, "ELF string table does not map uniquely to file bytes")
    start, length = ranges[0]
    def string(offset):
        require(offset < length, "ELF dependency string offset is outside table")
        data = get(start + offset, min(256, length - offset))
        require(b"\0" in data, "ELF dependency string is unterminated or oversized")
        value = data.split(b"\0", 1)[0].decode("ascii", errors="strict")
        require(NAME.fullmatch(value) and value not in (".", ".."), "unsafe ELF dependency name")
        return value
    names = [string(offset) for offset in needed]
    require(len(names) == len(set(names)), "duplicate ELF dependency name")
    return {"interpreter": interp, "needed": names, "soname": string(strings[14]) if 14 in strings else None}


def validate_listing(output, root, metadata, sonames, prefix, loader):
    """Require the complete transitive startup set and exact sealed-prefix paths.

    glibc displays the explicitly invoked loader's argv[0], not its memfd path.
    That one display label is bound separately to the sealed executable fd.
    """
    require(type(output) is bytes and 0 < len(output) <= 65536, "loader listing exceeds bound")
    expected, visited, queue = {Path(loader).name}, set(), [root]
    for path in queue:
        if path in visited: continue
        require(path in metadata and len(visited) <= MAX_FILES + 6, "loader closure is incomplete or oversized")
        visited.add(path)
        for name in metadata[path]["needed"]:
            require(name in sonames, "loader dependency lacks a sealed resolution")
            expected.add(name)
            queue.append(sonames[name])
    observed, virtual = set(), 0
    for line in output.decode("utf-8", errors="strict").splitlines():
        line = line.strip()
        if re.fullmatch(r"linux-vdso\.so\.1 \(0x[0-9a-f]+\)", line):
            virtual += 1
            continue
        if re.fullmatch(re.escape(loader) + r" \(0x[0-9a-f]+\)", line):
            require(metadata[root]["interpreter"] is None, "unexpected direct loader display")
            name = Path(loader).name
        else:
            match = re.fullmatch(r"([^\s]+) => (/[^\s]+) \(0x[0-9a-f]+\)", line)
            require(match is not None, "unrecognized or missing loader dependency")
            name, path = match.groups()
            if name == INTERPRETER:
                name = Path(loader).name
                require(metadata[root]["interpreter"] == INTERPRETER and path == loader, "explicit loader display differs")
            else:
                require(NAME.fullmatch(name) and path == prefix + "/" + name, "loader resolved outside sealed library prefix")
        require(name not in observed, "duplicate loader dependency")
        observed.add(name)
    require(virtual == 1 and observed == expected, "loader resolution does not equal transitive ELF dependencies")
    return sorted(observed)


class RuntimeInventory:
    """Borrow a live compiler phase; retain newly observed loader/library bytes.

    Newly observed libraries are NOT retroactively pinned by the original
    preflight. The private injections are for small real-file test fixtures.
    Enter and exit inside the compiler phase, under its caller's host/lock lease.
    """
    def __init__(self, phase, *, _tool_factory=None, _library_root=None, _loader=None, _include_plugin=True,
                 _preload=Path("/etc/ld.so.preload")):
        require(type(phase) in (compiler.CompilerExecution, build_tools.BuildToolExecution) or _tool_factory is not None,
                "runtime requires a live compiler or build-tool phase")
        self.phase = phase
        self._build_tools = type(phase) is build_tools.BuildToolExecution
        self.library_root = LIBRARY_ROOT if _library_root is None else Path(_library_root)
        self.loader_path = LOADER if _loader is None else Path(_loader)
        self.preload = Path(_preload)
        self._factory = builder._StreamedTool if _tool_factory is None else _tool_factory
        self._plugin = _include_plugin and not self._build_tools
        self._stack, self._state = ExitStack(), "new"
        self._files, self._metadata, self._links, self._roots, self._lists = {}, {}, {}, {}, {}
        self._total = 0
        self._commands = []

    def _add(self, path):
        path = Path(path)
        self._guard.add_file_path(path)
        target = path.resolve(strict=True)
        self._links[str(path)] = str(target)
        if str(target) in self._files: return str(target)
        require(len(self._files) < MAX_FILES, "runtime file count exceeds bound")
        tool = self._factory(target, maximum_bytes=64 << 20, permitted_modes=(0o644, 0o755))
        self._stack.callback(tool.close)
        value = os.fstat(tool.fd)
        self._total += value.st_size
        require(self._total <= MAX_BYTES, "runtime file bytes exceed bound")
        require(all((value.st_dev, value.st_ino) != (os.fstat(other.fd).st_dev, os.fstat(other.fd).st_ino)
                    for other in (*self._files.values(), *self.phase._tools.values())), "runtime files alias other owned inputs")
        self._files[str(target)] = tool
        return str(target)

    def __enter__(self):
        require(self._state == "new", "runtime inventory cannot be reused")
        self._state = "entering"
        try:
            self.phase.validate_current()
            require(self.library_root.resolve(strict=True) == self.library_root and self.loader_path.resolve(strict=True) == self.loader_path,
                    "runtime roots are not canonical")
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 runtime inputs"))
            self._guard.add_exact_directory_entries(self.preload.parent, (self.preload.name,))
            require(not os.path.lexists(self.preload), "global loader preload configuration exists")
            self._roots = {role: str(tool.path) for role, tool in self.phase._tools.items()}
            self.loader_key = self._add(self.loader_path)
            if self._plugin:
                path = Path(self.phase.pins["collect2"]["path"]).parent / "liblto_plugin.so"
                self._roots["gcc-link-plugin"] = self._add(path)
            queue = list(self.phase._tools.values()) + list(self._files.values())
            mappings = {}
            for tool in queue:
                key = str(tool.path)
                if key in self._metadata: continue
                fd = tool.executable_descriptor
                try:
                    record = elf_startup(lambda offset, count: os.pread(fd, count, offset), os.fstat(fd).st_size)
                except (builder.host.PreflightError, ValueError) as error:
                    raise builder.host.PreflightError(f"runtime ELF {tool.path}: {error}") from error
                self._metadata[key] = record
                for name in record["needed"]:
                    path = self.library_root / name
                    target = path.resolve(strict=True)
                    require(target.is_relative_to(self.library_root), "runtime dependency escapes system library root")
                    key = self._add(path)
                    require(name not in mappings or mappings[name] == key, "conflicting SONAME resolution")
                    mappings[name] = key
                    queue.append(self._files[key])
            self._sonames = mappings
            require(all(self._metadata[key]["soname"] == name for name, key in mappings.items()),
                    "runtime SONAME differs from dependency name")
            name = "v19-libs-" + secrets.token_hex(16)
            os.mkdir(name, mode=0o700, dir_fd=self.phase._parent_fd)
            self.root = self.phase.parent / name
            fd = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=self.phase._parent_fd)
            self._stack.callback(os.close, fd)
            self.phase.validate_current()
            links = {name: f"/proc/self/fd/{self._files[key].executable_descriptor}" for name, key in mappings.items()}
            for name, target in links.items(): os.symlink(target, name, dir_fd=fd)
            os.fchmod(fd, 0o500)
            self._guard.verify()
            self._guard.add_directory_path(self.root)
            self._guard._add_watch(self.root, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR | provenance.IN_ATTRIB, None)
            self.prefix = self._stack.enter_context(provenance._retain_exact_symlink_directory(self.root, links, "v19 sealed library prefix"))
            require(provenance._stable_fields(os.fstat(fd)) == provenance._stable_fields(os.fstat(self.prefix.descriptor)),
                    "runtime prefix was replaced")
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def validate_current(self):
        require(self._state == "held", "runtime inventory is not held")
        try:
            self.phase.validate_current()
            self._guard.verify()
            require(not os.path.lexists(self.preload), "global loader preload configuration appeared")
            for path, target in self._links.items(): require(str(Path(path).resolve(strict=True)) == target, "runtime resolution changed")
            for tool in self._files.values(): tool.validate_current()
            self.prefix.verify()
            require(not os.get_inheritable(self.prefix.descriptor), "runtime prefix descriptor inheritance changed")
            self._guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def list_dependencies(self, role):
        """Inspect one sealed root with an explicit sealed loader; no tool job."""
        self.validate_current()
        try:
            require(type(role) is str and role in self._roots, "unknown runtime root")
            require(role not in self._lists, "runtime root was already inspected")
            target = self.phase._tools[role] if role in self.phase._tools else self._files[self._roots[role]]
            loader = self._files[self.loader_key]
            argv = [str(self.loader_path), "--list", "--inhibit-cache", "--glibc-hwcaps-mask", "",
                    "--library-path", f"/proc/self/fd/{self.prefix.descriptor}", f"/proc/self/fd/{target.executable_descriptor}"]
            output = provenance._run(argv, "v19 sealed loader dependency listing", maximum_bytes=65536, timeout=30,
                executable_descriptor=loader.executable_descriptor,
                inherited_descriptors=(self.prefix.descriptor, target.executable_descriptor,
                                       *(tool.executable_descriptor for tool in self._files.values())),
                environment_overrides=builder.ENVIRONMENT)
            self.validate_current()
            resolved = validate_listing(output, self._roots[role], self._metadata, self._sonames,
                                        f"/proc/self/fd/{self.prefix.descriptor}", str(self.loader_path))
            self._lists[role] = {"argv": argv, "stdout": output.decode("utf-8", errors="strict"),
                                 "stdout_sha256": builder.hashlib.sha256(output).hexdigest(),
                                 "resolved_dependencies": resolved, "sealed_prefix_resolution_verified": True}
            return output
        except BaseException:
            self._state = "failed"
            raise

    def run_tool(self, role, argv, *, input_descriptors=(), configuration=None, _runner=None):
        """Run one declared build root through the sealed loader/startup runtime.

        The caller still owns the exact recipe and inputs. Extra input fds must
        already be fully sealed; this is not dlopen, nested-tool or complete data routing.
        A same-phase configuration owner redirects CMake's GnuTLS policy only.
        Compiler jobs must use RuntimeDispatch instead.
        """
        self.validate_current()
        record = None
        try:
            require(self._build_tools, "direct runtime jobs require the build-tool profile")
            logical = self.phase.arguments(role, argv)
            require(configuration is None or (type(configuration) is build_tools.CMakePriorityConfiguration and
                    configuration.phase is self.phase and role == "cmake"), "build job configuration owner differs")
            environment = dict(builder.ENVIRONMENT)
            configuration_fds = ()
            if configuration is not None:
                configuration_before = configuration.record()
                environment.update(configuration_before["environment"])
                configuration_fds = (configuration.descriptor(),)
            require(type(input_descriptors) is tuple and all(type(fd) is int and 3 <= fd < 65536 for fd in input_descriptors) and
                    len(input_descriptors) <= 64 and len(set(input_descriptors)) == len(input_descriptors),
                    "build input descriptors differ")
            def input_fields():
                fields = []
                for fd in input_descriptors:
                    require(not os.get_inheritable(fd) and
                            builder.fcntl.fcntl(fd, getattr(builder.fcntl, "F_GET_SEALS", 1034)) == 15,
                            "build input descriptor is not privately sealed")
                    fields.append(provenance._stable_fields(os.fstat(fd)))
                return fields
            inputs = input_fields()
            target, loader = self.phase._tools[role], self._files[self.loader_key]
            loader_argv = [str(self.loader_path), "--inhibit-cache", "--glibc-hwcaps-mask", "", "--library-path",
                f"/proc/self/fd/{self.prefix.descriptor}", "--argv0", logical[0],
                f"/proc/self/fd/{target.executable_descriptor}", *logical[1:]]
            require(len(loader_argv) <= 512, "effective build tool arguments exceed bound")
            record = {"role": role, "logical_argv": logical, "loader_argv": loader_argv,
                      "input_descriptors": list(input_descriptors), "status": "running"}
            if configuration is not None: record["configuration"] = configuration_before
            self._commands.append(record)
            inherited = (self.prefix.descriptor, *(tool.executable_descriptor for tool in self.phase._tools.values()),
                         *(tool.executable_descriptor for tool in self._files.values()), *input_descriptors, *configuration_fds)
            output = (provenance._run if _runner is None else _runner)(loader_argv, "v19 sealed build tool " + role,
                maximum_bytes=1 << 20, timeout=600, executable_descriptor=loader.executable_descriptor,
                inherited_descriptors=tuple(dict.fromkeys(inherited)), environment_overrides=environment)
            require(type(output) is bytes and len(output) <= 1 << 20, "build tool output exceeds bound")
            self.validate_current()
            require(input_fields() == inputs, "build input descriptor identity changed")
            if configuration is not None:
                require(configuration.record() == configuration_before, "CMake policy ownership changed during job")
            record.update(status="exit-zero", stdout_sha256=builder.hashlib.sha256(output).hexdigest())
            return output
        except BaseException as error:
            if record is not None: record.update(status="failed", failure=type(error).__name__ + ": " + str(error))
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-runtime-inventory/v2" if self._build_tools else "leopard2-v19-runtime-inventory/v1",
            **({"root_profile": "build-tools", "commands": self._commands} if self._build_tools else {}), "roots": self._roots,
            "metadata": self._metadata, "sonames": self._sonames, "loader": self.loader_key,
            "files": {path: {**tool.executable_record(), "source_mode": os.fstat(tool.fd).st_mode}
                      for path, tool in self._files.items()}, "bytes": self._total, "loader_listings": self._lists,
            "library_descriptor_mappings": {name: f"/proc/self/fd/{self._files[key].executable_descriptor}"
                                            for name, key in self._sonames.items()},
            "sealed_library_prefix": str(self.root), "libraries_pinned_by_original_preflight": False,
            "compiler_data_owned": False, "full_runtime_execution_owned": False, "fresh_build_recipe_integrated": False,
            "atomic_snapshot": False, "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise builder.host.PreflightError("failed runtime inventory cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
