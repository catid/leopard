#!/usr/bin/python3
"""Bounded sealed GCC header view; leopard-79h.38.5.4.8.2.2.2.3.

Owns a caller-pinned, observed header set, not arbitrary compiler data or every
possible include. Callers still retain source/recipe ownership and must prove
actual header reads and complete object equality for each qualified recipe.
"""
from __future__ import annotations
from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import re
import secrets
import stat

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_runtime_inventory.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("compiler header dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_headers_runtime", dependency)
runtime = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtime)
builder, provenance, require = runtime.builder, runtime.provenance, runtime.require

CPP_INCLUDE_ROOTS = ("/usr/include/c++/13", "/usr/include/x86_64-linux-gnu/c++/13",
                     "/usr/include/c++/13/backward", "/usr/lib/gcc/x86_64-linux-gnu/13/include",
                     "/usr/local/include", "/usr/include/x86_64-linux-gnu", "/usr/include")
C_INCLUDE_ROOTS = CPP_INCLUDE_ROOTS[3:]
C_PREDEFINITION_HEADER = "/usr/include/stdc-predef.h"
MAX_FILES, MAX_DIRECTORIES, MAX_FILE_BYTES, MAX_TOTAL_BYTES = 512, 128, 2 << 20, 16 << 20


def checked_path(value):
    require(type(value) is str and len(value) < 4096 and len(Path(value).parts) <= 32 and
            re.fullmatch(r"/[A-Za-z0-9_./+-]+", value), "header path exceeds qualified profile")
    builder.host.canonical_path(value, nonroot=True)
    return Path(value)


def header_pins(rows, include_roots, language="c++"):
    require(type(language) is str and language in ("c", "c++"), "unsupported header language")
    expected = C_INCLUDE_ROOTS if language == "c" else CPP_INCLUDE_ROOTS
    require(type(include_roots) is tuple and include_roots == expected,
            "header include order is not the qualified language profile")
    roots = tuple(checked_path(value) for value in include_roots)
    require(type(rows) is list and 0 < len(rows) <= MAX_FILES, "header count exceeds bound")
    pins = {}
    for row in rows:
        require(type(row) is dict and set(row) == {"path", "sha256", "size"} and
                type(row["size"]) is int and 0 < row["size"] <= MAX_FILE_BYTES and
                type(row["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]),
                "header pin differs")
        path = checked_path(row["path"])
        require(path not in pins and any(path != root and path.is_relative_to(root) for root in roots),
                "duplicate or out-of-profile header")
        pins[path] = copy.deepcopy(row)
    require(sum(row["size"] for row in pins.values()) <= MAX_TOTAL_BYTES, "header set exceeds byte bound")
    return pins


class _PinnedInputView:
    """Shared descriptor/namespace ownership, with caller-qualified mappings.

    Header and linker policies supply their own pins, limits and arguments;
    this private primitive does not confer recipe or search-closure authority.
    """
    def __init__(self, inventory, pins, source_roots, aliases, *, limits, root_prefix, _file_factory=None):
        require(type(inventory) is runtime.RuntimeInventory or _file_factory is not None,
                "input view requires a live runtime inventory")
        self.inventory, self.phase = inventory, inventory.phase
        self._pins, self._source_roots = copy.deepcopy(pins), tuple(source_roots)
        self._max_files, self._max_directories, self._max_file_bytes = limits
        self._root_prefix = root_prefix
        require(type(aliases) is dict and 0 < len(aliases) <= self._max_files * 4 and
                all(isinstance(path, Path) and not path.is_absolute() and str(path) != "." and
                    checked_path("/" + str(path)).relative_to("/") == path and target in self._pins
                    for path, target in aliases.items()) and set(aliases.values()) == set(self._pins),
                "input view alias inventory differs")
        self._aliases = copy.deepcopy(aliases)
        self._factory = builder._StreamedTool if _file_factory is None else _file_factory
        self._stack, self._state = ExitStack(), "new"
        self._files, self._source_dirs, self._view_dirs, self._mappings = {}, {}, {}, {}

    def __enter__(self):
        require(self._state == "new", "header owner cannot be reused")
        self._state = "entering"
        try:
            self.inventory.validate_current()
            self._source_guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 header inputs"))
            source_directories = {parent for path in self._pins for parent in path.parents}
            source_directories.update(Path(value) for value in self._source_roots)
            source_directories.update(parent for value in self._source_roots for parent in Path(value).parents)
            require(len(source_directories) <= self._max_directories, "input source directory count exceeds bound")
            for path in sorted(source_directories, key=str):
                require(path.resolve(strict=True) == path, "header directory alias is outside profile")
                if path != Path("/"): self._source_guard.add_directory_path(path)
                self._source_guard._add_watch(path, self._source_guard._DIRECTORY_MASK |
                    provenance.IN_ONLYDIR | provenance.IN_ATTRIB, {b""})
                fd = builder.host.LinuxReader.open_directory(str(path))
                self._stack.callback(os.close, fd)
                self._source_dirs[path] = fd, builder.streamed._directory_identity(os.fstat(fd))
            identities = set()
            for path, pin in self._pins.items():
                file = self._factory(path, maximum_bytes=self._max_file_bytes, permitted_modes=(0o644, 0o755),
                                     _guard=self._source_guard)
                self._stack.callback(file.close)
                self._files[path] = file
                value = os.fstat(file.fd)
                require((file.sha256, value.st_size) == (pin["sha256"], pin["size"]), "header differs from pin")
                inode = value.st_dev, value.st_ino
                require(inode not in identities, "header pins alias one inode")
                identities.add(inode)
                require(3 <= file.executable_descriptor < 65536, "header descriptor exceeds profile")
            name = self._root_prefix + secrets.token_hex(16)
            os.mkdir(name, mode=0o700, dir_fd=self.phase._parent_fd)
            self.root = self.phase.parent / name
            fd = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=self.phase._parent_fd)
            self._stack.callback(os.close, fd)
            self._view_dirs[Path(".")] = fd, None
            self.phase.validate_current()
            require(self.root.resolve(strict=True) == self.root and
                    builder.streamed._directory_identity(os.fstat(fd)) ==
                    builder.streamed._directory_identity(self.root.lstat()), "header root was redirected")
            directories = {parent for relative in self._aliases for parent in relative.parents if parent != Path(".")}
            # Include-root parents are needed even when a root has no headers.
            directories.update(Path(value).relative_to("/") for value in self._source_roots)
            directories.update(parent.relative_to("/") for value in self._source_roots
                               for parent in Path(value).parents if parent != Path("/"))
            require(len(directories) + 1 <= self._max_directories, "input view directory count exceeds bound")
            for relative in sorted(directories, key=lambda value: (len(value.parts), str(value))):
                parent_fd = self._view_dirs[relative.parent][0]
                os.mkdir(relative.name, mode=0o700, dir_fd=parent_fd)
                fd = os.open(relative.name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=parent_fd)
                self._stack.callback(os.close, fd)
                self._view_dirs[relative] = fd, None
            for relative, path in self._aliases.items():
                file = self._files[path]
                target = f"/proc/self/fd/{file.executable_descriptor}"
                os.symlink(target, relative.name, dir_fd=self._view_dirs[relative.parent][0])
                self._mappings[relative] = target
            for relative, (fd, _) in self._view_dirs.items():
                os.fchmod(fd, 0o500)
                self._view_dirs[relative] = fd, provenance._stable_fields(os.fstat(fd))
            self._view_guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 header view"))
            self._view_guard.add_directory_path(self.root)
            for relative in self._view_dirs:
                self._view_guard._add_watch(self.root / relative, self._view_guard._DIRECTORY_MASK |
                    provenance.IN_ONLYDIR | provenance.IN_ATTRIB, None)
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def validate_current(self):
        require(self._state == "held", "header owner is not held")
        try:
            self.inventory.validate_current()
            self._source_guard.verify()
            self._view_guard.verify()
            for path, (fd, identity) in self._source_dirs.items():
                require(not os.get_inheritable(fd) and path.resolve(strict=True) == path and
                        builder.streamed._directory_identity(os.fstat(fd)) == identity ==
                        builder.streamed._directory_identity(path.lstat()), "header source directory changed")
            for file in self._files.values(): file.validate_current()
            require(self.root.resolve(strict=True) == self.root, "header view pathname changed")
            for relative, (fd, fields) in self._view_dirs.items():
                require(not os.get_inheritable(fd) and provenance._stable_fields(os.fstat(fd)) == fields ==
                        provenance._stable_fields((self.root / relative).lstat()), "header view directory changed")
                expected_dirs = {str(path.name) for path in self._view_dirs if path != Path(".") and path.parent == relative}
                expected_links = {path.name: target for path, target in self._mappings.items() if path.parent == relative}
                with os.scandir(fd) as entries:
                    observed = {}
                    for entry in entries:
                        require(len(observed) < len(self._aliases) + self._max_directories, "input view entry count exceeds bound")
                        observed[entry.name] = (os.readlink(entry.name, dir_fd=fd) if entry.is_symlink() else
                                               "directory" if entry.is_dir(follow_symlinks=False) else "unexpected")
                require(observed == {**dict.fromkeys(expected_dirs, "directory"), **expected_links}, "header view inventory changed")
            self._source_guard.verify()
            self._view_guard.verify()
        except BaseException:
            self._state = "failed"
            raise

    def descriptors(self):
        self.validate_current()
        return tuple(file.executable_descriptor for file in self._files.values())

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise builder.host.PreflightError("failed input view cannot complete")
        finally:
            self._state = "closed"
            self._stack.close()


class CompilerSourceInputs(_PinnedInputView):
    """Seal one pinned translation unit and its declared quoted includes.

    Quoted includes must be under the source's directory. Rewriting that one
    source argument lets its relative includes resolve inside the same view.
    Absolute includes and complete source/read closure are not claimed.
    """
    def __init__(self, inventory, source, quoted_headers, *, origin="system", _file_factory=None):
        require(type(inventory) is runtime.RuntimeInventory or _file_factory is not None,
                "source inputs require a live runtime inventory")
        require(type(origin) is str and origin in ("system", "generated"), "source input origin differs")
        require(inventory.phase.language in ("c", "c++"), "source input language differs")
        require(type(quoted_headers) is list and len(quoted_headers) < 128, "quoted header count exceeds bound")
        pins = {}
        for row in [source, *quoted_headers]:
            require(type(row) is dict and set(row) == {"path", "sha256", "size"} and
                    type(row["size"]) is int and 0 < row["size"] <= MAX_FILE_BYTES and
                    type(row["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]), "source input pin differs")
            path = checked_path(row["path"])
            require(path not in pins, "duplicate source input")
            pins[path] = copy.deepcopy(row)
        self.source = checked_path(source["path"])
        suffixes = (".c",) if inventory.phase.language == "c" else (".cpp", ".cxx", ".cc", ".C")
        require(self.source.suffix in suffixes and all(path.is_relative_to(self.source.parent) for path in pins),
                "source suffix or quoted header location differs")
        require(sum(row["size"] for row in pins.values()) <= MAX_TOTAL_BYTES, "source input bytes exceed bound")
        self.origin = origin
        trusted_owner = (0, 0) if origin == "system" else (os.getuid(), os.getgid())
        def factory(path, **kwargs):
            return builder._StreamedTool(path, _trusted_owner=trusted_owner, **kwargs)
        super().__init__(inventory, pins, (str(self.source.parent),),
            {path.relative_to("/"): path for path in pins},
            limits=(128, MAX_DIRECTORIES, MAX_FILE_BYTES), root_prefix="v19-source-inputs-",
            _file_factory=factory if _file_factory is None else _file_factory)

    def arguments(self, argv):
        self.validate_current()
        try:
            require(type(argv) is list and 0 < len(argv) <= 480 and argv[0] == self.phase.logical_driver and
                    all(type(value) is str and "\0" not in value for value in argv), "source arguments differ")
            count = argv.count(str(self.source))
            require(count <= 1, "source argument is duplicated")
            source_suffixes = (".c", ".cpp", ".cxx", ".cc", ".C")
            if not count:
                require("-c" not in argv and not any(value.endswith(source_suffixes) for value in argv[1:]),
                        "compile request does not use the pinned source")
                return list(argv)  # A separate object-only link still borrows this owner.
            require(argv.count("-c") <= 1 and argv.count("-o") == 1 and "-S" not in argv and "-E" not in argv and
                    not any(value != str(self.source) and value.endswith(source_suffixes) for value in argv[1:]),
                    "source profile requires one explicit translation unit")
            selected = str(self.root / str(self.source).lstrip("/"))
            return [argv[0], f"-ffile-prefix-map={self.root}=",
                    *[selected if value == str(self.source) else value for value in argv[1:]]]
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-compiler-source-inputs/v1", "root": str(self.root),
            "source": str(self.source), "origin": self.origin, "language": self.phase.language,
            "files": [dict(pin, descriptor=self._files[path].executable_descriptor,
                seals=self._files[path].executable_record()["seals"]) for path, pin in self._pins.items()],
            "input_bytes": sum(row["size"] for row in self._pins.values()),
            "declared_source_inputs_sealed": True, "full_source_read_closure_owned": False,
            "source_identity_owned": False, "fresh_build_recipe_integrated": False, "benchmark_executed": False})


class CompilerHeaders(_PinnedInputView):
    """Qualified C/C++ include policies over a bounded, sealed input view.

    The declared set is sealed; absolute includes, optional specs and complete
    positive/negative search equivalence remain separate obligations.
    """
    def __init__(self, inventory, pins, include_roots=None, *, _file_factory=None):
        require(inventory.phase.language in ("c", "c++"), "header profile requires C or C++")
        self.include_roots = ((C_INCLUDE_ROOTS if inventory.phase.language == "c" else CPP_INCLUDE_ROOTS)
                              if include_roots is None else include_roots)
        selected = header_pins(pins, self.include_roots, inventory.phase.language)
        require(Path(C_PREDEFINITION_HEADER) in selected,
                "header profile requires the implicit predefinition header")
        super().__init__(inventory, selected, self.include_roots, {path.relative_to("/"): path for path in selected},
            limits=(MAX_FILES, MAX_DIRECTORIES, MAX_FILE_BYTES), root_prefix="v19-headers-", _file_factory=_file_factory)

    def arguments(self, argv, *, compile_only=True):
        self.validate_current()
        try:
            require(type(compile_only) is bool, "header compilation mode is not boolean")
            require(type(argv) is list and len(argv) <= 480 and argv and argv[0] == self.phase.logical_driver and
                    all(type(value) is str and "\0" not in value for value in argv) and
                    argv.count("-c") == int(compile_only) and argv.count("-o") == 1 and
                    "-S" not in argv and "-E" not in argv, "header profile requires one explicit compile")
            require(not any(value.startswith(("-isystem", "-idirafter", "-iquote", "-iprefix", "-iwithprefix",
                "-isysroot", "--sysroot", "-include", "-imacros", "-nostdinc", "-Wp,", "-Xpreprocessor",
                "-fno-canonical-system-headers", "-fpreprocessed", "-E")) for value in argv[1:]),
                "header arguments override qualified include routing")
            extra = ["-nostdinc"] + (["-nostdinc++"] if self.phase.language == "c++" else [])
            extra += [f"-ffile-prefix-map={self.root}="]
            for path in self.include_roots: extra += ["-isystem", str(self.root / path.lstrip("/"))]
            # -nostdinc also suppresses GCC's implicit predefinition header.
            # Restore its contents explicitly from the retained view; output
            # equality alone may miss lost feature macros in a CMake probe.
            # This applies to C++ as well as C: omp.h alone does not restore it.
            extra += ["-include", str(self.root / C_PREDEFINITION_HEADER.lstrip("/"))]
            return [argv[0], *extra, *argv[1:]]
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-compiler-headers/v3" if self.phase.language == "c++" else "leopard2-v19-compiler-headers/v2",
            **({"implicit_predefinition_sealed": True} if self.phase.language == "c++" else {}),
            "root": str(self.root), "language": self.phase.language,
            "include_roots": list(self.include_roots), "files": [dict(pin, descriptor=self._files[path].executable_descriptor,
                seals=self._files[path].executable_record()["seals"]) for path, pin in self._pins.items()],
            "header_bytes": sum(row["size"] for row in self._pins.values()),
            "source_directories": len(self._source_dirs), "view_directories": len(self._view_dirs),
            "shared_source_guard": True, "declared_headers_sealed": True,
            "compiler_data_owned": False, "full_header_read_closure_owned": False,
            "negative_search_closure_owned": False, "fresh_build_recipe_integrated": False,
            "live_acquisition_armed": False, "benchmark_executed": False})
