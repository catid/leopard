#!/usr/bin/python3
"""Qualified sealed startup/link-data policy; leopard-79h.38.5.4.8.2.2.2.3.

Preserves original script bytes in a private sysroot, not arbitrary linker
scripts or complete host search closure. No benchmark/acquisition integration.
"""
from __future__ import annotations
import copy
import importlib.util
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_compiler_headers.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("linker input dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_link_input_headers", dependency)
header_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(header_module)
runtime, builder, provenance, require = header_module.runtime, header_module.builder, header_module.provenance, header_module.require
GCC_ROOT = "/usr/lib/gcc/x86_64-linux-gnu/13/"
SYSTEM_ROOT = "/usr/lib/x86_64-linux-gnu/"
LINK_INPUT_PATHS = tuple(GCC_ROOT + name for name in ("crtbeginS.o", "crtendS.o", "libgcc.a", "libgcc_s.so")) + tuple(
    SYSTEM_ROOT + name for name in ("Scrt1.o", "crti.o", "crtn.o", "ld-linux-x86-64.so.2", "libc.so", "libc.so.6",
        "libc_nonshared.a", "libgcc_s.so.1", "libgomp.so.1.0.0", "libm.so", "libm.so.6", "libmvec.so.1",
        "libpthread.a", "libstdc++.so.6.0.33"))
OPENMP_LINK_INPUT_PATHS = tuple(GCC_ROOT + name for name in
    ("libgomp.spec", "crtoffloadbegin.o", "crtoffloadend.o"))
MAX_FILES, MAX_DIRECTORIES, MAX_FILE_BYTES, MAX_TOTAL_BYTES = 64, 64, 8 << 20, 32 << 20


def link_pins(rows, *, openmp=False):
    require(type(openmp) is bool, "OpenMP link selection is not boolean")
    expected = LINK_INPUT_PATHS + (OPENMP_LINK_INPUT_PATHS if openmp else ())
    require(type(rows) is list and len(rows) == len(expected) <= MAX_FILES, "link input count differs")
    pins = {}
    for row in rows:
        require(type(row) is dict and set(row) == {"path", "sha256", "size"} and
                type(row["size"]) is int and 0 < row["size"] <= MAX_FILE_BYTES and
                type(row["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]), "link data pin differs")
        path = header_module.checked_path(row["path"])
        require(path not in pins and str(path) in expected, "duplicate or unqualified link input")
        pins[path] = copy.deepcopy(row)
    require(set(map(str, pins)) == set(expected) and sum(row["size"] for row in pins.values()) <= MAX_TOTAL_BYTES,
            "link input coverage or byte total differs")
    return pins


def link_aliases(pins):
    result = {}
    for path in pins:
        names = {path.relative_to("/"), Path("gcc-prefix") / path.name}
        if str(path).startswith(SYSTEM_ROOT): names.add(Path("lib/x86_64-linux-gnu") / path.name)
        if path.name == "ld-linux-x86-64.so.2": names.add(Path("lib64/ld-linux-x86-64.so.2"))
        if path.name == "libstdc++.so.6.0.33": names.add(Path("gcc-prefix/libstdc++.so"))
        if path.name == "libgomp.so.1.0.0": names.add(Path("gcc-prefix/libgomp.so"))
        for relative in sorted(names):
            require(relative not in result, "duplicate link data alias")
            result[relative] = path
    return result


class LinkerInputs(header_module._PinnedInputView):
    """Retain GCC13 link data with explicit OpenMP and C++ configuration profiles.

    The original scripts remain byte-for-byte unchanged. The private GCC prefix
    supplies startup files and libraries; its enclosing sysroot supplies the
    absolute paths inside those scripts. The helper prefix always comes first.
    This inventory is newly observed, not part of the original preflight pins.
    """
    def __init__(self, inventory, pins, *, openmp=False, cpp_configuration=False, _file_factory=None):
        require(inventory.phase.language in ("c", "c++"), "link input profile requires C or C++")
        require(type(openmp) is bool, "OpenMP link profile requires an explicit boolean selection")
        require(type(cpp_configuration) is bool, "C++ configuration link profile requires an explicit boolean selection")
        require(not cpp_configuration or (inventory.phase.language == "c++" and not openmp),
                "C++ configuration link profile requires plain C++ without OpenMP")
        self.openmp = openmp
        self.cpp_configuration = cpp_configuration
        selected = link_pins(pins, openmp=openmp)
        super().__init__(inventory, selected, (), link_aliases(selected),
            limits=(MAX_FILES, MAX_DIRECTORIES, MAX_FILE_BYTES), root_prefix="v19-link-inputs-", _file_factory=_file_factory)

    def arguments(self, argv, helper_descriptor):
        self.validate_current()
        try:
            # Validate caller arguments before inserting the one allowed data
            # -B prefix; it contains no helper executable roles.
            effective = runtime.compiler.driver_arguments(argv, self.phase.logical_driver, helper_descriptor)
            require(len(argv) <= 480 and argv.count("-o") == 1 and
                    all(value not in argv for value in ("-c", "-E", "-S")), "link profile requires one explicit link")
            require(not any(value.startswith(("--sysroot", "-isysroot", "-Wl,", "-Xlinker", "-L", "-l",
                "-nostdlib", "-nodefaultlibs", "-nostartfiles", "-static", "-shared", "-r", "-T", "-u"))
                for value in argv[1:]), "link arguments override qualified input routing")
            replacements = {GCC_ROOT + "libgomp.so": "libgomp.so", SYSTEM_ROOT + "libpthread.a": "libpthread.a"}
            options = argv[1:]
            if self.openmp:
                require(argv.count("-fopenmp") == 1, "OpenMP link profile requires exactly one -fopenmp")
                options = [value for value in options if value != "-fopenmp"]
            require(not any(value.startswith(("-fopenmp", "-fno-openmp", "-foffload", "-fno-offload",
                        "-fopenacc", "-pthread", "-pg")) or value == "-p" for value in options),
                    "link flags differ from selected default or OpenMP profile")
            if self.phase.language == "c++" and not self.openmp and not self.cpp_configuration:
                require(all(argv.count(path) == 1 for path in replacements), "qualified explicit link libraries differ")
            else:
                require(not any(path in argv for path in replacements), "probe profile uses implicit libraries only")
            prefix = self.root / "gcc-prefix"
            effective = [effective[0], effective[1], "-B" + str(prefix) + "/", "--sysroot=" + str(self.root), *effective[2:]]
            return [str(prefix / replacements[value]) if value in replacements else value for value in effective]
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": ("leopard2-v19-linker-inputs/v4" if self.cpp_configuration else
                "leopard2-v19-linker-inputs/v3" if self.openmp else "leopard2-v19-linker-inputs/v2"),
            **({"openmp_link_enabled": True} if self.openmp else {}),
            **({"cpp_configuration_link_enabled": True} if self.cpp_configuration else {}),
            "root": str(self.root), "language": self.phase.language,
            "prefix": str(self.root / "gcc-prefix"), "files": [dict(pin, descriptor=self._files[path].executable_descriptor,
                seals=self._files[path].executable_record()["seals"]) for path, pin in self._pins.items()],
            "mappings": {str(path): target for path, target in self._mappings.items()},
            "input_bytes": sum(row["size"] for row in self._pins.values()), "source_directories": len(self._source_dirs),
            "view_directories": len(self._view_dirs), "declared_link_inputs_sealed": True, "original_script_bytes_preserved": True,
            "compiler_data_owned": False, "full_link_read_closure_owned": False, "negative_search_closure_owned": False,
            "fresh_build_recipe_integrated": False, "live_acquisition_armed": False, "benchmark_executed": False})
