#!/usr/bin/env python3
"""Live generation-3 K65 arming boundary.

This module is the only generation-3 component allowed to touch host identity,
sealed build authorities, executable files, an evidence lane, or a child
process.  It deliberately has no CLI timing mode and exposes one production
entry point.  Only the live capability returned after an atomically durable
attempt directory may launch the fixed sequential campaign.

The exact-source boundary assumes a trusted CPython interpreter and standard
library state at bootstrap.  It closes and continuously revalidates the
reviewed repository controller graph, not arbitrary preseeded stdlib objects
or hostile code already executing inside the caller's process.  The trusted
caller also does not asynchronously mutate process-global state such as umask
while an acquisition or retained campaign method is executing.
"""

from __future__ import annotations

import _frozen_importlib
import _frozen_importlib_external
import copy
import ctypes
import errno
import fcntl
import functools
import hashlib
import importlib.machinery
import importlib.util
import os
from pathlib import Path
import platform
import re
import resource
import secrets
import selectors
import signal
import stat
import struct
import subprocess
import sys
import threading
import time
import types
from typing import Any, Callable, Mapping, NoReturn, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
TOOLS = REPO_ROOT / "tools"
MAIN_COMPARE = HERE.parent / "main_compare"


_EXECUTING_MODULE_CODE = sys._getframe().f_code
_BOOTSTRAP_MARKER = globals().get("__k65_gen3_exact_source_bootstrap__")
if (importlib.util.spec_from_file_location is not
        _frozen_importlib_external.spec_from_file_location or
        importlib.util.module_from_spec is not
        _frozen_importlib.module_from_spec or
        importlib.machinery.PathFinder is not
        _frozen_importlib_external.PathFinder):
    raise RuntimeError("Python import helpers were polluted before bootstrap")
_ORIGINAL_SPEC_FROM_FILE_LOCATION = \
    _frozen_importlib_external.spec_from_file_location
_ORIGINAL_MODULE_FROM_SPEC = _frozen_importlib.module_from_spec
_MODULE_SPEC_TYPE = _frozen_importlib.ModuleSpec
_MODULE_TYPE = type(sys)
_BOOTSTRAP_MODULE_REGISTRY = sys.modules
if type(_BOOTSTRAP_MODULE_REGISTRY) is not dict or any(
        type(alias) is not str
        for alias in tuple(_BOOTSTRAP_MODULE_REGISTRY.keys())):
    raise RuntimeError("Python module registry has unsafe keys")
_EXACT_IMPORT_LOCK = threading.RLock()
_FORK_BARRIER = threading.RLock()
_BOOTSTRAP_PID = os.getpid()
_ACTIVE_CAMPAIGNS: dict[int, Any] = {}
_CONTROLLED_FORK_PERMIT: list[tuple[int, int, int] | None] = [None]
_LOADED_CONTROLLER_SOURCE_SHA256: dict[Path, str] = {}
_EXACT_CONTROLLER_MODULES: dict[str, tuple[Any, Path, str]] = {}
_EXACT_CONTROLLER_LOADS: list[tuple[str, Any, Path, str]] = []
_MISSING_MODULE_ALIAS = object()
_LIVE_ARMER_MODULE = _BOOTSTRAP_MODULE_REGISTRY.get(__name__)
_MAX_CONTROLLER_SOURCE_BYTES = 16 * 1024 * 1024
_REQUIRED_CONTROLLER_RELATIVE_PATHS = (
    "experiments/leopard2/main_compare/pair_qualification_contract.py",
    "experiments/leopard2/main_compare/pair_qualification_acquire.py",
    "experiments/leopard2/main_compare/pair_qualification_bridge_contract.py",
    "experiments/leopard2/main_compare/pair_qualification_verify.py",
    "experiments/leopard2/main_compare/git_capture.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_preregistration.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_exact_source_bootstrap.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_gen3_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_execution_contract.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_gen3_acquire.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k66r16_b64_tail_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k9r5_b1024_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k5r5_b64_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/run_t8_two_block_abba.py",
    "experiments/leopard2/main_compare/run_abba.py",
    "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py",
    "tools/leopard2_build_provenance.py",
    "tools/leopard2_exact_main_baseline_acquire.py",
    "tools/leopard2_exact_main_baseline.py",
    "tools/leopard2_exact_main_baseline_record.py",
    "tools/leopard2_exact_main_baseline_verifier.py",
)
_REQUIRED_CONTROLLER_PATHS = frozenset(
    (REPO_ROOT / relative).resolve(strict=True)
    for relative in _REQUIRED_CONTROLLER_RELATIVE_PATHS)


def _read_exact_source_bytes(path: Path) -> tuple[bytes, str]:
    expected = path.resolve(strict=True)
    descriptor = os.open(
        expected, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        pathname = os.stat(expected, follow_symlinks=False)
        if not (stat.S_ISREG(before.st_mode) and
                before.st_uid == os.geteuid() and before.st_nlink == 1 and
                0 < before.st_size <= _MAX_CONTROLLER_SOURCE_BYTES and
                (before.st_dev, before.st_ino) ==
                    (pathname.st_dev, pathname.st_ino)):
            raise RuntimeError(f"controller source is unsafe: {expected}")
        chunks: list[bytes] = []
        offset = 0
        while offset < before.st_size:
            block = os.pread(
                descriptor, min(1 << 20, before.st_size - offset), offset)
            if not block:
                raise RuntimeError(
                    f"controller source read made no progress: {expected}")
            chunks.append(block)
            offset += len(block)
        data = b"".join(chunks)
        after = os.fstat(descriptor)
        pathname_after = os.stat(expected, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_mode, value.st_uid,
            value.st_nlink, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns)
        if not (len(data) == before.st_size and
                identity(before) == identity(after) == identity(pathname_after)):
            raise RuntimeError(
                f"controller source changed while loading: {expected}")
        return data, hashlib.sha256(data).hexdigest()
    finally:
        os.close(descriptor)


class _ExactSourceLoader:
    def __init__(self, name: str, path: Path) -> None:
        self.name = name
        self.path = path.resolve(strict=True)

    def create_module(self, unused_specification: Any) -> None:
        return None

    def get_filename(self, unused_name: str) -> str:
        return str(self.path)

    def is_package(self, unused_name: str) -> bool:
        return False

    def exec_module(self, module: Any) -> None:
        if (globals().get("_frozen_controller_ledger") is not None or
                type(_EXACT_CONTROLLER_MODULES) is not dict or
                type(_EXACT_CONTROLLER_LOADS) is not list):
            raise RuntimeError(
                "exact controller ledger is frozen against further loads")
        if type(module) is not _MODULE_TYPE:
            raise RuntimeError("exact controller module type changed")
        specification = vars(module).get("__spec__")
        if (type(specification) is not _MODULE_SPEC_TYPE or
                specification.name != self.name or
                specification.loader is not self or
                specification.origin != str(self.path) or
                specification.submodule_search_locations is not None):
            raise RuntimeError("exact controller module specification changed")
        data, digest = _read_exact_source_bytes(self.path)
        module.__file__ = str(self.path)
        module.__cached__ = None
        module.__exact_source_sha256__ = digest
        code = compile(
            data, str(self.path), "exec", dont_inherit=True,
            optimize=sys.flags.optimize)
        exec(code, module.__dict__)
        _EXACT_CONTROLLER_MODULES[self.name] = (module, self.path, digest)
        _EXACT_CONTROLLER_LOADS.append(
            (self.name, module, self.path, digest))
        _LOADED_CONTROLLER_SOURCE_SHA256[self.path] = digest


def _new_exact_specification(name: str, path: Path) -> Any:
    if type(name) is not str or not name:
        raise RuntimeError("exact controller module name is invalid")
    expected = path.resolve(strict=True)
    loader = _ExactSourceLoader(name, expected)
    specification = _MODULE_SPEC_TYPE(
        name, loader, origin=str(expected), is_package=False)
    specification._set_fileattr = True
    specification.cached = None
    if (type(specification) is not _MODULE_SPEC_TYPE or
            specification.name != name or
            specification.loader is not loader or
            specification.origin != str(expected) or
            specification.submodule_search_locations is not None):
        raise RuntimeError("cannot construct exact controller specification")
    return specification


def _exact_module_from_spec(specification: Any) -> Any:
    if (type(specification) is not _MODULE_SPEC_TYPE or
            type(specification.name) is not str or
            type(specification.loader) is not _ExactSourceLoader or
            specification.loader.name != specification.name or
            specification.origin != str(specification.loader.path) or
            specification.submodule_search_locations is not None):
        raise RuntimeError("non-exact controller specification was rejected")
    module = _MODULE_TYPE(specification.name)
    module.__loader__ = specification.loader
    module.__package__ = specification.parent
    module.__spec__ = specification
    module.__file__ = str(specification.loader.path)
    module.__cached__ = None
    return module


class _ExactRepositoryFinder:
    def find_spec(
        self, fullname: str, path: Any = None, target: Any = None,
    ) -> Any:
        specification = importlib.machinery.PathFinder.find_spec(
            fullname, path, target)
        if specification is None or specification.origin is None:
            return specification
        try:
            resolved = Path(specification.origin).resolve(strict=True)
        except (OSError, TypeError, ValueError):
            return specification
        if not (resolved.is_relative_to(REPO_ROOT) and
                resolved.suffix == ".py"):
            return specification
        if resolved not in _REQUIRED_CONTROLLER_PATHS:
            raise RuntimeError(
                f"unbound repository Python dependency: {resolved}")
        return _new_exact_specification(fullname, resolved)


_EXACT_REPOSITORY_FINDER = _ExactRepositoryFinder()


def _exact_spec_from_file_location(
    name: str, location: Any, *arguments: Any, **keywords: Any,
) -> Any:
    try:
        resolved = Path(location).resolve(strict=True)
        in_repository = resolved.is_relative_to(REPO_ROOT)
    except (OSError, TypeError, ValueError):
        in_repository = False
        resolved = Path("/")
    if (not in_repository or resolved.suffix != ".py" or
            resolved not in _REQUIRED_CONTROLLER_PATHS):
        raise RuntimeError(f"unbound Python dependency: {location}")
    if (arguments or set(keywords) - {
            "loader", "submodule_search_locations"} or
            keywords.get("loader") is not None or
            keywords.get("submodule_search_locations") is not None):
        raise RuntimeError("custom controller module specification rejected")
    return _new_exact_specification(name, resolved)


class _ExactImportScope:
    def __enter__(self) -> None:
        _EXACT_IMPORT_LOCK.acquire()
        try:
            self.alias_bindings: list[tuple[str, Any, Any]] = []
            self.dynamic_alias_bindings: list[tuple[str, Any]] = []
            self.modules_object: dict[str, Any] | None = None
            self.modules_snapshot: dict[str, Any] = {}
            self.meta_path_object: list[Any] | None = None
            self.meta_path_snapshot: tuple[Any, ...] = ()
            self.previous_spec: Any = _MISSING_MODULE_ALIAS
            self.previous_module: Any = _MISSING_MODULE_ALIAS
            self.live_surface_bindings: tuple[tuple[Any, str, Any], ...] = ()
            self.live_surface_getter: Any = None
            self.controller_graph_getter: Any = None
            self.controller_graph: tuple[Any, ...] = ()
            self.child_environment_getter: Any = None
            self.child_environment: Any = None
            self.live_execution_getter: Any = None
            self.live_execution_graph: tuple[Any, ...] = ()
            modules_object = sys.modules
            meta_path_object = sys.meta_path
            if type(modules_object) is not dict or \
                    type(meta_path_object) is not list:
                raise RuntimeError("Python import registries were replaced")
            self.modules_object = modules_object
            registry_keys = tuple(modules_object.keys())
            if any(type(alias) is not str for alias in registry_keys):
                raise RuntimeError(
                    "Python module registry contains a non-string key")
            self.modules_snapshot = dict(modules_object)
            self.meta_path_object = meta_path_object
            self.meta_path_snapshot = tuple(meta_path_object)
            self.previous_spec = importlib.util.spec_from_file_location
            self.previous_module = importlib.util.module_from_spec
            if not (self.previous_spec is
                        _ORIGINAL_SPEC_FROM_FILE_LOCATION or
                    self.previous_spec is _exact_spec_from_file_location) or \
                    not (self.previous_module is
                            _ORIGINAL_MODULE_FROM_SPEC or
                         self.previous_module is _exact_module_from_spec):
                raise RuntimeError("Python source loader was replaced")
            surface_getter = globals().get("_frozen_live_public_surface")
            if surface_getter is not None:
                self.live_surface_getter = surface_getter
                surface_bindings = surface_getter()
                if (_LIVE_PUBLIC_CALLABLE_IDENTITIES is not
                        surface_bindings):
                    raise RuntimeError(
                        "live public callable authority binding changed")
                self.live_surface_bindings = surface_bindings
                changed_surface = [
                    (owner, name, expected)
                    for owner, name, expected in surface_bindings
                    if vars(owner).get(name, _MISSING_MODULE_ALIAS) is not
                        expected
                ]
                if changed_surface:
                    for owner, name, expected in changed_surface:
                        setattr(owner, name, expected)
                    raise RuntimeError("live public callable identity changed")
                graph_getter = globals().get(
                    "_frozen_controller_execution_graph")
                if graph_getter is None:
                    raise RuntimeError(
                        "frozen controller graph getter is absent")
                graph = graph_getter()
                if type(graph) is not tuple or len(graph) != 5:
                    raise RuntimeError("frozen controller graph differs")
                (armer_modules, controller_module_edges,
                 executable_edges, class_member_edges,
                 function_states) = graph
                if (_ARMER_MODULE_REFERENCE_IDENTITIES is not armer_modules or
                        _CONTROLLER_MODULE_REFERENCE_IDENTITIES is not
                            controller_module_edges or
                        _CONTROLLER_EXECUTABLE_REFERENCE_IDENTITIES is not
                            executable_edges or
                        _CONTROLLER_CLASS_MEMBER_REFERENCE_IDENTITIES is not
                            class_member_edges or
                        _CONTROLLER_FUNCTION_EXECUTION_STATES is not
                            function_states):
                    raise RuntimeError(
                        "frozen controller graph binding changed")
                for name, expected in armer_modules:
                    if globals().get(name) is not expected:
                        raise RuntimeError(
                            f"live armer module reference changed: {name}")
                for owner, name, expected in controller_module_edges:
                    if vars(owner).get(name) is not expected:
                        raise RuntimeError(
                            "controller module reference changed: "
                            f"{owner.__name__}.{name}")
                _validate_function_execution_states(
                    function_states, "controller")
                self.controller_graph_getter = graph_getter
                self.controller_graph = graph
                environment_getter = globals().get(
                    "_frozen_child_environment")
                if environment_getter is None:
                    raise RuntimeError(
                        "frozen child environment getter is absent")
                frozen_environment = environment_getter()
                if (CHILD_ENVIRONMENT is not frozen_environment or
                        type(frozen_environment) is not
                            types.MappingProxyType):
                    raise RuntimeError(
                        "frozen child environment binding changed")
                self.child_environment_getter = environment_getter
                self.child_environment = frozen_environment
                live_execution_getter = globals().get(
                    "_frozen_live_execution_graph")
                if live_execution_getter is None:
                    raise RuntimeError(
                        "frozen live execution graph getter is absent")
                live_execution_graph = live_execution_getter()
                if (type(live_execution_graph) is not tuple or
                        len(live_execution_graph) != 3):
                    raise RuntimeError(
                        "frozen live execution graph differs")
                (live_executables, live_class_members,
                 live_function_states) = live_execution_graph
                if (_LIVE_ARMER_EXECUTABLE_REFERENCE_IDENTITIES is not
                        live_executables or
                        _LIVE_ARMER_CLASS_MEMBER_REFERENCE_IDENTITIES is not
                            live_class_members or
                        _LIVE_ARMER_FUNCTION_EXECUTION_STATES is not
                            live_function_states):
                    raise RuntimeError(
                        "frozen live execution graph binding changed")
                for owner, name, expected in live_executables:
                    if vars(owner).get(name, _MISSING_MODULE_ALIAS) is not \
                            expected:
                        raise RuntimeError(
                            f"live armer executable changed: {name}")
                for owner, name, expected in live_class_members:
                    if vars(owner).get(name, _MISSING_MODULE_ALIAS) is not \
                            expected:
                        raise RuntimeError(
                            "live armer class member changed: "
                            f"{owner.__qualname__}.{name}")
                _validate_function_execution_states(
                    live_function_states, "live armer")
                self.live_execution_getter = live_execution_getter
                self.live_execution_graph = live_execution_graph
            ledger_getter = globals().get("_frozen_controller_ledger")
            if ledger_getter is None:
                if (type(_EXACT_CONTROLLER_LOADS) is not list or
                        type(_EXACT_CONTROLLER_MODULES) is not dict):
                    raise RuntimeError(
                        "controller ledger froze without its authority")
                self.controller_ledger_frozen = False
                self.controller_loads = _EXACT_CONTROLLER_LOADS
                self.controller_modules = _EXACT_CONTROLLER_MODULES
            else:
                frozen_loads, frozen_modules = ledger_getter()
                if (type(frozen_loads) is not tuple or
                        type(frozen_modules) is not types.MappingProxyType or
                        _EXACT_CONTROLLER_LOADS is not frozen_loads or
                        _EXACT_CONTROLLER_MODULES is not frozen_modules):
                    raise RuntimeError(
                        "frozen controller ledger binding changed")
                self.controller_ledger_frozen = True
                self.controller_loads = frozen_loads
                self.controller_modules = frozen_modules
            if _LIVE_ARMER_MODULE is None:
                raise RuntimeError("live armer module identity is absent")
            protected_aliases: dict[str, Any] = {
                __name__: _LIVE_ARMER_MODULE,
            }
            protected_aliases.update({
                alias: module
                for alias, (module, unused_path, unused_digest) in tuple(
                    self.controller_modules.items())
            })
            for alias, module in tuple(modules_object.items()):
                if alias in protected_aliases or module is None:
                    continue
                if (type(module) is not _MODULE_TYPE or
                        _module_file_marker_is_unsafe(module) or
                        _controller_module_path(module) is not None):
                    protected_aliases[alias] = _MISSING_MODULE_ALIAS
            for alias, module in protected_aliases.items():
                previous_module = modules_object.get(
                    alias, _MISSING_MODULE_ALIAS)
                self.alias_bindings.append(
                    (alias, module, previous_module))
                if module is _MISSING_MODULE_ALIAS:
                    modules_object.pop(alias, None)
                elif previous_module is not module:
                    modules_object[alias] = module
            meta_path_object[:] = [
                _EXACT_REPOSITORY_FINDER,
                *(entry for entry in meta_path_object
                  if entry is not _EXACT_REPOSITORY_FINDER),
            ]
            importlib.util.spec_from_file_location = \
                _exact_spec_from_file_location
            importlib.util.module_from_spec = _exact_module_from_spec
        except BaseException:
            try:
                previous_spec = getattr(
                    self, "previous_spec", _MISSING_MODULE_ALIAS)
                if previous_spec is not _MISSING_MODULE_ALIAS:
                    importlib.util.spec_from_file_location = previous_spec
                previous_module = getattr(
                    self, "previous_module", _MISSING_MODULE_ALIAS)
                if previous_module is not _MISSING_MODULE_ALIAS:
                    importlib.util.module_from_spec = previous_module
                meta_path_object = getattr(
                    self, "meta_path_object", None)
                if meta_path_object is not None:
                    sys.meta_path = meta_path_object
                    meta_path_object[:] = getattr(
                        self, "meta_path_snapshot", ())
                modules_object = getattr(self, "modules_object", None)
                if modules_object is not None:
                    sys.modules = modules_object
                    for alias, unused_module, previous_module in reversed(
                            getattr(self, "alias_bindings", ())):
                        if previous_module is _MISSING_MODULE_ALIAS:
                            modules_object.pop(alias, None)
                        else:
                            modules_object[alias] = previous_module
            finally:
                _EXACT_IMPORT_LOCK.release()
            raise

    def __exit__(self, unused_type: Any, unused_value: Any,
                 unused_traceback: Any) -> None:
        error: BaseException | None = None
        force_registry_rebuild = False
        cleanup_errors: list[BaseException] = []
        try:
            if importlib.util.spec_from_file_location is not \
                    _exact_spec_from_file_location:
                error = RuntimeError("Python source loader changed in exact scope")
            if (importlib.util.module_from_spec is not
                    _exact_module_from_spec and error is None):
                error = RuntimeError(
                    "Python module constructor changed in exact scope")
            if (sys.meta_path is not self.meta_path_object or
                    not self.meta_path_object or
                    self.meta_path_object[0] is not
                        _EXACT_REPOSITORY_FINDER or
                    sum(entry is _EXACT_REPOSITORY_FINDER
                        for entry in self.meta_path_object) != 1):
                if error is None:
                    error = RuntimeError(
                        "exact repository finder identity or order changed")
            if sys.modules is not self.modules_object and error is None:
                error = RuntimeError("Python module registry changed in scope")
            if (self.controller_ledger_frozen and
                    (_EXACT_CONTROLLER_LOADS is not self.controller_loads or
                     _EXACT_CONTROLLER_MODULES is not
                        self.controller_modules) and
                    error is None):
                error = RuntimeError(
                    "frozen controller ledger binding changed in scope")
            if self.live_surface_getter is not None:
                if (globals().get("_frozen_live_public_surface") is not
                        self.live_surface_getter or
                        _LIVE_PUBLIC_CALLABLE_IDENTITIES is not
                        self.live_surface_bindings):
                    if error is None:
                        error = RuntimeError(
                            "live public callable authority changed in scope")
                changed_surface = [
                    (owner, name, expected)
                    for owner, name, expected in self.live_surface_bindings
                    if vars(owner).get(name, _MISSING_MODULE_ALIAS) is not
                        expected
                ]
                if changed_surface:
                    for owner, name, expected in changed_surface:
                        setattr(owner, name, expected)
                    if error is None:
                        error = RuntimeError(
                            "live public callable identity changed in scope")
            if self.controller_graph_getter is not None:
                current_graph = self.controller_graph_getter()
                if (globals().get("_frozen_controller_execution_graph") is not
                        self.controller_graph_getter or
                        type(current_graph) is not tuple or
                        len(current_graph) != len(self.controller_graph) or
                        any(current is not retained for current, retained in
                            zip(current_graph, self.controller_graph))):
                    if error is None:
                        error = RuntimeError(
                            "frozen controller graph changed in scope")
                else:
                    (armer_modules, controller_module_edges, unused_a,
                     unused_b, function_states) = self.controller_graph
                    for name, expected in armer_modules:
                        if globals().get(name) is not expected and error is None:
                            error = RuntimeError(
                                f"live armer module reference changed: {name}")
                    for owner, name, expected in controller_module_edges:
                        if (vars(owner).get(name) is not expected and
                                error is None):
                            error = RuntimeError(
                                "controller module reference changed: "
                                f"{owner.__name__}.{name}")
                    if error is None:
                        try:
                            _validate_function_execution_states(
                                function_states, "controller")
                        except BaseException as state_error:
                            error = state_error
            if self.child_environment_getter is not None and (
                    globals().get("_frozen_child_environment") is not
                        self.child_environment_getter or
                    CHILD_ENVIRONMENT is not self.child_environment or
                    self.child_environment_getter() is not
                        self.child_environment) and error is None:
                error = RuntimeError(
                    "frozen child environment changed in scope")
            if self.live_execution_getter is not None:
                current_live_graph = self.live_execution_getter()
                if (globals().get("_frozen_live_execution_graph") is not
                        self.live_execution_getter or
                        type(current_live_graph) is not tuple or
                        len(current_live_graph) !=
                            len(self.live_execution_graph) or
                        any(current is not retained
                            for current, retained in zip(
                                current_live_graph,
                                self.live_execution_graph))):
                    if error is None:
                        error = RuntimeError(
                            "frozen live execution graph changed in scope")
                else:
                    (live_executables, live_class_members,
                     live_function_states) = self.live_execution_graph
                    for owner, name, expected in live_executables:
                        if (vars(owner).get(
                                name, _MISSING_MODULE_ALIAS) is not expected and
                                error is None):
                            error = RuntimeError(
                                f"live armer executable changed: {name}")
                    for owner, name, expected in live_class_members:
                        if (vars(owner).get(
                                name, _MISSING_MODULE_ALIAS) is not expected and
                                error is None):
                            error = RuntimeError(
                                "live armer class member changed: "
                                f"{owner.__qualname__}.{name}")
                    if error is None:
                        try:
                            _validate_function_execution_states(
                                live_function_states, "live armer")
                        except BaseException as state_error:
                            error = state_error
            registry_keys = tuple(self.modules_object.keys())
            if any(type(alias) is not str for alias in registry_keys):
                force_registry_rebuild = True
                if error is None:
                    error = RuntimeError(
                        "Python module registry gained a non-string key")
            elif any(
                    module is not None and (
                        type(module) is not _MODULE_TYPE or
                        _module_file_marker_is_unsafe(module))
                    for module in self.modules_object.values()):
                force_registry_rebuild = True
                if error is None:
                    error = RuntimeError(
                        "Python module registry gained a non-plain module")
            else:
                for alias, module, unused_previous in self.alias_bindings:
                    current = self.modules_object.get(
                        alias, _MISSING_MODULE_ALIAS)
                    if current is not module and error is None:
                        label = (
                            "unreviewed controller alias appeared in scope"
                            if module is _MISSING_MODULE_ALIAS else
                            "exact controller alias changed in scope")
                        error = RuntimeError(f"{label}: {alias}")
                protected_names = {
                    alias for alias, unused_module, unused_previous
                    in self.alias_bindings
                }
                for alias, module in tuple(self.modules_object.items()):
                    if (alias in protected_names or module is None or
                            _controller_module_path(module) is None):
                        continue
                    if not self.controller_ledger_frozen:
                        ledger = self.controller_modules.get(alias)
                        if (ledger is not None and ledger[0] is module and
                                _controller_module_path(module) == ledger[1] and
                                vars(module).get(
                                    "__exact_source_sha256__") == ledger[2] and
                                any(load_alias == alias and
                                    load_module is module and
                                    load_path == ledger[1] and
                                    load_digest == ledger[2]
                                    for load_alias, load_module, load_path,
                                        load_digest in self.controller_loads)):
                            continue
                    previous_module = self.modules_snapshot.get(
                        alias, _MISSING_MODULE_ALIAS)
                    self.dynamic_alias_bindings.append(
                        (alias, previous_module))
                    if error is None:
                        error = RuntimeError(
                            f"new controller alias appeared in scope: {alias}")
        except BaseException as validation_error:
            force_registry_rebuild = True
            if error is None:
                error = validation_error
        finally:
            try:
                try:
                    importlib.util.spec_from_file_location = \
                        self.previous_spec
                    importlib.util.module_from_spec = self.previous_module
                except BaseException as cleanup_error:
                    cleanup_errors.append(cleanup_error)
                try:
                    sys.meta_path = self.meta_path_object
                    self.meta_path_object[:] = self.meta_path_snapshot
                except BaseException as cleanup_error:
                    cleanup_errors.append(cleanup_error)
                try:
                    sys.modules = self.modules_object
                    if force_registry_rebuild:
                        self.modules_object.clear()
                        self.modules_object.update(self.modules_snapshot)
                    else:
                        for alias, previous_module in reversed(
                                self.dynamic_alias_bindings):
                            if previous_module is _MISSING_MODULE_ALIAS:
                                self.modules_object.pop(alias, None)
                            else:
                                self.modules_object[alias] = previous_module
                        for alias, unused_module, previous_module in reversed(
                                self.alias_bindings):
                            if previous_module is _MISSING_MODULE_ALIAS:
                                self.modules_object.pop(alias, None)
                            else:
                                self.modules_object[alias] = previous_module
                except BaseException as cleanup_error:
                    cleanup_errors.append(cleanup_error)
                    try:
                        sys.modules = self.modules_object
                        self.modules_object.clear()
                        self.modules_object.update(self.modules_snapshot)
                    except BaseException as rebuild_error:
                        cleanup_errors.append(rebuild_error)
            finally:
                _EXACT_IMPORT_LOCK.release()
        if cleanup_errors:
            cleanup_summary = "; ".join(repr(item) for item in cleanup_errors)
            if error is None:
                error = RuntimeError(
                    "exact import scope cleanup failed: " + cleanup_summary)
            else:
                error.add_note(
                    "exact import scope cleanup also failed: " +
                    cleanup_summary)
        if error is not None and unused_type is None:
            raise error
        if error is not None and isinstance(unused_value, BaseException):
            unused_value.add_note(f"exact import scope also failed: {error!r}")


def _with_exact_source_imports(function: Callable[..., Any]) \
        -> Callable[..., Any]:
    @functools.wraps(function)
    def wrapped(*arguments: Any, **keywords: Any) -> Any:
        # A post-fork child can inherit both the process-global exact-import
        # RLock and campaign locks in their locked state.  Every foreign call
        # is terminal before either lock; calls without a campaign authority
        # are rejected as well.
        campaign_type = globals().get("ArmedCampaign")
        if os.getpid() != _BOOTSTRAP_PID:
            if (arguments and isinstance(campaign_type, type) and
                    isinstance(arguments[0], campaign_type)):
                arguments[0]._reject_foreign_process()
            error_type = globals().get("ArmingError", RuntimeError)
            raise error_type(
                "generation-3 live authority cannot be entered after fork")
        _FORK_BARRIER.acquire()
        try:
            result: Any = _MISSING_MODULE_ALIAS
            body_succeeded = False
            try:
                with _ExactImportScope():
                    result = function(*arguments, **keywords)
                    body_succeeded = True
            except BaseException as scope_error:
                # A successful acquisition transfers every retained capability
                # to its returned ArmedCampaign before this context exits.  If
                # the exit gate then detects mutation, revoke it before the
                # fork-serialization lease is released.
                campaign_type = globals().get("ArmedCampaign")
                campaign = None
                if body_succeeded and isinstance(campaign_type, type):
                    if isinstance(result, campaign_type):
                        campaign = result
                    elif (arguments and
                          isinstance(arguments[0], campaign_type)):
                        campaign = arguments[0]
                if campaign is not None:
                    try:
                        campaign.close()
                    except BaseException as cleanup_error:
                        scope_error.add_note(
                            "ArmedCampaign cleanup after exact-scope failure "
                            f"also failed: {cleanup_error!r}")
                raise
            campaign_type = globals().get("ArmedCampaign")
            if (body_succeeded and isinstance(campaign_type, type) and
                    isinstance(result, campaign_type)):
                try:
                    _register_active_campaign(result)
                except BaseException as registration_error:
                    try:
                        result.close()
                    except BaseException as cleanup_error:
                        registration_error.add_note(
                            "ArmedCampaign cleanup after registry failure "
                            f"also failed: {cleanup_error!r}")
                    raise
            return result
        finally:
            _FORK_BARRIER.release()
    return wrapped


def _controller_module_path(module: Any) -> Path | None:
    if type(module) is not _MODULE_TYPE:
        return None
    raw_path = vars(module).get("__file__")
    if raw_path is None:
        return None
    try:
        if type(raw_path) is str:
            normalized = raw_path
        elif type(raw_path) is bytes:
            normalized = os.fsdecode(raw_path)
        elif type(raw_path) is type(HERE):
            normalized = str(raw_path)
        else:
            return None
        if not normalized:
            return None
        resolved = Path(normalized).resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError):
        return None
    return resolved if resolved in _REQUIRED_CONTROLLER_PATHS else None


def _module_file_marker_is_unsafe(module: Any) -> bool:
    if type(module) is not _MODULE_TYPE:
        return True
    raw_path = vars(module).get("__file__")
    return (raw_path is not None and
            type(raw_path) not in (str, bytes, type(HERE)))


def _purge_preexisting_controller_aliases() -> None:
    module_registry = sys.modules
    if type(module_registry) is not dict:
        raise RuntimeError("Python module registry was replaced")
    aliases = tuple(module_registry.keys())
    if any(type(alias) is not str for alias in aliases):
        raise RuntimeError("Python module registry contains a non-string key")
    for alias, module in tuple(module_registry.items()):
        if (module is None or
                (alias == __name__ and module is _LIVE_ARMER_MODULE)):
            continue
        path = _controller_module_path(module)
        if path is not None:
            if module_registry.get(alias) is module:
                del module_registry[alias]


def _load(name: str, path: Path) -> Any:
    if globals().get("_frozen_controller_ledger") is not None:
        raise RuntimeError(
            "exact controller ledger is frozen against further loads")
    expected = path.resolve(strict=True)
    loaded = sys.modules.get(name)
    if loaded is not None:
        data, digest = _read_exact_source_bytes(expected)
        del data
        if (type(loaded) is not _MODULE_TYPE or
                _controller_module_path(loaded) != expected or
                vars(loaded).get("__exact_source_sha256__") != digest):
            raise RuntimeError(f"{name} came from another path")
        ledger = _EXACT_CONTROLLER_MODULES.get(name)
        if ledger is None or ledger != (loaded, expected, digest):
            raise RuntimeError(f"{name} bypassed the exact-source ledger")
        _LOADED_CONTROLLER_SOURCE_SHA256[expected] = digest
        return loaded
    with _ExactImportScope():
        specification = _new_exact_specification(name, expected)
        loader = specification.loader
        module = _exact_module_from_spec(specification)
        if (type(loader) is not _ExactSourceLoader or
                specification.loader is not loader or
                type(module) is not _MODULE_TYPE or
                module.__spec__ is not specification or
                module.__loader__ is not loader):
            raise RuntimeError(f"cannot construct exact module for {expected}")
        sys.modules[name] = module
        try:
            loader.exec_module(module)
        except BaseException:
            sys.modules.pop(name, None)
            raise
    if _controller_module_path(module) != expected:
        raise RuntimeError(f"{name} resolved to another path")
    return module


def _record_self_source_identity() -> None:
    path = Path(__file__).resolve(strict=True)
    data, digest = _read_exact_source_bytes(path)
    del data
    marker = _BOOTSTRAP_MARKER
    if (type(marker) is not dict or set(marker) != {
            "code", "source_path", "source_sha256", "bootstrap_path",
            "bootstrap_sha256", "bootstrap_direct_source"} or
            _LIVE_ARMER_MODULE is None or
            type(_LIVE_ARMER_MODULE) is not _MODULE_TYPE or
            sys.modules.get(__name__) is not _LIVE_ARMER_MODULE or
            vars(_LIVE_ARMER_MODULE) is not globals() or
            marker["code"] is not _EXECUTING_MODULE_CODE or
            marker["source_path"] != path or
            marker["source_sha256"] != digest or
            marker["bootstrap_direct_source"] is not True):
        raise RuntimeError(
            "generation-3 armer requires the exact-source bootstrap")
    bootstrap_path = Path(marker["bootstrap_path"]).resolve(strict=True)
    bootstrap_data, bootstrap_digest = _read_exact_source_bytes(bootstrap_path)
    del bootstrap_data
    if (bootstrap_path !=
            (HERE / "k65_gen3_exact_source_bootstrap.py").resolve(strict=True) or
            marker["bootstrap_sha256"] != bootstrap_digest):
        raise RuntimeError("generation-3 exact-source bootstrap changed")
    globals()["__exact_source_sha256__"] = digest
    _LOADED_CONTROLLER_SOURCE_SHA256[path] = digest
    _LOADED_CONTROLLER_SOURCE_SHA256[bootstrap_path] = bootstrap_digest


_record_self_source_identity()
_purge_preexisting_controller_aliases()


prereg = _load("leopard2_k65_gen3_prereg_for_live_armer",
               HERE / "k65_gen3_preregistration.py")
if tuple(prereg.REQUIRED_CONTROLLER_PATHS) != \
        _REQUIRED_CONTROLLER_RELATIVE_PATHS:
    raise RuntimeError(
        "bootstrap controller closure differs from preregistration")
execution = _load("leopard2_k65_gen3_execution_for_live_armer",
                  HERE / "k65_gen3_execution_contract.py")
plan_contract = _load("leopard2_k65_gen3_plan_for_live_armer",
                      HERE / "run_k65r65_b64_packed_terminal_gen3_abba.py")
build_provenance = _load("leopard2_build_provenance",
                         TOOLS / "leopard2_build_provenance.py")
git_capture = _load("git_capture", MAIN_COMPARE / "git_capture.py")
balanced_evidence_common = _load(
    "balanced_evidence_common",
    HERE.parent / "decoder_dispatch" / "balanced_evidence_common.py")
result_contract = _load(
    "leopard2_k65_result_contract_for_gen3_live_armer",
    HERE / "run_k65r65_b64_packed_terminal_abba.py")
pair_acquire = _load("leopard2_pair_acquire_for_k65_gen3_live_armer",
                     MAIN_COMPARE / "pair_qualification_acquire.py")
pair_bridge = _load("leopard2_pair_bridge_for_k65_gen3_live_armer",
                    MAIN_COMPARE / "pair_qualification_bridge_contract.py")
pair_verify = _load("leopard2_pair_verify_for_k65_gen3_live_armer",
                    MAIN_COMPARE / "pair_qualification_verify.py")
identity_contract = _load("leopard2_exact_main_baseline",
                          TOOLS / "leopard2_exact_main_baseline.py")
baseline_acquire = _load("leopard2_exact_main_acquire_for_k65_gen3_live_armer",
                         TOOLS / "leopard2_exact_main_baseline_acquire.py")
authority_record = _load("leopard2_exact_main_baseline_record",
                         TOOLS / "leopard2_exact_main_baseline_record.py")
authority_verifier = _load(
    "leopard2_exact_main_baseline_verifier",
    TOOLS / "leopard2_exact_main_baseline_verifier.py")
contract = prereg.contract


def _capture_function_execution_states(
    roots: Sequence[Any],
) -> tuple[tuple[Any, ...], ...]:
    """Snapshot every Python function reachable through reviewed callables."""
    observed: list[tuple[Any, ...]] = []
    seen: set[int] = set()

    def visit(value: Any) -> None:
        if isinstance(value, (staticmethod, classmethod)):
            visit(value.__func__)
            return
        if isinstance(value, property):
            for accessor in (value.fget, value.fset, value.fdel):
                if accessor is not None:
                    visit(accessor)
            return
        if isinstance(value, type):
            for member in vars(value).values():
                visit(member)
            return
        if type(value) is not types.FunctionType or id(value) in seen:
            return
        seen.add(id(value))
        defaults = value.__defaults__
        kwdefaults = value.__kwdefaults__
        closure = value.__closure__
        cells: list[tuple[Any, Any]] = []
        for cell in closure or ():
            try:
                content = cell.cell_contents
            except ValueError:
                content = _MISSING_MODULE_ALIAS
            cells.append((cell, content))
        observed.append((
            value, value.__code__, defaults,
            tuple(defaults or ()), kwdefaults,
            tuple(sorted((kwdefaults or {}).items())),
            value.__globals__, closure, tuple(cells),
        ))
        for default in defaults or ():
            visit(default)
        for default in (kwdefaults or {}).values():
            visit(default)
        for unused_cell, content in cells:
            visit(content)

    for root in roots:
        visit(root)
    return tuple(observed)


def _validate_function_execution_states(
    states: Sequence[tuple[Any, ...]], label: str,
) -> None:
    for index, state in enumerate(states):
        _require(type(state) is tuple and len(state) == 9,
                 f"{label} function state {index} is malformed")
        (function, code, defaults, default_values, kwdefaults,
         kwdefault_items, global_namespace, closure, cells) = state
        _require(
            type(function) is types.FunctionType and
            function.__code__ is code and
            function.__defaults__ is defaults and
            len(function.__defaults__ or ()) == len(default_values) and
            all(current is expected for current, expected in
                zip(function.__defaults__ or (), default_values)) and
            function.__kwdefaults__ is kwdefaults and
            set((function.__kwdefaults__ or {}).keys()) ==
                {name for name, unused in kwdefault_items} and
            all((function.__kwdefaults__ or {}).get(name) is expected
                for name, expected in kwdefault_items) and
            function.__globals__ is global_namespace and
            function.__closure__ is closure and
            len(function.__closure__ or ()) == len(cells),
            f"{label} function execution state changed")
        for current, (expected_cell, expected_content) in zip(
                function.__closure__ or (), cells):
            try:
                current_content = current.cell_contents
            except ValueError:
                current_content = _MISSING_MODULE_ALIAS
            _require(current is expected_cell and
                     current_content is expected_content,
                     f"{label} function closure changed")


# Freeze every module-valued import edge after the exact-source graph has
# finished loading.  A path check alone is insufficient here: Python code can
# replace an executed dependency with a module that has no ``__file__`` (or a
# forged one), while functions in the owner continue to resolve that global at
# call time.  Identity joins cover both reviewed controller edges and stdlib
# modules on which those controllers depend.
_ARMER_MODULE_REFERENCE_IDENTITIES = tuple(
    (name, value)
    for name, value in globals().items()
    if isinstance(value, types.ModuleType)
)
_CONTROLLER_MODULE_REFERENCE_IDENTITIES = tuple(
    (module, name, value)
    for module in dict.fromkeys(
        loaded_module for _, loaded_module, _, _ in _EXACT_CONTROLLER_LOADS)
    for name, value in vars(module).items()
    if isinstance(value, types.ModuleType)
)
_CONTROLLER_EXECUTABLE_REFERENCE_IDENTITIES = tuple(
    (module, name, value)
    for module in dict.fromkeys(
        loaded_module for _, loaded_module, _, _ in _EXACT_CONTROLLER_LOADS)
    for name, value in vars(module).items()
    if ((type(value) is types.FunctionType and
         value.__globals__ is vars(module)) or
        (isinstance(value, type) and
         getattr(value, "__module__", None) == module.__name__))
)
_CONTROLLER_CLASS_MEMBER_REFERENCE_IDENTITIES = tuple(
    (class_value, member_name, member)
    for module in dict.fromkeys(
        loaded_module for _, loaded_module, _, _ in _EXACT_CONTROLLER_LOADS)
    for class_value in vars(module).values()
    if (isinstance(class_value, type) and
        getattr(class_value, "__module__", None) == module.__name__)
    for member_name, member in vars(class_value).items()
    if (callable(member) or
        isinstance(member, (staticmethod, classmethod, property)))
)
_CONTROLLER_FUNCTION_EXECUTION_STATES = _capture_function_execution_states((
    *(value for unused_owner, unused_name, value in
      _CONTROLLER_EXECUTABLE_REFERENCE_IDENTITIES),
    *(value for unused_owner, unused_name, value in
      _CONTROLLER_CLASS_MEMBER_REFERENCE_IDENTITIES),
))


def _frozen_controller_execution_graph(
    armer_modules: tuple[tuple[str, Any], ...] =
        _ARMER_MODULE_REFERENCE_IDENTITIES,
    controller_modules: tuple[tuple[Any, str, Any], ...] =
        _CONTROLLER_MODULE_REFERENCE_IDENTITIES,
    executables: tuple[tuple[Any, str, Any], ...] =
        _CONTROLLER_EXECUTABLE_REFERENCE_IDENTITIES,
    class_members: tuple[tuple[Any, str, Any], ...] =
        _CONTROLLER_CLASS_MEMBER_REFERENCE_IDENTITIES,
    function_states: tuple[tuple[Any, ...], ...] =
        _CONTROLLER_FUNCTION_EXECUTION_STATES,
) -> tuple[tuple[tuple[str, Any], ...],
           tuple[tuple[Any, str, Any], ...],
           tuple[tuple[Any, str, Any], ...],
           tuple[tuple[Any, str, Any], ...],
           tuple[tuple[Any, ...], ...]]:
    """Return bootstrap-captured executable identities after ledger freeze."""
    return (armer_modules, controller_modules, executables, class_members,
            function_states)


# No controller source is loaded after bootstrap.  Convert the construction
# ledgers themselves into immutable authority before exposing either public
# workflow entry point; production scope validation must never authenticate a
# module merely because a caller appended a plausible-looking ledger row.
_EXACT_CONTROLLER_LOADS = tuple(_EXACT_CONTROLLER_LOADS)
_EXACT_CONTROLLER_MODULES = types.MappingProxyType(
    dict(_EXACT_CONTROLLER_MODULES))


def _frozen_controller_ledger(
    loads: tuple[tuple[str, Any, Path, str], ...] = _EXACT_CONTROLLER_LOADS,
    modules: Any = _EXACT_CONTROLLER_MODULES,
) -> tuple[tuple[tuple[str, Any, Path, str], ...], Any]:
    """Return bootstrap-captured identities, independent of global rebinding."""
    return loads, modules

GENERATION = 3
PATH_VARIANT_RELATIVE = \
    "artifacts/path-variant/leopard_main_benchmark"
PATH_VARIANT_RAW_SHA256 = \
    "0baae845bbf30d2b3b213c02501a31c2d15fd125965d898f7a100fa6d0ede46d"
PATH_VARIANT_SIZE = 1_175_456
CANONICAL_RAW_SHA256 = prereg.EXACT_MAIN_EXECUTABLE_SHA256
NORMALIZED_COMBINED_SHA256 = \
    "ddfef166af6c1dafd989019f87694526693623fa4ea2aa9e4d74f97c012fa093"
AUTHORITY_RECORD_FILE_SHA256 = \
    "0dba0d41057536397e4b37d0dad5387e66e33be76ce939dc39fcf3c373226fcc"
AUTHORITY_LEDGER_SHA256 = \
    "20ca2823aed2330474cba186b0029e1cfd39dfced3b6c76fa966dcac1ea3fae2"
SEALED_EXECUTABLE_PROTOCOL = "linux-sealed-memfd-executable/v1"
PUBLICATION_RULE = prereg.OUTPUT_LANE_PUBLICATION_RULE
MAX_EXECUTABLE_BYTES = 64 * 1024 * 1024
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_AUTHORITY_JSON_BYTES = 128 * 1024 * 1024
MAX_RUNTIME_CLOSURE_FILES = 64
AUTHORITATIVE_LOCK_PATH = "/tmp/leopard-gf8-authoritative.lock"
CANDIDATE_AUTHORITY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-candidate-authority/v1"
AUTHORITY_BUNDLE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-authority-bundle/v2"
LANE_BINDING_SCHEMA = \
    prereg.OUTPUT_LANE_TOPOLOGY_SCHEMA
PREARMED_STATE_ENTRY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-prearmed-state-entry/v1"
ATTEMPT_MANIFEST_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-attempt-manifest/v2"
JOURNAL_INTENT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-launch-intent/v1"
JOURNAL_RESULT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-child-result/v1"
JOURNAL_COMPLETE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-campaign-complete/v1"
CAMPAIGN_TRANSCRIPT_ENTRY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-transcript-entry/v1"
CAMPAIGN_TRANSCRIPT_ALLOCATION_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-transcript-allocation/v1"
LAUNCH_CONTEXT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-launch-context/v1"
ATTEMPT_NAME = re.compile(r"attempt-([0-9]{3})")
STAGING_NAME = re.compile(r"\.attempt-([0-9]{3})\.staging-([0-9]+)-([0-9a-f]{32})")
JOURNAL_INTENT_NAME = re.compile(r"intent-([0-9]{4})\.json")
JOURNAL_RESULT_NAME = re.compile(r"result-([0-9]{4})\.json")
JOURNAL_STAGING_NAME = re.compile(
    r"\.(intent-[0-9]{4}\.json|result-[0-9]{4}\.json|complete\.json)"
    r"\.staging-([0-9]+)-([0-9a-f]{32})")
ATTEMPTS_DIRECTORY = "attempts"
LANE_BINDING_FILE = "lane.json"
ARMED_FILE = "armed.json"
AUTHORITY_BUNDLE_FILE = "authority-bundle.json"
ATTEMPT_MANIFEST_FILE = "sidecar-manifest.json"
JOURNAL_DIRECTORY = "journal"
LAUNCH_CONSUMED_FILE = "launch-consumed.json"
CAMPAIGN_CONSUMED_NAME = re.compile(
    r"attempt-([0-9]{3})-transcript\.jsonl")
CAMPAIGN_PROGRESS_DIRECTORY_NAME = re.compile(
    r"attempt-([0-9]{3})-journal-checkpoints")
JOURNAL_CHECKPOINT_NAME = re.compile(
    r"checkpoint-([0-9]{4})-"
    r"(armed|intent-[0-9]{4}|result-[0-9]{4}|complete)\.marker")
PREARMED_HISTORY_NAME = re.compile(
    r"prearmed-([0-9]{4})-history\.jsonl")
PREARMED_BOUNDARY_NAME = re.compile(
    r"prearmed-([0-9]{4})-"
    r"(preregistered|qualifying|qualified|bridging|bridged|arming|presampling)"
    r"-reached\.marker")
PREARMED_BOUNDARY_STATES = (
    "PREREGISTERED", "QUALIFYING", "QUALIFIED", "BRIDGING", "BRIDGED",
    "ARMING", "PRESAMPLING",
)
PREARMED_STATE_TRANSITIONS = types.MappingProxyType({
    "INIT": "PREREGISTERED",
    "PREREGISTERED": "QUALIFYING",
    "QUALIFYING": "QUALIFIED",
    "QUALIFIED": "BRIDGING",
    "BRIDGING": "BRIDGED",
    "BRIDGED": "ARMING",
    "ARMING": "PRESAMPLING",
})
CHILD_ENVIRONMENT = types.MappingProxyType({
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
})
MAX_CHILD_STDOUT_BYTES = 16 * 1024 * 1024
MAX_CHILD_STDERR_BYTES = 8 * 1024 * 1024
MAX_CAMPAIGN_TRANSCRIPT_BYTES = 4 * 1024 * 1024
MAX_PREARMED_HISTORY_BYTES = 1 * 1024 * 1024
_ARMED_CAMPAIGN_TOKEN = object()


def _frozen_child_environment(
    environment: Mapping[str, str] = CHILD_ENVIRONMENT,
) -> Mapping[str, str]:
    """Return the immutable launch environment captured at source load."""
    return environment

MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
MFD_EXEC = getattr(os, "MFD_EXEC", 0x0010)
F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
REQUIRED_SEALS = F_SEAL_SEAL | F_SEAL_SHRINK | F_SEAL_GROW | F_SEAL_WRITE

IN_MODIFY = 0x00000002
IN_ATTRIB = 0x00000004
IN_CLOSE_WRITE = 0x00000008
IN_MOVED_FROM = 0x00000040
IN_MOVED_TO = 0x00000080
IN_CREATE = 0x00000100
IN_DELETE = 0x00000200
IN_DELETE_SELF = 0x00000400
IN_MOVE_SELF = 0x00000800
IN_UNMOUNT = 0x00002000
IN_Q_OVERFLOW = 0x00004000
IN_IGNORED = 0x00008000
IN_MASK_ADD = 0x20000000
RUNTIME_MUTATION_MASK = (
    IN_MODIFY | IN_ATTRIB | IN_CLOSE_WRITE | IN_MOVED_FROM | IN_MOVED_TO |
    IN_CREATE | IN_DELETE | IN_DELETE_SELF | IN_MOVE_SELF | IN_UNMOUNT |
    IN_Q_OVERFLOW)
RLIMIT_NAMES = (
    "RLIMIT_AS", "RLIMIT_CORE", "RLIMIT_CPU", "RLIMIT_DATA", "RLIMIT_FSIZE",
    "RLIMIT_MEMLOCK", "RLIMIT_MSGQUEUE", "RLIMIT_NICE", "RLIMIT_NOFILE",
    "RLIMIT_NPROC", "RLIMIT_RSS", "RLIMIT_RTPRIO", "RLIMIT_RTTIME",
    "RLIMIT_SIGPENDING", "RLIMIT_STACK",
)


class ArmingError(ValueError):
    """Live arming or launch authority was not exact and fail-closed."""


def _fail(message: str) -> NoReturn:
    raise ArmingError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_fd(descriptor: int, size: int) -> str:
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        block = os.pread(descriptor, min(1 << 20, size - offset), offset)
        _require(bool(block), "descriptor hash read made no progress")
        digest.update(block)
        offset += len(block)
    _require(not os.pread(descriptor, 1, offset),
             "descriptor grew during hashing")
    return digest.hexdigest()


def _read_regular_nofollow(path: Path, maximum: int, label: str) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        metadata = os.fstat(descriptor)
        _require(stat.S_ISREG(metadata.st_mode) and
                 0 < metadata.st_size <= maximum,
                 f"{label} is not one bounded regular file")
        chunks: list[bytes] = []
        offset = 0
        while offset < metadata.st_size:
            block = os.pread(descriptor, min(1 << 20, metadata.st_size - offset),
                             offset)
            _require(bool(block), f"{label} read made no progress")
            chunks.append(block)
            offset += len(block)
        return b"".join(chunks)
    except ArmingError:
        raise
    except OSError as error:
        raise ArmingError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _bounded_text_file(path: Path, maximum: int, label: str) -> str:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        metadata = os.fstat(descriptor)
        _require(stat.S_ISREG(metadata.st_mode),
                 f"{label} is not one regular file")
        chunks: list[bytes] = []
        total = 0
        while True:
            block = os.read(descriptor, min(65536, maximum + 1 - total))
            if not block:
                break
            chunks.append(block)
            total += len(block)
            _require(total <= maximum, f"{label} exceeds its byte bound")
        text = b"".join(chunks).decode("utf-8")
    except UnicodeError as error:
        raise ArmingError(f"{label} is not UTF-8") from error
    except OSError as error:
        raise ArmingError(f"cannot read {label}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    value = text.strip()
    _require(value != "" and "\x00" not in value,
             f"{label} is empty or contains NUL")
    return value


def _capture_process_umask() -> int:
    status = _bounded_text_file(
        Path("/proc/self/status"), 1 << 20, "process status")
    umask_values = [
        line.split(":", 1)[1].strip()
        for line in status.splitlines()
        if line.startswith("Umask:")
    ]
    _require(len(umask_values) == 1 and
             re.fullmatch(r"[0-7]{4}", umask_values[0]) is not None,
             "process umask is unavailable or malformed")
    return int(umask_values[0], 8)


def _validate_owner_creation_umask(label: str) -> None:
    _require((_capture_process_umask() & 0o700) == 0,
             f"{label} parent umask masks required owner permissions")


def _capture_parent_launch_state() -> dict[str, Any]:
    """Capture the parent attributes relevant to enacting fixed child state."""
    parameter = os.sched_getparam(0)
    limits: list[dict[str, Any]] = []
    for name in RLIMIT_NAMES:
        soft, hard = resource.getrlimit(getattr(resource, name))
        limits.append({"name": name, "soft": soft, "hard": hard})
    return {
        "nice": os.getpriority(os.PRIO_PROCESS, 0),
        "scheduler_policy": os.sched_getscheduler(0),
        "scheduler_priority": parameter.sched_priority,
        "umask": _capture_process_umask(),
        "rlimits": limits,
    }


def _validate_launch_context(value: Any) -> dict[str, Any]:
    try:
        record = prereg.validate_launch_context_contract(value)
    except Exception as error:
        raise ArmingError("generation-3 launch context was rejected") from error
    _require(record["schema"] == LAUNCH_CONTEXT_SCHEMA and
             record["scheduler_policy"] == os.SCHED_OTHER and
             tuple(item["name"] for item in record["rlimits"]) == RLIMIT_NAMES,
             "generation-3 launch context is unsupported on this host")
    for item in record["rlimits"]:
        _require(type(getattr(resource, item["name"], None)) is int,
                 f"generation-3 resource limit {item['name']} is unsupported")
    return record


def _validate_launch_context_current(
    expected: Mapping[str, Any], label: str,
) -> None:
    retained = _validate_launch_context(expected)
    current = _capture_parent_launch_state()
    _require(current["nice"] == retained["nice"] and
             current["scheduler_policy"] == retained["scheduler_policy"] and
             current["scheduler_priority"] == retained["scheduler_priority"],
             f"{label} parent scheduling context changed")
    for actual, required in zip(current["rlimits"], retained["rlimits"]):
        required_hard = required["hard"]
        actual_hard = actual["hard"]
        _require(actual["name"] == required["name"] and
                 (actual_hard == resource.RLIM_INFINITY or
                  (required_hard != resource.RLIM_INFINITY and
                   actual_hard >= required_hard)),
                 f"{label} cannot enact {required['name']}")


def current_host_facts(
    *, topology_sha256: str, allowed_cpus: Sequence[int],
) -> dict[str, Any]:
    """Capture stable host authority plus the current boot/affinity instance."""
    machine_id = _bounded_text_file(Path("/etc/machine-id"), 4096, "machine ID")
    boot_id = _bounded_text_file(
        Path("/proc/sys/kernel/random/boot_id"), 4096, "boot ID")
    cpuinfo = _bounded_text_file(Path("/proc/cpuinfo"), 16 << 20, "CPU info")
    models = {
        line.split(":", 1)[1].strip()
        for line in cpuinfo.splitlines()
        if line.lower().startswith("model name") and ":" in line
    }
    _require(len(models) == 1, "host CPU model is absent or heterogeneous")
    online = pair_acquire.parse_cpu_list(
        _bounded_text_file(Path(pair_acquire.ONLINE_CPUS_PATH),
                           pair_acquire.MAX_HOST_FILE_BYTES,
                           "online CPU list").encode("ascii"),
        "online CPU list")
    allowed = sorted(os.sched_getaffinity(0))
    _require(list(allowed_cpus) == allowed,
             "live allowed CPU set differs from qualification acquisition")
    ticks = os.sysconf("SC_CLK_TCK")
    authority = execution.host_authority_record(
        machine_id_sha256=_sha256_bytes(machine_id.encode("utf-8")),
        hostname=platform.node(), architecture=platform.machine(),
        cpu_model=next(iter(models)))
    instance = execution.host_instance_record(
        authority=authority, boot_id=boot_id,
        kernel_release=platform.release(), online_cpus=online,
        allowed_cpus=allowed, clock_ticks_per_second=ticks,
        topology_sha256=topology_sha256)
    return {"authority": authority, "instance": instance}


def _capture_live_host_instance(
    registration: Mapping[str, Any],
) -> dict[str, Any]:
    """Capture the production host through the non-injectable system reader."""
    qualification = registration["qualification"]
    _require(qualification["policy_evaluation_order"] == ["pair-and-domain"] and
             len(qualification["policies"]) == 1,
             "live host capture requires the frozen Track-A policy")
    reader = pair_acquire.SystemHostReader()
    allowed = sorted(reader.allowed_cpus())
    topology = pair_acquire._acquire_topology(
        reader, allowed,
        qualification["policies"][0]["candidate_primary_cpus"])
    host = current_host_facts(
        topology_sha256=contract.canonical_sha256(topology),
        allowed_cpus=allowed)
    _require(contract.exact_json_equal(
        _prereg_host_projection(host["authority"]),
        registration["host_authority"]),
        "live host differs from the ratified host authority")
    return host["instance"]


def _validate_runtime_file(entry: Mapping[str, Any], label: str) -> None:
    raw_path = entry.get("path")
    _require(type(raw_path) is str and Path(raw_path).is_absolute(),
             f"{label} path is not absolute")
    path = Path(raw_path)
    try:
        resolved = path.resolve(strict=True)
        data = _read_regular_nofollow(
            resolved, MAX_EXECUTABLE_BYTES, label)
        metadata = os.stat(resolved, follow_symlinks=False)
    except (OSError, ArmingError) as error:
        raise ArmingError(f"{label} cannot be revalidated") from error
    _require(len(data) == entry.get("size") and
             _sha256_bytes(data) == entry.get("sha256"),
             f"{label} bytes changed")
    if "mode" in entry:
        _require(stat.S_IMODE(metadata.st_mode) == entry["mode"],
                 f"{label} mode changed")


def _validate_runtime_closures_current(
    candidate_runtime: Mapping[str, Any],
    exact_main_record: Mapping[str, Any],
) -> None:
    for entry, label in _runtime_file_entries(
            candidate_runtime, exact_main_record):
        _validate_runtime_file(entry, label)


@_with_exact_source_imports
def _validate_controller_bindings_current(
    registration: Mapping[str, Any],
    _graph_getter: Callable[..., Any] = _frozen_controller_execution_graph,
    _environment_getter: Callable[..., Any] = _frozen_child_environment,
) -> None:
    controller_loads, controller_modules = _frozen_controller_ledger()
    _require(globals().get("_frozen_controller_execution_graph") is
             _graph_getter,
             "frozen controller graph getter binding changed")
    (armer_module_references, controller_module_references,
     executable_references, class_member_references,
     function_states) = _graph_getter()
    _require(_EXACT_CONTROLLER_LOADS is controller_loads and
             _EXACT_CONTROLLER_MODULES is controller_modules,
             "frozen controller ledger binding changed")
    _require(_CONTROLLER_EXECUTABLE_REFERENCE_IDENTITIES is
                executable_references and
             _CONTROLLER_CLASS_MEMBER_REFERENCE_IDENTITIES is
                class_member_references,
             "frozen controller execution graph binding changed")
    _require(_ARMER_MODULE_REFERENCE_IDENTITIES is armer_module_references and
             _CONTROLLER_MODULE_REFERENCE_IDENTITIES is
                controller_module_references,
             "frozen controller module graph binding changed")
    _require(_CONTROLLER_FUNCTION_EXECUTION_STATES is function_states,
             "frozen controller function-state binding changed")
    _validate_function_execution_states(function_states, "controller")
    _require(globals().get("_frozen_child_environment") is
                _environment_getter and
             CHILD_ENVIRONMENT is _environment_getter() and
             type(CHILD_ENVIRONMENT) is types.MappingProxyType,
             "frozen child environment binding changed")
    observed = prereg.current_controller_bindings(REPO_ROOT)
    _require(contract.exact_json_equal(
        observed, registration["controller_bindings"]),
        "reviewed controller bytes changed")
    binding_by_path = {
        (REPO_ROOT / binding["path"]).resolve(strict=True): binding["sha256"]
        for binding in observed
    }
    _require(tuple(binding["path"] for binding in observed) ==
             tuple(sorted(_REQUIRED_CONTROLLER_RELATIVE_PATHS)) and
             set(binding_by_path) == _REQUIRED_CONTROLLER_PATHS,
             "reviewed controller closure differs from bootstrap closure")
    for binding in observed:
        path = (REPO_ROOT / binding["path"]).resolve(strict=True)
        _require(_LOADED_CONTROLLER_SOURCE_SHA256.get(path) ==
                 binding["sha256"],
                 f"loaded controller bytes differ: {binding['path']}")

    self_module = _LIVE_ARMER_MODULE
    self_path = Path(__file__).resolve(strict=True)
    bootstrap_path = (HERE / "k65_gen3_exact_source_bootstrap.py").resolve(
        strict=True)
    _require(self_module is not None and
             sys.modules.get(__name__) is self_module,
             "live armer module alias was replaced")
    for name, expected in armer_module_references:
        _require(globals().get(name) is expected,
                 f"live armer module reference changed: {name}")
    frozen_controller_references = {
        (id(owner), name): expected
        for owner, name, expected in controller_module_references
    }
    for owner, name, expected in controller_module_references:
        _require(vars(owner).get(name) is expected,
                 f"controller module reference changed: "
                 f"{owner.__name__}.{name}")
    for owner, name, expected in executable_references:
        _require(vars(owner).get(name, _MISSING_MODULE_ALIAS) is expected,
                 f"controller executable binding changed: "
                 f"{owner.__name__}.{name}")
    for owner, name, expected in class_member_references:
        _require(vars(owner).get(name, _MISSING_MODULE_ALIAS) is expected,
                 f"controller class binding changed: "
                 f"{owner.__module__}.{owner.__qualname__}.{name}")
    ledger_objects: dict[int, tuple[Any, Path, str]] = {}
    ledger_globals: dict[int, tuple[Any, Path, str]] = {}
    loaded_paths: set[Path] = set()
    modules_by_alias: dict[str, set[int]] = {}
    for alias, module, path, digest in controller_loads:
        entry = (module, path, digest)
        _require(path in _REQUIRED_CONTROLLER_PATHS and
                 vars(module).get("__exact_source_sha256__") == digest and
                 _controller_module_path(module) == path and
                 digest == binding_by_path.get(path),
                 f"exact controller alias changed: {alias}")
        previous = ledger_objects.get(id(module))
        _require(previous is None or previous == entry,
                 f"controller object identity collision: {alias}")
        ledger_objects[id(module)] = entry
        ledger_globals[id(module.__dict__)] = entry
        modules_by_alias.setdefault(alias, set()).add(id(module))
        loaded_paths.add(path)
    for alias, object_ids in modules_by_alias.items():
        current = sys.modules.get(alias)
        _require(current is None or id(current) in object_ids,
                 f"exact controller alias was replaced: {alias}")

    special_paths = {self_path, bootstrap_path}
    _require(loaded_paths | special_paths == _REQUIRED_CONTROLLER_PATHS,
             "not every reviewed controller was exact-source loaded")
    for alias, module in tuple(sys.modules.items()):
        if module is None:
            continue
        path = _controller_module_path(module)
        if path is None:
            continue
        if module is self_module:
            _require(alias == __name__ and path == self_path and
                     vars(module).get("__exact_source_sha256__") ==
                        binding_by_path[self_path],
                     "live armer bootstrap identity changed")
            continue
        _require(id(module) in modules_by_alias.get(alias, set()),
                 f"controller alias escaped exact-source ledger: {alias}")

    def validate_function(value: Any, label: str) -> None:
        if not isinstance(value, types.FunctionType):
            return
        try:
            source_path = Path(value.__code__.co_filename).resolve(strict=True)
        except (OSError, RuntimeError, ValueError):
            return
        if source_path not in _REQUIRED_CONTROLLER_PATHS:
            return
        if value.__globals__ is globals():
            _require(source_path == self_path,
                     f"controller function source changed: {label}")
            return
        entry = ledger_globals.get(id(value.__globals__))
        _require(entry is not None and entry[1] == source_path,
                 f"controller function escaped exact-source ledger: {label}")

    roots = [
        prereg, execution, plan_contract, build_provenance, git_capture,
        balanced_evidence_common, result_contract, pair_acquire, pair_bridge,
        pair_verify, identity_contract, baseline_acquire, authority_record,
        authority_verifier,
    ]
    roots.extend(entry[1] for entry in controller_loads)
    seen: set[int] = set()
    for module in roots:
        if id(module) in seen:
            continue
        seen.add(id(module))
        for name, value in vars(module).items():
            label = f"{module.__name__}.{name}"
            if isinstance(value, types.ModuleType):
                expected = frozen_controller_references.get(
                    (id(module), name))
                _require(expected is value,
                         f"controller module reference escaped frozen graph: "
                         f"{label}")
                path = _controller_module_path(expected)
                if id(expected) in ledger_objects:
                    _require(path == ledger_objects[id(expected)][1],
                             f"controller module reference escaped ledger: "
                             f"{label}")
                else:
                    _require(path is None,
                             f"unledgered controller module reference: "
                             f"{label}")
            validate_function(value, label)
            if isinstance(value, type):
                for member_name, member in vars(value).items():
                    if isinstance(member, (staticmethod, classmethod)):
                        member = member.__func__
                    validate_function(member, f"{label}.{member_name}")

    _require(getattr(git_capture, "_build_provenance", None) is
             build_provenance,
             "Git capture does not use the bound build-provenance module")


def _prereg_host_projection(host_authority: Mapping[str, Any]) -> dict[str, Any]:
    validated = execution.validate_host_authority(host_authority)
    return {
        name: validated[name]
        for name in ("machine_id_sha256", "hostname", "architecture", "cpu_model")
    }


def bind_exact_main_path_variant(authority_lane: Path) -> dict[str, Any]:
    """Replay the sealed authority and bind its non-pooled normalized variant."""
    _require(isinstance(authority_lane, Path) and not authority_lane.is_symlink(),
             "exact-main authority lane is not a direct path")
    try:
        root = authority_lane.resolve(strict=True)
    except OSError as error:
        raise ArmingError("exact-main authority lane is absent") from error
    _require(root == authority_lane and root.is_dir(),
             "exact-main authority lane is not canonical")
    try:
        with authority_verifier.read_sealed_tree(root) as sealed_tree:
            verdict = authority_verifier.verify_authority_lane(sealed_tree)
            record_data = sealed_tree.read_file(
                "baseline-authority.json", maximum_bytes=MAX_JSON_BYTES)
            data = sealed_tree.read_file(
                PATH_VARIANT_RELATIVE, maximum_bytes=MAX_EXECUTABLE_BYTES)
            sealed_tree.reverify()
    except Exception as error:
        raise ArmingError("sealed exact-main authority did not replay") from error
    _require(
        type(verdict) is dict and
        verdict.get("schema") == authority_verifier.VERIFIER_SCHEMA and
        verdict.get("outcome") == "promoted_authority" and
        verdict.get("promoted") is True and
        verdict.get("record", {}).get("record_sha256") ==
            prereg.EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
        verdict.get("seal", {}).get("sha256sums_sha256") ==
            AUTHORITY_LEDGER_SHA256,
        "sealed exact-main authority verdict differs")
    _require(_sha256_bytes(record_data) == AUTHORITY_RECORD_FILE_SHA256,
             "exact-main authority record file hash differs")
    try:
        record = authority_record.load_baseline_authority_record(record_data)
    except Exception as error:
        raise ArmingError("exact-main authority record is invalid") from error
    build = record.get("builds", {}).get("path_variant", {})
    normalized = record.get("identity", {}).get("path_variant")
    relative = build.get("executable", {}).get("retained_relative_path")
    _require(
        record.get("record_sha256") == prereg.EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
        relative == PATH_VARIANT_RELATIVE and
        build.get("executable", {}).get("sha256") == PATH_VARIANT_RAW_SHA256 and
        build.get("executable", {}).get("size") == PATH_VARIANT_SIZE and
        record.get("identity", {}).get("normalized_match") is True and
        record.get("identity", {}).get("path_variant_raw_bytes_differ") is True,
        "path-variant authority semantics differ")
    _require(len(data) == PATH_VARIANT_SIZE and
             _sha256_bytes(data) == PATH_VARIANT_RAW_SHA256 and
             PATH_VARIANT_RAW_SHA256 != CANONICAL_RAW_SHA256,
             "path-variant raw executable differs")
    try:
        replayed = identity_contract.verify_normalized_code_identity_against_elf_bytes(
            data, normalized, roots=build["roots"])
    except Exception as error:
        raise ArmingError("path-variant normalized identity did not replay") from error
    _require(replayed.get("combined_sha256") == NORMALIZED_COMBINED_SHA256,
             "path-variant normalized combined identity differs")
    return {
        "root": root,
        "artifact_relative_path": PATH_VARIANT_RELATIVE,
        "artifact_data": data,
        "verdict": verdict,
        "record": record,
        "record_data": record_data,
    }


_CANDIDATE_AUTHORITY_KEYS = frozenset((
    "schema", "generation", "status", "source", "candidate",
    "build_closure", "controller_bindings_sha256", "inventory",
))
_CANDIDATE_SOURCE_KEYS = frozenset((
    "commit", "tree", "git_capture_sha256", "source_archive_sha256",
))
_CANDIDATE_EXECUTABLE_KEYS = frozenset((
    "relative_path", "sha256", "size",
))
_CANDIDATE_BUILD_KEYS = frozenset((
    "build_provenance_sha256", "reproducible_build_core_sha256",
    "reproducible_build_proof_sha256", "attestation_header_sha256",
    "runtime_closure_sha256",
))
_INVENTORY_KEYS = frozenset(("relative_path", "sha256", "size"))
_CANDIDATE_PAYLOAD_PATHS = (
    "artifacts/bench_leopard2",
    "build/benchmark-source-attestation.h",
    "build/build-provenance.json",
    "build/reproducible-build-core.json",
    "build/reproducible-build-proof.json",
    "runtime/runtime-closure.json",
    "source/candidate-source.tar",
    "source/git-capture.json",
)


def _exact_mapping(value: Any, keys: frozenset[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys,
             f"{label} is not an exact object")
    return value


def _canonical_json_document(data: bytes, label: str) -> dict[str, Any]:
    value = contract.strict_json_loads(data, label)
    _require(type(value) is dict and data == contract.canonical_json_bytes(value),
             f"{label} bytes are not one canonical object")
    return value


def _hex_digest(value: Any, label: str) -> str:
    _require(type(value) is str and
             re.fullmatch(r"[0-9a-f]{64}", value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _candidate_runtime_closure(value: Any) -> dict[str, Any]:
    record = _exact_mapping(
        value, frozenset(("schema", "dependencies", "launchers")),
        "candidate runtime closure")
    _require(record["schema"] ==
             "leopard2-k65-gen3-runtime-closure/v1" and
             type(record["dependencies"]) is list and
             type(record["launchers"]) is list and
             len(record["launchers"]) == 2 and
             2 <= len(record["dependencies"]) + len(record["launchers"]) <=
                MAX_RUNTIME_CLOSURE_FILES,
             "candidate runtime closure metadata differs")
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(record["dependencies"] + record["launchers"]):
        item = _exact_mapping(
            raw, frozenset(("path", "sha256", "size", "mode", "role")),
            f"candidate runtime entry {index}")
        path = Path(str(item["path"]))
        _require(path.is_absolute() and "\0" not in str(path) and
                 type(item["size"]) is int and item["size"] > 0 and
                 type(item["mode"]) is int and 0 <= item["mode"] <= 0o7777 and
                 type(item["role"]) is str and item["role"],
                 f"candidate runtime entry {index} is malformed")
        normalized.append({
            "path": str(path),
            "sha256": _hex_digest(
                item["sha256"], f"candidate runtime entry {index} hash"),
            "size": item["size"], "mode": item["mode"],
            "role": item["role"],
        })
    split = len(record["dependencies"])
    _require(normalized == sorted(
        normalized, key=lambda item: (item["role"], item["path"])) and
        len({item["path"] for item in normalized}) == len(normalized) and
        len({(item["role"], item["path"]) for item in normalized}) ==
            len(normalized) and
        {item["path"] for item in normalized[split:]} ==
            {"/usr/bin/prlimit", "/usr/bin/taskset"},
        "candidate runtime closure is not canonical, bounded, and unique")
    return {
        "schema": record["schema"],
        "dependencies": normalized[:split],
        "launchers": normalized[split:],
    }


def _runtime_file_entries(
    candidate_runtime: Mapping[str, Any], exact_main_record: Mapping[str, Any],
) -> list[tuple[Mapping[str, Any], str]]:
    candidate = _candidate_runtime_closure(candidate_runtime)
    entries: list[tuple[Mapping[str, Any], str]] = [
        (entry, f"candidate runtime entry {index}")
        for index, entry in enumerate(
            candidate["dependencies"] + candidate["launchers"])
    ]
    runtime = exact_main_record.get("runtime_closure", {})
    records = runtime.get("records")
    _require(type(records) is list and len(records) == 3,
             "exact-main runtime closure differs")
    path_variant = next(
        (record for record in records
         if type(record) is dict and record.get("role") == "path_variant"),
        None)
    _require(type(path_variant) is dict and
             path_variant.get("executable_sha256") == PATH_VARIANT_RAW_SHA256 and
             type(path_variant.get("dependencies")) is list and
             0 < len(path_variant["dependencies"]) <=
                MAX_RUNTIME_CLOSURE_FILES,
             "exact-main path-variant runtime closure is absent")
    observed_exact_paths: set[str] = set()
    for index, entry in enumerate(path_variant["dependencies"]):
        _require(type(entry) is dict,
                 "exact-main runtime entry is not an object")
        if entry.get("kind") == "virtual":
            _require(entry.get("path") is None and entry.get("sha256") is None and
                     entry.get("size") is None,
                     "exact-main virtual runtime entry differs")
            continue
        _require(entry.get("kind") == "file",
                 "exact-main runtime entry kind differs")
        _require(type(entry.get("path")) is str and
                 entry["path"] not in observed_exact_paths,
                 "exact-main runtime entries are not unique")
        observed_exact_paths.add(entry["path"])
        entries.append((entry, f"exact-main runtime entry {index}"))
    return entries


class _RetainedRuntimeGuard:
    """Retain exact runtime inodes and detect path mutation across exec.

    Hashes are established once while watches are live.  Thereafter retained
    descriptors, complete path-component watches, inode/metadata replay, and
    fail-closed queue handling bridge the final check-to-exec window without
    re-reading megabytes of shared libraries into benchmark caches.
    """

    def __init__(self) -> None:
        self._inotify_fd = -1
        self._files: list[dict[str, Any]] = []
        self._directory_targets: dict[int, set[str]] = {}
        self._file_watches: set[int] = set()
        self._closed = False

    @classmethod
    def acquire(
        cls, candidate_runtime: Mapping[str, Any],
        exact_main_record: Mapping[str, Any],
    ) -> "_RetainedRuntimeGuard":
        guard = cls()
        try:
            guard._open_inotify()
            entries = _runtime_file_entries(
                candidate_runtime, exact_main_record)
            for entry, unused_label in entries:
                guard._watch_components(Path(str(entry["path"])))
            for entry, label in entries:
                guard._retain_file(entry, label)
            guard.validate_current("initial runtime closure")
            return guard
        except BaseException:
            guard.close()
            raise

    def _open_inotify(self) -> None:
        libc = ctypes.CDLL(None, use_errno=True)
        init = libc.inotify_init1
        init.argtypes = [ctypes.c_int]
        init.restype = ctypes.c_int
        descriptor = init(os.O_NONBLOCK | os.O_CLOEXEC)
        if descriptor < 0:
            error = ctypes.get_errno()
            raise ArmingError("runtime mutation guard could not be created") \
                from OSError(error, os.strerror(error))
        self._inotify_fd = descriptor

    def _add_watch(self, path: Path, mask: int) -> int:
        libc = ctypes.CDLL(None, use_errno=True)
        add = libc.inotify_add_watch
        add.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
        add.restype = ctypes.c_int
        watch = add(self._inotify_fd, os.fsencode(path), mask | IN_MASK_ADD)
        if watch < 0:
            error = ctypes.get_errno()
            raise ArmingError(f"runtime mutation watch failed for {path}") \
                from OSError(error, os.strerror(error))
        return watch

    def _watch_components(self, path: Path) -> None:
        _require(path.is_absolute(), "runtime watch path is not absolute")
        parent = Path("/")
        for component in path.parts[1:]:
            watch = self._add_watch(parent, RUNTIME_MUTATION_MASK)
            self._directory_targets.setdefault(watch, set()).add(component)
            parent = parent / component

    def _retain_file(self, entry: Mapping[str, Any], label: str) -> None:
        raw_path = Path(str(entry["path"]))
        try:
            resolved = raw_path.resolve(strict=True)
            self._watch_components(resolved)
            descriptor = os.open(
                resolved, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        except OSError as error:
            raise ArmingError(f"{label} could not be retained") from error
        try:
            metadata = os.fstat(descriptor)
            path_metadata = os.stat(resolved, follow_symlinks=False)
            _require(stat.S_ISREG(metadata.st_mode) and
                     0 < metadata.st_size <= MAX_EXECUTABLE_BYTES and
                     metadata.st_size == entry.get("size") and
                     _sha256_fd(descriptor, metadata.st_size) ==
                        entry.get("sha256") and
                     (metadata.st_dev, metadata.st_ino) ==
                        (path_metadata.st_dev, path_metadata.st_ino),
                     f"{label} retained bytes or inode differ")
            if "mode" in entry:
                _require(stat.S_IMODE(metadata.st_mode) == entry["mode"],
                         f"{label} retained mode differs")
            file_watch = self._add_watch(resolved, RUNTIME_MUTATION_MASK)
            self._file_watches.add(file_watch)
            snapshot = (
                metadata.st_dev, metadata.st_ino, metadata.st_mode,
                metadata.st_uid, metadata.st_gid, metadata.st_nlink,
                metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
            )
            self._files.append({
                "descriptor": descriptor, "raw_path": raw_path,
                "resolved_path": resolved, "snapshot": snapshot,
                "label": label,
            })
            descriptor = -1
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    def _reject_mutation_events(self, label: str) -> None:
        while True:
            try:
                data = os.read(self._inotify_fd, 64 * 1024)
            except BlockingIOError:
                return
            except OSError as error:
                raise ArmingError(
                    f"{label}: runtime mutation guard could not be read") from error
            _require(bool(data), f"{label}: runtime mutation guard closed")
            offset = 0
            while offset < len(data):
                _require(offset + 16 <= len(data),
                         f"{label}: truncated runtime mutation event")
                watch, mask, unused_cookie, name_size = struct.unpack_from(
                    "iIII", data, offset)
                end = offset + 16 + name_size
                _require(end <= len(data),
                         f"{label}: truncated runtime mutation name")
                name = data[offset + 16:end].split(b"\0", 1)[0]
                offset = end
                if mask & (IN_Q_OVERFLOW | IN_IGNORED):
                    _fail(f"{label}: runtime mutation watch lost authority")
                if watch in self._file_watches:
                    _fail(f"{label}: retained runtime file was mutated")
                targets = self._directory_targets.get(watch)
                if targets is None:
                    _fail(f"{label}: runtime mutation watch is unknown")
                decoded = os.fsdecode(name)
                if not decoded or decoded in targets:
                    _fail(f"{label}: retained runtime path was mutated")

    def validate_current(self, label: str) -> None:
        _require(not self._closed and self._inotify_fd >= 0,
                 f"{label}: runtime mutation guard is closed")
        for retained in self._files:
            descriptor = retained["descriptor"]
            metadata = os.fstat(descriptor)
            observed = (
                metadata.st_dev, metadata.st_ino, metadata.st_mode,
                metadata.st_uid, metadata.st_gid, metadata.st_nlink,
                metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
            )
            try:
                resolved = retained["raw_path"].resolve(strict=True)
                path_metadata = os.stat(
                    retained["resolved_path"], follow_symlinks=False)
            except OSError as error:
                raise ArmingError(
                    f"{label}: retained runtime path is unavailable") from error
            _require(observed == retained["snapshot"] and
                     resolved == retained["resolved_path"] and
                     (metadata.st_dev, metadata.st_ino) ==
                        (path_metadata.st_dev, path_metadata.st_ino),
                     f"{label}: retained runtime identity changed")
        self._reject_mutation_events(label)

    def close(self) -> None:
        if self._closed and not self._files and self._inotify_fd < 0:
            return
        self._closed = True
        while self._files:
            retained = self._files[-1]
            descriptor = retained["descriptor"]
            try:
                os.close(descriptor)
            except OSError:
                pass
            retained["descriptor"] = -1
            self._files.pop()
        if self._inotify_fd >= 0:
            descriptor = self._inotify_fd
            try:
                os.close(descriptor)
            except OSError:
                pass
            self._inotify_fd = -1


def _acquire_runtime_guard(
    candidate_runtime: Mapping[str, Any], exact_main_record: Mapping[str, Any],
) -> _RetainedRuntimeGuard:
    return _RetainedRuntimeGuard.acquire(candidate_runtime, exact_main_record)


class _RetainedArchiveBytes:
    """Minimal read-only tree surface for the audited archive verifier."""

    def __init__(self, data: bytes) -> None:
        self._data = data

    def size(self, name: str) -> int:
        _require(name == "candidate-source.tar",
                 "candidate source archive name differs")
        return len(self._data)

    def pread(self, name: str, offset: int, size: int) -> bytes:
        _require(name == "candidate-source.tar" and
                 type(offset) is int and offset >= 0 and
                 type(size) is int and size >= 0 and
                 offset + size <= len(self._data),
                 "candidate source archive read is out of bounds")
        return self._data[offset:offset + size]


def _validate_candidate_source_archive(
    data: bytes, validated_git: Mapping[str, Any], commit: str,
) -> None:
    """Prove canonical Git-archive members equal the retained Git blobs."""
    tracked_list = validated_git.get("tracked_files")
    _require(type(tracked_list) is list and tracked_list,
             "candidate Git capture has no tracked-file closure")
    tracked: dict[str, Mapping[str, Any]] = {}
    for index, item in enumerate(tracked_list):
        _require(type(item) is dict and type(item.get("path")) is str and
                 item["path"] not in tracked,
                 f"candidate tracked file {index} is malformed")
        tracked[item["path"]] = item
    try:
        archive = _RetainedArchiveBytes(data)
        members = authority_verifier._verify_source_archive(
            archive, "candidate-source.tar",
            expected_prefix="candidate-source/", expected_commit=commit)
        authority_verifier._verify_archive_git_binding(
            members, tracked, "candidate source")
    except ArmingError:
        raise
    except Exception as error:
        raise ArmingError(
            "candidate source archive did not match its Git capture") from error


def bind_candidate_authority_lane(
    authority_lane: Path, preregistration: Mapping[str, Any],
) -> dict[str, Any]:
    """Replay a preregistration-pinned candidate source/build authority."""
    registration = prereg.validate_preregistration(
        preregistration, verify_files=True)
    _require(isinstance(authority_lane, Path) and
             authority_lane.is_absolute() and not authority_lane.is_symlink(),
             "candidate authority lane is not a direct absolute path")
    try:
        root = authority_lane.resolve(strict=True)
    except OSError as error:
        raise ArmingError("candidate authority lane is absent") from error
    _require(root == authority_lane and root.is_dir(),
             "candidate authority lane is not canonical")
    candidate_binding = registration["candidate_executable"]
    try:
        with authority_verifier.read_sealed_tree(root) as sealed_tree:
            authority_verifier.verify_tree_metadata(sealed_tree)
            ledger = authority_verifier.verify_sha256sums(sealed_tree)
            record_data = sealed_tree.read_file(
                "candidate-authority.json", maximum_bytes=MAX_JSON_BYTES)
            record = _canonical_json_document(
                record_data, "candidate authority record")
            _exact_mapping(
                record, _CANDIDATE_AUTHORITY_KEYS,
                "candidate authority record")
            source = _exact_mapping(
                record["source"], _CANDIDATE_SOURCE_KEYS,
                "candidate authority source")
            candidate = _exact_mapping(
                record["candidate"], _CANDIDATE_EXECUTABLE_KEYS,
                "candidate authority executable")
            build = _exact_mapping(
                record["build_closure"], _CANDIDATE_BUILD_KEYS,
                "candidate authority build closure")
            inventory = record["inventory"]
            _require(type(inventory) is list and
                     len(inventory) == len(_CANDIDATE_PAYLOAD_PATHS),
                     "candidate authority inventory length differs")
            normalized_inventory: list[dict[str, Any]] = []
            for index, raw in enumerate(inventory):
                item = _exact_mapping(
                    raw, _INVENTORY_KEYS,
                    f"candidate authority inventory {index}")
                path = item["relative_path"]
                _require(type(path) is str and
                         path in _CANDIDATE_PAYLOAD_PATHS and
                         type(item["size"]) is int and item["size"] > 0,
                         f"candidate authority inventory {index} differs")
                normalized_inventory.append({
                    "relative_path": path,
                    "sha256": _hex_digest(
                        item["sha256"],
                        f"candidate authority inventory {index} hash"),
                    "size": item["size"],
                })
            normalized_inventory.sort(key=lambda item: item["relative_path"])
            _require(
                [item["relative_path"] for item in normalized_inventory] ==
                    list(_CANDIDATE_PAYLOAD_PATHS) and
                inventory == normalized_inventory,
                "candidate authority inventory is not exact and canonical")
            expected_files = set(_CANDIDATE_PAYLOAD_PATHS) | {
                "candidate-authority.json", "SHA256SUMS",
                "TREE-METADATA.json",
            }
            _require(set(sealed_tree.files) == expected_files,
                     "candidate authority file set differs")
            payload: dict[str, bytes] = {}
            for item in normalized_inventory:
                data = sealed_tree.read_file(
                    item["relative_path"],
                    maximum_bytes=(MAX_EXECUTABLE_BYTES
                                   if item["relative_path"] ==
                                      "artifacts/bench_leopard2"
                                   else MAX_AUTHORITY_JSON_BYTES))
                _require(len(data) == item["size"] and
                         _sha256_bytes(data) == item["sha256"],
                         f"candidate authority payload differs: "
                         f"{item['relative_path']}")
                payload[item["relative_path"]] = data
            sealed_tree.reverify()
    except ArmingError:
        raise
    except Exception as error:
        raise ArmingError("sealed candidate authority did not replay") from error
    source_expected = registration["candidate_source"]
    _require(
        record["schema"] == CANDIDATE_AUTHORITY_SCHEMA and
        record["generation"] == GENERATION and
        record["status"] == "promoted_authority" and
        _sha256_bytes(record_data) ==
            candidate_binding["authority_record_sha256"] and
        ledger["sha256"] == candidate_binding["authority_ledger_sha256"] and
        source["commit"] == source_expected["commit"] and
        source["tree"] == source_expected["tree"] and
        candidate["relative_path"] == "artifacts/bench_leopard2" and
        candidate["sha256"] == candidate_binding["sha256"] and
        candidate["size"] == candidate_binding["size"] and
        build["build_provenance_sha256"] ==
            candidate_binding["build_provenance_sha256"] and
        build["reproducible_build_core_sha256"] ==
            candidate_binding["reproducible_build_core_sha256"] and
        record["controller_bindings_sha256"] ==
            contract.canonical_sha256(registration["controller_bindings"]),
        "candidate authority does not match preregistration")
    inventory_by_path = {
        item["relative_path"]: item for item in normalized_inventory}
    joins = {
        "source/git-capture.json": source["git_capture_sha256"],
        "source/candidate-source.tar": source["source_archive_sha256"],
        "build/build-provenance.json": build["build_provenance_sha256"],
        "build/reproducible-build-core.json":
            build["reproducible_build_core_sha256"],
        "build/reproducible-build-proof.json":
            build["reproducible_build_proof_sha256"],
        "build/benchmark-source-attestation.h":
            build["attestation_header_sha256"],
        "runtime/runtime-closure.json": build["runtime_closure_sha256"],
        "artifacts/bench_leopard2": candidate["sha256"],
    }
    _require(all(inventory_by_path[path]["sha256"] == digest
                 for path, digest in joins.items()),
             "candidate authority record-to-inventory joins differ")
    git_record = _canonical_json_document(
        payload["source/git-capture.json"], "candidate Git capture")
    try:
        validated_git = git_capture.validate_git_capture(
            git_record, str(REPO_ROOT), source["commit"],
            require_detached=False)
    except Exception as error:
        raise ArmingError("candidate Git capture did not replay") from error
    _require(validated_git["tree"] == source["tree"],
             "candidate Git capture tree differs")
    _validate_candidate_source_archive(
        payload["source/candidate-source.tar"], validated_git,
        source["commit"])
    provenance = _canonical_json_document(
        payload["build/build-provenance.json"],
        "candidate build provenance")
    proof = _canonical_json_document(
        payload["build/reproducible-build-proof.json"],
        "candidate reproducible-build proof")
    core = _canonical_json_document(
        payload["build/reproducible-build-core.json"],
        "candidate reproducible-build core")
    try:
        build_provenance.validate_reproducible_build_proof(
            proof, provenance, label="K65 candidate")
        expected_core = build_provenance.compare_reproducible_builds(
            provenance, provenance)
        expected_header = \
            build_provenance._canonical_replay_attestation_header_bytes(
                provenance)
    except Exception as error:
        raise ArmingError("candidate reproducible-build proof did not replay") \
            from error
    _require(contract.exact_json_equal(core, expected_core) and
             provenance.get("tracked_source_manifest", {}).get("git", {}).get(
                 "commit") == source["commit"] and
             provenance.get("tracked_source_manifest", {}).get("git", {}).get(
                 "tree") == source["tree"] and
             provenance.get("tracked_source_manifest", {}).get("git", {}).get(
                 "dirty") is False and
             provenance.get("executable", {}).get("sha256") ==
                 candidate["sha256"] and
             provenance.get("executable", {}).get("size") ==
                 candidate["size"] and
             payload["build/benchmark-source-attestation.h"] == expected_header,
             "candidate source/build/executable closure differs")
    runtime = _candidate_runtime_closure(_canonical_json_document(
        payload["runtime/runtime-closure.json"],
        "candidate runtime closure"))
    artifact_data = payload["artifacts/bench_leopard2"]
    _require(artifact_data.startswith(b"\x7fELF"),
             "candidate authority artifact is not ELF")
    return {
        "root": root,
        "record": record,
        "record_data": record_data,
        "ledger": ledger,
        "artifact_data": artifact_data,
        "runtime_closure": runtime,
        "source_authority": execution.candidate_source_authority_record(
            registration,
            build_provenance_sha256=build["build_provenance_sha256"],
            reproducible_build_core_sha256=
                build["reproducible_build_core_sha256"],
            authority_record_sha256=_sha256_bytes(record_data),
            authority_ledger_sha256=ledger["sha256"],
            host_authority={
                "schema": execution.HOST_AUTHORITY_SCHEMA,
                **registration["host_authority"],
            }),
    }


def capture_sealed_executable_bytes(
    data: bytes, *, expected_sha256: str, expected_size: int, label: str,
    authority_relative_path: str,
) -> tuple[int, dict[str, Any]]:
    """Install already fd-anchored authority bytes into an immutable memfd."""
    _require(type(data) is bytes and len(data) == expected_size and
             0 < len(data) <= MAX_EXECUTABLE_BYTES and
             _sha256_bytes(data) == expected_sha256 and
             data.startswith(b"\x7fELF"),
             f"{label} authority bytes differ")
    snapshot = -1
    try:
        snapshot = _create_memfd(
            "leopard2-k65-gen3-" + label.replace(" ", "-"))
        _write_all(snapshot, data)
        os.fchmod(snapshot, 0o500)
        fcntl.fcntl(snapshot, F_ADD_SEALS, REQUIRED_SEALS)
        retained = sealed_descriptor_identity(
            snapshot, label, capture_nonce=secrets.token_hex(32))
        _require(retained["raw_sha256"] == expected_sha256 and
                 retained["size"] == expected_size,
                 f"{label} sealed snapshot differs")
        descriptor = snapshot
        snapshot = -1
        return descriptor, {
            "source": {
                "authority_relative_path": authority_relative_path,
                "size": expected_size,
                "raw_sha256": expected_sha256,
            },
            "snapshot": retained,
        }
    finally:
        if snapshot >= 0:
            os.close(snapshot)


def _create_memfd(name: str) -> int:
    creator = getattr(os, "memfd_create", None)
    _require(sys.platform.startswith("linux") and creator is not None,
             "sealed executable snapshots require Linux memfd_create")
    try:
        return creator(name, MFD_CLOEXEC | MFD_ALLOW_SEALING | MFD_EXEC)
    except OSError as error:
        if error.errno != errno.EINVAL:
            raise
    return creator(name, MFD_CLOEXEC | MFD_ALLOW_SEALING)


def _process_start_ticks() -> int:
    value = _bounded_text_file(Path("/proc/self/stat"), 1 << 20,
                               "current process stat")
    close = value.rfind(")")
    _require(close > 0, "current process stat is malformed")
    fields = value[close + 2:].split()
    _require(len(fields) >= 20 and fields[19].isdigit(),
             "current process start identity is malformed")
    return int(fields[19])


def sealed_descriptor_identity(
    descriptor: int, label: str, *, capture_nonce: str,
) -> dict[str, Any]:
    _require(type(descriptor) is int and descriptor >= 0,
             f"{label} descriptor is invalid")
    try:
        metadata = os.fstat(descriptor)
        seals = fcntl.fcntl(descriptor, F_GET_SEALS)
    except OSError as error:
        raise ArmingError(f"cannot inspect {label} descriptor") from error
    _require(stat.S_ISREG(metadata.st_mode) and
             0 < metadata.st_size <= MAX_EXECUTABLE_BYTES and
             metadata.st_mode & 0o111 and
             seals & REQUIRED_SEALS == REQUIRED_SEALS,
             f"{label} snapshot is not a sealed executable")
    _require(type(capture_nonce) is str and
             re.fullmatch(r"[0-9a-f]{64}", capture_nonce) is not None,
             f"{label} capture nonce is invalid")
    prefix = os.pread(descriptor, 4, 0)
    _require(prefix == b"\x7fELF", f"{label} snapshot is not ELF")
    return {
        "protocol": SEALED_EXECUTABLE_PROTOCOL,
        "boot_id": _bounded_text_file(
            Path("/proc/sys/kernel/random/boot_id"), 4096, "boot ID"),
        "session_nonce": capture_nonce,
        "creator_pid": os.getpid(),
        "creator_start_ticks": _process_start_ticks(),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "gid": metadata.st_gid,
        "size": metadata.st_size,
        "mode": stat.S_IMODE(metadata.st_mode),
        "raw_sha256": _sha256_fd(descriptor, metadata.st_size),
        "seals": seals,
        "elf": True,
    }


def revalidate_sealed_descriptor(
    descriptor: int, expected: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    _require(type(expected) is dict and
             type(expected.get("session_nonce")) is str,
             f"{label} retained descriptor identity is malformed")
    observed = sealed_descriptor_identity(
        descriptor, label, capture_nonce=expected["session_nonce"])
    _require(contract.exact_json_equal(observed, expected),
             f"{label} retained descriptor identity changed")
    return observed


def _canonical_lane(path: Path) -> Path:
    _require(isinstance(path, Path) and path.is_absolute() and
             not path.is_symlink(), "arming lane is not a direct absolute path")
    try:
        lane = path.resolve(strict=True)
        metadata = os.stat(lane, follow_symlinks=False)
    except OSError as error:
        raise ArmingError("arming lane is absent") from error
    _require(lane == path and stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             stat.S_IMODE(metadata.st_mode) & 0o077 == 0,
             "arming lane is not one canonical owner-private directory")
    return lane


def _validate_lane_authority_separation(
    lane: Path, candidate_authority: Mapping[str, Any],
    exact_main: Mapping[str, Any],
) -> Path:
    lane_root = _canonical_lane(lane)
    candidate_root = candidate_authority.get("root")
    exact_main_root = exact_main.get("root")
    _require(
        isinstance(candidate_root, Path) and candidate_root.is_absolute() and
        candidate_root == candidate_root.resolve(strict=False) and
        isinstance(exact_main_root, Path) and exact_main_root.is_absolute() and
        exact_main_root == exact_main_root.resolve(strict=False),
        "bound authority roots are not canonical absolute paths")

    def overlaps(first: Path, second: Path) -> bool:
        return (first == second or first in second.parents or
                second in first.parents)

    _require(not any(overlaps(first, second) for first, second in (
        (lane_root, candidate_root),
        (lane_root, exact_main_root),
        (candidate_root, exact_main_root),
    )), "output and authority lanes overlap")
    return lane_root


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        count = os.write(descriptor, content[offset:])
        _require(type(count) is int and count > 0,
                 "ARMED publication write made no progress")
        offset += count


def _read_exact_fd(descriptor: int, maximum: int, label: str) -> bytes:
    metadata = os.fstat(descriptor)
    _require(stat.S_ISREG(metadata.st_mode) and
             0 < metadata.st_size <= maximum,
             f"{label} is not one bounded regular file")
    chunks: list[bytes] = []
    offset = 0
    while offset < metadata.st_size:
        block = os.pread(
            descriptor, min(1 << 20, metadata.st_size - offset), offset)
        _require(bool(block), f"{label} read made no progress")
        chunks.append(block)
        offset += len(block)
    return b"".join(chunks)


def _read_immutable_file_at(
    directory_fd: int, name: str, maximum: int, label: str,
) -> tuple[int, bytes]:
    descriptor = -1
    try:
        descriptor = os.open(
            name, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        metadata = os.fstat(descriptor)
        path_metadata = os.stat(
            name, dir_fd=directory_fd, follow_symlinks=False)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) == 0o400 and
                 (metadata.st_dev, metadata.st_ino) ==
                    (path_metadata.st_dev, path_metadata.st_ino),
                 f"{label} is not one immutable owner record")
        stable = (metadata.st_dev, metadata.st_ino, metadata.st_mode,
                  metadata.st_uid, metadata.st_gid, metadata.st_nlink,
                  metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns)
        data = _read_exact_fd(descriptor, maximum, label)
        after = os.fstat(descriptor)
        _require((after.st_dev, after.st_ino, after.st_mode,
                  after.st_uid, after.st_gid, after.st_nlink,
                  after.st_size, after.st_mtime_ns, after.st_ctime_ns) == stable,
                 f"{label} changed while read")
        retained = descriptor
        descriptor = -1
        return retained, data
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _write_immutable_file_at(
    directory_fd: int, name: str, content: bytes, label: str,
) -> dict[str, Any]:
    _require(type(content) is bytes and 0 < len(content) <= MAX_JSON_BYTES,
             f"{label} content exceeds its byte bound")
    descriptor = -1
    try:
        descriptor = os.open(
            name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW, 0o600, dir_fd=directory_fd)
        _write_all(descriptor, content)
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        metadata = os.fstat(descriptor)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) == 0o400 and
                 metadata.st_size == len(content),
                 f"{label} immutable metadata differs")
        return {"name": name, "sha256": _sha256_bytes(content),
                "size": len(content)}
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _renameat2_noreplace(
    directory_fd: int, source_name: str, destination_name: str,
) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    _require(renameat2 is not None,
             "atomic attempt publication requires Linux renameat2")
    renameat2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint)
    renameat2.restype = ctypes.c_int
    result = renameat2(
        directory_fd, os.fsencode(source_name), directory_fd,
        os.fsencode(destination_name), 1)
    if result != 0:
        error_number = ctypes.get_errno()
        if error_number == errno.EEXIST:
            raise ArmingError("atomic publication target already exists")
        raise ArmingError(
            f"cannot atomically publish authority record: "
            f"{os.strerror(error_number)}")


def _write_atomic_journal_file_at(
    journal_fd: int, name: str, content: bytes, label: str,
) -> None:
    _require(
        (name == "complete.json" or
         JOURNAL_INTENT_NAME.fullmatch(name) is not None or
         JOURNAL_RESULT_NAME.fullmatch(name) is not None) and
        type(content) is bytes and 0 < len(content) <= MAX_JSON_BYTES,
        f"{label} publication input differs")
    staging_name = f".{name}.staging-{os.getpid()}-{secrets.token_hex(16)}"
    _require(JOURNAL_STAGING_NAME.fullmatch(staging_name) is not None,
             f"{label} staging name differs")
    staging_fd = -1
    published = False
    try:
        staging_fd = os.open(
            staging_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW,
            0o600, dir_fd=journal_fd)
        _write_all(staging_fd, content)
        os.fsync(staging_fd)
        os.fchmod(staging_fd, 0o400)
        os.fsync(staging_fd)
        metadata = os.fstat(staging_fd)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) == 0o400 and
                 metadata.st_size == len(content),
                 f"{label} staging metadata differs")
        _renameat2_noreplace(journal_fd, staging_name, name)
        published = True
        os.fsync(journal_fd)
        final_fd, final_data = _read_immutable_file_at(
            journal_fd, name, MAX_JSON_BYTES, label)
        try:
            _require(final_data == content,
                     f"{label} bytes changed during atomic publication")
        finally:
            os.close(final_fd)
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors: list[BaseException] = []
        if staging_fd >= 0:
            try:
                os.close(staging_fd)
            except BaseException as error:
                cleanup_errors.append(error)
        if not published:
            try:
                os.unlink(staging_name, dir_fd=journal_fd)
                os.fsync(journal_fd)
            except FileNotFoundError:
                pass
            except BaseException as error:
                cleanup_errors.append(error)
        if cleanup_errors and active_error is None:
            raise ArmingError(f"{label} staging cleanup failed") \
                from cleanup_errors[0]


def _lane_binding_record(
    lane: Path, registration: Mapping[str, Any], *, lane_nonce: str,
    arming_lock: Mapping[str, Any], attempts_directory: Mapping[str, Any],
    campaign_markers: Sequence[Mapping[str, Any]],
    campaign_progress_directories: Sequence[Mapping[str, Any]],
    prearmed_history_markers: Sequence[Mapping[str, Any]],
    prearmed_boundary_markers: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    _require(re.fullmatch(r"[0-9a-f]{64}", lane_nonce) is not None,
             "generation-3 lane nonce differs")
    return {
        "schema": LANE_BINDING_SCHEMA,
        "generation": GENERATION,
        "lane_nonce": lane_nonce,
        "absolute_path": str(lane),
        "budgets": {
            "setup_invalid": registration["budgets"]["setup_invalid"],
            "environment_rejected":
                registration["budgets"]["environment_rejected"],
            "evidence_attempts":
                registration["budgets"]["evidence_attempts"],
        },
        "arming_lock": _validate_lane_lock_binding(arming_lock),
        "attempts_directory": _validate_lane_attempts_binding(
            attempts_directory),
        "campaign_markers": [
            _validate_campaign_marker_binding(marker)
            for marker in campaign_markers
        ],
        "campaign_progress_directories": [
            _validate_campaign_progress_directory_binding(directory)
            for directory in campaign_progress_directories
        ],
        "prearmed_history_markers": [
            _validate_prearmed_history_marker_binding(marker)
            for marker in prearmed_history_markers
        ],
        "prearmed_boundary_markers": [
            _validate_prearmed_boundary_marker_binding(marker)
            for marker in prearmed_boundary_markers
        ],
        "publication_rule": PUBLICATION_RULE,
    }


def _prearmed_history_limit(registration: Mapping[str, Any]) -> int:
    budgets = registration["budgets"]
    return (budgets["setup_invalid"] + budgets["environment_rejected"] +
            budgets["evidence_attempts"])


def _validate_lane_binding(
    value: Any, lane: Path, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], candidate_authority: Mapping[str, Any],
    exact_main: Mapping[str, Any],
) -> dict[str, Any]:
    _exact_mapping(value, frozenset((
        "schema", "generation", "lane_nonce", "absolute_path",
        "budgets", "arming_lock", "attempts_directory", "campaign_markers",
        "campaign_progress_directories",
        "prearmed_history_markers", "prearmed_boundary_markers",
        "publication_rule",
    )), "generation-3 lane binding")
    _require(contract.exact_json_equal(value["budgets"], {
        "setup_invalid": registration["budgets"]["setup_invalid"],
        "environment_rejected":
            registration["budgets"]["environment_rejected"],
        "evidence_attempts": registration["budgets"]["evidence_attempts"],
    }), "generation-3 lane budget topology differs")
    arming_lock = _validate_lane_lock_binding(value["arming_lock"])
    attempts_directory = _validate_lane_attempts_binding(
        value["attempts_directory"])
    _require(type(value["campaign_markers"]) is list and
             len(value["campaign_markers"]) ==
                registration["budgets"]["evidence_attempts"],
             "generation-3 campaign marker inventory differs")
    markers = [
        _validate_campaign_marker_binding(marker)
        for marker in value["campaign_markers"]
    ]
    _require(
        [marker["name"] for marker in markers] == [
            f"attempt-{index:03d}-transcript.jsonl"
            for index in range(1, len(markers) + 1)
        ],
        "generation-3 campaign markers are not the exact slot inventory")
    _require(type(value["campaign_progress_directories"]) is list and
             len(value["campaign_progress_directories"]) ==
                registration["budgets"]["evidence_attempts"],
             "generation-3 campaign progress inventory differs")
    progress_directories = [
        _validate_campaign_progress_directory_binding(directory)
        for directory in value["campaign_progress_directories"]
    ]
    _require(
        [directory["name"] for directory in progress_directories] == [
            f"attempt-{index:03d}-journal-checkpoints"
            for index in range(1, len(progress_directories) + 1)
        ],
        "generation-3 campaign progress directories are not the exact slot "
        "inventory")
    history_limit = _prearmed_history_limit(registration)
    _require(type(value["prearmed_history_markers"]) is list and
             len(value["prearmed_history_markers"]) == history_limit,
             "generation-3 pre-ARMED history marker inventory differs")
    history_markers = [
        _validate_prearmed_history_marker_binding(marker)
        for marker in value["prearmed_history_markers"]
    ]
    _require(
        [marker["name"] for marker in history_markers] == [
            f"prearmed-{index:04d}-history.jsonl"
            for index in range(1, history_limit + 1)
        ],
        "generation-3 pre-ARMED histories are not the exact slot inventory")
    _require(type(value["prearmed_boundary_markers"]) is list and
             len(value["prearmed_boundary_markers"]) ==
                history_limit * len(PREARMED_BOUNDARY_STATES),
             "generation-3 pre-ARMED boundary inventory differs")
    boundary_markers = [
        _validate_prearmed_boundary_marker_binding(marker)
        for marker in value["prearmed_boundary_markers"]
    ]
    _require(
        [marker["name"] for marker in boundary_markers] == [
            f"prearmed-{index:04d}-{state.lower()}-reached.marker"
            for index in range(1, history_limit + 1)
            for state in PREARMED_BOUNDARY_STATES
        ],
        "generation-3 pre-ARMED boundaries are not the exact slot inventory")
    expected = _lane_binding_record(
        lane, registration, lane_nonce=value["lane_nonce"],
        arming_lock=arming_lock, attempts_directory=attempts_directory,
        campaign_markers=markers,
        campaign_progress_directories=progress_directories,
        prearmed_history_markers=history_markers,
        prearmed_boundary_markers=boundary_markers)
    _require(contract.exact_json_equal(value, expected),
             "generation-3 lane binding differs")
    manifest_binding = registration["output_lane"]["lane_manifest"]
    manifest_data = contract.canonical_json_bytes(expected)
    _require(len(manifest_data) == manifest_binding["size"] and
             _sha256_bytes(manifest_data) == manifest_binding["sha256"],
             "generation-3 lane binding differs from preregistered bytes")
    return copy.deepcopy(expected)


def _validate_attempt_staging_inventory(
    attempts_fd: int, name: str, *, remove_complete: bool = True,
) -> None:
    staging_fd = -1
    journal_fd = -1
    try:
        try:
            staging_fd = os.open(
                name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=attempts_fd)
        except OSError as error:
            raise ArmingError(
                "stale attempt staging directory cannot be opened") from error
        metadata = os.fstat(staging_fd)
        _require(stat.S_ISDIR(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 stat.S_IMODE(metadata.st_mode) in (0o700, 0o500),
                 "stale attempt staging directory is unsafe")
        children = sorted(os.listdir(staging_fd))
        _require(set(children) <= {
            AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
            JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE},
            "stale attempt staging directory contains an unsafe entry")
        _require(remove_complete or not (
            stat.S_IMODE(metadata.st_mode) == 0o500 and
            set(children) == {
                AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
                JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE}),
            "complete attempt staging cannot be classified as unpublished")
        for child in children:
            child_metadata = os.stat(
                child, dir_fd=staging_fd, follow_symlinks=False)
            if child == JOURNAL_DIRECTORY:
                _require(stat.S_ISDIR(child_metadata.st_mode) and
                         child_metadata.st_uid == os.geteuid() and
                         stat.S_IMODE(child_metadata.st_mode) == 0o700,
                         "stale attempt journal directory is unsafe")
                journal_fd = os.open(
                    child,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC |
                    os.O_NOFOLLOW,
                    dir_fd=staging_fd)
                _require(not os.listdir(journal_fd),
                         "stale unpublished journal is not empty")
                os.close(journal_fd)
                journal_fd = -1
                continue
            _require(stat.S_ISREG(child_metadata.st_mode) and
                     child_metadata.st_uid == os.geteuid() and
                     child_metadata.st_nlink == 1 and
                     stat.S_IMODE(child_metadata.st_mode) in (0o600, 0o400) and
                     0 <= child_metadata.st_size <= MAX_JSON_BYTES,
                     "stale attempt staging file is unsafe")
            if child == LAUNCH_CONSUMED_FILE:
                _require(child_metadata.st_size == 0,
                         "stale unpublished launch marker is consumed")
    finally:
        if journal_fd >= 0:
            os.close(journal_fd)
        if staging_fd >= 0:
            os.close(staging_fd)


def _preflight_lane_inventory(
    lane_fd: int, evidence_attempt_limit: int,
    prearmed_attempt_limit: int, *, lock_required: bool,
) -> None:
    _require(type(evidence_attempt_limit) is int and
             evidence_attempt_limit > 0,
             "generation-3 evidence-attempt limit is invalid")
    _require(type(prearmed_attempt_limit) is int and
             prearmed_attempt_limit > evidence_attempt_limit,
             "generation-3 pre-ARMED attempt limit is invalid")
    names = set(os.listdir(lane_fd))
    if not lock_required:
        _require(not names,
                 "new generation-3 lane must be empty before lock creation")
        return
    _require(".arming.lock" in names,
             "existing generation-3 lane has no arming lock")
    for name in names:
        marker = CAMPAIGN_CONSUMED_NAME.fullmatch(name)
        progress = CAMPAIGN_PROGRESS_DIRECTORY_NAME.fullmatch(name)
        history = PREARMED_HISTORY_NAME.fullmatch(name)
        boundary = PREARMED_BOUNDARY_NAME.fullmatch(name)
        allowed = (
            name in {".arming.lock", LANE_BINDING_FILE, ATTEMPTS_DIRECTORY} or
            (marker is not None and
             1 <= int(marker.group(1)) <= evidence_attempt_limit) or
            (progress is not None and
             1 <= int(progress.group(1)) <= evidence_attempt_limit) or
            (history is not None and
             1 <= int(history.group(1)) <= prearmed_attempt_limit) or
            (boundary is not None and
             1 <= int(boundary.group(1)) <= prearmed_attempt_limit)
        )
        _require(allowed,
                 "generation-3 lane contains an unrecognized entry")
    lane_binding_mode: int | None = None
    if LANE_BINDING_FILE in names:
        metadata = os.stat(
            LANE_BINDING_FILE, dir_fd=lane_fd, follow_symlinks=False)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) in (0o400, 0o600) and
                 0 <= metadata.st_size <= MAX_JSON_BYTES,
                 "generation-3 lane binding is unsafe")
        lane_binding_mode = stat.S_IMODE(metadata.st_mode)
    for name in sorted(
            entry for entry in names
            if CAMPAIGN_CONSUMED_NAME.fullmatch(entry) is not None):
        metadata = os.stat(name, dir_fd=lane_fd, follow_symlinks=False)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) in (0o400, 0o600) and
                 0 <= metadata.st_size <= MAX_CAMPAIGN_TRANSCRIPT_BYTES and
                 (lane_binding_mode == 0o400 or metadata.st_size == 0),
                 "generation-3 campaign transcript slot is unsafe")
    for name in sorted(
            entry for entry in names
            if CAMPAIGN_PROGRESS_DIRECTORY_NAME.fullmatch(entry) is not None):
        descriptor = os.open(
            name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        try:
            metadata = os.fstat(descriptor)
            children = sorted(os.listdir(descriptor))
            _require(stat.S_ISDIR(metadata.st_mode) and
                     metadata.st_uid == os.geteuid() and
                     metadata.st_nlink == 2 and
                     stat.S_IMODE(metadata.st_mode) == 0o500 and
                     len(children) ==
                        2 * prereg.OUTPUT_LANE_CHILD_PROCESS_COUNT + 2 and
                     all(JOURNAL_CHECKPOINT_NAME.fullmatch(child) is not None
                         for child in children),
                     "generation-3 campaign checkpoint directory is unsafe")
        finally:
            os.close(descriptor)
    for name in sorted(
            entry for entry in names
            if PREARMED_HISTORY_NAME.fullmatch(entry) is not None):
        metadata = os.stat(name, dir_fd=lane_fd, follow_symlinks=False)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) in (0o400, 0o600) and
                 0 <= metadata.st_size <= MAX_PREARMED_HISTORY_BYTES and
                 (lane_binding_mode == 0o400 or metadata.st_size == 0),
                 "generation-3 pre-ARMED history slot is unsafe")
    for name in sorted(
            entry for entry in names
            if PREARMED_BOUNDARY_NAME.fullmatch(entry) is not None):
        metadata = os.stat(name, dir_fd=lane_fd, follow_symlinks=False)
        _require(stat.S_ISREG(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 metadata.st_nlink == 1 and
                 stat.S_IMODE(metadata.st_mode) in (0o400, 0o600) and
                 metadata.st_size == 0,
                 "generation-3 pre-ARMED boundary slot is unsafe")
    if ATTEMPTS_DIRECTORY not in names:
        return
    attempts_fd = -1
    try:
        try:
            attempts_fd = os.open(
                ATTEMPTS_DIRECTORY,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
        except OSError as error:
            raise ArmingError(
                "generation-3 attempts namespace cannot be preflighted") \
                from error
        metadata = os.fstat(attempts_fd)
        _require(stat.S_ISDIR(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 stat.S_IMODE(metadata.st_mode) in (0o500, 0o700),
                 "generation-3 attempts directory is unsafe")
        children = sorted(os.listdir(attempts_fd))
        final_matches = [
            ATTEMPT_NAME.fullmatch(child) for child in children
            if ATTEMPT_NAME.fullmatch(child) is not None]
        staging_matches = [
            STAGING_NAME.fullmatch(child) for child in children
            if STAGING_NAME.fullmatch(child) is not None]
        _require(len(final_matches) + len(staging_matches) == len(children),
            "generation-3 attempts namespace contains an unknown entry")
        final_indexes = sorted(int(match.group(1)) for match in final_matches)
        _require(
            final_indexes == list(range(1, len(final_indexes) + 1)) and
            len(final_indexes) <= evidence_attempt_limit,
            "generation-3 final attempts are not a gapless in-budget prefix")
        _require(len(staging_matches) <= 1,
                 "generation-3 attempts namespace has multiple staging "
                 "remnants")
        if staging_matches:
            staging_index = int(staging_matches[0].group(1))
            _require(staging_index == len(final_indexes) + 1 and
                     staging_index <= evidence_attempt_limit,
                     "generation-3 staging remnant is not the unique next "
                     "in-budget attempt")
        for child, match in zip(children, (
                ATTEMPT_NAME.fullmatch(name) for name in children)):
            if match is None:
                continue
            child_metadata = os.stat(
                child, dir_fd=attempts_fd, follow_symlinks=False)
            _require(stat.S_ISDIR(child_metadata.st_mode) and
                     child_metadata.st_uid == os.geteuid() and
                     stat.S_IMODE(child_metadata.st_mode) == 0o500,
                     "generation-3 final attempt directory is unsafe")
        for child in children:
            if STAGING_NAME.fullmatch(child) is not None:
                _validate_attempt_staging_inventory(attempts_fd, child)
    finally:
        if attempts_fd >= 0:
            os.close(attempts_fd)


def _validate_open_lane_lock_identity(
    lane: Path, lane_fd: int, lock_fd: int, *, allowed_lane_modes: frozenset[int],
) -> None:
    try:
        retained_lane = os.fstat(lane_fd)
        current_lane = os.lstat(lane)
        retained_lock = os.fstat(lock_fd)
        current_lock = os.stat(
            ".arming.lock", dir_fd=lane_fd, follow_symlinks=False)
    except OSError as error:
        raise ArmingError(
            "generation-3 lane or lock identity cannot be revalidated") \
            from error
    _require(
        stat.S_ISDIR(retained_lane.st_mode) and
        retained_lane.st_uid == os.geteuid() and
        stat.S_IMODE(retained_lane.st_mode) in allowed_lane_modes and
        (retained_lane.st_dev, retained_lane.st_ino) ==
            (current_lane.st_dev, current_lane.st_ino) and
        stat.S_ISREG(retained_lock.st_mode) and
        retained_lock.st_uid == os.geteuid() and
        retained_lock.st_nlink == 1 and
        stat.S_IMODE(retained_lock.st_mode) == 0o600 and
        (retained_lock.st_dev, retained_lock.st_ino) ==
            (current_lock.st_dev, current_lock.st_ino),
        "generation-3 lane or lock changed while acquiring authority")


def _open_lane_and_lock(
    lane_path: Path, lane_fd: int, registration: Mapping[str, Any], *,
    evidence_attempt_limit: int, prearmed_attempt_limit: int,
) -> tuple[Path, int, int]:
    lock_fd = -1
    try:
        lane = _canonical_lane(lane_path)
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane, lane_fd, registration)
        lane_metadata = os.fstat(lane_fd)
        path_metadata = os.lstat(lane)
        _require((lane_metadata.st_dev, lane_metadata.st_ino) ==
                 (path_metadata.st_dev, path_metadata.st_ino),
                 "generation-3 lane changed while opening")
        initial_names = set(os.listdir(lane_fd))
        lock_exists = ".arming.lock" in initial_names
        _require(lock_exists,
                 "pre-ratified generation-3 lane has no arming lock")
        _preflight_lane_inventory(
            lane_fd, evidence_attempt_limit, prearmed_attempt_limit,
            lock_required=True)
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane, lane_fd, registration)
        try:
            lock_fd = os.open(
                ".arming.lock", os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
        except OSError as error:
            raise ArmingError("generation-3 lane lock cannot be opened") \
                from error
        lock_metadata = os.fstat(lock_fd)
        _require(stat.S_ISREG(lock_metadata.st_mode) and
                 lock_metadata.st_uid == os.geteuid() and
                 lock_metadata.st_nlink == 1 and
                 stat.S_IMODE(lock_metadata.st_mode) == 0o600,
                 "generation-3 lane lock is unsafe")
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane, lane_fd, registration)
        _validate_open_lane_lock_identity(
            lane, lane_fd, lock_fd, allowed_lane_modes=frozenset((0o500,)))
        _preflight_lane_inventory(
            lane_fd, evidence_attempt_limit, prearmed_attempt_limit,
            lock_required=True)
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane, lane_fd, registration)
        _validate_open_lane_lock_identity(
            lane, lane_fd, lock_fd, allowed_lane_modes=frozenset((0o500,)))
        return lane, lane_fd, lock_fd
    except BaseException:
        if lock_fd >= 0:
            os.close(lock_fd)
        os.close(lane_fd)
        raise


def _validate_lane_capability(
    lane: Path, lane_fd: int, lock_fd: int,
    registration: Mapping[str, Any],
) -> None:
    try:
        retained_lane = os.fstat(lane_fd)
        path_lane = os.lstat(lane)
        retained_lock = os.fstat(lock_fd)
        path_lock = os.stat(
            ".arming.lock", dir_fd=lane_fd, follow_symlinks=False)
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        raise ArmingError("generation-3 lane capability is no longer held") \
            from error
    _require(
        stat.S_ISDIR(retained_lane.st_mode) and
        retained_lane.st_uid == os.geteuid() and
        stat.S_IMODE(retained_lane.st_mode) == 0o500 and
        (retained_lane.st_dev, retained_lane.st_ino) ==
            (path_lane.st_dev, path_lane.st_ino) and
        stat.S_ISREG(retained_lock.st_mode) and
        retained_lock.st_uid == os.geteuid() and
        retained_lock.st_nlink == 1 and
        stat.S_IMODE(retained_lock.st_mode) == 0o600 and
        (retained_lock.st_dev, retained_lock.st_ino) ==
            (path_lock.st_dev, path_lock.st_ino),
        "generation-3 lane or lock inode changed")
    prereg.validate_output_lane_descriptor_identity_for_preregistration(
        lane, lane_fd, registration)


def _attempts_binding_record(metadata: os.stat_result) -> dict[str, Any]:
    _require(stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             stat.S_IMODE(metadata.st_mode) == 0o500,
             "generation-3 attempts namespace is not sealed")
    return {
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o500,
    }


def _validate_attempts_capability(
    lane_fd: int, attempts_fd: int, binding: Mapping[str, Any],
    attempt_fd: int, attempt_name: str,
) -> None:
    expected = _exact_mapping(
        binding, frozenset(("device", "inode", "uid", "mode")),
        "retained attempts namespace binding")
    try:
        retained = os.fstat(attempts_fd)
        pathname = os.stat(
            ATTEMPTS_DIRECTORY, dir_fd=lane_fd, follow_symlinks=False)
        retained_attempt = os.fstat(attempt_fd)
        pathname_attempt = os.stat(
            attempt_name, dir_fd=attempts_fd, follow_symlinks=False)
    except OSError as error:
        raise ArmingError("retained attempts namespace changed") from error
    observed = _attempts_binding_record(retained)
    _require(contract.exact_json_equal(observed, expected) and
             stat.S_ISDIR(pathname.st_mode) and
             pathname.st_uid == expected["uid"] and
             stat.S_IMODE(pathname.st_mode) == expected["mode"] and
             (retained.st_dev, retained.st_ino) ==
                (pathname.st_dev, pathname.st_ino) and
             (retained_attempt.st_dev, retained_attempt.st_ino) ==
                (pathname_attempt.st_dev, pathname_attempt.st_ino),
             "retained attempts namespace or attempt inode changed")


def _validate_retained_attempt_capability(
    attempt_path: Path, attempt_fd: int, journal_fd: int,
    journal_binding: Mapping[str, Any], launch_marker_fd: int,
    launch_marker_binding: Mapping[str, Any], armed_fd: int,
    bundle_fd: int, manifest_fd: int, armed_data: bytes,
    bundle_data: bytes, manifest_data: bytes,
) -> None:
    try:
        retained = os.fstat(attempt_fd)
        pathname = os.lstat(attempt_path)
        armed_status = os.fstat(armed_fd)
        armed_path_status = os.stat(
            ARMED_FILE, dir_fd=attempt_fd, follow_symlinks=False)
        bundle_status = os.fstat(bundle_fd)
        bundle_path_status = os.stat(
            AUTHORITY_BUNDLE_FILE, dir_fd=attempt_fd,
            follow_symlinks=False)
        manifest_status = os.fstat(manifest_fd)
        manifest_path_status = os.stat(
            ATTEMPT_MANIFEST_FILE, dir_fd=attempt_fd,
            follow_symlinks=False)
    except OSError as error:
        raise ArmingError("retained attempt capability changed") from error
    _require(
        stat.S_ISDIR(retained.st_mode) and retained.st_uid == os.geteuid() and
        stat.S_IMODE(retained.st_mode) == 0o500 and
        set(os.listdir(attempt_fd)) == {
            AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
            JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE} and
        (retained.st_dev, retained.st_ino) ==
            (pathname.st_dev, pathname.st_ino) and
        all(stat.S_ISREG(status.st_mode) and
            status.st_uid == os.geteuid() and status.st_nlink == 1 and
            stat.S_IMODE(status.st_mode) == 0o400
            for status in (armed_status, bundle_status, manifest_status)) and
        (armed_status.st_dev, armed_status.st_ino) ==
            (armed_path_status.st_dev, armed_path_status.st_ino) and
        (bundle_status.st_dev, bundle_status.st_ino) ==
            (bundle_path_status.st_dev, bundle_path_status.st_ino) and
        (manifest_status.st_dev, manifest_status.st_ino) ==
            (manifest_path_status.st_dev, manifest_path_status.st_ino) and
        _read_exact_fd(armed_fd, MAX_JSON_BYTES, "retained ARMED record") ==
            armed_data and
        _read_exact_fd(bundle_fd, MAX_JSON_BYTES, "retained authority bundle") ==
            bundle_data and
        _read_exact_fd(manifest_fd, MAX_JSON_BYTES, "retained manifest") ==
            manifest_data,
        "retained attempt directory or authority bytes changed")
    _validate_journal_capability(
        attempt_fd, journal_fd, journal_binding,
        "retained execution journal")
    _validate_launch_marker_capability(
        attempt_fd, launch_marker_fd, launch_marker_binding,
        "retained launch-consumed marker")


def _validate_preallocated_lane_lock(
    lane_fd: int, lock_fd: int, binding_value: Mapping[str, Any],
) -> None:
    binding = _validate_lane_lock_binding(binding_value)
    retained = os.fstat(lock_fd)
    pathname = os.stat(
        binding["name"], dir_fd=lane_fd, follow_symlinks=False)
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        stat.S_IMODE(retained.st_mode) == binding["mode"] and
        retained.st_size == 0 and
        retained.st_mtime_ns == binding["initial_mtime_ns"] and
        retained.st_ctime_ns == binding["initial_ctime_ns"] and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino),
        "preallocated generation-3 lane lock changed")


def _validate_preallocated_attempts_directory(
    lane: Path, lane_fd: int, attempts_fd: int,
    binding_value: Mapping[str, Any],
) -> None:
    binding = _validate_lane_attempts_binding(binding_value)
    retained = os.fstat(attempts_fd)
    pathname = os.stat(
        binding["name"], dir_fd=lane_fd, follow_symlinks=False)
    retained_handle, retained_mount_id = \
        prereg._capture_output_lane_file_handle(
            lane_fd, os.fsencode(binding["name"]),
            "retained preallocated attempts directory")
    pathname_handle, pathname_mount_id = \
        prereg._capture_output_lane_file_handle(
            -100,
            os.fsencode(lane / binding["name"]),
            "path preallocated attempts directory")
    _require(
        stat.S_ISDIR(retained.st_mode) and
        retained.st_uid == binding["uid"] and
        stat.S_IMODE(retained.st_mode) in (0o500, 0o700) and
        retained.st_nlink >= binding["initial_link_count"] and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        retained_mount_id == pathname_mount_id and
        contract.exact_json_equal(retained_handle, binding["file_handle"]) and
        contract.exact_json_equal(pathname_handle, binding["file_handle"]),
        "preallocated generation-3 attempts directory changed")


def _expected_lane_inventory(
    binding: Mapping[str, Any],
) -> set[str]:
    return {
        ".arming.lock", LANE_BINDING_FILE, ATTEMPTS_DIRECTORY,
        *(marker["name"] for marker in binding["campaign_markers"]),
        *(directory["name"] for directory in
          binding["campaign_progress_directories"]),
        *(marker["name"] for marker in
          binding["prearmed_history_markers"]),
        *(marker["name"] for marker in
          binding["prearmed_boundary_markers"]),
    }


def _initialize_lane_locked(
    lane: Path, lane_fd: int, lock_fd: int, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], candidate_authority: Mapping[str, Any],
    exact_main: Mapping[str, Any],
) -> dict[str, Any]:
    """Replay the preregistered immutable topology without root mutation."""
    _validate_open_lane_lock_identity(
        lane, lane_fd, lock_fd, allowed_lane_modes=frozenset((0o500,)))
    prereg.validate_output_lane_descriptor_identity_for_preregistration(
        lane, lane_fd, registration)
    limit = registration["budgets"]["evidence_attempts"]
    history_limit = _prearmed_history_limit(registration)
    _preflight_lane_inventory(
        lane_fd, limit, history_limit, lock_required=True)
    binding_fd, binding_data = _read_immutable_file_at(
        lane_fd, LANE_BINDING_FILE, MAX_JSON_BYTES,
        "pre-ratified generation-3 lane binding")
    try:
        binding = _validate_lane_binding(
            _canonical_json_document(
                binding_data, "pre-ratified generation-3 lane binding"),
            lane, registration, logical_plan, candidate_authority, exact_main)
    finally:
        os.close(binding_fd)
    _require(set(os.listdir(lane_fd)) == _expected_lane_inventory(binding),
             "pre-ratified generation-3 lane inventory differs")
    _validate_preallocated_lane_lock(
        lane_fd, lock_fd, binding["arming_lock"])
    attempts_fd = os.open(
        ATTEMPTS_DIRECTORY,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=lane_fd)
    try:
        _validate_preallocated_attempts_directory(
            lane, lane_fd, attempts_fd, binding["attempts_directory"])
        unused_names, prior, prepared = _read_prior_attempts_before_recovery(
            attempts_fd, lane=lane, lane_fd=lane_fd,
            registration=registration, logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=binding)
        del unused_names
        authority = _read_budget_authority_before_recovery(
            lane_fd, registration=registration, logical_plan=logical_plan,
            lane_binding=binding, prior_attempts=[
                *prior,
                *([] if prepared is None else [prepared[1]]),
            ])
        _recover_prearmed_histories(
            lane_fd, authority, registration=registration,
            logical_plan=logical_plan, lane_binding=binding)
    finally:
        os.close(attempts_fd)
    prereg.validate_output_lane_descriptor_identity_for_preregistration(
        lane, lane_fd, registration)
    _validate_preallocated_lane_lock(
        lane_fd, lock_fd, binding["arming_lock"])
    return binding


def _qualification_evidence_document(
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "policy_a_scan": copy.deepcopy(evidence["policy_a_scan"]),
        "policy_b_scan": copy.deepcopy(evidence["policy_b_scan"]),
        "track_selection": copy.deepcopy(evidence["track_selection"]),
        "expected_frozen_pair": copy.deepcopy(evidence["expected_frozen_pair"]),
        "acquisition": contract.strict_json_loads(
            evidence["acquisition_data"], "qualification acquisition"),
        "bridge": contract.strict_json_loads(
            evidence["bridge_data"], "qualification bridge"),
        "independent_verdict": contract.strict_json_loads(
            evidence["independent_verdict_data"],
            "independent qualification verdict"),
    }


def _qualification_evidence_wire(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "policy_a_scan", "policy_b_scan", "track_selection",
        "expected_frozen_pair", "acquisition", "bridge",
        "independent_verdict",
    )), "retained qualification evidence")
    return {
        "policy_a_scan": copy.deepcopy(record["policy_a_scan"]),
        "policy_b_scan": copy.deepcopy(record["policy_b_scan"]),
        "track_selection": copy.deepcopy(record["track_selection"]),
        "expected_frozen_pair": copy.deepcopy(record["expected_frozen_pair"]),
        "acquisition_data": contract.canonical_json_bytes(record["acquisition"]),
        "bridge_data": contract.canonical_json_bytes(record["bridge"]),
        "independent_verdict_data": contract.canonical_json_bytes(
            record["independent_verdict"]),
    }


def _authority_bundle_record(
    *, lane_binding: Mapping[str, Any], preregistration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], host_instance: Mapping[str, Any],
    launch_context: Mapping[str, Any],
    candidate_authority_record: Mapping[str, Any],
    exact_main_authority_record: Mapping[str, Any],
    exact_main_verifier: Mapping[str, Any], artifact_bundle: Mapping[str, Any],
    qualification_binding: Mapping[str, Any],
    qualification_evidence: Mapping[str, Any],
    candidate_descriptor: Mapping[str, Any],
    main_descriptor: Mapping[str, Any], budget_commit: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": AUTHORITY_BUNDLE_SCHEMA,
        "lane_binding": copy.deepcopy(lane_binding),
        "preregistration": copy.deepcopy(preregistration),
        "logical_plan": copy.deepcopy(logical_plan),
        "host_instance": copy.deepcopy(host_instance),
        "launch_context": copy.deepcopy(
            _validate_launch_context(launch_context)),
        "candidate_authority_record": copy.deepcopy(
            candidate_authority_record),
        "exact_main_authority_record": copy.deepcopy(
            exact_main_authority_record),
        "exact_main_verifier": copy.deepcopy(exact_main_verifier),
        "artifact_bundle": copy.deepcopy(artifact_bundle),
        "qualification_binding": copy.deepcopy(qualification_binding),
        "qualification_evidence": _qualification_evidence_document(
            qualification_evidence),
        "budget_commit": copy.deepcopy(budget_commit),
        "descriptor_binding": {
            "candidate_control": copy.deepcopy(candidate_descriptor),
            "main": copy.deepcopy(main_descriptor),
        },
    }


def _bundle_components(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "schema", "lane_binding", "preregistration", "logical_plan",
        "host_instance", "launch_context", "candidate_authority_record",
        "exact_main_authority_record", "exact_main_verifier",
        "artifact_bundle", "qualification_binding", "qualification_evidence",
        "budget_commit", "descriptor_binding",
    )), "generation-3 authority bundle")
    _require(record["schema"] == AUTHORITY_BUNDLE_SCHEMA,
             "generation-3 authority bundle schema differs")
    descriptors = _exact_mapping(
        record["descriptor_binding"], frozenset(("candidate_control", "main")),
        "generation-3 descriptor binding")
    launch_context = _validate_launch_context(record["launch_context"])
    return {**record, "descriptor_binding": descriptors,
            "launch_context": launch_context,
            "qualification_evidence_wire": _qualification_evidence_wire(
                record["qualification_evidence"])}


def _validate_descriptor_bundle(
    descriptors: Mapping[str, Any], artifacts: Mapping[str, Any],
) -> None:
    roles = artifacts["roles"]
    keys = frozenset((
        "protocol", "boot_id", "session_nonce", "creator_pid",
        "creator_start_ticks", "device", "inode", "uid", "gid", "size",
        "mode", "raw_sha256", "seals", "elf",
    ))

    def validate(value: Any, label: str) -> dict[str, Any]:
        identity = _exact_mapping(
            value, keys, f"retained {label} descriptor binding")
        _require(
            identity["protocol"] == SEALED_EXECUTABLE_PROTOCOL and
            type(identity["boot_id"]) is str and
            re.fullmatch(
                r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-"
                r"[0-9a-f]{4}-[0-9a-f]{12}",
                identity["boot_id"]) is not None and
            type(identity["session_nonce"]) is str and
            re.fullmatch(r"[0-9a-f]{64}", identity["session_nonce"])
                is not None and
            type(identity["creator_pid"]) is int and
            identity["creator_pid"] > 0 and
            type(identity["creator_start_ticks"]) is int and
            identity["creator_start_ticks"] > 0 and
            type(identity["device"]) is int and identity["device"] >= 0 and
            type(identity["inode"]) is int and identity["inode"] > 0 and
            type(identity["uid"]) is int and
            identity["uid"] == os.geteuid() and
            type(identity["gid"]) is int and identity["gid"] >= 0 and
            type(identity["size"]) is int and
            0 < identity["size"] <= MAX_EXECUTABLE_BYTES and
            type(identity["mode"]) is int and identity["mode"] == 0o500 and
            type(identity["raw_sha256"]) is str and
            re.fullmatch(r"[0-9a-f]{64}", identity["raw_sha256"])
                is not None and
            type(identity["seals"]) is int and
            identity["seals"] == REQUIRED_SEALS and
            identity["elf"] is True,
            f"retained {label} descriptor binding is malformed")
        return identity

    candidate = validate(descriptors["candidate_control"], "candidate")
    main = validate(descriptors["main"], "main")
    candidate_handle = contract.canonical_sha256(candidate)
    main_handle = contract.canonical_sha256(main)
    _require(
        roles["candidate"]["handle_id"] ==
            roles["control"]["handle_id"] == candidate_handle and
        roles["candidate"]["handle_device"] == candidate["device"] and
        roles["candidate"]["handle_inode"] == candidate["inode"] and
        roles["main"]["handle_id"] == main_handle and
        roles["main"]["handle_device"] == main["device"] and
        roles["main"]["handle_inode"] == main["inode"] and
        candidate["raw_sha256"] == roles["candidate"]["raw_sha256"] and
        candidate["size"] == roles["candidate"]["size"] ==
            roles["control"]["size"] and
        main["raw_sha256"] == roles["main"]["raw_sha256"] and
        main["size"] == roles["main"]["size"],
        "retained descriptor and artifact-role identities differ")


def _journal_binding_record(metadata: os.stat_result) -> dict[str, Any]:
    _require(stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             metadata.st_nlink == 2 and
             stat.S_IMODE(metadata.st_mode) == 0o700,
             "execution journal metadata is unsafe")
    return {
        "name": JOURNAL_DIRECTORY,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o700,
        "link_count": 2,
    }


def _launch_marker_binding_record(metadata: os.stat_result) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o400 and
             metadata.st_size == 0,
             "launch-consumed marker metadata is unsafe")
    return {
        "name": LAUNCH_CONSUMED_FILE,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o400,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _validate_launch_marker_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "launch-consumed marker binding")
    _require(record["name"] == LAUNCH_CONSUMED_FILE and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             type(record["mode"]) is int and record["mode"] == 0o400 and
             type(record["link_count"]) is int and
             record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "launch-consumed marker binding is malformed")
    return copy.deepcopy(record)


def _validate_recoverable_launch_marker_capability(
    attempt_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    binding = _validate_launch_marker_binding(expected)
    try:
        retained = os.fstat(marker_fd)
        pathname = os.stat(
            LAUNCH_CONSUMED_FILE, dir_fd=attempt_fd,
            follow_symlinks=False)
    except OSError as error:
        raise ArmingError(f"{label} changed") from error
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        stat.S_IMODE(retained.st_mode) in (binding["mode"], 0o600) and
        0 <= retained.st_size <= MAX_JSON_BYTES and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISREG(pathname.st_mode) and pathname.st_uid == binding["uid"] and
        pathname.st_nlink == 1 and
        stat.S_IMODE(pathname.st_mode) == stat.S_IMODE(retained.st_mode) and
        pathname.st_size == retained.st_size and
        (retained.st_size != 0 or
         (retained.st_mtime_ns == binding["initial_mtime_ns"] and
          retained.st_ctime_ns == binding["initial_ctime_ns"])),
        f"{label} inode or monotone state changed")
    return retained


def _validate_launch_marker_capability(
    attempt_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    metadata = _validate_recoverable_launch_marker_capability(
        attempt_fd, marker_fd, expected, label)
    _require(stat.S_IMODE(metadata.st_mode) == 0o400,
             f"{label} writable mode differs")
    return metadata


def _reconcile_launch_marker_prefix(
    attempt_fd: int, binding: Mapping[str, Any], *,
    expected_data: bytes | None, label: str,
) -> None:
    probe = os.open(
        LAUNCH_CONSUMED_FILE,
        os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=attempt_fd)
    writable = -1
    try:
        metadata = _validate_recoverable_launch_marker_capability(
            attempt_fd, probe, binding, label)
        observed = b"" if metadata.st_size == 0 else _read_exact_fd(
            probe, MAX_JSON_BYTES, label)
        if expected_data is None:
            _require(not observed and
                     stat.S_IMODE(metadata.st_mode) == binding["mode"],
                     f"{label} was consumed without a durable journal intent")
            # A pristine marker is itself the monotone no-launch witness.
            # Do not issue a same-mode chmod/fsync during replay: chmod would
            # advance its kernel-only ctime and erase the distinction between
            # untouched and coherently truncated post-launch state.
            os.close(probe)
            probe = -1
            return
        else:
            _require(expected_data.startswith(observed),
                     f"{label} is not an exact durable-intent prefix")
            if len(observed) < len(expected_data):
                os.fchmod(probe, 0o600)
                os.fsync(probe)
                writable = os.open(
                    LAUNCH_CONSUMED_FILE,
                    os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
                    dir_fd=attempt_fd)
                current = os.fstat(writable)
                _require(
                    (current.st_dev, current.st_ino) ==
                        (binding["device"], binding["inode"]) and
                    current.st_size == len(observed),
                    f"{label} changed while opening for recovery")
                offset = len(observed)
                while offset < len(expected_data):
                    count = os.pwrite(
                        writable, expected_data[offset:], offset)
                    _require(count > 0,
                             f"{label} prefix recovery made no progress")
                    offset += count
                os.fsync(writable)
        active = writable if writable >= 0 else probe
        os.fchmod(active, 0o400)
        os.fsync(active)
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors: list[BaseException] = []
        for descriptor in (writable, probe):
            if descriptor < 0:
                continue
            try:
                os.fchmod(descriptor, 0o400)
            except BaseException as error:
                cleanup_errors.append(error)
            try:
                os.fsync(descriptor)
            except BaseException as error:
                cleanup_errors.append(error)
            try:
                os.close(descriptor)
            except BaseException as error:
                cleanup_errors.append(error)
        if cleanup_errors and active_error is None:
            raise ArmingError(f"{label} recovery cleanup failed") \
                from cleanup_errors[0]
    descriptor = os.open(
        LAUNCH_CONSUMED_FILE,
        os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=attempt_fd)
    try:
        metadata = _validate_launch_marker_capability(
            attempt_fd, descriptor, binding, label)
        observed = b"" if metadata.st_size == 0 else _read_exact_fd(
            descriptor, MAX_JSON_BYTES, label)
        _require(observed == (b"" if expected_data is None else expected_data),
                 f"{label} recovery bytes differ")
    finally:
        os.close(descriptor)


def _campaign_marker_binding_record(
    metadata: os.stat_result, name: str,
) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o400 and
             metadata.st_size == 0 and
             CAMPAIGN_CONSUMED_NAME.fullmatch(name) is not None,
             "campaign-consumed marker metadata is unsafe")
    return {
        "name": name,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o400,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _validate_lane_lock_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "preallocated lane-lock binding")
    _require(record["name"] == ".arming.lock" and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             record["mode"] == 0o600 and record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "preallocated lane-lock binding differs")
    return copy.deepcopy(record)


def _validate_lane_attempts_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "initial_link_count",
        "file_handle",
    )), "preallocated attempts-directory binding")
    _require(record["name"] == ATTEMPTS_DIRECTORY and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             record["mode"] == 0o500 and record["initial_link_count"] == 2,
             "preallocated attempts-directory binding differs")
    record["file_handle"] = prereg._validate_output_lane_file_handle(
        record["file_handle"])
    return copy.deepcopy(record)


def _validate_campaign_progress_marker_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "campaign journal-progress marker binding")
    _require(type(record["name"]) is str and
             JOURNAL_CHECKPOINT_NAME.fullmatch(record["name"]) is not None and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             record["mode"] == 0o400 and record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "campaign journal-progress marker binding differs")
    return copy.deepcopy(record)


def _validate_campaign_progress_directory_binding(
    value: Any,
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns", "file_handle", "markers",
    )), "campaign journal-checkpoint directory binding")
    _require(type(record["name"]) is str and
             CAMPAIGN_PROGRESS_DIRECTORY_NAME.fullmatch(record["name"])
                is not None and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             record["mode"] == 0o500 and record["link_count"] == 2 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0 and
             type(record["markers"]) is list and
             len(record["markers"]) ==
                2 * prereg.OUTPUT_LANE_CHILD_PROCESS_COUNT + 2,
             "campaign journal-checkpoint directory binding differs")
    record["file_handle"] = prereg._validate_output_lane_file_handle(
        record["file_handle"])
    markers = [
        _validate_campaign_progress_marker_binding(marker)
        for marker in record["markers"]
    ]
    _require(
        [marker["name"] for marker in markers] == [
            prereg._journal_checkpoint_name(sequence)
            for sequence in range(len(markers))
        ],
        "campaign journal checkpoints are not the exact ordered inventory")
    record["markers"] = markers
    return copy.deepcopy(record)


def _read_campaign_progress_marker(
    directory_fd: int, marker_fd: int, binding_value: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    binding = _validate_campaign_progress_marker_binding(binding_value)
    retained = os.fstat(marker_fd)
    pathname = os.stat(
        binding["name"], dir_fd=directory_fd, follow_symlinks=False)
    mode = stat.S_IMODE(retained.st_mode)
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        mode in (0o400, 0o600) and
        0 <= retained.st_size <= MAX_JSON_BYTES and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISREG(pathname.st_mode) and
        pathname.st_uid == binding["uid"] and pathname.st_nlink == 1 and
        stat.S_IMODE(pathname.st_mode) == mode and
        pathname.st_size == retained.st_size and
        pathname.st_mtime_ns == retained.st_mtime_ns and
        pathname.st_ctime_ns == retained.st_ctime_ns,
        f"{label} inode or monotone state changed")
    pristine = (
        mode == 0o400 and retained.st_size == 0 and
        retained.st_mtime_ns == binding["initial_mtime_ns"] and
        retained.st_ctime_ns == binding["initial_ctime_ns"])
    reached = not pristine
    _require(pristine != reached,
             f"{label} monotone state is ambiguous")
    _require(
        not reached or retained.st_ctime_ns != binding["initial_ctime_ns"],
        f"{label} changed without advancing its inode authority")
    _require(
        not (reached and mode == 0o400) or retained.st_size > 0,
        f"{label} is a sealed empty or erased checkpoint")
    data = b"" if retained.st_size == 0 else _read_exact_fd(
        marker_fd, MAX_JSON_BYTES, label)
    after = os.fstat(marker_fd)
    _require(
        (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
         after.st_nlink, after.st_size, after.st_mtime_ns,
         after.st_ctime_ns) ==
        (retained.st_dev, retained.st_ino, retained.st_mode, retained.st_uid,
         retained.st_nlink, retained.st_size, retained.st_mtime_ns,
         retained.st_ctime_ns),
        f"{label} changed while its payload was read")
    return {
        "binding": binding,
        "pristine": pristine,
        "reached": reached,
        "mode": mode,
        "ctime_ns": retained.st_ctime_ns,
        "data": data,
    }


def _flip_campaign_progress_marker(
    directory_fd: int, marker_fd: int, binding_value: Mapping[str, Any],
    payload: bytes, label: str,
) -> dict[str, Any]:
    _require(type(payload) is bytes and 0 < len(payload) <= MAX_JSON_BYTES,
             f"{label} checkpoint payload is invalid")
    binding = _validate_campaign_progress_marker_binding(binding_value)
    before = _read_campaign_progress_marker(
        directory_fd, marker_fd, binding, f"{label} before flip")
    if before["mode"] == 0o400 and before["reached"]:
        _require(before["data"] == payload,
                 f"{label} sealed checkpoint payload differs")
        return before
    if before["pristine"]:
        os.fchmod(marker_fd, 0o600)
        os.fsync(marker_fd)
    in_progress = _read_campaign_progress_marker(
        directory_fd, marker_fd, binding, f"{label} during flip")
    _require(
        in_progress["reached"] and in_progress["mode"] == 0o600 and
        payload.startswith(in_progress["data"]),
        f"{label} is not a recoverable exact payload prefix")
    writable_fd = -1
    try:
        writable_fd = os.open(
            binding["name"], os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        writable = _read_campaign_progress_marker(
            directory_fd, writable_fd, binding,
            f"{label} writable checkpoint")
        _require(
            writable["mode"] == 0o600 and
            writable["data"] == in_progress["data"],
            f"{label} changed while opening its writable capability")
        offset = len(writable["data"])
        while offset < len(payload):
            written = os.pwrite(writable_fd, payload[offset:], offset)
            _require(written > 0,
                     f"{label} checkpoint write made no progress")
            offset += written
        os.fsync(writable_fd)
        committed = _read_campaign_progress_marker(
            directory_fd, writable_fd, binding,
            f"{label} committed checkpoint")
        _require(committed["mode"] == 0o600 and
                 committed["data"] == payload,
                 f"{label} checkpoint payload did not commit exactly")
        os.fchmod(writable_fd, 0o400)
        os.fsync(writable_fd)
        after = _read_campaign_progress_marker(
            directory_fd, writable_fd, binding, f"{label} after flip")
        _require(after["reached"] and after["mode"] == 0o400 and
                 after["data"] == payload,
                 f"{label} did not reach its exact sealed checkpoint")
        return after
    finally:
        if writable_fd >= 0:
            os.close(writable_fd)


def _seal_campaign_progress_marker(
    directory_fd: int, marker_fd: int, binding_value: Mapping[str, Any],
    payload: bytes, label: str,
) -> dict[str, Any]:
    before = _read_campaign_progress_marker(
        directory_fd, marker_fd, binding_value, f"{label} before recovery")
    _require(before["reached"], f"{label} is not reached")
    return _flip_campaign_progress_marker(
        directory_fd, marker_fd, binding_value, payload, label)


def _open_campaign_progress_directory(
    lane: Path, lane_fd: int, binding_value: Mapping[str, Any], label: str,
) -> tuple[int, dict[str, Any]]:
    binding = _validate_campaign_progress_directory_binding(binding_value)
    descriptor = -1
    try:
        descriptor = os.open(
            binding["name"],
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        retained = os.fstat(descriptor)
        pathname = os.stat(
            binding["name"], dir_fd=lane_fd, follow_symlinks=False)
        retained_handle, retained_mount_id = \
            prereg._capture_output_lane_file_handle(
                lane_fd, os.fsencode(binding["name"]),
                f"retained {label}")
        pathname_handle, pathname_mount_id = \
            prereg._capture_output_lane_file_handle(
                -100, os.fsencode(lane / binding["name"]), f"path {label}")
        _require(
            stat.S_ISDIR(retained.st_mode) and
            retained.st_uid == binding["uid"] and
            stat.S_IMODE(retained.st_mode) == binding["mode"] and
            retained.st_nlink == binding["link_count"] and
            retained.st_mtime_ns == binding["initial_mtime_ns"] and
            retained.st_ctime_ns == binding["initial_ctime_ns"] and
            (retained.st_dev, retained.st_ino) ==
                (binding["device"], binding["inode"]) ==
                (pathname.st_dev, pathname.st_ino) and
            retained_mount_id == pathname_mount_id and
            contract.exact_json_equal(
                retained_handle, binding["file_handle"]) and
            contract.exact_json_equal(
                pathname_handle, binding["file_handle"]) and
            set(os.listdir(descriptor)) == {
                marker["name"] for marker in binding["markers"]},
            f"{label} differs from its preregistered topology")
        retained_descriptor = descriptor
        descriptor = -1
        return retained_descriptor, binding
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _validate_campaign_progress_directory_capability(
    lane: Path, lane_fd: int, directory_fd: int,
    binding_value: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    binding = _validate_campaign_progress_directory_binding(binding_value)
    retained = os.fstat(directory_fd)
    pathname = os.stat(
        binding["name"], dir_fd=lane_fd, follow_symlinks=False)
    retained_handle, retained_mount_id = \
        prereg._capture_output_lane_file_handle(
            lane_fd, os.fsencode(binding["name"]), f"retained {label}")
    pathname_handle, pathname_mount_id = \
        prereg._capture_output_lane_file_handle(
            -100, os.fsencode(lane / binding["name"]), f"path {label}")
    _require(
        stat.S_ISDIR(retained.st_mode) and
        stat.S_ISDIR(pathname.st_mode) and
        (retained.st_dev, retained.st_ino) ==
            (pathname.st_dev, pathname.st_ino) ==
            (binding["device"], binding["inode"]) and
        retained.st_uid == pathname.st_uid == binding["uid"] and
        stat.S_IMODE(retained.st_mode) ==
            stat.S_IMODE(pathname.st_mode) == binding["mode"] and
        retained.st_nlink == pathname.st_nlink == binding["link_count"] and
        retained.st_mtime_ns == pathname.st_mtime_ns ==
            binding["initial_mtime_ns"] and
        retained.st_ctime_ns == pathname.st_ctime_ns ==
            binding["initial_ctime_ns"] and
        retained_mount_id == pathname_mount_id and
        contract.exact_json_equal(retained_handle, binding["file_handle"]) and
        contract.exact_json_equal(pathname_handle, binding["file_handle"]),
        f"{label} retained capability differs")
    return binding


def _read_campaign_checkpoint_frontier(
    directory_fd: int, binding_value: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    binding = _validate_campaign_progress_directory_binding(binding_value)
    observed: list[dict[str, Any]] = []
    for marker_binding in binding["markers"]:
        marker_fd = os.open(
            marker_binding["name"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        try:
            marker = _read_campaign_progress_marker(
                directory_fd, marker_fd, marker_binding,
                f"{label} {marker_binding['name']}")
        finally:
            os.close(marker_fd)
        observed.append(marker)
    reached = [marker["reached"] for marker in observed]
    frontier = 0
    while frontier < len(reached) and reached[frontier]:
        frontier += 1
    _require(not any(reached[frontier:]),
             f"{label} checkpoints are not one exact prefix")
    _require(all(marker["mode"] == 0o400
                 for marker in observed[:max(0, frontier - 1)]),
             f"{label} has an impossible writable predecessor")
    return {"binding": binding, "markers": observed, "frontier": frontier}


def _flip_campaign_checkpoint(
    directory_fd: int, directory_binding: Mapping[str, Any], sequence: int,
    payload: bytes, label: str,
) -> dict[str, Any]:
    binding = _validate_campaign_progress_directory_binding(
        directory_binding)
    _require(type(sequence) is int and 0 <= sequence < len(binding["markers"]),
             f"{label} checkpoint sequence is invalid")
    if sequence > 0:
        prior_binding = binding["markers"][sequence - 1]
        prior_fd = os.open(
            prior_binding["name"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        try:
            prior = _read_campaign_progress_marker(
                directory_fd, prior_fd, prior_binding,
                f"{label} prior checkpoint")
            _require(prior["reached"] and prior["mode"] == 0o400,
                     f"{label} prior checkpoint is not sealed")
        finally:
            os.close(prior_fd)
    marker_binding = binding["markers"][sequence]
    marker_fd = os.open(
        marker_binding["name"],
        os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=directory_fd)
    try:
        result = _flip_campaign_progress_marker(
            directory_fd, marker_fd, marker_binding, payload, label)
    finally:
        os.close(marker_fd)
    return result


def _prearmed_history_marker_binding_record(
    metadata: os.stat_result, name: str,
) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o400 and
             metadata.st_size == 0 and
             PREARMED_HISTORY_NAME.fullmatch(name) is not None,
             "pre-ARMED history marker metadata is unsafe")
    return {
        "name": name,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o400,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _prearmed_boundary_marker_binding_record(
    metadata: os.stat_result, name: str,
) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o400 and
             metadata.st_size == 0 and
             PREARMED_BOUNDARY_NAME.fullmatch(name) is not None,
             "pre-ARMED boundary marker metadata is unsafe")
    return {
        "name": name,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o400,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _validate_prearmed_history_marker_binding(
    value: Any,
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "pre-ARMED history marker binding")
    _require(type(record["name"]) is str and
             PREARMED_HISTORY_NAME.fullmatch(record["name"]) is not None and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and
             record["uid"] == os.geteuid() and
             type(record["mode"]) is int and record["mode"] == 0o400 and
             type(record["link_count"]) is int and
             record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "pre-ARMED history marker binding is malformed")
    return copy.deepcopy(record)


def _validate_prearmed_boundary_marker_binding(
    value: Any,
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "pre-ARMED boundary marker binding")
    _require(type(record["name"]) is str and
             PREARMED_BOUNDARY_NAME.fullmatch(record["name"]) is not None and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and
             record["uid"] == os.geteuid() and
             type(record["mode"]) is int and record["mode"] == 0o400 and
             type(record["link_count"]) is int and
             record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "pre-ARMED boundary marker binding is malformed")
    return copy.deepcopy(record)


def _prearmed_boundary_bindings_for_acquisition(
    lane_binding: Mapping[str, Any], acquisition_index: int,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    _require(type(acquisition_index) is int and acquisition_index > 0,
             "pre-ARMED boundary acquisition index is invalid")
    offset = (acquisition_index - 1) * len(PREARMED_BOUNDARY_STATES)
    values = lane_binding["prearmed_boundary_markers"][
        offset:offset + len(PREARMED_BOUNDARY_STATES)]
    _require(len(values) == len(PREARMED_BOUNDARY_STATES),
             "pre-ARMED boundary slot is missing")
    bindings = tuple(
        _validate_prearmed_boundary_marker_binding(value)
        for value in values)
    expected_names = tuple(
        f"prearmed-{acquisition_index:04d}-{state.lower()}-reached.marker"
        for state in PREARMED_BOUNDARY_STATES)
    _require(tuple(binding["name"] for binding in bindings) == expected_names,
             "pre-ARMED boundary slot differs")
    return bindings  # type: ignore[return-value]


def _read_prearmed_boundary_marker(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    binding = _validate_prearmed_boundary_marker_binding(expected)
    try:
        retained = os.fstat(marker_fd)
        pathname = os.stat(
            binding["name"], dir_fd=lane_fd, follow_symlinks=False)
    except OSError as error:
        raise ArmingError(f"{label} changed") from error
    mode = stat.S_IMODE(retained.st_mode)
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        mode in (binding["mode"], 0o600) and retained.st_size == 0 and
        retained.st_mtime_ns == binding["initial_mtime_ns"] and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISREG(pathname.st_mode) and
        pathname.st_uid == binding["uid"] and pathname.st_nlink == 1 and
        stat.S_IMODE(pathname.st_mode) == mode and pathname.st_size == 0 and
        pathname.st_mtime_ns == retained.st_mtime_ns and
        pathname.st_ctime_ns == retained.st_ctime_ns,
        f"{label} inode or monotone state changed")
    pristine = (
        mode == binding["mode"] and
        retained.st_ctime_ns == binding["initial_ctime_ns"])
    reached = retained.st_ctime_ns != binding["initial_ctime_ns"]
    _require(pristine != reached,
             f"{label} monotone state is ambiguous")
    match = PREARMED_BOUNDARY_NAME.fullmatch(binding["name"])
    _require(match is not None, f"{label} name differs")
    return {
        "binding": binding,
        "acquisition_index": int(match.group(1)),
        "state": match.group(2).upper(),
        "reached": reached,
        "mode": mode,
        "mtime_ns": retained.st_mtime_ns,
        "ctime_ns": retained.st_ctime_ns,
    }


def _read_prearmed_boundary_frontier(
    lane_fd: int, lane_binding: Mapping[str, Any], acquisition_index: int,
) -> dict[str, Any]:
    observed: list[dict[str, Any]] = []
    for binding in _prearmed_boundary_bindings_for_acquisition(
            lane_binding, acquisition_index):
        descriptor = -1
        try:
            descriptor = os.open(
                binding["name"],
                os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
            marker = _read_prearmed_boundary_marker(
                lane_fd, descriptor, binding,
                f"pre-ARMED boundary {binding['name']}")
        except OSError as error:
            raise ArmingError(
                f"pre-ARMED boundary {binding['name']} cannot be opened") \
                from error
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        observed.append(marker)
    reached = [marker["reached"] for marker in observed]
    frontier = 0
    while frontier < len(reached) and reached[frontier]:
        frontier += 1
    _require(not any(reached[frontier:]),
             "pre-ARMED boundary markers are not one exact prefix")
    _require(all(
        marker["mode"] == 0o400
        for marker in observed[:max(0, frontier - 1)]),
        "pre-ARMED boundary prefix has an impossible writable predecessor")
    _require(
        tuple(marker["state"] for marker in observed) ==
            PREARMED_BOUNDARY_STATES,
        "pre-ARMED boundary marker states differ")
    return {
        "markers": observed,
        "frontier": frontier,
        "states": ["INIT", *PREARMED_BOUNDARY_STATES[:frontier]],
    }


def _flip_prearmed_boundary_marker(
    lane_fd: int, binding_value: Mapping[str, Any], expected_state: str,
) -> dict[str, Any]:
    binding = _validate_prearmed_boundary_marker_binding(binding_value)
    match = PREARMED_BOUNDARY_NAME.fullmatch(binding["name"])
    _require(match is not None and match.group(2).upper() == expected_state,
             "pre-ARMED boundary transition differs")
    descriptor = os.open(
        binding["name"], os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=lane_fd)
    try:
        before = _read_prearmed_boundary_marker(
            lane_fd, descriptor, binding,
            f"pristine {expected_state} boundary marker")
        _require(not before["reached"] and before["mode"] == 0o400,
                 f"{expected_state} boundary was already reached")
        os.fchmod(descriptor, 0o600)
        os.fsync(descriptor)
        in_progress = _read_prearmed_boundary_marker(
            lane_fd, descriptor, binding,
            f"in-progress {expected_state} boundary marker")
        _require(in_progress["reached"] and in_progress["mode"] == 0o600,
                 f"{expected_state} boundary did not become durable")
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        after = _read_prearmed_boundary_marker(
            lane_fd, descriptor, binding,
            f"sealed {expected_state} boundary marker")
        _require(after["reached"] and after["mode"] == 0o400,
                 f"{expected_state} boundary did not seal")
        return after
    finally:
        os.close(descriptor)


def _validate_recoverable_prearmed_history_marker(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    binding = _validate_prearmed_history_marker_binding(expected)
    try:
        retained = os.fstat(marker_fd)
        pathname = os.stat(
            binding["name"], dir_fd=lane_fd, follow_symlinks=False)
    except OSError as error:
        raise ArmingError(f"{label} changed") from error
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        stat.S_IMODE(retained.st_mode) in (binding["mode"], 0o600) and
        0 <= retained.st_size <= MAX_PREARMED_HISTORY_BYTES and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISREG(pathname.st_mode) and
        pathname.st_uid == binding["uid"] and pathname.st_nlink == 1 and
        stat.S_IMODE(pathname.st_mode) == stat.S_IMODE(retained.st_mode) and
        pathname.st_size == retained.st_size and
        (retained.st_size != 0 or
         stat.S_IMODE(retained.st_mode) == 0o600 or
         (retained.st_mtime_ns == binding["initial_mtime_ns"] and
          retained.st_ctime_ns == binding["initial_ctime_ns"])),
        f"{label} inode or monotone state changed")
    return retained


def _validate_campaign_marker_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
        "initial_mtime_ns", "initial_ctime_ns",
    )), "campaign-consumed marker binding")
    _require(type(record["name"]) is str and
             CAMPAIGN_CONSUMED_NAME.fullmatch(record["name"]) is not None and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             type(record["mode"]) is int and record["mode"] == 0o400 and
             type(record["link_count"]) is int and
             record["link_count"] == 1 and
             type(record["initial_mtime_ns"]) is int and
             record["initial_mtime_ns"] > 0 and
             type(record["initial_ctime_ns"]) is int and
             record["initial_ctime_ns"] > 0,
             "campaign-consumed marker binding is malformed")
    return copy.deepcopy(record)


def _validate_recoverable_published_campaign_marker(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    binding = _validate_campaign_marker_binding(expected)
    try:
        retained = os.fstat(marker_fd)
        pathname = os.stat(
            binding["name"], dir_fd=lane_fd,
            follow_symlinks=False)
    except OSError as error:
        raise ArmingError(f"{label} changed") from error
    _require(
        stat.S_ISREG(retained.st_mode) and
        retained.st_uid == binding["uid"] and retained.st_nlink == 1 and
        stat.S_IMODE(retained.st_mode) in (binding["mode"], 0o600) and
        0 <= retained.st_size <= MAX_CAMPAIGN_TRANSCRIPT_BYTES and
        (retained.st_dev, retained.st_ino) ==
            (binding["device"], binding["inode"]) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISREG(pathname.st_mode) and pathname.st_uid == binding["uid"] and
        pathname.st_nlink == 1 and
        stat.S_IMODE(pathname.st_mode) == stat.S_IMODE(retained.st_mode) and
        pathname.st_size == retained.st_size and
        (retained.st_size != 0 or
         (retained.st_mtime_ns == binding["initial_mtime_ns"] and
          retained.st_ctime_ns >= binding["initial_ctime_ns"])),
        f"{label} inode or monotone state changed")
    return retained


def _validate_campaign_marker_capability(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    binding = _validate_campaign_marker_binding(expected)
    retained = _validate_recoverable_published_campaign_marker(
        lane_fd, marker_fd, binding, label)
    _require(stat.S_IMODE(retained.st_mode) == binding["mode"],
             f"{label} inode or monotone state changed")
    return retained


def _recover_published_campaign_marker_mode(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> os.stat_result:
    binding = _validate_campaign_marker_binding(expected)
    metadata = _validate_recoverable_published_campaign_marker(
        lane_fd, marker_fd, binding, label)
    if stat.S_IMODE(metadata.st_mode) == 0o600:
        os.fchmod(marker_fd, binding["mode"])
        os.fsync(marker_fd)
    return _validate_campaign_marker_capability(
        lane_fd, marker_fd, binding, label)


def _validate_unallocated_campaign_marker(
    lane_fd: int, marker_fd: int, expected: Mapping[str, Any], label: str,
) -> None:
    binding = _validate_campaign_marker_binding(expected)
    metadata = _validate_campaign_marker_capability(
        lane_fd, marker_fd, binding, label)
    _require(metadata.st_size == 0 and
             metadata.st_mtime_ns == binding["initial_mtime_ns"] and
             metadata.st_ctime_ns == binding["initial_ctime_ns"],
             f"{label} was allocated without a published attempt")


def _prearmed_state_entry_record(
    *, acquisition_index: int, sequence: int, from_state: str, to_state: str,
    prior_entry_sha256: str | None, preregistration_sha256: str,
    plan_sha256: str, lane_binding_sha256: str,
) -> dict[str, Any]:
    _require(type(acquisition_index) is int and
             1 <= acquisition_index <= 3000 and
             type(sequence) is int and 1 <= sequence <= 7 and
             PREARMED_STATE_TRANSITIONS.get(from_state) == to_state and
             (prior_entry_sha256 is None if sequence == 1 else
              re.fullmatch(r"[0-9a-f]{64}", prior_entry_sha256 or "")
                  is not None),
             "pre-ARMED state transition differs")
    return {
        "schema": PREARMED_STATE_ENTRY_SCHEMA,
        "generation": GENERATION,
        "acquisition_index": acquisition_index,
        "sequence": sequence,
        "from_state": from_state,
        "to_state": to_state,
        "prior_entry_sha256": prior_entry_sha256,
        "preregistration_sha256": _hex_digest(
            preregistration_sha256, "pre-ARMED preregistration hash"),
        "plan_sha256": _hex_digest(
            plan_sha256, "pre-ARMED plan hash"),
        "lane_binding_sha256": _hex_digest(
            lane_binding_sha256, "pre-ARMED lane-binding hash"),
    }


def _read_prearmed_history(
    lane_fd: int, marker_fd: int, binding_value: Mapping[str, Any], *,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    lane_binding: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    binding = _validate_prearmed_history_marker_binding(binding_value)
    before = _validate_recoverable_prearmed_history_marker(
        lane_fd, marker_fd, binding, label)
    data = b"" if before.st_size == 0 else _read_exact_fd(
        marker_fd, MAX_PREARMED_HISTORY_BYTES, label)
    after = _validate_recoverable_prearmed_history_marker(
        lane_fd, marker_fd, binding, label)
    _require(
        len(data) == before.st_size == after.st_size and
        (before.st_mode, before.st_mtime_ns, before.st_ctime_ns) ==
            (after.st_mode, after.st_mtime_ns, after.st_ctime_ns),
        f"{label} changed while read")
    mode = stat.S_IMODE(before.st_mode)
    blocks = data.splitlines(keepends=True)
    partial = b""
    if blocks and not blocks[-1].endswith(b"\n"):
        partial = blocks.pop()
    _require(not partial or mode == 0o600,
             f"{label} sealed an incomplete state entry")
    complete_data = b"".join(blocks)
    states = ["INIT"]
    prior: str | None = None
    preregistration_sha256 = contract.canonical_sha256(registration)
    plan_sha256 = contract.canonical_sha256(logical_plan)
    lane_binding_sha256 = contract.canonical_sha256(lane_binding)
    match = PREARMED_HISTORY_NAME.fullmatch(binding["name"])
    _require(match is not None, f"{label} name differs")
    acquisition_index = int(match.group(1))
    for sequence, block in enumerate(blocks, 1):
        value = _exact_mapping(
            _canonical_json_document(
                block, f"{label} state entry {sequence}"),
            frozenset((
                "schema", "generation", "acquisition_index", "sequence",
                "from_state", "to_state", "prior_entry_sha256",
                "preregistration_sha256", "plan_sha256",
                "lane_binding_sha256",
            )), f"{label} state entry {sequence}")
        expected = _prearmed_state_entry_record(
            acquisition_index=acquisition_index, sequence=sequence,
            from_state=states[-1], to_state=value["to_state"],
            prior_entry_sha256=prior,
            preregistration_sha256=preregistration_sha256,
            plan_sha256=plan_sha256,
            lane_binding_sha256=lane_binding_sha256)
        _require(contract.exact_json_equal(value, expected) and
                 block == contract.canonical_json_bytes(expected),
                 f"{label} state entry {sequence} differs")
        states.append(expected["to_state"])
        prior = _sha256_bytes(block)
    if partial:
        next_state = PREARMED_STATE_TRANSITIONS.get(states[-1])
        _require(next_state is not None,
                 f"{label} extends terminal pre-ARMED state")
        expected_partial_block = contract.canonical_json_bytes(
            _prearmed_state_entry_record(
                acquisition_index=acquisition_index,
                sequence=len(states), from_state=states[-1],
                to_state=next_state, prior_entry_sha256=prior,
                preregistration_sha256=preregistration_sha256,
                plan_sha256=plan_sha256,
                lane_binding_sha256=lane_binding_sha256))
        _require(expected_partial_block.startswith(partial),
                 f"{label} has a non-canonical partial state entry")
    pristine_unused = (
        mode == 0o400 and not data and
        before.st_mtime_ns == binding["initial_mtime_ns"] and
        before.st_ctime_ns == binding["initial_ctime_ns"])
    used = not pristine_unused
    _require(used or pristine_unused, f"{label} unused state differs")
    return {
        "binding": binding,
        "acquisition_index": acquisition_index,
        "states": states,
        "used": used,
        "mode": mode,
        "data": data,
        "complete_data": complete_data,
        "complete_size": len(complete_data),
        "partial": partial,
        "prior_entry_sha256": prior,
    }


def _open_prearmed_history_for_allocation(
    lane_fd: int, binding_value: Mapping[str, Any], *,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    lane_binding: Mapping[str, Any],
) -> int:
    binding = _validate_prearmed_history_marker_binding(binding_value)
    probe = os.open(
        binding["name"], os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=lane_fd)
    descriptor = -1
    try:
        history = _read_prearmed_history(
            lane_fd, probe, binding, registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding,
            label="unused pre-ARMED history slot")
        _require(not history["used"],
                 "pre-ARMED history slot was already allocated")
        os.fchmod(probe, 0o600)
        os.fsync(probe)
        descriptor = os.open(
            binding["name"], os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        allocated = _read_prearmed_history(
            lane_fd, descriptor, binding, registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding,
            label="allocated pre-ARMED history slot")
        _require(allocated["used"] and allocated["mode"] == 0o600 and
                 allocated["data"] == b"",
                 "pre-ARMED history allocation did not reach its fixed point")
        retained = descriptor
        descriptor = -1
        return retained
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(probe)


def _append_prearmed_state(
    lane_fd: int, marker_fd: int, binding_value: Mapping[str, Any],
    to_state: str, *, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], lane_binding: Mapping[str, Any],
) -> dict[str, Any]:
    binding = _validate_prearmed_history_marker_binding(binding_value)
    before = _read_prearmed_history(
        lane_fd, marker_fd, binding, registration=registration,
        logical_plan=logical_plan, lane_binding=lane_binding,
        label="active pre-ARMED history")
    _require(before["used"] and before["mode"] == 0o600 and
             not before["partial"],
             "active pre-ARMED history is not appendable")
    entry = _prearmed_state_entry_record(
        acquisition_index=before["acquisition_index"],
        sequence=len(before["states"]), from_state=before["states"][-1],
        to_state=to_state,
        prior_entry_sha256=before["prior_entry_sha256"],
        preregistration_sha256=contract.canonical_sha256(registration),
        plan_sha256=contract.canonical_sha256(logical_plan),
        lane_binding_sha256=contract.canonical_sha256(lane_binding))
    block = contract.canonical_json_bytes(entry)
    offset = len(before["data"])
    while offset < len(before["data"]) + len(block):
        count = os.pwrite(
            marker_fd, block[offset - len(before["data"]):], offset)
        _require(count > 0, "pre-ARMED history append made no progress")
        offset += count
    os.fsync(marker_fd)
    after = _read_prearmed_history(
        lane_fd, marker_fd, binding, registration=registration,
        logical_plan=logical_plan, lane_binding=lane_binding,
        label="appended pre-ARMED history")
    _require(after["data"] == before["data"] + block and
             after["states"][-1] == to_state and not after["partial"],
             "pre-ARMED history append differs")
    return after


def _seal_prearmed_history(
    lane_fd: int, marker_fd: int, binding_value: Mapping[str, Any], *,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    lane_binding: Mapping[str, Any],
) -> dict[str, Any]:
    before = _read_prearmed_history(
        lane_fd, marker_fd, binding_value, registration=registration,
        logical_plan=logical_plan, lane_binding=lane_binding,
        label="pre-ARMED history before sealing")
    _require(before["used"] and before["mode"] == 0o600 and
             not before["partial"],
             "pre-ARMED history cannot be sealed")
    os.fchmod(marker_fd, 0o400)
    os.fsync(marker_fd)
    after = _read_prearmed_history(
        lane_fd, marker_fd, binding_value, registration=registration,
        logical_plan=logical_plan, lane_binding=lane_binding,
        label="sealed pre-ARMED history")
    _require(after["mode"] == 0o400 and after["data"] == before["data"],
             "sealed pre-ARMED history differs")
    return after


def _prearmed_history_reference(
    history: Mapping[str, Any],
) -> dict[str, Any]:
    _require(history["used"] and not history["partial"] and
             history["states"][-1] == "PRESAMPLING" and
             history["complete_data"] == history["data"],
             "evidence history did not reach exact PRESAMPLING")
    return {
        "marker": copy.deepcopy(history["binding"]),
        "acquisition_index": history["acquisition_index"],
        "history_sha256": _sha256_bytes(history["complete_data"]),
        "history_size": history["complete_size"],
        "highest_state": "PRESAMPLING",
    }


def _validate_prearmed_history_reference(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "marker", "acquisition_index", "history_sha256", "history_size",
        "highest_state",
    )), "pre-ARMED evidence history reference")
    marker = _validate_prearmed_history_marker_binding(record["marker"])
    match = PREARMED_HISTORY_NAME.fullmatch(marker["name"])
    _require(match is not None and
             type(record["acquisition_index"]) is int and
             record["acquisition_index"] == int(match.group(1)) and
             type(record["history_size"]) is int and
             0 < record["history_size"] <= MAX_PREARMED_HISTORY_BYTES and
             record["highest_state"] == "PRESAMPLING",
             "pre-ARMED evidence history reference differs")
    return {
        "marker": marker,
        "acquisition_index": record["acquisition_index"],
        "history_sha256": _hex_digest(
            record["history_sha256"], "pre-ARMED history hash"),
        "history_size": record["history_size"],
        "highest_state": "PRESAMPLING",
    }


def _budget_commit_record(
    history_reference: Mapping[str, Any],
    prospective_ledger: Mapping[str, Any],
) -> dict[str, Any]:
    reference = _validate_prearmed_history_reference(history_reference)
    return {
        "history": reference,
        "prospective_ledger": copy.deepcopy(prospective_ledger),
        "prospective_ledger_sha256": contract.canonical_sha256(
            prospective_ledger),
    }


def _validate_budget_commit(
    value: Any, registration: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "history", "prospective_ledger", "prospective_ledger_sha256",
    )), "generation-3 budget commit")
    history = _validate_prearmed_history_reference(record["history"])
    ledger = plan_contract.validate_budget_ledger(
        record["prospective_ledger"], registration)
    _require(record["prospective_ledger_sha256"] ==
             contract.canonical_sha256(ledger),
             "generation-3 prospective budget-ledger hash differs")
    return _budget_commit_record(history, ledger)


def _open_campaign_marker_for_append(
    lane_fd: int, binding: Mapping[str, Any], *,
    associated_with_published_attempt: bool = False,
) -> int:
    """Open one unused prebound slot writable, then immediately reseal it."""
    expected = _validate_campaign_marker_binding(binding)
    probe = -1
    descriptor = -1
    try:
        probe = os.open(
            expected["name"], os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        if associated_with_published_attempt:
            metadata = _recover_published_campaign_marker_mode(
                lane_fd, probe, expected,
                "published attempt campaign transcript slot")
        else:
            _validate_unallocated_campaign_marker(
                lane_fd, probe, expected, "unused campaign transcript slot")
            metadata = os.fstat(probe)
        retained_size = metadata.st_size
        os.fchmod(probe, 0o600)
        os.fsync(probe)
        descriptor = os.open(
            expected["name"], os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        writable = os.fstat(descriptor)
        _require((writable.st_dev, writable.st_ino) ==
                 (expected["device"], expected["inode"]) and
                 writable.st_size == retained_size,
                 "campaign transcript slot changed while opening")
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        _validate_campaign_marker_capability(
            lane_fd, descriptor, expected,
            "opened campaign transcript slot")
        os.fchmod(probe, 0o400)
        os.fsync(probe)
        os.close(probe)
        probe = -1
        retained = descriptor
        descriptor = -1
        return retained
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors: list[BaseException] = []
        for current in (descriptor, probe):
            if current < 0:
                continue
            try:
                os.fchmod(current, 0o400)
            except BaseException as error:
                cleanup_errors.append(error)
            try:
                os.fsync(current)
            except BaseException as error:
                cleanup_errors.append(error)
            try:
                os.close(current)
            except BaseException as error:
                cleanup_errors.append(error)
        if cleanup_errors and active_error is None:
            raise ArmingError(
                "campaign transcript writable capability cleanup failed") \
                from cleanup_errors[0]


def _campaign_transcript_allocation_record(
    armed: Mapping[str, Any],
) -> dict[str, Any]:
    _require(type(armed.get("evidence_attempt")) is int and
             armed["evidence_attempt"] > 0 and
             type(armed.get("authority_bundle_sha256")) is str and
             type(armed.get("attempt_manifest_sha256")) is str and
             type(armed.get("lane_binding_sha256")) is str and
             type(armed.get("host_instance_sha256")) is str and
             type(armed.get("selected_pair")) is dict,
             "campaign transcript allocation binding is malformed")
    return {
        "schema": CAMPAIGN_TRANSCRIPT_ALLOCATION_SCHEMA,
        "sequence": 0,
        "evidence_attempt": armed["evidence_attempt"],
        "armed_sha256": contract.canonical_sha256(armed),
        "authority_bundle_sha256": _hex_digest(
            armed["authority_bundle_sha256"],
            "transcript authority-bundle hash"),
        "attempt_manifest_sha256": _hex_digest(
            armed["attempt_manifest_sha256"],
            "transcript attempt-manifest hash"),
        "lane_binding_sha256": _hex_digest(
            armed["lane_binding_sha256"], "transcript lane-binding hash"),
        "host_instance_sha256": _hex_digest(
            armed["host_instance_sha256"], "transcript host-instance hash"),
        "selected_pair_sha256": contract.canonical_sha256(
            armed["selected_pair"]),
    }


def _validate_campaign_transcript_allocation(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "schema", "sequence", "evidence_attempt", "armed_sha256",
        "authority_bundle_sha256", "attempt_manifest_sha256",
        "lane_binding_sha256", "host_instance_sha256",
        "selected_pair_sha256",
    )), "campaign transcript allocation")
    _require(record["schema"] == CAMPAIGN_TRANSCRIPT_ALLOCATION_SCHEMA and
             record["sequence"] == 0 and
             type(record["evidence_attempt"]) is int and
             record["evidence_attempt"] > 0,
             "campaign transcript allocation metadata differs")
    for name in (
            "armed_sha256", "authority_bundle_sha256",
            "attempt_manifest_sha256", "lane_binding_sha256",
            "host_instance_sha256", "selected_pair_sha256"):
        _hex_digest(record[name], f"campaign transcript allocation {name}")
    return copy.deepcopy(record)


def _campaign_transcript_binding_sha256(
    armed: Mapping[str, Any],
) -> str:
    return contract.canonical_sha256(
        _campaign_transcript_allocation_record(armed))


def _campaign_transcript_entry_record(
    *, sequence: int, name: str, data_sha256: str,
    prior_entry_sha256: str | None, attempt_binding_sha256: str,
) -> dict[str, Any]:
    _require(type(sequence) is int and sequence >= 1 and
             (name == "complete.json" or
              JOURNAL_INTENT_NAME.fullmatch(name) is not None or
              JOURNAL_RESULT_NAME.fullmatch(name) is not None),
             "campaign transcript entry position is malformed")
    return {
        "schema": CAMPAIGN_TRANSCRIPT_ENTRY_SCHEMA,
        "sequence": sequence,
        "name": name,
        "data_sha256": _hex_digest(
            data_sha256, "campaign transcript data hash"),
        "prior_entry_sha256": (
            None if prior_entry_sha256 is None else
            _hex_digest(prior_entry_sha256,
                        "campaign transcript prior-entry hash")),
        "attempt_binding_sha256": _hex_digest(
            attempt_binding_sha256,
            "campaign transcript attempt-binding hash"),
    }


def _expected_campaign_transcript_blocks(
    armed: Mapping[str, Any],
    journal_records: Sequence[tuple[str, str]],
) -> list[bytes]:
    allocation = _campaign_transcript_allocation_record(armed)
    blocks = [contract.canonical_json_bytes(allocation)]
    prior = contract.canonical_sha256(allocation)
    attempt_binding_sha256 = prior
    for sequence, (name, data_sha256) in enumerate(journal_records, 1):
        entry = _campaign_transcript_entry_record(
            sequence=sequence, name=name, data_sha256=data_sha256,
            prior_entry_sha256=prior,
            attempt_binding_sha256=attempt_binding_sha256)
        blocks.append(contract.canonical_json_bytes(entry))
        prior = contract.canonical_sha256(entry)
    _require(sum(map(len, blocks)) <= MAX_CAMPAIGN_TRANSCRIPT_BYTES,
             "expected campaign transcript exceeds its byte bound")
    return blocks


def _expected_campaign_transcript_bytes(
    armed: Mapping[str, Any],
    journal_records: Sequence[tuple[str, str]],
) -> bytes:
    result = b"".join(_expected_campaign_transcript_blocks(
        armed, journal_records))
    _require(len(result) <= MAX_CAMPAIGN_TRANSCRIPT_BYTES,
             "expected campaign transcript exceeds its byte bound")
    return result


def _read_campaign_transcript_bytes(
    lane_fd: int, marker_fd: int, binding: Mapping[str, Any], label: str, *,
    recoverable: bool = False,
) -> bytes:
    validator = (
        _validate_recoverable_published_campaign_marker
        if recoverable else _validate_campaign_marker_capability)
    before = validator(lane_fd, marker_fd, binding, label)
    data = b"" if before.st_size == 0 else _read_exact_fd(
        marker_fd, MAX_CAMPAIGN_TRANSCRIPT_BYTES, label)
    after = validator(lane_fd, marker_fd, binding, label)
    _require(
        len(data) == before.st_size == after.st_size and
        (before.st_mode, before.st_mtime_ns, before.st_ctime_ns) ==
            (after.st_mode, after.st_mtime_ns, after.st_ctime_ns),
        f"{label} changed while it was read")
    return data


def _parse_campaign_transcript_bytes(
    data: bytes, label: str,
) -> dict[str, Any]:
    if not data:
        return {"allocation": None, "journal_entries": []}
    blocks = data.splitlines(keepends=True)
    _require(bool(blocks) and b"".join(blocks) == data and
             all(block.endswith(b"\n") for block in blocks),
             f"{label} framing differs")
    allocation = _validate_campaign_transcript_allocation(
        _canonical_json_document(blocks[0], f"{label} allocation"))
    entries: list[dict[str, Any]] = []
    prior = contract.canonical_sha256(allocation)
    for sequence, block in enumerate(blocks[1:], 1):
        value = _exact_mapping(
            _canonical_json_document(
                block, f"{label} entry {sequence}"),
            frozenset((
                "schema", "sequence", "name", "data_sha256",
                "prior_entry_sha256", "attempt_binding_sha256",
            )), f"{label} entry {sequence}")
        _require(value["schema"] == CAMPAIGN_TRANSCRIPT_ENTRY_SCHEMA,
                 f"{label} entry schema differs")
        expected = _campaign_transcript_entry_record(
            sequence=sequence, name=value["name"],
            data_sha256=value["data_sha256"], prior_entry_sha256=prior,
            attempt_binding_sha256=value["attempt_binding_sha256"])
        _require(contract.exact_json_equal(value, expected),
                 f"{label} entry chain differs")
        entries.append(expected)
        prior = contract.canonical_sha256(expected)
    return {"allocation": allocation, "journal_entries": entries}


def _read_campaign_transcript(
    lane_fd: int, marker_fd: int, binding: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    return _parse_campaign_transcript_bytes(
        _read_campaign_transcript_bytes(
            lane_fd, marker_fd, binding, label),
        label)


def _validate_campaign_transcript(
    lane_fd: int, marker_fd: int, binding: Mapping[str, Any], *,
    journal_records: Sequence[tuple[str, str]],
    armed: Mapping[str, Any], label: str,
    allow_unallocated: bool = False,
) -> dict[str, Any]:
    observed = _read_campaign_transcript(
        lane_fd, marker_fd, binding, label)
    if observed["allocation"] is None:
        _require(allow_unallocated and not journal_records,
                 f"{label} lacks its durable ARMED allocation")
        return observed
    allocation = _campaign_transcript_allocation_record(armed)
    _require(contract.exact_json_equal(
        observed["allocation"], allocation),
        f"{label} allocation does not bind its published attempt")
    expected: list[dict[str, Any]] = []
    prior = contract.canonical_sha256(allocation)
    attempt_binding_sha256 = prior
    _require(len(observed["journal_entries"]) == len(journal_records),
             f"{label} journal length differs")
    for sequence, ((name, data_sha256), observed_entry) in enumerate(
            zip(journal_records, observed["journal_entries"]), 1):
        entry = _campaign_transcript_entry_record(
            sequence=sequence, name=name, data_sha256=data_sha256,
            prior_entry_sha256=prior,
            attempt_binding_sha256=attempt_binding_sha256)
        expected.append(entry)
        prior = contract.canonical_sha256(entry)
    _require(contract.exact_json_equal(
        observed["journal_entries"], expected),
             f"{label} does not exactly authenticate the journal")
    return {"allocation": allocation, "journal_entries": expected}


def _validate_campaign_checkpoint_payloads(
    directory_fd: int, directory_binding: Mapping[str, Any], *,
    armed: Mapping[str, Any],
    journal_records: Sequence[tuple[str, str]], label: str,
    allow_recoverable_tail: bool = False,
) -> dict[str, Any]:
    checkpoint = _read_campaign_checkpoint_frontier(
        directory_fd, directory_binding, label)
    blocks = _expected_campaign_transcript_blocks(armed, journal_records)
    frontier = checkpoint["frontier"]
    if allow_recoverable_tail:
        _require(
            frontier in (len(blocks), len(blocks) - 1) and frontier >= 1,
            f"{label} is not the exact journal/checkpoint crash frontier")
    else:
        _require(frontier == len(blocks),
                 f"{label} checkpoint length differs from the journal")
    journal_ahead = frontier == len(blocks) - 1
    for sequence, marker in enumerate(checkpoint["markers"][:frontier]):
        expected = blocks[sequence]
        if marker["mode"] == 0o400:
            _require(marker["data"] == expected,
                     f"{label} checkpoint {sequence} payload differs")
        else:
            _require(
                allow_recoverable_tail and not journal_ahead and
                sequence == frontier - 1 and
                expected.startswith(marker["data"]),
                f"{label} checkpoint {sequence} is not an exact recoverable "
                "payload prefix")
    return {
        "binding": checkpoint["binding"],
        "markers": checkpoint["markers"],
        "frontier": frontier,
        "blocks": blocks,
        "payload_prefix": b"".join(
            marker["data"] for marker in checkpoint["markers"][:frontier]),
    }


def _validate_campaign_checkpoint_payload_at(
    directory_fd: int, directory_binding: Mapping[str, Any], sequence: int,
    payload: bytes | None, label: str,
) -> dict[str, Any]:
    binding = _validate_campaign_progress_directory_binding(
        directory_binding)
    _require(type(sequence) is int and 0 <= sequence < len(binding["markers"]),
             f"{label} checkpoint sequence is invalid")
    marker_binding = binding["markers"][sequence]
    marker_fd = os.open(
        marker_binding["name"],
        os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=directory_fd)
    try:
        marker = _read_campaign_progress_marker(
            directory_fd, marker_fd, marker_binding, label)
    finally:
        os.close(marker_fd)
    if payload is None:
        _require(marker["pristine"] and marker["mode"] == 0o400 and
                 marker["data"] == b"",
                 f"{label} future checkpoint is not pristine")
    else:
        _require(type(payload) is bytes and payload and
                 marker["reached"] and marker["mode"] == 0o400 and
                 marker["data"] == payload,
                 f"{label} committed checkpoint payload differs")
    return marker


def _validate_live_campaign_checkpoint_tail(
    directory_fd: int, directory_binding: Mapping[str, Any], *,
    armed: Mapping[str, Any],
    journal_records: Sequence[tuple[str, str]], label: str,
) -> None:
    binding = _validate_campaign_progress_directory_binding(
        directory_binding)
    blocks = _expected_campaign_transcript_blocks(armed, journal_records)
    # Revalidate every reached payload before each child boundary.  The
    # current checkpoint's hash chain authenticates journal contents, but an
    # older preregistered checkpoint inode can still be independently erased
    # or rewritten in place without changing the directory metadata.
    for sequence, block in enumerate(blocks):
        _validate_campaign_checkpoint_payload_at(
            directory_fd, binding, sequence, block,
            f"{label} reached {sequence}")
    if len(blocks) < len(binding["markers"]):
        _validate_campaign_checkpoint_payload_at(
            directory_fd, binding, len(blocks), None,
            f"{label} next")


def _reconcile_campaign_transcript_prefix(
    lane_fd: int, binding: Mapping[str, Any], *,
    armed: Mapping[str, Any],
    journal_records: Sequence[tuple[str, str]], label: str,
    progress_directory_fd: int,
    progress_directory_binding: Mapping[str, Any],
) -> None:
    checkpoint = _validate_campaign_checkpoint_payloads(
        progress_directory_fd, progress_directory_binding,
        armed=armed, journal_records=journal_records,
        label=f"{label} checkpoints", allow_recoverable_tail=True)
    probe = os.open(
        binding["name"], os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=lane_fd)
    try:
        _recover_published_campaign_marker_mode(
            lane_fd, probe, binding, label)
        observed_data = _read_campaign_transcript_bytes(
            lane_fd, probe, binding, label)
    finally:
        os.close(probe)
    checkpoint_prefix = checkpoint["payload_prefix"]
    _require(checkpoint_prefix.startswith(observed_data),
             f"{label} advanced beyond or differs from its checkpoints")

    if checkpoint["frontier"] < len(checkpoint["blocks"]):
        sequence = checkpoint["frontier"]
        _flip_campaign_checkpoint(
            progress_directory_fd, progress_directory_binding, sequence,
            checkpoint["blocks"][sequence],
            f"{label} recovered durable journal record")
    else:
        sequence = checkpoint["frontier"] - 1
        marker = checkpoint["markers"][sequence]
        if marker["mode"] == 0o600:
            marker_binding = checkpoint["binding"]["markers"][sequence]
            marker_fd = os.open(
                marker_binding["name"],
                os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=progress_directory_fd)
            try:
                _seal_campaign_progress_marker(
                    progress_directory_fd, marker_fd, marker_binding,
                    checkpoint["blocks"][sequence],
                    f"{label} recovered partial checkpoint")
            finally:
                os.close(marker_fd)

    expected_data = b"".join(checkpoint["blocks"])
    _require(expected_data.startswith(observed_data),
             f"{label} is not the exact checkpointed transcript prefix")
    descriptor = _open_campaign_marker_for_append(
        lane_fd, binding, associated_with_published_attempt=True)
    try:
        current_data = _read_campaign_transcript_bytes(
            lane_fd, descriptor, binding, label)
        _require(current_data == observed_data,
                 f"{label} changed before transcript recovery")
        offset = len(current_data)
        while offset < len(expected_data):
            count = os.pwrite(
                descriptor, expected_data[offset:], offset)
            _require(count > 0,
                     f"{label} recovery write made no progress")
            offset += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    descriptor = os.open(
        binding["name"], os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=lane_fd)
    try:
        _validate_campaign_transcript(
            lane_fd, descriptor, binding, journal_records=journal_records,
            armed=armed, label=label)
    finally:
        os.close(descriptor)
    _validate_campaign_checkpoint_payloads(
        progress_directory_fd, progress_directory_binding,
        armed=armed, journal_records=journal_records,
        label=f"{label} recovered checkpoints")


def _append_campaign_transcript_allocation(
    lane_fd: int, marker_fd: int, binding: Mapping[str, Any],
    armed: Mapping[str, Any],
) -> None:
    data = contract.canonical_json_bytes(
        _campaign_transcript_allocation_record(armed))
    observed = _read_campaign_transcript_bytes(
        lane_fd, marker_fd, binding,
        "campaign transcript before ARMED allocation")
    _require(len(observed) <= len(data) and data.startswith(observed),
             "campaign transcript is not an exact ARMED allocation prefix")
    written = len(observed)
    while written < len(data):
        count = os.pwrite(marker_fd, data[written:], written)
        _require(count > 0,
                 "campaign transcript allocation write made no progress")
        written += count
    os.fsync(marker_fd)
    _validate_campaign_transcript(
        lane_fd, marker_fd, binding, journal_records=[], armed=armed,
        label="campaign transcript after ARMED allocation")


def _validate_journal_binding(value: Any) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "name", "device", "inode", "uid", "mode", "link_count",
    )), "execution journal binding")
    _require(record["name"] == JOURNAL_DIRECTORY and
             type(record["device"]) is int and record["device"] >= 0 and
             type(record["inode"]) is int and record["inode"] > 0 and
             type(record["uid"]) is int and record["uid"] == os.geteuid() and
             type(record["mode"]) is int and record["mode"] == 0o700 and
             type(record["link_count"]) is int and
             record["link_count"] == 2,
             "execution journal binding is malformed")
    return copy.deepcopy(record)


def _validate_journal_capability(
    attempt_fd: int, journal_fd: int, expected: Mapping[str, Any], label: str,
) -> None:
    binding = _validate_journal_binding(expected)
    try:
        retained = os.fstat(journal_fd)
        pathname = os.stat(
            JOURNAL_DIRECTORY, dir_fd=attempt_fd, follow_symlinks=False)
    except OSError as error:
        raise ArmingError(f"{label} changed") from error
    _require(
        contract.exact_json_equal(
            _journal_binding_record(retained), binding) and
        (retained.st_dev, retained.st_ino) ==
            (pathname.st_dev, pathname.st_ino) and
        stat.S_ISDIR(pathname.st_mode) and
        pathname.st_uid == os.geteuid() and
        pathname.st_nlink == 2 and
        stat.S_IMODE(pathname.st_mode) == 0o700,
        f"{label} inode changed")


def _attempt_manifest_record(
    authority_bundle_data: bytes, journal_binding: Mapping[str, Any],
    launch_marker_binding: Mapping[str, Any],
    campaign_marker_binding: Mapping[str, Any],
    campaign_progress_directory_binding: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": ATTEMPT_MANIFEST_SCHEMA,
        "files": [{
            "name": AUTHORITY_BUNDLE_FILE,
            "sha256": _sha256_bytes(authority_bundle_data),
            "size": len(authority_bundle_data),
        }],
        "journal": _validate_journal_binding(journal_binding),
        "launch_marker": _validate_launch_marker_binding(
            launch_marker_binding),
        "campaign_marker": _validate_campaign_marker_binding(
            campaign_marker_binding),
        "campaign_progress_directory":
            _validate_campaign_progress_directory_binding(
                campaign_progress_directory_binding),
    }


def _validate_attempt_manifest(
    value: Any, authority_bundle_data: bytes,
    journal_metadata: os.stat_result,
) -> dict[str, Any]:
    record = _exact_mapping(
        value, frozenset((
            "schema", "files", "journal", "launch_marker",
            "campaign_marker", "campaign_progress_directory",
        )),
        "attempt manifest")
    expected = _attempt_manifest_record(
        authority_bundle_data, _journal_binding_record(journal_metadata),
        _validate_launch_marker_binding(record["launch_marker"]),
        _validate_campaign_marker_binding(record["campaign_marker"]),
        _validate_campaign_progress_directory_binding(
            record["campaign_progress_directory"]))
    _require(contract.exact_json_equal(record, expected),
             "attempt manifest differs from its fixed point")
    return copy.deepcopy(expected)


def _journal_command_template(
    child: Mapping[str, Any], selected_pair: Mapping[str, Any],
) -> list[str]:
    return [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(selected_pair["benchmark_cpu"]),
        "<sealed-executable>", *child["argv_tail"],
    ]


def _journal_intent_record(
    *, child: Mapping[str, Any], selected_pair: Mapping[str, Any],
    artifact_handle_id: str, evidence_attempt: int,
    prior_journal_sha256: str | None, started_monotonic_ns: int,
    launch_context: Mapping[str, Any],
) -> dict[str, Any]:
    _require(type(started_monotonic_ns) is int and started_monotonic_ns > 0,
             "journal launch timestamp is invalid")
    return {
        "schema": JOURNAL_INTENT_SCHEMA,
        "evidence_attempt": evidence_attempt,
        "child_index": child["index"],
        "logical_child_sha256": contract.canonical_sha256(child),
        "prior_journal_sha256": prior_journal_sha256,
        "implementation": child["implementation"],
        "selected_pair": copy.deepcopy(selected_pair),
        "artifact_handle_id": artifact_handle_id,
        "command_template": _journal_command_template(child, selected_pair),
        "environment_sha256": contract.canonical_sha256(
            dict(_frozen_child_environment())),
        "launch_context_sha256": contract.canonical_sha256(
            _validate_launch_context(launch_context)),
        "cwd": "/",
        "stdin": "DEVNULL",
        "stdout_limit_bytes": MAX_CHILD_STDOUT_BYTES,
        "stderr_limit_bytes": MAX_CHILD_STDERR_BYTES,
        "timeout_seconds": child["timeout_budget"]["timeout_seconds"],
        "creator_pid": os.getpid(),
        "creator_start_ticks": _process_start_ticks(),
        "started_monotonic_ns": started_monotonic_ns,
    }


def _validate_journal_intent(
    value: Any, *, child: Mapping[str, Any],
    selected_pair: Mapping[str, Any], artifact_handle_id: str,
    evidence_attempt: int, prior_journal_sha256: str | None,
    launch_context: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "schema", "evidence_attempt", "child_index", "logical_child_sha256",
        "prior_journal_sha256", "implementation", "selected_pair",
        "artifact_handle_id", "command_template", "environment_sha256",
        "launch_context_sha256", "cwd", "stdin", "stdout_limit_bytes",
        "stderr_limit_bytes",
        "timeout_seconds", "creator_pid", "creator_start_ticks",
        "started_monotonic_ns",
    )), "launch intent")
    _require(type(record["creator_pid"]) is int and
             record["creator_pid"] > 0 and
             type(record["creator_start_ticks"]) is int and
             record["creator_start_ticks"] > 0,
             "launch intent creator identity differs")
    expected = _journal_intent_record(
        child=child, selected_pair=selected_pair,
        artifact_handle_id=artifact_handle_id,
        evidence_attempt=evidence_attempt,
        prior_journal_sha256=prior_journal_sha256,
        started_monotonic_ns=record["started_monotonic_ns"],
        launch_context=launch_context)
    expected["creator_pid"] = record["creator_pid"]
    expected["creator_start_ticks"] = record["creator_start_ticks"]
    _require(contract.exact_json_equal(record, expected),
             "launch intent differs from the frozen child")
    return copy.deepcopy(expected)


def _validate_launch_marker_state(
    *, attempt_fd: int, marker_fd: int,
    marker_binding: Mapping[str, Any], journal_fd: int,
    lane_fd: int, campaign_marker_fd: int,
    campaign_marker_binding: Mapping[str, Any],
    logical_plan: Mapping[str, Any], artifacts: Mapping[str, Any],
    armed: Mapping[str, Any], launch_context: Mapping[str, Any], label: str,
    require_allocation: bool = True,
) -> bool:
    metadata = _validate_launch_marker_capability(
        attempt_fd, marker_fd, marker_binding, label)
    campaign_metadata = _validate_campaign_marker_capability(
        lane_fd, campaign_marker_fd, campaign_marker_binding,
        f"{label} lane slot")
    transcript = _read_campaign_transcript(
        lane_fd, campaign_marker_fd, campaign_marker_binding,
        f"{label} lane transcript")
    allocated = transcript["allocation"] is not None
    launched = bool(transcript["journal_entries"])
    _require((allocated or not require_allocation) and
             (campaign_metadata.st_size > 0) == allocated and
             (metadata.st_size > 0) == launched,
             f"{label} allocation or launch state differs")
    if metadata.st_size == 0:
        return False
    marker_data = _read_exact_fd(marker_fd, MAX_JSON_BYTES, label)
    child = logical_plan["child_plans"][0]
    role = child["implementation"]
    marker_intent = _validate_journal_intent(
        _canonical_json_document(marker_data, label), child=child,
        selected_pair=armed["selected_pair"],
        artifact_handle_id=artifacts["roles"][role]["handle_id"],
        evidence_attempt=armed["evidence_attempt"],
        prior_journal_sha256=None, launch_context=launch_context)
    try:
        intent_fd, intent_data = _read_immutable_file_at(
            journal_fd, "intent-0000.json", MAX_JSON_BYTES,
            f"{label} first journal intent")
    except OSError as error:
        raise ArmingError(
            f"{label} has no retained first journal intent") from error
    try:
        _require(marker_data == intent_data and
                 transcript["journal_entries"][0]["name"] ==
                    "intent-0000.json" and
                 transcript["journal_entries"][0]["data_sha256"] ==
                    _sha256_bytes(intent_data) and
                 contract.exact_json_equal(
                     marker_intent,
                     _canonical_json_document(
                         intent_data, f"{label} first journal intent")),
                 f"{label} and first journal intent differ")
    finally:
        os.close(intent_fd)
    return True


def _journal_result_record(
    *, child_index: int, intent_sha256: str, outcome: str,
    finished_monotonic_ns: int, elapsed_ns: int,
    returncode: int | None, stdout_sha256: str, stderr_sha256: str,
    payload: Mapping[str, Any] | None, normalized: Mapping[str, Any] | None,
    error: str | None,
) -> dict[str, Any]:
    _require(outcome in {
        "accepted", "timeout", "nonzero", "output-limit", "output-rejected"},
        "child result outcome differs")
    return {
        "schema": JOURNAL_RESULT_SCHEMA,
        "child_index": child_index,
        "intent_sha256": _hex_digest(intent_sha256, "launch intent hash"),
        "outcome": outcome,
        "finished_monotonic_ns": finished_monotonic_ns,
        "elapsed_ns": elapsed_ns,
        "returncode": returncode,
        "stdout_sha256": _hex_digest(stdout_sha256, "child stdout hash"),
        "stderr_sha256": _hex_digest(stderr_sha256, "child stderr hash"),
        "payload": copy.deepcopy(payload),
        "normalized": copy.deepcopy(normalized),
        "error": error,
    }


def _normalize_gen3_result(
    implementation: str, payload: Mapping[str, Any], cell: Mapping[str, Any],
    registration: Mapping[str, Any],
) -> dict[str, Any]:
    """Normalize one result while attributing the bytes Gen3 really execs."""
    normalized = result_contract.validate_result(
        implementation, payload, cell,
        registration["candidate_source"]["commit"],
        registration["candidate_source"]["tree"],
        registration["campaign"]["iterations_per_child"],
        registration["campaign"]["warmup_per_child"])
    if implementation == "main":
        attribution = normalized.get("exact_main_authority_attribution")
        _require(
            type(attribution) is dict and set(attribution) == {
                "pure_avx2", "main_source_commit",
                "authority_record_sha256", "executable_sha256"} and
            attribution["pure_avx2"] is True and
            attribution["main_source_commit"] ==
                result_contract.EXACT_MAIN_SOURCE_COMMIT and
            attribution["authority_record_sha256"] ==
                prereg.EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
            attribution["executable_sha256"] == CANONICAL_RAW_SHA256 and
            PATH_VARIANT_RAW_SHA256 != CANONICAL_RAW_SHA256,
            "shared exact-main normalization attribution changed")
        attribution["executable_sha256"] = PATH_VARIANT_RAW_SHA256
    return normalized


def _validate_journal_result(
    value: Any, *, child: Mapping[str, Any], intent_sha256: str,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_mapping(value, frozenset((
        "schema", "child_index", "intent_sha256", "outcome",
        "finished_monotonic_ns", "elapsed_ns", "returncode",
        "stdout_sha256", "stderr_sha256", "payload", "normalized", "error",
    )), "child result")
    _require(record["schema"] == JOURNAL_RESULT_SCHEMA and
             record["child_index"] == child["index"] and
             record["intent_sha256"] == intent_sha256 and
             type(record["finished_monotonic_ns"]) is int and
             record["finished_monotonic_ns"] > 0 and
             type(record["elapsed_ns"]) is int and record["elapsed_ns"] >= 0 and
             type(record["stdout_sha256"]) is str and
             type(record["stderr_sha256"]) is str,
             "child result metadata differs")
    _hex_digest(record["stdout_sha256"], "child stdout hash")
    _hex_digest(record["stderr_sha256"], "child stderr hash")
    if record["outcome"] == "accepted":
        _require(record["returncode"] == 0 and type(record["payload"]) is dict and
                 type(record["normalized"]) is dict and record["error"] is None,
                 "accepted child result is incomplete")
        cell = logical_plan["cells"][child["cell_index"]]
        expected_normalized = _normalize_gen3_result(
            child["implementation"], record["payload"], cell,
            registration)
        _require(contract.exact_json_equal(
            record["normalized"], expected_normalized),
            "accepted child normalization changed")
    else:
        _require(record["outcome"] in {
            "timeout", "nonzero", "output-limit", "output-rejected"} and
            record["payload"] is None and record["normalized"] is None and
            type(record["error"]) is str and 0 < len(record["error"]) <= 4096,
            "rejected child result is malformed")
        if record["outcome"] == "nonzero":
            _require(type(record["returncode"]) is int and
                     record["returncode"] != 0,
                     "nonzero child result lacks its return code")
        else:
            _require(record["returncode"] is None or
                     type(record["returncode"]) is int,
                     "rejected child return code differs")
    return copy.deepcopy(record)


def _validated_journal_staging_files(
    journal_fd: int, names: Sequence[str],
) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    for name in sorted(names):
        match = JOURNAL_STAGING_NAME.fullmatch(name)
        if match is None:
            continue
        descriptor = os.open(
            name, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=journal_fd)
        try:
            metadata = os.fstat(descriptor)
            pathname = os.stat(
                name, dir_fd=journal_fd, follow_symlinks=False)
            _require(
                stat.S_ISREG(metadata.st_mode) and
                metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
                stat.S_IMODE(metadata.st_mode) in (0o400, 0o600) and
                0 <= metadata.st_size <= MAX_JSON_BYTES and
                (metadata.st_dev, metadata.st_ino) ==
                    (pathname.st_dev, pathname.st_ino),
                "execution journal staging entry is unsafe")
        finally:
            os.close(descriptor)
        result.append((name, match.group(1)))
    _require(len(result) <= 1,
             "execution journal has conflicting staging entries")
    return result


def _validate_journal_directory(
    attempt_fd: int, journal_fd: int, *,
    journal_binding: Mapping[str, Any], registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], artifacts: Mapping[str, Any],
    armed: Mapping[str, Any], launch_context: Mapping[str, Any],
    cleanup_staging: bool = True,
) -> tuple[int, str | None, bool, list[tuple[str, str]]]:
    _validate_journal_capability(
        attempt_fd, journal_fd, journal_binding, "execution journal")
    try:
        metadata = os.fstat(journal_fd)
        _require(stat.S_ISDIR(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 stat.S_IMODE(metadata.st_mode) == 0o700,
                 "execution journal directory is unsafe")
        all_names = set(os.listdir(journal_fd))
        staging = _validated_journal_staging_files(
            journal_fd, sorted(all_names))
        staging_names = {name for name, unused_target in staging}
        names = all_names - staging_names
        allowed_names = {"complete.json"}
        for index in range(len(logical_plan["child_plans"])):
            allowed_names.add(f"intent-{index:04d}.json")
            allowed_names.add(f"result-{index:04d}.json")
        _require(names <= allowed_names,
                 "execution journal contains an unrecognized entry")
        chain: str | None = None
        completed = 0
        trailing_intent = False
        terminal_failure = False
        journal_records: list[tuple[str, str]] = []
        roles = artifacts["roles"]
        for index, child in enumerate(logical_plan["child_plans"]):
            intent_name = f"intent-{index:04d}.json"
            result_name = f"result-{index:04d}.json"
            has_intent = intent_name in names
            has_result = result_name in names
            _require(not has_result or has_intent,
                     "execution journal result has no launch intent")
            if not has_intent:
                _require(not has_result, "execution journal has a gap")
                break
            _require(not trailing_intent,
                     "execution journal continued after an incomplete child")
            intent_fd, intent_data = _read_immutable_file_at(
                journal_fd, intent_name, MAX_JSON_BYTES,
                f"journal intent {index}")
            os.close(intent_fd)
            intent = _validate_journal_intent(
                _canonical_json_document(intent_data, f"journal intent {index}"),
                child=child, selected_pair=armed["selected_pair"],
                artifact_handle_id=roles[child["implementation"]]["handle_id"],
                evidence_attempt=armed["evidence_attempt"],
                prior_journal_sha256=chain,
                launch_context=launch_context)
            intent_sha = contract.canonical_sha256(intent)
            _require(intent_sha == _sha256_bytes(intent_data),
                     "journal intent canonical hash differs")
            journal_records.append((intent_name, intent_sha))
            if not has_result:
                chain = intent_sha
                trailing_intent = True
                break
            result_fd, result_data = _read_immutable_file_at(
                journal_fd, result_name, MAX_JSON_BYTES,
                f"journal result {index}")
            os.close(result_fd)
            result = _validate_journal_result(
                _canonical_json_document(result_data, f"journal result {index}"),
                child=child, intent_sha256=intent_sha,
                registration=registration, logical_plan=logical_plan)
            chain = contract.canonical_sha256(result)
            _require(chain == _sha256_bytes(result_data),
                     "journal result canonical hash differs")
            journal_records.append((result_name, chain))
            completed += 1
            if result["outcome"] != "accepted":
                terminal_failure = True
                _require(index == len(logical_plan["child_plans"]) - 1 or
                         not any(
                            f"intent-{later:04d}.json" in names
                            for later in range(index + 1,
                                               len(logical_plan["child_plans"]))),
                         "execution journal continued after terminal failure")
                break
        expected_entries = 2 * completed + (1 if trailing_intent else 0)
        observed_entries = len(names - {"complete.json"})
        _require(observed_entries == expected_entries,
                 "execution journal is not one gapless prefix")
        complete = "complete.json" in names
        if complete:
            complete_fd, complete_data = _read_immutable_file_at(
                journal_fd, "complete.json", MAX_JSON_BYTES,
                "campaign-complete record")
            os.close(complete_fd)
            record = _exact_mapping(
                _canonical_json_document(
                    complete_data, "campaign-complete record"),
                frozenset(("schema", "evidence_attempt", "child_count",
                           "prior_journal_sha256", "completed_monotonic_ns")),
                "campaign-complete record")
            _require(record["schema"] == JOURNAL_COMPLETE_SCHEMA and
                     record["evidence_attempt"] == armed["evidence_attempt"] and
                     record["child_count"] == len(logical_plan["child_plans"]) and
                     record["prior_journal_sha256"] == chain and
                     type(record["completed_monotonic_ns"]) is int and
                     record["completed_monotonic_ns"] > 0 and
                     completed == len(logical_plan["child_plans"]) and
                     not trailing_intent and not terminal_failure,
                     "campaign-complete record differs")
            chain = contract.canonical_sha256(record)
            _require(chain == _sha256_bytes(complete_data),
                     "campaign-complete canonical hash differs")
            journal_records.append(("complete.json", chain))
        next_name: str | None
        if complete or terminal_failure:
            next_name = None
        elif trailing_intent:
            next_name = f"result-{completed:04d}.json"
        elif completed < len(logical_plan["child_plans"]):
            next_name = f"intent-{completed:04d}.json"
        else:
            next_name = "complete.json"
        if staging:
            staging_name, target_name = staging[0]
            _require(next_name is not None and target_name == next_name and
                     target_name not in names,
                     "execution journal staging target is not the unique next "
                     "record")
            if cleanup_staging:
                os.unlink(staging_name, dir_fd=journal_fd)
                os.fsync(journal_fd)
        _validate_journal_capability(
            attempt_fd, journal_fd, journal_binding, "execution journal")
        return completed, chain, complete, journal_records
    except OSError as error:
        raise ArmingError("execution journal could not be replayed") from error


def _validate_authority_bundle_and_armed(
    bundle_value: Any, armed_value: Any, *, lane: Path,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    candidate_authority: Mapping[str, Any], exact_main: Mapping[str, Any],
    lane_binding: Mapping[str, Any], authority_bundle_data: bytes,
    attempt_manifest_data: bytes,
) -> tuple[dict[str, Any], dict[str, Any]]:
    bundle = _bundle_components(bundle_value)
    _require(contract.exact_json_equal(bundle["lane_binding"], lane_binding),
             "authority bundle lane binding differs")
    retained_registration = prereg.validate_preregistration(
        bundle["preregistration"], verify_files=True)
    _require(contract.exact_json_equal(retained_registration, registration),
             "authority bundle preregistration differs")
    retained_plan = plan_contract.validate_campaign_plan(
        bundle["logical_plan"], retained_registration)
    _require(contract.exact_json_equal(retained_plan, logical_plan),
             "authority bundle logical plan differs")
    budget_commit = _validate_budget_commit(
        bundle["budget_commit"], retained_registration)
    history_reference = budget_commit["history"]
    _require(
        1 <= history_reference["acquisition_index"] <=
            len(lane_binding["prearmed_history_markers"]) and
        contract.exact_json_equal(
            history_reference["marker"],
            lane_binding["prearmed_history_markers"]
                        [history_reference["acquisition_index"] - 1]),
        "authority bundle pre-ARMED marker differs from its lane")
    _require(contract.exact_json_equal(
        bundle["launch_context"], retained_registration["launch_context"]),
        "authority bundle launch context differs from preregistration")
    host = execution.validate_host_instance(bundle["host_instance"])
    artifacts = execution.validate_artifact_bundle(
        bundle["artifact_bundle"], retained_registration,
        host_authority=host["authority"])
    qualification = execution.validate_qualification_binding(
        bundle["qualification_binding"], retained_registration, host,
        evidence=bundle["qualification_evidence_wire"])
    _validate_descriptor_bundle(bundle["descriptor_binding"], artifacts)
    _require(
        contract.exact_json_equal(
            bundle["candidate_authority_record"],
            candidate_authority["record"]) and
        contract.exact_json_equal(
            bundle["exact_main_authority_record"], exact_main["record"]) and
        contract.exact_json_equal(
            bundle["exact_main_verifier"], exact_main["verdict"]) and
        artifacts["candidate_source_authority"]["authority_record_sha256"] ==
            _sha256_bytes(candidate_authority["record_data"]) and
        artifacts["candidate_source_authority"]["authority_ledger_sha256"] ==
            candidate_authority["ledger"]["sha256"] and
        artifacts["exact_main_authority"]["verifier_verdict_sha256"] ==
            exact_main["verdict"]["verdict_sha256"],
        "authority bundle candidate or exact-main authority differs")
    armed = execution.validate_armed_record(
        armed_value, retained_registration, retained_plan, host, artifacts,
        qualification,
        qualification_evidence=bundle["qualification_evidence_wire"])
    _require(
        armed["authority_bundle_sha256"] ==
            _sha256_bytes(authority_bundle_data) and
        armed["attempt_manifest_sha256"] ==
            _sha256_bytes(attempt_manifest_data) and
        armed["lane_binding_sha256"] ==
            contract.canonical_sha256(lane_binding) and
        budget_commit["prospective_ledger"]["evidence_attempts_used"] ==
            armed["evidence_attempt"] and
        contract.exact_json_equal(
            budget_commit["prospective_ledger"]["frozen_pair"],
            armed["selected_pair"]),
        "ARMED record does not bind its retained attempt directory")
    return copy.deepcopy(bundle_value), armed


def _clean_staging_directories(
    attempts_fd: int, *, remove_complete: bool = True,
) -> None:
    changed = False
    names = sorted(
        name for name in os.listdir(attempts_fd)
        if STAGING_NAME.fullmatch(name) is not None)
    # No recovery mutation is permitted until every remnant is known-safe.
    # The acquisition preflight admits at most one, but this helper also has
    # direct test coverage and must remain atomic if called independently.
    for name in names:
        _validate_attempt_staging_inventory(
            attempts_fd, name, remove_complete=remove_complete)
    for name in names:
        staging_fd = os.open(
            name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC |
            os.O_NOFOLLOW, dir_fd=attempts_fd)
        try:
            metadata = os.fstat(staging_fd)
            _require(stat.S_ISDIR(metadata.st_mode) and
                     metadata.st_uid == os.geteuid() and
                     stat.S_IMODE(metadata.st_mode) in (0o700, 0o500),
                     "stale attempt staging directory is unsafe")
            children = sorted(os.listdir(staging_fd))
            _require(set(children) <= {
                AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
                JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE},
                "stale attempt staging directory contains an unsafe entry")
            _require(remove_complete or not (
                stat.S_IMODE(metadata.st_mode) == 0o500 and
                set(children) == {
                    AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
                    JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE}),
                "complete attempt staging cannot be classified as unpublished")
            # Validate the complete crash remnant before removing any part of
            # it.  A crash may leave a just-created output at mode 0600 (and
            # with zero or partial contents), or the staging directory may
            # already have reached its final read-only mode.  Both are safe
            # unpublished states; an unexpected object poisons the lane
            # without partially erasing the evidence that made it unsafe.
            for child in children:
                child_metadata = os.stat(
                    child, dir_fd=staging_fd, follow_symlinks=False)
                if child == JOURNAL_DIRECTORY:
                    _require(stat.S_ISDIR(child_metadata.st_mode) and
                             child_metadata.st_uid == os.geteuid() and
                             stat.S_IMODE(child_metadata.st_mode) == 0o700,
                             "stale attempt journal directory is unsafe")
                    journal_fd = os.open(
                        child, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC |
                        os.O_NOFOLLOW, dir_fd=staging_fd)
                    try:
                        _require(not os.listdir(journal_fd),
                                 "stale unpublished journal is not empty")
                    finally:
                        os.close(journal_fd)
                    continue
                _require(stat.S_ISREG(child_metadata.st_mode) and
                         child_metadata.st_uid == os.geteuid() and
                         child_metadata.st_nlink == 1 and
                         stat.S_IMODE(child_metadata.st_mode) in (0o600, 0o400)
                         and 0 <= child_metadata.st_size <= MAX_JSON_BYTES,
                         "stale attempt staging file is unsafe")
                if child == LAUNCH_CONSUMED_FILE:
                    _require(child_metadata.st_size == 0,
                             "stale unpublished launch marker is consumed")
            os.fchmod(staging_fd, 0o700)
            os.fsync(staging_fd)
            if JOURNAL_DIRECTORY in children:
                os.rmdir(JOURNAL_DIRECTORY, dir_fd=staging_fd)
            for child in children:
                if child != JOURNAL_DIRECTORY:
                    os.unlink(child, dir_fd=staging_fd)
            os.fsync(staging_fd)
        finally:
            os.close(staging_fd)
        os.rmdir(name, dir_fd=attempts_fd)
        changed = True
    if changed:
        os.fsync(attempts_fd)


def _read_attempt_directory(
    attempts_fd: int, name: str, *, lane: Path, lane_fd: int,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    candidate_authority: Mapping[str, Any], exact_main: Mapping[str, Any],
    lane_binding: Mapping[str, Any], require_allocation: bool = True,
    read_only: bool = False,
) -> tuple[
    dict[str, Any], dict[str, Any],
    tuple[int, str | None, bool, list[tuple[str, str]]],
]:
    try:
        directory_fd = os.open(
            name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=attempts_fd)
    except OSError as error:
        raise ArmingError(f"retained {name} directory cannot be opened") \
            from error
    descriptors: list[int] = []
    try:
        metadata = os.fstat(directory_fd)
        _require(stat.S_ISDIR(metadata.st_mode) and
                 metadata.st_uid == os.geteuid() and
                 stat.S_IMODE(metadata.st_mode) == 0o500 and
                 set(os.listdir(directory_fd)) == {
                    AUTHORITY_BUNDLE_FILE, ATTEMPT_MANIFEST_FILE, ARMED_FILE,
                    JOURNAL_DIRECTORY, LAUNCH_CONSUMED_FILE},
                 f"retained {name} directory is not exact and immutable")
        bundle_fd, bundle_data = _read_immutable_file_at(
            directory_fd, AUTHORITY_BUNDLE_FILE, MAX_JSON_BYTES,
            f"retained {name} authority bundle")
        descriptors.append(bundle_fd)
        manifest_fd, manifest_data = _read_immutable_file_at(
            directory_fd, ATTEMPT_MANIFEST_FILE, MAX_JSON_BYTES,
            f"retained {name} manifest")
        descriptors.append(manifest_fd)
        armed_fd, armed_data = _read_immutable_file_at(
            directory_fd, ARMED_FILE, MAX_JSON_BYTES,
            f"retained {name} ARMED record")
        descriptors.append(armed_fd)
        marker_fd = os.open(
            LAUNCH_CONSUMED_FILE,
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        descriptors.append(marker_fd)
        journal_fd = os.open(
            JOURNAL_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory_fd)
        descriptors.append(journal_fd)
        journal_metadata = os.fstat(journal_fd)
        manifest = _validate_attempt_manifest(
            _canonical_json_document(manifest_data, f"retained {name} manifest"),
            bundle_data, journal_metadata)
        campaign_marker_binding = manifest["campaign_marker"]
        campaign_marker_fd = os.open(
            campaign_marker_binding["name"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=lane_fd)
        descriptors.append(campaign_marker_fd)
        _validate_journal_capability(
            directory_fd, journal_fd, manifest["journal"],
            f"retained {name} journal")
        bundle, armed = _validate_authority_bundle_and_armed(
            _canonical_json_document(bundle_data, f"retained {name} bundle"),
            _canonical_json_document(armed_data, f"retained {name} ARMED"),
            lane=lane, registration=registration, logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding, authority_bundle_data=bundle_data,
            attempt_manifest_data=manifest_data)
        attempt_index = armed["evidence_attempt"]
        _require(type(attempt_index) is int and
                 1 <= attempt_index <= len(lane_binding["campaign_markers"]) and
                 contract.exact_json_equal(
                    campaign_marker_binding,
                    lane_binding["campaign_markers"][attempt_index - 1]),
                 f"retained {name} transcript slot is not its lane anchor")
        progress_binding = _validate_campaign_progress_directory_binding(
            manifest["campaign_progress_directory"])
        _require(
            contract.exact_json_equal(
                progress_binding,
                lane_binding["campaign_progress_directories"]
                            [attempt_index - 1]),
            f"retained {name} checkpoint directory is not its lane anchor")
        progress_fd, observed_progress_binding = \
            _open_campaign_progress_directory(
                lane, lane_fd, progress_binding,
                f"retained {name} campaign checkpoints")
        descriptors.append(progress_fd)
        _require(contract.exact_json_equal(
            observed_progress_binding, progress_binding),
            f"retained {name} checkpoint directory binding changed")
        if read_only and not require_allocation:
            _validate_recoverable_published_campaign_marker(
                lane_fd, campaign_marker_fd, campaign_marker_binding,
                f"retained {name} campaign transcript")
        elif read_only:
            _validate_campaign_marker_capability(
                lane_fd, campaign_marker_fd, campaign_marker_binding,
                f"retained {name} campaign transcript")
        else:
            _recover_published_campaign_marker_mode(
                lane_fd, campaign_marker_fd, campaign_marker_binding,
                f"retained {name} campaign transcript")
        journal_status = _validate_journal_directory(
            directory_fd, journal_fd, journal_binding=manifest["journal"],
            registration=registration,
            logical_plan=logical_plan, artifacts=bundle["artifact_bundle"],
            armed=armed, launch_context=bundle["launch_context"],
            cleanup_staging=require_allocation and not read_only)
        if require_allocation and read_only:
            _require(not any(
                JOURNAL_STAGING_NAME.fullmatch(entry) is not None
                for entry in os.listdir(journal_fd)),
                f"retained {name} journal contains a staging remnant")
        transcript_label = f"retained {name} campaign transcript"
        transcript_data = _read_campaign_transcript_bytes(
            lane_fd, campaign_marker_fd, campaign_marker_binding,
            transcript_label,
            recoverable=read_only and not require_allocation)
        expected_transcript = _expected_campaign_transcript_bytes(
            armed, journal_status[3])
        checkpoint = _read_campaign_checkpoint_frontier(
            progress_fd, progress_binding,
            f"retained {name} campaign checkpoints")
        if checkpoint["frontier"] == 0:
            _require(
                STAGING_NAME.fullmatch(name) is not None and
                not journal_status[3] and not transcript_data,
                f"retained {name} bypassed its ARMED checkpoint")
        else:
            checkpoint = _validate_campaign_checkpoint_payloads(
                progress_fd, progress_binding, armed=armed,
                journal_records=journal_status[3],
                label=f"retained {name} campaign checkpoints",
                allow_recoverable_tail=not require_allocation)
            checkpoint_prefix = checkpoint["payload_prefix"]
            _require(checkpoint_prefix.startswith(transcript_data),
                     f"{transcript_label} advanced beyond or differs from "
                     "its exact checkpoints")
        _require(expected_transcript.startswith(transcript_data),
                 f"{transcript_label} is not an exact durable-journal prefix")
        marker_label = f"retained {name} launch-consumed marker"
        marker_validator = (
            _validate_launch_marker_capability
            if read_only and require_allocation else
            _validate_recoverable_launch_marker_capability)
        marker_metadata = marker_validator(
            directory_fd, marker_fd, manifest["launch_marker"], marker_label)
        marker_data = b"" if marker_metadata.st_size == 0 else _read_exact_fd(
            marker_fd, MAX_JSON_BYTES, marker_label)
        expected_marker_data: bytes | None = None
        if journal_status[3]:
            _require(journal_status[3][0][0] == "intent-0000.json",
                     f"retained {name} journal does not begin at child zero")
            first_fd, expected_marker_data = _read_immutable_file_at(
                journal_fd, "intent-0000.json", MAX_JSON_BYTES,
                f"retained {name} first journal intent")
            os.close(first_fd)
            _require(_sha256_bytes(expected_marker_data) ==
                     journal_status[3][0][1],
                     f"retained {name} first journal intent hash changed")
        _require(
            (not marker_data if expected_marker_data is None else
             expected_marker_data.startswith(marker_data)),
            f"{marker_label} is not an exact durable-intent prefix")
        if not require_allocation:
            if transcript_data == expected_transcript and not read_only:
                _validate_campaign_transcript(
                    lane_fd, campaign_marker_fd, campaign_marker_binding,
                    journal_records=journal_status[3], armed=armed,
                    label=transcript_label)
            marker_consumed = expected_marker_data is not None
        else:
            _validate_campaign_transcript(
                lane_fd, campaign_marker_fd, campaign_marker_binding,
                journal_records=journal_status[3], armed=armed,
                label=transcript_label)
            marker_consumed = _validate_launch_marker_state(
                attempt_fd=directory_fd, marker_fd=marker_fd,
                marker_binding=manifest["launch_marker"],
                journal_fd=journal_fd, lane_fd=lane_fd,
                campaign_marker_fd=campaign_marker_fd,
                campaign_marker_binding=campaign_marker_binding,
                logical_plan=logical_plan,
                artifacts=bundle["artifact_bundle"], armed=armed,
                launch_context=bundle["launch_context"],
                label=f"retained {name} launch-consumed marker")
        if not marker_consumed:
            _require(journal_status[1] is None,
                     f"retained {name} journal bypassed its launch marker")
        return bundle, armed, journal_status
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
        os.close(directory_fd)


def _reconcile_attempt_execution_prefixes(
    attempts_fd: int, name: str, *, lane: Path, lane_fd: int,
    lane_binding: Mapping[str, Any], bundle: Mapping[str, Any],
    armed: Mapping[str, Any], registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any],
    journal_status: tuple[
        int, str | None, bool, list[tuple[str, str]]],
) -> None:
    attempt_fd = journal_fd = bundle_fd = manifest_fd = progress_fd = -1
    try:
        attempt_fd = os.open(
            name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=attempts_fd)
        bundle_fd, bundle_data = _read_immutable_file_at(
            attempt_fd, AUTHORITY_BUNDLE_FILE, MAX_JSON_BYTES,
            f"retained {name} authority bundle during recovery")
        manifest_fd, manifest_data = _read_immutable_file_at(
            attempt_fd, ATTEMPT_MANIFEST_FILE, MAX_JSON_BYTES,
            f"retained {name} manifest during recovery")
        journal_fd = os.open(
            JOURNAL_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=attempt_fd)
        manifest = _validate_attempt_manifest(
            _canonical_json_document(
                manifest_data, f"retained {name} manifest during recovery"),
            bundle_data, os.fstat(journal_fd))
        cleaned_status = _validate_journal_directory(
            attempt_fd, journal_fd, journal_binding=manifest["journal"],
            registration=registration, logical_plan=logical_plan,
            artifacts=bundle["artifact_bundle"], armed=armed,
            launch_context=bundle["launch_context"], cleanup_staging=True)
        _require(cleaned_status == journal_status,
                 f"retained {name} journal changed during staging cleanup")
        marker_binding = _validate_campaign_marker_binding(
            lane_binding["campaign_markers"][armed["evidence_attempt"] - 1])
        _require(contract.exact_json_equal(
            manifest["campaign_marker"], marker_binding),
            f"retained {name} recovery transcript binding differs")
        progress_binding = _validate_campaign_progress_directory_binding(
            lane_binding["campaign_progress_directories"]
                        [armed["evidence_attempt"] - 1])
        _require(contract.exact_json_equal(
            manifest["campaign_progress_directory"], progress_binding),
            f"retained {name} recovery checkpoint binding differs")
        progress_fd, observed_progress_binding = \
            _open_campaign_progress_directory(
                lane, lane_fd, progress_binding,
                f"retained {name} recovery checkpoints")
        _require(contract.exact_json_equal(
            observed_progress_binding, progress_binding),
            f"retained {name} recovery checkpoint directory changed")
        _reconcile_campaign_transcript_prefix(
            lane_fd, marker_binding, armed=armed,
            journal_records=journal_status[3],
            label=f"retained {name} campaign transcript recovery",
            progress_directory_fd=progress_fd,
            progress_directory_binding=progress_binding)
        expected_marker_data: bytes | None = None
        if journal_status[3]:
            _require(journal_status[3][0][0] == "intent-0000.json",
                     f"retained {name} recovery journal begins after child zero")
            first_fd, expected_marker_data = _read_immutable_file_at(
                journal_fd, "intent-0000.json", MAX_JSON_BYTES,
                f"retained {name} first intent during recovery")
            try:
                _require(_sha256_bytes(expected_marker_data) ==
                         journal_status[3][0][1],
                         f"retained {name} recovery intent hash changed")
            finally:
                os.close(first_fd)
        _reconcile_launch_marker_prefix(
            attempt_fd, manifest["launch_marker"],
            expected_data=expected_marker_data,
            label=f"retained {name} launch marker recovery")
    finally:
        for descriptor in (
                progress_fd, journal_fd, manifest_fd, bundle_fd, attempt_fd):
            if descriptor >= 0:
                os.close(descriptor)


def _validate_prior_attempt_invariants(
    attempts: Sequence[tuple[
        Mapping[str, Any], Mapping[str, Any],
        tuple[int, str | None, bool, list[tuple[str, str]]],
    ]],
    *, reject_consumed: bool = True,
) -> None:
    prior_records: list[Mapping[str, Any]] = []
    for index, (bundle, armed, journal_status) in enumerate(attempts, 1):
        if reject_consumed:
            _require(
                journal_status[1] is None,
                "a durable launch intent permanently consumes the campaign")
        elif journal_status[1] is not None:
            _require(index == len(attempts),
                     "an attempt followed a durably consumed campaign")
        expected_prior = None if index == 1 else contract.canonical_sha256(
            prior_records[-1])
        _require(armed["evidence_attempt"] == index and
                 armed["prior_armed_chain_sha256"] == expected_prior,
                 f"attempt {index} ARMED chain differs")
        expected_frozen_pair = (
            None if index == 1 else prior_records[0]["selected_pair"])
        _require(contract.exact_json_equal(
            bundle["qualification_evidence"]["expected_frozen_pair"],
            expected_frozen_pair),
            f"attempt {index} qualification frozen-pair input differs")
        if prior_records:
            first = prior_records[0]
            for field in (
                "preregistration_sha256", "plan_sha256",
                "host_instance_sha256", "controller_bindings_sha256",
                "candidate_source", "selected_pair", "lane_binding_sha256",
            ):
                _require(contract.exact_json_equal(armed[field], first[field]),
                         f"attempt {index} changed frozen {field}")
            first_bundle = attempts[0][0]
            for field in (
                "candidate_authority_record", "exact_main_authority_record",
                "exact_main_verifier",
            ):
                _require(contract.exact_json_equal(
                    bundle[field], first_bundle[field]),
                    f"attempt {index} changed frozen {field}")
        prior_records.append(armed)


def _read_prior_attempts_before_recovery(
    attempts_fd: int, *, lane: Path, lane_fd: int,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    candidate_authority: Mapping[str, Any], exact_main: Mapping[str, Any],
    lane_binding: Mapping[str, Any], require_exact: bool = False,
) -> tuple[
    list[str],
    list[tuple[
        dict[str, Any], dict[str, Any],
        tuple[int, str | None, bool, list[tuple[str, str]]],
    ]],
    tuple[
        str,
        tuple[
            dict[str, Any], dict[str, Any],
            tuple[int, str | None, bool, list[tuple[str, str]]],
        ],
    ] | None,
]:
    metadata = os.fstat(attempts_fd)
    _require(stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and
             stat.S_IMODE(metadata.st_mode) in (0o500, 0o700),
             "generation-3 attempts directory is unsafe")
    all_names = sorted(os.listdir(attempts_fd))
    names = [
        name for name in all_names if ATTEMPT_NAME.fullmatch(name) is not None]
    staging_names = [
        name for name in all_names if STAGING_NAME.fullmatch(name) is not None]
    _require(len(names) + len(staging_names) == len(all_names),
             "attempts directory contains an unrecognized entry")
    indexes = [int(ATTEMPT_NAME.fullmatch(name).group(1)) for name in names]
    limit = len(lane_binding["campaign_markers"])
    _require(indexes == list(range(1, len(indexes) + 1)) and
             len(indexes) <= limit,
             "ARMED attempt directories are not a gapless in-budget prefix")
    _require(len(staging_names) <= 1,
             "attempts directory has conflicting staging remnants")
    _require(not require_exact or not staging_names,
             "live attempts inventory contains a staging remnant")
    if staging_names:
        staging_match = STAGING_NAME.fullmatch(staging_names[0])
        _require(staging_match is not None and
                 int(staging_match.group(1)) == len(names) + 1 and
                 len(names) + 1 <= limit,
                 "attempt staging is not the unique next in-budget attempt")
        _validate_attempt_staging_inventory(attempts_fd, staging_names[0])
    preliminary = [
        _read_attempt_directory(
            attempts_fd, name, lane=lane, registration=registration,
            lane_fd=lane_fd, logical_plan=logical_plan,
            candidate_authority=candidate_authority,
            exact_main=exact_main, lane_binding=lane_binding,
            require_allocation=require_exact, read_only=True)
        for name in names
    ]
    prepared: tuple[
        str,
        tuple[
            dict[str, Any], dict[str, Any],
            tuple[int, str | None, bool, list[tuple[str, str]]],
        ],
    ] | None = None
    next_index = len(names) + 1
    for index, progress_value in enumerate(
            lane_binding["campaign_progress_directories"], 1):
        if index <= len(names):
            continue
        progress_binding = _validate_campaign_progress_directory_binding(
            progress_value)
        progress_fd = -1
        try:
            progress_fd, observed_progress_binding = \
                _open_campaign_progress_directory(
                    lane, lane_fd, progress_binding,
                    f"unallocated attempt {index} campaign checkpoints")
            _require(contract.exact_json_equal(
                observed_progress_binding, progress_binding),
                f"unallocated attempt {index} checkpoint binding changed")
            checkpoint = _read_campaign_checkpoint_frontier(
                progress_fd, progress_binding,
                f"unallocated attempt {index} campaign checkpoints")
        finally:
            if progress_fd >= 0:
                os.close(progress_fd)
        if index == next_index and staging_names and \
                checkpoint["frontier"] > 0:
            _require(checkpoint["frontier"] == 1,
                     "prepared attempt advanced beyond its ARMED checkpoint")
            prepared_attempt = _read_attempt_directory(
                attempts_fd, staging_names[0], lane=lane,
                registration=registration, lane_fd=lane_fd,
                logical_plan=logical_plan,
                candidate_authority=candidate_authority,
                exact_main=exact_main, lane_binding=lane_binding,
                require_allocation=False, read_only=True)
            _require(
                prepared_attempt[1]["evidence_attempt"] == next_index and
                not prepared_attempt[2][3],
                "prepared attempt is not the exact empty next attempt")
            prepared = (staging_names[0], prepared_attempt)
        else:
            _require(checkpoint["frontier"] == 0,
                     f"unallocated attempt {index} has a reached checkpoint")
    combined = [
        *preliminary,
        *([] if prepared is None else [prepared[1]]),
    ]
    _validate_prior_attempt_invariants(combined, reject_consumed=False)
    allocated_count = len(names) + (prepared is not None)
    for binding in lane_binding["campaign_markers"][allocated_count:]:
        expected = _validate_campaign_marker_binding(binding)
        descriptor = -1
        try:
            descriptor = os.open(
                expected["name"],
                os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
            _validate_unallocated_campaign_marker(
                lane_fd, descriptor, expected,
                "unallocated generation-3 campaign transcript slot")
        except OSError as error:
            raise ArmingError(
                "unallocated generation-3 campaign transcript slot cannot "
                "be opened") from error
        finally:
            if descriptor >= 0:
                os.close(descriptor)
    return names, preliminary, prepared


def _read_budget_authority_before_recovery(
    lane_fd: int, *, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], lane_binding: Mapping[str, Any],
    prior_attempts: Sequence[tuple[
        Mapping[str, Any], Mapping[str, Any],
        tuple[int, str | None, bool, list[tuple[str, str]]],
    ]],
    require_exact: bool = False,
) -> dict[str, Any]:
    _require(type(require_exact) is bool,
             "pre-ARMED history replay policy is invalid")
    observed: list[dict[str, Any]] = []
    unused_seen = False
    for binding_value in lane_binding["prearmed_history_markers"]:
        binding = _validate_prearmed_history_marker_binding(binding_value)
        descriptor = -1
        try:
            descriptor = os.open(
                binding["name"],
                os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
            history = _read_prearmed_history(
                lane_fd, descriptor, binding, registration=registration,
                logical_plan=logical_plan, lane_binding=lane_binding,
                label=f"pre-ARMED history {binding['name']}")
        except OSError as error:
            raise ArmingError(
                f"pre-ARMED history {binding['name']} cannot be opened") \
                from error
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        boundary = _read_prearmed_boundary_frontier(
            lane_fd, lane_binding, history["acquisition_index"])
        actual_frontier = len(history["states"]) - 1
        _require(actual_frontier <= boundary["frontier"],
                 "pre-ARMED history advanced without its durable boundary")
        if require_exact:
            _require(
                history["mode"] == 0o400 and not history["partial"] and
                actual_frontier == boundary["frontier"] and
                history["states"] == boundary["states"] and
                all(marker["mode"] == 0o400
                    for marker in boundary["markers"]),
                "live pre-ARMED history differs from its exact sealed "
                "boundary")
        history["boundary"] = boundary
        history["effective_states"] = list(boundary["states"])
        consumed = history["used"] or boundary["frontier"] > 0
        history["inode_used"] = history["used"]
        history["used"] = consumed
        if consumed:
            _require(not unused_seen,
                     "pre-ARMED histories are not one allocated prefix")
            observed.append(history)
        else:
            _require(boundary["frontier"] == 0 and
                     all(not marker["reached"] and marker["mode"] == 0o400
                         for marker in boundary["markers"]),
                     "unused pre-ARMED slot has a reached boundary")
            unused_seen = True

    links: dict[int, tuple[dict[str, Any], Mapping[str, Any], int]] = {}
    linked_indexes: list[int] = []
    for evidence_index, (bundle, armed, unused_journal) in enumerate(
            prior_attempts, 1):
        commit = _validate_budget_commit(
            bundle["budget_commit"], registration)
        acquisition_index = commit["history"]["acquisition_index"]
        _require(acquisition_index not in links and
                 armed["evidence_attempt"] == evidence_index,
                 "evidence attempts do not uniquely link budget histories")
        links[acquisition_index] = (commit, armed, evidence_index)
        linked_indexes.append(acquisition_index)
    _require(linked_indexes == sorted(linked_indexes),
             "evidence attempts reorder pre-ARMED histories")

    counters = {"setup": 0, "environment": 0, "evidence": 0}
    histories: list[dict[str, Any]] = []
    environment_states = {
        "QUALIFYING", "QUALIFIED", "BRIDGING", "BRIDGED", "PRESAMPLING",
    }
    for history in observed:
        acquisition_index = history["acquisition_index"]
        link = links.get(acquisition_index)
        states = list(history["effective_states"])
        selected_pair: Mapping[str, Any] | None = None
        if link is None:
            lane_class = (
                "environment" if states[-1] in environment_states else
                "setup")
        else:
            commit, armed, evidence_index = link
            reference = commit["history"]
            _require(history["boundary"]["frontier"] ==
                        len(PREARMED_BOUNDARY_STATES) and
                     history["states"][-1] == "PRESAMPLING" and
                     states[-1] == "PRESAMPLING" and
                     reference["history_size"] == history["complete_size"] and
                     reference["history_sha256"] ==
                        _sha256_bytes(history["complete_data"]) and
                     contract.exact_json_equal(
                         reference["marker"], history["binding"]),
                     "evidence attempt pre-ARMED history bytes differ")
            states.append("ARMED")
            lane_class = "evidence"
            selected_pair = armed["selected_pair"]
            _require(evidence_index == counters["evidence"] + 1,
                     "evidence budget history index differs")
        counters[lane_class] += 1
        history_record = plan_contract.attempt_state_history_record(
            lane_class=lane_class, lane_index=counters[lane_class],
            states=states, selected_pair=selected_pair)
        histories.append(history_record)
        prefix_ledger = plan_contract.budget_ledger_record(
            registration, histories)
        if link is not None:
            _require(contract.exact_json_equal(
                prefix_ledger, link[0]["prospective_ledger"]),
                "evidence attempt prospective budget ledger differs")
    _require(set(links) <= {
        history["acquisition_index"] for history in observed},
        "evidence attempt links an unused pre-ARMED history")
    ledger = plan_contract.budget_ledger_record(registration, histories)
    _require(ledger["evidence_attempts_used"] == len(prior_attempts),
             "durable evidence and budget-ledger counts differ")
    return {
        "histories": histories,
        "ledger": ledger,
        "observed": observed,
        "next_acquisition_index": len(observed) + 1,
    }


def _recover_prearmed_histories(
    lane_fd: int, authority: Mapping[str, Any], *,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    lane_binding: Mapping[str, Any],
) -> None:
    # The complete authority was validated above before this first recovery
    # write.  Marker frontier is the monotone authority; history is repaired
    # only by appending its unique canonical suffix to that frontier.
    for expected in authority["observed"]:
        for marker_expected in expected["boundary"]["markers"]:
            marker_binding = marker_expected["binding"]
            marker_fd = os.open(
                marker_binding["name"],
                os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=lane_fd)
            try:
                current_marker = _read_prearmed_boundary_marker(
                    lane_fd, marker_fd, marker_binding,
                    "pre-ARMED boundary before recovery")
                _require(
                    current_marker["reached"] == marker_expected["reached"] and
                    current_marker["mode"] == marker_expected["mode"] and
                    current_marker["ctime_ns"] == marker_expected["ctime_ns"],
                    "pre-ARMED boundary changed before recovery")
                if current_marker["mode"] == 0o600:
                    _require(current_marker["reached"],
                             "pristine pre-ARMED boundary became writable")
                    os.fchmod(marker_fd, 0o400)
                    os.fsync(marker_fd)
                    sealed_marker = _read_prearmed_boundary_marker(
                        lane_fd, marker_fd, marker_binding,
                        "pre-ARMED boundary after recovery")
                    _require(sealed_marker["reached"] and
                             sealed_marker["mode"] == 0o400,
                             "pre-ARMED boundary recovery did not seal")
            finally:
                os.close(marker_fd)

        binding = expected["binding"]
        descriptor = -1
        try:
            if not expected["inode_used"]:
                descriptor = _open_prearmed_history_for_allocation(
                    lane_fd, binding, registration=registration,
                    logical_plan=logical_plan, lane_binding=lane_binding)
            else:
                probe = os.open(
                    binding["name"],
                    os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
                    dir_fd=lane_fd)
                try:
                    current_probe = _read_prearmed_history(
                        lane_fd, probe, binding, registration=registration,
                        logical_plan=logical_plan, lane_binding=lane_binding,
                        label="pre-ARMED history before recovery")
                    _require(current_probe["data"] == expected["data"] and
                             current_probe["mode"] == expected["mode"],
                             "pre-ARMED history changed before recovery")
                    if current_probe["mode"] == 0o400:
                        os.fchmod(probe, 0o600)
                        os.fsync(probe)
                finally:
                    os.close(probe)
                descriptor = os.open(
                    binding["name"],
                    os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW,
                    dir_fd=lane_fd)
            current = _read_prearmed_history(
                lane_fd, descriptor, binding, registration=registration,
                logical_plan=logical_plan, lane_binding=lane_binding,
                label="recoverable pre-ARMED history")
            _require(current["data"] == expected["data"] and
                     current["complete_size"] == expected["complete_size"] and
                     current["mode"] == 0o600,
                     "pre-ARMED history changed before recovery")
            if current["partial"]:
                os.ftruncate(descriptor, current["complete_size"])
                os.fsync(descriptor)
                current = _read_prearmed_history(
                    lane_fd, descriptor, binding, registration=registration,
                    logical_plan=logical_plan, lane_binding=lane_binding,
                    label="truncated pre-ARMED history before suffix recovery")
            target_frontier = expected["boundary"]["frontier"]
            while len(current["states"]) - 1 < target_frontier:
                next_state = PREARMED_BOUNDARY_STATES[
                    len(current["states"]) - 1]
                current = _append_prearmed_state(
                    lane_fd, descriptor, binding, next_state,
                    registration=registration, logical_plan=logical_plan,
                    lane_binding=lane_binding)
            os.fchmod(descriptor, 0o400)
            os.fsync(descriptor)
            sealed = _read_prearmed_history(
                lane_fd, descriptor, binding, registration=registration,
                logical_plan=logical_plan, lane_binding=lane_binding,
                label="recovered pre-ARMED history")
            _require(sealed["used"] and sealed["mode"] == 0o400 and
                     not sealed["partial"] and
                     sealed["data"] == current["complete_data"] and
                     sealed["states"] == expected["effective_states"],
                     "pre-ARMED history recovery differs")
        finally:
            if descriptor >= 0:
                os.close(descriptor)


def _require_budget_capacity(ledger: Mapping[str, Any]) -> None:
    _require(ledger["setup_invalid_used"] < ledger["setup_invalid_limit"],
             "setup-invalid budget is exhausted")
    _require(ledger["environment_rejected_used"] <
             ledger["environment_rejected_limit"],
             "environment-rejected budget is exhausted")
    _require(ledger["evidence_attempts_used"] <
             ledger["evidence_attempts_limit"],
             "evidence-attempt budget is exhausted")


def _recover_prepared_attempt_staging(
    attempts_fd: int, staging_name: str, *, lane: Path, lane_fd: int,
    registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
    candidate_authority: Mapping[str, Any], exact_main: Mapping[str, Any],
    lane_binding: Mapping[str, Any],
    expected: tuple[
        Mapping[str, Any], Mapping[str, Any],
        tuple[int, str | None, bool, list[tuple[str, str]]],
    ],
) -> str:
    observed = _read_attempt_directory(
        attempts_fd, staging_name, lane=lane, lane_fd=lane_fd,
        registration=registration, logical_plan=logical_plan,
        candidate_authority=candidate_authority, exact_main=exact_main,
        lane_binding=lane_binding, require_allocation=False, read_only=True)
    _require(
        all(contract.exact_json_equal(first, second)
            for first, second in zip(observed[:2], expected[:2])) and
        observed[2] == expected[2] and not observed[2][3],
        "prepared attempt staging changed before recovery")
    armed = observed[1]
    attempt_index = armed["evidence_attempt"]
    marker_binding = _validate_campaign_marker_binding(
        lane_binding["campaign_markers"][attempt_index - 1])
    progress_binding = _validate_campaign_progress_directory_binding(
        lane_binding["campaign_progress_directories"][attempt_index - 1])
    progress_fd = -1
    try:
        progress_fd, observed_progress_binding = \
            _open_campaign_progress_directory(
                lane, lane_fd, progress_binding,
                "prepared attempt recovery checkpoints")
        _require(contract.exact_json_equal(
            observed_progress_binding, progress_binding),
            "prepared attempt recovery checkpoint binding changed")
        _reconcile_campaign_transcript_prefix(
            lane_fd, marker_binding, armed=armed, journal_records=[],
            label="prepared attempt campaign transcript recovery",
            progress_directory_fd=progress_fd,
            progress_directory_binding=progress_binding)
    finally:
        if progress_fd >= 0:
            os.close(progress_fd)
    final_name = f"attempt-{attempt_index:03d}"
    _require(final_name not in os.listdir(attempts_fd),
             "prepared attempt publication target appeared")
    os.fchmod(attempts_fd, 0o700)
    os.fsync(attempts_fd)
    _renameat2_noreplace(attempts_fd, staging_name, final_name)
    os.fsync(attempts_fd)
    os.fchmod(attempts_fd, 0o500)
    os.fsync(attempts_fd)
    recovered = _read_attempt_directory(
        attempts_fd, final_name, lane=lane, lane_fd=lane_fd,
        registration=registration, logical_plan=logical_plan,
        candidate_authority=candidate_authority, exact_main=exact_main,
        lane_binding=lane_binding)
    _require(
        all(contract.exact_json_equal(first, second)
            for first, second in zip(recovered[:2], expected[:2])) and
        recovered[2] == expected[2],
        "recovered prepared attempt differs")
    return final_name


def _prospective_budget_commit(
    baseline_authority: Mapping[str, Any],
    active_history: Mapping[str, Any], selected_pair: Mapping[str, Any],
    registration: Mapping[str, Any],
) -> dict[str, Any]:
    reference = _prearmed_history_reference(active_history)
    evidence_index = (
        baseline_authority["ledger"]["evidence_attempts_used"] + 1)
    history = plan_contract.attempt_state_history_record(
        lane_class="evidence", lane_index=evidence_index,
        states=(*active_history["states"], "ARMED"),
        selected_pair=selected_pair)
    prospective = plan_contract.budget_ledger_record(
        registration, [*baseline_authority["histories"], history])
    return _budget_commit_record(reference, prospective)


def _load_prior_attempts(
    lane: Path, lane_fd: int, *, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], candidate_authority: Mapping[str, Any],
    exact_main: Mapping[str, Any], lane_binding: Mapping[str, Any],
) -> list[tuple[
        dict[str, Any], dict[str, Any],
        tuple[int, str | None, bool, list[tuple[str, str]]]]]:
    try:
        attempts_fd = os.open(
            ATTEMPTS_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
    except OSError as error:
        raise ArmingError(
            "generation-3 attempts namespace cannot be replayed") from error
    recovery_started = False
    try:
        names, preliminary, prepared = _read_prior_attempts_before_recovery(
            attempts_fd, lane=lane, lane_fd=lane_fd,
            registration=registration, logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding)
        recovery_started = True
        if prepared is not None:
            staging_name, prepared_attempt = prepared
            final_name = _recover_prepared_attempt_staging(
                attempts_fd, staging_name, lane=lane, lane_fd=lane_fd,
                registration=registration, logical_plan=logical_plan,
                candidate_authority=candidate_authority,
                exact_main=exact_main, lane_binding=lane_binding,
                expected=prepared_attempt)
            names.append(final_name)
            preliminary.append(prepared_attempt)
        if stat.S_IMODE(os.fstat(attempts_fd).st_mode) != 0o700:
            os.fchmod(attempts_fd, 0o700)
            os.fsync(attempts_fd)
        _clean_staging_directories(attempts_fd)
        os.fchmod(attempts_fd, 0o500)
        os.fsync(attempts_fd)
        for index, (bundle, armed, journal_status) in enumerate(
                preliminary, 1):
            _reconcile_attempt_execution_prefixes(
                attempts_fd, names[index - 1], lane=lane, lane_fd=lane_fd,
                lane_binding=lane_binding, bundle=bundle, armed=armed,
                registration=registration, logical_plan=logical_plan,
                journal_status=journal_status)
        result = [
            _read_attempt_directory(
                attempts_fd, name, lane=lane, registration=registration,
                lane_fd=lane_fd,
                logical_plan=logical_plan,
                candidate_authority=candidate_authority,
                exact_main=exact_main, lane_binding=lane_binding)
            for name in names
        ]
        _validate_prior_attempt_invariants(result)
    finally:
        try:
            if recovery_started:
                os.fchmod(attempts_fd, 0o500)
                os.fsync(attempts_fd)
        finally:
            os.close(attempts_fd)
    return result


def _publish_attempt_directory(
    lane: Path, lane_fd: int, *, registration: Mapping[str, Any],
    logical_plan: Mapping[str, Any], candidate_authority: Mapping[str, Any],
    exact_main: Mapping[str, Any], lane_binding: Mapping[str, Any],
    authority_bundle: Mapping[str, Any],
    record_factory: Callable[[int, str | None, str, str, int], Mapping[str, Any]],
    durability_validator: Callable[[], None],
) -> tuple[Path, dict[str, Any], int, int, int, int, int, int, int, int]:
    _validate_owner_creation_umask(
        "before generation-3 attempt publication")
    prior = _load_prior_attempts(
        lane, lane_fd, registration=registration, logical_plan=logical_plan,
        candidate_authority=candidate_authority, exact_main=exact_main,
        lane_binding=lane_binding)
    attempt = len(prior) + 1
    _require(attempt <= registration["budgets"]["evidence_attempts"],
             "generation-3 evidence-attempt budget is exhausted")
    prior_sha = contract.canonical_sha256(prior[-1][1]) if prior else None
    if prior:
        expected_pair = prior[0][1]["selected_pair"]
        observed_pair = authority_bundle["qualification_binding"]["selected_pair"]
        _require(contract.exact_json_equal(observed_pair, expected_pair),
                 "new attempt changed the frozen CPU pair")
    attempts_fd = -1
    staging_fd = -1
    final_fd = journal_fd = marker_fd = campaign_marker_fd = -1
    campaign_progress_fd = -1
    armed_fd = bundle_fd = manifest_fd = -1
    published = False
    allocation_prepared = False
    transferred = False
    try:
        attempts_fd = os.open(
            ATTEMPTS_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        attempts_metadata = os.fstat(attempts_fd)
        _require(stat.S_ISDIR(attempts_metadata.st_mode) and
                 attempts_metadata.st_uid == os.geteuid() and
                 stat.S_IMODE(attempts_metadata.st_mode) == 0o500,
                 "attempt publication namespace is not sealed")
        staging_name = (
            f".attempt-{attempt:03d}.staging-{os.getpid()}-"
            f"{secrets.token_hex(16)}")
        final_name = f"attempt-{attempt:03d}"
        os.fchmod(attempts_fd, 0o700)
        os.fsync(attempts_fd)
        os.mkdir(staging_name, 0o700, dir_fd=attempts_fd)
        os.fchmod(attempts_fd, 0o500)
        os.fsync(attempts_fd)
        staging_fd = os.open(
            staging_name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=attempts_fd)
        staging_metadata = os.fstat(staging_fd)
        os.mkdir(JOURNAL_DIRECTORY, 0o700, dir_fd=staging_fd)
        journal_fd = os.open(
            JOURNAL_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=staging_fd)
        journal_binding = _journal_binding_record(os.fstat(journal_fd))
        _require(not os.listdir(journal_fd),
                 "new execution journal is not empty")
        marker_fd = os.open(
            LAUNCH_CONSUMED_FILE,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW, 0o600, dir_fd=staging_fd)
        os.fsync(marker_fd)
        os.fchmod(marker_fd, 0o400)
        os.fsync(marker_fd)
        launch_marker_binding = _launch_marker_binding_record(
            os.fstat(marker_fd))
        campaign_marker_binding = _validate_campaign_marker_binding(
            lane_binding["campaign_markers"][attempt - 1])
        campaign_marker_fd = os.open(
            campaign_marker_binding["name"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=lane_fd)
        _validate_unallocated_campaign_marker(
            lane_fd, campaign_marker_fd, campaign_marker_binding,
            "attempt campaign transcript slot")
        campaign_progress_binding = \
            _validate_campaign_progress_directory_binding(
                lane_binding["campaign_progress_directories"][attempt - 1])
        campaign_progress_fd, observed_progress_binding = \
            _open_campaign_progress_directory(
                lane, lane_fd, campaign_progress_binding,
                "attempt campaign checkpoints")
        _require(contract.exact_json_equal(
            observed_progress_binding, campaign_progress_binding),
            "attempt campaign checkpoints differ")
        _require(_read_campaign_checkpoint_frontier(
            campaign_progress_fd, campaign_progress_binding,
            "unused attempt campaign checkpoints")["frontier"] == 0,
            "attempt campaign checkpoints were already consumed")
        bundle_data = contract.canonical_json_bytes(authority_bundle)
        manifest = _attempt_manifest_record(
            bundle_data, journal_binding, launch_marker_binding,
            campaign_marker_binding, campaign_progress_binding)
        manifest_data = contract.canonical_json_bytes(manifest)
        armed_ns = time.monotonic_ns()
        armed = dict(record_factory(
            attempt, prior_sha, _sha256_bytes(bundle_data),
            _sha256_bytes(manifest_data), armed_ns))
        _validate_authority_bundle_and_armed(
            authority_bundle, armed, lane=lane, registration=registration,
            logical_plan=logical_plan, candidate_authority=candidate_authority,
            exact_main=exact_main, lane_binding=lane_binding,
            authority_bundle_data=bundle_data,
            attempt_manifest_data=manifest_data)
        durability_validator()
        _write_immutable_file_at(
            staging_fd, AUTHORITY_BUNDLE_FILE, bundle_data,
            "staged authority bundle")
        _write_immutable_file_at(
            staging_fd, ATTEMPT_MANIFEST_FILE, manifest_data,
            "staged attempt manifest")
        armed_data = contract.canonical_json_bytes(armed)
        _write_immutable_file_at(
            staging_fd, ARMED_FILE, armed_data, "staged ARMED record")
        os.fsync(journal_fd)
        os.fsync(staging_fd)
        os.fchmod(staging_fd, 0o500)
        os.fsync(staging_fd)
        durability_validator()
        _read_attempt_directory(
            attempts_fd, staging_name, lane=lane, registration=registration,
            lane_fd=lane_fd,
            logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding, require_allocation=False)
        allocation_prepared = True
        allocation_data = contract.canonical_json_bytes(
            _campaign_transcript_allocation_record(armed))
        _flip_campaign_checkpoint(
            campaign_progress_fd, campaign_progress_binding, 0,
            allocation_data, "prepared ARMED allocation")
        os.close(campaign_marker_fd)
        campaign_marker_fd = -1
        campaign_marker_fd = _open_campaign_marker_for_append(
            lane_fd, campaign_marker_binding,
            associated_with_published_attempt=True)
        _append_campaign_transcript_allocation(
            lane_fd, campaign_marker_fd, campaign_marker_binding, armed)
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane, lane_fd, registration)
        os.fchmod(attempts_fd, 0o700)
        os.fsync(attempts_fd)
        _renameat2_noreplace(attempts_fd, staging_name, final_name)
        published = True
        os.fsync(attempts_fd)
        os.fchmod(attempts_fd, 0o500)
        os.fsync(attempts_fd)
        durability_validator()
        _read_attempt_directory(
            attempts_fd, final_name, lane=lane, registration=registration,
            lane_fd=lane_fd,
            logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding, require_allocation=False)
        final_fd = os.open(
            final_name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=attempts_fd)
        final_metadata = os.fstat(final_fd)
        _require((final_metadata.st_dev, final_metadata.st_ino) ==
                 (staging_metadata.st_dev, staging_metadata.st_ino),
                 "published attempt inode differs from its staging inode")
        armed_fd, observed_armed = _read_immutable_file_at(
            final_fd, ARMED_FILE, MAX_JSON_BYTES, "live ARMED record")
        bundle_fd, observed_bundle = _read_immutable_file_at(
            final_fd, AUTHORITY_BUNDLE_FILE, MAX_JSON_BYTES,
            "live authority bundle")
        manifest_fd, observed_manifest = _read_immutable_file_at(
            final_fd, ATTEMPT_MANIFEST_FILE, MAX_JSON_BYTES,
            "live attempt manifest")
        _validate_journal_capability(
            final_fd, journal_fd, journal_binding,
            "published execution journal")
        _validate_launch_marker_capability(
            final_fd, marker_fd, launch_marker_binding,
            "published launch-consumed marker")
        _validate_campaign_marker_capability(
            lane_fd, campaign_marker_fd, campaign_marker_binding,
            "published campaign-consumed marker")
        _require(not os.listdir(journal_fd),
                 "new execution journal is not empty and owner-private")
        _require(observed_armed == armed_data and
                 observed_bundle == bundle_data and
                 observed_manifest == manifest_data,
                 "published attempt bytes changed")
        _read_attempt_directory(
            attempts_fd, final_name, lane=lane, registration=registration,
            lane_fd=lane_fd, logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding)
        durability_validator()
        os.fchmod(attempts_fd, 0o500)
        os.fsync(attempts_fd)
        _require(stat.S_IMODE(os.fstat(attempts_fd).st_mode) == 0o500,
                 "published attempts namespace did not seal")
        result = (
            lane / ATTEMPTS_DIRECTORY / final_name, copy.deepcopy(armed),
            final_fd, journal_fd, marker_fd, campaign_marker_fd,
            campaign_progress_fd,
            armed_fd, bundle_fd, manifest_fd)
        os.close(staging_fd)
        staging_fd = -1
        os.close(attempts_fd)
        attempts_fd = -1
        campaign_progress_fd = -1
        transferred = True
        return result
    finally:
        active_error = sys.exc_info()[1]
        cleanup_errors: list[BaseException] = []
        if not transferred:
            if staging_fd >= 0:
                try:
                    os.close(staging_fd)
                except BaseException as error:
                    cleanup_errors.append(error)
            for descriptor in (
                    manifest_fd, bundle_fd, armed_fd, campaign_marker_fd,
                    campaign_progress_fd, marker_fd, journal_fd, final_fd):
                if descriptor >= 0:
                    try:
                        os.close(descriptor)
                    except BaseException as error:
                        cleanup_errors.append(error)
            if attempts_fd >= 0:
                if not published and not allocation_prepared:
                    try:
                        os.fchmod(attempts_fd, 0o700)
                        os.fsync(attempts_fd)
                        _clean_staging_directories(
                            attempts_fd, remove_complete=True)
                    except BaseException as error:
                        cleanup_errors.append(error)
                try:
                    os.fchmod(attempts_fd, 0o500)
                except BaseException as error:
                    cleanup_errors.append(error)
                try:
                    os.fsync(attempts_fd)
                except BaseException as error:
                    cleanup_errors.append(error)
                try:
                    os.close(attempts_fd)
                except BaseException as error:
                    cleanup_errors.append(error)
            if cleanup_errors and active_error is None:
                raise ArmingError(
                    "attempt publication capability cleanup failed") \
                    from cleanup_errors[0]


def _acquire_track_a_qualification(
    preregistration: Mapping[str, Any], *,
    expected_frozen_pair: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Acquire one production Track-A chain with no injectable host reader."""
    registration = prereg.validate_preregistration(
        preregistration, verify_files=True)
    qualification = registration["qualification"]
    _require(qualification["track_b_permitted"] is False and
             qualification["policy_evaluation_order"] == ["pair-and-domain"] and
             len(qualification["policies"]) == 1,
             "live armer currently requires preregistered Track A only")
    host_reader = pair_acquire.SystemHostReader()
    policy = qualification["policies"][0]
    geometry = qualification["geometry"]
    allowed_before = sorted(host_reader.allowed_cpus())
    topology_before = pair_acquire._acquire_topology(
        host_reader, allowed_before, policy["candidate_primary_cpus"])
    host_before = current_host_facts(
        topology_sha256=contract.canonical_sha256(topology_before),
        allowed_cpus=allowed_before)
    _require(contract.exact_json_equal(
        _prereg_host_projection(host_before["authority"]),
        registration["host_authority"]),
        "live host differs from the ratified host authority")
    acquisition = pair_acquire.acquire_pair_qualification(
        host_reader, policy_value=policy,
        window_count=geometry["scan_window_count"],
        nominal_window_ns=geometry["scan_nominal_window_ns"],
        frozen_pair_from_prior_attempt=expected_frozen_pair)
    selection = prereg.qualification_track_selection_record(
        registration, policy_a_scan=acquisition["scan"], policy_b_scan=None,
        expected_frozen_pair=expected_frozen_pair)
    _require(selection["selected_track"] == "pair-and-domain" and
             selection["selected_pair"] is not None,
             "Track A did not select one qualified pair")
    bridge_geometry = {
        "minimum_window_count": geometry["bridge_minimum_window_count"],
        "maximum_window_count": geometry["bridge_maximum_window_count"],
        "nominal_window_ns": geometry["bridge_nominal_window_ns"],
        "maximum_handoff_elapsed_ns": geometry["maximum_handoff_elapsed_ns"],
    }
    scan_tail = acquisition["scan"]["windows"][-1]["after"]
    observed_cpus = sorted(int(cpu) for cpu in scan_tail["cpus"])
    snapshots = [copy.deepcopy(scan_tail)]
    for unused_index in range(bridge_geometry["minimum_window_count"]):
        host_reader.sleep_ns(bridge_geometry["nominal_window_ns"])
        snapshots.append(pair_acquire._capture_snapshot(
            host_reader, observed_cpus))
    bridge = pair_bridge.pair_qualification_bridge_record(
        acquisition, expected_policy=policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_acquisition_window_count=geometry["scan_window_count"],
        expected_acquisition_nominal_window_ns=
            geometry["scan_nominal_window_ns"],
        expected_bridge_geometry=bridge_geometry, snapshots=snapshots)
    acquisition_data = contract.canonical_json_bytes(acquisition)
    bridge_data = contract.canonical_json_bytes(bridge)
    verdict = pair_verify.require_accepted_pair_qualification_bundle(
        acquisition_data, bridge_data, expected_policy=policy,
        expected_policy_sha256=qualification["policy_sha256s"][0],
        expected_frozen_pair=expected_frozen_pair,
        expected_acquisition_window_count=geometry["scan_window_count"],
        expected_acquisition_nominal_window_ns=
            geometry["scan_nominal_window_ns"],
        expected_bridge_geometry=bridge_geometry)
    verdict_data = contract.canonical_json_bytes(verdict)
    host_after = current_host_facts(
        topology_sha256=acquisition["topology_after_sha256"],
        allowed_cpus=acquisition["allowed_cpu_set_after_scan"])
    _require(contract.exact_json_equal(host_before, host_after),
             "host instance changed across qualification and bridge")
    evidence = {
        "policy_a_scan": copy.deepcopy(acquisition["scan"]),
        "policy_b_scan": None,
        "track_selection": selection,
        "expected_frozen_pair": copy.deepcopy(expected_frozen_pair),
        "acquisition_data": acquisition_data,
        "bridge_data": bridge_data,
        "independent_verdict_data": verdict_data,
    }
    binding = execution.qualification_binding_record(
        registration, host_after["instance"], evidence=evidence)
    return {
        "host_instance": host_after["instance"],
        "qualification_evidence": evidence,
        "qualification_binding": binding,
    }


PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37


def _child_subreaper_state() -> bool:
    value = ctypes.c_int()
    libc = ctypes.CDLL(None, use_errno=True)
    result = libc.prctl(
        ctypes.c_int(PR_GET_CHILD_SUBREAPER), ctypes.byref(value),
        ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        error = ctypes.get_errno()
        raise ArmingError("child-subreaper state could not be read") \
            from OSError(error, os.strerror(error))
    return value.value != 0


def _set_child_subreaper(enabled: bool) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    result = libc.prctl(
        ctypes.c_int(PR_SET_CHILD_SUBREAPER), ctypes.c_ulong(int(enabled)),
        ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        error = ctypes.get_errno()
        raise ArmingError("child-subreaper state could not be changed") \
            from OSError(error, os.strerror(error))
    _require(_child_subreaper_state() is enabled,
             "child-subreaper state did not reach its requested value")


def _direct_child_pids() -> list[int]:
    tasks = os.listdir("/proc/self/task")
    _require(tasks == [str(os.getpid())] or
             (len(tasks) == 1 and tasks[0] == str(os.getpid())),
             "child supervisor requires a single-threaded authority process")
    descriptor = os.open(
        f"/proc/self/task/{os.getpid()}/children",
        os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        data = os.read(descriptor, 1 << 20)
        _require(not os.read(descriptor, 1),
                 "direct-child inventory exceeded its byte bound")
    finally:
        os.close(descriptor)
    try:
        children = [int(item) for item in data.decode("ascii").split()]
    except (UnicodeDecodeError, ValueError) as error:
        raise ArmingError("direct-child inventory is malformed") from error
    _require(all(child > 0 for child in children) and
             len(children) == len(set(children)),
             "direct-child inventory is not one unique PID set")
    return sorted(children)


class _ChildTreeScope:
    """Make descendants observable and freeze parent signal semantics."""

    def __init__(self, launch_context: Mapping[str, Any]) -> None:
        self._launch_context = _validate_launch_context(launch_context)
        self._previous = False
        self._previous_mask: set[signal.Signals] = set()
        self._previous_pending: set[signal.Signals] = set()
        self._mask_active = False
        self._entered = False

    @staticmethod
    def _blockable_signals() -> set[signal.Signals]:
        return {
            member for member in signal.valid_signals()
            if member not in {signal.SIGKILL, signal.SIGSTOP}
        }

    def child_preexec(self) -> None:
        """Install the exact preregistered child state before any exec."""
        context = self._launch_context
        blockable = self._blockable_signals()
        for member in sorted(blockable, key=int):
            signal.signal(member, signal.SIG_DFL)
        signal.pthread_sigmask(signal.SIG_SETMASK, set())
        for limit in context["rlimits"]:
            resource.setrlimit(
                getattr(resource, limit["name"]),
                (limit["soft"], limit["hard"]))
        os.sched_setscheduler(
            0, context["scheduler_policy"],
            os.sched_param(context["scheduler_priority"]))
        os.setpriority(os.PRIO_PROCESS, 0, context["nice"])
        os.umask(context["umask"])

        _require(not signal.pthread_sigmask(signal.SIG_BLOCK, set()) and
                 all(signal.getsignal(member) is signal.SIG_DFL
                     for member in blockable),
                 "child signal context did not reach its fixed point")
        _require(os.getpriority(os.PRIO_PROCESS, 0) == context["nice"] and
                 os.sched_getscheduler(0) == context["scheduler_policy"] and
                 os.sched_getparam(0).sched_priority ==
                    context["scheduler_priority"] and
                 all(resource.getrlimit(getattr(resource, limit["name"])) ==
                     (limit["soft"], limit["hard"])
                     for limit in context["rlimits"]),
                 "child process context did not reach its fixed point")
        observed_umask = os.umask(context["umask"])
        _require(observed_umask == context["umask"],
                 "child umask did not reach its fixed point")

    def _restore_signal_authority(self) -> None:
        if not self._mask_active:
            return
        restore_error: BaseException | None = None
        try:
            _require(signal.getsignal(signal.SIGCHLD) is signal.SIG_DFL,
                     "SIGCHLD changed while child signals were blocked")
        except BaseException as error:
            restore_error = error
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, self._previous_mask)
        except BaseException as error:
            if restore_error is None:
                restore_error = error
        finally:
            self._mask_active = False
        try:
            if signal.getsignal(signal.SIGCHLD) is not signal.SIG_DFL:
                signal.signal(signal.SIGCHLD, signal.SIG_DFL)
                _fail("SIGCHLD changed while restoring child signal authority")
        except BaseException as error:
            if restore_error is None:
                restore_error = error
        if restore_error is not None:
            raise ArmingError("child signal authority changed") \
                from restore_error

    def __enter__(self) -> "_ChildTreeScope":
        _require(signal.getsignal(signal.SIGCHLD) is signal.SIG_DFL,
                 "child supervisor requires the default SIGCHLD disposition")
        _require(not _direct_child_pids(),
                 "child supervisor began with an existing child process")
        self._previous = _child_subreaper_state()
        try:
            self._previous_mask = set(signal.pthread_sigmask(
                signal.SIG_BLOCK, self._blockable_signals()))
            self._mask_active = True
            self._previous_pending = set(signal.sigpending())
            _require(signal.getsignal(signal.SIGCHLD) is signal.SIG_DFL,
                     "SIGCHLD changed while entering child authority")
            if not self._previous:
                _set_child_subreaper(True)
            self._entered = True
            _require(not _direct_child_pids(),
                     "child supervisor raced with another child process")
            return self
        except BaseException as entry_error:
            cleanup_error: BaseException | None = None
            try:
                if not self._previous and _child_subreaper_state():
                    _set_child_subreaper(False)
            except BaseException as error:
                cleanup_error = error
            finally:
                try:
                    self._restore_signal_authority()
                except BaseException as error:
                    if cleanup_error is None:
                        cleanup_error = error
                finally:
                    self._entered = False
            if cleanup_error is not None:
                raise ArmingError(
                    "child authority entry rollback failed") from cleanup_error
            raise entry_error

    def close(self) -> None:
        if not self._entered:
            return
        cleanup_error: BaseException | None = None
        try:
            _require(not _direct_child_pids(),
                     "child supervisor retained an unreaped child process")
            unexpected = set(signal.sigpending()) - self._previous_pending - {
                signal.SIGCHLD}
            _require(not unexpected,
                     "child signal authority observed an unexpected parent signal")
        except BaseException as error:
            cleanup_error = error
        finally:
            try:
                if not self._previous:
                    _set_child_subreaper(False)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
            finally:
                try:
                    self._restore_signal_authority()
                except BaseException as error:
                    if cleanup_error is None:
                        cleanup_error = error
                finally:
                    self._entered = False
        if cleanup_error is not None:
            raise ArmingError(
                f"child authority scope cleanup failed: {cleanup_error}") \
                from cleanup_error


def _wait_pidfd_exit(pidfd: int, deadline: float) -> bool:
    watcher = selectors.DefaultSelector()
    try:
        watcher.register(pidfd, selectors.EVENT_READ)
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return False
            if watcher.select(min(0.25, remaining)):
                return True
    finally:
        watcher.close()


def _reap_residual_descendants(leader_pid: int, deadline: float) -> bool:
    observed = False
    while True:
        residual = [pid for pid in _direct_child_pids() if pid != leader_pid]
        if not residual:
            return observed
        observed = True
        for pid in residual:
            pidfd = -1
            try:
                pidfd = os.pidfd_open(pid, 0)
                try:
                    signal.pidfd_send_signal(pidfd, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _require(_wait_pidfd_exit(pidfd, deadline),
                         "escaped child did not exit after SIGKILL")
                waited, unused_status = os.waitpid(pid, 0)
                _require(waited == pid,
                         "escaped child could not be reaped exactly")
            except ChildProcessError as error:
                raise ArmingError(
                    "escaped child left the subreaper authority") from error
            finally:
                if pidfd >= 0:
                    os.close(pidfd)


def _kill_process_group(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def _run_child_process(
    command: Sequence[str], *, descriptor: int, timeout_seconds: int,
    expected_launch_context: Mapping[str, Any],
) -> dict[str, Any]:
    """Run one fixed child while bounding both captured streams."""
    launch_context = _validate_launch_context(expected_launch_context)
    _validate_launch_context_current(
        launch_context, "before child supervisor")
    started_ns = time.monotonic_ns()
    tree_scope = _ChildTreeScope(launch_context)
    process: subprocess.Popen[bytes] | None = None
    pidfd = -1
    selector: selectors.BaseSelector | None = None
    stdout_fd = stderr_fd = -1
    streams: dict[int, tuple[Any, int]] = {}
    buffers: dict[int, bytearray] = {}
    open_streams: set[int] = set()
    outcome: str | None = None
    error_text: str | None = None
    returncode: int | None = None
    deadline = time.monotonic() + timeout_seconds
    termination_deadline: float | None = None
    leader_exited = False
    group_terminated = False
    leader_reaped = False
    residual_descendant = False

    def terminate_group() -> None:
        nonlocal group_terminated, termination_deadline
        _require(process is not None,
                 "child process does not exist for group termination")
        if not group_terminated:
            _kill_process_group(process)
            group_terminated = True
        if termination_deadline is None:
            termination_deadline = time.monotonic() + 5.0

    def handle_leader_event() -> None:
        nonlocal leader_exited
        _require(selector is not None,
                 "child selector does not exist for leader event")
        if leader_exited:
            return
        leader_exited = True
        try:
            selector.unregister(pidfd)
        except KeyError:
            pass
        # The pidfd is readable before wait() reaps the leader, so its numeric
        # PID/PGID cannot yet be reused.  Terminate any descendants now; never
        # issue a numeric kill after the successful reap.
        terminate_group()

    tree_scope.__enter__()
    try:
        controlled_fork_permit = _arm_controlled_process_fork(descriptor)
        try:
            try:
                process = subprocess.Popen(
                    list(command), env=dict(_frozen_child_environment()),
                    cwd="/", stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    close_fds=True, pass_fds=(descriptor,),
                    start_new_session=True,
                    preexec_fn=tree_scope.child_preexec)
            finally:
                _clear_controlled_process_fork(controlled_fork_permit)
        except OSError as error:
            return {
                "outcome": "output-rejected", "returncode": None,
                "stdout": b"", "stderr": b"",
                "elapsed_ns": time.monotonic_ns() - started_ns,
                "error": f"child process creation failed: {error}"[:4096],
            }
        selector = selectors.DefaultSelector()
        _require(process.stdout is not None and process.stderr is not None,
                 "child pipes were not created")
        pidfd = os.pidfd_open(process.pid, 0)
        stdout_fd = process.stdout.fileno()
        stderr_fd = process.stderr.fileno()
        streams = {
            stdout_fd: (process.stdout, MAX_CHILD_STDOUT_BYTES),
            stderr_fd: (process.stderr, MAX_CHILD_STDERR_BYTES),
        }
        buffers = {descriptor_number: bytearray()
                   for descriptor_number in streams}
        open_streams = set(streams)
        for descriptor_number in streams:
            os.set_blocking(descriptor_number, False)
            selector.register(descriptor_number, selectors.EVENT_READ)
        selector.register(pidfd, selectors.EVENT_READ)
        while open_streams:
            now = time.monotonic()
            if now >= deadline and termination_deadline is None:
                if outcome is None:
                    outcome = "timeout"
                    error_text = f"child exceeded {timeout_seconds} seconds"
                terminate_group()
            if (termination_deadline is not None and
                    now >= termination_deadline):
                raise ArmingError(
                    "child pipes remained open after process-group termination")
            wait_deadline = (termination_deadline
                             if termination_deadline is not None else deadline)
            wait_seconds = max(0.001, min(0.25, wait_deadline - now))
            events = selector.select(wait_seconds)
            for key, unused_mask in events:
                descriptor_number = int(key.fd)
                if descriptor_number == pidfd:
                    handle_leader_event()
                    continue
                try:
                    block = os.read(descriptor_number, 65536)
                except BlockingIOError:
                    continue
                if not block:
                    selector.unregister(descriptor_number)
                    open_streams.remove(descriptor_number)
                    continue
                buffer = buffers[descriptor_number]
                limit = streams[descriptor_number][1]
                if len(buffer) + len(block) > limit:
                    retained = max(0, limit - len(buffer))
                    buffer.extend(block[:retained])
                    if outcome is None:
                        outcome = "output-limit"
                        error_text = "child output exceeded its frozen byte limit"
                        terminate_group()
                else:
                    buffer.extend(block)
        while not leader_exited:
            now = time.monotonic()
            wait_deadline = (termination_deadline
                             if termination_deadline is not None else deadline)
            if now >= wait_deadline:
                if termination_deadline is not None:
                    raise ArmingError(
                        "child leader did not exit after process-group termination")
                if outcome is None:
                    outcome = "timeout"
                    error_text = f"child exceeded {timeout_seconds} seconds"
                terminate_group()
                continue
            for key, unused_mask in selector.select(
                    max(0.001, min(0.25, wait_deadline - now))):
                if int(key.fd) == pidfd:
                    handle_leader_event()
        _require(group_terminated,
                 "child process group was not terminated before reap")
        residual_descendant = _reap_residual_descendants(
            process.pid, time.monotonic() + 5.0)
        returncode = process.wait(timeout=5)
        leader_reaped = True
    finally:
        # Any exceptional post-Popen path still terminates the fresh session
        # before reaping.  Once wait() succeeds, no numeric kill is permitted:
        # that PID/PGID may immediately be reused by an unrelated process.
        cleanup_error: BaseException | None = None
        if process is not None and not group_terminated:
            try:
                terminate_group()
            except BaseException as error:
                cleanup_error = error
        if process is not None and pidfd >= 0 and not leader_exited:
            try:
                leader_exited = _wait_pidfd_exit(
                    pidfd, time.monotonic() + 5.0)
                _require(leader_exited,
                         "child leader did not exit during cleanup")
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if process is not None and leader_exited:
            try:
                residual_descendant = (
                    _reap_residual_descendants(
                        process.pid, time.monotonic() + 5.0) or
                    residual_descendant)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if process is not None and not leader_reaped:
            try:
                process.wait(timeout=5)
                leader_reaped = True
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if process is not None and leader_reaped:
            try:
                residual_descendant = (
                    _reap_residual_descendants(
                        process.pid, time.monotonic() + 5.0) or
                    residual_descendant)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if selector is not None:
            try:
                selector.close()
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if process is not None and process.stdout is not None:
            try:
                process.stdout.close()
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if process is not None and process.stderr is not None:
            try:
                process.stderr.close()
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if pidfd >= 0:
            try:
                os.close(pidfd)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        try:
            tree_scope.close()
        except BaseException as error:
            if cleanup_error is None:
                cleanup_error = error
        try:
            _validate_launch_context_current(
                launch_context, "after child supervisor")
        except BaseException as error:
            if cleanup_error is None:
                cleanup_error = error
        if cleanup_error is not None:
            raise ArmingError(
                f"child process tree could not be contained: {cleanup_error}") \
                from cleanup_error
    elapsed_ns = time.monotonic_ns() - started_ns
    stdout = bytes(buffers[stdout_fd])
    stderr = bytes(buffers[stderr_fd])
    if residual_descendant and outcome is None:
        outcome = "output-rejected"
        error_text = "child left a residual descendant process"
    if outcome is None and returncode != 0:
        outcome = "nonzero"
        error_text = (
            f"child exited {returncode}: " +
            stderr.decode("utf-8", errors="replace")[-1000:])[:4096]
    if outcome is None:
        outcome = "accepted"
    return {
        "outcome": outcome, "returncode": returncode,
        "stdout": stdout, "stderr": stderr, "elapsed_ns": elapsed_ns,
        "error": error_text,
    }


def _actual_child_command(
    child: Mapping[str, Any], selected_pair: Mapping[str, Any], descriptor: int,
) -> list[str]:
    return [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(selected_pair["benchmark_cpu"]),
        f"/proc/self/fd/{descriptor}", *child["argv_tail"],
    ]


@_with_exact_source_imports
def acquire_and_arm(
    *, lane: Path, preregistration_bytes: bytes,
    candidate_authority_lane: Path, exact_main_authority_lane: Path,
) -> "ArmedCampaign":
    """Create the sole live Gen3 authority and retain every capability.

    The signature intentionally offers no reader, evidence, pair, executable,
    attempt, launcher, fault, or subprocess injection surface.
    """
    _require(type(preregistration_bytes) is bytes and
             0 < len(preregistration_bytes) <= MAX_JSON_BYTES,
             "preregistration bytes are absent or oversized")
    value = contract.strict_json_loads(
        preregistration_bytes, "generation-3 preregistration")
    _require(type(value) is dict and
             preregistration_bytes == contract.canonical_json_bytes(value),
             "preregistration bytes are not canonical")

    canonical_lock: Any | None = None
    canonical_held = False
    lane_fd = lane_lock_fd = candidate_fd = main_fd = -1
    attempts_fd = attempt_fd = journal_fd = marker_fd = -1
    armed_fd = bundle_fd = manifest_fd = -1
    lane_binding_fd = campaign_marker_fd = -1
    campaign_progress_fd = -1
    prearmed_history_fd = -1
    runtime_guard: _RetainedRuntimeGuard | None = None
    transferred = False
    try:
        prelock_registration = prereg.validate_preregistration(
            value, verify_files=False)
        lane_path, lane_fd = \
            prereg.open_output_lane_identity_for_preregistration(
                lane, prelock_registration)
        _validate_owner_creation_umask(
            "before canonical build/test lock acquisition")
        canonical_lock = baseline_acquire.CanonicalFileLock(
            AUTHORITATIVE_LOCK_PATH)
        try:
            canonical_lock.__enter__()
            canonical_held = True
        except Exception as error:
            raise ArmingError("canonical build/test lock could not be acquired") \
                from error
        registration = prereg.validate_preregistration(
            value, verify_files=True)
        _require(contract.exact_json_equal(
            registration, prelock_registration),
            "preregistration changed across canonical-lock acquisition")
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            lane_path, lane_fd, registration)
        _validate_controller_bindings_current(registration)
        logical_plan = plan_contract.campaign_plan_record(
            preregistration=registration)
        candidate_authority = bind_candidate_authority_lane(
            candidate_authority_lane, registration)
        exact_main = bind_exact_main_path_variant(exact_main_authority_lane)
        lane_path = _validate_lane_authority_separation(
            lane_path, candidate_authority, exact_main)
        runtime_guard = _acquire_runtime_guard(
            candidate_authority["runtime_closure"], exact_main["record"])
        candidate_fd, candidate_identity = capture_sealed_executable_bytes(
            candidate_authority["artifact_data"],
            expected_sha256=registration["candidate_executable"]["sha256"],
            expected_size=registration["candidate_executable"]["size"],
            label="candidate-control",
            authority_relative_path="artifacts/bench_leopard2")
        main_fd, main_identity = capture_sealed_executable_bytes(
            exact_main["artifact_data"], expected_sha256=PATH_VARIANT_RAW_SHA256,
            expected_size=PATH_VARIANT_SIZE, label="exact-main-path-variant",
            authority_relative_path=PATH_VARIANT_RELATIVE)
        retained_lane_fd = lane_fd
        lane_fd = -1
        lane_path, lane_fd, lane_lock_fd = _open_lane_and_lock(
            lane_path, retained_lane_fd, registration,
            evidence_attempt_limit=registration["budgets"][
                "evidence_attempts"],
            prearmed_attempt_limit=_prearmed_history_limit(registration))
        lane_binding = _initialize_lane_locked(
            lane_path, lane_fd, lane_lock_fd, registration, logical_plan,
            candidate_authority, exact_main)
        prior = _load_prior_attempts(
            lane_path, lane_fd, registration=registration,
            logical_plan=logical_plan, candidate_authority=candidate_authority,
            exact_main=exact_main, lane_binding=lane_binding)
        budget_authority = _read_budget_authority_before_recovery(
            lane_fd, registration=registration, logical_plan=logical_plan,
            lane_binding=lane_binding, prior_attempts=prior)
        _require_budget_capacity(budget_authority["ledger"])
        acquisition_index = budget_authority["next_acquisition_index"]
        _require(acquisition_index <=
                 len(lane_binding["prearmed_history_markers"]),
                 "generation-3 pre-ARMED history slots are exhausted")
        prearmed_binding = lane_binding["prearmed_history_markers"] \
            [acquisition_index - 1]
        prearmed_boundaries = _prearmed_boundary_bindings_for_acquisition(
            lane_binding, acquisition_index)
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[0], "PREREGISTERED")
        prearmed_history_fd = _open_prearmed_history_for_allocation(
            lane_fd, prearmed_binding, registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "PREREGISTERED", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        launch_context = _validate_launch_context(
            registration["launch_context"])
        _validate_launch_context_current(
            launch_context, "before generation-3 qualification")
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[1], "QUALIFYING")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "QUALIFYING", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        expected_pair = prior[0][1]["selected_pair"] if prior else None
        if prior:
            _require(contract.exact_json_equal(
                _capture_live_host_instance(registration),
                prior[0][0]["host_instance"]),
                "generation-3 retry changed the full frozen host instance")
            _require(contract.exact_json_equal(
                launch_context, prior[0][0]["launch_context"]),
                "generation-3 retry changed the launch context")

        qualification = _acquire_track_a_qualification(
            registration, expected_frozen_pair=expected_pair)
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[2], "QUALIFIED")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "QUALIFIED", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        host = execution.validate_host_instance(qualification["host_instance"])
        evidence = qualification["qualification_evidence"]
        binding = execution.validate_qualification_binding(
            qualification["qualification_binding"], registration, host,
            evidence=evidence)
        _require(contract.exact_json_equal(
            evidence["expected_frozen_pair"], expected_pair),
            "fresh qualification did not use the replay-derived frozen pair")
        if prior:
            _require(contract.exact_json_equal(
                host, prior[0][0]["host_instance"]),
                "generation-3 retry changed the full frozen host instance")
        _require(contract.exact_json_equal(
            binding["selected_pair"],
            expected_pair if expected_pair is not None
            else binding["selected_pair"]),
            "generation-3 retry changed the frozen pair")
        for field in (
                "acquisition_data", "bridge_data",
                "independent_verdict_data"):
            prereg.reject_denylisted_evidence_hash(
                _sha256_bytes(evidence[field]), registration)
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[3], "BRIDGING")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "BRIDGING", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        bridge = contract.strict_json_loads(
            evidence["bridge_data"], "fresh qualification bridge")
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[4], "BRIDGED")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "BRIDGED", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[5], "ARMING")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "ARMING", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)

        candidate_snapshot = candidate_identity["snapshot"]
        main_snapshot = main_identity["snapshot"]
        candidate_handle = contract.canonical_sha256(candidate_snapshot)
        main_handle = contract.canonical_sha256(main_snapshot)
        roles = {
            "candidate": {
                "role": "candidate",
                "raw_sha256": candidate_snapshot["raw_sha256"],
                "size": candidate_snapshot["size"],
                "launch_protocol": execution.CANDIDATE_LAUNCH_PROTOCOL,
                "handle_id": candidate_handle,
                "handle_device": candidate_snapshot["device"],
                "handle_inode": candidate_snapshot["inode"],
            },
            "control": {
                "role": "control",
                "raw_sha256": candidate_snapshot["raw_sha256"],
                "size": candidate_snapshot["size"],
                "launch_protocol": execution.CONTROL_LAUNCH_PROTOCOL,
                "handle_id": candidate_handle,
                "handle_device": candidate_snapshot["device"],
                "handle_inode": candidate_snapshot["inode"],
            },
            "main": {
                "role": "main", "raw_sha256": main_snapshot["raw_sha256"],
                "size": main_snapshot["size"],
                "launch_protocol": execution.EXACT_MAIN_LAUNCH_PROTOCOL,
                "handle_id": main_handle,
                "handle_device": main_snapshot["device"],
                "handle_inode": main_snapshot["inode"],
            },
        }
        exact_main_authority = execution.exact_main_authority_record(
            verifier_verdict_sha256=exact_main["verdict"]["verdict_sha256"])
        artifacts = execution.artifact_bundle_record(
            registration, roles=roles,
            exact_main_authority=exact_main_authority,
            candidate_source_authority=candidate_authority["source_authority"],
            host_authority=host["authority"])
        _flip_prearmed_boundary_marker(
            lane_fd, prearmed_boundaries[6], "PRESAMPLING")
        active_history = _append_prearmed_state(
            lane_fd, prearmed_history_fd, prearmed_binding,
            "PRESAMPLING", registration=registration,
            logical_plan=logical_plan, lane_binding=lane_binding)
        budget_commit = _prospective_budget_commit(
            budget_authority, active_history, binding["selected_pair"],
            registration)
        authority_bundle = _authority_bundle_record(
            lane_binding=lane_binding, preregistration=registration,
            logical_plan=logical_plan, host_instance=host,
            launch_context=launch_context,
            candidate_authority_record=candidate_authority["record"],
            exact_main_authority_record=exact_main["record"],
            exact_main_verifier=exact_main["verdict"],
            artifact_bundle=artifacts, qualification_binding=binding,
            qualification_evidence=evidence,
            candidate_descriptor=candidate_snapshot,
            main_descriptor=main_snapshot, budget_commit=budget_commit)

        def durability_validator() -> None:
            now = time.monotonic_ns()
            _require(bridge["bridge_finished_monotonic_ns"] <= now <=
                     bridge["bridge_deadline_monotonic_ns"],
                     "qualification-to-ARMED handoff deadline expired")
            try:
                canonical_lock.validate_current()
            except Exception as error:
                raise ArmingError("canonical build/test lock changed") from error
            _validate_lane_capability(
                lane_path, lane_fd, lane_lock_fd, registration)
            retained_history = _read_prearmed_history(
                lane_fd, prearmed_history_fd, prearmed_binding,
                registration=registration, logical_plan=logical_plan,
                lane_binding=lane_binding,
                label="durable pre-ARMED evidence history")
            retained_frontier = _read_prearmed_boundary_frontier(
                lane_fd, lane_binding, acquisition_index)
            _require(
                retained_frontier["frontier"] ==
                    len(PREARMED_BOUNDARY_STATES) and
                all(marker["mode"] == 0o400
                    for marker in retained_frontier["markers"]) and
                retained_history["states"][-1] == "PRESAMPLING" and
                not retained_history["partial"] and
                budget_commit["history"]["history_sha256"] ==
                    _sha256_bytes(retained_history["data"]) and
                budget_commit["history"]["history_size"] ==
                    len(retained_history["data"]),
                "durable pre-ARMED evidence history changed")
            _require(contract.exact_json_equal(
                _capture_live_host_instance(registration), host),
                "live host changed before durable ARMED")
            _validate_launch_context_current(
                launch_context, "before durable ARMED")
            _validate_controller_bindings_current(registration)
            runtime_guard.validate_current("durable ARMED runtime closure")
            revalidate_sealed_descriptor(
                candidate_fd, candidate_snapshot, "candidate/control")
            revalidate_sealed_descriptor(
                main_fd, main_snapshot, "exact-main")

        def record_factory(
            attempt: int, prior_sha: str | None, bundle_sha: str,
            manifest_sha: str, armed_ns: int,
        ) -> Mapping[str, Any]:
            return execution.armed_record(
                registration, logical_plan, host, artifacts, binding,
                qualification_evidence=evidence, evidence_attempt=attempt,
                prior_armed_chain_sha256=prior_sha,
                authority_bundle_sha256=bundle_sha,
                attempt_manifest_sha256=manifest_sha,
                lane_binding_sha256=contract.canonical_sha256(lane_binding),
                armed_monotonic_ns=armed_ns)

        (attempt_path, armed, attempt_fd, journal_fd, marker_fd,
         campaign_marker_fd, campaign_progress_fd, armed_fd, bundle_fd,
         manifest_fd) = _publish_attempt_directory(
            lane_path, lane_fd, registration=registration,
            logical_plan=logical_plan, candidate_authority=candidate_authority,
            exact_main=exact_main, lane_binding=lane_binding,
            authority_bundle=authority_bundle, record_factory=record_factory,
            durability_validator=durability_validator)
        sealed_history = _seal_prearmed_history(
            lane_fd, prearmed_history_fd, prearmed_binding,
            registration=registration, logical_plan=logical_plan,
            lane_binding=lane_binding)
        _require(
            budget_commit["history"]["history_sha256"] ==
                _sha256_bytes(sealed_history["data"]) and
            budget_commit["history"]["history_size"] ==
                len(sealed_history["data"]),
            "sealed pre-ARMED evidence history differs from its commit")
        attempts_fd = os.open(
            ATTEMPTS_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=lane_fd)
        attempts_binding = _attempts_binding_record(os.fstat(attempts_fd))
        _validate_attempts_capability(
            lane_fd, attempts_fd, attempts_binding,
            attempt_fd, attempt_path.name)
        lane_binding_fd, lane_binding_data = _read_immutable_file_at(
            lane_fd, LANE_BINDING_FILE, MAX_JSON_BYTES,
            "retained generation-3 lane binding")
        _require(lane_binding_data == contract.canonical_json_bytes(lane_binding),
                 "retained lane binding changed after publication")
        durability_validator()
        os.close(prearmed_history_fd)
        prearmed_history_fd = -1
        campaign = ArmedCampaign(
            _token=_ARMED_CAMPAIGN_TOKEN, attempt_path=attempt_path,
            record=armed, canonical_lock=canonical_lock,
            runtime_guard=runtime_guard, lane=lane_path,
            lane_fd=lane_fd, lane_lock_fd=lane_lock_fd,
            lane_binding_fd=lane_binding_fd,
            campaign_marker_fd=campaign_marker_fd,
            campaign_progress_fd=campaign_progress_fd,
            candidate_control_fd=candidate_fd, main_fd=main_fd,
            attempts_fd=attempts_fd, attempt_fd=attempt_fd,
            journal_fd=journal_fd,
            launch_marker_fd=marker_fd,
            armed_fd=armed_fd, bundle_fd=bundle_fd, manifest_fd=manifest_fd,
            registration=registration, logical_plan=logical_plan,
            candidate_authority=candidate_authority, exact_main=exact_main,
            lane_binding=lane_binding, authority_bundle=authority_bundle,
            candidate_snapshot=candidate_snapshot, main_snapshot=main_snapshot)
        campaign._revalidate_authority()
        transferred = True
        canonical_held = False
        canonical_lock = None
        runtime_guard = None
        lane_fd = lane_lock_fd = lane_binding_fd = campaign_marker_fd = -1
        campaign_progress_fd = -1
        candidate_fd = main_fd = -1
        attempts_fd = attempt_fd = journal_fd = marker_fd = -1
        armed_fd = bundle_fd = manifest_fd = -1
        return campaign
    finally:
        if not transferred:
            active_error = sys.exc_info()[1]
            cleanup_errors: list[BaseException] = []
            if prearmed_history_fd >= 0:
                try:
                    cleanup_history = _read_prearmed_history(
                        lane_fd, prearmed_history_fd, prearmed_binding,
                        registration=registration,
                        logical_plan=logical_plan,
                        lane_binding=lane_binding,
                        label="failed acquisition pre-ARMED history")
                    if (cleanup_history["mode"] == 0o600 and
                            cleanup_history["complete_size"] > 0 and
                            not cleanup_history["partial"]):
                        _seal_prearmed_history(
                            lane_fd, prearmed_history_fd, prearmed_binding,
                            registration=registration,
                            logical_plan=logical_plan,
                            lane_binding=lane_binding)
                except BaseException as error:
                    if active_error is None:
                        cleanup_errors.append(error)
            for descriptor in (
                    manifest_fd, bundle_fd, armed_fd, marker_fd, journal_fd,
                    attempt_fd, attempts_fd, campaign_marker_fd,
                    campaign_progress_fd,
                    prearmed_history_fd, main_fd, candidate_fd,
                    lane_binding_fd):
                if descriptor >= 0:
                    try:
                        os.close(descriptor)
                    except BaseException as error:
                        cleanup_errors.append(error)
            if lane_lock_fd >= 0:
                try:
                    fcntl.flock(lane_lock_fd, fcntl.LOCK_UN)
                except BaseException as error:
                    cleanup_errors.append(error)
                try:
                    os.close(lane_lock_fd)
                except BaseException as error:
                    cleanup_errors.append(error)
            if lane_fd >= 0:
                try:
                    os.close(lane_fd)
                except BaseException as error:
                    cleanup_errors.append(error)
            if runtime_guard is not None:
                try:
                    runtime_guard.close()
                except BaseException as error:
                    cleanup_errors.append(error)
            if canonical_held and canonical_lock is not None:
                try:
                    canonical_lock.__exit__(None, None, None)
                except BaseException as error:
                    cleanup_errors.append(error)
            if cleanup_errors and active_error is None:
                raise ArmingError(
                    "acquisition capability cleanup failed") \
                    from cleanup_errors[0]


class ArmedCampaign:
    """Opaque, process-bound authority for one exact sequential campaign."""

    __slots__ = (
        "attempt_path", "record", "preregistration", "logical_plan",
        "_canonical_lock", "_runtime_guard", "_lane", "_lane_fd",
        "_lane_lock_fd", "_lane_binding_fd", "_campaign_marker_fd",
        "_campaign_progress_fd", "_campaign_progress_binding",
        "_attempts_fd", "_attempt_fd", "_journal_fd",
        "_launch_marker_fd", "_armed_fd", "_bundle_fd", "_manifest_fd",
        "_candidate_fd", "_main_fd", "_candidate_authority",
        "_exact_main", "_lane_binding", "_authority_bundle",
        "_armed_data", "_bundle_data", "_journal_binding",
        "_launch_marker_binding", "_campaign_marker_binding",
        "_manifest_data", "_lane_binding_data", "_attempts_binding",
        "_creator_pid", "_creator_start_ticks", "_state", "_closed",
        "_cleanup_complete",
        "_journal_chain", "_journal_records",
        "_transcript_attempt_binding", "_campaign_marker_writable",
        "_operation_lock", "_operation_owner",
    )

    def __init__(
        self, *, _token: object, attempt_path: Path,
        record: Mapping[str, Any], canonical_lock: Any,
        runtime_guard: _RetainedRuntimeGuard, lane: Path,
        lane_fd: int, lane_lock_fd: int, lane_binding_fd: int,
        campaign_marker_fd: int, campaign_progress_fd: int,
        candidate_control_fd: int, main_fd: int, attempts_fd: int,
        attempt_fd: int,
        journal_fd: int, launch_marker_fd: int, armed_fd: int,
        bundle_fd: int, manifest_fd: int,
        registration: Mapping[str, Any], logical_plan: Mapping[str, Any],
        candidate_authority: Mapping[str, Any], exact_main: Mapping[str, Any],
        lane_binding: Mapping[str, Any], authority_bundle: Mapping[str, Any],
        candidate_snapshot: Mapping[str, Any],
        main_snapshot: Mapping[str, Any],
    ) -> None:
        _require(_token is _ARMED_CAMPAIGN_TOKEN,
                 "ArmedCampaign is created only by acquire_and_arm")
        self.attempt_path = attempt_path
        self.record = copy.deepcopy(record)
        self.preregistration = copy.deepcopy(registration)
        self.logical_plan = copy.deepcopy(logical_plan)
        self._canonical_lock = canonical_lock
        self._runtime_guard = runtime_guard
        self._lane = lane
        self._lane_fd = lane_fd
        self._lane_lock_fd = lane_lock_fd
        self._lane_binding_fd = lane_binding_fd
        self._campaign_marker_fd = campaign_marker_fd
        self._campaign_progress_fd = campaign_progress_fd
        self._attempts_fd = attempts_fd
        self._attempt_fd = attempt_fd
        self._journal_fd = journal_fd
        self._launch_marker_fd = launch_marker_fd
        self._armed_fd = armed_fd
        self._bundle_fd = bundle_fd
        self._manifest_fd = manifest_fd
        self._candidate_fd = candidate_control_fd
        self._main_fd = main_fd
        self._candidate_authority = copy.deepcopy(candidate_authority)
        self._exact_main = copy.deepcopy(exact_main)
        # Avoid duplicating multi-megabyte executable bytes in the retained
        # records; the sealed memfds are the executable capabilities.
        self._candidate_authority.pop("artifact_data", None)
        self._exact_main.pop("artifact_data", None)
        self._lane_binding = copy.deepcopy(lane_binding)
        self._authority_bundle = copy.deepcopy(authority_bundle)
        descriptors = self._authority_bundle["descriptor_binding"]
        _validate_descriptor_bundle(
            descriptors, self._authority_bundle["artifact_bundle"])
        _require(
            contract.exact_json_equal(
                candidate_snapshot, descriptors["candidate_control"]) and
            contract.exact_json_equal(main_snapshot, descriptors["main"]),
            "retained executable snapshots differ from the authority bundle")
        self._armed_data = contract.canonical_json_bytes(self.record)
        self._bundle_data = contract.canonical_json_bytes(authority_bundle)
        self._journal_binding = _journal_binding_record(
            os.fstat(self._journal_fd))
        self._launch_marker_binding = _launch_marker_binding_record(
            os.fstat(self._launch_marker_fd))
        self._campaign_marker_binding = _validate_campaign_marker_binding(
            self._lane_binding["campaign_markers"]
                              [self.record["evidence_attempt"] - 1])
        self._campaign_progress_binding = \
            _validate_campaign_progress_directory_binding(
                self._lane_binding["campaign_progress_directories"]
                                  [self.record["evidence_attempt"] - 1])
        _validate_campaign_checkpoint_payloads(
            self._campaign_progress_fd, self._campaign_progress_binding,
            armed=self.record, journal_records=[],
            label="new campaign checkpoints")
        _validate_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding, journal_records=[],
            armed=self.record, label="new campaign transcript slot")
        self._manifest_data = contract.canonical_json_bytes(
            _attempt_manifest_record(
                self._bundle_data, self._journal_binding,
                self._launch_marker_binding,
                self._campaign_marker_binding,
                self._campaign_progress_binding))
        self._lane_binding_data = contract.canonical_json_bytes(lane_binding)
        self._attempts_binding = _attempts_binding_record(
            os.fstat(self._attempts_fd))
        self._creator_pid = os.getpid()
        self._creator_start_ticks = _process_start_ticks()
        self._state = "ready"
        self._closed = False
        self._cleanup_complete = False
        self._operation_lock = threading.Lock()
        self._operation_owner: int | None = None
        self._journal_chain: str | None = None
        self._journal_records: list[tuple[str, str]] = []
        self._transcript_attempt_binding = \
            _campaign_transcript_binding_sha256(self.record)
        self._campaign_marker_writable = (
            fcntl.fcntl(self._campaign_marker_fd, fcntl.F_GETFL) &
            os.O_ACCMODE) == os.O_RDWR
        _require(self._campaign_marker_writable,
                 "new campaign transcript is not retained writable")
        _require(candidate_control_fd != main_fd and
                 (candidate_snapshot["device"], candidate_snapshot["inode"]) !=
                 (main_snapshot["device"], main_snapshot["inode"]),
                 "candidate/control and main do not hold distinct memfds")

    def __enter__(self) -> "ArmedCampaign":
        _require(not self._closed, "ArmedCampaign is closed")
        return self

    def __exit__(self, unused_type: Any, unused_value: Any,
                 unused_traceback: Any) -> None:
        self.close()

    def close(self) -> None:
        # Never touch a threading.Lock inherited from another process: it may
        # have been owned by a parent thread that does not exist in this child.
        if os.getpid() != self._creator_pid:
            self._reject_foreign_process()
        with _FORK_BARRIER:
            _require(self._operation_owner != threading.get_ident(),
                     "ArmedCampaign cannot close from its running operation")
            with self._operation_lock:
                self._close_locked()

    def _reject_foreign_process(self) -> NoReturn:
        _require(os.getpid() != self._creator_pid,
                 "foreign-process rejection requires a forked campaign")
        # Registered at-fork cleanup normally reaches this object first.  Keep
        # the fallback serialized as well so a nested fork cannot snapshot a
        # partially discarded foreign authority.
        with _FORK_BARRIER:
            self._close_locked()
        # A parent close may already have set ``_closed`` when the fork took
        # its memory snapshot.  Cleanup can then be a no-op, but this branch
        # must still be terminal before any inherited Python lock is touched.
        raise ArmingError(
            "forked process has no live campaign authority")

    def _close_locked(self) -> None:
        errors: list[BaseException] = []
        foreign_process = True
        try:
            foreign_process = (
                os.getpid() != self._creator_pid or
                _process_start_ticks() != self._creator_start_ticks)
        except BaseException as error:
            errors.append(error)
        if self._cleanup_complete:
            _unregister_active_campaign(self)
            if foreign_process:
                raise ArmingError(
                    "forked process has no live campaign authority")
            return
        self._closed = True
        descriptor_attributes = (
            "_main_fd", "_candidate_fd", "_manifest_fd", "_bundle_fd",
            "_armed_fd", "_launch_marker_fd", "_journal_fd",
            "_campaign_marker_fd", "_campaign_progress_fd",
            "_attempt_fd", "_attempts_fd",
            "_lane_binding_fd",
        )
        for attribute in descriptor_attributes:
            descriptor = getattr(self, attribute)
            if descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError as error:
                    if not (foreign_process and error.errno == errno.EBADF):
                        errors.append(error)
                finally:
                    setattr(self, attribute, -1)
        if self._lane_lock_fd >= 0:
            descriptor = self._lane_lock_fd
            if not foreign_process:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                except OSError as error:
                    errors.append(error)
            try:
                os.close(descriptor)
            except OSError as error:
                if not (foreign_process and error.errno == errno.EBADF):
                    errors.append(error)
            finally:
                self._lane_lock_fd = -1
        if self._lane_fd >= 0:
            descriptor = self._lane_fd
            try:
                os.close(descriptor)
            except OSError as error:
                if not (foreign_process and error.errno == errno.EBADF):
                    errors.append(error)
            finally:
                self._lane_fd = -1
        try:
            self._runtime_guard.close()
        except BaseException as error:
            errors.append(error)
        if foreign_process:
            descriptor = getattr(self._canonical_lock, "descriptor", -1)
            if type(descriptor) is int and descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError as error:
                    if error.errno != errno.EBADF:
                        errors.append(error)
                try:
                    self._canonical_lock.descriptor = -1
                    self._canonical_lock.identity = None
                except BaseException as error:
                    errors.append(error)
        else:
            try:
                self._canonical_lock.__exit__(None, None, None)
            except BaseException as error:
                errors.append(error)
        if not errors:
            try:
                _unregister_active_campaign(self)
            except BaseException as error:
                errors.append(error)
        self._cleanup_complete = not errors
        if errors:
            raise ArmingError("ArmedCampaign capability cleanup failed") \
                from errors[0]
        if foreign_process:
            raise ArmingError(
                "forked process discarded inherited campaign capabilities "
                "without unlocking parent authority")

    def _revalidate_authority(self) -> None:
        _require(not self._closed and
                 os.getpid() == self._creator_pid and
                 _process_start_ticks() == self._creator_start_ticks,
                 "ArmedCampaign is closed, forked, or process-replaced")
        try:
            self._canonical_lock.validate_current()
        except Exception as error:
            raise ArmingError("canonical build/test lock changed") from error
        _validate_lane_capability(
            self._lane, self._lane_fd, self._lane_lock_fd,
            self.preregistration)
        _validate_preallocated_lane_lock(
            self._lane_fd, self._lane_lock_fd,
            self._lane_binding["arming_lock"])
        lane_binding_status = os.fstat(self._lane_binding_fd)
        lane_binding_path = os.stat(
            LANE_BINDING_FILE, dir_fd=self._lane_fd, follow_symlinks=False)
        _require((lane_binding_status.st_dev, lane_binding_status.st_ino) ==
                 (lane_binding_path.st_dev, lane_binding_path.st_ino) and
                 _read_exact_fd(
                     self._lane_binding_fd, MAX_JSON_BYTES,
                     "retained lane binding") == self._lane_binding_data,
                 "retained lane binding changed")
        prereg.validate_output_lane_descriptor_identity_for_preregistration(
            self._lane, self._lane_fd, self.preregistration)
        _require(
            set(os.listdir(self._lane_fd)) ==
                _expected_lane_inventory(self._lane_binding),
            "live generation-3 lane inventory differs")
        _validate_attempts_capability(
            self._lane_fd, self._attempts_fd, self._attempts_binding,
            self._attempt_fd, self.attempt_path.name)
        _validate_retained_attempt_capability(
            self.attempt_path, self._attempt_fd, self._journal_fd,
            self._journal_binding, self._launch_marker_fd,
            self._launch_marker_binding, self._armed_fd,
            self._bundle_fd, self._manifest_fd, self._armed_data,
            self._bundle_data, self._manifest_data)
        bundle, armed = _validate_authority_bundle_and_armed(
            self._authority_bundle, self.record, lane=self._lane,
            registration=self.preregistration,
            logical_plan=self.logical_plan,
            candidate_authority=self._candidate_authority,
            exact_main=self._exact_main, lane_binding=self._lane_binding,
            authority_bundle_data=self._bundle_data,
            attempt_manifest_data=self._manifest_data)
        _require(contract.exact_json_equal(bundle, self._authority_bundle) and
                 contract.exact_json_equal(armed, self.record),
                 "retained authority bundle replay changed")
        artifacts = bundle["artifact_bundle"]
        descriptors = bundle["descriptor_binding"]
        _validate_descriptor_bundle(descriptors, artifacts)
        revalidate_sealed_descriptor(
            self._candidate_fd, descriptors["candidate_control"],
            "candidate/control")
        revalidate_sealed_descriptor(
            self._main_fd, descriptors["main"], "exact-main")
        _require(contract.exact_json_equal(
            _capture_live_host_instance(self.preregistration),
            bundle["host_instance"]),
            "live host changed before or after a campaign child")
        _validate_launch_context_current(
            bundle["launch_context"],
            "before or after a generation-3 campaign child")
        _validate_controller_bindings_current(self.preregistration)
        self._runtime_guard.validate_current("live campaign runtime closure")
        completed, chain, complete, journal_records = \
            _validate_journal_directory(
            self._attempt_fd, self._journal_fd,
            journal_binding=self._journal_binding,
            registration=self.preregistration,
            logical_plan=self.logical_plan, artifacts=artifacts,
            armed=self.record, launch_context=bundle["launch_context"],
            cleanup_staging=False)
        _require(not any(
            JOURNAL_STAGING_NAME.fullmatch(entry) is not None
            for entry in os.listdir(self._journal_fd)),
            "live execution journal contains a staging remnant")
        _validate_campaign_progress_directory_capability(
            self._lane, self._lane_fd, self._campaign_progress_fd,
            self._campaign_progress_binding,
            "retained campaign checkpoint directory")
        _validate_live_campaign_checkpoint_tail(
            self._campaign_progress_fd, self._campaign_progress_binding,
            armed=self.record, journal_records=journal_records,
            label="retained campaign checkpoints")
        _validate_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            journal_records=journal_records,
            armed=self.record,
            label="retained campaign transcript")
        marker_consumed = _validate_launch_marker_state(
            attempt_fd=self._attempt_fd, marker_fd=self._launch_marker_fd,
            marker_binding=self._launch_marker_binding,
            journal_fd=self._journal_fd, logical_plan=self.logical_plan,
            lane_fd=self._lane_fd,
            campaign_marker_fd=self._campaign_marker_fd,
            campaign_marker_binding=self._campaign_marker_binding,
            artifacts=artifacts, armed=self.record,
            launch_context=bundle["launch_context"],
            label="retained launch-consumed marker")
        _require(marker_consumed == (chain is not None),
                 "launch-consumed marker and journal state differ")
        _require(chain == self._journal_chain and
                 journal_records == self._journal_records and
                 (self._state == "complete") == complete and
                 (self._state != "ready" or completed == 0),
                 "execution journal state differs from the live campaign")
        names, replayed, prepared = _read_prior_attempts_before_recovery(
            self._attempts_fd, lane=self._lane, lane_fd=self._lane_fd,
            registration=self.preregistration,
            logical_plan=self.logical_plan,
            candidate_authority=self._candidate_authority,
            exact_main=self._exact_main, lane_binding=self._lane_binding,
            require_exact=True)
        _require(
            prepared is None and
            len(names) == len(replayed) == self.record["evidence_attempt"] and
            contract.exact_json_equal(replayed[-1][0], self._authority_bundle)
            and contract.exact_json_equal(replayed[-1][1], self.record) and
            replayed[-1][2] == (completed, chain, complete, journal_records),
            "live campaign differs from the exact retained attempt replay")
        budget_authority = _read_budget_authority_before_recovery(
            self._lane_fd, registration=self.preregistration,
            logical_plan=self.logical_plan, lane_binding=self._lane_binding,
            prior_attempts=replayed, require_exact=True)
        budget_commit = _validate_budget_commit(
            self._authority_bundle["budget_commit"], self.preregistration)
        _require(
            contract.exact_json_equal(
                budget_authority["ledger"],
                budget_commit["prospective_ledger"]) and
            budget_authority["ledger"]["evidence_attempts_used"] ==
                self.record["evidence_attempt"],
            "live campaign budget histories differ from exact replay")

    def _selected_executable_capability(
        self, role: str,
    ) -> tuple[int, Mapping[str, Any]]:
        descriptors = self._authority_bundle["descriptor_binding"]
        if role in ("candidate", "control"):
            descriptor = self._candidate_fd
            expected = descriptors["candidate_control"]
            label = "selected candidate/control"
        else:
            _require(role == "main", "logical child implementation differs")
            descriptor = self._main_fd
            expected = descriptors["main"]
            label = "selected exact-main"
        revalidate_sealed_descriptor(descriptor, expected, label)
        return descriptor, expected

    def _consume_launch_authority(self, intent_data: bytes) -> None:
        _require(type(intent_data) is bytes and
                 0 < len(intent_data) <= MAX_JSON_BYTES,
                 "launch-consumed intent bytes are invalid")
        metadata = _validate_launch_marker_capability(
            self._attempt_fd, self._launch_marker_fd,
            self._launch_marker_binding,
            "launch-consumed marker before first child")
        _require(metadata.st_size == 0 and
                 set(os.listdir(self._journal_fd)) == {"intent-0000.json"},
                 "first-child journal was not durably published exactly once")
        intent_fd, published_intent = _read_immutable_file_at(
            self._journal_fd, "intent-0000.json", MAX_JSON_BYTES,
            "first-child durable journal intent")
        try:
            _require(published_intent == intent_data,
                     "first-child durable journal intent bytes differ")
        finally:
            os.close(intent_fd)
        campaign_metadata = _validate_campaign_marker_capability(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            "campaign-consumed marker before first child")
        _require(campaign_metadata.st_size > 0,
                 "campaign transcript was not committed before launch")
        transcript = _read_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            "campaign transcript before first child")
        entries = transcript["journal_entries"]
        _require(transcript["allocation"] is not None and
                 len(entries) == 1 and
                 entries[0]["name"] == "intent-0000.json" and
                 entries[0]["data_sha256"] == _sha256_bytes(intent_data) and
                 entries[0]["attempt_binding_sha256"] ==
                    self._transcript_attempt_binding,
                 "campaign transcript did not commit the first intent")

        def commit(descriptor: int, label: str) -> None:
            offset = 0
            while offset < len(intent_data):
                written = os.pwrite(descriptor, intent_data[offset:], offset)
                _require(written > 0, f"{label} write made no progress")
                offset += written
            os.fsync(descriptor)

        commit(self._launch_marker_fd, "launch-consumed marker")
        metadata = _validate_launch_marker_capability(
            self._attempt_fd, self._launch_marker_fd,
            self._launch_marker_binding,
            "launch-consumed marker after first-child commitment")
        _require(metadata.st_size == len(intent_data) and
                 _read_exact_fd(
                     self._launch_marker_fd, MAX_JSON_BYTES,
                     "launch-consumed marker") == intent_data,
                 "launch-consumed marker bytes differ after commitment")

    def _append_campaign_transcript(
        self, name: str, data: bytes,
    ) -> None:
        _require(self._campaign_marker_writable,
                 "campaign transcript writable capability is absent")
        prospective_records = [
            *self._journal_records, (name, _sha256_bytes(data)),
        ]
        _validate_live_campaign_checkpoint_tail(
            self._campaign_progress_fd, self._campaign_progress_binding,
            armed=self.record, journal_records=prospective_records,
            label="campaign checkpoints before transcript append")
        _validate_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            journal_records=self._journal_records,
            armed=self.record,
            label="campaign transcript before append")
        existing = _read_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            "campaign transcript before append")
        _require(existing["allocation"] is not None,
                 "campaign transcript lacks its ARMED allocation")
        journal_entries = existing["journal_entries"]
        prior = contract.canonical_sha256(
            journal_entries[-1] if journal_entries
            else existing["allocation"])
        entry = _campaign_transcript_entry_record(
            sequence=len(journal_entries) + 1, name=name,
            data_sha256=_sha256_bytes(data), prior_entry_sha256=prior,
            attempt_binding_sha256=self._transcript_attempt_binding)
        entry_data = contract.canonical_json_bytes(entry)
        offset = os.fstat(self._campaign_marker_fd).st_size
        _require(offset + len(entry_data) <= MAX_CAMPAIGN_TRANSCRIPT_BYTES,
                 "campaign transcript exceeds its byte bound")
        written = 0
        while written < len(entry_data):
            count = os.pwrite(
                self._campaign_marker_fd, entry_data[written:],
                offset + written)
            _require(count > 0,
                     "campaign transcript append made no progress")
            written += count
        os.fsync(self._campaign_marker_fd)
        _validate_campaign_transcript(
            self._lane_fd, self._campaign_marker_fd,
            self._campaign_marker_binding,
            journal_records=prospective_records,
            armed=self.record,
            label="campaign transcript after append")

    def _write_journal(
        self, name: str, value: Mapping[str, Any], *,
        consume_launch: bool = False,
    ) -> bytes:
        _validate_journal_capability(
            self._attempt_fd, self._journal_fd, self._journal_binding,
            "execution journal before write")
        _validate_live_campaign_checkpoint_tail(
            self._campaign_progress_fd, self._campaign_progress_binding,
            armed=self.record, journal_records=self._journal_records,
            label="execution journal checkpoints before write")
        data = contract.canonical_json_bytes(value)
        _write_atomic_journal_file_at(
            self._journal_fd, name, data, f"execution journal {name}")
        os.fsync(self._journal_fd)
        _validate_journal_capability(
            self._attempt_fd, self._journal_fd, self._journal_binding,
            "execution journal after durable publication")
        prospective_records = [
            *self._journal_records, (name, _sha256_bytes(data)),
        ]
        checkpoint_blocks = _expected_campaign_transcript_blocks(
            self.record, prospective_records)
        checkpoint_sequence = len(prospective_records)
        _require(checkpoint_sequence == len(self._journal_records) + 1,
                 "journal checkpoint sequence differs")
        _flip_campaign_checkpoint(
            self._campaign_progress_fd, self._campaign_progress_binding,
            checkpoint_sequence, checkpoint_blocks[checkpoint_sequence],
            f"execution journal {name}")
        self._append_campaign_transcript(name, data)
        if consume_launch:
            self._consume_launch_authority(data)
        _validate_journal_capability(
            self._attempt_fd, self._journal_fd, self._journal_binding,
            "execution journal after write")
        self._journal_records.append((name, _sha256_bytes(data)))
        return data

    @_with_exact_source_imports
    def run_all(self) -> dict[str, Any]:
        """Run logical children 0..1649 once, in order, with no retry knobs."""
        if os.getpid() != self._creator_pid:
            self._reject_foreign_process()
        _require(self._operation_owner != threading.get_ident(),
                 "ArmedCampaign.run_all is not reentrant")
        self._operation_lock.acquire()
        self._operation_owner = threading.get_ident()
        try:
            return ArmedCampaign._run_all_locked(self)
        finally:
            self._operation_owner = None
            self._operation_lock.release()

    def _run_all_locked(self) -> dict[str, Any]:
        _require(self._state == "ready" and not self._closed,
                 "ArmedCampaign.run_all may be called exactly once")
        _validate_owner_creation_umask(
            "before generation-3 journal publication")
        self._state = "running"
        try:
            for index, child in enumerate(self.logical_plan["child_plans"]):
                _require(child["index"] == index,
                         "logical child schedule is not gapless")
                self._revalidate_authority()
                role = child["implementation"]
                descriptor, unused_snapshot = \
                    self._selected_executable_capability(role)
                started_ns = time.monotonic_ns()
                intent = _journal_intent_record(
                    child=child, selected_pair=self.record["selected_pair"],
                    artifact_handle_id=
                        self._authority_bundle["artifact_bundle"]["roles"]
                                              [role]["handle_id"],
                    evidence_attempt=self.record["evidence_attempt"],
                    prior_journal_sha256=self._journal_chain,
                    started_monotonic_ns=started_ns,
                    launch_context=self._authority_bundle["launch_context"])
                intended_data = contract.canonical_json_bytes(intent)
                intent_data = self._write_journal(
                    f"intent-{index:04d}.json", intent,
                    consume_launch=(index == 0))
                _require(intent_data == intended_data,
                         "journal intent bytes changed during publication")
                self._journal_chain = _sha256_bytes(intent_data)
                self._revalidate_authority()
                descriptor, expected_snapshot = \
                    self._selected_executable_capability(role)
                command = _actual_child_command(
                    child, self.record["selected_pair"], descriptor)
                revalidate_sealed_descriptor(
                    descriptor, expected_snapshot,
                    "immediate pre-launch executable")
                process_result = _run_child_process(
                    command, descriptor=descriptor,
                    timeout_seconds=child["timeout_budget"]["timeout_seconds"],
                    expected_launch_context=
                        self._authority_bundle["launch_context"])
                _validate_launch_context_current(
                    self._authority_bundle["launch_context"],
                    f"immediately after child {index}")
                self._runtime_guard.validate_current(
                    f"runtime closure after child {index}")
                stdout = process_result["stdout"]
                stderr = process_result["stderr"]
                outcome = process_result["outcome"]
                payload: Mapping[str, Any] | None = None
                normalized: Mapping[str, Any] | None = None
                error_text = process_result["error"]
                if outcome == "accepted":
                    try:
                        parsed = contract.strict_json_loads(
                            stdout, f"child {index} output")
                        _require(type(parsed) is dict,
                                 f"child {index} output is not an object")
                        cell = self.logical_plan["cells"][child["cell_index"]]
                        normalized_value = _normalize_gen3_result(
                            role, parsed, cell, self.preregistration)
                        payload = parsed
                        normalized = normalized_value
                    except Exception as error:
                        outcome = "output-rejected"
                        error_text = f"child output rejected: {error}"[:4096]
                result = _journal_result_record(
                    child_index=index, intent_sha256=self._journal_chain,
                    outcome=outcome,
                    finished_monotonic_ns=time.monotonic_ns(),
                    elapsed_ns=process_result["elapsed_ns"],
                    returncode=process_result["returncode"],
                    stdout_sha256=_sha256_bytes(stdout),
                    stderr_sha256=_sha256_bytes(stderr), payload=payload,
                    normalized=normalized, error=error_text)
                result_data = self._write_journal(
                    f"result-{index:04d}.json", result)
                self._journal_chain = _sha256_bytes(result_data)
                if outcome != "accepted":
                    self._state = "failed"
                    self._revalidate_authority()
                    raise ArmingError(
                        f"generation-3 child {index} terminated: {error_text}")
                self._revalidate_authority()
            complete = {
                "schema": JOURNAL_COMPLETE_SCHEMA,
                "evidence_attempt": self.record["evidence_attempt"],
                "child_count": len(self.logical_plan["child_plans"]),
                "prior_journal_sha256": self._journal_chain,
                "completed_monotonic_ns": time.monotonic_ns(),
            }
            complete_data = self._write_journal("complete.json", complete)
            self._journal_chain = _sha256_bytes(complete_data)
            self._state = "complete"
            self._revalidate_authority()
            return {
                "state": "complete",
                "evidence_attempt": self.record["evidence_attempt"],
                "child_count": len(self.logical_plan["child_plans"]),
                "journal_chain_sha256": self._journal_chain,
            }
        except BaseException:
            if self._state == "running":
                self._state = "failed"
            raise


def main() -> int:
    print(
        "K65 generation-3 live armer has no standalone timing CLI; "
        "use the reviewed acquisition workflow",
        file=sys.stderr)
    return 2


def _register_active_campaign(
    campaign: ArmedCampaign,
    registry: dict[int, Any] = _ACTIVE_CAMPAIGNS,
) -> None:
    """Retain every live capability until explicit or fork-child cleanup."""
    _require(type(registry) is dict,
             "active campaign registry type changed")
    key = id(campaign)
    current = registry.get(key, _MISSING_MODULE_ALIAS)
    _require(current is _MISSING_MODULE_ALIAS or current is campaign,
             "active campaign registry identity collision")
    registry[key] = campaign


def _unregister_active_campaign(
    campaign: ArmedCampaign,
    registry: dict[int, Any] = _ACTIVE_CAMPAIGNS,
) -> None:
    """Forget a campaign only after all of its capabilities are revoked."""
    _require(type(registry) is dict,
             "active campaign registry type changed")
    key = id(campaign)
    current = registry.get(key, _MISSING_MODULE_ALIAS)
    _require(current is _MISSING_MODULE_ALIAS or current is campaign,
             "active campaign registry entry changed")
    if current is campaign:
        del registry[key]


def _discard_inherited_campaigns_after_fork(
    registry: dict[int, Any] = _ACTIVE_CAMPAIGNS,
    campaign_type: type = ArmedCampaign,
    close_method: Callable[[ArmedCampaign], None] = ArmedCampaign._close_locked,
) -> bool:
    """Revoke the child copy of every process-bound campaign capability."""
    if type(registry) is not dict:
        return False
    campaigns = tuple(registry.items())
    valid = all(
        type(key) is int and key == id(campaign) and
        isinstance(campaign, campaign_type)
        for key, campaign in campaigns)
    for unused_key, campaign in campaigns:
        if not isinstance(campaign, campaign_type):
            continue
        try:
            close_method(campaign)
        except ArmingError:
            pass
        except BaseException:
            valid = False
        if not campaign._cleanup_complete:
            valid = False
    return valid and not registry


def _arm_controlled_process_fork(
    descriptor: int,
    registry: dict[int, Any] = _ACTIVE_CAMPAIGNS,
    holder: list[tuple[int, int, int] | None] = _CONTROLLED_FORK_PERMIT,
) -> tuple[int, int, int] | None:
    """Permit one Popen fork to carry a selected campaign fd into exec."""
    _require(type(descriptor) is int and descriptor >= 0,
             "controlled process fork descriptor is invalid")
    _require(type(registry) is dict and type(holder) is list and
             len(holder) == 1 and holder[0] is None,
             "controlled process fork state changed")
    owners = [
        campaign for campaign in registry.values()
        if isinstance(campaign, ArmedCampaign) and descriptor in {
            campaign._candidate_fd, campaign._main_fd,
        }
    ]
    _require(len(owners) <= 1,
             "controlled process fork descriptor has ambiguous authority")
    if not owners:
        return None
    _require(os.getpid() == _BOOTSTRAP_PID and
             not owners[0]._closed and not owners[0]._cleanup_complete,
             "controlled process fork campaign is not live")
    permit = (os.getpid(), threading.get_ident(), descriptor)
    holder[0] = permit
    return permit


def _clear_controlled_process_fork(
    permit: tuple[int, int, int] | None,
    holder: list[tuple[int, int, int] | None] = _CONTROLLED_FORK_PERMIT,
) -> None:
    """Clear the parent copy of a one-shot process-fork permit."""
    _require(type(holder) is list and len(holder) == 1,
             "controlled process fork holder changed")
    if permit is None:
        _require(holder[0] is None,
                 "unexpected controlled process fork permit appeared")
        return
    _require(holder[0] is None or holder[0] == permit,
             "controlled process fork permit changed")
    holder[0] = None


def _consume_controlled_process_fork_in_child(
    registry: dict[int, Any] = _ACTIVE_CAMPAIGNS,
    holder: list[tuple[int, int, int] | None] = _CONTROLLED_FORK_PERMIT,
) -> bool:
    """Consume only the exact same-thread Popen permit copied by fork."""
    if type(holder) is not list or len(holder) != 1:
        return False
    permit = holder[0]
    holder[0] = None
    if (type(permit) is not tuple or len(permit) != 3 or
            type(permit[0]) is not int or type(permit[1]) is not int or
            type(permit[2]) is not int):
        return False
    parent_pid, parent_thread, descriptor = permit
    if (parent_pid != _BOOTSTRAP_PID or os.getpid() == parent_pid or
            os.getppid() != parent_pid or
            threading.get_ident() != parent_thread or descriptor < 0 or
            type(registry) is not dict):
        return False
    owners = [
        campaign for campaign in registry.values()
        if isinstance(campaign, ArmedCampaign) and descriptor in {
            campaign._candidate_fd, campaign._main_fd,
        }
    ]
    if len(owners) != 1:
        return False
    try:
        os.fstat(descriptor)
    except OSError:
        return False
    return True


def _before_process_fork(
    barrier: threading.RLock = _FORK_BARRIER,
) -> None:
    """Wait until no Gen3 authority operation is between durable states."""
    barrier.acquire()


def _after_process_fork_in_parent(
    barrier: threading.RLock = _FORK_BARRIER,
) -> None:
    """Release the prepare lease in the unchanged parent process."""
    barrier.release()


def _after_process_fork_in_child(
    barrier: threading.RLock = _FORK_BARRIER,
    discard: Callable[[], bool] = _discard_inherited_campaigns_after_fork,
    consume_permit: Callable[[], bool] =
        _consume_controlled_process_fork_in_child,
    fail_stop: Callable[[int], NoReturn] = os._exit,
) -> None:
    """Revoke inherited authority before any child Python code can run."""
    cleanup_succeeded = False
    try:
        try:
            cleanup_succeeded = consume_permit() or discard()
        except BaseException:
            cleanup_succeeded = False
    finally:
        try:
            barrier.release()
        except BaseException:
            cleanup_succeeded = False
    if not cleanup_succeeded:
        # Returning from a failed child cleanup would expose process-bound
        # descriptors to arbitrary post-fork code.  Exit before os.fork can
        # return in the child instead.
        fail_stop(125)


_require(hasattr(os, "register_at_fork"),
         "generation-3 authority requires Python at-fork serialization")
os.register_at_fork(
    before=_before_process_fork,
    after_in_parent=_after_process_fork_in_parent,
    after_in_child=_after_process_fork_in_child)


_LIVE_ARMER_EXECUTABLE_REFERENCE_IDENTITIES = tuple(
    (_LIVE_ARMER_MODULE, name, value)
    for name, value in globals().items()
    if ((type(value) is types.FunctionType and
         value.__globals__ is globals()) or
        (isinstance(value, type) and
         getattr(value, "__module__", None) == __name__))
)
_LIVE_ARMER_CLASS_MEMBER_REFERENCE_IDENTITIES = tuple(
    (class_value, member_name, member)
    for unused_owner, unused_name, class_value in
        _LIVE_ARMER_EXECUTABLE_REFERENCE_IDENTITIES
    if isinstance(class_value, type)
    for member_name, member in vars(class_value).items()
    if (callable(member) or
        isinstance(member, (staticmethod, classmethod, property)))
)
_LIVE_ARMER_FUNCTION_EXECUTION_STATES = _capture_function_execution_states((
    *(value for unused_owner, unused_name, value in
      _LIVE_ARMER_EXECUTABLE_REFERENCE_IDENTITIES),
    *(value for unused_owner, unused_name, value in
      _LIVE_ARMER_CLASS_MEMBER_REFERENCE_IDENTITIES),
))


def _frozen_live_execution_graph(
    executables: tuple[tuple[Any, str, Any], ...] =
        _LIVE_ARMER_EXECUTABLE_REFERENCE_IDENTITIES,
    class_members: tuple[tuple[Any, str, Any], ...] =
        _LIVE_ARMER_CLASS_MEMBER_REFERENCE_IDENTITIES,
    function_states: tuple[tuple[Any, ...], ...] =
        _LIVE_ARMER_FUNCTION_EXECUTION_STATES,
) -> tuple[tuple[tuple[Any, str, Any], ...],
           tuple[tuple[Any, str, Any], ...],
           tuple[tuple[Any, ...], ...]]:
    """Return the post-definition live-armer callable identity graph."""
    return executables, class_members, function_states


_LIVE_PUBLIC_CALLABLE_IDENTITIES = (
    (_LIVE_ARMER_MODULE, "acquire_and_arm", acquire_and_arm),
    (_LIVE_ARMER_MODULE, "ArmedCampaign", ArmedCampaign),
    (_LIVE_ARMER_MODULE, "ArmingError", ArmingError),
    (_LIVE_ARMER_MODULE, "_frozen_controller_execution_graph",
     _frozen_controller_execution_graph),
    (_LIVE_ARMER_MODULE, "_frozen_child_environment",
     _frozen_child_environment),
    (_LIVE_ARMER_MODULE, "_frozen_live_execution_graph",
     _frozen_live_execution_graph),
    (ArmedCampaign, "run_all", ArmedCampaign.run_all),
)


def _frozen_live_public_surface(
    bindings: tuple[tuple[Any, str, Any], ...] =
        _LIVE_PUBLIC_CALLABLE_IDENTITIES,
) -> tuple[tuple[Any, str, Any], ...]:
    """Return post-bootstrap public callable identities for scope gates."""
    return bindings


__all__ = (
    "ArmedCampaign", "ArmingError", "acquire_and_arm",
)


if __name__ == "__main__":
    raise SystemExit(main())
