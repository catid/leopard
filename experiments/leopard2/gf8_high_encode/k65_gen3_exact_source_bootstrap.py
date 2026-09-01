#!/usr/bin/env python3
"""Direct-source bootstrap for the generation-3 K65 live armer.

The live armer is deliberately not importable through Python's ordinary source
loader because a timestamp-valid ``.pyc`` cannot prove which exact source bytes
were executed.  Obtain :func:`load_live_armer` by executing this file directly
from source (for example with :func:`runpy.run_path`).  The function performs a
stable descriptor read, compiles those bytes without a bytecode cache, and
passes the exact code object and byte digest into the armer before execution.

Trust boundary: the caller supplies a trusted CPython interpreter and standard
library state.  This bootstrap attests the reviewed repository-source closure;
it is not an interpreter or arbitrary preseeded-stdlib attestation mechanism.
"""

from __future__ import annotations

import _frozen_importlib
import hashlib
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any


if (type(sys._getframe) is not type(len) or
        sys._getframe.__self__ is not sys):
    raise RuntimeError("Python sys module alias was replaced")
_MODULE_TYPE = type(sys._getframe.__self__)
if type(sys) is not _MODULE_TYPE or type(_frozen_importlib) is not _MODULE_TYPE:
    raise RuntimeError("Python bootstrap module alias was replaced")
_EXECUTING_BOOTSTRAP_CODE = sys._getframe().f_code
_MODULE_SPEC_TYPE = type(sys.__spec__)
_FUNCTION_TYPE = type(lambda: None)
if _frozen_importlib.ModuleSpec is not _MODULE_SPEC_TYPE:
    raise RuntimeError("Python module specification type was replaced")
HERE = Path(__file__).resolve().parent
LIVE_ARMER_PATH = HERE / "run_k65r65_b64_packed_terminal_gen3_acquire.py"
MAX_SOURCE_BYTES = 16 * 1024 * 1024
MARKER_NAME = "__k65_gen3_exact_source_bootstrap__"
# Exact-byte self-digest trust anchor.
_BOOTSTRAP_SELF_SHA256 = \
    "0b768924c5929c023d1b2c133f3d2a205c2a07d78c6da0921bdb5c4d07eff3bd"
_BOOTSTRAP_DIGEST_PATTERN = re.compile(
    rb'(?m)^_BOOTSTRAP_SELF_SHA256 = \\\n'
    rb'    "([0-9a-f]{64})"$')


def _validated_module_registry() -> dict[str, Any]:
    registry = sys.modules
    if type(registry) is not dict:
        raise RuntimeError("Python module registry was replaced")
    aliases = tuple(registry.keys())
    if any(type(alias) is not str for alias in aliases):
        raise RuntimeError("Python module registry contains a non-string key")
    return registry


_validated_module_registry()


def _read_exact_source_bytes(path: Path) -> tuple[bytes, str]:
    expected = path.resolve(strict=True)
    descriptor = os.open(
        expected, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        pathname = os.stat(expected, follow_symlinks=False)
        if not (stat.S_ISREG(before.st_mode) and
                before.st_uid == os.geteuid() and before.st_nlink == 1 and
                0 < before.st_size <= MAX_SOURCE_BYTES and
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


def _validate_direct_source_execution() -> str:
    path = Path(__file__).resolve(strict=True)
    source, digest = _read_exact_source_bytes(path)
    matches = list(_BOOTSTRAP_DIGEST_PATTERN.finditer(source))
    if len(matches) != 1 or \
            matches[0].group(1).decode("ascii") != _BOOTSTRAP_SELF_SHA256:
        raise RuntimeError("K65 generation-3 bootstrap self-digest is malformed")
    normalized = (
        source[:matches[0].start(1)] + b"0" * 64 +
        source[matches[0].end(1):])
    if hashlib.sha256(normalized).hexdigest() != _BOOTSTRAP_SELF_SHA256:
        raise RuntimeError("K65 generation-3 bootstrap exact bytes changed")
    current = compile(
        source, str(path), "exec", dont_inherit=True,
        optimize=sys.flags.optimize)
    if (globals().get("__cached__") is not None or
            globals().get("__spec__") is not None or
            current != _EXECUTING_BOOTSTRAP_CODE):
        raise RuntimeError(
            "K65 generation-3 bootstrap must be executed directly from "
            "current source")
    return digest


_DIRECT_SOURCE_BOOTSTRAP_SHA256 = _validate_direct_source_execution()


class _LiveArmerSourceLoader:
    def __init__(
        self, name: str, path: Path, source: bytes, source_sha256: str,
        bootstrap_path: Path, bootstrap_sha256: str,
    ) -> None:
        self.name = name
        self.path = path
        self.source = source
        self.source_sha256 = source_sha256
        self.bootstrap_path = bootstrap_path
        self.bootstrap_sha256 = bootstrap_sha256

    def create_module(self, unused_specification: Any) -> None:
        return None

    def exec_module(self, module: Any) -> None:
        specification = getattr(module, "__spec__", None)
        if (type(module) is not _MODULE_TYPE or
                type(specification) is not _MODULE_SPEC_TYPE or
                specification.name != self.name or
                specification.loader is not self or
                specification.origin != str(self.path) or
                specification.submodule_search_locations is not None):
            raise RuntimeError("live armer module specification changed")
        code = compile(
            self.source, str(self.path), "exec", dont_inherit=True,
            optimize=sys.flags.optimize)
        module.__dict__[MARKER_NAME] = {
            "code": code,
            "source_path": self.path,
            "source_sha256": self.source_sha256,
            "bootstrap_path": self.bootstrap_path,
            "bootstrap_sha256": self.bootstrap_sha256,
            "bootstrap_direct_source": True,
        }
        exec(code, module.__dict__)


def _validate_live_armer_surface(module: Any, source_path: Path) -> None:
    namespace = module.__dict__
    acquire = namespace.get("acquire_and_arm")
    campaign_type = namespace.get("ArmedCampaign")
    run_all = (vars(campaign_type).get("run_all")
               if type(campaign_type) is type else None)
    functions = (acquire, getattr(acquire, "__wrapped__", None), run_all,
                 getattr(run_all, "__wrapped__", None))
    if (type(namespace) is not dict or
            type(acquire) is not _FUNCTION_TYPE or
            type(campaign_type) is not type or
            type(run_all) is not _FUNCTION_TYPE or
            namespace.get("ArmingError") is None or
            type(namespace.get("ArmingError")) is not type or
            namespace.get("__all__") != (
                "ArmedCampaign", "ArmingError", "acquire_and_arm") or
            any(type(function) is not _FUNCTION_TYPE or
                function.__globals__ is not namespace or
                Path(function.__code__.co_filename).resolve(strict=True) !=
                    source_path
                for function in functions)):
        raise RuntimeError("live armer callable surface changed during bootstrap")


def load_live_armer(
    name: str = "leopard2_k65_gen3_live_armer",
) -> Any:
    """Load the live armer from exact current source, never from ``.pyc``."""
    module_registry = _validated_module_registry()
    if not (type(name) is str and name and name not in module_registry):
        raise RuntimeError("live armer module name is invalid or already loaded")
    bootstrap_sha256 = _validate_direct_source_execution()
    if bootstrap_sha256 != _DIRECT_SOURCE_BOOTSTRAP_SHA256:
        raise RuntimeError(
            "K65 generation-3 bootstrap changed after direct execution")
    bootstrap_path = Path(__file__).resolve(strict=True)
    source_path = LIVE_ARMER_PATH.resolve(strict=True)
    source, source_sha256 = _read_exact_source_bytes(source_path)
    loader = _LiveArmerSourceLoader(
        name, source_path, source, source_sha256,
        bootstrap_path, bootstrap_sha256)
    specification = _MODULE_SPEC_TYPE(
        name, loader, origin=str(source_path), is_package=False)
    specification._set_fileattr = True
    specification.cached = None
    module = _MODULE_TYPE(name)
    module.__loader__ = loader
    module.__package__ = specification.parent
    module.__spec__ = specification
    module.__file__ = str(source_path)
    module.__cached__ = None
    if (type(specification) is not _MODULE_SPEC_TYPE or
            specification.loader is not loader or
            type(module) is not _MODULE_TYPE or
            module.__spec__ is not specification or
            module.__loader__ is not loader):
        raise RuntimeError("cannot construct exact-source live armer module")
    module_registry[name] = module
    try:
        loader.exec_module(module)
    except BaseException:
        module_registry.pop(name, None)
        raise
    if (type(module) is not _MODULE_TYPE or
            sys.modules is not module_registry or
            module_registry.get(name) is not module):
        module_registry.pop(name, None)
        raise RuntimeError("live armer module identity changed during bootstrap")
    try:
        _validate_live_armer_surface(module, source_path)
    except BaseException:
        module_registry.pop(name, None)
        raise
    return module


__all__ = ("load_live_armer",)
