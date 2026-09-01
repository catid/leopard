#!/usr/bin/env python3
"""Construct and independently replay K65 generation-3 candidate authority.

This module is deliberately offline with respect to benchmark execution.  It
validates the exact sealed source/build/runtime payload consumed by the live
generation-3 armer.  Its only ELF process inspection is bounded ``ldd`` replay
through immutable descriptors; it never enters the benchmark program or any
timed mode.  The acquisition controller uses the CLI in a fresh isolated
Python process before publishing an authority lane.
"""

from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import selectors
import signal
import stat
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence


_REPOSITORY_PYCACHE_PREFIX = "/dev/null"
if sys.pycache_prefix not in (None, _REPOSITORY_PYCACHE_PREFIX):
    raise RuntimeError("candidate verifier Python cache prefix is unsafe")
sys.pycache_prefix = _REPOSITORY_PYCACHE_PREFIX
sys.dont_write_bytecode = True


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
TOOLS = REPO_ROOT / "tools"
MAIN_COMPARE = HERE.parent / "main_compare"


def _inside_repository(path: Path) -> bool:
    try:
        path.relative_to(REPO_ROOT)
        return True
    except ValueError:
        return False


# No repository directory remains an import search surface.  Every repository
# dependency below is loaded from its one canonical pathname instead.
_initial_search_path = list(sys.path)
_safe_search_path: list[str] = []
for _entry in sys.path:
    try:
        _resolved_entry = Path(_entry or os.getcwd()).resolve(strict=True)
    except OSError:
        _safe_search_path.append(_entry)
        continue
    if not _inside_repository(_resolved_entry):
        _safe_search_path.append(_entry)
sys.path[:] = _safe_search_path
del _entry, _resolved_entry, _safe_search_path


def _load_repository_module(
    module_name: str, expected_path: Path, label: str,
) -> Any:
    """Execute one repository module from its exact canonical pathname."""
    expected = expected_path.resolve(strict=True)
    if expected != expected_path or not expected.is_file():
        raise RuntimeError(f"{label} path is not one canonical file")
    specification = importlib.util.spec_from_file_location(
        module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {label} from its exact path")
    module = importlib.util.module_from_spec(specification)
    previous = sys.modules.get(module_name)
    original_search_path = list(sys.path)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
        observed = Path(
            str(getattr(module, "__file__", ""))).resolve(strict=True)
        if observed != expected:
            raise RuntimeError(f"{label} resolved outside this source tree")
    except BaseException:
        if previous is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    finally:
        # Legacy modules may temporarily prepend a sibling directory.  Do not
        # allow that directory to become a later or lazy import surface.
        sys.path[:] = original_search_path
    return module


# Preinstall exact dependencies under the public names used by legacy
# transitive imports, then load their consumers under private names.
_identity_contract = _load_repository_module(
    "leopard2_exact_main_baseline",
    TOOLS / "leopard2_exact_main_baseline.py",
    "exact-main identity contract")
_record_contract = _load_repository_module(
    "leopard2_exact_main_baseline_record",
    TOOLS / "leopard2_exact_main_baseline_record.py",
    "exact-main record contract")
build_provenance = _load_repository_module(
    "leopard2_build_provenance",
    TOOLS / "leopard2_build_provenance.py",
    "build-provenance module")
sealed_verifier = _load_repository_module(
    "leopard2_k65_gen3_exact_main_verifier",
    TOOLS / "leopard2_exact_main_baseline_verifier.py",
    "sealed-tree verifier module")
git_capture = _load_repository_module(
    "leopard2_k65_gen3_git_capture",
    MAIN_COMPARE / "git_capture.py", "Git-capture module")
_pair_for_bridge = _load_repository_module(
    "leopard2_pair_qualification_contract_for_bridge",
    MAIN_COMPARE / "pair_qualification_contract.py",
    "bridge pair-qualification contract")
_bridge_for_prereg = _load_repository_module(
    "leopard2_bridge_contract_for_k65_gen3_prereg",
    MAIN_COMPARE / "pair_qualification_bridge_contract.py",
    "preregistration bridge contract")
_pair_for_prereg = _load_repository_module(
    "leopard2_pair_qualification_contract_for_k65_gen3_prereg",
    MAIN_COMPARE / "pair_qualification_contract.py",
    "preregistration pair-qualification contract")
prereg = _load_repository_module(
    "leopard2_k65_gen3_preregistration",
    HERE / "k65_gen3_preregistration.py",
    "generation-3 preregistration module")
if __name__ != "__main__":
    # A library import must not remove paths from its caller's search list.  The
    # exact module aliases intentionally remain installed for transitive users.
    sys.path[:] = _initial_search_path
del _initial_search_path


contract = prereg.contract
GENERATION = 3
CANDIDATE_AUTHORITY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-candidate-authority/v1"
VERDICT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-candidate-verdict/v1"
RUNTIME_CLOSURE_SCHEMA = "leopard2-k65-gen3-runtime-closure/v1"
RUNTIME_DISCOVERY_SCHEMA = \
    "leopard2-k65-gen3-runtime-discovery-transcript/v1"
TERMINAL_PATH = "candidate-authority.json"
MAX_EXECUTABLE_BYTES = 64 * 1024 * 1024
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_AUTHORITY_PAYLOAD_BYTES = 128 * 1024 * 1024
MAX_RUNTIME_CLOSURE_FILES = 64
MAX_LDD_OUTPUT_BYTES = 4 * 1024 * 1024
LDD_TIMEOUT_SECONDS = 120
LDD_PATH = "/usr/bin/ldd"
LDD_INTERPRETER_PATH = "/bin/bash"
LDD_DISCOVERY_ROLE_PREFIX = "discovery:"
LDD_INTERPRETER_ROLE = "discovery-interpreter"
PAYLOAD_PATHS = (
    "artifacts/bench_leopard2",
    "build/benchmark-source-attestation.h",
    "build/build-provenance.json",
    "build/reproducible-build-core.json",
    "build/reproducible-build-proof.json",
    "runtime/runtime-closure.json",
    "source/candidate-source.tar",
    "source/git-capture.json",
)
SEALED_PATHS = frozenset(PAYLOAD_PATHS) | {
    TERMINAL_PATH, "SHA256SUMS", "TREE-METADATA.json"}
_JSON_PAYLOAD_PATHS = frozenset((
    "build/build-provenance.json",
    "build/reproducible-build-core.json",
    "build/reproducible-build-proof.json",
    "runtime/runtime-closure.json",
    "source/git-capture.json",
))

_RECORD_KEYS = frozenset((
    "schema", "generation", "status", "source", "candidate",
    "build_closure", "controller_bindings_sha256", "inventory",
))
_SOURCE_KEYS = frozenset((
    "commit", "tree", "git_capture_sha256", "source_archive_sha256",
))
_CANDIDATE_KEYS = frozenset(("relative_path", "sha256", "size"))
_BUILD_KEYS = frozenset((
    "build_provenance_sha256", "reproducible_build_core_sha256",
    "reproducible_build_proof_sha256", "attestation_header_sha256",
    "runtime_closure_sha256",
))
_INVENTORY_KEYS = frozenset(("relative_path", "sha256", "size"))
_RUNTIME_ENTRY_KEYS = frozenset(("path", "sha256", "size", "mode", "role"))
_VERDICT_KEYS = frozenset((
    "schema", "outcome", "promoted", "root", "record", "seal", "source",
    "candidate", "build_closure", "controller_bindings_sha256",
    "inventory_file_count", "runtime_file_count", "runtime_discovery",
))
_RUNTIME_DISCOVERY_VERDICT_KEYS = frozenset((
    "schema", "transcript_sha256", "ldd_sha256", "interpreter_sha256",
    "file_row_count", "virtual_row_count", "replay_count",
))
_HEX40 = re.compile(r"[0-9a-f]{40}")
_HEX256 = re.compile(r"[0-9a-f]{64}")
_LDD_ADDRESS = r"\(0x[0-9A-Fa-f]+\)"
_LDD_MAPPED = re.compile(
    rf"^[ \t]*(?P<soname>[^ /\t=]+)[ \t]+=>[ \t]+"
    rf"(?P<path>/[^ \t]+)[ \t]+{_LDD_ADDRESS}[ \t]*$")
_LDD_DIRECT = re.compile(
    rf"^[ \t]*(?P<path>/[^ \t]+)[ \t]+{_LDD_ADDRESS}[ \t]*$")
_LDD_VIRTUAL = re.compile(
    rf"^[ \t]*(?P<soname>linux-(?:vdso(?:32|64)?|gate)\.so\.[0-9]+)"
    rf"[ \t]+{_LDD_ADDRESS}[ \t]*$")


class CandidateAuthorityError(RuntimeError):
    """Candidate authority bytes or their sealed tree are invalid."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CandidateAuthorityError(message)


def _exact_mapping(
    value: Any, keys: frozenset[str], label: str,
) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys,
             f"{label} is not an exact object")
    return value


def _hex_digest(value: Any, label: str) -> str:
    _require(type(value) is str and _HEX256.fullmatch(value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_document(data: bytes, label: str) -> dict[str, Any]:
    try:
        value = contract.strict_json_loads(data, label)
        _require(type(value) is dict and
                 data == contract.canonical_json_bytes(value),
                 f"{label} bytes are not one canonical object")
        return value
    except CandidateAuthorityError:
        raise
    except Exception as error:
        raise CandidateAuthorityError(f"{label} is invalid") from error


def payload_identity(data: bytes) -> dict[str, Any]:
    """Return the canonical byte identity accepted by the record constructor."""
    _require(type(data) is bytes and 0 < len(data) <=
             MAX_AUTHORITY_PAYLOAD_BYTES,
             "candidate payload is not bounded non-empty bytes")
    return {"size": len(data), "sha256": _sha256(data)}


def _payload_maximum(path: str) -> int:
    if path == "artifacts/bench_leopard2":
        return MAX_EXECUTABLE_BYTES
    if path in _JSON_PAYLOAD_PATHS:
        return contract.MAX_JSON_BYTES
    return MAX_AUTHORITY_PAYLOAD_BYTES


def normalize_runtime_discovery(raw: bytes) -> tuple[dict[str, Any], ...]:
    """Normalize a bounded C-locale ldd transcript without losing rows."""
    _require(type(raw) is bytes and 0 < len(raw) <= MAX_LDD_OUTPUT_BYTES,
             "runtime discovery output is not bounded non-empty bytes")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as error:
        raise CandidateAuthorityError(
            "runtime discovery output is not strict ASCII") from error
    _require(text.endswith("\n") and "\r" not in text and "\0" not in text,
             "runtime discovery output is not LF-framed")
    lines = text[:-1].split("\n")
    _require(0 < len(lines) <= MAX_RUNTIME_CLOSURE_FILES and
             all(line.strip(" \t") for line in lines),
             "runtime discovery row count differs")
    rows: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        match = _LDD_MAPPED.fullmatch(line)
        if match is not None:
            rows.append({
                "soname": match.group("soname"), "kind": "file",
                "path": match.group("path"),
            })
            continue
        match = _LDD_DIRECT.fullmatch(line)
        if match is not None:
            path = match.group("path")
            rows.append({
                "soname": path.rsplit("/", 1)[-1], "kind": "file",
                "path": path,
            })
            continue
        match = _LDD_VIRTUAL.fullmatch(line)
        if match is not None:
            rows.append({
                "soname": match.group("soname"), "kind": "virtual",
                "path": None,
            })
            continue
        raise CandidateAuthorityError(
            f"runtime discovery row {index} is unsupported")
    rows.sort(key=lambda item: item["soname"])
    sonames = [item["soname"] for item in rows]
    paths = [item["path"] for item in rows if item["kind"] == "file"]
    _require(
        sonames == sorted(set(sonames)) and
        len(paths) == len(set(paths)) and
        all(type(path) is str and path.startswith("/") and
            os.path.normpath(path) == path and "\0" not in path
            for path in paths),
        "runtime discovery rows are not canonical and unique")
    return tuple(copy.deepcopy(rows))


def runtime_discovery_sha256(
    rows: Sequence[Mapping[str, Any]],
) -> str:
    """Hash the complete normalized transcript retained in the role field."""
    _require(type(rows) in (list, tuple) and bool(rows),
             "runtime discovery transcript is not a sequence")
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(rows):
        row = _exact_mapping(
            raw, frozenset(("soname", "kind", "path")),
            f"runtime discovery row {index}")
        soname = row["soname"]
        _require(type(soname) is str and 0 < len(soname) <= 255 and
                 "/" not in soname and "\t" not in soname and
                 all(0x21 <= ord(character) <= 0x7e
                     for character in soname),
                 f"runtime discovery row {index} soname differs")
        if row["kind"] == "virtual":
            _require(row["path"] is None,
                     f"runtime discovery row {index} virtual path differs")
            normalized.append({
                "soname": soname, "kind": "virtual", "path": None})
        else:
            path = row["path"]
            _require(row["kind"] == "file" and type(path) is str and
                     path.startswith("/") and os.path.normpath(path) == path and
                     "\0" not in path,
                     f"runtime discovery row {index} file path differs")
            normalized.append({
                "soname": soname, "kind": "file", "path": path})
    _require(normalized == sorted(normalized, key=lambda item: item["soname"]) and
             len({item["soname"] for item in normalized}) == len(normalized) and
             len({item["path"] for item in normalized
                  if item["kind"] == "file"}) ==
                sum(item["kind"] == "file" for item in normalized),
             "runtime discovery transcript is not canonical and unique")
    return contract.canonical_sha256({
        "schema": RUNTIME_DISCOVERY_SCHEMA, "rows": normalized})


def validate_runtime_closure(value: Any) -> dict[str, Any]:
    """Validate the runtime closure identically to the live Gen3 consumer."""
    record = _exact_mapping(
        value, frozenset(("schema", "dependencies", "launchers")),
        "candidate runtime closure")
    _require(
        record["schema"] == RUNTIME_CLOSURE_SCHEMA and
        type(record["dependencies"]) is list and
        type(record["launchers"]) is list and
        len(record["launchers"]) == 2 and
        2 <= len(record["dependencies"]) + len(record["launchers"]) <=
            MAX_RUNTIME_CLOSURE_FILES,
        "candidate runtime closure metadata differs")
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(
            record["dependencies"] + record["launchers"]):
        item = _exact_mapping(
            raw, _RUNTIME_ENTRY_KEYS,
            f"candidate runtime entry {index}")
        path = Path(str(item["path"]))
        _require(
            path.is_absolute() and "\0" not in str(path) and
            type(item["size"]) is int and item["size"] > 0 and
            type(item["mode"]) is int and 0 <= item["mode"] <= 0o7777 and
            type(item["role"]) is str and bool(item["role"]),
            f"candidate runtime entry {index} is malformed")
        normalized.append({
            "path": str(path),
            "sha256": _hex_digest(
                item["sha256"], f"candidate runtime entry {index} hash"),
            "size": item["size"],
            "mode": item["mode"],
            "role": item["role"],
        })
    split = len(record["dependencies"])
    dependencies = normalized[:split]
    launchers = normalized[split:]
    discovery = [
        item for item in dependencies
        if item["role"].startswith(LDD_DISCOVERY_ROLE_PREFIX)
    ]
    interpreters = [
        item for item in dependencies
        if item["role"] == LDD_INTERPRETER_ROLE
    ]
    ordinary = [
        item for item in dependencies if item not in discovery + interpreters
    ]
    _require(
        normalized == sorted(
            normalized, key=lambda item: (item["role"], item["path"])) and
        len({item["path"] for item in normalized}) == len(normalized) and
        len({(item["role"], item["path"]) for item in normalized}) ==
            len(normalized) and
        len(discovery) == 1 and discovery[0]["path"] == LDD_PATH and
        _HEX256.fullmatch(
            discovery[0]["role"][len(LDD_DISCOVERY_ROLE_PREFIX):]) is not None and
        len(interpreters) == 1 and
        interpreters[0]["path"] == LDD_INTERPRETER_PATH and
        all(item["role"] == "dependency" for item in ordinary) and
        {item["path"] for item in launchers} ==
            {"/usr/bin/prlimit", "/usr/bin/taskset"} and
        all(item["role"] == "launcher" for item in launchers),
        "candidate runtime closure is not canonical, bounded, and unique")
    return {
        "schema": record["schema"],
        "dependencies": dependencies,
        "launchers": launchers,
    }


def validate_candidate_authority_record(value: Any) -> dict[str, Any]:
    """Validate the exact terminal shape and every record/inventory join."""
    record = _exact_mapping(
        value, _RECORD_KEYS, "candidate authority record")
    source = _exact_mapping(
        record["source"], _SOURCE_KEYS, "candidate authority source")
    candidate = _exact_mapping(
        record["candidate"], _CANDIDATE_KEYS,
        "candidate authority executable")
    build = _exact_mapping(
        record["build_closure"], _BUILD_KEYS,
        "candidate authority build closure")
    _require(
        record["schema"] == CANDIDATE_AUTHORITY_SCHEMA and
        type(record["generation"]) is int and
        record["generation"] == GENERATION and
        type(record["status"]) is str and
        record["status"] == "promoted_authority" and
        type(source["commit"]) is str and
        _HEX40.fullmatch(source["commit"]) is not None and
        type(source["tree"]) is str and
        _HEX40.fullmatch(source["tree"]) is not None and
        candidate["relative_path"] == "artifacts/bench_leopard2" and
        type(candidate["size"]) is int and
        0 < candidate["size"] <= MAX_EXECUTABLE_BYTES,
        "candidate authority identity fields differ")
    for label, digest in (
            ("candidate executable", candidate["sha256"]),
            ("Git capture", source["git_capture_sha256"]),
            ("source archive", source["source_archive_sha256"]),
            ("build provenance", build["build_provenance_sha256"]),
            ("reproducible-build core",
             build["reproducible_build_core_sha256"]),
            ("reproducible-build proof",
             build["reproducible_build_proof_sha256"]),
            ("attestation header", build["attestation_header_sha256"]),
            ("runtime closure", build["runtime_closure_sha256"]),
            ("controller bindings", record["controller_bindings_sha256"])):
        _hex_digest(digest, label)
    inventory = record["inventory"]
    _require(type(inventory) is list and len(inventory) == len(PAYLOAD_PATHS),
             "candidate authority inventory length differs")
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(inventory):
        item = _exact_mapping(
            raw, _INVENTORY_KEYS,
            f"candidate authority inventory {index}")
        _require(
            type(item["relative_path"]) is str and
            item["relative_path"] in PAYLOAD_PATHS and
            type(item["size"]) is int and item["size"] > 0 and
            item["size"] <= _payload_maximum(item["relative_path"]),
            f"candidate authority inventory {index} differs")
        normalized.append({
            "relative_path": item["relative_path"],
            "sha256": _hex_digest(
                item["sha256"],
                f"candidate authority inventory {index} hash"),
            "size": item["size"],
        })
    normalized.sort(key=lambda item: item["relative_path"])
    _require(
        [item["relative_path"] for item in normalized] ==
            list(PAYLOAD_PATHS) and inventory == normalized,
        "candidate authority inventory is not exact and canonical")
    by_path = {item["relative_path"]: item for item in normalized}
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
    _require(all(by_path[path]["sha256"] == digest
                 for path, digest in joins.items()),
             "candidate authority record-to-inventory joins differ")
    _require(by_path["artifacts/bench_leopard2"]["size"] ==
             candidate["size"],
             "candidate authority executable size join differs")
    return copy.deepcopy(record)


def candidate_authority_record(
    *,
    source_commit: str,
    source_tree: str,
    controller_bindings_sha256: str,
    payload_identities: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Construct the terminal solely from the exact payload identities."""
    _require(type(payload_identities) is dict and
             set(payload_identities) == set(PAYLOAD_PATHS),
             "candidate payload identity set differs")
    normalized: dict[str, dict[str, Any]] = {}
    for path in PAYLOAD_PATHS:
        identity = _exact_mapping(
            payload_identities[path], frozenset(("size", "sha256")),
            f"candidate payload identity {path}")
        _require(type(identity["size"]) is int and identity["size"] > 0,
                 f"candidate payload identity {path} has invalid size")
        _require(identity["size"] <= _payload_maximum(path),
                 f"candidate payload identity {path} exceeds its bound")
        normalized[path] = {
            "size": identity["size"],
            "sha256": _hex_digest(
                identity["sha256"], f"candidate payload identity {path}"),
        }
    inventory = [
        {"relative_path": path, **normalized[path]}
        for path in PAYLOAD_PATHS
    ]
    record = {
        "schema": CANDIDATE_AUTHORITY_SCHEMA,
        "generation": GENERATION,
        "status": "promoted_authority",
        "source": {
            "commit": source_commit,
            "tree": source_tree,
            "git_capture_sha256":
                normalized["source/git-capture.json"]["sha256"],
            "source_archive_sha256":
                normalized["source/candidate-source.tar"]["sha256"],
        },
        "candidate": {
            "relative_path": "artifacts/bench_leopard2",
            "sha256": normalized["artifacts/bench_leopard2"]["sha256"],
            "size": normalized["artifacts/bench_leopard2"]["size"],
        },
        "build_closure": {
            "build_provenance_sha256":
                normalized["build/build-provenance.json"]["sha256"],
            "reproducible_build_core_sha256":
                normalized["build/reproducible-build-core.json"]["sha256"],
            "reproducible_build_proof_sha256":
                normalized["build/reproducible-build-proof.json"]["sha256"],
            "attestation_header_sha256":
                normalized["build/benchmark-source-attestation.h"]["sha256"],
            "runtime_closure_sha256":
                normalized["runtime/runtime-closure.json"]["sha256"],
        },
        "controller_bindings_sha256": controller_bindings_sha256,
        "inventory": inventory,
    }
    return validate_candidate_authority_record(record)


class _ArchiveBytes:
    """Minimal byte-range surface consumed by the audited tar verifier."""

    def __init__(self, data: bytes) -> None:
        self.data = data

    def size(self, name: str) -> int:
        _require(name == "candidate-source.tar",
                 "candidate source archive name differs")
        return len(self.data)

    def pread(self, name: str, offset: int, size: int) -> bytes:
        _require(
            name == "candidate-source.tar" and type(offset) is int and
            type(size) is int and offset >= 0 and size >= 0 and
            offset + size <= len(self.data),
            "candidate source archive read is out of bounds")
        return self.data[offset:offset + size]


def controller_bindings_sha256(source_root: Path | str) -> str:
    """Hash the live armer's exact reviewed controller binding list."""
    root = Path(source_root).resolve(strict=True)
    try:
        bindings = prereg.current_controller_bindings(root)
        return contract.canonical_sha256(bindings)
    except Exception as error:
        raise CandidateAuthorityError(
            "current controller bindings could not be captured") from error


def validate_candidate_authority_verdict(
    value: Any,
    *,
    expected_root: Path | str | None = None,
    expected_record_sha256: str | None = None,
    expected_ledger_sha256: str | None = None,
    expected_controller_bindings_sha256: str | None = None,
) -> dict[str, Any]:
    """Validate the exact fresh-process receipt consumed by the producer."""
    verdict = _exact_mapping(value, _VERDICT_KEYS, "candidate verdict")
    record = _exact_mapping(
        verdict["record"], frozenset(("relative_path", "schema", "sha256")),
        "candidate verdict record")
    seal = _exact_mapping(
        verdict["seal"], frozenset((
            "sha256sums_sha256", "checksum_line_count",
            "metadata_entry_count")), "candidate verdict seal")
    source = _exact_mapping(
        verdict["source"], frozenset(("commit", "tree", "detached")),
        "candidate verdict source")
    candidate = _exact_mapping(
        verdict["candidate"], _CANDIDATE_KEYS,
        "candidate verdict executable")
    build = _exact_mapping(
        verdict["build_closure"], _BUILD_KEYS,
        "candidate verdict build closure")
    discovery = _exact_mapping(
        verdict["runtime_discovery"], _RUNTIME_DISCOVERY_VERDICT_KEYS,
        "candidate verdict runtime discovery")
    try:
        verdict_root = Path(verdict["root"])
        canonical_verdict_root = verdict_root.resolve(strict=True)
    except (OSError, TypeError, ValueError) as error:
        raise CandidateAuthorityError(
            "candidate verdict lane root is unavailable") from error
    _require(
        verdict["schema"] == VERDICT_SCHEMA and
        verdict["outcome"] == "promoted_authority" and
        verdict["promoted"] is True and
        type(verdict["root"]) is str and verdict_root.is_absolute() and
        verdict_root == canonical_verdict_root and
        record["relative_path"] == TERMINAL_PATH and
        record["schema"] == CANDIDATE_AUTHORITY_SCHEMA and
        type(seal["checksum_line_count"]) is int and
        seal["checksum_line_count"] == len(SEALED_PATHS) - 1 and
        type(seal["metadata_entry_count"]) is int and
        seal["metadata_entry_count"] >= len(SEALED_PATHS) and
        type(source["commit"]) is str and
        _HEX40.fullmatch(source["commit"]) is not None and
        type(source["tree"]) is str and
        _HEX40.fullmatch(source["tree"]) is not None and
        type(source["detached"]) is bool and
        candidate["relative_path"] == "artifacts/bench_leopard2" and
        type(candidate["size"]) is int and
        0 < candidate["size"] <= MAX_EXECUTABLE_BYTES and
        verdict["inventory_file_count"] == len(PAYLOAD_PATHS) and
        type(verdict["runtime_file_count"]) is int and
        4 <= verdict["runtime_file_count"] <= MAX_RUNTIME_CLOSURE_FILES and
        discovery["schema"] == RUNTIME_DISCOVERY_SCHEMA and
        type(discovery["file_row_count"]) is int and
        1 <= discovery["file_row_count"] <= MAX_RUNTIME_CLOSURE_FILES and
        type(discovery["virtual_row_count"]) is int and
        0 <= discovery["virtual_row_count"] <= MAX_RUNTIME_CLOSURE_FILES and
        discovery["file_row_count"] + discovery["virtual_row_count"] <=
            MAX_RUNTIME_CLOSURE_FILES and
        discovery["replay_count"] == 2,
        "candidate verdict semantics differ")
    for label, digest in (
            ("record", record["sha256"]),
            ("ledger", seal["sha256sums_sha256"]),
            ("candidate", candidate["sha256"]),
            ("controller", verdict["controller_bindings_sha256"]),
            ("runtime transcript", discovery["transcript_sha256"]),
            ("runtime ldd", discovery["ldd_sha256"]),
            ("runtime interpreter", discovery["interpreter_sha256"]),
            *((label, digest) for label, digest in build.items())):
        _hex_digest(digest, f"candidate verdict {label}")
    if expected_root is not None:
        root = Path(os.path.abspath(os.fspath(expected_root))).resolve(strict=True)
        _require(verdict["root"] == str(root),
                 "candidate verdict names another lane root")
    for observed, expected, label in (
            (record["sha256"], expected_record_sha256, "record"),
            (seal["sha256sums_sha256"], expected_ledger_sha256, "ledger"),
            (verdict["controller_bindings_sha256"],
             expected_controller_bindings_sha256, "controller bindings")):
        if expected is not None:
            _hex_digest(expected, f"expected candidate verdict {label}")
            _require(observed == expected,
                     f"candidate verdict {label} differs")
    return copy.deepcopy(verdict)


def _verify_source_archive(
    data: bytes, git_record: Mapping[str, Any], commit: str,
) -> None:
    tracked_list = git_record.get("tracked_files")
    _require(type(tracked_list) is list and bool(tracked_list),
             "candidate Git capture has no tracked-file closure")
    tracked: dict[str, Mapping[str, Any]] = {}
    for index, item in enumerate(tracked_list):
        _require(type(item) is dict and type(item.get("path")) is str and
                 item["path"] not in tracked,
                 f"candidate tracked file {index} is malformed")
        tracked[item["path"]] = item
    try:
        members = sealed_verifier._verify_source_archive(
            _ArchiveBytes(data), "candidate-source.tar",
            expected_prefix="candidate-source/", expected_commit=commit)
        sealed_verifier._verify_archive_git_binding(
            members, tracked, "candidate source")
    except CandidateAuthorityError:
        raise
    except Exception as error:
        raise CandidateAuthorityError(
            "candidate source archive did not match its Git capture") \
            from error


def _verify_runtime_files(runtime: Mapping[str, Any]) -> None:
    """Rehash every current runtime path without executing the candidate."""
    for index, entry in enumerate(
            runtime["dependencies"] + runtime["launchers"]):
        raw = Path(entry["path"])
        descriptor = -1
        try:
            resolved = raw.resolve(strict=True)
            descriptor = os.open(
                resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
            before = os.fstat(descriptor)
            pathname = os.stat(resolved, follow_symlinks=False)
            _require(
                before.st_size == entry["size"] and
                before.st_size <= MAX_EXECUTABLE_BYTES and
                (before.st_dev, before.st_ino) ==
                    (pathname.st_dev, pathname.st_ino) and
                (before.st_mode & 0o170000) == 0o100000 and
                before.st_mode & 0o7777 == entry["mode"],
                f"candidate runtime entry {index} metadata differs")
            digest = hashlib.sha256()
            remaining = before.st_size
            while remaining:
                chunk = os.read(descriptor, min(1 << 20, remaining))
                _require(bool(chunk),
                         f"candidate runtime entry {index} ended while read")
                digest.update(chunk)
                remaining -= len(chunk)
            _require(os.read(descriptor, 1) == b"" and
                     digest.hexdigest() == entry["sha256"],
                     f"candidate runtime entry {index} bytes differ")
            after = os.fstat(descriptor)
            pathname_after = os.stat(resolved, follow_symlinks=False)
            fingerprint = (
                before.st_dev, before.st_ino, before.st_mode, before.st_uid,
                before.st_gid, before.st_nlink, before.st_size,
                before.st_mtime_ns, before.st_ctime_ns,
            )
            _require(
                (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
                 after.st_gid, after.st_nlink, after.st_size,
                 after.st_mtime_ns, after.st_ctime_ns) ==
                    fingerprint and
                (pathname_after.st_dev, pathname_after.st_ino,
                 pathname_after.st_mode, pathname_after.st_uid,
                 pathname_after.st_gid, pathname_after.st_nlink,
                 pathname_after.st_size, pathname_after.st_mtime_ns,
                 pathname_after.st_ctime_ns) == fingerprint and
                raw.resolve(strict=True) == resolved,
                f"candidate runtime entry {index} changed while verified")
        except CandidateAuthorityError:
            raise
        except OSError as error:
            raise CandidateAuthorityError(
                f"candidate runtime entry {index} could not be verified") \
                from error
        finally:
            if descriptor >= 0:
                os.close(descriptor)


def _retain_runtime_entry(
    entry: Mapping[str, Any], label: str,
) -> tuple[int, Path, Path, tuple[int, ...]]:
    """Open and hash one runtime program while retaining its exact inode."""
    raw = Path(str(entry["path"]))
    descriptor = -1
    try:
        resolved = raw.resolve(strict=True)
        descriptor = os.open(
            resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        before = os.fstat(descriptor)
        pathname = os.stat(resolved, follow_symlinks=False)
        fingerprint = (
            before.st_dev, before.st_ino, before.st_mode, before.st_uid,
            before.st_gid, before.st_nlink, before.st_size,
            before.st_mtime_ns, before.st_ctime_ns,
        )
        _require(
            stat.S_ISREG(before.st_mode) and
            before.st_size == entry["size"] and
            stat.S_IMODE(before.st_mode) == entry["mode"] and
            (pathname.st_dev, pathname.st_ino, pathname.st_mode,
             pathname.st_uid, pathname.st_gid, pathname.st_nlink,
             pathname.st_size, pathname.st_mtime_ns,
             pathname.st_ctime_ns) == fingerprint,
            f"{label} metadata differs")
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            chunk = os.pread(
                descriptor, min(1 << 20, before.st_size - offset), offset)
            _require(bool(chunk), f"{label} ended while retained")
            digest.update(chunk)
            offset += len(chunk)
        _require(digest.hexdigest() == entry["sha256"],
                 f"{label} bytes differ")
        retained = (descriptor, raw, resolved, fingerprint)
        _verify_retained_runtime_entry(retained, label)
        descriptor = -1
        return retained
    except CandidateAuthorityError:
        raise
    except OSError as error:
        raise CandidateAuthorityError(f"{label} could not be retained") \
            from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _verify_retained_runtime_entry(
    retained: tuple[int, Path, Path, tuple[int, ...]], label: str,
) -> None:
    descriptor, raw, resolved, fingerprint = retained
    try:
        after = os.fstat(descriptor)
        pathname = os.stat(resolved, follow_symlinks=False)
        _require(
            (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
             after.st_gid, after.st_nlink, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns) == fingerprint and
            (pathname.st_dev, pathname.st_ino, pathname.st_mode,
             pathname.st_uid, pathname.st_gid, pathname.st_nlink,
             pathname.st_size, pathname.st_mtime_ns,
             pathname.st_ctime_ns) == fingerprint and
            raw.resolve(strict=True) == resolved,
            f"{label} changed while retained")
    except CandidateAuthorityError:
        raise
    except OSError as error:
        raise CandidateAuthorityError(f"{label} could not be rechecked") \
            from error


def _bounded_ldd_output(
    executable_descriptor: int,
    ldd_descriptor: int,
    interpreter_descriptor: int,
) -> bytes:
    """Run the sealed ldd tool on one inherited immutable executable fd."""
    descriptors = (
        executable_descriptor, ldd_descriptor, interpreter_descriptor)
    _require(all(type(item) is int and item >= 0 for item in descriptors) and
             len(set(descriptors)) == len(descriptors),
             "runtime discovery descriptors are invalid")
    process: subprocess.Popen[bytes] | None = None
    selector = selectors.DefaultSelector()
    output = bytearray()
    diagnostics = bytearray()
    failure: str | None = None
    try:
        process = subprocess.Popen(
            [f"/proc/self/fd/{interpreter_descriptor}",
             f"/proc/self/fd/{ldd_descriptor}",
             f"/proc/self/fd/{executable_descriptor}"],
            cwd=str(REPO_ROOT),
            env={
                "HOME": "/nonexistent", "LANG": "C", "LC_ALL": "C",
                "PATH": "/usr/bin:/bin", "TZ": "UTC",
            },
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, close_fds=True,
            pass_fds=descriptors, start_new_session=True)
        _require(process.stdout is not None and process.stderr is not None,
                 "runtime discovery pipes were not created")
        streams = {
            process.stdout.fileno(): output,
            process.stderr.fileno(): diagnostics,
        }
        for stream in (process.stdout, process.stderr):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ)
        deadline = time.monotonic() + LDD_TIMEOUT_SECONDS
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = "runtime discovery exceeded its timeout"
                break
            for key, unused_event in selector.select(min(remaining, 0.1)):
                descriptor = key.fileobj.fileno()
                try:
                    chunk = os.read(descriptor, 1 << 16)
                except BlockingIOError:
                    continue
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                streams[descriptor].extend(chunk)
                if len(output) + len(diagnostics) > MAX_LDD_OUTPUT_BYTES:
                    failure = "runtime discovery exceeded its output bound"
                    break
            if failure is not None:
                break
        if failure is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = "runtime discovery exceeded its timeout"
            else:
                try:
                    status = process.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    failure = "runtime discovery exceeded its timeout"
                else:
                    _require(status == 0 and diagnostics == b"" and output,
                             "runtime discovery failed or emitted diagnostics")
        if failure is not None:
            raise CandidateAuthorityError(failure)
        return bytes(output)
    except CandidateAuthorityError:
        raise
    except (OSError, ValueError) as error:
        raise CandidateAuthorityError(
            "runtime discovery process could not run") from error
    finally:
        selector.close()
        if process is not None:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=5)
                except (OSError, subprocess.TimeoutExpired):
                    pass
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()


def _verify_runtime_discovery(
    runtime: Mapping[str, Any], artifact: bytes,
) -> dict[str, Any]:
    """Replay the sealed transcript twice without entering benchmark mode."""
    dependencies = runtime["dependencies"]
    discovery = next(
        item for item in dependencies
        if item["role"].startswith(LDD_DISCOVERY_ROLE_PREFIX))
    interpreter = next(
        item for item in dependencies
        if item["role"] == LDD_INTERPRETER_ROLE)
    expected_digest = discovery["role"][len(LDD_DISCOVERY_ROLE_PREFIX):]
    expected_paths = sorted(
        item["path"] for item in dependencies
        if item["role"] == "dependency")
    descriptor = -1
    ldd_retained: tuple[int, Path, Path, tuple[int, ...]] | None = None
    interpreter_retained: \
        tuple[int, Path, Path, tuple[int, ...]] | None = None
    try:
        ldd_retained = _retain_runtime_entry(
            discovery, "runtime discovery tool")
        interpreter_retained = _retain_runtime_entry(
            interpreter, "runtime discovery interpreter")
        flags = os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING
        descriptor = os.memfd_create("leopard-k65-runtime-discovery", flags)
        offset = 0
        while offset < len(artifact):
            written = os.write(descriptor, artifact[offset:offset + (1 << 20)])
            _require(written > 0,
                     "runtime discovery artifact write made no progress")
            offset += written
        os.fchmod(descriptor, 0o500)
        seals = (fcntl.F_SEAL_SEAL | fcntl.F_SEAL_SHRINK |
                 fcntl.F_SEAL_GROW | fcntl.F_SEAL_WRITE)
        fcntl.fcntl(descriptor, fcntl.F_ADD_SEALS, seals)
        first = normalize_runtime_discovery(_bounded_ldd_output(
            descriptor, ldd_retained[0], interpreter_retained[0]))
        second = normalize_runtime_discovery(_bounded_ldd_output(
            descriptor, ldd_retained[0], interpreter_retained[0]))
        _verify_retained_runtime_entry(
            ldd_retained, "runtime discovery tool")
        _verify_retained_runtime_entry(
            interpreter_retained, "runtime discovery interpreter")
    except CandidateAuthorityError:
        raise
    except (AttributeError, OSError) as error:
        raise CandidateAuthorityError(
            "immutable runtime discovery artifact could not be prepared") \
            from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if ldd_retained is not None:
            os.close(ldd_retained[0])
        if interpreter_retained is not None:
            os.close(interpreter_retained[0])
    _require(
        contract.exact_json_equal(list(first), list(second)) and
        runtime_discovery_sha256(first) == expected_digest and
        sorted(item["path"] for item in first
               if item["kind"] == "file") == expected_paths,
        "independent runtime discovery differs from its sealed transcript")
    return {
        "schema": RUNTIME_DISCOVERY_SCHEMA,
        "transcript_sha256": expected_digest,
        "ldd_sha256": discovery["sha256"],
        "interpreter_sha256": interpreter["sha256"],
        "file_row_count": sum(item["kind"] == "file" for item in first),
        "virtual_row_count": sum(
            item["kind"] == "virtual" for item in first),
        "replay_count": 2,
    }


def verify_candidate_authority_lane(
    lane_root: Path | str,
    source_root: Path | str,
    *,
    expected_controller_bindings_sha256: str | None = None,
    require_detached: bool = False,
    verify_current_source: bool = False,
) -> dict[str, Any]:
    """Replay one sealed lane without entering the candidate benchmark."""
    requested_root = Path(os.path.abspath(os.fspath(lane_root)))
    source_path = Path(source_root).resolve(strict=True)
    _require(requested_root.is_absolute() and not requested_root.is_symlink(),
             "candidate authority lane is not one direct absolute path")
    root = requested_root.resolve(strict=True)
    _require(root == requested_root and source_path == REPO_ROOT,
             "candidate authority lane/source path is not canonical here")
    _require(type(require_detached) is bool and
             type(verify_current_source) is bool,
             "candidate source replay policy is invalid")
    if expected_controller_bindings_sha256 is not None:
        _hex_digest(
            expected_controller_bindings_sha256,
            "expected controller bindings")
    try:
        with sealed_verifier.read_sealed_tree(root) as tree:
            metadata = sealed_verifier.verify_tree_metadata(tree)
            ledger = sealed_verifier.verify_sha256sums(tree)
            _require(set(tree.files) == SEALED_PATHS,
                     "candidate authority file set differs")
            record_data = tree.read_file(
                TERMINAL_PATH, maximum_bytes=MAX_JSON_BYTES)
            record = validate_candidate_authority_record(
                _canonical_document(record_data, "candidate authority record"))
            inventory = record["inventory"]
            payload: dict[str, bytes] = {}
            for item in inventory:
                path = item["relative_path"]
                maximum = _payload_maximum(path)
                data = tree.read_file(path, maximum_bytes=maximum)
                _require(len(data) == item["size"] and
                         _sha256(data) == item["sha256"],
                         f"candidate authority payload differs: {path}")
                payload[path] = data
            tree.reverify()
    except CandidateAuthorityError:
        raise
    except Exception as error:
        raise CandidateAuthorityError(
            "sealed candidate authority did not replay") from error

    current_controller_sha256 = controller_bindings_sha256(source_path)
    _require(record["controller_bindings_sha256"] ==
             current_controller_sha256,
             "candidate authority differs from current controller bindings")
    if expected_controller_bindings_sha256 is not None:
        _require(record["controller_bindings_sha256"] ==
                 expected_controller_bindings_sha256,
                 "candidate authority controller bindings differ")
    git_record = _canonical_document(
        payload["source/git-capture.json"], "candidate Git capture")
    try:
        validated_git = git_capture.validate_git_capture(
            git_record, str(source_path), record["source"]["commit"],
            require_detached=require_detached)
    except Exception as error:
        raise CandidateAuthorityError(
            "candidate Git capture did not replay") from error
    _require(validated_git["tree"] == record["source"]["tree"],
             "candidate Git capture tree differs")
    if verify_current_source:
        try:
            current = git_capture.capture_git_identity(
                source_path, record["source"]["commit"],
                require_detached=require_detached)
        except Exception as error:
            raise CandidateAuthorityError(
                "current candidate source could not be recaptured") from error
        _require(contract.exact_json_equal(current, validated_git),
                 "current candidate source differs from its sealed capture")
    _verify_source_archive(
        payload["source/candidate-source.tar"], validated_git,
        record["source"]["commit"])

    provenance = _canonical_document(
        payload["build/build-provenance.json"],
        "candidate build provenance")
    proof = _canonical_document(
        payload["build/reproducible-build-proof.json"],
        "candidate reproducible-build proof")
    core = _canonical_document(
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
        raise CandidateAuthorityError(
            "candidate reproducible-build proof did not replay") from error
    candidate = record["candidate"]
    manifest_git = provenance.get("tracked_source_manifest", {}).get("git", {})
    _require(
        contract.exact_json_equal(core, expected_core) and
        manifest_git.get("commit") == record["source"]["commit"] and
        manifest_git.get("tree") == record["source"]["tree"] and
        manifest_git.get("dirty") is False and
        provenance.get("executable", {}).get("sha256") ==
            candidate["sha256"] and
        provenance.get("executable", {}).get("size") == candidate["size"] and
        payload["build/benchmark-source-attestation.h"] == expected_header,
        "candidate source/build/executable closure differs")
    runtime = validate_runtime_closure(_canonical_document(
        payload["runtime/runtime-closure.json"],
        "candidate runtime closure"))
    _verify_runtime_files(runtime)
    artifact = payload["artifacts/bench_leopard2"]
    _require(artifact.startswith(b"\x7fELF"),
             "candidate authority artifact is not ELF")
    runtime_discovery = _verify_runtime_discovery(runtime, artifact)
    _verify_runtime_files(runtime)
    return {
        "schema": VERDICT_SCHEMA,
        "outcome": "promoted_authority",
        "promoted": True,
        "root": str(root),
        "record": {
            "relative_path": TERMINAL_PATH,
            "schema": record["schema"],
            "sha256": _sha256(record_data),
        },
        "seal": {
            "sha256sums_sha256": ledger["sha256"],
            "checksum_line_count": ledger["line_count"],
            "metadata_entry_count": len(metadata["entries"]),
        },
        "source": {
            "commit": record["source"]["commit"],
            "tree": record["source"]["tree"],
            "detached": validated_git["detached"],
        },
        "candidate": copy.deepcopy(candidate),
        "build_closure": copy.deepcopy(record["build_closure"]),
        "controller_bindings_sha256":
            record["controller_bindings_sha256"],
        "inventory_file_count": len(record["inventory"]),
        "runtime_file_count":
            len(runtime["dependencies"]) + len(runtime["launchers"]),
        "runtime_discovery": runtime_discovery,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="independently replay a sealed K65 Gen3 candidate lane")
    parser.add_argument("lane", type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--controller-bindings-sha256", required=True)
    arguments = parser.parse_args(argv)
    try:
        verdict = verify_candidate_authority_lane(
            arguments.lane, arguments.source_root,
            expected_controller_bindings_sha256=
                arguments.controller_bindings_sha256,
            require_detached=True, verify_current_source=True)
        sys.stdout.buffer.write(contract.canonical_json_bytes(verdict))
        return 0
    except Exception as error:
        message = " ".join(str(error).replace("\r", " ").replace(
            "\n", " ").split())[:4096]
        sys.stderr.write((message or type(error).__name__) + "\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
