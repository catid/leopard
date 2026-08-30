#!/usr/bin/env python3
"""Pure authority records for the reproducible exact-main baseline.

This module validates already-acquired JSON facts.  It performs no file I/O,
host discovery, process execution, or artifact parsing.  The acquisition and
offline-verification layers are separate and must rederive every file-backed
claim before a record may be sealed or trusted.  Failure producers must
normalize signal deaths to shell status ``128 + signal`` and flatten diagnostic
messages to one display-safe line before constructing a failure record.
"""

from __future__ import annotations

import copy
from datetime import datetime, timezone
import re
from typing import Any, Mapping, NoReturn, Sequence

from leopard2_exact_main_baseline import (
    EMPTY_CONTENT_SHA256,
    ExactMainBaselineError,
    canonical_json_bytes,
    canonical_sha256,
    exact_json_equal,
    strict_json_loads,
    validate_normalized_code_identity,
)


AUTHORITY_SCHEMA = "leopard2-gf8-exact-main-pure-avx2-baseline/v1"
ACQUISITION_FAILURE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline-acquisition-failure/v1"
VERIFICATION_FAILURE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline-verification-failure/v1"
VERIFIER_VERDICT_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline-verification/v1"
BUILD_PROFILE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-build-profile/v1"
RUNTIME_CLOSURE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-runtime-closure/v1"
ATTESTATION_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-attestation/v1"
PROMOTION_GATE = \
    "same-path-bytes-path-variant-normalized-zero-census/v1"
SEAL_PROTOCOL = "owner-only-tree-sha256sums/v1"
CANONICAL_LDD_NORMALIZATION = "canonical-ldd-C-v1"
BUILD_CLOSURE_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-build-closure/v1"
BENCHMARK_SCHEMA = "leopard-main-benchmark-v1"
CTEST_SUMMARY_LINE = "100% tests passed, 0 tests failed out of 1"

BEAD_ID = "leopard-79h.38.5.12.1"
BASELINE_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
BASELINE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
BASELINE_SSE2NEON_COMMIT = "cad518a93b326f0f644b7972d488d04eaa2b0475"
MINIMUM_HARNESS_COMMIT = \
    "f909ff767802c48640951d13b1d2597292c57fa5"
HISTORICAL_EXECUTABLE_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"
HISTORICAL_ARCHIVE_SHA256 = \
    "49250fb7528df955e57e96f84b309cf0c26f203122d112f3f501336c0eb131d0"

MAX_TEXT_LENGTH = 4096
MAX_ARGUMENTS = 128
MAX_DEPENDENCIES = 128
MAX_RETAINED_FILES = 256
MAX_CLOSURE_FILES = 1 << 16
MAX_FILE_BYTES = (1 << 63) - 1
MAX_CPU_COUNT = 4096
MAX_CANONICAL_LDD_BYTES = 1 << 20

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_OID = re.compile(r"[0-9a-f]{40}")
_TIMESTAMP = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z")
_ROLE = re.compile(r"[a-z][a-z0-9_-]{0,127}")

BUILD_ROLES = ("canonical_first", "canonical_second", "path_variant")
STAGE_NAMES = (
    "source_acquisition",
    "canonical_first_build",
    "canonical_second_build",
    "path_variant_build",
    "independent_verification",
    "seal",
)
FAILURE_STAGES = STAGE_NAMES[:-1]
TOOL_ROLES = (
    "archiver", "cmake", "compiler", "ctest", "git", "ldd", "linker",
    "make", "python", "ranlib",
)
SUBTOOL_ROLES = ("assembler", "cc1plus", "collect2")
VERSION_ROLES = tuple(sorted(TOOL_ROLES + SUBTOOL_ROLES))
ADAPTER_PATHS = (
    "experiments/leopard2/main_compare/CMakeLists.txt",
    "experiments/leopard2/main_compare/legacy_main_benchmark.cpp",
    "experiments/leopard2/main_compare/test_legacy_main_benchmark.py",
)
BASELINE_CPP_SOURCES = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)

AUTHORITY_KEYS = frozenset((
    "schema", "created_utc", "bead", "host", "lane", "source",
    "adapter", "toolchain", "configure", "builds", "runtime_closure",
    "attestation", "identity", "promotion", "superseded_references",
    "record_sha256",
))
HOST_KEYS = frozenset((
    "hostname", "kernel", "architecture", "cpu_model", "online_cpus",
    "sc_clk_tck",
))
LANE_KEYS = frozenset((
    "root", "attempt", "attempt_budget", "record_relative_path",
    "seal_protocol", "stages",
))
STAGE_KEYS = frozenset(("name", "status", "log"))
SOURCE_KEYS = frozenset(("baseline", "adapter_repository"))
BASELINE_SOURCE_KEYS = frozenset((
    "commit", "tree", "submodule", "git_capture", "archive",
))
ADAPTER_SOURCE_KEYS = frozenset((
    "commit", "tree", "clean", "git_capture", "archive",
))
SUBMODULE_KEYS = frozenset(("path", "commit"))
ARCHIVE_KEYS = frozenset((
    "relative_path", "prefix", "size", "sha256", "replay_sha256",
    "replay_identical",
))
FILE_KEYS = frozenset(("relative_path", "size", "sha256"))
TEXT_KEYS = frozenset(("relative_path", "size", "sha256"))
ADAPTER_KEYS = frozenset((
    "minimum_harness_commit", "files", "files_sha256",
))
ADAPTER_FILE_KEYS = frozenset((
    "path", "git_blob_sha1", "size", "sha256",
))
TOOLCHAIN_KEYS = frozenset(("tools", "subtools", "versions"))
TOOL_KEYS = frozenset((
    "role", "path", "resolved_path", "size", "mode", "sha256",
))
VERSION_KEYS = frozenset((
    "role", "argv", "stdout", "stderr", "exit_status",
))
PROFILE_KEYS = frozenset((
    "schema", "name", "generator", "build_type", "language",
    "release_flags", "isa_flags", "warning_flags", "openmp",
    "pure_avx2", "adapter_definitions", "environment",
    "cache_requirements", "build_job_count",
))
BUILDS_KEYS = frozenset(BUILD_ROLES)
BUILD_KEYS = frozenset((
    "role", "roots", "configure_argv", "build_argv", "configure_log",
    "build_log", "cmake_cache", "compile_commands", "executable",
    "archive", "closure",
))
ROOTS_KEYS = frozenset((
    "adapter_source_root", "baseline_source_root", "build_root",
))
BUILD_ARTIFACT_KEYS = frozenset((
    "name", "build_relative_path", "retained_relative_path", "size",
    "sha256",
))
CLOSURE_KEYS = frozenset((
    "relative_path", "size", "sha256", "file_count",
))
BUILD_CLOSURE_KEYS = frozenset((
    "schema", "role", "build_root", "files", "file_count",
))
CLOSURE_FILE_KEYS = frozenset(("relative_path", "size", "sha256"))
RUNTIME_KEYS = frozenset(("schema", "normalization", "records"))
RUNTIME_RECORD_KEYS = frozenset((
    "role", "executable_sha256", "canonical_ldd_output", "dependencies",
))
DEPENDENCY_KEYS = frozenset((
    "soname", "kind", "path", "size", "sha256",
))
ATTESTATION_KEYS = frozenset(("schema", "test_controller", "records"))
ATTESTATION_RECORD_KEYS = frozenset((
    "role", "argv", "stdout", "stderr", "exit_status",
    "reported_schema", "main_source_commit", "pure_avx2", "round_trip",
    "ctest",
))
CTEST_KEYS = frozenset((
    "argv", "stdout", "stderr", "exit_status", "passed", "failed",
))
IDENTITY_KEYS = frozenset((
    "canonical_first", "canonical_second", "path_variant",
    "combined_sha256", "canonical_raw_bytes_identical",
    "path_variant_raw_bytes_differ", "normalized_match",
))
PROMOTION_KEYS = frozenset((
    "gate", "same_path_executable_bytes_identical",
    "same_path_archive_bytes_identical",
    "path_variant_raw_executable_differs", "path_variant_normalized_match",
    "selected_section_census_zero", "pure_avx2_attested",
    "historical_references_non_authoritative", "promoted",
))
REFERENCE_KEYS = frozenset(("label", "value", "authority"))
FAILURE_KEYS = frozenset((
    "schema", "created_utc", "bead", "lane", "stage", "error",
    "retained_files", "authority_record_sha256", "promoted",
    "record_sha256",
))
ERROR_KEYS = frozenset(("kind", "message", "exit_status"))


class ExactMainBaselineRecordError(ExactMainBaselineError):
    """An exact-main authority or failure record violates its contract."""


def _fail(message: str) -> NoReturn:
    raise ExactMainBaselineRecordError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _dict(value: Any, keys: frozenset[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label} is not an exact JSON object")
    _require(set(value) == keys, f"{label} has an unexpected key set")
    return value


def _integer(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int, f"{label} is not an exact integer")
    _require(minimum <= value <= maximum,
             f"{label} is outside its structural bound")
    return value


def _text(value: Any, label: str) -> str:
    _require(type(value) is str and 0 < len(value) <= MAX_TEXT_LENGTH,
             f"{label} is not a bounded non-empty string")
    _require(all(
        ord(character) >= 0x20 and
        not (0x7F <= ord(character) <= 0x9F) and
        not (0xD800 <= ord(character) <= 0xDFFF)
        for character in value
    ), f"{label} contains a non-display-safe character")
    return value


def _token(value: Any, label: str) -> str:
    token = _text(value, label)
    _require(_ROLE.fullmatch(token) is not None,
             f"{label} is not a canonical token")
    return token


def _sha256(value: Any, label: str) -> str:
    _require(type(value) is str and _SHA256.fullmatch(value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _git_oid(value: Any, label: str) -> str:
    _require(type(value) is str and _GIT_OID.fullmatch(value) is not None,
             f"{label} is not a lowercase Git object ID")
    return value


def _byte_identity(size_value: Any, digest_value: Any, label: str,
                   *, minimum_size: int) -> tuple[int, str]:
    size = _integer(size_value, minimum_size, MAX_FILE_BYTES, f"{label} size")
    digest = _sha256(digest_value, f"{label} SHA-256")
    _require((size == 0) == (digest == EMPTY_CONTENT_SHA256),
             f"{label} size and empty-content SHA-256 disagree")
    return size, digest


def _timestamp(value: Any, label: str) -> str:
    _require(type(value) is str and _TIMESTAMP.fullmatch(value) is not None,
             f"{label} is not a canonical UTC timestamp")
    try:
        parsed = datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as error:
        raise ExactMainBaselineRecordError(
            f"{label} is not a real UTC timestamp") from error
    _require(parsed.replace(tzinfo=timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ") == value,
        f"{label} is not a canonical UTC timestamp")
    return value


def _absolute_path(value: Any, label: str) -> str:
    # Keep this language identical to the sibling identity contract's private
    # _absolute_posix_path: build roots are joined to its census roots exactly.
    path = _text(value, label)
    _require(all(0x21 <= ord(character) <= 0x7E for character in path),
             f"{label} is not a portable path")
    _require(path.startswith("/") and path != "/" and "//" not in path,
             f"{label} is not a canonical absolute POSIX path")
    parts = path.split("/")[1:]
    _require(all(part not in ("", ".", "..") for part in parts),
             f"{label} is not a canonical absolute POSIX path")
    return path


def _relative_path(value: Any, label: str) -> str:
    path = _text(value, label)
    _require(all(0x21 <= ord(character) <= 0x7E for character in path),
             f"{label} is not a portable path")
    _require(not path.startswith("/") and not path.endswith("/") and
             "//" not in path,
             f"{label} is not a canonical relative POSIX path")
    _require(all(part not in ("", ".", "..") for part in path.split("/")),
             f"{label} is not a canonical relative POSIX path")
    return path


def _require_independent_paths(paths: Sequence[str], label: str) -> None:
    _require(len(paths) == len(set(paths)), f"{label} contain duplicates")
    _require(all(
        not (left + "/").startswith(right + "/") and
        not (right + "/").startswith(left + "/")
        for index, left in enumerate(paths)
        for right in paths[index + 1:]
    ), f"{label} overlap by containment")


def _arguments(value: Any, label: str) -> list[str]:
    _require(type(value) is list and 0 < len(value) <= MAX_ARGUMENTS,
             f"{label} is not a bounded argument list")
    return [_text(argument, f"{label} argument {index}")
            for index, argument in enumerate(value)]


def _text_identity(value: Any, label: str,
                   *, expected_path: str) -> dict[str, Any]:
    record = _dict(value, TEXT_KEYS, label)
    path = _relative_path(record["relative_path"], f"{label} path")
    _require(path == expected_path, f"{label} path changed")
    size, digest = _byte_identity(
        record["size"], record["sha256"], label, minimum_size=0)
    return {"relative_path": path, "size": size, "sha256": digest}


def _file(value: Any, label: str, *,
          expected_path: str | None = None,
          minimum_size: int = 1) -> dict[str, Any]:
    record = _dict(value, FILE_KEYS, label)
    path = _relative_path(record["relative_path"], f"{label} path")
    if expected_path is not None:
        _require(path == expected_path, f"{label} path changed")
    size, digest = _byte_identity(
        record["size"], record["sha256"], label,
        minimum_size=minimum_size)
    return {
        "relative_path": path,
        "size": size,
        "sha256": digest,
    }


def _archive(value: Any, label: str, *, expected_prefix: str,
             expected_path: str) -> dict[str, Any]:
    archive = _dict(value, ARCHIVE_KEYS, label)
    prefix = _text(archive["prefix"], f"{label} prefix")
    _require(prefix == expected_prefix,
             f"{label} prefix changed")
    relative_path = _relative_path(
        archive["relative_path"], f"{label} path")
    _require(relative_path == expected_path, f"{label} path changed")
    size, digest = _byte_identity(
        archive["size"], archive["sha256"], label, minimum_size=1)
    replay_digest = _sha256(
        archive["replay_sha256"], f"{label} replay SHA-256")
    _require(type(archive["replay_identical"]) is bool and
             archive["replay_identical"] is True and
             replay_digest == digest,
             f"{label} was not independently replayed byte-identically")
    return {
        "relative_path": relative_path,
        "prefix": prefix,
        "size": size,
        "sha256": digest,
        "replay_sha256": replay_digest,
        "replay_identical": True,
    }


def _stage_records(value: Any, *, successful: bool) -> list[dict[str, Any]]:
    _require(type(value) is list and value,
             "baseline lane stages are not a non-empty list")
    maximum = len(STAGE_NAMES) if successful else len(FAILURE_STAGES)
    _require(len(value) <= maximum, "baseline lane has too many stages")
    records = []
    for index, raw in enumerate(value):
        stage = _dict(raw, STAGE_KEYS, f"baseline lane stage {index}")
        expected_name = STAGE_NAMES[index]
        _require(stage["name"] == expected_name and
                 type(stage["name"]) is str,
                 "baseline lane stage order changed")
        expected_status = (
            "complete" if successful or index < len(value) - 1 else "failed")
        _require(stage["status"] == expected_status and
                 type(stage["status"]) is str,
                 "baseline lane stage status changed")
        records.append({
            "name": expected_name,
            "status": expected_status,
            "log": _file(
                stage["log"], f"baseline lane {expected_name} log",
                expected_path=f"logs/{index:02d}-{expected_name}.log"),
        })
    _require(not successful or len(records) == len(STAGE_NAMES),
             "promoted baseline lane is incomplete")
    return records


def _lane(value: Any, *, successful: bool) -> dict[str, Any]:
    lane = _dict(value, LANE_KEYS, "baseline lane")
    record_path = _relative_path(
        lane["record_relative_path"], "baseline lane record path")
    expected_path = "baseline-authority.json" if successful else "FAILED.json"
    _require(record_path == expected_path, "baseline lane record path changed")
    _require(lane["attempt_budget"] == 3 and
             type(lane["attempt_budget"]) is int,
             "baseline attempt budget changed")
    attempt = _integer(
        lane["attempt"], 1, 3, "baseline lane attempt")
    return {
        "root": _absolute_path(lane["root"], "baseline lane root"),
        "attempt": attempt,
        "attempt_budget": 3,
        "record_relative_path": record_path,
        "seal_protocol": _fixed_text(
            lane["seal_protocol"], SEAL_PROTOCOL,
            "baseline lane seal protocol"),
        "stages": _stage_records(lane["stages"], successful=successful),
    }


def _fixed_text(value: Any, expected: str, label: str) -> str:
    _require(type(value) is str and value == expected, f"{label} changed")
    return expected


def _host(value: Any) -> dict[str, Any]:
    host = _dict(value, HOST_KEYS, "baseline host")
    cpus_value = host["online_cpus"]
    _require(type(cpus_value) is list and
             1 <= len(cpus_value) <= MAX_CPU_COUNT,
             "baseline host online CPU list is invalid")
    cpus = [_integer(cpu, 0, (1 << 31) - 1,
                     f"baseline host online CPU {index}")
            for index, cpu in enumerate(cpus_value)]
    _require(cpus == sorted(set(cpus)),
             "baseline host online CPUs are not sorted and unique")
    return {
        "hostname": _text(host["hostname"], "baseline host name"),
        "kernel": _text(host["kernel"], "baseline host kernel"),
        "architecture": _fixed_text(
            host["architecture"], "x86_64", "baseline host architecture"),
        "cpu_model": _text(host["cpu_model"], "baseline host CPU model"),
        "online_cpus": cpus,
        "sc_clk_tck": _integer(
            host["sc_clk_tck"], 1, 1_000_000,
            "baseline host SC_CLK_TCK"),
    }


def _source(value: Any) -> dict[str, Any]:
    source = _dict(value, SOURCE_KEYS, "baseline source")
    baseline = _dict(
        source["baseline"], BASELINE_SOURCE_KEYS, "Leopard1 source")
    _require(baseline["commit"] == BASELINE_COMMIT and
             type(baseline["commit"]) is str,
             "Leopard1 source commit changed")
    _require(baseline["tree"] == BASELINE_TREE and
             type(baseline["tree"]) is str,
             "Leopard1 source tree changed")
    submodule = _dict(
        baseline["submodule"], SUBMODULE_KEYS, "Leopard1 submodule")
    submodule_path = _fixed_text(
        submodule["path"], "sse2neon", "Leopard1 sse2neon path")
    submodule_commit = _fixed_text(
        submodule["commit"], BASELINE_SSE2NEON_COMMIT,
        "Leopard1 sse2neon commit")
    adapter = _dict(
        source["adapter_repository"], ADAPTER_SOURCE_KEYS,
        "adapter repository source")
    _require(type(adapter["clean"]) is bool and adapter["clean"] is True,
             "adapter repository was not clean")
    result = {
        "baseline": {
            "commit": BASELINE_COMMIT,
            "tree": BASELINE_TREE,
            "submodule": {
                "path": submodule_path,
                "commit": submodule_commit,
            },
            "git_capture": _file(
                baseline["git_capture"], "Leopard1 git capture",
                expected_path="source/leopard-main-git-capture.json"),
            "archive": _archive(
                baseline["archive"], "Leopard1 source archive",
                expected_prefix="leopard-main-source/",
                expected_path="source/leopard-main-source.tar"),
        },
        "adapter_repository": {
            "commit": _git_oid(
                adapter["commit"], "adapter repository commit"),
            "tree": _git_oid(adapter["tree"], "adapter repository tree"),
            "clean": True,
            "git_capture": _file(
                adapter["git_capture"], "adapter repository git capture",
                expected_path="source/adapter-git-capture.json"),
            "archive": _archive(
                adapter["archive"], "adapter source archive",
                expected_prefix="leopard2-adapter-source/",
                expected_path="source/leopard2-adapter-source.tar"),
        },
    }
    return result


def _adapter(value: Any) -> dict[str, Any]:
    adapter = _dict(value, ADAPTER_KEYS, "exact-main adapter")
    _require(adapter["minimum_harness_commit"] == MINIMUM_HARNESS_COMMIT and
             type(adapter["minimum_harness_commit"]) is str,
             "exact-main minimum harness commit changed")
    values = adapter["files"]
    _require(type(values) is list and len(values) == len(ADAPTER_PATHS),
             "exact-main adapter file inventory changed")
    files = []
    for index, raw in enumerate(values):
        record = _dict(raw, ADAPTER_FILE_KEYS,
                       f"exact-main adapter file {index}")
        path = _relative_path(
            record["path"], f"exact-main adapter file {index} path")
        _require(path == ADAPTER_PATHS[index],
                 "exact-main adapter file order changed")
        size, digest = _byte_identity(
            record["size"], record["sha256"],
            f"exact-main adapter file {index}", minimum_size=1)
        files.append({
            "path": path,
            "git_blob_sha1": _git_oid(
                record["git_blob_sha1"],
                f"exact-main adapter file {index} Git blob"),
            "size": size,
            "sha256": digest,
        })
    blob_identities: dict[str, tuple[int, str]] = {}
    identity_blobs: dict[tuple[int, str], str] = {}
    for record in files:
        blob_oid = record["git_blob_sha1"]
        blob_identity = (record["size"], record["sha256"])
        if blob_oid in blob_identities:
            _require(blob_identities[blob_oid] == blob_identity,
                     "exact-main adapter files reuse one Git blob object ID "
                     "at two identities")
        else:
            blob_identities[blob_oid] = blob_identity
        if blob_identity in identity_blobs:
            _require(identity_blobs[blob_identity] == blob_oid,
                     "exact-main adapter files reuse one byte identity with "
                     "two Git blob object IDs")
        else:
            identity_blobs[blob_identity] = blob_oid
    digest = _sha256(adapter["files_sha256"],
                     "exact-main adapter inventory SHA-256")
    _require(digest == canonical_sha256(files),
             "exact-main adapter inventory SHA-256 is inconsistent")
    return {
        "minimum_harness_commit": MINIMUM_HARNESS_COMMIT,
        "files": files,
        "files_sha256": digest,
    }


def _tool(value: Any, label: str, expected_role: str) -> dict[str, Any]:
    tool = _dict(value, TOOL_KEYS, label)
    _require(tool["role"] == expected_role and
             type(tool["role"]) is str,
             f"{label} role changed")
    mode = _integer(tool["mode"], 0, 0o7777, f"{label} mode")
    _require(mode & 0o111 != 0, f"{label} is not executable")
    size, digest = _byte_identity(
        tool["size"], tool["sha256"], label, minimum_size=1)
    return {
        "role": expected_role,
        "path": _absolute_path(tool["path"], f"{label} path"),
        "resolved_path": _absolute_path(
            tool["resolved_path"], f"{label} resolved path"),
        "size": size,
        "mode": mode,
        "sha256": digest,
    }


def _toolchain(value: Any) -> dict[str, Any]:
    toolchain = _dict(value, TOOLCHAIN_KEYS, "baseline toolchain")
    tools_value = toolchain["tools"]
    subtools_value = toolchain["subtools"]
    versions_value = toolchain["versions"]
    _require(type(tools_value) is list and len(tools_value) == len(TOOL_ROLES),
             "baseline tool inventory changed")
    _require(type(subtools_value) is list and
             len(subtools_value) == len(SUBTOOL_ROLES),
             "baseline compiler-subtool inventory changed")
    tools = [_tool(raw, f"baseline tool {role}", role)
             for role, raw in zip(TOOL_ROLES, tools_value)]
    subtools = [_tool(raw, f"baseline compiler subtool {role}", role)
                for role, raw in zip(SUBTOOL_ROLES, subtools_value)]
    resolved = {
        record["role"]: record["resolved_path"]
        for record in tools + subtools
    }
    _require(type(versions_value) is list and
             len(versions_value) == len(VERSION_ROLES),
             "baseline tool version inventory changed")
    versions = []
    for index, (role, raw) in enumerate(zip(VERSION_ROLES, versions_value)):
        version = _dict(raw, VERSION_KEYS, f"baseline version {index}")
        _require(version["role"] == role and type(version["role"]) is str,
                 "baseline tool version order changed")
        argv = _arguments(version["argv"], f"baseline {role} version argv")
        _require(argv[0] == resolved[role],
                 f"baseline {role} version did not use the retained tool")
        _require(version["exit_status"] == 0 and
                 type(version["exit_status"]) is int,
                 f"baseline {role} version command failed")
        versions.append({
            "role": role,
            "argv": argv,
            "stdout": _text_identity(
                version["stdout"], f"baseline {role} version stdout",
                expected_path=f"toolchain/versions/{role}.stdout"),
            "stderr": _text_identity(
                version["stderr"], f"baseline {role} version stderr",
                expected_path=f"toolchain/versions/{role}.stderr"),
            "exit_status": 0,
        })
    return {"tools": tools, "subtools": subtools, "versions": versions}


def exact_main_build_profile() -> dict[str, Any]:
    """Return the detached immutable-by-validation v1 build profile."""
    return {
        "schema": BUILD_PROFILE_SCHEMA,
        "name": "canonical-pure-avx2-exact-main/v1",
        "generator": "Unix Makefiles",
        "build_type": "Release",
        "language": "gnu++11",
        "release_flags": ["-g", "-O0", "-O3"],
        "isa_flags": [
            "-march=x86-64", "-mtune=generic", "-mavx2", "-mno-avx512f",
        ],
        "warning_flags": ["-Wall", "-Wextra"],
        "openmp": True,
        "pure_avx2": True,
        "adapter_definitions": [
            f'LEOPARD_MAIN_SOURCE_COMMIT="{BASELINE_COMMIT}"',
            "LEO_MAIN_PURE_AVX2_PROFILE=1",
        ],
        "environment": [
            {"name": "GIT_ATTR_NOSYSTEM", "value": "1"},
            {"name": "GIT_CONFIG_GLOBAL", "value": "/dev/null"},
            {"name": "GIT_CONFIG_NOSYSTEM", "value": "1"},
            {"name": "GIT_CONFIG_SYSTEM", "value": "/dev/null"},
            {"name": "GIT_NO_REPLACE_OBJECTS", "value": "1"},
            {"name": "GIT_OPTIONAL_LOCKS", "value": "0"},
            {"name": "LANG", "value": "C"},
            {"name": "LC_ALL", "value": "C"},
            {"name": "PATH", "value": "/usr/bin:/bin"},
            {"name": "TZ", "value": "UTC"},
        ],
        "cache_requirements": [
            {"name": "CMAKE_BUILD_TYPE", "type": "STRING",
             "value": "Release"},
            {"name": "CMAKE_CXX_FLAGS_RELEASE", "type": "STRING",
             "value": "-g -O0 -O3"},
            {"name": "LEOPARD_MAIN_SOURCE_DIR", "type": "PATH",
             "value": "${BASELINE_SOURCE_ROOT}"},
            {"name": "LEO_MAIN_PURE_AVX2", "type": "BOOL",
             "value": "ON"},
        ],
        "build_job_count": 1,
    }


def exact_main_build_cache_requirements(
    roots: Mapping[str, Any],
) -> tuple[dict[str, str], ...]:
    """Expand the fixed cache requirements for one validated build root set."""
    roots_value = _dict(roots, ROOTS_KEYS, "exact-main build roots")
    canonical_roots = {
        key: _absolute_path(roots_value[key], f"exact-main {key}")
        for key in ("adapter_source_root", "baseline_source_root", "build_root")
    }
    _require_independent_paths(
        list(canonical_roots.values()), "exact-main build roots")
    requirements = copy.deepcopy(
        exact_main_build_profile()["cache_requirements"])
    for requirement in requirements:
        if requirement["value"] == "${BASELINE_SOURCE_ROOT}":
            requirement["value"] = canonical_roots["baseline_source_root"]
    return tuple(requirements)


def _validated_build_roots(roots: Mapping[str, Any]) -> dict[str, str]:
    roots_value = _dict(roots, ROOTS_KEYS, "exact-main build roots")
    canonical = {
        key: _absolute_path(roots_value[key], f"exact-main {key}")
        for key in ("adapter_source_root", "baseline_source_root", "build_root")
    }
    _require_independent_paths(list(canonical.values()),
                               "exact-main build roots")
    return canonical


def _validated_tool_path(value: Any, label: str) -> str:
    return _absolute_path(value, label)


def exact_main_configure_argv(
    *, cmake: str, compiler: str, roots: Mapping[str, Any],
) -> list[str]:
    """Return the frozen CMake configure argv for one exact-main build."""
    canonical = _validated_build_roots(roots)
    cmake_path = _validated_tool_path(cmake, "exact-main CMake path")
    compiler_path = _validated_tool_path(
        compiler, "exact-main compiler path")
    return [
        cmake_path,
        "-S", canonical["adapter_source_root"] +
        "/experiments/leopard2/main_compare",
        "-B", canonical["build_root"],
        "-G", "Unix Makefiles",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DLEO_MAIN_PURE_AVX2=ON",
        "-DLEOPARD_MAIN_SOURCE_DIR=" + canonical["baseline_source_root"],
        "-DCMAKE_CXX_COMPILER=" + compiler_path,
    ]


def exact_main_build_argv(
    *, cmake: str, roots: Mapping[str, Any],
) -> list[str]:
    """Return the frozen serial CMake build argv."""
    canonical = _validated_build_roots(roots)
    cmake_path = _validated_tool_path(cmake, "exact-main CMake path")
    return [
        cmake_path, "--build", canonical["build_root"],
        "--target", "leopard_main_benchmark", "--", "-j1",
    ]


def exact_main_benchmark_argv(*, executable_path: str) -> list[str]:
    """Return the fixed correctness-only benchmark attestation argv."""
    executable = _absolute_path(
        executable_path, "exact-main benchmark executable path")
    return [
        executable,
        "--k", "8", "--r", "4", "--bytes", "64", "--loss", "1",
        "--batch", "2", "--reuse", "1", "--iterations", "2",
        "--warmup", "1", "--threads", "1", "--seed", "7",
        "--json", "-",
    ]


def exact_main_ctest_argv(
    *, ctest: str, build_root: str,
) -> list[str]:
    """Return the frozen single-smoke-test CTest argv."""
    ctest_path = _validated_tool_path(ctest, "exact-main CTest path")
    root = _absolute_path(build_root, "exact-main CTest build root")
    return [
        ctest_path, "--test-dir", root,
        "--output-on-failure", "-R", "^leopard_main_benchmark_smoke$",
    ]


def exact_main_compile_sources(
    roots: Mapping[str, Any],
) -> tuple[str, ...]:
    """Return the exact source inventory in deterministic compile order."""
    canonical = _validated_build_roots(roots)
    return (
        canonical["adapter_source_root"] +
        "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp",
        *(canonical["baseline_source_root"] + "/" + source
          for source in BASELINE_CPP_SOURCES),
    )


def exact_main_compile_command_argv(
    *,
    roots: Mapping[str, Any],
    compiler: str,
    source: str,
) -> tuple[list[str], str]:
    """Return the exact compile argv and CMake-relative output for one TU."""
    canonical = _validated_build_roots(roots)
    compiler_path = _validated_tool_path(
        compiler, "exact-main compiler path")
    sources = exact_main_compile_sources(canonical)
    _require(type(source) is str and source in sources,
             "exact-main compile source is outside the frozen inventory")
    adapter_source = sources[0]
    is_adapter = source == adapter_source
    profile = exact_main_build_profile()
    output = (
        "CMakeFiles/leopard_main_benchmark.dir/"
        "legacy_main_benchmark.cpp.o" if is_adapter else
        "CMakeFiles/leopard_main_exact.dir" + source + ".o")
    argv = [
        compiler_path,
        *(["-D" + definition
           for definition in profile["adapter_definitions"]]
          if is_adapter else []),
        "-I" + canonical["baseline_source_root"],
        *profile["release_flags"],
        "-std=" + profile["language"],
        *profile["isa_flags"],
        *profile["warning_flags"],
        *(["-fopenmp"] if profile["openmp"] is True else []),
        "-o", output,
        "-c", source,
    ]
    return argv, output


def parse_cmake_cache(content: bytes) -> dict[str, tuple[str, str]]:
    """Parse the strict LF-framed CMake cache grammar used by v1."""
    _require(type(content) is bytes and len(content) <= MAX_FILE_BYTES,
             "CMakeCache.txt is not bounded bytes")
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ExactMainBaselineRecordError(
            "CMakeCache.txt is not strict UTF-8") from error
    _require(text.endswith("\n") and "\r" not in text and "\0" not in text,
             "CMakeCache.txt is not canonical LF text")
    result: dict[str, tuple[str, str]] = {}
    for line in text[:-1].split("\n"):
        if not line or line.startswith("//") or line.startswith("#"):
            continue
        head, separator, value = line.partition("=")
        name, typed, entry_type = head.partition(":")
        _require(bool(separator) and bool(typed) and name and entry_type and
                 name not in result,
                 "CMakeCache.txt contains a malformed or duplicate entry")
        result[name] = (entry_type, value)
    return result


def validate_cmake_cache(
    content: bytes, roots: Mapping[str, Any],
) -> dict[str, tuple[str, str]]:
    """Validate and return the exact required cache projection."""
    cache = parse_cmake_cache(content)
    _require(cache.get("CMAKE_CXX_FLAGS") == ("STRING", ""),
             "CMake cache injects global C++ flags")
    for requirement in exact_main_build_cache_requirements(roots):
        _require(cache.get(requirement["name"]) ==
                 (requirement["type"], requirement["value"]),
                 f"CMake cache changed {requirement['name']}")
    return copy.deepcopy(cache)


def _compile_argv(entry: Mapping[str, Any], label: str) -> list[str]:
    if set(entry) == {"directory", "file", "output", "arguments"}:
        _require(type(entry["arguments"]) is list and
                 all(type(argument) is str for argument in entry["arguments"]),
                 f"{label} arguments are invalid")
        return list(entry["arguments"])
    _require(set(entry) == {"directory", "file", "output", "command"} and
             type(entry["command"]) is str,
             f"{label} has an unexpected shape")
    command = entry["command"]
    arguments: list[str] = []
    current: list[str] = []
    state = "plain"
    quoted = False
    index = 0
    while index < len(command):
        character = command[index]
        if state == "plain":
            if character in " \t\n":
                if current or quoted:
                    arguments.append("".join(current))
                    current = []
                    quoted = False
            elif character == "'":
                state = "single"
                quoted = True
            elif character == '"':
                state = "double"
                quoted = True
            elif character == "\\":
                index += 1
                _require(index < len(command),
                         f"{label} command cannot be parsed")
                current.append(command[index])
            else:
                current.append(character)
        elif state == "single":
            if character == "'":
                state = "plain"
            else:
                current.append(character)
        else:
            if character == '"':
                state = "plain"
            elif character == "\\":
                index += 1
                _require(index < len(command),
                         f"{label} command cannot be parsed")
                escaped = command[index]
                if escaped in ('"', "\\", "$", "`"):
                    current.append(escaped)
                elif escaped == "\n":
                    pass
                else:
                    current.extend(("\\", escaped))
            else:
                current.append(character)
        index += 1
    _require(state == "plain", f"{label} command cannot be parsed")
    if current or quoted:
        arguments.append("".join(current))
    _require(bool(arguments), f"{label} command cannot be parsed")
    return arguments


def validate_compile_commands(
    value: Any,
    *,
    roots: Mapping[str, Any],
    compiler: str,
    profile: Mapping[str, Any],
) -> tuple[str, ...]:
    """Validate the exact five-TU compile closure for one build."""
    canonical_roots = _validated_build_roots(roots)
    compiler_path = _validated_tool_path(
        compiler, "exact-main compiler path")
    _require(type(profile) is dict and
             exact_json_equal(profile, exact_main_build_profile()),
             "exact-main compile profile changed")
    _require(type(value) is list and 1 <= len(value) <= 4096,
             "exact-main compile command inventory is invalid")
    expected_sources = set(exact_main_compile_sources(canonical_roots))
    observed_sources: set[str] = set()
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(value):
        _require(type(raw) is dict,
                 f"exact-main compile command {index} is invalid")
        argv = _compile_argv(raw, f"exact-main compile command {index}")
        source = raw.get("file")
        _require(type(source) is str and source in expected_sources and
                 source not in observed_sources,
                 f"exact-main compile command {index} differs from the profile")
        expected_argv, expected_output = exact_main_compile_command_argv(
            roots=canonical_roots, compiler=compiler_path, source=source)
        _require(raw.get("directory") == canonical_roots["build_root"] and
                 raw.get("output") == expected_output and
                 argv == expected_argv,
                 f"exact-main compile command {index} differs from the profile")
        observed_sources.add(source)
        normalized.append(copy.deepcopy(raw))
    _require(observed_sources == expected_sources,
             "exact-main compile closure changed its exact source set")
    return tuple(raw["file"] for raw in normalized)


def validate_build_closure(
    value: Any, *, role: str, build: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate a canonical build-tree census and its two artifact joins."""
    _require(role in BUILD_ROLES and type(role) is str,
             "build closure role is invalid")
    _require(type(build) is dict and "roots" in build and
             "closure" in build and "executable" in build and
             "archive" in build,
             f"{role} build is missing closure claims")
    closure = _dict(value, BUILD_CLOSURE_KEYS, f"{role} build closure")
    build_root = _validated_build_roots(build["roots"])["build_root"]
    closure_claim = _dict(
        build["closure"], CLOSURE_KEYS, f"{role} build closure claim")
    _require(closure["schema"] == BUILD_CLOSURE_SCHEMA and
             closure["role"] == role and
             closure["build_root"] == build_root and
             type(closure["file_count"]) is int and
             1 <= closure["file_count"] <= MAX_CLOSURE_FILES and
             closure["file_count"] == closure_claim["file_count"] and
             type(closure["files"]) is list and
             len(closure["files"]) == closure["file_count"],
             f"{role} build closure header changed")
    files: list[dict[str, Any]] = []
    for index, raw in enumerate(closure["files"]):
        item = _dict(raw, CLOSURE_FILE_KEYS,
                     f"{role} closure file {index}")
        relative_path = _relative_path(
            item["relative_path"], f"{role} closure file {index} path")
        size, digest = _byte_identity(
            item["size"], item["sha256"],
            f"{role} closure file {index}", minimum_size=0)
        files.append({"relative_path": relative_path,
                      "size": size, "sha256": digest})
    paths = [item["relative_path"] for item in files]
    _require(paths == sorted(set(paths)),
             f"{role} build closure paths are not sorted and unique")
    by_path = {item["relative_path"]: item for item in files}
    for artifact_name in ("executable", "archive"):
        artifact = build[artifact_name]
        path = artifact["build_relative_path"]
        _require(path in by_path and
                 by_path[path]["size"] == artifact["size"] and
                 by_path[path]["sha256"] == artifact["sha256"],
                 f"{role} build closure does not bind its {artifact_name}")
    return {
        "schema": BUILD_CLOSURE_SCHEMA,
        "role": role,
        "build_root": build_root,
        "files": files,
        "file_count": len(files),
    }


def _argv_options(argv: Sequence[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    index = 1
    while index < len(argv):
        option = argv[index]
        _require(option.startswith("--") and index + 1 < len(argv),
                 "attestation argv is malformed")
        _require(option not in result, "attestation argv repeats an option")
        result[option] = argv[index + 1]
        index += 2
    return result


def validate_attestation_stdout(
    value: Any, *, argv: Sequence[str], reported_schema: str,
) -> None:
    """Validate correctness/profile claims in benchmark JSON, not timing."""
    arguments = _arguments(argv, "attestation argv")
    schema = _fixed_text(
        reported_schema, BENCHMARK_SCHEMA, "benchmark reported schema")
    _require(type(value) is dict and
             value.get("schema") == schema and
             value.get("schema") == BENCHMARK_SCHEMA and
             type(value.get("build")) is dict and
             value["build"].get("main_source_commit") == BASELINE_COMMIT and
             value["build"].get("pure_avx2") is True and
             type(value.get("correctness")) is dict and
             value["correctness"].get("round_trip") is True and
             type(value.get("parameters")) is dict,
             "retained benchmark attestation claims changed")
    options = _argv_options(arguments)
    try:
        expected = {
            "K": int(options["--k"]),
            "R": int(options["--r"]),
            "shard_bytes": int(options["--bytes"]),
            "loss_count": int(options["--loss"]),
            "batch": int(options["--batch"]),
            "reuse": int(options["--reuse"]),
            "iterations": int(options["--iterations"]),
            "warmup": int(options["--warmup"]),
            "thread_count": int(options["--threads"]),
            "seed": int(options["--seed"]),
        }
    except (KeyError, ValueError) as error:
        raise ExactMainBaselineRecordError(
            "attestation argv is malformed") from error
    _require(all(type(value["parameters"].get(key)) is int and
                 value["parameters"].get(key) == expected_value
                 for key, expected_value in expected.items()),
             "retained benchmark attestation parameters changed")
    return None


def validate_exact_main_build(
    value: Any, *, role: str, tools: Mapping[str, Any],
) -> dict[str, Any]:
    """Public pure validator for one exact-main build record."""
    _require(type(tools) is dict and
             set(("cmake", "compiler")).issubset(tools),
             "exact-main build tool map is incomplete")
    return _build(value, role, {
        key: _validated_tool_path(raw, f"exact-main {key} path")
        for key, raw in tools.items()
    })


def _profile(value: Any) -> dict[str, Any]:
    _dict(value, PROFILE_KEYS, "exact-main build profile")
    expected = exact_main_build_profile()
    _require(exact_json_equal(value, expected),
             "exact-main build profile changed")
    return copy.deepcopy(expected)


def _build_artifact(value: Any, label: str, *, expected_name: str,
                    expected_build_path: str, expected_retained_path: str
                    ) -> dict[str, Any]:
    artifact = _dict(value, BUILD_ARTIFACT_KEYS, label)
    _require(artifact["name"] == expected_name and
             type(artifact["name"]) is str,
             f"{label} name changed")
    build_path = _relative_path(
        artifact["build_relative_path"], f"{label} build path")
    retained_path = _relative_path(
        artifact["retained_relative_path"], f"{label} retained path")
    _require(build_path == expected_build_path,
             f"{label} build path changed")
    _require(retained_path == expected_retained_path,
             f"{label} retained path changed")
    size, digest = _byte_identity(
        artifact["size"], artifact["sha256"], label, minimum_size=1)
    return {
        "name": expected_name,
        "build_relative_path": build_path,
        "retained_relative_path": retained_path,
        "size": size,
        "sha256": digest,
    }


def _build(value: Any, role: str,
           tools: Mapping[str, str]) -> dict[str, Any]:
    build = _dict(value, BUILD_KEYS, f"{role} build")
    _require(build["role"] == role and type(build["role"]) is str,
             f"{role} build role changed")
    roots_value = _dict(build["roots"], ROOTS_KEYS, f"{role} build roots")
    roots = {
        key: _absolute_path(roots_value[key], f"{role} {key}")
        for key in ("adapter_source_root", "baseline_source_root", "build_root")
    }
    _require_independent_paths(list(roots.values()), f"{role} build roots")
    configure_argv = _arguments(
        build["configure_argv"], f"{role} configure argv")
    expected_configure = exact_main_configure_argv(
        roots=roots, cmake=tools["cmake"], compiler=tools["compiler"])
    _require(configure_argv == expected_configure,
             f"{role} configure command changed")
    build_argv = _arguments(build["build_argv"], f"{role} build argv")
    expected_build = exact_main_build_argv(
        roots=roots, cmake=tools["cmake"])
    _require(build_argv == expected_build, f"{role} build command changed")
    role_path = role.replace("_", "-")
    executable = _build_artifact(
        build["executable"], f"{role} executable",
        expected_name="leopard_main_benchmark",
        expected_build_path="leopard_main_benchmark",
        expected_retained_path=(
            f"artifacts/{role_path}/leopard_main_benchmark"))
    archive = _build_artifact(
        build["archive"], f"{role} archive",
        expected_name="libleopard_main_exact.a",
        expected_build_path="libleopard_main_exact.a",
        expected_retained_path=(
            f"artifacts/{role_path}/libleopard_main_exact.a"))
    closure_value = _dict(
        build["closure"], CLOSURE_KEYS, f"{role} build closure")
    closure_size, closure_digest = _byte_identity(
        closure_value["size"], closure_value["sha256"],
        f"{role} build closure", minimum_size=1)
    closure = {
        "relative_path": _relative_path(
            closure_value["relative_path"], f"{role} closure path"),
        "size": closure_size,
        "sha256": closure_digest,
        "file_count": _integer(
            closure_value["file_count"], 1, MAX_CLOSURE_FILES,
            f"{role} closure file count"),
    }
    _require(closure["relative_path"] ==
             f"builds/{role_path}/build-closure.json",
             f"{role} closure path changed")
    configure_log = _file(
        build["configure_log"], f"{role} configure log")
    build_log = _file(build["build_log"], f"{role} build log")
    cmake_cache = _file(build["cmake_cache"], f"{role} CMake cache")
    compile_commands = _file(
        build["compile_commands"], f"{role} compile commands")
    expected_paths = {
        "configure_log": f"builds/{role_path}/configure.log",
        "build_log": f"builds/{role_path}/build.log",
        "cmake_cache": f"builds/{role_path}/CMakeCache.txt",
        "compile_commands": f"builds/{role_path}/compile_commands.json",
    }
    for name, artifact in (
            ("configure_log", configure_log),
            ("build_log", build_log),
            ("cmake_cache", cmake_cache),
            ("compile_commands", compile_commands)):
        _require(artifact["relative_path"] == expected_paths[name],
                 f"{role} {name} path changed")
    result = {
        "role": role,
        "roots": roots,
        "configure_argv": configure_argv,
        "build_argv": build_argv,
        "configure_log": configure_log,
        "build_log": build_log,
        "cmake_cache": cmake_cache,
        "compile_commands": compile_commands,
        "executable": executable,
        "archive": archive,
        "closure": closure,
    }
    retained_paths = [
        result[key]["relative_path"]
        for key in ("configure_log", "build_log", "cmake_cache",
                    "compile_commands")
    ] + [
        executable["retained_relative_path"],
        archive["retained_relative_path"],
        closure["relative_path"],
    ]
    _require(len(set(retained_paths)) == len(retained_paths),
             f"{role} retained paths are not unique")
    return result


def _builds(value: Any, toolchain: Mapping[str, Any], lane_root: str
            ) -> dict[str, dict[str, Any]]:
    records = _dict(value, BUILDS_KEYS, "exact-main builds")
    tool_paths = {
        record["role"]: record["resolved_path"]
        for record in toolchain["tools"]
    }
    builds = {
        role: _build(records[role], role, tool_paths) for role in BUILD_ROLES
    }
    canonical_roots = builds["canonical_first"]["roots"]
    _require(exact_json_equal(
        canonical_roots, builds["canonical_second"]["roots"]),
        "same-path canonical builds did not use the same roots")
    variant_roots = builds["path_variant"]["roots"]
    _require_independent_paths(
        [lane_root] + list(canonical_roots.values()) +
        list(variant_roots.values()),
        "lane and canonical/path-variant build roots",
    )
    return builds


def _require_unique_retained_paths(paths: Sequence[str], label: str) -> None:
    _require_independent_paths(paths, f"{label} retained paths")


def _require_digest_size_function(
    byte_identities: Sequence[tuple[str, int]],
    label: str,
) -> None:
    observed: dict[str, int] = {}
    for digest, size in byte_identities:
        if digest in observed:
            _require(observed[digest] == size,
                     f"{label} reuse one SHA-256 at two byte lengths")
        else:
            observed[digest] = size


def _require_digest_claim_function(
    claims: Sequence[tuple[str, str, Any]],
    label: str,
) -> None:
    observed: dict[tuple[str, str], Any] = {}
    for digest, claim_name, claim_value in claims:
        key = (digest, claim_name)
        if key in observed:
            _require(exact_json_equal(observed[key], claim_value),
                     f"{label} reuse one SHA-256 with two "
                     "content-derived identities")
        else:
            observed[key] = copy.deepcopy(claim_value)


def _require_host_file_function(
    host_files: Sequence[tuple[str, int, str]],
    label: str,
) -> None:
    observed: dict[str, tuple[int, str]] = {}
    for path, size, digest in host_files:
        identity = (size, digest)
        if path in observed:
            _require(observed[path] == identity,
                     f"{label} reuse one absolute path with two identities")
        else:
            observed[path] = identity


def _require_tool_path_function(
    tools: Sequence[Mapping[str, Any]],
    label: str,
) -> None:
    resolutions: dict[str, str] = {}
    modes: dict[str, int] = {}
    for tool in tools:
        path = tool["path"]
        resolved_path = tool["resolved_path"]
        if path in resolutions:
            _require(resolutions[path] == resolved_path,
                     f"{label} resolve one tool path two ways")
        else:
            resolutions[path] = resolved_path
        if resolved_path in modes:
            _require(modes[resolved_path] == tool["mode"],
                     f"{label} give one resolved tool path two modes")
        else:
            modes[resolved_path] = tool["mode"]


def _require_git_object_consistency(
    source: Mapping[str, Any],
    adapter: Mapping[str, Any],
    label: str,
) -> None:
    object_types: dict[str, str] = {}

    def add_object(object_id: str, object_type: str) -> None:
        if object_id in object_types:
            _require(object_types[object_id] == object_type,
                     f"{label} reuse one Git object ID as two object types")
        else:
            object_types[object_id] = object_type

    commit_ids = (
        source["baseline"]["commit"],
        source["baseline"]["submodule"]["commit"],
        source["adapter_repository"]["commit"],
        adapter["minimum_harness_commit"],
    )
    tree_ids = (
        source["baseline"]["tree"],
        source["adapter_repository"]["tree"],
    )
    for object_id in commit_ids:
        add_object(object_id, "commit")
    for object_id in tree_ids:
        add_object(object_id, "tree")
    for adapter_file in adapter["files"]:
        add_object(adapter_file["git_blob_sha1"], "blob")

    commit_trees: dict[str, str] = {}
    for commit_id, tree_id in (
        (source["baseline"]["commit"], source["baseline"]["tree"]),
        (source["adapter_repository"]["commit"],
         source["adapter_repository"]["tree"]),
    ):
        if commit_id in commit_trees:
            _require(commit_trees[commit_id] == tree_id,
                     f"{label} reuse one commit object ID with two trees")
        else:
            commit_trees[commit_id] = tree_id


def _retained_inventory_entry(
    record: Mapping[str, Any],
    *,
    path_key: str = "relative_path",
) -> dict[str, Any]:
    return {
        "relative_path": record[path_key],
        "size": record["size"],
        "sha256": record["sha256"],
    }


def _authority_retained_inventory_from_validated(
    lane: Mapping[str, Any],
    source: Mapping[str, Any],
    toolchain: Mapping[str, Any],
    builds: Mapping[str, Mapping[str, Any]],
    runtime: Mapping[str, Any],
    attestation: Mapping[str, Any],
) -> tuple[dict[str, Any], ...]:
    inventory: list[dict[str, Any]] = [{
        "relative_path": lane["record_relative_path"],
        "size": None,
        "sha256": None,
    }]
    inventory.extend(
        _retained_inventory_entry(stage["log"])
        for stage in lane["stages"])
    for source_role in ("baseline", "adapter_repository"):
        inventory.extend((
            _retained_inventory_entry(source[source_role]["git_capture"]),
            _retained_inventory_entry(source[source_role]["archive"]),
        ))
    for version in toolchain["versions"]:
        inventory.extend((
            _retained_inventory_entry(version["stdout"]),
            _retained_inventory_entry(version["stderr"]),
        ))
    for role in BUILD_ROLES:
        build = builds[role]
        inventory.extend(
            _retained_inventory_entry(build[key])
            for key in ("configure_log", "build_log", "cmake_cache",
                        "compile_commands"))
        inventory.extend((
            _retained_inventory_entry(
                build["executable"], path_key="retained_relative_path"),
            _retained_inventory_entry(
                build["archive"], path_key="retained_relative_path"),
            _retained_inventory_entry(build["closure"]),
        ))
    inventory.extend(
        _retained_inventory_entry(record["canonical_ldd_output"])
        for record in runtime["records"])
    inventory.append(_retained_inventory_entry(attestation["test_controller"]))
    for record in attestation["records"]:
        inventory.extend((
            _retained_inventory_entry(record["stdout"]),
            _retained_inventory_entry(record["stderr"]),
            _retained_inventory_entry(record["ctest"]["stdout"]),
            _retained_inventory_entry(record["ctest"]["stderr"]),
        ))
    inventory.sort(key=lambda item: item["relative_path"])
    return tuple(copy.deepcopy(inventory))


def _authority_evidence_inventory(
    lane: Mapping[str, Any],
    source: Mapping[str, Any],
    adapter: Mapping[str, Any],
    toolchain: Mapping[str, Any],
    builds: Mapping[str, Mapping[str, Any]],
    runtime: Mapping[str, Any],
    attestation: Mapping[str, Any],
    identity: Mapping[str, Any],
) -> None:
    retained_inventory = _authority_retained_inventory_from_validated(
        lane, source, toolchain, builds, runtime, attestation)
    paths = [item["relative_path"] for item in retained_inventory]
    byte_identities: list[tuple[str, int]] = [
        (item["sha256"], item["size"])
        for item in retained_inventory
        if item["sha256"] is not None and item["size"] is not None
    ]
    byte_claims: list[tuple[str, str, Any]] = []
    host_files: list[tuple[str, int, str]] = []

    def add_bytes(record: Mapping[str, Any],
                  digest_key: str = "sha256") -> None:
        byte_identities.append((record[digest_key], record["size"]))

    for source_role in ("baseline", "adapter_repository"):
        add_bytes(source[source_role]["archive"], "replay_sha256")
        byte_claims.append((
            source[source_role]["archive"]["sha256"], "archive_prefix",
            source[source_role]["archive"]["prefix"],
        ))
    for adapter_file in adapter["files"]:
        add_bytes(adapter_file)
    for tool in toolchain["tools"] + toolchain["subtools"]:
        add_bytes(tool)
        host_files.append((
            tool["resolved_path"], tool["size"], tool["sha256"]))
    for role in BUILD_ROLES:
        build = builds[role]
        byte_claims.append((
            build["closure"]["sha256"], "closure_file_count",
            build["closure"]["file_count"],
        ))
    for runtime_record in runtime["records"]:
        for dependency in runtime_record["dependencies"]:
            if dependency["kind"] == "file":
                add_bytes(dependency)
                host_files.append((
                    dependency["path"], dependency["size"],
                    dependency["sha256"]))
    add_bytes(attestation["test_controller"])
    for role in BUILD_ROLES:
        normalized = identity[role]
        add_bytes(normalized["artifact"])
        for section in normalized["sections"]:
            if section["content_sha256"] is not None:
                byte_identities.append(
                    (section["content_sha256"], section["size"]))
    _require_unique_retained_paths(paths, "baseline authority record")
    _require_digest_size_function(
        byte_identities, "baseline authority record")
    _require_digest_claim_function(
        byte_claims, "baseline authority record")
    _require_host_file_function(host_files, "baseline authority record")
    _require_tool_path_function(
        toolchain["tools"] + toolchain["subtools"],
        "baseline authority record",
    )
    _require_git_object_consistency(
        source, adapter, "baseline authority record")


def _dependency(value: Any, label: str) -> dict[str, Any]:
    dependency = _dict(value, DEPENDENCY_KEYS, label)
    soname = _text(dependency["soname"], f"{label} soname")
    kind = _token(dependency["kind"], f"{label} kind")
    _require(kind in ("file", "virtual"), f"{label} kind changed")
    if kind == "virtual":
        _require(dependency["path"] is None and
                 dependency["size"] is None and
                 dependency["sha256"] is None,
                 f"{label} virtual dependency has file claims")
        return {"soname": soname, "kind": kind, "path": None,
                "size": None, "sha256": None}
    size, digest = _byte_identity(
        dependency["size"], dependency["sha256"], label, minimum_size=1)
    return {
        "soname": soname,
        "kind": kind,
        "path": _absolute_path(dependency["path"], f"{label} path"),
        "size": size,
        "sha256": digest,
    }


def parse_canonical_ldd_output(data: bytes) -> tuple[dict[str, Any], ...]:
    """Parse the frozen canonical-ldd-C-v1 byte grammar.

    File-backed rows are ``SONAME<TAB>file<TAB>/absolute/path<LF>`` and
    virtual rows are ``SONAME<TAB>virtual<LF>``.  Rows are sorted uniquely by
    SONAME.  File sizes and hashes remain separate record claims because the
    normalized ldd text intentionally contains no mutable host identities.
    """
    _require(type(data) is bytes and 0 < len(data) <= MAX_CANONICAL_LDD_BYTES,
             "canonical ldd output is not bounded non-empty bytes")
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as error:
        raise ExactMainBaselineRecordError(
            "canonical ldd output is not strict ASCII") from error
    _require(text.endswith("\n") and "\r" not in text,
             "canonical ldd output is not LF-framed")
    lines = text[:-1].split("\n")
    _require(0 < len(lines) <= MAX_DEPENDENCIES and all(lines),
             "canonical ldd output row count is invalid")
    records: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        fields = line.split("\t")
        _require(len(fields) in (2, 3),
                 f"canonical ldd row {index} is malformed")
        soname = _text(fields[0], f"canonical ldd row {index} soname")
        _require(all(0x21 <= ord(character) <= 0x7E
                     for character in soname) and "/" not in soname,
                 f"canonical ldd row {index} soname is not portable")
        if len(fields) == 2:
            _require(fields[1] == "virtual",
                     f"canonical ldd row {index} kind changed")
            records.append({"soname": soname, "kind": "virtual", "path": None})
        else:
            _require(fields[1] == "file",
                     f"canonical ldd row {index} kind changed")
            records.append({
                "soname": soname,
                "kind": "file",
                "path": _absolute_path(
                    fields[2], f"canonical ldd row {index} path"),
            })
    sonames = [record["soname"] for record in records]
    _require(sonames == sorted(set(sonames)),
             "canonical ldd rows are not sorted and unique")
    reconstructed = "".join(
        (f'{record["soname"]}\tvirtual\n'
         if record["kind"] == "virtual" else
         f'{record["soname"]}\tfile\t{record["path"]}\n')
        for record in records
    ).encode("ascii")
    _require(reconstructed == data,
             "canonical ldd output is not byte-canonical")
    return tuple(copy.deepcopy(records))


def _runtime_closure(value: Any,
                     builds: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    closure = _dict(value, RUNTIME_KEYS, "exact-main runtime closure")
    _fixed_text(closure["schema"], RUNTIME_CLOSURE_SCHEMA,
                "runtime closure schema")
    normalization = _fixed_text(
        closure["normalization"], CANONICAL_LDD_NORMALIZATION,
        "runtime closure normalization")
    values = closure["records"]
    _require(type(values) is list and len(values) == len(BUILD_ROLES),
             "runtime closure build inventory changed")
    records = []
    for index, (role, raw) in enumerate(zip(BUILD_ROLES, values)):
        record = _dict(
            raw, RUNTIME_RECORD_KEYS, f"runtime closure record {index}")
        _require(record["role"] == role and type(record["role"]) is str,
                 "runtime closure record order changed")
        executable_sha = _sha256(
            record["executable_sha256"],
            f"{role} runtime closure executable SHA-256")
        _require(executable_sha == builds[role]["executable"]["sha256"],
                 f"{role} runtime closure executable changed")
        dependencies_value = record["dependencies"]
        _require(type(dependencies_value) is list and
                 1 <= len(dependencies_value) <= MAX_DEPENDENCIES,
                 f"{role} runtime dependency count is invalid")
        dependencies = [
            _dependency(dependency, f"{role} dependency {dep_index}")
            for dep_index, dependency in enumerate(dependencies_value)
        ]
        sonames = [dependency["soname"] for dependency in dependencies]
        _require(sonames == sorted(set(sonames)),
                 f"{role} runtime dependencies are not sorted and unique")
        file_paths = [
            dependency["path"] for dependency in dependencies
            if dependency["kind"] == "file"
        ]
        _require(len(file_paths) == len(set(file_paths)),
                 f"{role} runtime dependency paths are not unique")
        ldd_output = _text_identity(
            record["canonical_ldd_output"],
            f"{role} canonical ldd output",
            expected_path=f"runtime/{role.replace('_', '-')}/ldd.txt")
        _require(ldd_output["size"] > 0,
                 f"{role} attested no runtime closure output")
        records.append({
            "role": role,
            "executable_sha256": executable_sha,
            "canonical_ldd_output": ldd_output,
            "dependencies": dependencies,
        })
    _require(all(
        exact_json_equal(records[0]["dependencies"], record["dependencies"])
        and records[0]["canonical_ldd_output"]["size"] ==
        record["canonical_ldd_output"]["size"]
        and records[0]["canonical_ldd_output"]["sha256"] ==
        record["canonical_ldd_output"]["sha256"]
        for record in records[1:]
    ), "canonical runtime closure changed between builds")
    return {
        "schema": RUNTIME_CLOSURE_SCHEMA,
        "normalization": normalization,
        "records": records,
    }


def _attestation(value: Any, builds: Mapping[str, Mapping[str, Any]],
                 toolchain: Mapping[str, Any],
                 adapter: Mapping[str, Any]) -> dict[str, Any]:
    attestation = _dict(value, ATTESTATION_KEYS, "exact-main attestation")
    _fixed_text(attestation["schema"], ATTESTATION_SCHEMA,
                "exact-main attestation schema")
    tools = {record["role"]: record["resolved_path"]
             for record in toolchain["tools"]}
    values = attestation["records"]
    _require(type(values) is list and len(values) == len(BUILD_ROLES),
             "exact-main attestation build inventory changed")
    records = []
    for index, (role, raw) in enumerate(zip(BUILD_ROLES, values)):
        record = _dict(
            raw, ATTESTATION_RECORD_KEYS, f"attestation record {index}")
        _require(record["role"] == role and type(record["role"]) is str,
                 "attestation record order changed")
        build = builds[role]
        expected_argv = exact_main_benchmark_argv(executable_path=(
            build["roots"]["build_root"] + "/" +
            build["executable"]["build_relative_path"]))
        argv = _arguments(record["argv"], f"{role} attestation argv")
        _require(argv == expected_argv, f"{role} attestation argv changed")
        _require(record["exit_status"] == 0 and
                 type(record["exit_status"]) is int,
                 f"{role} attestation command failed")
        _require(record["reported_schema"] == BENCHMARK_SCHEMA and
                 type(record["reported_schema"]) is str,
                 f"{role} benchmark schema changed")
        _require(record["main_source_commit"] == BASELINE_COMMIT and
                 type(record["main_source_commit"]) is str,
                 f"{role} attestation source commit changed")
        _require(type(record["pure_avx2"]) is bool and
                 record["pure_avx2"] is True,
                 f"{role} did not attest pure AVX2")
        _require(type(record["round_trip"]) is bool and
                 record["round_trip"] is True,
                 f"{role} did not attest correctness")
        stdout = _text_identity(
            record["stdout"], f"{role} attestation stdout",
            expected_path=(
                f"attestations/{role.replace('_', '-')}/benchmark.stdout.json"))
        stderr = _text_identity(
            record["stderr"], f"{role} attestation stderr",
            expected_path=(
                f"attestations/{role.replace('_', '-')}/benchmark.stderr"))
        _require(stdout["size"] > 0,
                 f"{role} attested no benchmark JSON")
        ctest_raw = _dict(record["ctest"], CTEST_KEYS,
                          f"{role} CTest attestation")
        ctest_argv = _arguments(
            ctest_raw["argv"], f"{role} CTest argv")
        expected_ctest = exact_main_ctest_argv(
            ctest=tools["ctest"], build_root=build["roots"]["build_root"])
        _require(ctest_argv == expected_ctest,
                 f"{role} CTest command changed")
        _require(ctest_raw["exit_status"] == 0 and
                 type(ctest_raw["exit_status"]) is int and
                 ctest_raw["passed"] == 1 and
                 type(ctest_raw["passed"]) is int and
                 ctest_raw["failed"] == 0 and
                 type(ctest_raw["failed"]) is int,
                 f"{role} CTest did not pass exactly one smoke test")
        ctest_stdout = _text_identity(
            ctest_raw["stdout"], f"{role} CTest stdout",
            expected_path=(
                f"attestations/{role.replace('_', '-')}/ctest.stdout.log"))
        ctest_stderr = _text_identity(
            ctest_raw["stderr"], f"{role} CTest stderr",
            expected_path=(
                f"attestations/{role.replace('_', '-')}/ctest.stderr.log"))
        _require(ctest_stdout["size"] > 0,
                 f"{role} CTest emitted no success record")
        records.append({
            "role": role,
            "argv": argv,
            "stdout": stdout,
            "stderr": stderr,
            "exit_status": 0,
            "reported_schema": BENCHMARK_SCHEMA,
            "main_source_commit": BASELINE_COMMIT,
            "pure_avx2": True,
            "round_trip": True,
            "ctest": {
                "argv": ctest_argv,
                "stdout": ctest_stdout,
                "stderr": ctest_stderr,
                "exit_status": 0,
                "passed": 1,
                "failed": 0,
            },
        })
    controller = _file(
        attestation["test_controller"],
        "exact-main attestation controller")
    adapter_controller = adapter["files"][2]
    _require(
        controller["relative_path"] ==
        "controllers/test_legacy_main_benchmark.py" and
        controller["size"] == adapter_controller["size"] and
        controller["sha256"] == adapter_controller["sha256"],
        "exact-main attestation controller differs from the adapter source",
    )
    return {
        "schema": ATTESTATION_SCHEMA,
        "test_controller": controller,
        "records": records,
    }


def _identity(value: Any, builds: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    identity = _dict(value, IDENTITY_KEYS, "exact-main identity group")
    records = {}
    for role in BUILD_ROLES:
        try:
            record = validate_normalized_code_identity(identity[role])
        except ExactMainBaselineError as error:
            raise ExactMainBaselineRecordError(
                f"{role} normalized identity is invalid") from error
        artifact = record["artifact"]
        _require(artifact["size"] == builds[role]["executable"]["size"] and
                 artifact["sha256"] == builds[role]["executable"]["sha256"],
                 f"{role} normalized identity is not bound to its executable")
        census_roots = {
            item["role"]: item["path"]
            for item in record["path_string_census"]["roots"]
        }
        _require(exact_json_equal(census_roots, builds[role]["roots"]),
                 f"{role} normalized identity roots changed")
        records[role] = record
    canonical_equal = exact_json_equal(
        records["canonical_first"], records["canonical_second"])
    variant_differs = (
        records["path_variant"]["artifact"]["sha256"] !=
        records["canonical_first"]["artifact"]["sha256"])
    combined = _sha256(identity["combined_sha256"],
                       "exact-main combined identity SHA-256")
    normalized_match = all(
        records[role]["combined_sha256"] == combined for role in BUILD_ROLES)
    _require(type(identity["canonical_raw_bytes_identical"]) is bool and
             identity["canonical_raw_bytes_identical"] is canonical_equal,
             "same-path normalized raw identity claim is inconsistent")
    _require(type(identity["path_variant_raw_bytes_differ"]) is bool and
             identity["path_variant_raw_bytes_differ"] is variant_differs,
             "path-variant raw identity claim is inconsistent")
    _require(type(identity["normalized_match"]) is bool and
             identity["normalized_match"] is normalized_match,
             "normalized identity match claim is inconsistent")
    _require(canonical_equal,
             "same-path normalized identities were not byte-identical")
    _require(variant_differs,
             "path-variant raw executable did not differ")
    _require(normalized_match,
             "path-variant normalized identities did not match")
    return {
        **records,
        "combined_sha256": combined,
        "canonical_raw_bytes_identical": True,
        "path_variant_raw_bytes_differ": True,
        "normalized_match": True,
    }


def _census_zero(identity: Mapping[str, Any]) -> bool:
    return all(
        root["occurrences"] == 0 and all(
            row["occurrences"] == 0 for row in root["sections"])
        for role in BUILD_ROLES
        for root in identity[role]["path_string_census"]["roots"]
    )


def superseded_historical_references() -> list[dict[str, Any]]:
    """Return the fixed historical hashes with explicit non-authority."""
    return [
        {
            "label": "historical_archive_sha256",
            "value": HISTORICAL_ARCHIVE_SHA256,
            "authority": False,
        },
        {
            "label": "historical_executable_sha256",
            "value": HISTORICAL_EXECUTABLE_SHA256,
            "authority": False,
        },
    ]


def _references(value: Any) -> list[dict[str, Any]]:
    expected = superseded_historical_references()
    _require(type(value) is list and len(value) == len(expected),
             "superseded historical reference inventory changed")
    for index, raw in enumerate(value):
        _dict(raw, REFERENCE_KEYS, f"historical reference {index}")
    _require(exact_json_equal(value, expected),
             "historical hashes were relabeled as authority")
    return expected


def _promotion(value: Any, *, builds: Mapping[str, Mapping[str, Any]],
               identity: Mapping[str, Any], attestation: Mapping[str, Any],
               references: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    promotion = _dict(value, PROMOTION_KEYS, "baseline promotion")
    executable_equal = (
        builds["canonical_first"]["executable"]["size"] ==
        builds["canonical_second"]["executable"]["size"] and
        builds["canonical_first"]["executable"]["sha256"] ==
        builds["canonical_second"]["executable"]["sha256"])
    archive_equal = (
        builds["canonical_first"]["archive"]["size"] ==
        builds["canonical_second"]["archive"]["size"] and
        builds["canonical_first"]["archive"]["sha256"] ==
        builds["canonical_second"]["archive"]["sha256"])
    variant_differs = (
        builds["path_variant"]["executable"]["sha256"] !=
        builds["canonical_first"]["executable"]["sha256"])
    normalized_match = identity["normalized_match"] is True
    census_zero = _census_zero(identity)
    attested = all(
        record["pure_avx2"] is True and
        record["main_source_commit"] == BASELINE_COMMIT and
        record["round_trip"] is True and
        record["ctest"]["passed"] == 1 and record["ctest"]["failed"] == 0
        for record in attestation["records"]
    )
    non_authority = all(reference["authority"] is False
                        for reference in references)
    expected = {
        "gate": PROMOTION_GATE,
        "same_path_executable_bytes_identical": executable_equal,
        "same_path_archive_bytes_identical": archive_equal,
        "path_variant_raw_executable_differs": variant_differs,
        "path_variant_normalized_match": normalized_match,
        "selected_section_census_zero": census_zero,
        "pure_avx2_attested": attested,
        "historical_references_non_authoritative": non_authority,
        "promoted": all((
            executable_equal, archive_equal, variant_differs,
            normalized_match, census_zero, attested, non_authority,
        )),
    }
    _require(exact_json_equal(promotion, expected),
             "baseline promotion claims are inconsistent")
    _require(expected["promoted"],
             "baseline authority record did not satisfy every promotion gate")
    return copy.deepcopy(expected)


def _record_sha256(value: Mapping[str, Any]) -> str:
    try:
        return canonical_sha256({
            key: copy.deepcopy(field) for key, field in value.items()
            if key != "record_sha256"
        })
    except ExactMainBaselineError as error:
        raise ExactMainBaselineRecordError(
            "baseline record is not canonical finite JSON") from error


def _successful_promotion() -> dict[str, Any]:
    return {
        "gate": PROMOTION_GATE,
        "same_path_executable_bytes_identical": True,
        "same_path_archive_bytes_identical": True,
        "path_variant_raw_executable_differs": True,
        "path_variant_normalized_match": True,
        "selected_section_census_zero": True,
        "pure_avx2_attested": True,
        "historical_references_non_authoritative": True,
        "promoted": True,
    }


def validate_baseline_authority_record(value: Any) -> dict[str, Any]:
    """Validate and detach one promoted exact-main authority record."""
    record = _dict(value, AUTHORITY_KEYS, "baseline authority record")
    _fixed_text(record["schema"], AUTHORITY_SCHEMA,
                "baseline authority schema")
    created = _timestamp(record["created_utc"], "baseline creation time")
    _fixed_text(record["bead"], BEAD_ID, "baseline authority bead")
    host = _host(record["host"])
    lane = _lane(record["lane"], successful=True)
    source = _source(record["source"])
    adapter = _adapter(record["adapter"])
    toolchain = _toolchain(record["toolchain"])
    configure = _profile(record["configure"])
    builds = _builds(record["builds"], toolchain, lane["root"])
    runtime = _runtime_closure(record["runtime_closure"], builds)
    attestation = _attestation(
        record["attestation"], builds, toolchain, adapter)
    identity = _identity(record["identity"], builds)
    _authority_evidence_inventory(
        lane, source, adapter, toolchain, builds, runtime, attestation,
        identity)
    references = _references(record["superseded_references"])
    promotion = _promotion(
        record["promotion"], builds=builds, identity=identity,
        attestation=attestation, references=references)
    canonical = {
        "schema": AUTHORITY_SCHEMA,
        "created_utc": created,
        "bead": BEAD_ID,
        "host": host,
        "lane": lane,
        "source": source,
        "adapter": adapter,
        "toolchain": toolchain,
        "configure": configure,
        "builds": builds,
        "runtime_closure": runtime,
        "attestation": attestation,
        "identity": identity,
        "promotion": promotion,
        "superseded_references": references,
        "record_sha256": _sha256(
            record["record_sha256"], "baseline authority record SHA-256"),
    }
    _require(canonical["record_sha256"] == _record_sha256(canonical),
             "baseline authority record SHA-256 is inconsistent")
    return copy.deepcopy(canonical)


def baseline_authority_record(
    *,
    created_utc: str,
    host: dict[str, Any],
    lane: dict[str, Any],
    source: dict[str, Any],
    adapter: dict[str, Any],
    toolchain: dict[str, Any],
    builds: dict[str, Any],
    runtime_closure: dict[str, Any],
    attestation: dict[str, Any],
    identity: dict[str, Any],
) -> dict[str, Any]:
    """Construct a promoted v1 authority record from acquired facts."""
    value = {
        "schema": AUTHORITY_SCHEMA,
        "created_utc": copy.deepcopy(created_utc),
        "bead": BEAD_ID,
        "host": copy.deepcopy(host),
        "lane": copy.deepcopy(lane),
        "source": copy.deepcopy(source),
        "adapter": copy.deepcopy(adapter),
        "toolchain": copy.deepcopy(toolchain),
        "configure": exact_main_build_profile(),
        "builds": copy.deepcopy(builds),
        "runtime_closure": copy.deepcopy(runtime_closure),
        "attestation": copy.deepcopy(attestation),
        "identity": copy.deepcopy(identity),
        "promotion": _successful_promotion(),
        "superseded_references": superseded_historical_references(),
        "record_sha256": "0" * 64,
    }
    value["record_sha256"] = _record_sha256(value)
    return validate_baseline_authority_record(value)


def load_baseline_authority_record(data: bytes) -> dict[str, Any]:
    """Strictly load one promoted exact-main authority record."""
    try:
        value = strict_json_loads(data, "exact-main authority record JSON")
    except ExactMainBaselineError as error:
        raise ExactMainBaselineRecordError(
            "exact-main authority record is not strict JSON") from error
    return validate_baseline_authority_record(value)


def authority_retained_inventory(
    record: Any,
) -> tuple[dict[str, Any], ...]:
    """Return the one validated authority-lane file inventory.

    The terminal record necessarily carries ``None`` for its own byte identity;
    canonical terminal bytes and the outer checksum ledger close that entry.
    """
    validated = validate_baseline_authority_record(record)
    return _authority_retained_inventory_from_validated(
        validated["lane"], validated["source"], validated["toolchain"],
        validated["builds"], validated["runtime_closure"],
        validated["attestation"],
    )


def _failure_retained_inventory_from_validated(
    lane: Mapping[str, Any],
    retained: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    inventory: list[dict[str, Any]] = [{
        "relative_path": lane["record_relative_path"],
        "size": None,
        "sha256": None,
    }]
    inventory.extend(
        _retained_inventory_entry(stage["log"])
        for stage in lane["stages"])
    inventory.extend(_retained_inventory_entry(item) for item in retained)
    inventory.sort(key=lambda item: item["relative_path"])
    return tuple(copy.deepcopy(inventory))


def _failure_record(value: Any, schema: str) -> dict[str, Any]:
    record = _dict(value, FAILURE_KEYS, "baseline failure record")
    _fixed_text(record["schema"], schema, "baseline failure schema")
    created = _timestamp(record["created_utc"], "baseline failure time")
    _fixed_text(record["bead"], BEAD_ID, "baseline failure bead")
    lane = _lane(record["lane"], successful=False)
    stage = _fixed_text(
        record["stage"], lane["stages"][-1]["name"],
        "baseline failure stage")
    error_value = _dict(record["error"], ERROR_KEYS,
                        "baseline failure error")
    error = {
        "kind": _token(error_value["kind"], "baseline failure kind"),
        "message": _text(error_value["message"], "baseline failure message"),
        "exit_status": _integer(
            error_value["exit_status"], 1, 255,
            "baseline failure exit status"),
    }
    retained_value = record["retained_files"]
    _require(type(retained_value) is list and
             len(retained_value) <= MAX_RETAINED_FILES,
             "baseline retained failure file inventory is invalid")
    retained = [_file(
                    item, f"baseline retained failure file {index}",
                    minimum_size=0)
                for index, item in enumerate(retained_value)]
    paths = [item["relative_path"] for item in retained]
    _require(paths == sorted(set(paths)),
             "baseline retained failure files are not sorted and unique")
    _require("baseline-authority.json" not in paths,
             "baseline failure retained a promoted terminal name")
    future_stage_logs = {
        f"logs/{index:02d}-{name}.log"
        for index, name in enumerate(STAGE_NAMES)
        if index >= len(lane["stages"])
    }
    _require(not future_stage_logs.intersection(paths),
             "baseline failure retained a log for a stage that never ran")
    failure_inventory = _failure_retained_inventory_from_validated(
        lane, retained)
    _require_unique_retained_paths(
        [item["relative_path"] for item in failure_inventory],
        "baseline failure record")
    failure_byte_records = (
        [stage_record["log"] for stage_record in lane["stages"]] + retained)
    _require_digest_size_function(
        [(item["sha256"], item["size"]) for item in failure_byte_records],
        "baseline failure record",
    )
    authority_sha = record["authority_record_sha256"]
    if schema == ACQUISITION_FAILURE_SCHEMA:
        _require(stage in STAGE_NAMES[:4],
                 "acquisition failure occurred after authority construction")
        _require(authority_sha is None,
                 "acquisition failure cites nonexistent authority")
        authority_sha = None
    else:
        _require(stage == "independent_verification",
                 "verification failure did not occur during verification")
        authority_sha = _sha256(
            authority_sha, "failed authority record SHA-256")
    _require(type(record["promoted"]) is bool and
             record["promoted"] is False,
             "baseline failure was relabeled as promoted")
    canonical = {
        "schema": schema,
        "created_utc": created,
        "bead": BEAD_ID,
        "lane": lane,
        "stage": stage,
        "error": error,
        "retained_files": retained,
        "authority_record_sha256": authority_sha,
        "promoted": False,
        "record_sha256": _sha256(
            record["record_sha256"], "baseline failure record SHA-256"),
    }
    _require(canonical["record_sha256"] == _record_sha256(canonical),
             "baseline failure record SHA-256 is inconsistent")
    return copy.deepcopy(canonical)


def validate_baseline_failure_record(value: Any) -> dict[str, Any]:
    """Validate either exact failure variant without schema coercion."""
    _require(type(value) is dict, "baseline failure is not a JSON object")
    schema = value.get("schema")
    _require(type(schema) is str and schema in (
        ACQUISITION_FAILURE_SCHEMA, VERIFICATION_FAILURE_SCHEMA),
        "baseline failure schema changed")
    return _failure_record(value, schema)


def failure_retained_inventory(
    record: Any,
) -> tuple[dict[str, Any], ...]:
    """Return the one validated failure-lane file inventory."""
    validated = validate_baseline_failure_record(record)
    return _failure_retained_inventory_from_validated(
        validated["lane"], validated["retained_files"])


def _construct_failure(
    *,
    schema: str,
    created_utc: str,
    lane: dict[str, Any],
    stage: str,
    error: dict[str, Any],
    retained_files: list[dict[str, Any]],
    authority_record_sha256: str | None,
) -> dict[str, Any]:
    value = {
        "schema": schema,
        "created_utc": copy.deepcopy(created_utc),
        "bead": BEAD_ID,
        "lane": copy.deepcopy(lane),
        "stage": copy.deepcopy(stage),
        "error": copy.deepcopy(error),
        "retained_files": copy.deepcopy(retained_files),
        "authority_record_sha256": copy.deepcopy(authority_record_sha256),
        "promoted": False,
        "record_sha256": "0" * 64,
    }
    value["record_sha256"] = _record_sha256(value)
    return validate_baseline_failure_record(value)


def baseline_acquisition_failure_record(
    *,
    created_utc: str,
    lane: dict[str, Any],
    stage: str,
    error: dict[str, Any],
    retained_files: list[dict[str, Any]],
) -> dict[str, Any]:
    return _construct_failure(
        schema=ACQUISITION_FAILURE_SCHEMA,
        created_utc=created_utc,
        lane=lane,
        stage=stage,
        error=error,
        retained_files=retained_files,
        authority_record_sha256=None,
    )


def baseline_verification_failure_record(
    *,
    created_utc: str,
    lane: dict[str, Any],
    stage: str,
    error: dict[str, Any],
    retained_files: list[dict[str, Any]],
    authority_record_sha256: str,
) -> dict[str, Any]:
    return _construct_failure(
        schema=VERIFICATION_FAILURE_SCHEMA,
        created_utc=created_utc,
        lane=lane,
        stage=stage,
        error=error,
        retained_files=retained_files,
        authority_record_sha256=authority_record_sha256,
    )


def load_baseline_failure_record(data: bytes) -> dict[str, Any]:
    try:
        value = strict_json_loads(data, "exact-main failure record JSON")
    except ExactMainBaselineError as error:
        raise ExactMainBaselineRecordError(
            "exact-main failure record is not strict JSON") from error
    return validate_baseline_failure_record(value)


__all__ = (
    "ACQUISITION_FAILURE_SCHEMA",
    "ADAPTER_PATHS",
    "ATTESTATION_SCHEMA",
    "AUTHORITY_SCHEMA",
    "BASELINE_COMMIT",
    "BASELINE_CPP_SOURCES",
    "BASELINE_SSE2NEON_COMMIT",
    "BASELINE_TREE",
    "BENCHMARK_SCHEMA",
    "BEAD_ID",
    "BUILD_CLOSURE_KEYS",
    "BUILD_CLOSURE_SCHEMA",
    "BUILD_PROFILE_SCHEMA",
    "BUILD_ROLES",
    "CANONICAL_LDD_NORMALIZATION",
    "CLOSURE_FILE_KEYS",
    "CTEST_SUMMARY_LINE",
    "ExactMainBaselineRecordError",
    "FAILURE_STAGES",
    "HISTORICAL_ARCHIVE_SHA256",
    "HISTORICAL_EXECUTABLE_SHA256",
    "MINIMUM_HARNESS_COMMIT",
    "MAX_CLOSURE_FILES",
    "MAX_CANONICAL_LDD_BYTES",
    "PROMOTION_GATE",
    "RUNTIME_CLOSURE_SCHEMA",
    "SEAL_PROTOCOL",
    "STAGE_NAMES",
    "SUBTOOL_ROLES",
    "TOOL_ROLES",
    "VERIFICATION_FAILURE_SCHEMA",
    "VERIFIER_VERDICT_SCHEMA",
    "VERSION_ROLES",
    "baseline_acquisition_failure_record",
    "baseline_authority_record",
    "baseline_verification_failure_record",
    "authority_retained_inventory",
    "canonical_json_bytes",
    "exact_main_benchmark_argv",
    "exact_main_build_argv",
    "exact_main_build_cache_requirements",
    "exact_main_build_profile",
    "exact_main_compile_command_argv",
    "exact_main_compile_sources",
    "exact_main_configure_argv",
    "exact_main_ctest_argv",
    "failure_retained_inventory",
    "load_baseline_authority_record",
    "load_baseline_failure_record",
    "parse_canonical_ldd_output",
    "parse_cmake_cache",
    "superseded_historical_references",
    "validate_attestation_stdout",
    "validate_baseline_authority_record",
    "validate_baseline_failure_record",
    "validate_build_closure",
    "validate_cmake_cache",
    "validate_compile_commands",
    "validate_exact_main_build",
)
