#!/usr/bin/env python3
"""Run and replay one authoritative, full-matrix C7 AVX2 benchmark.

This tool never builds C7.  It binds an already-built standalone executable and
its core archive to an exact core commit beneath a clean tooling revision,
validates the final v4 A/B build provenance live, holds both a coordinator
reservation and the shared per-user physical-core lease, pins the child to one
SMT thread, and rejects the run when the sibling records any non-idle work.
``verify`` and ``verify-failure`` are portable: they replay canonical retained
bytes and every derived value without requiring the original source, build,
executable, or reservation paths.  A verified failure remains a failure.
``self-test`` is synthetic and never executes benchmark timing.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import stat
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any, Mapping, Sequence


RAW_SCHEMA = "leopard2-c7-authoritative-raw/v2"
MANIFEST_SCHEMA = "leopard2-c7-authoritative-manifest/v2"
FAILURE_SCHEMA = "leopard2-c7-authoritative-failure/v3"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
ISOLATION_SCHEMA = "leopard2-c7-authoritative-isolation/v1"
BUILD_MANIFEST_SCHEMA = "leopard2-c7-build-run-manifest/v4"
BUILD_ATTESTATION_SCHEMA = "leopard2-c7-authoritative-build-attestation/v1"
CHILD_SCHEMA = "leopard2-c7-exact-low/v1"
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")
SAMPLE_COUNT = 7
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_LOG_BYTES = 1024 * 1024
SOURCE_RELATIVE = Path("experiments/leopard2/non_power_of_two/c7/c7_exact_low.cpp")
RUNNER_RELATIVE = Path("experiments/leopard2/non_power_of_two/c7/run_authoritative.py")
BUILD_RUNNER_RELATIVE = Path(
    "experiments/leopard2/non_power_of_two/c7/run_matrix.py")
BUILD_VALIDATOR_RELATIVE = Path(
    "experiments/leopard2/non_power_of_two/c7/validate_evidence.py")
BUILD_NAMES = ("scalar", "ssse3", "avx2", "auto", "asan-ubsan")
BUILD_SCOPE = (
    "correctness plus one affinity-selected non-authoritative harness "
    "smoke; no promotion timing")
CORE_SOURCE_CLOSURE = (
    "CMakeLists.txt", "Leopard2Backend.cpp", "Leopard2Backend.h",
    "Leopard2BackendAVX2.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h", "Leopard2Dispatch.h", "Leopard2Plan.cpp",
    "Leopard2Plan.h", "LeopardCommon.cpp", "LeopardCommon.h",
    "LeopardFF16.cpp", "LeopardFF16.h", "LeopardFF8.cpp", "LeopardFF8.h",
    "cmake/leopardConfig.cmake.in", "leopard.cpp", "leopard.h",
    "leopard2.cpp", "leopard2.h",
)
CHILD_ENVIRONMENT = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1", "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin", "TZ": "UTC",
}
VALIDATOR_ENVIRONMENT = {
    **CHILD_ENVIRONMENT, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
}
CPU_STAT_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal",
)
CPU_STAT_IDLE_FIELDS = ("idle", "iowait")

EXPECTED_PROFILE = {
    "family": 3, "version": 1, "coordinate_map": 1,
    "systematic": "0..K-1", "parity": "K..K+R-1",
    "production_enabled": False,
}
EXPECTED_CORRECTNESS = {
    "gf8_cases": 9, "gf16_cases": 5, "coefficients": 118717,
    "gf16_vandermonde_coefficients": 1500, "encode_executions": 117,
    "encode_symbol_comparisons": 1030423,
    "subset_encode_executions": 117, "decode_plans": 403,
    "decode_executions": 403, "decode_symbol_comparisons": 272487,
    "maximum_loss_plans": 117, "unavailable_parity_plans": 175,
    "no_loss_null_calls": 14, "parity_rebuilds": 403,
    "odd_gf16_rejections": 10, "overlap_rejections": 59,
    "parity_output_overlap_rejections": 13,
    "restored_output_overlap_rejections": 12,
    "restored_input_overlap_rejections": 20,
    "selected_parity_null_rejections": 14, "survivor_null_rejections": 6,
    "atomic_rejection_bytes_checked": 61570,
    "read_only_input_alias_calls": 13,
    "read_only_input_alias_symbol_comparisons": 2139,
    "decode_read_only_input_alias_calls": 117,
    "decode_read_only_input_alias_symbol_comparisons": 6025,
    "hot_path_allocations": 0, "digest_fnv64": "0xec4179e9f2776a58",
}

# K, R, bytes, batch, losses, exact field, padded field.
EXPECTED_CELLS = (
    (3, 253, 64, 8, 3, 1, 2),
    (3, 253, 1024, 8, 3, 1, 2),
    (3, 253, 65536, 1, 3, 1, 2),
    (16, 240, 1024, 8, 8, 1, 1),
    (64, 192, 65536, 1, 8, 1, 1),
    (100, 156, 1024, 8, 8, 1, 2),
    (127, 129, 65536, 1, 8, 1, 2),
    (128, 128, 1024, 8, 8, 1, 1),
    (192, 64, 1024, 8, 8, 1, 2),
    (248, 8, 65536, 1, 8, 1, 2),
    (129, 100, 1024, 8, 8, 2, 2),
    (1000, 200, 1024, 1, 8, 2, 2),
)
CELL_KEYS = {
    "K", "R", "batch", "bytes", "exact_coefficients", "exact_decode",
    "exact_decode_samples_us", "exact_decode_setup",
    "exact_decode_setup_samples_us", "exact_decode_terms", "exact_encode",
    "exact_encode_samples_us", "exact_field", "exact_setup",
    "exact_setup_samples_us", "losses", "padded_decode",
    "padded_decode_samples_us", "padded_decode_scratch",
    "padded_decode_setup", "padded_decode_setup_samples_us", "padded_encode",
    "padded_encode_samples_us", "padded_encode_scratch", "padded_field",
    "padded_setup", "padded_setup_samples_us",
}
SUMMARY_SAMPLE_PAIRS = (
    ("exact_setup", "exact_setup_samples_us"),
    ("padded_setup", "padded_setup_samples_us"),
    ("exact_decode_setup", "exact_decode_setup_samples_us"),
    ("padded_decode_setup", "padded_decode_setup_samples_us"),
    ("exact_encode", "exact_encode_samples_us"),
    ("padded_encode", "padded_encode_samples_us"),
    ("exact_decode", "exact_decode_samples_us"),
    ("padded_decode", "padded_decode_samples_us"),
)


class EvidenceError(RuntimeError):
    """The run or evidence is not authoritative."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def typed_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            typed_equal(left[key], right[key]) for key in left)
    if isinstance(left, (list, tuple)):
        return len(left) == len(right) and all(
            typed_equal(a, b) for a, b in zip(left, right))
    return bool(left == right)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def validate_utc(value: object, label: str) -> str:
    require(isinstance(value, str) and value.endswith("Z"),
            f"{label} is not canonical UTC")
    try:
        parsed = dt.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise EvidenceError(f"{label} is not valid UTC") from error
    require(parsed.tzinfo is not None and parsed.utcoffset() == dt.timedelta(0),
            f"{label} is not UTC")
    return value


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"),
                          allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def signed(value: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(dict(value))
    require("digest" not in result, "digest field already exists")
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def verify_signature(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    payload = copy.deepcopy(value)
    digest = payload.pop("digest", None)
    require(isinstance(digest, str) and SHA256_RE.fullmatch(digest) is not None,
            f"{label} has no canonical digest")
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            f"{label} digest mismatch")
    return value


def strict_json(data: bytes, label: str) -> Any:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in items:
            if key in result:
                raise EvidenceError(f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise EvidenceError(f"{label} contains non-finite number {value}")

    def finite_float(value: str) -> float:
        parsed = float(value)
        if not math.isfinite(parsed):
            raise EvidenceError(f"{label} contains non-finite number {value}")
        return parsed

    try:
        return json.loads(data.decode("utf-8"), object_pairs_hook=pairs,
                          parse_constant=reject_constant, parse_float=finite_float)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise EvidenceError(f"invalid {label}: {error}") from error


def canonical_pretty_bytes(value: object) -> bytes:
    try:
        return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) +
                "\n").encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical pretty JSON: {error}") from error


def read_bounded(path: Path, maximum_bytes: int = MAX_ARTIFACT_BYTES) -> bytes:
    with path.open("rb") as stream:
        data = stream.read(maximum_bytes + 1)
    require(len(data) <= maximum_bytes, f"{path.name} exceeds the evidence bound")
    return data


def write_exclusive(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace {path}") from error


def write_json_exclusive(path: Path, value: object) -> None:
    write_exclusive(path, canonical_bytes(value) + b"\n")


def run_checked(arguments: Sequence[str], cwd: Path | None = None) -> bytes:
    completed = subprocess.run(list(arguments), cwd=None if cwd is None else str(cwd),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               check=False, timeout=30, env=CHILD_ENVIRONMENT)
    if completed.returncode:
        raise EvidenceError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}: "
            + completed.stderr.decode("utf-8", errors="replace").strip())
    return completed.stdout


def file_identity(path: Path, kind: str) -> dict[str, Any]:
    path = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        info = os.fstat(descriptor)
        require(stat.S_ISREG(info.st_mode), f"{kind} is not a regular file")
        require(0 < info.st_size <= MAX_ARTIFACT_BYTES,
                f"{kind} size is outside the evidence bound")
        mode = stat.S_IMODE(info.st_mode)
        if kind in ("runner", "executable", "taskset", "python"):
            require(mode & 0o111 != 0, f"{kind} is not executable")
        digest = hashlib.sha256()
        prefix = b""
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            if len(prefix) < 8:
                prefix += chunk[:8 - len(prefix)]
            digest.update(chunk)
        if kind == "archive":
            require(prefix == b"!<arch>\n", "C7 archive is not ar format")
        return {"kind": kind, "name": path.name, "size": info.st_size,
                "mode": mode, "sha256": digest.hexdigest()}
    finally:
        os.close(descriptor)


def committed_file_identity(root: Path, commit: str, relative: Path,
                            kind: str) -> dict[str, Any]:
    path = (root / relative).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise EvidenceError(f"{kind} escapes the source checkout") from error
    record = file_identity(path, kind)
    committed = run_checked(("git", "-C", str(root), "show",
                             f"{commit}:{relative.as_posix()}"))
    require(record["sha256"] == sha256_bytes(committed),
            f"{kind} differs from committed source")
    tree_entry = run_checked(("git", "-C", str(root), "ls-tree", "-z", commit,
                              "--", relative.as_posix()))
    require(tree_entry.endswith(b"\0") and tree_entry.count(b"\0") == 1 and
            b"\t" in tree_entry, f"{kind} committed mode is absent")
    header, retained_path = tree_entry[:-1].split(b"\t", 1)
    fields = header.split()
    require(len(fields) == 3 and fields[1] == b"blob" and
            fields[0] in (b"100644", b"100755") and
            retained_path == relative.as_posix().encode("utf-8"),
            f"{kind} committed mode is invalid")
    require(bool(record["mode"] & 0o111) == (fields[0] == b"100755"),
            f"{kind} executable mode differs from committed source")
    record["path"] = relative.as_posix()
    return record


def git_identity(root: Path, tooling_commit: str,
                 core_commit: str) -> dict[str, Any]:
    root = root.resolve(strict=True)
    require(GIT_SHA_RE.fullmatch(tooling_commit) is not None and
            GIT_SHA_RE.fullmatch(core_commit) is not None,
            "expected tooling/core commit is not full SHA-1")
    head = run_checked(("git", "-C", str(root), "rev-parse", "HEAD")).decode().strip()
    require(head == tooling_commit,
            f"source HEAD {head} differs from tooling commit {tooling_commit}")
    status_text = run_checked(("git", "-C", str(root), "status", "--porcelain=v1",
                               "--untracked-files=all")).decode()
    require(status_text == "", f"source checkout is not clean: {status_text!r}")
    ancestor = subprocess.run(
        ["git", "-C", str(root), "merge-base", "--is-ancestor",
         core_commit, tooling_commit], stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, check=False, timeout=30,
        env=CHILD_ENVIRONMENT)
    require(ancestor.returncode == 0,
            "core commit is not an ancestor of tooling commit")
    tooling_tree = run_checked(("git", "-C", str(root), "rev-parse",
                                f"{tooling_commit}^{{tree}}")).decode().strip()
    core_tree = run_checked(("git", "-C", str(root), "rev-parse",
                             f"{core_commit}^{{tree}}")).decode().strip()
    listing = run_checked(("git", "-C", str(root), "ls-tree", "-r", "-z", "HEAD"))
    return {"tooling_commit": head, "tooling_tree": tooling_tree,
            "core_commit": core_commit, "core_tree": core_tree,
            "core_is_ancestor": True, "clean": True,
            "tracked_tree_sha256": sha256_bytes(listing)}


def manifest_artifact(value: object, expected_path: Path,
                      label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"bytes", "path", "sha256"},
            f"{label} build artifact schema changed")
    require(type(value["bytes"]) is int and 0 < value["bytes"] <= MAX_ARTIFACT_BYTES and
            value["path"] == expected_path.as_posix() and
            isinstance(value["sha256"], str) and
            SHA256_RE.fullmatch(value["sha256"]) is not None,
            f"{label} build artifact identity changed")
    return value


def artifact_matches_identity(record: Mapping[str, Any],
                              identity: Mapping[str, Any], label: str) -> None:
    require(record["bytes"] == identity["size"] and
            record["sha256"] == identity["sha256"] and
            Path(record["path"]).name == identity["name"],
            f"{label} differs from the supplied exact file")


def resolve_build_artifact(root: Path, record: Mapping[str, Any], label: str) -> Path:
    relative = record["path"]
    require(isinstance(relative, str) and relative and "\\" not in relative and
            ":" not in relative and not os.path.isabs(relative),
            f"{label} build path is not canonical")
    pure = Path(relative)
    require(pure.as_posix() == relative and
            all(part not in ("", ".", "..") for part in pure.parts),
            f"{label} build path is not canonical")
    candidate = root.joinpath(*pure.parts).resolve(strict=True)
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise EvidenceError(f"{label} build path escapes the checkout") from error
    return candidate


def derive_build_attestation(
    value: object, manifest_identity: Mapping[str, Any],
    expected_tooling_commit: str, expected_core_commit: str,
    source_identity: Mapping[str, Any], build_runner_identity: Mapping[str, Any],
    build_validator_identity: Mapping[str, Any], archive_identity: Mapping[str, Any],
    executable_identity: Mapping[str, Any], *, source_root: Path | None = None,
    archive_path: Path | None = None, executable_path: Path | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict), "C7 build manifest is not an object")
    require(value.get("schema") == BUILD_MANIFEST_SCHEMA and
            value.get("status") == "pass" and value.get("scope") == BUILD_SCOPE and
            value.get("tooling_git_sha") == expected_tooling_commit and
            value.get("core_git_sha") == expected_core_commit,
            "C7 build manifest version, status, scope, or Git identity changed")
    source = manifest_artifact(value.get("source"), SOURCE_RELATIVE, "C7 source")
    build_runner = manifest_artifact(
        value.get("runner"), BUILD_RUNNER_RELATIVE, "C7 build runner")
    build_validator = manifest_artifact(
        value.get("validator"), BUILD_VALIDATOR_RELATIVE, "C7 build validator")
    artifact_matches_identity(source, source_identity, "C7 source")
    artifact_matches_identity(build_runner, build_runner_identity, "C7 build runner")
    artifact_matches_identity(
        build_validator, build_validator_identity, "C7 build validator")

    builds = value.get("builds")
    require(isinstance(builds, list) and
            [item.get("name") if isinstance(item, dict) else None for item in builds] ==
            list(BUILD_NAMES), "C7 build matrix changed")
    avx2 = builds[2]
    require(avx2.get("name") == "avx2" and avx2.get("backend") == "avx2" and
            avx2.get("sanitizer") is False,
            "C7 authoritative build is not the non-sanitized AVX2 build")
    build_dir = avx2.get("build_dir")
    require(isinstance(build_dir, str) and
            build_dir.endswith("/core-avx2") and "\\" not in build_dir and
            ":" not in build_dir and not os.path.isabs(build_dir) and
            Path(build_dir).as_posix() == build_dir and
            all(part not in ("", ".", "..") for part in Path(build_dir).parts),
            "C7 AVX2 build directory changed")
    library = manifest_artifact(
        avx2.get("library"), Path(build_dir) / "liblibleopard.a", "C7 AVX2 archive")
    executable_record = manifest_artifact(
        avx2.get("executable"), Path(build_dir).parent / "c7-avx2",
        "C7 AVX2 executable")
    artifact_matches_identity(library, archive_identity, "C7 AVX2 archive")
    artifact_matches_identity(
        executable_record, executable_identity, "C7 AVX2 executable")
    if source_root is not None:
        require(archive_path is not None and executable_path is not None,
                "local C7 build paths are absent")
        require(resolve_build_artifact(source_root, library, "C7 AVX2 archive") ==
                archive_path.resolve(strict=True) and
                resolve_build_artifact(source_root, executable_record,
                                       "C7 AVX2 executable") ==
                executable_path.resolve(strict=True),
                "supplied C7 AVX2 paths differ from the verified build manifest")

    compile_argv = avx2.get("compile_argv")
    require(isinstance(compile_argv, list) and
            all(isinstance(item, str) for item in compile_argv) and
            compile_argv.count("-O2") == 1 and compile_argv.count("-fopenmp") == 1 and
            not any(re.fullmatch(r"-O(?:0|1|3|g|fast|s|z)", item)
                    for item in compile_argv) and
            not any(item.startswith("-fsanitize") or
                    item == "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1"
                    for item in compile_argv) and
            source["path"] in compile_argv and library["path"] in compile_argv and
            "-o" in compile_argv and
            compile_argv.index("-o") + 1 < len(compile_argv) and
            compile_argv[compile_argv.index("-o") + 1] == executable_record["path"],
            "C7 AVX2 exact compile command is not optimized and non-sanitized")
    instrumentation = avx2.get("instrumentation")
    zero_counts = {"asan_lines": 0, "ubsan_lines": 0}
    require(isinstance(instrumentation, dict) and
            instrumentation.get("required_compile_macro") is False and
            typed_equal(instrumentation.get("executable_counts"), zero_counts) and
            typed_equal(instrumentation.get("core_archive_counts"), zero_counts),
            "C7 AVX2 sanitizer attribution changed")
    reproducibility = value.get("reproducibility")
    require(isinstance(reproducibility, dict) and
            isinstance(reproducibility.get("comparison"), dict) and
            reproducibility["comparison"].get("status") == "pass" and
            isinstance(reproducibility.get("fingerprints"), dict) and
            typed_equal(reproducibility["fingerprints"].get("avx2"), {
                "library_sha256": library["sha256"],
                "executable_sha256": executable_record["sha256"],
            }), "C7 v4 A/B reproducibility attestation is absent or changed")
    require(manifest_identity.get("kind") == "build-manifest" and
            type(manifest_identity.get("size")) is int and
            manifest_identity["size"] > 0 and
            isinstance(manifest_identity.get("sha256"), str) and
            SHA256_RE.fullmatch(manifest_identity["sha256"]) is not None,
            "C7 build manifest file identity changed")
    return {
        "schema": BUILD_ATTESTATION_SCHEMA,
        "validation": {
            "checkout_head_required": True,
            "mode": "live",
            "validator_stdout": "C7 evidence validation passed (live)\n",
        },
        "manifest": {"bytes": manifest_identity["size"],
                     "schema": BUILD_MANIFEST_SCHEMA,
                     "sha256": manifest_identity["sha256"]},
        "status": "pass", "scope": BUILD_SCOPE,
        "tooling_commit": expected_tooling_commit,
        "core_commit": expected_core_commit,
        "comparison_status": "pass",
        "source": copy.deepcopy(source),
        "build_runner": copy.deepcopy(build_runner),
        "build_validator": copy.deepcopy(build_validator),
        "avx2": {
            "name": "avx2", "backend": "avx2", "sanitizer": False,
            "library": copy.deepcopy(library),
            "executable": copy.deepcopy(executable_record),
            "compile_argv_sha256": sha256_bytes(canonical_bytes(compile_argv)),
            "optimization": "-O2", "openmp": "-fopenmp",
            "instrumentation": {
                "required_compile_macro": False,
                "executable_counts": copy.deepcopy(zero_counts),
                "core_archive_counts": copy.deepcopy(zero_counts),
            },
        },
    }


def run_build_validator(source_root: Path, build_manifest: Path) -> None:
    validator = (source_root / BUILD_VALIDATOR_RELATIVE).resolve(strict=True)
    completed = subprocess.run(
        [sys.executable, str(validator), str(build_manifest),
         "--source-root", str(source_root), "--evidence-root", str(source_root),
         "--live", "--require-checkout-head"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        timeout=300, env=VALIDATOR_ENVIRONMENT)
    require(completed.returncode == 0,
            "C7 v4 live build validation failed: " +
            completed.stderr.decode("utf-8", errors="replace").strip())
    require(completed.stdout == b"C7 evidence validation passed (live)\n" and
            completed.stderr == b"", "C7 v4 live validator output changed")


def verified_build_attestation(
    source_root: Path, build_manifest: Path, expected_tooling_commit: str,
    expected_core_commit: str, source_identity: Mapping[str, Any],
    build_runner_identity: Mapping[str, Any],
    build_validator_identity: Mapping[str, Any], archive_identity: Mapping[str, Any],
    executable_identity: Mapping[str, Any], archive_path: Path, executable_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    build_manifest = build_manifest.resolve(strict=True)
    try:
        build_manifest.relative_to(source_root.resolve())
    except ValueError as error:
        raise EvidenceError("C7 v4 build manifest is outside the source checkout") from error
    before = file_identity(build_manifest, "build-manifest")
    run_build_validator(source_root, build_manifest)
    after = file_identity(build_manifest, "build-manifest")
    require(typed_equal(before, after), "C7 build manifest changed during live validation")
    data = read_bounded(build_manifest)
    require(len(data) == before["size"] and sha256_bytes(data) == before["sha256"],
            "C7 build manifest changed after live validation")
    parsed = strict_json(data, "C7 v4 build manifest")
    require(data == canonical_pretty_bytes(parsed),
            "C7 v4 build manifest is not canonical pretty JSON")
    attestation = derive_build_attestation(
        parsed, before, expected_tooling_commit, expected_core_commit,
        source_identity, build_runner_identity, build_validator_identity,
        archive_identity, executable_identity, source_root=source_root,
        archive_path=archive_path, executable_path=executable_path)
    return before, attestation


def input_snapshot(source_root: Path, expected_tooling_commit: str,
                   expected_core_commit: str, archive: Path,
                   executable: Path, taskset: Path,
                   build_manifest: Path) -> dict[str, Any]:
    source_root = source_root.resolve(strict=True)
    runner = Path(__file__).resolve(strict=True)
    require(runner == (source_root / RUNNER_RELATIVE).resolve(strict=True),
            "runner is not the committed runner in --source-root")
    core_source_closure = [
        committed_file_identity(
            source_root, expected_core_commit, Path(relative), "core-source")
        for relative in CORE_SOURCE_CLOSURE
    ]
    git = git_identity(source_root, expected_tooling_commit, expected_core_commit)
    runner = committed_file_identity(
        source_root, expected_tooling_commit, RUNNER_RELATIVE, "runner")
    source = committed_file_identity(
        source_root, expected_tooling_commit, SOURCE_RELATIVE, "source")
    build_runner = committed_file_identity(
        source_root, expected_tooling_commit, BUILD_RUNNER_RELATIVE, "build-runner")
    build_validator = committed_file_identity(
        source_root, expected_tooling_commit, BUILD_VALIDATOR_RELATIVE,
        "build-validator")
    archive_identity = file_identity(archive, "archive")
    executable_identity = file_identity(executable, "executable")
    build_manifest_identity, build_attestation = verified_build_attestation(
        source_root, build_manifest, expected_tooling_commit, expected_core_commit,
        source, build_runner, build_validator, archive_identity, executable_identity,
        archive, executable)
    result = {
        "git": git, "runner": runner, "source": source,
        "build_runner": build_runner, "build_validator": build_validator,
        "build_manifest": build_manifest_identity,
        "build_attestation": build_attestation,
        "archive": archive_identity,
        "executable": executable_identity,
        "taskset": file_identity(taskset, "taskset"),
        "python": file_identity(Path(sys.executable), "python"),
        "core_source_closure": core_source_closure,
    }
    result["binding_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_input_snapshot(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "git", "runner", "source", "build_runner", "build_validator",
        "build_manifest", "build_attestation", "archive", "executable",
        "taskset", "python", "core_source_closure", "binding_sha256"},
        "input snapshot schema changed")
    payload = {key: value[key] for key in value if key != "binding_sha256"}
    require(isinstance(value["binding_sha256"], str) and
            SHA256_RE.fullmatch(value["binding_sha256"]) is not None and
            value["binding_sha256"] == sha256_bytes(canonical_bytes(payload)),
            "input snapshot binding changed")
    git = value["git"]
    require(isinstance(git, dict) and set(git) == {
        "tooling_commit", "tooling_tree", "core_commit", "core_tree",
        "core_is_ancestor", "clean", "tracked_tree_sha256"} and
        all(isinstance(git[key], str) and
            GIT_SHA_RE.fullmatch(git[key]) is not None for key in
            ("tooling_commit", "tooling_tree", "core_commit", "core_tree")) and
        git["core_is_ancestor"] is True and
        git["clean"] is True and
        isinstance(git["tracked_tree_sha256"], str) and
        SHA256_RE.fullmatch(git["tracked_tree_sha256"]) is not None,
        "retained Git identity is invalid")
    expected_kinds = {
        "runner": "runner", "source": "source",
        "build_runner": "build-runner", "build_validator": "build-validator",
        "build_manifest": "build-manifest", "archive": "archive",
        "executable": "executable", "taskset": "taskset", "python": "python",
    }
    committed_paths = {
        "runner": RUNNER_RELATIVE, "source": SOURCE_RELATIVE,
        "build_runner": BUILD_RUNNER_RELATIVE,
        "build_validator": BUILD_VALIDATOR_RELATIVE,
    }
    for key in ("runner", "source", "build_runner", "build_validator",
                "build_manifest", "archive", "executable", "taskset", "python"):
        record = value[key]
        expected_fields = {"kind", "name", "size", "mode", "sha256"}
        if key in committed_paths:
            expected_fields.add("path")
        require(isinstance(record, dict) and set(record) == expected_fields and
            record["kind"] == expected_kinds[key] and
            isinstance(record["name"], str) and record["name"] and
            Path(record["name"]).name == record["name"] and
            type(record["size"]) is int and record["size"] > 0 and
            type(record["mode"]) is int and 0 <= record["mode"] <= 0o7777 and
            isinstance(record["sha256"], str) and
            SHA256_RE.fullmatch(record["sha256"]) is not None,
            f"retained {key} identity is invalid")
        if key in ("runner", "executable", "taskset", "python"):
            require(record["mode"] & 0o111 != 0,
                    f"retained {key} is not executable")
        if key in ("build_runner", "build_validator"):
            require(record["mode"] & 0o111 == 0,
                    f"retained {key} differs from committed nonexecutable mode")
        if key in committed_paths:
            expected_path = committed_paths[key]
            require(record["path"] == expected_path.as_posix() and
                    record["name"] == expected_path.name,
                    f"retained {key} committed path changed")
    attestation = value["build_attestation"]
    require(isinstance(attestation, dict) and set(attestation) == {
        "schema", "validation", "manifest", "status", "scope",
        "tooling_commit", "core_commit", "comparison_status", "source",
        "build_runner", "build_validator", "avx2"} and
        attestation["schema"] == BUILD_ATTESTATION_SCHEMA and
        typed_equal(attestation["validation"], {
            "checkout_head_required": True, "mode": "live",
            "validator_stdout": "C7 evidence validation passed (live)\n"}) and
        attestation["status"] == "pass" and attestation["scope"] == BUILD_SCOPE and
        attestation["tooling_commit"] == git["tooling_commit"] and
        attestation["core_commit"] == git["core_commit"] and
        attestation["comparison_status"] == "pass",
        "retained C7 v4 build attestation changed")
    manifest = attestation["manifest"]
    require(isinstance(manifest, dict) and set(manifest) == {
        "bytes", "schema", "sha256"} and
        manifest["schema"] == BUILD_MANIFEST_SCHEMA and
        manifest["bytes"] == value["build_manifest"]["size"] and
        manifest["sha256"] == value["build_manifest"]["sha256"],
        "retained C7 build manifest attestation changed")
    require(typed_equal(attestation["source"], {
        "bytes": value["source"]["size"], "path": SOURCE_RELATIVE.as_posix(),
        "sha256": value["source"]["sha256"]}) and
        typed_equal(attestation["build_runner"], {
            "bytes": value["build_runner"]["size"],
            "path": BUILD_RUNNER_RELATIVE.as_posix(),
            "sha256": value["build_runner"]["sha256"]}) and
        typed_equal(attestation["build_validator"], {
            "bytes": value["build_validator"]["size"],
            "path": BUILD_VALIDATOR_RELATIVE.as_posix(),
            "sha256": value["build_validator"]["sha256"]}),
        "retained build tooling attestation changed")
    avx2 = attestation["avx2"]
    require(isinstance(avx2, dict) and set(avx2) == {
        "name", "backend", "sanitizer", "library", "executable",
        "compile_argv_sha256", "optimization", "openmp", "instrumentation"} and
        avx2["name"] == "avx2" and avx2["backend"] == "avx2" and
        avx2["sanitizer"] is False and avx2["optimization"] == "-O2" and
        avx2["openmp"] == "-fopenmp" and
        isinstance(avx2["compile_argv_sha256"], str) and
        SHA256_RE.fullmatch(avx2["compile_argv_sha256"]) is not None and
        isinstance(avx2["library"], dict) and
        set(avx2["library"]) == {"bytes", "path", "sha256"} and
        avx2["library"]["bytes"] == value["archive"]["size"] and
        avx2["library"]["sha256"] == value["archive"]["sha256"] and
        isinstance(avx2["executable"], dict) and
        set(avx2["executable"]) == {"bytes", "path", "sha256"} and
        avx2["executable"]["bytes"] == value["executable"]["size"] and
        avx2["executable"]["sha256"] == value["executable"]["sha256"] and
        typed_equal(avx2["instrumentation"], {
            "required_compile_macro": False,
            "executable_counts": {"asan_lines": 0, "ubsan_lines": 0},
            "core_archive_counts": {"asan_lines": 0, "ubsan_lines": 0},
        }), "retained non-sanitized AVX2 build attestation changed")
    closure = value["core_source_closure"]
    require(isinstance(closure, list) and len(closure) == len(CORE_SOURCE_CLOSURE),
            "core source closure length changed")
    for record, expected in zip(closure, CORE_SOURCE_CLOSURE):
        expected_path = Path(expected)
        require(isinstance(record, dict) and set(record) == {
            "kind", "name", "path", "size", "mode", "sha256"} and
            record["kind"] == "core-source" and
            record["path"] == expected_path.as_posix() and
            record["name"] == expected_path.name and
            type(record["size"]) is int and record["size"] > 0 and
            type(record["mode"]) is int and 0 <= record["mode"] <= 0o7777 and
            isinstance(record["sha256"], str) and
            SHA256_RE.fullmatch(record["sha256"]) is not None,
            f"core source closure changed at {expected}")
    return value


def parse_cpu_list(text: str) -> set[int]:
    result: set[int] = set()
    for item in text.strip().split(","):
        if "-" in item:
            first, last = map(int, item.split("-", 1))
            require(first <= last, "invalid CPU list")
            result.update(range(first, last + 1))
        elif item:
            result.add(int(item))
    return result


def read_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="ascii").strip()
    except FileNotFoundError:
        return None


def cpu_record(cpu: int) -> dict[str, Any]:
    root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
    topology = {name: read_optional(root / "topology" / name) for name in
                ("core_id", "physical_package_id", "die_id",
                 "thread_siblings_list", "core_siblings_list")}
    frequency = {name: read_optional(root / "cpufreq" / name) for name in
                 ("scaling_driver", "scaling_governor",
                  "energy_performance_preference", "scaling_min_freq",
                  "scaling_max_freq", "cpuinfo_min_freq", "cpuinfo_max_freq")}
    return {"cpu": cpu, "online": read_optional(root / "online"),
            "topology": topology, "frequency_policy": frequency}


def host_identity(cpu: int, sibling: int, allowed: set[int]) -> dict[str, Any]:
    uname = platform.uname()
    return {
        "system": {"system": uname.system, "release": uname.release,
                   "version": uname.version, "machine": uname.machine,
                   "python": platform.python_version()},
        "allowed_at_launch": sorted(allowed),
        "timing_cpu": cpu_record(cpu), "sibling_cpu": cpu_record(sibling),
        "turbo": {str(path): read_optional(path) for path in (
            Path("/sys/devices/system/cpu/intel_pstate/no_turbo"),
            Path("/sys/devices/system/cpu/amd_pstate/status"),
            Path("/sys/devices/system/cpu/cpufreq/boost"))},
    }


def validate_host(value: object, cpu: int, sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "system", "allowed_at_launch", "timing_cpu", "sibling_cpu", "turbo"},
        "host identity schema changed")
    system = value["system"]
    require(isinstance(system, dict) and set(system) == {
        "system", "release", "version", "machine", "python"} and
        all(isinstance(item, str) and item for item in system.values()),
        "host system identity is invalid")
    allowed = value["allowed_at_launch"]
    require(isinstance(allowed, list) and
            all(type(item) is int and item >= 0 for item in allowed) and
            allowed == sorted(set(allowed)) and cpu in allowed and sibling in allowed,
            "host launch affinity is invalid")
    topology_keys = {
        "core_id", "physical_package_id", "die_id",
        "thread_siblings_list", "core_siblings_list",
    }
    frequency_keys = {
        "scaling_driver", "scaling_governor", "energy_performance_preference",
        "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
        "cpuinfo_max_freq",
    }
    for key, expected_cpu in (("timing_cpu", cpu), ("sibling_cpu", sibling)):
        record = value[key]
        require(isinstance(record, dict) and set(record) == {
            "cpu", "online", "topology", "frequency_policy"} and
            record["cpu"] == expected_cpu and
            (record["online"] is None or isinstance(record["online"], str)) and
            isinstance(record["topology"], dict) and
            set(record["topology"]) == topology_keys and
            isinstance(record["frequency_policy"], dict) and
            set(record["frequency_policy"]) == frequency_keys and
            all(item is None or isinstance(item, str)
                for item in (*record["topology"].values(),
                             *record["frequency_policy"].values())),
            f"host {key} identity is invalid")
        siblings = record["topology"]["thread_siblings_list"]
        require(isinstance(siblings, str) and
                parse_cpu_list(siblings) == {cpu, sibling},
                f"host {key} SMT topology changed")
    turbo_keys = {
        "/sys/devices/system/cpu/intel_pstate/no_turbo",
        "/sys/devices/system/cpu/amd_pstate/status",
        "/sys/devices/system/cpu/cpufreq/boost",
    }
    require(isinstance(value["turbo"], dict) and
            set(value["turbo"]) == turbo_keys and
            all(item is None or isinstance(item, str)
                for item in value["turbo"].values()),
            "host turbo identity is invalid")
    return value


def validate_topology(cpu: int, sibling: int) -> tuple[set[int], set[int]]:
    require(type(cpu) is int and type(sibling) is int and
            cpu >= 0 and sibling >= 0 and cpu != sibling, "invalid CPU pair")
    allowed = set(os.sched_getaffinity(0))
    require(cpu in allowed and sibling in allowed,
            "timing CPU and sibling must both be in launch affinity")
    text = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list") \
        .read_text(encoding="ascii")
    require(parse_cpu_list(text) == {cpu, sibling}, "requested CPUs are not one SMT pair")
    housekeeping = allowed - {cpu, sibling}
    require(housekeeping, "no housekeeping CPU remains")
    return allowed, housekeeping


def cpu_stat_snapshot(cpu: int) -> dict[str, Any]:
    require(type(cpu) is int and cpu >= 0, "CPU stat identity is invalid")
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            tokens = line.split()
            require(tokens[0] == f"cpu{cpu}" and
                    len(tokens) >= 1 + len(CPU_STAT_FIELDS),
                    f"CPU {cpu} has an incomplete /proc/stat record")
            try:
                values = [int(item) for item in
                          tokens[1:1 + len(CPU_STAT_FIELDS)]]
            except ValueError as error:
                raise EvidenceError(
                    f"CPU {cpu} has a non-integer /proc/stat record") from error
            require(all(value >= 0 for value in values),
                    f"CPU {cpu} has a negative /proc/stat counter")
            fields = dict(zip(CPU_STAT_FIELDS, values))
            return {"cpu": cpu, "fields": fields,
                    "total_jiffies": sum(values)}
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def validate_cpu_stat_snapshot(value: object, expected_cpu: int,
                               label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "cpu", "fields", "total_jiffies"} and
        type(value["cpu"]) is int and value["cpu"] == expected_cpu and
        isinstance(value["fields"], dict) and
        set(value["fields"]) == set(CPU_STAT_FIELDS) and
        all(type(value["fields"][name]) is int and value["fields"][name] >= 0
            for name in CPU_STAT_FIELDS) and
        type(value["total_jiffies"]) is int and
        value["total_jiffies"] == sum(value["fields"].values()),
        f"{label} CPU stat snapshot is invalid")
    return value


def cpu_stat_delta(before: Mapping[str, Any],
                   after: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(before, dict) and isinstance(after, dict) and
            type(before.get("cpu")) is int and type(after.get("cpu")) is int and
            before.get("cpu") == after.get("cpu"),
            "CPU stat snapshots refer to different CPUs")
    expected_cpu = before["cpu"]
    before = validate_cpu_stat_snapshot(before, expected_cpu, "before")
    after = validate_cpu_stat_snapshot(after, expected_cpu, "after")
    before_fields, after_fields = before["fields"], after["fields"]
    deltas: dict[str, int] = {}
    for name in CPU_STAT_FIELDS:
        first, last = before_fields[name], after_fields[name]
        require(type(first) is int and first >= 0 and
                type(last) is int and last >= first,
                f"CPU stat counter {name} moved backwards")
        deltas[name] = last - first
    idle = sum(deltas[name] for name in CPU_STAT_IDLE_FIELDS)
    nonidle = sum(value for name, value in deltas.items()
                  if name not in CPU_STAT_IDLE_FIELDS)
    return {"cpu": before["cpu"], "fields": deltas,
            "idle_jiffies": idle, "nonidle_jiffies": nonidle,
            "total_jiffies": idle + nonidle}


def pair_lease_payload(cpu: int, sibling: int,
                       uid: int | None = None) -> dict[str, Any]:
    retained_uid = os.getuid() if uid is None else uid
    require(type(retained_uid) is int and retained_uid >= 0,
            "CPU pair lease UID is invalid")
    return {"cpus": sorted((cpu, sibling)), "schema": PAIR_LEASE_SCHEMA,
            "uid": retained_uid}


def pair_lease_name(cpu: int, sibling: int, uid: int | None = None) -> str:
    retained_uid = os.getuid() if uid is None else uid
    first, second = sorted((cpu, sibling))
    return f"leopard2-cpu-pair-{retained_uid}-{first}-{second}.lock"


def pair_lease_runtime_root(uid: int | None = None) -> Path:
    retained_uid = os.getuid() if uid is None else uid
    return Path("/run/user") / str(retained_uid)


def pair_lease_directory(uid: int | None = None) -> Path:
    return pair_lease_runtime_root(uid) / "leopard2-cpu-leases"


class PairLease:
    """Serialize Leopard2 evidence runners by one physical CPU pair."""

    def __init__(self, cpu: int, sibling: int,
                 runtime_root: Path | None = None):
        require(type(cpu) is int and type(sibling) is int and
                cpu >= 0 and sibling >= 0 and cpu != sibling,
                "pair lease requires two distinct non-negative CPUs")
        self.cpu = cpu
        self.sibling = sibling
        self.runtime_root = (pair_lease_runtime_root() if runtime_root is None
                             else runtime_root)
        self.root = self.runtime_root / "leopard2-cpu-leases"
        self.path = self.root / pair_lease_name(cpu, sibling)
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None

    def _validate_runtime_root(self) -> os.stat_result:
        metadata = os.lstat(self.runtime_root)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "CPU lease runtime root is not an owned mode-0700 directory")
        return metadata

    def _validate_directory(self) -> os.stat_result:
        self._validate_runtime_root()
        metadata = os.lstat(self.root)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "CPU pair lease directory is not an owned mode-0700 directory")
        return metadata

    def validate_current(self) -> None:
        require(self.descriptor is not None and self.identity is not None,
                "CPU pair lease is not held")
        directory = self._validate_directory()
        descriptor = os.fstat(self.descriptor)
        path = os.lstat(self.path)
        require(stat.S_ISREG(descriptor.st_mode) and
                descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) == (path.st_dev, path.st_ino) and
                (descriptor.st_dev, descriptor.st_ino) ==
                    (self.identity["device"], self.identity["inode"]),
                "CPU pair lease path was replaced or its metadata changed")
        require((directory.st_dev, directory.st_ino) ==
                (self.identity["directory_device"],
                 self.identity["directory_inode"]),
                "CPU pair lease directory was replaced")
        expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        retained = os.read(self.descriptor, 4096)
        require(retained == expected, "CPU pair lease contents changed while held")

    def __enter__(self) -> dict[str, Any]:
        self._validate_runtime_root()
        try:
            self.root.mkdir(mode=0o700)
            os.chmod(self.root, 0o700)
        except FileExistsError:
            pass
        self._validate_directory()
        flags = os.O_RDWR | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_NOFOLLOW", 0)
        try:
            self.descriptor = os.open(
                self.path, flags | os.O_CREAT | os.O_EXCL, 0o600)
            os.fchmod(self.descriptor, 0o600)
        except FileExistsError:
            self.descriptor = os.open(self.path, flags)
        except OSError as error:
            raise EvidenceError(f"cannot open CPU pair lease: {error}") from error
        try:
            metadata = os.fstat(self.descriptor)
            require(stat.S_ISREG(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and metadata.st_nlink == 1 and
                    stat.S_IMODE(metadata.st_mode) == 0o600,
                    "CPU pair lease file has unsafe ownership, type, links, or permissions")
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError(
                    f"physical CPU pair is already leased: {self.path}") from error
            expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
            os.lseek(self.descriptor, 0, os.SEEK_SET)
            retained = os.read(self.descriptor, 4096)
            if not retained:
                require(os.write(self.descriptor, expected) == len(expected),
                        "short write to CPU pair lease")
                os.fsync(self.descriptor)
                retained = expected
            require(retained == expected,
                    "CPU pair lease has unexpected or noncanonical contents")
            directory = self._validate_directory()
            self.identity = {
                "device": metadata.st_dev,
                "directory_device": directory.st_dev,
                "directory_inode": directory.st_ino,
                "inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_pair_wide",
                "path": str(self.path.resolve()),
                "payload": pair_lease_payload(self.cpu, self.sibling),
                "sha256": sha256_bytes(retained),
            }
            self.validate_current()
            return self.identity
        except Exception:
            if self.descriptor is not None:
                try:
                    fcntl.flock(self.descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(self.descriptor)
                    self.descriptor = None
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.descriptor is not None:
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None


def validate_pair_lease_identity(value: object, cpu: int,
                                 sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "device", "directory_device", "directory_inode", "inode", "lock",
        "path", "payload", "sha256"}, "CPU pair lease identity is incomplete")
    payload = value["payload"]
    require(isinstance(payload, dict) and set(payload) == {"cpus", "schema", "uid"} and
            payload["schema"] == PAIR_LEASE_SCHEMA and
            typed_equal(payload["cpus"], sorted((cpu, sibling))) and
            type(payload["uid"]) is int and payload["uid"] >= 0,
            "CPU pair lease payload does not match the campaign")
    expected_path = pair_lease_directory(payload["uid"]) / \
        pair_lease_name(cpu, sibling, payload["uid"])
    require(value["path"] == str(expected_path),
            "CPU pair lease path does not identify the reserved pair")
    require(all(type(value[name]) is int and value[name] >= 0 for name in
                ("device", "directory_device", "directory_inode", "inode")),
            "CPU pair lease filesystem identity is invalid")
    require(value["lock"] == "exclusive_nonblocking_pair_wide" and
            value["sha256"] == sha256_bytes(canonical_bytes(payload)),
            "CPU pair lease lock or digest is invalid")
    return value


def isolation_record(cpu: int, sibling: int, pair_lease: Mapping[str, Any],
                     before_cpu: Mapping[str, Any], after_cpu: Mapping[str, Any],
                     before_sibling: Mapping[str, Any],
                     after_sibling: Mapping[str, Any], started_ns: int,
                     ended_ns: int) -> dict[str, Any]:
    before_cpu = validate_cpu_stat_snapshot(
        before_cpu, cpu, "before timing CPU")
    after_cpu = validate_cpu_stat_snapshot(after_cpu, cpu, "after timing CPU")
    before_sibling = validate_cpu_stat_snapshot(
        before_sibling, sibling, "before sibling CPU")
    after_sibling = validate_cpu_stat_snapshot(
        after_sibling, sibling, "after sibling CPU")
    require(type(started_ns) is int and type(ended_ns) is int,
            "isolation monotonic timestamps are invalid")
    benchmark_delta = cpu_stat_delta(before_cpu, after_cpu)
    sibling_delta = cpu_stat_delta(before_sibling, after_sibling)
    accepted = (ended_ns > started_ns and benchmark_delta["nonidle_jiffies"] > 0 and
                sibling_delta["total_jiffies"] > 0 and
                sibling_delta["nonidle_jiffies"] == 0)
    return {
        "schema": ISOLATION_SCHEMA, "accepted": accepted,
        "timing_cpu": cpu, "sibling_cpu": sibling,
        "before": {"timing_cpu": dict(before_cpu),
                   "sibling_cpu": dict(before_sibling),
                   "monotonic_ns": started_ns},
        "after": {"timing_cpu": dict(after_cpu),
                  "sibling_cpu": dict(after_sibling),
                  "monotonic_ns": ended_ns},
        "delta": {"timing_cpu": benchmark_delta,
                  "sibling_cpu": sibling_delta},
        "duration_ns": ended_ns - started_ns,
        "pair_lease": dict(pair_lease),
        "policy": {"counter_source": "/proc/stat",
                   "idle_fields": list(CPU_STAT_IDLE_FIELDS),
                   "nonidle_fields": [name for name in CPU_STAT_FIELDS
                                      if name not in CPU_STAT_IDLE_FIELDS],
                   "sibling_max_nonidle_jiffies": 0,
                   "timing_min_nonidle_jiffies": 1,
                   "sibling_min_total_jiffies": 1},
    }


def validate_isolation(value: object, cpu: int, sibling: int, *,
                       require_accepted: bool = True) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "accepted", "timing_cpu", "sibling_cpu", "before", "after",
        "delta", "duration_ns", "pair_lease", "policy"},
        "isolation schema changed")
    require(value["schema"] == ISOLATION_SCHEMA and
            value["timing_cpu"] == cpu and value["sibling_cpu"] == sibling,
            "isolation CPU pair changed")
    validate_pair_lease_identity(value["pair_lease"], cpu, sibling)
    before, after = value["before"], value["after"]
    require(isinstance(before, dict) and isinstance(after, dict) and
            set(before) == {"timing_cpu", "sibling_cpu", "monotonic_ns"} and
            set(after) == {"timing_cpu", "sibling_cpu", "monotonic_ns"} and
            type(before["monotonic_ns"]) is int and
            type(after["monotonic_ns"]) is int,
            "isolation snapshots changed")
    validate_cpu_stat_snapshot(before["timing_cpu"], cpu,
                               "before timing CPU")
    validate_cpu_stat_snapshot(after["timing_cpu"], cpu,
                               "after timing CPU")
    validate_cpu_stat_snapshot(before["sibling_cpu"], sibling,
                               "before sibling CPU")
    validate_cpu_stat_snapshot(after["sibling_cpu"], sibling,
                               "after sibling CPU")
    expected = isolation_record(
        cpu, sibling, value["pair_lease"], before["timing_cpu"],
        after["timing_cpu"], before["sibling_cpu"], after["sibling_cpu"],
        before["monotonic_ns"], after["monotonic_ns"])
    require(typed_equal(value, expected), "isolation derivation changed")
    if require_accepted:
        require(value["accepted"] is True,
                "timing CPU was inactive or SMT sibling was active during timing")
    return value


def parse_reservation(raw: bytes, cpu: int, sibling: int) -> dict[str, Any]:
    payload = strict_json(raw, "CPU reservation")
    require(isinstance(payload, dict) and set(payload) == {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema", "status"},
        "CPU reservation fields changed")
    require(raw == canonical_bytes(payload),
            "CPU reservation is not canonical JSON without trailing newline")
    require(payload["schema"] == RESERVATION_SCHEMA and payload["status"] == "held",
            "CPU reservation is stale or not held")
    require(type(payload["benchmark_cpu"]) is int and
            type(payload["reserved_sibling"]) is int and
            payload["benchmark_cpu"] == cpu and payload["reserved_sibling"] == sibling,
            "CPU reservation pair differs")
    require(isinstance(payload["owner"], str) and payload["owner"].strip(),
            "CPU reservation owner is empty")
    require(isinstance(payload["nonce"], str) and 8 <= len(payload["nonce"]) <= 256,
            "CPU reservation nonce is invalid")
    return payload


class Reservation:
    def __init__(self, path: Path, cpu: int, sibling: int):
        self.path = path.resolve(strict=True)
        self.cpu = cpu
        self.sibling = sibling
        self.fd: int | None = None
        self.raw = b""
        self.identity: dict[str, Any] | None = None
        self.device: int | None = None
        self.inode: int | None = None

    def __enter__(self) -> dict[str, Any]:
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_NOFOLLOW", 0)
        self.fd = os.open(self.path, flags)
        try:
            descriptor = os.fstat(self.fd)
            current = os.lstat(self.path)
            require(stat.S_ISREG(descriptor.st_mode) and
                    descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                    stat.S_IMODE(descriptor.st_mode) & 0o022 == 0 and
                    (descriptor.st_dev, descriptor.st_ino) ==
                    (current.st_dev, current.st_ino),
                    "CPU reservation has unsafe ownership, type, links, or permissions")
            self.device, self.inode = descriptor.st_dev, descriptor.st_ino
            fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            os.close(self.fd)
            self.fd = None
            raise EvidenceError("CPU reservation is already locked") from error
        try:
            self.raw = os.read(self.fd, MAX_LOG_BYTES + 1)
            require(len(self.raw) <= MAX_LOG_BYTES,
                    "CPU reservation is unexpectedly large")
            payload = parse_reservation(self.raw, self.cpu, self.sibling)
            self.identity = {"schema": RESERVATION_SCHEMA, "bytes": len(self.raw),
                             "sha256": sha256_bytes(self.raw), "payload": payload,
                             "lock": "exclusive_nonblocking"}
            return self.identity
        except Exception:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            os.close(self.fd)
            self.fd = None
            raise

    def validate_current(self) -> None:
        require(self.fd is not None and self.identity is not None,
                "CPU reservation lock was lost")
        descriptor = os.fstat(self.fd)
        current = os.lstat(self.path)
        require((descriptor.st_dev, descriptor.st_ino) ==
                (current.st_dev, current.st_ino) ==
                (self.device, self.inode) and descriptor.st_uid == os.getuid() and
                descriptor.st_nlink == 1 and stat.S_ISREG(descriptor.st_mode) and
                stat.S_IMODE(descriptor.st_mode) & 0o022 == 0,
                "CPU reservation path was replaced or its metadata changed")
        os.lseek(self.fd, 0, os.SEEK_SET)
        raw = os.read(self.fd, MAX_LOG_BYTES + 1)
        require(raw == self.raw, "CPU reservation changed while locked")
        parse_reservation(raw, self.cpu, self.sibling)

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.fd is not None:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            os.close(self.fd)
            self.fd = None


def validate_reservation_record(value: object, cpu: int, sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "bytes", "sha256", "payload", "lock"} and
        value["schema"] == RESERVATION_SCHEMA and
        value["lock"] == "exclusive_nonblocking", "reservation record changed")
    require(type(value["bytes"]) is int and 0 < value["bytes"] <= MAX_LOG_BYTES and
            isinstance(value["sha256"], str) and
            SHA256_RE.fullmatch(value["sha256"]) is not None,
            "reservation byte identity is invalid")
    raw = canonical_bytes(value["payload"])
    require(value["bytes"] == len(raw) and value["sha256"] == sha256_bytes(raw),
            "reservation bytes or checksum changed")
    require(parse_reservation(raw, cpu, sibling) == value["payload"],
            "reservation semantics changed")
    return value


def finite(value: object, label: str) -> float:
    require(type(value) in (int, float), f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result > 0, f"{label} is not finite positive")
    return result


def validate_summary(cell: Mapping[str, Any], summary_key: str,
                     samples_key: str) -> dict[str, float]:
    samples_value = cell.get(samples_key)
    require(isinstance(samples_value, list) and len(samples_value) == SAMPLE_COUNT,
            f"{samples_key} does not contain seven samples")
    samples = [finite(value, samples_key) for value in samples_value]
    median = statistics.median(samples)
    mad = statistics.median(abs(value - median) for value in samples)
    summary = cell.get(summary_key)
    require(isinstance(summary, dict) and set(summary) == {"median_us", "mad_us"},
            f"{summary_key} schema changed")
    require(type(summary["median_us"]) in (int, float) and
            type(summary["mad_us"]) in (int, float) and
            float(summary["median_us"]) == median and float(summary["mad_us"]) == mad,
            f"{summary_key} differs from raw samples")
    return {"median_us": median, "mad_us": mad}


def validate_child_result(value: object, cpu: int, inputs: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "affinity", "allocation_tracking", "benchmarks", "core_git_sha",
        "correctness", "library_sha256", "omp_dynamic", "omp_num_threads",
        "production_constructor_rejected", "profile", "requested_backend",
        "runtime_backend", "sanitizer", "sanitizer_features", "schema",
        "source_sha256", "status", "timing_scope"}, "C7 result schema changed")
    require(value["schema"] == CHILD_SCHEMA and value["status"] == "pass" and
            typed_equal(value["profile"], EXPECTED_PROFILE) and
            value["production_constructor_rejected"] is True and
            value["timing_scope"] == "candidate-authoritative-requires-runner-manifest" and
            value["requested_backend"] == "avx2" and
            value["runtime_backend"] == "avx2" and
            typed_equal(value["affinity"], [cpu]) and
            value["omp_num_threads"] == "1" and value["omp_dynamic"] == "FALSE" and
            value["sanitizer"] == "none" and
            typed_equal(value["sanitizer_features"], {"address": False, "undefined": False}) and
            value["allocation_tracking"] == "global-new" and
            typed_equal(value["correctness"], EXPECTED_CORRECTNESS),
            "C7 result identity, backend, affinity, or correctness changed")
    require(value["source_sha256"] == inputs["source"]["sha256"] and
            value["library_sha256"] == inputs["archive"]["sha256"] and
            value["core_git_sha"] == inputs["git"]["core_commit"],
            "C7 embedded source/archive/Git identity changed")
    cells = value["benchmarks"]
    require(isinstance(cells, list) and len(cells) == len(EXPECTED_CELLS),
            "C7 result does not contain the full 12-cell matrix")
    normalized = []
    for index, (cell, expected) in enumerate(zip(cells, EXPECTED_CELLS)):
        require(isinstance(cell, dict) and set(cell) == CELL_KEYS,
                f"C7 cell {index} schema changed")
        k, r, byte_count, batch, losses, exact_field, padded_field = expected
        geometry = [cell[key] for key in
                    ("K", "R", "bytes", "batch", "losses", "exact_field", "padded_field")]
        require(typed_equal(geometry, list(expected)), f"C7 cell {index} order/geometry changed")
        require(type(cell["exact_coefficients"]) is int and
                cell["exact_coefficients"] == k * r and
                type(cell["exact_decode_terms"]) is int and
                cell["exact_decode_terms"] == k * losses,
                f"C7 cell {index} exact accounting changed")
        require(type(cell["padded_encode_scratch"]) is int and
                cell["padded_encode_scratch"] >= 0 and
                type(cell["padded_decode_scratch"]) is int and
                cell["padded_decode_scratch"] >= 0,
                f"C7 cell {index} scratch accounting is invalid")
        summaries = {name: validate_summary(cell, name, samples)
                     for name, samples in SUMMARY_SAMPLE_PAIRS}
        normalized.append({"geometry": list(expected), "summaries": summaries})
    return {"schema": "leopard2-c7-authoritative-validation/v1",
            "cell_count": len(cells), "sample_count_per_metric": SAMPLE_COUNT,
            "correctness_checksum": value["correctness"]["digest_fnv64"],
            "cells": normalized}


def artifact_relative_path(relative: object) -> Path:
    require(isinstance(relative, str) and relative and
            "\\" not in relative and ":" not in relative and
            not os.path.isabs(relative),
            "artifact path is not relative")
    relative_path = Path(relative)
    require(relative_path.as_posix() == relative,
            "artifact path is not canonical")
    parts = relative_path.parts
    require(all(part not in ("", ".", "..") for part in parts),
            "artifact path is not canonical")
    return relative_path


def safe_artifact(root: Path, relative: object) -> Path:
    relative_path = artifact_relative_path(relative)
    path = root.joinpath(*relative_path.parts).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise EvidenceError("artifact path escapes evidence root") from error
    return path


def artifact_record(root: Path, path: Path) -> dict[str, Any]:
    size = path.stat().st_size
    require(type(size) is int and 0 <= size <= MAX_ARTIFACT_BYTES,
            "retained artifact is unexpectedly large")
    return {"path": path.relative_to(root).as_posix(), "size": size,
            "sha256": sha256_file(path)}


def validate_artifact_record(record: object, label: str,
                             maximum_bytes: int = MAX_ARTIFACT_BYTES) -> dict[str, Any]:
    require(isinstance(record, dict) and set(record) == {"path", "size", "sha256"},
            f"{label} artifact record changed")
    require(type(record["size"]) is int and 0 <= record["size"] <= maximum_bytes and
            isinstance(record["sha256"], str) and
            SHA256_RE.fullmatch(record["sha256"]) is not None,
            f"{label} artifact identity is invalid")
    artifact_relative_path(record["path"])
    return record


def read_artifact(root: Path, record: object, label: str,
                  maximum_bytes: int = MAX_ARTIFACT_BYTES) -> bytes:
    record = validate_artifact_record(record, label, maximum_bytes)
    path = safe_artifact(root, record["path"])
    require(path.is_file(), f"missing retained {label}")
    with path.open("rb") as stream:
        data = stream.read(maximum_bytes + 1)
    require(record["size"] == len(data) and record["sha256"] == sha256_bytes(data),
            f"retained {label} checksum changed")
    return data


def retain_build_provenance(root: Path, build_manifest: Path,
                            inputs: Mapping[str, Any]) -> dict[str, Any]:
    data = read_bounded(build_manifest.resolve(strict=True))
    identity = inputs["build_manifest"]
    require(len(data) == identity["size"] and
            sha256_bytes(data) == identity["sha256"],
            "C7 build manifest changed before retention")
    path = root / "provenance/build-run-manifest-v4.json"
    write_exclusive(path, data)
    return {"manifest": artifact_record(root, path),
            "attestation": copy.deepcopy(inputs["build_attestation"])}


def validate_build_provenance(value: object, root: Path,
                              inputs: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"manifest", "attestation"},
            "retained C7 build provenance schema changed")
    data = read_artifact(root, value["manifest"], "C7 v4 build manifest")
    require(data == canonical_pretty_bytes(
        strict_json(data, "retained C7 v4 build manifest")),
        "retained C7 v4 build manifest is not canonical pretty JSON")
    require(value["manifest"]["size"] == inputs["build_manifest"]["size"] and
            value["manifest"]["sha256"] == inputs["build_manifest"]["sha256"],
            "retained C7 build manifest differs from the live-validated input")
    parsed = strict_json(data, "retained C7 v4 build manifest")
    expected = derive_build_attestation(
        parsed, inputs["build_manifest"], inputs["git"]["tooling_commit"],
        inputs["git"]["core_commit"], inputs["source"], inputs["build_runner"],
        inputs["build_validator"], inputs["archive"], inputs["executable"])
    require(typed_equal(value["attestation"], expected) and
            typed_equal(inputs["build_attestation"], expected),
            "retained C7 v4 build attestation changed")
    return value


def expected_stderr() -> bytes:
    return b"".join(f"C7 benchmark {index}/12\n".encode("ascii")
                    for index in range(1, 13))


def validate_request(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "backend", "cell_geometry", "child_environment", "command", "cpu",
        "sibling", "timeout_seconds"} and value["backend"] == "avx2" and
        typed_equal(value["cell_geometry"],
                    [list(cell) for cell in EXPECTED_CELLS]) and
        typed_equal(value["child_environment"], CHILD_ENVIRONMENT),
        "authoritative request changed")
    cpu, sibling = value["cpu"], value["sibling"]
    require(type(cpu) is int and type(sibling) is int and
            cpu >= 0 and sibling >= 0 and cpu != sibling and
            type(value["timeout_seconds"]) in (int, float) and
            math.isfinite(value["timeout_seconds"]) and
            value["timeout_seconds"] > 0,
            "authoritative request CPU or timeout is invalid")
    require(typed_equal(value["command"], [
        "${TASKSET}", "-c", str(cpu), "${C7_EXECUTABLE}",
        "--backend", "avx2", "${RESULT_JSON}"]),
        "authoritative child command changed")
    return value


def validate_raw(value: object, root: Path, check_files: bool = True,
                 result_override: object | None = None) -> dict[str, Any]:
    raw = verify_signature(value, "C7 raw evidence")
    require(set(raw) == {"schema", "created_utc", "request", "inputs_before",
                         "inputs_after", "host_before", "host_after", "reservation",
                         "isolation", "build_provenance", "child",
                         "validated_output", "digest"} and
            raw["schema"] == RAW_SCHEMA, "raw evidence schema changed")
    validate_utc(raw["created_utc"], "raw creation time")
    request = validate_request(raw["request"])
    cpu, sibling = request["cpu"], request["sibling"]
    inputs = validate_input_snapshot(raw["inputs_before"])
    require(typed_equal(raw["inputs_after"], inputs), "input identity changed during run")
    validate_build_provenance(raw["build_provenance"], root, inputs)
    validate_host(raw["host_before"], cpu, sibling)
    validate_host(raw["host_after"], cpu, sibling)
    require(typed_equal(raw["host_before"], raw["host_after"]),
            "host topology/frequency policy changed during run")
    validate_reservation_record(raw["reservation"], cpu, sibling)
    validate_isolation(raw["isolation"], cpu, sibling)
    child = raw["child"]
    require(isinstance(child, dict) and set(child) == {
        "command", "environment", "returncode", "timed_out", "duration_ns",
        "stdout", "stderr", "result"} and child["command"] == request["command"] and
        typed_equal(child["environment"], CHILD_ENVIRONMENT) and
        type(child["returncode"]) is int and child["returncode"] == 0 and
        child["timed_out"] is False and
        type(child["duration_ns"]) is int and child["duration_ns"] > 0 and
        child["duration_ns"] == raw["isolation"]["duration_ns"],
        "C7 child failed, timed out, or changed invocation")
    if check_files:
        stdout = read_artifact(root, child["stdout"], "stdout", MAX_LOG_BYTES)
        stderr = read_artifact(root, child["stderr"], "stderr", MAX_LOG_BYTES)
        result_bytes = read_artifact(root, child["result"], "result")
        require(stdout == b"", "C7 child unexpectedly wrote stdout")
        require(stderr == expected_stderr(), "C7 child stderr progress changed")
        parsed = strict_json(result_bytes, "C7 result")
    else:
        require(isinstance(result_override, dict), "parsed C7 result is absent")
        parsed = result_override
    normalized = validate_child_result(parsed, cpu, inputs)
    require(typed_equal(raw["validated_output"], normalized),
            "retained C7 validation summary changed")
    return normalized


def validate_manifest(value: object, root: Path) -> dict[str, Any]:
    manifest = verify_signature(value, "C7 manifest")
    require(set(manifest) == {"schema", "created_utc", "valid", "raw", "request",
                              "inputs", "reservation", "isolation", "build_provenance",
                              "validated_output", "digest"} and
            manifest["schema"] == MANIFEST_SCHEMA and manifest["valid"] is True,
            "C7 manifest schema changed")
    validate_utc(manifest["created_utc"], "manifest creation time")
    raw_bytes = read_artifact(root, manifest["raw"], "raw")
    raw = strict_json(raw_bytes, "C7 raw evidence")
    require(raw_bytes == canonical_bytes(raw) + b"\n",
            "raw evidence is not canonical JSON")
    normalized = validate_raw(raw, root, check_files=True)
    require(manifest["created_utc"] == raw["created_utc"],
            "manifest creation time differs from raw evidence")
    for key, expected in (("request", raw["request"]),
                          ("inputs", raw["inputs_before"]),
                          ("reservation", raw["reservation"]),
                          ("isolation", raw["isolation"]),
                          ("build_provenance", raw["build_provenance"]),
                          ("validated_output", normalized)):
        require(typed_equal(manifest[key], expected), f"manifest {key} differs from raw")
    return normalized


def validate_failure_child(value: object, request: Mapping[str, Any],
                           root: Path | None, check_files: bool) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "command", "environment", "returncode", "timed_out", "duration_ns",
        "stdout", "stderr", "result"},
        "failure child record changed")
    require(typed_equal(value["command"], request["command"]) and
            typed_equal(value["environment"], CHILD_ENVIRONMENT) and
            type(value["returncode"]) is int and
            type(value["timed_out"]) is bool and
            type(value["duration_ns"]) is int and value["duration_ns"] > 0,
            "failure child invocation or outcome is invalid")
    validate_artifact_record(value["stdout"], "failure stdout", MAX_LOG_BYTES)
    validate_artifact_record(value["stderr"], "failure stderr", MAX_LOG_BYTES)
    require(value["result"] is None or isinstance(value["result"], dict),
            "failure result artifact record changed")
    if value["result"] is not None:
        validate_artifact_record(value["result"], "failure result")
    if check_files:
        require(root is not None, "failure artifact root is absent")
        read_artifact(root, value["stdout"], "failure stdout", MAX_LOG_BYTES)
        read_artifact(root, value["stderr"], "failure stderr", MAX_LOG_BYTES)
        if value["result"] is not None:
            read_artifact(root, value["result"], "failure result")
    return value


def classify_failure(child: Mapping[str, Any] | None,
                     error_type: str, error: str) -> str:
    if child is None:
        return "validation"
    if (child["timed_out"] is True and child["returncode"] == 124 and
            error_type == "EvidenceError" and error == "C7 child timed out"):
        return "timeout"
    if (child["timed_out"] is False and child["returncode"] != 0 and
            error_type == "EvidenceError" and
            error == f"C7 child exited {child['returncode']}"):
        return "child-exit"
    return "validation"


def validate_failure(value: object, root: Path | None = None,
                     check_files: bool = False) -> dict[str, Any]:
    failure = verify_signature(value, "C7 failure evidence")
    require(isinstance(failure, dict) and set(failure) == {
        "schema", "status", "failure_kind", "created_utc", "error_type",
        "error", "traceback", "request", "inputs_before", "inputs_after",
        "host_before", "host_after", "reservation", "pair_lease", "isolation",
        "build_provenance", "child", "digest"} and
        failure["schema"] == FAILURE_SCHEMA and failure["status"] == "failed" and
        failure["failure_kind"] in ("validation", "timeout", "child-exit") and
        isinstance(failure["error_type"], str) and failure["error_type"] and
        isinstance(failure["error"], str) and failure["error"] and
        isinstance(failure["traceback"], str) and failure["traceback"],
        "failure evidence schema changed")
    validate_utc(failure["created_utc"], "failure creation time")

    request = failure["request"]
    if request is not None:
        request = validate_request(request)
    inputs_before = failure["inputs_before"]
    if inputs_before is not None:
        inputs_before = validate_input_snapshot(inputs_before)
    inputs_after = failure["inputs_after"]
    if inputs_after is not None:
        require(inputs_before is not None,
                "failure inputs-after has no inputs-before context")
        validate_input_snapshot(inputs_after)

    for key in ("host_before", "host_after"):
        host = failure[key]
        if host is not None:
            require(request is not None,
                    f"failure {key.replace('_', '-')} has no request context")
            validate_host(host, request["cpu"], request["sibling"])
    reservation = failure["reservation"]
    if reservation is not None:
        require(request is not None,
                "failure reservation has no request context")
        validate_reservation_record(
            reservation, request["cpu"], request["sibling"])
    pair_lease = failure["pair_lease"]
    if pair_lease is not None:
        require(request is not None,
                "failure pair lease has no request context")
        validate_pair_lease_identity(
            pair_lease, request["cpu"], request["sibling"])

    provenance = failure["build_provenance"]
    if provenance is not None:
        require(inputs_before is not None and isinstance(provenance, dict) and
                set(provenance) == {"manifest", "attestation"},
                "failure build provenance has no input context or changed schema")
        validate_artifact_record(
            provenance["manifest"], "failure C7 v4 build manifest")
        require(typed_equal(provenance["attestation"],
                            inputs_before["build_attestation"]) and
                provenance["manifest"]["size"] ==
                    inputs_before["build_manifest"]["size"] and
                provenance["manifest"]["sha256"] ==
                    inputs_before["build_manifest"]["sha256"],
                "failure build provenance differs from its input attestation")
        if check_files:
            require(root is not None, "failure provenance root is absent")
            validate_build_provenance(provenance, root, inputs_before)

    isolation = failure["isolation"]
    if isolation is not None:
        require(request is not None,
                "failure isolation has no request context")
        validate_isolation(isolation, request["cpu"], request["sibling"],
                           require_accepted=False)
        require(pair_lease is not None and
                typed_equal(isolation["pair_lease"], pair_lease),
                "failure isolation differs from the held pair lease")

    child = failure["child"]
    if child is not None:
        require(request is not None,
                "failure child has no request context")
        child = validate_failure_child(child, request, root, check_files)
        require(inputs_before is not None and
                failure["host_before"] is not None and reservation is not None and
                pair_lease is not None and provenance is not None,
                "failure child is missing launch provenance")
        if isolation is not None:
            require(child["duration_ns"] == isolation["duration_ns"],
                    "failure child and isolation durations differ")
    if isolation is not None:
        require(child is not None,
                "failure isolation has no child context")

    expected_kind = classify_failure(
        child, failure["error_type"], failure["error"])
    require(failure["failure_kind"] == expected_kind,
            "failure kind differs from the child outcome")
    if expected_kind == "timeout":
        require(child is not None and child["returncode"] == 124 and
                failure["error_type"] == "EvidenceError" and
                failure["error"] == "C7 child timed out",
                "timeout failure status or error relationship changed")
    elif expected_kind == "child-exit":
        require(child is not None and child["timed_out"] is False and
                child["returncode"] != 0 and
                failure["error_type"] == "EvidenceError" and
                failure["error"] ==
                    f"C7 child exited {child['returncode']}",
                "child-exit failure status or error relationship changed")
    else:
        require(failure["error"] != "C7 child timed out" and
                re.fullmatch(r"C7 child exited -?[0-9]+", failure["error"]) is None,
                "validation failure uses a child-outcome error")
    return failure


def run_campaign(options: argparse.Namespace) -> int:
    output = options.output.resolve()
    require(not output.exists(), f"output already exists: {output}")
    output.mkdir(parents=True)
    child_record: dict[str, Any] | None = None
    isolation: dict[str, Any] | None = None
    build_provenance: dict[str, Any] | None = None
    request: dict[str, Any] | None = None
    inputs_before: dict[str, Any] | None = None
    inputs_after: dict[str, Any] | None = None
    host_before: dict[str, Any] | None = None
    host_after: dict[str, Any] | None = None
    reservation_record: dict[str, Any] | None = None
    pair_lease_record: dict[str, Any] | None = None
    original_affinity: set[int] | None = None
    try:
        taskset = Path("/usr/bin/taskset").resolve(strict=True)
        executable = options.executable.resolve(strict=True)
        archive = options.archive.resolve(strict=True)
        build_manifest = options.build_manifest.resolve(strict=True)
        source_root = options.source_root.resolve(strict=True)
        require(type(options.timeout) in (int, float) and
                math.isfinite(options.timeout) and options.timeout > 0,
                "--timeout must be positive and finite")
        candidate_request = {
            "backend": "avx2",
            "cell_geometry": [list(cell) for cell in EXPECTED_CELLS],
            "child_environment": dict(CHILD_ENVIRONMENT), "cpu": options.cpu,
            "sibling": options.sibling, "timeout_seconds": options.timeout,
            "command": ["${TASKSET}", "-c", str(options.cpu),
                        "${C7_EXECUTABLE}", "--backend", "avx2",
                        "${RESULT_JSON}"],
        }
        validate_request(candidate_request)
        request = candidate_request
        original_affinity = set(os.sched_getaffinity(0))
        allowed, housekeeping = validate_topology(options.cpu, options.sibling)
        host_before = host_identity(options.cpu, options.sibling, allowed)
        reservation_guard = Reservation(
            options.reservation_file, options.cpu, options.sibling)
        pair_guard = PairLease(options.cpu, options.sibling)
        with reservation_guard as reservation_record, \
                pair_guard as pair_lease_record:
            os.sched_setaffinity(0, housekeeping)
            inputs_before = input_snapshot(
                source_root, options.expected_tooling_commit,
                options.expected_core_commit, archive, executable, taskset,
                build_manifest)
            build_provenance = retain_build_provenance(
                output, build_manifest, inputs_before)
            result_path = output / "child/result.json"
            stdout_path = output / "child/stdout.bin"
            stderr_path = output / "child/stderr.bin"
            result_path.parent.mkdir(parents=True, exist_ok=True)
            before_cpu = cpu_stat_snapshot(options.cpu)
            before_sibling = cpu_stat_snapshot(options.sibling)
            started = time.monotonic_ns()
            timed_out = False
            command = [str(taskset), "-c", str(options.cpu), str(executable),
                       "--backend", "avx2", str(result_path)]
            try:
                completed = subprocess.run(command, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE, check=False,
                                           timeout=options.timeout,
                                           env=dict(CHILD_ENVIRONMENT))
            except subprocess.TimeoutExpired as error:
                timed_out = True
                completed = subprocess.CompletedProcess(
                    command, 124, error.stdout or b"", error.stderr or b"")
            ended = time.monotonic_ns()
            after_cpu = cpu_stat_snapshot(options.cpu)
            after_sibling = cpu_stat_snapshot(options.sibling)
            write_exclusive(stdout_path, completed.stdout)
            write_exclusive(stderr_path, completed.stderr)
            result_record = (artifact_record(output, result_path)
                             if result_path.is_file() else None)
            child_record = {
                "command": request["command"], "environment": dict(CHILD_ENVIRONMENT),
                "returncode": completed.returncode, "timed_out": timed_out,
                "duration_ns": ended - started,
                "stdout": artifact_record(output, stdout_path),
                "stderr": artifact_record(output, stderr_path),
                "result": result_record,
            }
            isolation = isolation_record(
                options.cpu, options.sibling, pair_lease_record,
                before_cpu, after_cpu, before_sibling, after_sibling,
                started, ended)
            pair_guard.validate_current()
            reservation_guard.validate_current()
            require(not timed_out, "C7 child timed out")
            require(completed.returncode == 0, f"C7 child exited {completed.returncode}")
            require(result_record is not None, "C7 child did not write result JSON")
            result_bytes = read_artifact(output, result_record, "result")
            parsed = strict_json(result_bytes, "C7 result")
            validate_isolation(isolation, options.cpu, options.sibling)
            validate_reservation_record(
                reservation_record, options.cpu, options.sibling)
            inputs_after = input_snapshot(
                source_root, options.expected_tooling_commit,
                options.expected_core_commit, archive, executable, taskset,
                build_manifest)
            require(typed_equal(inputs_after, inputs_before),
                    "source/archive/executable changed during C7 run")
            host_after = host_identity(options.cpu, options.sibling, allowed)
            require(typed_equal(host_after, host_before),
                    "host topology/frequency policy changed during C7 run")
            normalized = validate_child_result(parsed, options.cpu, inputs_before)
            created = utc_now()
            raw = signed({
                "schema": RAW_SCHEMA, "created_utc": created, "request": request,
                "inputs_before": inputs_before, "inputs_after": inputs_after,
                "host_before": host_before, "host_after": host_after,
                "reservation": reservation_record, "isolation": isolation,
                "build_provenance": build_provenance,
                "child": child_record, "validated_output": normalized,
            })
            validate_raw(raw, output, check_files=True)
            raw_path = output / "raw.json"
            write_json_exclusive(raw_path, raw)
            manifest = signed({
                "schema": MANIFEST_SCHEMA, "created_utc": created, "valid": True,
                "raw": artifact_record(output, raw_path), "request": request,
                "inputs": inputs_before, "reservation": reservation_record,
                "isolation": isolation, "build_provenance": build_provenance,
                "validated_output": normalized,
            })
            validate_manifest(manifest, output)
            reservation_guard.validate_current()
            pair_guard.validate_current()
            write_json_exclusive(output / "manifest.json", manifest)
    except Exception as error:
        error_type = type(error).__name__
        error_text = str(error)
        failure = signed({
            "schema": FAILURE_SCHEMA, "status": "failed",
            "failure_kind": classify_failure(
                child_record, error_type, error_text),
            "created_utc": utc_now(), "error_type": error_type,
            "error": error_text, "traceback": traceback.format_exc(),
            "request": request, "inputs_before": inputs_before,
            "inputs_after": inputs_after, "host_before": host_before,
            "host_after": host_after, "reservation": reservation_record,
            "pair_lease": pair_lease_record, "isolation": isolation,
            "build_provenance": build_provenance,
            "child": child_record,
        })
        validate_failure(failure, output, check_files=True)
        failure_path = output / "failure.json"
        if not failure_path.exists():
            write_json_exclusive(failure_path, failure)
        raise
    finally:
        if original_affinity is not None:
            os.sched_setaffinity(0, original_affinity)
    print(output / "manifest.json")
    return 0


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    root = manifest_path.parent
    manifest_bytes = manifest_path.read_bytes()
    manifest = strict_json(manifest_bytes, "C7 manifest")
    require(manifest_bytes == canonical_bytes(manifest) + b"\n",
            "manifest is not canonical JSON")
    normalized = validate_manifest(manifest, root)
    print(json.dumps({"status": "PASS", "cells": normalized["cell_count"],
                      "manifest_sha256": sha256_file(manifest_path)}, sort_keys=True))
    return 0


def verify_failure_campaign(options: argparse.Namespace) -> int:
    failure_path = options.failure.resolve(strict=True)
    root = failure_path.parent
    failure_bytes = failure_path.read_bytes()
    failure = strict_json(failure_bytes, "C7 failure evidence")
    require(failure_bytes == canonical_bytes(failure) + b"\n",
            "failure evidence is not canonical JSON")
    validated = validate_failure(failure, root, check_files=True)
    print(json.dumps({
        "status": "VERIFIED_FAILURE",
        "failure_kind": validated["failure_kind"],
        "failure_sha256": sha256_file(failure_path),
    }, sort_keys=True))
    return 0


def synthetic_result(cpu: int, inputs: Mapping[str, Any]) -> dict[str, Any]:
    cells = []
    for expected in EXPECTED_CELLS:
        k, r, byte_count, batch, losses, exact_field, padded_field = expected
        cell: dict[str, Any] = {
            "K": k, "R": r, "bytes": byte_count, "batch": batch,
            "losses": losses, "exact_field": exact_field,
            "padded_field": padded_field, "exact_coefficients": k * r,
            "exact_decode_terms": k * losses, "padded_encode_scratch": 64,
            "padded_decode_scratch": 128,
        }
        for index, (summary, samples) in enumerate(SUMMARY_SAMPLE_PAIRS, 1):
            values = [float(index)] * SAMPLE_COUNT
            cell[samples] = values
            cell[summary] = {"median_us": float(index), "mad_us": 0.0}
        cells.append(cell)
    return {"schema": CHILD_SCHEMA, "status": "pass",
            "profile": copy.deepcopy(EXPECTED_PROFILE),
            "production_constructor_rejected": True,
            "timing_scope": "candidate-authoritative-requires-runner-manifest",
            "requested_backend": "avx2", "runtime_backend": "avx2",
            "affinity": [cpu], "omp_num_threads": "1", "omp_dynamic": "FALSE",
            "source_sha256": inputs["source"]["sha256"],
            "core_git_sha": inputs["git"]["core_commit"],
            "library_sha256": inputs["archive"]["sha256"], "sanitizer": "none",
            "sanitizer_features": {"address": False, "undefined": False},
            "allocation_tracking": "global-new",
            "correctness": copy.deepcopy(EXPECTED_CORRECTNESS), "benchmarks": cells}


def synthetic_build_manifest(inputs: Mapping[str, Any]) -> dict[str, Any]:
    def artifact(key: str, path: Path) -> dict[str, Any]:
        return {"bytes": inputs[key]["size"], "path": path.as_posix(),
                "sha256": inputs[key]["sha256"]}

    build_dir = Path(".research/fixture/build/core-avx2")
    library = {"bytes": inputs["archive"]["size"],
               "path": (build_dir / "liblibleopard.a").as_posix(),
               "sha256": inputs["archive"]["sha256"]}
    executable = {"bytes": inputs["executable"]["size"],
                  "path": (build_dir.parent / "c7-avx2").as_posix(),
                  "sha256": inputs["executable"]["sha256"]}
    builds: list[dict[str, Any]] = []
    for name in BUILD_NAMES:
        if name == "avx2":
            builds.append({
                "name": "avx2", "backend": "avx2", "sanitizer": False,
                "build_dir": build_dir.as_posix(), "library": library,
                "executable": executable,
                "compile_argv": [
                    "/usr/bin/g++", "-O2", SOURCE_RELATIVE.as_posix(),
                    library["path"], "-pthread", "-fopenmp", "-o",
                    executable["path"],
                ],
                "instrumentation": {
                    "required_compile_macro": False,
                    "executable_counts": {"asan_lines": 0, "ubsan_lines": 0},
                    "core_archive_counts": {"asan_lines": 0, "ubsan_lines": 0},
                },
            })
        else:
            builds.append({"name": name})
    return {
        "schema": BUILD_MANIFEST_SCHEMA, "status": "pass", "scope": BUILD_SCOPE,
        "tooling_git_sha": inputs["git"]["tooling_commit"],
        "core_git_sha": inputs["git"]["core_commit"],
        "source": artifact("source", SOURCE_RELATIVE),
        "runner": artifact("build_runner", BUILD_RUNNER_RELATIVE),
        "validator": artifact("build_validator", BUILD_VALIDATOR_RELATIVE),
        "builds": builds,
        "reproducibility": {
            "comparison": {"status": "pass"},
            "fingerprints": {"avx2": {
                "library_sha256": library["sha256"],
                "executable_sha256": executable["sha256"],
            }},
        },
    }


def synthetic_inputs() -> dict[str, Any]:
    def identity(kind: str, name: str, mode: int = 0o644) -> dict[str, Any]:
        return {"kind": kind, "name": name, "size": 1,
                "mode": mode, "sha256": "a" * 64}

    result = {
        "git": {"tooling_commit": "b" * 40, "tooling_tree": "c" * 40,
                "core_commit": "e" * 40, "core_tree": "f" * 40,
                "core_is_ancestor": True, "clean": True,
                "tracked_tree_sha256": "d" * 64},
        "runner": {**identity("runner", RUNNER_RELATIVE.name, 0o755),
                   "path": RUNNER_RELATIVE.as_posix()},
        "source": {**identity("source", SOURCE_RELATIVE.name),
                   "path": SOURCE_RELATIVE.as_posix()},
        "build_runner": {**identity("build-runner", BUILD_RUNNER_RELATIVE.name),
                         "path": BUILD_RUNNER_RELATIVE.as_posix()},
        "build_validator": {
            **identity("build-validator", BUILD_VALIDATOR_RELATIVE.name),
            "path": BUILD_VALIDATOR_RELATIVE.as_posix()},
        "archive": identity("archive", "liblibleopard.a"),
        "executable": identity("executable", "c7-avx2", 0o755),
        "taskset": identity("taskset", "taskset", 0o755),
        "python": identity("python", "python3", 0o755),
        "core_source_closure": [
            {**identity("core-source", Path(path).name), "path": path}
            for path in CORE_SOURCE_CLOSURE
        ],
    }
    build_data = canonical_pretty_bytes(synthetic_build_manifest(result))
    result["build_manifest"] = {
        "kind": "build-manifest", "name": "build-run-manifest.json",
        "size": len(build_data), "mode": 0o644, "sha256": sha256_bytes(build_data),
    }
    result["build_attestation"] = derive_build_attestation(
        synthetic_build_manifest(result), result["build_manifest"],
        result["git"]["tooling_commit"], result["git"]["core_commit"],
        result["source"], result["build_runner"], result["build_validator"],
        result["archive"], result["executable"])
    result["binding_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def synthetic_cpu_stat(cpu: int, *, user: int, idle: int) -> dict[str, Any]:
    fields = {"user": user, "nice": 0, "system": 0, "idle": idle,
              "iowait": 0, "irq": 0, "softirq": 0, "steal": 0}
    return {"cpu": cpu, "fields": fields, "total_jiffies": sum(fields.values())}


def synthetic_pair_lease(cpu: int, sibling: int) -> dict[str, Any]:
    payload = pair_lease_payload(cpu, sibling, uid=1000)
    return {"device": 1, "directory_device": 1, "directory_inode": 2,
            "inode": 3, "lock": "exclusive_nonblocking_pair_wide",
            "path": str(pair_lease_directory(1000) /
                        pair_lease_name(cpu, sibling, 1000)),
            "payload": payload, "sha256": sha256_bytes(canonical_bytes(payload))}


def synthetic_host(cpu: int, sibling: int) -> dict[str, Any]:
    topology = {
        "core_id": "0", "physical_package_id": "0", "die_id": "0",
        "thread_siblings_list": f"{cpu},{sibling}",
        "core_siblings_list": f"{cpu},{sibling},2",
    }
    frequency = {
        "scaling_driver": "fixture", "scaling_governor": "performance",
        "energy_performance_preference": "performance", "scaling_min_freq": "1",
        "scaling_max_freq": "2", "cpuinfo_min_freq": "1",
        "cpuinfo_max_freq": "2",
    }
    def record(index: int) -> dict[str, Any]:
        return {"cpu": index, "online": "1", "topology": copy.deepcopy(topology),
                "frequency_policy": copy.deepcopy(frequency)}
    return {
        "system": {"system": "Linux", "release": "fixture",
                   "version": "fixture", "machine": "x86_64", "python": "3.11"},
        "allowed_at_launch": sorted((cpu, sibling, 2)),
        "timing_cpu": record(cpu), "sibling_cpu": record(sibling),
        "turbo": {
            "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
            "/sys/devices/system/cpu/amd_pstate/status": None,
            "/sys/devices/system/cpu/cpufreq/boost": None,
        },
    }


def synthetic_bundle(root: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cpu, sibling = 0, 1
    inputs = synthetic_inputs()
    result = synthetic_result(cpu, inputs)
    result_path = root / "child/result.json"
    stdout_path = root / "child/stdout.bin"
    stderr_path = root / "child/stderr.bin"
    write_json_exclusive(result_path, result)
    write_exclusive(stdout_path, b"")
    write_exclusive(stderr_path, expected_stderr())
    retained_build_manifest = root / "provenance/build-run-manifest-v4.json"
    write_exclusive(retained_build_manifest,
                    canonical_pretty_bytes(synthetic_build_manifest(inputs)))
    build_provenance = {
        "manifest": artifact_record(root, retained_build_manifest),
        "attestation": copy.deepcopy(inputs["build_attestation"]),
    }
    payload = {"benchmark_cpu": cpu, "nonce": "fixture-nonce", "owner": "self-test",
               "reserved_sibling": sibling, "schema": RESERVATION_SCHEMA,
               "status": "held"}
    reservation_bytes = canonical_bytes(payload)
    reservation = {"schema": RESERVATION_SCHEMA, "bytes": len(reservation_bytes),
                   "sha256": sha256_bytes(reservation_bytes), "payload": payload,
                   "lock": "exclusive_nonblocking"}
    isolation = isolation_record(
        cpu, sibling, synthetic_pair_lease(cpu, sibling),
        synthetic_cpu_stat(cpu, user=0, idle=0),
        synthetic_cpu_stat(cpu, user=1, idle=0),
        synthetic_cpu_stat(sibling, user=0, idle=0),
        synthetic_cpu_stat(sibling, user=0, idle=1), 1, 2)
    request = {
        "backend": "avx2", "cell_geometry": [list(cell) for cell in EXPECTED_CELLS],
        "child_environment": copy.deepcopy(CHILD_ENVIRONMENT), "cpu": cpu,
        "sibling": sibling, "timeout_seconds": 10.0,
        "command": ["${TASKSET}", "-c", "0", "${C7_EXECUTABLE}",
                    "--backend", "avx2", "${RESULT_JSON}"],
    }
    normalized = validate_child_result(result, cpu, inputs)
    created = "2026-07-17T00:00:00Z"
    raw = signed({
        "schema": RAW_SCHEMA, "created_utc": created, "request": request,
        "inputs_before": inputs, "inputs_after": copy.deepcopy(inputs),
        "host_before": synthetic_host(cpu, sibling),
        "host_after": synthetic_host(cpu, sibling), "reservation": reservation,
        "isolation": isolation, "build_provenance": build_provenance,
        "child": {"command": request["command"],
                  "environment": copy.deepcopy(CHILD_ENVIRONMENT),
                  "returncode": 0, "timed_out": False, "duration_ns": 1,
                  "stdout": artifact_record(root, stdout_path),
                  "stderr": artifact_record(root, stderr_path),
                  "result": artifact_record(root, result_path)},
        "validated_output": normalized,
    })
    raw_path = root / "raw.json"
    write_json_exclusive(raw_path, raw)
    manifest = signed({
        "schema": MANIFEST_SCHEMA, "created_utc": created, "valid": True,
        "raw": artifact_record(root, raw_path), "request": request,
        "inputs": inputs, "reservation": reservation, "isolation": isolation,
        "build_provenance": build_provenance, "validated_output": normalized,
    })
    return manifest, raw, result


def self_test() -> int:
    inputs = synthetic_inputs()
    validate_input_snapshot(inputs)
    result = synthetic_result(0, inputs)
    validate_child_result(result, 0, inputs)
    rejected = 0

    def expect_rejected(function: Any, label: str) -> None:
        nonlocal rejected
        try:
            function()
        except EvidenceError:
            rejected += 1
        else:
            raise EvidenceError(f"{label} mutation was accepted")

    for label, mutate in (
        ("summary", lambda value: value["benchmarks"][0]["exact_encode"].update(
            median_us=99.0)),
        ("sample count", lambda value: value["benchmarks"][0][
            "exact_encode_samples_us"].pop()),
        ("cell order", lambda value: value["benchmarks"].reverse()),
        ("backend", lambda value: value.update(runtime_backend="scalar")),
        ("embedded source", lambda value: value.update(source_sha256="f" * 64)),
        ("Boolean scratch", lambda value: value["benchmarks"][0].update(
            padded_encode_scratch=False)),
    ):
        changed = copy.deepcopy(result)
        mutate(changed)
        expect_rejected(lambda changed=changed: validate_child_result(
            changed, 0, inputs), label)

    stale = {"benchmark_cpu": 0, "nonce": "fixture-nonce", "owner": "self-test",
             "reserved_sibling": 1, "schema": RESERVATION_SCHEMA,
             "status": "released"}
    expect_rejected(lambda: parse_reservation(canonical_bytes(stale), 0, 1),
                    "stale reservation")

    busy = isolation_record(
        0, 1, synthetic_pair_lease(0, 1),
        synthetic_cpu_stat(0, user=0, idle=0),
        synthetic_cpu_stat(0, user=1, idle=0),
        synthetic_cpu_stat(1, user=0, idle=0),
        synthetic_cpu_stat(1, user=1, idle=1), 1, 2)
    expect_rejected(lambda: validate_isolation(busy, 0, 1), "busy sibling")

    with tempfile.TemporaryDirectory(prefix="leopard2-c7-authoritative-") as directory:
        root = Path(directory)
        manifest, raw, _fixture_result = synthetic_bundle(root)
        validate_manifest(manifest, root)
        changed_raw = copy.deepcopy(raw)
        changed_raw["request"]["cell_geometry"].reverse()
        changed_raw = signed({key: value for key, value in changed_raw.items()
                              if key != "digest"})
        expect_rejected(lambda: validate_raw(changed_raw, root), "raw schedule")
        changed_raw = copy.deepcopy(raw)
        changed_raw["host_after"]["timing_cpu"]["frequency_policy"][
            "scaling_governor"] = "powersave"
        changed_raw = signed({key: value for key, value in changed_raw.items()
                              if key != "digest"})
        expect_rejected(lambda: validate_raw(changed_raw, root), "host drift")
        result_path = root / "child/result.json"
        original = result_path.read_bytes()
        result_path.write_bytes(b"{}\n")
        expect_rejected(lambda: validate_manifest(manifest, root), "artifact bytes")
        result_path.write_bytes(original)
        changed_manifest = copy.deepcopy(manifest)
        changed_manifest["validated_output"]["cell_count"] = 11
        changed_manifest = signed({key: value for key, value in changed_manifest.items()
                                   if key != "digest"})
        expect_rejected(lambda: validate_manifest(
            changed_manifest, root), "manifest summary")

    print(json.dumps({"status": "PASS", "cells": len(EXPECTED_CELLS),
                      "samples_per_metric": SAMPLE_COUNT,
                      "mutations_rejected": rejected}, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="execute one fixed 12-cell authoritative AVX2 campaign")
    run.add_argument("--source-root", required=True, type=Path,
                     help="clean checkout containing the committed C7 tools")
    run.add_argument("--expected-tooling-commit", required=True,
                     help="full SHA-1 that must equal clean checkout HEAD")
    run.add_argument("--expected-core-commit", required=True,
                     help="full ancestor SHA-1 embedded in the C7 executable")
    run.add_argument("--archive", required=True, type=Path,
                     help="exact AVX2 liblibleopard.a linked into the harness")
    run.add_argument("--executable", required=True, type=Path,
                     help="exact non-sanitized AVX2 c7_exact_low executable")
    run.add_argument("--build-manifest", required=True, type=Path,
                     help=("final v4 A/B build manifest validated live and bound "
                           "to the supplied AVX2 archive/executable"))
    run.add_argument("--reservation-file", required=True, type=Path,
                     help="canonical held v1 CPU-pair reservation to lock")
    run.add_argument("--output", required=True, type=Path,
                     help="new directory for raw, manifest, logs, or failure")
    run.add_argument("--cpu", required=True, type=int,
                     help="allowed logical CPU used for all timing")
    run.add_argument("--sibling", required=True, type=int,
                     help="the timing CPU's allowed SMT sibling to reserve")
    run.add_argument("--timeout", type=float, default=1800.0,
                     help="positive child timeout in seconds (default: 1800)")
    run.set_defaults(function=run_campaign)
    verify = commands.add_parser(
        "verify", help="portably replay a retained canonical evidence bundle")
    verify.add_argument("--manifest", required=True, type=Path,
                        help="retained manifest.json")
    verify.set_defaults(function=verify_campaign)
    verify_failure = commands.add_parser(
        "verify-failure",
        help="portably authenticate a failed campaign without treating it as success")
    verify_failure.add_argument("--failure", required=True, type=Path,
                                help="retained failure.json")
    verify_failure.set_defaults(function=verify_failure_campaign)
    test = commands.add_parser(
        "self-test", help="run synthetic validation and mutation tests only")
    test.set_defaults(function=lambda unused: self_test())
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    try:
        return int(options.function(options))
    except (EvidenceError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"C7 authoritative evidence error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
