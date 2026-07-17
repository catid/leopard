#!/usr/bin/env python3
"""Portable and optional live validation for C7 v3/v4 manifests."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import pathlib
import posixpath
import re
import shlex
import shutil
import statistics
import subprocess
import tarfile
import tempfile
import zlib
from typing import Any, Iterable

from run_matrix import (
    BACKENDS,
    BUILD_NAMES,
    C7_BINARY_MAX_BYTES,
    EXPECTED_ARCHIVE_MEMBER_COUNTS,
    EXPECTED_ARCHIVE_SANITIZER_COUNTS,
    EXPECTED_EXECUTABLE_SANITIZER_COUNTS,
    EXPECTED_SOURCE_CLOSURE,
    EXPECTED_TOOLING_CLOSURE,
    LEGACY_EXPECTED_ARCHIVE_MEMBER_COUNTS,
    LEGACY_EXPECTED_ARCHIVE_SANITIZER_COUNTS,
    LEGACY_EXPECTED_EXECUTABLE_SANITIZER_COUNTS,
    LEGACY_MANIFEST_SCHEMA,
    LEGACY_PEER_ATTESTATION_SCHEMA,
    MANIFEST_SCHEMA,
    NORMALIZATION_SCHEMA,
    NORMALIZATION_TOKEN,
    PEER_ATTESTATION_SCHEMA,
    PEER_BUNDLE_MAX_ARCHIVE_BYTES,
    PEER_BUNDLE_MAX_MEMBER_BYTES,
    PEER_BUNDLE_MAX_MEMBERS,
    PEER_BUNDLE_MAX_TOTAL_BYTES,
    PEER_BUNDLE_SCHEMA,
    PREFIX_MAP_TARGET,
    PROGRAM_ROLES,
    canonical_peer_tar,
    peer_portable_artifact_records,
    typed_equal,
)


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")
CORE_SOURCES = {
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF16.cpp",
    "LeopardFF8.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
}
ARCHIVE_MEMBERS = tuple(EXPECTED_ARCHIVE_MEMBER_COUNTS)
EXPECTED_CORRECTNESS = {
    "gf8_cases": 9,
    "gf16_cases": 5,
    "coefficients": 118717,
    "gf16_vandermonde_coefficients": 1500,
    "encode_executions": 117,
    "encode_symbol_comparisons": 1030423,
    "subset_encode_executions": 117,
    "decode_plans": 403,
    "decode_executions": 403,
    "decode_symbol_comparisons": 272487,
    "maximum_loss_plans": 117,
    "unavailable_parity_plans": 175,
    "no_loss_null_calls": 14,
    "parity_rebuilds": 403,
    "odd_gf16_rejections": 10,
    "overlap_rejections": 59,
    "parity_output_overlap_rejections": 13,
    "restored_output_overlap_rejections": 12,
    "restored_input_overlap_rejections": 20,
    "selected_parity_null_rejections": 14,
    "survivor_null_rejections": 6,
    "atomic_rejection_bytes_checked": 61570,
    "read_only_input_alias_calls": 13,
    "read_only_input_alias_symbol_comparisons": 2139,
    "decode_read_only_input_alias_calls": 117,
    "decode_read_only_input_alias_symbol_comparisons": 6025,
    "detached_plan_executions": 14,
    "detached_plan_symbol_comparisons": 3598,
    "owned_backend_codec_executions": 2,
    "detached_backend_plan_executions": 2,
    "owned_backend_symbol_comparisons": 196,
    "concurrent_backend_contexts": 3,
    "concurrent_backend_executions": 384,
    "concurrent_backend_trace_calls": 30720,
    "concurrent_wire_digest_comparisons": 4,
    "exhaustive_small_plans": 163,
    "exhaustive_small_executions": 4720,
    "exhaustive_small_symbol_comparisons": 8192,
    "malformed_plan_rejections": 102,
    "hot_path_allocations": 0,
    "digest_fnv64": "0x329cac84bfd9f27",
}
LEGACY_EXPECTED_CORRECTNESS = dict(EXPECTED_CORRECTNESS)
for _new_correctness_key in (
    "detached_plan_executions", "detached_plan_symbol_comparisons",
    "owned_backend_codec_executions", "detached_backend_plan_executions",
    "owned_backend_symbol_comparisons",
    "concurrent_backend_contexts", "concurrent_backend_executions",
    "concurrent_backend_trace_calls", "concurrent_wire_digest_comparisons",
    "exhaustive_small_plans", "exhaustive_small_executions",
    "exhaustive_small_symbol_comparisons", "malformed_plan_rejections",
):
    del LEGACY_EXPECTED_CORRECTNESS[_new_correctness_key]
LEGACY_EXPECTED_CORRECTNESS["digest_fnv64"] = "0xec4179e9f2776a58"
EXPECTED_TRACED_BACKENDS = {
    "scalar": 1, "ssse3": 2, "avx2": 3, "auto": 3,
}
EXPECTED_PROFILE = {
    "family": 3, "version": 1, "coordinate_map": 1,
    "systematic": "0..K-1", "parity": "K..K+R-1",
    "production_enabled": False,
}
EXPECTED_RUNTIME = {
    "scalar": "scalar", "ssse3": "ssse3", "avx2": "avx2", "auto": "avx2",
}
ABSOLUTE_PROJECT_PATH = re.compile(
    r"(?<![A-Za-z0-9_.$}{-])(?:[A-Za-z]:[\\/]|/)"
    r"(?:[A-Za-z0-9_.+@~-]+[\\/])+(?:"
    r"CMakeLists\.txt|cmake[\\/]leopardConfig\.cmake\.in|"
    r"experiments[\\/]leopard2[\\/]|leopard2?\.(?:cpp|h)|"
    r"Leopard(?:2|Common|FF)[A-Za-z0-9_]*\.(?:cpp|h))")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def canonical_pretty_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def strict_json_loads(document: str | bytes, label: str) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label} contains a duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> Any:
        raise ValueError(f"{label} contains a non-standard JSON value: {value}")

    def finite_float(value: str) -> float:
        parsed = float(value)
        if not math.isfinite(parsed):
            raise ValueError(f"{label} contains a non-finite JSON number")
        return parsed

    try:
        return json.loads(
            document, object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant, parse_float=finite_float)
    except (json.JSONDecodeError, RecursionError) as error:
        raise ValueError(f"{label} JSON is invalid") from error


def validate_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        raise ValueError(f"{label} is not canonical SHA-256")
    return value


def _walk_dicts(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_dicts(child)


def validate_comparison(
    comparison: object, fingerprints: dict[str, Any],
    normalized_text_records: int,
) -> None:
    if type(normalized_text_records) is not int or normalized_text_records < 0:
        raise ValueError("normalized text record count is invalid")
    if typed_equal(comparison, {"status": "not-run"}):
        return
    required = {
        "build_names", "checkout_roots_scanned", "current_scan",
        "fingerprints_sha256", "peer_attestation", "peer_manifest_sha256",
        "peer_evidence_bundle", "peer_manifest", "peer_scan", "status",
    }
    if not isinstance(comparison, dict) or set(comparison) != required:
        raise ValueError("A/B comparison attestation schema changed")
    if (comparison["status"] != "pass" or
            not typed_equal(comparison["build_names"], list(BUILD_NAMES)) or
            type(comparison["checkout_roots_scanned"]) is not int or
            comparison["checkout_roots_scanned"] != 2):
        raise ValueError("A/B comparison attestation changed")
    validate_sha(comparison["peer_manifest_sha256"],
                 "A/B peer manifest hash")
    for key, filename in (
        ("peer_manifest", "peer-manifest.json"),
        ("peer_evidence_bundle", "peer-evidence-bundle.tar.gz"),
        ("peer_attestation", "peer-reproducibility-attestation.json"),
    ):
        record = comparison[key]
        if (not isinstance(record, dict) or
                pathlib.PurePosixPath(str(record.get("path", ""))).name !=
                filename):
            raise ValueError(f"A/B {key} artifact identity changed")
    if comparison["peer_manifest"]["sha256"] != comparison[
            "peer_manifest_sha256"]:
        raise ValueError("A/B peer manifest hashes disagree")
    if comparison["fingerprints_sha256"] != canonical_json_sha256(fingerprints):
        raise ValueError("A/B fingerprint attestation changed")
    expected_scan = {
        "normalized_text_records": normalized_text_records,
        "archives": len(BUILD_NAMES),
        "executables": len(BUILD_NAMES),
    }
    if (not typed_equal(comparison["current_scan"], expected_scan) or
            not typed_equal(comparison["peer_scan"], expected_scan)):
        raise ValueError("A/B root-byte scan attestation changed")


def manifest_program_records(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "taskset": manifest["taskset"],
        "builds": {
            build["name"]: {role: build[role] for role in PROGRAM_ROLES}
            for build in manifest["builds"]
        },
    }


def manifest_binary_records(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        build["name"]: {
            key: build[key] for key in ("library", "executable")
        }
        for build in manifest["builds"]
    }


def _canonical_relative_path(path: object, label: str) -> str:
    if not isinstance(path, str) or not path or "\\" in path or ":" in path:
        raise ValueError(f"{label} path is not canonical")
    pure = pathlib.PurePosixPath(path)
    if (pure.is_absolute() or pure.as_posix() != path or
            any(part in ("", ".", "..") for part in pure.parts)):
        raise ValueError(f"{label} path is not canonical")
    return path


def _read_file_once(path: pathlib.Path, maximum_bytes: int) -> bytes:
    """Read at most one bounded snapshot from an already-contained path."""
    with path.open("rb") as stream:
        return stream.read(maximum_bytes + 1)


def read_verified_artifact_bytes(
    record: object, label: str, evidence_root: pathlib.Path, *,
    maximum_bytes: int, minimum_bytes: int = 1,
) -> bytes:
    """Read once, then authenticate and return that exact byte string."""
    if not isinstance(record, dict) or set(record) != {
            "bytes", "path", "sha256"}:
        raise ValueError(f"{label} artifact schema changed")
    size = record["bytes"]
    if (type(size) is not int or type(minimum_bytes) is not int or
            not 0 <= minimum_bytes <= size <= maximum_bytes):
        raise ValueError(f"{label} byte count is invalid")
    expected = validate_sha(record["sha256"], f"{label} hash")
    try:
        path = resolve_path(record["path"], evidence_root)
        contents = _read_file_once(path, maximum_bytes)
    except (OSError, ValueError) as error:
        raise ValueError(f"{label} retained artifact cannot be read") from error
    if (len(contents) != size or
            hashlib.sha256(contents).hexdigest() != expected):
        raise ValueError(f"{label} bytes disagree with the manifest")
    return contents


class VerifiedTextArtifact:
    """Path identity plus immutable, authenticated retained text bytes."""

    def __init__(self, path: pathlib.Path, contents: bytes, text: str) -> None:
        self.path = path
        self.contents = contents
        self.text = text

    @property
    def name(self) -> str:
        return self.path.name

    def read_bytes(self) -> bytes:
        return self.contents

    def read_text(self, encoding: str = "utf-8") -> str:
        if encoding.lower().replace("_", "-") != "utf-8":
            raise ValueError("verified text is available only as UTF-8")
        return self.text


class VerifiedBinaryArtifact:
    """Path identity plus immutable, authenticated binary bytes."""

    def __init__(self, path: pathlib.Path, contents: bytes) -> None:
        self.path = path
        self.contents = contents

    @property
    def name(self) -> str:
        return self.path.name


def _peer_bundle_records(peer: dict[str, Any]) -> dict[str, dict[str, Any]]:
    try:
        records = peer_portable_artifact_records(peer)
    except (KeyError, TypeError, RuntimeError) as error:
        raise ValueError("peer portable artifact graph is incomplete") from error
    expected: dict[str, dict[str, Any]] = {}
    total = 0
    for index, record in enumerate(records):
        if not isinstance(record, dict) or set(record) != {
                "bytes", "path", "sha256", "source_root_tokens"}:
            raise ValueError("peer portable artifact record schema changed")
        path = _canonical_relative_path(
            record["path"], f"peer portable artifact {index}")
        size = record["bytes"]
        if (path in expected or type(size) is not int or size < 0 or
                size > PEER_BUNDLE_MAX_MEMBER_BYTES):
            raise ValueError("peer portable artifact size/path is invalid")
        validate_sha(record["sha256"], f"peer portable artifact {index}")
        expected[path] = record
        total += size
    if (len(expected) > PEER_BUNDLE_MAX_MEMBERS or
            total > PEER_BUNDLE_MAX_TOTAL_BYTES):
        raise ValueError("peer portable evidence exceeds bundle limits")
    return expected


def _read_peer_bundle(
    contents: bytes, expected: dict[str, dict[str, Any]],
) -> dict[str, bytes]:
    if not contents or len(contents) > PEER_BUNDLE_MAX_ARCHIVE_BYTES:
        raise ValueError("peer bundle compressed size is invalid")
    if contents[:10] != b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff":
        raise ValueError("peer bundle gzip metadata is noncanonical")
    # Bound expansion before tar parsing.  Canonical tar padding is small, but
    # leave one MiB for the index and block records while retaining a hard
    # global ceiling for malicious manifests.
    expected_payload = sum(record["bytes"] for record in expected.values())
    expansion_limit = min(
        PEER_BUNDLE_MAX_ARCHIVE_BYTES,
        expected_payload + (len(expected) + 8) * 2048 + 1024 * 1024)
    try:
        decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        raw = decompressor.decompress(contents, expansion_limit + 1)
    except zlib.error as error:
        raise ValueError("peer bundle gzip stream is invalid") from error
    if len(raw) > expansion_limit:
        raise ValueError("peer bundle exceeds its decompression limit")
    if (not decompressor.eof or decompressor.unused_data or
            decompressor.unconsumed_tail):
        raise ValueError("peer bundle must contain exactly one gzip member")
    members: dict[str, bytes] = {}
    total = 0
    try:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as archive:
            infos = archive.getmembers()
            if len(infos) > PEER_BUNDLE_MAX_MEMBERS + 1:
                raise ValueError("peer bundle member count exceeds its limit")
            if [info.name for info in infos] != sorted(info.name for info in infos):
                raise ValueError("peer bundle member order is noncanonical")
            for info in infos:
                name = _canonical_relative_path(info.name, "peer bundle member")
                if (not info.isreg() or name in members or info.size < 0 or
                        info.size > PEER_BUNDLE_MAX_MEMBER_BYTES or
                        info.mode != 0o444 or info.mtime != 0 or
                        info.uid != 0 or info.gid != 0 or
                        info.uname != "" or info.gname != "" or
                        info.pax_headers):
                    raise ValueError("peer bundle contains an unsafe member")
                total += info.size
                if total > PEER_BUNDLE_MAX_TOTAL_BYTES + 1024 * 1024:
                    raise ValueError("peer bundle payload exceeds its limit")
                extracted = archive.extractfile(info)
                if extracted is None:
                    raise ValueError("peer bundle member cannot be read")
                payload = extracted.read(info.size + 1)
                if len(payload) != info.size:
                    raise ValueError("peer bundle member size changed")
                members[name] = payload
    except (tarfile.TarError, OSError) as error:
        raise ValueError("peer bundle tar stream is invalid") from error

    expected_names = {"index.json", *{
        f"files/{path}" for path in expected
    }}
    if set(members) != expected_names:
        raise ValueError("peer bundle member set changed")
    try:
        index = strict_json_loads(members["index.json"], "peer bundle index")
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError("peer bundle index is invalid") from error
    expected_index = {
        "schema": PEER_BUNDLE_SCHEMA,
        "files": [
            {
                "path": path, "bytes": record["bytes"],
                "sha256": record["sha256"],
            }
            for path, record in sorted(expected.items())
        ],
    }
    if not typed_equal(index, expected_index):
        raise ValueError("peer bundle index disagrees with peer manifest")
    files = {path: members[f"files/{path}"] for path in expected}
    for path, record in expected.items():
        payload = files[path]
        if (len(payload) != record["bytes"] or
                hashlib.sha256(payload).hexdigest() != record["sha256"]):
            raise ValueError("peer bundle bytes disagree with peer manifest")
    try:
        canonical_tar = canonical_peer_tar(files)
    except RuntimeError as error:
        raise ValueError("peer bundle cannot be canonically reconstructed") from error
    if raw != canonical_tar:
        raise ValueError("peer bundle tar encoding or metadata is noncanonical")
    return files


def _materialize_peer_bundle(
    files: dict[str, bytes], destination: pathlib.Path,
) -> None:
    destination.mkdir(mode=0o700)
    for path, contents in files.items():
        target = destination.joinpath(*pathlib.PurePosixPath(path).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(contents)
        target.chmod(0o444)


def validate_peer_attestation(
    comparison: dict[str, Any], manifest: dict[str, Any],
    source_root: pathlib.Path, evidence_root: pathlib.Path, *, live: bool,
) -> None:
    if typed_equal(comparison, {"status": "not-run"}):
        return
    peer_manifest_record = comparison["peer_manifest"]
    peer_manifest_bytes = read_verified_artifact_bytes(
        peer_manifest_record, "A/B peer manifest", evidence_root,
        maximum_bytes=16 * 1024 * 1024)
    try:
        peer = strict_json_loads(peer_manifest_bytes, "A/B peer manifest")
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError("A/B peer manifest JSON is invalid") from error
    if not isinstance(peer, dict):
        raise ValueError("A/B peer manifest is not an object")
    if peer_manifest_bytes != canonical_pretty_json_bytes(peer):
        raise ValueError("A/B peer manifest JSON is noncanonical")
    if not typed_equal(
            peer.get("reproducibility", {}).get("comparison"),
            {"status": "not-run"}):
        raise ValueError("A/B peer manifest is not an initial independent run")
    legacy = manifest.get("schema") == LEGACY_MANIFEST_SCHEMA
    identity_keys = ("core_git_sha",) if legacy else (
        "tooling_git_sha", "core_git_sha")
    try:
        identity_mismatch = (
            not typed_equal(peer.get("schema"), manifest.get("schema")) or
            any(not typed_equal(peer.get(key), manifest.get(key))
                for key in identity_keys) or
            any(not typed_equal(peer.get(key), manifest.get(key)) for key in (
                "source", "runner", "validator", "normalization")) or
            not typed_equal(peer.get("taskset"), manifest.get("taskset")) or
            not typed_equal(
                manifest_program_records(peer),
                manifest_program_records(manifest)) or
            not typed_equal(
                manifest_binary_records(peer),
                manifest_binary_records(manifest)) or
            not typed_equal(
                peer.get("reproducibility", {}).get("fingerprints"),
                manifest.get("reproducibility", {}).get("fingerprints")))
    except (KeyError, TypeError) as error:
        raise ValueError("A/B peer identity graph is incomplete") from error
    if identity_mismatch:
        raise ValueError("A/B peer does not bind the same source/program outputs")

    bundle_record = comparison["peer_evidence_bundle"]
    bundle_bytes = read_verified_artifact_bytes(
        bundle_record, "A/B peer evidence bundle", evidence_root,
        maximum_bytes=PEER_BUNDLE_MAX_ARCHIVE_BYTES)
    expected_bundle_records = _peer_bundle_records(peer)
    files = _read_peer_bundle(bundle_bytes, expected_bundle_records)
    with tempfile.TemporaryDirectory(prefix="c7-peer-replay-") as directory:
        peer_evidence_root = pathlib.Path(directory) / "evidence"
        _materialize_peer_bundle(files, peer_evidence_root)
        # The historical peer's exact tools and binary inode bytes were checked
        # during capture.  Retained replay is intentionally portable; --live
        # validates the current independent outputs, whose records are required
        # above to equal the peer exactly.
        validate_manifest(
            peer, source_root=source_root, evidence_root=peer_evidence_root,
            live=False, require_checkout_head=False)

    record = comparison["peer_attestation"]
    attestation_bytes = read_verified_artifact_bytes(
        record, "A/B peer attestation", evidence_root,
        maximum_bytes=1024 * 1024)
    if pathlib.PurePosixPath(record["path"]).name != (
            "peer-reproducibility-attestation.json"):
        raise ValueError("A/B peer attestation path changed")
    report = strict_json_loads(attestation_bytes, "A/B peer attestation")
    if attestation_bytes != canonical_pretty_json_bytes(report):
        raise ValueError("A/B peer attestation JSON is noncanonical")
    required = {
        "binary_artifacts", "checks", "core_git_sha", "fingerprints",
        "normalized_text_records_sha256", "peer_evidence_bundle",
        "peer_manifest_artifact",
        "program_records_sha256", "root_scan", "runs_sha256", "schema",
        "source_closures_sha256", "status", "tooling",
    }
    if not legacy:
        required.add("tooling_git_sha")
    if set(report) != required:
        raise ValueError("A/B peer attestation schema changed")
    builds = {build["name"]: build for build in manifest["builds"]}
    expected_binaries = manifest_binary_records(manifest)
    expected_closures = {
        name: builds[name]["source_closure"] for name in BUILD_NAMES
    }
    expected_attestation_schema = (
        LEGACY_PEER_ATTESTATION_SCHEMA if legacy else PEER_ATTESTATION_SCHEMA)
    expected_checks = {
        "portable_semantics": "pass",
        "live_tools_and_outputs": "pass",
        ("git_toplevel_and_head" if legacy else
         "clean_tooling_head_and_core_ancestor"): "pass",
        "program_identity_match": "pass",
        "binary_and_text_root_scan": "pass",
    }
    if (report["schema"] != expected_attestation_schema or
            report["status"] != "pass" or
            report["core_git_sha"] != manifest["core_git_sha"] or
            (not legacy and report["tooling_git_sha"] !=
             manifest["tooling_git_sha"]) or
            not typed_equal(report["tooling"], {
                key: manifest[key] for key in ("source", "runner", "validator")
            }) or not typed_equal(
                report["fingerprints"],
                manifest["reproducibility"]["fingerprints"]) or
            not typed_equal(report["binary_artifacts"], expected_binaries) or
            report["program_records_sha256"] != canonical_json_sha256(
                manifest_program_records(manifest)) or
            report["source_closures_sha256"] != canonical_json_sha256(
                expected_closures) or not typed_equal(
                    report["root_scan"], comparison["peer_scan"]) or
            not typed_equal(report["checks"], expected_checks) or
            not typed_equal(
                report["peer_manifest_artifact"], peer_manifest_record) or
            not typed_equal(report["peer_evidence_bundle"], bundle_record) or
            report["normalized_text_records_sha256"] != canonical_json_sha256(
                sorted(
                    (dict(item) for item in peer_portable_artifact_records(peer)),
                    key=lambda item: item["path"])) or
            report["runs_sha256"] != canonical_json_sha256(peer["runs"])):
        raise ValueError("A/B peer attestation contents changed")
    serialized = json.dumps(report, sort_keys=True)
    if (str(source_root.resolve()) in serialized or
            str(evidence_root.resolve()) in serialized or
            ABSOLUTE_PROJECT_PATH.search(serialized)):
        raise ValueError("A/B peer attestation leaked a checkout path")


def validate_git_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not GIT_SHA_RE.fullmatch(value):
        raise ValueError(f"{label} is not a canonical git SHA")
    return value


def resolve_path(path_text: str, source_root: pathlib.Path) -> pathlib.Path:
    if (not isinstance(path_text, str) or not path_text or
            "\\" in path_text or ":" in path_text):
        raise ValueError("artifact path is not canonical checkout-relative POSIX")
    pure = pathlib.PurePosixPath(path_text)
    if (pure.is_absolute() or pure.as_posix() != path_text or
            any(part in ("", ".", "..") for part in pure.parts)):
        raise ValueError("artifact path is not canonical checkout-relative POSIX")
    root = source_root.resolve()
    candidate = root.joinpath(*pure.parts)
    resolved = candidate.resolve(strict=False)
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ValueError("artifact path escapes the source checkout") from error
    return candidate


def validate_artifact(
    record: object, label: str, source_root: pathlib.Path, *, required: bool,
    check_if_present: bool = True,
) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "bytes", "path", "sha256"}:
        raise ValueError(f"{label} artifact schema changed")
    path_text = record["path"]
    if type(record["bytes"]) is not int or record["bytes"] < 0:
        raise ValueError(f"{label} byte count is invalid")
    expected = validate_sha(record["sha256"], f"{label} hash")
    try:
        path = resolve_path(path_text, source_root)
    except ValueError as error:
        raise ValueError(f"{label} path is invalid: {error}") from error
    if required and not path.is_file():
        raise ValueError(f"{label} retained artifact is missing")
    if path.is_file() and check_if_present:
        if path.stat().st_size != record["bytes"] or sha256(path) != expected:
            raise ValueError(f"{label} bytes disagree with the manifest")
    return path


def validate_normalized_text(
    record: object, label: str, artifact_root: pathlib.Path, *,
    require_token: bool, checkout_root: pathlib.Path | None = None,
) -> VerifiedTextArtifact:
    if not isinstance(record, dict) or set(record) != {
            "bytes", "path", "sha256", "source_root_tokens"}:
        raise ValueError(f"{label} normalized-artifact schema changed")
    token_count = record["source_root_tokens"]
    if type(token_count) is not int or token_count < 0:
        raise ValueError(f"{label} token count is invalid")
    generic = {key: record[key] for key in ("bytes", "path", "sha256")}
    contents = read_verified_artifact_bytes(
        generic, label, artifact_root,
        maximum_bytes=PEER_BUNDLE_MAX_MEMBER_BYTES, minimum_bytes=0)
    try:
        path = resolve_path(record["path"], artifact_root)
        text = contents.decode("utf-8", errors="strict")
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError(f"{label} retained text is invalid") from error
    if text.count(NORMALIZATION_TOKEN) != token_count:
        raise ValueError(f"{label} normalization token count changed")
    if require_token and token_count == 0:
        raise ValueError(f"{label} lacks its required normalization token")
    checkout_root = artifact_root if checkout_root is None else checkout_root
    if str(checkout_root.resolve()) in text or ABSOLUTE_PROJECT_PATH.search(text):
        raise ValueError(f"{label} leaked an absolute checkout path")
    return VerifiedTextArtifact(path, contents, text)


def validate_program_record(record: object, label: str) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "path", "sha256", "version"}:
        raise ValueError(f"{label} program schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not pathlib.Path(path_text).is_absolute():
        raise ValueError(f"{label} path is not absolute")
    validate_sha(record["sha256"], f"{label} hash")
    if not isinstance(record["version"], str) or not record["version"].strip() or (
            not record["version"].endswith("\n")):
        raise ValueError(f"{label} version output is not canonical")
    return pathlib.Path(path_text)


def validate_program_live(record: dict, label: str) -> pathlib.Path:
    path = validate_program_record(record, label)
    if not path.is_file() or sha256(path) != record["sha256"]:
        raise ValueError(f"{label} exact executable is unavailable or changed")
    completed = subprocess.run(
        [str(path), "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.stdout != record["version"]:
        raise ValueError(f"{label} exact version output changed")
    return path


def validate_argv(
    argv: object, token_count: object, label: str, *, require_token: bool,
    source_root: pathlib.Path = ROOT,
) -> list[str]:
    if not isinstance(argv, list) or not all(isinstance(item, str) for item in argv):
        raise ValueError(f"{label} argv is not an exact string array")
    if type(token_count) is not int or token_count < 0:
        raise ValueError(f"{label} argv token count is invalid")
    text = "\n".join(argv)
    if text.count(NORMALIZATION_TOKEN) != token_count:
        raise ValueError(f"{label} argv token count changed")
    if require_token and token_count == 0:
        raise ValueError(f"{label} argv lacks its normalization token")
    if str(source_root.resolve()) in text or ABSOLUTE_PROJECT_PATH.search(text):
        raise ValueError(f"{label} argv leaked an absolute checkout path")
    return argv


def validate_git_artifact(
    commit: str, record: dict, expected_relative: str, label: str,
    source_root: pathlib.Path,
) -> None:
    if record["path"] != expected_relative:
        raise ValueError(f"{label} repository path changed")
    completed = subprocess.run(
        ["git", "show", f"{commit}:{expected_relative}"], cwd=source_root,
        check=False, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if (completed.returncode != 0 or
            len(completed.stdout) != record["bytes"] or
            hashlib.sha256(completed.stdout).hexdigest() != record["sha256"]):
        raise ValueError(f"{label} is not bound to its declared commit")


def parse_cache(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text.splitlines():
        match = re.fullmatch(r"([^#/:=][^:=]*):([^=]+)=(.*)", line)
        if match:
            key, value = match.group(1), match.group(3)
            if key in values:
                raise ValueError(f"duplicate CMake cache key: {key}")
            values[key] = value
    return values


def normalize_checkout_text(
    text: str, source_root: pathlib.Path,
) -> str:
    pattern = re.compile(
        re.escape(str(source_root.resolve())) + r"(?=\Z|[/\\=\s\"'])")
    return pattern.sub(NORMALIZATION_TOKEN, text)


def dependency_source_closure(
    entries: object, build_dir: str, source_root: pathlib.Path,
    evidence_root: pathlib.Path, retained_directory: str, *, live: bool,
) -> tuple[str, ...]:
    if not isinstance(entries, list) or len(entries) != len(ARCHIVE_MEMBERS):
        raise ValueError("dependency-file matrix changed")
    closure = {"CMakeLists.txt", "cmake/leopardConfig.cmake.in"}
    build_paths: list[str] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or set(entry) != {"build_path", "retained"}:
            raise ValueError("dependency-file record schema changed")
        build_path = entry["build_path"]
        if not isinstance(build_path, str):
            raise ValueError("dependency build path is invalid")
        try:
            raw_path = resolve_path(build_path, evidence_root)
        except ValueError as error:
            raise ValueError("dependency build path escapes evidence root") from error
        expected_retained_path = (
            pathlib.PurePosixPath(retained_directory) /
            pathlib.PurePosixPath(build_path).name).as_posix()
        if (not isinstance(entry["retained"], dict) or
                entry["retained"].get("path") != expected_retained_path):
            raise ValueError("retained dependency path or build role changed")
        retained = validate_normalized_text(
            entry["retained"], f"dependency file {index}", evidence_root,
            require_token=False, checkout_root=source_root)
        if retained.name != pathlib.PurePosixPath(build_path).name:
            raise ValueError("retained dependency identity changed")
        retained_text = retained.read_text(encoding="utf-8")
        if live:
            if not raw_path.is_file():
                raise ValueError("live dependency file is missing")
            live_text = raw_path.read_text(encoding="utf-8", errors="strict")
            if normalize_checkout_text(live_text, source_root) != retained_text:
                raise ValueError("live dependency file differs from retained bytes")
        flattened = retained_text.replace("\\\n", " ")
        target, separator, dependency_text = flattened.partition(":")
        expected_target = build_path[len(build_dir) + 1:-2]
        if separator != ":" or target != expected_target:
            raise ValueError("retained dependency file is malformed")
        for token in shlex.split(dependency_text):
            if token.startswith(f"{NORMALIZATION_TOKEN}/"):
                candidate = token[len(NORMALIZATION_TOKEN) + 1:]
            elif token.startswith("/"):
                continue
            else:
                candidate = posixpath.normpath(
                    posixpath.join(build_dir, token))
            if (candidate == build_dir or
                    candidate.startswith(f"{build_dir}/")):
                continue
            if (candidate.startswith("../") or candidate == ".." or
                    candidate.startswith("/") or "\\" in candidate):
                raise ValueError("dependency source path is not canonical")
            closure.add(candidate)
        build_paths.append(build_path)
    expected_build_paths = sorted(
        f"{build_dir}/CMakeFiles/{target}/{member}.d"
        for member in ARCHIVE_MEMBERS
        for target in [
            "leopard2_backend_avx2.dir" if member ==
            "Leopard2BackendAVX2.cpp.o" else
            "leopard2_backend_ssse3.dir" if member ==
            "Leopard2BackendSSSE3.cpp.o" else "libleopard.dir"
        ])
    if build_paths != expected_build_paths:
        raise ValueError("dependency-file build paths or identities changed")
    return tuple(sorted(closure))


def validate_source_closure_paths(closure: object) -> list[str]:
    if not isinstance(closure, list) or len(closure) != len(
            EXPECTED_SOURCE_CLOSURE):
        raise ValueError("core source closure is incomplete")
    if not all(isinstance(entry, dict) for entry in closure):
        raise ValueError("core source closure record is malformed")
    paths = [entry.get("path") for entry in closure]
    if paths != list(EXPECTED_SOURCE_CLOSURE):
        raise ValueError("source closure exact set/order changed")
    return paths


def compiler_identity(record: dict, family: str) -> tuple[str, str]:
    first = record["version"].splitlines()[0]
    if family == "GNU":
        match = re.search(r"\b([0-9]+(?:\.[0-9]+){1,3})\Z", first)
    else:
        match = re.search(r"clang version ([0-9]+(?:\.[0-9]+){1,3})", first)
    if not match:
        raise ValueError("recorded compiler version is not understood")
    return family, match.group(1)


def validate_repository(
    source_root: pathlib.Path, tooling_sha: str, core_sha: str, *,
    require_checkout_head: bool,
) -> pathlib.Path:
    root = source_root.resolve()
    top = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
        check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if top.returncode != 0 or pathlib.Path(top.stdout.strip()).resolve() != root:
        raise ValueError("source root is not the exact Git worktree top level")
    ancestor = subprocess.run(
        ["git", "-C", str(root), "merge-base", "--is-ancestor", core_sha,
         tooling_sha], check=False, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)
    if ancestor.returncode != 0:
        raise ValueError("manifest core commit is not an ancestor of tooling")
    if require_checkout_head:
        head = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"], check=False,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if head.returncode != 0 or head.stdout.strip() != tooling_sha:
            raise ValueError("checkout HEAD differs from its tooling commit")
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain=v1",
             "--untracked-files=all"], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if status.returncode != 0 or status.stdout:
            raise ValueError("tooling checkout is not clean")
    return root


def validate_symbol_scan(
    text: str, target: str, *, archive: bool,
    expected_counts: dict[str, int],
    expected_members: dict[str, dict[str, int]],
) -> None:
    if text and not text.endswith("\n"):
        raise ValueError("symbol scan is not canonical line output")
    counts = {"asan_lines": 0, "ubsan_lines": 0}
    members = {name: {"asan_lines": 0, "ubsan_lines": 0}
               for name in expected_members}
    line_pattern = re.compile(
        r"(?:(?:[0-9a-f]{16}) |(?: {17}))([A-Za-z?]) (\S+)\Z")
    for line in text.splitlines():
        prefix = f"{target}:"
        if not line.startswith(prefix):
            raise ValueError("symbol scan target prefix changed")
        body = line[len(prefix):]
        member = None
        if archive:
            if ":" not in body:
                raise ValueError("archive scan omitted its member")
            member, body = body.split(":", 1)
            if member not in members:
                raise ValueError("archive scan named an unknown member")
        match = line_pattern.fullmatch(body)
        if not match:
            raise ValueError("symbol scan line is not canonical nm output")
        symbol = match.group(2)
        asan = "__asan_" in symbol
        ubsan = "__ubsan_" in symbol
        if not asan and not ubsan:
            raise ValueError("symbol scan contains a non-sanitizer symbol")
        counts["asan_lines"] += asan
        counts["ubsan_lines"] += ubsan
        if member is not None:
            members[member]["asan_lines"] += asan
            members[member]["ubsan_lines"] += ubsan
    if (not typed_equal(counts, expected_counts) or
            not typed_equal(members, expected_members)):
        raise ValueError("sanitizer counts or archive attribution changed")
    if expected_counts["asan_lines"]:
        for symbol in (
            "__asan_init", "__asan_report_load1", "__asan_report_store1",
            "__ubsan_handle_pointer_overflow",
            "__ubsan_handle_type_mismatch_v1",
        ):
            if symbol not in text:
                raise ValueError("sanitizer symbol family is incomplete")


def live_nm(
    record: dict, target: VerifiedBinaryArtifact,
    retained: VerifiedTextArtifact,
) -> None:
    nm = validate_program_live(record, "nm")
    with tempfile.TemporaryDirectory(prefix="c7-nm-snapshot-") as directory:
        snapshot = pathlib.Path(directory) / target.name
        snapshot.write_bytes(target.contents)
        snapshot.chmod(0o400)
        completed = subprocess.run(
            [str(nm), "--print-file-name", snapshot.name], cwd=snapshot.parent,
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completed.stderr:
        raise ValueError("live nm emitted stderr")
    lines = [line for line in completed.stdout.splitlines()
             if b"__asan_" in line or b"__ubsan_" in line]
    actual = b"\n".join(lines) + (b"\n" if lines else b"")
    if actual != retained.read_bytes():
        raise ValueError("live nm output differs from retained bytes")


def validate_cpp_result(
    data: dict[str, Any], backend: str, timing_scope: str,
    affinity: list[int], sanitizer: bool, *, legacy: bool = False,
) -> None:
    expected_correctness = (LEGACY_EXPECTED_CORRECTNESS if legacy else
                            dict(EXPECTED_CORRECTNESS))
    if not legacy:
        traced_backends = EXPECTED_TRACED_BACKENDS[backend]
        expected_correctness["concurrent_backend_contexts"] = traced_backends
        expected_correctness["concurrent_backend_executions"] = (
            traced_backends * 128)
        expected_correctness["concurrent_backend_trace_calls"] = (
            traced_backends * 10240)
        expected_correctness["concurrent_wire_digest_comparisons"] = (
            (traced_backends - 1) * 2)
    required = {
        "affinity", "allocation_tracking", "benchmarks", "core_git_sha",
        "correctness", "library_sha256", "omp_dynamic", "omp_num_threads",
        "production_constructor_rejected", "profile", "requested_backend",
        "runtime_backend", "sanitizer", "sanitizer_features", "schema",
        "source_sha256", "status", "timing_scope",
    }
    if not legacy:
        required.add("exact_byte_backend")
    if set(data) != required:
        raise ValueError("unexpected C7 child result keys")
    if (data["schema"] != "leopard2-c7-exact-low/v1" or
            data["status"] != "pass" or
            not typed_equal(data["profile"], EXPECTED_PROFILE) or
            data["production_constructor_rejected"] is not True or
            data["timing_scope"] != timing_scope or
            not typed_equal(data["affinity"], affinity) or
            data["omp_num_threads"] != "1" or data["omp_dynamic"] != "FALSE" or
            not typed_equal(data["correctness"], expected_correctness) or
            not isinstance(data["benchmarks"], list)):
        raise ValueError("C7 child identity or correctness result changed")
    validate_sha(data["source_sha256"], "C7 child source hash")
    validate_sha(data["library_sha256"], "C7 child library hash")
    validate_git_sha(data["core_git_sha"], "C7 child core commit")
    if sanitizer:
        if (data["requested_backend"] != "auto" or
                data["runtime_backend"] != EXPECTED_RUNTIME["auto"] or
                (not legacy and data["exact_byte_backend"] !=
                 EXPECTED_RUNTIME["auto"]) or
                data["sanitizer"] != "asan-ubsan" or
                data["allocation_tracking"] != "disabled-for-sanitizer" or
                not typed_equal(data["sanitizer_features"], {
                    "address": True, "undefined": True})):
            raise ValueError("sanitizer child provenance changed")
    elif (data["requested_backend"] != backend or
            data["runtime_backend"] != EXPECTED_RUNTIME[backend] or
            (not legacy and data["exact_byte_backend"] !=
             EXPECTED_RUNTIME[backend]) or
            data["sanitizer"] != "none" or
            data["allocation_tracking"] != "global-new" or
            not typed_equal(data["sanitizer_features"], {
                "address": False, "undefined": False})):
        raise ValueError("ordinary child backend provenance changed")
    if (timing_scope == "none-correctness-only" and
            not typed_equal(data["benchmarks"], [])):
        raise ValueError("correctness child unexpectedly contains timing")


def validate_smoke(
    data: dict[str, Any], affinity: list[int], *, legacy: bool = False,
) -> None:
    validate_cpp_result(
        data, "auto", "non-authoritative-smoke", affinity, False,
        legacy=legacy)
    cells = data["benchmarks"]
    if len(cells) != 1 or not isinstance(cells[0], dict):
        raise ValueError("smoke must contain exactly one benchmark cell")
    cell = cells[0]
    required = {
        "K", "R", "batch", "bytes", "exact_coefficients", "exact_decode",
        "exact_decode_samples_us", "exact_decode_setup",
        "exact_decode_setup_samples_us", "exact_decode_terms", "exact_encode",
        "exact_encode_samples_us", "exact_field", "exact_setup",
        "exact_setup_samples_us", "losses", "padded_decode",
        "padded_decode_samples_us", "padded_decode_scratch",
        "padded_decode_setup", "padded_decode_setup_samples_us",
        "padded_encode", "padded_encode_samples_us", "padded_encode_scratch",
        "padded_field", "padded_setup", "padded_setup_samples_us",
    }
    if set(cell) != required:
        raise ValueError("smoke benchmark cell schema changed")
    if not typed_equal(
            [cell[key] for key in ("K", "R", "bytes", "batch", "losses")],
            [3, 253, 64, 8, 3]):
        raise ValueError("smoke benchmark geometry changed")
    if not typed_equal({
        "exact_field": cell["exact_field"],
        "padded_field": cell["padded_field"],
        "exact_coefficients": cell["exact_coefficients"],
        "exact_decode_terms": cell["exact_decode_terms"],
    }, {
        "exact_field": 1,
        "padded_field": 2,
        "exact_coefficients": 759,
        "exact_decode_terms": 9,
    }):
        raise ValueError("smoke field or exact accounting changed")
    for key in ("padded_encode_scratch", "padded_decode_scratch"):
        if type(cell[key]) is not int or cell[key] < 0:
            raise ValueError("smoke scratch accounting is invalid")
    pairs = (
        ("exact_setup", "exact_setup_samples_us"),
        ("padded_setup", "padded_setup_samples_us"),
        ("exact_decode_setup", "exact_decode_setup_samples_us"),
        ("padded_decode_setup", "padded_decode_setup_samples_us"),
        ("exact_encode", "exact_encode_samples_us"),
        ("padded_encode", "padded_encode_samples_us"),
        ("exact_decode", "exact_decode_samples_us"),
        ("padded_decode", "padded_decode_samples_us"),
    )
    for summary_key, samples_key in pairs:
        samples = cell[samples_key]
        if (not isinstance(samples, list) or len(samples) != 7 or
                not all(type(value) in (int, float) and math.isfinite(value) and
                        value > 0 for value in samples)):
            raise ValueError("smoke raw sample schema changed")
        summary = cell[summary_key]
        median = statistics.median(samples)
        mad = statistics.median(abs(value - median) for value in samples)
        # Historical v3 evidence can serialize an exact zero MAD as integer 0
        # even though recomputing it from floating samples produces 0.0.  That
        # is a legitimate measurement representation difference, but bool is
        # never a numeric measurement: reject False/True explicitly instead of
        # relying on Python's bool-as-int equality.
        if (type(summary) is not dict or
                set(summary) != {"median_us", "mad_us"} or
                any(type(summary[key]) not in (int, float) for key in summary) or
                summary["median_us"] != median or
                summary["mad_us"] != mad):
            raise ValueError("smoke summary differs from raw samples")


def validate_manifest(
    data: dict[str, Any], *, source_root: pathlib.Path = ROOT,
    evidence_root: pathlib.Path | None = None, live: bool = False,
    require_checkout_head: bool = False,
) -> None:
    schema = data.get("schema")
    legacy = schema == LEGACY_MANIFEST_SCHEMA
    if schema not in (LEGACY_MANIFEST_SCHEMA, MANIFEST_SCHEMA):
        raise ValueError("manifest schema changed")
    required_top = {
        "builds", "core_git_sha", "normalization", "reproducibility",
        "runner", "runs", "schema", "scope", "source", "status",
        "taskset", "validator",
    }
    if not legacy:
        required_top.add("tooling_git_sha")
    if set(data) != required_top:
        raise ValueError("manifest keys changed")
    if (data["status"] != "pass" or data["scope"] !=
            "correctness plus one affinity-selected non-authoritative harness "
            "smoke; no promotion timing"):
        raise ValueError("manifest status or scope changed")
    if not typed_equal(data["normalization"], {
        "schema": NORMALIZATION_SCHEMA,
        "token": NORMALIZATION_TOKEN,
        "operation": "replace exact source-root prefix only",
    }):
        raise ValueError("normalization identity changed")
    source_root = source_root.resolve()
    evidence_root = (
        source_root if evidence_root is None else evidence_root.resolve())
    serialized = json.dumps(data, sort_keys=True)
    if str(source_root.resolve()) in serialized or ABSOLUTE_PROJECT_PATH.search(
            serialized):
        raise ValueError("manifest leaked an absolute checkout path")
    core_sha = validate_git_sha(data["core_git_sha"], "manifest core commit")
    tooling_sha = (core_sha if legacy else validate_git_sha(
        data["tooling_git_sha"], "manifest tooling commit"))
    validate_repository(
        source_root, tooling_sha, core_sha,
        require_checkout_head=require_checkout_head)
    for key, relative, label in (
        ("source", EXPECTED_TOOLING_CLOSURE[0], "C7 source"),
        ("runner", EXPECTED_TOOLING_CLOSURE[1], "C7 runner"),
        ("validator", EXPECTED_TOOLING_CLOSURE[2], "C7 validator"),
    ):
        # The current checkout contains v4 tooling, so historical v3 tooling
        # is authenticated from its immutable Git object rather than requiring
        # those old bytes to remain materialized at the same worktree path.
        validate_artifact(
            data[key], label, source_root, required=not legacy,
            check_if_present=not legacy)
        validate_git_artifact(
            tooling_sha, data[key], relative, label, source_root)
    taskset = validate_program_record(data["taskset"], "taskset")
    if live:
        validate_program_live(data["taskset"], "taskset")

    builds = data["builds"]
    if (not isinstance(builds, list) or not typed_equal(
            [item.get("name") for item in builds], list(BUILD_NAMES))):
        raise ValueError("build matrix changed")
    by_name: dict[str, dict] = {}
    build_root: str | None = None
    program_roles = (
        "ar", "c_compiler", "cmake", "cmake_linker", "compiler",
        "cxx_compiler", "gmake", "launcher_python", "link_driver", "nm", "ranlib",
        "standalone_linker",
    )
    for build in builds:
        required = {
            "ar", "argv_source_root_tokens", "backend", "build_argv",
            "build_dir", "build_stderr", "build_stdout", "c_compiler",
            "cmake", "cmake_cache", "cmake_linker", "compile_argv",
            "compile_stderr", "compile_stdout", "compiler",
            "configure_argv", "configure_stderr", "configure_stdout",
            "cxx_compiler", "dependency_files", "executable", "gmake",
            "instrumentation", "jobs", "launcher_python", "library",
            "link_driver", "name", "nm", "prefix_map_flags", "ranlib",
            "sanitizer", "source_closure", "standalone_linker",
        }
        if set(build) != required:
            raise ValueError("build record schema changed")
        name = build["name"]
        sanitizer = name == "asan-ubsan"
        if build["sanitizer"] is not sanitizer or build["backend"] != (
                "auto" if sanitizer else name):
            raise ValueError("backend/sanitizer label changed")
        records = {role: validate_program_record(build[role], f"{name} {role}")
                   for role in program_roles}
        if live:
            for role in program_roles:
                validate_program_live(build[role], f"{name} {role}")
            completed = subprocess.run(
                [str(records["compiler"]), "-print-prog-name=ld"],
                check=True, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            selected = pathlib.Path(completed.stdout.strip())
            if not selected.is_absolute():
                found = shutil.which(str(selected))
                if not found:
                    raise ValueError("live compiler-selected linker is unavailable")
                selected = pathlib.Path(found)
            if selected.resolve() != records["standalone_linker"].resolve():
                raise ValueError("live compiler-selected linker role changed")
        if (not typed_equal(build["compiler"], build["cxx_compiler"]) or
                not typed_equal(
                    build["link_driver"], build["cxx_compiler"])):
            raise ValueError("standalone compiler/link-driver role changed")
        token_counts = build["argv_source_root_tokens"]
        if not isinstance(token_counts, dict) or set(token_counts) != {
                "configure", "build", "compile"}:
            raise ValueError("build argv normalization record changed")
        configure = validate_argv(
            build["configure_argv"], token_counts["configure"],
            f"{name} configure", require_token=True, source_root=source_root)
        build_argv = validate_argv(
            build["build_argv"], token_counts["build"],
            f"{name} build", require_token=True, source_root=source_root)
        compile_argv = validate_argv(
            build["compile_argv"], token_counts["compile"],
            f"{name} compile", require_token=True, source_root=source_root)
        jobs = build["jobs"]
        if type(jobs) is not int or not 1 <= jobs <= 8:
            raise ValueError("build job count is not a typed value in 1..8")
        flags = build["prefix_map_flags"]
        mandatory = [
            f"-ffile-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}",
            f"-fdebug-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}",
        ]
        optional = f"-fmacro-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}"
        if not any(typed_equal(flags, expected) for expected in (
                mandatory, [*mandatory, optional])):
            raise ValueError("prefix-map flag set changed")
        build_dir = build["build_dir"]
        if (not isinstance(build_dir, str) or pathlib.Path(build_dir).is_absolute() or
                not build_dir.endswith(f"/core-{name}")):
            raise ValueError("build directory is not portable or canonical")
        try:
            resolve_path(build_dir, evidence_root)
        except ValueError as error:
            raise ValueError("build directory escapes the source checkout") from error
        candidate_root = build_dir.rsplit("/", 1)[0]
        if build_root is None:
            build_root = candidate_root
        elif candidate_root != build_root:
            raise ValueError("builds do not share one output root")
        expected_library = f"{build_dir}/liblibleopard.a"
        expected_executable = f"{build_root}/c7-{name}"
        if (build["library"]["path"] != expected_library or
                build["executable"]["path"] != expected_executable):
            raise ValueError("build output path changed")
        library = validate_artifact(
            build["library"], f"{name} archive", evidence_root, required=False,
            check_if_present=False)
        executable = validate_artifact(
            build["executable"], f"{name} executable", evidence_root,
            required=False, check_if_present=False)
        library_snapshot = executable_snapshot = None
        if live:
            library_snapshot = VerifiedBinaryArtifact(
                library, read_verified_artifact_bytes(
                    build["library"], f"{name} archive", evidence_root,
                    maximum_bytes=C7_BINARY_MAX_BYTES))
            executable_snapshot = VerifiedBinaryArtifact(
                executable, read_verified_artifact_bytes(
                    build["executable"], f"{name} executable", evidence_root,
                    maximum_bytes=C7_BINARY_MAX_BYTES))
        log_paths: dict[str, VerifiedTextArtifact] = {}
        for label in (
            "configure_stdout", "configure_stderr", "build_stdout",
            "build_stderr", "compile_stdout", "compile_stderr", "cmake_cache",
        ):
            log_paths[label] = validate_normalized_text(
                build[label], f"{name} {label}", evidence_root,
                require_token=label in {
                    "configure_stdout", "build_stdout", "compile_stderr",
                    "cmake_cache"}, checkout_root=source_root)

        prefix_text = " ".join(flags)
        core_flags = prefix_text + (
            " -fsanitize=address,undefined -fno-omit-frame-pointer"
            if sanitizer else "")
        linker_flags = (
            "-fsanitize=address,undefined -fno-omit-frame-pointer"
            if sanitizer else "")
        expected_configure = [
            str(records["cmake"]), "-S", NORMALIZATION_TOKEN, "-B",
            f"{NORMALIZATION_TOKEN}/{build_dir}", "-G", "Unix Makefiles",
            f"-DCMAKE_BUILD_TYPE={'Debug' if sanitizer else 'Release'}",
            f"-DCMAKE_C_COMPILER={records['c_compiler']}",
            f"-DCMAKE_CXX_COMPILER={records['cxx_compiler']}",
            f"-DLEO2_BACKEND_VARIANT={build['backend']}",
            "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=OFF",
            "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
            f"-DENABLE_OPENMP={'OFF' if sanitizer else 'ON'}",
            f"-DCMAKE_C_FLAGS={core_flags}",
            f"-DCMAKE_CXX_FLAGS={core_flags}",
            f"-DCMAKE_EXE_LINKER_FLAGS={linker_flags}",
            f"-DCMAKE_CXX_COMPILER_LAUNCHER={records['launcher_python']};"
            f"{NORMALIZATION_TOKEN}/{data['runner']['path']};--compiler-launch",
        ]
        expected_build = [
            str(records["cmake"]), "--build",
            f"{NORMALIZATION_TOKEN}/{build_dir}", "--target", "libleopard",
            "--verbose", "--", f"-j{jobs}",
        ]
        expected_compile_prefix = [
            str(records["compiler"]), "-v", "-Wl,-v", "-std=c++11", "-g",
            "-Wall", "-Wextra", "-Wpedantic", "-Werror", *flags,
            f"-I{NORMALIZATION_TOKEN}",
            f'-DLEO2_C7_SOURCE_SHA256="{data["source"]["sha256"]}"',
            f'-DLEO2_C7_CORE_GIT_SHA="{core_sha}"',
            f'-DLEO2_C7_LIBRARY_SHA256="{build["library"]["sha256"]}"',
        ]
        if sanitizer:
            expected_compile_prefix += [
                "-O1", "-fsanitize=address,undefined",
                "-fno-omit-frame-pointer",
                "-DLEO2_C7_DISABLE_GLOBAL_NEW_TRACKING=1",
                '-DLEO2_C7_SANITIZER_MODE="asan-ubsan"',
                "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1",
            ]
        else:
            expected_compile_prefix += ["-O2"]
        expected_compile_prefix += [
            data["source"]["path"], expected_library, "-pthread",
        ]
        if not sanitizer:
            expected_compile_prefix += ["-fopenmp"]
        expected_compile_prefix += [
            "-o", expected_executable]
        if (not typed_equal(configure, expected_configure) or
                not typed_equal(build_argv, expected_build) or
                not typed_equal(compile_argv, expected_compile_prefix)):
            raise ValueError("exact normalized build command changed")

        cache = parse_cache(log_paths["cmake_cache"].read_text(encoding="utf-8"))
        expected_cache = {
            "CMAKE_AR": str(records["ar"]),
            "CMAKE_BUILD_TYPE": "Debug" if sanitizer else "Release",
            "CMAKE_COMMAND": str(records["cmake"]),
            "CMAKE_CXX_COMPILER": str(records["cxx_compiler"]),
            "CMAKE_CXX_COMPILER_LAUNCHER": (
                f"{records['launcher_python']};{NORMALIZATION_TOKEN}/"
                f"{data['runner']['path']};--compiler-launch"),
            "CMAKE_CXX_FLAGS": core_flags,
            "CMAKE_C_COMPILER": str(records["c_compiler"]),
            "CMAKE_C_FLAGS": core_flags,
            "CMAKE_EXE_LINKER_FLAGS": linker_flags,
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_HOME_DIRECTORY": NORMALIZATION_TOKEN,
            "CMAKE_LINKER": str(records["cmake_linker"]),
            "CMAKE_MAKE_PROGRAM": str(records["gmake"]),
            "CMAKE_RANLIB": str(records["ranlib"]),
            "ENABLE_OPENMP": "OFF" if sanitizer else "ON",
            "LEO2_BACKEND_VARIANT": build["backend"],
            "LEO2_BUILD_BENCHMARKS": "OFF",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
            "leopard_BINARY_DIR": f"{NORMALIZATION_TOKEN}/{build_dir}",
        }
        for key, value in expected_cache.items():
            if cache.get(key) != value:
                raise ValueError(f"CMake cache role changed: {key}")
        configure_text = log_paths["configure_stdout"].read_text(encoding="utf-8")
        family = "Clang" if sanitizer else "GNU"
        _, c_version = compiler_identity(build["c_compiler"], family)
        _, cxx_version = compiler_identity(build["cxx_compiler"], family)
        required_configure = {
            f"-- The C compiler identification is {family} {c_version}",
            f"-- The CXX compiler identification is {family} {cxx_version}",
            f"-- Check for working C compiler: {records['c_compiler']} - skipped",
            f"-- Check for working CXX compiler: {records['cxx_compiler']} - skipped",
            "-- Detecting C compile features - done",
            "-- Detecting CXX compile features - done",
            "-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success",
            "-- Found Threads: TRUE  ",
            "-- Performing Test LEO2_FLAG_MSSSE3 - Success",
            "-- Performing Test LEO2_FLAG_MAVX2 - Success",
            f"-- Build files have been written to: {NORMALIZATION_TOKEN}/{build_dir}",
        }
        if not required_configure.issubset(set(configure_text.splitlines())):
            raise ValueError("configure log lacks required semantic events")
        if log_paths["configure_stderr"].read_bytes() or log_paths[
                "build_stderr"].read_bytes():
            raise ValueError("configure or core-build stderr is nonempty")
        build_text = log_paths["build_stdout"].read_text(encoding="utf-8")
        built_sources = set(re.findall(
            r"Building CXX object [^\n ]*/([^/ ]+\.cpp)\.o", build_text))
        if built_sources != CORE_SOURCES:
            raise ValueError("core build source set changed")
        for required in (
            str(records["cmake"]), str(records["cxx_compiler"]),
            str(records["gmake"]), str(records["ar"]), str(records["ranlib"]),
            str(records["launcher_python"]),
            f"{NORMALIZATION_TOKEN}/{data['runner']['path']}", *flags,
            "Built target leopard2_backend_ssse3",
            "Built target leopard2_backend_avx2",
            "Linking CXX static library liblibleopard.a", "Built target libleopard",
        ):
            if required not in build_text:
                raise ValueError("core build log lost an exact role or flag")
        compile_stdout = log_paths["compile_stdout"].read_text(encoding="utf-8")
        compile_stderr = log_paths["compile_stderr"].read_text(encoding="utf-8")
        linker_first = build["standalone_linker"]["version"].splitlines()[0]
        if compile_stdout != f"{linker_first}\n":
            raise ValueError("standalone linker output changed")
        for required in (
            "c7_exact_low.cpp", "liblibleopard.a", "-std=c++11",
            str(records["standalone_linker"]), NORMALIZATION_TOKEN,
            PREFIX_MAP_TARGET,
        ):
            if required not in compile_stderr:
                raise ValueError("standalone compile closure changed")

        instrumentation = build["instrumentation"]
        required_instrumentation = {
            "archive_members", "core_archive_counts",
            "core_archive_member_counts", "core_archive_symbol_scan",
            "executable_counts", "executable_symbol_scan",
            "required_compile_macro",
        }
        if set(instrumentation) != required_instrumentation:
            raise ValueError("instrumentation schema changed")
        zero_counts = {"asan_lines": 0, "ubsan_lines": 0}
        sanitizer_executable_proof = (
            LEGACY_EXPECTED_EXECUTABLE_SANITIZER_COUNTS if legacy else
            EXPECTED_EXECUTABLE_SANITIZER_COUNTS)
        sanitizer_archive_proof = (
            LEGACY_EXPECTED_ARCHIVE_SANITIZER_COUNTS if legacy else
            EXPECTED_ARCHIVE_SANITIZER_COUNTS)
        sanitizer_member_proof = (
            LEGACY_EXPECTED_ARCHIVE_MEMBER_COUNTS if legacy else
            EXPECTED_ARCHIVE_MEMBER_COUNTS)
        expected_executable_counts = (
            sanitizer_executable_proof if sanitizer else zero_counts)
        expected_archive_counts = (
            sanitizer_archive_proof if sanitizer else zero_counts)
        expected_members = (
            sanitizer_member_proof if sanitizer else
            {member: dict(zero_counts) for member in ARCHIVE_MEMBERS})
        if (instrumentation["required_compile_macro"] is not sanitizer or
                not typed_equal(
                    instrumentation["archive_members"], list(ARCHIVE_MEMBERS)) or
                not typed_equal(
                    instrumentation["executable_counts"],
                    expected_executable_counts) or
                not typed_equal(
                    instrumentation["core_archive_counts"],
                    expected_archive_counts) or
                not typed_equal(
                    instrumentation["core_archive_member_counts"],
                    expected_members)):
            raise ValueError("sanitizer summary or member attribution changed")
        executable_scan = validate_normalized_text(
            instrumentation["executable_symbol_scan"],
            f"{name} executable scan", evidence_root, require_token=False,
            checkout_root=source_root)
        archive_scan = validate_normalized_text(
            instrumentation["core_archive_symbol_scan"],
            f"{name} archive scan", evidence_root, require_token=False,
            checkout_root=source_root)
        validate_symbol_scan(
            executable_scan.read_text(encoding="utf-8"), executable.name,
            archive=False, expected_counts=expected_executable_counts,
            expected_members={})
        validate_symbol_scan(
            archive_scan.read_text(encoding="utf-8"), library.name,
            archive=True, expected_counts=expected_archive_counts,
            expected_members=expected_members)
        if live:
            if library_snapshot is None or executable_snapshot is None:
                raise AssertionError("live binary snapshots were not captured")
            live_nm(build["nm"], executable_snapshot, executable_scan)
            live_nm(build["nm"], library_snapshot, archive_scan)

        closure = build["source_closure"]
        validate_source_closure_paths(closure)
        for index, entry in enumerate(closure):
            validate_artifact(
                entry, f"{name} source closure {index}", source_root,
                required=not legacy, check_if_present=not legacy)
            relative = entry["path"]
            if pathlib.Path(relative).is_absolute():
                raise ValueError("source closure path is not portable")
            validate_git_artifact(
                core_sha, entry, relative, "source closure", source_root)
        derived_closure = dependency_source_closure(
            build["dependency_files"], build_dir, source_root, evidence_root,
            (pathlib.PurePosixPath(build["configure_stdout"]["path"]).parent /
             f"{name}-dependencies").as_posix(), live=live)
        if derived_closure != EXPECTED_SOURCE_CLOSURE:
            raise ValueError("dependency files do not reproduce source closure")
        if by_name and not typed_equal(
                closure, next(iter(by_name.values()))["source_closure"]):
            raise ValueError("equivalent builds have different source closures")
        by_name[name] = build

    reproducibility = data["reproducibility"]
    if not isinstance(reproducibility, dict) or set(reproducibility) != {
            "comparison", "fingerprints", "prefix_map_target"} or reproducibility[
            "prefix_map_target"] != PREFIX_MAP_TARGET:
        raise ValueError("reproducibility record changed")
    expected_fingerprints = {
        name: {
            "library_sha256": by_name[name]["library"]["sha256"],
            "executable_sha256": by_name[name]["executable"]["sha256"],
        }
        for name in BUILD_NAMES
    }
    if not typed_equal(reproducibility["fingerprints"], expected_fingerprints):
        raise ValueError("reproducibility fingerprints changed")
    normalized_text_record_count = sum(
        1 for item in _walk_dicts(data)
        if set(item) == {
            "bytes", "path", "sha256", "source_root_tokens"})
    validate_comparison(
        reproducibility["comparison"], expected_fingerprints,
        normalized_text_record_count)
    validate_peer_attestation(
        reproducibility["comparison"], data, source_root, evidence_root,
        live=live)

    runs = data["runs"]
    expected_runs = [*BUILD_NAMES, "smoke-nonauthoritative"]
    if (not isinstance(runs, list) or not typed_equal(
            [run.get("name") for run in runs], expected_runs)):
        raise ValueError("run matrix changed")
    correctness_cpus: list[int] = []
    for run in runs:
        if set(run) != {
            "argv", "argv_source_root_tokens", "build", "environment", "kind",
            "name", "observed_affinity", "requested_cpu", "result", "stderr",
            "stdout",
        }:
            raise ValueError("run record schema changed")
        name = run["name"]
        smoke = name == "smoke-nonauthoritative"
        build_name = "auto" if smoke else name
        if (run["build"] != build_name or run["kind"] != (
                "non-authoritative-smoke" if smoke else "correctness")):
            raise ValueError("run build or kind changed")
        cpu = run["requested_cpu"]
        if (type(cpu) is not int or cpu < 0 or
                not typed_equal(run["observed_affinity"], [cpu])):
            raise ValueError("run CPU or observed affinity changed")
        if not smoke:
            correctness_cpus.append(cpu)
        expected_environment = {
            "LC_ALL": "C", "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
        }
        if build_name == "asan-ubsan":
            expected_environment.update({
                "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
                "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
            })
        if not typed_equal(run["environment"], expected_environment):
            raise ValueError("run environment changed")
        argv = validate_argv(
            run["argv"], run["argv_source_root_tokens"], f"{name} run",
            require_token=True, source_root=source_root)
        result = validate_normalized_text(
            run["result"], f"{name} result", evidence_root,
            require_token=False, checkout_root=source_root)
        stdout = validate_normalized_text(
            run["stdout"], f"{name} stdout", evidence_root,
            require_token=False, checkout_root=source_root)
        stderr = validate_normalized_text(
            run["stderr"], f"{name} stderr", evidence_root,
            require_token=False, checkout_root=source_root)
        expected_argv = [
            str(taskset), "-c", str(cpu),
            f"{NORMALIZATION_TOKEN}/{by_name[build_name]['executable']['path']}",
            "--backend", "auto" if build_name == "asan-ubsan" else
            by_name[build_name]["backend"],
            (run["result"]["path"] if pathlib.Path(run["result"]["path"]).is_absolute()
             else f"{NORMALIZATION_TOKEN}/{run['result']['path']}"),
            "--benchmark-smoke" if smoke else "--correctness-only",
        ]
        if not typed_equal(argv, expected_argv):
            raise ValueError("normalized run command changed")
        if stdout.read_bytes() or stderr.read_text(encoding="utf-8") != (
                "C7 benchmark 1/1\n" if smoke else ""):
            raise ValueError("run log semantics changed")
        child = strict_json_loads(
            result.read_text(encoding="utf-8"), f"{name} child result")
        build = by_name[build_name]
        if (child.get("source_sha256") != data["source"]["sha256"] or
                child.get("core_git_sha") != core_sha or
                child.get("library_sha256") != build["library"]["sha256"]):
            raise ValueError("run result is not bound to build inputs")
        if smoke:
            validate_smoke(child, [cpu], legacy=legacy)
        else:
            validate_cpp_result(
                child, "auto" if build_name == "asan-ubsan" else build_name,
                "none-correctness-only", [cpu], build_name == "asan-ubsan",
                legacy=legacy)
    if len(correctness_cpus) != 5 or len(set(correctness_cpus)) != 5:
        raise ValueError("correctness runs did not use five distinct CPUs")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument(
        "--source-root", type=pathlib.Path, default=ROOT,
        help="checkout root containing the manifest's retained artifacts")
    parser.add_argument(
        "--evidence-root", type=pathlib.Path,
        help="optional immutable root for build/log/result artifacts")
    parser.add_argument(
        "--live", action="store_true",
        help="require exact recorded tools/build outputs and replay nm scans")
    parser.add_argument(
        "--require-checkout-head", action="store_true",
        help=("require a clean source-root HEAD equal to the manifest tooling "
              "commit"))
    arguments = parser.parse_args()
    data = strict_json_loads(
        arguments.manifest.read_text(encoding="utf-8"), "C7 manifest")
    validate_manifest(
        data, source_root=arguments.source_root,
        evidence_root=arguments.evidence_root, live=arguments.live,
        require_checkout_head=arguments.require_checkout_head)
    print("C7 evidence validation passed (live)" if arguments.live else
          "C7 evidence validation passed (portable)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
