#!/usr/bin/env python3
"""Shared, fail-closed helpers for balanced forced-decoder evidence."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import stat
import subprocess
import sys
from typing import Any, Sequence


MATRIX_SCHEMA = "leopard2-balanced-forced-matrix/v1"
RUN_SCHEMA = "leopard2-balanced-forced-abba/v2"
CHECKPOINT_SCHEMA = "leopard2-balanced-forced-checkpoint/v1"
RECORD_SCHEMA = "leopard2-balanced-forced-record/v1"
SUMMARY_SCHEMA = "leopard2-balanced-forced-summary/v2"
BENCHMARK_SCHEMA = "leopard2-benchmark-v3"
ROUND_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
MODES = ("generic", "materialized", "tiled")
BACKENDS = ("scalar", "ssse3", "avx2")
MODE_SELECTORS = {
    "generic": ("--force-generic",),
    "materialized": ("--force-specialized", "--force-materialized"),
    "tiled": ("--force-specialized", "--force-tiled"),
}
MODE_PARAMETERS = {
    "generic": {
        "force_generic_decode": True,
        "force_specialized_decode": False,
        "force_tiled_decode": False,
        "force_materialized_decode": False,
    },
    "materialized": {
        "force_generic_decode": False,
        "force_specialized_decode": True,
        "force_tiled_decode": False,
        "force_materialized_decode": True,
    },
    "tiled": {
        "force_generic_decode": False,
        "force_specialized_decode": True,
        "force_tiled_decode": True,
        "force_materialized_decode": False,
    },
}
CPU_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal")
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
RUNNER_RELATIVE = "experiments/leopard2/decoder_dispatch/run_balanced_abba.py"
ANALYZER_RELATIVE = "experiments/leopard2/decoder_dispatch/analyze_balanced.py"
COMMON_RELATIVE = "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
LEASE_RELATIVE = "tools/leopard2_jerasure_compare.py"
EXECUTABLE_LINK_VALUE_FREE_FLAGS = frozenset({
    "-Wall", "-Wextra", "-Wpedantic", "-fopenmp", "-g", "-DNDEBUG",
    "-O0", "-O1", "-O2", "-O3", "-Os", "-Oz",
})
EXTERNAL_LINK_INPUT_ROLES = {
    "libgomp.so": ("openmp_runtime_shared", "shared_library"),
    "libpthread.a": ("pthread_support_archive", "archive"),
}
EXTERNAL_LINK_INPUT_ORDER = (
    "openmp_runtime_shared", "pthread_support_archive",
)
MAX_EXTERNAL_LINK_INPUT_BYTES = 64 * 1024 * 1024
EFFECTIVE_FLAG_ALLOWLISTS = {
    "release": frozenset({"-g", "-O0", "-O3", "-DNDEBUG"}),
    "link": frozenset({
        "-Wall", "-Wextra", "-Wpedantic", "-fopenmp", "-g", "-O0",
        "-O3", "-DNDEBUG", "-o",
    }),
    "compile": frozenset({
        "-Wall", "-Wextra", "-Wpedantic", "-fopenmp", "-g", "-O0",
        "-O3", "-DNDEBUG", "-std=gnu++11", "-march=native",
        "-mssse3", "-mno-avx", "-mavx2", "-mno-avx512f",
        "-mavx512f", "-mavx512bw", "-mavx512vl",
        "-mprefer-vector-width=256", "-falign-functions=64", "-o", "-c",
    }),
}


class EvidenceError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def validate_effective_flags(
    tokens: Sequence[str], label: str, policy: str,
) -> None:
    """Apply one context-specific allowlist and require effective final -O3."""
    require(isinstance(tokens, Sequence) and
            not isinstance(tokens, (str, bytes)) and
            all(isinstance(token, str) for token in tokens),
            f"{label} flag stream is invalid")
    allowed = EFFECTIVE_FLAG_ALLOWLISTS.get(policy)
    require(allowed is not None, f"{label} effective-flag policy is unknown")
    flags = (list(tokens) if policy == "release" else
             [token for token in tokens if token.startswith("-")])
    unknown = [
        token for token in flags
        if token not in allowed
    ]
    require(not unknown,
            f"{label} contains noncanonical or ambiguous flags: {unknown}")
    optimizations = [token for token in tokens if token.startswith("-O")]
    require(optimizations and optimizations[-1] == "-O3",
            f"{label} final optimization flag is not -O3: {optimizations}")


def validate_external_link_operand_path(
    operand: object, role: str, label: str,
) -> str:
    """Accept only the two canonical explicit GCC OpenMP support inputs."""
    require(isinstance(operand, str) and operand and
            Path(operand).is_absolute() and
            os.path.normpath(operand) == operand and
            not any(character.isspace() for character in operand) and
            "@" not in operand and "\\" not in operand,
            f"{label} {role} operand path is not canonical")
    if role == "openmp_runtime_shared":
        require(re.fullmatch(
                    r"/usr/lib/gcc/x86_64-linux-gnu/"
                    r"[0-9]+(?:\.[0-9]+)*/libgomp\.so",
                    operand) is not None,
                f"{label} OpenMP runtime is outside the canonical GCC root")
    elif role == "pthread_support_archive":
        require(re.fullmatch(
                    r"/(?:usr/)?lib(?:64|/x86_64-linux-gnu)/libpthread\.a",
                    operand) is not None,
                f"{label} pthread archive is outside a canonical system root")
    else:
        raise EvidenceError(f"{label} external-link role is unknown: {role!r}")
    return operand


def _stat_identity(metadata: os.stat_result) -> tuple[int, ...]:
    """Return every mutable identity field used by an external snapshot."""
    return (
        metadata.st_mode, metadata.st_nlink, metadata.st_size,
        metadata.st_mtime_ns, metadata.st_ctime_ns,
        metadata.st_dev, metadata.st_ino,
    )


def _resolve_external_lexical_path(
    path: Path, label: str,
) -> tuple[Path, list[dict[str, Any]], tuple[tuple[object, ...], ...]]:
    """Resolve and inventory every symlink traversed by an absolute operand.

    ``Path.resolve`` returns only the terminal path, which is insufficient for
    provenance: either the public ``libgomp.so`` link or a link named by one of
    its targets can be retargeted without changing the terminal inode already
    open by the validator.  This small realpath walk retains the complete
    lexical chain and a private inode/mode snapshot used for the after-read
    comparison.
    """
    require(path.is_absolute() and os.path.normpath(str(path)) == str(path),
            f"{label} lexical path is not canonical")
    pending = list(path.parts[1:])
    resolved_parts: list[str] = []
    public_chain: list[dict[str, Any]] = []
    private_chain: list[tuple[object, ...]] = []
    traversals = 0
    while pending:
        component = pending.pop(0)
        if component in {"", "."}:
            continue
        if component == "..":
            require(resolved_parts, f"{label} lexical path escapes root")
            resolved_parts.pop()
            continue
        candidate = Path("/", *resolved_parts, component)
        metadata = os.lstat(candidate)
        if stat.S_ISLNK(metadata.st_mode):
            traversals += 1
            require(traversals <= 40 and metadata.st_nlink == 1,
                    f"{label} lexical symlink chain is unsafe")
            target = os.readlink(candidate)
            require(target and "\x00" not in target,
                    f"{label} lexical symlink target is invalid")
            public_chain.append({
                "path": str(candidate),
                "target": target,
                "mode": stat.S_IMODE(metadata.st_mode),
            })
            private_chain.append((
                str(candidate), target, *_stat_identity(metadata)))
            target_path = Path(target)
            if target_path.is_absolute():
                resolved_parts = []
                target_components = list(target_path.parts[1:])
            else:
                target_components = list(target_path.parts)
            pending = target_components + pending
            continue
        require(not pending or stat.S_ISDIR(metadata.st_mode),
                f"{label} lexical parent is not a directory")
        resolved_parts.append(component)
    resolved = Path("/", *resolved_parts)
    require(resolved.is_absolute(), f"{label} did not resolve absolutely")
    return resolved, public_chain, tuple(private_chain)


def current_external_file_snapshot(
    path: Path, label: str,
) -> tuple[os.stat_result, str, bytes, Path, list[dict[str, Any]]]:
    """Return a private immutable snapshot bound to the full lexical chain.

    The exact bytes are retained until all digest and file-format consumers
    have finished.  Therefore an in-place rewrite after the descriptor read
    cannot turn a stale digest into an accepted description of newly consumed
    bytes; the closure consumes this immutable ``bytes`` value, then verifies
    the lexical path and descriptor identity again.
    """
    resolved, chain_before, private_chain_before = \
        _resolve_external_lexical_path(path, label)
    before = os.lstat(resolved)
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
            0 < before.st_size <= MAX_EXTERNAL_LINK_INPUT_BYTES,
            f"{label} is not a bounded single-link regular file")
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        initial = os.fstat(descriptor)
        path_initial = os.lstat(resolved)
        require(stat.S_ISREG(initial.st_mode) and initial.st_nlink == 1 and
                (initial.st_dev, initial.st_ino) ==
                (before.st_dev, before.st_ino) ==
                (path_initial.st_dev, path_initial.st_ino) and
                _stat_identity(initial) == _stat_identity(before) ==
                _stat_identity(path_initial) and
                0 < initial.st_size <= MAX_EXTERNAL_LINK_INPUT_BYTES,
                f"{label} changed before its identity read")
        blocks: list[bytes] = []
        retained = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
            retained += len(block)
            require(retained <= MAX_EXTERNAL_LINK_INPUT_BYTES,
                    f"{label} exceeds its identity bound")
        snapshot = b"".join(blocks)
        require(len(snapshot) == retained,
                f"{label} immutable snapshot length differs")
        final = os.fstat(descriptor)
        resolved_after, chain_after, private_chain_after = \
            _resolve_external_lexical_path(path, label)
        path_final = os.lstat(resolved_after)
        descriptor_final = os.fstat(descriptor)
        require(retained == initial.st_size and
                stat.S_ISREG(final.st_mode) and final.st_nlink == 1 and
                stat.S_ISREG(path_final.st_mode) and path_final.st_nlink == 1 and
                resolved_after == resolved and
                chain_after == chain_before and
                private_chain_after == private_chain_before and
                (final.st_dev, final.st_ino) ==
                (initial.st_dev, initial.st_ino) ==
                (path_final.st_dev, path_final.st_ino) and
                _stat_identity(final) == _stat_identity(initial) ==
                _stat_identity(path_final) == _stat_identity(descriptor_final),
                f"{label} changed during its identity read")
        return (
            initial, hashlib.sha256(snapshot).hexdigest(), snapshot,
            resolved, chain_before,
        )
    finally:
        os.close(descriptor)


def validate_external_link_input_shape(
    value: object, label: str,
) -> list[dict[str, Any]]:
    """Bind grammar operands to one ordered identity record per input."""
    require(isinstance(value, list) and
            len(value) == len(EXTERNAL_LINK_INPUT_ORDER),
            f"{label} external-link identity closure differs")
    records: list[dict[str, Any]] = []
    for index, expected_role in enumerate(EXTERNAL_LINK_INPUT_ORDER):
        record = value[index]
        require(isinstance(record, dict) and set(record) == {
                    "operand", "role", "artifact",
                    "lexical_symlink_chain"} and
                record.get("role") == expected_role and
                isinstance(record.get("artifact"), dict),
                f"{label} external-link identity {index} shape differs")
        operand = validate_external_link_operand_path(
            record.get("operand"), expected_role, label)
        _, expected_kind = EXTERNAL_LINK_INPUT_ROLES[Path(operand).name]
        artifact = record["artifact"]
        require(artifact.get("kind") == expected_kind and
                isinstance(artifact.get("path"), str) and
                Path(artifact["path"]).is_absolute() and
                type(artifact.get("size")) is int and artifact["size"] > 0 and
                isinstance(artifact.get("sha256"), str) and
                re.fullmatch(r"[0-9a-f]{64}", artifact["sha256"]) is not None,
                f"{label} external artifact identity {index} is incomplete")
        chain = record.get("lexical_symlink_chain")
        require(isinstance(chain, list) and
                all(isinstance(link, dict) and set(link) == {
                        "path", "target", "mode"} and
                    isinstance(link.get("path"), str) and
                    Path(link["path"]).is_absolute() and
                    os.path.normpath(link["path"]) == link["path"] and
                    isinstance(link.get("target"), str) and link["target"] and
                    type(link.get("mode")) is int and
                    0 <= link["mode"] <= 0o7777
                    for link in chain),
                f"{label} external lexical symlink chain {index} is invalid")
        try:
            metadata, current_sha256, snapshot, resolved, current_chain = \
                current_external_file_snapshot(
                    Path(operand), f"{label} external operand {index}")
        except (OSError, RuntimeError) as error:
            raise EvidenceError(
                f"{label} external operand {operand!r} does not resolve: {error}") \
                from error
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                0 < metadata.st_size <= MAX_EXTERNAL_LINK_INPUT_BYTES and
                artifact["path"] == str(resolved) and
                artifact["size"] == metadata.st_size and
                artifact.get("mode") == stat.S_IMODE(metadata.st_mode) and
                artifact["sha256"] == current_sha256 and
                chain == current_chain and
                ("mtime_ns" not in artifact or
                 artifact["mtime_ns"] == metadata.st_mtime_ns),
                f"{label} external operand does not match its current resolved identity")
        require((expected_kind == "archive" and
                 snapshot.startswith(b"!<arch>\n")) or
                (expected_kind == "shared_library" and
                 snapshot.startswith(b"\x7fELF")),
                f"{label} external operand has the wrong file format")
        if expected_kind == "archive":
            require(artifact["path"] == operand and artifact["size"] >= 8 and
                    not chain,
                    f"{label} pthread operand does not bind its exact archive")
        else:
            require(re.fullmatch(
                        r"/(?:usr/)?lib(?:64|/x86_64-linux-gnu)/"
                        r"libgomp\.so(?:\.[0-9]+)+",
                        artifact["path"]) is not None,
                    f"{label} OpenMP operand resolves outside its runtime root")
            require(chain and chain[0]["path"] == operand,
                    f"{label} OpenMP operand does not bind its lexical symlink")
        require(
                ((expected_kind == "archive" and
                  Path(artifact["path"]).name == "libpthread.a") or
                 (expected_kind == "shared_library" and
                  re.fullmatch(r"libgomp\.so(?:\.[0-9]+)+",
                               Path(artifact["path"]).name) is not None)),
                f"{label} external operand resolves to a different library")
        records.append(record)
    operands = [record["operand"] for record in records]
    resolved = [record["artifact"]["path"] for record in records]
    require(len(operands) == len(set(operands)) and
            len(resolved) == len(set(resolved)),
            f"{label} external-link inputs are duplicated or aliased")
    return records


def validate_executable_link_semantics(
    tokens: Sequence[str], *, compiler_invocation: str, archive_name: str,
    executable_name: str, benchmark_object: str,
    external_link_inputs: object, label: str,
) -> None:
    """Consume every executable-link token with one fail-closed production."""
    require(isinstance(tokens, Sequence) and
            not isinstance(tokens, (str, bytes)) and tokens and
            all(isinstance(token, str) and token for token in tokens),
            f"{label} token stream is invalid")
    require(isinstance(compiler_invocation, str) and compiler_invocation and
            isinstance(archive_name, str) and Path(archive_name).name == archive_name and
            isinstance(executable_name, str) and
            Path(executable_name).name == executable_name and
            isinstance(benchmark_object, str) and benchmark_object.endswith(".o") and
            not Path(benchmark_object).is_absolute() and
            "\\" not in benchmark_object and "@" not in benchmark_object,
            f"{label} expected semantic closure is invalid")
    validate_effective_flags(tokens, label, "link")
    require(isinstance(external_link_inputs, list) and
            len(external_link_inputs) == len(EXTERNAL_LINK_INPUT_ORDER),
            f"{label} external-link grammar closure differs")
    external_operands: list[str] = []
    for index, expected_role in enumerate(EXTERNAL_LINK_INPUT_ORDER):
        record = external_link_inputs[index]
        require(isinstance(record, dict) and
                record.get("role") == expected_role,
                f"{label} external-link grammar role {index} differs")
        external_operands.append(validate_external_link_operand_path(
            record.get("operand"), expected_role, label))
    require(tokens[0] == compiler_invocation,
            f"{label} compiler invocation differs")
    seen_archive = 0
    seen_object = 0
    seen_output = 0
    seen_external: list[str] = []
    index = 1
    while index < len(tokens):
        token = tokens[index]
        if token in EXECUTABLE_LINK_VALUE_FREE_FLAGS:
            index += 1
            continue
        if token == "-o":
            require(index + 1 < len(tokens) and
                    tokens[index + 1] == executable_name,
                    f"{label} output operand differs")
            seen_output += 1
            index += 2
            continue
        if token == archive_name:
            seen_archive += 1
            index += 1
            continue
        if token == benchmark_object:
            seen_object += 1
            index += 1
            continue
        if token in external_operands:
            seen_external.append(token)
            index += 1
            continue
        raise EvidenceError(
            f"{label} contains an unrecognized or controlling token: {token!r}")
    require(seen_archive == 1 and seen_object == 1 and seen_output == 1 and
            seen_external == external_operands,
            f"{label} declared input/output closure differs")
    # Snapshot the external operands only after every grammar token has been
    # consumed.  File-format and digest checks use immutable byte snapshots,
    # and the snapshot helper re-resolves the complete lexical chain after its
    # read, so there is no later path-based consumer that could accept stale
    # same-inode bytes based on timestamps alone.
    validate_external_link_input_shape(external_link_inputs, label)


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise EvidenceError(f"cannot read {label} {path}: {error}") from error
    require(isinstance(value, dict), f"{label} is not a JSON object")
    return value, raw


def require_keys(value: object, keys: set[str], label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == keys,
            f"{label} keys changed: expected {sorted(keys)}, got {sorted(value)}")
    return value


def strict_int(value: object, label: str, minimum: int = 0) -> int:
    require(type(value) is int and value >= minimum,
            f"{label} is not an integer >= {minimum}")
    return value


def hex_digest(value: object, label: str, length: int = 64) -> str:
    require(isinstance(value, str) and len(value) == length and value == value.lower(),
            f"{label} is not lowercase {length}-digit hexadecimal")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{label} is not hexadecimal") from error
    return value


def finite(value: object, label: str, minimum: float = 0.0) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result >= minimum,
            f"{label} is not finite and >= {minimum}")
    return result


def checked_output(arguments: list[str]) -> str:
    try:
        result = subprocess.run(
            arguments, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"command failed: {arguments!r}: {error}") from error
    require(not result.stderr, f"command emitted stderr: {arguments!r}")
    return result.stdout.strip()


def file_identity(path: Path, root: Path | None, label: str) -> dict[str, Any]:
    path = path.absolute()
    require(path.is_file(), f"missing {label}: {path}")
    descriptor = path.stat()
    relative: str | None = None
    if root is not None:
        try:
            relative = str(path.resolve().relative_to(root.resolve()))
        except ValueError:
            relative = None
    return {
        "path": str(path),
        "realpath": str(path.resolve()),
        "relative_path": relative,
        "sha256": sha256_file(path),
        "size": descriptor.st_size,
        "mode": stat.S_IMODE(descriptor.st_mode),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
        "mtime_ns": descriptor.st_mtime_ns,
    }


def git_identity(root: Path, expected_commit: str) -> dict[str, Any]:
    root = root.resolve()
    hex_digest(expected_commit, "expected commit", 40)
    head = checked_output(["git", "-C", str(root), "rev-parse", "HEAD"])
    tree = checked_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"])
    expected_tree = checked_output(
        ["git", "-C", str(root), "rev-parse", expected_commit + "^{tree}"])
    status = checked_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=all"])
    require(head == expected_commit, f"source HEAD {head} != expected {expected_commit}")
    require(tree == expected_tree, "source HEAD tree differs from expected commit tree")
    require(not status, f"source tree is dirty: {status!r}")
    return {
        "root": str(root), "head": head, "tree": tree,
        "expected_commit": expected_commit, "status": "clean",
        "status_sha256": EMPTY_SHA256,
    }


def find_build_root(binary: Path) -> Path:
    for candidate in (binary.parent, *binary.parents):
        if (candidate / "CMakeCache.txt").is_file():
            return candidate
    raise EvidenceError(f"cannot find CMakeCache.txt above {binary}")


def refresh_build(binary: Path) -> list[str]:
    root = find_build_root(binary)
    command = ["cmake", "--build", str(root), "--target", "bench_leopard2", "-j", "1"]
    try:
        completed = subprocess.run(
            command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"benchmark build refresh failed: {error}") from error
    return command


def _unique_file(candidates: list[Path], label: str) -> Path:
    unique = sorted({item.resolve() for item in candidates if item.is_file()})
    require(len(unique) == 1,
            f"expected exactly one {label}, found {[str(item) for item in unique]}")
    return unique[0]


def build_identity(source_root: Path, binary_relative: str) -> dict[str, Any]:
    source_root = source_root.resolve()
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "benchmark binary path must be source-root-relative")
    binary = source_root / relative
    require(binary.is_file() and os.access(binary, os.X_OK),
            f"benchmark executable is absent or not executable: {binary}")
    build_root = find_build_root(binary)
    cache = build_root / "CMakeCache.txt"
    homes = [line.split("=", 1)[1] for line in cache.read_text(
        encoding="utf-8", errors="strict").splitlines()
             if line.startswith("CMAKE_HOME_DIRECTORY:INTERNAL=")]
    require(len(homes) == 1 and Path(homes[0]).resolve() == source_root,
            "CMake cache does not name the clean source root")
    benchmark_source = source_root / "bench/leopard2/benchmark.cpp"
    decoder_source = source_root / "leopard2.cpp"
    dispatch_source = source_root / "Leopard2Dispatch.h"
    benchmark_object = _unique_file([
        item for item in build_root.rglob("benchmark.cpp.o")
        if "bench_leopard2.dir" in str(item.parent)
    ] + [
        item for item in build_root.rglob("benchmark.cpp.obj")
        if "bench_leopard2.dir" in str(item.parent)
    ], "bench_leopard2 benchmark object")
    decoder_object = _unique_file([
        item for item in build_root.rglob("leopard2.cpp.o")
        if "leopard.dir" in str(item.parent) and "test_hooks" not in str(item)
    ] + [
        item for item in build_root.rglob("leopard2.cpp.obj")
        if "leopard.dir" in str(item.parent) and "test_hooks" not in str(item)
    ], "production leopard2 object")
    archive_names = ("libleopard.a", "leopard.lib")
    archive = _unique_file(
        [item for name in archive_names for item in build_root.rglob(name)],
        "production Leopard archive")
    graph_candidates = [build_root / "build.ninja", build_root / "Makefile"]
    graph = _unique_file(graph_candidates, "CMake build graph")
    identities = {
        "benchmark": file_identity(benchmark_source, source_root, "benchmark source"),
        "decoder": file_identity(decoder_source, source_root, "decoder source"),
        "dispatch": file_identity(dispatch_source, source_root, "dispatch source"),
    }
    objects = {
        "benchmark": file_identity(benchmark_object, source_root, "benchmark object"),
        "decoder": file_identity(decoder_object, source_root, "decoder object"),
    }
    archive_identity = file_identity(archive, source_root, "production archive")
    binary_identity = file_identity(binary, source_root, "benchmark executable")
    require(objects["benchmark"]["mtime_ns"] <= binary_identity["mtime_ns"],
            "benchmark executable predates its benchmark object")
    require(objects["decoder"]["mtime_ns"] <= archive_identity["mtime_ns"] <=
            binary_identity["mtime_ns"], "executable/archive/object timestamps are stale")
    return {
        "root": str(build_root),
        "cmake_home": str(source_root),
        "cache": file_identity(cache, source_root, "CMake cache"),
        "graph": file_identity(graph, source_root, "CMake build graph"),
        "sources": identities,
        "objects": objects,
        "archive": archive_identity,
        "binary": binary_identity,
    }


def _external_file_identity(path: Path, label: str) -> dict[str, Any]:
    return file_identity(path, None, label)


def stable_cpuinfo(cpus: list[int]) -> list[dict[str, str]]:
    """Retain CPU identity fields while excluding frequency/load telemetry."""
    wanted = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "bugs", "address sizes",
        "CPU implementer", "CPU architecture", "CPU variant", "CPU part",
        "CPU revision", "Features",
    }
    blocks = Path("/proc/cpuinfo").read_text(
        encoding="utf-8", errors="strict").strip().split("\n\n")
    by_processor: dict[int, dict[str, str]] = {}
    for block in blocks:
        record: dict[str, str] = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            key, value = (item.strip() for item in line.split(":", 1))
            if key in wanted:
                record[key] = value
        try:
            processor = int(record.get("processor", "-1"))
        except ValueError:
            continue
        by_processor[processor] = record
    result = []
    for cpu in cpus:
        require(cpu in by_processor, f"CPU {cpu} is absent from stable /proc/cpuinfo")
        result.append(by_processor[cpu])
    return result


def runtime_identity(binary: Path, cpu: int, sibling: int,
                     affinity: set[int]) -> dict[str, Any]:
    require(cpu in affinity and sibling in affinity and cpu != sibling,
            "runtime CPU pair is outside affinity or not distinct")
    records = []
    expected = sorted((cpu, sibling))
    for logical in expected:
        topology = Path(f"/sys/devices/system/cpu/cpu{logical}/topology")
        try:
            sibling_text = (topology / "thread_siblings_list").read_text(
                encoding="ascii").strip()
            core_id = (topology / "core_id").read_text(encoding="ascii").strip()
            package_id = (topology / "physical_package_id").read_text(
                encoding="ascii").strip()
        except OSError as error:
            raise EvidenceError(f"cannot read CPU {logical} topology: {error}") from error
        parsed: list[int] = []
        for part in sibling_text.split(","):
            bounds = part.split("-", 1)
            start = int(bounds[0])
            stop = int(bounds[-1])
            parsed.extend(range(start, stop + 1))
        require(sorted(set(parsed)) == expected,
                f"CPU {logical} is not in the exact SMT pair {expected}")
        records.append({
            "cpu": logical, "thread_siblings_list": sibling_text,
            "thread_siblings": sorted(set(parsed)), "core_id": core_id,
            "physical_package_id": package_id,
        })
    require(records[0]["core_id"] == records[1]["core_id"] and
            records[0]["physical_package_id"] == records[1]["physical_package_id"],
            "CPU pair does not share one physical core")
    try:
        ldd = subprocess.run(
            ["ldd", str(binary)], check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"cannot inspect benchmark runtime libraries: {error}") from error
    require(not ldd.stderr, "ldd emitted stderr")
    library_paths = []
    for line in ldd.stdout.splitlines():
        words = line.strip().split()
        candidate = None
        if len(words) >= 3 and words[1] == "=>" and words[2].startswith("/"):
            candidate = words[2]
        elif words and words[0].startswith("/"):
            candidate = words[0]
        if candidate:
            library_paths.append(Path(candidate))
    libraries = [_external_file_identity(path, "runtime library")
                 for path in sorted(set(library_paths))]
    normalized_ldd = re.sub(r"\s+\(0x[0-9a-fA-F]+\)", "", ldd.stdout)
    uname = platform.uname()
    environment = {key: os.environ.get(key) for key in (
        "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT", "OMP_NUM_THREADS",
        "OMP_DYNAMIC", "OMP_PROC_BIND")}
    require(not environment["LD_PRELOAD"] and not environment["LD_AUDIT"],
            "LD_PRELOAD and LD_AUDIT must be unset for authoritative evidence")
    return {
        "platform": {
            "system": uname.system, "node": uname.node, "release": uname.release,
            "version": uname.version, "machine": uname.machine,
        },
        "python": {
            "version": sys.version, "implementation": platform.python_implementation(),
            "executable": _external_file_identity(
                Path(sys.executable), "Python executable"),
            "byteorder": sys.byteorder,
        },
        "affinity": sorted(affinity),
        "cpu": cpu,
        "sibling": sibling,
        "topology": records,
        "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
        "cpuinfo": stable_cpuinfo(expected),
        "runtime_libraries": libraries,
        "ldd_normalized_sha256": sha256_bytes(normalized_ldd.encode("utf-8")),
        "environment": environment,
    }


def normalize_matrix(value: object) -> tuple[str, list[dict[str, Any]]]:
    matrix = require_keys(value, {"schema", "name", "cases"}, "matrix")
    require(matrix["schema"] == MATRIX_SCHEMA, "unexpected matrix schema")
    name = matrix["name"]
    require(isinstance(name, str) and name, "matrix name is empty")
    raw_cases = matrix["cases"]
    require(isinstance(raw_cases, list) and raw_cases, "matrix cases are empty")
    cases = []
    names = set()
    keys = {
        "name", "K", "R", "profile", "field", "backend", "shard_bytes",
        "loss_count", "batch", "reuse", "iterations", "warmup", "seed",
        "control_mode", "candidate_mode",
    }
    for index, raw in enumerate(raw_cases):
        case = require_keys(raw, keys, f"matrix case {index}")
        case_name = case["name"]
        require(isinstance(case_name, str) and case_name and case_name not in names,
                f"matrix case {index} name is empty or duplicated")
        names.add(case_name)
        for key in ("K", "R", "shard_bytes", "loss_count", "batch", "reuse",
                    "iterations", "seed"):
            strict_int(case[key], f"case {case_name} {key}", 1)
        strict_int(case["warmup"], f"case {case_name} warmup")
        require(case["profile"] == "legacy_high_v1" and case["field"] == "gf8",
                f"case {case_name} is outside balanced legacy-high GF8")
        k = case["K"]
        r = case["R"]
        require(k == r and 5 <= k <= 128 and case["loss_count"] == k,
                f"case {case_name} is not balanced full-original recovery")
        side = 1 << (r - 1).bit_length()
        parent = 1 << (k + side - 1).bit_length()
        require(side >= 8 and parent == 2 * side and parent <= 256,
                f"case {case_name} has an invalid balanced GF8 parent")
        require(case["backend"] in BACKENDS,
                f"case {case_name} uses an unsupported backend")
        control = case["control_mode"]
        candidate = case["candidate_mode"]
        require(control in MODES and candidate in MODES and control != candidate,
                f"case {case_name} has invalid or identical forced modes")
        normalized = dict(case)
        normalized["padded_side"] = side
        normalized["parent_count"] = parent
        cases.append(normalized)
    return name, cases


def role_mode(case: dict[str, Any], role: str) -> str:
    require(role in ("control", "candidate"), f"unknown role: {role}")
    return case[role + "_mode"]


def benchmark_command(binary: Path, case: dict[str, Any], output: Path,
                      cpu: int, role: str) -> list[str]:
    mode = role_mode(case, role)
    return [
        "/usr/bin/taskset", "-c", str(cpu),
        "/usr/bin/env", "OMP_NUM_THREADS=1", "OMP_DYNAMIC=false",
        "OMP_PROC_BIND=close", str(binary.absolute()),
        "--k", str(case["K"]), "--r", str(case["R"]), "--profile", "high",
        "--field", "gf8", "--backend", case["backend"],
        "--bytes", str(case["shard_bytes"]), "--loss", str(case["loss_count"]),
        "--batch", str(case["batch"]), "--reuse", str(case["reuse"]),
        "--iterations", str(case["iterations"]), "--warmup", str(case["warmup"]),
        "--threads", "1", "--seed", str(case["seed"]),
        *MODE_SELECTORS[mode], "--skip-legacy", "--retain-samples",
        "--report-decode-path",
        "--json", str(output.absolute()),
    ]


def cpu_snapshot(cpu: int) -> dict[str, int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            values = [int(item) for item in line.split()[1:1 + len(CPU_FIELDS)]]
            require(len(values) == len(CPU_FIELDS), f"CPU {cpu} stat row is incomplete")
            return dict(zip(CPU_FIELDS, values))
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def cpu_delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    require(set(before) == set(CPU_FIELDS) and set(after) == set(CPU_FIELDS),
            "CPU snapshot fields changed")
    result = {key: after[key] - before[key] for key in CPU_FIELDS}
    require(all(value >= 0 for value in result.values()), "CPU counters moved backwards")
    result["idle_total"] = result["idle"] + result["iowait"]
    result["nonidle_total"] = sum(
        result[key] for key in CPU_FIELDS if key not in {"idle", "iowait"})
    result["total"] = result["idle_total"] + result["nonidle_total"]
    return result


def isolation(before_cpu: dict[str, int], after_cpu: dict[str, int],
              before_sibling: dict[str, int], after_sibling: dict[str, int]) -> dict[str, Any]:
    timed = cpu_delta(before_cpu, after_cpu)
    sibling = cpu_delta(before_sibling, after_sibling)
    accepted = (timed["nonidle_total"] >= 1 and sibling["total"] >= 1 and
                sibling["nonidle_total"] == 0)
    return {
        "policy": {
            "source": "/proc/stat", "timed_min_nonidle_jiffies": 1,
            "sibling_min_total_jiffies": 1, "sibling_max_nonidle_jiffies": 0,
        },
        "accepted": accepted,
        "timed_before": before_cpu, "timed_after": after_cpu, "timed_delta": timed,
        "sibling_before": before_sibling, "sibling_after": after_sibling,
        "sibling_delta": sibling,
    }
