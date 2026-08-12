#!/usr/bin/env python3
"""Reproducible Leopard2 / Leopard-main / Wirehair performance atlas.

The campaign is deliberately single-process and single-core.  It owns the
project's canonical benchmark lock, freezes all three executables, writes one
JSON record per invocation, and can resume after interruption without timing
an executable that has changed underneath it.  Plotting uses only the Python
standard library so the checked-in SVGs can be regenerated on a minimal host.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import datetime as dt
import fcntl
import hashlib
import html
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence


MANIFEST_SCHEMA = "leopard2-performance-atlas-manifest/v1"
SUMMARY_SCHEMA = "leopard2-performance-atlas-summary/v1"
RAW_BUNDLE_SCHEMA = "leopard2-performance-atlas-raw-bundle/v1"
RUNNER_VERSION = 1
R_COUNT = 32
K_MAX = 224
BLOCK_BYTES = (64, 1024, 4096, 1024 * 1024)
CODECS = ("leopard2", "leopard1", "wirehair")
CODEC_LABELS = {
    "leopard2": "Leopard2",
    "leopard1": "Leopard main",
    "wirehair": "Wirehair (shipping)",
}
CODEC_COLORS = {
    "leopard2": "#cf222e",
    "leopard1": "#0969da",
    "wirehair": "#1a7f37",
}
LOSS_LABELS = ("one", "two", "ten_percent", "full")
LOSS_DISPLAY = {
    "one": "1 random source erasure",
    "two": "2 random source erasures",
    "ten_percent": "10% random source erasures",
    "full": "maximum random source erasures",
}
CANONICAL_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
DEFAULT_AS_LIMIT = 2 * 1024 * 1024 * 1024
EXPECTED_MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXPECTED_WIREHAIR_COMMIT = "067ca7cdb66aed424ec23f97557429bf791c6f0c"
EXPECTED_WIREHAIR_PROFILE_ID = 0x4D241359DB07BB07
MAIN_ARCHIVE_SOURCE_RELATIVE = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)
WIREHAIR_ARCHIVE_SOURCE_RELATIVE = {
    "wirehair.cpp": "wirehair.cpp",
    "gf256.cpp": "gf256.cpp",
    "WirehairCodec.cpp": "WirehairCodec.cpp",
    "WirehairTools.cpp": "WirehairTools.cpp",
    "WirehairV2GF16.cpp": "codec/WirehairV2GF16.cpp",
    "WirehairV2Codec.cpp": "codec/WirehairV2Codec.cpp",
    "WirehairV2Peel.cpp": "codec/WirehairV2Peel.cpp",
    "WirehairV2Plan.cpp": "codec/WirehairV2Plan.cpp",
    "WirehairV2Policy.cpp": "codec/WirehairV2Policy.cpp",
    "WirehairV2Precode.cpp": "codec/WirehairV2Precode.cpp",
    "WirehairV2PrecodeDecode.cpp": "codec/WirehairV2PrecodeDecode.cpp",
    "WirehairV2PrecodeEncode.cpp": "codec/WirehairV2PrecodeEncode.cpp",
    "WirehairV2Profile.cpp": "codec/WirehairV2Profile.cpp",
    "WirehairV2Seeds.cpp": "codec/WirehairV2Seeds.cpp",
    "WirehairV2Solve.cpp": "codec/WirehairV2Solve.cpp",
}


class AtlasError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AtlasError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(str(path.parent), os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_bytes(path, (json.dumps(value, indent=2, sort_keys=True) +
                              "\n").encode("utf-8"))


def read_json(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as source:
            return json.load(source)
    except (OSError, json.JSONDecodeError) as error:
        raise AtlasError(f"cannot read valid JSON {path}: {error}") from error


def git_output(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *arguments], check=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if completed.returncode != 0:
        raise AtlasError(
            f"git {' '.join(arguments)} failed in {root}: "
            f"{completed.stderr.strip()}")
    return completed.stdout.strip()


def source_identity(root: Path, expected_commit: str | None = None) -> dict[str, Any]:
    root = root.resolve()
    commit = git_output(root, "rev-parse", "HEAD")
    if expected_commit is not None:
        require(commit == expected_commit,
                f"{root} is {commit}, expected {expected_commit}")
    status = git_output(root, "status", "--porcelain=v1", "--untracked-files=all")
    require(not status, f"tracked source tree is dirty: {root}")
    tree = git_output(root, "rev-parse", "HEAD^{tree}")
    return {"root": str(root), "commit": commit, "tree": tree,
            "tracked_clean": True}


def executable_identity(path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    stat = path.stat()
    require(path.is_file(), f"not a regular executable: {path}")
    require(os.access(path, os.X_OK), f"not executable: {path}")
    return {
        "path": str(path),
        "size": stat.st_size,
        "sha256": sha256_file(path),
        "mode": stat.st_mode & 0o7777,
    }


def data_file_identity(path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    stat = path.stat()
    require(path.is_file(), f"not a regular file: {path}")
    return {"path": str(path), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": sha256_file(path), "mode": stat.st_mode & 0o7777}


def verify_identity(identity: Mapping[str, Any]) -> None:
    current = executable_identity(Path(str(identity["path"])))
    for key in ("size", "sha256"):
        require(current[key] == identity[key],
                f"executable identity changed for {identity['path']}: {key}")


def load_build_provenance_module(source_root: Path) -> Any:
    path = source_root / "tools" / "leopard2_build_provenance.py"
    require(path.is_file(), f"Leopard2 provenance module is missing: {path}")
    name = "leopard2_atlas_build_provenance"
    specification = importlib.util.spec_from_file_location(name, path)
    require(specification is not None and specification.loader is not None,
            "cannot load Leopard2 build-provenance module")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def capture_leopard2_build_closure(build_root: Path, source_root: Path,
                                   executable: Path) -> dict[str, Any]:
    module = load_build_provenance_module(source_root.resolve())
    try:
        closure = module.candidate_build_provenance(
            build_root.resolve(), source_root.resolve(), executable.resolve(),
            "bench_leopard2_allk")
    except Exception as error:
        raise AtlasError(f"Leopard2 production build closure failed: {error}") from error
    require(isinstance(closure, dict) and closure,
            "Leopard2 production build closure is empty")
    return closure


def compile_tokens(row: Mapping[str, Any]) -> list[str]:
    arguments = row.get("arguments")
    if isinstance(arguments, list) and all(isinstance(item, str)
                                           for item in arguments):
        return list(arguments)
    command = row.get("command")
    require(isinstance(command, str), "compile-command row lacks arguments")
    return shlex.split(command)


def cmake_cache_entries(path: Path) -> dict[str, str]:
    cache: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        key_and_type, value = line.split("=", 1)
        cache[key_and_type.split(":", 1)[0]] = value
    return cache


def resolve_command_path(value: str, directory: Path,
                         description: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute() or "/" in value:
        if not candidate.is_absolute():
            candidate = directory / candidate
        try:
            return candidate.resolve(strict=True)
        except OSError as error:
            raise AtlasError(
                f"{description} does not resolve to a file: {value}: "
                f"{error}") from error
    found = shutil.which(value)
    require(found is not None, f"cannot resolve {description}: {value}")
    return Path(found).resolve(strict=True)


def compile_row_paths(row: Mapping[str, Any], output_root: Path,
                      label: str) -> tuple[Path, Path]:
    directory_value = row.get("directory")
    require(isinstance(directory_value, str),
            f"{label} compile-command row lacks a directory")
    directory = Path(directory_value).resolve(strict=True)
    source_value = row.get("file")
    output_value = row.get("output")
    require(isinstance(source_value, str) and isinstance(output_value, str),
            f"{label} compile-command row lacks source/output identity")
    source = Path(source_value)
    output = Path(output_value)
    if not source.is_absolute():
        source = directory / source
    if not output.is_absolute():
        # CMake's compile database reports ``output`` relative to the
        # top-level build tree even when ``directory`` is a nested binary
        # directory.  The command's -o operand remains relative to
        # ``directory``.  Bind both representations to the same object so a
        # plausible but stale output field cannot attest a different file.
        output = output_root / output
    resolved_output = output.resolve(strict=True)
    tokens = compile_tokens(row)
    output_operands = [tokens[index + 1]
                       for index, token in enumerate(tokens[:-1])
                       if token == "-o"]
    require(len(output_operands) == 1,
            f"{label} compile command must contain exactly one -o operand")
    command_output = Path(output_operands[0])
    if not command_output.is_absolute():
        command_output = directory / command_output
    require(command_output.resolve(strict=True) == resolved_output,
            f"{label} compile-command output field differs from -o operand")
    return source.resolve(strict=True), resolved_output


def recipe_object_paths(recipe: str, working_directory: Path,
                        label: str) -> list[Path]:
    objects: list[Path] = []
    for token in shlex.split(recipe):
        if not token.endswith((".o", ".obj")):
            continue
        path = Path(token)
        if not path.is_absolute():
            path = working_directory / path
        try:
            objects.append(path.resolve(strict=True))
        except OSError as error:
            raise AtlasError(
                f"{label} recipe references a missing object {token}: "
                f"{error}") from error
    return objects


def run_binary_capture(command: Sequence[str], label: str) -> bytes:
    try:
        completed = subprocess.run(
            list(command), check=False, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise AtlasError(f"{label} failed: {error}") from error
    require(completed.returncode == 0,
            f"{label} failed with {completed.returncode}: "
            f"{completed.stderr.decode('utf-8', errors='replace')[-2000:]}")
    return completed.stdout


def capture_external_avx2_compile_build(
        build_root: Path, executable: Path, archive_relative: Path,
        archive_link_relative: Path,
        archive_sources: Mapping[str, Path],
        executable_sources: Mapping[str, Path],
        cache_contract: Mapping[str, str], archive_target_fragment: str,
        executable_target_fragment: str, label: str) -> dict[str, Any]:
    build = build_root.resolve(strict=True)
    expected_executable = (build / executable.name).resolve(strict=True)
    require(executable.resolve(strict=True) == expected_executable,
            f"{label} executable is not the declared build target")
    commands_path = build / "compile_commands.json"
    cache_path = build / "CMakeCache.txt"
    executable_link_path = build / f"CMakeFiles/{executable.name}.dir/link.txt"
    archive_link_path = build / archive_link_relative
    archive_path = build / archive_relative
    for path in (commands_path, cache_path, executable_link_path,
                 archive_link_path, archive_path):
        require(path.is_file(), f"{label} build evidence is missing: {path}")
    rows = read_json(commands_path)
    require(isinstance(rows, list) and rows,
            f"{label} compile-command graph is empty")
    expected_by_role = {
        "archive": {name: path.resolve(strict=True)
                    for name, path in archive_sources.items()},
        "executable": {name: path.resolve(strict=True)
                       for name, path in executable_sources.items()},
    }
    require(set(expected_by_role["archive"]).isdisjoint(
            expected_by_role["executable"]),
            f"{label} source basenames are ambiguous across targets")
    observed: dict[str, dict[str, dict[str, Any]]] = {
        "archive": {}, "executable": {}}
    required_flags = {
        "-march=x86-64", "-mtune=generic", "-mavx2",
        "-mno-avx512f",
    }
    forbidden_prefixes = ("-march=native", "-mtune=native", "-mavx512")
    for row in rows:
        require(isinstance(row, dict) and isinstance(row.get("file"), str),
                f"{label} compile-command row is malformed")
        output = row.get("output")
        require(isinstance(output, str),
                f"{label} compile-command row lacks an output identity")
        normalized_output = output.replace("\\", "/")
        if archive_target_fragment in normalized_output:
            role = "archive"
        elif executable_target_fragment in normalized_output:
            role = "executable"
        else:
            continue
        source_path, object_path = compile_row_paths(row, build, label)
        source_name = source_path.name
        expected_source = expected_by_role[role].get(source_name)
        require(expected_source is not None,
                f"{label} compile graph has an unexpected {role} source: "
                f"{source_path}")
        require(source_path == expected_source,
                f"{label} {source_name} resolves to {source_path}, expected "
                f"the declared source {expected_source}")
        require(source_name not in observed[role],
                f"{label} compile graph duplicates {role} source "
                f"{source_name}")
        tokens = compile_tokens(row)
        missing = required_flags - set(tokens)
        require(not missing,
                f"{label} {row['file']} lacks the required global AVX2 flags: "
                f"{sorted(missing)}")
        require(not any(token.startswith(forbidden_prefixes) and
                        not token.startswith("-mno-avx512") for token in tokens),
                f"{label} {row['file']} enables a forbidden ISA")
        source_identity_value = data_file_identity(source_path)
        object_identity_value = data_file_identity(object_path)
        require(object_identity_value["mtime_ns"] >=
                source_identity_value["mtime_ns"],
                f"{label} object predates its source: {source_path}")
        compile_driver = resolve_command_path(
            tokens[0], Path(str(row["directory"])),
            f"{label} compile driver for {source_name}")
        observed[role][source_name] = {
            "source": source_identity_value,
            "object": object_identity_value,
            "compile_driver": executable_identity(compile_driver),
            "compile_command_sha256": sha256_bytes(
                canonical_json(tokens).encode("utf-8")),
        }
    for role in ("archive", "executable"):
        require(set(observed[role]) == set(expected_by_role[role]),
                f"{label} {role} compile graph source set differs: missing="
                f"{sorted(set(expected_by_role[role]) - set(observed[role]))}, "
                f"extra={sorted(set(observed[role]) - set(expected_by_role[role]))}")

    cache = cmake_cache_entries(cache_path)
    for key, value in cache_contract.items():
        require(cache.get(key) == value,
                f"{label} cache {key} differs: {cache.get(key)!r}")
    compiler_value = cache.get("CMAKE_CXX_COMPILER")
    archiver_value = cache.get("CMAKE_AR")
    ranlib_value = cache.get("CMAKE_RANLIB")
    require(bool(compiler_value) and bool(archiver_value) and bool(ranlib_value),
            f"{label} cache omits compiler/archive tool identities")
    compiler_path = resolve_command_path(
        str(compiler_value), build, f"{label} C++ compiler")
    archiver_path = resolve_command_path(
        str(archiver_value), build, f"{label} archiver")
    ranlib_path = resolve_command_path(
        str(ranlib_value), build, f"{label} ranlib")
    for role_records in observed.values():
        for source_name, record in role_records.items():
            require(record["compile_driver"]["path"] == str(compiler_path),
                    f"{label} {source_name} was compiled by "
                    f"{record['compile_driver']['path']}, not cached compiler "
                    f"{compiler_path}")

    archive_link_text = archive_link_path.read_text(encoding="utf-8")
    require(Path(str(archiver_value)).name in archive_link_text and
            Path(str(ranlib_value)).name in archive_link_text,
            f"{label} archive recipe omits cached archive tools")
    archive_recipe_objects = recipe_object_paths(
        archive_link_text, archive_path.parent, f"{label} archive")
    expected_archive_objects = [
        Path(record["object"]["path"])
        for record in observed["archive"].values()]
    require(len(archive_recipe_objects) == len(expected_archive_objects) and
            set(archive_recipe_objects) == set(expected_archive_objects),
            f"{label} archive recipe object graph differs")

    executable_link_text = executable_link_path.read_text(encoding="utf-8")
    require(Path(str(compiler_value)).name in executable_link_text,
            f"{label} executable recipe omits cached compiler driver")
    executable_recipe_objects = recipe_object_paths(
        executable_link_text, expected_executable.parent,
        f"{label} executable")
    expected_executable_objects = [
        Path(record["object"]["path"])
        for record in observed["executable"].values()]
    require(len(executable_recipe_objects) == len(expected_executable_objects)
            and set(executable_recipe_objects) ==
            set(expected_executable_objects),
            f"{label} executable recipe object graph differs")
    require(str(archive_path) in executable_link_text or
            archive_path.name in executable_link_text,
            f"{label} executable link recipe omits its codec archive")

    member_output = run_binary_capture(
        [str(archiver_path), "t", str(archive_path)],
        f"{label} archive member listing")
    members = [line for line in member_output.decode(
        "utf-8", errors="strict").splitlines() if line]
    require(len(members) == len(set(members)),
            f"{label} archive contains duplicate member names")
    expected_members = {path.name: path for path in expected_archive_objects}
    require(len(expected_members) == len(expected_archive_objects),
            f"{label} archive object basenames are ambiguous")
    require(set(members) == set(expected_members),
            f"{label} archive member set differs: missing="
            f"{sorted(set(expected_members) - set(members))}, extra="
            f"{sorted(set(members) - set(expected_members))}")
    archive_members: list[dict[str, Any]] = []
    for member in members:
        member_bytes = run_binary_capture(
            [str(archiver_path), "p", str(archive_path), member],
            f"{label} archive member {member}")
        object_path = expected_members[member]
        object_hash = sha256_file(object_path)
        member_hash = sha256_bytes(member_bytes)
        require(member_hash == object_hash,
                f"{label} archive member {member} differs from its "
                f"compiled object")
        archive_members.append({
            "name": member, "size": len(member_bytes),
            "sha256": member_hash,
            "object_path": str(object_path),
        })
    archive_identity_value = data_file_identity(archive_path)
    require(all(archive_identity_value["mtime_ns"] >=
                record["object"]["mtime_ns"]
                for record in observed["archive"].values()),
            f"{label} archive predates one or more compiled objects")
    executable_identity_value = executable_identity(executable)
    executable_mtime = executable.resolve(strict=True).stat().st_mtime_ns
    require(executable_mtime >= archive_identity_value["mtime_ns"] and
            all(executable_mtime >= record["object"]["mtime_ns"]
                for record in observed["executable"].values()),
            f"{label} executable predates one or more link inputs")
    return {
        "build_root": str(build),
        "compile_commands": data_file_identity(commands_path),
        "cmake_cache": data_file_identity(cache_path),
        "archive_link_recipe": data_file_identity(archive_link_path),
        "executable_link_recipe": data_file_identity(executable_link_path),
        "archive": archive_identity_value,
        "executable": executable_identity_value,
        "source_object_graph": observed,
        "archive_members": archive_members,
        "tools": {
            "cxx_compiler": executable_identity(compiler_path),
            "archiver": executable_identity(archiver_path),
            "ranlib": executable_identity(ranlib_path),
        },
        "required_flags": sorted(required_flags),
        "cache_contract": dict(cache_contract),
    }


def freeze_executable(source: Path, destination: Path,
                      required_sha256: str) -> dict[str, Any]:
    source_identity_value = executable_identity(source)
    require(source_identity_value["sha256"] == required_sha256,
            f"SHA-256 mismatch for {source}: got "
            f"{source_identity_value['sha256']}, expected {required_sha256}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        existing = executable_identity(destination)
        require(existing["sha256"] == required_sha256,
                f"frozen executable collision: {destination}")
    else:
        temporary = destination.with_name(destination.name + ".tmp")
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()
        shutil.copyfile(source, temporary)
        os.chmod(temporary, 0o555)
        require(sha256_file(temporary) == required_sha256,
                f"copy changed executable bytes: {source}")
        os.replace(temporary, destination)
    frozen = executable_identity(destination)
    require(frozen["mode"] & 0o222 == 0,
            f"frozen executable remains writable: {destination}")
    return frozen


def k_values() -> tuple[int, ...]:
    values = set(range(1, K_MAX + 1, 2))
    power = 1
    while power <= K_MAX:
        values.add(power)
        power *= 2
    values.add(K_MAX)
    result = tuple(sorted(values))
    require(len(result) == 120, "unexpected all-K sample count")
    return result


def conceptual_loss(k: int, label: str) -> int | None:
    if label == "one":
        return 1
    if label == "two":
        return 2 if k >= 2 else None
    if label == "ten_percent":
        return min(R_COUNT, max(1, (k + 9) // 10))
    if label == "full":
        return min(k, R_COUNT)
    raise AtlasError(f"unknown loss label: {label}")


def loss_map(k: int) -> dict[int, list[str]]:
    by_count: dict[int, list[str]] = {}
    for label in LOSS_LABELS:
        count = conceptual_loss(k, label)
        if count is not None:
            by_count.setdefault(count, []).append(label)
    return by_count


def seed_for(k: int, block_bytes: int, losses: int) -> int:
    # Keep a K/loss pattern identical across shard-size panels.  Source bytes
    # still differ because each benchmark expands this seed across a different
    # total message length.
    del block_bytes
    material = f"leopard2-atlas-v1:{k}:{R_COUNT}:{losses}".encode()
    value = int.from_bytes(hashlib.sha256(material).digest()[:8], "big")
    return value or 1


def repetition_policy(k: int, block_bytes: int) -> tuple[int, int, int]:
    working_bytes = max(k, R_COUNT) * block_bytes
    target = 8 * 1024 * 1024
    reuse = max(1, min(4096, (target + working_bytes - 1) // working_bytes))
    if block_bytes >= 1024 * 1024:
        return reuse, 5, 1
    return reuse, 9, 2


def available_codecs(k: int) -> tuple[str, ...]:
    result = ["leopard2"]
    if k >= R_COUNT:
        result.append("leopard1")
    if k >= 2:
        result.append("wirehair")
    return tuple(result)


def ceil_power_of_two(value: int) -> int:
    require(isinstance(value, int) and not isinstance(value, bool) and
            value > 0, "power-of-two input is invalid")
    return 1 << (value - 1).bit_length()


def expected_geometry(k: int, profile: str) -> tuple[int, int]:
    """Return (parent_count, padded_side) for the frozen atlas profile."""
    if profile == "low_v1":
        padded = ceil_power_of_two(k)
        return ceil_power_of_two(padded + R_COUNT), padded
    require(profile == "legacy_high_v1", "unexpected atlas profile")
    padded = ceil_power_of_two(R_COUNT)
    return ceil_power_of_two(k + padded), padded


def build_manifest() -> dict[str, Any]:
    cells: list[dict[str, Any]] = []
    for block_bytes in BLOCK_BYTES:
        for k in k_values():
            reuse, iterations, warmup = repetition_policy(k, block_bytes)
            for losses, labels in sorted(loss_map(k).items()):
                cell_id = f"k{k:03d}_b{block_bytes:07d}_l{losses:02d}"
                cells.append({
                    "id": cell_id,
                    "K": k,
                    "R": R_COUNT,
                    "shard_bytes": block_bytes,
                    "loss_count": losses,
                    "loss_labels": labels,
                    "seed": seed_for(k, block_bytes, losses),
                    "reuse": reuse,
                    "iterations": iterations,
                    "warmup": warmup,
                    "available_codecs": list(available_codecs(k)),
                })
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "runner_version": RUNNER_VERSION,
        "matrix": {
            "K_policy": "every odd K, every power of two, and endpoint 224",
            "K_values": list(k_values()),
            "R": R_COUNT,
            "shard_bytes": list(BLOCK_BYTES),
            "loss_labels": dict(LOSS_DISPLAY),
            "loss_formulas": {
                "one": "1",
                "two": "2 when K >= 2",
                "ten_percent": "min(R,max(1,ceil(K/10)))",
                "full": "min(K,R)",
            },
            "missing_coordinate_class": "original/source shards only",
            "unique_cell_count": len(cells),
        },
        "cells": cells,
    }
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == MANIFEST_SCHEMA,
            "unexpected manifest schema")
    matrix = manifest.get("matrix")
    cells = manifest.get("cells")
    require(isinstance(matrix, dict) and isinstance(cells, list),
            "manifest lacks matrix or cells")
    require(matrix.get("K_values") == list(k_values()),
            "manifest K values drifted")
    require(matrix.get("R") == R_COUNT,
            "manifest R drifted")
    require(matrix.get("shard_bytes") == list(BLOCK_BYTES),
            "manifest block sizes drifted")
    require(matrix.get("unique_cell_count") == len(cells),
            "manifest cell count is inconsistent")
    expected = build_expected_cells_without_recursion()
    observed = {cell.get("id"): cell for cell in cells
                if isinstance(cell, dict)}
    require(len(observed) == len(cells), "duplicate or malformed cell IDs")
    require(observed == expected, "manifest cells differ from canonical matrix")


def build_expected_cells_without_recursion() -> dict[str, dict[str, Any]]:
    expected: dict[str, dict[str, Any]] = {}
    for block_bytes in BLOCK_BYTES:
        for k in k_values():
            reuse, iterations, warmup = repetition_policy(k, block_bytes)
            for losses, labels in sorted(loss_map(k).items()):
                cell_id = f"k{k:03d}_b{block_bytes:07d}_l{losses:02d}"
                expected[cell_id] = {
                    "id": cell_id, "K": k, "R": R_COUNT,
                    "shard_bytes": block_bytes, "loss_count": losses,
                    "loss_labels": labels,
                    "seed": seed_for(k, block_bytes, losses),
                    "reuse": reuse, "iterations": iterations,
                    "warmup": warmup,
                    "available_codecs": list(available_codecs(k)),
                }
    return expected


def allowed_cpus() -> list[int]:
    if hasattr(os, "sched_getaffinity"):
        cpus = sorted(os.sched_getaffinity(0))
    else:
        cpus = list(range(os.cpu_count() or 1))
    require(bool(cpus), "process has no allowed CPUs")
    return cpus


def physical_cpu(cpus: Sequence[int]) -> int:
    seen: set[tuple[str, str]] = set()
    for cpu in cpus:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        try:
            package = (topology / "physical_package_id").read_text().strip()
            core = (topology / "core_id").read_text().strip()
        except OSError:
            return cpu
        key = package, core
        if key not in seen:
            return cpu
        seen.add(key)
    return cpus[0]


@contextlib.contextmanager
def canonical_lock() -> Iterator[dict[str, Any]]:
    CANONICAL_LOCK.parent.mkdir(parents=True, exist_ok=True)
    with CANONICAL_LOCK.open("a+b") as lock:
        print(f"waiting for canonical lock {CANONICAL_LOCK}", flush=True)
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        before = os.fstat(lock.fileno())
        identity = {"path": str(CANONICAL_LOCK), "device": before.st_dev,
                    "inode": before.st_ino}
        print("acquired canonical lock", flush=True)
        try:
            yield identity
        finally:
            after = os.fstat(lock.fileno())
            require((after.st_dev, after.st_ino) ==
                    (before.st_dev, before.st_ino),
                    "canonical lock identity changed")
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def codec_order(cell_id: str) -> tuple[str, ...]:
    # Filter first so two-codec low-rate cells remain balanced rather than
    # inheriting a biased projection of a three-codec rotation.
    cell_k = int(cell_id[1:4])
    available = available_codecs(cell_k)
    parts = cell_id.split("_")
    require(len(parts) == 3 and parts[1].startswith("b") and
            parts[2].startswith("l"), "cell ID cannot drive codec rotation")
    block_index = BLOCK_BYTES.index(int(parts[1][1:]))
    losses = int(parts[2][1:])
    offset = (cell_k + block_index + losses) % len(available)
    return available[offset:] + available[:offset]


def command_for(codec: str, executable: Path,
                cell: Mapping[str, Any]) -> list[str]:
    common = [
        str(executable), "--k", str(cell["K"]),
        "--r", str(cell["R"]), "--bytes", str(cell["shard_bytes"]),
        "--loss", str(cell["loss_count"]), "--reuse", str(cell["reuse"]),
        "--iterations", str(cell["iterations"]),
        "--warmup", str(cell["warmup"]), "--seed", str(cell["seed"]),
        "--json", "-",
    ]
    if codec == "leopard2":
        return common + [
            "--batch", "1", "--threads", "1", "--profile", "auto",
            "--field", "auto", "--backend", "auto", "--skip-legacy",
            "--retain-samples", "--measure-one-shot-encode",
            "--measure-one-shot-decode",
        ]
    if codec == "leopard1":
        return common + ["--batch", "1", "--threads", "1"]
    if codec == "wirehair":
        return common
    raise AtlasError(f"unknown codec: {codec}")


def run_invocation(codec: str, executable: Path, cell: Mapping[str, Any],
                   cpu: int, as_limit: int, timeout: int) -> tuple[dict[str, Any], str]:
    command = ["/usr/bin/taskset", "-c", str(cpu), "/usr/bin/prlimit",
               f"--as={as_limit}", "--"] + command_for(codec, executable, cell)
    environment = dict(os.environ)
    environment.update({
        "OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "true",
        "OMP_PLACES": "cores", "MALLOC_ARENA_MAX": "1",
    })
    completed = subprocess.run(
        command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=environment, timeout=timeout)
    if completed.returncode != 0:
        raise AtlasError(
            f"{codec} failed for {cell['id']} with {completed.returncode}: "
            f"{completed.stderr[-2000:]}")
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise AtlasError(
            f"{codec} emitted invalid JSON for {cell['id']}: {error}") from error
    validate_payload(codec, payload, cell)
    return payload, completed.stderr


def attest_leopard2(executable: Path, cpu: int, as_limit: int,
                    timeout: int, expected_source: Mapping[str, Any]) -> dict[str, Any]:
    command = [
        "/usr/bin/taskset", "-c", str(cpu), "/usr/bin/prlimit",
        f"--as={as_limit}", "--", str(executable),
        "--k", "32", "--r", "32", "--bytes", "64", "--loss", "1",
        "--batch", "1", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--threads", "1", "--seed", "1",
        "--profile", "auto", "--field", "auto", "--backend", "auto",
        "--skip-legacy", "--retain-samples", "--attest-source", "--json", "-",
    ]
    completed = subprocess.run(
        command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, timeout=timeout,
        env={**os.environ, "OMP_NUM_THREADS": "1", "MALLOC_ARENA_MAX": "1"})
    require(completed.returncode == 0,
            "Leopard2 source attestation failed: " + completed.stderr[-2000:])
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise AtlasError(f"Leopard2 attestation emitted invalid JSON: {error}") from error
    require(payload.get("schema") == "leopard2-benchmark-v5",
            "Leopard2 attestation schema differs")
    build = payload.get("build")
    require(isinstance(build, dict), "Leopard2 attestation lacks build identity")
    require(build.get("source_commit") == expected_source["commit"],
            "Leopard2 executable commit does not match the declared source")
    require(build.get("source_tree") == expected_source["tree"],
            "Leopard2 executable tree does not match the declared source")
    require(build.get("source_tracked_dirty") is False,
            "Leopard2 executable was built from tracked-dirty source")
    require(payload.get("correctness", {}).get("leopard2_round_trip") is True,
            "Leopard2 attestation preflight did not round-trip")
    return payload


def finite_positive(value: Any, name: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{name} is not numeric")
    converted = float(value)
    require(math.isfinite(converted) and converted > 0.0,
            f"{name} is not finite and positive")
    return converted


def finite_nonnegative(value: Any, name: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{name} is not numeric")
    converted = float(value)
    require(math.isfinite(converted) and converted >= 0.0,
            f"{name} is not finite and nonnegative")
    return converted


def validate_metric(metric: Any, name: str, iterations: int,
                    median_key: str) -> None:
    require(isinstance(metric, dict), f"{name} is not an object")
    finite_positive(metric.get(median_key), f"{name}.{median_key}")
    sample_key = ("samples_us_per_batch_call" if
                  median_key.endswith("per_batch_call") else "samples_us")
    samples = metric.get(sample_key)
    require(isinstance(samples, list) and len(samples) == iterations,
            f"{name}.{sample_key} has the wrong length")
    sample_values = [finite_positive(value, f"{name}.{sample_key}[{index}]")
                     for index, value in enumerate(samples)]
    median = finite_positive(metric.get(median_key), f"{name}.{median_key}")
    expected_median = statistics.median(sample_values)
    require(math.isclose(median, expected_median, rel_tol=1e-9,
                         abs_tol=0.000002),
            f"{name}.{median_key} differs from retained samples")
    min_key = ("minimum_us_per_batch_call" if
               median_key.endswith("per_batch_call") else "minimum_us")
    max_key = ("maximum_us_per_batch_call" if
               median_key.endswith("per_batch_call") else "maximum_us")
    mad_key = ("mad_us_per_batch_call" if
               median_key.endswith("per_batch_call") else "mad_us")
    if min_key in metric or max_key in metric or mad_key in metric:
        require(math.isclose(finite_positive(metric.get(min_key),
                                               f"{name}.{min_key}"),
                             min(sample_values),
                             rel_tol=1e-9, abs_tol=0.000002),
                f"{name}.{min_key} differs from retained samples")
        require(math.isclose(finite_positive(metric.get(max_key),
                                               f"{name}.{max_key}"),
                             max(sample_values),
                             rel_tol=1e-9, abs_tol=0.000002),
                f"{name}.{max_key} differs from retained samples")
        expected_mad = statistics.median(
            [abs(value - expected_median) for value in sample_values])
        require(math.isclose(finite_nonnegative(metric.get(mad_key),
                                                  f"{name}.{mad_key}"),
                             expected_mad,
                             rel_tol=1e-9, abs_tol=0.000002),
                f"{name}.{mad_key} differs from retained samples")


def validate_payload(codec: str, payload: Any,
                     cell: Mapping[str, Any]) -> None:
    require(isinstance(payload, dict), f"{codec} result is not an object")
    parameters = payload.get("parameters")
    require(isinstance(parameters, dict), f"{codec} lacks parameters")
    exact = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["shard_bytes"],
        "loss_count": cell["loss_count"], "reuse": cell["reuse"],
        "iterations": cell["iterations"], "warmup": cell["warmup"],
        "seed": cell["seed"],
    }
    for key, value in exact.items():
        require(parameters.get(key) == value,
                f"{codec} parameter {key} differs for {cell['id']}")
    if codec in {"leopard2", "leopard1", "wirehair"}:
        require(parameters.get("batch") == 1 and
                parameters.get("thread_count") == 1,
                f"{codec} atlas invocation is not one stripe/one thread")
    missing = parameters.get("missing_original_indices")
    require(isinstance(missing, list) and len(missing) == cell["loss_count"],
            f"{codec} missing-index count differs for {cell['id']}")
    require(missing == sorted(set(missing)),
            f"{codec} missing indices are not sorted and unique")
    require(all(isinstance(index, int) and 0 <= index < cell["K"]
                for index in missing), f"{codec} missing index out of range")
    digests = payload.get("workload_digests")
    require(isinstance(digests, dict) and
            digests.get("algorithm") == "fnv1a64",
            f"{codec} lacks FNV workload digests")
    digest_keys = ["original_data", "recovered_originals"]
    digest_keys.append("generated_repair" if codec == "wirehair" else
                       "transmitted_parity")
    for key in digest_keys:
        value = digests.get(key)
        require(isinstance(value, str) and len(value) == 16 and
                all(character in "0123456789abcdef" for character in value),
                f"{codec} has invalid {key} digest")
    metrics = payload.get("metrics")
    require(isinstance(metrics, dict), f"{codec} lacks metrics")
    iterations = int(cell["iterations"])
    if codec == "leopard2":
        require(payload.get("schema") == "leopard2-benchmark-v9",
                "unexpected Leopard2 benchmark schema")
        correctness = payload.get("correctness")
        require(isinstance(correctness, dict) and
                correctness.get("leopard2_round_trip") is True,
                "Leopard2 round trip failed")
        resolved = payload.get("resolved")
        require(isinstance(resolved, dict), "Leopard2 lacks resolution")
        expected_profile = "low_v1" if cell["K"] < R_COUNT else "legacy_high_v1"
        expected_parent, expected_padded = expected_geometry(
            int(cell["K"]), expected_profile)
        require(parameters.get("requested_profile") == "auto" and
                parameters.get("requested_field") == "auto" and
                parameters.get("requested_backend") == "auto" and
                all(parameters.get(name) is False for name in (
                    "force_generic_decode", "force_specialized_decode",
                    "force_tiled_decode", "force_materialized_decode")) and
                parameters.get("skip_legacy") is True and
                parameters.get("retain_samples") is True and
                parameters.get("measure_one_shot_encode") is True and
                parameters.get("measure_one_shot_decode") is True,
                "Leopard2 request/timing policy differs from the atlas")
        require(resolved.get("profile") == expected_profile,
                f"Leopard2 AUTO profile differs for {cell['id']}")
        require(resolved.get("field") == "gf8", "atlas must remain in GF8")
        require(resolved.get("backend") == "avx2",
                "atlas requires the ISA-matched AVX2 Leopard2 backend")
        require(resolved.get("thread_count") == 1 and
                resolved.get("parent_count") == expected_parent and
                resolved.get("padded_side") == expected_padded,
                f"Leopard2 resolved geometry differs for {cell['id']}")
        require(correctness.get("legacy_comparison") is None,
                "Leopard2 atlas unexpectedly timed/compared its legacy path")
        memory = payload.get("memory")
        require(isinstance(memory, dict) and
                memory.get("encode_scratch_bytes_batch") ==
                    memory.get("encode_scratch_bytes_per_stripe") and
                memory.get("decode_scratch_bytes_batch") ==
                    memory.get("decode_scratch_bytes_per_stripe"),
                "Leopard2 batch-one scratch accounting differs")
        for name in ("encode_execution", "one_shot_encode",
                     "decode_execution", "one_shot_decode_including_setup"):
            validate_metric(metrics.get(name), f"leopard2.{name}", iterations,
                            "median_us_per_batch_call")
        for name in ("codec_setup", "decode_plan_setup"):
            validate_metric(metrics.get(name), f"leopard2.{name}", iterations,
                            "median_us")
        require(metrics.get("rate_semantics") ==
                "offered_received counts all non-null shard pointers supplied; "
                "a plan may read a deterministic subset",
                "Leopard2 rate semantics changed")
    elif codec == "leopard1":
        require(payload.get("schema") == "leopard-main-benchmark-v1",
                "unexpected Leopard-main benchmark schema")
        correctness = payload.get("correctness")
        require(isinstance(correctness, dict) and
                correctness.get("round_trip") is True,
                "Leopard-main round trip failed")
        require(payload.get("build", {}).get("main_source_commit") ==
                EXPECTED_MAIN_COMMIT, "Leopard-main source commit differs")
        require(payload.get("build", {}).get("pure_avx2") is True,
                "Leopard-main adapter is not attested as pure AVX2")
        expected_parent, expected_padded = expected_geometry(
            int(cell["K"]), "legacy_high_v1")
        resolved = payload.get("resolved")
        require(parameters.get("logical_shard_bytes") == cell["shard_bytes"] and
                isinstance(resolved, dict) and
                resolved.get("profile") == "legacy_high_v1" and
                resolved.get("field") == "gf8" and
                resolved.get("thread_count") == 1 and
                resolved.get("parent_count") == expected_parent and
                resolved.get("padded_side") == expected_padded and
                resolved.get("padded_application_bytes") is False and
                resolved.get("padding_policy") == "zero suffix per shard",
                f"Leopard-main geometry/timing shape differs for {cell['id']}")
        encode_work = 2 * expected_padded
        memory = payload.get("memory")
        require(isinstance(memory, dict) and
                memory.get("alignment") == 64 and
                memory.get("encode_work_shards_per_stripe") == encode_work and
                memory.get("decode_work_shards_per_stripe") == expected_parent and
                memory.get("encode_work_bytes_batch") ==
                    encode_work * cell["shard_bytes"] and
                memory.get("decode_work_bytes_batch") ==
                    expected_parent * cell["shard_bytes"],
                "Leopard-main work geometry differs")
        validate_metric(metrics.get("encode_execution"),
                        "leopard1.encode_execution", iterations,
                        "median_us_per_batch_call")
        validate_metric(metrics.get("decode_including_setup"),
                        "leopard1.decode_including_setup", iterations,
                        "median_us_per_batch_call")
        require(metrics.get("codec_setup") is None and
                metrics.get("decode_timing_includes_setup") is True and
                metrics.get("rate_semantics") ==
                    "offered_received counts every non-null original and all R "
                    "supplied parity shards",
                "Leopard-main timing semantics changed")
    elif codec == "wirehair":
        require(payload.get("schema") ==
                "leopard2-performance-atlas-wirehair-v1/v2",
                "unexpected Wirehair adapter schema")
        require(payload.get("correctness", {}).get("round_trip") is True,
                "Wirehair round trip failed")
        require(payload.get("build", {}).get("wirehair_source_commit") ==
                EXPECTED_WIREHAIR_COMMIT, "Wirehair source commit differs")
        build = payload.get("build")
        require(isinstance(build, dict) and
                build.get("pure_avx2") is False and
                build.get("isa_claim") == "wirehair_v1_compact_path_avx2" and
                build.get("target_qualified_avx512_helpers_present") is True and
                build.get("wirehair_v1_wide_xor_forced_off") is True and
                build.get("runtime_wide_xor_enabled") is False and
                build.get("measured_path_avx512") is False,
                "Wirehair measured-path ISA contract differs")
        active_features = build.get("active_x86_features")
        require(isinstance(active_features, dict) and
                active_features.get("avx2") is True and
                active_features.get("gfni") is False,
                "Wirehair active feature contract differs")
        require(build.get("wirehair_abi_version") == 2 and
                build.get("wire_profile_version") == 1 and
                build.get("wire_profile_id") ==
                    EXPECTED_WIREHAIR_PROFILE_ID,
                "Wirehair ABI/wire-profile identity differs")
        require(payload.get("path_semantics") == {
                    "codec": "shipping_wirehair_v1",
                    "threading": "single_caller_thread",
                    "wide_xor": "forced_off_on_benchmark_thread",
                    "avx512_target_helpers":
                        "present_but_unreachable_from_measured_v1_compact_path",
                }, "Wirehair execution-path semantics changed")
        require(payload.get("timing_semantics") == {
                    "message_precode_setup":
                        "fresh wirehair_encoder_create_ex; no repair emission",
                    "encode_execution":
                        "reuse one message-precode encoder; emit exactly R "
                        "repairs; exclude encoder creation",
                    "encode_including_setup":
                        "fresh wirehair_encoder_create_ex then emit exactly R "
                        "repairs",
                    "decode_including_setup":
                        "fresh decoder; ingest surviving systematic blocks then "
                        "ascending repair blocks through solve; recover missing "
                        "originals",
                }, "Wirehair timing semantics changed")
        for name in ("message_precode_setup", "encode_execution",
                     "encode_including_setup", "decode_including_setup"):
            validate_metric(metrics.get(name), f"wirehair.{name}", iterations,
                            "median_us")
        decode_input = payload.get("decode_input")
        require(isinstance(decode_input, dict), "Wirehair lacks decode input")
        consumed = decode_input.get("repair_blocks_consumed")
        extra = decode_input.get("extra_repair_blocks")
        require(isinstance(consumed, int) and consumed >= cell["loss_count"],
                "Wirehair consumed too few repair blocks")
        require(consumed <= cell["loss_count"] + 32 and
                decode_input.get("surviving_systematic_blocks") ==
                    cell["K"] - cell["loss_count"] and
                decode_input.get("arrival_order") ==
                    "surviving_systematic_ascending_then_repair_ascending",
                "Wirehair decode input geometry differs")
        require(extra == consumed - cell["loss_count"],
                "Wirehair repair overhead is inconsistent")
    else:
        raise AtlasError(f"unknown codec {codec}")


def compare_cell_payloads(payloads: Mapping[str, Mapping[str, Any]],
                          cell: Mapping[str, Any]) -> None:
    reference = payloads["leopard2"]
    ref_parameters = reference["parameters"]
    ref_digests = reference["workload_digests"]
    for codec, payload in payloads.items():
        require(payload["parameters"]["missing_original_indices"] ==
                ref_parameters["missing_original_indices"],
                f"{codec} loss pattern differs for {cell['id']}")
        require(payload["workload_digests"]["original_data"] ==
                ref_digests["original_data"],
                f"{codec} original data differs for {cell['id']}")
        require(payload["workload_digests"]["recovered_originals"] ==
                ref_digests["recovered_originals"],
                f"{codec} recovered data differs for {cell['id']}")
    if "leopard1" in payloads:
        require(payloads["leopard1"]["workload_digests"]["transmitted_parity"] ==
                ref_digests["transmitted_parity"],
                f"Leopard2 parity differs from exact main for {cell['id']}")


def capture_host(cpu: int, lock_identity: Mapping[str, Any],
                 as_limit: int) -> dict[str, Any]:
    def optional(command: Sequence[str]) -> str | None:
        try:
            completed = subprocess.run(command, check=False,
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                text=True, timeout=10)
            return completed.stdout.strip() if completed.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            return None
    return {
        "captured_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hostname": platform.node(), "platform": platform.platform(),
        "kernel": platform.release(), "python": platform.python_version(),
        "allowed_cpus": allowed_cpus(), "benchmark_cpu": cpu,
        "isa_comparison": (
            "Leopard2/main pure AVX2; Wirehair-v1 compact AVX2 path with "
            "target-qualified AVX-512 helpers present but wide-XOR forced off"),
        "address_space_limit_bytes": as_limit,
        "canonical_lock": dict(lock_identity),
        "lscpu": optional(["lscpu", "--json"]),
        "numactl_hardware": optional(["numactl", "--hardware"]),
    }


def resume_host_identity(host: Mapping[str, Any]) -> dict[str, Any]:
    """Return only stable host facts; lscpu includes live frequency fields."""
    lock = host.get("canonical_lock")
    lock_path = lock.get("path") if isinstance(lock, dict) else None
    return {
        key: host.get(key) for key in (
            "hostname", "platform", "kernel", "python", "allowed_cpus",
            "benchmark_cpu", "address_space_limit_bytes", "isa_comparison")
    } | {"canonical_lock_path": lock_path}


def reproduction_command(args: argparse.Namespace) -> str:
    command = [
        sys.executable, str(Path(__file__).resolve()), "all",
        "--output", str(args.output.resolve()),
        "--leopard2", str(args.leopard2.resolve()),
        "--leopard2-build-root", str(args.leopard2_build_root.resolve()),
        "--leopard2-sha256", args.leopard2_sha256,
        "--leopard1", str(args.leopard1.resolve()),
        "--leopard1-build-root", str(args.leopard1_build_root.resolve()),
        "--leopard1-sha256", args.leopard1_sha256,
        "--wirehair", str(args.wirehair.resolve()),
        "--wirehair-build-root", str(args.wirehair_build_root.resolve()),
        "--wirehair-sha256", args.wirehair_sha256,
        "--leopard-source", str(args.leopard_source.resolve()),
        "--main-source", str(args.main_source.resolve()),
        "--wirehair-source", str(args.wirehair_source.resolve()),
        "--address-space-limit", str(args.address_space_limit),
        "--timeout", str(args.timeout),
    ]
    if args.cpu is not None:
        command += ["--cpu", str(args.cpu)]
    if args.max_cells is not None:
        command += ["--max-cells", str(args.max_cells)]
    return shlex.join(command)


def load_cell_results(raw_root: Path, cell: Mapping[str, Any]) -> dict[str, Any]:
    payloads: dict[str, Any] = {}
    for codec in cell["available_codecs"]:
        path = raw_root / cell["id"] / f"{codec}.json"
        require(path.exists(), f"missing raw result: {path}")
        payload = read_json(path)
        validate_payload(codec, payload, cell)
        payloads[codec] = payload
    compare_cell_payloads(payloads, cell)
    return payloads


def run_campaign(args: argparse.Namespace, manifest: Mapping[str, Any]) -> None:
    output_root = args.output.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "manifest.json"
    if manifest_path.exists():
        require(read_json(manifest_path) == manifest,
                "existing manifest differs from this runner")
    else:
        atomic_write_json(manifest_path, manifest)

    with canonical_lock() as lock_identity:
        cpus = allowed_cpus()
        cpu = args.cpu if args.cpu is not None else physical_cpu(cpus)
        require(cpu in cpus, f"CPU {cpu} is not in the allowed affinity set")
        leopard2_closure = capture_leopard2_build_closure(
            args.leopard2_build_root, args.leopard_source, args.leopard2)
        leopard1_closure = capture_external_avx2_compile_build(
            args.leopard1_build_root, args.leopard1,
            Path("libleopard_main_exact.a"),
            Path("CMakeFiles/leopard_main_exact.dir/link.txt"),
            {name: args.main_source / name
             for name in MAIN_ARCHIVE_SOURCE_RELATIVE},
            {"legacy_main_benchmark.cpp": args.leopard_source /
             "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"},
            {"LEO_MAIN_PURE_AVX2": "ON", "CMAKE_BUILD_TYPE": "Release"},
            "CMakeFiles/leopard_main_exact.dir/",
            "CMakeFiles/leopard_main_benchmark.dir/",
            "Leopard main")
        wirehair_closure = capture_external_avx2_compile_build(
            args.wirehair_build_root, args.wirehair,
            Path("wirehair-source/libwirehair.a"),
            Path("wirehair-source/CMakeFiles/wirehair.dir/link.txt"),
            {name: args.wirehair_source / relative
             for name, relative in WIREHAIR_ARCHIVE_SOURCE_RELATIVE.items()},
            {"wirehair_v1_benchmark.cpp": args.leopard_source /
             "experiments/leopard2/performance_atlas/"
             "wirehair_v1_benchmark.cpp"},
            {"MARCH_NATIVE": "OFF", "BUILD_TESTS": "OFF",
             "BUILD_CODEC_V2": "OFF", "CMAKE_BUILD_TYPE": "Release"},
            "wirehair-source/CMakeFiles/wirehair.dir/",
            "CMakeFiles/wirehair_v1_benchmark.dir/",
            "Wirehair")
        frozen_root = output_root / "frozen"
        paths = {
            "leopard2": args.leopard2,
            "leopard1": args.leopard1,
            "wirehair": args.wirehair,
        }
        hashes = {
            "leopard2": args.leopard2_sha256,
            "leopard1": args.leopard1_sha256,
            "wirehair": args.wirehair_sha256,
        }
        frozen: dict[str, dict[str, Any]] = {}
        for codec in CODECS:
            destination = frozen_root / f"{codec}-{hashes[codec][:16]}"
            frozen[codec] = freeze_executable(paths[codec], destination,
                                                hashes[codec])
        sources = {
            "leopard2": source_identity(args.leopard_source),
            "leopard1": source_identity(args.main_source,
                                           EXPECTED_MAIN_COMMIT),
            "wirehair": source_identity(args.wirehair_source,
                                           EXPECTED_WIREHAIR_COMMIT),
        }
        attestation = attest_leopard2(
            Path(frozen["leopard2"]["path"]), cpu,
            args.address_space_limit, args.timeout, sources["leopard2"])
        metadata = {
            "schema": "leopard2-performance-atlas-run-metadata/v1",
            "host": capture_host(cpu, lock_identity, args.address_space_limit),
            "executables": frozen,
            "sources": sources,
            "leopard2_embedded_source_attestation": attestation["build"],
            "build_closures": {
                "leopard2": leopard2_closure,
                "leopard1": leopard1_closure,
                "wirehair": wirehair_closure,
            },
            "runner": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256_file(Path(__file__).resolve()),
                "version": RUNNER_VERSION,
            },
            "reproduction_command": reproduction_command(args),
        }
        metadata_path = output_root / "run_metadata.json"
        if metadata_path.exists():
            old = read_json(metadata_path)
            # The capture time is deliberately not part of resume identity.
            require(resume_host_identity(old.get("host", {})) ==
                    resume_host_identity(metadata["host"]) and
                    old.get("executables") == metadata["executables"] and
                    old.get("sources") == metadata["sources"] and
                    old.get("leopard2_embedded_source_attestation") ==
                    metadata["leopard2_embedded_source_attestation"] and
                    old.get("build_closures") == metadata["build_closures"] and
                    old.get("runner") == metadata["runner"] and
                    old.get("reproduction_command") ==
                    metadata["reproduction_command"],
                    "existing run metadata differs from current campaign")
        else:
            atomic_write_json(metadata_path, metadata)
        reproduce_path = output_root / "REPRODUCE.txt"
        reproduce = metadata["reproduction_command"] + "\n"
        if reproduce_path.exists():
            require(reproduce_path.read_text(encoding="utf-8") == reproduce,
                    "existing reproduction command differs")
        else:
            atomic_write_bytes(reproduce_path, reproduce.encode("utf-8"))

        cells = list(manifest["cells"])
        if args.max_cells is not None:
            cells = cells[:args.max_cells]
        raw_root = output_root / "raw"
        total = sum(len(cell["available_codecs"]) for cell in cells)
        complete = 0
        started = time.monotonic()
        for cell in cells:
            cell_root = raw_root / cell["id"]
            cell_root.mkdir(parents=True, exist_ok=True)
            for codec in codec_order(cell["id"]):
                if codec not in cell["available_codecs"]:
                    continue
                result_path = cell_root / f"{codec}.json"
                stderr_path = cell_root / f"{codec}.stderr.txt"
                if result_path.exists():
                    validate_payload(codec, read_json(result_path), cell)
                    complete += 1
                    continue
                for identity in frozen.values():
                    verify_identity(identity)
                print(f"[{complete + 1}/{total}] {cell['id']} {codec}",
                      flush=True)
                payload, stderr = run_invocation(
                    codec, Path(frozen[codec]["path"]), cell, cpu,
                    args.address_space_limit, args.timeout)
                atomic_write_json(result_path, payload)
                atomic_write_bytes(stderr_path, stderr.encode("utf-8"))
                complete += 1
            load_cell_results(raw_root, cell)
            elapsed = time.monotonic() - started
            rate = complete / elapsed if elapsed > 0 else 0.0
            eta = (total - complete) / rate if rate > 0 else math.inf
            print(f"validated {cell['id']}; {complete}/{total}, "
                  f"elapsed={elapsed:.1f}s eta={eta:.1f}s", flush=True)
        for identity in frozen.values():
            verify_identity(identity)


def metric_median(metric: Mapping[str, Any], setup: bool = False) -> float:
    key = "median_us" if setup else "median_us_per_batch_call"
    return finite_positive(metric.get(key), key)


def normalize_result(codec: str, payload: Mapping[str, Any],
                     cell: Mapping[str, Any]) -> dict[str, Any]:
    metrics = payload["metrics"]
    if codec == "leopard2":
        codec_setup = metric_median(metrics["codec_setup"], True)
        encode_execution = metric_median(metrics["encode_execution"])
        plan_setup = metric_median(metrics["decode_plan_setup"], True)
        decode_execution = metric_median(metrics["decode_execution"])
        decode_first = metric_median(
            metrics["one_shot_decode_including_setup"])
        result = {
            "encode_execution_us": encode_execution,
            "decode_first_us": decode_first,
            "decode_reused_us": decode_execution,
            "codec_setup_us": codec_setup,
            "decode_plan_setup_us": plan_setup,
            "extra_repair_blocks": 0,
            "profile": payload["resolved"]["profile"],
            "field": payload["resolved"]["field"],
            "backend": payload["resolved"]["backend"],
            "encode_scratch_bytes": payload["memory"][
                "encode_scratch_bytes_per_stripe"],
            "decode_scratch_bytes": payload["memory"][
                "decode_scratch_bytes_per_stripe"],
        }
    elif codec == "leopard1":
        encode = metric_median(metrics["encode_execution"])
        decode = metric_median(metrics["decode_including_setup"])
        result = {
            "encode_execution_us": encode,
            "decode_first_us": decode, "decode_reused_us": decode,
            "codec_setup_us": None, "decode_plan_setup_us": None,
            "extra_repair_blocks": 0,
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "legacy runtime dispatch",
            "encode_scratch_bytes": payload["memory"][
                "encode_work_bytes_batch"],
            "decode_scratch_bytes": payload["memory"][
                "decode_work_bytes_batch"],
        }
    elif codec == "wirehair":
        complete_encode = metric_median(
            metrics["encode_including_setup"], True)
        result = {
            # Wirehair create consumes the message and performs its precode;
            # subsequent wirehair_encode calls only emit repair rows.  The
            # complete create+emit measurement is the comparable encode.
            "encode_execution_us": complete_encode,
            "wirehair_repair_emission_us": metric_median(
                metrics["encode_execution"], True),
            "decode_first_us": metric_median(
                metrics["decode_including_setup"], True),
            "decode_reused_us": None,
            "codec_setup_us": metric_median(
                metrics["message_precode_setup"], True),
            "decode_plan_setup_us": None,
            "extra_repair_blocks": payload["decode_input"][
                "extra_repair_blocks"],
            "profile": "shipping fountain codec", "field": "gf256 internals",
            "backend": "wirehair runtime",
            "encode_scratch_bytes": None, "decode_scratch_bytes": None,
        }
    else:
        raise AtlasError(f"unknown codec {codec}")
    k = cell["K"]
    block_bytes = cell["shard_bytes"]
    losses = cell["loss_count"]
    source_bytes = k * block_bytes
    parity_bytes = R_COUNT * block_bytes
    repair_bytes = losses * block_bytes
    received_shards = (k - losses +
        (int(payload["decode_input"]["repair_blocks_consumed"])
         if codec == "wirehair" else R_COUNT))
    received_bytes = received_shards * block_bytes
    encode_us = result["encode_execution_us"]
    result["encode_execution_message_GBps"] = (
        source_bytes / (encode_us * 1000.0))
    result["encode_execution_output_GBps"] = (
        parity_bytes / (encode_us * 1000.0))
    if codec == "wirehair":
        emission_us = result["wirehair_repair_emission_us"]
        result["wirehair_repair_emission_output_GBps"] = (
            parity_bytes / (emission_us * 1000.0))
    for prefix, usec in (
            ("decode_first", result["decode_first_us"]),
            ("decode_reused", result["decode_reused_us"])):
        if usec is None:
            result[prefix + "_message_GBps"] = None
            result[prefix + "_repaired_GBps"] = None
            result[prefix + "_received_GBps"] = None
        else:
            result[prefix + "_message_GBps"] = source_bytes / (usec * 1000.0)
            result[prefix + "_repaired_GBps"] = repair_bytes / (usec * 1000.0)
            result[prefix + "_received_GBps"] = received_bytes / (usec * 1000.0)
    result["decode_received_shards"] = received_shards
    return result


def build_summary(output_root: Path, manifest: Mapping[str, Any],
                  allow_partial: bool = False) -> dict[str, Any]:
    raw_root = output_root / "raw"
    rows: list[dict[str, Any]] = []
    complete_cells = 0
    for cell in manifest["cells"]:
        available = [raw_root / cell["id"] / f"{codec}.json"
                     for codec in cell["available_codecs"]]
        if not all(path.exists() for path in available):
            if allow_partial:
                continue
            missing = [str(path) for path in available if not path.exists()]
            raise AtlasError("incomplete campaign: " + ", ".join(missing))
        payloads = load_cell_results(raw_root, cell)
        for codec, payload in payloads.items():
            row = {
                "cell_id": cell["id"], "codec": codec,
                "K": cell["K"], "R": cell["R"],
                "shard_bytes": cell["shard_bytes"],
                "loss_count": cell["loss_count"],
                "loss_labels": list(cell["loss_labels"]),
                "seed": cell["seed"],
                "missing_original_indices": payload["parameters"][
                    "missing_original_indices"],
                "original_digest": payload["workload_digests"]["original_data"],
                "recovered_digest": payload["workload_digests"][
                    "recovered_originals"],
            }
            row.update(normalize_result(codec, payload, cell))
            rows.append(row)
        complete_cells += 1
    summary = {
        "schema": SUMMARY_SCHEMA,
        "manifest_sha256": sha256_bytes(canonical_json(manifest).encode()),
        "run_metadata_sha256": sha256_file(output_root / "run_metadata.json"),
        "complete": complete_cells == len(manifest["cells"]),
        "complete_cell_count": complete_cells,
        "expected_cell_count": len(manifest["cells"]),
        "rows": rows,
    }
    validate_summary(summary, manifest, require_complete=not allow_partial)
    return summary


def validate_summary(summary: Mapping[str, Any], manifest: Mapping[str, Any],
                     require_complete: bool = True) -> None:
    require(summary.get("schema") == SUMMARY_SCHEMA,
            "unexpected summary schema")
    expected_manifest_hash = sha256_bytes(canonical_json(manifest).encode())
    require(summary.get("manifest_sha256") == expected_manifest_hash,
            "summary is bound to a different manifest")
    metadata_hash = summary.get("run_metadata_sha256")
    require(isinstance(metadata_hash, str) and len(metadata_hash) == 64 and
            all(character in "0123456789abcdef" for character in metadata_hash),
            "summary lacks a run-metadata identity")
    if require_complete:
        require(summary.get("complete") is True, "summary is incomplete")
    require(summary.get("expected_cell_count") == len(manifest["cells"]),
            "summary expected-cell count differs")
    rows = summary.get("rows")
    require(isinstance(rows, list), "summary rows are missing")
    cells = {cell["id"]: cell for cell in manifest["cells"]}
    seen: set[tuple[str, str]] = set()
    represented: dict[str, set[str]] = {}
    for row in rows:
        require(isinstance(row, dict), "summary row is not an object")
        key = row.get("cell_id"), row.get("codec")
        require(key not in seen, f"duplicate summary row {key}")
        seen.add(key)
        cell_id, codec = key
        require(isinstance(cell_id, str) and cell_id in cells,
                "summary references an unknown cell")
        cell = cells[cell_id]
        require(codec in cell["available_codecs"],
                "summary codec is unavailable for its cell")
        represented.setdefault(cell_id, set()).add(codec)
        exact = {
            "K": cell["K"], "R": cell["R"],
            "shard_bytes": cell["shard_bytes"],
            "loss_count": cell["loss_count"],
            "loss_labels": cell["loss_labels"], "seed": cell["seed"],
        }
        for name, value in exact.items():
            require(row.get(name) == value,
                    f"summary {cell_id}/{codec} {name} differs")
        for name in ("encode_execution_us", "decode_first_us"):
            finite_positive(row.get(name), f"summary.{name}")
        reused = row.get("decode_reused_us")
        if codec == "wirehair":
            require(reused is None, "Wirehair cannot claim a reusable plan")
        else:
            finite_positive(reused, "summary.decode_reused_us")
        k_bytes = cell["K"] * cell["shard_bytes"]
        parity_bytes = cell["R"] * cell["shard_bytes"]
        repaired_bytes = cell["loss_count"] * cell["shard_bytes"]

        def rate(name: str, byte_count: int, usec_name: str) -> None:
            actual = finite_positive(row.get(name), f"summary.{name}")
            expected = byte_count / (float(row[usec_name]) * 1000.0)
            require(math.isclose(actual, expected, rel_tol=1e-12,
                                 abs_tol=1e-15),
                    f"summary {cell_id}/{codec} has forged {name}")

        rate("encode_execution_message_GBps", k_bytes,
             "encode_execution_us")
        rate("encode_execution_output_GBps", parity_bytes,
             "encode_execution_us")
        rate("decode_first_message_GBps", k_bytes, "decode_first_us")
        rate("decode_first_repaired_GBps", repaired_bytes,
             "decode_first_us")
        received_shards = row.get("decode_received_shards")
        require(isinstance(received_shards, int) and received_shards > 0,
                "summary decode-received shard count is invalid")
        received_bytes = received_shards * cell["shard_bytes"]
        rate("decode_first_received_GBps", received_bytes,
             "decode_first_us")
        if codec != "wirehair":
            rate("decode_reused_message_GBps", k_bytes,
                 "decode_reused_us")
            rate("decode_reused_repaired_GBps", repaired_bytes,
                 "decode_reused_us")
            rate("decode_reused_received_GBps", received_bytes,
                 "decode_reused_us")
        else:
            require(row.get("decode_reused_message_GBps") is None and
                    row.get("decode_reused_repaired_GBps") is None and
                    row.get("decode_reused_received_GBps") is None,
                    "Wirehair reusable rates must remain unavailable")
            emission = finite_positive(row.get("wirehair_repair_emission_us"),
                                       "summary.wirehair_repair_emission_us")
            expected_emission_rate = parity_bytes / (emission * 1000.0)
            emission_rate = finite_positive(row.get(
                "wirehair_repair_emission_output_GBps"),
                "summary.wirehair_repair_emission_output_GBps")
            require(math.isclose(emission_rate,
                        expected_emission_rate, rel_tol=1e-12,
                        abs_tol=1e-15),
                    "Wirehair marginal emission rate is forged")
    for cell_id, codecs in represented.items():
        require(codecs == set(cells[cell_id]["available_codecs"]),
                f"summary cell {cell_id} is only partially represented")
    complete_count = len(represented)
    require(summary.get("complete_cell_count") == complete_count,
            "summary complete-cell count differs")
    is_complete = complete_count == len(cells)
    require(summary.get("complete") is is_complete,
            "summary completeness flag differs")
    if require_complete:
        expected_rows = {(cell["id"], codec) for cell in manifest["cells"]
                         for codec in cell["available_codecs"]}
        require(seen == expected_rows, "complete summary row set differs")


def write_raw_bundle(output_root: Path, manifest: Mapping[str, Any],
                     destination: Path | None = None) -> Path:
    records: list[dict[str, Any]] = []
    raw_root = output_root / "raw"
    for cell in manifest["cells"]:
        for codec, payload in load_cell_results(raw_root, cell).items():
            records.append({"cell_id": cell["id"], "codec": codec,
                            "payload": payload})
    metadata_hash = sha256_file(output_root / "run_metadata.json")
    value = {"schema": RAW_BUNDLE_SCHEMA,
             "manifest_sha256": sha256_bytes(canonical_json(manifest).encode()),
             "run_metadata_sha256": metadata_hash,
             "record_count": len(records), "records": records}
    validate_raw_bundle(value, manifest, metadata_hash)
    path = destination or output_root / "raw_bundle.json"
    atomic_write_json(path, value)
    return path


def validate_raw_bundle(bundle: Mapping[str, Any], manifest: Mapping[str, Any],
                        run_metadata_sha256: str) -> None:
    require(bundle.get("schema") == RAW_BUNDLE_SCHEMA,
            "raw bundle schema differs")
    require(bundle.get("manifest_sha256") ==
            sha256_bytes(canonical_json(manifest).encode()),
            "raw bundle manifest identity differs")
    require(bundle.get("run_metadata_sha256") == run_metadata_sha256,
            "raw bundle run-metadata identity differs")
    records = bundle.get("records")
    require(isinstance(records, list) and
            bundle.get("record_count") == len(records),
            "raw bundle record count differs")
    cells = {cell["id"]: cell for cell in manifest["cells"]}
    expected = {(cell["id"], codec) for cell in manifest["cells"]
                for codec in cell["available_codecs"]}
    seen: set[tuple[str, str]] = set()
    per_cell: dict[str, dict[str, Any]] = {}
    for record in records:
        require(isinstance(record, dict), "raw bundle record is malformed")
        cell_id, codec = record.get("cell_id"), record.get("codec")
        require(isinstance(cell_id, str) and cell_id in cells and
                codec in cells[cell_id]["available_codecs"],
                "raw bundle record identifies an unavailable cell/codec")
        key = cell_id, codec
        require(key not in seen, "raw bundle record is duplicated")
        payload = record.get("payload")
        validate_payload(codec, payload, cells[cell_id])
        seen.add(key)
        per_cell.setdefault(cell_id, {})[codec] = payload
    require(seen == expected, "raw bundle record set differs")
    for cell_id, payloads in per_cell.items():
        compare_cell_payloads(payloads, cells[cell_id])


def format_bytes(value: int) -> str:
    if value == 1024 * 1024:
        return "1 MiB"
    if value % 1024 == 0:
        return f"{value // 1024} KiB"
    return f"{value} B"


def svg_text(x: float, y: float, text_value: str, *, size: int = 12,
             anchor: str = "start", weight: str = "normal",
             color: str = "#24292f", rotate: int | None = None) -> str:
    transform = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
            f'text-anchor="{anchor}" font-weight="{weight}" fill="{color}"'
            f'{transform}>{html.escape(text_value)}</text>')


def nice_log_ticks(low: float, high: float) -> list[float]:
    if low <= 0 or high <= 0:
        return []
    start = math.floor(math.log10(low))
    stop = math.ceil(math.log10(high))
    ticks: list[float] = []
    for exponent in range(start, stop + 1):
        for multiplier in (1, 2, 5):
            value = multiplier * (10 ** exponent)
            if low * 0.999 <= value <= high * 1.001:
                ticks.append(value)
    # Narrow panels can sit wholly between adjacent 1/2/5 ticks.  Boundaries
    # guarantee a labeled axis even then; near-duplicates are removed.
    ticks.extend((low, high))
    ordered = sorted(ticks)
    result: list[float] = []
    for value in ordered:
        if not result or not math.isclose(value, result[-1], rel_tol=1e-6):
            result.append(value)
    return result


def format_axis(value: float) -> str:
    if value >= 1000:
        return f"{value:.0f}"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2g}"
    if value >= 0.01:
        return f"{value:.2f}"
    return f"{value:.1e}"


@dataclasses.dataclass(frozen=True)
class PlotSpec:
    filename: str
    title: str
    subtitle: str
    metric: str
    y_label: str
    loss_label: str | None = None
    codecs: tuple[str, ...] = CODECS
    baseline_codec: str | None = None
    logarithmic_y: bool = True


def index_rows(summary: Mapping[str, Any]) -> dict[tuple[int, int, str], dict[str, Any]]:
    return {(row["K"], row["shard_bytes"], row["codec"]): row
            for row in summary["rows"]}


def rows_for_label(summary: Mapping[str, Any], loss_label: str | None,
                   block_bytes: int) -> list[dict[str, Any]]:
    rows = []
    for row in summary["rows"]:
        if row["shard_bytes"] != block_bytes:
            continue
        if loss_label is None:
            # Encoding is loss-independent.  Use the canonical one-loss cell.
            if "one" not in row["loss_labels"]:
                continue
        elif loss_label not in row["loss_labels"]:
            continue
        rows.append(row)
    return rows


def render_plot(summary: Mapping[str, Any], spec: PlotSpec, path: Path) -> None:
    width, height = 1400, 900
    margin_x, margin_top, gap_x, gap_y = 78, 105, 58, 60
    panel_width = (width - 2 * margin_x - gap_x) / 2
    panel_height = (height - margin_top - 65 - gap_y) / 2
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text { font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }</style>',
        svg_text(width / 2, 35, spec.title, size=23, anchor="middle", weight="600"),
        svg_text(width / 2, 61, spec.subtitle, size=13, anchor="middle", color="#57606a"),
    ]
    legend_x = width / 2 - (len(spec.codecs) * 145) / 2
    for index, codec in enumerate(spec.codecs):
        x = legend_x + index * 145
        elements.append(f'<line x1="{x:.1f}" y1="82" x2="{x + 28:.1f}" y2="82" '
                        f'stroke="{CODEC_COLORS[codec]}" stroke-width="3"/>')
        elements.append(svg_text(x + 35, 86, CODEC_LABELS[codec], size=12))

    for panel_index, block_bytes in enumerate(BLOCK_BYTES):
        column = panel_index % 2
        row_index = panel_index // 2
        left = margin_x + column * (panel_width + gap_x)
        top = margin_top + row_index * (panel_height + gap_y)
        plot_left, plot_top = left + 65, top + 30
        plot_width, plot_height = panel_width - 82, panel_height - 70
        panel_rows = rows_for_label(summary, spec.loss_label, block_bytes)
        series: dict[str, list[tuple[int, float]]] = {}
        for codec in spec.codecs:
            values: list[tuple[int, float]] = []
            for item in panel_rows:
                if item["codec"] != codec:
                    continue
                value = item.get(spec.metric)
                if value is None:
                    continue
                if spec.baseline_codec is not None:
                    baseline = next((candidate for candidate in panel_rows
                        if candidate["codec"] == spec.baseline_codec and
                        candidate["K"] == item["K"] and
                        candidate["loss_count"] == item["loss_count"]), None)
                    if baseline is None:
                        continue
                    baseline_value = baseline.get(spec.metric)
                    if baseline_value is None:
                        continue
                    # Time metrics are lower-is-better; rates are higher-is-better.
                    if spec.metric.endswith("_us"):
                        value = float(baseline_value) / float(value)
                    else:
                        value = float(value) / float(baseline_value)
                value = float(value)
                if math.isfinite(value) and value > 0:
                    values.append((item["K"], value))
            # A conceptual loss label can share a raw cell with another label,
            # but still contributes one point per K and codec.
            series[codec] = sorted(set(values))
        all_values = [value for values in series.values() for _, value in values]
        require(all_values, f"plot {spec.filename} panel {block_bytes} has no data")
        low, high = min(all_values), max(all_values)
        if spec.baseline_codec is not None:
            low, high = min(low, 1.0), max(high, 1.0)
        if spec.logarithmic_y:
            low = 10 ** (math.floor(math.log10(low) * 5) / 5)
            high = 10 ** (math.ceil(math.log10(high) * 5) / 5)
            if high <= low:
                high = low * 1.25
            y_map = lambda value: plot_top + plot_height * (
                math.log(high) - math.log(value)) / (math.log(high) - math.log(low))
            ticks = nice_log_ticks(low, high)
        else:
            span = high - low
            low = max(0.0, low - span * 0.08)
            high += max(span * 0.08, high * 0.01)
            y_map = lambda value: plot_top + plot_height * (high - value) / (high - low)
            ticks = [low + (high - low) * index / 5 for index in range(6)]
        x_map = lambda value: plot_left + plot_width * math.log2(value) / math.log2(K_MAX)
        elements.append(f'<rect x="{left:.1f}" y="{top:.1f}" width="{panel_width:.1f}" '
                        f'height="{panel_height:.1f}" rx="6" fill="#f6f8fa" stroke="#d0d7de"/>')
        elements.append(svg_text(left + panel_width / 2, top + 20,
                                 format_bytes(block_bytes), size=14,
                                 anchor="middle", weight="600"))
        for tick in ticks:
            y = y_map(tick)
            if plot_top - 1 <= y <= plot_top + plot_height + 1:
                elements.append(f'<line x1="{plot_left:.1f}" y1="{y:.1f}" '
                    f'x2="{plot_left + plot_width:.1f}" y2="{y:.1f}" '
                    'stroke="#d8dee4" stroke-width="1"/>')
                elements.append(svg_text(plot_left - 7, y + 4, format_axis(tick),
                                         size=10, anchor="end", color="#57606a"))
        x_ticks = (1, 2, 4, 8, 16, 32, 64, 128, K_MAX)
        for tick in x_ticks:
            x = x_map(tick)
            elements.append(f'<line x1="{x:.1f}" y1="{plot_top:.1f}" '
                f'x2="{x:.1f}" y2="{plot_top + plot_height:.1f}" '
                'stroke="#eaeef2" stroke-width="1"/>')
            elements.append(svg_text(x, plot_top + plot_height + 17, str(tick),
                                     size=10, anchor="middle", color="#57606a"))
        if low <= 1.0 <= high and spec.baseline_codec is not None:
            y = y_map(1.0)
            elements.append(f'<line x1="{plot_left:.1f}" y1="{y:.1f}" '
                f'x2="{plot_left + plot_width:.1f}" y2="{y:.1f}" '
                'stroke="#6e7781" stroke-width="1.5" stroke-dasharray="5 4"/>')
        for codec, values in series.items():
            if not values:
                continue
            points = " ".join(f"{x_map(k):.1f},{y_map(value):.1f}"
                              for k, value in values)
            elements.append(f'<polyline points="{points}" fill="none" '
                f'stroke="{CODEC_COLORS[codec]}" stroke-width="2.2" '
                'stroke-linejoin="round" stroke-linecap="round"/>')
            for k, value in values:
                if k & (k - 1) == 0 or k == K_MAX:
                    elements.append(f'<circle cx="{x_map(k):.1f}" cy="{y_map(value):.1f}" '
                        f'r="2.8" fill="#fff" stroke="{CODEC_COLORS[codec]}" '
                        'stroke-width="1.6"/>')
        elements.append(svg_text(plot_left + plot_width / 2,
                                 plot_top + plot_height + 36, "K (log₂ scale)",
                                 size=11, anchor="middle"))
        elements.append(svg_text(left + 14, plot_top + plot_height / 2,
                                 spec.y_label, size=11, anchor="middle", rotate=-90))
    elements.append(svg_text(width / 2, height - 16,
        "Lines include every odd K; hollow markers identify powers of two and K=224. Missing codec cells are not connected.",
        size=11, anchor="middle", color="#57606a"))
    elements.append("</svg>\n")
    atomic_write_bytes(path, "\n".join(elements).encode("utf-8"))


def plot_specs() -> list[PlotSpec]:
    common = "R=32; every odd K plus powers of two through K=224; pinned single core"
    specs = [
        PlotSpec("encode_execution_message_gbps.svg", "Full-message encode throughput",
                 common, "encode_execution_message_GBps", "message GB/s"),
        PlotSpec("encode_execution_output_gbps.svg", "Full-message parity/repair output throughput",
                 common + "; exactly 32 output shards generated",
                 "encode_execution_output_GBps", "output GB/s"),
        PlotSpec("encode_execution_latency_us.svg", "Full-message encode latency",
                 common, "encode_execution_us", "microseconds"),
        PlotSpec("encode_speedup_vs_leopard1.svg", "Leopard2 encode speedup over Leopard main",
                 common + "; values above 1 are faster",
                 "encode_execution_us", "speedup (×)", codecs=("leopard2",),
                 baseline_codec="leopard1"),
        PlotSpec("encode_speedup_vs_wirehair.svg", "Leopard2 encode speedup over Wirehair",
                 common + "; values above 1 are faster",
                 "encode_execution_us", "speedup (×)", codecs=("leopard2",),
                 baseline_codec="wirehair"),
        PlotSpec("wirehair_marginal_repair_output_gbps.svg",
                 "Wirehair marginal repair-row emission",
                 common + "; excludes message-dependent encoder creation/precode",
                 "wirehair_repair_emission_output_GBps", "repair output GB/s",
                 codecs=("wirehair",)),
        PlotSpec("wirehair_marginal_repair_latency_us.svg",
                 "Wirehair marginal repair-row emission latency",
                 common + "; excludes message-dependent encoder creation/precode",
                 "wirehair_repair_emission_us", "microseconds",
                 codecs=("wirehair",)),
    ]
    for loss_label in LOSS_LABELS:
        display = LOSS_DISPLAY[loss_label]
        suffix = loss_label
        specs.extend([
            PlotSpec(f"decode_{suffix}_first_message_gbps.svg",
                     f"Decode throughput: {display}",
                     common + "; setup included", "decode_first_message_GBps",
                     "message-equivalent GB/s", loss_label),
            PlotSpec(f"decode_{suffix}_first_repaired_gbps.svg",
                     f"Repaired-output throughput: {display}",
                     common + "; setup included", "decode_first_repaired_GBps",
                     "repaired GB/s", loss_label),
            PlotSpec(f"decode_{suffix}_first_received_gbps.svg",
                     f"Received-input throughput: {display}",
                     common + "; Leopard offers all R parity; Wirehair counts consumed repairs",
                     "decode_first_received_GBps", "received-input GB/s",
                     loss_label),
            PlotSpec(f"decode_{suffix}_first_latency_us.svg",
                     f"Decode latency: {display}", common + "; setup included",
                     "decode_first_us", "microseconds", loss_label),
            PlotSpec(f"decode_{suffix}_reused_message_gbps.svg",
                     f"Repeated decode-call throughput: {display}",
                     common + "; Leopard2 reuses a plan; Leopard main repeats setup; Wirehair omitted",
                     "decode_reused_message_GBps", "message-equivalent GB/s",
                     loss_label, codecs=("leopard2", "leopard1")),
            PlotSpec(f"decode_{suffix}_speedup_vs_leopard1.svg",
                     f"Leopard2 decode speedup over Leopard main: {display}",
                     common + "; setup included; values above 1 are faster",
                     "decode_first_us", "speedup (×)", loss_label,
                     codecs=("leopard2",), baseline_codec="leopard1"),
            PlotSpec(f"decode_{suffix}_speedup_vs_wirehair.svg",
                     f"Leopard2 decode speedup over Wirehair: {display}",
                     common + "; setup included; values above 1 are faster",
                     "decode_first_us", "speedup (×)", loss_label,
                     codecs=("leopard2",), baseline_codec="wirehair"),
        ])
    return specs


def render_overhead_plot(summary: Mapping[str, Any], loss_label: str,
                         path: Path) -> None:
    # Reuse the generic renderer with a one-off derived metric copied into rows.
    copied = json.loads(json.dumps(summary))
    for row in copied["rows"]:
        value = row.get("extra_repair_blocks")
        row["extra_repair_blocks_plot"] = (float(value) + 0.05
            if row["codec"] == "wirehair" else None)
    render_plot(copied, PlotSpec(
        path.name, f"Wirehair repair overhead: {LOSS_DISPLAY[loss_label]}",
        "Displayed value is extra blocks + 0.05 so zero remains visible on a log axis",
        "extra_repair_blocks_plot", "extra repair blocks + 0.05", loss_label,
        codecs=("wirehair",)), path)


def generate_plots(summary: Mapping[str, Any], output: Path) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for spec in plot_specs():
        path = output / spec.filename
        render_plot(summary, spec, path)
        paths.append(path)
    for loss_label in LOSS_LABELS:
        path = output / f"wirehair_overhead_{loss_label}.svg"
        render_overhead_plot(summary, loss_label, path)
        paths.append(path)
    return paths


def geometric_median(values: Sequence[float]) -> float | None:
    positive = [value for value in values if value > 0 and math.isfinite(value)]
    if not positive:
        return None
    return math.exp(statistics.median([math.log(value) for value in positive]))


def speedup_statistics(summary: Mapping[str, Any], baseline: str,
                       metric: str, loss_label: str | None = None) -> dict[str, Any]:
    rows = summary["rows"]
    lookup = {(row["K"], row["shard_bytes"], row["loss_count"], row["codec"]): row
              for row in rows}
    ratios: list[tuple[float, dict[str, Any]]] = []
    for row in rows:
        if row["codec"] != "leopard2":
            continue
        if loss_label is not None and loss_label not in row["loss_labels"]:
            continue
        if loss_label is None and "one" not in row["loss_labels"]:
            continue
        other = lookup.get((row["K"], row["shard_bytes"],
                            row["loss_count"], baseline))
        if other is None or row.get(metric) is None or other.get(metric) is None:
            continue
        if metric.endswith("_us"):
            ratio = float(other[metric]) / float(row[metric])
        else:
            ratio = float(row[metric]) / float(other[metric])
        ratios.append((ratio, {
            "K": row["K"], "shard_bytes": row["shard_bytes"],
            "loss_count": row["loss_count"], "speedup": ratio,
        }))
    require(ratios, f"no comparable ratios for {baseline} {metric}")
    ratios.sort(key=lambda item: item[0])
    return {
        "comparable_cell_count": len(ratios),
        "median_speedup": geometric_median([ratio for ratio, _ in ratios]),
        "faster_cell_fraction": sum(ratio > 1.0 for ratio, _ in ratios) / len(ratios),
        "slowest": ratios[0][1], "fastest": ratios[-1][1],
    }


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    output = ["| " + " | ".join(headers) + " |",
              "| " + " | ".join("---" for _ in headers) + " |"]
    output.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(output)


def generate_performance_readme(summary: Mapping[str, Any],
                                manifest: Mapping[str, Any],
                                run_metadata: Mapping[str, Any],
                                plot_directory: Path,
                                destination: Path,
                                summary_path: Path,
                                raw_bundle: Path) -> None:
    encode_rows: list[list[str]] = []
    for baseline in ("leopard1", "wirehair"):
        steady = speedup_statistics(summary, baseline, "encode_execution_us")
        encode_rows.append([
            CODEC_LABELS[baseline],
            str(steady["comparable_cell_count"]),
            f"{steady['median_speedup']:.3f}×",
            f"{100 * steady['faster_cell_fraction']:.1f}%",
        ])
    decode_rows: list[list[str]] = []
    for loss_label in LOSS_LABELS:
        for baseline in ("leopard1", "wirehair"):
            stats = speedup_statistics(
                summary, baseline, "decode_first_us", loss_label)
            decode_rows.append([
                LOSS_DISPLAY[loss_label], CODEC_LABELS[baseline],
                str(stats["comparable_cell_count"]),
                f"{stats['median_speedup']:.3f}×",
                f"{100 * stats['faster_cell_fraction']:.1f}%",
                f"K={stats['slowest']['K']}, {format_bytes(stats['slowest']['shard_bytes'])}, "
                f"{stats['slowest']['speedup']:.3f}×",
            ])
    source = run_metadata["sources"]
    binaries = run_metadata["executables"]
    host = run_metadata["host"]
    relative_plot = os.path.relpath(plot_directory, destination.parent)
    relative_raw = os.path.relpath(raw_bundle, destination.parent)
    relative_summary = os.path.relpath(summary_path, destination.parent)
    text = f"""# Leopard2 performance atlas

This atlas compares Leopard2 against the exact Leopard `main` baseline and the
shipping Wirehair codec over **120 K values**, **four shard sizes**, and **four
source-erasure regimes**. It contains {len(manifest['cells']):,} unique workload
cells and {len(summary['rows']):,} validated codec results. All graphs are
generated from the checked-in machine-readable evidence; unavailable cells are
left blank rather than interpolated.

> Values above 1× in speedup graphs mean Leopard2 is faster. These are
> single-core results on one recorded host, not universal performance claims.

## Headline graphs

### Encoding

![Full-message encode throughput]({relative_plot}/encode_execution_message_gbps.svg)

![Leopard2 encode speedup over Leopard main]({relative_plot}/encode_speedup_vs_leopard1.svg)

![Leopard2 encode speedup over Wirehair]({relative_plot}/encode_speedup_vs_wirehair.svg)

{markdown_table(
    ['Baseline', 'Comparable cells', 'Median full-message speedup',
     'Cells faster'], encode_rows)}

### Decode with setup included

![One-erasure decode throughput]({relative_plot}/decode_one_first_message_gbps.svg)

![Two-erasure decode throughput]({relative_plot}/decode_two_first_message_gbps.svg)

![Ten-percent decode throughput]({relative_plot}/decode_ten_percent_first_message_gbps.svg)

![Maximum-loss decode throughput]({relative_plot}/decode_full_first_message_gbps.svg)

{markdown_table(
    ['Loss regime', 'Baseline', 'Comparable cells', 'Median speedup',
     'Cells faster', 'Worst observed cell'], decode_rows)}

## Complete graph index

Encoding:

- [Full-message throughput]({relative_plot}/encode_execution_message_gbps.svg)
- [Parity/repair-output throughput]({relative_plot}/encode_execution_output_gbps.svg)
- [Full-message latency]({relative_plot}/encode_execution_latency_us.svg)
- [Speedup over Leopard main]({relative_plot}/encode_speedup_vs_leopard1.svg)
- [Speedup over Wirehair]({relative_plot}/encode_speedup_vs_wirehair.svg)
- [Wirehair marginal repair-output throughput (precode excluded)]({relative_plot}/wirehair_marginal_repair_output_gbps.svg)
- [Wirehair marginal repair-output latency (precode excluded)]({relative_plot}/wirehair_marginal_repair_latency_us.svg)

Each decode regime has setup-inclusive message throughput, received-input and
repaired-output throughput, latency, repeated-call/Leopard2-plan throughput,
Leopard-main speedup, Wirehair speedup, and Wirehair overhead plots:

{chr(10).join(
    f'- **{LOSS_DISPLAY[label]}:** '
    f'[setup-inclusive message throughput]({relative_plot}/decode_{label}_first_message_gbps.svg), '
    f'[repaired throughput]({relative_plot}/decode_{label}_first_repaired_gbps.svg), '
    f'[received-input throughput]({relative_plot}/decode_{label}_first_received_gbps.svg), '
    f'[latency]({relative_plot}/decode_{label}_first_latency_us.svg), '
    f'[repeated call / Leopard2 plan reuse]({relative_plot}/decode_{label}_reused_message_gbps.svg), '
    f'[vs Leopard main]({relative_plot}/decode_{label}_speedup_vs_leopard1.svg), '
    f'[vs Wirehair]({relative_plot}/decode_{label}_speedup_vs_wirehair.svg), '
    f'[Wirehair overhead]({relative_plot}/wirehair_overhead_{label}.svg)'
    for label in LOSS_LABELS)}

## Methodology and comparability

- **Code shape:** R=32. K includes every odd integer from 1 through 223,
  every power of two through 128, and endpoint 224. The maximum is deliberately
  224 so Leopard2's high-profile GF8 parent remains N=256.
- **Profiles:** Leopard2 AUTO selects low-v1 for K<32 and legacy-high-v1 for
  K≥32. The selected field is GF8 for every cell. Exact Leopard main is defined
  only for K≥32 because its public API requires R≤K. Shipping Wirehair is
  defined for K≥2, so K=1 contains Leopard2 alone.
- **Shard sizes:** 64 B, 1 KiB, 4 KiB, and 1 MiB.
- **Erasure patterns:** deterministic random original/source erasures. “10%” is
  min(32,max(1,ceil(K/10))); “full” is min(K,32). Equivalent low-K patterns are
  timed once and tagged with every applicable label.
- **Encoding:** every codec generates exactly 32 parity/repair blocks. Message
  throughput uses K×shard-bytes; output throughput uses 32×shard-bytes.
  Leopard2 and Leopard main reuse code-dependent state and perform one complete
  stripe encode. Wirehair encoder creation is message-dependent and performs
  its matrix solve/precode, so its full-message value is the jointly timed
  create-plus-32-repair emission. Its much smaller marginal repair-row emission
  is a separately labeled diagnostic, never the headline comparison. No
  synthetic "first use" number is formed by adding separately measured medians.
- **Decode:** setup-inclusive plots compare Leopard2's public one-shot decode,
  exact Leopard main's decode (which includes its setup), and a fresh Wirehair
  decoder plus internal allocation, ordered surviving-source/repair ingestion,
  matrix solve, and recovery. Reusable-plan plots compare
  Leopard2 execution with exact main; Wirehair has no equivalent immutable
  erasure-pattern plan and is omitted there. Exact Leopard main is shown for
  context in those repeated-call plots, but each of its calls still includes
  locator/setup work; only Leopard2 reuses a byte-independent plan.
- **Decode rate denominators:** “source-message-equivalent” is K×shard-bytes
  divided by latency and normalizes one logical message across loss counts.
  Repaired-output uses missing-count×shard-bytes. Separate received-input plots
  use every non-null shard offered to Leopard (surviving originals plus all R
  parity shards) and the source/repair blocks actually consumed by Wirehair.
- **Wirehair caveat:** Wirehair is a fountain code, not the same systematic MDS
  Reed–Solomon wire profile. It may require repair overhead beyond the number of
  missing sources; that overhead is graphed separately. Its parity bytes are not
  expected to match Leopard.
- **Buffer shape:** every atlas message is exactly K full blocks; this harness
  does not exercise Wirehair's partial final-message-block behavior.
- **Correctness:** every invocation round-trips. Leopard2 parity matches exact
  Leopard main byte-for-byte wherever main is defined. All codecs use identical
  source bytes, missing-original indices, and recovered-original fingerprints.
- **ISA:** Leopard2 and exact Leopard main use an attested AVX2 ceiling, and
  Leopard2 AUTO must resolve to AVX2 in every cell. The pinned shipping
  Wirehair artifact contains target-qualified AVX-512 helpers that cannot be
  compiled out at this revision. Those helpers are reachable only through its
  opt-in thread-wide XOR path; the adapter forces that path off before and after
  every measured phase and records `measured_path_avx512=false`. Thus the
  measured Wirehair-v1 compact path is AVX2, while its artifact is not falsely
  described as AVX2-only.
- **Isolation:** a single pinned physical CPU, one thread, the canonical project
  lock, immutable executable copies, and a 2 GiB address-space cap. Runs are
  serialized to avoid both timing contamination and prior OOM failure modes.
- **Statistics:** each timing value is the median of {manifest['cells'][0]['iterations']}
  samples for small/medium blocks and 5 samples at 1 MiB, after warmup. Reuse is
  chosen deterministically to target about 8 MiB of work per sample and is
  reported in each raw result. Summary-table speedups are log-medians: the
  exponential of the median log-ratio across comparable cells.

## Evidence identity

Host `{host['hostname']}`, kernel `{host['kernel']}`, pinned CPU
`{host['benchmark_cpu']}`, allowed CPU set `{host['allowed_cpus']}`.

{markdown_table(
    ['Component', 'Source commit', 'Source tree', 'Executable SHA-256'], [
      ['Leopard2', source['leopard2']['commit'], source['leopard2']['tree'],
       binaries['leopard2']['sha256']],
      ['Leopard main', source['leopard1']['commit'], source['leopard1']['tree'],
       binaries['leopard1']['sha256']],
      ['Wirehair shipping codec', source['wirehair']['commit'],
       source['wirehair']['tree'], binaries['wirehair']['sha256']],
    ])}

Evidence files:

- [normalized summary]({relative_summary})
- [complete raw bundle]({relative_raw})
- `manifest.json` — exact deterministic workload matrix
- `run_metadata.json` — host, source, lock, runner, and executable identities

## Reproduction

The complete command, including the three required executable SHA-256 values,
is retained in `REPRODUCE.txt` next to this README. The core workflow is:

    python3 experiments/leopard2/performance_atlas/test_generate_atlas.py -v
    python3 experiments/leopard2/performance_atlas/generate_atlas.py all \\
      --output <result-directory> \\
      --leopard2-build-root <clean-build> \\
      --leopard2 <clean-build>/bench_leopard2_allk --leopard2-sha256 <sha256> \\
      --leopard1-build-root <main-build> \\
      --leopard1 <main-build>/leopard_main_benchmark --leopard1-sha256 <sha256> \\
      --wirehair-build-root <wirehair-build> \\
      --wirehair <wirehair-build>/wirehair_v1_benchmark --wirehair-sha256 <sha256> \\
      --leopard-source <clean-leopard2-source> \\
      --main-source <exact-main-source> \\
      --wirehair-source <pinned-wirehair-source>

The runner is resumable: each invocation is an atomic JSON file, and every
existing record is fully revalidated before it is reused.
"""
    atomic_write_bytes(destination, text.encode("utf-8"))


def self_test() -> None:
    values = k_values()
    require(values[0] == 1 and values[-1] == 224 and len(values) == 120,
            "K-grid self-test failed")
    require(all(k % 2 == 1 or k & (k - 1) == 0 or k == 224
                for k in values), "K-grid admitted an unintended even value")
    require(loss_map(1) == {1: ["one", "ten_percent", "full"]},
            "K=1 loss de-duplication failed")
    require(loss_map(2) == {1: ["one", "ten_percent"],
                            2: ["two", "full"]},
            "K=2 loss de-duplication failed")
    require(conceptual_loss(101, "ten_percent") == 11,
            "10-percent rounding failed")
    require(conceptual_loss(224, "full") == 32,
            "full-loss clipping failed")
    require(available_codecs(1) == ("leopard2",),
            "K=1 availability failed")
    require(available_codecs(31) == ("leopard2", "wirehair"),
            "K=31 availability failed")
    require(available_codecs(32) == CODECS,
            "K=32 availability failed")
    manifest = build_manifest()
    require(manifest["matrix"]["unique_cell_count"] ==
            len(build_expected_cells_without_recursion()),
            "manifest count failed")
    mutated = json.loads(json.dumps(manifest))
    mutated["cells"][0]["seed"] ^= 1
    try:
        validate_manifest(mutated)
    except AtlasError:
        pass
    else:
        raise AtlasError("manifest validator accepted a changed seed")
    require(codec_order("k032_b0000064_l02") ==
            codec_order("k032_b0000064_l02"),
            "codec order is not deterministic")
    print(f"PASS: atlas self-test ({len(values)} K values, "
          f"{len(manifest['cells'])} unique cells)")


def parse_args(arguments: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("manifest", "run", "summarize",
                                         "plot", "all", "self-test"))
    parser.add_argument("--output", type=Path, default=Path(
        ".research/leopard2/performance-atlas-run"))
    parser.add_argument("--leopard2", type=Path)
    parser.add_argument("--leopard2-build-root", type=Path)
    parser.add_argument("--leopard1", type=Path)
    parser.add_argument("--leopard1-build-root", type=Path)
    parser.add_argument("--wirehair", type=Path)
    parser.add_argument("--wirehair-build-root", type=Path)
    parser.add_argument("--leopard2-sha256")
    parser.add_argument("--leopard1-sha256")
    parser.add_argument("--wirehair-sha256")
    parser.add_argument("--leopard-source", type=Path,
                        default=Path.cwd())
    parser.add_argument("--main-source", type=Path)
    parser.add_argument("--wirehair-source", type=Path,
        default=Path(".research/leopard2/wirehair-v1"))
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--address-space-limit", type=int,
                        default=DEFAULT_AS_LIMIT)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--max-cells", type=int,
                        help="bounded smoke/resume run; never authoritative")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--raw-bundle", type=Path)
    parser.add_argument("--plots", type=Path)
    parser.add_argument("--performance-readme", type=Path)
    args = parser.parse_args(arguments)
    if args.mode in ("run", "all"):
        required = ("leopard2", "leopard2_build_root",
                    "leopard1", "leopard1_build_root",
                    "wirehair", "wirehair_build_root",
                    "leopard2_sha256", "leopard1_sha256", "wirehair_sha256",
                    "main_source")
        for name in required:
            require(getattr(args, name) is not None,
                    f"--{name.replace('_', '-')} is required for {args.mode}")
        for name in ("leopard2_sha256", "leopard1_sha256", "wirehair_sha256"):
            value = getattr(args, name)
            require(len(value) == 64 and all(character in "0123456789abcdef"
                    for character in value), f"--{name.replace('_', '-')} is invalid")
    require(args.address_space_limit >= 256 * 1024 * 1024,
            "address-space limit is implausibly small")
    require(args.timeout > 0, "timeout must be positive")
    return args


def main(arguments: Sequence[str]) -> int:
    try:
        args = parse_args(arguments)
        if args.mode == "self-test":
            self_test()
            return 0
        manifest = build_manifest()
        if args.mode == "manifest":
            atomic_write_json(args.output, manifest)
            return 0
        if args.mode in ("run", "all"):
            run_campaign(args, manifest)
        if args.mode in ("summarize", "all"):
            summary = build_summary(args.output.resolve(), manifest,
                                    allow_partial=args.allow_partial)
            summary_path = args.summary or args.output / "summary.json"
            atomic_write_json(summary_path, summary)
            if summary["complete"]:
                write_raw_bundle(args.output.resolve(), manifest,
                                 args.raw_bundle)
        if args.mode in ("plot", "all"):
            summary_path = args.summary or args.output / "summary.json"
            summary = read_json(summary_path)
            validate_summary(summary, manifest,
                             require_complete=not args.allow_partial)
            plots = args.plots or args.output / "plots"
            generated = generate_plots(summary, plots)
            print(f"generated {len(generated)} SVG plots in {plots}")
            if summary.get("complete") is True:
                metadata = read_json(args.output / "run_metadata.json")
                raw_bundle = args.raw_bundle or args.output / "raw_bundle.json"
                require(raw_bundle.exists(), "complete raw bundle is missing")
                readme = args.performance_readme or args.output / "README_PERFORMANCE.md"
                generate_performance_readme(
                    summary, manifest, metadata, plots, readme, summary_path,
                    raw_bundle)
        return 0
    except (AtlasError, OSError, subprocess.SubprocessError) as error:
        print(f"performance atlas: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
