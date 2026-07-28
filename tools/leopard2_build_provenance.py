#!/usr/bin/env python3
"""Fail-closed provenance for production Leopard2 benchmark builds.

The benchmark's embedded source attestation covers only the benchmark
translation unit.  It cannot prove which objects were archived into Leopard2.
This module binds a benchmark to its exact Release CMake cache, compiler
commands, source/object pairs, archive recipe/members, link recipe, and output
bytes.  The final-map and all-K runners both use the same closure.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import stat
import subprocess
import tempfile
from typing import Any, Mapping, Sequence


MAX_FILE_BYTES = 256 * 1024 * 1024
MAX_METADATA_BYTES = 16 * 1024 * 1024

CORE_LIBRARY_SOURCES = {
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp",
    "LeopardCommon.cpp",
}

GIT_ENVIRONMENT = {
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
    "LANG": "C",
    "LC_ALL": "C",
    "PATH": "/usr/bin:/bin",
}


class BuildProvenanceError(RuntimeError):
    """The candidate build is not valid benchmark evidence."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BuildProvenanceError(message)


def _stable_fields(status: os.stat_result) -> tuple[int, ...]:
    return (
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
    )


def file_snapshot(
    path: Path | str, label: str, *, maximum_bytes: int = MAX_FILE_BYTES,
) -> tuple[dict[str, Any], bytes]:
    """Read one regular file once and bind the returned bytes to its pathname."""
    resolved = Path(path).resolve(strict=True)
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode), f"{label} is not a regular file")
        require(0 <= before.st_size <= maximum_bytes,
                f"{label} exceeds its retained byte bound")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            require(bool(block), f"{label} ended before its recorded size")
            chunks.append(block)
            remaining -= len(block)
        require(not os.read(descriptor, 1),
                f"{label} grew beyond its recorded size while read")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    require(_stable_fields(before) == _stable_fields(after),
            f"{label} changed while read")
    path_status = resolved.stat()
    require(_stable_fields(after) == _stable_fields(path_status),
            f"{label} pathname changed while read")
    content = b"".join(chunks)
    identity = {
        "path": str(resolved),
        "sha256": hashlib.sha256(content).hexdigest(),
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": after.st_mode,
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
    }
    return identity, content


def file_identity(
    path: Path | str, label: str, *, maximum_bytes: int = MAX_FILE_BYTES,
) -> dict[str, Any]:
    return file_snapshot(path, label, maximum_bytes=maximum_bytes)[0]


def _run(
    command: Sequence[str], label: str, *, maximum_bytes: int = 4 << 20,
    timeout: float = 120,
) -> bytes:
    completed = subprocess.run(
        list(command), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=GIT_ENVIRONMENT, timeout=timeout, check=False)
    require(completed.returncode == 0,
            f"{label} failed with rc={completed.returncode}: "
            f"{completed.stderr.decode(errors='replace').strip()}")
    require(len(completed.stdout) <= maximum_bytes,
            f"{label} output exceeds its retained byte bound")
    return completed.stdout


def parse_cmake_cache(content: bytes) -> dict[str, str]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError("CMakeCache.txt is not strict UTF-8") from error
    result: dict[str, str] = {}
    for line in text.splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name_and_type, value = line.split("=", 1)
        if ":" not in name_and_type:
            continue
        name, unused_type = name_and_type.split(":", 1)
        require(name not in result,
                f"CMakeCache.txt contains duplicate key {name}")
        result[name] = value
    return result


def _cmake_true(value: str | None) -> bool:
    return (value or "").upper() in {"1", "ON", "TRUE", "YES", "Y"}


def _cmake_false(value: str | None) -> bool:
    return (value or "").upper() in {
        "0", "OFF", "FALSE", "NO", "N", "IGNORE", "NOTFOUND", "",
    } or (value or "").upper().endswith("-NOTFOUND")


def _resolve_build_operand(build: Path, directory: Path, value: str) -> Path:
    operand = Path(value)
    if not operand.is_absolute():
        operand = directory / operand
    resolved = operand.resolve(strict=False)
    require(resolved.is_relative_to(build),
            f"build operand escapes the candidate build: {value}")
    return resolved


def _recipe_objects(content: bytes, build: Path, label: str) -> list[Path]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError(f"{label} is not strict UTF-8") from error
    objects: list[Path] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError as error:
            raise BuildProvenanceError(f"cannot parse {label}: {error}") from error
        require(tokens and all("@" not in token for token in tokens),
                f"{label} contains an empty command or response file")
        for token in tokens:
            if token.endswith((".o", ".obj")):
                objects.append(_resolve_build_operand(build, build, token))
    require(objects and len(objects) == len(set(objects)),
            f"{label} object closure is empty or contains duplicates")
    return objects


def _compile_tokens(entry: Mapping[str, Any]) -> list[str]:
    has_arguments = "arguments" in entry
    has_command = "command" in entry
    require(has_arguments != has_command,
            "compile command must have exactly one argv representation")
    if has_arguments:
        value = entry["arguments"]
        require(isinstance(value, list) and value and
                all(isinstance(item, str) and item for item in value),
                "compile command arguments are malformed")
        return list(value)
    value = entry["command"]
    require(isinstance(value, str) and value,
            "compile command string is malformed")
    try:
        return shlex.split(value, posix=True)
    except ValueError as error:
        raise BuildProvenanceError(
            f"cannot parse compile command: {error}") from error


def _compile_output(
    entry: Mapping[str, Any], tokens: Sequence[str], build: Path,
) -> Path:
    directory_value = entry.get("directory")
    require(isinstance(directory_value, str) and directory_value,
            "compile command has no directory")
    directory = Path(directory_value).resolve(strict=True)
    require(directory == build,
            "compile command did not run at the exact candidate build root")
    output_positions = [index for index, token in enumerate(tokens)
                        if token == "-o"]
    require(len(output_positions) == 1 and
            output_positions[0] + 1 < len(tokens),
            "compile command has no unique output operand")
    output_value = tokens[output_positions[0] + 1]
    retained_output = entry.get("output")
    require(isinstance(retained_output, str) and
            retained_output == output_value,
            "compile command output metadata differs from its argv")
    return _resolve_build_operand(build, directory, output_value)


def _validate_compile_flags(tokens: Sequence[str], source: Path) -> str:
    require(tokens.count("-c") == 1 and "-O3" in tokens,
            f"compile command is not an optimized compile for {source.name}")
    require(all("@" not in token for token in tokens),
            f"compile command uses an unbound response file for {source.name}")
    forbidden_common = (
        "-march=native", "-mtune=native", "-Ofast", "-ffast-math",
        "-flto", "-fprofile-generate", "-fprofile-use",
        "-fsanitize=", "-O0", "-O1", "-O2", "-Os", "-Og",
    )
    require(not any(token == item or token.startswith(item)
                    for token in tokens for item in forbidden_common),
            f"compile command has a non-canonical performance flag for "
            f"{source.name}")
    positive_avx512 = [token for token in tokens
                       if token.startswith("-mavx512")]
    require(not positive_avx512,
            f"pure-AVX2 build compiled {source.name} with AVX-512")
    name = source.name
    if name.startswith("Leopard2BackendGFNI"):
        # The GFNI member is the AVX2 algorithms recompiled with the affine
        # kernels.  Its data path stays VEX 256-bit, so it is admissible in a
        # pure-AVX2 build, but it must carry the affine flag explicitly and
        # remain EVEX-free.
        require("-mavx2" in tokens and "-mgfni" in tokens and
                "-mno-avx512f" in tokens,
                f"GFNI source lacks exact AVX2/GFNI/no-AVX512 flags: {name}")
        return "avx2-gfni-no-avx512"
    require("-mgfni" not in tokens,
            f"non-GFNI source compiled with the affine ISA: {name}")
    if name.startswith("Leopard2BackendAVX2"):
        require("-mavx2" in tokens and "-mno-avx512f" in tokens,
                f"AVX2 source lacks exact AVX2/no-AVX512 flags: {name}")
        return "avx2-no-avx512"
    if name == "Leopard2BackendSSSE3.cpp":
        require("-mssse3" in tokens and "-mno-avx" in tokens and
                not any(token.startswith("-mavx") for token in tokens),
                "SSSE3 source lacks exact SSSE3/no-AVX flags")
        return "ssse3-no-avx"
    require(not any(token.startswith(("-mavx", "-mssse3"))
                    for token in tokens),
            f"portable source has a target-local ISA flag: {name}")
    return "portable-core"


def _tracked_files(source: Path) -> set[Path]:
    git = Path("/usr/bin/git").resolve(strict=True)
    raw = _run((str(git), "-C", str(source), "ls-files", "-z"),
               "candidate tracked-file inventory", maximum_bytes=32 << 20)
    result = set()
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            relative = item.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise BuildProvenanceError(
                "candidate tracked path is not strict UTF-8") from error
        path = (source / relative).resolve(strict=True)
        require(path.is_relative_to(source),
                "candidate tracked path escapes the source root")
        result.add(path)
    require(result, "candidate tracked-file inventory is empty")
    return result


def _expected_library_sources(
    source: Path, cache: Mapping[str, str], tracked: set[Path],
) -> set[Path]:
    names = set(CORE_LIBRARY_SOURCES)
    names.add("LeopardFF8.cpp")
    if _cmake_true(cache.get("LEOPARD_ENABLE_GF16")):
        names.add("LeopardFF16.cpp")
    names.add("Leopard2BackendSSSE3.cpp")
    avx2_sources = {
        path.name for path in tracked
        if path.parent == source and
        path.name.startswith("Leopard2BackendAVX2") and
        path.suffix == ".cpp"
    }
    require(avx2_sources,
            "candidate source contains no production AVX2 translation unit")
    names.update(avx2_sources)
    # The GFNI member joins the archive whenever its compiler probe passes.
    # It contains no AVX-512 encoding, so a pure-AVX2 candidate retains it.
    if _cmake_true(cache.get("LEO2_FLAG_MGFNI")):
        names.add("Leopard2BackendGFNI.cpp")
    expected = {(source / name).resolve(strict=True) for name in names}
    require(expected.issubset(tracked),
            "candidate build expects a library source not tracked at HEAD")
    return expected


def candidate_build_provenance(
    build_root: Path | str,
    source_root: Path | str,
    executable: Path | str,
    executable_target: str,
) -> dict[str, Any]:
    """Capture and validate one pure-AVX2 production benchmark build."""
    build = Path(build_root).resolve(strict=True)
    source = Path(source_root).resolve(strict=True)
    require(build.is_dir(), "candidate build root is not a directory")
    require(source.is_dir(), "candidate source root is not a directory")
    require(executable_target in {"bench_leopard2", "bench_leopard2_allk"},
            "unsupported candidate benchmark target")

    expected_executable = (build / executable_target).resolve(strict=True)
    requested_executable = Path(executable).resolve(strict=True)
    require(requested_executable == expected_executable,
            "candidate executable is not the declared CMake build target")
    archive = (build / "libleopard.a").resolve(strict=True)
    cache_path = build / "CMakeCache.txt"
    commands_path = build / "compile_commands.json"
    archive_recipe_path = build / "CMakeFiles/leopard.dir/link.txt"
    executable_recipe_path = \
        build / f"CMakeFiles/{executable_target}.dir/link.txt"

    cache_identity, cache_bytes = file_snapshot(
        cache_path, "candidate CMake cache", maximum_bytes=MAX_METADATA_BYTES)
    cache = parse_cmake_cache(cache_bytes)
    required_exact = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "ENABLE_OPENMP": "ON",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
    }
    for name, expected in required_exact.items():
        require(cache.get(name) == expected,
                f"candidate CMake cache {name} is {cache.get(name)!r}, "
                f"expected {expected!r}")
    if executable_target == "bench_leopard2_allk":
        require(cache.get("LEO2_BUILD_ALLK_DIAGNOSTIC") == "ON",
                "all-K candidate target was not explicitly enabled")
    variant = cache.get("LEO2_BACKEND_VARIANT")
    require(variant in {"auto", "avx2"},
            "candidate backend variant is not production auto or AVX2")
    require(_cmake_true(cache.get("LEO2_FLAG_MAVX2")) and
            _cmake_true(cache.get("LEO2_FLAG_MNO_AVX512F")),
            "candidate compiler lacks the required pure-AVX2 flags")
    for name in ("LEO2_FLAG_MAVX512F", "LEO2_FLAG_MAVX512BW",
                 "LEO2_FLAG_MAVX512VL"):
        require(_cmake_false(cache.get(name)),
                f"candidate pure-AVX2 build did not disable {name}")
    home = Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve(strict=True)
    require(home == source,
            "candidate CMake cache points at another source root")
    global_flags = " ".join((cache.get("CMAKE_CXX_FLAGS", ""),
                             cache.get("CMAKE_CXX_FLAGS_RELEASE", "")))
    require(cache.get("CMAKE_CXX_FLAGS") == "" and
            cache.get("CMAKE_CXX_FLAGS_RELEASE") == "-O3 -DNDEBUG",
            "candidate CMake flags differ from the canonical Release profile")
    for name in ("CMAKE_CXX_COMPILER_LAUNCHER", "CMAKE_CXX_COMPILER_ARG1",
                 "CMAKE_CXX_COMPILER_TARGET", "CMAKE_SYSROOT",
                 "CMAKE_TOOLCHAIN_FILE"):
        require(not cache.get(name),
                f"candidate CMake cache uses unsupported {name}")
    require(_cmake_false(cache.get("CMAKE_POSITION_INDEPENDENT_CODE")),
            "candidate CMake cache changes position-independent code policy")
    require(not any(flag in global_flags for flag in (
                "-march=native", "-mtune=native", "-mavx", "-mssse3",
                "-mavx512")),
            "candidate CMake cache raises the global ISA floor")

    compiler_invocation = cache.get("CMAKE_CXX_COMPILER", "")
    require(bool(compiler_invocation),
            "candidate CMake cache has no C++ compiler")
    compiler = Path(compiler_invocation).resolve(strict=True)
    compiler_identity = file_identity(compiler, "candidate C++ compiler")
    compiler_version = _run((str(compiler), "--version"),
                            "candidate compiler version")

    commands_identity, commands_bytes = file_snapshot(
        commands_path, "candidate compile commands",
        maximum_bytes=MAX_METADATA_BYTES)
    try:
        entries = json.loads(commands_bytes.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise BuildProvenanceError(
            f"candidate compile_commands.json is invalid: {error}") from error
    require(isinstance(entries, list) and entries,
            "candidate compile command closure is empty")
    by_output: dict[Path, tuple[Mapping[str, Any], list[str], Path]] = {}
    for entry in entries:
        require(isinstance(entry, dict) and
                isinstance(entry.get("file"), str),
                "candidate compile command entry is malformed")
        tokens = _compile_tokens(entry)
        require(tokens and
                Path(tokens[0]).resolve(strict=True) == compiler,
                "candidate compile command uses another compiler")
        source_operand = Path(entry["file"]).resolve(strict=True)
        compile_positions = [index for index, token in enumerate(tokens)
                             if token == "-c"]
        require(len(compile_positions) == 1 and
                compile_positions[0] + 1 < len(tokens) and
                Path(tokens[compile_positions[0] + 1]).resolve(strict=True) ==
                source_operand,
                "candidate compile command source metadata differs from argv")
        output = _compile_output(entry, tokens, build)
        require(output not in by_output,
                f"duplicate candidate compile output: {output}")
        by_output[output] = (entry, tokens, source_operand)

    archive_recipe_identity, archive_recipe_bytes = file_snapshot(
        archive_recipe_path, "candidate archive link recipe",
        maximum_bytes=MAX_METADATA_BYTES)
    executable_recipe_identity, executable_recipe_bytes = file_snapshot(
        executable_recipe_path, "candidate executable link recipe",
        maximum_bytes=MAX_METADATA_BYTES)
    archive_objects = _recipe_objects(
        archive_recipe_bytes, build, "candidate archive link recipe")
    executable_objects = _recipe_objects(
        executable_recipe_bytes, build, "candidate executable link recipe")
    require(len(executable_objects) == 1,
            "candidate executable does not have one benchmark object")

    tracked = _tracked_files(source)
    expected_library_sources = _expected_library_sources(
        source, cache, tracked)
    closure_records: list[dict[str, Any]] = []
    archive_sources: set[Path] = set()
    for role, objects in (("archive", archive_objects),
                          ("benchmark", executable_objects)):
        for obj in objects:
            require(obj in by_output,
                    f"{role} object has no exact compile command: {obj}")
            entry, tokens, source_operand = by_output[obj]
            require(source_operand.is_relative_to(source) and
                    source_operand in tracked,
                    f"{role} object source is not tracked at candidate HEAD")
            source_identity = file_identity(
                source_operand, f"candidate {role} source")
            object_identity = file_identity(obj, f"candidate {role} object")
            require(object_identity["mtime_ns"] >= source_identity["mtime_ns"],
                    f"candidate object predates source: {source_operand.name}")
            flag_profile = _validate_compile_flags(tokens, source_operand)
            if role == "archive":
                archive_sources.add(source_operand)
            else:
                require(source_operand ==
                        (source / "bench/leopard2/benchmark.cpp").resolve(
                            strict=True),
                        "candidate executable was not compiled from the exact "
                        "Leopard2 benchmark source")
            closure_records.append({
                "role": role,
                "source": source_identity,
                "object": object_identity,
                "compile_entry": {
                    "directory": entry["directory"],
                    "file": entry["file"],
                    "output": entry["output"],
                    "arguments": tokens,
                },
                "flag_profile": flag_profile,
            })
    require(archive_sources == expected_library_sources,
            "candidate archive source closure differs from the production "
            "pure-AVX2 source set")

    archive_identity = file_identity(archive, "candidate production archive")
    executable_identity = file_identity(
        expected_executable, "candidate benchmark executable")
    require(stat.S_ISREG(executable_identity["mode"]) and
            executable_identity["mode"] & 0o111,
            "candidate benchmark is not executable")
    archive_records = [record for record in closure_records
                       if record["role"] == "archive"]
    benchmark_record = next(record for record in closure_records
                            if record["role"] == "benchmark")
    require(all(archive_identity["mtime_ns"] >=
                record["object"]["mtime_ns"] for record in archive_records),
            "candidate archive predates one of its exact objects")
    require(executable_identity["mtime_ns"] >= archive_identity["mtime_ns"] and
            executable_identity["mtime_ns"] >=
            benchmark_record["object"]["mtime_ns"],
            "candidate benchmark predates one of its exact link inputs")

    ar_invocation = cache.get("CMAKE_AR", "")
    require(bool(ar_invocation), "candidate CMake cache has no archiver")
    archiver = Path(ar_invocation).resolve(strict=True)
    archiver_identity = file_identity(archiver, "candidate archiver")
    members = _run((str(archiver), "t", str(archive)),
                   "candidate archive inventory").decode(
                       "utf-8", errors="strict").splitlines()
    expected_members = [path.name for path in archive_objects]
    require(members == expected_members and len(members) == len(set(members)),
            "candidate archive members differ from its exact object recipe")
    object_identity_by_name = {
        Path(record["object"]["path"]).name: record["object"]
        for record in archive_records
    }
    require(set(object_identity_by_name) == set(members),
            "candidate archive object basenames are ambiguous")
    archive_member_identities = []
    for member in members:
        content = _run(
            (str(archiver), "p", str(archive), member),
            f"candidate archive member {member}",
            maximum_bytes=MAX_FILE_BYTES)
        identity = {
            "member": member, "size": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
        expected_object = object_identity_by_name[member]
        require(identity["size"] == expected_object["size"] and
                identity["sha256"] == expected_object["sha256"],
                f"candidate archive member bytes differ from object: {member}")
        archive_member_identities.append(identity)

    try:
        executable_recipe_text = executable_recipe_bytes.decode(
            "utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError(
            "candidate executable recipe is not strict UTF-8") from error
    link_tokens = []
    for line in executable_recipe_text.splitlines():
        if line.strip():
            link_tokens.extend(shlex.split(line, posix=True))
    require(str(archive) in link_tokens or archive.name in link_tokens,
            "candidate executable recipe does not link the validated archive")
    require(not any("@" in token for token in link_tokens),
            "candidate executable recipe uses an unbound response file")

    return {
        "schema": "leopard2-production-build-closure/v1",
        "build_root": str(build),
        "source_root": str(source),
        "executable_target": executable_target,
        "validated_cache": {
            **required_exact,
            "LEO2_BUILD_ALLK_DIAGNOSTIC": cache.get(
                "LEO2_BUILD_ALLK_DIAGNOSTIC"),
            "LEO2_BUILD_TESTS": cache.get("LEO2_BUILD_TESTS"),
            "ENABLE_OPENMP": cache.get("ENABLE_OPENMP"),
            "LEO2_BACKEND_VARIANT": variant,
            "LEOPARD_ENABLE_GF16": cache.get("LEOPARD_ENABLE_GF16"),
            "LEO2_FLAG_MAVX2": cache.get("LEO2_FLAG_MAVX2"),
            "LEO2_FLAG_MNO_AVX512F": cache.get("LEO2_FLAG_MNO_AVX512F"),
            "LEO2_FLAG_MAVX512F": cache.get("LEO2_FLAG_MAVX512F"),
            "LEO2_FLAG_MAVX512BW": cache.get("LEO2_FLAG_MAVX512BW"),
            "LEO2_FLAG_MAVX512VL": cache.get("LEO2_FLAG_MAVX512VL"),
            "CMAKE_CXX_FLAGS": cache.get("CMAKE_CXX_FLAGS"),
            "CMAKE_CXX_FLAGS_RELEASE": cache.get("CMAKE_CXX_FLAGS_RELEASE"),
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
        },
        "cmake_cache": cache_identity,
        "compile_commands": commands_identity,
        "archive_link_recipe": archive_recipe_identity,
        "executable_link_recipe": executable_recipe_identity,
        "compiler": compiler_identity,
        "compiler_version_sha256": hashlib.sha256(compiler_version).hexdigest(),
        "archiver": archiver_identity,
        "archive_members": members,
        "archive_member_identities": archive_member_identities,
        "archive": archive_identity,
        "executable": executable_identity,
        "source_object_compile_closure": sorted(
            closure_records, key=lambda record: (
                record["role"], record["source"]["path"],
                record["object"]["path"])),
    }


def compare_reproducible_builds(
    candidate: Mapping[str, Any], rebuilt: Mapping[str, Any],
) -> dict[str, Any]:
    """Require a clean rebuild to reproduce every linked project byte."""
    require(candidate.get("schema") ==
            "leopard2-production-build-closure/v1" and
            rebuilt.get("schema") ==
            "leopard2-production-build-closure/v1",
            "candidate/rebuild provenance schema differs")
    require(candidate.get("source_root") == rebuilt.get("source_root") and
            candidate.get("executable_target") ==
            rebuilt.get("executable_target"),
            "candidate/rebuild source or target differs")
    source = Path(str(candidate["source_root"]))

    def closure_map(record: Mapping[str, Any]) -> dict[tuple[str, str], Any]:
        result = {}
        closure = record.get("source_object_compile_closure")
        require(isinstance(closure, list) and closure,
                "candidate/rebuild compile closure is empty")
        for item in closure:
            require(isinstance(item, dict) and
                    isinstance(item.get("role"), str) and
                    isinstance(item.get("source"), dict) and
                    isinstance(item.get("object"), dict),
                    "candidate/rebuild compile closure is malformed")
            source_path = Path(str(item["source"].get("path", "")))
            require(source_path.is_relative_to(source),
                    "candidate/rebuild source escapes the bound source root")
            key = (item["role"], source_path.relative_to(source).as_posix())
            require(key not in result,
                    "candidate/rebuild compile closure contains duplicates")
            result[key] = item
        return result

    candidate_closure = closure_map(candidate)
    rebuilt_closure = closure_map(rebuilt)
    require(set(candidate_closure) == set(rebuilt_closure),
            "clean rebuild source/object closure differs")
    object_comparisons = []
    for key in sorted(candidate_closure):
        original = candidate_closure[key]
        reproduction = rebuilt_closure[key]
        require(original["source"]["sha256"] ==
                reproduction["source"]["sha256"] and
                original["source"]["size"] ==
                reproduction["source"]["size"],
                f"clean rebuild source bytes differ: {key[1]}")
        require(original["flag_profile"] == reproduction["flag_profile"],
                f"clean rebuild ISA flag profile differs: {key[1]}")
        require(original["object"]["sha256"] ==
                reproduction["object"]["sha256"] and
                original["object"]["size"] ==
                reproduction["object"]["size"],
                f"clean rebuild object bytes differ: {key[1]}")
        object_comparisons.append({
            "role": key[0], "source": key[1],
            "sha256": original["object"]["sha256"],
            "size": original["object"]["size"],
            "flag_profile": original["flag_profile"],
        })
    require(candidate.get("archive_member_identities") ==
            rebuilt.get("archive_member_identities"),
            "clean rebuild archive member bytes differ")
    for artifact in ("archive", "executable"):
        original = candidate.get(artifact)
        reproduction = rebuilt.get(artifact)
        require(isinstance(original, dict) and
                isinstance(reproduction, dict) and
                original.get("sha256") == reproduction.get("sha256") and
                original.get("size") == reproduction.get("size"),
                f"clean rebuild {artifact} bytes differ")
    require(candidate.get("compiler", {}).get("sha256") ==
            rebuilt.get("compiler", {}).get("sha256") and
            candidate.get("compiler_version_sha256") ==
            rebuilt.get("compiler_version_sha256"),
            "clean rebuild compiler identity differs")
    return {
        "schema": "leopard2-reproducible-build-proof/v1",
        "method": "runner-owned-empty-directory-configure-build-byte-compare",
        "source_root": str(source),
        "executable_target": candidate["executable_target"],
        "compiler_sha256": candidate["compiler"]["sha256"],
        "objects": object_comparisons,
        "archive_members": candidate["archive_member_identities"],
        "archive_sha256": candidate["archive"]["sha256"],
        "executable_sha256": candidate["executable"]["sha256"],
    }


def verify_reproducible_candidate_build(
    candidate: Mapping[str, Any], *, jobs: int | None = None,
) -> dict[str, Any]:
    """Reconfigure/build in an empty directory and compare all linked bytes."""
    source = Path(str(candidate.get("source_root", ""))).resolve(strict=True)
    target = str(candidate.get("executable_target", ""))
    cache = candidate.get("validated_cache")
    require(isinstance(cache, dict),
            "candidate build has no retained CMake semantics")
    cmake = Path("/usr/bin/cmake").resolve(strict=True)
    if jobs is None:
        jobs = min(128, len(os.sched_getaffinity(0)))
    require(type(jobs) is int and 1 <= jobs <= 128,
            "clean rebuild job count must be in 1..128")

    temporary = Path(tempfile.mkdtemp(prefix="leopard2-proof-build-"))
    try:
        configure = [
            str(cmake), "-S", str(source), "-B", str(temporary),
            "-G", "Unix Makefiles",
        ]
        cmake_values = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_CXX_FLAGS": cache.get("CMAKE_CXX_FLAGS"),
            "CMAKE_CXX_FLAGS_RELEASE": cache.get("CMAKE_CXX_FLAGS_RELEASE"),
            "ENABLE_OPENMP": cache.get("ENABLE_OPENMP"),
            "LEO2_BACKEND_VARIANT": cache.get("LEO2_BACKEND_VARIANT"),
            "LEO2_BUILD_ALLK_DIAGNOSTIC": cache.get(
                "LEO2_BUILD_ALLK_DIAGNOSTIC"),
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": cache.get("LEO2_BUILD_TESTS"),
            "LEO2_ENABLE_CUDA": "OFF",
            "LEOPARD_ENABLE_GF8": "ON",
            "LEOPARD_ENABLE_GF16": cache.get("LEOPARD_ENABLE_GF16"),
            "LEO2_FLAG_MAVX2": cache.get("LEO2_FLAG_MAVX2"),
            "LEO2_FLAG_MNO_AVX512F": cache.get("LEO2_FLAG_MNO_AVX512F"),
            "LEO2_FLAG_MAVX512F": "FALSE",
            "LEO2_FLAG_MAVX512BW": "FALSE",
            "LEO2_FLAG_MAVX512VL": "FALSE",
        }
        require(all(value is not None for value in cmake_values.values()),
                "candidate build lacks a CMake value needed for clean rebuild")
        configure.extend(
            f"-D{name}={value}" for name, value in cmake_values.items())
        _run(configure, "runner-owned candidate configure",
             maximum_bytes=16 << 20, timeout=300)
        _run((str(cmake), "--build", str(temporary), "--target", target,
              "--parallel", str(jobs)),
             "runner-owned candidate build", maximum_bytes=32 << 20,
             timeout=1800)
        rebuilt = candidate_build_provenance(
            temporary, source, temporary / target, target)
        return compare_reproducible_builds(candidate, rebuilt)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
