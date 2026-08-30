#!/usr/bin/env python3
"""Linux adversarial tests for the sealed exact-main baseline verifier."""

from __future__ import annotations

import ast
import base64
import copy
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import pathlib
import subprocess
import sys
import tarfile
import tempfile
import time
import unittest
from unittest import mock


TOOLS = pathlib.Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
import leopard2_exact_main_baseline as identity_contract  # noqa: E402
import leopard2_exact_main_baseline_record as record_contract  # noqa: E402
import leopard2_exact_main_baseline_verifier as verifier  # noqa: E402


def load_module(name: str, path: pathlib.Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


identity_fixtures = load_module(
    "leopard2_exact_main_baseline_fixture_helpers",
    TOOLS / "leopard2_exact_main_baseline_test.py",
)
record_fixtures = load_module(
    "leopard2_exact_main_baseline_record_fixture_helpers",
    TOOLS / "leopard2_exact_main_baseline_record_test.py",
)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def byte_identity(path: str, content: bytes) -> dict:
    return {"relative_path": path, "size": len(content), "sha256": sha256(content)}


def artifact_identity(
    name: str,
    build_path: str,
    retained_path: str,
    content: bytes,
) -> dict:
    return {
        "name": name,
        "build_relative_path": build_path,
        "retained_relative_path": retained_path,
        "size": len(content),
        "sha256": sha256(content),
    }


def roots(role: str) -> dict[str, str]:
    family = "canonical" if role != "path_variant" else "path-variant"
    base = f"/tmp/leopard-exact-main-{family}"
    return {
        "adapter_source_root": base + "/adapter",
        "baseline_source_root": base + "/baseline",
        "build_root": base + "/build",
    }


def _tar_info(name: str, *, directory: bool, size: int = 0) -> tarfile.TarInfo:
    info = tarfile.TarInfo(name)
    info.type = tarfile.DIRTYPE if directory else tarfile.REGTYPE
    info.mode = 0o775 if directory else 0o664
    info.uid = 0
    info.gid = 0
    info.uname = "root"
    info.gname = "root"
    info.mtime = 0
    info.size = 0 if directory else size
    return info


def source_tar(
    prefix: str,
    commit: str,
    members: dict[str, bytes],
    *,
    extra_directories: tuple[str, ...] = (),
) -> bytes:
    stream = io.BytesIO()
    with tarfile.open(
            fileobj=stream, mode="w", format=tarfile.PAX_FORMAT,
            pax_headers={"comment": commit}) as archive:
        archive.addfile(_tar_info(prefix[:-1], directory=True))
        directories: set[str] = set(extra_directories)
        for relative in sorted(members):
            components = relative.split("/")[:-1]
            for index in range(1, len(components) + 1):
                directories.add("/".join(components[:index]))
        for directory in sorted(directories):
            archive.addfile(_tar_info(prefix + directory, directory=True))
        for relative, content in sorted(members.items()):
            archive.addfile(
                _tar_info(prefix + relative, directory=False, size=len(content)),
                io.BytesIO(content),
            )
    return stream.getvalue()


def oversized_pax_archive() -> bytes:
    original = source_tar("fixture/", "1" * 40, {"file": b"x"})
    header = bytearray(original[:512])
    size = verifier.MAX_TAR_PAX_BYTES + 1
    header[124:136] = f"{size:011o}\0".encode("ascii")
    header[148:156] = b" " * 8
    header[148:156] = f"{sum(header):06o}\0 ".encode("ascii")
    padded = (size + 511) // 512 * 512
    return bytes(header) + b"\0" * padded + b"\0" * 1024


def git_object_id(kind: str, content: bytes) -> str:
    return hashlib.sha1(
        f"{kind} {len(content)}\0".encode("ascii") + content,
        usedforsecurity=False,
    ).hexdigest()


def git_object_identity(kind: str, content: bytes) -> dict:
    return {
        "encoding": "base64",
        "size": len(content),
        "sha256": sha256(content),
        "object_id": git_object_id(kind, content),
        "base64": base64.b64encode(content).decode("ascii"),
    }


def git_capture_digest(value: object) -> str:
    content = identity_contract.canonical_json_bytes(value)
    if not content.endswith(b"\n"):
        raise AssertionError("canonical fixture JSON lost its delimiter")
    return sha256(content[:-1])


def git_capture(
    path: str,
    members: dict[str, bytes],
    *,
    submodule: tuple[str, str] | None,
) -> tuple[dict, str, str]:
    leaves: dict[str, tuple[str, str, str]] = {
        relative: ("100644", "blob", git_object_id("blob", content))
        for relative, content in members.items()
    }
    if submodule is not None:
        relative, commit = submodule
        leaves[relative] = ("160000", "commit", commit)

    root: dict[str, object] = {}
    for relative, leaf in leaves.items():
        node = root
        parts = relative.split("/")
        for component in parts[:-1]:
            node = node.setdefault(component, {})  # type: ignore[assignment]
        node[parts[-1]] = leaf

    objects: dict[str, bytes] = {}

    def encode_tree(node: dict[str, object]) -> str:
        entries: list[tuple[bytes, bytes]] = []
        for name, value in node.items():
            if isinstance(value, dict):
                object_id = encode_tree(value)
                mode = "40000"
                sort_key = name.encode("utf-8") + b"/"
            else:
                mode, _kind, object_id = value
                sort_key = name.encode("utf-8")
            encoded = (mode.encode("ascii") + b" " + name.encode("utf-8") +
                       b"\0" + bytes.fromhex(object_id))
            entries.append((sort_key, encoded))
        content = b"".join(encoded for _key, encoded in sorted(entries))
        object_id = git_object_id("tree", content)
        objects[object_id] = content
        return object_id

    tree = encode_tree(root)
    commit_content = (
        f"tree {tree}\n"
        "author Fixture <fixture@example.invalid> 0 +0000\n"
        "committer Fixture <fixture@example.invalid> 0 +0000\n"
        "\nexact-main verifier fixture\n"
    ).encode("ascii")
    head = git_object_id("commit", commit_content)
    tracked_files = []
    submodules = []
    for relative in sorted(leaves):
        mode, git_type, object_id = leaves[relative]
        if git_type == "commit":
            nested = {
                "schema": "leopard2-git-source-capture/v2",
                "path": path + "/" + relative,
                "head": object_id,
                "tracked_status": "clean",
            }
            nested_sha = git_capture_digest(nested)
            tracked_files.append({
                "path": relative,
                "git_mode": mode,
                "git_type": git_type,
                "object_id": object_id,
                "kind": "submodule",
                "identity_sha256": nested_sha,
            })
            submodules.append({
                "path": relative,
                "object_id": object_id,
                "identity_sha256": nested_sha,
                "identity": nested,
            })
        else:
            tracked_files.append({
                "path": relative,
                "git_mode": mode,
                "git_type": git_type,
                "object_id": object_id,
                "kind": "regular",
            })
    listing = b"".join(
        (f"{item['git_mode']} {item['git_type']} {item['object_id']}\t"
         f"{item['path']}\0").encode("utf-8")
        for item in tracked_files)
    stage = b"".join(
        (f"{item['git_mode']} {item['object_id']} 0\t{item['path']}\0").
        encode("utf-8") for item in tracked_files)
    flags = b"".join(
        f"H {item['path']}\0".encode("utf-8")
        for item in tracked_files)
    capture = {
        "schema": "leopard2-git-source-capture/v2",
        "path": path,
        "head": head,
        "tree": tree,
        "detached": True,
        "head_ref": None,
        "superproject_worktree": None,
        "tracked_tree_listing_sha256": sha256(listing),
        "tracked_status": "clean",
        "commit_object": git_object_identity("commit", commit_content),
        "tree_objects": [
            git_object_identity("tree", objects[object_id])
            for object_id in sorted(objects)
        ],
        "git_executable": {},
        "git_metadata": {},
        "worktree_guard_policy": "fixture",
        "config": {"size": 0, "sha256": sha256(b"")},
        "index": {
            "entry_count": len(tracked_files),
            "stage": {"size": len(stage), "sha256": sha256(stage)},
            "flags_v": {"size": len(flags), "sha256": sha256(flags)},
            "flags_f": {"size": len(flags), "sha256": sha256(flags)},
        },
        "tracked_files": tracked_files,
        "tracked_files_sha256": git_capture_digest(tracked_files),
        "submodules": submodules,
    }
    return capture, head, tree


_REAL_BASELINE_COMMIT = record_contract.BASELINE_COMMIT
_REAL_BASELINE_TREE = record_contract.BASELINE_TREE
_TEST_BASELINE_MEMBERS = {
    name: f"fixture exact-main source {name}\n".encode("ascii")
    for name in ("leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
                 "LeopardFF16.cpp")
}
_unused_capture, _TEST_BASELINE_COMMIT, _TEST_BASELINE_TREE = git_capture(
    "/tmp/leopard-exact-main-canonical/baseline",
    _TEST_BASELINE_MEMBERS,
    submodule=("sse2neon", record_contract.BASELINE_SSE2NEON_COMMIT),
)


def verifier_cli_command(root: pathlib.Path) -> list[str]:
    code = (
        "import sys;"
        f"sys.path.insert(0,{str(TOOLS)!r});"
        "import leopard2_exact_main_baseline_record as r;"
        f"r.BASELINE_COMMIT={_TEST_BASELINE_COMMIT!r};"
        f"r.BASELINE_TREE={_TEST_BASELINE_TREE!r};"
        "import leopard2_exact_main_baseline_verifier as v;"
        f"v.BASELINE_COMMIT={_TEST_BASELINE_COMMIT!r};"
        f"v.BASELINE_TREE={_TEST_BASELINE_TREE!r};"
        "raise SystemExit(v.main([sys.argv[1]]))"
    )
    return [sys.executable, "-I", "-S", "-B", "-c", code, str(root)]


def ldd_bytes() -> bytes:
    return (
        b"ld-linux-x86-64.so.2\tfile\t"
        b"/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2\n"
        b"libc.so.6\tfile\t/usr/lib/x86_64-linux-gnu/libc.so.6\n"
        b"linux-vdso.so.1\tvirtual\n"
    )


def dependencies() -> list[dict]:
    return [
        {
            "soname": "ld-linux-x86-64.so.2", "kind": "file",
            "path": "/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2",
            "size": 234000, "sha256": record_fixtures.digest("ld-linux"),
        },
        {
            "soname": "libc.so.6", "kind": "file",
            "path": "/usr/lib/x86_64-linux-gnu/libc.so.6",
            "size": 2000000, "sha256": record_fixtures.digest("libc"),
        },
        {
            "soname": "linux-vdso.so.1", "kind": "virtual",
            "path": None, "size": None, "sha256": None,
        },
    ]


def cache_bytes(role_roots: dict[str, str]) -> bytes:
    return (
        "# exact-main fixture\n"
        "CMAKE_BUILD_TYPE:STRING=Release\n"
        "CMAKE_CXX_FLAGS:STRING=\n"
        "CMAKE_CXX_FLAGS_RELEASE:STRING=-g -O0 -O3\n"
        f"LEOPARD_MAIN_SOURCE_DIR:PATH={role_roots['baseline_source_root']}\n"
        "LEO_MAIN_PURE_AVX2:BOOL=ON\n"
    ).encode("utf-8")


def compile_commands_bytes(
    role_roots: dict[str, str],
    *,
    form: str = "arguments",
    escape_definition_quotes: bool = True,
) -> bytes:
    definitions = [
        f'-DLEOPARD_MAIN_SOURCE_COMMIT="{record_contract.BASELINE_COMMIT}"',
        "-DLEO_MAIN_PURE_AVX2_PROFILE=1",
    ]
    sources = [
        role_roots["adapter_source_root"] +
            "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp",
        *[role_roots["baseline_source_root"] + "/" + name for name in (
            "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
            "LeopardFF16.cpp")],
    ]
    value = []
    for index, source in enumerate(sources):
        output = (
            "CMakeFiles/leopard_main_benchmark.dir/"
            "legacy_main_benchmark.cpp.o" if index == 0 else
            "CMakeFiles/leopard_main_exact.dir" + source + ".o")
        arguments = [record_fixtures.TOOL_PATHS["compiler"]]
        if index == 0:
            arguments.extend(definitions)
        arguments.extend((
            "-I" + role_roots["baseline_source_root"],
            "-g", "-O0", "-O3", "-std=gnu++11",
            "-march=x86-64", "-mtune=generic", "-mavx2", "-mno-avx512f",
            "-Wall", "-Wextra", "-fopenmp",
            "-o", output, "-c", source,
        ))
        entry = {
            "directory": role_roots["build_root"],
            "file": source,
            "output": output,
        }
        if form == "arguments":
            entry["arguments"] = arguments
        elif form == "command":
            command_arguments = []
            for argument in arguments:
                if escape_definition_quotes and argument.startswith("-D"):
                    argument = argument.replace('"', r'\"')
                command_arguments.append(argument)
            entry["command"] = " ".join(command_arguments)
        else:
            raise AssertionError("unsupported compile-command fixture form")
        value.append(entry)
    return identity_contract.canonical_json_bytes(value)


def benchmark_stdout() -> bytes:
    value = {
        "schema": "leopard-main-benchmark-v1",
        "build": {
            "main_source_commit": record_contract.BASELINE_COMMIT,
            "pure_avx2": True,
            "cplusplus": 201103,
        },
        "parameters": {
            "K": 8, "R": 4, "shard_bytes": 64,
            "logical_shard_bytes": 64, "loss_count": 1,
            "missing_original_indices": [0], "batch": 2, "reuse": 1,
            "iterations": 2, "warmup": 1, "thread_count": 1, "seed": 7,
        },
        "resolved": {},
        "correctness": {"round_trip": True,
                        "logical_prefix_fingerprinted": True},
        "workload_digests": {},
        "memory": {},
        "timing": {},
    }
    return identity_contract.canonical_json_bytes(value)


def _version_bytes(role: str) -> tuple[bytes, bytes]:
    return f"{role} fixture version\n".encode("ascii"), b""


def _stage_bytes(name: str) -> bytes:
    return f"stage {name} complete\n".encode("ascii")


def build_closure_bytes(
    role: str,
    role_roots: dict[str, str],
    executable: bytes,
    archive: bytes,
) -> bytes:
    files = [
        {"relative_path": "leopard_main_benchmark",
         "size": len(executable), "sha256": sha256(executable)},
        {"relative_path": "libleopard_main_exact.a",
         "size": len(archive), "sha256": sha256(archive)},
    ]
    files.sort(key=lambda item: item["relative_path"])
    return identity_contract.canonical_json_bytes({
        "schema": verifier.BUILD_CLOSURE_SCHEMA,
        "role": role,
        "build_root": role_roots["build_root"],
        "files": files,
        "file_count": len(files),
    })


def _write_files(root: pathlib.Path, files: dict[str, bytes]) -> None:
    for relative, content in sorted(files.items()):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)


def _chmod_tree(root: pathlib.Path, *, sealed: bool) -> None:
    directory_mode = 0o500 if sealed else 0o700
    file_mode = 0o400 if sealed else 0o600
    for directory, subdirectories, filenames in os.walk(root, topdown=True):
        os.chmod(directory, directory_mode)
        for name in subdirectories:
            os.chmod(pathlib.Path(directory) / name, directory_mode)
        for name in filenames:
            os.chmod(pathlib.Path(directory) / name, file_mode)
    os.chmod(root, directory_mode)


def _metadata_document(root: pathlib.Path) -> dict:
    entries = []
    for directory, subdirectories, filenames in os.walk(root, topdown=True):
        subdirectories.sort()
        filenames.sort()
        directory_path = pathlib.Path(directory)
        relative_directory = directory_path.relative_to(root).as_posix()
        relative_directory = "." if relative_directory == "." else relative_directory
        status = os.lstat(directory_path)
        entries.append({
            "gid": status.st_gid, "mode": f"{stat_mode(status):04o}",
            "nlink": status.st_nlink, "path": relative_directory,
            "type": "directory", "uid": status.st_uid,
        })
        for filename in filenames:
            path = directory_path / filename
            relative = path.relative_to(root).as_posix()
            if relative == "TREE-METADATA.json":
                continue
            status = os.lstat(path)
            entries.append({
                "gid": status.st_gid, "mode": f"{stat_mode(status):04o}",
                "nlink": status.st_nlink, "path": relative,
                "type": "file", "uid": status.st_uid,
            })
    entries.sort(key=lambda item: item["path"])
    root_status = os.lstat(root)
    return {
        "entries": entries,
        "excluded_paths": ["TREE-METADATA.json"],
        "final_mode_policy": "observed mode with all write bits removed",
        "root": ".",
        "schema": verifier.TREE_METADATA_SCHEMA,
        "self_policy": {
            "gid": root_status.st_gid, "mode": "0400", "nlink": 1,
            "sha256_binding":
                "exactly one ./TREE-METADATA.json checksum entry",
            "type": "file", "uid": root_status.st_uid,
        },
        "uid_gid_policy": {
            "gid": root_status.st_gid,
            "rule": "every retained node has the invoking effective uid and gid",
            "uid": root_status.st_uid,
        },
    }


def stat_mode(status: os.stat_result) -> int:
    return status.st_mode & 0o7777


def seal_lane(root: pathlib.Path) -> None:
    metadata_path = root / "TREE-METADATA.json"
    sums_path = root / "SHA256SUMS"
    if not metadata_path.exists():
        metadata_path.write_bytes(b"")
    if not sums_path.exists():
        sums_path.write_bytes(b"")
    _chmod_tree(root, sealed=True)
    metadata = identity_contract.canonical_json_bytes(_metadata_document(root))
    os.chmod(metadata_path, 0o600)
    metadata_path.write_bytes(metadata)
    os.chmod(metadata_path, 0o400)
    paths = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS")
    sums = b"".join(
        f"{sha256(path.read_bytes())}  ./{path.relative_to(root).as_posix()}\n".
        encode("ascii") for path in paths)
    os.chmod(sums_path, 0o600)
    sums_path.write_bytes(sums)
    os.chmod(sums_path, 0o400)


def reseal_lane(root: pathlib.Path) -> None:
    _chmod_tree(root, sealed=False)
    metadata = root / "TREE-METADATA.json"
    sums = root / "SHA256SUMS"
    if metadata.exists():
        metadata.unlink()
    if sums.exists():
        sums.unlink()
    seal_lane(root)


def _resign_record(record: dict) -> dict:
    value = copy.deepcopy(record)
    value["record_sha256"] = identity_contract.canonical_sha256({
        key: copy.deepcopy(field) for key, field in value.items()
        if key != "record_sha256"
    })
    return value


def _replace_claimed_file(
    lane_root: pathlib.Path,
    claim: dict,
    content: bytes,
) -> None:
    (lane_root / claim["relative_path"]).write_bytes(content)
    claim["size"] = len(content)
    claim["sha256"] = sha256(content)


def authority_files_and_record() -> tuple[dict[str, bytes], dict]:
    files: dict[str, bytes] = {}
    canonical_elf = identity_fixtures.synthetic_elf(
        rodata=b"portable allocatable data", debug=b"canonical debug bytes")
    variant_elf = identity_fixtures.synthetic_elf(
        rodata=b"portable allocatable data", debug=b"variant debug bytes differ")
    canonical_archive = b"!<arch>\ncanonical exact-main archive\n"
    variant_archive = b"!<arch>\npath-variant exact-main archive\n"
    executables = {
        "canonical_first": canonical_elf,
        "canonical_second": canonical_elf,
        "path_variant": variant_elf,
    }
    archives = {
        "canonical_first": canonical_archive,
        "canonical_second": canonical_archive,
        "path_variant": variant_archive,
    }

    adapter_contents = {
        path: f"fixture adapter {path}\n".encode("ascii")
        for path in record_contract.ADAPTER_PATHS
    }
    canonical_roots = roots("canonical_first")
    baseline_capture, baseline_commit, baseline_tree = git_capture(
        canonical_roots["baseline_source_root"], _TEST_BASELINE_MEMBERS,
        submodule=("sse2neon", record_contract.BASELINE_SSE2NEON_COMMIT))
    adapter_capture, adapter_commit, adapter_tree = git_capture(
        canonical_roots["adapter_source_root"], adapter_contents,
        submodule=("sse2neon", record_contract.BASELINE_SSE2NEON_COMMIT))
    if (baseline_commit, baseline_tree) != (
            record_contract.BASELINE_COMMIT, record_contract.BASELINE_TREE):
        raise AssertionError("synthetic baseline trust anchors drifted")
    adapter_archive_content = source_tar(
        "leopard2-adapter-source/", adapter_commit, adapter_contents,
        extra_directories=("sse2neon",))
    baseline_archive_content = source_tar(
        "leopard-main-source/", record_contract.BASELINE_COMMIT,
        _TEST_BASELINE_MEMBERS, extra_directories=("sse2neon",))
    files["source/leopard2-adapter-source.tar"] = adapter_archive_content
    files["source/leopard-main-source.tar"] = baseline_archive_content

    files["source/leopard-main-git-capture.json"] = \
        identity_contract.canonical_json_bytes(baseline_capture)
    files["source/adapter-git-capture.json"] = \
        identity_contract.canonical_json_bytes(adapter_capture)

    for index, name in enumerate(record_contract.STAGE_NAMES):
        files[f"logs/{index:02d}-{name}.log"] = _stage_bytes(name)

    toolchain = record_fixtures.toolchain_fixture()
    for role in record_contract.VERSION_ROLES:
        stdout, stderr = _version_bytes(role)
        stdout_path = f"toolchain/versions/{role}.stdout"
        stderr_path = f"toolchain/versions/{role}.stderr"
        files[stdout_path] = stdout
        files[stderr_path] = stderr
        version = next(item for item in toolchain["versions"]
                       if item["role"] == role)
        version["stdout"] = byte_identity(stdout_path, stdout)
        version["stderr"] = byte_identity(stderr_path, stderr)

    builds: dict[str, dict] = {}
    identities: dict[str, dict] = {}
    runtime_records = []
    attestation_records = []
    ldd = ldd_bytes()
    benchmark = benchmark_stdout()
    for role in record_contract.BUILD_ROLES:
        role_path = role.replace("_", "-")
        role_roots = roots(role)
        configure_log = f"{role} configure complete\n".encode("ascii")
        build_log = f"{role} build complete\n".encode("ascii")
        cache = cache_bytes(role_roots)
        commands = compile_commands_bytes(role_roots)
        closure = build_closure_bytes(
            role, role_roots, executables[role], archives[role])
        build_paths = {
            "configure_log": f"builds/{role_path}/configure.log",
            "build_log": f"builds/{role_path}/build.log",
            "cmake_cache": f"builds/{role_path}/CMakeCache.txt",
            "compile_commands": f"builds/{role_path}/compile_commands.json",
            "closure": f"builds/{role_path}/build-closure.json",
        }
        files[build_paths["configure_log"]] = configure_log
        files[build_paths["build_log"]] = build_log
        files[build_paths["cmake_cache"]] = cache
        files[build_paths["compile_commands"]] = commands
        files[build_paths["closure"]] = closure
        executable_path = f"artifacts/{role_path}/leopard_main_benchmark"
        archive_path = f"artifacts/{role_path}/libleopard_main_exact.a"
        files[executable_path] = executables[role]
        files[archive_path] = archives[role]
        builds[role] = {
            "role": role,
            "roots": role_roots,
            "configure_argv": [
                record_fixtures.TOOL_PATHS["cmake"],
                "-S", role_roots["adapter_source_root"] +
                    "/experiments/leopard2/main_compare",
                "-B", role_roots["build_root"], "-G", "Unix Makefiles",
                "-DCMAKE_BUILD_TYPE=Release", "-DLEO_MAIN_PURE_AVX2=ON",
                "-DLEOPARD_MAIN_SOURCE_DIR=" +
                    role_roots["baseline_source_root"],
                "-DCMAKE_CXX_COMPILER=" +
                    record_fixtures.TOOL_PATHS["compiler"],
            ],
            "build_argv": [
                record_fixtures.TOOL_PATHS["cmake"], "--build",
                role_roots["build_root"], "--target",
                "leopard_main_benchmark", "--", "-j1",
            ],
            "configure_log": byte_identity(
                build_paths["configure_log"], configure_log),
            "build_log": byte_identity(build_paths["build_log"], build_log),
            "cmake_cache": byte_identity(build_paths["cmake_cache"], cache),
            "compile_commands": byte_identity(
                build_paths["compile_commands"], commands),
            "executable": artifact_identity(
                "leopard_main_benchmark", "leopard_main_benchmark",
                executable_path, executables[role]),
            "archive": artifact_identity(
                "libleopard_main_exact.a", "libleopard_main_exact.a",
                archive_path, archives[role]),
            "closure": {
                **byte_identity(build_paths["closure"], closure),
                "file_count": 2,
            },
        }
        identities[role] = \
            identity_contract.normalized_code_identity_from_elf_bytes(
                executables[role], roots=role_roots)

        ldd_path = f"runtime/{role_path}/ldd.txt"
        files[ldd_path] = ldd
        runtime_records.append({
            "role": role,
            "executable_sha256": sha256(executables[role]),
            "canonical_ldd_output": byte_identity(ldd_path, ldd),
            "dependencies": dependencies(),
        })

        stdout_path = f"attestations/{role_path}/benchmark.stdout.json"
        stderr_path = f"attestations/{role_path}/benchmark.stderr"
        ctest_stdout_path = f"attestations/{role_path}/ctest.stdout.log"
        ctest_stderr_path = f"attestations/{role_path}/ctest.stderr.log"
        ctest_stdout = b"100% tests passed, 0 tests failed out of 1\n"
        files[stdout_path] = benchmark
        files[stderr_path] = b""
        files[ctest_stdout_path] = ctest_stdout
        files[ctest_stderr_path] = b""
        executable_build_path = role_roots["build_root"] + \
            "/leopard_main_benchmark"
        attestation_records.append({
            "role": role,
            "argv": [
                executable_build_path,
                "--k", "8", "--r", "4", "--bytes", "64", "--loss", "1",
                "--batch", "2", "--reuse", "1", "--iterations", "2",
                "--warmup", "1", "--threads", "1", "--seed", "7",
                "--json", "-",
            ],
            "stdout": byte_identity(stdout_path, benchmark),
            "stderr": byte_identity(stderr_path, b""),
            "exit_status": 0,
            "reported_schema": "leopard-main-benchmark-v1",
            "main_source_commit": record_contract.BASELINE_COMMIT,
            "pure_avx2": True,
            "round_trip": True,
            "ctest": {
                "argv": [
                    record_fixtures.TOOL_PATHS["ctest"], "--test-dir",
                    role_roots["build_root"], "--output-on-failure", "-R",
                    "^leopard_main_benchmark_smoke$",
                ],
                "stdout": byte_identity(ctest_stdout_path, ctest_stdout),
                "stderr": byte_identity(ctest_stderr_path, b""),
                "exit_status": 0, "passed": 1, "failed": 0,
            },
        })

    controller_path = "controllers/test_legacy_main_benchmark.py"
    controller_bytes = adapter_contents[record_contract.ADAPTER_PATHS[2]]
    files[controller_path] = controller_bytes
    adapter_files = [{
        "path": path,
        "git_blob_sha1": hashlib.sha1(
            f"blob {len(content)}\0".encode("ascii") + content,
            usedforsecurity=False).hexdigest(),
        "size": len(content),
        "sha256": sha256(content),
    } for path, content in sorted(adapter_contents.items())]

    lane = {
        "root": "/tmp/leopard-exact-main-baseline-a1",
        "attempt": 1, "attempt_budget": 3,
        "record_relative_path": "baseline-authority.json",
        "seal_protocol": record_contract.SEAL_PROTOCOL,
        "stages": [{
            "name": name, "status": "complete",
            "log": byte_identity(f"logs/{index:02d}-{name}.log", files[
                f"logs/{index:02d}-{name}.log"]),
        } for index, name in enumerate(record_contract.STAGE_NAMES)],
    }
    source = {
        "baseline": {
            "commit": record_contract.BASELINE_COMMIT,
            "tree": record_contract.BASELINE_TREE,
            "submodule": {"path": "sse2neon",
                          "commit": record_contract.BASELINE_SSE2NEON_COMMIT},
            "git_capture": byte_identity(
                "source/leopard-main-git-capture.json",
                files["source/leopard-main-git-capture.json"]),
            "archive": {
                **byte_identity("source/leopard-main-source.tar",
                                baseline_archive_content),
                "prefix": "leopard-main-source/",
                "replay_sha256": sha256(baseline_archive_content),
                "replay_identical": True,
            },
        },
        "adapter_repository": {
            "commit": adapter_commit, "tree": adapter_tree, "clean": True,
            "git_capture": byte_identity(
                "source/adapter-git-capture.json",
                files["source/adapter-git-capture.json"]),
            "archive": {
                **byte_identity("source/leopard2-adapter-source.tar",
                                adapter_archive_content),
                "prefix": "leopard2-adapter-source/",
                "replay_sha256": sha256(adapter_archive_content),
                "replay_identical": True,
            },
        },
    }
    adapter = {
        "minimum_harness_commit": record_contract.MINIMUM_HARNESS_COMMIT,
        "files": adapter_files,
        "files_sha256": identity_contract.canonical_sha256(adapter_files),
    }
    combined = identities["canonical_first"]["combined_sha256"]
    identity = {
        **identities,
        "combined_sha256": combined,
        "canonical_raw_bytes_identical": True,
        "path_variant_raw_bytes_differ": True,
        "normalized_match": True,
    }
    attestation = {
        "schema": record_contract.ATTESTATION_SCHEMA,
        "test_controller": byte_identity(controller_path, controller_bytes),
        "records": attestation_records,
    }
    record = record_contract.baseline_authority_record(
        created_utc="2026-08-29T23:00:00Z",
        host=record_fixtures.host_fixture(), lane=lane, source=source,
        adapter=adapter, toolchain=toolchain, builds=builds,
        runtime_closure={
            "schema": record_contract.RUNTIME_CLOSURE_SCHEMA,
            "normalization": record_contract.CANONICAL_LDD_NORMALIZATION,
            "records": runtime_records,
        },
        attestation=attestation, identity=identity,
    )
    files["baseline-authority.json"] = \
        identity_contract.canonical_json_bytes(record)
    return files, record


def failure_files_and_record(
    *, verification_failure: bool = False, retained: bool = True,
) -> tuple[dict[str, bytes], dict]:
    final_index = 4 if verification_failure else 2
    stage_names = list(record_contract.STAGE_NAMES[:final_index + 1])
    files = {
        f"logs/{index:02d}-{name}.log":
            f"{name} {'failed' if index == final_index else 'complete'}\n".
            encode("ascii")
        for index, name in enumerate(stage_names)
    }
    retained_records = []
    if retained:
        files["diagnostics/compiler.stderr"] = b"compiler diagnostic\n"
        retained_records.append(byte_identity(
            "diagnostics/compiler.stderr", files["diagnostics/compiler.stderr"]))
    lane = {
        "root": "/tmp/leopard-exact-main-baseline-a1",
        "attempt": 1, "attempt_budget": 3,
        "record_relative_path": "FAILED.json",
        "seal_protocol": record_contract.SEAL_PROTOCOL,
        "stages": [{
            "name": name,
            "status": "failed" if index == final_index else "complete",
            "log": byte_identity(
                f"logs/{index:02d}-{name}.log",
                files[f"logs/{index:02d}-{name}.log"]),
        } for index, name in enumerate(stage_names)],
    }
    if verification_failure:
        record = record_contract.baseline_verification_failure_record(
            created_utc="2026-08-29T23:01:00Z", lane=lane,
            stage=stage_names[-1],
            error={"kind": "verification_error", "message": "fixture failed",
                   "exit_status": 1},
            retained_files=retained_records,
            authority_record_sha256="a" * 64,
        )
    else:
        record = record_contract.baseline_acquisition_failure_record(
            created_utc="2026-08-29T23:01:00Z", lane=lane,
            stage=stage_names[-1],
            error={"kind": "build_error", "message": "fixture failed",
                   "exit_status": 1},
            retained_files=retained_records,
        )
    files["FAILED.json"] = identity_contract.canonical_json_bytes(record)
    return files, record


class Lane:
    def __init__(self, *, kind: str = "authority", retained: bool = True):
        self._anchors = contextlib.ExitStack()
        for module in (record_contract, verifier):
            self._anchors.enter_context(mock.patch.object(
                module, "BASELINE_COMMIT", _TEST_BASELINE_COMMIT))
            self._anchors.enter_context(mock.patch.object(
                module, "BASELINE_TREE", _TEST_BASELINE_TREE))
        try:
            self.temporary = tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-verifier-test.")
            self.root = pathlib.Path(self.temporary.name) / "lane"
            self.root.mkdir(mode=0o700)
            if kind == "authority":
                self.files, self.record = authority_files_and_record()
            else:
                self.files, self.record = failure_files_and_record(
                    verification_failure=kind == "verification_failure",
                    retained=retained)
            _write_files(self.root, self.files)
            seal_lane(self.root)
        except BaseException:
            self._anchors.close()
            raise

    def __enter__(self) -> "Lane":
        return self

    def __exit__(self, _kind: object, _value: object, _traceback: object) -> None:
        try:
            _chmod_tree(self.root, sealed=False)
            self.temporary.cleanup()
        finally:
            self._anchors.close()

    def writable(self) -> None:
        _chmod_tree(self.root, sealed=False)

    def reseal(self) -> None:
        reseal_lane(self.root)


class ExactMainBaselineVerifierTest(unittest.TestCase):
    def assertInvalid(self, root: pathlib.Path) -> None:  # noqa: N802
        with self.assertRaises(verifier.SealedLaneError):
            verifier.verify_sealed_lane(root)

    def test_authority_and_failure_verdicts(self) -> None:
        self.assertEqual(verifier.BASELINE_COMMIT, _REAL_BASELINE_COMMIT)
        self.assertEqual(verifier.BASELINE_TREE, _REAL_BASELINE_TREE)
        with Lane() as lane:
            capture = json.loads((
                lane.root / "source/adapter-git-capture.json").read_bytes())
            self.assertEqual(
                capture["tracked_files_sha256"],
                git_capture_digest(capture["tracked_files"]))
            self.assertNotEqual(
                capture["tracked_files_sha256"],
                identity_contract.canonical_sha256(capture["tracked_files"]))
            verdict = verifier.verify_sealed_lane(lane.root)
            self.assertEqual(verdict["outcome"], "promoted_authority")
            self.assertTrue(verdict["promoted"])
            self.assertEqual(verdict["seal"]["file_count"], 76)
            self.assertEqual(verdict["seal"]["checksum_line_count"], 75)
            self.assertEqual(verdict["seal"]["directory_count"], 22)
            self.assertEqual(verdict["producer_attested"], sorted((
                "/attestation/records/*/{exit_status,stderr/content}",
                "/attestation/records/*/stdout/uninterpreted-fields",
                "/attestation/records/*/ctest/"
                    "{exit_status,failed,passed,stdout/content,stderr/content}",
                "/builds/*/{configure_argv,build_argv,configure_log/content,"
                    "build_log/content}",
                "/builds/*/closure/files/<non-artifact>/{size,sha256}",
                "/created_utc",
                "/host",
                "/lane/{attempt,attempt_budget,root}",
                "/lane/stages/*/log/content",
                "/runtime_closure/records/*/dependencies/*/{sha256,size}",
                "/source/*/archive/{replay_identical,replay_sha256}",
                "/source/*/git_capture/{config,git_executable,git_metadata,"
                    "superproject_worktree,worktree_guard_policy}",
                "/source/*/git_capture/{detached,head_ref,tracked_status}",
                "/source/*/git_capture/submodules/*/identity/interior",
                "/toolchain/subtools",
                "/toolchain/tools",
                "/toolchain/versions",
            )))
            self.assertEqual(verdict["verdict_sha256"],
                             identity_contract.canonical_sha256({
                                 key: copy.deepcopy(value)
                                 for key, value in verdict.items()
                                 if key != "verdict_sha256"
                             }))
        for kind, retained in (
                ("acquisition_failure", True),
                ("acquisition_failure", False),
                ("verification_failure", True)):
            with self.subTest(kind=kind, retained=retained), \
                    Lane(kind=kind, retained=retained) as lane:
                verdict = verifier.verify_sealed_lane(lane.root)
                self.assertEqual(verdict["outcome"], "verified_failure")
                self.assertFalse(verdict["promoted"])

    def test_shape_link_and_mode_adversaries(self) -> None:
        mutations = ("missing", "extra_file", "extra_directory", "hardlink",
                     "symlink", "fifo", "file_mode", "directory_mode",
                     "root_mode")
        for mutation in mutations:
            with self.subTest(mutation=mutation), Lane() as lane:
                lane.writable()
                victim = lane.root / "logs/00-source_acquisition.log"
                if mutation == "missing":
                    victim.unlink()
                elif mutation == "extra_file":
                    (lane.root / "extra.bin").write_bytes(b"extra")
                elif mutation == "extra_directory":
                    (lane.root / "empty-extra").mkdir()
                elif mutation == "hardlink":
                    os.link(victim, lane.root / "hardlink")
                elif mutation == "symlink":
                    victim.unlink()
                    os.symlink("../source/leopard-main-source.tar", victim)
                elif mutation == "fifo":
                    victim.unlink()
                    os.mkfifo(victim)
                elif mutation == "file_mode":
                    os.chmod(victim, 0o444)
                elif mutation == "directory_mode":
                    os.chmod(victim.parent, 0o550)
                else:
                    os.chmod(lane.root, 0o550)
                if mutation not in ("file_mode", "directory_mode", "root_mode"):
                    _chmod_tree(lane.root, sealed=True)
                self.assertInvalid(lane.root)

    def test_root_identity_and_resource_bounds(self) -> None:
        with Lane() as lane, mock.patch.object(
                os, "geteuid", return_value=os.geteuid() + 1):
            self.assertInvalid(lane.root)
        with Lane() as lane, mock.patch.object(
                os, "getegid", return_value=os.getegid() + 1):
            self.assertInvalid(lane.root)
        for name, value in (
                ("MAX_TREE_DEPTH", 1), ("MAX_TREE_NODES", 2),
                ("MAX_SEALED_FILE_BYTES", 1),
                ("MAX_SEALED_TOTAL_BYTES", 1),
                ("MAX_PARSED_JSON_BYTES", 1), ("MAX_TAR_MEMBERS", 1),
                ("MAX_GIT_COMMIT_BYTES", 1), ("MAX_GIT_OBJECT_BYTES", 1),
                ("MAX_GIT_TREE_OBJECTS", 0),
                ("MAX_GIT_TREE_TOTAL_BYTES", 1),
                ("MAX_GIT_TRACKED_FILES", 1)):
            with self.subTest(bound=name), Lane() as lane, \
                    mock.patch.object(verifier, name, value):
                self.assertInvalid(lane.root)

        with Lane() as lane:
            real_stat = os.stat

            class ForeignOwner:
                def __init__(self, status):
                    self._status = status
                    self.st_uid = status.st_uid + 1

                def __getattr__(self, name):
                    return getattr(self._status, name)

            def foreign_node(path, *args, **kwargs):
                status = real_stat(path, *args, **kwargs)
                return ForeignOwner(status) if path == "logs" else status

            with mock.patch.object(os, "stat", side_effect=foreign_node), \
                    self.assertRaises(verifier.SealedLaneError) as caught:
                verifier.verify_sealed_lane(lane.root)
            self.assertIn("has another owner", str(caught.exception))

        with Lane() as lane:
            real_scandir = os.scandir
            observed = {"count": 0}

            class CountingScandir:
                def __init__(self, path):
                    self._context = real_scandir(path)

                def __enter__(self):
                    iterator = self._context.__enter__()

                    class CountingIterator:
                        def __iter__(self):
                            return self

                        def __next__(self):
                            value = next(iterator)
                            observed["count"] += 1
                            return value

                    return CountingIterator()

                def __exit__(self, *arguments):
                    return self._context.__exit__(*arguments)

            with mock.patch.object(verifier, "MAX_TREE_NODES", 2), \
                    mock.patch.object(
                        os, "scandir", side_effect=CountingScandir):
                self.assertInvalid(lane.root)
            self.assertLessEqual(observed["count"], 3)

    def test_checksum_and_metadata_adversaries(self) -> None:
        mutations = ("file_byte", "ledger_order", "ledger_separator",
                     "ledger_crlf", "metadata_pretty", "metadata_schema",
                     "metadata_entry")
        for mutation in mutations:
            with self.subTest(mutation=mutation), Lane() as lane:
                lane.writable()
                if mutation == "file_byte":
                    path = lane.root / "logs/00-source_acquisition.log"
                    path.write_bytes(path.read_bytes() + b"changed")
                elif mutation.startswith("ledger_"):
                    path = lane.root / "SHA256SUMS"
                    content = path.read_bytes()
                    if mutation == "ledger_order":
                        lines = content.splitlines(keepends=True)
                        content = b"".join(reversed(lines))
                    elif mutation == "ledger_separator":
                        content = content.replace(b"  ./", b" ./", 1)
                    else:
                        content = content.replace(b"\n", b"\r\n")
                    path.write_bytes(content)
                else:
                    path = lane.root / "TREE-METADATA.json"
                    value = json.loads(path.read_bytes())
                    if mutation == "metadata_pretty":
                        path.write_text(json.dumps(value, indent=2) + "\n")
                    elif mutation == "metadata_schema":
                        value["schema"] += "-changed"
                        path.write_bytes(identity_contract.canonical_json_bytes(value))
                    else:
                        value["entries"] = value["entries"][:-1]
                        path.write_bytes(identity_contract.canonical_json_bytes(value))
                _chmod_tree(lane.root, sealed=True)
                self.assertInvalid(lane.root)

    def test_terminal_and_seal_presence_are_exact(self) -> None:
        mutations = ("missing_sums", "missing_metadata", "both_terminals",
                     "neither_terminal", "unclaimed_failure_file")
        for mutation in mutations:
            kind = "acquisition_failure" if \
                mutation == "unclaimed_failure_file" else "authority"
            with self.subTest(mutation=mutation), Lane(kind=kind) as lane:
                lane.writable()
                if mutation == "missing_sums":
                    (lane.root / "SHA256SUMS").unlink()
                elif mutation == "missing_metadata":
                    (lane.root / "TREE-METADATA.json").unlink()
                elif mutation == "both_terminals":
                    _unused_files, failure = failure_files_and_record()
                    (lane.root / "FAILED.json").write_bytes(
                        identity_contract.canonical_json_bytes(failure))
                elif mutation == "neither_terminal":
                    (lane.root / "baseline-authority.json").unlink()
                else:
                    (lane.root / "diagnostics/unclaimed.txt").parent.mkdir(
                        parents=True, exist_ok=True)
                    (lane.root / "diagnostics/unclaimed.txt").write_bytes(
                        b"unclaimed\n")
                if mutation in ("both_terminals", "unclaimed_failure_file"):
                    lane.reseal()
                else:
                    _chmod_tree(lane.root, sealed=True)
                self.assertInvalid(lane.root)

    def test_canonical_record_and_semantic_adversaries(self) -> None:
        with Lane() as lane:
            lane.writable()
            path = lane.root / "baseline-authority.json"
            path.write_text(json.dumps(lane.record, indent=2) + "\n")
            lane.reseal()
            self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            role = "path_variant"
            executable_path = record["builds"][role]["executable"][
                "retained_relative_path"]
            artifact = bytearray((lane.root / executable_path).read_bytes())
            text_offset = artifact.find(b"\x90\xc3")
            self.assertGreaterEqual(text_offset, 0)
            artifact[text_offset] ^= 1
            changed_bytes = bytes(artifact)
            (lane.root / executable_path).write_bytes(changed_bytes)
            changed_sha = sha256(changed_bytes)
            record["builds"][role]["executable"]["sha256"] = changed_sha
            record["builds"][role]["executable"]["size"] = len(changed_bytes)
            record["runtime_closure"]["records"][2][
                "executable_sha256"] = changed_sha
            old_identity = record["identity"][role]
            record["identity"][role] = \
                identity_contract.normalized_code_identity_record(
                    artifact={"size": len(changed_bytes), "sha256": changed_sha},
                    sections=copy.deepcopy(old_identity["sections"]),
                    path_string_census=copy.deepcopy(
                        old_identity["path_string_census"]),
                )
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            changed = ldd_bytes().replace(b"/libc.so.6", b"/libz.so.6")
            for runtime_record in record["runtime_closure"]["records"]:
                _replace_claimed_file(
                    lane.root, runtime_record["canonical_ldd_output"], changed)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            changed = b"not-a-canonical-runtime-closure"
            for runtime_record in record["runtime_closure"]["records"]:
                _replace_claimed_file(
                    lane.root, runtime_record["canonical_ldd_output"], changed)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

    def test_archive_and_attestation_semantic_adversaries(self) -> None:
        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            changed = oversized_pax_archive()
            claim = record["source"]["adapter_repository"]["archive"]
            _replace_claimed_file(lane.root, claim, changed)
            claim["replay_sha256"] = sha256(changed)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            with mock.patch.object(
                    verifier, "_tar_member_bytes",
                    side_effect=AssertionError("oversized pax was read")):
                self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            path = "source/leopard2-adapter-source.tar"
            changed = bytearray((lane.root / path).read_bytes())
            changed[148] ^= 1
            changed_bytes = bytes(changed)
            (lane.root / path).write_bytes(changed_bytes)
            archive = record["source"]["adapter_repository"]["archive"]
            archive["size"] = len(changed_bytes)
            archive["sha256"] = sha256(changed_bytes)
            archive["replay_sha256"] = sha256(changed_bytes)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            stream = io.BytesIO()
            with tarfile.open(
                    fileobj=stream, mode="w", format=tarfile.PAX_FORMAT,
                    pax_headers={"comment": record["source"]
                                 ["adapter_repository"]["commit"]}) as archive:
                archive.addfile(_tar_info(
                    "leopard2-adapter-source", directory=True))
                link = tarfile.TarInfo(
                    "leopard2-adapter-source/forbidden-link")
                link.type = tarfile.SYMTYPE
                link.linkname = "/etc/passwd"
                link.mode = 0o777
                link.mtime = 0
                archive.addfile(link)
                for relative, expected in zip(
                        record_contract.ADAPTER_PATHS,
                        record["adapter"]["files"]):
                    content = (lane.root /
                               "controllers/test_legacy_main_benchmark.py").\
                        read_bytes() if relative == record_contract.ADAPTER_PATHS[2] \
                        else f"fixture adapter {relative}\n".encode("ascii")
                    self.assertEqual(len(content), expected["size"])
                    archive.addfile(_tar_info(
                        "leopard2-adapter-source/" + relative,
                        directory=False, size=len(content)), io.BytesIO(content))
            changed = stream.getvalue()
            path = "source/leopard2-adapter-source.tar"
            (lane.root / path).write_bytes(changed)
            archive_claim = record["source"]["adapter_repository"]["archive"]
            archive_claim["size"] = len(changed)
            archive_claim["sha256"] = sha256(changed)
            archive_claim["replay_sha256"] = sha256(changed)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

    def test_toolchain_flag_overrides_and_command_form(self) -> None:
        for mutation in (
                "isa_override", "definition_override",
                "definition_alias_override", "definition_passthrough_override",
                "cache_flag_injection", "cache_flag_removed"):
            with self.subTest(mutation=mutation), Lane() as lane:
                lane.writable()
                record = copy.deepcopy(lane.record)
                for role in record_contract.BUILD_ROLES:
                    build = record["builds"][role]
                    command_claim = build["compile_commands"]
                    commands = json.loads(
                        (lane.root / command_claim["relative_path"]).read_bytes())
                    if mutation == "isa_override":
                        for command in commands:
                            command["arguments"].append("-mavx512f")
                    elif mutation == "definition_override":
                        commands[0]["arguments"].extend((
                            '-DLEOPARD_MAIN_SOURCE_COMMIT="' + "0" * 40 + '"',
                            "-DLEO_MAIN_PURE_AVX2_PROFILE=0",
                        ))
                    elif mutation == "definition_alias_override":
                        commands[0]["arguments"].append(
                            "--define-macro=LEO_MAIN_PURE_AVX2_PROFILE=0")
                    elif mutation == "definition_passthrough_override":
                        commands[0]["arguments"].append(
                            "-Wp,-DLEO_MAIN_PURE_AVX2_PROFILE=0")
                    _replace_claimed_file(
                        lane.root, command_claim,
                        identity_contract.canonical_json_bytes(commands))
                    if mutation in ("cache_flag_injection",
                                    "cache_flag_removed"):
                        cache_claim = build["cmake_cache"]
                        cache = (lane.root / cache_claim["relative_path"]).read_bytes()
                        replacement = b"" if mutation == "cache_flag_removed" \
                            else b"CMAKE_CXX_FLAGS:STRING=-mavx512f\n"
                        cache = cache.replace(
                            b"CMAKE_CXX_FLAGS:STRING=\n", replacement)
                        _replace_claimed_file(lane.root, cache_claim, cache)
                record = _resign_record(record)
                (lane.root / "baseline-authority.json").write_bytes(
                    identity_contract.canonical_json_bytes(record))
                lane.reseal()
                self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            for role in record_contract.BUILD_ROLES:
                build = record["builds"][role]
                _replace_claimed_file(
                    lane.root, build["compile_commands"],
                    compile_commands_bytes(build["roots"], form="command"))
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertTrue(verifier.verify_sealed_lane(lane.root)["promoted"])

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            for role in record_contract.BUILD_ROLES:
                build = record["builds"][role]
                _replace_claimed_file(
                    lane.root, build["compile_commands"],
                    compile_commands_bytes(
                        build["roots"], form="command",
                        escape_definition_quotes=False))
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        injections = (
            ("@/tmp/isa.rsp",),
            ("-imacros", "/tmp/override.h"),
            ("-include", "/tmp/override.h"),
            ("-ffast-math",),
            ("-Xclang", "-target-feature", "+avx512f"),
            ("-ULEO_MAIN_PURE_AVX2_PROFILE",),
            ("-U", "LEO_MAIN_PURE_AVX2_PROFILE"),
            ("--undefine-macro=LEO_MAIN_PURE_AVX2_PROFILE",),
            ("--undefine-macro", "LEO_MAIN_PURE_AVX2_PROFILE"),
            ("-D", "LEO_MAIN_PURE_AVX2_PROFILE=0"),
            ("-Xpreprocessor", "-DLEO_MAIN_PURE_AVX2_PROFILE=0"),
            ("-Xpreprocessor", "-ULEO_MAIN_PURE_AVX2_PROFILE"),
        )
        for form in ("arguments", "command"):
            for injection in injections:
                with self.subTest(form=form, injection=injection), Lane() as lane:
                    lane.writable()
                    record = copy.deepcopy(lane.record)
                    for role in record_contract.BUILD_ROLES:
                        build = record["builds"][role]
                        commands = json.loads(compile_commands_bytes(
                            build["roots"], form=form))
                        if form == "arguments":
                            commands[0]["arguments"].extend(injection)
                        else:
                            commands[0]["command"] += " " + " ".join(injection)
                        _replace_claimed_file(
                            lane.root, build["compile_commands"],
                            identity_contract.canonical_json_bytes(commands))
                    record = _resign_record(record)
                    (lane.root / "baseline-authority.json").write_bytes(
                        identity_contract.canonical_json_bytes(record))
                    lane.reseal()
                    self.assertInvalid(lane.root)

    def test_git_capture_and_source_archive_closure_adversaries(self) -> None:
        def empty_tree_dag(depth: int):
            empty_tree = b""
            child = git_object_id("tree", empty_tree)
            tree_objects = {child: empty_tree}
            for _level in range(depth):
                content = b"".join((
                    b"40000 a\0" + bytes.fromhex(child),
                    b"40000 b\0" + bytes.fromhex(child),
                ))
                child = git_object_id("tree", content)
                tree_objects[child] = content
            return child, tree_objects

        positive_root, positive_objects = empty_tree_dag(15)
        self.assertEqual(verifier._flatten_git_tree_objects(
            positive_root, positive_objects), [])
        child, tree_objects = empty_tree_dag(22)
        started = time.monotonic()
        with self.assertRaises(verifier.SealedLaneError) as caught:
            verifier._flatten_git_tree_objects(child, tree_objects)
        self.assertIn("expansion exceeds its bound", str(caught.exception))
        self.assertLess(time.monotonic() - started, 5.0)

        for mutation in (
                "commit_object", "tree_object", "tracked_digest",
                "submodule_object", "submodule_digest"):
            with self.subTest(mutation=mutation), Lane() as lane:
                lane.writable()
                record = copy.deepcopy(lane.record)
                claim = record["source"]["adapter_repository"]["git_capture"]
                capture = json.loads(
                    (lane.root / claim["relative_path"]).read_bytes())
                if mutation == "commit_object":
                    capture["commit_object"]["sha256"] = "f" * 64
                elif mutation == "tree_object":
                    capture["tree_objects"][0]["object_id"] = "f" * 40
                elif mutation == "tracked_digest":
                    capture["tracked_files_sha256"] = "f" * 64
                elif mutation == "submodule_object":
                    capture["submodules"][0]["object_id"] = "f" * 40
                else:
                    capture["submodules"][0]["identity_sha256"] = "f" * 64
                _replace_claimed_file(
                    lane.root, claim,
                    identity_contract.canonical_json_bytes(capture))
                record = _resign_record(record)
                (lane.root / "baseline-authority.json").write_bytes(
                    identity_contract.canonical_json_bytes(record))
                lane.reseal()
                self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            members = dict(_TEST_BASELINE_MEMBERS)
            members["leopard.cpp"] += b"forged source\n"
            archive = source_tar(
                record["source"]["baseline"]["archive"]["prefix"],
                record["source"]["baseline"]["commit"], members)
            claim = record["source"]["baseline"]["archive"]
            _replace_claimed_file(lane.root, claim, archive)
            claim["replay_sha256"] = sha256(archive)
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            members = {
                path: f"fixture adapter {path}\n".encode("ascii")
                for path in record_contract.ADAPTER_PATHS
            }
            changed_path = record_contract.ADAPTER_PATHS[0]
            members[changed_path] += b"forged source\n"
            source = record["source"]["adapter_repository"]
            archive = source_tar(
                source["archive"]["prefix"], source["commit"], members)
            _replace_claimed_file(lane.root, source["archive"], archive)
            source["archive"]["replay_sha256"] = sha256(archive)
            for item in record["adapter"]["files"]:
                if item["path"] == changed_path:
                    item["size"] = len(members[changed_path])
                    item["sha256"] = sha256(members[changed_path])
                    item["git_blob_sha1"] = git_object_id(
                        "blob", members[changed_path])
            record["adapter"]["files_sha256"] = \
                identity_contract.canonical_sha256(record["adapter"]["files"])
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

    def test_build_git_and_census_semantic_adversaries(self) -> None:
        changes = ("cache", "compile", "closure", "git_capture")
        for change in changes:
            with self.subTest(change=change), Lane() as lane:
                lane.writable()
                record = copy.deepcopy(lane.record)
                build = record["builds"]["canonical_first"]
                if change == "cache":
                    path = build["cmake_cache"]["relative_path"]
                    content = (lane.root / path).read_bytes().replace(
                        b"CMAKE_BUILD_TYPE:STRING=Release",
                        b"CMAKE_BUILD_TYPE:STRING=Debug")
                    claim = build["cmake_cache"]
                elif change == "compile":
                    path = build["compile_commands"]["relative_path"]
                    commands = json.loads((lane.root / path).read_bytes())
                    commands.pop()
                    content = identity_contract.canonical_json_bytes(commands)
                    claim = build["compile_commands"]
                elif change == "closure":
                    path = build["closure"]["relative_path"]
                    closure = json.loads((lane.root / path).read_bytes())
                    closure["role"] = "path_variant"
                    content = identity_contract.canonical_json_bytes(closure)
                    claim = build["closure"]
                else:
                    path = record["source"]["adapter_repository"][
                        "git_capture"]["relative_path"]
                    capture = json.loads((lane.root / path).read_bytes())
                    capture["tree"] = "3" * 40
                    content = identity_contract.canonical_json_bytes(capture)
                    claim = record["source"]["adapter_repository"]["git_capture"]
                (lane.root / path).write_bytes(content)
                claim["size"] = len(content)
                claim["sha256"] = sha256(content)
                record = _resign_record(record)
                (lane.root / "baseline-authority.json").write_bytes(
                    identity_contract.canonical_json_bytes(record))
                lane.reseal()
                self.assertInvalid(lane.root)

        with Lane() as lane:
            lane.writable()
            record = copy.deepcopy(lane.record)
            role = "path_variant"
            build = record["builds"][role]
            contaminated = identity_fixtures.synthetic_elf(
                rodata=build["roots"]["build_root"].encode("ascii"),
                debug=b"variant debug bytes differ")
            path = build["executable"]["retained_relative_path"]
            (lane.root / path).write_bytes(contaminated)
            build["executable"]["size"] = len(contaminated)
            build["executable"]["sha256"] = sha256(contaminated)
            record["runtime_closure"]["records"][2]["executable_sha256"] = \
                sha256(contaminated)
            record["identity"][role] = \
                identity_contract.normalized_code_identity_from_elf_bytes(
                    contaminated, roots=build["roots"])
            record["promotion"]["selected_section_census_zero"] = False
            record["promotion"]["promoted"] = False
            record = _resign_record(record)
            (lane.root / "baseline-authority.json").write_bytes(
                identity_contract.canonical_json_bytes(record))
            lane.reseal()
            self.assertInvalid(lane.root)

        for mutation in ("pure_avx2", "parameters"):
            with self.subTest(attestation=mutation), Lane() as lane:
                lane.writable()
                record = copy.deepcopy(lane.record)
                path = "attestations/canonical-first/benchmark.stdout.json"
                value = json.loads((lane.root / path).read_bytes())
                if mutation == "pure_avx2":
                    value["build"]["pure_avx2"] = False
                else:
                    value["parameters"]["K"] += 1
                changed = identity_contract.canonical_json_bytes(value)
                (lane.root / path).write_bytes(changed)
                stdout = record["attestation"]["records"][0]["stdout"]
                stdout["size"] = len(changed)
                stdout["sha256"] = sha256(changed)
                record = _resign_record(record)
                (lane.root / "baseline-authority.json").write_bytes(
                    identity_contract.canonical_json_bytes(record))
                lane.reseal()
                self.assertInvalid(lane.root)

    def test_cli_exit_contract(self) -> None:
        for kind, expected in (("authority", 0), ("acquisition_failure", 3),
                               ("verification_failure", 3)):
            with self.subTest(kind=kind), Lane(kind=kind) as lane:
                completed = subprocess.run(
                    verifier_cli_command(lane.root), check=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(completed.returncode, expected)
                self.assertTrue(completed.stdout)
                self.assertEqual(completed.stderr, b"")
                verdict = identity_contract.strict_json_loads(
                    completed.stdout, "verifier CLI verdict")
                self.assertEqual(verdict["promoted"], expected == 0)
        completed = subprocess.run([
            sys.executable, "-I", "-S", "-B",
            str(TOOLS / "leopard2_exact_main_baseline_verifier.py"), "--help",
        ], check=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 2)
        self.assertEqual(completed.stdout, b"")
        self.assertTrue(completed.stderr.startswith(b"usage:"))
        with Lane() as lane:
            lane.writable()
            (lane.root / "logs/00-source_acquisition.log").write_bytes(b"bad")
            _chmod_tree(lane.root, sealed=True)
            completed = subprocess.run(
                verifier_cli_command(lane.root), check=False,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(completed.returncode, 1)
            self.assertEqual(completed.stdout, b"")
            self.assertTrue(completed.stderr.startswith(
                b"invalid sealed exact-main baseline:"))
        with Lane() as lane:
            reader, writer = os.pipe()
            os.close(reader)
            try:
                process = subprocess.Popen(
                    verifier_cli_command(lane.root), stdout=writer,
                    stderr=subprocess.PIPE, close_fds=True)
            finally:
                os.close(writer)
            _stdout, stderr = process.communicate()
            self.assertEqual(process.returncode, 1)
            self.assertTrue(stderr.startswith(
                b"invalid sealed exact-main baseline:"))

    def test_production_imports_and_host_access_are_read_only(self) -> None:
        source = (TOOLS / "leopard2_exact_main_baseline_verifier.py").read_text()
        parsed = ast.parse(source)
        imports = {
            alias.name.split(".")[0]
            for node in ast.walk(parsed)
            if isinstance(node, ast.Import)
            for alias in node.names
        } | {
            (node.module or "").split(".")[0]
            for node in ast.walk(parsed)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertEqual(imports, {
            "__future__", "base64", "binascii", "copy", "hashlib",
            "importlib", "os", "shlex", "stat", "sys", "typing",
            "leopard2_exact_main_baseline",
            "leopard2_exact_main_baseline_record",
        })
        self.assertFalse(any(
            isinstance(node, ast.Call) and
            isinstance(node.func, ast.Name) and node.func.id == "__import__"
            for node in ast.walk(parsed)))
        real_open = os.open

        def guarded_open(path, flags, *args, **kwargs):
            write_flags = (os.O_WRONLY | os.O_RDWR | os.O_CREAT |
                           os.O_TRUNC | os.O_APPEND)
            if flags & write_flags:
                raise AssertionError("write")
            return real_open(path, flags, *args, **kwargs)

        with Lane(kind="acquisition_failure") as lane, \
             mock.patch.object(os, "chmod", side_effect=AssertionError("write")), \
             mock.patch.object(os, "open", side_effect=guarded_open), \
             mock.patch.object(os, "write", side_effect=AssertionError("write")), \
             mock.patch.object(os, "link", side_effect=AssertionError("write")), \
             mock.patch.object(os, "symlink", side_effect=AssertionError("write")), \
             mock.patch.object(os, "unlink", side_effect=AssertionError("write")), \
             mock.patch.object(os, "remove", side_effect=AssertionError("write")), \
             mock.patch.object(os, "rmdir", side_effect=AssertionError("write")), \
             mock.patch.object(os, "rename", side_effect=AssertionError("write")), \
             mock.patch.object(os, "replace", side_effect=AssertionError("write")), \
             mock.patch.object(os, "truncate", side_effect=AssertionError("write")), \
             mock.patch.object(os, "mkdir", side_effect=AssertionError("write")), \
             mock.patch.object(os, "system", side_effect=AssertionError("process")), \
             mock.patch.object(os, "popen", side_effect=AssertionError("process")), \
             mock.patch.object(os, "posix_spawn", side_effect=AssertionError("process")), \
             mock.patch.object(os, "fork", side_effect=AssertionError("process")), \
             mock.patch.object(subprocess, "run", side_effect=AssertionError("process")), \
             mock.patch.object(subprocess, "Popen", side_effect=AssertionError("process")):
            verdict = verifier.verify_sealed_lane(lane.root)
            self.assertEqual(verdict["outcome"], "verified_failure")


if __name__ == "__main__":
    unittest.main()
