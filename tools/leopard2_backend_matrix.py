#!/usr/bin/env python3
"""Build and compare deterministic Leopard2 x86 backend variants.

The runner uses isolated build trees, pins each test process to an allowed CPU,
and writes resumable, machine-readable results.  It is intentionally standard
library only so that backend verification adds no runtime dependency.
"""

from __future__ import print_function

import argparse
import collections
import concurrent.futures
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCHEMA = "leopard2-backend-matrix/v1"
VARIANTS = ("auto", "scalar", "ssse3", "avx2")
COMPARE_TESTS = (
    "direct_oracle",
    "backend_ops",
    "context_backends",
    "legacy_golden",
    "api",
    "public_api_contract",
    "random",
    "locator",
    "active_lch",
    "gf16_tails",
    "gf16_padded_odd",
    "gf16_legacy_encoder_matrix",
    "low_gf16_direct_rows",
    "decode_high_acceptance",
    "decode_low_acceptance",
    "decode_plan_schedule",
    "direct_encode",
    "arbitrary_counts_acceptance",
    "max_counts",
    "encode_concurrency",
    "codec_options_abi",
    "direct_repair",
    "boundaries",
    "transform_differential",
    "fuzz_smoke",
)
BUILD_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE", "CMAKE_GENERATOR",
    "CMAKE_C_FLAGS", "CMAKE_C_FLAGS_RELEASE",
    "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
    "CMAKE_STATIC_LINKER_FLAGS", "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT", "LEO2_BUILD_TESTS",
    "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS", "LEO2_ENABLE_CUDA",
)

EXPECTED_COMPILE_SOURCE_COUNTS = {
    "Leopard2Backend.cpp": 1,
    "Leopard2BackendAVX2.cpp": 1,
    "Leopard2BackendSSSE3.cpp": 1,
    "Leopard2BackendScalar.cpp": 1,
    "Leopard2CpuFeatures.cpp": 1,
    "Leopard2Plan.cpp": 1,
    "LeopardCommon.cpp": 1,
    "LeopardFF16.cpp": 1,
    "LeopardFF8.cpp": 1,
    "leopard.cpp": 1,
    "leopard2.cpp": 1,
    "tests/leopard2/direct_oracle.cpp": 13,
    "tests/leopard2/direct_repair.cpp": 1,
    "tests/leopard2/fuzz_api.cpp": 1,
    "tests/leopard2/fuzz_replay.cpp": 1,
    "tests/leopard2/test_active_lch.cpp": 1,
    "tests/leopard2/test_api.cpp": 1,
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp": 1,
    "tests/leopard2/test_backend_failures.cpp": 1,
    "tests/leopard2/test_backend_ops.cpp": 1,
    "tests/leopard2/test_context_backends.cpp": 1,
    "tests/leopard2/test_boundaries.cpp": 1,
    "tests/leopard2/test_codec_options_abi.c": 1,
    "tests/leopard2/test_decode_high_acceptance.cpp": 1,
    "tests/leopard2/test_decode_low_acceptance.cpp": 1,
    "tests/leopard2/test_decode_plan_schedule.cpp": 1,
    "tests/leopard2/test_direct_encode.cpp": 1,
    "tests/leopard2/test_direct_oracle.cpp": 1,
    "tests/leopard2/test_direct_repair.cpp": 1,
    "tests/leopard2/test_encode_concurrency.cpp": 1,
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp": 1,
    "tests/leopard2/test_gf16_padded_odd.cpp": 1,
    "tests/leopard2/test_gf16_tails.cpp": 1,
    "tests/leopard2/test_legacy_golden.cpp": 1,
    "tests/leopard2/test_locator.cpp": 1,
    "tests/leopard2/test_low_gf16_direct_rows.cpp": 1,
    "tests/leopard2/test_max_counts.cpp": 1,
    "tests/leopard2/test_public_api_contract.cpp": 1,
    "tests/leopard2/test_random.cpp": 1,
    "tests/leopard2/test_transform_differential.cpp": 1,
}
SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "Leopard2CpuFeatures.cpp",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "Leopard2Direct.h",
    "Leopard2Dispatch.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "tests/leopard2/test_legacy_golden.cpp",
    "tests/leopard2/test_backend_failures.cpp",
    "tests/leopard2/test_backend_ops.cpp",
    "tests/leopard2/test_context_backends.cpp",
    "tests/leopard2/legacy_golden_vectors.h",
    "tests/leopard2/test_api.cpp",
    "tests/leopard2/test_public_api_contract.cpp",
    "tests/leopard2/test_random.cpp",
    "tests/leopard2/test_locator.cpp",
    "tests/leopard2/test_boundaries.cpp",
    "tests/leopard2/test_active_lch.cpp",
    "tests/leopard2/test_gf16_tails.cpp",
    "tests/leopard2/test_gf16_padded_odd.cpp",
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp",
    "tests/leopard2/test_low_gf16_direct_rows.cpp",
    "tests/leopard2/test_decode_high_acceptance.cpp",
    "tests/leopard2/test_decode_low_acceptance.cpp",
    "tests/leopard2/test_decode_plan_schedule.cpp",
    "tests/leopard2/test_direct_encode.cpp",
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp",
    "tests/leopard2/test_max_counts.cpp",
    "tests/leopard2/test_encode_concurrency.cpp",
    "tests/leopard2/test_codec_options_abi.c",
    "tests/leopard2/test_transform_differential.cpp",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
    "tests/leopard2/test_direct_oracle.cpp",
    "tests/leopard2/direct_repair.cpp",
    "tests/leopard2/direct_repair.h",
    "tests/leopard2/test_direct_repair.cpp",
    "tests/leopard2/fuzz_api.cpp",
    "tests/leopard2/fuzz_replay.cpp",
    "tests/cmake/test_cuda_optional.cmake",
    "cmake/leopardConfig.cmake.in",
    "tools/check_leopard2_portable_isa.sh",
    "tools/leopard2_backend_matrix.py",
)


class MatrixError(Exception):
    """An actionable matrix configuration or execution error."""


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical_bytes(value))


def normalized_output(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def named_ctest_executed(stdout, test_name, stderr=b""):
    output = normalized_output(stdout + stderr)
    return b"Test #" in output and test_name.encode("ascii") in output


def portable_ctest_executed(stdout, stderr=b""):
    return named_ctest_executed(stdout, "leopard2_portable_isa", stderr)


def cuda_ctest_executed(stdout, stderr=b""):
    return named_ctest_executed(stdout, "leopard2_cuda_optional", stderr)


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, str(path))
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def compact_cpu_list(cpus):
    values = sorted(set(int(cpu) for cpu in cpus))
    if not values:
        return ""
    result = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        result.append(str(first) if first == previous else "{}-{}".format(first, previous))
        first = previous = value
    result.append(str(first) if first == previous else "{}-{}".format(first, previous))
    return ",".join(result)


def allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            cpus = sorted(os.sched_getaffinity(0))
            if cpus:
                return cpus
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def cpu_flags_for(allowed):
    try:
        text = Path("/proc/cpuinfo").read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return []
    wanted = set(allowed)
    fallback = []
    for block in text.split("\n\n"):
        fields = {}
        for line in block.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key.strip().lower()] = value.strip()
        flags = fields.get("flags", fields.get("features", "")).split()
        if flags and not fallback:
            fallback = flags
        try:
            processor = int(fields.get("processor", "-1"))
        except ValueError:
            processor = -1
        if processor in wanted:
            return sorted(set(flags))
    return sorted(set(fallback))


def compiler_identity(compiler):
    resolved = shutil.which(compiler)
    if not resolved:
        raise MatrixError("C++ compiler not found: {}".format(compiler))
    completed = subprocess.run(
        [resolved, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output = normalized_output(completed.stdout + completed.stderr).decode(
        "utf-8", errors="replace"
    ).strip()
    if completed.returncode != 0:
        raise MatrixError("cannot query compiler {}: {}".format(resolved, output))
    return {
        "executable": str(Path(resolved).resolve()),
        "binary_sha256": digest_bytes(Path(resolved).read_bytes()),
        "version": output,
        "version_sha256": digest_bytes(output.encode("utf-8")),
    }


def variant_flags(variant):
    if variant == "ssse3":
        return ["-mssse3", "-mno-avx"]
    if variant == "avx2":
        return ["-mavx2", "-mno-avx512f"]
    return []


def compiler_accepts(compiler, flags):
    if not flags:
        return True, ""
    completed = subprocess.run(
        [compiler, "-x", "c++", "-std=c++11", "-fsyntax-only", "-"] + flags,
        input=b"int main() { return 0; }\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    message = normalized_output(completed.stdout + completed.stderr).decode(
        "utf-8", errors="replace"
    ).strip()
    return completed.returncode == 0, message


def availability(variant, machine, compiler):
    if variant == "auto":
        return True, ""
    architecture = machine["architecture"].lower()
    if architecture not in ("x86_64", "amd64", "i386", "i486", "i586", "i686"):
        return False, "forced variants are x86-only on this implementation"
    if variant == "scalar":
        return True, ""
    required = "avx2" if variant == "avx2" else "ssse3"
    if required not in machine["cpu_flags"]:
        return False, "host CPU does not advertise {}".format(required)
    supported, message = compiler_accepts(compiler["executable"], variant_flags(variant))
    if not supported:
        detail = ": " + message if message else ""
        return False, "compiler rejects {}{}".format(" ".join(variant_flags(variant)), detail)
    return True, ""


def source_fingerprint(source):
    files = {}
    for relative in SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise MatrixError("required matrix input is missing: {}".format(path))
        files[relative] = digest_bytes(path.read_bytes())
    return {"digest": digest_value(files), "files": files}


def executable_path(build, name):
    direct = build / name
    if direct.is_file():
        return direct
    release = build / "Release" / (name + (".exe" if os.name == "nt" else ""))
    if release.is_file():
        return release
    windows = build / (name + ".exe")
    if windows.is_file():
        return windows
    raise MatrixError("built executable not found: {}".format(name))


def read_cache_variant(build):
    cache = build / "CMakeCache.txt"
    try:
        lines = cache.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise MatrixError("cannot read {}: {}".format(cache, error))
    prefix = "LEO2_BACKEND_VARIANT:STRING="
    for line in lines:
        if line.startswith(prefix):
            return line[len(prefix):]
    raise MatrixError("LEO2_BACKEND_VARIANT missing from {}".format(cache))


def seal_command(record):
    payload = dict(record)
    payload.pop("command_sha256", None)
    record["command_sha256"] = digest_value(payload)
    return record


def run_command(label, argv, cwd, log_dir, timeout, environment=None, hash_output=False):
    timed_out = False
    try:
        completed = subprocess.run(
            [str(item) for item in argv],
            cwd=str(cwd),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        returncode = completed.returncode
        stdout = normalized_output(completed.stdout)
        stderr = normalized_output(completed.stderr)
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = 124
        stdout = normalized_output(error.stdout or b"")
        stderr = normalized_output(error.stderr or b"")
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_name = label + ".stdout.log"
    stderr_name = label + ".stderr.log"
    (log_dir / stdout_name).write_bytes(stdout)
    (log_dir / stderr_name).write_bytes(stderr)
    record = {
        "argv": [str(item) for item in argv],
        "cwd": str(cwd),
        "label": label,
        "returncode": returncode,
        "stderr_log": stderr_name,
        "stdout_log": stdout_name,
        "timed_out": timed_out,
    }
    record["stderr_sha256"] = digest_bytes(stderr)
    record["stdout_sha256"] = digest_bytes(stdout)
    return record


def prepare_fresh_directory(path):
    path = Path(path)
    if path.exists():
        if path.is_symlink() or not path.is_dir():
            raise MatrixError("refusing non-directory build path: {}".format(path))
        shutil.rmtree(str(path))
    path.mkdir(parents=True, exist_ok=False)
    return True


def parse_cache(path):
    values = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        values[key] = value
    return values


def normalized_build_identity(build, source, c_compiler, compiler, cmake, ctest):
    cache_path = build / "CMakeCache.txt"
    compile_path = build / "compile_commands.json"
    if not cache_path.is_file() or not compile_path.is_file():
        raise MatrixError("fresh build omitted cache or compile commands")
    cache = parse_cache(cache_path)
    normalized_cache = {}
    for key in BUILD_CACHE_KEYS:
        if key not in cache:
            raise MatrixError("fresh cache omits {}".format(key))
        normalized_cache[key] = cache[key]
    compiler_cache = {
        "CMAKE_C_COMPILER": (c_compiler, "@tool/cc"),
        "CMAKE_CXX_COMPILER": (compiler, "@tool/cxx"),
    }
    for key, (identity_record, tag) in compiler_cache.items():
        cached = Path(cache.get(key, "")).resolve()
        if cached != Path(identity_record["executable"]).resolve():
            raise MatrixError("fresh cache/compiler mismatch: {}".format(key))
        normalized_cache[key] = tag
    commands = json.loads(compile_path.read_text(encoding="utf-8"))
    if not isinstance(commands, list) or not commands:
        raise MatrixError("compile_commands.json is empty or malformed")
    normalized_commands = []
    seen = set()
    replacements = ((str(source.resolve()), "@source"),
                    (str(build.resolve()), "@build"),
                    (str(Path(c_compiler["executable"]).resolve()), "@tool/cc"),
                    (str(Path(compiler["executable"]).resolve()), "@tool/cxx"))
    for entry in commands:
        path = Path(entry.get("file", ""))
        if not path.is_absolute():
            path = Path(entry.get("directory", build)) / path
        path = path.resolve()
        try:
            relative = path.relative_to(source.resolve()).as_posix()
        except ValueError:
            raise MatrixError("compile source escapes source root: {}".format(path))
        raw = entry.get("arguments")
        if raw is None:
            raw = shlex.split(entry.get("command", ""))
        argv = [str(value) for value in raw]
        if not argv:
            raise MatrixError("empty compile command: {}".format(relative))
        expected_compiler = c_compiler if relative.endswith(".c") else compiler
        executable = Path(argv[0])
        if not executable.is_absolute():
            resolved = shutil.which(argv[0], path="/usr/bin:/bin")
            if not resolved:
                raise MatrixError("compile executable not found: {}".format(argv[0]))
            executable = Path(resolved)
        if executable.resolve() != Path(expected_compiler["executable"]).resolve():
            raise MatrixError("compile language/tool mismatch: {}".format(relative))
        normalized_argv = []
        for argument in argv:
            for actual, logical in replacements:
                argument = argument.replace(actual, logical)
            normalized_argv.append(argument)
        normalized_argv[0] = "@tool/cc" if relative.endswith(".c") else "@tool/cxx"
        source_tag = "@source/" + relative
        if normalized_argv.count(source_tag) != 1:
            raise MatrixError("compile command source mismatch: {}".format(relative))
        language = "C" if relative.endswith(".c") else "CXX"
        identity = (relative, tuple(normalized_argv))
        if identity in seen:
            raise MatrixError("duplicate compile command: {}".format(relative))
        seen.add(identity)
        normalized_commands.append({"file": relative, "language": language,
                                    "argv": normalized_argv})
    normalized_commands.sort(key=lambda value: (value["file"], value["argv"]))
    counts = collections.Counter(value["file"] for value in normalized_commands)
    if dict(counts) != EXPECTED_COMPILE_SOURCE_COUNTS:
        raise MatrixError("compile-command source multiset mismatch")

    def identity(path, name):
        resolved = Path(path).resolve()
        completed = subprocess.run(
            [str(resolved), "--version"], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        raw = normalized_output(completed.stdout + completed.stderr)
        if completed.returncode != 0 or not raw:
            raise MatrixError("cannot identify {}".format(name))
        return {"basename": resolved.name,
                "binary_sha256": digest_bytes(resolved.read_bytes()),
                "version_sha256": digest_bytes(raw)}

    result = {
        "cache": normalized_cache,
        "cache_sha256": digest_value(normalized_cache),
        "compile_commands": normalized_commands,
        "compile_commands_sha256": digest_value(normalized_commands),
        "test_executables": {},
        "tools": {
            "cmake": identity(cmake, "cmake"),
            "ctest": identity(ctest, "ctest"),
            "cxx": {"basename": Path(compiler["executable"]).name,
                    "binary_sha256": compiler["binary_sha256"],
                    "version_sha256": compiler["version_sha256"]},
            "cc": {"basename": Path(c_compiler["executable"]).name,
                   "binary_sha256": c_compiler["binary_sha256"],
                   "version_sha256": c_compiler["version_sha256"]},
        },
    }
    result["digest"] = digest_value(result)
    return result


def run_variant(context, variant, index):
    result_dir = context["result_dir"] / variant
    result_path = context["result_dir"] / (variant + ".json")
    build = context["build_root"] / variant
    available, reason = availability(variant, context["machine"], context["compiler"])
    environment = {
        "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
        "OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin",
    }
    if variant != "auto":
        environment["LEO2_EXPECT_BACKEND"] = variant
    identity_input = {
        "c_compiler": context["c_compiler"],
        "compiler": context["compiler"],
        "generator": context["generator"],
        "jobs_per_variant": context["jobs_per_variant"],
        "machine": context["machine"],
        "source": context["source_fingerprint"],
        "variant": variant,
        "environment": environment,
    }
    configuration_id = digest_value(identity_input)

    if context["resume"] and result_path.is_file():
        try:
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            if (previous.get("configuration_id") == configuration_id and
                    previous.get("status") in ("passed", "unavailable")):
                return previous
        except (OSError, UnicodeError, ValueError):
            pass

    base = {
        "configuration_id": configuration_id,
        "resumed": False,
        "schema": SCHEMA,
        "source_fingerprint": context["source_fingerprint"]["digest"],
        "variant": variant,
    }
    if not available:
        base.update({"commands": [], "reason": reason, "status": "unavailable", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    # A source fingerprint is not enough if an older build tree can retain an
    # object whose mtime is newer than the checkout.  Every non-resumed run is
    # therefore configured from an actually empty per-variant directory.
    prepare_fresh_directory(build)
    if result_dir.exists():
        shutil.rmtree(str(result_dir))
    result_dir.mkdir(parents=True, exist_ok=False)
    commands = []
    base["build_environment"] = dict(environment)
    configure = [
        context["cmake"], "-S", context["source"], "-B", build,
        "-G", context["generator"],
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_COMPILER={}".format(context["c_compiler"]["executable"]),
        "-DCMAKE_CXX_COMPILER={}".format(context["compiler"]["executable"]),
        "-DLEO2_BACKEND_VARIANT={}".format(variant),
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_BUILD_BENCHMARKS=OFF",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
    ]
    base["fresh_build"] = {
        "configured_from_empty": True,
        "identity_sha256": digest_value({
            "configuration_id": configuration_id,
            "configure_argv": [str(value) for value in configure],
            "environment": environment,
        }),
    }
    command = run_command(
        "configure", configure, context["source"], result_dir,
        context["timeout"], environment
    )
    seal_command(command)
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "configure failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    base["build_identity"] = normalized_build_identity(
        build, context["source"], context["c_compiler"], context["compiler"],
        context["cmake"], context["ctest"])

    selected = read_cache_variant(build)
    if selected != variant:
        base.update({
            "commands": commands,
            "reason": "CMake selected {!r}, expected {!r}".format(selected, variant),
            "status": "failed",
            "tests": {},
        })
        atomic_write_json(result_path, base)
        return base

    targets = [
        "leopard2_direct_oracle_test",
        "leopard2_backend_ops_test", "leopard2_backend_failures_test",
        "leopard2_context_backends_test",
        "leopard2_legacy_golden_test", "leopard2_api_test",
        "leopard2_public_api_contract_test",
        "leopard2_random_test", "leopard2_locator_test",
        "leopard2_active_lch_test", "leopard2_gf16_tails_test",
        "leopard2_gf16_padded_odd_test",
        "leopard2_gf16_legacy_encoder_matrix_test",
        "leopard2_low_gf16_direct_rows_test",
        "leopard2_decode_high_acceptance_test",
        "leopard2_decode_low_acceptance_test",
        "leopard2_decode_plan_schedule_test",
        "leopard2_direct_encode_test",
        "leopard2_arbitrary_counts_acceptance_test",
        "leopard2_max_counts_test",
        "leopard2_encode_concurrency_test", "leopard2_codec_options_abi_test",
        "leopard2_direct_repair_test", "leopard2_boundaries_test",
        "leopard2_transform_differential_test", "leopard2_fuzz_smoke",
    ]
    build_command = [
        context["cmake"], "--build", build, "--config", "Release",
        "-j", str(context["jobs_per_variant"]), "--target",
    ] + targets
    command = run_command(
        "build", build_command, context["source"], result_dir,
        context["timeout"], environment
    )
    seal_command(command)
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "build failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    test_specs = {
        "direct_oracle": ("leopard2_direct_oracle_test", []),
        "backend_ops": ("leopard2_backend_ops_test", []),
        "context_backends": ("leopard2_context_backends_test", []),
        "legacy_golden": ("leopard2_legacy_golden_test", []),
        "api": ("leopard2_api_test", []),
        "public_api_contract": ("leopard2_public_api_contract_test", []),
        "random": ("leopard2_random_test", [
            "--seed", "0x4c656f7061726432", "--cases", "64", "--threads", "1"
        ]),
        "locator": ("leopard2_locator_test", []),
        "active_lch": ("leopard2_active_lch_test", []),
        "gf16_tails": ("leopard2_gf16_tails_test", []),
        "gf16_padded_odd": ("leopard2_gf16_padded_odd_test", []),
        "gf16_legacy_encoder_matrix": (
            "leopard2_gf16_legacy_encoder_matrix_test", []),
        "low_gf16_direct_rows": ("leopard2_low_gf16_direct_rows_test", []),
        "decode_high_acceptance": (
            "leopard2_decode_high_acceptance_test", []),
        "decode_low_acceptance": (
            "leopard2_decode_low_acceptance_test", []),
        "decode_plan_schedule": ("leopard2_decode_plan_schedule_test", []),
        "direct_encode": ("leopard2_direct_encode_test", []),
        "arbitrary_counts_acceptance": (
            "leopard2_arbitrary_counts_acceptance_test", []),
        "max_counts": ("leopard2_max_counts_test", []),
        "encode_concurrency": ("leopard2_encode_concurrency_test", []),
        "codec_options_abi": ("leopard2_codec_options_abi_test", []),
        "direct_repair": ("leopard2_direct_repair_test", []),
        "boundaries": ("leopard2_boundaries_test", []),
        "transform_differential": ("leopard2_transform_differential_test", []),
        "fuzz_smoke": ("leopard2_fuzz_smoke", []),
    }
    tests = {}
    pin_cpu = context["allowed_cpus"][index % len(context["allowed_cpus"])]
    for name in COMPARE_TESTS:
        target, arguments = test_specs[name]
        executable = executable_path(build, target)
        argv = [str(executable)] + arguments
        if context["taskset"]:
            argv = [context["taskset"], "-c", str(pin_cpu)] + argv
        test_environment = environment
        if name == "context_backends":
            test_environment = dict(environment)
            test_environment["OMP_DYNAMIC"] = "FALSE"
            test_environment["OMP_NUM_THREADS"] = "4"
        command = run_command(
            "test_" + name, argv, context["source"], result_dir,
            context["timeout"], test_environment, hash_output=True
        )
        command["executable_sha256"] = digest_bytes(executable.read_bytes())
        base["build_identity"]["test_executables"][name] = {
            "path": "@build/" + executable.relative_to(build).as_posix(),
            "sha256": command["executable_sha256"],
        }
        seal_command(command)
        tests[name] = command
        commands.append(command)
        if command["returncode"] != 0:
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": "{} test failed".format(name),
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    failure_command = [
        context["ctest"], "--test-dir", build, "-C", "Release",
        "-R", "^leopard2_backend_failure_", "--output-on-failure",
    ]
    command = run_command(
        "test_backend_failures", failure_command, context["source"],
        result_dir, context["timeout"], environment, hash_output=True
    )
    tests["backend_failures"] = command
    commands.append(command)
    failure_was_run = named_ctest_executed(
        (result_dir / command["stdout_log"]).read_bytes(),
        "leopard2_backend_failure_",
        (result_dir / command["stderr_log"]).read_bytes(),
    )
    command["ctest_executed"] = failure_was_run
    seal_command(command)
    if command["returncode"] != 0 or not failure_was_run:
        base.update({
            "commands": commands,
            "pin_cpu": pin_cpu,
            "reason": (
                "backend failure matrix failed" if command["returncode"] != 0
                else "backend failure tests were not registered or executed"
            ),
            "selected_cache_variant": selected,
            "status": "failed",
            "tests": tests,
        })
        atomic_write_json(result_path, base)
        return base

    if variant in VARIANTS:
        portable_command = [
            context["ctest"], "--test-dir", build, "-C", "Release",
            "-R", "^leopard2_portable_isa$", "--output-on-failure",
        ]
        command = run_command(
            "test_portable_isa", portable_command, context["source"], result_dir,
            context["timeout"], environment, hash_output=True
        )
        tests["portable_isa"] = command
        commands.append(command)
        portable_was_run = portable_ctest_executed(
            (result_dir / command["stdout_log"]).read_bytes(),
            (result_dir / command["stderr_log"]).read_bytes(),
        )
        command["ctest_executed"] = portable_was_run
        seal_command(command)
        if command["returncode"] != 0 or not portable_was_run:
            if command["returncode"] == 0:
                reason = (
                    "portable-ISA test was not registered or executed; install "
                    "objdump or llvm-objdump and a POSIX sh, then reconfigure"
                )
            else:
                reason = "portable-ISA test failed"
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": reason,
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    if variant == "auto":
        cuda_command = [
            context["ctest"], "--test-dir", build, "-C", "Release",
            "-R", "^leopard2_cuda_optional$", "--output-on-failure",
        ]
        command = run_command(
            "test_cuda_optional", cuda_command, context["source"], result_dir,
            context["timeout"], environment, hash_output=True
        )
        tests["cuda_optional"] = command
        commands.append(command)
        cuda_was_run = cuda_ctest_executed(
            (result_dir / command["stdout_log"]).read_bytes(),
            (result_dir / command["stderr_log"]).read_bytes(),
        )
        command["ctest_executed"] = cuda_was_run
        seal_command(command)
        if command["returncode"] != 0 or not cuda_was_run:
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": (
                    "optional-CUDA test failed" if command["returncode"] != 0
                    else "optional-CUDA test was not registered or executed"
                ),
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    identity_payload = dict(base["build_identity"])
    identity_payload.pop("digest", None)
    base["build_identity"]["digest"] = digest_value(identity_payload)
    base.update({
        "commands": commands,
        "expected_runtime_backend": None if variant == "auto" else variant,
        "pin_cpu": pin_cpu,
        "reason": "",
        "selected_cache_variant": selected,
        "status": "passed",
        "tests": tests,
    })
    atomic_write_json(result_path, base)
    return base


def compare_results(results):
    passed = {result["variant"]: result for result in results if result["status"] == "passed"}
    mismatches = []
    if "auto" not in passed:
        return [{"reason": "auto variant did not pass"}]
    for variant in sorted(passed):
        if variant == "auto":
            continue
        for test in COMPARE_TESTS:
            for stream in ("stdout", "stderr"):
                key = stream + "_sha256"
                expected = passed["auto"]["tests"][test][key]
                actual = passed[variant]["tests"][test][key]
                if actual != expected:
                    mismatches.append({
                        "actual": actual,
                        "expected": expected,
                        "stream": stream,
                        "test": test,
                        "variant": variant,
                    })
    return mismatches


def matrix_run(arguments):
    source = Path(arguments.source).resolve()
    build_root = Path(arguments.build_root)
    if not build_root.is_absolute():
        build_root = source / build_root
    result_dir = Path(arguments.result_dir)
    if not result_dir.is_absolute():
        result_dir = source / result_dir
    requested = []
    for value in arguments.variants.split(","):
        variant = value.strip().lower()
        if variant and variant not in requested:
            requested.append(variant)
    invalid = [variant for variant in requested if variant not in VARIANTS]
    if invalid or not requested:
        raise MatrixError("variants must be a comma-separated subset of {}".format(",".join(VARIANTS)))

    cpus = allowed_cpus()
    jobs = min(max(1, arguments.jobs), len(cpus), 128)
    workers = min(max(1, arguments.variant_workers), len(requested), jobs)
    c_compiler = compiler_identity(arguments.c_compiler)
    compiler = compiler_identity(arguments.compiler)
    machine = {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_flags": cpu_flags_for(cpus),
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
    }
    source_state = source_fingerprint(source)
    cmake = shutil.which(arguments.cmake)
    ctest = shutil.which(arguments.ctest)
    if not cmake or not ctest:
        raise MatrixError("cmake and ctest must both be available")
    context = {
        "allowed_cpus": cpus,
        "build_root": build_root,
        "cmake": str(Path(cmake).resolve()),
        "c_compiler": c_compiler,
        "compiler": compiler,
        "ctest": str(Path(ctest).resolve()),
        "generator": arguments.generator,
        "jobs_per_variant": max(1, jobs // workers),
        "machine": machine,
        "result_dir": result_dir,
        "resume": not arguments.no_resume,
        "source": source,
        "source_fingerprint": source_state,
        "taskset": shutil.which("taskset") if os.name == "posix" else None,
        "timeout": arguments.timeout,
    }
    print(
        "backend matrix: variants={} workers={} jobs/variant={} cpus={}".format(
            ",".join(requested), workers, context["jobs_per_variant"],
            machine["allowed_cpu_list"]
        ),
        flush=True,
    )
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(run_variant, context, variant, index): variant
            for index, variant in enumerate(requested)
        }
        for future in concurrent.futures.as_completed(futures):
            variant = futures[future]
            result = future.result()
            results.append(result)
            print("{}: {}".format(variant, result["status"]), flush=True)
    results.sort(key=lambda item: VARIANTS.index(item["variant"]))

    final_state = source_fingerprint(source)
    source_changed = final_state["digest"] != source_state["digest"]
    mismatches = compare_results(results)
    failed = [result["variant"] for result in results if result["status"] == "failed"]
    summary = {
        "c_compiler": c_compiler,
        "compiler": compiler,
        "generator": arguments.generator,
        "jobs": jobs,
        "jobs_per_variant": context["jobs_per_variant"],
        "machine": machine,
        "mismatches": mismatches,
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "source_fingerprint": source_state,
        "status": "failed" if failed or mismatches or source_changed else "passed",
        "variant_workers": workers,
        "variants": results,
    }
    atomic_write_json(result_dir / "matrix.json", summary)
    print("matrix: {} ({})".format(summary["status"], result_dir / "matrix.json"), flush=True)
    return 0 if summary["status"] == "passed" else 1


def self_test():
    assert compact_cpu_list([3, 2, 1, 7, 9, 8]) == "1-3,7-9"
    assert compact_cpu_list([]) == ""
    assert normalized_output(b"a\r\nb\rc\n") == b"a\nb\nc\n"
    assert portable_ctest_executed(
        b"1/1 Test #1: leopard2_portable_isa ... Passed\n"
    )
    assert not portable_ctest_executed(b"No tests were found!!!\n")
    assert cuda_ctest_executed(
        b"1/1 Test #2: leopard2_cuda_optional ... Passed\n"
    )
    assert not cuda_ctest_executed(b"No tests were found!!!\n")
    assert digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2})
    assert variant_flags("scalar") == []
    assert variant_flags("ssse3") == ["-mssse3", "-mno-avx"]
    assert availability(
        "scalar",
        {"architecture": "x86_64", "cpu_flags": []},
        {"executable": "not-needed-for-empty-flags"},
    ) == (True, "")
    with tempfile.TemporaryDirectory(prefix="leo2-backend-self-test-") as directory:
        path = Path(directory) / "result.json"
        value = {"z": [3, 2, 1], "a": "stable"}
        atomic_write_json(path, value)
        assert json.loads(path.read_text(encoding="utf-8")) == value
        stale = Path(directory) / "stale-build"
        stale.mkdir()
        (stale / "copied-object.o").write_bytes(b"stale")
        assert prepare_fresh_directory(stale)
        assert list(stale.iterdir()) == []
        timeout_record = run_command(
            "timeout", [sys.executable, "-c", "import time; time.sleep(1)"],
            directory, Path(directory), 0.01
        )
        assert timeout_record["returncode"] == 124
        assert timeout_record["timed_out"] is True
    print("leopard2 backend matrix self-test passed")
    return 0


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command")
    run = subparsers.add_parser("run", help="build and compare backend variants")
    default_source = str(Path(__file__).resolve().parents[1])
    run.add_argument("--source", default=default_source)
    run.add_argument("--build-root", default="build/leopard2-backend-matrix")
    run.add_argument("--result-dir", default="results/leopard2/backend-matrix")
    run.add_argument("--variants", default=",".join(VARIANTS))
    run.add_argument("--jobs", type=int, default=min(os.cpu_count() or 1, 128))
    run.add_argument("--variant-workers", type=int, default=len(VARIANTS))
    run.add_argument("--timeout", type=int, default=900)
    run.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    run.add_argument("--c-compiler", default=os.environ.get("CC", "cc"))
    run.add_argument("--cmake", default="cmake")
    run.add_argument("--ctest", default="ctest")
    run.add_argument("--generator", default="Ninja" if shutil.which("ninja") else "Unix Makefiles")
    run.add_argument("--no-resume", action="store_true")
    subparsers.add_parser("self-test", help="run deterministic utility tests")
    return result


def main():
    arguments = parser().parse_args()
    try:
        if arguments.command == "self-test":
            return self_test()
        if arguments.command == "run":
            if arguments.jobs <= 0 or arguments.variant_workers <= 0 or arguments.timeout <= 0:
                raise MatrixError("jobs, variant-workers, and timeout must be positive")
            return matrix_run(arguments)
        parser().print_help()
        return 2
    except (MatrixError, OSError, subprocess.SubprocessError) as error:
        print("backend matrix error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
