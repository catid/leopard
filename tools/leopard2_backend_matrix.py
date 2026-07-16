#!/usr/bin/env python3
"""Build and compare deterministic Leopard2 x86 backend variants.

The runner uses isolated build trees, pins each test process to an allowed CPU,
and writes resumable, machine-readable results.  It is intentionally standard
library only so that backend verification adds no runtime dependency.
"""

from __future__ import print_function

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCHEMA = "leopard2-backend-matrix/v1"
VARIANTS = ("auto", "scalar", "ssse3", "avx2")
COMPARE_TESTS = (
    "legacy_golden",
    "api",
    "random",
    "active_lch",
    "gf16_padded_odd",
    "gf16_legacy_encoder_matrix",
    "low_gf16_direct_rows",
    "decode_high_acceptance",
    "decode_low_acceptance",
    "max_counts",
    "encode_concurrency",
    "codec_options_abi",
    "transform_differential",
)
SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "Leopard2Dispatch.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "tests/leopard2/test_legacy_golden.cpp",
    "tests/leopard2/legacy_golden_vectors.h",
    "tests/leopard2/test_api.cpp",
    "tests/leopard2/test_random.cpp",
    "tests/leopard2/test_boundaries.cpp",
    "tests/leopard2/test_active_lch.cpp",
    "tests/leopard2/test_gf16_padded_odd.cpp",
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp",
    "tests/leopard2/test_low_gf16_direct_rows.cpp",
    "tests/leopard2/test_decode_high_acceptance.cpp",
    "tests/leopard2/test_decode_low_acceptance.cpp",
    "tests/leopard2/test_max_counts.cpp",
    "tests/leopard2/test_encode_concurrency.cpp",
    "tests/leopard2/test_codec_options_abi.c",
    "tests/leopard2/test_transform_differential.cpp",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
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
        "version": output,
        "version_sha256": digest_bytes(output.encode("utf-8")),
    }


def variant_flags(variant):
    if variant in ("scalar", "ssse3"):
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
    if hash_output:
        record["stderr_sha256"] = digest_bytes(stderr)
        record["stdout_sha256"] = digest_bytes(stdout)
    return record


def run_variant(context, variant, index):
    result_dir = context["result_dir"] / variant
    result_path = context["result_dir"] / (variant + ".json")
    build = context["build_root"] / variant
    available, reason = availability(variant, context["machine"], context["compiler"])
    identity_input = {
        "compiler": context["compiler"],
        "generator": context["generator"],
        "jobs_per_variant": context["jobs_per_variant"],
        "machine": context["machine"],
        "source": context["source_fingerprint"],
        "variant": variant,
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

    build.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    commands = []
    environment = os.environ.copy()
    environment["CMAKE_BUILD_PARALLEL_LEVEL"] = str(context["jobs_per_variant"])
    environment["OMP_NUM_THREADS"] = "1"
    if variant == "auto":
        environment.pop("LEO2_EXPECT_BACKEND", None)
    else:
        environment["LEO2_EXPECT_BACKEND"] = variant
    configure = [
        context["cmake"], "-S", context["source"], "-B", build,
        "-G", context["generator"],
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_CXX_COMPILER={}".format(context["compiler"]["executable"]),
        "-DLEO2_BACKEND_VARIANT={}".format(variant),
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_BUILD_BENCHMARKS=OFF",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
    ]
    command = run_command(
        "configure", configure, context["source"], result_dir,
        context["timeout"], environment
    )
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "configure failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

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
        "leopard2_legacy_golden_test", "leopard2_api_test",
        "leopard2_random_test", "leopard2_active_lch_test",
        "leopard2_gf16_padded_odd_test",
        "leopard2_gf16_legacy_encoder_matrix_test",
        "leopard2_low_gf16_direct_rows_test",
        "leopard2_decode_high_acceptance_test",
        "leopard2_decode_low_acceptance_test", "leopard2_max_counts_test",
        "leopard2_encode_concurrency_test", "leopard2_codec_options_abi_test",
        "leopard2_transform_differential_test",
    ]
    build_command = [
        context["cmake"], "--build", build, "--config", "Release",
        "-j", str(context["jobs_per_variant"]), "--target",
    ] + targets
    command = run_command(
        "build", build_command, context["source"], result_dir,
        context["timeout"], environment
    )
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "build failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    test_specs = {
        "legacy_golden": ("leopard2_legacy_golden_test", []),
        "api": ("leopard2_api_test", []),
        "random": ("leopard2_random_test", [
            "--seed", "0x4c656f7061726432", "--cases", "64", "--threads", "1"
        ]),
        "active_lch": ("leopard2_active_lch_test", []),
        "gf16_padded_odd": ("leopard2_gf16_padded_odd_test", []),
        "gf16_legacy_encoder_matrix": (
            "leopard2_gf16_legacy_encoder_matrix_test", []),
        "low_gf16_direct_rows": ("leopard2_low_gf16_direct_rows_test", []),
        "decode_high_acceptance": (
            "leopard2_decode_high_acceptance_test", []),
        "decode_low_acceptance": (
            "leopard2_decode_low_acceptance_test", []),
        "max_counts": ("leopard2_max_counts_test", []),
        "encode_concurrency": ("leopard2_encode_concurrency_test", []),
        "codec_options_abi": ("leopard2_codec_options_abi_test", []),
        "transform_differential": ("leopard2_transform_differential_test", []),
    }
    tests = {}
    pin_cpu = context["allowed_cpus"][index % len(context["allowed_cpus"])]
    for name in COMPARE_TESTS:
        target, arguments = test_specs[name]
        executable = executable_path(build, target)
        argv = [str(executable)] + arguments
        if context["taskset"]:
            argv = [context["taskset"], "-c", str(pin_cpu)] + argv
        command = run_command(
            "test_" + name, argv, context["source"], result_dir,
            context["timeout"], environment, hash_output=True
        )
        command["executable_sha256"] = digest_bytes(executable.read_bytes())
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
        if command["returncode"] != 0:
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": "optional-CUDA test failed",
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

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
        "compiler": compiler,
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
    assert digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2})
    assert variant_flags("scalar") == ["-mssse3", "-mno-avx"]
    with tempfile.TemporaryDirectory(prefix="leo2-backend-self-test-") as directory:
        path = Path(directory) / "result.json"
        value = {"z": [3, 2, 1], "a": "stable"}
        atomic_write_json(path, value)
        assert json.loads(path.read_text(encoding="utf-8")) == value
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
