#!/usr/bin/env python3
"""Run resumable exact-mask versus mature-prefix encoder experiments.

The C++ cell benchmark contains both paths in one binary and measures their
setup-inclusive pair in independent ABBA rounds. This runner adds deterministic
manifests, source/object/executable/runtime identity,
resumable per-cell artifacts, CPU pinning, and fail-closed authority rules.
It never configures or builds Leopard: a pinned measurement must consume a
clean, already-built executable whose embedded Git SHA matches the source tree.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import errno
import hashlib
import json
import math
import os
import platform
import shlex
import shutil
import socket
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

try:
    import fcntl
except ImportError:  # pragma: no cover - pinned evidence is Linux-only.
    fcntl = None


SCHEMA = "leopard2-sparse-encode-crossover/v2"
JOB_SCHEMA = "leopard2-sparse-encode-crossover-job/v2"
BENCHMARK_SCHEMA = "leopard2-sparse-encode-benchmark-v4"
ATTESTATION_SCHEMA = "leopard2-benchmark-isolation-attestation/v1"
LINK_SIDECAR_SCHEMA = "leopard2-sparse-encode-link-sidecar/v1"
CELL_SCHEMA = "leopard2-sparse-encode-cells/v2"
SOURCE_FILES = (
    "CMakeLists.txt",
    "cmake/WriteSparseEncodeEvidenceSidecar.cmake",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX512.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h",
    "Leopard2Dispatch.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "bench/leopard2/sparse_encode_benchmark.cpp",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
    "tests/cmake/test_sparse_encode_benchmark_registration.cmake",
    "tools/leopard2_sparse_encode_benchmark_json_test.py",
    "tools/leopard2_sparse_encode_crossover.py",
)
KNOWN_BACKENDS = ("scalar", "ssse3", "avx2", "avx512", "auto")


class CrossoverError(RuntimeError):
    pass


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_value(value: Any) -> str:
    return digest_bytes(canonical_bytes(value))


def atomic_write(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(value)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write(
        path,
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True).encode("utf-8")
        + b"\n",
    )


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError(f"cannot read {path}: {error}") from error


def parse_csv_unsigned(text: str, name: str, maximum: int) -> list[int]:
    result: list[int] = []
    for raw in text.split(","):
        item = raw.strip()
        if not item or not item.isdigit():
            raise CrossoverError(f"{name} must contain unsigned decimal integers")
        value = int(item)
        if value <= 0 or value > maximum:
            raise CrossoverError(f"{name} values must be in 1..{maximum}")
        if value in result:
            raise CrossoverError(f"{name} contains duplicate {value}")
        result.append(value)
    if not result:
        raise CrossoverError(f"{name} is empty")
    return sorted(result)


def parse_backends(text: str) -> list[str]:
    result = [item.strip().lower() for item in text.split(",") if item.strip()]
    if not result or len(result) != len(set(result)):
        raise CrossoverError("backends must be a nonempty unique list")
    invalid = sorted(set(result) - set(KNOWN_BACKENDS))
    if invalid:
        raise CrossoverError(f"unsupported backends: {','.join(invalid)}")
    return result


def validate_resume_policy(mode: str, no_resume: bool) -> None:
    if mode == "pinned" and not no_resume:
        raise CrossoverError(
            "pinned evidence requires --no-resume so every retained timing "
            "belongs to the measured isolation interval"
        )


def allowed_cpus() -> list[int]:
    if hasattr(os, "sched_getaffinity"):
        try:
            cpus = sorted(os.sched_getaffinity(0))
            if cpus:
                return cpus
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def compact_cpu_list(cpus: list[int]) -> str:
    values = sorted(set(cpus))
    if not values:
        return ""
    ranges: list[str] = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(first) if first == previous else f"{first}-{previous}")
        first = previous = value
    ranges.append(str(first) if first == previous else f"{first}-{previous}")
    return ",".join(ranges)


def parse_cpu_list(text: str) -> list[int]:
    values: list[int] = []
    for item in text.split(","):
        bounds = item.strip().split("-")
        if len(bounds) == 1 and bounds[0].isdigit():
            first = last = int(bounds[0])
        elif len(bounds) == 2 and all(value.isdigit() for value in bounds):
            first, last = map(int, bounds)
            if last < first:
                raise CrossoverError("descending CPU topology range")
        else:
            raise CrossoverError("invalid CPU topology list")
        values.extend(range(first, last + 1))
    if not values or len(values) != len(set(values)):
        raise CrossoverError("CPU topology list is empty or duplicated")
    return sorted(values)


def exact_smt_pair(cpu: int, allowed: list[int]) -> dict[str, Any]:
    records = []
    root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    sibling_text = read_optional(root / "thread_siblings_list")
    if sibling_text is None:
        raise CrossoverError("pinned evidence requires Linux SMT topology")
    siblings = parse_cpu_list(sibling_text)
    if len(siblings) != 2 or cpu not in siblings:
        raise CrossoverError("pinned evidence requires an exact two-thread SMT pair")
    if any(value not in allowed for value in siblings):
        raise CrossoverError("the complete SMT pair is outside process affinity")
    if not set(allowed).difference(siblings):
        raise CrossoverError("no housekeeping CPU remains outside the timed SMT pair")
    for logical_cpu in siblings:
        topology = Path(f"/sys/devices/system/cpu/cpu{logical_cpu}/topology")
        text = read_optional(topology / "thread_siblings_list")
        core = read_optional(topology / "core_id")
        package = read_optional(topology / "physical_package_id")
        if text is None or core is None or package is None or parse_cpu_list(text) != siblings:
            raise CrossoverError("SMT siblings do not report mutual exact topology")
        records.append({
            "cpu": logical_cpu, "core_id": core,
            "physical_package_id": package,
            "thread_siblings_list": text,
        })
    if records[0]["core_id"] != records[1]["core_id"] or (
        records[0]["physical_package_id"] != records[1]["physical_package_id"]
    ):
        raise CrossoverError("SMT pair disagrees on core/package identity")
    return {"cpus": siblings, "records": records}


def cpu_stat(cpu: int) -> dict[str, Any]:
    prefix = f"cpu{cpu} "
    text_value = read_optional(Path("/proc/stat"))
    if text_value is None:
        raise CrossoverError("cannot read /proc/stat isolation counters")
    for line in text_value.splitlines():
        if line.startswith(prefix):
            tokens = line.split()
            if len(tokens) < 9:
                raise CrossoverError("CPU scheduler counter record is incomplete")
            try:
                fields = [int(value) for value in tokens[1:9]]
            except ValueError as error:
                raise CrossoverError("CPU scheduler counters are not integers") from error
            return {"cpu": cpu, "fields": fields}
    raise CrossoverError(f"CPU {cpu} is absent from /proc/stat")


def cpu_stat_delta(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    for label, snapshot in (("before", before), ("after", after)):
        fields = snapshot.get("fields") if isinstance(snapshot, dict) else None
        if (not isinstance(snapshot, dict) or not isinstance(snapshot.get("cpu"), int)
                or not isinstance(fields, list) or len(fields) != 8
                or any(not isinstance(value, int) or isinstance(value, bool)
                       or value < 0 for value in fields)):
            raise CrossoverError(f"{label} CPU scheduler snapshot is malformed")
    if before.get("cpu") != after.get("cpu"):
        raise CrossoverError("CPU scheduler snapshots use different CPUs")
    differences = [last - first for first, last in zip(before["fields"], after["fields"])]
    if any(value < 0 for value in differences):
        raise CrossoverError("CPU scheduler counter moved backwards")
    # /proc/stat fields 3 and 4 (zero based) are idle and iowait.
    nonidle = sum(value for index, value in enumerate(differences) if index not in (3, 4))
    return {"cpu": before["cpu"], "fields": differences,
            "nonidle_jiffies": nonidle, "total_jiffies": sum(differences)}


def sibling_delta_is_idle(delta: dict[str, Any]) -> bool:
    return (
        isinstance(delta, dict)
        and isinstance(delta.get("total_jiffies"), int)
        and delta["total_jiffies"] > 0
        and delta.get("nonidle_jiffies") == 0
    )


class CpuPairLease:
    def __init__(self, cpus: list[int], root: Path | None = None):
        if len(cpus) != 2 or cpus != sorted(set(cpus)):
            raise CrossoverError("CPU pair lease requires two distinct sorted CPUs")
        self.cpus = cpus
        self.root = root or (Path("/run/user") / str(os.getuid()))
        self.directory = self.root / "leopard2-cpu-leases"
        self.path = self.directory / (
            f"leopard2-cpu-pair-{os.getuid()}-{cpus[0]}-{cpus[1]}.lock"
        )
        self.descriptor: int | None = None
        self.kernel_socket: socket.socket | None = None

    def payload(self) -> dict[str, Any]:
        return {
            "cpus": self.cpus, "schema": "leopard2-cpu-pair-lease/v1",
            "uid": os.getuid(),
        }

    def kernel_name(self) -> bytes:
        material = canonical_bytes({
            "cpus": self.cpus, "root": os.path.abspath(self.directory),
            "schema": "leopard2-cpu-pair-lease/v1", "uid": os.getuid(),
        })
        return b"\0leopard2-pair-v1-" + hashlib.sha256(material).hexdigest()[
            :40
        ].encode("ascii")

    def __enter__(self) -> "CpuPairLease":
        if fcntl is None or not sys.platform.startswith("linux"):
            raise CrossoverError("protected CPU pair leases require Linux fcntl")
        self.kernel_socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.kernel_socket.set_inheritable(False)
        try:
            self.kernel_socket.bind(self.kernel_name())
        except OSError as error:
            self.close()
            if error.errno == errno.EADDRINUSE:
                raise CrossoverError("physical CPU pair already has a kernel lease") from error
            raise CrossoverError(f"cannot bind CPU pair kernel lease: {error}") from error
        try:
            root_stat = self.root.stat()
        except OSError as error:
            self.close()
            raise CrossoverError(f"cannot inspect per-user runtime directory: {error}") from error
        if (not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_uid != os.getuid()
                or stat.S_IMODE(root_stat.st_mode) != 0o700):
            self.close()
            raise CrossoverError("per-user runtime directory is not private and owned")
        self.directory.mkdir(mode=0o700, exist_ok=True)
        if self.directory.is_symlink():
            self.close()
            raise CrossoverError("CPU lease directory must not be a symlink")
        directory_stat = self.directory.stat()
        if (directory_stat.st_uid != os.getuid()
                or stat.S_IMODE(directory_stat.st_mode) != 0o700):
            self.close()
            raise CrossoverError("CPU lease directory is not private and owned")
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        try:
            self.descriptor = os.open(self.path, flags, 0o600)
        except OSError as error:
            self.close()
            raise CrossoverError(f"cannot open CPU pair lease: {error}") from error
        os.fchmod(self.descriptor, 0o600)
        file_stat = os.fstat(self.descriptor)
        if (not stat.S_ISREG(file_stat.st_mode) or file_stat.st_uid != os.getuid()
                or file_stat.st_nlink != 1):
            self.close()
            raise CrossoverError("CPU pair lease file has unsafe identity")
        try:
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            self.close()
            raise CrossoverError("physical CPU pair is already leased") from error
        except OSError as error:
            self.close()
            raise CrossoverError(f"cannot lock physical CPU pair: {error}") from error
        payload = canonical_bytes(self.payload())
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        retained = os.read(self.descriptor, max(4096, len(payload) + 1))
        if not retained:
            if os.write(self.descriptor, payload) != len(payload):
                self.close()
                raise CrossoverError("short write to CPU pair lease")
            os.fsync(self.descriptor)
        elif retained != payload:
            self.close()
            raise CrossoverError("CPU pair lease has noncanonical contents")
        return self

    def identity(self) -> dict[str, Any]:
        if (self.descriptor is None or self.kernel_socket is None
                or self.kernel_socket.fileno() < 0
                or self.kernel_socket.getsockname() != self.kernel_name()):
            raise CrossoverError("CPU pair lease is not held")
        file_stat = os.fstat(self.descriptor)
        path_stat = self.path.stat()
        if ((file_stat.st_dev, file_stat.st_ino) !=
                (path_stat.st_dev, path_stat.st_ino)):
            raise CrossoverError("CPU pair lease path identity changed")
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        if os.read(self.descriptor, 4096) != canonical_bytes(self.payload()):
            raise CrossoverError("CPU pair lease payload changed")
        return {
            "cpus": self.cpus, "device": file_stat.st_dev,
            "inode": file_stat.st_ino, "path": str(self.path.resolve()),
            "mechanism": "abstract-af-unix-bind-plus-fcntl-flock",
            "payload": self.payload(),
            "kernel_name_sha256": hashlib.sha256(self.kernel_name()).hexdigest(),
        }

    def close(self) -> None:
        if self.descriptor is not None:
            descriptor, self.descriptor = self.descriptor, None
            if fcntl is not None:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            os.close(descriptor)
        if self.kernel_socket is not None:
            kernel_socket, self.kernel_socket = self.kernel_socket, None
            kernel_socket.close()

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self.close()


def read_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def machine_identity(cpus: list[int]) -> dict[str, Any]:
    uname = platform.uname()
    model = None
    cpuinfo = read_optional(Path("/proc/cpuinfo"))
    if cpuinfo:
        for line in cpuinfo.splitlines():
            if line.lower().startswith("model name") and ":" in line:
                model = line.split(":", 1)[1].strip()
                break
    return {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_model": model,
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "uname": {
            "machine": uname.machine,
            "node": uname.node,
            "release": uname.release,
            "system": uname.system,
            "version": uname.version,
        },
    }


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    stat_result = resolved.stat()
    if not resolved.is_file():
        raise CrossoverError(f"build input is not a regular file: {resolved}")
    return {
        "path": str(resolved),
        "sha256": digest_bytes(resolved.read_bytes()),
        "size": stat_result.st_size,
        "mtime_ns": stat_result.st_mtime_ns,
    }


def parse_link_sidecar(path: Path) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeError) as error:
        raise CrossoverError(f"cannot read build-time link sidecar: {error}") from error
    values: dict[str, str] = {}
    for line in lines:
        if line.count("=") != 1:
            raise CrossoverError("build-time link sidecar is malformed")
        key, value = line.split("=", 1)
        if not key or not value or key in values:
            raise CrossoverError("build-time link sidecar has invalid fields")
        values[key] = value
    expected = {
        "schema", "executable_sha256", "executable_size",
        "production_archive_sha256", "production_archive_size",
        "benchmark_object_sha256", "benchmark_object_size",
        "oracle_object_sha256", "oracle_object_size", "link_recipe_kind",
        "benchmark_link_recipe_sha256", "benchmark_link_recipe_size",
        "production_link_recipe_sha256", "production_link_recipe_size",
    }
    if set(values) != expected or values.get("schema") != LINK_SIDECAR_SCHEMA:
        raise CrossoverError("build-time link sidecar has an unknown schema")
    for prefix in (
        "executable", "production_archive", "benchmark_object", "oracle_object"
    ):
        digest = values[f"{prefix}_sha256"]
        size = values[f"{prefix}_size"]
        if (len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
                or not size.isdigit()):
            raise CrossoverError("build-time link sidecar has invalid artifact identity")
    if values["link_recipe_kind"] not in (
        "cmake-link-txt-v1", "ninja-tool-commands-v1", "missing"
    ):
        raise CrossoverError("build-time link sidecar has invalid recipe kind")
    for prefix in ("benchmark_link_recipe", "production_link_recipe"):
        digest = values[f"{prefix}_sha256"]
        size = values[f"{prefix}_size"]
        if (digest, size) == ("missing", "missing"):
            continue
        if (len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
                or not size.isdigit()):
            raise CrossoverError("build-time link sidecar has invalid recipe identity")
    return values


def recipe_file_tokens(recipe: bytes, build_root: Path) -> list[Path]:
    try:
        tokens = shlex.split(recipe.decode("utf-8"), posix=True)
    except (UnicodeError, ValueError) as error:
        raise CrossoverError(f"cannot parse CMake link recipe: {error}") from error
    result: list[Path] = []
    for token in tokens:
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = build_root / candidate
        try:
            resolved = candidate.resolve(strict=True)
        except OSError:
            continue
        if resolved.is_file():
            result.append(resolved)
    return result


def capture_link_recipes(
    fields: dict[str, str], build_root: Path,
    benchmark_recipe: Path, production_recipes: list[Path],
) -> tuple[str, bytes, bytes, dict[str, Any]]:
    if benchmark_recipe.is_file() and len(production_recipes) == 1:
        production_recipe = production_recipes[0]
        benchmark_bytes = benchmark_recipe.read_bytes()
        production_bytes = production_recipe.read_bytes()
        return (
            "cmake-link-txt-v1", benchmark_bytes, production_bytes,
            {
                "benchmark": file_identity(benchmark_recipe),
                "production_archive": file_identity(production_recipe),
            },
        )
    generator = fields.get("CMAKE_GENERATOR", "")
    make_program = fields.get("CMAKE_MAKE_PROGRAM", "")
    if "Ninja" not in generator or not make_program:
        return "missing", b"", b"", {}
    program = Path(make_program).resolve(strict=True)

    def commands(target: str) -> bytes:
        completed = subprocess.run(
            [str(program), "-C", str(build_root), "-t", "commands", target],
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
        )
        if completed.returncode != 0 or not completed.stdout:
            raise CrossoverError(
                f"cannot capture Ninja command graph for {target}: "
                + completed.stderr.decode("utf-8", errors="replace").strip()
            )
        return completed.stdout

    return (
        "ninja-tool-commands-v1",
        commands("bench_leopard2_sparse_encode"),
        commands("leopard"),
        {"build_program": file_identity(program)},
    )


def build_metadata(executable: Path, source: Path) -> dict[str, Any]:
    executable = executable.resolve(strict=True)
    cache = None
    for directory in (executable.parent, *list(executable.parents)[1:4]):
        candidate = directory / "CMakeCache.txt"
        if candidate.is_file():
            cache = candidate
            break
    if cache is None:
        raise CrossoverError("CMakeCache.txt was not found near the executable")
    data = cache.read_bytes()
    fields: dict[str, str] = {}
    for raw in data.decode("utf-8", errors="replace").splitlines():
        if raw.startswith("//") or raw.startswith("#") or "=" not in raw:
            continue
        key_type, value = raw.split("=", 1)
        key = key_type.split(":", 1)[0]
        if key in (
            "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
            "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_HOME_DIRECTORY",
            "CMAKE_GENERATOR", "CMAKE_MAKE_PROGRAM",
            "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS",
            "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
        ):
            fields[key] = value
    home = fields.get("CMAKE_HOME_DIRECTORY")
    if not home or Path(home).resolve() != source.resolve():
        raise CrossoverError(
            "benchmark build cache belongs to a different source root"
        )
    build_root = cache.parent.resolve()
    try:
        executable.relative_to(build_root)
    except ValueError as error:
        raise CrossoverError("benchmark executable is outside its build root") from error
    benchmark_target_dir = (
        build_root / "CMakeFiles" / "leopard2_sparse_encode_benchmark_object.dir"
    )
    oracle_target_dir = (
        build_root / "CMakeFiles" / "leopard2_sparse_encode_oracle_object.dir"
    )
    objects = sorted(
        list(benchmark_target_dir.rglob("sparse_encode_benchmark.cpp.o"))
        + list(benchmark_target_dir.rglob("sparse_encode_benchmark.cpp.obj"))
    )
    if len(objects) != 1:
        raise CrossoverError(
            "expected exactly one sparse benchmark object in the current build root"
        )
    oracle_objects = sorted(
        list(oracle_target_dir.rglob("direct_oracle.cpp.o"))
        + list(oracle_target_dir.rglob("direct_oracle.cpp.obj"))
    )
    if len(oracle_objects) != 1:
        raise CrossoverError(
            "expected exactly one direct-oracle object in the current build root"
        )
    direct_archives = [
        build_root / "libleopard.a", build_root / "leopard.lib",
        build_root / "Release" / "leopard.lib",
    ]
    archives = [path for path in direct_archives if path.is_file()]
    if len(archives) != 1:
        raise CrossoverError(
            "expected exactly one production leopard archive in the current build root"
        )
    archive = archives[0].resolve(strict=True)
    benchmark_recipe = (
        build_root / "CMakeFiles" / "bench_leopard2_sparse_encode.dir" / "link.txt"
    )
    production_recipe_candidates = [
        build_root / "CMakeFiles" / "leopard.dir" / "link.txt",
        build_root / "CMakeFiles" / "libleopard.dir" / "link.txt",
    ]
    production_recipes = [path for path in production_recipe_candidates if path.is_file()]
    sidecar_path = Path(str(executable) + ".leopard2-evidence")
    if not sidecar_path.is_file():
        raise CrossoverError("build-time sparse evidence sidecar is missing")
    sidecar = parse_link_sidecar(sidecar_path)
    executable_identity = file_identity(executable)
    archive_identity = file_identity(archive)
    benchmark_objects = [file_identity(objects[0]), file_identity(oracle_objects[0])]
    if (sidecar["executable_sha256"] != executable_identity["sha256"]
            or int(sidecar["executable_size"]) != executable_identity["size"]
            or sidecar["production_archive_sha256"] != archive_identity["sha256"]
            or int(sidecar["production_archive_size"]) != archive_identity["size"]
            or sidecar["benchmark_object_sha256"] != benchmark_objects[0]["sha256"]
            or int(sidecar["benchmark_object_size"]) != benchmark_objects[0]["size"]
            or sidecar["oracle_object_sha256"] != benchmark_objects[1]["sha256"]
            or int(sidecar["oracle_object_size"]) != benchmark_objects[1]["size"]):
        raise CrossoverError("build artifacts differ from their post-link sidecar")

    recipe_kind, benchmark_recipe_bytes, production_recipe_bytes, \
        recipe_identities = capture_link_recipes(
            fields, build_root, benchmark_recipe, production_recipes
        )
    recipes_attested = recipe_kind != "missing"
    production_objects: list[dict[str, Any]] = []
    if recipes_attested:
        if (sidecar["link_recipe_kind"] != recipe_kind
                or sidecar["benchmark_link_recipe_sha256"] !=
                    digest_bytes(benchmark_recipe_bytes)
                or int(sidecar["benchmark_link_recipe_size"]) !=
                    len(benchmark_recipe_bytes)
                or sidecar["production_link_recipe_sha256"] !=
                    digest_bytes(production_recipe_bytes)
                or int(sidecar["production_link_recipe_size"]) !=
                    len(production_recipe_bytes)):
            raise CrossoverError("link graph differs from its post-link sidecar")

        def command_segments(recipe: bytes) -> list[list[Path]]:
            segments: list[list[Path]] = []
            for line in recipe.splitlines():
                for segment in line.split(b" && "):
                    tokens = recipe_file_tokens(segment, build_root)
                    if tokens:
                        segments.append(tokens)
            return segments

        required_benchmark_inputs = [objects[0].resolve(), oracle_objects[0].resolve(), archive]
        benchmark_matches = [
            tokens for tokens in command_segments(benchmark_recipe_bytes)
            if executable in tokens
            and all(tokens.count(path) == 1 for path in required_benchmark_inputs)
        ]
        if len(benchmark_matches) != 1:
            raise CrossoverError(
                "benchmark link recipe does not bind its two objects and production archive"
            )
        production_matches = []
        for tokens in command_segments(production_recipe_bytes):
            object_paths = sorted(set(
                path for path in tokens if path.suffix in (".o", ".obj")
            ))
            if tokens.count(archive) == 1 and object_paths:
                production_matches.append((tokens, object_paths))
        if len(production_matches) != 1:
            raise CrossoverError(
                "production archive recipe does not bind its object graph and output"
            )
        object_paths = production_matches[0][1]
        production_objects = [file_identity(path) for path in object_paths]
        if archive_identity["mtime_ns"] < max(
            identity["mtime_ns"] for identity in production_objects
        ):
            raise CrossoverError("production archive is older than a recipe-bound object")
    elif sidecar["link_recipe_kind"] != "missing" or any(
        sidecar[name] != "missing" for name in (
        "benchmark_link_recipe_sha256", "benchmark_link_recipe_size",
        "production_link_recipe_sha256", "production_link_recipe_size",
    )):
        raise CrossoverError("post-link sidecar names a link recipe that is now missing")

    newest_direct_input = max(
        [archive_identity["mtime_ns"]]
        + [identity["mtime_ns"] for identity in benchmark_objects]
    )
    if executable_identity["mtime_ns"] < newest_direct_input:
        raise CrossoverError("benchmark executable is older than a direct link input")
    return {
        "build_root": str(build_root),
        "cmake_cache": {
            "path": str(cache.resolve()), "sha256": digest_bytes(data),
            "selected_fields": fields,
        },
        "executable": executable_identity,
        "link_graph_attested": recipes_attested,
        "link_recipes": recipe_identities,
        "post_link_sidecar": {
            "file": file_identity(sidecar_path),
            "values": sidecar,
        },
        "linked_inputs": {
            "benchmark_object": benchmark_objects[0],
            "direct_oracle_object": benchmark_objects[1],
            "production_archive": archive_identity,
            "production_objects": production_objects,
        },
    }


def source_fingerprint(source: Path) -> dict[str, Any]:
    files: dict[str, str] = {}
    missing: list[str] = []
    for relative in SOURCE_FILES:
        path = source / relative
        if path.is_file():
            files[relative] = digest_bytes(path.read_bytes())
        else:
            missing.append(relative)
    if missing:
        raise CrossoverError(f"required source files are missing: {','.join(missing)}")
    return {"digest": digest_value(files), "files": files}


def git_state(source: Path) -> dict[str, Any]:
    def git(*arguments: str) -> str:
        completed = subprocess.run(
            ["git", *arguments], cwd=source, check=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30,
        )
        if completed.returncode != 0:
            raise CrossoverError(
                f"git {' '.join(arguments)} failed: {completed.stderr.strip()}"
            )
        return completed.stdout.strip()

    sha = git("rev-parse", "HEAD")
    dirty_text = git("status", "--porcelain", "--untracked-files=normal")
    if len(sha) != 40 or any(character not in "0123456789abcdef" for character in sha):
        raise CrossoverError("Git HEAD is not a full lowercase SHA-1")
    return {"dirty": bool(dirty_text), "sha": sha}


def base_cells() -> list[dict[str, Any]]:
    return [
        {"profile": "high", "field": "gf8", "K": 192, "R": 64,
         "cell_kind": "sparse_candidate", "mask_name": "edge_sparse",
         "requested_parity": [0, 31, 63]},
        {"profile": "high", "field": "gf8", "K": 192, "R": 64,
         "cell_kind": "sparse_candidate", "mask_name": "scattered_sparse",
         "requested_parity": [3, 15, 39, 63]},
        {"profile": "high", "field": "gf8", "K": 192, "R": 64,
         "cell_kind": "prefix_neighbor", "mask_name": "dense_prefix_neighbor",
         "requested_parity": list(range(64))},
        {"profile": "low", "field": "gf8", "K": 32, "R": 192,
         "cell_kind": "sparse_candidate", "mask_name": "edge_sparse",
         "requested_parity": [0, 31, 32, 95, 191]},
        {"profile": "low", "field": "gf8", "K": 32, "R": 192,
         "cell_kind": "sparse_candidate", "mask_name": "scattered_sparse",
         "requested_parity": [3, 23, 47, 111, 191]},
        {"profile": "low", "field": "gf8", "K": 32, "R": 192,
         "cell_kind": "prefix_neighbor", "mask_name": "dense_prefix_neighbor",
         "requested_parity": list(range(32))},
        {"profile": "high", "field": "gf16", "K": 1000, "R": 200,
         "cell_kind": "sparse_candidate", "mask_name": "edge_sparse",
         "requested_parity": [0, 63, 127, 199]},
        {"profile": "high", "field": "gf16", "K": 1000, "R": 200,
         "cell_kind": "sparse_candidate", "mask_name": "scattered_sparse",
         "requested_parity": [3, 31, 95, 159, 199]},
        {"profile": "high", "field": "gf16", "K": 1000, "R": 200,
         "cell_kind": "prefix_neighbor", "mask_name": "dense_prefix_neighbor",
         "requested_parity": list(range(200))},
        {"profile": "low", "field": "gf16", "K": 128, "R": 896,
         "cell_kind": "sparse_candidate", "mask_name": "edge_sparse",
         "requested_parity": [0, 127, 128, 383, 895]},
        {"profile": "low", "field": "gf16", "K": 128, "R": 896,
         "cell_kind": "sparse_candidate", "mask_name": "scattered_sparse",
         "requested_parity": [7, 63, 135, 255, 519, 895]},
        {"profile": "low", "field": "gf16", "K": 128, "R": 896,
         "cell_kind": "prefix_neighbor", "mask_name": "dense_prefix_neighbor",
         "requested_parity": list(range(128))},
    ]


def make_cells(backends: list[str], shard_bytes: list[int]) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for base in base_cells():
        for backend in sorted(backends):
            for byte_count in sorted(shard_bytes):
                cell = dict(base)
                cell.update(backend=backend, shard_bytes=byte_count)
                cells.append(cell)
    return sorted(
        cells,
        key=lambda cell: (
            cell["field"], cell["profile"], cell["K"], cell["R"],
            cell["mask_name"], cell["backend"], cell["shard_bytes"],
        ),
    )


def validate_cells(value: Any) -> list[dict[str, Any]]:
    """Validate a complete, explicitly targeted benchmark cell manifest."""
    if not isinstance(value, dict) or set(value) != {"schema", "cells"} or (
        value.get("schema") != CELL_SCHEMA
    ):
        raise CrossoverError("cell manifest has an unknown schema or fields")
    raw_cells = value.get("cells")
    if not isinstance(raw_cells, list) or not raw_cells:
        raise CrossoverError("cell manifest contains no cells")

    cells: list[dict[str, Any]] = []
    expected_keys = {
        "profile", "field", "K", "R", "cell_kind", "mask_name",
        "requested_parity", "backend", "shard_bytes",
    }
    for index, raw in enumerate(raw_cells):
        if not isinstance(raw, dict) or set(raw) != expected_keys:
            raise CrossoverError(
                f"cell {index} has unknown or missing fields"
            )
        profile = raw["profile"]
        field = raw["field"]
        backend = raw["backend"]
        cell_kind = raw["cell_kind"]
        k = raw["K"]
        r = raw["R"]
        byte_count = raw["shard_bytes"]
        if profile not in ("high", "low") or field not in ("gf8", "gf16"):
            raise CrossoverError(f"cell {index} has an invalid profile or field")
        if backend not in KNOWN_BACKENDS:
            raise CrossoverError(f"cell {index} has an invalid backend")
        if cell_kind not in ("sparse_candidate", "prefix_neighbor"):
            raise CrossoverError(f"cell {index} has an invalid inference role")
        if not isinstance(k, int) or isinstance(k, bool) or k <= 0 or (
            not isinstance(r, int) or isinstance(r, bool) or r <= 0
        ):
            raise CrossoverError(f"cell {index} has invalid K or R")
        if not isinstance(byte_count, int) or isinstance(byte_count, bool) or (
            byte_count <= 0 or byte_count > 1 << 30
        ):
            raise CrossoverError(f"cell {index} has invalid shard bytes")
        if field == "gf16" and byte_count % 2 != 0:
            raise CrossoverError(f"cell {index} has an odd GF16 byte count")
        if not isinstance(raw["mask_name"], str) or not raw["mask_name"].strip():
            raise CrossoverError(f"cell {index} has no mask name")

        requested_raw = raw["requested_parity"]
        if requested_raw == "all":
            requested = list(range(r))
        elif isinstance(requested_raw, dict) and set(requested_raw) == {"prefix"}:
            prefix_count = requested_raw["prefix"]
            if not isinstance(prefix_count, int) or isinstance(prefix_count, bool) or (
                prefix_count <= 0 or prefix_count > r
            ):
                raise CrossoverError(
                    f"cell {index} has an invalid requested parity prefix"
                )
            requested = list(range(prefix_count))
        else:
            requested = requested_raw
        if not isinstance(requested, list) or not requested or any(
            not isinstance(item, int) or isinstance(item, bool)
            or item < 0 or item >= r for item in requested
        ) or requested != sorted(set(requested)):
            raise CrossoverError(f"cell {index} has invalid requested parity")
        if cell_kind == "prefix_neighbor" and requested != list(
            range(requested[-1] + 1)
        ):
            raise CrossoverError(
                f"cell {index} labels a non-prefix mask as a prefix neighbor"
            )
        if cell_kind == "sparse_candidate" and requested == list(
            range(requested[-1] + 1)
        ):
            raise CrossoverError(
                f"cell {index} labels a prefix mask as a sparse candidate"
            )

        padded_side_value = r if profile == "high" else k
        padded_side = 1
        while padded_side < padded_side_value:
            padded_side <<= 1
        field_order = 256 if field == "gf8" else 65536
        parent_span = k + padded_side if profile == "high" else padded_side + r
        if parent_span > field_order:
            raise CrossoverError(f"cell {index} does not fit its field")

        cell = dict(raw)
        cell["requested_parity"] = requested
        cells.append(cell)

    cells.sort(key=lambda cell: (
        cell["field"], cell["profile"], cell["K"], cell["R"],
        cell["mask_name"], cell["backend"], cell["shard_bytes"],
    ))
    identities = [canonical_bytes(cell) for cell in cells]
    if len(identities) != len(set(identities)):
        raise CrossoverError("cell manifest contains duplicate cells")

    groups: dict[tuple[Any, ...], set[str]] = {}
    seen_masks: set[tuple[Any, ...]] = set()
    seen_requested: set[tuple[Any, ...]] = set()
    for cell in cells:
        group = (
            cell["profile"], cell["field"], cell["K"], cell["R"],
            cell["backend"], cell["shard_bytes"],
        )
        groups.setdefault(group, set()).add(cell["cell_kind"])
        mask_identity = group + (cell["mask_name"],)
        if mask_identity in seen_masks:
            raise CrossoverError(
                "cell manifest reuses a mask name within one inference group"
            )
        seen_masks.add(mask_identity)
        requested_identity = group + (tuple(cell["requested_parity"]),)
        if requested_identity in seen_requested:
            raise CrossoverError(
                "cell manifest duplicates a requested mask within one inference group"
            )
        seen_requested.add(requested_identity)
    required_roles = {"sparse_candidate", "prefix_neighbor"}
    if any(roles != required_roles for roles in groups.values()):
        raise CrossoverError(
            "every cell-manifest inference group requires sparse and prefix roles"
        )
    return cells


def load_cell_manifest(path: Path) -> tuple[list[dict[str, Any]], dict[str, str]]:
    try:
        payload = path.read_bytes()
        value = json.loads(payload.decode("utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError(f"cannot read {path}: {error}") from error
    cells = validate_cells(value)
    return cells, {"path": str(path), "sha256": digest_bytes(payload)}


def load_cells(path: Path) -> list[dict[str, Any]]:
    return load_cell_manifest(path)[0]


def validate_cell_manifest_identity(manifest: dict[str, Any]) -> None:
    identity = manifest.get("settings", {}).get("cell_manifest")
    if identity is None:
        return
    if not isinstance(identity, dict) or set(identity) != {"path", "sha256"}:
        raise CrossoverError("cell manifest identity is malformed")
    path_text = identity.get("path")
    expected_sha = identity.get("sha256")
    if not isinstance(path_text, str) or not Path(path_text).is_absolute() or (
        not isinstance(expected_sha, str) or len(expected_sha) != 64
    ):
        raise CrossoverError("cell manifest identity is malformed")
    path = Path(path_text)
    cells, current = load_cell_manifest(path)
    if current != identity:
        raise CrossoverError("cell manifest changed before an invocation")
    # Hashing alone must not let a later parser/schema change reinterpret cells.
    if canonical_bytes(cells) != canonical_bytes(
        [job.get("cell") for job in manifest.get("jobs", [])]
    ):
        raise CrossoverError("cell manifest cells differ from the run manifest")


def mask_text(indices: list[int]) -> str:
    return ",".join(str(value) for value in indices)


def benchmark_command(
    executable: Path, cell: dict[str, Any], settings: dict[str, Any]
) -> list[str]:
    return [
        str(executable),
        "--profile", cell["profile"],
        "--field", cell["field"],
        "--k", str(cell["K"]),
        "--r", str(cell["R"]),
        "--bytes", str(cell["shard_bytes"]),
        "--requested-parity", mask_text(cell["requested_parity"]),
        "--backend", cell["backend"],
        "--iterations", str(settings["iterations"]),
        "--rounds", str(settings["rounds"]),
        "--warmups", str(settings["warmups"]),
        "--setup-iterations", str(settings["setup_iterations"]),
        "--reuse", mask_text(settings["reuse"]),
        "--memory-mib", str(settings["memory_mib"]),
        "--seed", str(settings["seed"]),
    ]


def validate_metric(
    metric: Any, samples: int, name: str, require_positive: bool = True
) -> None:
    if not isinstance(metric, dict) or not isinstance(metric.get("samples"), list):
        raise CrossoverError(f"missing metric {name}")
    values = metric["samples"]
    if len(values) != samples or not all(
        isinstance(value, (int, float)) and math.isfinite(value)
        and (value > 0 if require_positive else True) for value in values
    ):
        raise CrossoverError(f"invalid samples for {name}")
    median = statistics.median(values)
    if abs(float(metric.get("median", -1)) - median) > 0.002:
        raise CrossoverError(f"wrong median for {name}")
    if abs(float(metric.get("minimum", -1)) - min(values)) > 0.002 or (
        abs(float(metric.get("maximum", -1)) - max(values)) > 0.002
    ):
        raise CrossoverError(f"wrong extrema for {name}")
    expected_mad = statistics.median(
        abs(float(value) - median) for value in values
    )
    if not isinstance(metric.get("mad"), (int, float)) or (
        metric["mad"] < 0 or abs(float(metric["mad"]) - expected_mad) > 0.002
    ):
        raise CrossoverError(f"wrong MAD for {name}")


def validate_benchmark_result(
    result: Any,
    cell: dict[str, Any],
    settings: dict[str, Any],
    expected_git: dict[str, Any],
    expected_executable: Path | None = None,
) -> None:
    if not isinstance(result, dict) or result.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    if result.get("authoritative") is not False:
        raise CrossoverError("raw benchmark improperly claimed authority")
    build = result.get("build")
    if not isinstance(build, dict):
        raise CrossoverError("benchmark omitted build identity")
    if build.get("source_git_sha") != expected_git["sha"]:
        raise CrossoverError("benchmark binary source SHA differs from runner source")
    if build.get("source_dirty") != int(expected_git["dirty"]):
        raise CrossoverError("benchmark binary dirty marker differs from runner source")
    if build.get("library_test_hooks") is not False:
        raise CrossoverError(
            "benchmark library contains test-hook instrumentation; "
            "link bench_leopard2_sparse_encode against production leopard"
        )
    if not isinstance(build.get("compiler"), str) or not build["compiler"]:
        raise CrossoverError("benchmark binary omitted compiler identity")
    if not isinstance(build.get("compiler_version"), str) or not build["compiler_version"]:
        raise CrossoverError("benchmark binary omitted compiler version")
    if not isinstance(build.get("cplusplus"), int):
        raise CrossoverError("benchmark binary omitted C++ language identity")
    runtime = result.get("runtime")
    if not isinstance(runtime, dict):
        raise CrossoverError("benchmark omitted runtime identity")
    runtime_path = runtime.get("executable_path")
    runtime_cpus = runtime.get("allowed_cpus")
    linux_runtime = sys.platform.startswith("linux")
    if runtime.get("linux_procfs_affinity_attested") is not linux_runtime:
        raise CrossoverError("benchmark runtime-attestation platform differs")
    if linux_runtime and (
        not isinstance(runtime_path, str) or not Path(runtime_path).is_absolute()
    ):
        raise CrossoverError("benchmark emitted an invalid runtime executable path")
    if linux_runtime and expected_executable is not None:
        try:
            same_executable = Path(runtime_path).samefile(expected_executable)
        except OSError as error:
            raise CrossoverError(f"cannot compare runtime executable: {error}") from error
        if not same_executable:
            raise CrossoverError("runtime executable differs from the manifest")
    if not linux_runtime and (runtime_path != "" or runtime_cpus != []):
        raise CrossoverError("non-Linux runtime identity did not fail closed")
    if runtime_cpus != settings.get("runtime_allowed_cpus"):
        raise CrossoverError("benchmark runtime affinity differs from the manifest")
    parameters = result.get("parameters")
    expected_parameters = {
        "profile": cell["profile"], "field": cell["field"],
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["shard_bytes"],
        "requested_parity": cell["requested_parity"],
        "requested_backend": cell["backend"],
        "iterations": settings["iterations"], "rounds": settings["rounds"],
        "warmups": settings["warmups"],
        "setup_iterations": settings["setup_iterations"],
        "reuse": settings["reuse"], "memory_mib": settings["memory_mib"],
        "seed": settings["seed"],
    }
    if not isinstance(parameters, dict):
        raise CrossoverError("benchmark omitted parameters")
    for name, value in expected_parameters.items():
        if parameters.get(name) != value:
            raise CrossoverError(f"benchmark parameter differs: {name}")
    resolved = result.get("resolved")
    if not isinstance(resolved, dict) or not isinstance(resolved.get("backend"), str):
        raise CrossoverError("benchmark omitted resolved backend")
    if cell["backend"] != "auto" and resolved["backend"] != cell["backend"]:
        raise CrossoverError("benchmark resolved a different backend")
    plan = result.get("plan")
    if not isinstance(plan, dict):
        raise CrossoverError("benchmark omitted plan accounting")
    if not (
        isinstance(plan.get("retained_butterflies"), int)
        and isinstance(plan.get("prefix_butterflies"), int)
        and isinstance(plan.get("full_butterflies"), int)
        and 0 < plan["retained_butterflies"] <= plan["full_butterflies"]
        and 0 < plan["prefix_butterflies"] <= plan["full_butterflies"]
    ):
        raise CrossoverError("benchmark emitted invalid transform accounting")
    if cell.get("cell_kind") == "sparse_candidate" and not (
        plan["retained_butterflies"] < plan["prefix_butterflies"]
    ):
        raise CrossoverError("sparse candidate did not reduce prefix work")
    if cell.get("cell_kind") not in ("sparse_candidate", "prefix_neighbor"):
        raise CrossoverError("cell has an unknown inference role")
    correctness = result.get("correctness")
    if not isinstance(correctness, dict) or (
        correctness.get("exact_prefix_parity_match") is not True
    ):
        raise CrossoverError("exact parity differs from mature prefix parity")
    for name in (
        "direct_generator_parity_match", "encode_decode_recovery_match",
        "input_immutable", "allocation_guards_match", "post_timing_recheck_match",
    ):
        if correctness.get(name) is not True:
            raise CrossoverError(f"benchmark correctness closure failed: {name}")
    expected_symbols = len(cell["requested_parity"]) * min(
        8, cell["shard_bytes"] if cell["field"] == "gf8"
        else cell["shard_bytes"] // 2
    )
    if correctness.get("direct_generator_symbols_checked") != expected_symbols:
        raise CrossoverError("benchmark checked the wrong direct-generator sample set")
    canary = correctness.get("unrequested_output_canary")
    if not isinstance(canary, dict) or canary.get("match") is not True or (
        canary.get("applicable") is not (cell["profile"] == "low")
    ):
        raise CrossoverError("benchmark output-canary evidence is invalid")
    for name in (
        "direct_generator_sample_sha256", "input_sha256", "parity_sha256",
        "recovery_source_sha256", "recovered_original_sha256",
    ):
        value = correctness.get(name)
        if not isinstance(value, str) or len(value) != 64 or any(
            character not in "0123456789abcdef" for character in value
        ):
            raise CrossoverError(f"benchmark emitted an invalid {name}")
    if correctness["recovery_source_sha256"] != correctness["recovered_original_sha256"]:
        raise CrossoverError("recovery SHA-256 differs from the missing source shard")
    digest = correctness.get("digest_fnv1a64")
    if not isinstance(digest, str) or len(digest) != 18 or not digest.startswith("0x"):
        raise CrossoverError("benchmark emitted an invalid parity digest")
    metrics = result.get("metrics")
    if not isinstance(metrics, dict):
        raise CrossoverError("benchmark omitted metrics")
    primary = result.get("primary_abba")
    if not isinstance(primary, dict) or (
        primary.get("design") != "pairwise_prefix_call_local_ABBA"
        or primary.get("order_policy") != "alternating_ABBA_BAAB"
    ):
        raise CrossoverError("benchmark primary comparison is not pairwise ABBA")
    rounds = primary.get("rounds")
    if not isinstance(rounds, list) or len(rounds) != settings["rounds"]:
        raise CrossoverError("benchmark emitted the wrong ABBA round count")
    round_gains = []
    for index, round_value in enumerate(rounds):
        if not isinstance(round_value, dict) or round_value.get("round") != index:
            raise CrossoverError("benchmark emitted malformed ABBA round identity")
        observations = round_value.get("observations_ns")
        order = round_value.get("order")
        expected_order = (["prefix", "call_local", "call_local", "prefix"]
                          if index % 2 == 0 else
                          ["call_local", "prefix", "prefix", "call_local"])
        if order != expected_order or not isinstance(observations, list) or len(observations) != 4 or not all(
            isinstance(value, (int, float)) and value > 0 for value in observations
        ):
            raise CrossoverError("benchmark emitted malformed ABBA observations")
        prefix_values = [float(value) for value, form in zip(observations, order)
                         if form == "prefix"]
        call_values = [float(value) for value, form in zip(observations, order)
                       if form == "call_local"]
        prefix_round = statistics.median(prefix_values)
        call_round = statistics.median(call_values)
        log_contrast = (statistics.mean(map(math.log, prefix_values)) -
                        statistics.mean(map(math.log, call_values)))
        gain = math.expm1(log_contrast) * 100.0
        if abs(float(round_value.get("prefix_median_ns", -1)) - prefix_round) > 0.002:
            raise CrossoverError("benchmark emitted wrong round prefix median")
        if abs(float(round_value.get("call_local_median_ns", -1)) - call_round) > 0.002:
            raise CrossoverError("benchmark emitted wrong round call-local median")
        if abs(float(round_value.get("log_contrast", -999)) - log_contrast) > 0.000002:
            raise CrossoverError("benchmark emitted wrong paired log contrast")
        if abs(float(round_value.get("candidate_gain_percent", -1)) - gain) > 0.002:
            raise CrossoverError("benchmark emitted wrong round gain")
        round_gains.append(gain)
    validate_metric(primary.get("round_gain_percent"), settings["rounds"],
                    "round_gain_percent", False)
    reported_round_gains = primary["round_gain_percent"]["samples"]
    if any(abs(float(reported) - expected) > 0.002
           for reported, expected in zip(reported_round_gains, round_gains)):
        raise CrossoverError("benchmark round-gain summary differs from ABBA rounds")
    names = (
        "schedule_setup_ns", "prefix_execution_ns",
        "exact_prepared_execution_ns", "exact_call_local_total_ns",
    )
    for name in names:
        validate_metric(metrics.get(name), settings["rounds"], name)
    amortized = metrics.get("amortized_exact")
    if not isinstance(amortized, list) or (
        [row.get("reuse") for row in amortized] != settings["reuse"]
    ):
        raise CrossoverError("benchmark emitted wrong reuse rows")
    setup = float(metrics["schedule_setup_ns"]["median"])
    execution = float(metrics["exact_prepared_execution_ns"]["median"])
    prefix = float(metrics["prefix_execution_ns"]["median"])
    call_local = float(metrics["exact_call_local_total_ns"]["median"])
    if abs(float(metrics.get("prefix_over_exact_prepared", -1)) - prefix / execution) > 0.002:
        raise CrossoverError("benchmark emitted wrong prepared-exact ratio")
    if abs(float(metrics.get("prefix_over_exact_call_local", -1)) - prefix / call_local) > 0.002:
        raise CrossoverError("benchmark emitted wrong call-local ratio")
    for row in amortized:
        modeled = execution + setup / row["reuse"]
        if abs(float(row.get("modeled_ns", -1)) - modeled) > 0.002:
            raise CrossoverError("benchmark emitted wrong amortized total")
        if abs(float(row.get("prefix_over_modeled_exact", -1)) - prefix / modeled) > 0.002:
            raise CrossoverError("benchmark emitted wrong amortized ratio")


def load_attestation(path: Path, cpu: int, reserved_cpus: list[int]) -> dict[str, Any]:
    value = read_json(path)
    if not isinstance(value, dict) or value.get("schema") != ATTESTATION_SCHEMA:
        raise CrossoverError("isolation attestation has an unknown schema")
    if value.get("cpu") != cpu:
        raise CrossoverError("isolation attestation names a different CPU")
    if value.get("reserved_cpus") != reserved_cpus:
        raise CrossoverError("isolation attestation omits the exact SMT sibling set")
    if value.get("smt_sibling_idle") is not True or (
        value.get("competing_work_idle") is not True
    ):
        raise CrossoverError("isolation attestation does not prove an idle host")
    if not isinstance(value.get("operator"), str) or not value["operator"].strip():
        raise CrossoverError("isolation attestation omits the operator")
    if not isinstance(value.get("timestamp_utc"), str) or not value["timestamp_utc"].strip():
        raise CrossoverError("isolation attestation omits the timestamp")
    return value


def make_manifest(
    source: Path,
    executable: Path,
    mode: str,
    cells: list[dict[str, Any]],
    settings: dict[str, Any],
    pin_cpu: int | None,
    attestation: dict[str, Any] | None,
    protected_isolation: dict[str, Any] | None = None,
    allowed_override: list[int] | None = None,
) -> dict[str, Any]:
    cpus = list(allowed_override) if allowed_override is not None else allowed_cpus()
    fingerprint = source_fingerprint(source)
    state = git_state(source)
    executable_sha = digest_bytes(executable.read_bytes())
    build = build_metadata(executable, source)
    machine = machine_identity(cpus)
    identity = {
        "cells": cells,
        "executable_sha256": executable_sha,
        "git": state,
        "machine": machine,
        "mode": mode,
        "pin_cpu": pin_cpu,
        "settings": settings,
        "source_fingerprint": fingerprint,
        "attestation": attestation,
        "protected_isolation": protected_isolation,
        "build_metadata": build,
    }
    configuration_id = digest_value(identity)
    jobs = []
    for cell in cells:
        job_identity = {
            "cell": cell,
            "configuration_id": configuration_id,
            "executable_sha256": executable_sha,
        }
        jobs.append({
            "cell": cell,
            "configuration_id": configuration_id,
            "job_id": digest_value(job_identity)[:20],
        })
    return {
        "attestation": attestation,
        "authoritative_requested": mode == "pinned",
        "build_metadata": build,
        "configuration_id": configuration_id,
        "executable": str(executable),
        "executable_sha256": executable_sha,
        "git": state,
        "jobs": jobs,
        "machine": machine,
        "mode": mode,
        "pin_cpu": pin_cpu,
        "protected_isolation": protected_isolation,
        "schema": SCHEMA,
        "settings": settings,
        "source": str(source),
        "source_fingerprint": fingerprint,
    }


def safe_artifact(result_dir: Path, relative: str, label: str) -> Path:
    if Path(relative).is_absolute():
        raise CrossoverError(f"{label} path is absolute")
    root = result_dir.resolve()
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as error:
        raise CrossoverError(f"{label} path escapes the result directory") from error
    return path


def expected_job_command(manifest: dict[str, Any], expected: dict[str, Any]) -> list[str]:
    executable = Path(manifest["executable"])
    command = benchmark_command(executable, expected["cell"], manifest["settings"])
    if manifest["mode"] == "pinned":
        command = [
            manifest["settings"]["taskset"], "-c", str(manifest["pin_cpu"]),
            *command,
        ]
    return command


def validate_static_identity(manifest: dict[str, Any]) -> None:
    source = Path(manifest["source"])
    executable = Path(manifest["executable"])
    if source_fingerprint(source) != manifest.get("source_fingerprint"):
        raise CrossoverError("source fingerprint changed before an invocation")
    if git_state(source) != manifest.get("git"):
        raise CrossoverError("Git identity changed before an invocation")
    if not executable.is_file() or (
        digest_bytes(executable.read_bytes()) != manifest.get("executable_sha256")
    ):
        raise CrossoverError("executable identity changed before an invocation")
    if build_metadata(executable, source) != manifest.get("build_metadata"):
        raise CrossoverError("object/archive identity changed before an invocation")
    validate_cell_manifest_identity(manifest)


def validate_job_artifacts(
    result_dir: Path,
    job: dict[str, Any],
    expected: dict[str, Any],
    manifest: dict[str, Any],
) -> None:
    expected_keys = {
        "benchmark", "cell", "command", "configuration_id", "error",
        "job_id", "returncode", "schema", "status", "stderr_path",
        "stderr_sha256", "stdout_path", "stdout_sha256",
    }
    if not isinstance(job, dict) or set(job) != expected_keys:
        raise CrossoverError("job artifact keys changed")
    if job.get("schema") != JOB_SCHEMA or job.get("job_id") != expected["job_id"]:
        raise CrossoverError("job artifact has stale identity")
    if job.get("configuration_id") != manifest["configuration_id"]:
        raise CrossoverError("job artifact belongs to another manifest")
    if canonical_bytes(job.get("cell")) != canonical_bytes(expected["cell"]):
        raise CrossoverError("job artifact cell differs from the manifest")
    if job.get("command") != expected_job_command(manifest, expected):
        raise CrossoverError("job command differs from the manifest")
    if job.get("status") != "passed" or job.get("returncode") != 0 or (
        job.get("error") is not None
    ):
        raise CrossoverError("only passed jobs can be resumed")
    for stream in ("stdout", "stderr"):
        relative = job.get(f"{stream}_path")
        expected_digest = job.get(f"{stream}_sha256")
        if not isinstance(relative, str) or not isinstance(expected_digest, str):
            raise CrossoverError("job artifact omits retained streams")
        expected_relative = f"jobs/{expected['job_id']}.{stream}"
        if relative != expected_relative:
            raise CrossoverError(f"job {stream} path differs from its identity")
        path = safe_artifact(result_dir, relative, f"job {stream}")
        if not path.is_file() or digest_bytes(path.read_bytes()) != expected_digest:
            raise CrossoverError(f"job {stream} artifact changed")
    stdout = safe_artifact(result_dir, job["stdout_path"], "job stdout")
    try:
        retained_benchmark = json.loads(stdout.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError(f"cannot parse retained benchmark stdout: {error}") from error
    if canonical_bytes(retained_benchmark) != canonical_bytes(job.get("benchmark")):
        raise CrossoverError("job JSON snapshot differs from retained stdout")
    validate_benchmark_result(
        retained_benchmark, expected["cell"], manifest["settings"], manifest["git"],
        Path(manifest["executable"]),
    )


def run_job(
    job: dict[str, Any], manifest: dict[str, Any], result_dir: Path,
    executable: Path, taskset: str | None, timeout: int, resume: bool,
) -> dict[str, Any]:
    validate_static_identity(manifest)
    job_dir = result_dir / "jobs"
    result_path = job_dir / f"{job['job_id']}.json"
    if resume and result_path.is_file():
        existing = read_json(result_path)
        try:
            validate_job_artifacts(result_dir, existing, job, manifest)
        except CrossoverError:
            pass
        else:
            validate_static_identity(manifest)
            return existing
    command = benchmark_command(executable, job["cell"], manifest["settings"])
    if taskset is not None:
        command = [taskset, "-c", str(manifest["pin_cpu"]), *command]
    environment = os.environ.copy()
    environment.update(OMP_DYNAMIC="FALSE", OMP_NUM_THREADS="1")
    invocation_error = None
    try:
        completed = subprocess.run(
            command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=timeout, env=environment,
        )
        returncode = completed.returncode
        stdout_bytes = completed.stdout
        stderr_bytes = completed.stderr
    except subprocess.TimeoutExpired as caught:
        returncode = None
        stdout_bytes = caught.stdout or b""
        stderr_bytes = caught.stderr or b""
        invocation_error = f"benchmark timed out after {timeout} seconds"
    except OSError as caught:
        returncode = None
        stdout_bytes = b""
        stderr_bytes = str(caught).encode("utf-8", errors="replace")
        invocation_error = f"cannot execute benchmark: {caught}"
    validate_static_identity(manifest)
    stdout_relative = f"jobs/{job['job_id']}.stdout"
    stderr_relative = f"jobs/{job['job_id']}.stderr"
    atomic_write(result_dir / stdout_relative, stdout_bytes)
    atomic_write(result_dir / stderr_relative, stderr_bytes)
    benchmark = None
    error = None
    status = "failed"
    if returncode == 0:
        try:
            benchmark = json.loads(stdout_bytes.decode("utf-8"))
            validate_benchmark_result(
                benchmark, job["cell"], manifest["settings"], manifest["git"], executable
            )
            status = "passed"
        except (UnicodeError, ValueError, CrossoverError) as caught:
            error = str(caught)
    elif invocation_error is not None:
        error = invocation_error
    else:
        error = f"benchmark exited {returncode}"
    result = {
        "benchmark": benchmark,
        "cell": job["cell"],
        "command": command,
        "configuration_id": manifest["configuration_id"],
        "error": error,
        "job_id": job["job_id"],
        "returncode": returncode,
        "schema": JOB_SCHEMA,
        "status": status,
        "stderr_path": stderr_relative,
        "stderr_sha256": digest_bytes(stderr_bytes),
        "stdout_path": stdout_relative,
        "stdout_sha256": digest_bytes(stdout_bytes),
    }
    atomic_write_json(result_path, result)
    return result


def paired_log_inference(rounds: list[dict[str, Any]]) -> dict[str, Any]:
    contrasts = []
    for value in rounds:
        order = value.get("order")
        observations = value.get("observations_ns")
        if not isinstance(order, list) or not isinstance(observations, list):
            raise CrossoverError("paired-log inference requires raw ABBA observations")
        prefix = [float(sample) for sample, form in zip(observations, order)
                  if form == "prefix"]
        call_local = [float(sample) for sample, form in zip(observations, order)
                      if form == "call_local"]
        if len(prefix) != 2 or len(call_local) != 2:
            raise CrossoverError("paired-log inference found an unpaired round")
        contrasts.append(
            statistics.mean(map(math.log, prefix)) -
            statistics.mean(map(math.log, call_local))
        )
    if len(contrasts) != 3:
        raise CrossoverError("promotion inference requires exactly three rounds")
    mean = statistics.mean(contrasts)
    standard_error = statistics.stdev(contrasts) / math.sqrt(len(contrasts))
    # Two-sided 95% Student-t critical value for df=3-1.
    critical = 4.302652729911275
    lower = mean - critical * standard_error
    upper = mean + critical * standard_error
    return {
        "confidence": 0.95,
        "degrees_of_freedom": 2,
        "log_contrast_mean": mean,
        "log_contrast_standard_error": standard_error,
        "log_contrast_student_t_interval": [lower, upper],
        "speedup_geometric_mean": math.exp(mean),
        "speedup_student_t_interval": [math.exp(lower), math.exp(upper)],
        "gain_geometric_mean_percent": math.expm1(mean) * 100.0,
        "gain_student_t_interval_percent": [
            math.expm1(lower) * 100.0, math.expm1(upper) * 100.0
        ],
        "credible_gain_at_least_5_percent": math.expm1(lower) * 100.0 >= 5.0,
        "credible_regression_below_minus_2_percent":
            math.expm1(upper) * 100.0 < -2.0,
    }


def summarize_jobs(jobs: list[dict[str, Any]], reuse: list[int]) -> dict[str, Any]:
    passed = [job for job in jobs if job.get("status") == "passed"]
    complete = bool(jobs) and len(passed) == len(jobs)
    prepared: list[float] = []
    call_local: list[float] = []
    sparse_candidate: list[dict[str, Any]] = []
    prefix_neighbor: list[dict[str, Any]] = []
    amortized: dict[int, list[float]] = {value: [] for value in reuse}
    for job in passed:
        metrics = job["benchmark"]["metrics"]
        primary_rounds = job["benchmark"]["primary_abba"]["rounds"]
        round_gains = [
            float(round_value["candidate_gain_percent"])
            for round_value in primary_rounds
        ]
        inference_result = paired_log_inference(primary_rounds)
        prefix = float(metrics["prefix_execution_ns"]["median"])
        exact = float(metrics["exact_prepared_execution_ns"]["median"])
        setup = float(metrics["schedule_setup_ns"]["median"])
        prepared.append((prefix / exact - 1.0) * 100.0)
        call_local.append(statistics.median(round_gains))
        inference = {
            "all_rounds_positive": all(value > 0 for value in round_gains),
            "job_id": job["job_id"],
            "paired_log_student_t": inference_result,
            "round_gain_percent": round_gains,
            "round_gain_median_percent": statistics.median(round_gains),
        }
        if job["cell"].get("cell_kind") == "sparse_candidate":
            sparse_candidate.append(inference)
        else:
            prefix_neighbor.append(inference)
        for reuse_count in reuse:
            modeled = exact + setup / reuse_count
            amortized[reuse_count].append((prefix / modeled - 1.0) * 100.0)

    def summary(values: list[float]) -> dict[str, Any]:
        return {
            "cells": len(values),
            "descriptive_gain_max_percent": max(values) if values else None,
            "descriptive_gain_median_percent": statistics.median(values) if values else None,
            "descriptive_gain_min_percent": min(values) if values else None,
            "regressions": sum(value < 0 for value in values),
            "cells_at_descriptive_5_percent": sum(value >= 5.0 for value in values),
            "severe_regressions_below_minus_2_percent":
                sum(value < -2.0 for value in values),
        }

    candidate_wins = [
        value for value in sparse_candidate
        if value["paired_log_student_t"]["credible_gain_at_least_5_percent"]
    ]
    suspicious_neighbors = [
        value for value in prefix_neighbor
        if value["paired_log_student_t"]["credible_gain_at_least_5_percent"]
    ]
    credible_regressions = [
        value for value in sparse_candidate + prefix_neighbor
        if value["paired_log_student_t"]["credible_regression_below_minus_2_percent"]
    ]
    return {
        "call_local_total": summary(call_local),
        "round_level_inference": {
            "prefix_neighbors": prefix_neighbor,
            "prefix_neighbors_with_blanket_5_percent_gain": sum(
                value["paired_log_student_t"]["credible_gain_at_least_5_percent"]
                for value in prefix_neighbor
            ),
            "sparse_candidates": sparse_candidate,
            "sparse_candidates_with_lower_ci_at_least_5_percent": sum(
                value["paired_log_student_t"]["credible_gain_at_least_5_percent"]
                for value in sparse_candidate
            ),
        },
        "exact_prepared_execution": summary(prepared),
        "jobs_failed": len(jobs) - len(passed),
        "jobs_passed": len(passed),
        "jobs_total": len(jobs),
        "modeled_amortized": {
            str(value): summary(amortized[value]) for value in reuse
        },
        "promotion_gate": {
            "eligible_sparse_job_ids": [value["job_id"] for value in candidate_wins],
            "credible_regression_job_ids": [
                value["job_id"] for value in credible_regressions
            ],
            "suspicious_prefix_neighbor_job_ids": [
                value["job_id"] for value in suspicious_neighbors
            ],
            "passed": complete and bool(candidate_wins) and not suspicious_neighbors
                and not credible_regressions,
            "production_promotion_sufficient": False,
            "scope": "preliminary kernel evidence only; neighboring K/R and "
                "end-to-end public encoder gates remain external",
            "primary_only": "call-local exact versus prefix paired log contrast; "
                "prepared and modeled metrics excluded",
            "rule": "at least one sparse lower 95% Student-t bound >=5%; "
                "no prefix-neighbor lower bound >=5%; no upper bound below -2%",
        },
    }


def load_manifest(result_dir: Path) -> dict[str, Any]:
    manifest = read_json(result_dir / "manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise CrossoverError("manifest has an unknown schema")
    if not isinstance(manifest.get("jobs"), list) or not manifest.get("configuration_id"):
        raise CrossoverError("manifest is incomplete")
    expected_keys = {
        "attestation", "authoritative_requested", "build_metadata",
        "configuration_id", "executable", "executable_sha256", "git",
        "jobs", "machine", "mode", "pin_cpu", "protected_isolation",
        "schema", "settings", "source", "source_fingerprint",
    }
    if set(manifest) != expected_keys:
        raise CrossoverError("manifest keys changed")
    identity = {
        "cells": [job.get("cell") for job in manifest["jobs"]],
        "executable_sha256": manifest["executable_sha256"],
        "git": manifest["git"], "machine": manifest["machine"],
        "mode": manifest["mode"], "pin_cpu": manifest["pin_cpu"],
        "settings": manifest["settings"],
        "source_fingerprint": manifest["source_fingerprint"],
        "attestation": manifest["attestation"],
        "protected_isolation": manifest["protected_isolation"],
        "build_metadata": manifest["build_metadata"],
    }
    if digest_value(identity) != manifest["configuration_id"]:
        raise CrossoverError("manifest configuration digest is invalid")
    seen = set()
    for job in manifest["jobs"]:
        if not isinstance(job, dict) or set(job) != {
            "cell", "configuration_id", "job_id"
        } or job.get("configuration_id") != manifest["configuration_id"]:
            raise CrossoverError("manifest job identity is invalid")
        expected_job_id = digest_value({
            "cell": job["cell"],
            "configuration_id": manifest["configuration_id"],
            "executable_sha256": manifest["executable_sha256"],
        })[:20]
        if job.get("job_id") != expected_job_id or expected_job_id in seen:
            raise CrossoverError("manifest job digest is invalid or duplicated")
        seen.add(expected_job_id)
    return manifest


def load_jobs(result_dir: Path, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    expected = {job["job_id"]: job for job in manifest["jobs"]}
    job_dir = result_dir / "jobs"
    expected_names = {
        f"{job_id}{suffix}" for job_id in expected
        for suffix in (".json", ".stdout", ".stderr")
    }
    retained_names = {path.name for path in job_dir.iterdir()} if job_dir.is_dir() else set()
    if retained_names != expected_names:
        raise CrossoverError("retained job artifact set differs from manifest")
    paths = {path.stem: path for path in sorted(job_dir.glob("*.json"))}
    if set(paths) != set(expected):
        raise CrossoverError("job file set differs from manifest")
    jobs = []
    for job_id in sorted(expected):
        value = read_json(paths[job_id])
        validate_job_artifacts(result_dir, value, expected[job_id], manifest)
        jobs.append(value)
    return jobs


def require_compatible_directory(result_dir: Path, manifest: dict[str, Any]) -> None:
    path = result_dir / "manifest.json"
    if not path.exists():
        if (result_dir / "jobs").is_dir() and any((result_dir / "jobs").iterdir()):
            raise CrossoverError("result directory has jobs but no manifest")
        return
    previous = load_manifest(result_dir)
    if previous["configuration_id"] != manifest["configuration_id"]:
        raise CrossoverError("result directory belongs to a different configuration")


def _run_matrix_impl(
    arguments: argparse.Namespace,
    protected_isolation: dict[str, Any] | None = None,
    allowed_override: list[int] | None = None,
) -> int:
    source = Path(arguments.source).resolve()
    executable = Path(arguments.executable).resolve()
    result_dir = Path(arguments.result_dir).resolve()
    if not executable.is_file():
        raise CrossoverError(f"benchmark executable is missing: {executable}")
    cpus = list(allowed_override) if allowed_override is not None else allowed_cpus()
    pin_cpu: int | None = None
    taskset: str | None = None
    attestation = None
    state = git_state(source)
    validate_resume_policy(arguments.command, arguments.no_resume)
    if arguments.command == "pinned":
        if state["dirty"]:
            raise CrossoverError("pinned evidence requires a clean source tree")
        pin_cpu = arguments.cpu
        if pin_cpu is None or pin_cpu not in cpus:
            raise CrossoverError(
                f"pinned CPU must be in allowed affinity {compact_cpu_list(cpus)}"
            )
        if arguments.workers != 1:
            raise CrossoverError("pinned measurements require exactly one worker")
        taskset = shutil.which(arguments.taskset)
        if taskset is None:
            raise CrossoverError("pinned measurements require taskset")
        if not arguments.isolation_attestation:
            raise CrossoverError("pinned measurements require --isolation-attestation")
        topology = exact_smt_pair(pin_cpu, cpus)
        attestation = load_attestation(
            Path(arguments.isolation_attestation).resolve(), pin_cpu,
            topology["cpus"],
        )
    reuse = parse_csv_unsigned(arguments.reuse, "reuse", 1000000)
    if arguments.cell_manifest:
        if arguments.backends is not None or arguments.bytes is not None:
            raise CrossoverError(
                "--cell-manifest cannot be combined with --backends or --bytes"
            )
        cell_manifest = Path(arguments.cell_manifest).resolve()
        cells, cell_manifest_identity = load_cell_manifest(cell_manifest)
    else:
        shard_bytes = parse_csv_unsigned(
            arguments.default_bytes if arguments.bytes is None else arguments.bytes,
            "bytes", 1 << 30,
        )
        backends = parse_backends(
            "auto" if arguments.backends is None else arguments.backends
        )
        cells = make_cells(backends, shard_bytes)
        cell_manifest_identity = None
    if arguments.workers is None:
        arguments.workers = min(128, len(cpus), len(cells))
    settings = {
        "iterations": arguments.iterations,
        "memory_mib": arguments.memory_mib,
        "placement_policy": "taskset single CPU" if pin_cpu is not None
            else "inherited allowed affinity with independent worker processes",
        "reuse": reuse,
        "rounds": arguments.rounds,
        "runtime_allowed_cpus": (
            [pin_cpu] if pin_cpu is not None else cpus
        ) if sys.platform.startswith("linux") else [],
        "seed": arguments.seed,
        "setup_iterations": arguments.setup_iterations,
        "taskset": str(Path(taskset).resolve()) if taskset is not None else None,
        "timeout_seconds": arguments.timeout,
        "warmups": arguments.warmups,
        "workers": arguments.workers,
    }
    if cell_manifest_identity is not None:
        settings["cell_manifest"] = cell_manifest_identity
    manifest = make_manifest(
        source, executable, arguments.command, cells, settings,
        pin_cpu, attestation, protected_isolation, cpus,
    )
    if arguments.command == "pinned" and not manifest[
        "build_metadata"
    ].get("link_graph_attested"):
        raise CrossoverError(
            "pinned evidence requires a generator-specific link graph bound "
            "by the post-link sidecar"
        )
    require_compatible_directory(result_dir, manifest)
    atomic_write_json(result_dir / "manifest.json", manifest)
    expected_job_ids = {job["job_id"] for job in manifest["jobs"]}
    existing_job_ids = {
        path.stem for path in (result_dir / "jobs").glob("*.json")
    } if (result_dir / "jobs").is_dir() else set()
    extra_job_ids = sorted(existing_job_ids - expected_job_ids)
    if extra_job_ids:
        raise CrossoverError(
            "result directory contains stale job artifacts: "
            + ",".join(extra_job_ids)
        )
    if (result_dir / "jobs").is_dir():
        expected_names = {
            f"{job_id}{suffix}" for job_id in expected_job_ids
            for suffix in (".json", ".stdout", ".stderr")
        }
        extra_names = sorted(
            path.name for path in (result_dir / "jobs").iterdir()
            if path.name not in expected_names
        )
        if extra_names:
            raise CrossoverError(
                "result directory contains stale streams: " + ",".join(extra_names)
            )
    print(
        f"sparse encode crossover: mode={arguments.command} cells={len(cells)} "
        f"allowed={manifest['machine']['allowed_cpu_list']}"
        + (f" pinned={pin_cpu}" if pin_cpu is not None else ""),
        flush=True,
    )
    jobs: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=arguments.workers) as pool:
        futures = {
            pool.submit(
                run_job, job, manifest, result_dir, executable,
                taskset, arguments.timeout, not arguments.no_resume,
            ): job
            for job in manifest["jobs"]
        }
        for future in concurrent.futures.as_completed(futures):
            job = futures[future]
            result = future.result()
            jobs.append(result)
            print(f"{job['job_id']} {result['status']}", flush=True)
    jobs.sort(key=lambda value: value["job_id"])
    validate_cell_manifest_identity(manifest)
    fingerprint_after = source_fingerprint(source)
    git_after = git_state(source)
    executable_after = digest_bytes(executable.read_bytes())
    build_after = build_metadata(executable, source)
    source_changed = (
        fingerprint_after["digest"] != manifest["source_fingerprint"]["digest"]
        or git_after != manifest["git"]
    )
    executable_changed = executable_after != manifest["executable_sha256"]
    build_changed = build_after != manifest["build_metadata"]
    failed = any(job["status"] != "passed" for job in jobs)
    # The outer pinned wrapper elevates authority only after the post-run sibling
    # scheduler counters and held pair lease validate. Intermediate artifacts
    # therefore remain fail-closed if the run is interrupted.
    authoritative = False
    status = "passed" if (
        not failed and not source_changed and not executable_changed and not build_changed
    ) else "failed"
    summary = summarize_jobs(jobs, reuse)
    matrix = {
        "attestation": attestation,
        "authoritative": authoritative,
        "build_inputs_changed_during_run": build_changed,
        "executable_changed_during_run": executable_changed,
        "jobs": jobs,
        "manifest_configuration_id": manifest["configuration_id"],
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "status": status,
        "summary": summary,
    }
    atomic_write_json(result_dir / "matrix.json", matrix)
    atomic_write_json(result_dir / "summary.json", summary)
    print(
        f"matrix: {status}; passed={summary['jobs_passed']}/{summary['jobs_total']}; "
        f"authoritative={str(authoritative).lower()} ({result_dir / 'matrix.json'})"
    )
    return 0 if status == "passed" else 1


def run_matrix(arguments: argparse.Namespace) -> int:
    if arguments.command != "pinned":
        return _run_matrix_impl(arguments)
    if arguments.cpu is None:
        raise CrossoverError("pinned measurements require --cpu")
    if not hasattr(os, "sched_setaffinity"):
        raise CrossoverError("pinned isolation requires Linux scheduler affinity")
    original_affinity = allowed_cpus()
    topology = exact_smt_pair(arguments.cpu, original_affinity)
    sibling = next(value for value in topology["cpus"] if value != arguments.cpu)
    housekeeping = sorted(set(original_affinity).difference(topology["cpus"]))
    with CpuPairLease(topology["cpus"]) as lease:
        protected = {"pair_lease": lease.identity(), "topology": topology}
        os.sched_setaffinity(0, housekeeping)
        try:
            before = cpu_stat(sibling)
            before_ns = time.monotonic_ns()
            result = _run_matrix_impl(
                arguments, protected, original_affinity
            )
            after_ns = time.monotonic_ns()
            after = cpu_stat(sibling)
            if lease.identity() != protected["pair_lease"]:
                raise CrossoverError("CPU pair lease identity changed during the run")
            if exact_smt_pair(arguments.cpu, original_affinity) != topology:
                raise CrossoverError("SMT topology changed during the run")
            delta = cpu_stat_delta(before, after)
            accepted = (
                result == 0 and after_ns > before_ns
                and sibling_delta_is_idle(delta)
            )
            isolation = {
                "accepted": accepted,
                "benchmark_cpu": arguments.cpu,
                "reserved_sibling": sibling,
                "runner_housekeeping_cpus": housekeeping,
                "before_monotonic_ns": before_ns,
                "after_monotonic_ns": after_ns,
                "before_sibling_stat": before,
                "after_sibling_stat": after,
                "sibling_delta": delta,
                "maximum_sibling_nonidle_jiffies": 0,
                "protected_identity": protected,
            }
            result_dir = Path(arguments.result_dir).resolve()
            manifest = load_manifest(result_dir)
            if manifest.get("protected_isolation") != protected:
                raise CrossoverError("manifest lost the held CPU pair identity")
            matrix = read_json(result_dir / "matrix.json")
            matrix["isolation"] = isolation
            matrix["authoritative"] = bool(
                accepted and matrix.get("status") == "passed"
            )
            if not accepted:
                matrix["status"] = "failed"
            atomic_write_json(result_dir / "matrix.json", matrix)
            print(
                "pinned isolation: "
                f"accepted={str(accepted).lower()} sibling_nonidle_jiffies="
                f"{delta['nonidle_jiffies']} authoritative="
                f"{str(matrix['authoritative']).lower()}",
                flush=True,
            )
            return 0 if matrix["status"] == "passed" else 1
        finally:
            os.sched_setaffinity(0, original_affinity)


def analyze(arguments: argparse.Namespace) -> int:
    result_dir = Path(arguments.result_dir).resolve()
    manifest = load_manifest(result_dir)
    source = Path(manifest.get("source", "")).resolve()
    executable = Path(manifest.get("executable", "")).resolve()
    if source_fingerprint(source) != manifest.get("source_fingerprint"):
        raise CrossoverError("current source fingerprint differs from the manifest")
    if git_state(source) != manifest.get("git"):
        raise CrossoverError("current Git identity differs from the manifest")
    if not executable.is_file() or (
        digest_bytes(executable.read_bytes()) != manifest.get("executable_sha256")
    ):
        raise CrossoverError("current executable identity differs from the manifest")
    if build_metadata(executable, source) != manifest.get("build_metadata"):
        raise CrossoverError("current object/archive identity differs from the manifest")
    validate_cell_manifest_identity(manifest)
    jobs = load_jobs(result_dir, manifest)
    matrix = read_json(result_dir / "matrix.json")
    if not isinstance(matrix, dict) or matrix.get("schema") != SCHEMA:
        raise CrossoverError("matrix has an unknown schema")
    expected_matrix_keys = {
        "attestation", "authoritative", "build_inputs_changed_during_run",
        "executable_changed_during_run", "jobs", "manifest_configuration_id",
        "schema", "source_changed_during_run", "status", "summary",
    }
    if manifest.get("mode") == "pinned":
        expected_matrix_keys.add("isolation")
    if set(matrix) != expected_matrix_keys or (
        matrix.get("attestation") != manifest.get("attestation")
    ):
        raise CrossoverError("matrix identity keys differ from the manifest")
    if matrix.get("manifest_configuration_id") != manifest["configuration_id"]:
        raise CrossoverError("matrix belongs to another manifest")
    if canonical_bytes(matrix.get("jobs")) != canonical_bytes(jobs):
        raise CrossoverError("matrix job snapshot differs from retained jobs")
    if (matrix.get("source_changed_during_run") is not False
            or matrix.get("executable_changed_during_run") is not False
            or matrix.get("build_inputs_changed_during_run") is not False
            or matrix.get("status") != "passed"):
        raise CrossoverError("matrix retained a failed identity/status derivation")
    summary = summarize_jobs(jobs, manifest["settings"]["reuse"])
    if canonical_bytes(summary) != canonical_bytes(matrix.get("summary")):
        raise CrossoverError("matrix summary differs from validated jobs")
    if canonical_bytes(read_json(result_dir / "summary.json")) != canonical_bytes(summary):
        raise CrossoverError("retained summary differs from validated jobs")
    if manifest.get("mode") == "pinned":
        isolation = matrix.get("isolation")
        if not isinstance(isolation, dict) or (
            isolation.get("protected_identity") != manifest.get("protected_isolation")
            or isolation.get("maximum_sibling_nonidle_jiffies") != 0
        ):
            raise CrossoverError("pinned matrix isolation identity is invalid")
        retained_pair = manifest["protected_isolation"]["topology"]["cpus"]
        expected_housekeeping = sorted(set(parse_cpu_list(
            manifest["machine"]["allowed_cpu_list"]
        )).difference(retained_pair))
        if isolation.get("runner_housekeeping_cpus") != expected_housekeeping:
            raise CrossoverError("runner housekeeping affinity is invalid")
        expected_delta = cpu_stat_delta(
            isolation.get("before_sibling_stat", {}),
            isolation.get("after_sibling_stat", {}),
        )
        accepted = (
            isolation.get("after_monotonic_ns", 0)
            > isolation.get("before_monotonic_ns", 0)
            and sibling_delta_is_idle(expected_delta)
        )
        if (isolation.get("sibling_delta") != expected_delta
                or isolation.get("accepted") is not accepted
                or matrix.get("authoritative") is not (
                    accepted and matrix.get("status") == "passed"
                )):
            raise CrossoverError("pinned matrix isolation derivation is invalid")
        if not accepted or matrix.get("authoritative") is not True:
            raise CrossoverError("pinned matrix did not retain authoritative isolation")
    elif matrix.get("authoritative") is not False:
        raise CrossoverError("screen matrix improperly claimed authority")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if matrix.get("status") == "passed" else 1


def self_test() -> int:
    if parse_backends("avx512") != ["avx512"]:
        raise CrossoverError("explicit AVX-512 backend parsing failed")
    defaults = build_parser().parse_args([
        "screen", "--executable", "benchmark", "--result-dir", "results",
    ])
    if defaults.backends is not None or defaults.bytes is not None or (
        defaults.default_bytes != "64,1024"
    ):
        raise CrossoverError("sparse runner default backend changed")
    cells_a = make_cells(["avx2", "scalar"], [1024, 64])
    cells_b = make_cells(["scalar", "avx2"], [64, 1024])
    if canonical_bytes(cells_a) != canonical_bytes(cells_b) or len(cells_a) != 48:
        raise CrossoverError("cell generation is not deterministic")
    custom_value = {
        "schema": CELL_SCHEMA,
        "cells": [
            {
                "profile": "low", "field": "gf16", "K": 128, "R": 896,
                "cell_kind": "sparse_candidate", "mask_name": "edge_sparse",
                "requested_parity": [0, 127, 128, 383, 895],
                "backend": "avx2", "shard_bytes": 1024,
            },
            {
                "profile": "low", "field": "gf16", "K": 128, "R": 896,
                "cell_kind": "prefix_neighbor", "mask_name": "dense_prefix",
                "requested_parity": {"prefix": 128},
                "backend": "avx2", "shard_bytes": 1024,
            },
        ],
    }
    custom_cells = validate_cells(custom_value)
    if len(custom_cells) != 2 or custom_cells[0]["requested_parity"] != list(
        range(128)
    ):
        raise CrossoverError("custom cell manifest normalization failed")

    def reject_cell_mutation(mutator: Any) -> None:
        changed = json.loads(json.dumps(custom_value))
        mutator(changed)
        try:
            validate_cells(changed)
        except CrossoverError:
            return
        raise CrossoverError("malformed custom cell manifest was accepted")

    reject_cell_mutation(lambda value: value.update({"schema": "unknown"}))
    reject_cell_mutation(lambda value: value["cells"][0].update({"extra": True}))
    reject_cell_mutation(lambda value: value["cells"][0].update(
        {"requested_parity": [896]}))
    reject_cell_mutation(lambda value: value["cells"][0].update(
        {"shard_bytes": 1023}))
    reject_cell_mutation(lambda value: value["cells"][0].update(
        {"requested_parity": [0, 1]}))
    reject_cell_mutation(lambda value: value["cells"][0].update(
        {"field": "gf8"}))
    reject_cell_mutation(lambda value: value["cells"].append(value["cells"][0]))
    reject_cell_mutation(lambda value: value["cells"].pop())
    reject_cell_mutation(lambda value: value["cells"][1].update(
        {"mask_name": "edge_sparse"}))
    reject_cell_mutation(lambda value: value["cells"].append({
        **value["cells"][0], "mask_name": "edge_sparse_alias",
    }))
    if parse_csv_unsigned("64,1,8", "reuse", 100) != [1, 8, 64]:
        raise CrossoverError("numeric list normalization failed")
    validate_resume_policy("screen", False)
    validate_resume_policy("pinned", True)
    try:
        validate_resume_policy("pinned", False)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("pinned timing reuse was accepted")
    idle_snapshot = {"cpu": 8, "fields": [0] * 8}
    idle_after = {"cpu": 8, "fields": [0, 0, 0, 10, 0, 0, 0, 0]}
    one_busy_jiffy = {"cpu": 8, "fields": [1, 0, 0, 10, 0, 0, 0, 0]}
    if not sibling_delta_is_idle(cpu_stat_delta(idle_snapshot, idle_after)):
        raise CrossoverError("idle sibling scheduler evidence was rejected")
    if sibling_delta_is_idle(cpu_stat_delta(idle_snapshot, one_busy_jiffy)):
        raise CrossoverError("one busy sibling jiffy was accepted")
    inference_probe = paired_log_inference([
        {"order": ["prefix", "call_local", "call_local", "prefix"],
         "observations_ns": [ratio, 1.0, 1.0, ratio]}
        for ratio in (1.10, 1.20, 1.30)
    ])
    if (inference_probe["degrees_of_freedom"] != 2
            or inference_probe["confidence"] != 0.95
            or inference_probe["speedup_student_t_interval"][0]
            >= inference_probe["speedup_geometric_mean"]
            or inference_probe["speedup_student_t_interval"][1]
            <= inference_probe["speedup_geometric_mean"]):
        raise CrossoverError("paired-log Student-t inference is malformed")
    with tempfile.TemporaryDirectory(prefix="leopard2-sparse-crossover-") as directory:
        root = Path(directory)
        for relative in SOURCE_FILES:
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(relative + "\n", encoding="utf-8")
        before = source_fingerprint(root)
        (root / SOURCE_FILES[0]).write_text("changed\n", encoding="utf-8")
        after = source_fingerprint(root)
        if before["digest"] == after["digest"]:
            raise CrossoverError("source fingerprint ignored a mutation")
        cell_path = (root / "cells.json").resolve()
        atomic_write_json(cell_path, custom_value)
        loaded_cells, cell_identity = load_cell_manifest(cell_path)
        cell_bound_manifest = {
            "settings": {"cell_manifest": cell_identity},
            "jobs": [{"cell": cell} for cell in loaded_cells],
        }
        validate_cell_manifest_identity(cell_bound_manifest)
        changed_cells = json.loads(json.dumps(custom_value))
        changed_cells["cells"][0]["mask_name"] = "changed"
        atomic_write_json(cell_path, changed_cells)
        try:
            validate_cell_manifest_identity(cell_bound_manifest)
        except CrossoverError:
            pass
        else:
            raise CrossoverError("changed cell manifest retained stale authority")
        attestation = {
            "schema": ATTESTATION_SCHEMA, "cpu": 7,
            "reserved_cpus": [7, 8],
            "smt_sibling_idle": True, "competing_work_idle": True,
            "operator": "self-test", "timestamp_utc": "2000-01-01T00:00:00Z",
        }
        path = root / "attestation.json"
        atomic_write_json(path, attestation)
        if load_attestation(path, 7, [7, 8]) != attestation:
            raise CrossoverError("attestation round trip failed")
        attestation["smt_sibling_idle"] = False
        atomic_write_json(path, attestation)
        try:
            load_attestation(path, 7, [7, 8])
        except CrossoverError:
            pass
        else:
            raise CrossoverError("unsafe isolation attestation was accepted")
        with CpuPairLease([7, 8], root=root) as held:
            identity = held.identity()
            if identity["payload"]["cpus"] != [7, 8]:
                raise CrossoverError("CPU pair lease identity is malformed")
            try:
                with CpuPairLease([7, 8], root=root):
                    pass
            except CrossoverError:
                pass
            else:
                raise CrossoverError("overlapping CPU pair lease was accepted")
        with CpuPairLease([7, 8], root=root) as reacquired:
            reacquired.identity()

    with tempfile.TemporaryDirectory(prefix="leopard2-sparse-build-id-") as directory:
        fixture = Path(directory).resolve()
        fixture_source = fixture / "source"
        fixture_source.mkdir()
        def make_build(name: str) -> dict[str, Path]:
            fixture_build = fixture / name
            benchmark_target = fixture_build / (
                "CMakeFiles/leopard2_sparse_encode_benchmark_object.dir/"
                "bench/leopard2"
            )
            oracle_target = fixture_build / (
                "CMakeFiles/leopard2_sparse_encode_oracle_object.dir/tests/leopard2"
            )
            production_target = fixture_build / "CMakeFiles/leopard.dir"
            link_target = fixture_build / "CMakeFiles/bench_leopard2_sparse_encode.dir"
            for path in (benchmark_target, oracle_target, production_target, link_target):
                path.mkdir(parents=True)
            paths = {
                "executable": fixture_build / "bench_leopard2_sparse_encode",
                "archive": fixture_build / "libleopard.a",
                "benchmark_object": benchmark_target /
                    "sparse_encode_benchmark.cpp.o",
                "oracle_object": oracle_target / "direct_oracle.cpp.o",
                "production_object": production_target / "leopard2.cpp.o",
                "benchmark_recipe": link_target / "link.txt",
                "production_recipe": production_target / "link.txt",
            }
            paths["benchmark_object"].write_bytes(b"benchmark object")
            paths["oracle_object"].write_bytes(b"oracle object")
            paths["production_object"].write_bytes(b"production object")
            paths["archive"].write_bytes(b"archive")
            paths["executable"].write_bytes(b"executable")
            paths["benchmark_recipe"].write_text(
                "cc {} {} -o {} {}\n".format(
                    paths["benchmark_object"], paths["oracle_object"],
                    paths["executable"], paths["archive"]
                ), encoding="utf-8",
            )
            paths["production_recipe"].write_text(
                "ar qc {} {}\n".format(
                    paths["archive"], paths["production_object"]
                ), encoding="utf-8",
            )
            (fixture_build / "CMakeCache.txt").write_text(
                "CMAKE_HOME_DIRECTORY:INTERNAL=" + str(fixture_source) + "\n"
                "CMAKE_BUILD_TYPE:STRING=Release\n"
                "CMAKE_GENERATOR:INTERNAL=Unix Makefiles\n"
                "LEO2_BUILD_BENCHMARKS:BOOL=ON\n",
                encoding="utf-8",
            )
            sidecar = Path(str(paths["executable"]) + ".leopard2-evidence")
            identities = {
                key: file_identity(paths[key])
                for key in (
                    "executable", "archive", "benchmark_object", "oracle_object",
                    "benchmark_recipe", "production_recipe",
                )
            }
            sidecar.write_text(
                "schema=" + LINK_SIDECAR_SCHEMA + "\n"
                "executable_sha256=" + identities["executable"]["sha256"] + "\n"
                "executable_size=" + str(identities["executable"]["size"]) + "\n"
                "production_archive_sha256=" + identities["archive"]["sha256"] + "\n"
                "production_archive_size=" + str(identities["archive"]["size"]) + "\n"
                "benchmark_object_sha256=" +
                    identities["benchmark_object"]["sha256"] + "\n"
                "benchmark_object_size=" +
                    str(identities["benchmark_object"]["size"]) + "\n"
                "oracle_object_sha256=" + identities["oracle_object"]["sha256"] + "\n"
                "oracle_object_size=" + str(identities["oracle_object"]["size"]) + "\n"
                "link_recipe_kind=cmake-link-txt-v1\n"
                "benchmark_link_recipe_sha256=" +
                    identities["benchmark_recipe"]["sha256"] + "\n"
                "benchmark_link_recipe_size=" +
                    str(identities["benchmark_recipe"]["size"]) + "\n"
                "production_link_recipe_sha256=" +
                    identities["production_recipe"]["sha256"] + "\n"
                "production_link_recipe_size=" +
                    str(identities["production_recipe"]["size"]) + "\n",
                encoding="ascii",
            )
            paths["sidecar"] = sidecar
            return paths

        pristine_paths = make_build("pristine")
        build_metadata(pristine_paths["executable"], fixture_source)

        def reject_build_mutation(name: str, key: str, value: bytes) -> None:
            paths = make_build(name)
            paths[key].write_bytes(value)
            try:
                build_metadata(paths["executable"], fixture_source)
            except CrossoverError:
                return
            raise CrossoverError(f"build identity accepted a mutated {key}")

        reject_build_mutation("bad-oracle", "oracle_object", b"mutated oracle")
        reject_build_mutation("bad-benchmark", "benchmark_object", b"mutated benchmark")
        reject_build_mutation("bad-archive", "archive", b"unrelated archive")
        reject_build_mutation("bad-executable", "executable", b"manual relink")
        reject_build_mutation("bad-link", "benchmark_recipe", b"cc unrelated.a\n")
        reject_build_mutation("bad-archive-link", "production_recipe", b"ar unrelated.a\n")
        try:
            build_metadata(pristine_paths["executable"], fixture / "other-source")
        except CrossoverError:
            pass
        else:
            raise CrossoverError("build from another source root was accepted")

    with tempfile.TemporaryDirectory(prefix="leopard2-sparse-manifest-") as directory:
        manifest_root = Path(directory)
        manifest_cell = dict(base_cells()[0])
        identity = {
            "cells": [manifest_cell], "executable_sha256": "1" * 64,
            "git": {"dirty": False, "sha": "a" * 40},
            "machine": {"allowed_cpu_list": "0-2"}, "mode": "screen",
            "pin_cpu": None, "settings": {"rounds": 3},
            "source_fingerprint": {"digest": "2" * 64},
            "attestation": None, "protected_isolation": None,
            "build_metadata": {"fixture": True},
        }
        configuration_id = digest_value(identity)
        job_id = digest_value({
            "cell": manifest_cell, "configuration_id": configuration_id,
            "executable_sha256": "1" * 64,
        })[:20]
        manifest_fixture = {
            "attestation": None, "authoritative_requested": False,
            "build_metadata": identity["build_metadata"],
            "configuration_id": configuration_id,
            "executable": "/fixture/benchmark", "executable_sha256": "1" * 64,
            "git": identity["git"],
            "jobs": [{"cell": manifest_cell, "configuration_id": configuration_id,
                      "job_id": job_id}],
            "machine": identity["machine"], "mode": "screen", "pin_cpu": None,
            "protected_isolation": None, "schema": SCHEMA,
            "settings": identity["settings"], "source": "/fixture/source",
            "source_fingerprint": identity["source_fingerprint"],
        }
        atomic_write_json(manifest_root / "manifest.json", manifest_fixture)
        load_manifest(manifest_root)
        manifest_fixture["settings"]["rounds"] = 4
        atomic_write_json(manifest_root / "manifest.json", manifest_fixture)
        try:
            load_manifest(manifest_root)
        except CrossoverError:
            pass
        else:
            raise CrossoverError("manifest mutation retained a stale digest")

    cell = dict(base_cells()[0])
    cell.update(backend="scalar", shard_bytes=64)
    settings = {
        "iterations": 2, "memory_mib": 16, "reuse": [1, 8, 64],
        "rounds": 3, "runtime_allowed_cpus": [0], "seed": 7,
        "setup_iterations": 2, "warmups": 1,
    }
    state = {"dirty": False, "sha": "a" * 40}
    setup = 10.0
    prefix = 30.0
    exact = 20.0
    benchmark = {
        "schema": BENCHMARK_SCHEMA,
        "authoritative": False,
        "build": {
            "source_git_sha": state["sha"], "source_dirty": 0,
            "library_test_hooks": False,
            "compiler": "self-test", "compiler_version": "1", "cplusplus": 201103,
        },
        "runtime": {"linux_procfs_affinity_attested": True,
                    "executable_path": "/self-test/benchmark",
                    "allowed_cpus": [0]},
        "parameters": {
            "profile": cell["profile"], "field": cell["field"],
            "K": cell["K"], "R": cell["R"], "shard_bytes": 64,
            "requested_parity": cell["requested_parity"],
            "requested_backend": "scalar", "iterations": 2, "rounds": 3,
            "warmups": 1, "setup_iterations": 2, "reuse": [1, 8, 64],
            "memory_mib": 16, "seed": 7,
        },
        "resolved": {"backend": "scalar"},
        "plan": {
            "full_butterflies": 100, "prefix_butterflies": 80,
            "retained_butterflies": 50,
        },
        "correctness": {
            "exact_prefix_parity_match": True,
            "direct_generator_parity_match": True,
            "direct_generator_symbols_checked": 24,
            "direct_generator_sample_sha256": "1" * 64,
            "encode_decode_recovery_match": True,
            "input_immutable": True,
            "allocation_guards_match": True,
            "unrequested_output_canary": {"applicable": False, "match": True},
            "post_timing_recheck_match": True,
            "input_sha256": "2" * 64,
            "parity_sha256": "3" * 64,
            "recovery_source_sha256": "4" * 64,
            "recovered_original_sha256": "4" * 64,
            "digest_fnv1a64": "0x0123456789abcdef",
        },
        "primary_abba": {
            "design": "pairwise_prefix_call_local_ABBA",
            "order_policy": "alternating_ABBA_BAAB",
            "rounds": [
                {
                    "round": index,
                    "order": (["prefix", "call_local", "call_local", "prefix"]
                              if index % 2 == 0 else
                              ["call_local", "prefix", "prefix", "call_local"]),
                    "observations_ns": ([30, 25, 25, 30]
                                        if index % 2 == 0 else [25, 30, 30, 25]),
                    "prefix_median_ns": 30,
                    "call_local_median_ns": 25,
                    "log_contrast": math.log(30) - math.log(25),
                    "candidate_gain_percent": 20,
                }
                for index in range(3)
            ],
            "round_gain_percent": {
                "median": 20, "mad": 0, "minimum": 20, "maximum": 20,
                "samples": [20, 20, 20],
            },
        },
        "metrics": {
            "schedule_setup_ns": {
                "median": setup, "mad": 0, "minimum": 10, "maximum": 10,
                "samples": [10, 10, 10],
            },
            "prefix_execution_ns": {
                "median": prefix, "mad": 0, "minimum": 30, "maximum": 30,
                "samples": [30, 30, 30],
            },
            "exact_prepared_execution_ns": {
                "median": exact, "mad": 0, "minimum": 20, "maximum": 20,
                "samples": [20, 20, 20],
            },
            "exact_call_local_total_ns": {
                "median": 25, "mad": 0, "minimum": 25, "maximum": 25,
                "samples": [25, 25, 25],
            },
            "prefix_over_exact_prepared": prefix / exact,
            "prefix_over_exact_call_local": prefix / 25,
            "amortized_exact": [
                {"reuse": reuse, "modeled_ns": exact + setup / reuse,
                 "prefix_over_modeled_exact": prefix / (exact + setup / reuse)}
                for reuse in (1, 8, 64)
            ],
        },
    }
    validate_benchmark_result(benchmark, cell, settings, state)
    benchmark["build"]["source_git_sha"] = "b" * 40
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("mismatched binary source SHA was accepted")
    benchmark["build"]["source_git_sha"] = state["sha"]
    benchmark["build"]["source_dirty"] = 1
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("mismatched binary dirty marker was accepted")
    benchmark["build"]["source_dirty"] = 0
    benchmark["build"]["library_test_hooks"] = True
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("instrumented benchmark archive was accepted")
    del benchmark["build"]["library_test_hooks"]
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("benchmark without hook identity was accepted")
    benchmark["build"]["library_test_hooks"] = False
    benchmark["metrics"]["prefix_execution_ns"]["mad"] = 1
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("malformed benchmark MAD was accepted")
    benchmark["metrics"]["prefix_execution_ns"]["mad"] = 0

    with tempfile.TemporaryDirectory(prefix="leopard2-sparse-resume-") as directory:
        result_root = Path(directory).resolve()
        executable = result_root / "bench_leopard2_sparse_encode"
        executable.write_bytes(b"self-test executable")
        executable.chmod(0o700)
        resume_settings = dict(settings)
        resume_settings["taskset"] = None
        resume_manifest = {
            "configuration_id": "resume-self-test",
            "executable": str(executable), "git": state,
            "mode": "screen", "pin_cpu": None,
            "settings": resume_settings,
        }
        expected = {
            "cell": cell, "configuration_id": "resume-self-test",
            "job_id": "resume-job",
        }
        retained = json.loads(json.dumps(benchmark))
        retained["runtime"]["executable_path"] = str(executable)
        stdout_path = result_root / "jobs/resume-job.stdout"
        stderr_path = result_root / "jobs/resume-job.stderr"
        atomic_write_json(stdout_path, retained)
        atomic_write(stderr_path, b"")
        pristine = {
            "benchmark": retained, "cell": cell,
            "command": expected_job_command(resume_manifest, expected),
            "configuration_id": "resume-self-test", "error": None,
            "job_id": "resume-job", "returncode": 0,
            "schema": JOB_SCHEMA, "status": "passed",
            "stderr_path": "jobs/resume-job.stderr",
            "stderr_sha256": digest_bytes(b""),
            "stdout_path": "jobs/resume-job.stdout",
            "stdout_sha256": digest_bytes(stdout_path.read_bytes()),
        }
        validate_job_artifacts(result_root, pristine, expected, resume_manifest)

        def reject_resume_mutation(mutator: Any) -> None:
            changed = json.loads(json.dumps(pristine))
            mutator(changed)
            try:
                validate_job_artifacts(result_root, changed, expected, resume_manifest)
            except CrossoverError:
                return
            raise CrossoverError("resume mutation was accepted")

        reject_resume_mutation(lambda value: value["command"].append("--tampered"))
        reject_resume_mutation(lambda value: value.update(
            {"stdout_path": "../resume-job.stdout"}))
        reject_resume_mutation(lambda value: value.update({"stdout_sha256": "0" * 64}))
        reject_resume_mutation(lambda value: value.update({"status": "failed"}))
        reject_resume_mutation(lambda value: value["cell"].update({"K": cell["K"] + 1}))
        reject_resume_mutation(lambda value: value["benchmark"]["runtime"].update(
            {"allowed_cpus": [1]}))

        mutated_stdout = json.loads(json.dumps(retained))
        mutated_stdout["runtime"]["allowed_cpus"] = [1]
        atomic_write_json(stdout_path, mutated_stdout)
        changed = json.loads(json.dumps(pristine))
        changed["benchmark"] = mutated_stdout
        changed["stdout_sha256"] = digest_bytes(stdout_path.read_bytes())
        try:
            validate_job_artifacts(result_root, changed, expected, resume_manifest)
        except CrossoverError:
            pass
        else:
            raise CrossoverError("rehashing a runtime-affinity mutation was accepted")

    print(
        "PASS sparse encode crossover self-test cells=48 "
        "cell_mutations=10 identity_mutations=24"
    )
    return 0


def positive(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("screen", "pinned"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--source", default=".")
        sub.add_argument("--executable", required=True)
        sub.add_argument("--result-dir", required=True)
        sub.add_argument(
            "--cell-manifest",
            help="complete v2 cell matrix; replaces built-in --backends/--bytes",
        )
        sub.add_argument("--backends")
        sub.add_argument("--bytes")
        sub.set_defaults(
            default_bytes="64,1024" if command == "screen"
            else "64,1024,65536,262144"
        )
        sub.add_argument("--reuse", default="1,8,64")
        sub.add_argument("--iterations", type=positive, default=8 if command == "screen" else 64)
        sub.add_argument("--rounds", type=positive, default=3)
        sub.add_argument("--warmups", type=positive, default=1 if command == "screen" else 4)
        sub.add_argument(
            "--setup-iterations", type=positive,
            default=8 if command == "screen" else 64,
        )
        sub.add_argument("--memory-mib", type=positive, default=768)
        sub.add_argument("--seed", type=positive, default=0x535041525345454E)
        sub.add_argument("--workers", type=positive, default=None if command == "screen" else 1)
        sub.add_argument("--timeout", type=positive, default=180)
        sub.add_argument("--no-resume", action="store_true")
        sub.add_argument("--cpu", type=int)
        sub.add_argument("--taskset", default="taskset")
        sub.add_argument("--isolation-attestation")
    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("--result-dir", required=True)
    subparsers.add_parser("self-test")
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.command in ("screen", "pinned"):
        if arguments.rounds != 3:
            raise CrossoverError("the promotion design requires exactly 3 rounds")
        if arguments.command == "screen" and (
            arguments.cpu is not None or arguments.isolation_attestation is not None
        ):
            raise CrossoverError("screen mode does not accept pinned-only options")
        return run_matrix(arguments)
    if arguments.command == "analyze":
        return analyze(arguments)
    return self_test()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CrossoverError, subprocess.TimeoutExpired) as error:
        print(f"sparse encode crossover: {error}", file=sys.stderr)
        raise SystemExit(1)
