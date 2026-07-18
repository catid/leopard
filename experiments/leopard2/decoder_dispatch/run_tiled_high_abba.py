#!/usr/bin/env python3
"""Collect same-binary ABBA evidence for materialized versus tiled high decode.

Current campaigns intentionally execute one clean source commit and one binary.
The control role adds ``--force-materialized`` and the candidate role adds
``--force-tiled``.  Every other argument and workload identity is identical.
The retained manifest binds the source tree, executable, CMake cache, runner,
CPU-pair lease implementation, exact argv, raw artifacts, and CPU deltas.

Historical two-worktree v1 manifests are replayed only by
``analyze_tiled_high.py``; this runner never emits the historical schema.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import time
from typing import Any


SCHEMA = "leopard2-tiling-followup/v2"
BENCHMARK_SEED = 424242
ROUND_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
ROLE_MODE = {"control": "materialized", "candidate": "tiled"}
MODE_SELECTOR = {
    "materialized": "--force-materialized",
    "tiled": "--force-tiled",
}
HEAVY_PROCESSES = {"ninja", "ctest", "cc1", "cc1plus", "bench_leopard2"}
CPU_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal")
RUNNER_RELATIVE = "experiments/leopard2/decoder_dispatch/run_tiled_high_abba.py"
LEASE_RELATIVE = "tools/leopard2_jerasure_compare.py"


class RunError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise RunError(message)


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_output(arguments: list[str]) -> str:
    result = subprocess.run(
        arguments, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True,
    )
    if result.stderr:
        fail("unexpected stderr from " + " ".join(arguments))
    return result.stdout.strip()


def require_full_oid(value: str, label: str) -> str:
    if len(value) != 40:
        fail(f"{label} is not a full Git object ID: {value!r}")
    try:
        int(value, 16)
    except ValueError:
        fail(f"{label} is not hexadecimal: {value!r}")
    return value.lower()


def git_identity(root: Path, expected_commit: str) -> dict[str, object]:
    root = root.resolve()
    expected_commit = require_full_oid(expected_commit, "expected commit")
    head = checked_output(["git", "-C", str(root), "rev-parse", "HEAD"])
    tree = checked_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"])
    expected_tree = checked_output(
        ["git", "-C", str(root), "rev-parse", expected_commit + "^{tree}"])
    status_text = checked_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=all"]
    )
    if head != expected_commit:
        fail(f"{root}: HEAD {head} != expected {expected_commit}")
    if tree != expected_tree:
        fail(f"{root}: HEAD tree {tree} != expected tree {expected_tree}")
    if status_text:
        fail(f"{root}: source is dirty: {status_text!r}")
    return {
        "root": str(root),
        "head": head,
        "tree": tree,
        "expected_commit": expected_commit,
        "status": "clean",
        "status_sha256": hashlib.sha256(b"").hexdigest(),
    }


def relative_to_root(path: Path, root: Path, label: str) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        fail(f"{label} is outside the source root: {path}")
    return ""


def file_identity(path: Path, root: Path, label: str) -> dict[str, object]:
    path = path.absolute()
    if not path.is_file():
        fail(f"missing {label}: {path}")
    descriptor = path.stat()
    return {
        "path": str(path),
        "realpath": str(path.resolve()),
        "relative_path": relative_to_root(path, root, label),
        "sha256": sha256_file(path),
        "size": descriptor.st_size,
        "mode": stat.S_IMODE(descriptor.st_mode),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
    }


def cmake_home(cache: Path) -> Path:
    prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
    matches = [line[len(prefix):] for line in cache.read_text(
        encoding="utf-8", errors="strict").splitlines() if line.startswith(prefix)]
    if len(matches) != 1 or not matches[0]:
        fail(f"{cache}: expected one CMAKE_HOME_DIRECTORY")
    return Path(matches[0]).resolve()


def build_identity(root: Path, binary_relative: str) -> dict[str, object]:
    root = root.resolve()
    relative = Path(binary_relative)
    if relative.is_absolute() or ".." in relative.parts:
        fail("--binary must be a source-root-relative path")
    binary = root / relative
    cache = binary.parent / "CMakeCache.txt"
    if not os.access(binary, os.X_OK):
        fail(f"benchmark executable is not executable: {binary}")
    if not cache.is_file():
        fail(f"missing CMake cache: {cache}")
    if cmake_home(cache) != root:
        fail(f"{cache}: CMAKE_HOME_DIRECTORY does not name {root}")
    return {
        "binary": file_identity(binary, root, "benchmark executable"),
        "cache": file_identity(cache, root, "CMake cache"),
        "cmake_home": str(root),
    }


def require_output_preserves_clean_source(output: Path, root: Path) -> None:
    try:
        output.relative_to(root)
    except ValueError:
        return
    result = subprocess.run(
        ["git", "-C", str(root), "check-ignore", "-q", str(output)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 1:
        fail("output inside the source root must be covered by .gitignore")
    if result.returncode != 0:
        fail(f"cannot verify output ignore policy ({result.returncode}): "
             f"{result.stderr.decode('utf-8', 'replace')}")


def active_heavy_processes() -> list[str]:
    found = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            command = (entry / "comm").read_text(encoding="utf-8").strip()
            arguments = (entry / "cmdline").read_bytes().replace(
                b"\0", b" ").decode("utf-8", "replace")
            project_python = command.startswith("python") and (
                "/home/catid/leopard" in arguments or "leopard2" in arguments)
            if command in HEAVY_PROCESSES or project_python:
                found.append(f"{entry.name} {command} {arguments}")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
    return found


def cpu_snapshot(cpu: int) -> dict[str, int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            values = [int(value) for value in line.split()[1:1 + len(CPU_FIELDS)]]
            if len(values) != len(CPU_FIELDS):
                fail(f"CPU {cpu} has an incomplete /proc/stat row")
            return dict(zip(CPU_FIELDS, values))
    fail(f"CPU {cpu} is absent from /proc/stat")
    return {}


def cpu_delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    if set(before) != set(CPU_FIELDS) or set(after) != set(CPU_FIELDS):
        fail("CPU snapshot fields changed")
    result = {key: after[key] - before[key] for key in CPU_FIELDS}
    if any(value < 0 for value in result.values()):
        fail("CPU counters moved backwards")
    result["idle_total"] = result["idle"] + result["iowait"]
    result["nonidle_total"] = sum(
        result[key] for key in CPU_FIELDS if key not in {"idle", "iowait"})
    result["total"] = result["idle_total"] + result["nonidle_total"]
    return result


def load_pair_lease(source: Path):
    specification = importlib.util.spec_from_file_location(
        "leopard2_tiled_high_pair_lease", source)
    if specification is None or specification.loader is None:
        fail(f"cannot import pair lease from {source}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    if not hasattr(module, "PairLease"):
        fail(f"{source} does not provide PairLease")
    return module.PairLease


def checkpoint_cases() -> tuple[tuple, ...]:
    cases = []
    for kib in (16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96):
        for loss in (1, 4, 8, 16, 32):
            cases.append((f"avx2-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "avx2", 1))
    for kib in (24, 32, 40, 48, 56, 64, 80, 96, 112):
        for loss in (1, 8, 16):
            cases.append((f"ssse3-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "ssse3", 1))
    for kib in (32, 48, 64, 80, 96):
        for loss in (1, 8, 16):
            cases.append((f"scalar-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "scalar", 1))
    for kib in (32, 48, 64, 80):
        for loss in (1, 8, 16):
            cases.append((f"avx2-t32-b{kib}k-l{loss}-batch2", 224, 32,
                          kib * 1024, 4, 4, loss, "avx2", 2))
    cases.append(("avx2-t32-b64k-l8-reuse1", 224, 32, 64 * 1024,
                  8, 1, 8, "avx2", 1))
    cases.append(("avx2-t32-b64k-l8-reuse64", 224, 32, 64 * 1024,
                  8, 64, 8, "avx2", 1))
    for recovery in (16, 64):
        original = 256 - recovery
        for kib in (32, 64, 96):
            cases.append((f"avx2-t{recovery}-b{kib}k-l8", original, recovery,
                          kib * 1024, 8, 8, 8, "avx2", 1))
    return tuple(cases)


def smoke_cases() -> tuple[tuple, ...]:
    return (
        ("avx2-t32-b32k-l8", 224, 32, 32 * 1024, 2, 2, 8, "avx2", 1),
        ("avx2-t32-b96k-l8", 224, 32, 96 * 1024, 2, 2, 8, "avx2", 1),
    )


def benchmark_command(
    executable: Path, case: tuple, output: Path, cpu: int,
    iterations: int, variant: str,
) -> list[str]:
    if variant not in ROLE_MODE:
        fail(f"unknown benchmark role: {variant}")
    _name, k, r, byte_count, warmup, reuse, loss, backend, batch = case
    selector = MODE_SELECTOR[ROLE_MODE[variant]]
    return [
        "/usr/bin/taskset", "-c", str(cpu),
        "/usr/bin/env", "OMP_NUM_THREADS=1", "OMP_DYNAMIC=false",
        "OMP_PROC_BIND=close", str(executable.absolute()),
        "--k", str(k), "--r", str(r), "--profile", "high",
        "--field", "gf8", "--backend", backend, "--bytes", str(byte_count),
        "--loss", str(loss), "--batch", str(batch), "--reuse", str(reuse),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(BENCHMARK_SEED),
        "--force-specialized", selector, "--skip-legacy", "--retain-samples",
        "--json", str(output.absolute()),
    ]


def artifact_identity(path: Path) -> dict[str, object]:
    if not path.is_file():
        fail(f"missing retained artifact: {path}")
    return {"name": path.name, "size": path.stat().st_size, "sha256": sha256_file(path)}


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    manifest = dict(manifest)
    manifest["content_sha256"] = canonical_sha256(manifest)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as output:
            json.dump(manifest, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root")
    parser.add_argument("--source-commit")
    parser.add_argument("--binary", default="build-audit/bench_leopard2")
    parser.add_argument("--lease-source", default="tools/leopard2_jerasure_compare.py")
    parser.add_argument("--output")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--sibling", type=int)
    parser.add_argument("--matrix", choices=("checkpoint", "smoke"), default="checkpoint")
    parser.add_argument("--iterations", type=int, default=11)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def self_test() -> None:
    cases = checkpoint_cases()
    if len(cases) != 117 or len({case[0] for case in cases}) != len(cases):
        fail("checkpoint matrix count or names changed")
    if len(smoke_cases()) != 2:
        fail("smoke matrix count changed")
    materialized = benchmark_command(
        Path("/tmp/bench"), cases[0], Path("/tmp/control.json"), 7, 11, "control")
    tiled = benchmark_command(
        Path("/tmp/bench"), cases[0], Path("/tmp/candidate.json"), 7, 11, "candidate")
    if materialized[7] != tiled[7]:
        fail("roles do not use the same binary argv[0]")
    if "--force-materialized" not in materialized or "--force-tiled" in materialized:
        fail("control selector mapping changed")
    if "--force-tiled" not in tiled or "--force-materialized" in tiled:
        fail("candidate selector mapping changed")
    common_control = [x for x in materialized if x != "--force-materialized" and
                      not x.endswith("control.json")]
    common_candidate = [x for x in tiled if x != "--force-tiled" and
                        not x.endswith("candidate.json")]
    if common_control != common_candidate:
        fail("role commands differ outside selector and JSON destination")
    if list(map(list, ROUND_ORDERS)) != [
            ["control", "candidate", "candidate", "control"],
            ["candidate", "control", "control", "candidate"],
            ["control", "candidate", "candidate", "control"]]:
        fail("round schedule changed")
    before = {key: 10 for key in CPU_FIELDS}
    after = {key: 11 for key in CPU_FIELDS}
    delta = cpu_delta(before, after)
    if delta["total"] != 8 or delta["nonidle_total"] != 6:
        fail("CPU delta accounting changed")
    print("tiled-high ABBA runner self-test passed: current v2, 117 cells, 1404 invocations")


def main() -> int:
    args = parse_arguments()
    if args.self_test:
        self_test()
        return 0
    required = ("source_root", "source_commit", "output", "cpu", "sibling")
    missing = ["--" + name.replace("_", "-") for name in required
               if getattr(args, name) is None]
    if missing:
        fail("missing required arguments: " + ", ".join(missing))
    if args.cpu == args.sibling or args.cpu < 0 or args.sibling < 0:
        fail("CPU and sibling must be distinct non-negative logical CPUs")
    if args.iterations < 3:
        fail("at least three internal timing samples are required")
    affinity = os.sched_getaffinity(0)
    if args.cpu not in affinity or args.sibling not in affinity:
        fail("requested CPU pair is outside process affinity")
    housekeeping = affinity - {args.cpu, args.sibling}
    if not housekeeping:
        fail("runner needs a housekeeping CPU outside the timed SMT pair")

    output = Path(args.output).absolute()
    if output.exists():
        fail(f"output already exists: {output}")
    source_root = Path(args.source_root).resolve()
    require_output_preserves_clean_source(output, source_root)
    source = git_identity(source_root, args.source_commit)
    build = build_identity(source_root, args.binary)
    binary = Path(build["binary"]["path"])
    runner = file_identity(Path(__file__), source_root, "runner")
    lease_path = Path(args.lease_source)
    if not lease_path.is_absolute():
        lease_path = source_root / lease_path
    lease_source = file_identity(lease_path, source_root, "pair-lease source")
    if runner["relative_path"] != RUNNER_RELATIVE:
        fail("runner is not executing from its canonical source path")
    if lease_source["relative_path"] != LEASE_RELATIVE:
        fail("pair-lease implementation is not the canonical source path")
    active = active_heavy_processes()
    if active:
        fail("heavy processes active before campaign: " + repr(active))

    cases = checkpoint_cases() if args.matrix == "checkpoint" else smoke_cases()
    output.mkdir(parents=True)
    raw = output / "raw"
    raw.mkdir()
    PairLease = load_pair_lease(Path(lease_source["path"]))
    records: list[dict[str, Any]] = []
    campaign_cpu_before = cpu_snapshot(args.cpu)
    campaign_sibling_before = cpu_snapshot(args.sibling)
    campaign_start_ns = time.time_ns()
    original_affinity = set(affinity)
    os.sched_setaffinity(0, housekeeping)
    try:
        with PairLease(args.cpu, args.sibling) as lease:
            lease_identity = lease.revalidate()
            lease_digest = canonical_sha256(lease_identity)
            sequence = 0
            for case in cases:
                name = case[0]
                for round_index, order in enumerate(ROUND_ORDERS):
                    for slot, variant in enumerate(order):
                        active = active_heavy_processes()
                        if active:
                            fail("heavy process appeared during campaign: " + repr(active))
                        current_lease = lease.revalidate()
                        if canonical_sha256(current_lease) != lease_digest:
                            fail("pair lease identity changed during campaign")
                        prefix = raw / f"{name}-round{round_index}-slot{slot}-{variant}"
                        json_path = prefix.with_suffix(".json")
                        stdout_path = prefix.with_suffix(".stdout")
                        stderr_path = prefix.with_suffix(".stderr")
                        command = benchmark_command(
                            binary, case, json_path, args.cpu, args.iterations, variant)
                        cpu_before = cpu_snapshot(args.cpu)
                        sibling_before = cpu_snapshot(args.sibling)
                        start_ns = time.time_ns()
                        completed = subprocess.run(
                            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        end_ns = time.time_ns()
                        cpu_after = cpu_snapshot(args.cpu)
                        sibling_after = cpu_snapshot(args.sibling)
                        stdout_path.write_bytes(completed.stdout)
                        stderr_path.write_bytes(completed.stderr)
                        if completed.returncode != 0:
                            fail(f"benchmark failed ({completed.returncode}): {command!r}")
                        if completed.stdout or completed.stderr or not json_path.is_file():
                            fail("benchmark emitted output or omitted its JSON result")
                        try:
                            raw_document = json.loads(json_path.read_text(encoding="utf-8"))
                            digests = raw_document["workload_digests"]
                        except (OSError, ValueError, KeyError, TypeError) as error:
                            fail(f"cannot read benchmark workload digests: {error}")
                        mode = ROLE_MODE[variant]
                        records.append({
                            "sequence": sequence,
                            "case": name,
                            "round": round_index,
                            "slot": slot,
                            "variant": variant,
                            "mode": mode,
                            "selector": MODE_SELECTOR[mode],
                            "command": command,
                            "command_sha256": canonical_sha256(command),
                            "benchmark_argv0": str(binary.absolute()),
                            "json_path_at_run": str(json_path.absolute()),
                            "iterations": args.iterations,
                            "seed": BENCHMARK_SEED,
                            "lease_identity_sha256": lease_digest,
                            "start_time_ns": start_ns,
                            "end_time_ns": end_ns,
                            "cpu_before": cpu_before,
                            "cpu_after": cpu_after,
                            "cpu_delta": cpu_delta(cpu_before, cpu_after),
                            "sibling_before": sibling_before,
                            "sibling_after": sibling_after,
                            "sibling_delta": cpu_delta(sibling_before, sibling_after),
                            "artifacts": {
                                "json": artifact_identity(json_path),
                                "stdout": artifact_identity(stdout_path),
                                "stderr": artifact_identity(stderr_path),
                            },
                            "workload_digests_sha256": canonical_sha256(digests),
                        })
                        sequence += 1
            final_lease = lease.revalidate()
            if canonical_sha256(final_lease) != lease_digest:
                fail("pair lease identity changed before release")
    finally:
        os.sched_setaffinity(0, original_affinity)

    active = active_heavy_processes()
    if active:
        fail("heavy processes active after campaign: " + repr(active))
    campaign_cpu_after = cpu_snapshot(args.cpu)
    campaign_sibling_after = cpu_snapshot(args.sibling)
    final = {
        "source": git_identity(source_root, args.source_commit),
        "build": build_identity(source_root, args.binary),
        "runner": file_identity(Path(__file__), source_root, "runner"),
        "lease_source": file_identity(Path(lease_source["path"]), source_root,
                                      "pair-lease source"),
    }
    if final != {"source": source, "build": build, "runner": runner,
                 "lease_source": lease_source}:
        fail("source/build/runner/lease identity changed during campaign")
    campaign_end_ns = time.time_ns()
    manifest: dict[str, object] = {
        "schema": SCHEMA,
        "status": "complete",
        "valid": True,
        "source": source,
        "build": build,
        "runner": runner,
        "lease_source": lease_source,
        "campaign": {
            "cpu": args.cpu,
            "sibling": args.sibling,
            "original_affinity": sorted(original_affinity),
            "housekeeping_affinity": sorted(housekeeping),
            "matrix": args.matrix,
            "iterations": args.iterations,
            "seed": BENCHMARK_SEED,
            "cases": [list(case) for case in cases],
            "round_orders": [list(order) for order in ROUND_ORDERS],
            "lease": lease_identity,
            "lease_identity_sha256": lease_digest,
            "start_time_ns": campaign_start_ns,
            "end_time_ns": campaign_end_ns,
            "cpu_before": campaign_cpu_before,
            "cpu_after": campaign_cpu_after,
            "cpu_delta": cpu_delta(campaign_cpu_before, campaign_cpu_after),
            "sibling_before": campaign_sibling_before,
            "sibling_after": campaign_sibling_after,
            "sibling_delta": cpu_delta(campaign_sibling_before,
                                       campaign_sibling_after),
        },
        "records": records,
        "final": final,
    }
    write_manifest(output / "manifest.json", manifest)
    print(output / "manifest.json")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"run_tiled_high_abba.py: {error}", file=sys.stderr)
        raise SystemExit(1)
