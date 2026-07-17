#!/usr/bin/env python3
"""Run a pinned ABBA comparison of materialized and tiled high decoders.

The two executables normally come from detached worktrees at the immediate
pre-tiling control and the integrated tiling endpoint.  The runner records
their exact source and binary identities, reserves one physical core and its
SMT sibling, rejects competing compiler/test/benchmark activity, and retains
every child result.  It does not decide promotion; analyze_tiled_high.py does.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time


SCHEMA = "leopard2-tiling-followup/v1"
ROUND_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
HEAVY_PROCESSES = {"ninja", "ctest", "cc1", "cc1plus", "bench_leopard2"}


class RunError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise RunError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_output(arguments: list[str]) -> str:
    result = subprocess.run(arguments, check=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if result.stderr:
        fail("unexpected stderr from " + " ".join(arguments))
    return result.stdout.strip()


def git_identity(root: Path, expected_commit: str) -> dict:
    if len(expected_commit) != 40:
        fail(f"expected commit is not a full Git object ID: {expected_commit!r}")
    try:
        int(expected_commit, 16)
    except ValueError:
        fail(f"expected commit is not hexadecimal: {expected_commit!r}")
    head = checked_output(["git", "-C", str(root), "rev-parse", "HEAD"])
    tree = checked_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"])
    status = checked_output([
        "git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"])
    if head != expected_commit:
        fail(f"{root}: HEAD {head} != expected {expected_commit}")
    if status:
        fail(f"{root}: tracked source is dirty")
    return {"root": str(root), "head": head, "tree": tree, "tracked_status": "clean"}


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


def cpu_stat(cpu: int) -> list[int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            return [int(value) for value in line.split()[1:]]
    fail(f"CPU {cpu} is absent from /proc/stat")
    return []


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
    # Target shape: N=256, T=32.  The dense single-stripe map brackets the
    # observed mid-size AVX2 crossover by byte count and requested outputs.
    for kib in (16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96):
        for loss in (1, 4, 8, 16, 32):
            cases.append((f"avx2-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "avx2", 1))
    # Lower backends are independently scoped rather than inheriting AVX2.
    for kib in (24, 32, 40, 48, 56, 64, 80, 96, 112):
        for loss in (1, 8, 16):
            cases.append((f"ssse3-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "ssse3", 1))
    for kib in (32, 48, 64, 80, 96):
        for loss in (1, 8, 16):
            cases.append((f"scalar-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "scalar", 1))
    # Batch execution has a different cache reuse pattern.  A future batch
    # dispatcher may retain tiled execution even where a single stripe does not.
    for kib in (32, 48, 64, 80):
        for loss in (1, 8, 16):
            cases.append((f"avx2-t32-b{kib}k-l{loss}-batch2", 224, 32,
                          kib * 1024, 4, 4, loss, "avx2", 2))
    # Execution is plan-reuse invariant, but these timing geometries verify the
    # benchmark/amortization boundary explicitly.
    cases.append(("avx2-t32-b64k-l8-reuse1", 224, 32, 64 * 1024,
                  8, 1, 8, "avx2", 1))
    cases.append(("avx2-t32-b64k-l8-reuse64", 224, 32, 64 * 1024,
                  8, 64, 8, "avx2", 1))
    # Neighboring transform sides prevent an over-broad T>=32 policy.
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


def benchmark_command(executable: Path, case: tuple, output: Path, cpu: int,
                      iterations: int) -> list[str]:
    _name, k, r, byte_count, warmup, reuse, loss, backend, batch = case
    return [
        "/usr/bin/taskset", "-c", str(cpu),
        "/usr/bin/env", "OMP_NUM_THREADS=1", "OMP_DYNAMIC=false",
        "OMP_PROC_BIND=close", str(executable),
        "--k", str(k), "--r", str(r), "--profile", "high",
        "--field", "gf8", "--backend", backend, "--bytes", str(byte_count),
        "--loss", str(loss), "--batch", str(batch), "--reuse", str(reuse),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", "424242", "--force-specialized",
        "--skip-legacy", "--retain-samples", "--json", str(output),
    ]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-root")
    parser.add_argument("--candidate-root")
    parser.add_argument("--control-commit")
    parser.add_argument("--candidate-commit")
    parser.add_argument("--control-binary", default="build-audit/bench_leopard2")
    parser.add_argument("--candidate-binary", default="build-audit/bench_leopard2")
    parser.add_argument("--lease-source", default="tools/leopard2_jerasure_compare.py")
    parser.add_argument("--output")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--sibling", type=int)
    parser.add_argument("--matrix", choices=("checkpoint", "smoke"),
                        default="checkpoint")
    parser.add_argument("--iterations", type=int, default=11)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def self_test() -> None:
    cases = checkpoint_cases()
    if len(cases) != 117 or len({case[0] for case in cases}) != len(cases):
        fail("checkpoint matrix count or names changed")
    if len(smoke_cases()) != 2:
        fail("smoke matrix count changed")
    command = benchmark_command(Path("/tmp/bench"), cases[0], Path("/tmp/out"), 7, 11)
    expected = ("--k", "--r", "--backend", "--force-specialized",
                "--skip-legacy", "--retain-samples")
    if any(value not in command for value in expected):
        fail("benchmark command omitted a required identity")
    if len(ROUND_ORDERS) != 3 or any(
            sorted(order) != ["candidate", "candidate", "control", "control"]
            for order in ROUND_ORDERS):
        fail("round schedule changed")
    print("tiled-high ABBA runner self-test passed: 117 cells, 1404 invocations")


def main() -> int:
    args = parse_arguments()
    if args.self_test:
        self_test()
        return 0
    required = ("control_root", "candidate_root", "control_commit",
                "candidate_commit", "output", "cpu", "sibling")
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

    output = Path(args.output).absolute()
    if output.exists():
        fail(f"output already exists: {output}")
    control_root = Path(args.control_root).absolute()
    candidate_root = Path(args.candidate_root).absolute()
    control_binary = control_root / args.control_binary
    candidate_binary = candidate_root / args.candidate_binary
    lease_source = Path(args.lease_source).absolute()
    for path, label in ((control_binary, "control executable"),
                        (candidate_binary, "candidate executable"),
                        (lease_source, "pair-lease source")):
        if not path.is_file():
            fail(f"missing {label}: {path}")

    source = {
        "control": git_identity(control_root, args.control_commit),
        "candidate": git_identity(candidate_root, args.candidate_commit),
    }
    identities = {
        "control_binary": sha256_file(control_binary),
        "candidate_binary": sha256_file(candidate_binary),
        "lease_source": sha256_file(lease_source),
    }
    for label, root in (("control", control_root), ("candidate", candidate_root)):
        cache = root / Path(args.control_binary if label == "control" else
                            args.candidate_binary).parent / "CMakeCache.txt"
        if not cache.is_file():
            fail(f"missing {label} CMake cache: {cache}")
        identities[label + "_cache"] = sha256_file(cache)
    active = active_heavy_processes()
    if active:
        fail("heavy processes active before campaign: " + repr(active))

    cases = checkpoint_cases() if args.matrix == "checkpoint" else smoke_cases()
    output.mkdir(parents=True)
    raw = output / "raw"
    raw.mkdir()
    start = {"cpu": cpu_stat(args.cpu), "sibling": cpu_stat(args.sibling),
             "loadavg": Path("/proc/loadavg").read_text(encoding="ascii").strip(),
             "time_ns": time.time_ns()}
    PairLease = load_pair_lease(lease_source)
    records = []
    with PairLease(args.cpu, args.sibling) as lease:
        lease_identity = lease.revalidate()
        for case in cases:
            name = case[0]
            for round_index, order in enumerate(ROUND_ORDERS):
                for slot, variant in enumerate(order):
                    active = active_heavy_processes()
                    if active:
                        fail("heavy process appeared during campaign: " + repr(active))
                    lease.revalidate()
                    executable = (control_binary if variant == "control"
                                  else candidate_binary)
                    prefix = raw / f"{name}-round{round_index}-slot{slot}-{variant}"
                    json_path = prefix.with_suffix(".json")
                    command = benchmark_command(
                        executable, case, json_path, args.cpu, args.iterations)
                    completed = subprocess.run(
                        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    prefix.with_suffix(".stdout").write_bytes(completed.stdout)
                    prefix.with_suffix(".stderr").write_bytes(completed.stderr)
                    if completed.returncode != 0:
                        fail(f"benchmark failed ({completed.returncode}): {command!r}")
                    if completed.stdout or completed.stderr or not json_path.is_file():
                        fail("benchmark emitted output or omitted its JSON result")
                    records.append({
                        "case": name, "round": round_index, "slot": slot,
                        "variant": variant, "command": command,
                        "json": json_path.name, "json_sha256": sha256_file(json_path),
                    })
    active = active_heavy_processes()
    if active:
        fail("heavy processes active after campaign: " + repr(active))
    end = {"cpu": cpu_stat(args.cpu), "sibling": cpu_stat(args.sibling),
           "loadavg": Path("/proc/loadavg").read_text(encoding="ascii").strip(),
           "time_ns": time.time_ns()}
    manifest = {
        "schema": SCHEMA, "valid": True,
        "control_commit": args.control_commit,
        "candidate_commit": args.candidate_commit,
        "source": source, "identities": identities,
        "runner_sha256": sha256_file(Path(__file__)),
        "cpu": args.cpu, "sibling": args.sibling,
        "lease": lease_identity, "start": start, "end": end,
        "matrix": args.matrix, "iterations": args.iterations,
        "cases": [list(case) for case in cases],
        "round_orders": [list(order) for order in ROUND_ORDERS],
        "records": records,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output / "manifest.json")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"run_tiled_high_abba.py: {error}", file=sys.stderr)
        raise SystemExit(1)
