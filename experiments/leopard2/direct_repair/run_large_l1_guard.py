#!/usr/bin/env python3
"""Memory-capped 16-MiB/reuse-64 guard for generalized GF8 AVX2 L=1.

Run this coordinator on a CPU outside the measured SMT pair and inside a
memory-capped cgroup.  It freezes the supplied executables, runs exactly one
benchmark process at a time, retains every attempt, rejects material sibling
activity, and emits paired-log-ratio confidence intervals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-general-l1-large-guard/v1"
SHARD_BYTES = 16 << 20
REUSE = 64
T95 = {4: 3.182446305284263, 6: 2.570581835636314, 10: 2.262157162740992}
ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
)
MAIN_ORDERS = (
    ("main", "candidate", "candidate", "main"),
    ("candidate", "main", "main", "candidate"),
)
CELLS = (
    {"id": "high-target", "K": 17, "R": 16, "profile": "high",
     "loss": 1, "role": "target", "compare_main": True},
    {"id": "low-target", "K": 17, "R": 31, "profile": "low",
     "loss": 1, "role": "target", "compare_main": False},
    {"id": "high-neighbor-loss2", "K": 17, "R": 16,
     "profile": "high", "loss": 2, "role": "neighbor",
     "compare_main": False},
    {"id": "low-neighbor-loss2", "K": 17, "R": 31,
     "profile": "low", "loss": 2, "role": "neighbor",
     "compare_main": False},
)
MAX_JSON_BYTES = 4 << 20
TIME_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


class GuardError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise GuardError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            block = source.read(1 << 20)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def write_json(path: Path, value: object) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) +
               "\n").encode("utf-8")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                         0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            require(written > 0, f"short write: {temporary}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, Any]:
    status = path.stat()
    require(0 < status.st_size <= MAX_JSON_BYTES,
            f"invalid benchmark JSON size: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"benchmark JSON is not an object: {path}")
    return value


def freeze(source: Path, destination: Path) -> dict[str, Any]:
    resolved = source.resolve(strict=True)
    require(resolved.is_file(), f"not a regular executable: {resolved}")
    with resolved.open("rb") as reader, destination.open("xb") as writer:
        shutil.copyfileobj(reader, writer, 1 << 20)
        writer.flush()
        os.fsync(writer.fileno())
    destination.chmod(0o555)
    return {
        "source": str(resolved),
        "snapshot": str(destination.resolve()),
        "size": destination.stat().st_size,
        "sha256": sha256(destination),
        "mode": oct(destination.stat().st_mode & 0o777),
    }


def cpu_nonidle(cpu: int) -> int:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            fields = [int(value) for value in line.split()[1:]]
            require(len(fields) >= 8, f"short /proc/stat row for CPU {cpu}")
            return sum(fields[index] for index in (0, 1, 2, 5, 6, 7))
    raise GuardError(f"CPU {cpu} is absent from /proc/stat")


def cgroup_state() -> dict[str, Any]:
    unified = None
    for line in Path("/proc/self/cgroup").read_text(encoding="ascii").splitlines():
        fields = line.split(":", 2)
        if len(fields) == 3 and fields[0] == "0" and fields[1] == "":
            unified = fields[2]
            break
    require(unified is not None, "unified cgroup path is unavailable")
    root = Path("/sys/fs/cgroup") / unified.lstrip("/")

    def read(name: str) -> str:
        return (root / name).read_text(encoding="ascii").strip()

    events = {}
    for line in read("memory.events").splitlines():
        name, value = line.split()
        events[name] = int(value)
    return {
        "path": str(root),
        "memory_max": read("memory.max"),
        "memory_swap_max": read("memory.swap.max"),
        "memory_peak": int(read("memory.peak")),
        "memory_events": events,
    }


def benchmark_arguments(
    executable: Path,
    implementation: str,
    cell: Mapping[str, Any],
    cpu: int,
    seed: int,
    json_path: Path,
    iterations: int,
    warmup: int,
) -> list[str]:
    common = [
        "env", "LANG=C", "LC_ALL=C", "OMP_DYNAMIC=FALSE",
        "OMP_NUM_THREADS=1", "taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(SHARD_BYTES), "--loss", str(cell["loss"]),
        "--batch", "1", "--reuse", str(REUSE),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(seed), "--json", str(json_path),
    ]
    if implementation == "main":
        common.extend(("--logical-bytes", str(SHARD_BYTES)))
        return common
    common.extend((
        "--profile", str(cell["profile"]), "--field", "gf8",
        "--backend", "avx2", "--skip-legacy", "--retain-samples",
        "--report-decode-path", "--report-direct-executor",
    ))
    if implementation == "control":
        common.append("--force-specialized")
    return common


def validate_result(
    implementation: str,
    cell: Mapping[str, Any],
    value: Mapping[str, Any],
) -> None:
    if implementation == "main":
        require(value["correctness"]["round_trip"], "main round trip failed")
        return
    require(value["correctness"]["leopard2_round_trip"],
            f"{implementation} round trip failed")
    resolved = value["resolved"]
    require(resolved["field"] == "gf8" and resolved["backend"] == "avx2",
            f"{implementation} selected an unexpected field/backend")
    path = resolved["selected_decode_path"]
    if cell["role"] == "target" and implementation == "candidate":
        require(path == "direct", f"target candidate did not select direct: {path}")
    elif implementation == "control":
        require(path != "direct", f"forced control selected direct: {path}")
    elif cell["role"] == "neighbor":
        require(path != "direct", f"neighbor candidate selected direct: {path}")


def run_one(
    binary: Path,
    implementation: str,
    cell: Mapping[str, Any],
    cpu: int,
    sibling: int,
    seed: int,
    raw: Path,
    stem: str,
    iterations: int,
    warmup: int,
    max_attempts: int,
    max_sibling_jiffies: int,
) -> dict[str, Any]:
    attempts = []
    for attempt in range(max_attempts):
        json_path = raw / f"{stem}-attempt{attempt}.json"
        time_path = raw / f"{stem}-attempt{attempt}.time"
        command = ["/usr/bin/time", "-v", "-o", str(time_path)]
        command.extend(benchmark_arguments(
            binary, implementation, cell, cpu, seed, json_path,
            iterations, warmup))
        sibling_before = cpu_nonidle(sibling)
        completed = subprocess.run(command, check=False)
        sibling_after = cpu_nonidle(sibling)
        require(completed.returncode == 0,
                f"benchmark failed ({completed.returncode}): {' '.join(command)}")
        value = load_json(json_path)
        validate_result(implementation, cell, value)
        match = TIME_RSS_RE.search(time_path.read_text(encoding="utf-8"))
        require(match is not None, f"maximum RSS is absent: {time_path}")
        sibling_delta = sibling_after - sibling_before
        record = {
            "attempt": attempt,
            "json": str(json_path),
            "time": str(time_path),
            "maximum_rss_kib": int(match.group(1)),
            "sibling_nonidle_jiffies": sibling_delta,
            "accepted": sibling_delta <= max_sibling_jiffies,
        }
        attempts.append(record)
        if record["accepted"]:
            return {
                "implementation": implementation,
                "result": value,
                "attempts": attempts,
                "accepted_attempt": attempt,
            }
    raise GuardError(f"all attempts contaminated: {cell['id']} {implementation}")


def paired_summary(rounds: Sequence[Mapping[str, Any]], metric: str) -> dict[str, Any]:
    logs = []
    pairs = []
    for round_value in rounds:
        values = round_value["invocations"]
        require(len(values) == 4, "paired round does not have four invocations")
        for first_index, second_index in ((0, 1), (2, 3)):
            first = values[first_index]
            second = values[second_index]
            if first["implementation"] == "candidate":
                candidate, baseline = first, second
            else:
                baseline, candidate = first, second
            require(candidate["implementation"] == "candidate",
                    "candidate order mismatch")
            require(baseline["implementation"] != "candidate",
                    "baseline order mismatch")
            baseline_value = metric_value(baseline["result"],
                                          baseline["implementation"], metric)
            candidate_value = metric_value(candidate["result"], "candidate", metric)
            ratio = baseline_value / candidate_value
            require(math.isfinite(ratio) and ratio > 0, "invalid timing ratio")
            logs.append(math.log(ratio))
            pairs.append({"baseline_us": baseline_value,
                          "candidate_us": candidate_value,
                          "speedup": ratio})
    require(len(logs) in T95, f"unsupported paired sample count: {len(logs)}")
    mean = statistics.mean(logs)
    error = T95[len(logs)] * statistics.stdev(logs) / math.sqrt(len(logs))
    return {
        "pairs": pairs,
        "geometric_mean_speedup": math.exp(mean),
        "confidence_level": 0.95,
        "confidence_interval_lower": math.exp(mean - error),
        "confidence_interval_upper": math.exp(mean + error),
    }


def metric_value(value: Mapping[str, Any], implementation: str, metric: str) -> float:
    metrics = value["metrics"]
    if implementation == "main":
        require(metric == "amortized", "main supports only amortized comparison")
        return float(metrics["decode_including_setup"]["median_us_per_batch_call"])
    if metric == "execution":
        return float(metrics["decode_execution"]["median_us_per_batch_call"])
    if metric == "amortized":
        return float(metrics["decode_amortized_at_reuse"]
                     ["derived_median_us_per_batch_call"])
    if metric == "setup":
        return float(metrics["decode_plan_setup"]["median_us"])
    raise GuardError(f"unknown metric: {metric}")


def run_comparison(
    name: str,
    cell: Mapping[str, Any],
    baseline: str,
    binaries: Mapping[str, Path],
    options: argparse.Namespace,
    raw: Path,
) -> dict[str, Any]:
    orders = MAIN_ORDERS if baseline == "main" else ORDERS
    rounds = []
    for round_index in range(options.rounds):
        order = orders[round_index % len(orders)]
        invocations = []
        for slot, implementation in enumerate(order):
            invocation = run_one(
                binaries[implementation], implementation, cell,
                options.cpu, options.sibling,
                980000 + list(item["id"] for item in CELLS).index(cell["id"]),
                raw, f"{cell['id']}-{name}-round{round_index}-slot{slot}",
                options.iterations, options.warmup, options.max_attempts,
                options.max_sibling_jiffies)
            invocations.append(invocation)
        digests = {
            json.dumps(item["result"]["workload_digests"], sort_keys=True)
            for item in invocations
        }
        require(len(digests) == 1,
                f"workload digest mismatch: {cell['id']} {name}")
        rounds.append({"round": round_index, "order": list(order),
                       "invocations": invocations})
    metrics = {"amortized": paired_summary(rounds, "amortized")}
    if baseline == "control":
        metrics["execution"] = paired_summary(rounds, "execution")
        metrics["setup"] = paired_summary(rounds, "setup")
    return {"name": name, "baseline": baseline, "rounds": rounds,
            "metrics": metrics}


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--expected-source-commit", required=True)
    parser.add_argument("--expected-source-tree", required=True)
    parser.add_argument("--expected-memory-max", required=True)
    parser.add_argument("--rounds", type=int, choices=(2, 3, 5), default=3)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--max-attempts", type=int, default=5)
    parser.add_argument("--max-sibling-jiffies", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    options = parse_args(argv)
    require(os.sched_getaffinity(0).isdisjoint({options.cpu, options.sibling}),
            "coordinator affinity overlaps the measured SMT pair")
    require(options.rounds * 2 in T95, "unsupported round count")
    require(options.iterations > 0 and options.warmup >= 0,
            "invalid iteration counts")
    options.output.mkdir(mode=0o700)
    raw = options.output / "raw"
    artifacts = options.output / "artifacts"
    raw.mkdir(mode=0o700)
    artifacts.mkdir(mode=0o700)
    frozen = {
        "candidate": freeze(options.candidate, artifacts / "candidate"),
        "main": freeze(options.main, artifacts / "main"),
    }
    binaries = {
        "candidate": Path(frozen["candidate"]["snapshot"]),
        "control": Path(frozen["candidate"]["snapshot"]),
        "main": Path(frozen["main"]["snapshot"]),
    }
    cgroup_before = cgroup_state()
    require(cgroup_before["memory_max"] == options.expected_memory_max,
            "runner is not inside the expected memory-capped cgroup")
    require(cgroup_before["memory_swap_max"] == "0",
            "runner cgroup permits swap")

    # Bind the frozen executable to the clean source identity before the
    # allocation-heavy campaign.
    attestation_path = raw / "source-attestation.json"
    attestation_command = benchmark_arguments(
        binaries["candidate"], "candidate", CELLS[0], options.cpu,
        979999, attestation_path, 1, 1)
    bytes_index = attestation_command.index("--bytes") + 1
    reuse_index = attestation_command.index("--reuse") + 1
    attestation_command[bytes_index] = "64"
    attestation_command[reuse_index] = "1"
    attestation_command.append("--attest-source")
    subprocess.run(attestation_command, check=True)
    attestation = load_json(attestation_path)
    build = attestation["build"]
    require(build["source_commit"] == options.expected_source_commit and
            build["source_tree"] == options.expected_source_tree and
            not build["source_tracked_dirty"],
            "candidate source attestation mismatch")

    results = []
    for cell in CELLS:
        comparisons = [run_comparison(
            "same-executable-transform", cell, "control", binaries,
            options, raw)]
        if cell["compare_main"]:
            comparisons.append(run_comparison(
                "exact-main", cell, "main", binaries, options, raw))
        results.append({"cell": cell, "comparisons": comparisons})

    cgroup_after = cgroup_state()
    for name in ("oom", "oom_kill"):
        require(cgroup_after["memory_events"].get(name, 0) ==
                cgroup_before["memory_events"].get(name, 0),
                f"cgroup {name} counter increased")
    require(sha256(binaries["candidate"]) == frozen["candidate"]["sha256"] and
            sha256(binaries["main"]) == frozen["main"]["sha256"],
            "frozen executable changed during campaign")
    maximum_rss = {}
    for implementation in ("candidate", "control", "main"):
        values = []
        for result in results:
            for comparison in result["comparisons"]:
                for round_value in comparison["rounds"]:
                    for invocation in round_value["invocations"]:
                        if invocation["implementation"] == implementation:
                            accepted = invocation["attempts"][
                                invocation["accepted_attempt"]]
                            values.append(accepted["maximum_rss_kib"])
        if values:
            maximum_rss[implementation] = max(values)
    summary = {
        "schema": SCHEMA,
        "parameters": {
            "shard_bytes": SHARD_BYTES,
            "reuse": REUSE,
            "rounds": options.rounds,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "cpu": options.cpu,
            "sibling": options.sibling,
            "coordinator_affinity": sorted(os.sched_getaffinity(0)),
            "max_sibling_nonidle_jiffies": options.max_sibling_jiffies,
        },
        "source_attestation": build,
        "frozen_binaries": frozen,
        "cgroup_before": cgroup_before,
        "cgroup_after": cgroup_after,
        "maximum_child_rss_kib": maximum_rss,
        "results": results,
    }
    write_json(options.output / "summary.json", summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except (GuardError, OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"large-L1 guard failed: {error}", file=sys.stderr)
        raise SystemExit(1)
