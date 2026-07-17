#!/usr/bin/env python3
"""Run independent, affinity-pinned Experiment W libFuzzer workers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import shutil
import signal
import stat
import subprocess
import sys
import time
from dataclasses import dataclass

sys.dont_write_bytecode = True


SCHEMA = "leopard2-code-identity-fuzz-campaign-v1"
SEED_STEP = 0x9E3779B9
MAX_LOG_BYTES = 16 * 1024 * 1024
FINDING_PATTERNS = (
    b"ERROR: AddressSanitizer",
    b"SUMMARY: UndefinedBehaviorSanitizer",
    b"runtime error:",
    b"Test unit written to",
)


@dataclass
class Worker:
    index: int
    cpu: int
    seed: int
    root: pathlib.Path
    process: subprocess.Popen[bytes]
    stdout_file: object
    stderr_file: object
    started: float
    timed_out: bool = False
    rss_limit_exceeded: bool = False
    maximum_external_rss_kib: int = 0


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def require_regular_file(path: pathlib.Path, executable: bool = False) -> None:
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"not a regular file: {path}")
    if executable and info.st_mode & 0o111 == 0:
        raise ValueError(f"not executable: {path}")


def require_empty_output(path: pathlib.Path) -> None:
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_dir():
            raise ValueError(f"results path is not a real directory: {path}")
        if any(path.iterdir()):
            raise ValueError(f"results directory is not empty: {path}")
    else:
        path.mkdir(parents=True)


def corpus_snapshot(path: pathlib.Path) -> dict[str, object]:
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"corpus is not a real directory: {path}")
    entries: dict[str, int] = {}
    total_bytes = 0
    for item in sorted(path.iterdir(), key=lambda value: value.name):
        require_regular_file(item)
        data = item.read_bytes()
        digest = sha256_bytes(data)
        entries[digest] = len(data)
        total_bytes += len(data)
    if not entries:
        raise ValueError(f"corpus is empty: {path}")
    encoded = json.dumps(
        sorted(entries.items()), separators=(",", ":")
    ).encode("ascii")
    return {
        "files": len(entries),
        "bytes": total_bytes,
        "content_sha256": sha256_bytes(encoded),
        "entries": entries,
    }


def copy_corpus(snapshot: dict[str, object], destination: pathlib.Path) -> None:
    destination.mkdir()
    entries = snapshot["entries"]
    if not isinstance(entries, dict):
        raise TypeError("corpus snapshot entries are malformed")
    source = snapshot["source"]
    if not isinstance(source, pathlib.Path):
        raise TypeError("corpus snapshot source is malformed")
    by_digest: dict[str, pathlib.Path] = {}
    for item in source.iterdir():
        if item.is_file() and not item.is_symlink():
            by_digest[sha256_file(item)] = item
    for digest in sorted(entries):
        shutil.copyfile(by_digest[digest], destination / digest)


def derive_seeds(base_seed: int, count: int) -> list[int]:
    seeds: list[int] = []
    used: set[int] = set()
    for index in range(count):
        seed = (base_seed + index * SEED_STEP) & 0xFFFFFFFF
        while seed == 0 or seed in used:
            seed = (seed + 1) & 0xFFFFFFFF
        used.add(seed)
        seeds.append(seed)
    return seeds


def git_identity(root: pathlib.Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="ascii", timeout=10,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain",
         "--untracked-files=normal"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", timeout=10,
    ).stdout)
    return commit, dirty


def pin_to_cpu(cpu: int) -> None:
    os.sched_setaffinity(0, {cpu})


def read_bounded(path: pathlib.Path) -> bytes:
    size = path.stat().st_size
    if size > MAX_LOG_BYTES:
        raise ValueError(f"worker log exceeds {MAX_LOG_BYTES} bytes: {path}")
    return path.read_bytes()


def linux_process_rss_kib(pid: int) -> int | None:
    """Read a live child's resident set without trusting fuzzer accounting."""
    try:
        status = pathlib.Path(f"/proc/{pid}/status").read_text(
            encoding="ascii", errors="strict"
        )
    except FileNotFoundError:
        return None
    for line in status.splitlines():
        if line.startswith("VmRSS:"):
            fields = line.split()
            if len(fields) != 3 or fields[2] != "kB":
                raise ValueError(f"malformed VmRSS for worker pid {pid}")
            return int(fields[1])
    return None


def kill_worker_group(worker: Worker) -> None:
    try:
        os.killpg(worker.process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    worker.process.wait(timeout=5)


def exact_int(pattern: bytes, data: bytes, label: str) -> int:
    matches = re.findall(pattern, data, re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f"expected one {label}, found {len(matches)}")
    return int(matches[0])


def worker_result(worker: Worker, requested_runs: int) -> dict[str, object]:
    stdout_path = worker.root / "stdout.txt"
    stderr_path = worker.root / "stderr.txt"
    stdout = read_bounded(stdout_path)
    stderr = read_bounded(stderr_path)
    combined = stdout + b"\n" + stderr
    for finding in FINDING_PATTERNS:
        if finding in combined:
            raise ValueError(
                f"worker {worker.index} contains finding marker {finding!r}"
            )
    if worker.timed_out:
        raise ValueError(f"worker {worker.index} exceeded its timeout")
    if worker.rss_limit_exceeded:
        raise ValueError(f"worker {worker.index} exceeded its external RSS cap")
    if worker.process.returncode != 0:
        raise ValueError(
            f"worker {worker.index} failed rc={worker.process.returncode}"
        )
    logged_seed = exact_int(rb"^INFO: Seed: ([0-9]+)$", combined, "seed")
    executed = exact_int(
        rb"^stat::number_of_executed_units:[ \t]+([0-9]+)$",
        combined, "execution count",
    )
    peak_rss = exact_int(
        rb"^stat::peak_rss_mb:[ \t]+([0-9]+)$", combined, "peak RSS"
    )
    if logged_seed != worker.seed:
        raise ValueError(
            f"worker {worker.index} logged seed {logged_seed}, expected {worker.seed}"
        )
    if executed != requested_runs:
        raise ValueError(
            f"worker {worker.index} executed {executed}, expected {requested_runs}"
        )
    artifacts = list((worker.root / "artifacts").iterdir())
    if artifacts:
        raise ValueError(f"worker {worker.index} retained crash artifacts")
    coverage = [int(value) for value in re.findall(rb"\bcov: ([0-9]+)", combined)]
    features = [int(value) for value in re.findall(rb"\bft: ([0-9]+)", combined)]
    if not coverage or not features:
        raise ValueError(f"worker {worker.index} has no coverage statistics")
    corpus = corpus_snapshot(worker.root / "corpus")
    corpus.pop("entries")
    seed_input_corpus = corpus_snapshot(worker.root / "seed-inputs")
    seed_input_corpus.pop("entries")
    return {
        "index": worker.index,
        "cpu": worker.cpu,
        "seed": worker.seed,
        "returncode": worker.process.returncode,
        "executed_units": executed,
        "maximum_edge_coverage": max(coverage),
        "maximum_feature_coverage": max(features),
        "peak_rss_mb": peak_rss,
        "maximum_external_rss_kib": worker.maximum_external_rss_kib,
        "stdout_sha256": sha256_bytes(stdout),
        "stderr_sha256": sha256_bytes(stderr),
        "corpus": corpus,
        "seed_input_corpus": seed_input_corpus,
        "seed_input_list_sha256": sha256_file(
            worker.root / "seed-inputs.txt"
        ),
        "artifacts": 0,
    }


def merge_corpora(workers: list[Worker], destination: pathlib.Path,
                  base: dict[str, object] | None = None) -> dict[str, object]:
    destination.mkdir()
    sources: list[pathlib.Path] = []
    if base is not None:
        source = base.get("source")
        if not isinstance(source, pathlib.Path):
            raise TypeError("base corpus source is malformed")
        sources.append(source)
    for worker in workers:
        sources.append(worker.root / "corpus")
    for source in sources:
        for item in sorted(source.iterdir(), key=lambda value: value.name):
            require_regular_file(item)
            data = item.read_bytes()
            digest = sha256_bytes(data)
            target = destination / digest
            if not target.exists():
                target.write_bytes(data)
    snapshot = corpus_snapshot(destination)
    snapshot.pop("entries")
    return snapshot


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fuzzer", required=True, type=pathlib.Path)
    parser.add_argument("--corpus", required=True, type=pathlib.Path)
    parser.add_argument("--results", required=True, type=pathlib.Path)
    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--runs", type=int, default=50000)
    parser.add_argument("--base-seed", type=int, default=1279609156)
    parser.add_argument("--max-len", type=int, default=65535)
    parser.add_argument("--rss-limit-mb", type=int, default=8192)
    parser.add_argument("--external-rss-limit-mb", type=int, default=10240)
    parser.add_argument("--timeout-seconds", type=float, default=600.0)
    arguments = parser.parse_args()

    fuzzer = arguments.fuzzer.resolve(strict=True)
    corpus_path = arguments.corpus.resolve(strict=True)
    results_path = arguments.results
    if not results_path.is_absolute():
        results_path = pathlib.Path.cwd() / results_path
    if any(character in str(results_path) for character in (",", "\n", "\r")):
        parser.error("results path cannot contain a comma or line break")
    require_regular_file(fuzzer, executable=True)
    allowed_cpus = sorted(os.sched_getaffinity(0))
    worker_count = arguments.workers or min(128, len(allowed_cpus))
    if (worker_count <= 0 or worker_count > len(allowed_cpus) or
            arguments.runs <= 0 or arguments.max_len <= 0 or
            arguments.max_len > 65535 or arguments.rss_limit_mb <= 0 or
            arguments.external_rss_limit_mb <= 0 or
            arguments.timeout_seconds <= 0 or
            not 0 <= arguments.base_seed <= 0xFFFFFFFF):
        parser.error("invalid campaign bounds or worker count")

    base = corpus_snapshot(corpus_path)
    base["source"] = corpus_path
    seeds = derive_seeds(arguments.base_seed, worker_count)
    source_root = pathlib.Path(__file__).resolve().parents[3]
    source_commit, source_dirty = git_identity(source_root)
    if source_dirty:
        raise RuntimeError("source tree is dirty; commit campaign code first")
    require_empty_output(results_path)

    environment = os.environ.copy()
    environment["ASAN_OPTIONS"] = (
        "detect_leaks=1:halt_on_error=1:allocator_may_return_null=0"
    )
    environment["UBSAN_OPTIONS"] = "halt_on_error=1:print_stacktrace=1"
    workers: list[Worker] = []
    try:
        for index in range(worker_count):
            root = results_path / f"worker-{index:03d}"
            root.mkdir()
            copy_corpus(base, root / "seed-inputs")
            (root / "corpus").mkdir()
            (root / "artifacts").mkdir()
            seed_inputs = sorted(
                (root / "seed-inputs").iterdir(), key=lambda value: value.name
            )
            (root / "seed-inputs.txt").write_text(
                ",".join(f"seed-inputs/{item.name}" for item in seed_inputs),
                encoding="ascii",
            )
            stdout_file = (root / "stdout.txt").open("wb")
            stderr_file = (root / "stderr.txt").open("wb")
            command = [
                str(fuzzer), "corpus",
                "-seed_inputs=@seed-inputs.txt",
                "-keep_seed=1",
                "-shuffle=0",
                "-reload=0",
                f"-runs={arguments.runs}",
                f"-seed={seeds[index]}",
                f"-max_len={arguments.max_len}",
                f"-rss_limit_mb={arguments.rss_limit_mb}",
                "-artifact_prefix=artifacts/",
                "-print_final_stats=1",
            ]
            try:
                process = subprocess.Popen(
                    command, stdin=subprocess.DEVNULL,
                    stdout=stdout_file, stderr=stderr_file,
                    env=environment, cwd=root, start_new_session=True,
                    preexec_fn=lambda cpu=allowed_cpus[index]: pin_to_cpu(cpu),
                )
            except BaseException:
                stdout_file.close()
                stderr_file.close()
                raise
            workers.append(Worker(
                index=index, cpu=allowed_cpus[index], seed=seeds[index],
                root=root, process=process, stdout_file=stdout_file,
                stderr_file=stderr_file, started=time.monotonic(),
            ))

        pending = set(range(worker_count))
        while pending:
            now = time.monotonic()
            for index in tuple(pending):
                worker = workers[index]
                if worker.process.poll() is not None:
                    pending.remove(index)
                    continue
                rss_kib = linux_process_rss_kib(worker.process.pid)
                if rss_kib is not None:
                    worker.maximum_external_rss_kib = max(
                        worker.maximum_external_rss_kib, rss_kib
                    )
                if (rss_kib is not None and
                        rss_kib > arguments.external_rss_limit_mb * 1024):
                    worker.rss_limit_exceeded = True
                    kill_worker_group(worker)
                    pending.remove(index)
                elif now - worker.started > arguments.timeout_seconds:
                    worker.timed_out = True
                    kill_worker_group(worker)
                    pending.remove(index)
            if pending:
                time.sleep(0.02)
    finally:
        for worker in workers:
            if worker.process.poll() is None:
                kill_worker_group(worker)
            worker.stdout_file.close()
            worker.stderr_file.close()

    worker_results = [worker_result(worker, arguments.runs)
                      for worker in workers]
    merged = merge_corpora(
        workers, results_path / "merged-corpus", base=base
    )
    base_public = dict(base)
    base_public.pop("entries")
    base_public.pop("source")
    manifest = {
        "schema": SCHEMA,
        "source_commit": source_commit,
        "source_dirty": source_dirty,
        "fuzzer_sha256": sha256_file(fuzzer),
        "request": {
            "workers": worker_count,
            "runs_per_worker": arguments.runs,
            "base_seed": arguments.base_seed,
            "seed_step": SEED_STEP,
            "max_input_bytes": arguments.max_len,
            "rss_limit_mb": arguments.rss_limit_mb,
            "external_rss_limit_mb": arguments.external_rss_limit_mb,
            "timeout_seconds": arguments.timeout_seconds,
        },
        "allowed_cpus": allowed_cpus,
        "selected_cpus": allowed_cpus[:worker_count],
        "base_corpus": base_public,
        "workers": worker_results,
        "aggregate": {
            "total_runs": worker_count * arguments.runs,
            "distinct_seeds": len({item["seed"] for item in worker_results}),
            "maximum_edge_coverage": max(
                int(item["maximum_edge_coverage"]) for item in worker_results
            ),
            "maximum_feature_coverage": max(
                int(item["maximum_feature_coverage"]) for item in worker_results
            ),
            "maximum_reported_worker_rss_mb": max(
                int(item["peak_rss_mb"]) for item in worker_results
            ),
            "maximum_externally_observed_worker_rss_kib": max(
                int(item["maximum_external_rss_kib"])
                for item in worker_results
            ),
            "crash_or_sanitizer_artifacts": 0,
            "merged_corpus": merged,
        },
    }
    manifest_path = results_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "manifest": str(manifest_path),
        "source_commit": source_commit,
        "workers": worker_count,
        "distinct_seeds": manifest["aggregate"]["distinct_seeds"],
        "total_runs": manifest["aggregate"]["total_runs"],
        "merged_corpus_files": merged["files"],
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
