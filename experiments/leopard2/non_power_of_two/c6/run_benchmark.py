#!/usr/bin/env python3
"""Run the authoritative C6 benchmark and bind its exact executable.

This runner performs no compilation.  It expects a prebuilt forced-AVX2 C6
executable, pins the child to one allowed CPU, forces OpenMP to one thread, and
writes a manifest only after the raw result passes its immediate provenance
checks.  Compilation and memory-intensive work must be stopped while it runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


SCHEMA = "leopard2-c6-benchmark-run/v1"
RAW_SCHEMA = "leopard2-c6-cpp/v1"
CORE_SHA = "48803c06fbd7a6802b4438af60e3104895938c9d"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    args = parser.parse_args()

    allowed = sorted(os.sched_getaffinity(0))
    if args.cpu not in allowed:
        raise SystemExit(f"CPU {args.cpu} is not in allowed affinity {allowed}")
    executable = args.executable.resolve(strict=True)
    library = args.library.resolve(strict=True)
    root = Path(__file__).resolve().parents[4]
    source = root / "experiments/leopard2/non_power_of_two/c6/c6_gf256.cpp"
    result = args.result.resolve()
    command = [str(executable), "--backend", "avx2", str(result)]
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"

    def pin_child() -> None:
        os.sched_setaffinity(0, {args.cpu})

    completed = subprocess.run(
        command, cwd=root, env=environment, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, preexec_fn=pin_child,
    )
    write(args.stdout, completed.stdout)
    write(args.stderr, completed.stderr)
    if completed.returncode != 0:
        raise SystemExit(
            f"benchmark failed with status {completed.returncode}; see {args.stderr}"
        )
    try:
        raw = json.loads(result.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SystemExit(f"cannot parse benchmark result: {error}") from error
    expected = {
        "schema": RAW_SCHEMA,
        "core_git_sha": CORE_SHA,
        "source_sha256": sha256(source),
        "library_sha256": sha256(library),
        "requested_backend": "avx2",
        "runtime_backend": "avx2",
        "sanitizer_mode": "none",
        "affinity": [args.cpu],
        "omp_num_threads": "1",
        "openmp_max_threads": 1,
    }
    for key, value in expected.items():
        if raw.get(key) != value:
            raise SystemExit(f"raw benchmark {key} differs: {raw.get(key)!r} != {value!r}")
    if not isinstance(raw.get("cells"), list) or len(raw["cells"]) != 50:
        raise SystemExit("raw benchmark does not contain the declared 50-cell matrix")

    manifest = {
        "schema": SCHEMA,
        "status": "pass",
        "command": [str(args.executable), "--backend", "avx2", str(args.result)],
        "cpu": args.cpu,
        "environment": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
        "core_git_sha": CORE_SHA,
        "runner_sha256": sha256(Path(__file__)),
        "source_sha256": sha256(source),
        "executable_sha256": sha256(executable),
        "library_sha256": sha256(library),
        "result_sha256": sha256(result),
        "stdout_sha256": sha256(args.stdout),
        "stderr_sha256": sha256(args.stderr),
        "result_path": str(args.result),
        "stdout_path": str(args.stdout),
        "stderr_path": str(args.stderr),
        "cell_count": 50,
    }
    write(args.manifest,
          (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
