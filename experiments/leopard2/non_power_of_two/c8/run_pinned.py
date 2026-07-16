#!/usr/bin/env python3
"""Run one C8 timing child pinned to a physical CPU and audit its sibling."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


def cpu_times(cpu: int) -> dict[str, int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            values = [int(value) for value in line.split()[1:]]
            values += [0] * (10 - len(values))
            names = ("user", "nice", "system", "idle", "iowait", "irq",
                     "softirq", "steal", "guest", "guest_nice")
            return dict(zip(names, values))
    raise ValueError(f"CPU {cpu} is absent from /proc/stat")


def delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    return {name: after[name] - before[name] for name in before}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--sibling", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args()
    command = arguments.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a child command is required after --")
    before = {str(cpu): cpu_times(cpu)
              for cpu in (arguments.cpu, arguments.sibling)}
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"})
    started_ns = time.time_ns()
    completed = subprocess.run(
        ["taskset", "-c", str(arguments.cpu), *command], env=environment,
        check=False,
    )
    ended_ns = time.time_ns()
    after = {str(cpu): cpu_times(cpu)
             for cpu in (arguments.cpu, arguments.sibling)}
    deltas = {key: delta(before[key], after[key]) for key in before}
    sibling_delta = deltas[str(arguments.sibling)]
    sibling_nonidle = sum(sibling_delta[name]
                          for name in ("user", "nice", "system", "irq",
                                      "softirq", "steal"))
    value = {
        "schema": "leopard2-c8-isolation-v1",
        "command": command,
        "timing_cpu": arguments.cpu,
        "sibling_cpu": arguments.sibling,
        "started_ns": started_ns,
        "ended_ns": ended_ns,
        "elapsed_seconds": (ended_ns - started_ns) / 1e9,
        "exit_code": completed.returncode,
        "before": before,
        "after": after,
        "delta": deltas,
        "sibling_nonidle_jiffies": sibling_nonidle,
        "sibling_idle": sibling_nonidle == 0,
        "environment": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
    }
    arguments.output.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8")
    if completed.returncode:
        return completed.returncode
    if sibling_nonidle:
        raise SystemExit("sibling CPU was not idle; discard the timing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
