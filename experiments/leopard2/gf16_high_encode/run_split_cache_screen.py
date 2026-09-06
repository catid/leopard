#!/usr/bin/env python3
"""Small, single-attempt same-source filter, not exact-main qualification."""
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def identity(record):
    return {key: value for key, value in record.items() if key != "samples_ns"}


def validate(record, cell, expected, measured):
    require(record.get("schema") == "leopard2-gf16-split-screen/v1", "schema")
    for key in ("k", "r", "bytes"):
        require(record.get(key) == cell[key], f"wrong {key}")
    require(record.get("cell") == cell["id"] and
            record.get("execution_route") == cell["route"], "route/cell")
    require(identity(record) == expected, "workload or scratch changed")
    samples = record.get("samples_ns")
    require(type(samples) is list and len(samples) == (21 if measured else 0),
            "sample count")
    require(all(type(value) is int and value > 0 for value in samples),
            "nonpositive/noninteger sample")


def analyze(plan, invocations):
    require(len(invocations) == 72, "incomplete attempt")
    cells = []
    for cell in plan["cells"]:
        ratios = []
        for round_id in range(3):
            group = [item for item in invocations if item["cell"] == cell["id"]
                     and item["round"] == round_id]
            require([item["slot"] for item in group] == [0, 1, 2, 3], "slots")
            require([item["variant"] for item in group] == plan["order"], "order")
            require(all(item["sibling_delta"] == 0 for item in group),
                    "sibling contamination")
            medians = [statistics.median(item["record"]["samples_ns"])
                       for item in group]
            ratios.append(math.sqrt(medians[0] * medians[3] /
                                    (medians[1] * medians[2])))
        cells.append({"cell": cell["id"], "role": cell["role"],
                      "round_ratios": ratios,
                      "off_over_on": math.exp(statistics.mean(
                          math.log(value) for value in ratios))})
    controls_ok = all(1 / 1.02 <= cell["off_over_on"] <= 1.02 for cell in cells
                      if cell["role"] == "unchanged_control")
    target_ok = any(cell["off_over_on"] >= 1.05 and
                    min(cell["round_ratios"]) > 1 for cell in cells
                    if cell["role"] == "target")
    neighbors_ok = all(cell["off_over_on"] >= 1 / 1.02 for cell in cells
                       if cell["role"] != "target")
    decision = ("inconclusive_controls" if not controls_ok else
                "continue_to_future_qualification" if target_ok and neighbors_ok
                else "reject_for_this_screen")
    return {"cells": cells, "decision": decision,
            "production_promotion": False, "exact_leopard1_claim": False}


def sibling_ticks(cpu):
    for line in Path("/proc/stat").read_text().splitlines():
        words = line.split()
        if words and words[0] == f"cpu{cpu}":
            values = [int(word) for word in words[1:]]
            require(len(values) >= 8, "short CPU counters")
            return sum(values[index] for index in (0, 1, 2, 5, 6, 7))
    raise ValueError("missing sibling CPU")


def run(bundle, output):
    output.mkdir(mode=0o700)  # Never overwrite an attempt or resume partial data.
    plan = json.loads((bundle / "split_cache_screen_plan.json").read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    require(plan["cpu"] == 4 and plan["sibling"] == 68 and
            plan["controller_cpu"] == 0 and plan["attempt_budget"] == 1,
            "unsupported frozen plan")
    require(Path("/sys/devices/system/cpu/cpu4/topology/thread_siblings_list")
            .read_text().strip() == "4,68", "topology changed")
    os.sched_setaffinity(0, {plan["controller_cpu"]})

    def verify():
        for name, expected in pins["files"].items():
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222, "unfrozen input")
            require(digest(path) == expected, f"changed input: {name}")
        require(pins["files"]["on.a"] == plan["archive_on_sha256"] and
                pins["files"]["off.a"] == plan["archive_off_sha256"], "archive pins")
        for variant in ("on", "off"):
            require(pins["files"][variant] == plan["executables"][variant],
                    "executable does not match preregistration")

    state = {"schema": "leopard2-gf16-split-screen-attempt/v1",
             "plan_sha256": digest(bundle / "split_cache_screen_plan.json"),
             "pins": pins, "preflight": [], "invocations": [], "complete": False}
    lock_fds = []
    try:
        root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(root.is_dir() and not root.is_symlink() and
                root.stat().st_uid == os.getuid() and
                root.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        pair = root / f"leopard2-cpu-pair-{os.getuid()}-4-68.lock"
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"), pair):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            lock_fds.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        environment = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
                       "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
                       "OMP_THREAD_LIMIT": "1"}

        def invoke(variant, cell, label, measured):
            stdout = output / f"{label}.stdout"
            stderr = output / f"{label}.stderr"
            command = ["/usr/bin/taskset", "-c", "4", "/usr/bin/prlimit",
                       "--cpu=30:30", "--fsize=1048576:1048576", "--",
                       str(bundle / variant), "--measure" if measured else
                       "--check", str(cell["id"])]
            before = sibling_ticks(plan["sibling"])
            with stdout.open("xb") as out, stderr.open("xb") as err:
                result = subprocess.run(command, stdout=out, stderr=err,
                                        env=environment, timeout=60, check=False)
            after = sibling_ticks(plan["sibling"])
            require(result.returncode == 0, f"child failed: {label}")
            require(stdout.stat().st_size <= 1048576, "oversized result")
            record = json.loads(stdout.read_text())
            verify()
            return record, after - before

        # New harness parity, scratch and route checks precede all clocks.
        expected = {}
        for cell in plan["cells"]:
            for variant in ("off", "on"):
                record, _ = invoke(variant, cell,
                                   f"check-{cell['id']}-{variant}", False)
                expected.setdefault(cell["id"], identity(record))
                validate(record, cell, expected[cell["id"]], False)
                state["preflight"].append(record)
        for cell in plan["cells"]:
            for round_id in range(3):
                for slot, variant in enumerate(plan["order"]):
                    record, delta = invoke(variant, cell,
                        f"cell-{cell['id']}-round-{round_id}-slot-{slot}", True)
                    state["invocations"].append({"cell": cell["id"],
                        "round": round_id, "slot": slot, "variant": variant,
                        "sibling_delta": delta, "record": record})
                    validate(record, cell, expected[cell["id"]], True)
                    require(delta == 0, "sibling contamination; attempt stopped")
                print(f"cell {cell['id']} round {round_id} retained", flush=True)
        verify()
        state["analysis"] = analyze(plan, state["invocations"])
        state["complete"] = True
    except Exception as error:
        state["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        (output / "attempt.json").write_text(json.dumps(state, indent=2) + "\n")
        for fd in reversed(lock_fds):
            os.close(fd)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: run_split_cache_screen.py frozen_bundle output")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
