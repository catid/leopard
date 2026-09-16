#!/usr/bin/env python3
"""One preregistered native release encode run; no resume, CPU search or retry."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import time

from check_encode import require, sha256, validate, BASELINE
from frozen_inputs import FrozenInputs, TIMING_ARTIFACTS
from timing_contract import (CANDIDATE, CELLS, GROUPS, SAMPLES, ROUNDS, MIN_GROUP_NS,
                             COMPARISONS, schedule, validate_timing, analyze)

SOURCE_NAMES = ("run_timing.py", "timing_contract.py", "check_encode.py", "frozen_inputs.py")
RUN_ARTIFACTS = TIMING_ARTIFACTS | {"expected.json", "qualification.json"}
RUNTIME_NAMES = frozenset({
    "/lib/x86_64-linux-gnu/libstdc++.so.6", "/lib/x86_64-linux-gnu/libgomp.so.1",
    "/lib/x86_64-linux-gnu/libgcc_s.so.1", "/lib/x86_64-linux-gnu/libc.so.6",
    "/lib/x86_64-linux-gnu/libm.so.6", "/lib64/ld-linux-x86-64.so.2",
    "/usr/bin/taskset", "/usr/bin/prlimit", "/usr/bin/python3.12"})


def host_identity():
    fields = dict(line.split(":", 1) for line in
                  Path("/proc/cpuinfo").read_text().split("\n\n")[0].splitlines() if ":" in line)
    fields = {key.strip(): value.strip() for key, value in fields.items()}
    return {**{key: fields[key] for key in
               ("vendor_id", "cpu family", "model", "model name", "microcode")},
            "hostname": os.uname().nodename, "kernel": os.uname().release}


def cpu_ticks(cpu):
    for line in Path("/proc/stat").read_text().splitlines():
        words = line.split()
        if words and words[0] == f"cpu{cpu}":
            values = [int(value) for value in words[1:]]
            require(len(values) >= 8, "short CPU counters")
            return sum(values[i] for i in (0, 1, 2, 5, 6, 7))
    raise ValueError("CPU counter missing")


def resource_snapshot():
    group = [line[3:] for line in Path("/proc/self/cgroup").read_text().splitlines()
             if line.startswith("0::")]
    require(len(group) == 1, "unified memory cgroup required")
    root = Path("/sys/fs/cgroup") / group[0].lstrip("/")
    events = dict(line.split() for line in (root / "memory.events").read_text().splitlines())
    result = {"cgroup": group[0], "events": {k: int(v) for k, v in events.items()}}
    for name in ("memory.max", "memory.peak", "memory.swap.max", "memory.swap.current"):
        result[name] = int((root / name).read_text())
    return result


def valid_resources(value):
    require(all(type(value[key]) is int and value[key] >= 0 for key in
                ("memory.max", "memory.peak", "memory.swap.max", "memory.swap.current")) and
            set(value["events"]) == {"low", "high", "max", "oom", "oom_kill", "oom_group_kill"} and
            all(type(number) is int for number in value["events"].values()), "resource types/inventory")
    require(value["memory.max"] == 268435456 and value["memory.peak"] <= 268435456 and
            value["memory.swap.max"] == value["memory.swap.current"] == 0 and
            all(number == 0 for number in value["events"].values()),
            "runtime resource gate failed")


def validate_plan(plan):
    require(plan["schema"] == "native-release-encode-plan/v1", "plan schema")
    exact = {"candidate_commit": CANDIDATE, "baseline_commit": BASELINE,
             "cpu": 52, "sibling": 116, "controller_cpu": 0, "passive_seconds": 10,
             "rounds": ROUNDS, "samples": SAMPLES, "attempt_budget": 1,
             "groups": list(GROUPS), "minimum_group_ns": MIN_GROUP_NS,
             "control_factor": 1.02, "warmup_calls": 4, "initial_checks": 2,
             "child_cpu_seconds": 180, "child_timeout_seconds": 180,
             "child_address_space_bytes": 268435456, "child_file_bytes": 1048576,
             "promotion": False}
    for key, expected in exact.items():
        require(type(plan.get(key)) is type(expected) and plan[key] == expected, "plan mismatch: " + key)
    require(json.dumps(plan["cells"]) == json.dumps([list(cell) for cell in CELLS]), "plan matrix changed")
    require(plan["comparisons"] == [[name, list(order)] for name, order in COMPARISONS],
            "plan comparison order changed")
    require(set(plan["source_sha256"]) == set(SOURCE_NAMES), "collector source inventory")
    require(set(plan["runtime_sha256"]) == RUNTIME_NAMES, "runtime inventory")
    require(set(plan["artifact_sha256"]) == RUN_ARTIFACTS, "plan artifact inventory")
    for digest in [plan["artifact_manifest_sha256"], *plan["source_sha256"].values(),
                   *plan["runtime_sha256"].values(), *plan["artifact_sha256"].values()]:
        require(type(digest) is str and re.fullmatch("[0-9a-f]{64}", digest), "invalid plan hash")
    for name in ("artifact_directory", "attempt_directory"):
        path = Path(plan[name])
        require(not path.is_absolute() and ".." not in path.parts and
                path.parts[:2] == (".research", "leopard-79h"), "unsafe campaign path")
    artifact, attempt = Path(plan["artifact_directory"]), Path(plan["attempt_directory"])
    require(not artifact.is_relative_to(attempt) and not attempt.is_relative_to(artifact),
            "attempt overlaps artifacts")


def collect(expected, launch, state, checkpoint):
    for cell in range(len(CELLS)):
        for implementation in ("main", "current"):
            label = f"check-{cell}-{implementation}"
            payload, observation = launch(implementation, cell, False, label)
            state["preflight"].append({"cell": cell, "implementation": implementation,
                                       "record": payload, **observation})
            checkpoint(state)
            validate(payload, cell, implementation, BASELINE if implementation == "main" else CANDIDATE)
            require(payload == expected[cell][implementation], "preflight does not match qualification")
    # The runtime launcher provides the single fixed passive gate at this point.
    launch(None, None, None, "passive")
    for cell, round_id, comparison, slot, implementation in schedule():
        label = f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}"
        payload, observation = launch(implementation, cell, True, label)
        state["invocations"].append({"cell": cell, "round": round_id, "comparison": comparison,
                                      "slot": slot, "implementation": implementation,
                                      "record": payload, **observation})
        checkpoint(state)
        require(observation["sibling_delta"] == 0 and observation["cpu_delta"] > 0,
                "CPU isolation gate failed")
        validate_timing(payload, cell, implementation, expected[cell][implementation])
    state["analysis"] = analyze(state["invocations"], expected)
    state["complete"] = True
    checkpoint(state)


def load_published_plan(plan_path, commit):
    require(type(commit) is str and re.fullmatch("[0-9a-f]{40}", commit), "invalid preregistration commit")
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    relative = plan_path.relative_to(root).as_posix()
    published = subprocess.check_output(["git", "show", f"{commit}:{relative}"])
    require(published == plan_path.read_bytes(), "plan is not the committed preregistration")
    subprocess.run(["git", "merge-base", "--is-ancestor", commit,
                    "origin/codex/native-release-evidence"], check=True)
    plan = json.loads(published)
    validate_plan(plan)
    return root, plan


def run(plan_path, commit):
    root, plan = load_published_plan(plan_path, commit)
    output = root / plan["attempt_directory"]
    bundle = root / plan["artifact_directory"]
    output.mkdir(mode=0o700)  # Fixed path: the sole attempt is now consumed.
    state = {"schema": "native-release-encode-attempt/v1", "preregistration_commit": commit,
             "plan_sha256": sha256(plan_path), "preflight": [], "invocations": [],
             "complete": False, "analysis": None}
    frozen = None
    locks = []

    def checkpoint(value):
        temporary = output / "attempt.tmp"
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        os.replace(temporary, output / "attempt.json")

    def verify():
        frozen.verify()
        require(host_identity() == plan["host"], "host identity changed")
        require(sha256(plan_path) == state["plan_sha256"], "plan changed")
        for name, digest in plan["source_sha256"].items():
            require(sha256(Path(__file__).parent / name) == digest, "collector source changed")
        for name, digest in plan["runtime_sha256"].items():
            require(sha256(Path(name)) == digest, "runtime file changed")
        for fd, path, identity in locks:
            info = path.lstat()
            require((info.st_dev, info.st_ino) == identity and
                    (os.fstat(fd).st_dev, os.fstat(fd).st_ino) == identity,
                    "lock identity changed")
        valid_resources(resource_snapshot())

    try:
        checkpoint(state)
        require(host_identity() == plan["host"], "host does not match preregistration")
        require({0, 52, 116} <= os.sched_getaffinity(0), "CPU pair not available to this process")
        for cpu in (52, 116):
            require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list").read_text().strip()
                    == "52,116", "physical pair changed")
        os.sched_setaffinity(0, {0})
        runtime_root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        info = runtime_root.lstat()
        require(stat.S_ISDIR(info.st_mode) and info.st_uid == os.getuid() and
                stat.S_IMODE(info.st_mode) == 0o700, "unsafe pair lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     runtime_root / f"leopard2-cpu-pair-{os.getuid()}-52-116.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            info = os.fstat(fd)
            locks.append((fd, path, (info.st_dev, info.st_ino)))
            require(stat.S_ISREG(info.st_mode) and info.st_uid == os.getuid(), "unsafe lock")
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        frozen = FrozenInputs(bundle, RUN_ARTIFACTS, plan["artifact_manifest_sha256"])
        require(frozen.pins == plan["artifact_sha256"], "plan artifact pins differ")
        verify()
        state["host"] = host_identity()
        state["artifact_directory_absolute"] = str(bundle)
        state["artifacts_before"] = dict(frozen.identities)
        state["resources_before"] = resource_snapshot()
        expected = json.loads((bundle / "expected.json").read_text())
        environment = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
                       "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1"}

        def launch(implementation, cell, measured, label):
            verify()
            if label == "passive":
                before = cpu_ticks(116)
                start = time.monotonic_ns()
                time.sleep(10)
                elapsed = time.monotonic_ns() - start
                after = cpu_ticks(116)
                state["passive"] = {"before": before, "after": after, "elapsed_ns": elapsed}
                checkpoint(state)
                require(before == after and elapsed >= 10000000000, "passive sibling gate failed")
                verify()
                return None
            command = ["/usr/bin/taskset", "-c", "52", "/usr/bin/prlimit", "--cpu=180:180",
                       "--as=268435456:268435456", "--fsize=1048576:1048576", "--",
                       frozen.executable(implementation + "-steady"),
                       "--measure" if measured else "--check", str(cell)]
            stdout, stderr = (output / (label + suffix) for suffix in (".stdout", ".stderr"))
            metadata = {"command": command, "stdout": stdout.name, "stderr": stderr.name,
                        "before": [cpu_ticks(52), cpu_ticks(116)]}
            metadata_path = output / (label + ".json")
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                child = subprocess.run(command, env=environment, stdout=out, stderr=err, timeout=180)
            metadata["after"] = [cpu_ticks(52), cpu_ticks(116)]
            metadata["exit_code"] = child.returncode
            metadata["stdout_sha256"] = sha256(stdout)
            metadata["stderr_sha256"] = sha256(stderr)
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            verify()
            require(child.returncode == 0 and stderr.stat().st_size == 0 and
                    stdout.stat().st_size <= 1048576, "benchmark child failed: " + label)
            record = json.loads(stdout.read_text())
            return record, {"sibling_delta": metadata["after"][1] - metadata["before"][1],
                            "cpu_delta": metadata["after"][0] - metadata["before"][0],
                            "metadata": metadata_path.name}

        def progress(value):
            checkpoint(value)
            count = len(value["invocations"])
            if count and count % 12 == 0:
                print(f"Retained {count}/288 timed invocations; no interim performance decision.", flush=True)

        collect(expected, launch, state, progress)
        verify()
        state["artifacts_after"] = dict(frozen.identities)
        state["resources_after"] = resource_snapshot()
    except BaseException as error:
        state["complete"] = False
        state["analysis"] = None
        state["failure"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        try:
            state["resources_final"] = resource_snapshot()
            checkpoint(state)
        finally:
            for fd, _, _ in reversed(locks):
                os.close(fd)


def replay(plan_path, output, bundle):
    plan = json.loads(plan_path.read_text())
    validate_plan(plan)
    state = json.loads((output / "attempt.json").read_text())
    require(state["schema"] == "native-release-encode-attempt/v1" and
            type(state["preregistration_commit"]) is str and
            re.fullmatch("[0-9a-f]{40}", state["preregistration_commit"]) and
            state["complete"] is True and "failure" not in state and
            state["plan_sha256"] == sha256(plan_path), "attempt is not complete")
    require(state["host"] == plan["host"] and len(state["preflight"]) == 16, "host/preflight mismatch")
    frozen = FrozenInputs(bundle, RUN_ARTIFACTS, plan["artifact_manifest_sha256"])
    require(frozen.pins == plan["artifact_sha256"], "replay artifact mismatch")
    qualified_expected = json.loads((bundle / "expected.json").read_text())
    require(state["artifacts_before"] == state["artifacts_after"] and
            set(state["artifacts_before"]) == RUN_ARTIFACTS | {"SHA256SUMS"},
            "recorded artifact metadata changed")
    original_bundle = Path(state["artifact_directory_absolute"])
    relative_bundle = Path(plan["artifact_directory"])
    require(original_bundle.parts[-len(relative_bundle.parts):] == relative_bundle.parts,
            "recorded artifact directory differs from plan")
    for name in ("resources_before", "resources_after", "resources_final"):
        valid_resources(state[name])
    passive = state["passive"]
    require(passive["before"] == passive["after"] and passive["elapsed_ns"] >= 10000000000,
            "passive interval failed")
    expected = [{} for _ in CELLS]
    for index, entry in enumerate(state["preflight"]):
        cell, implementation = index // 2, ("main", "current")[index % 2]
        require((entry["cell"], entry["implementation"]) == (cell, implementation), "preflight order")
        validate(entry["record"], cell, implementation, BASELINE if implementation == "main" else CANDIDATE)
        require(entry["record"] == qualified_expected[cell][implementation], "qualified preflight mismatch")
        expected[cell][implementation] = entry["record"]
    for index, entry in enumerate(state["preflight"] + state["invocations"]):
        label = (f"check-{entry['cell']}-{entry['implementation']}" if index < 16 else
                 f"cell-{entry['cell']}-round-{entry['round']}-{entry['comparison']}-slot-{entry['slot']}")
        require(entry["metadata"] == label + ".json", "observation label mismatch")
        metadata = json.loads((output / entry["metadata"]).read_text())
        require(metadata["stdout"] == label + ".stdout" and metadata["stderr"] == label + ".stderr",
                "raw output label mismatch")
        command = ["/usr/bin/taskset", "-c", "52", "/usr/bin/prlimit", "--cpu=180:180",
                   "--as=268435456:268435456", "--fsize=1048576:1048576", "--",
                   str(original_bundle / (entry["implementation"] + "-steady")),
                   "--check" if index < 16 else "--measure", str(entry["cell"])]
        require(metadata["command"] == command, "recorded child command mismatch")
        stdout, stderr = output / metadata["stdout"], output / metadata["stderr"]
        require(metadata["exit_code"] == 0 and stderr.stat().st_size == 0 and
                sha256(stdout) == metadata["stdout_sha256"] and
                sha256(stderr) == metadata["stderr_sha256"] and
                json.loads(stdout.read_text()) == entry["record"], "raw child mismatch")
        require(entry["sibling_delta"] == metadata["after"][1] - metadata["before"][1] and
                entry["cpu_delta"] == metadata["after"][0] - metadata["before"][0], "CPU observation mismatch")
        require(index < 16 or entry["cpu_delta"] > 0, "benchmark CPU did no work")
    result = analyze(state["invocations"], expected)
    require(result == state["analysis"], "analysis replay differs")
    frozen.verify()
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run", "replay"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--preregistration-commit")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bundle", type=Path)
    args = parser.parse_args()
    if args.mode == "run":
        require(args.preregistration_commit and args.output is None and args.bundle is None,
                "run requires commit and fixed plan paths")
        run(args.plan.resolve(), args.preregistration_commit)
    else:
        require(args.output is not None and args.bundle is not None, "replay requires retained output and bundle")
        replay(args.plan.resolve(), args.output.resolve(), args.bundle.resolve())
