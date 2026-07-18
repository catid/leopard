#!/usr/bin/env python3
"""Collect source-bound same-binary balanced forced-decoder ABBA evidence.

The matrix is external and versioned.  Every row compares two *explicitly*
forced modes selected from generic, materialized Algorithm 5, and tiled
Algorithm 5.  AUTO is never a role.  Each invocation retains a raw benchmark
JSON plus a signed record envelope that binds source, object, archive,
executable, matrix, argv, runtime, CPU lease, and strict sibling-idle evidence.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import balanced_evidence_common as common  # noqa: E402
sys.path.pop(0)


HEAVY_PROCESSES = {"ninja", "ctest", "cc1", "cc1plus", "bench_leopard2"}


def fail(message: str) -> None:
    raise common.EvidenceError(message)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    document = dict(value)
    document.pop("content_sha256", None)
    document["content_sha256"] = common.canonical_sha256(document)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(document, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def artifact(path: Path, label: str, empty: bool | None = None) -> dict[str, Any]:
    common.require(path.is_file(), f"missing {label}: {path}")
    size = path.stat().st_size
    if empty is True:
        common.require(size == 0, f"{label} is not empty")
    if empty is False:
        common.require(size > 0, f"{label} is empty")
    return {"name": path.name, "size": size, "sha256": common.sha256_file(path)}


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
                "run_balanced_abba.py" in arguments or
                "run_tiled_high_abba.py" in arguments or
                "main_compare/run_abba.py" in arguments)
            if command in HEAVY_PROCESSES or project_python:
                found.append(f"{entry.name} {command} {arguments}")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
    return found


def load_pair_lease(path: Path):
    specification = importlib.util.spec_from_file_location(
        "leopard2_balanced_pair_lease", path)
    if specification is None or specification.loader is None:
        fail(f"cannot import CPU pair lease from {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    if not hasattr(module, "PairLease"):
        fail(f"CPU pair lease module does not export PairLease: {path}")
    return module.PairLease


def output_is_ignored(path: Path, source_root: Path) -> None:
    try:
        path.relative_to(source_root)
    except ValueError:
        return
    result = subprocess.run(
        ["git", "-C", str(source_root), "check-ignore", "-q", str(path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    common.require(result.returncode == 0,
                   "output inside the source root must be ignored by Git")


def validate_raw_mode(path: Path, case: dict[str, Any], role: str) -> None:
    document, _ = common.load_json(path, "benchmark result")
    common.require(document.get("schema") == common.BENCHMARK_SCHEMA,
                   f"{path}: benchmark did not retain schema v3 path observations")
    common.require(set(document) == {
        "schema", "build", "parameters", "resolved", "correctness",
        "workload_digests", "memory", "metrics", "legacy",
    }, f"{path}: benchmark top-level shape changed")
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    correctness = document.get("correctness")
    common.require(isinstance(parameters, dict) and isinstance(resolved, dict) and
                   isinstance(correctness, dict), f"{path}: incomplete benchmark identity")
    common.require(set(parameters) == {
        "K", "R", "requested_profile", "requested_field", "requested_backend",
        "force_generic_decode", "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples",
        "report_decode_path", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations", "warmup",
        "thread_count", "seed",
    }, f"{path}: benchmark parameter shape changed")
    mode = common.role_mode(case, role)
    for key, expected in common.MODE_PARAMETERS[mode].items():
        actual = parameters.get(key)
        common.require(type(actual) is bool and actual is expected,
                       f"{path}: output does not prove exact forced {mode}: {key}")
    expected = {
        "K": case["K"], "R": case["R"],
        "requested_profile": "legacy_high_v1", "requested_field": "gf8",
        "requested_backend": case["backend"], "skip_legacy": True,
        "retain_samples": True, "report_decode_path": True,
        "shard_bytes": case["shard_bytes"],
        "loss_count": case["loss_count"], "batch": case["batch"],
        "reuse": case["reuse"], "iterations": case["iterations"],
        "warmup": case["warmup"], "thread_count": 1, "seed": case["seed"],
    }
    for key, wanted in expected.items():
        common.require(type(parameters.get(key)) is type(wanted) and
                       parameters.get(key) == wanted,
                       f"{path}: benchmark parameter {key} is not exact")
    missing = parameters.get("missing_original_indices")
    common.require(missing == list(range(case["K"])),
                   f"{path}: balanced full-loss coordinates are not canonical")
    resolved_expected = {
        "profile": "legacy_high_v1", "field": "gf8",
        "backend": case["backend"], "thread_count": 1,
        "parent_count": case["parent_count"], "padded_side": case["padded_side"],
        "selected_decode_path": mode,
        "selected_decode_rule": "forced_" + mode,
        "decode_required_work_slots": (
            case["parent_count"] if mode != "tiled" else
            2 * case["padded_side"] + case["loss_count"]),
        "decode_aligned_prefix_bytes": case["shard_bytes"] & ~63,
        "decode_tail_bytes": case["shard_bytes"] & 63,
        "decode_rounded_bytes": (case["shard_bytes"] + 63) & ~63,
        "decode_multi_item_batch": case["batch"] > 1,
    }
    for key, wanted in resolved_expected.items():
        common.require(type(resolved.get(key)) is type(wanted) and
                       resolved.get(key) == wanted,
                       f"{path}: resolved codec {key} is not exact")
    common.require(set(resolved) == set(resolved_expected),
                   f"{path}: resolved path shape changed")
    common.require(correctness.get("leopard2_round_trip") is True and
                   correctness.get("legacy_comparison") is None,
                   f"{path}: round-trip/skip-legacy contract failed")
    samples = document.get("metrics", {}).get("decode_execution", {}).get(
        "samples_us_per_batch_call")
    common.require(isinstance(samples, list) and len(samples) == case["iterations"] and
                   all(isinstance(item, (int, float)) and not isinstance(item, bool)
                       and item > 0 for item in samples),
                   f"{path}: raw decode observations are absent or invalid")


def envelope_prefix(case: dict[str, Any], round_index: int,
                    slot: int, role: str) -> str:
    return f"{case['name']}-round{round_index}-slot{slot}-{role}"


def validate_existing_envelope(
    path: Path, expected: dict[str, Any], raw_root: Path,
) -> dict[str, Any]:
    envelope, _ = common.load_json(path, "retained record envelope")
    digest = common.hex_digest(envelope.get("content_sha256"), "record content digest")
    unsigned = dict(envelope)
    del unsigned["content_sha256"]
    common.require(common.canonical_sha256(unsigned) == digest,
                   f"{path}: record content digest mismatch")
    for key, value in expected.items():
        common.require(envelope.get(key) == value,
                       f"{path}: retained record {key} differs")
    common.require(envelope.get("schema") == common.RECORD_SCHEMA and
                   envelope.get("status") == "complete",
                   f"{path}: record schema/status differs")
    isolation = envelope.get("isolation")
    common.require(isinstance(isolation, dict) and isolation.get("accepted") is True and
                   isolation.get("sibling_delta", {}).get("nonidle_total") == 0,
                   f"{path}: retained record lacks accepted idle-sibling evidence")
    artifacts = envelope.get("artifacts")
    common.require(isinstance(artifacts, dict) and
                   set(artifacts) == {"benchmark", "stdout", "stderr"},
                   f"{path}: retained artifact map differs")
    for key, empty in (("benchmark", False), ("stdout", True), ("stderr", True)):
        item = artifacts[key]
        common.require(isinstance(item, dict) and set(item) == {"name", "size", "sha256"},
                       f"{path}: retained {key} artifact identity differs")
        target = raw_root / item["name"]
        common.require(target.is_file() and target.stat().st_size == item["size"] and
                       common.sha256_file(target) == item["sha256"],
                       f"{path}: retained {key} bytes changed")
        common.require((item["size"] == 0) is empty,
                       f"{path}: retained {key} emptiness differs")
    return envelope


def expected_envelope_fields(
    sequence: int, case: dict[str, Any], round_index: int, slot: int, role: str,
    command: list[str], identity_digests: dict[str, str], lease_digest: str,
) -> dict[str, Any]:
    return {
        "sequence": sequence, "case": case["name"], "round": round_index,
        "slot": slot, "role": role, "mode": common.role_mode(case, role),
        "command": command, "command_sha256": common.canonical_sha256(command),
        "source_identity_sha256": identity_digests["source"],
        "build_identity_sha256": identity_digests["build"],
        "runner_identity_sha256": identity_digests["runner"],
        "common_identity_sha256": identity_digests["common"],
        "matrix_identity_sha256": identity_digests["matrix"],
        "lease_source_identity_sha256": identity_digests["lease_source"],
        "runtime_identity_sha256": identity_digests["runtime"],
        "lease_identity_sha256": lease_digest,
    }


def self_test() -> None:
    matrix = {
        "schema": common.MATRIX_SCHEMA, "name": "self-test",
        "cases": [],
    }
    pairs = (("generic", "materialized"), ("generic", "tiled"),
             ("materialized", "tiled"))
    for index, pair in enumerate(pairs):
        matrix["cases"].append({
            "name": f"pair-{index}", "K": 8, "R": 8,
            "profile": "legacy_high_v1", "field": "gf8", "backend": "scalar",
            "shard_bytes": 256, "loss_count": 8, "batch": 1, "reuse": 2,
            "iterations": 3, "warmup": 1, "seed": 100 + index,
            "control_mode": pair[0], "candidate_mode": pair[1],
        })
    name, cases = common.normalize_matrix(matrix)
    common.require(name == "self-test" and len(cases) == 3,
                   "external matrix normalization failed")
    for case in cases:
        for role in ("control", "candidate"):
            command = common.benchmark_command(
                Path("/tmp/bench"), case, Path("/tmp/result.json"), 7, role)
            mode = common.role_mode(case, role)
            for selector in common.MODE_SELECTORS[mode]:
                common.require(selector in command, "forced selector is absent")
            common.require("--json" in command and "--retain-samples" in command and
                           "--report-decode-path" in command,
                           "raw observation command contract changed")
    common.require([list(item) for item in common.ROUND_ORDERS] == [
        ["control", "candidate", "candidate", "control"],
        ["candidate", "control", "control", "candidate"],
        ["control", "candidate", "candidate", "control"],
    ], "independent clustered round order changed")
    before = {key: 10 for key in common.CPU_FIELDS}
    timed_after = dict(before)
    timed_after["user"] += 1
    sibling_after = dict(before)
    sibling_after["idle"] += 1
    common.require(common.isolation(before, timed_after, before, sibling_after)["accepted"],
                   "idle-sibling acceptance failed")
    sibling_after["system"] += 1
    common.require(not common.isolation(
        before, timed_after, before, sibling_after)["accepted"],
        "sibling interference was accepted")
    auto_case = json.loads(json.dumps(matrix))
    auto_case["cases"][0]["candidate_mode"] = "auto"
    try:
        common.normalize_matrix(auto_case)
    except common.EvidenceError:
        pass
    else:
        fail("AUTO was accepted as a forced benchmark role")
    print("balanced forced-path runner self-test passed: 3 modes, 3 clustered rounds")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root")
    parser.add_argument("--source-commit")
    parser.add_argument("--binary", default="build-release/bench_leopard2")
    parser.add_argument("--matrix")
    parser.add_argument("--output")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--sibling", type=int)
    parser.add_argument("--lease-source", default=common.LEASE_RELATIVE)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        return 0
    required = ("source_root", "source_commit", "matrix", "output", "cpu", "sibling")
    missing = ["--" + key.replace("_", "-") for key in required
               if getattr(args, key) is None]
    if missing:
        fail("missing required arguments: " + ", ".join(missing))
    common.require(args.cpu >= 0 and args.sibling >= 0 and args.cpu != args.sibling,
                   "CPU and sibling must be distinct non-negative CPUs")
    source_root = Path(args.source_root).resolve()
    output = Path(args.output).absolute()
    output_is_ignored(output, source_root)
    matrix_path = Path(args.matrix).absolute()
    matrix_document, _ = common.load_json(matrix_path, "external matrix")
    matrix_name, cases = common.normalize_matrix(matrix_document)
    source = common.git_identity(source_root, args.source_commit)
    binary = source_root / args.binary
    build_refresh_command = common.refresh_build(binary)
    source = common.git_identity(source_root, args.source_commit)
    build = common.build_identity(source_root, args.binary)
    runner = common.file_identity(Path(__file__), source_root, "runner")
    common_identity = common.file_identity(
        source_root / common.COMMON_RELATIVE, source_root, "common evidence module")
    matrix_identity = common.file_identity(matrix_path, source_root, "matrix")
    lease_path = Path(args.lease_source)
    if not lease_path.is_absolute():
        lease_path = source_root / lease_path
    lease_source = common.file_identity(lease_path, source_root, "CPU pair lease source")
    common.require(runner["relative_path"] == common.RUNNER_RELATIVE and
                   common_identity["relative_path"] == common.COMMON_RELATIVE and
                   lease_source["relative_path"] == common.LEASE_RELATIVE,
                   "runner/common/lease source is not canonical")
    affinity = set(os.sched_getaffinity(0))
    common.require(args.cpu in affinity and args.sibling in affinity,
                   "CPU pair is outside process affinity")
    housekeeping = affinity - {args.cpu, args.sibling}
    common.require(housekeeping, "no housekeeping CPU remains outside timed pair")
    runtime = common.runtime_identity(binary, args.cpu, args.sibling, affinity)
    identities = {
        "source": source, "build": build, "runner": runner,
        "common": common_identity, "matrix": matrix_identity,
        "lease_source": lease_source, "runtime": runtime,
    }
    identity_digests = {key: common.canonical_sha256(value)
                        for key, value in identities.items()}
    active = active_heavy_processes()
    common.require(not active, "heavy process active before campaign: " + repr(active))
    if output.exists():
        common.require(args.resume and output.is_dir(),
                       f"output exists; pass --resume only for an incomplete directory: {output}")
    else:
        common.require(not args.resume, "--resume requires an existing output directory")
        output.mkdir(parents=True)
    raw_root = output / "raw"
    raw_root.mkdir(exist_ok=True)
    manifest_path = output / "manifest.json"
    common.require(not manifest_path.exists(), "completed manifest already exists")
    checkpoint_path = output / "checkpoint.json"

    PairLease = load_pair_lease(Path(lease_source["path"]))
    original_affinity = set(affinity)
    os.sched_setaffinity(0, housekeeping)
    try:
        with PairLease(args.cpu, args.sibling) as lease:
            lease_identity = lease.revalidate()
            lease_digest = common.canonical_sha256(lease_identity)
            campaign = {
                "matrix_name": matrix_name,
                "cpu": args.cpu, "sibling": args.sibling,
                "original_affinity": sorted(original_affinity),
                "housekeeping_affinity": sorted(housekeeping),
                "round_orders": [list(order) for order in common.ROUND_ORDERS],
                "case_count": len(cases), "record_count": len(cases) * 12,
                "build_refresh_command": build_refresh_command,
                "lease": lease_identity, "lease_identity_sha256": lease_digest,
            }
            checkpoint_base: dict[str, Any] = {
                "schema": common.CHECKPOINT_SCHEMA, "status": "in_progress",
                **identities, "campaign": campaign,
            }
            completed: list[dict[str, Any]] = []
            if checkpoint_path.exists():
                retained, _ = common.load_json(checkpoint_path, "campaign checkpoint")
                retained_digest = common.hex_digest(
                    retained.get("content_sha256"), "checkpoint content digest")
                unsigned = dict(retained)
                del unsigned["content_sha256"]
                common.require(common.canonical_sha256(unsigned) == retained_digest,
                               "checkpoint content digest mismatch")
                retained_completed = unsigned.pop("completed_records", None)
                common.require(unsigned == checkpoint_base and
                               isinstance(retained_completed, list),
                               "checkpoint source/build/runtime/matrix identity changed")
                completed = retained_completed
            else:
                common.require(not args.resume, "resume checkpoint is absent")
                write_json(checkpoint_path, {**checkpoint_base, "completed_records": []})

            record_artifacts: list[dict[str, Any]] = []
            sequence = 0
            for case in cases:
                for round_index, order in enumerate(common.ROUND_ORDERS):
                    for slot, role in enumerate(order):
                        prefix = envelope_prefix(case, round_index, slot, role)
                        benchmark_path = raw_root / (prefix + ".benchmark.json")
                        stdout_path = raw_root / (prefix + ".stdout")
                        stderr_path = raw_root / (prefix + ".stderr")
                        envelope_path = raw_root / (prefix + ".record.json")
                        command = common.benchmark_command(
                            binary, case, benchmark_path, args.cpu, role)
                        expected = expected_envelope_fields(
                            sequence, case, round_index, slot, role, command,
                            identity_digests, lease_digest)
                        if envelope_path.exists():
                            common.require(args.resume, "record exists outside resume mode")
                            validate_existing_envelope(envelope_path, expected, raw_root)
                            validate_raw_mode(benchmark_path, case, role)
                        else:
                            common.require(not any(path.exists() for path in (
                                benchmark_path, stdout_path, stderr_path)),
                                f"partial record exists without envelope: {prefix}")
                            active = active_heavy_processes()
                            common.require(not active,
                                "heavy process appeared during campaign: " + repr(active))
                            current_lease = lease.revalidate()
                            common.require(common.canonical_sha256(current_lease) == lease_digest,
                                           "CPU pair lease changed during campaign")
                            timed_before = common.cpu_snapshot(args.cpu)
                            sibling_before = common.cpu_snapshot(args.sibling)
                            start_ns = time.time_ns()
                            completed_process = subprocess.run(
                                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            end_ns = time.time_ns()
                            timed_after = common.cpu_snapshot(args.cpu)
                            sibling_after = common.cpu_snapshot(args.sibling)
                            stdout_path.write_bytes(completed_process.stdout)
                            stderr_path.write_bytes(completed_process.stderr)
                            common.require(completed_process.returncode == 0,
                                           f"benchmark failed ({completed_process.returncode})")
                            common.require(not completed_process.stdout and
                                           not completed_process.stderr and
                                           benchmark_path.is_file(),
                                           "benchmark emitted streams or omitted JSON")
                            validate_raw_mode(benchmark_path, case, role)
                            isolation = common.isolation(
                                timed_before, timed_after, sibling_before, sibling_after)
                            envelope: dict[str, Any] = {
                                "schema": common.RECORD_SCHEMA, "status": "complete",
                                **expected, "start_time_ns": start_ns,
                                "end_time_ns": end_ns, "isolation": isolation,
                                "artifacts": {
                                    "benchmark": artifact(
                                        benchmark_path, "benchmark JSON", empty=False),
                                    "stdout": artifact(stdout_path, "stdout", empty=True),
                                    "stderr": artifact(stderr_path, "stderr", empty=True),
                                },
                            }
                            write_json(envelope_path, envelope)
                            common.require(isolation["accepted"],
                                           f"{prefix}: timed CPU/sibling isolation failed")
                        envelope_artifact = artifact(
                            envelope_path, "record envelope", empty=False)
                        record_artifacts.append(envelope_artifact)
                        sequence += 1
                        completed = record_artifacts.copy()
                        write_json(checkpoint_path, {
                            **checkpoint_base, "completed_records": completed})
            common.require(sequence == campaign["record_count"],
                           "campaign record count changed")
            final_lease = lease.revalidate()
            common.require(common.canonical_sha256(final_lease) == lease_digest,
                           "CPU pair lease changed before release")
    finally:
        os.sched_setaffinity(0, original_affinity)

    active = active_heavy_processes()
    common.require(not active, "heavy process active after campaign: " + repr(active))
    final_identities = {
        "source": common.git_identity(source_root, args.source_commit),
        "build": common.build_identity(source_root, args.binary),
        "runner": common.file_identity(Path(__file__), source_root, "runner"),
        "common": common.file_identity(
            source_root / common.COMMON_RELATIVE, source_root, "common evidence module"),
        "matrix": common.file_identity(matrix_path, source_root, "matrix"),
        "lease_source": common.file_identity(
            lease_path, source_root, "CPU pair lease source"),
        "runtime": common.runtime_identity(binary, args.cpu, args.sibling, original_affinity),
    }
    common.require(final_identities == identities,
                   "source/object/executable/matrix/runtime identity changed during campaign")
    manifest = {
        "schema": common.RUN_SCHEMA, "status": "complete", "valid": True,
        **identities, "campaign": campaign, "record_artifacts": record_artifacts,
        "final": final_identities,
    }
    write_json(manifest_path, manifest)
    checkpoint_path.unlink()
    print(manifest_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except common.EvidenceError as error:
        print(f"run_balanced_abba.py: {error}", file=sys.stderr)
        raise SystemExit(1)
