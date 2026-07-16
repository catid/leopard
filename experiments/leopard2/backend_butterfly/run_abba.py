#!/usr/bin/env python3
"""Run and verify fail-closed Leopard2 two-way-butterfly ABBA evidence."""

from __future__ import print_function

import argparse
import copy
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


SCHEMA = "leopard2-backend-butterfly-abba/v2"
SEQUENCES = (("A1", "baseline"), ("B1", "candidate"),
             ("B2", "candidate"), ("A2", "baseline"))
CELLS = (
    {
        "name": "high-gf8", "K": 240, "R": 16, "profile": "high",
        "resolved_profile": "legacy_high_v1", "field": "gf8",
        "bytes": 65536, "loss": 4,
    },
    {
        "name": "low-gf8", "K": 32, "R": 224, "profile": "low",
        "resolved_profile": "low_v1", "field": "gf8",
        "bytes": 65536, "loss": 16,
    },
    {
        "name": "high-gf16", "K": 1000, "R": 200, "profile": "high",
        "resolved_profile": "legacy_high_v1", "field": "gf16",
        "bytes": 16384, "loss": 8,
    },
    {
        "name": "low-gf16", "K": 128, "R": 1024, "profile": "low",
        "resolved_profile": "low_v1", "field": "gf16",
        "bytes": 16384, "loss": 16,
    },
)
SOURCE_FILES = (
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "LeopardFF8.cpp",
    "LeopardFF16.cpp",
    "tests/leopard2/test_backend_ops.cpp",
    "tools/check_leopard2_portable_isa.sh",
    "experiments/leopard2/backend_butterfly/run_abba.py",
)
ENTRY_PATTERN = re.compile(
    r"^(high-gf8|low-gf8|high-gf16|low-gf16)-"
    r"r([123])-(A1|B1|B2|A2)-(baseline|candidate)$"
)


class EvidenceError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise EvidenceError(message)


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path):
    return sha256_bytes(Path(path).read_bytes())


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    os.replace(str(temporary), str(path))


def parse_cpu_list(text):
    cpus = set()
    for part in text.strip().split(","):
        if not part:
            continue
        if "-" in part:
            first, last = (int(value) for value in part.split("-", 1))
            cpus.update(range(first, last + 1))
        else:
            cpus.add(int(part))
    return sorted(cpus)


def sibling_cpus(cpu):
    path = Path("/sys/devices/system/cpu/cpu{}/topology/thread_siblings_list".format(cpu))
    try:
        return parse_cpu_list(path.read_text(encoding="ascii"))
    except (OSError, ValueError) as error:
        raise EvidenceError("cannot read CPU sibling topology: {}".format(error))


def source_fingerprint(source_root):
    files = {}
    for relative in SOURCE_FILES:
        path = source_root / relative
        require(path.is_file(), "source fingerprint input missing: " + relative)
        files[relative] = sha256_file(path)
    return {"digest": sha256_bytes(canonical_bytes(files)), "files": files}


def require_provenance_identifiers(baseline_commit, matrix_fingerprint, matrix_sha):
    if baseline_commit == "self-test":
        require(matrix_fingerprint == "self-test" and matrix_sha == "self-test",
                "inconsistent self-test provenance")
        return
    require(re.fullmatch(r"[0-9a-f]{40}", baseline_commit) is not None,
            "baseline commit must be a full lowercase SHA-1")
    require(re.fullmatch(r"[0-9a-f]{64}", matrix_fingerprint) is not None,
            "matrix source fingerprint must be a lowercase SHA-256")
    require(re.fullmatch(r"[0-9a-f]{64}", matrix_sha) is not None,
            "matrix digest must be a lowercase SHA-256")


def expected_jobs():
    jobs = []
    for cell in CELLS:
        for round_number in (1, 2, 3):
            for sequence, build in SEQUENCES:
                name = "{}-r{}-{}-{}".format(
                    cell["name"], round_number, sequence, build)
                jobs.append((name, cell, round_number, sequence, build))
    return jobs


def benchmark_arguments(cell):
    return [
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--profile", cell["profile"], "--field", cell["field"],
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", "1", "--reuse", "8", "--iterations", "7",
        "--warmup", "3", "--threads", "1", "--backend", "auto",
        "--seed", "42", "--json", "-",
    ]


def check_raw(raw, cell, label):
    require(raw.get("schema") == "leopard2-benchmark-v1", label + " schema")
    parameters = raw.get("parameters", {})
    expected = {
        "K": cell["K"], "R": cell["R"],
        "requested_profile": cell["resolved_profile"],
        "requested_field": cell["field"], "shard_bytes": cell["bytes"],
        "loss_count": cell["loss"], "batch": 1, "reuse": 8,
        "iterations": 7, "warmup": 3, "thread_count": 1, "seed": 42,
    }
    for key, expected_value in expected.items():
        require(parameters.get(key) == expected_value,
                "{} parameter {} expected {!r}, got {!r}".format(
                    label, key, expected_value, parameters.get(key)))
    require(raw.get("resolved", {}).get("backend") == "avx2",
            label + " did not resolve AVX2")
    require(raw.get("correctness", {}).get("leopard2_round_trip") is True,
            label + " round trip failed")
    metrics = raw.get("metrics", {})
    encode = float(metrics["encode_execution"]["median_us_per_batch_call"])
    decode = float(metrics["decode_execution"]["median_us_per_batch_call"])
    require(encode > 0.0 and decode > 0.0, label + " invalid timing")
    return encode, decode


def summarize(entries, raw_by_name):
    summary = []
    for cell in CELLS:
        cell_result = {"name": cell["name"], "parameters": {
            "K": cell["K"], "R": cell["R"],
            "profile": cell["resolved_profile"], "field": cell["field"],
            "shard_bytes": cell["bytes"], "losses": cell["loss"],
        }}
        for metric_index, metric in ((0, "encode"), (1, "decode")):
            medians = {}
            result = {}
            for build in ("baseline", "candidate"):
                values = []
                for entry in entries:
                    if entry["cell"] == cell["name"] and entry["build"] == build:
                        values.append(raw_by_name[entry["name"]][metric_index])
                require(len(values) == 6,
                        "{} {} {} requires six invocations".format(
                            cell["name"], metric, build))
                median = statistics.median(values)
                mad = statistics.median(abs(value - median) for value in values)
                medians[build] = median
                result[build + "_median_us"] = median
                result[build + "_mad_us"] = mad
                result[build + "_minimum_us"] = min(values)
                result[build + "_maximum_us"] = max(values)
            speedup = (medians["baseline"] / medians["candidate"] - 1.0) * 100.0
            result["speedup_percent"] = speedup
            result["promotion_threshold_percent"] = 5.0
            result["promoted"] = speedup >= 5.0
            cell_result[metric] = result
        summary.append(cell_result)
    return summary


def stable_raw_digest(entries):
    evidence = []
    for entry in sorted(entries, key=lambda item: item["name"]):
        evidence.append({
            "name": entry["name"],
            "command_sha256": entry["command_sha256"],
            "binary_sha256": entry["binary_sha256"],
            "stdout_sha256": entry["stdout_sha256"],
            "stderr_sha256": entry["stderr_sha256"],
        })
    return sha256_bytes(canonical_bytes(evidence))


def validate_manifest(manifest_path, source_root, binaries=None):
    manifest_path = Path(manifest_path)
    output_directory = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest.get("schema") == SCHEMA, "unsupported manifest schema")
    require(manifest.get("status") == "passed", "campaign status is not passed")
    entries = manifest.get("entries", [])
    require(len(entries) == 48, "manifest must contain exactly 48 entries")
    expected = {name: (cell, round_number, sequence, build)
                for name, cell, round_number, sequence, build in expected_jobs()}
    seen = set()
    raw_by_name = {}
    provenance = manifest["provenance"]
    require_provenance_identifiers(
        provenance["baseline_source_commit"],
        provenance["candidate_matrix_source_fingerprint"],
        provenance["candidate_matrix_sha256"])
    execution = provenance["execution"]
    require(execution["enforced_affinity"] == [execution["benchmark_cpu"]],
            "manifest affinity is not a singleton benchmark CPU")
    require(execution["omp_environment"] == {
        "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1"},
        "OpenMP environment is not the required fixed single-thread state")
    require(execution["reserved_sibling"] in execution["thread_siblings"],
            "reserved CPU is not a topology sibling")
    require(execution["benchmark_cpu"] in execution["thread_siblings"],
            "benchmark CPU absent from sibling topology")
    require(execution["reserved_sibling"] != execution["benchmark_cpu"],
            "benchmark CPU cannot reserve itself as sibling")
    require(bool(execution["sibling_reservation_evidence"].strip()),
            "missing sibling reservation evidence")
    fingerprint = source_fingerprint(source_root)
    require(fingerprint == provenance["candidate_source"],
            "candidate source fingerprint mismatch")
    runner_path = source_root / "experiments/leopard2/backend_butterfly/run_abba.py"
    require(sha256_file(runner_path) == provenance["runner_sha256"],
            "runner hash mismatch")

    if binaries is not None:
        for build in ("baseline", "candidate"):
            path = Path(binaries[build]).resolve()
            require(path.is_file(), build + " binary missing")
            require(sha256_file(path) == provenance["binaries"][build]["sha256"],
                    build + " binary hash mismatch")

    for entry in entries:
        name = entry.get("name", "")
        require(name in expected, "unexpected or relabelled entry: " + name)
        require(name not in seen, "duplicate entry: " + name)
        seen.add(name)
        cell, round_number, sequence, build = expected[name]
        require((entry["cell"], entry["round"], entry["sequence"], entry["build"]) ==
                (cell["name"], round_number, sequence, build),
                "entry geometry mismatch: " + name)
        stdout_path = output_directory / entry["stdout_file"]
        stderr_path = output_directory / entry["stderr_file"]
        require(stdout_path.is_file() and stderr_path.is_file(),
                "raw output missing: " + name)
        stdout = stdout_path.read_bytes()
        stderr = stderr_path.read_bytes()
        require(sha256_bytes(stdout) == entry["stdout_sha256"],
                "stdout mutation: " + name)
        require(sha256_bytes(stderr) == entry["stderr_sha256"],
                "stderr mutation: " + name)
        binary = provenance["binaries"][build]
        require(entry["binary_sha256"] == binary["sha256"],
                "binary relabel: " + name)
        expected_argv = [binary["path"]] + benchmark_arguments(cell)
        require(entry["argv"] == expected_argv, "argv mismatch: " + name)
        command_record = {
            "affinity": execution["enforced_affinity"],
            "argv": expected_argv,
            "environment": execution["omp_environment"],
        }
        require(entry["command_sha256"] == sha256_bytes(canonical_bytes(command_record)),
                "command hash mismatch: " + name)
        raw = json.loads(stdout.decode("utf-8"))
        raw_by_name[name] = check_raw(raw, cell, name)
    require(seen == set(expected), "campaign geometry is incomplete")
    require(stable_raw_digest(entries) == manifest["raw_evidence_sha256"],
            "stable raw evidence digest mismatch")
    recomputed = summarize(entries, raw_by_name)
    require(canonical_bytes(recomputed) == canonical_bytes(manifest["summary"]),
            "promotion summary does not replay from raw output")
    for cell in recomputed:
        require(cell["encode"]["promoted"] and cell["decode"]["promoted"],
                cell["name"] + " misses the 5% promotion threshold")
    return manifest


def run_campaign(args, source_root):
    require(hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"),
            "Linux scheduler affinity APIs are required")
    initial_affinity = sorted(os.sched_getaffinity(0))
    require(args.cpu in initial_affinity, "benchmark CPU is not initially allowed")
    siblings = sibling_cpus(args.cpu)
    require(args.reserved_sibling in siblings,
            "reserved sibling does not share a core with benchmark CPU")
    require(args.reserved_sibling != args.cpu, "reserved sibling equals benchmark CPU")
    require(args.reservation_evidence.strip(), "reservation evidence is required")
    require_provenance_identifiers(
        args.baseline_commit, args.matrix_source_fingerprint, args.matrix_sha256)
    if args.baseline_commit != "self-test":
        completed = subprocess.run(
            ["git", "cat-file", "-e", args.baseline_commit + "^{commit}"],
            cwd=str(source_root), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False)
        require(completed.returncode == 0,
                "baseline commit is not present in the repository")
    output_directory = args.output.resolve()
    require(not output_directory.exists() or not any(output_directory.iterdir()),
            "output directory must be absent or empty")
    output_directory.mkdir(parents=True, exist_ok=True)

    binaries = {
        "baseline": Path(args.baseline).resolve(),
        "candidate": Path(args.candidate).resolve(),
    }
    expected_hashes = {
        "baseline": args.baseline_sha256,
        "candidate": args.candidate_sha256,
    }
    binary_records = {}
    for build, path in binaries.items():
        require(path.is_file(), build + " binary not found")
        digest = sha256_file(path)
        require(re.fullmatch(r"[0-9a-f]{64}", expected_hashes[build]) is not None,
                build + " expected binary hash is malformed")
        require(digest == expected_hashes[build], build + " binary hash mismatch before run")
        binary_records[build] = {"path": str(path), "sha256": digest}

    os.sched_setaffinity(0, {args.cpu})
    enforced_affinity = sorted(os.sched_getaffinity(0))
    require(enforced_affinity == [args.cpu], "failed to enforce singleton affinity")
    child_environment = os.environ.copy()
    child_environment["OMP_NUM_THREADS"] = "1"
    child_environment["OMP_DYNAMIC"] = "FALSE"
    omp_environment = {"OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1"}
    candidate_source = source_fingerprint(source_root)
    entries = []
    raw_by_name = {}
    campaign_start = time.time()
    for index, (name, cell, round_number, sequence, build) in enumerate(expected_jobs(), 1):
        require(sorted(os.sched_getaffinity(0)) == enforced_affinity,
                "runner affinity changed before " + name)
        argv = [str(binaries[build])] + benchmark_arguments(cell)
        command_record = {
            "affinity": enforced_affinity,
            "argv": argv,
            "environment": omp_environment,
        }
        print("[{}/48] {}".format(index, name), file=sys.stderr, flush=True)
        try:
            completed = subprocess.run(
                argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=child_environment, timeout=args.timeout, check=False)
        except subprocess.TimeoutExpired as error:
            raise EvidenceError("benchmark timeout: {} ({})".format(name, error))
        stdout_file = name + ".stdout.json"
        stderr_file = name + ".stderr.txt"
        (output_directory / stdout_file).write_bytes(completed.stdout)
        (output_directory / stderr_file).write_bytes(completed.stderr)
        require(completed.returncode == 0,
                "benchmark failed {} rc={}".format(name, completed.returncode))
        try:
            raw = json.loads(completed.stdout.decode("utf-8"))
        except (UnicodeError, ValueError) as error:
            raise EvidenceError("invalid benchmark JSON {}: {}".format(name, error))
        raw_by_name[name] = check_raw(raw, cell, name)
        entries.append({
            "name": name, "cell": cell["name"], "round": round_number,
            "sequence": sequence, "build": build, "argv": argv,
            "command_sha256": sha256_bytes(canonical_bytes(command_record)),
            "binary_sha256": binary_records[build]["sha256"],
            "stdout_file": stdout_file,
            "stdout_sha256": sha256_bytes(completed.stdout),
            "stderr_file": stderr_file,
            "stderr_sha256": sha256_bytes(completed.stderr),
            "returncode": completed.returncode,
        })
    campaign_end = time.time()
    require(sorted(os.sched_getaffinity(0)) == enforced_affinity,
            "runner affinity changed during campaign")
    for build, path in binaries.items():
        require(sha256_file(path) == binary_records[build]["sha256"],
                build + " binary changed during campaign")
    summary = summarize(entries, raw_by_name)
    manifest = {
        "schema": SCHEMA,
        "status": "passed",
        "provenance": {
            "runner_sha256": sha256_file(Path(__file__).resolve()),
            "baseline_source_commit": args.baseline_commit,
            "candidate_matrix_source_fingerprint": args.matrix_source_fingerprint,
            "candidate_matrix_sha256": args.matrix_sha256,
            "candidate_source": candidate_source,
            "binaries": binary_records,
            "execution": {
                "initial_affinity": initial_affinity,
                "enforced_affinity": enforced_affinity,
                "benchmark_cpu": args.cpu,
                "thread_siblings": siblings,
                "reserved_sibling": args.reserved_sibling,
                "sibling_reservation_evidence": args.reservation_evidence,
                "omp_environment": omp_environment,
                "timeout_seconds": args.timeout,
            },
        },
        "campaign": {
            "order_per_round": [item[0] for item in SEQUENCES],
            "rounds": 3, "entry_count": 48,
            "samples_per_invocation": 7, "warmups_per_invocation": 3,
            "reuse_per_sample": 8, "batch": 1, "threads": 1, "seed": 42,
            "started_unix": campaign_start, "finished_unix": campaign_end,
        },
        "entries": entries,
        "raw_evidence_sha256": stable_raw_digest(entries),
        "summary": summary,
    }
    atomic_json(output_directory / "manifest.json", manifest)
    validate_manifest(output_directory / "manifest.json", source_root, binaries)
    return manifest


def write_mock(path, factor):
    source = """#!/usr/bin/env python3
import json,sys
a=sys.argv[1:]
def v(name): return a[a.index(name)+1]
k=int(v('--k')); r=int(v('--r')); profile=v('--profile'); field=v('--field')
base=float(k+r+int(v('--loss'))+1)
factor=%.8f
out={'schema':'leopard2-benchmark-v1','parameters':{'K':k,'R':r,'requested_profile':'legacy_high_v1' if profile=='high' else 'low_v1','requested_field':field,'shard_bytes':int(v('--bytes')),'loss_count':int(v('--loss')),'batch':int(v('--batch')),'reuse':int(v('--reuse')),'iterations':int(v('--iterations')),'warmup':int(v('--warmup')),'thread_count':int(v('--threads')),'seed':int(v('--seed'))},'resolved':{'backend':'avx2'},'correctness':{'leopard2_round_trip':True},'metrics':{'encode_execution':{'median_us_per_batch_call':base*factor},'decode_execution':{'median_us_per_batch_call':base*2*factor}}}
print(json.dumps(out,sort_keys=True))
""" % factor
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def expect_failure(callback, label):
    try:
        callback()
    except (EvidenceError, KeyError, OSError, ValueError, json.JSONDecodeError):
        return
    raise EvidenceError("mutation was accepted: " + label)


def self_test(source_root):
    require(hasattr(os, "sched_getaffinity"), "affinity unavailable")
    cpu = sorted(os.sched_getaffinity(0))[0]
    siblings = sibling_cpus(cpu)
    reserved = next((item for item in siblings if item != cpu), None)
    require(reserved is not None, "self-test needs an SMT sibling")
    with tempfile.TemporaryDirectory(prefix="leo2-butterfly-runner-") as temporary:
        root = Path(temporary)
        baseline = root / "baseline.py"
        candidate = root / "candidate.py"
        write_mock(baseline, 1.0)
        write_mock(candidate, 0.8)
        args = argparse.Namespace(
            cpu=cpu, reserved_sibling=reserved,
            reservation_evidence="runner self-test reservation",
            output=root / "evidence", baseline=baseline, candidate=candidate,
            baseline_sha256=sha256_file(baseline),
            candidate_sha256=sha256_file(candidate),
            baseline_commit="self-test", matrix_source_fingerprint="self-test",
            matrix_sha256="self-test", timeout=30,
        )
        run_campaign(args, source_root)
        manifest_path = args.output / "manifest.json"
        binaries = {"baseline": baseline, "candidate": candidate}
        validate_manifest(manifest_path, source_root, binaries)

        first = args.output / "high-gf8-r1-A1-baseline.stdout.json"
        original = first.read_bytes()
        first.write_bytes(original + b" ")
        expect_failure(lambda: validate_manifest(manifest_path, source_root, binaries),
                       "raw stdout tamper")
        first.write_bytes(original)

        missing = first.with_suffix(first.suffix + ".missing")
        first.rename(missing)
        expect_failure(lambda: validate_manifest(manifest_path, source_root, binaries),
                       "missing raw output")
        missing.rename(first)

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        relabelled = copy.deepcopy(manifest)
        relabelled["entries"][0]["name"] = relabelled["entries"][1]["name"]
        relabel_path = args.output / "relabelled.json"
        atomic_json(relabel_path, relabelled)
        expect_failure(lambda: validate_manifest(relabel_path, source_root, binaries),
                       "duplicate/relabelled entry")

        removed = copy.deepcopy(manifest)
        removed["entries"].pop()
        removed_path = args.output / "removed.json"
        atomic_json(removed_path, removed)
        expect_failure(lambda: validate_manifest(removed_path, source_root, binaries),
                       "missing manifest entry")
    print("butterfly ABBA runner self-test passed: replay + raw/missing/relabel mutations")


def add_run_arguments(parser):
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--baseline-sha256", required=True)
    parser.add_argument("--candidate-sha256", required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--matrix-source-fingerprint", required=True)
    parser.add_argument("--matrix-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--reserved-sibling", type=int, required=True)
    parser.add_argument("--reservation-evidence", required=True)
    parser.add_argument("--timeout", type=int, default=120)


def main():
    source_root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")
    run_parser = subparsers.add_parser("run")
    add_run_arguments(run_parser)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--manifest", type=Path, required=True)
    verify_parser.add_argument("--baseline", type=Path, required=True)
    verify_parser.add_argument("--candidate", type=Path, required=True)
    subparsers.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "run":
            run_campaign(args, source_root)
            print("butterfly ABBA campaign passed: entries=48 metrics=8 threshold>=5%")
        elif args.command == "verify":
            validate_manifest(args.manifest, source_root, {
                "baseline": args.baseline, "candidate": args.candidate})
            print("butterfly ABBA evidence replay passed")
        elif args.command == "self-test":
            self_test(source_root)
        else:
            parser.error("a command is required")
    except (EvidenceError, KeyError, OSError, ValueError) as error:
        print("butterfly ABBA evidence failed: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
