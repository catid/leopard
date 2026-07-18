#!/usr/bin/env python3
"""Create and audit deterministic nested-thread-safe Leopard2 fuzz campaigns.

The campaign uses ``leopard2_lab.py`` for content-addressed executables,
stable seeds, affinity, memory limits, timeouts, resumable per-job results, and
deterministic merge.  Each sanitizer replay is deliberately a one-CPU,
one-thread job; the lab runner overrides inherited native-runtime defaults and
rejects observed oversubscription before the result can be accepted.
"""

from __future__ import print_function

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

import leopard2_lab as lab  # noqa: E402


CAMPAIGN_SCHEMA = "leopard2-fuzz-campaign/v1"
AUDIT_SCHEMA = "leopard2-fuzz-campaign-audit/v1"
TARGETS = ("api", "pruned")


class CampaignError(Exception):
    pass


def _positive_int(value, label, allow_zero=False):
    if (isinstance(value, bool) or not isinstance(value, int) or
            value < 0 or (value == 0 and not allow_zero)):
        raise CampaignError("{} must be {} integer".format(
            label, "a non-negative" if allow_zero else "a positive"))
    return value


def _campaign_spec(arguments):
    topology = lab.detect_topology()
    worker_count = min(len(topology["allowed_cpus"]), 128)
    seeds_per_target = (
        arguments.seeds_per_target
        if arguments.seeds_per_target is not None else worker_count)
    _positive_int(worker_count, "visible worker count")
    _positive_int(seeds_per_target, "seeds per target")
    _positive_int(arguments.iterations, "iterations")
    _positive_int(arguments.minimum_memory_mb, "minimum memory", True)
    _positive_int(arguments.rss_limit_mb, "RSS limit")
    if arguments.minimum_memory_mb > arguments.rss_limit_mb:
        raise CampaignError("minimum memory may not exceed the RSS limit")

    executables = {
        "api": str(Path(arguments.api_executable).resolve()),
        "pruned": str(Path(arguments.pruned_executable).resolve()),
    }
    for target, executable in executables.items():
        path = Path(executable)
        if not path.is_file() or not os.access(str(path), os.X_OK):
            raise CampaignError(
                "{} fuzz executable is missing or not executable: {}".format(
                    target, executable))

    jobs = []
    for target in TARGETS:
        for index in range(seeds_per_target):
            jobs.append({
                "id": "{}.{:03d}".format(target, index),
                "command": [
                    executables[target], "{seed}", str(arguments.iterations)],
                "env": {
                    "ASAN_OPTIONS": "abort_on_error=1:halt_on_error=1",
                    "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
                },
            })
    return {
        "schema": CAMPAIGN_SCHEMA,
        "root": str(Path(arguments.root).resolve()),
        "workers": worker_count,
        "base_seed": arguments.base_seed,
        "metadata": {
            "campaign": "api-and-pruned-asan-ubsan",
            "iterations_per_seed": arguments.iterations,
            "seeds_per_target": seeds_per_target,
            "targets": list(TARGETS),
        },
        "defaults": {
            "timeout_seconds": arguments.timeout_seconds,
            # ASan reserves a very large virtual shadow mapping, so RLIMIT_AS
            # is disabled and the lab's sampled process-session RSS cap is the
            # enforceable per-job bound.
            "memory_mb": 0,
            "minimum_memory_mb": arguments.minimum_memory_mb,
            "rss_limit_mb": arguments.rss_limit_mb,
            "cpu_count": 1,
            "thread_runtime": {
                "max_threads": 1,
                "allow_internal_team": False,
            },
        },
        "jobs": jobs,
    }


def create_manifest(arguments):
    specification = _campaign_spec(arguments)
    manifest = lab.build_manifest(specification)
    lab._atomic_write_json(arguments.output, manifest)
    print("wrote {} nested-safe fuzz jobs to {} (digest {})".format(
        len(manifest["jobs"]), arguments.output,
        manifest["manifest_digest"]))


def _load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise CampaignError("cannot read {}: {}".format(path, error))


def audit_campaign(arguments):
    manifest = lab.validate_manifest(_load_json(arguments.manifest))
    metadata = manifest.get("source_spec", {}).get("metadata", {})
    if (metadata.get("campaign") != "api-and-pruned-asan-ubsan" or
            metadata.get("targets") != list(TARGETS)):
        raise CampaignError("manifest is not an API/pruned sanitizer campaign")
    expected_per_target = metadata.get("seeds_per_target")
    _positive_int(expected_per_target, "manifest seeds per target")
    counts = {target: 0 for target in TARGETS}
    seeds = set()
    for job in manifest["jobs"]:
        target = job["id"].split(".", 1)[0]
        if target not in counts:
            raise CampaignError("unexpected campaign job {}".format(job["id"]))
        counts[target] += 1
        if job["seed"] in seeds:
            raise CampaignError("campaign seeds are not distinct")
        seeds.add(job["seed"])
        policy = job["thread_runtime"]
        if (len(job["cpu_set"]) != 1 or policy["max_threads"] != 1 or
                policy["allow_internal_team"] or
                policy["effective_env"].get("OMP_NUM_THREADS") != "1" or
                policy["effective_env"].get("OMP_DYNAMIC") != "FALSE"):
            raise CampaignError(
                "job {} is not a one-CPU, one-thread workload".format(
                    job["id"]))
        if job.get("rss_limit_mb", 0) <= 0:
            raise CampaignError("job {} has no RSS bound".format(job["id"]))
    if any(count != expected_per_target for count in counts.values()):
        raise CampaignError("campaign target seed counts are incomplete")

    destination = Path(arguments.output)
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior audit output {}: {}".format(
                destination, error))
    # Keep the lab's ordinary deterministic merge separate.  The requested
    # campaign audit artifact is written only after the stricter live-sample
    # gates below pass, so a failed audit cannot leave a newly forged-looking
    # "audited" output behind.
    merged = lab.merge_results(
        manifest, arguments.output_dir, allow_missing=False)
    if merged["summary"] != {"missing": 0, "success": len(manifest["jobs"])}:
        raise CampaignError(
            "campaign has non-success results: {}".format(merged["summary"]))
    for result in merged["results"]:
        observation = result["thread_runtime"]["observation"]
        # The generic lab runner seeds one launch observation after its
        # affinity-setting pre-exec hook succeeds so sub-millisecond commands
        # can still be resumed.  Fuzz evidence is held to a stronger gate: at
        # least one subsequent /proc session sample must observe affinity and
        # resident memory while the workload is live.
        if (observation["sample_count"] < 2 or
                observation["affinity_sample_count"] < 2 or
                observation["peak_thread_count"] > 1 or
                observation["peak_process_count"] > 1 or
                observation["oversubscribed"] or
                observation["outside_cpu_set"] or
                observation["rss_exceeded"]):
            raise CampaignError(
                "job {} violated its runtime allocation".format(
                    result["job_id"]))
    lab._atomic_write_json(destination, {
        "schema": AUDIT_SCHEMA,
        "manifest_digest": manifest["manifest_digest"],
        "job_count": len(manifest["jobs"]),
        "distinct_seed_count": len(seeds),
        "summary": merged["summary"],
        "merged_results": merged,
    })
    print("audited {} jobs, {} distinct seeds, one thread per CPU".format(
        len(manifest["jobs"]), len(seeds)))


def self_test():
    with tempfile.TemporaryDirectory(
            prefix="leopard2-fuzz-campaign-self-test-") as directory:
        arguments = argparse.Namespace(
            api_executable=sys.executable,
            pruned_executable=sys.executable,
            output=str(Path(directory) / "manifest.json"),
            root=directory,
            seeds_per_target=2,
            iterations=17,
            base_seed=19,
            timeout_seconds=5.0,
            minimum_memory_mb=1,
            rss_limit_mb=2,
        )
        first = lab.build_manifest(_campaign_spec(arguments))
        second = lab.build_manifest(_campaign_spec(arguments))
        if lab._canonical_json_bytes(first) != lab._canonical_json_bytes(second):
            raise CampaignError("campaign manifest is not deterministic")
        if len(first["jobs"]) != 4 or len({
                job["seed"] for job in first["jobs"]}) != 4:
            raise CampaignError("campaign did not assign distinct stable seeds")
        for job in first["jobs"]:
            if (len(job["cpu_set"]) != 1 or
                    job["thread_runtime"]["max_threads"] != 1 or
                    job["thread_runtime"]["effective_env"][
                        "OMP_NUM_THREADS"] != "1" or
                    job["rss_limit_mb"] != 2):
                raise CampaignError("campaign job is not nested-thread safe")
    print("leopard2_fuzz_campaign self-test: PASS")


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create", help="create a signed lab manifest")
    create.add_argument("--api-executable", required=True)
    create.add_argument("--pruned-executable", required=True)
    create.add_argument("--output", required=True)
    create.add_argument("--root", default=os.getcwd())
    create.add_argument("--seeds-per-target", type=int)
    create.add_argument("--iterations", type=int, default=8192)
    create.add_argument("--base-seed", type=int, default=0x79A1501)
    create.add_argument("--timeout-seconds", type=float, default=3600.0)
    create.add_argument("--minimum-memory-mb", type=int, default=512)
    create.add_argument("--rss-limit-mb", type=int, default=2048)

    audit = subparsers.add_parser("audit", help="validate and merge a campaign")
    audit.add_argument("--manifest", required=True)
    audit.add_argument("--output-dir", required=True)
    audit.add_argument("--output", required=True)
    subparsers.add_parser("self-test", help="verify deterministic safe specs")
    return parser


def main(argv=None):
    arguments = _parser().parse_args(argv)
    try:
        if arguments.command == "create":
            create_manifest(arguments)
        elif arguments.command == "audit":
            audit_campaign(arguments)
        else:
            self_test()
        return 0
    except (CampaignError, lab.LabError) as error:
        print("leopard2_fuzz_campaign: error: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
