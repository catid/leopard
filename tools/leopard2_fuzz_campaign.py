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
import math
import os
import sys
import tempfile
from pathlib import Path

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

import leopard2_lab as lab  # noqa: E402


CAMPAIGN_SCHEMA = "leopard2-fuzz-campaign/v2"
AUDIT_SCHEMA = "leopard2-fuzz-campaign-audit/v2"
TARGETS = ("api", "pruned")
CAMPAIGN_NAME = "api-and-pruned-asan-ubsan"
SANITIZER_ENVIRONMENT = {
    "ASAN_OPTIONS": "abort_on_error=1:halt_on_error=1",
    "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
}


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
    executable_identities = {
        target: lab._file_identity(executable)
        for target, executable in executables.items()
    }

    jobs = []
    for target in TARGETS:
        for index in range(seeds_per_target):
            jobs.append({
                "id": "{}.{:03d}".format(target, index),
                "command": [
                    executables[target], "{seed}", str(arguments.iterations)],
                "env": dict(SANITIZER_ENVIRONMENT),
            })
    return {
        "schema": CAMPAIGN_SCHEMA,
        "root": str(Path(arguments.root).resolve()),
        "workers": worker_count,
        "base_seed": arguments.base_seed,
        "metadata": {
            "campaign": CAMPAIGN_NAME,
            "iterations_per_seed": arguments.iterations,
            "seeds_per_target": seeds_per_target,
            "targets": list(TARGETS),
            "target_executables": executable_identities,
            "sanitizer_environment": dict(SANITIZER_ENVIRONMENT),
            "timeout_seconds": arguments.timeout_seconds,
            "memory_policy": {
                "address_space_limit_mb": 0,
                "minimum_memory_mb": arguments.minimum_memory_mb,
                "rss_limit_mb": arguments.rss_limit_mb,
            },
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
    _validate_campaign_manifest(manifest)
    lab._atomic_write_json(arguments.output, manifest)
    print("wrote {} nested-safe fuzz jobs to {} (digest {})".format(
        len(manifest["jobs"]), arguments.output,
        manifest["manifest_digest"]))


def _load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise CampaignError("cannot read {}: {}".format(path, error))


def _positive_number(value, label):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0):
        raise CampaignError("{} must be a positive number".format(label))
    return value


def _validate_executable_identity(value, label):
    if (not isinstance(value, dict) or
            set(value) != {"path", "sha256", "size_bytes"} or
            not isinstance(value.get("path"), str) or
            not Path(value["path"]).is_absolute() or
            not isinstance(value.get("sha256"), str) or
            len(value["sha256"]) != 64 or
            any(character not in "0123456789abcdef"
                for character in value["sha256"]) or
            isinstance(value.get("size_bytes"), bool) or
            not isinstance(value.get("size_bytes"), int) or
            value["size_bytes"] < 0):
        raise CampaignError("{} has an invalid executable identity".format(label))
    return value


def _validate_campaign_manifest(value):
    """Bind an audit to the exact canonical API/pruned sanitizer campaign."""
    manifest = lab.validate_manifest(value)
    source_spec = manifest.get("source_spec", {})
    if source_spec.get("schema") != CAMPAIGN_SCHEMA:
        raise CampaignError("manifest is not a versioned fuzz campaign")
    metadata = source_spec.get("metadata")
    expected_metadata_fields = {
        "campaign", "iterations_per_seed", "seeds_per_target", "targets",
        "target_executables", "sanitizer_environment", "timeout_seconds",
        "memory_policy",
    }
    if (not isinstance(metadata, dict) or
            set(metadata) != expected_metadata_fields or
            metadata.get("campaign") != CAMPAIGN_NAME or
            metadata.get("targets") != list(TARGETS) or
            metadata.get("sanitizer_environment") != SANITIZER_ENVIRONMENT):
        raise CampaignError("manifest campaign metadata is not canonical")

    seeds_per_target = metadata.get("seeds_per_target")
    iterations = metadata.get("iterations_per_seed")
    _positive_int(seeds_per_target, "manifest seeds per target")
    _positive_int(iterations, "manifest iterations per seed")
    timeout_seconds = _positive_number(
        metadata.get("timeout_seconds"), "manifest timeout")
    memory_policy = metadata.get("memory_policy")
    if (not isinstance(memory_policy, dict) or set(memory_policy) != {
            "address_space_limit_mb", "minimum_memory_mb", "rss_limit_mb"}):
        raise CampaignError("manifest campaign memory policy is not canonical")
    address_space_limit_mb = memory_policy.get("address_space_limit_mb")
    minimum_memory_mb = memory_policy.get("minimum_memory_mb")
    rss_limit_mb = memory_policy.get("rss_limit_mb")
    _positive_int(address_space_limit_mb, "manifest address-space limit", True)
    _positive_int(minimum_memory_mb, "manifest minimum memory", True)
    _positive_int(rss_limit_mb, "manifest RSS limit")
    if (address_space_limit_mb != 0 or
            minimum_memory_mb > rss_limit_mb):
        raise CampaignError("manifest campaign memory policy is invalid")

    executable_identities = metadata.get("target_executables")
    if (not isinstance(executable_identities, dict) or
            set(executable_identities) != set(TARGETS)):
        raise CampaignError("manifest campaign executable roles are invalid")
    for target in TARGETS:
        _validate_executable_identity(
            executable_identities[target], "{} target".format(target))

    base_seed = manifest.get("base_seed")
    _positive_int(base_seed, "manifest base seed", True)
    topology = manifest.get("topology")
    allowed_cpus = topology.get("allowed_cpus") if isinstance(topology, dict) else None
    if (not isinstance(allowed_cpus, list) or not allowed_cpus or
            not all(isinstance(cpu, int) and not isinstance(cpu, bool) and
                    cpu >= 0 for cpu in allowed_cpus) or
            allowed_cpus != sorted(set(allowed_cpus))):
        raise CampaignError("manifest campaign CPU topology is invalid")
    expected_workers = min(len(allowed_cpus), 128)
    if (manifest.get("workers") != expected_workers or
            manifest.get("cpu_policy") != "physical-first"):
        raise CampaignError("manifest campaign worker policy is not canonical")
    try:
        cpu_order = lab._cpu_order(topology, "physical-first")
    except (KeyError, TypeError, lab.LabError) as error:
        raise CampaignError(
            "manifest campaign CPU topology is invalid: {}".format(error))
    if (not all(isinstance(cpu, int) and not isinstance(cpu, bool) and
                cpu >= 0 for cpu in cpu_order) or
            len(cpu_order) != len(allowed_cpus) or
            set(cpu_order) != set(allowed_cpus)):
        raise CampaignError("manifest campaign CPU order is invalid")

    expected_roles = {}
    for target in TARGETS:
        for index in range(seeds_per_target):
            expected_roles["{}.{:03d}".format(target, index)] = target
    jobs = manifest["jobs"]
    if [job["id"] for job in jobs] != sorted(expected_roles):
        raise CampaignError("manifest campaign job IDs are incomplete or noncanonical")
    expected_job_fields = {
        "id", "command", "cwd", "env", "timeout_seconds", "memory_mb",
        "rss_limit_mb", "minimum_memory_mb", "cpu_set", "thread_runtime",
        "seed", "executable", "job_digest",
    }
    expected_thread_runtime = {
        "schema": lab.THREAD_POLICY_SCHEMA,
        "max_threads": 1,
        "allow_internal_team": False,
        "effective_env": lab._thread_environment(1),
    }
    for position, job in enumerate(jobs):
        job_id = job["id"]
        target = expected_roles[job_id]
        identity = executable_identities[target]
        expected_cpu_set = [cpu_order[position % len(cpu_order)]]
        if (set(job) != expected_job_fields or
                job.get("command") != [
                    identity["path"], "{seed}", str(iterations)] or
                job.get("cwd") != "." or
                job.get("env") != SANITIZER_ENVIRONMENT or
                job.get("timeout_seconds") != timeout_seconds or
                job.get("memory_mb") != address_space_limit_mb or
                job.get("minimum_memory_mb") != minimum_memory_mb or
                job.get("rss_limit_mb") != rss_limit_mb or
                job.get("cpu_set") != expected_cpu_set or
                job.get("thread_runtime") != expected_thread_runtime or
                job.get("seed") != lab._stable_seed(base_seed, job_id) or
                job.get("executable") != identity):
            raise CampaignError(
                "job {} is not a canonical {} sanitizer replay".format(
                    job_id, target))
    return manifest


def _prepare_audit_destination(arguments):
    """Clear stale output before validation without deleting input evidence."""
    destination = Path(arguments.output).absolute()
    resolved_destination = destination.resolve()
    resolved_manifest = Path(arguments.manifest).absolute().resolve()
    resolved_results = Path(arguments.output_dir).absolute().resolve()
    if resolved_destination == resolved_manifest:
        raise CampaignError("audit output may not overwrite its manifest")
    if (resolved_destination == resolved_results or
            resolved_results in resolved_destination.parents):
        raise CampaignError("audit output may not overwrite per-job evidence")
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior audit output {}: {}".format(
                destination, error))
    return destination


def audit_campaign(arguments):
    destination = _prepare_audit_destination(arguments)
    manifest = _validate_campaign_manifest(_load_json(arguments.manifest))
    metadata = manifest.get("source_spec", {}).get("metadata", {})
    expected_per_target = metadata.get("seeds_per_target")
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
        _validate_campaign_manifest(first)
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

        def clone(value):
            return json.loads(json.dumps(value))

        def resign(value, job=None):
            if job is not None:
                unsigned_job = dict(job)
                unsigned_job.pop("job_digest", None)
                job["job_digest"] = lab._digest(unsigned_job)
            unsigned_manifest = dict(value)
            unsigned_manifest.pop("manifest_digest", None)
            value["manifest_digest"] = lab._digest(unsigned_manifest)
            return value

        def expect_invalid(value, label):
            try:
                _validate_campaign_manifest(value)
            except (CampaignError, lab.LabError):
                return
            raise CampaignError("{} mutation was accepted".format(label))

        wrong_schema = clone(first)
        wrong_schema["source_spec"]["schema"] = CAMPAIGN_SCHEMA + ".forged"
        expect_invalid(resign(wrong_schema), "campaign schema")

        wrong_command = clone(first)
        wrong_command_job = wrong_command["jobs"][0]
        wrong_command_job["command"] = [
            sys.executable, "-c", "import time; time.sleep(1)"]
        expect_invalid(
            resign(wrong_command, wrong_command_job), "fuzzer command")

        wrong_environment = clone(first)
        wrong_environment_job = wrong_environment["jobs"][0]
        wrong_environment_job["env"]["ASAN_OPTIONS"] = "detect_leaks=0"
        expect_invalid(
            resign(wrong_environment, wrong_environment_job),
            "sanitizer environment")

        wrong_memory = clone(first)
        wrong_memory_job = wrong_memory["jobs"][0]
        wrong_memory_job["rss_limit_mb"] += 1
        expect_invalid(
            resign(wrong_memory, wrong_memory_job), "memory policy")

        wrong_identity = clone(first)
        wrong_identity["source_spec"]["metadata"]["target_executables"][
            "api"]["sha256"] = "0" * 64
        expect_invalid(resign(wrong_identity), "executable role identity")

        forged_metadata = clone(first["source_spec"]["metadata"])
        forged_metadata["seeds_per_target"] = 1
        forged_manifest = lab.build_manifest({
            "schema": CAMPAIGN_SCHEMA,
            "root": directory,
            "base_seed": arguments.base_seed,
            "metadata": forged_metadata,
            "defaults": {
                "timeout_seconds": arguments.timeout_seconds,
                "memory_mb": 0,
                "minimum_memory_mb": arguments.minimum_memory_mb,
                "rss_limit_mb": arguments.rss_limit_mb,
                "cpu_count": 1,
                "thread_runtime": {
                    "max_threads": 1,
                    "allow_internal_team": False,
                },
            },
            "jobs": [{
                "id": "{}.{:03d}".format(target, 0),
                "command": [
                    sys.executable, "-c", "import time; time.sleep(1)"],
                "env": dict(SANITIZER_ENVIRONMENT),
            } for target in TARGETS],
        })
        expect_invalid(forged_manifest, "non-fuzzer campaign")

        invalid_manifest_path = Path(directory) / "invalid-manifest.json"
        stale_audit_path = Path(directory) / "stale-audit.json"
        lab._atomic_write_json(invalid_manifest_path, {"schema": "invalid"})
        lab._atomic_write_json(stale_audit_path, {"stale": True})
        invalid_audit_arguments = argparse.Namespace(
            manifest=str(invalid_manifest_path),
            output_dir=str(Path(directory) / "unused-results"),
            output=str(stale_audit_path),
        )
        try:
            audit_campaign(invalid_audit_arguments)
        except (CampaignError, lab.LabError):
            pass
        else:
            raise CampaignError("invalid campaign audit unexpectedly passed")
        if stale_audit_path.exists():
            raise CampaignError("failed audit retained stale output")
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
