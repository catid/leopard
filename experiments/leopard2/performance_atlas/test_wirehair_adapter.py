#!/usr/bin/env python3
"""Bounded runtime and build-evidence smoke test for the Wirehair adapter."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


SCHEMA = "leopard2-performance-atlas-wirehair-v1/v2"
REQUIRED_FLAGS = {
    "-march=x86-64",
    "-mtune=generic",
    "-mavx2",
    "-mno-avx512f",
    "-mno-avx512bw",
    "-mno-avx512vl",
}
FORBIDDEN_FLAGS = {"-march=native", "-mtune=native"}
HEX64 = re.compile(r"^[0-9a-f]{16}$")


def fail(message):
    raise RuntimeError(message)


def command_tokens(entry):
    if "arguments" in entry:
        return list(entry["arguments"])
    if "command" in entry:
        return shlex.split(entry["command"])
    fail("compile command has neither 'arguments' nor 'command'")


def canonical(path):
    return Path(path).resolve(strict=True)


def inspect_compile_commands(path, source_dir):
    commands = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(commands, list) or not commands:
        fail("compile_commands.json is empty or malformed")

    atlas_source = canonical(Path(__file__).with_name(
        "wirehair_v1_benchmark.cpp"))
    required_sources = {
        atlas_source,
        canonical(source_dir / "wirehair.cpp"),
        canonical(source_dir / "gf256.cpp"),
    }
    observed_sources = set()
    adapter_tokens = None

    for entry in commands:
        try:
            source = canonical(entry["file"])
        except (KeyError, FileNotFoundError) as error:
            fail("compile command has an invalid source path: %s" % error)
        if source not in required_sources:
            continue
        tokens = command_tokens(entry)
        missing = sorted(REQUIRED_FLAGS.difference(tokens))
        if missing:
            fail("%s compile command lacks: %s" %
                 (source.name, ", ".join(missing)))
        forbidden = sorted(flag for flag in tokens
                           if flag in FORBIDDEN_FLAGS or
                           flag.startswith("-mavx512"))
        if forbidden:
            fail("%s compile command contains forbidden ISA flags: %s" %
                 (source.name, ", ".join(forbidden)))
        observed_sources.add(source)
        if source == atlas_source:
            adapter_tokens = tokens

    missing_sources = required_sources.difference(observed_sources)
    if missing_sources:
        fail("compile_commands.json lacks atlas sources: %s" %
             ", ".join(sorted(str(path) for path in missing_sources)))
    if not any(token == "-DLEO2_ATLAS_PURE_AVX2=0"
               for token in adapter_tokens):
        fail("adapter compile command does not bind pure_avx2=false")
    if not any(token == "-DLEO2_ATLAS_WIREHAIR_V1_COMPACT_PATH=1"
               for token in adapter_tokens):
        fail("adapter compile command does not force the compact v1 path")


def host_supports_avx2():
    if not sys.platform.startswith("linux"):
        return True
    try:
        cpuinfo = Path("/proc/cpuinfo").read_text(
            encoding="utf-8", errors="replace").lower()
    except OSError:
        return True
    return " avx2 " in " " + cpuinfo.replace("\n", " ") + " "


def require_exact_keys(value, expected, label):
    actual = set(value)
    if actual != set(expected):
        fail("%s keys differ: missing=%s extra=%s" % (
            label, sorted(set(expected) - actual),
            sorted(actual - set(expected))))


def require_number(value, label, allow_zero=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        fail("%s is not numeric" % label)
    if not math.isfinite(value) or value < 0 or (not allow_zero and value == 0):
        fail("%s is not a valid nonnegative timing/rate" % label)


def check_summary(summary, label, rate_names):
    expected = {
        "median_us", "mad_us", "minimum_us", "maximum_us", "samples_us",
    }.union(rate_names)
    require_exact_keys(summary, expected, label)
    samples = summary["samples_us"]
    if not isinstance(samples, list) or len(samples) != 1:
        fail("%s must contain exactly one bounded timing sample" % label)
    for name in ("median_us", "mad_us", "minimum_us", "maximum_us"):
        require_number(summary[name], "%s.%s" % (label, name))
    require_number(samples[0], "%s.samples_us[0]" % label)
    if abs(summary["median_us"] - samples[0]) > 1e-6:
        fail("%s median does not match its sole sample" % label)
    if summary["minimum_us"] > summary["median_us"] or \
            summary["median_us"] > summary["maximum_us"]:
        fail("%s timing order is inconsistent" % label)
    for name in rate_names:
        require_number(summary[name], "%s.%s" % (label, name))


def validate_payload(payload, expected_commit):
    require_exact_keys(payload, {
        "schema", "build", "parameters", "correctness",
        "workload_digests", "decode_input", "path_semantics",
        "timing_semantics", "metrics",
    }, "payload")
    if payload["schema"] != SCHEMA:
        fail("unexpected schema: %r" % payload["schema"])

    build = payload["build"]
    if build.get("wirehair_source_commit") != expected_commit:
        fail("Wirehair source commit is not bound in runtime output")
    if build.get("pure_avx2") is not False:
        fail("Wirehair must disclose target-qualified AVX-512 helpers")
    if build.get("wirehair_abi_version") != 2:
        fail("unexpected Wirehair ABI version")
    if build.get("wire_profile_version") != 1:
        fail("unexpected Wirehair wire-profile version")
    if build.get("isa_claim") != "wirehair_v1_compact_path_avx2":
        fail("unexpected Wirehair ISA claim")
    if build.get("target_qualified_avx512_helpers_present") is not True:
        fail("target-qualified AVX-512 helper disclosure is absent")
    if build.get("wirehair_v1_wide_xor_forced_off") is not True or \
            build.get("runtime_wide_xor_enabled") is not False or \
            build.get("measured_path_avx512") is not False:
        fail("the compact v1 measurement path is not bound")
    features = build.get("active_x86_features", {})
    if features.get("avx2") is not True:
        fail("runtime did not activate AVX2")

    parameters = payload["parameters"]
    expected_parameters = {
        "K": 8, "R": 4, "shard_bytes": 64, "loss_count": 2,
        "reuse": 1, "iterations": 1, "warmup": 0,
        "seed": 8675309, "batch": 1, "thread_count": 1,
    }
    for name, expected in expected_parameters.items():
        if parameters.get(name) != expected:
            fail("parameter %s is %r, expected %r" %
                 (name, parameters.get(name), expected))
    losses = parameters.get("missing_original_indices")
    if not isinstance(losses, list) or len(losses) != 2 or \
            losses != sorted(set(losses)) or \
            any(not isinstance(index, int) or index < 0 or index >= 8
                for index in losses):
        fail("missing-original indices are malformed")
    if payload["correctness"] != {"round_trip": True}:
        fail("Wirehair round trip did not pass")

    digests = payload["workload_digests"]
    require_exact_keys(digests, {
        "algorithm", "original_data", "generated_repair",
        "recovered_originals",
    }, "workload_digests")
    if digests["algorithm"] != "fnv1a64":
        fail("unexpected digest algorithm")
    for name in ("original_data", "generated_repair", "recovered_originals"):
        if not isinstance(digests[name], str) or not HEX64.match(digests[name]):
            fail("%s is not a canonical FNV-1a digest" % name)

    decode_input = payload["decode_input"]
    require_exact_keys(decode_input, {
        "surviving_systematic_blocks", "repair_blocks_consumed",
        "extra_repair_blocks", "arrival_order",
    }, "decode_input")
    repairs = decode_input["repair_blocks_consumed"]
    if decode_input["surviving_systematic_blocks"] != 6 or \
            not isinstance(repairs, int) or repairs < 2 or \
            decode_input["extra_repair_blocks"] != repairs - 2:
        fail("decode-input accounting is inconsistent")
    if decode_input["arrival_order"] != \
            "surviving_systematic_ascending_then_repair_ascending":
        fail("decode-input arrival order is not deterministic")

    if payload["path_semantics"] != {
            "codec": "shipping_wirehair_v1",
            "threading": "single_caller_thread",
            "wide_xor": "forced_off_on_benchmark_thread",
            "avx512_target_helpers":
                "present_but_unreachable_from_measured_v1_compact_path",
            }:
        fail("Wirehair path semantics are not fully bound")
    timing_semantics = payload["timing_semantics"]
    require_exact_keys(timing_semantics, {
        "message_precode_setup", "encode_execution",
        "encode_including_setup", "decode_including_setup",
    }, "timing_semantics")
    if any(not isinstance(value, str) or not value
           for value in timing_semantics.values()):
        fail("Wirehair timing semantics contain an empty description")

    metrics = payload["metrics"]
    require_exact_keys(metrics, {
        "message_precode_setup", "encode_execution",
        "encode_including_setup", "decode_including_setup",
    }, "metrics")
    check_summary(metrics["message_precode_setup"],
                  "metrics.message_precode_setup", set())
    encode_rates = {"message_equivalent_GB_per_s", "repair_output_GB_per_s"}
    check_summary(metrics["encode_execution"],
                  "metrics.encode_execution", encode_rates)
    check_summary(metrics["encode_including_setup"],
                  "metrics.encode_including_setup", encode_rates)
    check_summary(metrics["decode_including_setup"],
                  "metrics.decode_including_setup",
                  {"received_input_GB_per_s", "repaired_output_GB_per_s"})


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--compile-commands", type=Path, required=True)
    parser.add_argument("--wirehair-source", type=Path, required=True)
    parser.add_argument("--expected-commit", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    binary = canonical(args.binary)
    compile_commands = canonical(args.compile_commands)
    source_dir = canonical(args.wirehair_source)
    inspect_compile_commands(compile_commands, source_dir)

    if not host_supports_avx2():
        print("SKIP: host does not advertise AVX2", file=sys.stderr)
        return 77

    command = [
        str(binary), "--k", "8", "--r", "4", "--bytes", "64",
        "--loss", "2", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--seed", "8675309", "--json", "-",
    ]
    environment = dict(os.environ)
    environment.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        command, check=False, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, timeout=60, env=environment)
    if completed.returncode != 0:
        fail("adapter failed (%d): %s" %
             (completed.returncode, completed.stderr.strip()))
    if completed.stderr:
        fail("adapter wrote unexpected stderr: %s" % completed.stderr.strip())
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        fail("adapter did not emit valid JSON: %s" % error)
    validate_payload(payload, args.expected_commit)
    print("Wirehair atlas adapter compile/runtime smoke passed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.SubprocessError,
            json.JSONDecodeError) as error:
        print("Wirehair atlas adapter smoke: %s" % error, file=sys.stderr)
        sys.exit(1)
