#!/usr/bin/env python3
"""Evaluate portable versus machine-profile FF8 circuit corpora offline.

Inputs are two isolated Release builds, their generated circuit files, and
paired quick-benchmark JSONL repetitions.  The output is a compact checked-in
summary with input checksums, source/compiled tradeoffs, workload weighting,
and packed-drift-normalized codec timing.  It never changes runtime policy.
"""

from __future__ import print_function

import argparse
import collections
import hashlib
import json
import math
import re
import statistics
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as cost_model  # noqa: E402
import inspect_ff8xor_assembly as inspector  # noqa: E402


FAMILIES = ("multiply", "fft", "ifft")
BENCHMARK_SCHEMA = "leopard.ff8xor.benchmark.jsonl.v2"
HISTORICAL_BENCHMARK_SCHEMA = "leopard.ff8xor.benchmark.jsonl.v1"
SUPPORTED_BENCHMARK_SCHEMAS = (
    HISTORICAL_BENCHMARK_SCHEMA, BENCHMARK_SCHEMA)
BENCHMARK_BACKENDS = ("packed_ff8", "ff8xor_native")
BENCHMARK_METADATA_FIELDS = (
    "schema", "compiler", "build_type", "build_flags", "cpu", "simd",
    "operating_system", "affinity", "quick", "warmups", "iterations",
    "minimum_sample_usec", "transpose_included", "measurement_order",
    "counter_backend", "pmu_requested", "cache_coloring_requested",
)
BENCHMARK_PROVENANCE_FIELDS = ("circuit_checksum",)
BENCHMARK_V2_PROVENANCE_FIELDS = (
    "circuit_cost_profile_id", "circuit_cost_profile_checksum",
)
BENCHMARK_CASE_FIELDS = (
    "schedule_id", "operation", "k", "r", "buffer_bytes", "loss_count",
)
BENCHMARK_ROW_FIELDS = BENCHMARK_CASE_FIELDS + (
    "record", "backend", "transpose_included", "measurement_order",
    "warmups", "iterations", "calls_per_sample", "median_us",
)
BENCHMARK_REPRODUCIBILITY_FIELDS = BENCHMARK_CASE_FIELDS + (
    "transpose_included", "cache_coloring_applied", "measurement_order",
    "warmups", "iterations", "locator_shift",
)
BENCHMARK_V2_TRAFFIC_FIELDS = (
    "modeled_payload_bytes_scheduled", "modeled_payload_bytes_elided",
    "modeled_payload_bytes_adjusted",
    "materialization_deferred_zero_fills",
    "materialization_added_zero_fills",
    "materialization_butterflies_skipped",
    "materialization_butterflies_reduced",
    "materialization_xors_skipped",
    "materialization_xors_replaced_by_copies",
    "materialization_identity_operations_elided",
)
BENCHMARK_CROSS_BUILD_FIELDS = BENCHMARK_CASE_FIELDS + (
    "transpose_included", "cache_coloring_applied", "measurement_order",
    "warmups", "iterations",
)
AVX2_FAMILIES = tuple(family + "_avx2" for family in FAMILIES)
AVX512_FAMILIES = tuple(
    family + suffix for family in FAMILIES
    for suffix in ("_avx512vl", "_avx512zmm"))


def file_checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def geometric_mean(values):
    if not values or any(value <= 0 for value in values):
        raise ValueError("geometric mean requires positive values")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def parse_array(text, name):
    match = re.search(
        r"\b%s\[256\]\s*=\s*\{(.*?)\};" % re.escape(name), text, re.S)
    if not match:
        raise ValueError("generated circuit array is missing: %s" % name)
    values = [int(value) for value in re.findall(r"\d+", match.group(1))]
    if len(values) != 256:
        raise ValueError("generated circuit array has the wrong length")
    return values


def parse_circuits(text, family):
    width = 8 if family == "multiply" else 16
    class_name = {"multiply": "Multiply", "fft": "FFT", "ifft": "IFFT"}[family]
    pattern = r"struct %sCircuit<(\d+)>\s*\{(.*?)\n\};" % class_name
    result = {}
    for coefficient, body in re.findall(pattern, text, re.S):
        gates = []
        gate_pattern = r"([xy])(\d) = XorValue\(\1\2, ([xy])(\d)\);"
        for destination_letter, destination_number, source_letter, \
                source_number in re.findall(gate_pattern, body):
            destination = int(destination_number)
            source = int(source_number)
            if destination_letter == "y":
                destination += 8
            if source_letter == "y":
                source += 8
            gates.append((destination, source))
        result[int(coefficient)] = tuple(gates)
    if set(result) != set(range(256)):
        raise ValueError("generated %s circuits are incomplete" % family)
    for gates in result.values():
        cost_model.validate_gates(gates, width)
    return result


def circuit_source(path, profile):
    text = path.read_text(encoding="utf-8")
    checksum_match = re.search(
        r"Circuit checksum \(SHA-256\): ([0-9a-f]{64})", text)
    if not checksum_match:
        raise ValueError("generated circuit checksum is missing")
    profile_match = re.search(
        r"Selection cost profile: ([^\s]+) \(([0-9a-f]{64})\)", text)
    if not profile_match:
        raise ValueError("generated circuit cost-profile provenance is missing")
    result = {
        "file_sha256": file_checksum(path),
        "circuit_checksum_sha256": checksum_match.group(1),
        "cost_profile_id": profile_match.group(1),
        "cost_profile_checksum_sha256": profile_match.group(2),
        "families": {},
    }
    for family, prefix, width in (
            ("multiply", "Multiply", 8), ("fft", "FFT", 16),
            ("ifft", "IFFT", 16)):
        gates = parse_array(text, "k%sGateCounts" % prefix)
        depths = parse_array(text, "k%sDepths" % prefix)
        circuits = parse_circuits(text, family)
        result["families"][family] = {
            "gate_sum": sum(gates),
            "gate_min": min(gates),
            "gate_max": max(gates),
            "depth_sum": sum(depths),
            "depth_min": min(depths),
            "depth_max": max(depths),
            "circuits": circuits,
            "modeled_costs": [cost_model.score_metrics(
                cost_model.source_metrics(circuits[index], width, profile),
                profile) for index in range(256)],
        }
    return result


def compiled_functions(artifact):
    nm_tool = inspector.find_tool(None, ("nm", "llvm-nm"))
    objdump_tool = inspector.find_tool(None, ("objdump", "llvm-objdump"))
    nm_output = inspector.run_tool([
        nm_tool, "-S", "-C", "--defined-only", str(artifact)])
    sizes = inspector.parse_sizes(nm_output)
    disassembly = inspector.run_tool([
        objdump_tool, "-drwC", "--no-show-raw-insn", str(artifact)])
    functions = inspector.parse_disassembly(disassembly, sizes)
    selected = {}
    for function in functions:
        key = (function["family"], function["coefficient"])
        if (key not in selected or function["code_bytes"] >
                selected[key]["code_bytes"]):
            selected[key] = function
    summary = inspector.aggregate(list(selected.values()))
    violations = inspector.hard_violations(
        list(selected.values()), summary,
        set(FAMILIES + AVX2_FAMILIES + AVX512_FAMILIES))
    spills = inspector.spill_warnings(list(selected.values()))
    return selected, {
        "artifact_sha256": file_checksum(artifact),
        "strict_pass": not violations,
        "spill_check_pass": not spills,
        "hard_violation_count": len(violations),
        "spill_warning_count": len(spills),
        "families": {
            family: summary[family]["totals"]
            for family in sorted(summary)},
    }


def weights(profile_artifact):
    workload = profile_artifact["workload"]
    return {
        "multiply": {int(key): value for key, value in
                     workload["multiply_log_weight"].items()},
        "fft": {int(key): value for key, value in
                workload["fft_skew_weight"].items()},
        "ifft": {int(key): value for key, value in
                 workload["ifft_skew_weight"].items()},
    }


def weighted_source(source, family_weights):
    result = {}
    for family in FAMILIES:
        item = source["families"][family]
        coefficient_weights = family_weights[family]
        result[family] = {
            "gate_weighted": sum(
                len(item["circuits"][coefficient]) * weight
                for coefficient, weight in coefficient_weights.items()),
            "modeled_cost_weighted": sum(
                item["modeled_costs"][coefficient] * weight
                for coefficient, weight in coefficient_weights.items()),
        }
    return result


def weighted_compiled(functions, family_weights, suffix):
    result = {}
    for family in FAMILIES:
        coefficient_weights = family_weights[family]
        compiled_family = family + suffix
        xor2 = xor3 = code = 0
        for coefficient, weight in coefficient_weights.items():
            function = functions[(compiled_family, coefficient)]
            xor2 += function["counts"]["vector_xor2"] * weight
            xor3 += function["counts"]["vector_xor3"] * weight
            code += function["code_bytes"] * weight
        result[family] = {
            "xor2_weighted": xor2,
            "xor3_weighted": xor3,
            "code_bytes_weighted": code,
        }
    return result


def _positive_finite_number(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(value) and value > 0)


def _sha256_string(value):
    return (isinstance(value, str) and
            re.match(r"^[0-9a-f]{64}$", value) is not None)


def _next_power_of_two(value):
    result = 1
    while result < value:
        result <<= 1
    return result


def _benchmark_metadata(path, records):
    metadata_records = [record for record in records
                        if record.get("record") == "metadata"]
    if len(metadata_records) != 1:
        raise ValueError("%s must contain exactly one metadata row" % path)
    metadata = metadata_records[0]
    required_fields = BENCHMARK_METADATA_FIELDS
    if metadata.get("schema") == BENCHMARK_SCHEMA:
        required_fields += (BENCHMARK_PROVENANCE_FIELDS +
                            BENCHMARK_V2_PROVENANCE_FIELDS)
    missing = [field for field in required_fields
               if field not in metadata]
    if missing:
        raise ValueError("%s metadata is missing %s" %
                         (path, ", ".join(missing)))
    if metadata["schema"] not in SUPPORTED_BENCHMARK_SCHEMAS:
        raise ValueError("%s has an unsupported benchmark schema" % path)
    text_fields = (
        "compiler", "build_type", "build_flags", "cpu", "simd",
        "operating_system", "affinity", "measurement_order",
        "counter_backend",
    )
    if any(not isinstance(metadata[field], str) or not metadata[field]
           for field in text_fields):
        raise ValueError("%s has invalid benchmark environment metadata" %
                         path)
    if "ABBA for paired end-to-end rows" not in \
            metadata["measurement_order"]:
        raise ValueError("%s metadata does not describe paired codec rows" %
                         path)
    if (not isinstance(metadata["minimum_sample_usec"], int) or
            isinstance(metadata["minimum_sample_usec"], bool) or
            metadata["minimum_sample_usec"] < 1 or
            not isinstance(metadata["pmu_requested"], bool) or
            not isinstance(metadata["cache_coloring_requested"], bool)):
        raise ValueError("%s has invalid benchmark sampling metadata" % path)
    if metadata["schema"] == BENCHMARK_SCHEMA:
        if not _sha256_string(metadata["circuit_checksum"]):
            raise ValueError("%s has an invalid circuit checksum" % path)
        if (not isinstance(metadata["circuit_cost_profile_id"], str) or
                not metadata["circuit_cost_profile_id"]):
            raise ValueError("%s has an invalid circuit cost-profile id" % path)
        if not _sha256_string(
                metadata["circuit_cost_profile_checksum"]):
            raise ValueError(
                "%s has an invalid circuit cost-profile checksum" % path)
    if metadata["quick"] is not True:
        raise ValueError("%s is not a quick benchmark repetition" % path)
    if metadata["transpose_included"] is not False:
        raise ValueError("%s includes transpose time" % path)
    if (not isinstance(metadata["warmups"], int) or
            isinstance(metadata["warmups"], bool) or
            metadata["warmups"] < 1 or
            not isinstance(metadata["iterations"], int) or
            isinstance(metadata["iterations"], bool) or
            metadata["iterations"] < 1):
        raise ValueError("%s has invalid benchmark repetition counts" % path)
    return metadata


def _case_descriptor(record, fields):
    return tuple((field, record.get(field)) for field in fields)


def benchmark_rows(path):
    records = []
    for line_number, line in enumerate(path.read_text(
            encoding="utf-8").splitlines(), 1):
        if not line:
            continue
        try:
            record = json.loads(line)
        except ValueError as error:
            raise ValueError("%s:%d is not valid JSON: %s" %
                             (path, line_number, error))
        if not isinstance(record, dict):
            raise ValueError("%s:%d is not a JSON object" %
                             (path, line_number))
        records.append(record)
    metadata = _benchmark_metadata(path, records)
    circuit_records = [record for record in records
                       if record.get("record") == "circuit_metadata"]
    if (len(circuit_records) != 3 or
            {record.get("family") for record in circuit_records} !=
            set(FAMILIES)):
        raise ValueError(
            "%s must contain one circuit-metadata row per family" % path)
    if any(not _sha256_string(record.get("checksum"))
           for record in circuit_records):
        raise ValueError("%s has invalid family circuit metadata" % path)
    circuit_checksums = {record.get("checksum")
                         for record in circuit_records}
    if len(circuit_checksums) != 1:
        raise ValueError("%s has inconsistent circuit metadata" % path)
    metadata = dict(metadata)
    circuit_checksum = next(iter(circuit_checksums))
    if metadata["schema"] == BENCHMARK_SCHEMA:
        if metadata["circuit_checksum"] != circuit_checksum:
            raise ValueError(
                "%s metadata and family circuit checksums differ" % path)
        if any(
                not isinstance(record.get("cost_profile_id"), str) or
                not record.get("cost_profile_id") or
                not _sha256_string(record.get("cost_profile_checksum"))
                for record in circuit_records):
            raise ValueError("%s has invalid family cost-profile metadata" %
                             path)
        profile_rows = {
            (record.get("cost_profile_id"),
             record.get("cost_profile_checksum"))
            for record in circuit_records}
        expected_profile = {(
            metadata["circuit_cost_profile_id"],
            metadata["circuit_cost_profile_checksum"])}
        if profile_rows != expected_profile:
            raise ValueError(
                "%s metadata and family cost profiles differ" % path)
    else:
        # Historical v1 put the circuit checksum only on its three family
        # rows.  Normalize it into metadata so the same provenance check can
        # bind old inputs to the supplied generated source.
        metadata["circuit_checksum"] = circuit_checksum
    rows = {}
    descriptors = {}
    for record in records:
        if (record.get("record") != "benchmark" or
                record.get("backend") not in BENCHMARK_BACKENDS):
            continue
        required_fields = BENCHMARK_ROW_FIELDS + (
            BENCHMARK_V2_TRAFFIC_FIELDS
            if metadata["schema"] == BENCHMARK_SCHEMA else
            ("modeled_payload_bytes",))
        missing = [field for field in required_fields
                   if field not in record]
        if missing:
            raise ValueError("%s benchmark row is missing %s" %
                             (path, ", ".join(missing)))
        if record["transpose_included"] is not False:
            raise ValueError("%s contains a transpose-inclusive codec row" %
                             path)
        if (not isinstance(record["schedule_id"], str) or
                not record["schedule_id"]):
            raise ValueError("%s contains an empty schedule id" % path)
        if record["operation"] not in ("encode", "decode"):
            raise ValueError("%s contains an unsupported codec operation" %
                             path)
        for field in ("k", "r", "buffer_bytes", "loss_count", "warmups",
                      "iterations", "calls_per_sample"):
            if (not isinstance(record[field], int) or
                    isinstance(record[field], bool) or record[field] < 0):
                raise ValueError("%s has invalid %s in a codec row" %
                                 (path, field))
        if (record["k"] < 1 or record["r"] < 1 or
                record["buffer_bytes"] < 1 or record["warmups"] < 1 or
                record["iterations"] < 1 or record["calls_per_sample"] < 1 or
                not _positive_finite_number(record["median_us"])):
            raise ValueError("%s has invalid codec timing parameters" % path)
        m = _next_power_of_two(record["r"])
        if (record["r"] > record["k"] or m + record["k"] > 256):
            raise ValueError("%s contains a codec case outside FF8" % path)
        if record["buffer_bytes"] % 64 != 0:
            raise ValueError("%s contains a misaligned shard size" % path)
        if record["warmups"] != metadata["warmups"]:
            raise ValueError("%s codec warmups disagree with metadata" % path)
        # ABBA contributes two observations per requested measurement round.
        if record["iterations"] != metadata["iterations"] * 2:
            raise ValueError("%s codec iterations disagree with metadata" %
                             path)
        if record["measurement_order"] != "ABBA":
            raise ValueError("%s contains an unpaired codec measurement" % path)
        if record["operation"] == "encode":
            if record["loss_count"] != 0:
                raise ValueError("%s encode row has a loss count" % path)
            expected_schedule = "encode-k%d-r%d-b%d" % (
                record["k"], record["r"], record["buffer_bytes"])
        else:
            if not (1 <= record["loss_count"] <= record["r"]):
                raise ValueError("%s decode row has an invalid loss count" %
                                 path)
            expected_schedule = "decode-k%d-r%d-b%d-loss%d" % (
                record["k"], record["r"], record["buffer_bytes"],
                record["loss_count"])
        if record["schedule_id"] != expected_schedule:
            raise ValueError("%s has a noncanonical schedule id" % path)

        if (not isinstance(record.get("cache_coloring_applied"), bool) or
                (record["backend"] == "packed_ff8" and
                 record.get("locator_shift") is not None) or
                (record["backend"] == "ff8xor_native" and
                 record["operation"] == "encode" and
                 record.get("locator_shift") is not None) or
                (record["backend"] == "ff8xor_native" and
                 record["operation"] == "decode" and
                 (not isinstance(record.get("locator_shift"), int) or
                  isinstance(record.get("locator_shift"), bool) or
                  not 0 <= record["locator_shift"] < 255))):
            raise ValueError("%s has invalid codec execution metadata" % path)

        if metadata["schema"] == BENCHMARK_SCHEMA:
            traffic = [record[field]
                       for field in BENCHMARK_V2_TRAFFIC_FIELDS]
            if any(not isinstance(value, int) or isinstance(value, bool) or
                   value < 0 for value in traffic):
                raise ValueError("%s has invalid payload traffic metadata" %
                                 path)
            scheduled, elided, adjusted = traffic[:3]
            if adjusted != scheduled - elided:
                raise ValueError("%s has inconsistent payload traffic" % path)
            if (record["backend"] == "packed_ff8" and any(traffic)) or \
                    (record["backend"] == "ff8xor_native" and
                     scheduled == 0):
                raise ValueError("%s has traffic for the wrong backend" % path)
        else:
            modeled = record["modeled_payload_bytes"]
            if (not isinstance(modeled, int) or isinstance(modeled, bool) or
                    modeled < 0 or
                    (record["backend"] == "packed_ff8" and modeled != 0) or
                    (record["backend"] == "ff8xor_native" and modeled == 0)):
                raise ValueError("%s has invalid modeled payload traffic" %
                                 path)
        key = (record["schedule_id"], record["backend"])
        if key in rows:
            raise ValueError("%s contains a duplicate codec row %r" %
                             (path, key))
        rows[key] = record
        descriptor_fields = BENCHMARK_REPRODUCIBILITY_FIELDS + (
            BENCHMARK_V2_TRAFFIC_FIELDS
            if metadata["schema"] == BENCHMARK_SCHEMA else
            ("modeled_payload_bytes",))
        descriptors[key] = _case_descriptor(record, descriptor_fields)

    backend_schedules = {
        backend: {schedule for schedule, item_backend in rows
                  if item_backend == backend}
        for backend in BENCHMARK_BACKENDS}
    if not backend_schedules[BENCHMARK_BACKENDS[0]]:
        raise ValueError("%s contains no native codec benchmark rows" % path)
    if backend_schedules[BENCHMARK_BACKENDS[0]] != \
            backend_schedules[BENCHMARK_BACKENDS[1]]:
        raise ValueError(
            "%s does not contain identical packed/native schedule sets" % path)
    for schedule_id in backend_schedules[BENCHMARK_BACKENDS[0]]:
        packed = rows[(schedule_id, "packed_ff8")]
        native = rows[(schedule_id, "ff8xor_native")]
        if _case_descriptor(packed, BENCHMARK_CASE_FIELDS) != \
                _case_descriptor(native, BENCHMARK_CASE_FIELDS):
            raise ValueError("%s has mismatched backend case metadata for %s" %
                             (path, schedule_id))
    return metadata, rows, descriptors


def _validate_benchmark_provenance(path, metadata, expected):
    if expected is None:
        return
    if metadata["circuit_checksum"] != expected["circuit_checksum"]:
        raise ValueError("%s uses an unrelated generated circuit corpus" % path)
    if metadata["schema"] == BENCHMARK_SCHEMA:
        for field in BENCHMARK_V2_PROVENANCE_FIELDS:
            if metadata[field] != expected[field]:
                raise ValueError(
                    "%s uses unrelated circuit cost-profile provenance" % path)


def benchmark_summary(portable_paths, machine_paths,
                      portable_provenance=None, machine_provenance=None):
    if len(portable_paths) != len(machine_paths) or not portable_paths:
        raise ValueError("portable/machine benchmark repetitions must pair")
    if len(portable_paths) < 3:
        raise ValueError("at least three paired benchmark repetitions required")
    resolved = [path.resolve() for path in portable_paths + machine_paths]
    if len(set(resolved)) != len(resolved):
        raise ValueError("benchmark repetition inputs must be distinct files")
    by_case = collections.defaultdict(list)
    raw_by_case = collections.defaultdict(list)
    per_repetition = []
    operation_ratios = collections.defaultdict(list)
    metadata = None
    expected_portable_metadata = None
    expected_machine_metadata = None
    expected_keys = None
    expected_portable_descriptors = None
    expected_machine_descriptors = None
    for portable_path, machine_path in zip(portable_paths, machine_paths):
        portable_metadata, portable, portable_descriptors = benchmark_rows(
            portable_path)
        machine_metadata, machine, machine_descriptors = benchmark_rows(
            machine_path)
        _validate_benchmark_provenance(
            portable_path, portable_metadata, portable_provenance)
        _validate_benchmark_provenance(
            machine_path, machine_metadata, machine_provenance)
        portable_environment = _case_descriptor(
            portable_metadata, BENCHMARK_METADATA_FIELDS)
        machine_environment = _case_descriptor(
            machine_metadata, BENCHMARK_METADATA_FIELDS)
        if portable_environment != machine_environment:
            raise ValueError("benchmark environments differ across builds")
        if metadata is None:
            metadata = portable_metadata
            expected_portable_metadata = portable_metadata
            expected_machine_metadata = machine_metadata
        elif (portable_metadata != expected_portable_metadata or
              machine_metadata != expected_machine_metadata):
            raise ValueError("benchmark metadata differs between repetitions")
        if set(portable) != set(machine):
            raise ValueError(
                "portable and machine benchmark row sets differ")
        for key in portable:
            if _case_descriptor(portable[key], BENCHMARK_CROSS_BUILD_FIELDS) != \
                    _case_descriptor(machine[key],
                                     BENCHMARK_CROSS_BUILD_FIELDS):
                raise ValueError(
                    "portable and machine benchmark case metadata differs")
        if expected_keys is None:
            expected_keys = set(portable)
            expected_portable_descriptors = portable_descriptors
            expected_machine_descriptors = machine_descriptors
        elif (set(portable) != expected_keys or
              portable_descriptors != expected_portable_descriptors or
              machine_descriptors != expected_machine_descriptors):
            raise ValueError(
                "benchmark row sets differ between repetitions")
        keys = sorted({key[0] for key in expected_keys})
        ratios = []
        for schedule_id in keys:
            packed_portable = portable[(schedule_id, "packed_ff8")]
            xor_portable = portable[(schedule_id, "ff8xor_native")]
            packed_machine = machine[(schedule_id, "packed_ff8")]
            xor_machine = machine[(schedule_id, "ff8xor_native")]
            normalized = (
                (packed_machine["median_us"] / xor_machine["median_us"]) /
                (packed_portable["median_us"] / xor_portable["median_us"]))
            raw = xor_portable["median_us"] / xor_machine["median_us"]
            by_case[schedule_id].append(normalized)
            raw_by_case[schedule_id].append(raw)
            operation_ratios[xor_portable["operation"]].append(normalized)
            ratios.append(normalized)
        per_repetition.append(geometric_mean(ratios))

    case_medians = [statistics.median(values) for values in by_case.values()]
    result = {
        "mode": "quick",
        "normalization": (
            "(machine packed/machine ff8xor) / "
            "(portable packed/portable ff8xor); greater than one is faster"),
        "repetitions": len(portable_paths),
        "case_count": len(by_case),
        "observation_count": sum(len(values) for values in by_case.values()),
        "normalized_geomean_all_observations": geometric_mean([
            value for values in by_case.values() for value in values]),
        "normalized_geomean_case_medians": geometric_mean(case_medians),
        "normalized_case_median_wins": sum(value > 1 for value in case_medians),
        "per_repetition_geomean": per_repetition,
        "per_repetition_median": statistics.median(per_repetition),
        "raw_geomean_all_observations": geometric_mean([
            value for values in raw_by_case.values() for value in values]),
        "by_operation": {
            operation: {
                "observation_count": len(values),
                "normalized_geomean": geometric_mean(values),
            } for operation, values in sorted(operation_ratios.items())},
        "environment": {
            field: metadata.get(field) for field in (
                "compiler", "build_type", "build_flags", "cpu",
                "simd", "operating_system", "affinity")},
        "input_schema": metadata["schema"],
        "portable_input_sha256": [file_checksum(path)
                                  for path in portable_paths],
        "machine_input_sha256": [file_checksum(path)
                                 for path in machine_paths],
    }
    return result


def delta(portable, machine):
    if portable == 0:
        return None
    return (float(machine) / portable - 1.0) * 100.0


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--portable-circuits", type=Path, required=True)
    parser.add_argument("--machine-circuits", type=Path, required=True)
    parser.add_argument("--portable-artifact", type=Path, required=True)
    parser.add_argument("--machine-artifact", type=Path, required=True)
    parser.add_argument("--portable-benchmark", type=Path, action="append",
                        required=True)
    parser.add_argument("--machine-benchmark", type=Path, action="append",
                        required=True)
    parser.add_argument(
        "--profiles", type=Path,
        default=ROOT / "generated" / "FF8XorCostProfiles.json")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    profile_artifact = cost_model.load_profile_artifact(arguments.profiles)
    machine_profile = cost_model.find_profile(
        profile_artifact, "amd-zen5-gcc13-avx2-avx512-v1")
    family_weights = weights(profile_artifact)
    portable_source = circuit_source(arguments.portable_circuits, machine_profile)
    machine_source = circuit_source(arguments.machine_circuits, machine_profile)
    portable_functions, portable_assembly = compiled_functions(
        arguments.portable_artifact)
    machine_functions, machine_assembly = compiled_functions(
        arguments.machine_artifact)
    portable_weighted_source = weighted_source(portable_source, family_weights)
    machine_weighted_source = weighted_source(machine_source, family_weights)
    portable_avx2 = weighted_compiled(
        portable_functions, family_weights, "_avx2")
    machine_avx2 = weighted_compiled(
        machine_functions, family_weights, "_avx2")
    portable_avx512 = weighted_compiled(
        portable_functions, family_weights, "_avx512vl")
    machine_avx512 = weighted_compiled(
        machine_functions, family_weights, "_avx512vl")

    # Remove bulky per-coefficient programs from the checked-in summary.
    for source in (portable_source, machine_source):
        for family in FAMILIES:
            del source["families"][family]["circuits"]
            del source["families"][family]["modeled_costs"]

    comparison = {}
    for family in FAMILIES:
        comparison[family] = {
            "source_gate_weighted_delta_percent": delta(
                portable_weighted_source[family]["gate_weighted"],
                machine_weighted_source[family]["gate_weighted"]),
            "source_modeled_cost_weighted_delta_percent": delta(
                portable_weighted_source[family]["modeled_cost_weighted"],
                machine_weighted_source[family]["modeled_cost_weighted"]),
            "avx2_xor2_weighted_delta_percent": delta(
                portable_avx2[family]["xor2_weighted"],
                machine_avx2[family]["xor2_weighted"]),
            "avx2_code_bytes_weighted_delta_percent": delta(
                portable_avx2[family]["code_bytes_weighted"],
                machine_avx2[family]["code_bytes_weighted"]),
            "avx512vl_xor2_weighted_delta_percent": delta(
                portable_avx512[family]["xor2_weighted"],
                machine_avx512[family]["xor2_weighted"]),
            "avx512vl_xor3_weighted_delta_percent": delta(
                portable_avx512[family]["xor3_weighted"],
                machine_avx512[family]["xor3_weighted"]),
        }

    benchmark = benchmark_summary(
        arguments.portable_benchmark,
        arguments.machine_benchmark,
        {
            "circuit_checksum":
                portable_source["circuit_checksum_sha256"],
            "circuit_cost_profile_id": portable_source["cost_profile_id"],
            "circuit_cost_profile_checksum":
                portable_source["cost_profile_checksum_sha256"],
        },
        {
            "circuit_checksum":
                machine_source["circuit_checksum_sha256"],
            "circuit_cost_profile_id": machine_source["cost_profile_id"],
            "circuit_cost_profile_checksum":
                machine_source["cost_profile_checksum_sha256"],
        })
    result = {
        "schema": "leopard.ff8xor.cost-candidate-evaluation.v1",
        "machine_profile_checksum_sha256":
            machine_profile["checksum_sha256"],
        "workload_schedule_checksum_sha256":
            profile_artifact["workload"]["schedule_checksum_sha256"],
        "portable_source": portable_source,
        "machine_source": machine_source,
        "portable_assembly": portable_assembly,
        "machine_assembly": machine_assembly,
        "weighted_comparison": comparison,
        "benchmark": benchmark,
        "decision": {
            "enable_machine_profile_by_default": False,
            "measured_codec_speedup":
                benchmark["normalized_geomean_case_medians"],
            "reason": (
                "compiled instruction/code-size improvements did not produce "
                "a repeatable positive quick-codec result"),
        },
        "reproduction": {
            "generation": (
                "python3 tools/generate_ff8_xor_circuits.py --cost-profile "
                "amd-zen5-gcc13-avx2-avx512-v1 --output MACHINE.inl"),
            "build": "CMake Release, GCC 13.3, AVX2 and AVX-512 enabled",
            "benchmark": (
                "seven alternating process-order repetitions of "
                "bench_leopard_ff8xor --quick --json --abba --cpu 0"),
            "benchmark_input_schema": benchmark["input_schema"],
            "benchmark_input_note": (
                "The checked evidence captured benchmark schema v1 before "
                "traffic-accounting schema v2 landed; the evaluator accepts "
                "both known schemas and uses v2 fields for current binaries."
                if benchmark["input_schema"] ==
                HISTORICAL_BENCHMARK_SCHEMA else
                "Current traffic-accounting benchmark schema."),
        },
    }
    result["checksum_sha256"] = cost_model.checksum(result)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n",
        encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Evaluation checksum: %s" % result["checksum_sha256"])
    print("Normalized codec speedup: %.6f" %
          benchmark["normalized_geomean_case_medians"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
