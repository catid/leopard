#!/usr/bin/env python3
"""Fail-closed regression tests for FF8 cost-evidence evaluators."""

from __future__ import print_function

import copy
import json
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import benchmark_ff8xor_cost_model as calibration  # noqa: E402
import benchmark_ff8xor_cost_pairs as pairs  # noqa: E402
import evaluate_ff8xor_cost_candidate as evaluator  # noqa: E402


def require_raises(function, message):
    try:
        function()
    except ValueError:
        return
    raise AssertionError(message)


def require_runtime_raises(function, message):
    try:
        function()
    except RuntimeError:
        return
    raise AssertionError(message)


def metadata():
    return {
        "record": "metadata",
        "schema": evaluator.BENCHMARK_SCHEMA,
        "compiler": "test compiler",
        "build_type": "Release (NDEBUG)",
        "build_flags": "-O3 test flags",
        "cpu": "test cpu",
        "simd": "packed=AVX2; ff8xor=AVX2 XOR circuits",
        "circuit_checksum": "1" * 64,
        "circuit_cost_profile_id": "test-profile",
        "circuit_cost_profile_checksum": "2" * 64,
        "operating_system": "test os",
        "affinity": "pinned to logical CPU 0",
        "quick": True,
        "warmups": 1,
        "iterations": 3,
        "minimum_sample_usec": 250,
        "transpose_included": False,
        "counter_backend": "disabled",
        "pmu_requested": False,
        "cache_coloring_requested": False,
        "measurement_order": (
            "ABBA for paired end-to-end rows; other rows sequential"),
    }


def benchmark(schedule_id, backend, operation, loss_count, median_us):
    return {
        "record": "benchmark",
        "backend": backend,
        "operation": operation,
        "schedule_id": schedule_id,
        "k": 8,
        "r": 2,
        "buffer_bytes": 1024,
        "loss_count": loss_count,
        "transpose_included": False,
        "cache_coloring_applied": False,
        "measurement_order": "ABBA",
        "warmups": 1,
        "iterations": 6,
        "calls_per_sample": 17 if backend == "packed_ff8" else 11,
        "median_us": median_us,
        "modeled_payload_bytes_scheduled": (
            0 if backend == "packed_ff8" else 4096),
        "modeled_payload_bytes_elided": 0,
        "modeled_payload_bytes_adjusted": (
            0 if backend == "packed_ff8" else 4096),
        "materialization_deferred_zero_fills": 0,
        "materialization_added_zero_fills": 0,
        "materialization_butterflies_skipped": 0,
        "materialization_butterflies_reduced": 0,
        "materialization_xors_skipped": 0,
        "materialization_xors_replaced_by_copies": 0,
        "materialization_identity_operations_elided": 0,
        "locator_shift": (None if backend == "packed_ff8" or
                          operation == "encode" else 17),
    }


def records():
    result = [metadata()]
    for family in evaluator.FAMILIES:
        result.append({
            "record": "circuit_metadata",
            "family": family,
            "checksum": "1" * 64,
            "cost_profile_id": "test-profile",
            "cost_profile_checksum": "2" * 64,
        })
    for schedule_id, operation, loss_count in (
            ("encode-k8-r2-b1024", "encode", 0),
            ("decode-k8-r2-b1024-loss1", "decode", 1)):
        result.append(benchmark(
            schedule_id, "packed_ff8", operation, loss_count, 2.0))
        result.append(benchmark(
            schedule_id, "ff8xor_native", operation, loss_count, 4.0))
    return result


def write(path, rows):
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n"
                            for row in rows), encoding="utf-8")


def make_inputs(directory, portable_rows=None, machine_rows=None,
                repetitions=3):
    portable_rows = records() if portable_rows is None else portable_rows
    machine_rows = records() if machine_rows is None else machine_rows
    portable = []
    machine = []
    for repetition in range(repetitions):
        portable_path = directory / ("portable-%d.jsonl" % repetition)
        machine_path = directory / ("machine-%d.jsonl" % repetition)
        write(portable_path, portable_rows)
        write(machine_path, machine_rows)
        portable.append(portable_path)
        machine.append(machine_path)
    return portable, machine


def main():
    if "cpu >= CPU_SETSIZE" not in calibration.generated_source():
        raise AssertionError(
            "coefficient calibration can overrun cpu_set_t")
    empty_circuits = {
        family: {coefficient: () for coefficient in range(256)}
        for family in evaluator.FAMILIES}
    pair_source, unused_pairs = pairs.generated_source(
        empty_circuits, empty_circuits)
    if "cpu >= CPU_SETSIZE" not in pair_source:
        raise AssertionError(
            "equivalent-pair calibration can overrun cpu_set_t")

    with tempfile.TemporaryDirectory(prefix="leopard-cost-evaluator-test-") \
            as temporary:
        directory = Path(temporary)
        portable, machine = make_inputs(directory)
        portable_provenance = {
            "circuit_checksum": "1" * 64,
            "circuit_cost_profile_id": "test-profile",
            "circuit_cost_profile_checksum": "2" * 64,
        }
        machine_provenance = dict(portable_provenance)
        summary = evaluator.benchmark_summary(
            portable, machine, portable_provenance, machine_provenance)
        if summary["case_count"] != 2 or summary["observation_count"] != 6:
            raise AssertionError("complete paired rows were not summarized")

        historical = records()
        historical[0]["schema"] = evaluator.HISTORICAL_BENCHMARK_SCHEMA
        for field in ("circuit_checksum", "circuit_cost_profile_id",
                      "circuit_cost_profile_checksum"):
            del historical[0][field]
        for row in historical:
            if row.get("record") == "circuit_metadata":
                del row["cost_profile_id"]
                del row["cost_profile_checksum"]
            elif row.get("record") == "benchmark":
                row["modeled_payload_bytes"] = row.pop(
                    "modeled_payload_bytes_scheduled")
                for field in evaluator.BENCHMARK_V2_TRAFFIC_FIELDS[1:]:
                    del row[field]
        path = directory / "historical-v1.jsonl"
        write(path, historical)
        historical_metadata, unused_rows, unused_descriptors = \
            evaluator.benchmark_rows(path)
        if historical_metadata["circuit_checksum"] != "1" * 64:
            raise AssertionError(
                "historical family checksum was not normalized")

        different_locator = records()
        for row in different_locator:
            if (row.get("record") == "benchmark" and
                    row.get("backend") == "ff8xor_native" and
                    row.get("operation") == "decode"):
                row["locator_shift"] = 23
        p, m = make_inputs(directory, machine_rows=different_locator)
        evaluator.benchmark_summary(p, m)

        unrelated_circuit = records()
        unrelated_circuit[0]["circuit_checksum"] = "3" * 64
        p, m = make_inputs(directory, machine_rows=unrelated_circuit)
        require_raises(
            lambda: evaluator.benchmark_summary(
                p, m, portable_provenance, machine_provenance),
            "unrelated benchmark circuit corpus was accepted")

        unrelated_profile = records()
        unrelated_profile[0]["circuit_cost_profile_checksum"] = "4" * 64
        p, m = make_inputs(directory, machine_rows=unrelated_profile)
        require_raises(
            lambda: evaluator.benchmark_summary(
                p, m, portable_provenance, machine_provenance),
            "unrelated benchmark cost profile was accepted")

        missing_provenance = records()
        del missing_provenance[0]["circuit_cost_profile_id"]
        path = directory / "missing-provenance.jsonl"
        write(path, missing_provenance)
        require_raises(
            lambda: evaluator.benchmark_rows(path),
            "v2 benchmark without circuit provenance was accepted")

        malformed_provenance_cases = []
        for field, value in (
                ("circuit_checksum", None),
                ("circuit_checksum", "not-a-checksum"),
                ("circuit_cost_profile_id", None),
                ("circuit_cost_profile_id", ""),
                ("circuit_cost_profile_checksum", None),
                ("circuit_cost_profile_checksum", "not-a-checksum")):
            malformed = records()
            malformed[0][field] = value
            malformed_provenance_cases.append((field, malformed))
        for field, rows in malformed_provenance_cases:
            path = directory / ("malformed-metadata-%s.jsonl" % field)
            write(path, rows)
            require_raises(
                lambda path=path: evaluator.benchmark_rows(path),
                "malformed metadata %s was accepted" % field)

        for field, value in (
                ("compiler", None),
                ("build_type", ""),
                ("build_flags", None),
                ("cpu", ""),
                ("simd", None),
                ("operating_system", ""),
                ("affinity", None),
                ("measurement_order", "sequential"),
                ("counter_backend", ""),
                ("minimum_sample_usec", 0),
                ("minimum_sample_usec", True),
                ("pmu_requested", 0),
                ("cache_coloring_requested", None)):
            malformed = records()
            malformed[0][field] = value
            path = directory / ("malformed-environment-%s-%s.jsonl" %
                                (field, str(value)))
            write(path, malformed)
            require_raises(
                lambda path=path: evaluator.benchmark_rows(path),
                "malformed environment metadata %s was accepted" % field)

        for field, value in (
                ("checksum", None),
                ("checksum", "not-a-checksum"),
                ("cost_profile_id", None),
                ("cost_profile_id", ""),
                ("cost_profile_checksum", None),
                ("cost_profile_checksum", "not-a-checksum")):
            malformed = records()
            for row in malformed:
                if row.get("record") == "circuit_metadata":
                    row[field] = value
            path = directory / ("malformed-family-%s.jsonl" % field)
            write(path, malformed)
            require_raises(
                lambda path=path: evaluator.benchmark_rows(path),
                "malformed family %s was accepted" % field)

        missing_family_checksum = records()
        for row in missing_family_checksum:
            if row.get("record") == "circuit_metadata":
                del row["checksum"]
        path = directory / "missing-family-checksum.jsonl"
        write(path, missing_family_checksum)
        require_raises(
            lambda: evaluator.benchmark_rows(path),
            "missing family circuit checksum was accepted")

        semantic_cases = []

        invalid = records()
        for row in invalid:
            if row.get("record") == "benchmark":
                row["r"] = 9
                row["schedule_id"] = ("encode-k8-r9-b1024" if
                    row["operation"] == "encode" else
                    "decode-k8-r9-b1024-loss%d" % row["loss_count"])
        semantic_cases.append(("r-greater-than-k", invalid))

        invalid = records()
        for row in invalid:
            if row.get("record") == "benchmark":
                row["buffer_bytes"] = 65
                row["schedule_id"] = ("encode-k8-r2-b65" if
                    row["operation"] == "encode" else
                    "decode-k8-r2-b65-loss%d" % row["loss_count"])
        semantic_cases.append(("misaligned-buffer", invalid))

        invalid = records()
        for row in invalid:
            if (row.get("record") == "benchmark" and
                    row["operation"] == "encode"):
                row["loss_count"] = 1
        semantic_cases.append(("encode-loss", invalid))

        invalid = records()
        for row in invalid:
            if (row.get("record") == "benchmark" and
                    row["operation"] == "decode"):
                row["loss_count"] = 3
                row["schedule_id"] = "decode-k8-r2-b1024-loss3"
        semantic_cases.append(("excess-decode-loss", invalid))

        invalid = records()
        for row in invalid:
            if row.get("record") == "benchmark":
                row["schedule_id"] = "invented-schedule"
        semantic_cases.append(("noncanonical-schedule", invalid))

        invalid = records()
        for row in invalid:
            if row.get("record") == "benchmark":
                row["measurement_order"] = "sequential"
        semantic_cases.append(("unpaired-measurement", invalid))

        invalid = records()
        for row in invalid:
            if row.get("record") == "benchmark":
                row["iterations"] = 5
        semantic_cases.append(("iteration-count", invalid))

        invalid = records()
        for row in invalid:
            if (row.get("record") == "benchmark" and
                    row["backend"] == "ff8xor_native"):
                row["modeled_payload_bytes_adjusted"] += 1
        semantic_cases.append(("inconsistent-traffic", invalid))

        for label, rows in semantic_cases:
            path = directory / ("invalid-%s.jsonl" % label)
            write(path, rows)
            require_raises(
                lambda path=path: evaluator.benchmark_rows(path),
                "invalid %s benchmark was accepted" % label)

        missing = records()
        del missing[-1]
        p, m = make_inputs(directory, machine_rows=missing)
        require_raises(
            lambda: evaluator.benchmark_summary(p, m),
            "missing backend row was silently intersected away")

        renamed = records()
        for row in renamed:
            if row.get("schedule_id") == "encode-k8-r2-b1024":
                row["schedule_id"] = "renamed-schedule"
        p, m = make_inputs(directory, machine_rows=renamed)
        require_raises(
            lambda: evaluator.benchmark_summary(p, m),
            "renamed schedule was silently intersected away")

        mismatched_environment = records()
        mismatched_environment[0]["cpu"] = "different cpu"
        p, m = make_inputs(
            directory, machine_rows=mismatched_environment)
        require_raises(
            lambda: evaluator.benchmark_summary(p, m),
            "mismatched benchmark environment was accepted")

        duplicate = records()
        duplicate.append(copy.deepcopy(duplicate[-1]))
        path = directory / "duplicate.jsonl"
        write(path, duplicate)
        require_raises(
            lambda: evaluator.benchmark_rows(path),
            "duplicate codec row was accepted")

        bad_schema = records()
        bad_schema[0]["schema"] = "stale.schema"
        path = directory / "bad-schema.jsonl"
        write(path, bad_schema)
        require_raises(
            lambda: evaluator.benchmark_rows(path),
            "stale benchmark schema was accepted")

        p, m = make_inputs(directory, repetitions=2)
        require_raises(
            lambda: evaluator.benchmark_summary(p, m),
            "fewer than three benchmark repetitions were accepted")

    for arguments in (
            (1, 3, 0), (2, 2, 0), (2, 3, -1),
            (1 << 64, 3, 0), (2, 1 << 32, 0), (2, 3, 1 << 32),
            (2, 3, True), (2.0, 3, 0)):
        require_raises(
            lambda arguments=arguments: pairs.validate_options(*arguments),
            "invalid equivalent-pair benchmark options were accepted")
    pairs.validate_options(2, 3, 0)
    for arguments in (
            (0, 3, 0), (1, 2, 0), (1, 3, -1),
            (1 << 64, 3, 0), (1, 1 << 32, 0), (1, 3, 1 << 32),
            (1, 3, True), (1.0, 3, 0)):
        require_raises(
            lambda arguments=arguments:
                calibration.validate_options(*arguments),
            "invalid coefficient calibration options were accepted")
    calibration.validate_options(1, 3, 0)
    require_raises(
        lambda: pairs.validate_pairs([], {
            family: {} for family in evaluator.FAMILIES}),
        "empty equivalent-pair corpus was accepted")
    require_raises(
        lambda: pairs.validate_pairs([("multiply", 7)], {
            family: {} for family in evaluator.FAMILIES}),
        "zero-weight equivalent-pair corpus was accepted")

    pair_header = (
        "family,coefficient,portable_median_ns,machine_median_ns,"
        "portable_best_ns,machine_best_ns,portable_mad_ns,machine_mad_ns")
    require_runtime_raises(
        lambda: pairs.parse_timings(
            pair_header + "\nmultiply,7,nan,1,1,1,0,0\n",
            [("multiply", 7)]),
        "non-finite equivalent-pair timing was accepted")
    calibration_lines = [
        "checksum," + "0" * 64,
        "coefficient,median_ns,best_ns,mad_ns",
    ] + ["%d,%s,1,0" % (coefficient,
                         "nan" if coefficient == 7 else "1")
         for coefficient in range(256)]
    require_runtime_raises(
        lambda: calibration.parse_timings("\n".join(calibration_lines)),
        "non-finite coefficient calibration timing was accepted")
    print("FF8 XOR cost-evaluator robustness tests passed")


if __name__ == "__main__":
    main()
