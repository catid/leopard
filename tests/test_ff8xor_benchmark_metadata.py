#!/usr/bin/env python3
"""Check that finite benchmark rows describe their actual sampling order."""

from __future__ import print_function

import csv
import io
import json
import subprocess
import sys
from pathlib import Path


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def next_power_of_two(value):
    result = 1
    while result < value:
        result <<= 1
    return result


MASK64 = (1 << 64) - 1


def shuffled_indices(count, seed):
    indices = list(range(count))
    state = seed if seed else 1
    for i in range(count, 1, -1):
        value = state
        value ^= value >> 12
        value ^= (value << 25) & MASK64
        value ^= value >> 27
        state = value & MASK64
        random_value = (state * 2685821657736338717) & MASK64
        other = random_value % i
        indices[i - 1], indices[other] = indices[other], indices[i - 1]
    return indices


def is_power_of_two(value):
    return value != 0 and value & (value - 1) == 0


def transform_memory_units(count_truncated, count, skew_base):
    units = 0
    distance = 1
    while distance < count:
        span = distance * 2
        for range_start in range(0, count_truncated, span):
            sentinel = is_power_of_two(
                skew_base + range_start + distance)
            units += distance * (3 if sentinel else 4)
        distance *= 2
    return units


def error_fft_memory_units(count_truncated, initial_distance,
                           skew_base, missing_locations):
    units = 0
    distance = initial_distance
    while distance:
        span = distance * 2
        for range_start in range(0, count_truncated, span):
            needed = any(range_start <= location < range_start + span
                         for location in missing_locations)
            if needed:
                sentinel = is_power_of_two(
                    skew_base + range_start + distance)
                units += distance * (3 if sentinel else 4)
        distance //= 2
    return units


def fused_derivative_boundary_units(count):
    half_count = count // 2
    units = 0
    for q in range(half_count):
        units += 4  # Read A/R and write X/Y.
        bit = 1
        while bit < half_count:
            if not q & bit:
                units += 1
            bit *= 2
    return units


def modeled_fused_decode_bytes(original_count, recovery_count,
                               buffer_bytes, loss_count):
    m = next_power_of_two(recovery_count)
    n = next_power_of_two(m + original_count)
    seed = (0xff8c000000000000 ^ (original_count << 40) ^
            (recovery_count << 24) ^ buffer_bytes)
    loss_seed = seed ^ (loss_count << 8) ^ 0xd3c0de
    order = shuffled_indices(original_count, loss_seed)
    missing_locations = [m + order[i] for i in range(loss_count)]

    half_count = n // 2
    left_derivative_buffers = sum(
        ((index ^ (index - 1)) + 1) >> 1
        for index in range(1, half_count))
    units = 0
    units += (original_count + loss_count) * 2
    units += n - original_count
    units += transform_memory_units(m + original_count, n, 0)
    units += left_derivative_buffers * 3
    units += fused_derivative_boundary_units(n)
    units += error_fft_memory_units(
        m + original_count, n // 4, 0, missing_locations)
    return units * buffer_bytes


def main():
    require(len(sys.argv) == 2, "expected benchmark executable path")
    repository_root = Path(__file__).resolve().parents[1]
    profile_artifact = json.loads((repository_root / "generated" /
                                   "FF8XorCostProfiles.json").read_text(
                                       encoding="utf-8"))
    portable_profile = next(
        profile for profile in profile_artifact["profiles"]
        if profile["id"] == "portable-default-v1")
    schedule_artifact = json.loads((repository_root / "generated" /
                                    "FF8XorScheduleCorpus.json").read_text(
                                        encoding="utf-8"))
    expected_locator_shifts = {
        record["id"]: record["locator_shift"]
        for record in schedule_artifact["records"]
        if record["operation"] == "decode"
    }
    result = subprocess.run(
        [sys.argv[1], "--quick", "--json", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(result.returncode == 0,
            "quick benchmark failed (%d): %s" %
            (result.returncode, result.stderr))

    records = [json.loads(line) for line in result.stdout.splitlines()
               if line.strip()]
    metadata = [record for record in records
                if record.get("record") == "metadata"]
    require(len(metadata) == 1,
            "benchmark did not emit exactly one metadata record")
    require(metadata[0].get("schema") ==
            "leopard.ff8xor.benchmark.jsonl.v2",
            "benchmark traffic-accounting schema version changed")
    require("XOR-batch pairs ABBA" in metadata[0]["measurement_order"],
            "metadata omitted unconditional XOR-batch ABBA sampling")
    require("locator-selector pairs" in metadata[0]["measurement_order"],
            "metadata omitted unconditional locator-selector ABBA sampling")
    require(metadata[0].get("circuit_cost_profile_id") ==
            "portable-default-v1",
            "metadata omitted the generated circuit selection profile")
    require(metadata[0].get("circuit_cost_profile_checksum") ==
            portable_profile["checksum_sha256"],
            "metadata circuit profile checksum is stale")
    require(len(metadata[0].get("circuit_checksum", "")) == 64,
            "metadata omitted the circuit checksum")
    require(metadata[0].get("ff8xor_mode_requested") == "auto",
            "metadata omitted the requested default kernel mode")
    require(metadata[0].get("four_buffer_mode_requested") == "disabled",
            "metadata omitted the requested default four-buffer mode")
    circuit_metadata = [record for record in records
                        if record.get("record") == "circuit_metadata"]
    require(len(circuit_metadata) == 3,
            "benchmark did not emit all circuit metadata families")
    require(all(record.get("cost_profile_id") ==
                metadata[0]["circuit_cost_profile_id"] and
                record.get("cost_profile_checksum") ==
                metadata[0]["circuit_cost_profile_checksum"]
                for record in circuit_metadata),
            "circuit rows disagree with top-level profile provenance")

    expected_operations = {
        "xor2_sequential", "xor2_batched",
        "xor4_sequential", "xor4_batched",
    }
    rows = [record for record in records
            if record.get("record") == "microbenchmark" and
            record.get("operation") in expected_operations]
    observed_operations = {record["operation"] for record in rows}
    require(observed_operations == expected_operations,
            "quick benchmark XOR-batch row set changed: %r" %
            sorted(observed_operations))
    require(len(rows) == 4 * len(expected_operations),
            "quick benchmark did not cover four XOR-batch sizes")
    for record in rows:
        require(record.get("measurement_order") == "ABBA",
                "%s mislabeled paired sampling as %r" %
                (record["operation"], record.get("measurement_order")))
        require("paired ABBA" in record.get("note", ""),
                "%s row lost paired-sampling note" % record["operation"])

    selector_rows = [record for record in records
                     if record.get("record") == "microbenchmark" and
                     record.get("operation") == "locator_shift_select"]
    expected_selector_schedules = {
        "locator_shift_select-k8-r2-b0-loss1",
        "locator_shift_select-k16-r4-b0-loss4",
        "locator_shift_select-k32-r8-b0-loss4",
        "locator_shift_select-k64-r16-b0-loss8",
        "locator_shift_select-k128-r32-b0-loss16",
        "locator_shift_select-k128-r128-b0-loss128",
    }
    require(len(selector_rows) == 12,
            "quick benchmark did not emit six locator-selector ABBA pairs")
    selector_pairs = {}
    for record in selector_rows:
        require(record.get("measurement_order") == "ABBA" and
                "paired ABBA" in record.get("note", ""),
                "locator-selector row mislabeled paired sampling")
        require(record.get("locator_shift") is not None and
                0 <= int(record["locator_shift"]) < 255,
                "locator-selector row omitted its selected shift")
        require(record.get("schedule_id", "").endswith(
                    "-loss%d" % int(record["loss_count"])),
                "locator-selector schedule ID omitted its loss count")
        selector_pairs.setdefault(record["schedule_id"], []).append(record)
    require(set(selector_pairs) == expected_selector_schedules,
            "locator-selector schedule set changed: %r" %
            sorted(selector_pairs))
    for schedule_id, pair in selector_pairs.items():
        require({record["backend"] for record in pair} == {
                    "ff8xor_selector_reference",
                    "ff8xor_selector_rotated"} and
                len(pair) == 2 and
                len({int(record["locator_shift"]) for record in pair}) == 1,
                "locator-selector pair disagrees for %s" % schedule_id)

    traffic_fields = {
        "modeled_payload_bytes_scheduled",
        "modeled_payload_bytes_elided",
        "modeled_payload_bytes_adjusted",
        "materialization_deferred_zero_fills",
        "materialization_added_zero_fills",
        "materialization_butterflies_skipped",
        "materialization_butterflies_reduced",
        "materialization_xors_skipped",
        "materialization_xors_replaced_by_copies",
        "materialization_identity_operations_elided",
        "four_buffer_fused_units",
        "four_buffer_payload_bytes_elided",
    }
    end_to_end = [record for record in records
                  if record.get("record") == "benchmark" and
                  record.get("backend") in
                  ("packed_ff8", "ff8xor_native")]
    require(end_to_end, "quick benchmark emitted no end-to-end rows")
    saw_tracked_decode = False
    operation_fields = traffic_fields - {
        "modeled_payload_bytes_scheduled",
        "modeled_payload_bytes_elided",
        "modeled_payload_bytes_adjusted",
    }
    for record in end_to_end:
        missing = traffic_fields - set(record)
        require(not missing,
                "%s/%s row omitted traffic fields: %r" %
                (record["backend"], record["operation"], sorted(missing)))
        scheduled = record["modeled_payload_bytes_scheduled"]
        elided = record["modeled_payload_bytes_elided"]
        adjusted = record["modeled_payload_bytes_adjusted"]
        require(adjusted == scheduled - elided,
                "%s/%s traffic accounting is inconsistent" %
                (record["backend"], record["operation"]))
        if record["backend"] == "packed_ff8":
            require(scheduled == 0 and elided == 0 and adjusted == 0,
                    "packed row unexpectedly reported ff8xor traffic")
            require(all(record[field] == 0 for field in operation_fields),
                    "packed row unexpectedly reported materialization work")
            continue

        require(scheduled > 0,
                "native end-to-end row omitted scheduled traffic")
        if record["operation"] == "decode":
            expected = modeled_fused_decode_bytes(
                record["k"], record["r"], record["buffer_bytes"],
                record["loss_count"])
            require(scheduled == expected,
                    "native decode fused-traffic model is stale: "
                    "expected %d, observed %d for %s" %
                    (expected, scheduled, record["schedule_id"]))
            require(record["schedule_id"] in expected_locator_shifts and
                    record.get("locator_shift") ==
                    expected_locator_shifts[record["schedule_id"]],
                    "native decode production locator gather selected a "
                    "stale shift for %s: expected %r, observed %r" %
                    (record["schedule_id"],
                     expected_locator_shifts.get(record["schedule_id"]),
                     record.get("locator_shift")))
        full_chunk_encode = (record["operation"] == "encode" and
                             record["k"] %
                             next_power_of_two(record["r"]) == 0)
        if record["buffer_bytes"] < 65536 or full_chunk_encode:
            require(elided == 0,
                    "bypassed native row reported stale elided traffic")
            require(all(record[field] == 0 for field in operation_fields),
                    "bypassed native row reported stale materialization ops")
        elif record["operation"] == "decode" and elided > 0:
            saw_tracked_decode = True
    require(saw_tracked_decode,
            "quick benchmark did not expose any tracked decode elision")

    csv_result = subprocess.run(
        [sys.argv[1], "--quick", "--csv", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(csv_result.returncode == 0,
            "quick CSV benchmark failed (%d): %s" %
            (csv_result.returncode, csv_result.stderr))
    reader = csv.DictReader(io.StringIO(csv_result.stdout))
    require(reader.fieldnames is not None,
            "CSV benchmark emitted no header")
    require(not (traffic_fields - set(reader.fieldnames)),
            "CSV header omitted traffic fields: %r" %
            sorted(traffic_fields - set(reader.fieldnames)))
    profile_fields = {"cost_profile_id", "cost_profile_checksum"}
    require(not (profile_fields - set(reader.fieldnames)),
            "CSV header omitted circuit profile provenance")
    require("ff8xor_mode_requested" in reader.fieldnames,
            "CSV header omitted the requested kernel mode")
    require("four_buffer_mode_requested" in reader.fieldnames,
            "CSV header omitted the requested four-buffer mode")
    csv_records = list(reader)
    require(csv_records, "CSV benchmark emitted no records")
    require(all(record["ff8xor_mode_requested"] == "auto"
                for record in csv_records),
            "CSV rows disagree about the requested kernel mode")
    require(all(record["four_buffer_mode_requested"] == "disabled"
                for record in csv_records),
            "CSV rows disagree about the requested four-buffer mode")
    require(all(None not in record for record in csv_records),
            "CSV benchmark emitted a row wider than its header")
    require(all(len(record) == len(reader.fieldnames) and
                all(value is not None for value in record.values())
                for record in csv_records),
            "CSV benchmark schema is not rectangular")
    csv_circuit_metadata = [record for record in csv_records
                            if record.get("record") == "metadata" and
                            record.get("backend") == "generated"]
    require(len(csv_circuit_metadata) == 3 and
            all(record["cost_profile_id"] == "portable-default-v1" and
                record["cost_profile_checksum"] ==
                portable_profile["checksum_sha256"]
                for record in csv_circuit_metadata),
            "CSV circuit metadata omitted profile provenance")

    csv_selector_rows = [record for record in csv_records
                         if record.get("record") == "microbenchmark" and
                         record.get("operation") == "locator_shift_select"]
    csv_selector_pairs = {}
    for record in csv_selector_rows:
        require(record.get("measurement_order") == "ABBA" and
                "paired ABBA" in record.get("note", "") and
                record.get("locator_shift") not in ("", "null") and
                0 <= int(record["locator_shift"]) < 255 and
                record.get("schedule_id", "").endswith(
                    "-loss%s" % record.get("loss_count")),
                "CSV locator-selector paired metadata is incomplete")
        csv_selector_pairs.setdefault(
            record["schedule_id"], []).append(record)
    require(set(csv_selector_pairs) == expected_selector_schedules,
            "CSV locator-selector schedule set changed: %r" %
            sorted(csv_selector_pairs))
    for schedule_id, pair in csv_selector_pairs.items():
        require({record["backend"] for record in pair} == {
                    "ff8xor_selector_reference",
                    "ff8xor_selector_rotated"} and
                len(pair) == 2 and
                len({int(record["locator_shift"]) for record in pair}) == 1,
                "CSV locator-selector pair disagrees for %s" % schedule_id)

    csv_end_to_end = [record for record in csv_records
                      if record.get("record") == "benchmark" and
                      record.get("backend") in
                      ("packed_ff8", "ff8xor_native")]
    require(csv_end_to_end,
            "quick CSV benchmark emitted no end-to-end rows")
    csv_saw_tracked_decode = False
    for record in csv_end_to_end:
        scheduled = int(record["modeled_payload_bytes_scheduled"])
        elided = int(record["modeled_payload_bytes_elided"])
        adjusted = int(record["modeled_payload_bytes_adjusted"])
        require(adjusted == scheduled - elided,
                "CSV %s/%s traffic accounting is inconsistent" %
                (record["backend"], record["operation"]))
        operation_values = [int(record[field]) for field in operation_fields]
        if record["backend"] == "packed_ff8":
            require(scheduled == 0 and elided == 0 and adjusted == 0 and
                    all(value == 0 for value in operation_values),
                    "CSV packed row unexpectedly reported materialization")
            continue
        buffer_bytes = int(record["buffer_bytes"])
        require(scheduled > 0,
                "CSV native row omitted scheduled traffic")
        if record["operation"] == "decode":
            expected = modeled_fused_decode_bytes(
                int(record["k"]), int(record["r"]), buffer_bytes,
                int(record["loss_count"]))
            require(scheduled == expected,
                    "CSV native decode fused-traffic model is stale: "
                    "expected %d, observed %d for %s" %
                    (expected, scheduled, record["schedule_id"]))
            require(record["schedule_id"] in expected_locator_shifts and
                    int(record["locator_shift"]) ==
                    expected_locator_shifts[record["schedule_id"]],
                    "CSV native decode production locator gather selected a "
                    "stale shift for %s" % record["schedule_id"])
        full_chunk_encode = (record["operation"] == "encode" and
                             int(record["k"]) %
                             next_power_of_two(int(record["r"])) == 0)
        if buffer_bytes < 65536 or full_chunk_encode:
            require(elided == 0 and
                    all(value == 0 for value in operation_values),
                    "CSV bypass row reported stale materialization")
        elif record["operation"] == "decode" and elided > 0:
            csv_saw_tracked_decode = True
    require(csv_saw_tracked_decode,
            "quick CSV benchmark did not expose tracked decode elision")

    invalid_mode = subprocess.run(
        [sys.argv[1], "--ff8xor-mode", "not-a-mode", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(invalid_mode.returncode == 2 and
            "Invalid --ff8xor-mode value" in invalid_mode.stderr,
            "benchmark accepted or misreported an invalid kernel mode")

    invalid_four_buffer_mode = subprocess.run(
        [sys.argv[1], "--four-buffer-mode", "not-a-mode", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(invalid_four_buffer_mode.returncode == 2 and
            "Invalid --four-buffer-mode value" in
            invalid_four_buffer_mode.stderr,
            "benchmark accepted or misreported an invalid four-buffer mode")

    mislabeled_four_buffer_mode = subprocess.run(
        [sys.argv[1], "--four-buffer-mode", "xor2", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(mislabeled_four_buffer_mode.returncode == 2 and
            "requires --ff8xor-mode avx512zmm" in
            mislabeled_four_buffer_mode.stderr,
            "benchmark allowed a four-buffer label on a non-ZMM run")

    print("FF8 XOR benchmark metadata tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
