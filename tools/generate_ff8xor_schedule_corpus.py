#!/usr/bin/env python3
"""Generate deterministic FF8 XOR transform/frequency workloads.

The corpus is deliberately independent of the C++ payload implementation.  It
reproduces the scalar field metadata, transform schedule, locator construction,
locator-shift selection, and decoder ErrorBitfield pruning used by
LeopardFF8Xor.cpp.  It is intended for profile-guided circuit generation and
for comparing optimization experiments against one stable workload mix.
"""

from __future__ import print_function

import argparse
import collections
import difflib
import hashlib
import json
import sys
from pathlib import Path

import generate_ff8_xor_circuits as circuits


SCHEMA = "leopard.ff8xor.schedule-corpus.v1"
FULL_COUNTS = ((8, 2), (16, 4), (32, 8), (64, 16), (128, 32), (128, 128))
FULL_BUFFER_BYTES = (1024, 4096, 64 * 1024, 1024 * 1024)
MASK64 = (1 << 64) - 1


def next_power_of_two(value):
    result = 1
    while result < value:
        result <<= 1
    return result


def sub_mod(a, b):
    difference = (a - b) & 0xFFFFFFFF
    return (difference + (difference >> circuits.FIELD_BITS)) & 0xFF


def fwht_2(data, left, right):
    a = data[left]
    b = data[right]
    data[left] = circuits.add_mod(a, b)
    data[right] = sub_mod(a, b)


def fwht_4(data, offset, stride):
    stride2 = stride << 1
    indices = (offset, offset + stride, offset + stride2,
               offset + stride2 + stride)
    values = [data[index] for index in indices]
    fwht_2(values, 0, 1)
    fwht_2(values, 2, 3)
    fwht_2(values, 0, 2)
    fwht_2(values, 1, 3)
    for index, value in zip(indices, values):
        data[index] = value


def fwht(data, count, count_truncated):
    distance = 1
    distance4 = 4
    while distance4 <= count:
        for range_start in range(0, count_truncated, distance4):
            for index in range(range_start, range_start + distance):
                fwht_4(data, index, distance)
        distance = distance4
        distance4 <<= 2
    if distance < count:
        for index in range(distance):
            fwht_2(data, index, index + distance)


def make_fft_metadata():
    temporary = [1 << bit for bit in range(1, circuits.FIELD_BITS)]
    # Logical FFTSkew[i] is padded[i + 1].  padded[0] is the -1 sentinel.
    padded = [0] * (circuits.FIELD_ORDER + 1)
    padded[0] = circuits.FIELD_MODULUS

    for level in range(circuits.FIELD_BITS - 1):
        step = 1 << (level + 1)
        padded[(1 << level)] = 0
        for bit in range(level, circuits.FIELD_BITS - 1):
            span = 1 << (bit + 1)
            for index in range((1 << level) - 1, span, step):
                padded[index + span + 1] = padded[index + 1] ^ temporary[bit]

        value = circuits.scalar_multiply_log(
            temporary[level], circuits.LOG_LUT[temporary[level] ^ 1])
        temporary[level] = circuits.FIELD_MODULUS - circuits.LOG_LUT[value]

        for bit in range(level + 1, circuits.FIELD_BITS - 1):
            total = circuits.add_mod(
                circuits.LOG_LUT[temporary[bit] ^ 1], temporary[level])
            temporary[bit] = circuits.scalar_multiply_log(temporary[bit], total)

    for index in range(circuits.FIELD_MODULUS):
        padded[index + 1] = circuits.LOG_LUT[padded[index + 1]]

    log_walsh = list(circuits.LOG_LUT)
    log_walsh[0] = 0
    fwht(log_walsh, circuits.FIELD_ORDER, circuits.FIELD_ORDER)
    return tuple(padded), tuple(log_walsh)


FFT_SKEW_PADDED, LOG_WALSH = make_fft_metadata()


class Random(object):
    def __init__(self, seed):
        self.state = seed if seed else 1

    def next(self):
        value = self.state
        value ^= value >> 12
        value ^= (value << 25) & MASK64
        value ^= value >> 27
        self.state = value & MASK64
        return (self.state * 2685821657736338717) & MASK64


def shuffled_indices(count, seed):
    indices = list(range(count))
    random = Random(seed)
    for index in range(count, 1, -1):
        other = random.next() % index
        indices[index - 1], indices[other] = indices[other], indices[index - 1]
    return indices


def loss_counts(recovery_count):
    result = []
    for value in (1, min(4, recovery_count), max(1, recovery_count // 2),
                  recovery_count):
        if value not in result:
            result.append(value)
    return result


def histogram(values):
    counts = collections.Counter(values)
    return {str(key): counts[key] for key in sorted(counts)}


def tuple_histogram(values):
    counts = collections.Counter(values)
    return {"%d,%d,%d" % key: counts[key] for key in sorted(counts)}


def transform_invocations(count_truncated, count, skew_base, inverse,
                          missing_locations=None):
    """Return literal two-way invocations and complete two-layer tuples."""
    stages = []
    if inverse:
        distances = []
        distance = 1
        while distance < count:
            distances.append(distance)
            distance <<= 1
    else:
        distances = []
        distance = count >> 1
        while distance:
            distances.append(distance)
            distance >>= 1

    def needed(distance, range_start):
        if missing_locations is None:
            return True
        # This is equivalent to ErrorBitfield::IsNeeded for n <= 256: each
        # level expands every requested output to its complete butterfly block.
        return any(range_start <= location < range_start + (distance << 1)
                   for location in missing_locations)

    for distance in distances:
        stage = []
        span = distance << 1
        for range_start in range(0, count_truncated, span):
            if not needed(distance, range_start):
                continue
            skew = FFT_SKEW_PADDED[skew_base + range_start + distance]
            stage.extend([skew] * distance)
        stages.append((distance, stage))

    invocations = [skew for unused_distance, stage in stages for skew in stage]

    # A tuple describes the three coefficients of one complete four-buffer,
    # two-layer unit: pair 0/1, pair 2/3, then pair 0/2.  Count it once per
    # payload butterfly lane (distance times), matching whole-buffer calls.
    tuples = []
    low_distance = 1
    while (low_distance << 1) < count:
        block = low_distance << 2
        for range_start in range(0, count_truncated, block):
            second_range = range_start + (low_distance << 1)
            if second_range >= count_truncated:
                continue
            if (not needed(low_distance, range_start) or
                    not needed(low_distance, second_range) or
                    not needed(low_distance << 1, range_start)):
                continue
            tuple_value = (
                FFT_SKEW_PADDED[skew_base + range_start + low_distance],
                FFT_SKEW_PADDED[skew_base + second_range + low_distance],
                FFT_SKEW_PADDED[skew_base + range_start + (low_distance << 1)])
            tuples.extend([tuple_value] * low_distance)
        low_distance <<= 2

    return invocations, tuples


def locator_logs(original_count, recovery_count, m, original_missing,
                 recovery_available):
    locations = [0] * circuits.FIELD_ORDER
    for index in range(recovery_count):
        if index not in recovery_available:
            locations[index] = 1
    for index in range(recovery_count, m):
        locations[index] = 1
    for index in original_missing:
        locations[m + index] = 1

    fwht(locations, circuits.FIELD_ORDER, m + original_count)
    for index in range(circuits.FIELD_ORDER):
        locations[index] = (locations[index] * LOG_WALSH[index]) % 255
    fwht(locations, circuits.FIELD_ORDER, circuits.FIELD_ORDER)
    return locations


def canonical_log(value):
    return 0 if value == circuits.FIELD_MODULUS else value


def shifted_log(value, shift):
    result = canonical_log(value) + shift
    if result >= circuits.FIELD_MODULUS:
        result -= circuits.FIELD_MODULUS
    return result


def inverse_shifted_log(value, shift):
    shifted = shifted_log(value, shift)
    return 0 if shifted == 0 else circuits.FIELD_MODULUS - shifted


def select_locator_shift(locations, original_count, recovery_count, m,
                         original_missing, recovery_available, gate_counts,
                         depths):
    best = None
    for shift in range(circuits.FIELD_MODULUS):
        logs = []
        for index in sorted(recovery_available):
            logs.append(shifted_log(locations[index], shift))
        for index in range(original_count):
            if index in original_missing:
                logs.append(inverse_shifted_log(locations[m + index], shift))
            else:
                logs.append(shifted_log(locations[m + index], shift))
        score = (sum(gate_counts[value] for value in logs),
                 sum(depths[value] for value in logs), shift)
        if best is None or score < best:
            best = score
    return best[2]


def payload_model(counts, buffer_bytes):
    factors = {
        "copy": 2,
        "zero": 1,
        "xor": 3,
        "multiply": 2,
        "butterfly": 4,
        "butterfly_sentinel": 3,
    }
    units = sum(counts.get(name, 0) * factor
                for name, factor in factors.items())
    return {
        "whole_buffer_operations": dict(sorted(counts.items())),
        "bytes_per_shard_byte": units,
        "bytes": units * buffer_bytes,
        "model": ("loads+stores; copy=2 zero=1 xor=3 multiply=2 "
                  "butterfly=4 sentinel-butterfly=3"),
    }


def encode_record(original_count, recovery_count, buffer_bytes):
    m = next_power_of_two(recovery_count)
    n = next_power_of_two(m + original_count)
    ifft_skews = []
    fft_skews = []
    ifft_tuples = []
    chunks = 0
    zero_buffers = 0
    chunk_start = 0
    while chunk_start < original_count:
        chunk_count = min(m, original_count - chunk_start)
        skews, chunk_tuples = transform_invocations(
            chunk_count, m, m + chunk_start, True)
        ifft_skews.extend(skews)
        ifft_tuples.extend(chunk_tuples)
        zero_buffers += m - chunk_count
        chunks += 1
        chunk_start += m
    skews, final_tuples = transform_invocations(recovery_count, m, 0, False)
    fft_skews.extend(skews)
    fft_tuples = final_tuples

    operation_counts = {
        "copy": original_count,
        "zero": zero_buffers,
        "xor": m * (chunks - 1),
        "butterfly": sum(1 for value in ifft_skews + fft_skews
                           if value != circuits.FIELD_MODULUS),
        "butterfly_sentinel": sum(1 for value in ifft_skews + fft_skews
                                    if value == circuits.FIELD_MODULUS),
    }
    return {
        "id": "encode-k%d-r%d-b%d" % (
            original_count, recovery_count, buffer_bytes),
        "operation": "encode",
        "k": original_count,
        "r": recovery_count,
        "buffer_bytes": buffer_bytes,
        "m": m,
        "n": n,
        "loss_count": 0,
        "seed": "0x%016x" % ((0xFF8C000000000000 ^
                                (original_count << 40) ^
                                (recovery_count << 24) ^ buffer_bytes) & MASK64),
        "multiply_log_frequency": {},
        "ifft_skew_frequency": histogram(ifft_skews),
        "fft_skew_frequency": histogram(fft_skews),
        "ifft_four_buffer_tuple_frequency": tuple_histogram(ifft_tuples),
        "fft_four_buffer_tuple_frequency": tuple_histogram(fft_tuples),
        "four_buffer_tuple_frequency": tuple_histogram(
            ifft_tuples + fft_tuples),
        "payload_memory": payload_model(operation_counts, buffer_bytes),
    }


def decode_record(original_count, recovery_count, buffer_bytes, loss_count,
                  gate_counts, depths):
    m = next_power_of_two(recovery_count)
    n = next_power_of_two(m + original_count)
    seed = (0xFF8C000000000000 ^ (original_count << 40) ^
            (recovery_count << 24) ^ buffer_bytes) & MASK64
    loss_seed = seed ^ (loss_count << 8) ^ 0xD3C0DE
    original_order = shuffled_indices(original_count, loss_seed)
    recovery_order = shuffled_indices(
        recovery_count, loss_seed ^ 0x7265636F76657279)
    original_missing = frozenset(original_order[:loss_count])
    recovery_available = frozenset(recovery_order[:loss_count])

    locations = locator_logs(original_count, recovery_count, m,
                             original_missing, recovery_available)
    shift = select_locator_shift(
        locations, original_count, recovery_count, m, original_missing,
        recovery_available, gate_counts, depths)

    input_multiply_logs = []
    for index in sorted(recovery_available):
        input_multiply_logs.append(shifted_log(locations[index], shift))
    for index in range(original_count):
        if index not in original_missing:
            input_multiply_logs.append(
                shifted_log(locations[m + index], shift))
    recovery_multiply_logs = []
    for index in sorted(original_missing):
        recovery_multiply_logs.append(
            inverse_shifted_log(locations[m + index], shift))
    multiply_logs = input_multiply_logs + recovery_multiply_logs

    ifft_skews, ifft_tuples = transform_invocations(
        m + original_count, n, 0, True)
    missing_locations = tuple(m + index for index in sorted(original_missing))
    fft_skews, fft_tuples = transform_invocations(
        m + original_count, n, 0, False, missing_locations)
    derivative_buffers = 0
    for index in range(1, n):
        derivative_buffers += ((index ^ (index - 1)) + 1) >> 1

    present_inputs = len(recovery_available) + original_count - loss_count
    operation_counts = {
        "multiply": present_inputs + loss_count,
        "zero": n - present_inputs,
        "xor": derivative_buffers,
        "butterfly": sum(1 for value in ifft_skews + fft_skews
                           if value != circuits.FIELD_MODULUS),
        "butterfly_sentinel": sum(1 for value in ifft_skews + fft_skews
                                    if value == circuits.FIELD_MODULUS),
    }
    return {
        "id": "decode-k%d-r%d-b%d-loss%d" % (
            original_count, recovery_count, buffer_bytes, loss_count),
        "operation": "decode",
        "k": original_count,
        "r": recovery_count,
        "buffer_bytes": buffer_bytes,
        "m": m,
        "n": n,
        "loss_count": loss_count,
        "seed": "0x%016x" % seed,
        "loss_seed": "0x%016x" % (loss_seed & MASK64),
        "missing_original_indices": sorted(original_missing),
        "available_recovery_indices": sorted(recovery_available),
        "locator_shift": shift,
        "input_multiply_log_frequency": histogram(input_multiply_logs),
        "recovery_multiply_log_frequency": histogram(recovery_multiply_logs),
        "multiply_log_frequency": histogram(multiply_logs),
        "input_multiply_gate_total": sum(
            gate_counts[value] for value in input_multiply_logs),
        "recovery_multiply_gate_total": sum(
            gate_counts[value] for value in recovery_multiply_logs),
        "multiply_gate_total": sum(gate_counts[value] for value in multiply_logs),
        "input_multiply_depth_sum": sum(
            depths[value] for value in input_multiply_logs),
        "recovery_multiply_depth_sum": sum(
            depths[value] for value in recovery_multiply_logs),
        "multiply_depth_sum": sum(depths[value] for value in multiply_logs),
        "ifft_skew_frequency": histogram(ifft_skews),
        "fft_skew_frequency": histogram(fft_skews),
        "ifft_four_buffer_tuple_frequency": tuple_histogram(ifft_tuples),
        "fft_four_buffer_tuple_frequency": tuple_histogram(fft_tuples),
        "four_buffer_tuple_frequency": tuple_histogram(ifft_tuples + fft_tuples),
        "payload_memory": payload_model(operation_counts, buffer_bytes),
    }


def generate_corpus():
    circuit_data = circuits.build_circuits()
    multiply_circuits = circuit_data["multiply_circuits"]
    gate_counts = [len(gates) for gates in multiply_circuits]
    depths = [circuits.circuit_depth(gates, circuits.WIRE_COUNT_MULTIPLY)
              for gates in multiply_circuits]

    records = []
    for original_count, recovery_count in FULL_COUNTS:
        for buffer_bytes in FULL_BUFFER_BYTES:
            records.append(encode_record(
                original_count, recovery_count, buffer_bytes))
            for loss_count in loss_counts(recovery_count):
                records.append(decode_record(
                    original_count, recovery_count, buffer_bytes, loss_count,
                    gate_counts, depths))

    records_text = json.dumps(records, sort_keys=True, separators=(",", ":"))
    schedule_checksum = hashlib.sha256(
        records_text.encode("utf-8")).hexdigest()
    return {
        "schema": SCHEMA,
        "schedule_checksum_sha256": schedule_checksum,
        "circuit_checksum_sha256": circuits.circuit_checksum(circuit_data),
        "field_polynomial": "0x11d",
        "cantor_basis": list(circuits.CANTOR_BASIS),
        "error_bitfield_optimization": True,
        "benchmark_order": "counts, buffer bytes, encode, unique loss counts",
        "record_count": len(records),
        "records": records,
    }


def render_corpus(corpus):
    return json.dumps(corpus, sort_keys=True, indent=2) + "\n"


def check_output(output, expected):
    try:
        current = output.read_text(encoding="utf-8")
    except FileNotFoundError:
        print("schedule corpus is missing: %s" % output, file=sys.stderr)
        return False
    if current == expected:
        print("FF8 XOR schedule corpus is up to date: %s" % output)
        return True
    print("generated FF8 XOR schedule corpus is stale: %s" % output,
          file=sys.stderr)
    for line in difflib.unified_diff(
            current.splitlines(), expected.splitlines(),
            fromfile=str(output), tofile="regenerated", lineterm=""):
        print(line, file=sys.stderr)
    return False


def parse_arguments():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="fail when the checked-in corpus is stale")
    parser.add_argument(
        "--output", type=Path,
        default=root / "generated" / "FF8XorScheduleCorpus.json",
        help="output path (default: %(default)s)")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    corpus = generate_corpus()
    rendered = render_corpus(corpus)
    if arguments.check:
        return 0 if check_output(arguments.output, rendered) else 1
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered, encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Schedule checksum: %s" % corpus["schedule_checksum_sha256"])
    print("Records: %d" % corpus["record_count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
