#!/usr/bin/env python3
"""Verify the rejected decoder locator-boundary fusion experiment.

The checked evidence is deliberately negative: a bounded identity/zero-only
first-IFFT fast path was exact, but its 768 static specializations nearly
doubled the core circuit text and did not improve the paired codec benchmark.
This checker preserves the independently reproducible algebra, schedule
frequency, evidence checksum, and aggregate timing arithmetic without keeping
the rejected payload code in the library.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import generate_ff8_xor_circuits as circuits  # noqa: E402
import generate_ff8xor_schedule_corpus as schedule  # noqa: E402


SCHEMA = "leopard.ff8xor.locator-boundary-evaluation.v1"
DEFAULT_EVIDENCE = (
    ROOT / "tools" / "profiles" /
    "FF8XorLocatorBoundaryEvaluationZen5Gcc13.json")
GENERATED_CIRCUITS = ROOT / "generated" / "LeopardFF8XorCircuits.inl"
SCHEDULE_CORPUS = ROOT / "generated" / "FF8XorScheduleCorpus.json"


def canonical_json(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")


def checksum(value):
    value = dict(value)
    value.pop("checksum_sha256", None)
    return hashlib.sha256(canonical_json(value)).hexdigest()


def load_evidence(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema") != SCHEMA:
        raise ValueError("unsupported locator-boundary evidence schema")
    if value.get("checksum_sha256") != checksum(value):
        raise ValueError("locator-boundary evidence checksum mismatch")
    if value.get("decision", {}).get("result") != "rejected":
        raise ValueError("the measured boundary candidate must remain rejected")
    if value.get("candidate", {}).get("retained") is not False:
        raise ValueError("rejected boundary candidate is marked retained")
    return value


def generated_array(name):
    text = GENERATED_CIRCUITS.read_text(encoding="utf-8")
    match = re.search(
        r"static const uint8_t %s\[256\] = \{(.*?)\n\};" %
        re.escape(name), text, re.DOTALL)
    if match is None:
        raise ValueError("generated circuit metadata is missing %s" % name)
    values = [int(token) for token in re.findall(r"\d+", match.group(1))]
    if len(values) != 256:
        raise ValueError("generated %s has %d entries" % (name, len(values)))
    return values


def verify_boundary_maps():
    """Compare independent staged and direct fused linear maps."""
    multiply_matrices = tuple(
        circuits.multiplication_matrix(logarithm)
        for logarithm in range(circuits.FIELD_ORDER))
    checked = 0
    for skew in range(circuits.FIELD_ORDER):
        for presence in (1, 2, 3):
            def staged(state):
                masked = 0
                if presence & 1:
                    masked |= state & 0xff
                if presence & 2:
                    masked |= state & 0xff00
                return circuits.scalar_inverse_butterfly(
                    masked, skew, multiply_matrices)

            def direct(state):
                x = (state & 0xff) if presence & 1 else 0
                y = ((state >> 8) & 0xff) if presence & 2 else 0
                y ^= x
                if skew != circuits.FIELD_MODULUS:
                    x ^= circuits.independent_scalar_multiply_log(y, skew)
                return x | (y << 8)

            staged_rows = circuits.linear_map(
                circuits.WIRE_COUNT_BUTTERFLY, staged)
            direct_rows = circuits.linear_map(
                circuits.WIRE_COUNT_BUTTERFLY, direct)
            if staged_rows != direct_rows:
                raise AssertionError(
                    "boundary map mismatch for skew=%d presence=%d" %
                    (skew, presence))
            checked += circuits.WIRE_COUNT_BUTTERFLY
    return checked


def automatic_schedule_cardinality(gate_counts, depths):
    eligible_pairs = 0
    eligible_skews = set()
    decode_records = 0
    for original_count, recovery_count in schedule.FULL_COUNTS:
        for buffer_bytes in schedule.FULL_BUFFER_BYTES:
            for loss_count in schedule.loss_counts(recovery_count):
                m = schedule.next_power_of_two(recovery_count)
                n = schedule.next_power_of_two(m + original_count)
                seed = (
                    0xFF8C000000000000 ^ (original_count << 40) ^
                    (recovery_count << 24) ^ buffer_bytes) & schedule.MASK64
                loss_seed = seed ^ (loss_count << 8) ^ 0xD3C0DE
                original_order = schedule.shuffled_indices(
                    original_count, loss_seed)
                recovery_order = schedule.shuffled_indices(
                    recovery_count,
                    loss_seed ^ 0x7265636F76657279)
                original_missing = frozenset(
                    original_order[:loss_count])
                recovery_available = frozenset(
                    recovery_order[:loss_count])
                locations = schedule.locator_logs(
                    original_count, recovery_count, m,
                    original_missing, recovery_available)
                shift = schedule.select_locator_shift(
                    locations, original_count, recovery_count, m,
                    original_missing, recovery_available,
                    gate_counts, depths)

                present = [False] * n
                logarithms = [0] * n
                for index in recovery_available:
                    present[index] = True
                    logarithms[index] = schedule.shifted_log(
                        locations[index], shift)
                for index in range(original_count):
                    if index not in original_missing:
                        position = m + index
                        present[position] = True
                        logarithms[position] = schedule.shifted_log(
                            locations[position], shift)

                for range_start in range(0, n, 2):
                    if range_start >= m + original_count:
                        continue
                    x_present = present[range_start]
                    y_present = present[range_start + 1]
                    if not x_present and not y_present:
                        continue
                    skew = schedule.FFT_SKEW_PADDED[range_start + 1]
                    # Existing x-zero/y-present sentinel handling moves only
                    # 2B; an out-of-place two-output kernel moves 3B.
                    if not x_present and y_present and skew == 255:
                        continue
                    if ((not x_present or logarithms[range_start] == 0) and
                            (not y_present or
                             logarithms[range_start + 1] == 0)):
                        eligible_pairs += 1
                        eligible_skews.add(skew)
                decode_records += 1
    return decode_records, eligible_pairs, eligible_skews


def geometric_mean(values):
    if not values or any(value <= 0 for value in values):
        raise ValueError("ratios must be a nonempty positive sequence")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def verify_evidence(evidence):
    proof_count = verify_boundary_maps()
    expected_proof = 256 * 3 * 16
    if proof_count != expected_proof:
        raise AssertionError("unexpected boundary linear-map proof count")

    gate_counts = generated_array("kMultiplyGateCounts")
    depths = generated_array("kMultiplyDepths")
    records, pairs, skews = automatic_schedule_cardinality(
        gate_counts, depths)
    corpus = json.loads(SCHEDULE_CORPUS.read_text(encoding="utf-8"))
    if evidence["schedule_corpus"] != {
            "automatic_decode_records": records,
            "eligible_first_layer_pairs": pairs,
            "eligible_unique_skews": len(skews),
            "possible_first_layer_skews": len(set(
                schedule.FFT_SKEW_PADDED[index + 1]
                for index in range(0, 256, 2))),
            "schedule_checksum_sha256":
                corpus["schedule_checksum_sha256"]}:
        raise ValueError("locator-boundary schedule evidence is stale")

    candidate = evidence["candidate"]
    if candidate["dense_first_layer_map_cardinality"] != 255 * 255 * 256:
        raise ValueError("dense boundary-map cardinality is stale")
    if candidate["generated_presence_skew_specializations"] != 3 * 256:
        raise ValueError("bounded specialization count is stale")

    benchmark = evidence["benchmark"]
    if benchmark["eligible_ratio_aggregation"] != (
            "equal weight per retained benchmark row and independent "
            "repetition; not schedule-frequency weighted"):
        raise ValueError("benchmark ratio aggregation is ambiguous")
    if benchmark["measurement_order"] != (
            "paired A-B-B-A; each sample is the average of its two "
            "same-arm calls"):
        raise ValueError("benchmark does not describe its paired averages")
    if benchmark["allocation_timed"] or benchmark["transpose_included"]:
        raise ValueError("boundary benchmark included out-of-scope work")
    if benchmark["work_buffer_placement"] != (
            "same physical work buffers for both arms"):
        raise ValueError("boundary benchmark arms used different placements")

    eligible_ratios = []
    for row in benchmark["rows"]:
        ratios = row["ratios_staged_over_fused"]
        staged_usec = row["staged_median_usec"]
        fused_usec = row["fused_median_usec"]
        if (len(ratios) != benchmark["repetitions"] or
                len(staged_usec) != len(ratios) or
                len(fused_usec) != len(ratios)):
            raise ValueError("benchmark row repetition count mismatch")
        for index, recorded_ratio in enumerate(ratios):
            measured_ratio = staged_usec[index] / fused_usec[index]
            if abs(measured_ratio - recorded_ratio) > 7e-7:
                raise ValueError("benchmark row ratio is inconsistent")
        geometric_mean(ratios)
        if row["eligible_pairs"]:
            eligible_ratios.extend(ratios)
    measured = geometric_mean(eligible_ratios)
    recorded = evidence["decision"]["eligible_row_geometric_mean_ratio"]
    if abs(measured - recorded) > 5e-11:
        raise ValueError("eligible-row aggregate ratio is stale")

    small = benchmark["rows"][0]
    small_measured = geometric_mean(small["ratios_staged_over_fused"])
    small_recorded = evidence["decision"][
        "small_transform_geometric_mean_ratio"]
    if abs(small_measured - small_recorded) > 5e-7:
        raise ValueError("small-transform aggregate ratio is stale")

    size = evidence["code_size"]
    if (size["candidate_core_object_text_bytes"] -
            size["baseline_core_object_text_bytes"] !=
            size["delta_text_bytes"]):
        raise ValueError("boundary code-size delta is inconsistent")
    ratio = (float(size["candidate_core_object_text_bytes"]) /
             size["baseline_core_object_text_bytes"])
    if abs(ratio - size["ratio_candidate_over_baseline"]) > 5e-11:
        raise ValueError("boundary code-size ratio is inconsistent")

    return proof_count, records, pairs, len(skews), measured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="verify checked evidence (also the default)")
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    try:
        evidence = load_evidence(args.evidence)
        proof, records, pairs, skew_count, ratio = verify_evidence(evidence)
    except (AssertionError, KeyError, OSError, ValueError) as error:
        print("FF8 XOR locator-boundary evidence failed: %s" % error,
              file=sys.stderr)
        return 1
    print(
        "FF8 XOR locator-boundary evidence is valid: proof=%d "
        "decode_records=%d eligible_pairs=%d eligible_skews=%d "
        "eligible_ratio=%.9f decision=rejected" %
        (proof, records, pairs, skew_count, ratio))
    return 0


if __name__ == "__main__":
    sys.exit(main())
