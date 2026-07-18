#!/usr/bin/env python3
"""Generate deterministic portable and calibrated FF8 circuit-cost profiles."""

from __future__ import print_function

import argparse
import collections
import difflib
import hashlib
import itertools
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as model  # noqa: E402


SCHEDULE_CORPUS_SCHEMA = "leopard.ff8xor.schedule-corpus.v1"


def read_checked_json(path, expected_schema=None):
    data = json.loads(path.read_text(encoding="utf-8"))
    if expected_schema is not None and data.get("schema") != expected_schema:
        raise ValueError("unexpected schema in %s" % path)
    if data.get("checksum_sha256") != model.checksum(data):
        raise ValueError("checksum mismatch in %s" % path)
    return data


def read_schedule_corpus(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != SCHEDULE_CORPUS_SCHEMA:
        raise ValueError("unexpected schedule-corpus schema in %s" % path)
    records = data.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("schedule corpus in %s has no records" % path)
    if data.get("record_count") != len(records):
        raise ValueError("schedule-corpus record count mismatch in %s" % path)
    records_text = json.dumps(
        records, sort_keys=True, separators=(",", ":"))
    expected_checksum = hashlib.sha256(
        records_text.encode("utf-8")).hexdigest()
    if data.get("schedule_checksum_sha256") != expected_checksum:
        raise ValueError("schedule-corpus checksum mismatch in %s" % path)
    return data


def add_histogram(destination, histogram, scale=1):
    for key, value in histogram.items():
        destination[key] += int(value) * scale


def workload_model(corpus):
    multiply = collections.Counter()
    fft = collections.Counter()
    ifft = collections.Counter()
    tuples = collections.Counter()
    for record in corpus["records"]:
        # A whole-buffer invocation scales linearly with payload bytes.  The
        # corpus gives every benchmark row equal representation, while this
        # byte weight makes the 1 MiB rows contribute their actual work.
        scale = int(record["buffer_bytes"])
        add_histogram(multiply, record["multiply_log_frequency"], scale)
        add_histogram(fft, record["fft_skew_frequency"], scale)
        add_histogram(ifft, record["ifft_skew_frequency"], scale)
        add_histogram(tuples, record["four_buffer_tuple_frequency"], scale)

    def ordered(counter):
        return {key: counter[key] for key in sorted(
            counter, key=lambda value: tuple(
                int(part) for part in value.split(",")))}

    top_tuples = sorted(
        tuples.items(), key=lambda item: (-item[1], item[0]))[:32]
    return {
        "schedule_checksum_sha256": corpus["schedule_checksum_sha256"],
        "record_count": corpus["record_count"],
        "weight_unit": "whole-buffer invocation times buffer bytes",
        "multiply_log_weight": ordered(multiply),
        "fft_skew_weight": ordered(fft),
        "ifft_skew_weight": ordered(ifft),
        "four_buffer_tuple_total_weight": sum(tuples.values()),
        "four_buffer_tuple_unique_count": len(tuples),
        "top_four_buffer_tuples": [
            {"tuple": key, "weight": value} for key, value in top_tuples],
    }


def candidate_weights():
    """Yield the finite, checked-in deterministic offline tuning portfolio."""
    for literal, xor2, depth, code, icache, live in itertools.product(
            (1, 25, 50),
            (850, 900, 950, 1000),
            (400, 500, 600, 700),
            (5, 10, 15, 20),
            (50, 100, 150, 200),
            (1, 2)):
        yield {
            "literal_xor2_count": literal,
            "xor2_count": xor2,
            # Filled after selection from the independent XOR3 benchmark.
            "xor3_count": xor2,
            "dependency_depth": depth,
            "peak_live_wires": 20,
            "live_range_events": live,
            "spill_references": 500000,
            "loads": 50,
            "stores": 50,
            "code_bytes": code,
            "icache_lines": icache,
        }


def score_record(record, weights):
    # This must be the exact predictor available to generator circuit_key.
    # Post-compile metrics are retained as calibration diagnostics only.
    metrics = record["source_metrics"]
    return sum(metrics[name] * weights[name] for name in model.FEATURES)


def fit_machine_weights(evidence, xor3):
    records = [record for record in evidence["records"]
               if record["coefficient"] not in (0, 255)]
    measured = [record["timing"]["median_ns"] for record in records]
    raw_gates = [record["source_gate_count"] for record in records]
    raw_correlation = model.spearman(raw_gates, measured)

    best = None
    portfolio_count = 0
    for weights in candidate_weights():
        portfolio_count += 1
        predicted = [score_record(record, weights) for record in records]
        correlation = model.spearman(predicted, measured)
        key = (-correlation, tuple(weights[name] for name in model.FEATURES))
        if best is None or key < best[0]:
            best = (key, dict(weights), correlation)
    weights = best[1]

    # An XOR3 operation replaces two XOR2 operations.  The four-independent-
    # chain YMM result is the closest isolated issue-throughput measurement to
    # a wide generated circuit.  Round to an integer micro-cost.
    xor3_record = next(record for record in xor3["profiles"]
                       if record["width"] == "ymm" and
                       record["shape"] == "four-independent-chains")
    weights["xor3_count"] = int(round(
        2.0 * weights["xor2_count"] / xor3_record["speedup"]))
    predicted = [score_record(record, weights) for record in records]
    fitted_correlation = model.spearman(predicted, measured)
    compiled_xor_correlation = model.spearman(
        [record["compiled_metrics"]["xor2_count"] for record in records],
        measured)

    folds = []
    for fold in range(5):
        selected = [record for record in records
                    if record["coefficient"] % 5 == fold]
        selected_measured = [record["timing"]["median_ns"]
                             for record in selected]
        folds.append({
            "fold": fold,
            "record_count": len(selected),
            "raw_gate_spearman": model.spearman(
                [record["source_gate_count"] for record in selected],
                selected_measured),
            "profile_spearman": model.spearman(
                [score_record(record, weights) for record in selected],
                selected_measured),
        })

    return weights, {
        "portfolio_candidate_count": portfolio_count,
        "fit_records": len(records),
        "objective": "maximize Spearman rank correlation; lexical weight tie-break",
        "predictor_metrics": "source_metrics available before generation",
        "raw_gate_spearman": raw_correlation,
        "profile_spearman": fitted_correlation,
        "actual_compiled_xor_count_spearman_oracle":
            compiled_xor_correlation,
        "coefficient_modulo_five_slices": folds,
        "xor3_issue_weight_source": {
            "width": xor3_record["width"],
            "shape": xor3_record["shape"],
            "measured_speedup": xor3_record["speedup"],
        },
    }


def build_output(corpus, evidence, xor3, candidate, equivalent_pairs):
    if evidence["circuit_checksum_sha256"] != \
            corpus["circuit_checksum_sha256"]:
        raise ValueError("cost evidence and schedule corpus use different circuits")
    weights, fit = fit_machine_weights(evidence, xor3)
    portable = model.checked(model.PORTABLE_DEFAULT_PROFILE)
    machine = model.checked({
        "id": "amd-zen5-gcc13-avx2-avx512-v1",
        "version": model.PROFILE_VERSION,
        "target_isa": "avx2 plus avx512f/avx512vl XOR3",
        "selection_policy": "modeled-cost",
        "register_budget": 16,
        "xor3_supported": True,
        "weights": weights,
        "measurement": {
            "kind": "offline measured machine/compiler profile",
            "runtime_calibration_required": False,
            "cpu_model": evidence["cpu_model"],
            "compiler": evidence["compiler"],
            "compiler_flags": evidence["compiler_flags"],
            "circuit_calibration_checksum_sha256":
                evidence["checksum_sha256"],
            "xor3_calibration_checksum_sha256": xor3["checksum_sha256"],
            "fit": fit,
            "notes": (
                "Use only for explicit offline regeneration and revalidate "
                "compiled assembly plus end-to-end codec timing."),
        },
    })
    for label, artifact in (("candidate evaluation", candidate),
                            ("equivalent-pair calibration",
                             equivalent_pairs)):
        if artifact["machine_profile_checksum_sha256"] != \
                machine["checksum_sha256"]:
            raise ValueError("%s uses a different machine profile" % label)
        if artifact["workload_schedule_checksum_sha256"] != \
                corpus["schedule_checksum_sha256"]:
            raise ValueError("%s uses a different workload schedule" % label)

    codec_speedup = candidate["benchmark"][
        "normalized_geomean_case_medians"]
    codec_benefit_percent = (codec_speedup - 1.0) * 100.0
    pair_speedup = equivalent_pairs["workload_weighted_geomean_speedup"]
    result = {
        "schema": model.SCHEMA,
        "profiles": [portable, machine],
        "workload": workload_model(corpus),
        "selection_evaluation": {
            "default_profile": portable["id"],
            "machine_profile_enabled_by_default": False,
            "measured_end_to_end_benefit_percent": codec_benefit_percent,
            "measured_end_to_end_note": (
                "Seven alternating quick-benchmark repetitions (98 "
                "observations across 14 cases), normalized by the packed "
                "control, measured %.6fx.  This is a negative result and "
                "does not enable the machine-specific circuit corpus." %
                codec_speedup),
            "candidate_evaluation": {
                "checksum_sha256": candidate["checksum_sha256"],
                "compiler_assembly_strict_pass": (
                    candidate["portable_assembly"]["strict_pass"] and
                    candidate["machine_assembly"]["strict_pass"]),
                "compiler_assembly_spill_check_pass": (
                    candidate["portable_assembly"]["spill_check_pass"] and
                    candidate["machine_assembly"]["spill_check_pass"]),
                "weighted_comparison": candidate["weighted_comparison"],
                "benchmark": {
                    key: candidate["benchmark"][key] for key in (
                        "repetitions", "case_count", "observation_count",
                        "normalized_geomean_all_observations",
                        "normalized_geomean_case_medians",
                        "normalized_case_median_wins",
                        "per_repetition_median", "by_operation")},
            },
            "equivalent_pair_evaluation": {
                "checksum_sha256": equivalent_pairs["checksum_sha256"],
                "pair_count": equivalent_pairs["pair_count"],
                "median_win_count": equivalent_pairs["median_win_count"],
                "median_geomean_speedup":
                    equivalent_pairs["median_geomean_speedup"],
                "workload_weighted_geomean_speedup": pair_speedup,
                "raw_gate_saving_vs_measured_speedup_spearman":
                    equivalent_pairs[
                        "raw_gate_saving_vs_measured_speedup_spearman"],
                "predicted_saving_vs_measured_speedup_spearman":
                    equivalent_pairs[
                        "predicted_saving_vs_measured_speedup_spearman"],
                "by_family": equivalent_pairs["by_family"],
                "interpretation": (
                    "The same-process ABBA comparison validates equivalent "
                    "maps, but the selected alternatives were %.6fx "
                    "workload-weighted and the predictor correlation was "
                    "negative.  The modeled-cost candidate is rejected." %
                    pair_speedup),
            },
            "correlation_scope_caveat": (
                "The positive calibration correlation ranks different field "
                "coefficients; it does not prove schedule ranking among "
                "equivalent circuits for one map.  The direct within-map pair "
                "experiment is the relevant selection test and was negative."),
            "rejected_variants": [
                {
                    "variant": "raw source gate count",
                    "reason": "lower measured rank correlation",
                    "spearman": fit["raw_gate_spearman"],
                },
                {
                    "variant": "current machine-profile circuit corpus",
                    "reason": (
                        "compiled XOR/code-size improvements did not survive "
                        "isolated equivalent-pair or end-to-end codec timing"),
                    "codec_speedup": codec_speedup,
                    "equivalent_pair_workload_weighted_speedup": pair_speedup,
                    "equivalent_pair_predictor_spearman": equivalent_pairs[
                        "predicted_saving_vs_measured_speedup_spearman"],
                },
                {
                    "variant": "VPTERNLOG benefit in streaming kernels",
                    "reason": (
                        "checked-in stream ratios are approximately neutral; "
                        "register-only throughput is not end-to-end evidence"),
                },
            ],
        },
    }
    result["checksum_sha256"] = model.checksum(result)
    return result


def render(value):
    return json.dumps(value, sort_keys=True, indent=2) + "\n"


def check_output(path, expected):
    try:
        current = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        print("FF8 XOR cost profiles are missing: %s" % path, file=sys.stderr)
        return False
    if current == expected:
        print("FF8 XOR cost profiles are up to date: %s" % path)
        return True
    print("FF8 XOR cost profiles are stale: %s" % path, file=sys.stderr)
    for line in difflib.unified_diff(
            current.splitlines(), expected.splitlines(),
            fromfile=str(path), tofile="regenerated", lineterm=""):
        print(line, file=sys.stderr)
    return False


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--corpus", type=Path,
        default=ROOT / "generated" / "FF8XorScheduleCorpus.json")
    parser.add_argument(
        "--calibration", type=Path,
        default=ROOT / "tools" / "profiles" /
        "FF8XorCostCalibrationZen5Gcc13Avx2.json")
    parser.add_argument(
        "--xor3-calibration", type=Path,
        default=ROOT / "tools" / "profiles" /
        "FF8XorXor3CalibrationZen5Gcc13.json")
    parser.add_argument(
        "--candidate-evaluation", type=Path,
        default=ROOT / "tools" / "profiles" /
        "FF8XorCostCandidateEvaluationZen5Gcc13.json")
    parser.add_argument(
        "--equivalent-pair-calibration", type=Path,
        default=ROOT / "tools" / "profiles" /
        "FF8XorCostEquivalentPairCalibrationZen5Gcc13Avx2.json")
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "generated" / "FF8XorCostProfiles.json")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    corpus = read_schedule_corpus(arguments.corpus)
    evidence = read_checked_json(
        arguments.calibration, model.CALIBRATION_SCHEMA)
    xor3 = read_checked_json(arguments.xor3_calibration,
                             "leopard.ff8xor.xor3-calibration.v1")
    candidate = read_checked_json(
        arguments.candidate_evaluation,
        "leopard.ff8xor.cost-candidate-evaluation.v1")
    equivalent_pairs = read_checked_json(
        arguments.equivalent_pair_calibration,
        "leopard.ff8xor.equivalent-pair-calibration.v1")
    output = build_output(
        corpus, evidence, xor3, candidate, equivalent_pairs)
    text = render(output)
    if arguments.check:
        return 0 if check_output(arguments.output, text) else 1
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(text, encoding="utf-8")
    machine = output["profiles"][1]
    fit = machine["measurement"]["fit"]
    print("Generated %s" % arguments.output)
    print("Profile checksum: %s" % output["checksum_sha256"])
    print("Raw gates Spearman: %.6f" % fit["raw_gate_spearman"])
    print("ISA-aware Spearman: %.6f" % fit["profile_spearman"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
