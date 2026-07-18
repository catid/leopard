#!/usr/bin/env python3
"""Regression checks for checked-in FF8 ISA-aware cost profiles."""

from __future__ import print_function

import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as model  # noqa: E402
import generate_ff8xor_cost_profiles as generator  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def require_raises(function, message):
    try:
        function()
    except ValueError:
        return
    raise AssertionError(message)


def main():
    subprocess.check_call([
        sys.executable, str(ROOT / "tools" /
                            "generate_ff8xor_cost_profiles.py"), "--check"])
    artifact = json.loads((ROOT / "generated" /
                           "FF8XorCostProfiles.json").read_text(
                               encoding="utf-8"))
    require(artifact["schema"] == model.SCHEMA, "wrong profile schema")
    require(artifact["checksum_sha256"] == model.checksum(artifact),
            "top-level profile checksum mismatch")
    require(len(artifact["profiles"]) == 2,
            "portable and machine profiles are both required")
    profiles = {profile["id"]: profile for profile in artifact["profiles"]}
    portable = profiles["portable-default-v1"]
    machine = profiles["amd-zen5-gcc13-avx2-avx512-v1"]
    model.validate_profile(portable)
    model.validate_profile(machine)
    require(portable["weights"] ==
            model.PORTABLE_DEFAULT_PROFILE["weights"],
            "checked-in portable profile differs from generator default")
    require(not portable["measurement"]["runtime_calibration_required"] and
            not machine["measurement"]["runtime_calibration_required"],
            "a cost profile unexpectedly requires runtime calibration")
    require(not portable["xor3_supported"] and machine["xor3_supported"],
            "profile XOR3 capabilities are wrong")
    require(portable["selection_policy"] ==
            "gate-depth-lexical-then-modeled-cost" and
            machine["selection_policy"] == "modeled-cost",
            "portable/machine selection policies are wrong")

    fit = machine["measurement"]["fit"]
    require(fit["profile_spearman"] > fit["raw_gate_spearman"],
            "ISA-aware ranking did not beat raw gates")
    require(fit["profile_spearman"] - fit["raw_gate_spearman"] > 0.10,
            "rank-correlation improvement is too small")
    for fold in fit["coefficient_modulo_five_slices"]:
        require(fold["profile_spearman"] > fold["raw_gate_spearman"],
                "ISA-aware ranking lost a deterministic slice")

    workload = artifact["workload"]
    require(workload["record_count"] == 104,
            "cost workload does not cover the complete corpus")
    require(sum(workload["multiply_log_weight"].values()) > 0 and
            sum(workload["fft_skew_weight"].values()) > 0 and
            sum(workload["ifft_skew_weight"].values()) > 0,
            "cost workload omitted an operation family")
    require(workload["four_buffer_tuple_total_weight"] > 0 and
            workload["four_buffer_tuple_unique_count"] > 0,
            "cost workload omitted tuple frequency")

    corpus_path = ROOT / "generated" / "FF8XorScheduleCorpus.json"
    corpus = generator.read_schedule_corpus(corpus_path)
    require(corpus["record_count"] == len(corpus["records"]) == 104,
            "validated schedule corpus is incomplete")
    with tempfile.TemporaryDirectory(prefix="leopard-cost-corpus-test-") \
            as temporary:
        stale = json.loads(json.dumps(corpus))
        stale["record_count"] -= 1
        stale_path = Path(temporary) / "stale-count.json"
        stale_path.write_text(json.dumps(stale), encoding="utf-8")
        require_raises(
            lambda: generator.read_schedule_corpus(stale_path),
            "stale schedule record count was accepted")
        stale = json.loads(json.dumps(corpus))
        stale["records"][0]["buffer_bytes"] += 64
        stale_path = Path(temporary) / "stale-checksum.json"
        stale_path.write_text(json.dumps(stale), encoding="utf-8")
        require_raises(
            lambda: generator.read_schedule_corpus(stale_path),
            "stale schedule checksum was accepted")

    evaluation = artifact["selection_evaluation"]
    require(evaluation["default_profile"] == portable["id"] and
            not evaluation["machine_profile_enabled_by_default"],
            "machine-specific calibration became an implicit default")
    require(evaluation["measured_end_to_end_benefit_percent"] < 0,
            "rejected machine corpus lost its negative codec result")
    candidate = evaluation["candidate_evaluation"]
    require(candidate["compiler_assembly_strict_pass"] and
            candidate["compiler_assembly_spill_check_pass"],
            "candidate assembly evidence failed its strict checks")
    require(candidate["benchmark"]["repetitions"] == 7 and
            candidate["benchmark"]["case_count"] == 14 and
            candidate["benchmark"]["observation_count"] == 98 and
            candidate["benchmark"][
                "normalized_geomean_case_medians"] < 1,
            "candidate codec evidence is incomplete or no longer negative")
    pairs = evaluation["equivalent_pair_evaluation"]
    require(pairs["pair_count"] == 362 and
            pairs["median_win_count"] < pairs["pair_count"] / 2 and
            pairs["median_geomean_speedup"] < 1 and
            pairs["workload_weighted_geomean_speedup"] < 1,
            "equivalent-schedule rejection evidence changed")
    require(pairs["predicted_saving_vs_measured_speedup_spearman"] >
            pairs["raw_gate_saving_vs_measured_speedup_spearman"],
            "modeled ranking no longer improves even slightly on raw gates")
    require(pairs["predicted_saving_vs_measured_speedup_spearman"] < 0 and
            "does not prove" in evaluation["correlation_scope_caveat"],
            "within-map negative correlation is not candidly documented")
    require(len(evaluation["rejected_variants"]) >= 3,
            "rejected cost variants were not recorded")
    print("FF8 XOR cost-profile tests passed")


if __name__ == "__main__":
    main()
