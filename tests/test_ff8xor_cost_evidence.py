#!/usr/bin/env python3
"""Validate checked ISA-cost measurements and their provenance links."""

from __future__ import print_function

import json
import hashlib
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as model  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def checked(path, schema):
    artifact = json.loads(path.read_text(encoding="utf-8"))
    require(artifact.get("schema") == schema,
            "%s has the wrong schema" % path.name)
    require(artifact.get("checksum_sha256") == model.checksum(artifact),
            "%s has a stale checksum" % path.name)
    return artifact


def main():
    profiles = model.load_profile_artifact(
        ROOT / "generated" / "FF8XorCostProfiles.json")
    machine = model.find_profile(
        profiles, "amd-zen5-gcc13-avx2-avx512-v1")
    profile_evaluation = profiles["selection_evaluation"]
    calibration = checked(
        ROOT / "tools" / "profiles" /
        "FF8XorCostCalibrationZen5Gcc13Avx2.json",
        model.CALIBRATION_SCHEMA)
    xor3 = checked(
        ROOT / "tools" / "profiles" /
        "FF8XorXor3CalibrationZen5Gcc13.json",
        "leopard.ff8xor.xor3-calibration.v1")
    candidate = checked(
        ROOT / "tools" / "profiles" /
        "FF8XorCostCandidateEvaluationZen5Gcc13.json",
        "leopard.ff8xor.cost-candidate-evaluation.v1")
    pairs = checked(
        ROOT / "tools" / "profiles" /
        "FF8XorCostEquivalentPairCalibrationZen5Gcc13Avx2.json",
        "leopard.ff8xor.equivalent-pair-calibration.v1")

    require(len(calibration["records"]) == 256 and
            calibration["rounds"] >= 3 and
            calibration["circuit_iterations_per_sample"] > 1,
            "multiplier calibration is incomplete")
    require(calibration["checksum_sha256"] == machine["measurement"][
                "circuit_calibration_checksum_sha256"] and
            xor3["checksum_sha256"] == machine["measurement"][
                "xor3_calibration_checksum_sha256"],
            "machine profile is not bound to its raw calibrations")
    require(any(record["width"] == "ymm" and
                record["shape"] == "four-independent-chains" and
                record["speedup"] > 1
                for record in xor3["profiles"]),
            "XOR3 issue-throughput evidence is missing")

    for artifact in (candidate, pairs):
        require(artifact["machine_profile_checksum_sha256"] ==
                machine["checksum_sha256"],
                "measurement uses a different machine profile")
        require(artifact["workload_schedule_checksum_sha256"] ==
                profiles["workload"]["schedule_checksum_sha256"],
                "measurement uses a different schedule corpus")
    require(candidate["benchmark"]["repetitions"] == 7 and
            candidate["benchmark"]["case_count"] == 14 and
            candidate["benchmark"]["observation_count"] == 98,
            "candidate codec benchmark is incomplete")
    require(candidate["benchmark"]["input_schema"] ==
            "leopard.ff8xor.benchmark.jsonl.v1" and
            "before" in candidate["reproduction"]["benchmark_input_note"],
            "historical benchmark schema provenance is unclear")
    require(not candidate["decision"]["enable_machine_profile_by_default"] and
            candidate["decision"]["measured_codec_speedup"] < 1,
            "negative codec candidate was not rejected")
    for assembly in (candidate["portable_assembly"],
                     candidate["machine_assembly"]):
        require(assembly["strict_pass"] and assembly["spill_check_pass"] and
                assembly["hard_violation_count"] == 0 and
                assembly["spill_warning_count"] == 0,
                "candidate assembly failed strict/no-spill inspection")
    require(all(candidate["weighted_comparison"][family][
                "avx2_xor2_weighted_delta_percent"] < 0
                for family in ("multiply", "fft", "ifft")),
            "compiled AVX2 tradeoff evidence changed unexpectedly")

    require(pairs["pair_count"] == len(pairs["records"]) == 362 and
            pairs["rounds"] >= 3 and pairs["base_iterations"] > 1 and
            pairs["pinned_cpu"] >= 0,
            "equivalent-pair calibration is incomplete")
    require(pairs["semantic_validation"] ==
            "portable and machine outputs compared before timing",
            "equivalent circuits were not validated before timing")
    identities = {(record["family"], record["coefficient"])
                  for record in pairs["records"]}
    require(len(identities) == pairs["pair_count"],
            "equivalent-pair calibration contains duplicates")
    expected_counts = {"multiply": 51, "fft": 176, "ifft": 135}
    require({family: pairs["by_family"][family]["pair_count"]
             for family in expected_counts} == expected_counts,
            "equivalent-pair family coverage changed")
    require(all(record["workload_weight"] >= 0 and
                math.isfinite(record["measured_speedup"]) and
                record["measured_speedup"] > 0
                for record in pairs["records"]),
            "equivalent-pair record has invalid timing/weight")
    require(pairs["workload_weighted_geomean_speedup"] < 1 and
            pairs["median_geomean_speedup"] < 1 and
            pairs["predicted_saving_vs_measured_speedup_spearman"] < 0,
            "rejected equivalent-pair result is no longer negative")

    require(profile_evaluation["candidate_evaluation"]["checksum_sha256"] ==
            candidate["checksum_sha256"] and
            profile_evaluation["equivalent_pair_evaluation"][
                "checksum_sha256"] == pairs["checksum_sha256"],
            "generated profile artifact is not bound to raw evaluations")

    # Bind the retained timing records to the exact generated programs they
    # measured.  The rejected machine corpus is not built by default, so its
    # checked candidate summary is the durable second side of this link.
    portable_path = ROOT / "generated" / "LeopardFF8XorCircuits.inl"
    portable_text = portable_path.read_text(encoding="utf-8")
    portable_file_checksum = hashlib.sha256(
        portable_path.read_bytes()).hexdigest()
    portable_circuit_match = re.search(
        r"Circuit checksum \(SHA-256\): ([0-9a-f]{64})", portable_text)
    require(portable_circuit_match is not None,
            "current portable circuit checksum is missing")
    require(pairs["portable_circuit_file_sha256"] ==
            portable_file_checksum,
            "equivalent-pair timings use a stale portable circuit file")
    require(candidate["portable_source"]["circuit_checksum_sha256"] ==
            portable_circuit_match.group(1),
            "candidate evaluation uses a different portable circuit map")
    require(pairs["machine_circuit_file_sha256"] ==
            candidate["machine_source"]["file_sha256"],
            "machine pair timings and codec evaluation use different circuits")
    print("FF8 XOR cost-evidence provenance tests passed")


if __name__ == "__main__":
    main()
