#!/usr/bin/env python3
"""Finite independent tests for deterministic FF8 circuit costing."""

from __future__ import print_function

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as model  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    profile = model.checked(model.PORTABLE_DEFAULT_PROFILE)
    model.validate_profile(profile)
    require(profile["checksum_sha256"] == model.checksum(profile),
            "portable profile checksum is unstable")
    require(model.canonical_json({"b": 1, "a": 2}) == b'{"a":2,"b":1}',
            "canonical JSON ordering changed")

    bad = copy.deepcopy(profile)
    bad["weights"]["xor2_count"] += 1
    try:
        model.validate_profile(bad)
        raise AssertionError("stale profile checksum was accepted")
    except ValueError:
        pass

    for field, value in (
            ("id", None),
            ("id", "unsafe profile id"),
            ("id", "unsafe\"profile"),
            ("version", True),
            ("target_isa", ""),
            ("register_budget", True),
            ("xor3_supported", 1),
            ("weights", []),
            ("measurement", "not-an-object")):
        bad = copy.deepcopy(model.PORTABLE_DEFAULT_PROFILE)
        bad[field] = value
        try:
            model.validate_profile(bad)
            raise AssertionError("malformed profile %s was accepted" % field)
        except ValueError:
            pass
    bad = copy.deepcopy(model.PORTABLE_DEFAULT_PROFILE)
    bad["weights"]["xor2_count"] = True
    try:
        model.validate_profile(bad)
        raise AssertionError("boolean cost weight was accepted")
    except ValueError:
        pass

    for gates, width in (
            (((0, 1),), True),
            (((False, 1),), 4),
            (((0, True),), 4)):
        try:
            model.validate_gates(gates, width)
            raise AssertionError("boolean gate coordinate was accepted")
        except ValueError:
            pass

    gates = ((0, 1), (0, 2), (3, 0))
    metrics = model.source_metrics(gates, 4, profile)
    require(metrics["literal_xor2_count"] == 3,
            "wrong literal XOR2 count")
    require(metrics["xor2_count"] == 3,
            "wrong distinct-form XOR2 estimate")
    require(metrics["dependency_depth"] == 3, "wrong dependency depth")
    require(metrics["peak_live_wires"] == 4, "wrong peak liveness")
    require(metrics["loads"] == 4 and metrics["stores"] == 4,
            "wrong source memory model")
    require(metrics["code_bytes"] == 18 and metrics["icache_lines"] == 1,
            "wrong source code-footprint model")
    require(metrics["live_range_events"] > 0,
            "live-range model returned no work")

    census = {
        "code_bytes": 130,
        "icache_lines": 3,
        "counts": {
            "vector_xor2": 2,
            "vector_xor3": 1,
            "memory_loads": 7,
            "folded_vector_memory_reads": 2,
            "memory_stores": 8,
            "vector_stack_refs": 4,
            "loop_vector_stack_refs": 3,
        },
    }
    compiled = model.assembly_metrics(metrics, census)
    require(compiled["literal_xor2_count"] == 3 and
            compiled["xor2_count"] == 2 and compiled["xor3_count"] == 1,
            "assembly XOR counts were not applied")
    require(compiled["loads"] == 9 and compiled["stores"] == 8,
            "assembly memory counts were not applied")
    require(compiled["spill_references"] == 3,
            "loop spill count was not preferred")
    require(compiled["code_bytes"] == 130 and compiled["icache_lines"] == 3,
            "assembly footprint was not applied")
    boolean_metrics = dict(metrics)
    boolean_metrics["xor2_count"] = True
    try:
        model.validate_metrics(boolean_metrics)
        raise AssertionError("boolean circuit metric was accepted")
    except ValueError:
        pass

    shorter = ((0, 1),)
    longer = ((0, 1), (2, 3))
    require(model.circuit_key(shorter, 4, profile) <
            model.circuit_key(longer, 4, profile),
            "portable cost did not prefer fewer XORs")
    require(model.score_metrics(metrics, profile) > 0,
            "valid circuit received a non-positive score")

    require(abs(model.spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1) < 1e-12,
            "positive rank correlation is wrong")
    require(abs(model.spearman([1, 2, 3, 4], [40, 30, 20, 10]) + 1) < 1e-12,
            "negative rank correlation is wrong")
    require(model.tied_ranks([3, 1, 1, 2]) == [4.0, 1.5, 1.5, 3.0],
            "tie ranks are wrong")

    # JSON round trip is part of the checked-in profile contract.
    require(json.loads(json.dumps(profile, sort_keys=True)) == profile,
            "profile is not JSON round-trip stable")
    artifact = model.load_profile_artifact(
        ROOT / "generated" / "FF8XorCostProfiles.json")
    require(model.find_profile(artifact, "portable-default-v1")["id"] ==
            "portable-default-v1", "profile lookup returned the wrong entry")
    try:
        model.find_profile(artifact, "missing-profile")
        raise AssertionError("unknown profile id was accepted")
    except ValueError:
        pass
    print("FF8 XOR cost-model tests passed")


if __name__ == "__main__":
    main()
