#!/usr/bin/env python3
"""Deterministic ISA-aware costing for generated FF8 XOR circuits.

This module is deliberately independent of the C++ payload implementation.
It turns a literal CNOT/XOR schedule plus optional compiler census data into a
small set of auditable features.  The normal circuit generator uses the
checked-in portable profile; machine profiles are an explicit offline choice
and never trigger runtime calibration.
"""

from __future__ import print_function

import hashlib
import json
import math
import re
from pathlib import Path


SCHEMA = "leopard.ff8xor.cost-profiles.v1"
CALIBRATION_SCHEMA = "leopard.ff8xor.cost-calibration.v1"
PROFILE_VERSION = 1
ICACHE_LINE_BYTES = 64

FEATURES = (
    "literal_xor2_count",
    "xor2_count",
    "xor3_count",
    "dependency_depth",
    "peak_live_wires",
    "live_range_events",
    "spill_references",
    "loads",
    "stores",
    "code_bytes",
    "icache_lines",
)


# Integer micro-costs keep profile provenance bit-for-bit deterministic on every
# Python implementation.  For the portable profile, the historical
# gate/depth/lexical key is authoritative and deliberately precedes this score;
# only an explicit machine profile lets modeled costs drive circuit selection.
PORTABLE_DEFAULT_PROFILE = {
    "id": "portable-default-v1",
    "version": PROFILE_VERSION,
    "target_isa": "portable-xor2",
    "selection_policy": "gate-depth-lexical-then-modeled-cost",
    "register_budget": 8,
    "xor3_supported": False,
    "weights": {
        "literal_xor2_count": 1000000,
        "xor2_count": 100,
        "xor3_count": 4000000,
        "dependency_depth": 10000,
        "peak_live_wires": 100,
        "live_range_events": 10,
        "spill_references": 8000000,
        "loads": 100,
        "stores": 100,
        "code_bytes": 1,
        "icache_lines": 100,
    },
    "measurement": {
        "kind": "portable architectural default",
        "runtime_calibration_required": False,
        "notes": (
            "Historical gate/depth/lexical selection is authoritative; the "
            "modeled score is provenance-only. XOR3 is conservatively "
            "unavailable."),
    },
}


def canonical_json(value):
    """Return the stable byte representation used by every profile checksum."""
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")


def checksum(value):
    """Hash *value*, excluding a top-level checksum field if present."""
    if isinstance(value, dict) and "checksum_sha256" in value:
        value = dict(value)
        del value["checksum_sha256"]
    return hashlib.sha256(canonical_json(value)).hexdigest()


def checked(value):
    """Return a shallow copy carrying its deterministic SHA-256 checksum."""
    result = dict(value)
    result["checksum_sha256"] = checksum(result)
    return result


def validate_profile(profile):
    if not isinstance(profile, dict):
        raise ValueError("cost profile must be an object")
    for field in ("id", "version", "target_isa", "selection_policy",
                  "register_budget",
                  "xor3_supported", "weights", "measurement"):
        if field not in profile:
            raise ValueError("cost profile is missing %s" % field)
    if (not isinstance(profile["version"], int) or
            isinstance(profile["version"], bool) or
            profile["version"] != PROFILE_VERSION):
        raise ValueError("unsupported cost profile version")
    for field in ("id", "target_isa", "selection_policy"):
        if not isinstance(profile[field], str) or not profile[field]:
            raise ValueError("cost profile %s must be a nonempty string" %
                             field)
    if re.match(r"^[A-Za-z0-9_.+-]+$", profile["id"]) is None:
        raise ValueError("cost profile id contains unsafe characters")
    if profile["selection_policy"] not in (
            "gate-depth-lexical-then-modeled-cost", "modeled-cost"):
        raise ValueError("unsupported cost-profile selection policy")
    if (not isinstance(profile["register_budget"], int) or
            isinstance(profile["register_budget"], bool) or
            profile["register_budget"] <= 0):
        raise ValueError("register budget must be a positive integer")
    if not isinstance(profile["xor3_supported"], bool):
        raise ValueError("xor3_supported must be boolean")
    weights = profile["weights"]
    if not isinstance(weights, dict):
        raise ValueError("cost profile weights must be an object")
    if set(weights) != set(FEATURES):
        raise ValueError("cost profile weights do not match feature schema")
    if any(not isinstance(weights[name], int) or
           isinstance(weights[name], bool) or weights[name] < 0
           for name in FEATURES):
        raise ValueError("cost weights must be non-negative integers")
    if not isinstance(profile["measurement"], dict):
        raise ValueError("cost profile measurement must be an object")
    expected = profile.get("checksum_sha256")
    if expected is not None and expected != checksum(profile):
        raise ValueError("cost profile checksum mismatch")
    return profile


def load_profile_artifact(path):
    """Load and checksum a generated profile artifact."""
    path = Path(path)
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema") != SCHEMA:
        raise ValueError("unsupported cost-profile artifact schema")
    if artifact.get("checksum_sha256") != checksum(artifact):
        raise ValueError("cost-profile artifact checksum mismatch")
    profiles = artifact.get("profiles")
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("cost-profile artifact has no profiles")
    identifiers = set()
    for profile in profiles:
        validate_profile(profile)
        if profile["id"] in identifiers:
            raise ValueError("duplicate cost-profile id")
        identifiers.add(profile["id"])
    return artifact


def find_profile(artifact, profile_id):
    """Return one validated profile by stable id."""
    matches = [profile for profile in artifact["profiles"]
               if profile["id"] == profile_id]
    if len(matches) != 1:
        raise ValueError("unknown cost-profile id: %s" % profile_id)
    return validate_profile(matches[0])


def validate_gates(gates, width):
    if (not isinstance(width, int) or isinstance(width, bool) or
            width <= 0):
        raise ValueError("wire width must be a positive integer")
    result = tuple(tuple(gate) for gate in gates)
    for gate in result:
        if len(gate) != 2:
            raise ValueError("CNOT gate must contain destination and source")
        destination, source = gate
        if (not isinstance(destination, int) or
                isinstance(destination, bool) or
                not isinstance(source, int) or isinstance(source, bool) or
                not (0 <= destination < width) or
                not (0 <= source < width) or destination == source):
            raise ValueError("invalid CNOT gate")
    return result


def dependency_depth(gates, width):
    """Return the conservative named-wire dependency depth."""
    last_layer = [0] * width
    maximum = 0
    for destination, source in validate_gates(gates, width):
        layer = max(last_layer[destination], last_layer[source]) + 1
        last_layer[destination] = layer
        last_layer[source] = layer
        maximum = max(maximum, layer)
    return maximum


def live_range_events(gates, width):
    """Estimate SSA live-range pressure from literal wire uses.

    Every in-place assignment defines a new form for its destination.  Count
    the span to the next use/definition of each form.  This is not a register
    allocator simulation; it is a stable discriminator for schedules with the
    same instruction count and depth.
    """
    validated = validate_gates(gates, width)
    definitions = [0] * width
    total = 0
    for position, (destination, source) in enumerate(validated, 1):
        total += position - definitions[source]
        total += position - definitions[destination]
        definitions[destination] = position
    end = len(validated) + 1
    total += sum(end - definition for definition in definitions)
    return total


def distinct_linear_form_count(gates, width):
    """Estimate optimized XOR2 work as distinct non-input SSA forms.

    This is the same deterministic predictor used by the generator's earlier
    exact-synthesis promotion guard.  Unlike post-compile instruction counts,
    it is available for every candidate before selection and therefore
    transfers honestly into the explicit machine-profile CLI path.
    """
    wires = [1 << wire for wire in range(width)]
    forms = set(wires)
    for destination, source in validate_gates(gates, width):
        wires[destination] ^= wires[source]
        forms.add(wires[destination])
    return len(forms) - width


def source_metrics(gates, width, profile=None):
    """Return portable features for a literal XOR2 circuit."""
    profile = validate_profile(profile or PORTABLE_DEFAULT_PROFILE)
    validated = validate_gates(gates, width)
    code_bytes = 6 * len(validated)
    return {
        "literal_xor2_count": len(validated),
        "xor2_count": distinct_linear_form_count(validated, width),
        "xor3_count": 0,
        "dependency_depth": dependency_depth(validated, width),
        "peak_live_wires": width,
        "live_range_events": live_range_events(validated, width),
        "spill_references": max(0, width - profile["register_budget"]),
        "loads": width,
        "stores": width,
        "code_bytes": code_bytes,
        "icache_lines": ((code_bytes + ICACHE_LINE_BYTES - 1) //
                         ICACHE_LINE_BYTES),
    }


def assembly_metrics(source, assembly):
    """Overlay measured compiler shape on portable source metrics."""
    result = dict(source)
    counts = assembly.get("counts", {})
    result.update({
        "literal_xor2_count": result["literal_xor2_count"],
        "xor2_count": int(counts.get("vector_xor2", result["xor2_count"])),
        "xor3_count": int(counts.get("vector_xor3", 0)),
        "spill_references": int(counts.get(
            "loop_vector_stack_refs", counts.get(
                "vector_stack_refs", result["spill_references"]))),
        "loads": int(counts.get("memory_loads", result["loads"])) + int(
            counts.get("folded_vector_memory_reads", 0)),
        "stores": int(counts.get("memory_stores", result["stores"])),
        "code_bytes": int(assembly.get("code_bytes", result["code_bytes"])),
    })
    if "icache_lines" in assembly:
        result["icache_lines"] = int(assembly["icache_lines"])
    else:
        result["icache_lines"] = int(math.ceil(
            float(result["code_bytes"]) / ICACHE_LINE_BYTES))
    return result


def validate_metrics(metrics):
    if set(metrics) != set(FEATURES):
        raise ValueError("circuit metrics do not match feature schema")
    for name in FEATURES:
        if (not isinstance(metrics[name], int) or
                isinstance(metrics[name], bool) or metrics[name] < 0):
            raise ValueError("metric %s must be a non-negative integer" % name)
    return metrics


def score_metrics(metrics, profile=None):
    """Return a deterministic integer score; lower is better."""
    profile = validate_profile(profile or PORTABLE_DEFAULT_PROFILE)
    validate_metrics(metrics)
    return sum(metrics[name] * profile["weights"][name]
               for name in FEATURES)


def circuit_key(gates, width, profile=None):
    """Return an ISA-aware stable preference key for equivalent circuits."""
    validated = validate_gates(gates, width)
    profile = validate_profile(profile or PORTABLE_DEFAULT_PROFILE)
    metrics = source_metrics(validated, width, profile)
    score = score_metrics(metrics, profile)
    if profile["selection_policy"] == \
            "gate-depth-lexical-then-modeled-cost":
        # This is an intentional zero-output migration for the checked-in
        # portable backend.  The profile is consumed and its score remains in
        # the key/provenance, while the historical gate/depth/lexical safety
        # order remains byte-for-byte stable.  Explicit offline machine
        # profiles use modeled-cost first.
        return (metrics["literal_xor2_count"], metrics["dependency_depth"],
                validated, score)
    return (score, metrics["literal_xor2_count"],
            metrics["dependency_depth"], validated)


def tied_ranks(values):
    """Return one-based average ranks, with deterministic tie handling."""
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    result = [0.0] * len(values)
    begin = 0
    while begin < len(order):
        end = begin + 1
        while end < len(order) and values[order[end]] == values[order[begin]]:
            end += 1
        rank = (begin + 1 + end) / 2.0
        for offset in range(begin, end):
            result[order[offset]] = rank
        begin = end
    return result


def spearman(values_a, values_b):
    """Return Spearman rank correlation without third-party dependencies."""
    if len(values_a) != len(values_b) or len(values_a) < 2:
        raise ValueError("correlation requires equal vectors of length >= 2")
    ranks_a = tied_ranks(values_a)
    ranks_b = tied_ranks(values_b)
    mean_a = sum(ranks_a) / len(ranks_a)
    mean_b = sum(ranks_b) / len(ranks_b)
    covariance = sum((a - mean_a) * (b - mean_b)
                     for a, b in zip(ranks_a, ranks_b))
    variance_a = sum((a - mean_a) ** 2 for a in ranks_a)
    variance_b = sum((b - mean_b) ** 2 for b in ranks_b)
    if variance_a == 0 or variance_b == 0:
        return 0.0
    return covariance / math.sqrt(variance_a * variance_b)
