#!/usr/bin/env python3
"""Audit whether external erasure-code results are fairly comparable.

This tool does not benchmark third-party code.  It classifies every job in a
deterministically regenerated Leopard2 matrix against explicit provider
constraints so an unsupported cell cannot be silently resized, moved to
another field, or reported as a wire-compatible comparison.  The matrix command
emits aggregate status/reason/qualification counts; ``classify`` returns the
per-cell detail used to build those aggregates.  ISA-L cells supported by the
bounded, default-off single-thread adapter are explicitly called
``adapter-available-unmeasured``: this audit never turns availability into a
throughput claim.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence


SCHEMA = "leopard2-external-comparison-audit/v1"
RETRIEVED = "2026-07-16"
# The future adapter uses a deterministic 64-bit-host contract instead of
# silently inheriting sizeof(long) from whichever machine runs this audit.
# Jerasure requires region sizes to be longword multiples and pointers to be
# longword aligned.  Eight-byte regions are safe for the intended 64-bit
# comparison builds and keep offline classifications reproducible.
JERASURE_ADAPTER_REGION_BYTES = 8
SOURCES = {
    "isa-l": {
        "name": "Intel ISA-L",
        "url": "https://github.com/intel/isa-l",
        "commit": "e8cc5e87fc64b4da434f32bc1fa18184622a4998",
        "license": "BSD-3-Clause",
        "kind": "codec",
        "coordinate_equivalence": "different systematic GF(256) matrix code",
        "adapter": {
            "status": "bounded-default-off",
            "entrypoint": "tools/leopard2_isal_compare.py",
            "generator": "gf_gen_cauchy1_matrix",
            "maximum_threads": 1,
            "production_dependency": False,
        },
    },
    "jerasure": {
        "name": "Jerasure 2.0",
        "url": "https://github.com/ceph/jerasure",
        "commit": "de1739cc8483696506829b52e7fda4f6bb195e6a",
        "license": "BSD-3-Clause",
        "kind": "codec",
        "coordinate_equivalence": "different systematic Vandermonde matrix code",
        "runtime_dependency": "GF-Complete",
    },
    "fastecc": {
        "name": "FastECC",
        "url": "https://github.com/Bulat-Ziganshin/FastECC",
        "commit": "b8ca7db6bf5556185c96009b161e8aec82af734e",
        "license": "Apache-2.0",
        "kind": "codec-research",
        "coordinate_equivalence": (
            "prime-field NTT code with parity symbols wider than input symbols"),
    },
    "ecc-benchmark": {
        "name": "ECC-Benchmark",
        "url": "https://github.com/Bulat-Ziganshin/ECC-Benchmark",
        "commit": "c43d4290f8525351821f7b04791cee3bdfbaccdd",
        "license": "MIT",
        "kind": "comparison-harness",
        "coordinate_equivalence": "not a codec or wire profile",
    },
}


class AuditError(ValueError):
    """Invalid matrix or source-cache evidence."""


def _ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def _binary_field_for_length(length: int) -> str:
    if length <= 256:
        return "gf8"
    if length <= 65536:
        return "gf16"
    return "unsupported"


def _resolved_identity(cell: Mapping[str, object]) -> tuple[str, str, int]:
    k = int(cell["K"])
    r = int(cell["R"])
    requested = str(cell["requested_profile"])
    if requested in ("low", "low_v1"):
        profile = "low_v1"
    elif requested in ("high", "legacy_high_v1"):
        profile = "legacy_high_v1"
    elif requested == "auto":
        profile = "low_v1" if r > k else "legacy_high_v1"
    else:
        raise AuditError(f"unsupported requested profile {requested!r}")
    if profile == "low_v1":
        parent = _ceil_pow2(_ceil_pow2(k) + r)
    else:
        parent = _ceil_pow2(k + _ceil_pow2(r))
    field = _binary_field_for_length(parent)
    return profile, field, parent


def classify(provider: str, cell: Mapping[str, object]) -> dict:
    if provider not in SOURCES:
        raise AuditError(f"unknown provider {provider!r}")
    k = int(cell["K"])
    r = int(cell["R"])
    shard_bytes = int(cell["shard_bytes"])
    profile, field, parent = _resolved_identity(cell)
    base = {
        "provider": provider,
        "K": k,
        "R": r,
        "shard_bytes": shard_bytes,
        "leopard2_profile": profile,
        "leopard2_field": field,
        "leopard2_parent": parent,
        "provider_field": None,
        "wire_compatible": False,
    }
    if provider == "isa-l":
        base["provider_field"] = "gf8" if k + r <= 256 else "unsupported"
        reasons = []
        if k + r > 256:
            reasons.append("public K+R exceeds the GF(256) evaluation-set bound")
        if shard_bytes > 0x7FFFFFFF:
            reasons.append("shard length exceeds ISA-L's signed-int API")
        thread_count = int(cell.get("thread_count", 1))
        adapter_gaps = []
        if thread_count != 1:
            adapter_gaps.append(
                "the bounded adapter has no persistent-pool implementation for "
                "thread counts above one")
        status = ("excluded" if reasons else
                  "adapter-required" if adapter_gaps else
                  "adapter-available-unmeasured")
        base.update({
            "status": status,
            "reasons": reasons or adapter_gaps or [
                "the reviewed default-off adapter can execute this public workload, "
                "but this audit contains no measurement"],
            "comparison_scope": (
                "public payload and repaired-output throughput only; parity bytes "
                "and generator matrices differ"),
            "qualifications": ([
                "the adapter uses gf_gen_cauchy1_matrix; ISA-L's gf_gen_rs_matrix "
                "documentation does not guarantee every submatrix is invertible "
                "for many larger K,R pairs",
                "the ISA-L adapter is standalone and default-off; production "
                "Leopard has no ISA-L build or runtime dependency",
            ] + ([
                "ISA-L remains in GF(256) while dyadic parent inflation selects "
                "GF16 for Leopard2; report that field advantage explicitly",
            ] if field == "gf16" and k + r <= 256 else [])),
        })
        return base
    if provider == "jerasure":
        provider_field = _binary_field_for_length(k + r)
        base["provider_field"] = provider_field
        reasons = []
        if provider_field not in ("gf8", "gf16"):
            reasons.append(
                "public K+R exceeds Jerasure's selected GF8/GF16 evaluation-set bound")
        if shard_bytes % JERASURE_ADAPTER_REGION_BYTES:
            reasons.append(
                "shard length is not a multiple of the deterministic 8-byte "
                "Jerasure adapter region contract")
        qualifications = [
            "a reviewed adapter must allocate 8-byte-aligned regions; staging, "
            "padding, and copy bytes for other layouts must be charged explicitly",
        ]
        if provider_field == "gf8" and field == "gf16":
            qualifications.append(
                "Jerasure can remain in GF(256) while dyadic parent inflation "
                "selects GF16 for Leopard2; report that field advantage explicitly")
        base.update({
            "status": "excluded" if reasons else "adapter-required",
            "reasons": reasons or [
                "mathematically comparable public erasure workload, but the "
                "reviewed GF-Complete/Jerasure adapter is not present"],
            "comparison_scope": (
                "public payload and repaired-output throughput only; field/basis "
                "representation, coordinates, generator matrices, and parity bytes differ"),
            "qualifications": qualifications,
        })
        return base
    if provider == "fastecc":
        base.update({
            "status": "excluded",
            "reasons": [
                "FastECC uses prime-field NTT symbols rather than Leopard GF(2^m)",
                "stored parity sectors are wider than input sectors (for example "
                "4100 bytes for a 4096-byte source sector)",
                "the current matrix requires equal source/parity shard bytes",
            ],
            "comparison_scope": (
                "future application-payload study only, with parity expansion and "
                "transfer bytes charged explicitly"),
        })
        return base
    base.update({
        "status": "excluded",
        "reasons": [
            "ECC-Benchmark is a historical harness, not an independent codec",
            "its one/all-loss trials and timing boundaries do not match the signed "
            "Leopard2 loss/reuse/setup matrix",
            "its bundled Leopard adapter targets the legacy API and cannot establish "
            "an independent Leopard2 result",
        ],
        "comparison_scope": "methodology review only",
    })
    return base


def _matrix_module():
    tools_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(tools_dir))
    try:
        import leopard2_benchmark_matrix as matrix
    finally:
        sys.path.pop(0)
    return matrix


def audit_matrix(preset: str) -> dict:
    matrix = _matrix_module()
    specification = matrix.make_spec(
        "/not-executed/bench_leopard2", preset, 1, 3, 1)
    providers = {}
    for provider in sorted(SOURCES):
        status_counts: Counter[str] = Counter()
        reason_counts: Counter[str] = Counter()
        qualification_counts: Counter[str] = Counter()
        for job in specification["jobs"]:
            result = classify(provider, job["benchmark_cell"])
            status_counts[result["status"]] += 1
            for reason in result["reasons"]:
                reason_counts[reason] += 1
            for qualification in result.get("qualifications", []):
                qualification_counts[qualification] += 1
        providers[provider] = {
            "source": dict(SOURCES[provider], retrieved=RETRIEVED),
            "job_status_counts": dict(sorted(status_counts.items())),
            "reason_counts": dict(sorted(reason_counts.items())),
            "qualification_counts": dict(sorted(qualification_counts.items())),
        }
    return {
        "schema": SCHEMA,
        "preset": preset,
        "job_count": len(specification["jobs"]),
        "policy": {
            "wire_compatibility_requires_identical_parity": True,
            "eligible_is_not_measured": True,
            "silent_parameter_substitution": False,
            "setup_and_execution_must_be_separate": True,
            "output_bytes_must_be_equivalent_or_explicitly_charged": True,
        },
        "providers": providers,
    }


def audit_cache(cache: Path) -> dict:
    directories = {
        "isa-l": "isa-l",
        "jerasure": "jerasure",
        "fastecc": "FastECC",
        "ecc-benchmark": "ECC-Benchmark",
    }
    sources = {}
    for provider, directory in directories.items():
        checkout = cache / directory
        record = {"path": str(checkout), "expected_commit": SOURCES[provider]["commit"]}
        try:
            completed = subprocess.run(
                ["git", "-C", str(checkout), "rev-parse", "HEAD"],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=10.0)
        except (OSError, subprocess.SubprocessError) as error:
            record.update(status="unavailable", detail=str(error))
        else:
            actual = completed.stdout.strip()
            record.update(
                status=("verified" if completed.returncode == 0 and
                        actual == SOURCES[provider]["commit"] else "mismatch"),
                actual_commit=actual or None,
                detail=completed.stderr.strip() or None)
        sources[provider] = record
    return {
        "schema": "leopard2-external-source-cache/v1",
        "retrieved": RETRIEVED,
        "sources": sources,
        "host_dependencies": {
            "nasm": shutil.which("nasm"),
            "ignored_cache_nasm": _cached_nasm(cache),
            "pkg_config": shutil.which("pkg-config"),
            "gf_complete_version": _pkg_config_version("gf_complete"),
        },
    }


def _cached_nasm(cache: Path) -> dict | None:
    executable = cache / "toolchains" / "nasm-2.16.03-install" / "bin" / "nasm"
    if not executable.is_file():
        return None
    completed = subprocess.run(
        [str(executable), "-v"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    return {
        "path": str(executable),
        "status": "verified" if completed.returncode == 0 and
                  "NASM version 2.16.03" in completed.stdout else "mismatch",
        "version": completed.stdout.strip() or completed.stderr.strip(),
    }


def _pkg_config_version(package: str) -> str | None:
    executable = shutil.which("pkg-config")
    if executable is None:
        return None
    completed = subprocess.run(
        [executable, "--modversion", package], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False)
    return completed.stdout.strip() if completed.returncode == 0 else None


def self_test() -> None:
    matrix = audit_matrix("required")
    if matrix["job_count"] != 7134:
        raise AuditError("required matrix cardinality changed")
    providers = matrix["providers"]
    if providers["fastecc"]["job_status_counts"] != {"excluded": 7134}:
        raise AuditError("FastECC unequal-storage exclusion changed")
    if providers["ecc-benchmark"]["job_status_counts"] != {"excluded": 7134}:
        raise AuditError("ECC-Benchmark methodology exclusion changed")
    isa_counts = providers["isa-l"]["job_status_counts"]
    if (not isa_counts.get("adapter-available-unmeasured") or
            not isa_counts.get("adapter-required") or not isa_counts.get("excluded")):
        raise AuditError(
            "ISA-L matrix must distinguish adapter-ready, adapter-gap, and "
            "over-256 rows")
    jerasure_counts = providers["jerasure"]["job_status_counts"]
    if jerasure_counts != {"adapter-required": 7134}:
        raise AuditError("even-byte GF8/GF16 Jerasure eligibility changed")
    representative = {
        "K": 240, "R": 16, "requested_profile": "high", "shard_bytes": 4096,
    }
    if classify("isa-l", representative)["status"] != "adapter-available-unmeasured":
        raise AuditError("representative ISA-L cell was incorrectly excluded")
    representative.update(K=512, R=1536, requested_profile="low")
    if classify("isa-l", representative)["status"] != "excluded":
        raise AuditError("GF16 ISA-L cell was incorrectly admitted")
    representative.update(K=129, R=100, requested_profile="high")
    boundary = classify("isa-l", representative)
    if (boundary["status"] != "adapter-available-unmeasured" or
            boundary["leopard2_field"] != "gf16" or
            not boundary["qualifications"]):
        raise AuditError("GF8 field-boundary comparison was incorrectly excluded")
    if not any("gf_gen_cauchy1_matrix" in qualification
               for qualification in boundary["qualifications"]):
        raise AuditError("ISA-L MDS-safe generator requirement is missing")
    threaded = dict(representative, thread_count=8)
    if classify("isa-l", threaded)["status"] != "adapter-required":
        raise AuditError("bounded ISA-L adapter silently admitted multicore work")
    no_loss = dict(representative, loss_count=0)
    if classify("isa-l", no_loss)["status"] != "adapter-available-unmeasured":
        raise AuditError("bounded ISA-L adapter incorrectly rejects no-loss decode")
    jerasure_boundary = classify("jerasure", representative)
    if (jerasure_boundary["status"] != "adapter-required" or
            jerasure_boundary["provider_field"] != "gf8" or
            jerasure_boundary["leopard2_field"] != "gf16" or
            not jerasure_boundary["qualifications"]):
        raise AuditError("Jerasure field-boundary advantage was not disclosed")
    jerasure_gf8_unaligned = classify("jerasure", {
        "K": 8, "R": 8, "requested_profile": "high", "shard_bytes": 1,
    })
    jerasure_gf16_unaligned = classify("jerasure", {
        "K": 129, "R": 100, "requested_profile": "high", "shard_bytes": 2,
    })
    jerasure_aligned = classify("jerasure", {
        "K": 129, "R": 100, "requested_profile": "high", "shard_bytes": 8,
    })
    if (jerasure_gf8_unaligned["status"] != "excluded" or
            jerasure_gf16_unaligned["status"] != "excluded" or
            jerasure_aligned["status"] != "adapter-required" or
            not any("8-byte-aligned" in qualification
                    for qualification in jerasure_aligned["qualifications"])):
        raise AuditError("Jerasure deterministic region contract is not enforced")
    print(json.dumps({
        "isa_l": isa_counts,
        "jerasure": jerasure_counts,
        "required_jobs": matrix["job_count"],
        "status": "PASS",
    }, sort_keys=True))


def audit_isal_checkpoint(path: Path, correctness_path: Path) -> dict:
    tools_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(tools_dir))
    try:
        import leopard2_isal_compare as comparison
    finally:
        sys.path.pop(0)
    try:
        document = json.loads(path.read_text())
        correctness = json.loads(correctness_path.read_text())
        comparison.validate_checkpoint(document, correctness=correctness)
    except (OSError, json.JSONDecodeError, comparison.ComparisonError) as error:
        raise AuditError(f"invalid ISA-L checkpoint: {error}") from error
    return {
        "schema": "leopard2-external-isal-checkpoint-audit/v1",
        "checkpoint": str(path),
        "artifact_sha256": document["artifact_sha256"],
        "cells": len(document["cells"]),
        "results": len(document["results"]),
        "wire_compatible": document["method"]["wire_compatible"],
        "status": "verified-bounded-checkpoint",
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    matrix_parser = subparsers.add_parser("matrix")
    matrix_parser.add_argument(
        "--preset", choices=("smoke", "checkpoint", "balanced-crossover", "required"),
        default="required")
    cache_parser = subparsers.add_parser("cache")
    cache_parser.add_argument(
        "--path", default=".research/leopard2", type=Path)
    checkpoint_parser = subparsers.add_parser("isa-l-checkpoint")
    checkpoint_parser.add_argument("path", type=Path)
    checkpoint_parser.add_argument(
        "--correctness-artifact", type=Path, required=True)
    subparsers.add_parser("self-test")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "matrix":
            output = audit_matrix(arguments.preset)
        elif arguments.command == "cache":
            output = audit_cache(arguments.path)
        elif arguments.command == "isa-l-checkpoint":
            output = audit_isal_checkpoint(
                arguments.path, arguments.correctness_artifact)
        else:
            self_test()
            return 0
    except AuditError as error:
        parser.error(str(error))
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
