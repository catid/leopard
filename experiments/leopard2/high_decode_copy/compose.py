#!/usr/bin/env python3
"""Bind same-source Algorithm 5 A/B evidence to independent exact-main runs."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-high-decode-copy-composite/v1"


def load_module(name: str, path: Path) -> ModuleType:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


ROOT = Path(__file__).resolve().parents[3]
SUPPORT = load_module(
    "leopard2_high_copy_compose_exact_main",
    ROOT / "experiments/leopard2/main_compare/run_abba.py")
A_B = load_module(
    "leopard2_high_copy_compose_ab",
    ROOT / "experiments/leopard2/high_decode_copy/run_abba.py")
EvidenceError = SUPPORT.EvidenceError


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_signed_manifest(path: Path, schema: str, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise EvidenceError(f"cannot read {label}: {error}") from error
    SUPPORT.verify_signature(value, label)
    require(isinstance(value, dict) and value.get("schema") == schema and
            value.get("valid") is True,
            f"{label} is not valid evidence")
    return value


def load_raw(path: Path, manifest: Mapping[str, Any], label: str) -> dict[str, Any]:
    identity = manifest.get("raw")
    require(isinstance(identity, dict), f"{label} has no raw bundle")
    raw_path = SUPPORT.safe_evidence_path(path.parent, identity.get("path"))
    require(raw_path.is_file() and raw_path.stat().st_size == identity.get("size") and
            SUPPORT.sha256_file(raw_path) == identity.get("sha256"),
            f"{label} raw file identity differs")
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    SUPPORT.verify_signature(raw, f"{label} raw bundle")
    require(raw.get("digest") == identity.get("payload_digest"),
            f"{label} raw payload identity differs")
    return raw


def result_documents(records: Sequence[Mapping[str, Any]], cell: str,
                     *, mode: str | None = None,
                     implementation: str | None = None) -> list[Mapping[str, Any]]:
    require((mode is None) != (implementation is None),
            "cross-bundle record selector must name exactly one provider")
    cell_key = "cell" if mode is not None else "cell_id"
    result = []
    for record in records:
        if record.get(cell_key) != cell:
            continue
        if mode is not None and record.get("mode") != mode:
            continue
        if implementation is not None and record.get("implementation") != implementation:
            continue
        document = record.get("result")
        require(isinstance(document, dict), "cross-bundle record has no result document")
        result.append(document)
    return result


def cross_validate(
    ab_raw: Mapping[str, Any],
    exact_raws: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    matrix = ab_raw.get("matrix")
    require(matrix == A_B.load_matrix(ROOT / A_B.MATRIX_RELATIVE),
            "A/B raw bundle does not contain the canonical matrix")
    ab_source = ab_raw.get("identities_initial", {}).get("source", {})
    commit = ab_source.get("head")
    require(isinstance(commit, str) and len(commit) == 40,
            "A/B source commit is absent")
    ab_records = ab_raw.get("records")
    require(isinstance(ab_records, list), "A/B records are absent")
    expected_main = {
        workspace: {cell["id"] for cell in matrix["cells"]
                    if cell["workspace"] == workspace and
                    cell["exact_main_eligible"]}
        for workspace in ("materialized", "tiled")
    }
    output: dict[str, Any] = {}
    for workspace, raw in exact_raws.items():
        require(workspace in expected_main,
                "unexpected exact-main workspace bundle")
        campaign = raw.get("campaign")
        identities = raw.get("identities_initial")
        records = raw.get("invocations")
        require(isinstance(campaign, dict) and isinstance(identities, dict) and
                isinstance(records, list) and
                campaign.get("candidate_mode") == workspace,
                f"exact-main {workspace} bundle has the wrong candidate mode")
        candidate_source = identities.get("candidate_source", {})
        baseline_source = identities.get("baseline_source", {})
        require(candidate_source.get("head") == commit and
                baseline_source.get("head") == SUPPORT.MAIN_COMMIT and
                baseline_source.get("detached") is True,
                f"exact-main {workspace} bundle is not independent production evidence "
                "at the identical Leopard2 source commit")
        cells = campaign.get("cells")
        require(isinstance(cells, list) and
                {cell.get("identifier") for cell in cells
                 if isinstance(cell, dict)} == expected_main[workspace],
                f"exact-main {workspace} cells differ from the eligible matrix subset")
        for cell in matrix["cells"]:
            if cell["workspace"] != workspace or not cell["exact_main_eligible"]:
                continue
            identifier = cell["id"]
            ab_documents = result_documents(ab_records, identifier, mode="no-copy")
            main_candidates = result_documents(
                records, identifier, implementation="candidate")
            main_baselines = result_documents(
                records, identifier, implementation="baseline")
            require(len(ab_documents) == 6 and len(main_candidates) == 6 and
                    len(main_baselines) == 6,
                    f"{identifier} does not have six observations per compared provider")
            digest = ab_documents[0].get("workload_digests")
            missing = ab_documents[0].get("parameters", {}).get(
                "missing_original_indices")
            all_documents = ab_documents + main_candidates + main_baselines
            require(all(document.get("workload_digests") == digest and
                        document.get("parameters", {}).get(
                            "missing_original_indices") == missing
                        for document in all_documents),
                    f"{identifier} exact-main and A/B workload/path digests differ")
            require(all(document.get("correctness", {}).get(
                        "leopard2_round_trip") is True
                        for document in main_candidates) and
                    all(document.get("correctness", {}).get("round_trip") is True
                        for document in main_baselines),
                    f"{identifier} exact-main comparison did not round trip")
            require(all(
                        document.get("resolved", {}).get("backend") ==
                            cell["backend"] and
                        document.get("resolved", {}).get("field") ==
                            cell["field"]
                        for document in main_candidates),
                    f"{identifier} exact-main candidate did not execute the "
                    "same backend and field as the A/B role")
            output[identifier] = {
                "classification": "same-source-and-independent-exact-main",
                "workload_digests": digest,
                "missing_original_indices": missing,
                "exact_main_workspace": workspace,
            }
    require(set(output) == expected_main["materialized"] | expected_main["tiled"],
            "composite omitted an exact-main-eligible cell")
    for cell in matrix["cells"]:
        if cell["exact_main_eligible"]:
            continue
        require(cell["shard_bytes"] % 64 != 0 and cell["id"] not in output,
                "a non-64-byte tail was relabeled as exact-main evidence")
        documents = result_documents(ab_records, cell["id"], mode="no-copy")
        require(len(documents) == 6,
                "same-source-only tail lacks its six no-copy observations")
        output[cell["id"]] = {
            "classification": "same-source-only-tail",
            "reason": "exact Leopard main requires shard_bytes divisible by 64",
            "workload_digests": documents[0]["workload_digests"],
            "missing_original_indices":
                documents[0]["parameters"]["missing_original_indices"],
        }
    require(set(output) == {cell["id"] for cell in matrix["cells"]},
            "composite cell classification is incomplete")
    return output


def compose(options: argparse.Namespace) -> int:
    output = options.output.resolve()
    require(not output.exists(), "composite output already exists")
    # Re-run both verifiers against current source/build closures.  The hook
    # A/B campaign can never stand in for either production exact-main run.
    require(A_B.verify_campaign(SimpleNamespace(
        manifest=options.ab_manifest,
        no_current_input_check=False)) == 0,
        "same-source A/B verification failed")
    main_paths = {
        "materialized": options.exact_main_materialized.resolve(strict=True),
        "tiled": options.exact_main_tiled.resolve(strict=True),
    }
    for path in main_paths.values():
        require(SUPPORT.verify_campaign(SimpleNamespace(
            manifest=path, no_current_input_check=False)) == 0,
            "independent production exact-main verification failed")
    ab_path = options.ab_manifest.resolve(strict=True)
    ab_manifest = load_signed_manifest(
        ab_path, A_B.MANIFEST_SCHEMA, "same-source A/B manifest")
    ab_raw = load_raw(ab_path, ab_manifest, "same-source A/B manifest")
    exact_manifests = {}
    exact_raws = {}
    for workspace, path in main_paths.items():
        manifest = json.loads(path.read_text(encoding="utf-8"))
        SUPPORT.verify_signature(manifest, f"exact-main {workspace} manifest")
        require(manifest.get("valid") is True and
                manifest.get("schema") == SUPPORT.MANIFEST_SCHEMA,
                f"exact-main {workspace} manifest must use current v4 evidence")
        exact_manifests[workspace] = manifest
        exact_raws[workspace] = load_raw(
            path, manifest, f"exact-main {workspace} manifest")
    cells = cross_validate(ab_raw, exact_raws)
    payload = SUPPORT.signed({
        "schema": SCHEMA,
        "created_utc": SUPPORT.utc_now(),
        "valid": True,
        "validity_is_independent_of_speed": True,
        "source_commit": ab_raw["identities_initial"]["source"]["head"],
        "exact_main_commit": SUPPORT.MAIN_COMMIT,
        "same_source_manifest": SUPPORT.artifact_identity(ab_path, "source_file"),
        "exact_main_manifests": {
            workspace: SUPPORT.artifact_identity(path, "source_file")
            for workspace, path in main_paths.items()
        },
        "same_source_analysis": ab_manifest["analysis"],
        "exact_main_analysis": {
            workspace: exact_manifests[workspace]["analysis"]
            for workspace in ("materialized", "tiled")
        },
        "cells": cells,
        "tail_policy": "non-64-byte cells are same-source-only and never acquire "
                       "an inferred exact-main result",
    })
    SUPPORT.write_json_exclusive(output, payload)
    print(output)
    return 0


def synthetic_raws() -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    matrix = A_B.load_matrix(ROOT / A_B.MATRIX_RELATIVE)
    commit = "a" * 40
    ab_records = []
    for cell in matrix["cells"]:
        for mode in ("no-copy", "copy-fallback"):
            for _ in range(6):
                document = {
                    "workload_digests": {"algorithm": "fnv1a64",
                        "original_data": cell["id"], "transmitted_parity": "p",
                        "recovered_originals": "r"},
                    "parameters": {"missing_original_indices": [cell["losses"]]},
                }
                ab_records.append({"cell": cell["id"], "mode": mode,
                                   "result": document})
    ab = {"matrix": matrix, "identities_initial": {"source": {"head": commit}},
          "records": ab_records}
    exact = {}
    for workspace in ("materialized", "tiled"):
        cells = [cell for cell in matrix["cells"]
                 if cell["workspace"] == workspace and cell["exact_main_eligible"]]
        records = []
        for cell in cells:
            digest = next(record["result"]["workload_digests"] for record in ab_records
                          if record["cell"] == cell["id"])
            missing = [cell["losses"]]
            for implementation in ("candidate", "baseline"):
                for _ in range(6):
                    correctness = {"leopard2_round_trip": True} if \
                        implementation == "candidate" else {"round_trip": True}
                    resolved = ({"backend": cell["backend"],
                                 "field": cell["field"],
                                 "selected_decode_path": workspace}
                                if implementation == "candidate" else {})
                    records.append({"cell_id": cell["id"],
                                    "implementation": implementation,
                                    "result": {"workload_digests": digest,
                                        "parameters": {
                                            "missing_original_indices": missing},
                                        "correctness": correctness,
                                        "resolved": resolved}})
        exact[workspace] = {
            "campaign": {"candidate_mode": workspace,
                         "cells": [{"identifier": cell["id"]} for cell in cells]},
            "identities_initial": {
                "candidate_source": {"head": commit},
                "baseline_source": {"head": SUPPORT.MAIN_COMMIT,
                                    "detached": True}},
            "invocations": records,
        }
    return ab, exact


def self_test() -> None:
    ab, exact = synthetic_raws()
    result = cross_validate(ab, exact)
    require(sum(value["classification"] == "same-source-only-tail"
                for value in result.values()) == 4,
            "synthetic composite tail count changed")
    mutations = []
    wrong_commit = copy.deepcopy(exact)
    wrong_commit["materialized"]["identities_initial"]["candidate_source"]["head"] = "b" * 40
    mutations.append(wrong_commit)
    wrong_mode = copy.deepcopy(exact)
    wrong_mode["tiled"]["campaign"]["candidate_mode"] = "materialized"
    mutations.append(wrong_mode)
    wrong_digest = copy.deepcopy(exact)
    wrong_digest["materialized"]["invocations"][0]["result"][
        "workload_digests"]["recovered_originals"] = "changed"
    mutations.append(wrong_digest)
    tail_leak = copy.deepcopy(exact)
    tail_cell = next(cell for cell in ab["matrix"]["cells"]
                     if not cell["exact_main_eligible"] and
                     cell["workspace"] == "materialized")
    tail_leak["materialized"]["campaign"]["cells"].append(
        {"identifier": tail_cell["id"]})
    mutations.append(tail_leak)
    wrong_backend = copy.deepcopy(exact)
    next(record for record in wrong_backend["materialized"]["invocations"]
         if record["implementation"] == "candidate")["result"]["resolved"][
             "backend"] = "scalar"
    mutations.append(wrong_backend)
    for mutation in mutations:
        try:
            cross_validate(ab, mutation)
        except EvidenceError:
            continue
        raise EvidenceError(
            "adversarial commit/mode/digest/tail/backend mutation passed")
    print("high-decode copy composite self-test passed: 5 mutations rejected")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test")
    compose_parser = commands.add_parser("compose")
    compose_parser.add_argument("--ab-manifest", required=True, type=Path)
    compose_parser.add_argument("--exact-main-materialized", required=True, type=Path)
    compose_parser.add_argument("--exact-main-tiled", required=True, type=Path)
    compose_parser.add_argument("--output", required=True, type=Path)
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    return (self_test() or 0) if options.command == "self-test" else compose(options)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EvidenceError, OSError, ValueError) as error:
        print(f"high-decode copy composite error: {error}", file=sys.stderr)
        raise SystemExit(1)
