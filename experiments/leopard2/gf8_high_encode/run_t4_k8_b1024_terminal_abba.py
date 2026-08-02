#!/usr/bin/env python3
"""Qualify the exact GF8 K=8/R=3..4/B=1024 T=4 promotion.

This configured specialization of the hardened T=4 ABBA runner compares a
clean candidate with an executable-section-identical same-source control that
disables only K=8/R=3..4 T=4 classification.  The family-wide T=4 selector is
required to stay enabled in both binaries.  The first 42 cells and seeds are
preserved from the authoritative B=512 campaign, B=1024 supplies the two new
targets, and B=2048 supplies two inert same-path controls.

The candidate and control SHA-256 values pin the immutable binaries qualified
by this runner.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


BENCHMARK_CPU = 14
RESERVED_SIBLING = 30
K8_T4_MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L25g_k8r3r4_t4_terminal_modeE"
FAMILY_T4_MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L28g_high_t4_batch_binding_modeE"
MIN_CREDIBLE_TARGET_SPEEDUP = 1.05
MAX_NEIGHBOR_REGRESSION_FACTOR = 1.02
MIN_NEIGHBOR_SPEEDUP = 1.0 / MAX_NEIGHBOR_REGRESSION_FACTOR
MAX_SAME_PATH_SPEEDUP = MAX_NEIGHBOR_REGRESSION_FACTOR


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t4_packed_terminal_family_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t4_k8_b1024_hardened_abba_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load hardened ABBA base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t4-k8-b1024-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t4-k8-b1024-terminal-summary/v1"
BASE.BENCHMARK_CPU = BENCHMARK_CPU
BASE.RESERVED_SIBLING = RESERVED_SIBLING
BASE.MODE_SYMBOL = K8_T4_MODE_SYMBOL
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.CONTROL_EXTRA_ARGUMENTS = ("--disable-k8r3r4-t4-terminal",)
BASE.CONTROL_BUILD_MARKER = \
    "k8r3r4_t4_terminal_diagnostic_disabled"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v11"
BASE.AUXILIARY_MODE_EXPECTATIONS = {
    FAMILY_T4_MODE_SYMBOL: {
        "candidate": 1,
        "control": 1,
    },
}
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.git_capture.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.git_capture._build_provenance.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.link_common.__file__).resolve(),
)
BASE.EXPECTED_BINARY_SHA256 = {
    "candidate": "0" * 64,
    "control": "0" * 64,
    "main":
        "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910",
}
BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = False
BASE.REQUIRE_EQUAL_EXECUTABLE_PATH_LENGTHS = True


def campaign_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []

    def append(
        original_count: int,
        recovery_count: int,
        shard_bytes: int,
        role: str,
        selected: bool,
    ) -> None:
        index = len(cells)
        cells.append({
            "id": (
                f"{role}-k{original_count}-r{recovery_count}"
                f"-b{shard_bytes}-q1"
            ),
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": BASE.CELL_REUSE,
            "role": role,
            "candidate_selected": selected,
            # The diagnostic control disables only K=8 classification.  Every
            # K=4..7 regression cell therefore remains on the same terminal
            # in candidate and control.
            "control_selected": original_count < 8,
            "seed": 0x74370000 + index,
        })

    # Preserve the first 42 cells, including order, role, and seed, from the
    # authoritative B=512 campaign.
    for original_count in range(4, 8):
        for recovery_count in (3, 4):
            for shard_bytes in (64, 128, 256):
                append(original_count, recovery_count, shard_bytes,
                       "preexisting_regression", True)
    for recovery_count in (3, 4):
        append(4, recovery_count, 1024,
               "preexisting_regression", True)

    for original_count in range(4, 8):
        for recovery_count in (3, 4):
            append(original_count, recovery_count, 512,
                   "b512_regression", True)

    for recovery_count in (3, 4):
        for shard_bytes in (64, 128, 256):
            append(8, recovery_count, shard_bytes,
                   "k8_regression", True)
        append(8, recovery_count, 512, "k8_b512_target", True)

    # Re-role the final two cells from the B=512 campaign without changing
    # their dimensions, order, or deterministic seeds.
    for recovery_count in (3, 4):
        append(8, recovery_count, 1024, "k8_b1024_target", True)

    # Exact-size selection must not spill into the next dyadic byte count.
    for recovery_count in (3, 4):
        append(8, recovery_count, 2048, "k8_b2048_control", False)

    expected_roles = {
        "preexisting_regression": 26,
        "b512_regression": 8,
        "k8_regression": 6,
        "k8_b512_target": 2,
        "k8_b1024_target": 2,
        "k8_b2048_control": 2,
    }
    role_counts = {
        role: sum(cell["role"] == role for cell in cells)
        for role in expected_roles
    }
    BASE.require(
        len(cells) == 46 and
        len({cell["id"] for cell in cells}) == len(cells) and
        len({cell["seed"] for cell in cells}) == len(cells) and
        role_counts == expected_roles and
        sum(cell["candidate_selected"] for cell in cells) == 44 and
        sum(cell["control_selected"] for cell in cells) == 34 and
        all(cell["control_selected"] == (cell["K"] < 8)
            for cell in cells),
        "K8/B=1024 T=4 acceptance matrix is incomplete")
    return cells


def evaluate_promotion(
    analyses: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any]:
    """Apply the predeclared timing gates without rejecting negative data."""
    expected_roles = {
        "preexisting_regression": 26,
        "b512_regression": 8,
        "k8_regression": 6,
        "k8_b512_target": 2,
        "k8_b1024_target": 2,
        "k8_b2048_control": 2,
    }
    role_counts = {
        role: sum(item["cell"]["role"] == role for item in analyses)
        for role in expected_roles
    }
    BASE.require(len(analyses) == 46 and role_counts == expected_roles,
                 "promotion evaluator received the wrong matrix")

    failures: list[dict[str, Any]] = []

    def fail(
        item: Mapping[str, Any],
        contrast: str,
        statistic: str,
        actual: float,
        relation: str,
        threshold: float,
    ) -> None:
        failures.append({
            "cell": item["cell"]["id"],
            "role": item["cell"]["role"],
            "contrast": contrast,
            "statistic": statistic,
            "actual": actual,
            "required_relation": relation,
            "threshold": threshold,
        })

    for item in analyses:
        role = item["cell"]["role"]
        control = item["candidate_vs_control"]
        main = item["candidate_vs_main"]

        # The two new targets and all eight previously promoted K=8 cells
        # require a credible five-percent win over both exact main and the
        # K8-disabled same-source fallback.
        if role in {
                "k8_regression", "k8_b512_target", "k8_b1024_target"}:
            for contrast_name, contrast in (
                    ("candidate_vs_control", control),
                    ("candidate_vs_main", main)):
                lower = float(contrast["ci95"][0])
                if lower < MIN_CREDIBLE_TARGET_SPEEDUP:
                    fail(item, contrast_name, "ci95_lower", lower,
                         ">=", MIN_CREDIBLE_TARGET_SPEEDUP)
            continue

        # K=4..7 use the same production terminal in both binaries.  Preserve
        # a symmetric two-percent same-source band and permit at most a
        # two-percent regression relative to exact main.
        if role in {"preexisting_regression", "b512_regression"}:
            same_source = float(control["speedup"])
            if same_source < MIN_NEIGHBOR_SPEEDUP:
                fail(item, "candidate_vs_control", "speedup", same_source,
                     ">=", MIN_NEIGHBOR_SPEEDUP)
            if same_source > MAX_SAME_PATH_SPEEDUP:
                fail(item, "candidate_vs_control", "speedup", same_source,
                     "<=", MAX_SAME_PATH_SPEEDUP)
            exact_main = float(main["speedup"])
            if exact_main < MIN_NEIGHBOR_SPEEDUP:
                fail(item, "candidate_vs_main", "speedup", exact_main,
                     ">=", MIN_NEIGHBOR_SPEEDUP)
            continue

        # B=2048 is outside both K8 selectors.  Its exact-main result remains
        # descriptive, while candidate/control must stay within two percent.
        BASE.require(role == "k8_b2048_control",
                     "promotion evaluator encountered an unknown role")
        same_source = float(control["speedup"])
        if same_source < MIN_NEIGHBOR_SPEEDUP:
            fail(item, "candidate_vs_control", "speedup", same_source,
                 ">=", MIN_NEIGHBOR_SPEEDUP)
        if same_source > MAX_SAME_PATH_SPEEDUP:
            fail(item, "candidate_vs_control", "speedup", same_source,
                 "<=", MAX_SAME_PATH_SPEEDUP)

    return {
        "schema":
            "leopard2-gf8-t4-k8-b1024-promotion-evaluation/v1",
        "scope": "timing_only",
        "passed": not failures,
        "thresholds": {
            "minimum_target_ci95_lower": MIN_CREDIBLE_TARGET_SPEEDUP,
            "maximum_neighbor_regression_factor":
                MAX_NEIGHBOR_REGRESSION_FACTOR,
            "minimum_neighbor_speedup": MIN_NEIGHBOR_SPEEDUP,
            "maximum_same_path_speedup": MAX_SAME_PATH_SPEEDUP,
        },
        "role_counts": role_counts,
        "failure_count": len(failures),
        "failures": failures,
        "external_static_gate": {
            "required": True,
            "maximum_executable_growth_bytes": 256,
            "new_backend_body_allowed": False,
        },
    }


BASE.campaign_cells = campaign_cells
BASE.PROMOTION_EVALUATOR = evaluate_promotion


if __name__ == "__main__":
    raise SystemExit(BASE.main())
