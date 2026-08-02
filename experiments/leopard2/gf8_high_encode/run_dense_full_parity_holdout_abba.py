#!/usr/bin/env python3
"""Qualify predeclared dense full-parity GF8 intervals on one quiet core.

This is the isolated, promotion-grade follow-up to the 480-cell saturated
screen.  It measures exactly the 24 holdouts recorded by that screen: 17
boundary/interior targets and seven adjacent controls.  Each cell receives
nine independent, symmetric ABBA rounds.  Candidate and control are the
same-source mode-one/mode-two binaries and must have identical executable
sections plus normalized full-file provenance.  Exact Leopard main is
observed only for the four legal K=R cells and is deliberately informational:
the interval threshold is five percent over the same-source current path.

Passing this runner qualifies arithmetic and validation cost only.  The
candidate source already contains a default-inert selector, so a later
integration campaign must still compare the calibrated production build with
its exact parent to exclude layout-only regressions.
"""

from __future__ import annotations

import importlib.util
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Sequence


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "dense_full_parity_holdout_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-dense-full-parity-holdout-abba/v1"
BASE.SUMMARY_SCHEMA = \
    "leopard2-gf8-dense-full-parity-holdout-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L34g_dense_full_parity_one_block_modeE"
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.TARGET_ROLES = ("selected_control", "selected_main")
BASE.MAIN_TARGET_ROLES = ("selected_main", "adjacent_main")
BASE.TARGET_CONTROL_FLOOR = 1.05
# Exact main is an important reported baseline at K=R, but it is not the
# same-source interval-promotion denominator.  In particular, K=T/B64
# adjacent cells help explain why a count is excluded without converting
# unsupported K<R cells into a false Leopard-main comparison.
BASE.TARGET_MAIN_FLOOR = 0.0
BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = True
BASE.REQUIRE_STRICT_ISOLATION = True
BASE.REQUIRE_STABLE_LEASE_ANCHOR = True
BALANCED_MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L28g_balanced_b64_terminal_modeE"
BASE.AUXILIARY_MODE_EXPECTATIONS = {
    BALANCED_MODE_SYMBOL: {"candidate": 1, "control": 1},
}
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve().with_name("git_capture.py"),
    Path(BASE.MAIN_SUPPORT.__file__).resolve().parents[1] /
        "decoder_dispatch" / "balanced_evidence_common.py",
    Path(BASE.MAIN_SUPPORT.__file__).resolve().parents[3] /
        "tools" / "leopard2_build_provenance.py",
)

# Repeat the three balanced label rotations three times.  Every individual
# round is symmetric around its midpoint; the nine round-level contrasts are
# the independent observations used for the confidence interval.
BASE.TARGET_ORDER = BASE.TARGET_ORDER * 3
BASE.NEIGHBOR_ORDER = BASE.NEIGHBOR_ORDER * 3


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    critical_values = {
        3: BASE.T95_DF2,
        4: BASE.T95_DF3,
        9: 2.3060041350333704,
    }
    BASE.require(len(values) in critical_values,
                 "three, four, or nine independent contrasts are required")
    center = statistics.mean(values)
    half_width = critical_values[len(values)] * statistics.stdev(values) / \
        math.sqrt(len(values))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(values),
    }


BASE.confidence_interval = confidence_interval


# This table is copied from the committed saturated-screen declaration.  Its
# order is stable and intentionally not synthesized from measured outcomes.
HOLDOUTS = (
    (1, 16, 64, True, ("lower_boundary",)),
    (8, 16, 64, True, ("interior",)),
    (15, 16, 64, True, ("upper_boundary",)),
    (16, 16, 64, False, ("adjacent_above",)),
    (1, 16, 256, True, ("lower_boundary",)),
    (8, 16, 256, True, ("interior",)),
    (16, 16, 256, True, ("upper_boundary",)),
    (1, 32, 64, True, ("lower_boundary",)),
    (16, 32, 64, True, ("interior",)),
    (31, 32, 64, True, ("upper_boundary",)),
    (32, 32, 64, False, ("adjacent_above",)),
    (1, 32, 256, True,
     ("lower_boundary", "interior", "upper_boundary")),
    (2, 32, 256, False, ("adjacent_above",)),
    (1, 64, 64, True, ("lower_boundary",)),
    (32, 64, 64, True, ("interior",)),
    (63, 64, 64, True, ("upper_boundary",)),
    (64, 64, 64, False, ("adjacent_above",)),
    (39, 64, 256, False, ("adjacent_below",)),
    (40, 64, 256, True,
     ("lower_boundary", "interior", "upper_boundary")),
    (41, 64, 256, False, ("adjacent_above",)),
    (1, 128, 64, True, ("lower_boundary",)),
    (19, 128, 64, True, ("interior",)),
    (38, 128, 64, True, ("upper_boundary",)),
    (39, 128, 64, False, ("adjacent_above",)),
)


def cells() -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for index, (original_count, recovery_count, shard_bytes,
                selected, holdout_roles) in enumerate(HOLDOUTS):
        uses_main = original_count >= recovery_count
        role = ("selected_main" if uses_main else "selected_control") \
            if selected else \
            ("adjacent_main" if uses_main else "adjacent_control")
        role_text = "selected" if selected else "adjacent"
        output.append({
            "id": (f"{role_text}-k{original_count}-r{recovery_count}-"
                   f"b{shard_bytes}-q1"),
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "holdout_roles": list(holdout_roles),
            "expected_selected": selected,
            "seed": (0x44504600 + recovery_count * 0x10000 +
                     shard_bytes * 0x100 + original_count),
        })
    return output


BASE.cells = cells


def synthetic_rounds(order: Sequence[str]) -> list[dict[str, Any]]:
    digests = {
        "algorithm": "fnv1a64",
        "original_data": "0000000000000001",
        "transmitted_parity": "0000000000000002",
        "recovered_originals": "0000000000000003",
    }
    timings = {"candidate": 1.0, "control": 1.1, "main": 1.2}
    return [{
        "round": round_index,
        "order": list(round_order),
        "invocations": [{
            "implementation": label,
            "normalized": {
                "digests": dict(digests),
                "encode_us": timings[label],
            },
        } for label in round_order],
        "isolation": {"accepted": True},
    } for round_index, round_order in enumerate(order)]


def self_test() -> int:
    selected = [cell for cell in cells() if BASE.cell_is_target(cell)]
    adjacent = [cell for cell in cells() if not BASE.cell_is_target(cell)]
    main_cells = [cell for cell in cells() if BASE.cell_uses_main(cell)]
    BASE.require(len(cells()) == 24 and len(selected) == 17 and
                 len(adjacent) == 7 and len(main_cells) == 4,
                 "predeclared 24-cell holdout classification changed")
    BASE.require(len(BASE.TARGET_ORDER) == 9 and
                 len(BASE.NEIGHBOR_ORDER) == 9,
                 "holdout campaign is not nine-round")
    BASE.require(len({cell["seed"] for cell in cells()}) == 24,
                 "holdout seeds are not unique")
    for cell in cells():
        orders = BASE.TARGET_ORDER if BASE.cell_uses_main(cell) \
            else BASE.NEIGHBOR_ORDER
        analysis = BASE.analyze(cell, synthetic_rounds(orders))
        BASE.require(
            ("main_over_candidate" in analysis) ==
                BASE.cell_uses_main(cell),
            "exact-main classification changed")
        BASE.require(len(analysis["control_over_candidate"]
                         ["round_log_ratios"]) == 9,
                     "analysis did not retain nine contrasts")
    BASE.require(BASE.TARGET_CONTROL_FLOOR == 1.05 and
                 BASE.NEIGHBOR_FLOOR == 1.0 / 1.02 and
                 BASE.TARGET_MAIN_FLOOR == 0.0,
                 "predeclared promotion thresholds changed")
    print("dense full-parity holdout runner self-test passed")
    return 0


def main() -> int:
    if sys.argv[1:] == ["--self-test"]:
        return self_test()
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
