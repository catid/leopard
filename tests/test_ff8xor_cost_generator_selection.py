#!/usr/bin/env python3
"""Verify portable and explicit machine profiles drive generator selection."""

from __future__ import print_function

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))
import ff8_xor_cost_model as cost_model  # noqa: E402

SPEC = importlib.util.spec_from_file_location(
    "ff8xor_circuit_generator_cost_test",
    str(TOOLS / "generate_ff8_xor_circuits.py"))
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def profiles():
    artifact = cost_model.load_profile_artifact(
        ROOT / "generated" / "FF8XorCostProfiles.json")
    return (
        cost_model.find_profile(artifact, "portable-default-v1"),
        cost_model.find_profile(
            artifact, "amd-zen5-gcc13-avx2-avx512-v1"))


def check_selected(label, rows, portable_gates, machine_gates,
                   portable, machine):
    width = len(rows)
    require(GENERATOR.circuit_matrix(width, portable_gates) == rows,
            "%s portable selection changed the map" % label)
    require(GENERATOR.circuit_matrix(width, machine_gates) == rows,
            "%s machine selection changed the map" % label)
    require(GENERATOR.selection_key(machine_gates, width, machine) <=
            GENERATOR.selection_key(portable_gates, width, machine),
            "%s machine profile selected a higher-cost circuit" % label)
    require(GENERATOR.selection_key(portable_gates, width, portable) <=
            GENERATOR.selection_key(machine_gates, width, portable),
            "%s portable profile no longer preserves historical ordering" %
            label)
    require(portable_gates != machine_gates,
            "%s fixture no longer exercises distinct selections" % label)


def main():
    portable, machine = profiles()

    multiplier_rows = GENERATOR.multiplication_matrix(25)
    portable_multiplier, unused = \
        GENERATOR.synthesize_reversible_map_portfolio(
            multiplier_rows, portable)
    machine_multiplier, unused = \
        GENERATOR.synthesize_reversible_map_portfolio(
            multiplier_rows, machine)
    check_selected(
        "multiplier log 25", multiplier_rows,
        portable_multiplier, machine_multiplier, portable, machine)

    multiplication_matrices = tuple(
        GENERATOR.multiplication_matrix(log_multiplier)
        for log_multiplier in range(GENERATOR.FIELD_ORDER))
    for inverse in (False, True):
        function = ((lambda state:
                     GENERATOR.scalar_inverse_butterfly(
                         state, 4, multiplication_matrices))
                    if inverse else
                    (lambda state:
                     GENERATOR.scalar_forward_butterfly(
                         state, 4, multiplication_matrices)))
        rows = GENERATOR.linear_map(
            GENERATOR.WIRE_COUNT_BUTTERFLY, function)
        direct = GENERATOR.direct_butterfly_circuit(
            4, inverse, multiplication_matrices)
        portable_gates, unused = GENERATOR.choose_butterfly_circuit(
            rows, direct, portable)
        machine_gates, unused = GENERATOR.choose_butterfly_circuit(
            rows, direct, machine)
        check_selected(
            ("IFFT" if inverse else "FFT") + " skew 4",
            rows, portable_gates, machine_gates, portable, machine)

    print("FF8 XOR generator cost-profile selection tests passed")


if __name__ == "__main__":
    main()
