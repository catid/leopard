#!/usr/bin/env python3
"""Deterministic unit tests for the FF8 binary-linear synthesis portfolio."""

from __future__ import print_function

import importlib.util
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "ff8_circuit_generator",
    str(ROOT / "tools" / "generate_ff8_xor_circuits.py"))
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)

# Explicit row-map fixtures produced from three deterministic sequences of
# reversible row additions.  Keep these as data: constructing them at test
# time with either circuit evaluator would make the expected value circular.
REPRESENTATIVE_WIDE_MAPS = (
    ("wide-64-gates", (
        0x00201001, 0x080000A0, 0x20020084, 0x40211402,
        0x00000010, 0x00000022, 0x82004440, 0x08000080,
        0x00000100, 0x820446E0, 0xC60442E0, 0x00800810,
        0x0225302A, 0x00002000, 0x40014400, 0xC604C2F0,
        0x40010400, 0x20000004, 0x002C3421, 0x00281401,
        0x40900608, 0x82001001, 0x40400008, 0x40AC30A9,
        0x49014400, 0x020420A0, 0x86004440, 0x3A2738AE,
        0x322738AE, 0x30000802, 0xB06738A6, 0x82000000,
    )),
    ("wide-127-gates", (
        0x31DA0F61, 0x48F34DD7, 0x838124D9, 0xDBD242FB,
        0xC2018874, 0x18802861, 0x05C104CB, 0xC080ACBD,
        0x91202745, 0x587B5F13, 0xD2818A74, 0xC2012855,
        0x02013105, 0x40F0807D, 0xB852CDFD, 0xC2018034,
        0xC2002855, 0x72112219, 0xFD57269E, 0x00080868,
        0x75030146, 0xC0A82E6D, 0x0255CDF2, 0xC281245D,
        0x1A810000, 0x639720DD, 0x7306231F, 0x3A17205D,
        0x1080AE29, 0x20020004, 0xB012000C, 0xC5C104C7,
    )),
    ("wide-191-gates", (
        0x02000001, 0x4AB32173, 0x0BF66073, 0x1869381C,
        0xC45188D4, 0x51CFFDB6, 0xD3E97C65, 0x555B2C11,
        0x68B41E1C, 0x87566610, 0xE9B51650, 0x03A65122,
        0x191700F6, 0x445F2405, 0x441140DE, 0xD10EAD01,
        0x7DA55057, 0x75A45055, 0x0314AC5C, 0x45241584,
        0x18BD1C44, 0x9D4E61FC, 0x0D110818, 0xD27299D6,
        0xD8F85544, 0x4CC56831, 0xC6EC4B37, 0x08000002,
        0x0B1447DC, 0xCA31F9D6, 0x863171D6, 0x80043E6A,
    )),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def local_apply_circuit(state, gates, width):
    """Apply CNOTs without using the generator's circuit evaluator."""
    require(0 <= state < (1 << width),
            "circuit input does not fit in %d wires" % width)
    result = state
    for gate_index, (destination, source) in enumerate(gates):
        require(0 <= destination < width and 0 <= source < width,
                "gate %d is outside the %d-wire circuit" %
                (gate_index, width))
        require(destination != source,
                "gate %d uses the same source and destination" % gate_index)
        if (result >> source) & 1:
            result ^= 1 << destination
    return result


def local_circuit_matrix(width, gates):
    """Construct circuit rows from basis states using the local evaluator."""
    rows = [0] * width
    for input_wire in range(width):
        output = local_apply_circuit(1 << input_wire, gates, width)
        for output_wire in range(width):
            if (output >> output_wire) & 1:
                rows[output_wire] |= 1 << input_wire
    return tuple(rows)


def local_compiled_xor_form_cost(width, gates):
    """Independently count distinct non-input linear forms in a circuit."""
    wires = [1 << wire for wire in range(width)]
    forms = set(wires)
    for destination, source in gates:
        wires[destination] ^= wires[source]
        forms.add(wires[destination])
    return len(forms.difference(1 << wire for wire in range(width)))


def local_apply_matrix(rows, state):
    """Apply binary matrix rows with a test-local parity calculation."""
    width = len(rows)
    require(0 <= state < (1 << width),
            "matrix input does not fit in %d wires" % width)
    output = 0
    for output_wire, row in enumerate(rows):
        require(0 <= row < (1 << width),
                "matrix row %d does not fit in %d wires" %
                (output_wire, width))
        parity = bin(row & state).count("1") & 1
        output |= parity << output_wire
    return output


def require_circuit_matches_map(label, rows, gates, states):
    """Independently compare a circuit with explicit expected matrix rows."""
    width = len(rows)
    require(local_circuit_matrix(width, gates) == rows,
            "%s circuit basis matrix mismatch" % label)
    for state_index, state in enumerate(states):
        expected = local_apply_matrix(rows, state)
        observed = local_apply_circuit(state, gates, width)
        require(observed == expected,
                "%s state %d mismatch: input=%#x expected=%#x got=%#x" %
                (label, state_index, state, expected, observed))


def deterministic_map(width, seed, gate_count):
    generator = random.Random(seed)
    gates = []
    for unused in range(gate_count):
        destination = generator.randrange(width)
        source = generator.randrange(width - 1)
        if source >= destination:
            source += 1
        gates.append((destination, source))
    return local_circuit_matrix(width, gates)


def test_matrix_algebra():
    for width, seed, gate_count in ((8, 81, 31), (16, 161, 73),
                                    (32, 321, 127)):
        rows = deterministic_map(width, seed, gate_count)
        inverse = GENERATOR.inverse_matrix(rows)
        require(GENERATOR.transpose_matrix(
                    GENERATOR.transpose_matrix(rows)) == rows,
                "%dx%d transpose was not an involution" % (width, width))
        values = (0, (1 << width) - 1, 1,
                  random.Random(seed ^ 0x55AA).randrange(1 << width))
        for value in values:
            transformed = local_apply_matrix(rows, value)
            require(local_apply_matrix(inverse, transformed) == value,
                    "%dx%d inverse did not recover input" % (width, width))


def test_commuting_cancellation_and_exact_templates():
    commuting = ((0, 1), (2, 3), (0, 1))
    reduced = GENERATOR.cancel_commuting_duplicate_gates(commuting)
    require(reduced == ((2, 3),),
            "commuting duplicate CNOTs did not cancel")
    require(local_circuit_matrix(4, reduced) ==
            local_circuit_matrix(4, commuting),
            "commuting cancellation changed the map")

    # This conjugation cannot commute-cancel, but exact four-wire search can
    # replace the three gates with the equivalent two-gate template.
    template = ((0, 1), (1, 2), (0, 1))
    require(GENERATOR.cancel_commuting_duplicate_gates(template) == template,
            "noncommuting CNOTs were cancelled")
    optimized = GENERATOR.optimize_exact_windows(template, 4)
    require(len(optimized) == 2,
            "bounded exact search missed a three-to-two CNOT template")
    require(local_circuit_matrix(4, optimized) ==
            local_circuit_matrix(4, template),
            "bounded exact rewrite changed the map")

    expected_group_sizes = (1, 6, 168, 20160)
    observed = tuple(
        len(GENERATOR.exact_cnot_table(width)) for width in range(1, 5))
    require(observed == expected_group_sizes,
            "exact GL(n,2) table sizes changed: %r" % (observed,))


def test_bounded_exact_eight_wire_search():
    # Rebuild the complete ball so this test proves deterministic traversal,
    # not merely repeatability of a memoized answer from another fixture.
    GENERATOR._EXACT_CNOT_BALL_8 = None
    GENERATOR._EXACT_CNOT_SUFFIXES_8.clear()
    GENERATOR._EXACT_CNOT_SEARCH_CACHE_8.clear()
    first_ball = GENERATOR.exact_cnot_ball_8()
    depth_counts = tuple(
        sum(len(path) == depth for path in first_ball.values())
        for depth in range(GENERATOR.EXACT_CNOT_8_HALF_DEPTH + 1))
    require(depth_counts == (1, 56, 1904, 50316),
            "8x8 exact-search BFS ball changed: %r" % (depth_counts,))
    first_items = tuple(sorted(first_ball.items()))
    GENERATOR._EXACT_CNOT_BALL_8 = None
    GENERATOR._EXACT_CNOT_SUFFIXES_8.clear()
    second_items = tuple(sorted(GENERATOR.exact_cnot_ball_8().items()))
    require(first_items == second_items,
            "8x8 exact-search BFS traversal was nondeterministic")

    # This six-gate fixture touches all eight wires.  The first three gates
    # have an equivalent two-gate form, while the remaining destinations are
    # independent, so the exact answer has five gates.  More importantly, five
    # target rows differ from identity.  Since each CNOT changes only one
    # destination row, that supplies a test-local lower-bound proof that no
    # four-gate implementation exists.
    fixture = (
        (0, 1), (1, 2), (0, 1),
        (3, 4), (5, 6), (6, 7),
    )
    fixture_rows = local_circuit_matrix(8, fixture)
    require(set(wire for gate in fixture for wire in gate) == set(range(8)),
            "8x8 exact-search fixture does not touch every wire")
    changed_rows = sum(
        row != (1 << row_index)
        for row_index, row in enumerate(fixture_rows))
    require(changed_rows == 5,
            "8x8 exact-search fixture lower bound changed")
    require(GENERATOR.synthesize_reversible_map_exact_bounded_8(
                fixture_rows, 4) is None,
            "8x8 exact search accepted a circuit below its row lower bound")

    exact = GENERATOR.synthesize_reversible_map_exact_bounded_8(
        fixture_rows, 5)
    require(exact is not None and len(exact) == 5,
            "8x8 exact search missed its shortest five-gate circuit")
    require_circuit_matches_map(
        "exact-8x8", fixture_rows, exact, tuple(range(1 << 8)))
    exact_again = GENERATOR.synthesize_reversible_map_exact_bounded_8(
        fixture_rows, 5)
    require(exact_again == exact,
            "8x8 exact-search result was nondeterministic")
    require(GENERATOR.compiled_xor_form_cost(exact, 8) ==
            local_compiled_xor_form_cost(8, exact),
            "compiled XOR form-cost proxy disagreed with local model")

    optimized = GENERATOR.optimize_exact_8wire_windows(fixture, 8)
    require(len(optimized) == 5,
            "8x8 exact peephole failed to shorten a six-gate window")
    require_circuit_matches_map(
        "exact-8x8-window", fixture_rows, optimized,
        (0, 1, 0x55, 0xAA, 0xFF))
    require(GENERATOR.optimize_exact_8wire_windows(fixture, 8) == optimized,
            "8x8 exact peephole was nondeterministic")

    # Exercise every supported bound with circuits generated independently of
    # the production search.  Once the solver reports an optimum, querying one
    # gate below it must fail, and test-local symbolic execution proves the
    # returned program on every 8-bit state.
    exact_random = random.Random(0x8EAC7)
    for source_gate_count in range(
            GENERATOR.EXACT_CNOT_8_MAX_DEPTH + 1):
        source_gates = []
        for unused in range(source_gate_count):
            destination = exact_random.randrange(8)
            source = exact_random.randrange(7)
            if source >= destination:
                source += 1
            source_gates.append((destination, source))
        rows = local_circuit_matrix(8, source_gates)
        shortest = GENERATOR.synthesize_reversible_map_exact_bounded_8(
            rows, GENERATOR.EXACT_CNOT_8_MAX_DEPTH)
        require(shortest is not None and
                len(shortest) <= source_gate_count,
                "8x8 exact search missed a bounded random circuit")
        require_circuit_matches_map(
            "exact-8x8-random-%d" % source_gate_count,
            rows, shortest, tuple(range(1 << 8)))
        if shortest:
            require(GENERATOR.synthesize_reversible_map_exact_bounded_8(
                        rows, len(shortest) - 1) is None,
                    "8x8 exact search returned a non-minimum circuit")


def test_bidirectional_greedy_synthesis():
    # Every input/output-side move must update the explicitly tracked inverse
    # exactly.  The greedy cost depends on both matrices, so a one-sided update
    # bug could otherwise produce plausible but invalid choices.
    rows = deterministic_map(8, 0xB1D1, 53)
    inverse = GENERATOR.inverse_matrix(rows)
    for side in (GENERATOR.GREEDY_OUTPUT_SIDE,
                 GENERATOR.GREEDY_INPUT_SIDE):
        for destination in range(8):
            for source in range(8):
                if destination == source:
                    continue
                changed, changed_inverse = \
                    GENERATOR.apply_bidirectional_cnot(
                        rows, inverse, (side, destination, source))
                require(GENERATOR.inverse_matrix(changed) == changed_inverse,
                        "bidirectional CNOT inverse invariant failed")

    # Include both native circuit widths.  These fixed maps are deliberately
    # dense enough that greedy sharing beats plain row reduction, and running
    # each twice proves deterministic tie breaking.
    fixtures = (
        ("greedy-8", GENERATOR.multiplication_matrix(42)),
        ("greedy-16", deterministic_map(16, 0x1616, 89)),
    )
    for label, fixture in fixtures:
        width = len(fixture)
        GENERATOR._GREEDY_CANONICAL_CIRCUITS.clear()
        first = GENERATOR.synthesize_reversible_map_greedy(fixture)
        second = GENERATOR.synthesize_reversible_map_greedy(fixture)
        fixture_inverse = GENERATOR.inverse_matrix(fixture)
        inverse_first = GENERATOR.synthesize_reversible_map_greedy(
            fixture_inverse)
        baseline = GENERATOR.synthesize_reversible_map(fixture)
        require(first == second,
                "%s greedy synthesis was nondeterministic" % label)
        require(inverse_first == tuple(reversed(first)),
                "%s inverse did not reuse the canonical circuit" % label)
        require(GENERATOR.circuit_matrix(width, first) == fixture,
                "%s greedy synthesis changed the map" % label)
        require(len(first) < len(baseline),
                "%s greedy fixture no longer improves row reduction" % label)

        # Clear the memoization and request the inverse first.  The result must
        # be independent of call order, not merely stable while cached.
        GENERATOR._GREEDY_CANONICAL_CIRCUITS.clear()
        inverse_again = GENERATOR.synthesize_reversible_map_greedy(
            fixture_inverse)
        first_after_inverse = GENERATOR.synthesize_reversible_map_greedy(
            fixture)
        require((inverse_again, first_after_inverse) ==
                (inverse_first, first),
                "%s canonical synthesis depended on call order" % label)


def test_portfolio_determinism_and_quality():
    fixtures = []
    for coefficient in (0, 1, 2, 17, 42, 127, 254, 255):
        fixtures.append((
            "multiply-%d" % coefficient,
            GENERATOR.multiplication_matrix(coefficient)))
    fixtures.extend((
        ("random-8", deterministic_map(8, 0x8008, 47)),
        ("random-16", deterministic_map(16, 0x1616, 89)),
    ))

    for label, rows in fixtures:
        width = len(rows)
        baseline = GENERATOR.synthesize_reversible_map(rows)
        selected, variant = \
            GENERATOR.synthesize_reversible_map_portfolio(rows)
        selected_again, variant_again = \
            GENERATOR.synthesize_reversible_map_portfolio(rows)
        require((selected, variant) == (selected_again, variant_again),
                "%s portfolio selection was nondeterministic" % label)
        require(local_circuit_matrix(width, selected) == rows,
                "%s selected the wrong map" % label)
        require(GENERATOR.compiled_xor_form_cost(selected, width) ==
                local_compiled_xor_form_cost(width, selected),
                "%s compiled XOR form-cost proxy mismatch" % label)
        require(GENERATOR.circuit_key(selected, width) <=
                GENERATOR.circuit_key(baseline, width),
                "%s portfolio regressed the baseline cost" % label)
        require(bool(variant), "%s did not report its chosen variant" % label)


def test_representative_wide_maps():
    # Consume explicit maps and evaluate synthesized CNOTs entirely through
    # test-local code.  This prevents a shared generator evaluator defect from
    # making both circuit construction and validation agree on the wrong map.
    width = 32
    basis_states = tuple(1 << bit for bit in range(width))
    for fixture_index, (label, rows) in enumerate(REPRESENTATIVE_WIDE_MAPS):
        state_random = random.Random(0x51A7E000 + fixture_index)
        deterministic_states = (
            0,
            (1 << width) - 1,
            0xAAAAAAAA,
            0x55555555,
            0x80000001,
        ) + tuple(state_random.randrange(1 << width) for unused in range(64))
        validation_states = basis_states + deterministic_states

        baseline = GENERATOR.synthesize_reversible_map(rows)
        require_circuit_matches_map(
            label + "/baseline", rows, baseline, validation_states)

        orders = GENERATOR.synthesis_column_orders(rows)
        for order_name, order in (orders[0], orders[-1]):
            for pivot in (GENERATOR.PIVOT_MODES[0],
                          GENERATOR.PIVOT_MODES[-1]):
                for reverse in (False, True):
                    circuit = GENERATOR.synthesize_reversible_map_ordered(
                        rows, order, pivot, reverse)
                    variant = "%s/ordered/%s/%s/%s" % (
                        label, order_name, pivot,
                        "reverse" if reverse else "forward")
                    require_circuit_matches_map(
                        variant, rows, circuit, validation_states)

        # Select a complete portfolio circuit for every representative map,
        # then prove it on all 32 basis vectors and deterministic states.
        selected, variant = \
            GENERATOR.synthesize_reversible_map_portfolio(rows)
        require(bool(variant), "%s portfolio variant was empty" % label)
        require_circuit_matches_map(
            label + "/portfolio/" + variant,
            rows, selected, validation_states)


def main():
    test_matrix_algebra()
    test_commuting_cancellation_and_exact_templates()
    test_bounded_exact_eight_wire_search()
    test_bidirectional_greedy_synthesis()
    test_portfolio_determinism_and_quality()
    test_representative_wide_maps()
    print("FF8 XOR synthesis portfolio tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
