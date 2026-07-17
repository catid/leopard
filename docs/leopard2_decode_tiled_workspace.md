# Leopard2 side-sized decode workspace

This document describes the bounded production candidate tracked by
`leopard-79h.26.2.2`.  It changes execution storage and traversal only.  The
field, coordinate maps, locator values, normalization factors, versioned wire
profiles, and requested decoded bytes are unchanged.

## Scratch geometry

Let `B` be the physical shard byte count, `A=floor(B/64)*64`, and let `L` be
the number of missing originals in an immutable decode plan.  Define
`W=B` for an aligned shard and `W=max(A,64)` for a ragged shard.  Ignoring the
much smaller range and pointer arrays, data storage is:

| Execution path | Aligned bytes | Ragged bytes |
| --- | ---: | ---: |
| Low Algorithm 4 | `B*min(N,2P)` | `W*min(N,2P)+64*(K+R)` |
| High Algorithm 5 | `B*min(N,2T+L)` | `W*min(N,2T+L)+64*(K+R)` |
| Generic fallback | `B*N` | `W*N+64*(K+R)` |

The pattern-independent one-shot high query substitutes `R` for `L`.  Since no
valid erasure recovery can request more than `R` missing originals, the query
is conservative for every later pattern.  Low requested originals already live
in the final `P`-slot accumulator and require no separate retention slots.  If
the tiled count is greater than or equal to `N`, execution retains the regular
materialized specialized kernel, so balanced codes never pay more scratch just
to select the tiled traversal.  A ragged input stages only its final public
coordinate tiles.  The aligned prefix executes directly from caller inputs;
the same work slots are then reused for the one 64-byte tail execution.

The public contract test measures the slope between 64- and 128-byte aligned
shards.  Range and pointer metadata are constant across those two sizes, so the
scratch delta divided by 64 is exactly the full-shard slot count.  It checks
GF8 and GF16 low/high cases plus an `N`-slot forced-generic control.  A second
slope test compares 65 with 129 bytes (66 with 130 for GF16): the fixed
`64*(K+R)` staging term cancels, leaving exactly 64 bytes per work slot.

## Low-rate equivalence

The retained planned decoder first materialized every one of the `N/P` parent
blocks, inverse-transformed every nonempty block, then reduced blocks 1 onward
into block 0.  The tiled traversal preserves the same order:

1. Load and inverse-transform block 0 into the `P`-slot accumulator.
2. Apply the same active-parent formal-derivative term to that accumulator.
3. For each later nonempty block in increasing parent order, load it into the
   reusable `P`-slot tile, apply its shifted inverse transform and fixed block
   factor, then XOR it into the accumulator.
4. Apply the same pruned final transform and inverse-locator factors.

No later step reads a reduced parent block, so reusing its tile cannot change
the result.  The old materialized planned kernel remains a differential oracle.

## High-rate equivalence

The high decoder likewise consumes parent blocks only through their contribution
to the `T`-coefficient accumulator.  The tiled traversal loads and inversely
transforms one block at a time, XORs it into the accumulator in the original
block order, and reuses the second `T`-slot tile.  After the syndrome/evaluator
steps, each requested output block is evaluated into that same tile.  Its
requested coordinates are multiplied by the unchanged reveal factors and
copied to `L` retained kernel-layout slots.  Public scatter then handles a
possible compact GF16 tail exactly as before.

The old `N`-materializing planned high kernel remains available internally and
is compared byte-for-byte with the new requested outputs in the decode-plan
schedule test.

## Validation and promotion state

The combined side-sized and fixed-tail checkpoint has passed the full 50-test
Release suite; strict GCC 13 warning-as-error tests; focused Clang 18
ASan/UBSan and TSan tests; and an AArch64/SSE2NEON compile-only check.  It also
directly compares tiled GF8/GF16 kernels with the retained materialized planned
kernels on sparse inputs and completely empty blocks.  Public contract tests
execute 65/129-byte GF8 and 66/130-byte GF16 splits through low, high, and
forced-generic paths and verify both recovered bytes and the fixed-staging
scratch slope.  The decode-plan schedule target cannot be linked with TSan
because its pre-existing test-only global allocation replacements conflict
with the TSan runtime; the public plan, high/low acceptance, and shared-context
tests cover the concurrent production entry points under TSan instead.

This is not yet a performance promotion claim.  Production integration still
requires an independent source review and isolated paired measurements covering
setup, execution, cache effects, and neighboring high/low regimes.  Forced
generic remains the correctness fallback.
