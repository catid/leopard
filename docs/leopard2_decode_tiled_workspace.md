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

## Calibrated execution policy

The side-sized traversal is the production default.  One narrow deterministic
exception retains the regular `N`-slot high kernel for a single stripe when all
of the following are true:

- the profile is legacy-high GF8 with `K=224`, `R=T=32`, and `N=256`;
- one through eight originals are missing;
- rounded shard bytes are 24 through 64 KiB on AVX2, or 32 through 64 KiB on
  SSSE3; and
- the call contains exactly one batch item on AVX2.  SSSE3 retains the
  materialized exception for every batch size pending a dedicated batch sweep.

Scalar, two-or-more-item AVX2 batches, more than eight missing originals, the
neighboring `T=16` and `T=64` parents, and sizes outside those byte intervals
remain tiled.  An automatic plan reserves enough scratch for the materialized
exception; an AVX2 multi-item batch reuses that caller allocation while
executing the smaller tiled traversal.  This conservative query keeps immutable
plans usable for both single and batch calls without allocating in execution.

`LEO2_CODEC_FORCE_TILED_DECODE` and
`LEO2_CODEC_FORCE_MATERIALIZED_DECODE` are diagnostic flags for differential
tests and offline calibration.  They select only a workspace/traversal for the
same specialized decoder, never a field, profile, coordinate map, or wire
format.  They are mutually exclusive; forcing generic together with either is
invalid.  Forced tiled reserves the full `2P` or `2T+L` geometry even when that
diagnostic workspace is larger than `N`; the normal dispatcher retains the
materialized traversal in that case.

## Isolated performance evidence

The retained machine-readable summary is historical v1 evidence:
`experiments/leopard2/decoder_dispatch/results/tiled_high_amd_9950x3d.json`.
Its file SHA-256 is
`812332d2edee285b59adb0a45751b111df91bc6fbd1a2294f5ad94b4037d3283` and its
canonical content SHA-256 is
`7654e3d95f83edf854b9ae0f4f2cb2f48ef849e40d0205e396b614ad686d9cd0`.
The identity-recorded raw manifest is intentionally kept in the ignored research
cache; its SHA-256 is
`011b9c2815d320a558c4b802eacd4f29819dfa5d3c565f2978768a0744d87bf3`.

That v1 run compared clean detached control commit
`3c75ec7131b5eb36bb4e07cb8e40cc6dc1620703` with tiled candidate commit
`bc7162a16b93fae4e66179cb477e3532176d8405`.  It contains 117 cells and 1,404
benchmark invocations: three independent paired-log ABBA/BAAB rounds, 11
internal samples, a dedicated physical core (CPU 14), and an idle SMT sibling
(CPU 30).  The matrix covers AVX2 and SSSE3 `T=32` byte/loss boundaries,
scalar, batch size two, plan reuse 1 and 64, and AVX2 `T=16`/`T=64` neighbors.
Setup and execution are reported separately.  Candidate scratch fell by
46.82% to 84.34% across the matrix.

The 14 cells with a 95% confidence interval entirely worse than a 2%
regression were confined to the single-stripe `T=32` AVX2/SSSE3 region encoded
above.  The matrix also found 46 credible gains of at least 5%.  Scalar `T=32`
gained 6.90% to 11.74%; AVX2 `T=16` gained 6.05% to 11.49%; and AVX2 `T=64`
gained 1.92% to 18.27%.  At batch size two, AVX2 `T=32` tied or won throughout
and gained roughly 20% to 50% in the 64-to-80-KiB cells.  These results support
a region dispatcher rather than rolling back the bounded workspace globally.

Replay the retained v1 manifest with the version-aware analyzer using:

    python3 experiments/leopard2/decoder_dispatch/analyze_tiled_high.py \
        --manifest HISTORICAL_V1_MANIFEST.json \
        --output /tmp/tiled-high-summary.json --require-rounds 3

New authoritative collection uses the v2 same-binary contract.  Build the
benchmark from the clean commit to be measured, then run:

    python3 experiments/leopard2/decoder_dispatch/run_tiled_high_abba.py \
        --source-root CLEAN_WORKTREE \
        --source-commit FULL_40_HEX_HEAD \
        --binary build-audit/bench_leopard2 \
        --output OUTPUT_DIRECTORY --cpu 14 --sibling 30
    python3 experiments/leopard2/decoder_dispatch/analyze_tiled_high.py \
        --manifest OUTPUT_DIRECTORY/manifest.json \
        --output /tmp/tiled-high-v2-summary.json --require-rounds 3

V2 invokes the exact same executable for both roles.  Control adds only
`--force-materialized`; candidate adds only `--force-tiled`.  Benchmark JSON
must report the corresponding boolean pair (`false,true` and `true,false`,
respectively).  Both keys must be present and boolean in v2.  Historical/mixed
v1 replay accepts either both absent (the retained campaign) or both present as
boolean false; its retained commands remain selector-free.  Partial presence,
an active v1 selector, both active, or a v2 role/selector mismatch is rejected.

The v2 manifest binds the clean source HEAD and tree, benchmark and CMake-cache
bytes, runner and pair-lease source bytes, exact argv including executable,
mode, iteration count, seed, and JSON destination, retained raw/stdout/stderr
hashes, workload digests, lease identity, and per-invocation plus campaign CPU
and sibling counter deltas.  Analysis recomputes median, MAD, minimum, maximum,
execution rates, and amortized decode rates from retained samples.  The current
schema fails closed on unknown fields or versions.  `--source-root` on the
analyzer permits relocation to another clean checkout only when the Git tree,
binary, cache, runner, and lease-source bytes still match the manifest, and the
analyzer itself is executed from its canonical path in that clean tree.

The logical CPU pair must be replaced by an allowed physical-core/sibling pair
on another host.  The runner moves itself to housekeeping CPUs while each child
is pinned to the measured CPU and records both SMT siblings.  Absolute timings
are host-specific; source/binary identity, round-trip digests, scratch geometry,
raw samples, and evidence validation are the reproducible contract.

## Correctness validation

The combined side-sized and fixed-tail checkpoint passed the full 50-test
Release suite; strict GCC 13 warning-as-error tests; focused Clang 18
ASan/UBSan and TSan tests; and an AArch64/SSE2NEON compile-only check.  It also
directly compared tiled GF8/GF16 kernels with the retained materialized planned
kernels on sparse inputs and completely empty blocks.  Public contract tests
execute 65/129-byte GF8 and 66/130-byte GF16 splits through low, high, and
forced-generic paths and verify both recovered bytes and the fixed-staging
scratch slope.  The calibrated-dispatch regression adds boundary/policy checks,
forced tiled-versus-materialized byte equality, scratch selection, and the
multi-item AUTO route.  The decode-plan schedule target cannot be linked with
TSan because its pre-existing test-only global allocation replacements conflict
with the TSan runtime; the public plan, high/low acceptance, and shared-context
tests cover the concurrent production entry points under TSan instead.

Forced generic remains the correctness fallback.
