# Low-rate coefficient-copy elimination

Low-rate encoding first interpolates the K systematic values, padded with
shortened zeroes to P = ceil_pow2(K), into one immutable P-shard coefficient
block.  It then evaluates that same block on every parity coset that contains a
requested output.

Before this change, every nonempty parity block began with a separate loop that
copied all P coefficient shards into a second P-shard evaluation workspace.
The normal in-place forward transform then reread and rewrote that workspace.
For B-byte shards and Q active parity blocks, the copy alone added Q * P * B
bytes read and the same number written.

The encoder now starts each parity transform with an out-of-place forward
butterfly supplied by the selected immutable backend table:

- P = 2 uses one two-way butterfly;
- P >= 4 uses P / 4 fused four-way butterflies for the first two radix-2
  layers where the field/backend/byte-size policy has qualified that fusion;
- other GF16 sizes use P / 2 out-of-place two-way butterflies for the first
  layer, followed by the established in-place second-layer schedule; and
- the remaining layers continue in place from exactly the state reached by the
  original transform schedule.

Each primitive loads coefficients directly into registers, performs the first
layer or layers, and stores the transformed values to the evaluation workspace.
There is no copy-then-transform compatibility wrapper.  The coefficient block
remains unchanged and is reused for all parity cosets.  Scalar, SSSE3, and AVX2
provide the same operation contract; platforms without a native vector table
use the portable scalar implementation.  The GF16 choice deliberately reuses
the existing `UseFusedButterfly4` thresholds: 64 transform bytes on every
backend, 128 only on AVX2, and the split schedule otherwise.  This avoids
reintroducing the all-size GF16 fusion that prior isolated evidence rejected.

Complete dense parity blocks now also execute their final forward layer
out-of-place into the caller's disjoint parity buffers.  This removes the
following scratch-to-output copy without changing the intermediate transform.
Sparse blocks and a partial final parity block retain scratch evaluation plus
scatter because not every final butterfly has a complete destination set.
For a ragged shard tail the direct destination is the existing fixed 64-byte
parity staging slot; the required compact public-tail gather remains.

This changes memory traffic and arithmetic scheduling only.  It does not change
the low-v1 field, active-parent basis, coordinate order, shortening,
puncturing, requested-output semantics, or parity bytes.  Public validation
already proves that every non-null parity destination is disjoint from inputs,
other outputs, scratch, and protected pointer metadata.  The direct final
layer therefore satisfies the backend's stricter out-of-place alias contract.

## Correctness and instrumentation

Backend startup known-answer tests cover the two- and four-way primitives with
read-only source checks, guard bytes, unaligned addresses, zero lengths, vector
tails, zero-skew sentinels, and the GF16 full-tile and compact-tail layouts.
The public low-encode differential tests compare GF8 and GF16 results against
the independent direct systematic generator oracle for P = 1, P = 2, P = 8,
R > P, sparse requested outputs, partial final parity blocks, and
non-vector-aligned byte counts.

Test-only counters report the first-layer and direct-final two- and four-way
calls made by the selected context, plus the number of dense blocks written
directly.  Direct and explicit-context tests require the expected qualified or
split dispatch.  Dense GF8/GF16 cases compare every output against the direct
systematic generator oracle and prove the expected full-block route for both
aligned and ragged passes.  The operation-count self-test independently checks
both field implementations for the former whole-P `memcpy` loop.  Its mutation
check reintroduces that loop and must be rejected.

## Performance status

The coefficient-copy traffic eliminated by the original change is exact:
2 * Q * P * B bytes across the old read and write pass.  The direct-final path
additionally removes one P * B read and one P * B write for every complete
dense parity block.  End-to-end promotion still depends on isolated, paired
measurements across shard sizes and backends; timings gathered while other
project workers are active are diagnostic only.  This change does not reduce
the current 2P encode scratch geometry, which is a separate tiling and
workspace problem.

The current-HEAD port was checked on 2026-07-18 against its untouched
2fce390c control.  Both GCC 13.3 Release builds were pinned to logical CPU 14
of an AMD Ryzen 9 9950X3D with OMP_NUM_THREADS=1.  Each value below is the
pooled median of 18 retained A/B/B/A samples, with four encode calls per
sample.  Other project workers were active, so these are directional
diagnostics rather than the authoritative Leopard-main comparison.  Every
candidate codeword had the same transmitted-parity digest as its control.

| K,R | Profile/field | Shard bytes | Control us | Candidate us | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| 8,248 | low/GF8 | 65,536 | 710.656 | 578.388 | 22.87% |
| 32,224 | low/GF8 | 65,536 | 968.490 | 805.687 | 20.21% |
| 100,128 | low/GF8 | 65,536 | 1191.721 | 1085.574 | 9.78% |
| 128,128 | low/GF8 | 65,536 | 1198.961 | 1130.363 | 6.07% |
| 32,224 | low/GF8 | 1,048,576 | 39211.182 | 37430.401 | 4.76% |
| 64,448 | low/GF16 | 65,536 | 3113.338 | 2825.431 | 10.19% |
| 127,129 | low/GF16 | 65,536 | 1940.769 | 1783.008 | 8.85% |
| 256,768 | low/GF16 | 65,536 | 10094.469 | 9795.004 | 3.06% |
| 100,300 | low/GF16 | 1,048,576 | 126567.543 | 112769.282 | 12.24% |
| 240,16 | high/GF8 | 65,536 | 880.691 | 882.170 | -0.17% |

The unchanged high-profile neighbor remained within noise.  Low-profile gains
were positive in every sampled cell, from 3.06% to 22.87%; the 1 MiB GF8 cell
is below the default 5% promotion threshold and therefore needs the isolated
coordinator run before any universal performance claim.

The following 2026-07-17 measurements are deliberately non-authoritative.  They
ran while other project workers were active, pinned to logical CPU 13 of an AMD
Ryzen 9 9950X3D, with GCC 13.3 Release builds, OMP_NUM_THREADS=1, 64 KiB shards,
15 median/MAD samples and eight encode calls per sample.  Control was commit
4070e4e527935026fb87593567587558f0a08d51; the candidate was the otherwise
identical working tree containing this change.

| K,R | Field/backend | Control us | Candidate us | Candidate time delta |
| --- | --- | ---: | ---: | ---: |
| 8,248 | GF8/AVX2 | 717.427 | 571.578 | -20.33% |
| 8,248 | GF8/SSSE3 | 1243.562 | 1009.136 | -18.85% |
| 64,192 | GF8/AVX2 | 1032.355 | 896.673 | -13.14% |

An earlier GF16/AVX2 K=100,R=156 diagnostic measured 1910.594 us for
the control and 1808.125 us for the candidate (-5.36%), but that run exercised
the now-rejected all-size four-way first layer.  It is retained as historical
diagnostic evidence and does not describe the current gated GF16 source.

A longer 31-sample scalar diagnostic at K=8,R=248 measured 4967.924 us for the
control and 4996.797 us for the candidate (+0.58%).  This is below the project's
2% neighboring-regime regression threshold but should be revisited in the
isolated promotion run; these concurrent diagnostics do not establish either a
promotion or a scalar regression.
