# GF8 T=32 two-tail high encoder

This checkpoint qualifies a direct systematic-generator-column path for the
legacy-compatible GF8 high encoder at `K=34,R=32`.  It does not change the
field, active parent, coordinate map, shortening, puncturing, parity order, or
wire bytes.  The source base is
`dc3d09e6aa3af748474813d7947ff64c417868f9`.

For this shape, the first 32 sources already form one complete inverse
transform block.  The mature fallback stages the remaining two sources as a
ragged second block and transforms it.  The new path encodes the complete block
through the retained transform, then accumulates the two precomputed Lagrange
generator columns directly into parity space.  A pure-AVX2 callback handles two
sources and two outputs per call, rereading each source once per output pair.
This reduces direct-tail traffic from eight to six shard streams per output
pair (four to three per output) relative to a one-output fused callback.

## Production selector

The path is selected only when all of the following hold:

- backend is explicitly resolved to AVX2;
- field/profile are GF8 and `legacy_high_v1`;
- `T=32`, `K=34`, and `R=32`;
- all 32 transmitted parity outputs are requested;
- the execution tile is at least 16,384 bytes.

Sparse parity requests, `K=33`, `K=35`, `R=31`, smaller execution tiles, and
all other backends retain their previous path.  The existing `K=T+1` direct
column selector is unchanged.

## Performance evidence

Measurements used GCC 13.3.0 on an AMD Ryzen 9 9950X3D.  Both binaries were
generic-tuned pure AVX2 builds with AVX-512 disabled and zero EVEX instructions.
ABBA timings were pinned to CPU 7 while SMT sibling 23 was reserved.  Accepted
invocations observed zero sibling work; contaminated invocations in the
exact-main run were rejected and retried.  Ratios are control time divided by
candidate time, so values above one favor this change.

The source-identical boundary comparison used four ABBA rounds and eight
observations per cell:

| K | R | bytes | selector | encode speedup | 95% CI |
| ---: | ---: | ---: | :--- | ---: | :--- |
| 34 | 32 | 16,383 | off | 1.007x | 0.992-1.026 |
| 34 | 32 | 16,384 | on | 1.147x | 1.139-1.156 |
| 34 | 32 | 16,385 | on | 1.132x | 1.125-1.138 |
| 34 | 32 | 32,768 | on | 1.183x | 1.172-1.200 |
| 34 | 32 | 65,536 | on | 1.187x | 1.178-1.199 |

This discontinuity supports the conservative 16 KiB cutoff: the immediately
smaller cell is neutral, while the lower confidence bound at the selected
boundary exceeds the 1.05 promotion rule by a wide margin.

The independent Leopard1 baseline is exact `main` commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, built through the standalone
adapter.  Three ABBA rounds at 64 KiB produced:

| K | R | Leopard2/main encode speedup | 95% CI |
| ---: | ---: | ---: | :--- |
| 33 | 32 | 1.410x | 1.314-1.512 |
| 34 | 32 | 1.312x | 1.301-1.322 |
| 35 | 32 | 1.114x | 1.112-1.117 |
| 34 | 31 | 1.113x | 1.102-1.124 |

Input, transmitted-parity, and recovered-original digests matched exactly in
every comparison.  The final control, candidate, and exact-main executable
SHA-256 values are recorded in `summary.json`.

## Correctness and safety

The backend startup KAT covers byte counts `0,1,7,31,32,33,63,64,65,257`,
unaligned vector tails, guards, source immutability, and aliased read-only
source ranges.  The encoder test independently constructs both full-parent
systematic Lagrange columns and checks exact parity at 65, 4097, 16,383,
16,384, and 16,385 bytes.  It proves the selector boundary, checks caller-input
immutability, and verifies that a selected-size sparse output mask
`{0,15,31}` uses the mature fallback while preserving unrequested sentinels.

All 64 coefficients in the two T=32 columns are nonzero and unequal to one,
so zero/one multiplier specialization has no eligible operation.  The emitted
AVX2 vector loop has one GCC-generated stack spill/reload for the eighth nibble
table per 32-byte iteration.  A no-spill rewrite was not pursued because the
measured two-output kernel already improves whole encode by about 5-6% over
the one-output fused kernel and replacing one hot L1 load with table/address
reload pressure has no credible further 5% whole-codec bound.

The final static archive grows by 5,008 bytes.  The three affected production
objects add 3,302 text bytes, 8 data bytes, and 128 BSS bytes.  No AVX-512/EVEX
instructions occur in the candidate executable.

The final Release suite passed 96/96 tests.  Clang 18.1.3
ASan+UBSan+LSan passed 8/8 focused tests with leak detection enabled and
undefined-behavior recovery disabled; the independent direct generator-oracle
test completed under instrumentation.  Exact results are bound in
`summary.json`.

## Rejected extensions

The one-output two-source prototype was retained only as a stepping stone; the
two-output callback was consistently about 5-6% faster.  A three-tail direct
path was rejected: T=16 was 0.857-0.923x and T=32 was 0.926x, 0.956x, and
1.009x at 4 KiB, 16 KiB, and 64 KiB.  Because direct traffic worsens at larger
tail counts, this measured result closes `L>=3`.  T=8, T=16, and T=64 two-tail
regions also missed the promotion confidence threshold and remain on the
mature fallback.

## Reproduction

Build the pure-AVX2 candidate and run correctness tests:

    cmake -S . -B build/pure-avx2 -G Ninja \
      -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON \
      -DLEO2_FLAG_MAVX512F=FALSE -DLEO2_FLAG_MAVX512BW=FALSE \
      -DLEO2_FLAG_MAVX512VL=FALSE
    cmake --build build/pure-avx2 -j "$(nproc)"
    ctest --test-dir build/pure-avx2 -j "$(nproc)" --output-on-failure

The ignored raw manifests named in `summary.json` contain all retained samples,
digests, affinity observations, and rejected contaminated invocations.
