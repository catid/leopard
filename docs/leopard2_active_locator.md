# Leopard2 active-parent locator setup

Leopard2 now evaluates dense erasure locators on the active power-of-two parent
instead of transforming all 256 or 65,536 ambient-field entries.  The retained
full-field implementation remains the differential oracle; sparse patterns still
use direct products when their measured cost is lower.

## Algorithm and storage

In Cantor coordinates, field addition is coordinate XOR.  Locator logarithms are
therefore an XOR convolution of the erasure indicator with `log(omega_i)`, with
the zero/self term omitted to produce the locator derivative at an erased
coordinate.  An active `N`-point Walsh transform needs the normalization

    N^-1 = 2^m / N  modulo (2^m - 1).

The full derivation and its relationship to the legacy full-field transform are
in `docs/leopard2_math_and_sources.md`.

The scaled transform of the fixed log kernel is precomputed once for every
proper parent size.  Their packed lengths sum to `field_order - 2`: 254 bytes in
GF8 and 131,068 bytes in GF16.  Plan setup transforms its existing `N`-entry
locator output in place, so dense setup has no additional shard-independent
temporary.  The old production path used a 256-entry or 65,536-entry automatic
array regardless of active `N`.  The full-field array now exists only in the
retained reference oracle.

Permanent locator contributions are cached only for the sparse direct region.
This is complete for current plan construction: deterministic selection retains
exactly `K` of `K+R` transmitted coordinates, so the dynamic real-plus-virtual
erasure count is always exactly `R`, the value used by codec setup's crossover
decision.  Coordinate- and transform-domain dense cache prototypes both retain
both Walsh transforms, add an `N`-entry modular pass, and require `N` persistent
field symbols.  The derivation, operation model, and negative promotion result
are in `docs/leopard2_permanent_locator_cache.md`.

GNU `size` on otherwise equivalent Release `bench_leopard2` executables reports
131,328 additional BSS bytes (the packed tables plus alignment), 7,480 additional
text bytes, and no data-segment change.  This is a process-wide fixed footprint;
codec and decode-plan object sizes do not grow with the ambient field order.

## Deterministic direct/Walsh crossover

Pinned whole-locator sweeps selected conservative count-only thresholds.  A
count at or below the table uses direct products; the next count uses the active
Walsh path.  These are execution choices only and cannot change the wire profile.

| Field | Active parent N | Largest direct erasure count |
| --- | ---: | ---: |
| GF8 | 2-8 | N |
| GF8 | 16 | 8 |
| GF8 | 32 | 9 |
| GF8 | 64 | 8 |
| GF8 | 128 | 9 |
| GF8 | 256 | 7 |
| GF16 | 2-32 | N |
| GF16 | 64 | 34 |
| GF16 | 128 | 24 |
| GF16 | 256-512 | 16 |
| GF16 | 1,024-65,536 | 14 |

The benchmark target `bench_leopard2_locator` compares direct, active-Walsh, and
full-field-oracle results before timing either candidate.  It reports median and
MAD, retains every raw sample, alternates candidates in ABBA order, records the
compiler/source/affinity/OpenMP identity, and exposes the production dispatch
decision.  A parser-backed CTest validates this JSON contract.

## Correctness evidence

`leopard2_locator_test` covers:

- every erasure subset of every active parent in test-only GF(2^4), including
  the active-transform normalization;
- every GF8 erasure subset through `N=16` against direct field algebra;
- direct, active-Walsh, dispatched, permanent-contribution, and full-field
  results across every GF8 and GF16 parent size;
- sampled scalar `O(E)` sums at dense GF16 sizes where producing every direct
  coordinate would dominate the suite, including modulo-65,535 wrap;
- permanent-mask reuse with permanent coordinates both included in and excluded
  from the dynamic erasure bitmap;
- exhaustive disjoint permanent/dynamic GF(2^4) decompositions through `N=8`,
  deterministic `N=16` decompositions, and GF8/GF16 add/remove/swap/dense
  coordinate- versus transform-cache prototypes;
- the exact-`R` dynamic/virtual selection invariant across 9,647 valid small
  public presence patterns;
- the measured dispatcher boundary for every parent size; and
- 256 concurrent preparations of the same immutable pattern.

The focused test compares 9,894,532 locator entries.  The complete release test
configuration passed all 45 tests.  At `ctest -j30`, the tests shared the 30
allowed logical CPUs and the suite completed in 5.69 seconds; that saturated run
is correctness evidence, not a latency measurement.

The focused test also passes Clang 18 ASan+UBSan with leak detection and
halt-on-undefined-behavior enabled.  Clang TSan passes the 256-preparation
concurrency case with `ignore_noninstrumented_modules=1`, as required by the
host's non-instrumented `libomp`; without that documented runtime option TSan
reports a `libomp` mutex-initialization race rather than a Leopard2 access.

## Pinned setup results

Measurements used an AMD Ryzen 9 9950X3D, GCC 13.3, Release builds, scalar codec
backend, 64-byte shards, one codec thread, CPU 14 affinity, 101 samples, and nine
warmups.  Baseline commit `cd685f5` used the full-field dense fallback.  Setup is
reported separately from execution.

| Field/profile | K,R / active N | Missing originals | Baseline plan setup | Active-N plan setup | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| GF16 low | 128,128 / 256 | 4 | 14.751 us | 2.260 us | 6.53x |
| GF16 low | 64,192 / 256 | 32 | 22.351 us | 2.360 us | 9.47x |
| GF16 high | 200,50 / 512 | 25 | 15.180 us | 3.770 us | 4.03x |
| GF16 high | 1000,200 / 2,048 | 100 | 194.188 us | 12.820 us | 15.15x |
| GF16 high | 4096,512 / 8,192 | 256 | 211.793 us | 49.361 us | 4.29x |
| GF8 low control | 64,192 / 256 | 32 | 1.270 us | 1.260 us | 1.01x |

The GF8 full-field control is intentionally unchanged within measurement noise.
Fresh-process A/B/B/A timing of `leo_init()` (101 processes per run on the same
pinned CPU) measured 10.283 ms for the baseline and 10.297 ms with the packed
active tables, a 0.13% difference below the run-to-run MAD.  The plan-setup gain
does not come from moving a visible cost into startup.

Pinned end-to-end plan execution was unchanged: the representative high-profile
decode averaged 1,467.632 us on the baseline and 1,468.227 us with active-parent
setup (a 0.04% difference).  The optimization changes plan construction, not the
byte-heavy decode schedule.

The direct/active crossover used nine deterministic erasure layouts at every
candidate boundary, with 21 or 31 samples per layout and calls per sample scaled
by parent size.  The provisional single-layout thresholds were rejected after
this sweep exposed branch-pattern bias.  Final cutoffs select the direct path
through the last ambiguous region and the active path only in a multi-layout
winning region.  At every first active cell the active path won all nine frozen
layouts by at least 13% on the calibration host.  These deterministic built-in
cutoffs are deliberately more conservative than that host's median crossover;
other architectures should retain correctness but still need their own offline
performance calibration.  Raw machine-readable results are generated under
`.research/leopard2/locator-active/` and remain ignored because they are
machine-specific.  The frozen v2 manifest contains 360 runs (40 boundary cells
times nine layouts), identifies clean implementation commit `dbca4fa`, retains
every ABBA sample, and has SHA-256
`e31ee8c8757acf13de9efeb79d4380e6ef875c65febb6192073d7651ffbc0ac6`.
Earlier single-order manifests are exploratory and are not promotion evidence.
The initialization A/B/B/A manifest has SHA-256
`2c2678e4d5cf1f8f899154c55f4b7fc5ba38f5ba082b0c346e37c130447f08e4`.
The frozen run recorded `amd-pstate-epp`, the `powersave` scaling governor with
`performance` energy preference, and isolated CPU 14; no other benchmark job
ran concurrently.

Representative reproduction commands:

    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$(nproc)"
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
      ctest --test-dir build/release \
      -R 'leopard2_locator($|_benchmark_smoke)' --output-on-failure
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE taskset -c 14 \
      build/release/bench_leopard2_locator \
      --field gf16 --n 2048 --erasures 15 --calls 64 \
      --iterations 31 --warmup 3 --seed 7

## GF8 AVX2 cold-plan screen (negative, 2026-07-22)

A bounded follow-up tested whether the remaining GF8 active-parent work should
move into the AVX2 backend.  The prototype formed the union erasure bitmap with
AVX2 and evaluated both modulo-255 Walsh transforms plus their pointwise
log-kernel product in 16-bit lanes.  It was context-specific: scalar and SSSE3
contexts retained the scalar routine, and the AVX2 translation unit remained
compiled with `-mavx2 -mno-avx512f`.

The arithmetic idea worked.  A same-binary ABBA microbenchmark measured the
complete locator kernel at 2.81x to 3.60x the scalar active-parent routine for
`N=32,64,128,256`.  The prototype matched the scalar/full-field oracles for
every basis erasure and 1,024 dense patterns at every one of those sizes, plus
permanent/dynamic unions: 579,040 additional compared locator entries.  A
production round trip also matched, and disassembly found zero EVEX-prefixed
instructions in the AVX2 object.

It was not promoted because the public first-use result missed the gate.  Five
paired ABBA rounds covered `(K,R)=(16,8),(16,10),(128,64),(192,64)`, one and
maximum loss where valid, 2 and 4 KiB shards, and reuse counts 1, 8, and 64.
For shapes that actually selected Walsh setup, plan construction improved by
1.067x--1.564x, but decode including amortized setup ranged from 0.992x to
1.026x by cell.  At reuse one its median was 1.011x and its best cell was only
1.026x, below the required 1.05x.  Neighbor checks at `K=15,17,127,129,191`
and the direct/Walsh `R=9/10` boundary likewise produced no first-use gain over
1.017x.  The two transforms are only about 0.19 us after vectorization; byte
execution and the rest of immutable-plan construction dominate.  Roughly 300
lines of backend/API/test plumbing are not justified by that end-to-end result.

The frozen target and neighbor artifacts are:

- `experiments/leopard2/locator_construction/results/avx2_walsh_negative_20260722.json`
  (SHA-256 `0dfa61b0fa7d589f7c8b855e3d9982955dd0648aba174c49594df5516c508c42`);
- `experiments/leopard2/locator_construction/results/avx2_walsh_neighbors_20260722.json`
  (SHA-256 `302b3d8eaf2ed6583a972d8664488ad878d1c48b173e304902ad9ba2132b652b`).

Both identify the exact frozen control and candidate binaries.  The control is
clean source commit `f85970de031a6aeb74e437a34b5b92b6588e398f`; CPU 14 was
pinned with sibling 30 idle, and the global benchmark lock excluded concurrent
compilation and timing.  This is a negative production result for GF8 AVX2
cold-plan acceleration, not a rejection of the vector arithmetic itself on a
future workload whose setup dominates execution.
