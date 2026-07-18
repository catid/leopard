# Leopard2 low-decoder reveal/scatter fusion

## Scope

This optimization changes only the byte-heavy final-output pass of the
`LEO2_PROFILE_LOW_V1` Algorithm 4 decoder.  It does not change the field,
coordinate map, locator construction, transform schedule, codeword, or public
API.  Materialized, side-tiled, and C1-pruned Algorithm 4 executions use the
same path.

## Arithmetic equivalence

For each requested systematic coordinate `i`, Algorithm 4 leaves an evaluator
value `E[i]` after its final forward LCH transform.  The established decoder
reveals the systematic value by multiplying by the inverse-locator factor

    output[i] = E[i] * alpha^(q - locator_log[i])

where `q` is 255 for GF8 and 65,535 for GF16.  Previously this was performed
in place in scratch.  Because the fixed-multiply backend uses restricted source
and destination pointers, the in-place helper first copied each 64-byte tile to
a temporary buffer.  The public decode layer then copied the revealed scratch
shard to the caller output.

For complete 64-byte tiles, Leopard2 now leaves `E[i]` unrevealed in scratch
and invokes the same fixed-multiply backend with `E[i]` as the source and the
caller output as the distinct destination.  Buffer validation already rejects
output overlap with inputs, other outputs, and scratch.  Consequently the
backend's non-aliasing contract holds, the field operation and factor are
unchanged, and both the temporary copy and final scatter copy disappear.

## Eligibility and fallback

The fused path is deterministic and is selected only when all of the following
are true:

- the specialized low-profile decoder is selected;
- the pass covers one or more complete 64-byte kernel tiles; and
- ordinary decode overlap validation has succeeded.

A ragged public tail remains on the established path: it is staged to one
64-byte kernel tile, revealed in scratch, and gathered back to the exact public
length.  This is required because the padded tail is larger than the public
destination and GF16's compact tail representation can differ from the kernel
layout.  Generic and high-profile decoder behavior is unchanged.

## Test evidence

Test-only counters distinguish low-profile shards revealed directly into the
destination from shards revealed in scratch.  Focused API cases cover:

- GF8 and GF16;
- AUTO, forced materialized, and forced tiled Algorithm 4 execution;
- one missing original and a complete missing systematic side;
- C1-pruned and full final-output transforms;
- aligned, tail-only, and aligned-plus-ragged lengths; and
- deliberately unaligned caller output buffers.

The complete-tile counter must equal the number of requested originals for
each aligned pass.  Tail-only cases must record no direct reveal, and a ragged
tail following an aligned prefix must record both one direct and one scratch
reveal per requested original.

The implementation was also checked through the full Release suite, strict GCC
warnings, Clang AddressSanitizer plus UndefinedBehaviorSanitizer, diagnostic
scalar/SSSE3/AVX2 builds, GF8/GF16 configuration builds, aliasing tests, and a
no-OpenMP ThreadSanitizer concurrent-plan run.  On this host, ThreadSanitizer
occasionally failed during runtime startup with `unexpected memory mapping`;
rerunning the affected no-OpenMP binaries produced clean results.  An
OpenMP-enabled GCC ThreadSanitizer initially reported races between consecutive
GF16 encoder parallel regions, outside this decode change.  The follow-up audit
in `docs/leopard2_openmp_tsan_audit.md` reproduced the same report with two
standards-valid consecutive OpenMP loops, verified the implicit barriers and
disjoint iteration ranges, and found no codec race.  Clang ThreadSanitizer with
its documented non-instrumented-OpenMP-runtime setting passed the focused
encoder and concurrency programs.

## Performance gate

The authoritative comparison used production Release builds with
`LEO2_BUILD_TESTS=OFF`, so test-only counter atomics were absent from both
binaries.  The control was immediate parent `789bd0c`; the candidate was
unconditional-fusion commit `7ff4807`.  Whole-codec decode execution was pinned
to CPU 0 with its sibling CPU 16 idle and `OMP_NUM_THREADS=1`.  Plan setup was
reported separately and execution reused the plan.  Runs used counterbalanced
control-candidate-candidate-control ordering.

Representative speedup is positive below; values are pooled medians from the
isolated production run.

| Path | Backend | 128 B, loss 1 | 128 B, max loss | 1 KiB, loss 1 | 1 KiB, max loss |
| --- | --- | ---: | ---: | ---: | ---: |
| tiled GF8, K=64 R=192 | AVX2 | -0.79% | +9.77% | +3.66% | +29.23% |
| tiled GF8, K=64 R=192 | SSSE3 | +0.95% | +3.79% | +10.54% | +13.27% |
| tiled GF8, K=64 R=192 | scalar | +1.21% | +4.11% | +0.82% | +7.00% |

The scalar binary also showed an unrelated approximately 7% encode shift, so
its small decode deltas are more layout-sensitive than the SIMD results.

For AVX2 K=17, R=33, larger-shard results were:

| Path | 4 KiB, loss 1 | 4 KiB, max loss | 64 KiB, loss 1 | 64 KiB, max loss |
| --- | ---: | ---: | ---: | ---: |
| tiled GF8 | +3.17% | +27.14% | +2.28% | +26.39% |
| materialized GF8 | +1.97% | +26.12% | +1.96% | +24.17% |
| tiled GF16 | +1.51% | +22.44% | +2.01% | +23.33% |
| materialized GF16 | +2.02% | +22.08% | +1.77% | +21.89% |

At 64 bytes, K=64/R=192 maximum-loss tiled GF8 improved by 4.42% on
AVX2, 3.84% on scalar, and 0.87% on SSSE3.  The worst materialized-GF8
neighbor was -1.22%, within the 2% neighboring-regression gate, and tiled GF16
neighbors ranged from +0.24% to +3.26%.  A proposed 128-byte tiled-GF8 gate was
tested and rejected because it was consistently equal or slower than
unconditional fusion.

One stable anomaly remains: K=8/R=16, GF8, SSSE3, 64-byte maximum-loss forced
transform was 4.68% slower in the new binary.  The same shift appeared in AUTO,
tiled, materialized, and untouched forced-generic schedules (approximately
4-7%), while disabling this fusion was still slower.  This points to
sub-microsecond function/link-layout sensitivity or a missing tiny-code direct
dispatcher path, not the fused multiply-to-destination arithmetic.  It is a
separate follow-up rather than a reason to select the inferior fallback.

Two earlier timing sets are deliberately excluded: a tests-enabled screen had
candidate-only test-counter overhead, and a later diagnostic accidentally ran
two benchmark sessions on CPU 0.  Neither contributed to the tables or the
promotion decision.  The isolated production evidence promotes unconditional
fusion for every complete low-profile kernel tile; ragged tails keep the
correctness fallback described above.

The authoritative machine-readable bundles from this host are:

- `/tmp/leo2_low_reveal_prod64_latin_isolated_20260718_0818`
  (324 JSON files plus `summary.tsv`);
- `/tmp/leo2_low_reveal_prod_rep_isolated_20260718_0824`
  (320 JSON files plus `summary.tsv`); and
- `/tmp/leo2_low_reveal_k8_ssse3_modes_20260718_0821`
  (96 JSON files plus `summary.tsv`).

The production benchmark binaries can be reproduced for either source tree
with:

    taskset -c 0 cmake -S SOURCE -B SOURCE/build/low-reveal-prod \
      -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=ON
    taskset -c 0 cmake --build SOURCE/build/low-reveal-prod \
      -j1 --target bench_leopard2

A representative execution command is:

    taskset -c 0 env OMP_NUM_THREADS=1 bench_leopard2 \
      --k 64 --r 192 --profile low --field gf8 --backend avx2 \
      --bytes 1024 --loss 64 --reuse 128 --iterations 21 --warmup 5 \
      --threads 1 --force-specialized --force-tiled \
      --retain-samples --json result.json
