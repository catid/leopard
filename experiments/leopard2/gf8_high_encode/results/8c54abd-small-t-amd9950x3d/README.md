# GF8 tiny-redundancy high encoder

This checkpoint qualifies the direct-input AVX2 legacy-high encoder for
redundancy transform sizes `T=2` and `T=4`.  The arithmetic, active-parent
coordinates, field representation, parity order, shortening, and puncturing
are unchanged.  The implementation commit is `8c54abd`; the conservative
production selector and its correctness tests were frozen in descendant
`db725de35eb021357b88fc8cdf254722bac0ee71`.

The mature `T<=4` path copied every message shard into a `T`-shard temporary,
inverse-transformed it, and XORed it into the parity accumulator.  The new
callback writes the first complete block directly into the accumulator and
inverse-transforms later complete blocks directly into XOR accumulation.  Only
the final `K mod T` message shards, if any, use the upper scratch half.  For a
qualified pass of `B` bytes this removes `K-(K mod T)` shard copies, or
`2 * (K-(K mod T)) * B` bytes of copy read/write traffic.  Scratch size and wire
bytes do not change.

## Production selector

The callback is selected only for the AVX2 GF8 legacy-high transform path, a
contiguous requested parity prefix, and these measured regions:

- `T=2`, `16 <= K <= 64`;
- `T=4`, `K >= 16`;
- at least 64 KiB per execution pass, or at least 4 KiB when `K >= 64`.

Ragged 64-byte staging passes, sparse/holey output masks, other backends, `K=8`,
and `T=2` above `K=64` retain the mature fallback.  This deliberately excludes
the small cell that missed the promotion threshold and the large-`K` `T=2`
cell whose isolated timings were inconclusive.

## Exact Leopard-main result

The independent baseline is detached Leopard commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  The candidate is clean commit
`8c54abd6a74eadae4267b94f1f1c63c519bada59`, built without test hooks.  Every
accepted cell used launch affinity `0,5,21`, pinned timings to CPU 5, reserved
SMT sibling 21, acquired the canonical pair lease, and observed zero non-idle
jiffies on the sibling.  Each cell used three independent ABBA rounds, nine
samples per child, and exact matching input, parity, recovery, build, and source
identities.

Ratios are exact-main time divided by Leopard2 time; values above one favor
Leopard2.

| K | R | shard bytes | speedup | 95% CI |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 2 | 65,536 | 1.466x | 1.443-1.489 |
| 16 | 4 | 65,536 | 1.200x | 1.196-1.204 |
| 64 | 2 | 4,096 | 1.189x | 1.165-1.213 |
| 64 | 2 | 65,536 | 1.471x | 1.440-1.503 |
| 64 | 4 | 4,096 | 1.155x | 1.118-1.192 |
| 64 | 4 | 65,536 | 1.297x | 1.276-1.319 |
| 240 | 4 | 65,536 | 1.292x | 1.273-1.310 |

All accepted lower confidence bounds exceed the default 1.05 promotion gate.
The compact machine record in `summary.json` binds every full ignored manifest
by its internal digest and file SHA-256.  It also records the rejected `K=8,
R=4` diagnostic and both fail-closed `K=240,R=2` launches; neither is used as
positive evidence.

## Correctness and safety

Release tests compare the new route with the independent direct systematic
generator, exact legacy parity, scalar, SSSE3, AVX-512, and the retained AVX2
fallback.  They cover a 65-byte backend KAT and one-byte ragged execution tail,
unaligned buffers, a partial final message block, a T=4 parity prefix, a holey
output mask, source immutability, evidence-derived selector boundaries, and
concurrent immutable-codec use.  The backend startup KAT separately exercises
zero-skew sentinels and guards around every buffer.

Focused Release validation on the selector descendant passed the following
nine tests:

    leopard2_auto_encode_backend
    leopard2_context_backends
    leopard2_direct_encode
    leopard2_legacy_golden
    leopard2_api
    leopard2_portable_isa
    leopard2_backend_ops
    leopard2_transform_differential
    leopard2_encode_concurrency

Clang 18.1.3 ASan+UBSan passed 8/8 focused tests with leak detection and
non-recovering undefined behavior.  The direct generator-oracle case completed
under instrumentation rather than being replaced by a self-comparison.  Broad
fuzzing was intentionally outside this optimization-focused checkpoint.

## Reproduction

Build the production candidate from commit `8c54abd`:

    cmake -S . -B /tmp/leopard2-small-t-production \
      -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=ON -DLEO2_BUILD_FUZZERS=OFF \
      -DLEO2_ENABLE_CUDA=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /tmp/leopard2-small-t-production -j "$(nproc)"

Build exact main through the standalone adapter:

    cmake -S experiments/leopard2/main_compare \
      -B /tmp/leopard-main-small-t -DCMAKE_BUILD_TYPE=Release \
      -DLEOPARD_MAIN_SOURCE_DIR=/path/to/detached/exact-main
    cmake --build /tmp/leopard-main-small-t -j "$(nproc)"

Create a canonical held reservation for an allowed physical-core pair, then
run one cell at a time.  This example uses the 64-by-4 cell and must be adapted
to the allowed CPU set on the reproduction host:

    taskset -c 0,5,21 python3 \
      experiments/leopard2/main_compare/run_abba.py run \
      --baseline /tmp/leopard-main-small-t/leopard_main_benchmark \
      --candidate /tmp/leopard2-small-t-production/bench_leopard2 \
      --baseline-archive /tmp/leopard-main-small-t/libleopard_main_exact.a \
      --candidate-archive /tmp/leopard2-small-t-production/libleopard.a \
      --baseline-build-dir /tmp/leopard-main-small-t \
      --candidate-build-dir /tmp/leopard2-small-t-production \
      --baseline-source-root /path/to/detached/exact-main \
      --candidate-source-root "$PWD" --candidate-commit 8c54abd \
      --candidate-mode auto --reservation-file /path/to/reservation.json \
      --output /tmp/leopard2-small-t-k64-r4 --cpu 5 \
      --reserved-sibling 21 --cell k64-r4-64k:64:4:65536:1:43 \
      --reuse 8 --iterations 9 --warmup 2 --timeout 600

Verify the resulting full manifest while build and source inputs remain
available:

    python3 experiments/leopard2/main_compare/run_abba.py verify \
      --manifest /tmp/leopard2-small-t-k64-r4/manifest.json
