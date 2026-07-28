# GFNI affine multiplication inside the Leopard2 codec

This experiment is the codec-integration follow-up that
`docs/leopard2_gfni_affine.md` asked for.  The earlier experiment established
that `VGF2P8AFFINEQB` computes Leopard's GF8 fixed multiplication correctly and
is much faster than nibble tables on a standalone chain, but it did not put the
instruction behind the codec's dispatch or run whole-codec benchmarks.

This directory does that.  It is an evaluation, not a promotion: the production
integration requirements are listed at the end of
`docs/leopard2_gfni_codec.md`.

## What the candidate changes

`Leopard2BackendAVX2.cpp` gains a `LEO2_GFNI_VARIANT` compile-time variant.
When it is defined:

- every fixed multiplication in the vector kernels is one
  `_mm256_gf2p8affine_epi64_epi8` instead of two `vpshufb` plus the AND/shift/XOR
  around them;
- the GF16 16-by-16 bit multiplication matrix is applied as its four 8-by-8
  blocks over Leopard's split low-byte/high-byte block layout, so a 64-byte
  GF16 tile costs four affine operations and two XORs instead of four ANDs, two
  shifts, eight shuffles and eight XORs;
- the nibble-table storage shape is reused unchanged, with each 16-byte row
  holding one affine matrix duplicated so a 128-bit broadcast fills all four
  64-bit lanes.  No vector call site changes;
- sub-vector scalar tails evaluate the same stored matrix bitwise, so the
  variant needs no second table.

The vector data path stays 256 bits wide and the backend identity, dispatch,
selectors, and wire profile are untouched.

## Building the evaluation binaries

The variant is reached with compiler flags only; no CMake graph change is
required, so the build graph contract test still passes.

    cmake -S . -B /tmp/l2-stock -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF
    cmake --build /tmp/l2-stock -j "$(nproc)" --target bench_leopard2

    cmake -S . -B /tmp/l2-gfni -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
        -DCMAKE_CXX_FLAGS="-mgfni -DLEO2_GFNI_VARIANT=1"
    cmake --build /tmp/l2-gfni -j "$(nproc)" --target bench_leopard2

A tests-enabled GFNI build additionally fails `leopard2_portable_isa`, which is
correct: that audit forbids GFNI inside the object labelled as the AVX2 member.
A production integration must give GFNI its own runtime-qualified member with
its own allowed-instruction list rather than relaxing that audit.

## Scripts

- `kernel_screen.cpp` — standalone kernel comparison of the nibble and affine
  forms at 256 and 512 bits, over an IFFT-shaped butterfly sweep.  It derives
  the nibble tables and the affine operand from one shared GF(2)-linear map, so
  agreement between the kernels also pins the GFNI operand bit order on the
  host.  Build with
  `g++ -O3 -std=c++17 -march=x86-64-v3 kernel_screen.cpp -o kernel_screen`.
- `cross_build_differential.py` — randomized cross-build differential.  Runs
  random `(K, R, bytes, loss, field, backend, threads, batch)` shapes through
  both binaries and requires identical original/parity/recovered digests, an
  identical round-trip flag, and identical accept/reject behaviour.
- `backend_ab.py` — within-binary clustered ABBA between two `--backend`
  selections.
- `main_compare_screen.py` — clustered ABBA screen of a Leopard2 binary against
  the exact Leopard main adapter built by
  `experiments/leopard2/main_compare`.  This is a diagnostic screen, not
  promotion evidence: it does not hold the CPU-pair lease, does not verify
  build closure, and does not check sibling idleness.  Run it pinned to one
  logical CPU with the machine otherwise quiet; running many cells in parallel
  makes memory-bandwidth contention dominate the ratio.

## Retained results

See `results/` and `docs/leopard2_gfni_codec.md`.

## Additional scripts added 2026-07-24

- `binary_ab.py` — clustered-ABBA A/B between two Leopard2 benchmark binaries.
  Used for candidate-versus-`HEAD` comparisons where both sides are Leopard2.
- `fixed_cost_attribution.py` — encode size sweep against the exact-main
  adapter that separates fixed per-call cost from byte-proportional cost by
  least squares.  `results/fixed_cost_attribution.txt` retains the run that
  established that Leopard2's byte-proportional encode cost already equals or
  beats exact main in every measured shape, and that the whole remaining
  small-shard deficit is fixed per-call cost.

To attribute a specific fixed cost, build with the component stubbed out and
diff against the ordinary build; the encode-validation attribution used a
temporary early `return LEO2_SUCCESS` in `ValidateEncodeBuffers` behind a
measurement-only macro.  Do not leave such a macro in the tree.
