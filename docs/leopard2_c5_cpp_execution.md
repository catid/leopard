# Leopard2 C5 C++ dyadic-block execution checkpoint

## Outcome

The bounded C++ candidate is correct, default-off, and rejected for production
promotion.  It preserves the existing padded-parent wire map for every claimed
case, but no tested `(field,q)` region met the required credible 10% speedup
without neighboring regressions.  One isolated GF8 `q=65`, 1 MiB cell reached
a 14.48% conservative gain; the other ten `q=65` cells had a median 7.47%
regression and eight exceeded the 2% regression limit.  No dispatcher entry or
new profile is justified by that outlier.

The standalone candidate, fail-closed merger, unit tests, and retained results
are under:

    experiments/leopard2/c5_dyadic_cpp/

They are not referenced by `CMakeLists.txt`, the public API, the production
dispatcher, or a serialized code identifier.  A normal CMake configure/build
therefore compiles none of this experiment.

## Construction and compatibility

For a nonzero normalized-LCH coefficient prefix of length `q`, split

    q = b_0 + b_1 + ... + b_s

into the canonical decreasing powers of two.  Block `j` begins at aligned
offset `o_j`.  Because the set bits of `o_j` and `t < b_j` are disjoint,

    X_(o_j+t)(x) = X_o_j(x) X_t(x).

On an aligned output coset of size `b_j`, `X_o_j` is constant.  The executor
therefore copies the block's coefficients into a complete production
`b_j`-point shifted LCH transform, multiplies by the precomputed coset factor,
and XOR-accumulates into the existing parent outputs.  Factors zero and one
are skipped/specialized.  The field polynomial, Cantor basis, symbol layout,
coordinate order, normalizers, and transforms remain the legacy Leopard ones.

The strongest modeled special case is `q = 2^a + 1`.  Here `X_(2^a)` is zero
on the lower coset and one on the upper coset.  Instead of materializing a
one-shard transform and an output-wide accumulation, the executor XORs the
tail coefficient into coefficient zero before the upper shifted transform.
This is the implementable fused ceiling requested by the earlier scalar/GF16
checkpoint, not the previously rejected materialized route.

The candidate's block factors are derived with an independent legacy-field
implementation and the linearized subspace-polynomial recurrence.  Each
normalizer and active subspace value is cross-checked against production test
hooks.  Complete candidate outputs are compared byte-for-byte with one padded
production parent transform.  Selected output symbols are additionally
evaluated directly as

    y(omega_p) = sum_(i=0)^(q-1) c_i X_i(omega_p),

so the production transform is not its only oracle.

This checkpoint implements no exact coordinate profile.  The exact-prefix
Lagrange/LCH algebra and small-field MDS evidence in
`docs/leopard2_c5_dyadic_blocks.md` remains valid, but the C++ experiment makes
no exact-size or legacy-systematic-encoder claim.  Exact encoders remain C6-C8
work and require a new serialized code identity.

## Correctness and safety evidence

The deterministic projection is identical across the portable auto build and
the forced scalar, SSSE3, and AVX2 libraries:

| Evidence | Result |
| --- | ---: |
| Cases per backend | 144 |
| Backend case executions | 576 |
| ASan+UBSan case executions | 144 |
| GF8 prefix values | 15 |
| GF16 prefix values | 8, through `q=8191` |
| Valid shard sizes | 1, 17, 63, 64, 65, 129, 4095, 4096, 4097 bytes |
| Direct normalized-LCH symbol checks per build | 388 |
| Independent/production factor checks per build | 1,260,048 |
| Deterministic correctness digest | `0x56af7bbc2fb6a888` |
| Candidate/padded compared output | byte-identical, including zero-padded kernel tails |
| Guard regions | intact for input, output, and candidate scratch |
| ASan/UBSan | pass, leak detection enabled, no diagnostic |

Legacy kernels operate in 64-byte quanta.  For a nonmultiple tail the harness
allocates and zeros the next full quantum, compares both the valid bytes and
the padded bytes, and places a 64-byte canary after every physical shard.  This
is explicit experimental staging, not a claim that the old transform hook
accepts an unpadded allocation.

The C++ source builds warning-clean under GCC 13.3 and Clang 18 with
`-Wall -Wextra -Wpedantic -Werror`.  The checkpoint merger has eight fail-closed
unit tests covering matrix shape, odd ABBA sample counts, content bindings,
the promotion rule, and neighboring-regression rejection.

## Execution, setup, scratch, and traffic

The candidate plan is immutable and execution performs no allocation.  The
reported setup samples include factor derivation and the experiment's
independent consistency checks, so they are conservative diagnostic setup
times rather than a proposed production-plan cost.  Across the retained matrix
their medians range from 16.50 to 738.11 microseconds.

Caller-owned input and output shards are excluded from scratch.  In the fused
`2^a+1` cases the candidate needs only `2^a` pointer slots; the padded baseline
needs `2^(a+1)` pointer slots.  Thus representative candidate/baseline scratch
is 256/512 bytes at `q=33`, 1,024/2,048 bytes at `q=129`, and 4,096/8,192 bytes
at `q=513` on this 64-bit host.  General multi-block cases also need a reusable
largest-tail-block workspace and retain scalar nontrivial-factor arithmetic;
they are correctness evidence, not a SIMD promotion candidate.

The logical traffic model includes input staging, output initialization,
butterfly reads/writes, tail injection, scratch staging, factor application,
and accumulation.  For the fused benchmark geometries it predicts fixed
reductions of 10.87%, 9.57%, 8.52%, 7.67%, and 6.97% as `q` grows from 33 to
513.  The measured result shows why operation/traffic counts cannot promote a
kernel: the existing padded implementation's fused radix schedule and
single-call traversal offset those logical savings.

## Pinned AVX2 timing

The authoritative run used the forced AVX2 library on an AMD Ryzen 9 9950X3D,
one allowed physical CPU (`taskset -c 15`), `OMP_NUM_THREADS=1`, with sibling
CPU 31 idle.  The host exposed CPUs 0-31, 16 physical cores, one socket and one
NUMA node.  The readable governor was `powersave`; the sampled idle frequency
before the run was 600 MHz and boost remained firmware-controlled.  No other
memory-intensive experiment was active during the short pinned phase.

There are 54 alternating-order ABBA cells.  Shard sizes are 64 B, 1 KiB,
64 KiB, and 1 MiB; batch and plan-reuse counts are 1 and 8 where bounded memory
permits.  Setup is outside execution.  Each reported credible gain is

    100 * ((padded_median - padded_MAD) /
           (candidate_median + candidate_MAD) - 1).

| Field | q / parent | Cells | Modeled traffic reduction | Credible gain min / median / max | Cells >=10% | Disposition |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| GF8 | 33 / 64 | 11 | 10.87% | -21.54 / -11.55 / +1.65% | 0 | reject |
| GF8 | 65 / 128 | 11 | 9.57% | -19.09 / -7.47 / +14.48% | 1 | reject: isolated 1 MiB cell |
| GF8 | 129 / 256 | 11 | 8.52% | -17.32 / -12.96 / -2.80% | 0 | reject |
| GF16 | 257 / 512 | 11 | 7.67% | -14.82 / -7.92 / +4.25% | 0 | reject |
| GF16 | 513 / 1024 | 10 | 6.97% | -16.97 / -5.90 / +1.19% | 0 | reject |

Representative execution medians (microseconds per stripe) are:

| Field, q | Bytes | Batch / reuse | Candidate | Padded | Credible gain |
| --- | ---: | ---: | ---: | ---: | ---: |
| GF8, 33 | 1 KiB | 1 / 8 | 6.363 | 5.653 | -11.55% |
| GF8, 33 | 1 MiB | 1 / 1 | 8,327.9 | 8,934.5 | +1.65% |
| GF8, 65 | 1 KiB | 1 / 8 | 13.475 | 12.528 | -7.37% |
| GF8, 65 | 1 MiB | 1 / 1 | 19,907.6 | 23,250.8 | +14.48% |
| GF8, 129 | 1 KiB | 1 / 8 | 29.423 | 25.854 | -12.51% |
| GF8, 129 | 1 MiB | 1 / 1 | 51,304.0 | 50,239.6 | -2.80% |
| GF16, 257 | 1 KiB | 1 / 8 | 83.025 | 76.812 | -7.92% |
| GF16, 257 | 1 MiB | 1 / 1 | 119,574.1 | 131,236.8 | -0.69% |
| GF16, 513 | 1 KiB | 1 / 8 | 184.647 | 173.933 | -6.02% |

The 1 MiB `q=65` point is real in this run but does not form the required
region.  Promoting it would add a special executor that loses most small and
medium cells, duplicates C1/C2 scheduling machinery, and has no evidence on a
second host.  C5 therefore closes as a measured negative result.  If a future
C1/C2 executor can fuse block boundaries without the duplicate traversal, it
must be evaluated independently; this result must not be reinterpreted as a
production threshold.

## Reproduction

All commands below use at most eight build jobs.  The raw artifacts bind the
standalone source, baseline core revision, and each linked archive by SHA-256.

    SOURCE=experiments/leopard2/c5_dyadic_cpp/c5_dyadic.cpp
    SOURCE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
    CORE_SHA=37e852774bce6f9effb1acf1fcde99a758ecfe6e

    for variant in auto scalar ssse3 avx2; do
        cmake -S . -B "build/c5/$variant" -G Ninja \
            -DCMAKE_BUILD_TYPE=Release \
            -DLEO2_BUILD_TESTS=ON \
            -DLEO2_BUILD_BENCHMARKS=OFF \
            -DLEO2_BACKEND_VARIANT="$variant"
        cmake --build "build/c5/$variant" -j8 --target libleopard
        LIBRARY="build/c5/$variant/liblibleopard.a"
        LIBRARY_SHA="$(sha256sum "$LIBRARY" | awk '{print $1}')"
        g++ -std=c++11 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
            -DLEO2_ENABLE_TEST_HOOKS=1 \
            -DLEO2_C5_SOURCE_SHA256="\"$SOURCE_SHA\"" \
            -DLEO2_C5_CORE_GIT_SHA="\"$CORE_SHA\"" \
            -DLEO2_C5_LIBRARY_SHA256="\"$LIBRARY_SHA\"" \
            -DLEO2_C5_SANITIZER_MODE='"none"' \
            -I. "$SOURCE" "$LIBRARY" -fopenmp -pthread \
            -o "build/c5/c5-$variant"
        OMP_NUM_THREADS=1 "build/c5/c5-$variant" \
            --mode backend --backend-label "$variant" \
            --output "experiments/leopard2/c5_dyadic_cpp/results/$variant.json"
    done

Run the authoritative timing alone on an allowed physical CPU, changing 15 if
the affinity set differs:

    OMP_NUM_THREADS=1 taskset -c 15 build/c5/c5-avx2 \
        --mode all --backend-label avx2 \
        --output experiments/leopard2/c5_dyadic_cpp/results/benchmark.json

Build and run ASan+UBSan with matching instrumentation on both the library and
standalone source:

    cmake -S . -B build/c5/asan -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_C_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
        -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
        -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=OFF \
        -DLEO2_BACKEND_VARIANT=auto
    cmake --build build/c5/asan -j8 --target libleopard

Use the same compile command as above with the ASan archive, `-O1 -g`,
`-fsanitize=address,undefined -fno-omit-frame-pointer
-fno-sanitize-recover=all`, and `LEO2_C5_SANITIZER_MODE="asan-ubsan"`; then:

    OMP_NUM_THREADS=1 ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      build/c5/c5-asan-ubsan --mode correctness \
        --backend-label asan-ubsan \
        --output experiments/leopard2/c5_dyadic_cpp/results/asan-ubsan.json

Validate the merger itself and regenerate the checkpoint:

    PYTHONDONTWRITEBYTECODE=1 python3 -X dev \
        experiments/leopard2/c5_dyadic_cpp/test_checkpoint.py -v

    PYTHONDONTWRITEBYTECODE=1 python3 -X dev \
        experiments/leopard2/c5_dyadic_cpp/checkpoint.py \
        --backend auto=experiments/leopard2/c5_dyadic_cpp/results/auto.json \
        --backend scalar=experiments/leopard2/c5_dyadic_cpp/results/scalar.json \
        --backend ssse3=experiments/leopard2/c5_dyadic_cpp/results/ssse3.json \
        --backend avx2=experiments/leopard2/c5_dyadic_cpp/results/avx2.json \
        --library auto=build/c5/auto/liblibleopard.a \
        --library scalar=build/c5/scalar/liblibleopard.a \
        --library ssse3=build/c5/ssse3/liblibleopard.a \
        --library avx2=build/c5/avx2/liblibleopard.a \
        --sanitizer experiments/leopard2/c5_dyadic_cpp/results/asan-ubsan.json \
        --sanitizer-library build/c5/asan/liblibleopard.a \
        --benchmark experiments/leopard2/c5_dyadic_cpp/results/benchmark.json \
        --source "$SOURCE" --repository . \
        --output experiments/leopard2/c5_dyadic_cpp/results/checkpoint.json

Retained hashes are listed in
`experiments/leopard2/c5_dyadic_cpp/SHA256SUMS`.  Mathematical lineage and the
wire-identity distinction follow R15 and the derivations in
`docs/leopard2_math_and_sources.md`; no third-party implementation code was
copied.
