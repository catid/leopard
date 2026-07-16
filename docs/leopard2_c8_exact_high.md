# Leopard2 C8 exact high-rate parity solve

## Outcome

Experiment C8 is complete as a default-off research checkpoint.  It promotes
no public profile, production kernel, or dispatcher entry.

The result is **inconclusive with no production promotion**:

- precomputed dense exact-`R` execution is substantially faster than the
  padded legacy GF16 encoder for reused tiny-shard cells at `K=500, R=2..4`;
- the same candidate loses every measured GF8 cell, loses again at GF16
  `R=5`, and loses badly at `R=17` and `R=33`;
- its generator-table setup is orders of magnitude slower than legacy codec
  setup and must be amortized over many stripes;
- the coordinate map is a separately versioned candidate, not legacy wire
  format, and it has no production serialization or C9 exact-profile decoder;
- the partial/transposed LCH and suffix-set Tang--Han adaptations do not yet
  have proved executable schedules.

The narrow GF16 region is retained for C9/C10.  A future deterministic
dispatcher can revisit it only after the exact profile has a serialized code
identity and a matching decoder.  The broad dense encoder, factorized solve,
Newton execution, and dyadic Schur execution are rejected as production-wide
candidates.

All implementation and evidence is under:

    experiments/leopard2/non_power_of_two/c8/

The root `CMakeLists.txt`, installed targets, `leopard2.h`, stable profile
constructors, and production AUTO dispatcher do not reference this experiment.
`LEO2_PROFILE_EXACT_EXPERIMENTAL_V1` therefore continues to return
`LEO2_UNSUPPORTED` from the stable codec constructor.

## Candidate wire profile

For distinct field points identified by Leopard's public Cantor coordinates,
the candidate uses:

    parity evaluations:     omega_0 .. omega_(R-1)
    systematic evaluations: omega_R .. omega_(R+K-1)
    polynomial bound:       degree less than K
    transmitted role order: [ parity ][ systematic ]

The deterministic research name is
`exact_high_prefix_v1_candidate`, coordinate-map version 1.  This name is not
a stable serialized identifier.  GF8 is possible when `K+R <= 256`; otherwise
the experiment uses legacy-representation GF16 when `K+R <= 65,536`.

Let the message values be `m_i = f(omega_(R+i))`.  The parity values are

    p_j = f(omega_j), 0 <= j < R,

where `f` is the unique polynomial of degree less than `K` interpolating the
message evaluations.  The code is systematic and MDS because evaluation of a
degree-less-than-`K` polynomial at any `K` distinct points is nonsingular.
The GF(2^4) exhaustive gate below independently verifies that argument for every
small geometry and every transmitted `K`-coordinate subset.

This layout deliberately retains the legacy high profile's parity-before-data
role convention.  It does not silently reuse the low-profile or C6 exact-prefix
identity, whose application-visible role order and evaluation assignment differ.

## Exact `R` parity constraints

The scalar algebra derives the parity solve independently of the direct
Lagrange generator.

Let `n=K+R`, `x_i=omega_i`, and

    u_i = 1 / product_(ell != i) (x_i + x_ell).

For `0 <= a < R`, a generalized-RS dual check is

    sum_(i=0)^(n-1) u_i x_i^a c_i = 0.

Partitioning parity columns from message columns gives

    A p + B m = 0.

The field has characteristic two, so subtraction and addition coincide, and

    p = inverse(A) B m.

`algebra.py` constructs `A` and `B`, inverts the exact `R`-by-`R` matrix, and
requires `inverse(A) B` to equal generator rows built independently from
Lagrange cardinal polynomials.  The C++ executor precomputes the equivalent
Lagrange generator in `O(K^2 + KR)` field operations rather than paying the
less attractive generic dual-matrix setup.  The exact-`R` solve and the faster
setup must agree coefficient for coefficient before byte execution is tested.

The executable stores `K*R` fixed multipliers in the immutable experimental
plan.  Encoding performs no allocation and needs no shard-data scratch:

    for each parity j:
        parity[j] = sum_i generator[j][i] * message[i]

The production FF8/FF16 fixed-multiplier backends perform byte execution.  A
separate carryless-multiply/reduce field implementation with explicit legacy
Cantor-to-polynomial maps derives every expected coefficient and output symbol.
It calls neither the production logarithm tables nor the LCH transforms.

## Other C8 candidates

### Newton interpolation

The scalar experiment computes ordinary divided differences on the exact
message nodes `omega_R .. omega_(R+K-1)` and evaluates the Newton form at all
parity points.  It agrees with both the dual solve and direct cardinal
evaluation.  Its execution work is quadratic in `K` before the `KR` parity
evaluation and is rejected as a hot path.  It remains an independent oracle.

### Dyadic Schur complement

For non-power-of-two `R`, let `b` be the largest power of two below `R` and
partition the parity matrix as

    A = [ A11 A12 ]
        [ A21 A22 ].

The experiment constructs

    S = A22 + A21 inverse(A11) A12

and solves through `inverse(A11)` and `inverse(S)`.  It agrees with the full
solve for every exercised vector.  Its factorized execution adds `R^2`-scale
work after forming the message constraints, and its all-GF8 symbolic gain never
reaches one.  Precombining it collapses back to the same dense `R`-by-`K`
generator, so no distinct Schur executor is promoted.

### Partial/transposed LCH transforms

The all-count matrix includes an explicitly optimistic transposed-LCH lower
bound proportional to `(K+R) ceil(log2 R)`.  It is labeled non-executable in
every row.  No irregular exact-`R` transposed schedule has yet been derived
that preserves this coordinate profile and beats the regular `T` kernels after
memory traffic.  The lower bound prioritizes possible follow-up work; it is not
correctness or performance evidence and cannot promote a kernel.

### Tang--Han epsilon transform

Tang--Han Appendix-A Algorithm 8, already implemented and tested in C3b,
interpolates a prefix of an additive-coset evaluation order.  C8's known
message set is the suffix `R .. R+K-1` of the public prefix.  Consecutive
integer indices are not generally `shift XOR (0..K-1)`, so mechanically applying
the prefix algorithm would change the coordinate map.  C8 records the published
prefix method only as a differently mapped control/lower model.  A valid suffix
adaptation needs a new derivation; no numerical relabeling is presented as a
Tang--Han implementation.

### Regular padded `T` encoder

The baseline is the production legacy-high encoder with

    T = ceil_pow2(R)
    N = ceil_pow2(K + T).

Its regular block IFFTs and final `T`-point FFT remain very effective because
`T < 2R`, SIMD kernels are mature, and codec setup is small.  This baseline
uses a different code outside the identity condition below; timing comparisons
therefore compare versioned code choices, not interchangeable kernels.

## Legacy wire-identity gate

Legacy high interpolation fixes all systematic parent coordinates
`T .. N-1`, including the shortened zero tail.  C8 exact interpolation fixes
only `R .. R+K-1`.  The two definitions coincide when and only when the claimed
legacy parent has neither parity padding nor shortened systematic coordinates:

    R = T and K + R = N.

Equivalently, both `R` and `K+R` are powers of two.  In those full-parent cases
the coefficient matrices and encoded bytes must be identical.  Everywhere
else C8 claims a distinct wire profile.

The retained gates include:

- ten full-parent GF(2^4) geometries with exact coefficient identity;
- 75 GF(2^4) geometries with explicit coefficient-difference witnesses (the
  remaining small geometries have a legacy parent larger than the test field);
- 66 executable full-parent byte comparisons across the tested sizes;
- eight executable nonidentity geometries, each with a byte-difference witness;
- GF8/GF16 boundary and random vectors that enforce the same geometry rule.

Correct recovery or equal dimensions alone is never treated as legacy
compatibility.

## Correctness and safety evidence

The deterministic scalar checkpoint passed:

| Gate | Count |
| --- | ---: |
| GF(2^4) public `(K,R)` geometries, `K+R <= 16` | 120 |
| GF(2^4) transmitted `K`-coordinate reconstruction subsets | 131,038 |
| GF(2^4) solver vectors | 57,670 |
| GF(2^4) parity symbols compared | 403,993 |
| GF8/GF16 differential geometries | 19 |
| GF8 all-count model pairs, `K,R>0`, `K+R<=256` | 32,640 |
| Model method rows | 228,480 |

The standalone C++ executable passed identically under auto, forced scalar,
forced SSSE3, and forced AVX2 libraries:

| Gate | Count/result per backend |
| --- | ---: |
| Exact generator coefficients | 78,362 |
| Independently evaluated encoded bytes | 99,420 |
| Correctness digest | `0x4bd56b387de0757e` |
| Hot-path allocations | 0 |

GF8 sizes are `1,7,31,64,65,257` bytes.  GF16 sizes are
`2,6,62,64,66,130` bytes and exercise compact and complete ALTMAP tiles.
ASan+UBSan, including leak detection, passed the same complete C++ projection.
GCC 13.3 compiled the source with `-Wall -Wextra -Wpedantic -Werror`.
`clang++` was not installed on this host, so no Clang result is claimed.

This host exposes 32 logical CPUs, not the mission's nominal 128, and the
coordinator's allowed set contained 30 CPUs while evidence was generated.  The
all-GF8 model used a 30-process pool; the 19-case differential stage naturally
had only 19 independent jobs.  The four production-backend builds/runs were
launched concurrently, and the final CTest gate used all 30 allowed CPUs.
ASan/UBSan correctness and the authoritative timing intentionally used one
thread; sanitizer determinism and single-core benchmark isolation require it.

`test_checkpoint.py` contains nine fail-closed merger tests covering incomplete
matrices, raw-sample and derived-timing forgery, memory-accounting forgery,
backend/source/digest changes, sibling activity, duplicate bindings, and
content hashes.

## All-GF8 operation model

The deterministic architecture-neutral score is

    3*multiplications + XORs + loads + 2*stores + 4*irregular_operations.

It is not a cycle predictor.  Across all 32,640 valid GF8 public pairs:

| Modeled winner | Cells |
| --- | ---: |
| Production padded `T` LCH | 29,359 |
| Precomputed dense exact generator | 3,281 |

The dense exact model exceeds a 1.10 padded-over-exact score in 2,980 cells but
has a global mean gain of 0.616.  Factorized exact solve and dyadic Schur never
reach 1.10.  The transposed-LCH and Tang--Han rows are explicitly non-executable
lower/control models and are not included as promotable winners.

## Pinned AVX2 benchmark

The authoritative run used GCC 13.3 on an AMD Ryzen 9 9950X3D, Linux
6.8.0-134, one socket, one NUMA node.  The process was pinned to physical core
CPU 15 with `OMP_NUM_THREADS=1` and `OMP_DYNAMIC=FALSE`.  Sibling CPU 31
accumulated zero non-idle jiffies during the retained 0.412-second run.  The
governor was `powersave` with energy-performance preference `performance`.
Unprivileged hardware counters were unavailable (`perf_event_paranoid=4`).

Each of 27 cells uses two warmups, eleven samples, and alternating candidate/
baseline order.  All eleven execution samples and all seven setup samples are
retained in the machine-readable artifact and independently re-summarized by
the fail-closed checkpoint.  Setup and execution are separate.  The
conservative execution gain is

    (padded_median - padded_MAD) /
    (exact_median + exact_MAD) - 1.

The conservative amortized gain additionally spreads setup over each cell's
declared `batch*reuse` executions.

Representative cells are:

| Field, K/R | Bytes; batch/reuse | Exact us | Padded us | Execution gain | Amortized gain |
| --- | ---: | ---: | ---: | ---: | ---: |
| GF8, 31/3 | 1 KiB; 8/4 | 1.547 | 1.298 | -16.27% | -19.62% |
| GF8, 120/7 | 1 KiB; 8/4 | 12.891 | 5.995 | -54.85% | -58.40% |
| GF8, 224/31 | 1 KiB; 1/1 | 99.981 | 18.630 | -81.46% | -92.59% |
| GF16, 500/2 | 64 B; 8/4 | 5.279 | 77.080 | +1345.77% | +120.87% |
| GF16, 500/3 | 64 B; 8/4 | 7.962 | 54.145 | +573.67% | +31.10% |
| GF16, 500/4 | 64 B; 8/4 | 10.569 | 57.724 | +443.75% | +46.78% |
| GF16, 500/5 | 64 B; 8/4 | 13.890 | 39.502 | +182.32% | -16.46% |
| GF16, 1000/17 | 1 KiB; 1/1 | 455.545 | 117.201 | -74.47% | -97.52% |
| GF16, 1000/33 | 64 B; 1/1 | 226.422 | 31.170 | -86.51% | -98.74% |

All 16 GF8 execution cells regress, with conservative gains from -88.92% to
-16.27%.  Seven GF16 execution cells exceed +10%, but only the three 64-byte
`R=2..4` cells retain +10% after their declared setup amortization.  The same
`K=500` matrix includes `R=2,3,4,5` neighbors and both 64-byte reused and
1-KiB one-shot cells, so the `R=5` and one-shot setup crossover is explicit.

Exact setup medians span 2.09--156.10 microseconds in GF8 and
917.75--5,230.10 microseconds in GF16.  Padded setup spans 0.09--0.37 and
0.93--38.07 microseconds respectively.  The retained exact coefficient tables
span 62--6,944 bytes in GF8 and 2,000--66,000 bytes in GF16.

Execution scratch is zero for the dense exact candidate.  Public padded encoder
scratch spans 1,088--1,051,712 bytes in GF8 and 12,352--90,368 bytes in GF16
for this matrix.  That scratch advantage does not offset dense direct memory
traffic outside tiny `R`: the exact logical model reads and accumulates every
one of the `K*R` terms.

## Disposition by method

| Method | Disposition |
| --- | --- |
| Precomputed dense exact-`R` generator | Retain only as inconclusive GF16 `R=2..4` candidate; reject universal use. |
| Factorized `R`-by-`R` solve | Reject execution; precombine reduces to the dense generator. |
| Newton interpolation/evaluation | Reject execution; retain independent oracle. |
| Dyadic Schur complement | Reject execution; correct but adds dense joins. |
| Partial/transposed exact-`R` LCH | Inconclusive lower model; no executable schedule or promotion claim. |
| Tang--Han epsilon suffix adaptation | Inconclusive derivation; prefix implementation cannot be relabeled. |
| Production padded `T` LCH | Retain as production default and compatibility profile. |

No method enters AUTO.  No CPU-dependent choice changes parity format.

## C7 integration caveat

C7 was marked in progress in Beads, but no local worktree, branch, commit, or
remote ref containing its exact-low implementation existed at the C8 base
revision.  C8 therefore proceeds independently from the completed C0--C6
evidence.  This does not block C8's parity-first candidate, but before C9/C10
freeze a shared exact profile they must compare C7's eventual coordinate map,
serialized identity, and decoder expectations.  Similar mathematics is not a
license to merge code identifiers.

## Reproduction

Run scalar algebra and the all-GF8 model on every allowed CPU (30 were available
during the retained run because another benchmark reserved one sibling pair):

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/algebra.py self-test
    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/algebra.py analyze \
        --workers "$JOBS" \
        --output-dir experiments/leopard2/non_power_of_two/c8/results
    python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/algebra.py verify \
        --output-dir experiments/leopard2/non_power_of_two/c8/results

Build each normal backend library and the standalone source.  Replace `VARIANT`
with `auto`, `scalar`, `ssse3`, and `avx2`:

    cmake -S . -B "build/c8/$VARIANT" -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF \
      -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
      -DLEO2_BACKEND_VARIANT="$VARIANT"
    cmake --build "build/c8/$VARIANT" -j "$JOBS" --target libleopard

    SOURCE=experiments/leopard2/non_power_of_two/c8/c8_exact_high.cpp
    SOURCE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
    CORE_SHA="$(git rev-parse dfa69baab6f056a6b09f4548bc196ac39797294a)"
    LIBRARY="build/c8/$VARIANT/liblibleopard.a"
    LIBRARY_SHA="$(sha256sum "$LIBRARY" | awk '{print $1}')"
    g++ -std=c++11 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
      -DLEO2_C8_SOURCE_SHA256=\"$SOURCE_SHA\" \
      -DLEO2_C8_CORE_GIT_SHA=\"$CORE_SHA\" \
      -DLEO2_C8_LIBRARY_SHA256=\"$LIBRARY_SHA\" \
      -DLEO2_C8_SANITIZER_MODE=\"none\" \
      -I. "$SOURCE" "$LIBRARY" -fopenmp -pthread \
      -o "build/c8/c8-$VARIANT"
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE "build/c8/c8-$VARIANT" \
      --mode correctness --backend-label "$VARIANT" \
      --output \
        "experiments/leopard2/non_power_of_two/c8/results/$VARIANT.json"

Build both the library and experiment with ASan+UBSan; define
`LEO2_C8_DISABLE_GLOBAL_NEW_TRACKING=1` for the standalone sanitizer source and
use `LEO2_C8_SANITIZER_MODE="asan-ubsan"`.  Run with:

    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
    UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      build/c8/c8-asan-ubsan --mode correctness \
        --backend-label asan-ubsan \
        --output \
          experiments/leopard2/non_power_of_two/c8/results/asan-ubsan.json

Run timing only in a quiet window, changing the physical/sibling pair after
checking `lscpu -e=CPU,CORE,SOCKET,NODE`:

    python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/run_pinned.py \
      --cpu 15 --sibling 31 \
      --output experiments/leopard2/non_power_of_two/c8/results/isolation.json \
      -- build/c8/c8-avx2 --mode benchmark --backend-label avx2 \
      --output experiments/leopard2/non_power_of_two/c8/results/benchmark.json

Run merger tests and regenerate the checkpoint with the retained backend,
sanitizer, benchmark, isolation, source, and rebuilt-library paths:

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/test_checkpoint.py -v
    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c8/checkpoint.py \
        --algebra-dir experiments/leopard2/non_power_of_two/c8/results \
        --algebra-source experiments/leopard2/non_power_of_two/c8/algebra.py \
        --backend auto=experiments/leopard2/non_power_of_two/c8/results/auto.json \
        --backend scalar=experiments/leopard2/non_power_of_two/c8/results/scalar.json \
        --backend ssse3=experiments/leopard2/non_power_of_two/c8/results/ssse3.json \
        --backend avx2=experiments/leopard2/non_power_of_two/c8/results/avx2.json \
        --library auto=build/c8/auto/liblibleopard.a \
        --library scalar=build/c8/scalar/liblibleopard.a \
        --library ssse3=build/c8/ssse3/liblibleopard.a \
        --library avx2=build/c8/avx2/liblibleopard.a \
        --sanitizer experiments/leopard2/non_power_of_two/c8/results/asan-ubsan.json \
        --sanitizer-library build/c8/asan/liblibleopard.a \
        --benchmark experiments/leopard2/non_power_of_two/c8/results/benchmark.json \
        --benchmark-library build/c8/avx2/liblibleopard.a \
        --isolation experiments/leopard2/non_power_of_two/c8/results/isolation.json \
        --source experiments/leopard2/non_power_of_two/c8/c8_exact_high.cpp \
        --repository . \
        --output experiments/leopard2/non_power_of_two/c8/results/checkpoint.json

Finally verify all retained hashes:

    sha256sum -c experiments/leopard2/non_power_of_two/c8/SHA256SUMS

## Sources

The implementation is a clean local derivation.  No research implementation
code was copied.

- Leopard-RS: https://github.com/catid/leopard
- Tang and Han, *New Decoding of Reed-Solomon Codes Based on FFT and Modular
  Approach*: https://arxiv.org/abs/2207.11079 and
  https://arxiv.org/pdf/2207.11079
- Coxon, *Fast Transforms over Finite Fields of Characteristic Two*:
  https://arxiv.org/abs/1807.07785 and
  https://arxiv.org/pdf/1807.07785
- Lin, Chung, and Han, *Novel Polynomial Basis and Its Application to
  Reed-Solomon Erasure Codes*: https://arxiv.org/abs/1404.3458 and
  https://arxiv.org/pdf/1404.3458
- Yu, Lin, Hou, and Li, *Reed-Solomon Coding Algorithms Based on Reed-Muller
  Transform for Any Number of Parities*:
  https://ieeexplore.ieee.org/document/10086680/

The broader source and retrieval log remains in
`docs/leopard2_math_and_sources.md`.
