# Leopard2 C2 parent-preserving truncated transforms

Status: scalar correctness and a standalone C++ production-relevance checkpoint
are complete.  The candidate is **not promoted**: compact scratch is useful,
but the scalar ragged-boundary executor is much slower than the padded fused
SIMD transform.  This does not define an exact-length code, does not change a
coordinate set, and is not enabled by the default build.

The implementation is
`experiments/leopard2/non_power_of_two/c2/truncated_transform.py`.  It accepts a
power-of-two parent, an aligned additive-coset shift, an exact active-input
prefix or mask, and an exact requested-output prefix or mask.  Inputs and
outputs use compact mask order.  Boundary recursion uses compact maps, and the
executor creates a dense vector only for a complete dyadic subtree whose input
and output masks are both full.  It therefore never creates an `N`-entry dense
vector merely to represent an unused parent suffix.

Complete-subtree descriptors call the full scalar butterfly kernel in this
checkpoint.  They are the exact boundaries at which a later C++ implementation
can call Leopard's existing fused scalar/SIMD kernels; this Python program is
not itself an optimized-kernel implementation.

This is a wire-preserving computation of the existing parent transform.  The
mathematical padded zeros and unrequested coordinates still exist in the parent
code.  They are not transmitted or materialized by this scalar executor.  A
successful decode or matching evaluation is not evidence for a new exact-size
wire profile; no such compatibility claim is made here.

## Recurrences and active cosets

For a node of size `2h`, write its LCH coefficients as pairs `(a_i,b_i)`, and
let `m` be Leopard's existing skew factor for the node's active coset.  The
forward butterfly is

    u_i = a_i + m b_i
    v_i = a_i + (m + 1) b_i = b_i + u_i

The left child evaluates `u` at the current shift and the right child evaluates
`v` at `shift + h`.  Because every accepted shift is aligned to the parent,
integer `shift + h` is also Cantor-coordinate XOR by `h`; it is not ordinary
field addition in a hidden representation.

The exact inverse pair is

    b_i = u_i + v_i
    a_i = (m + 1) u_i + m v_i = u_i + m (u_i + v_i)

Plan setup propagates structural nonzero inputs forward and requested outputs
backward through these exact coefficients.  For the inverse, a requested
`a_i` needs `u_i` unless `m + 1` is zero and needs `v_i` unless `m` is zero;
a requested `b_i` needs both.  Thus zero and one skew factors are handled by
the derivation rather than test exceptions.

The scalar skew generator reconstructs the field elements stored logarithmically
in Leopard's `FFTSkew`.  It preserves the legacy GF8 polynomial, Cantor basis,
symbol representation, and coordinate order.  The same recurrence is tested
for every power-of-two active parent in GF(2^4), and through 256 coordinates in
legacy GF8.

## Relationship to Coxon's truncated transforms

R15 treats LCH/Lagrange conversion as a basis conversion and develops pruned,
cache-friendly Algorithms 3 and 4 with mixed intermediate representations.  C2
adopts the essential parent-embedding and recursive-truncation principles:

- retain the enclosing dyadic parent and discard only computations that cannot
  affect requested coordinates;
- decompose prefixes into complete dyadic regions plus a recursive ragged
  boundary; and
- keep shifted-coset normalization explicit at each child.

This checkpoint is not a line-by-line transcription of Coxon's mixed
Lagrange/LCH in-place algorithms, and it does not claim their auxiliary-space or
operation bounds.  Its purpose is narrower: establish an independently verified
compact recurrence in Leopard's precise coordinates before studying the more
specialized basis-conversion alternatives in C3.

Primary reference:

- Nicholas Coxon, *Fast Transforms over Finite Fields of Characteristic Two*,
  https://arxiv.org/abs/1807.07785 and
  https://arxiv.org/pdf/1807.07785, especially Sections 3.3 and 3.4.

The legacy coordinate and active-parent derivations remain documented in
`docs/leopard2_math_and_sources.md`; C2 changes none of those definitions.

## Independent correctness oracles

The primary oracle does not use the candidate skew graph.  It directly builds
each normalized LCH basis polynomial from subspace polynomials, evaluates that
basis at every shifted parent point, and forms the full evaluation matrix.  The
inverse oracle inverts that independently constructed matrix with GF Gaussian
elimination.  Equality on every permitted unit input proves equality of the two
GF-linear maps for all field-valued inputs.

A secondary padded oracle runs the complete legacy butterfly recurrence over an
explicit `N`-symbol vector.  It catches coordinate/schedule mistakes separately
from the direct polynomial construction.  Candidate plans also execute with an
allocation audit: any dense allocation whose local input or output mask is not
full fails the test.

The retained deterministic result is
`experiments/leopard2/non_power_of_two/c2/results/self_test.json`.  The final
SHA-256 values are recorded below after the reproduction commands.

The deterministic campaign covers:

| Check | Count |
| --- | ---: |
| Every prefix geometry, aligned coset, and direction in GF(2^4) | 1,374 plans |
| Every sparse input/output mask through N=8 in every aligned GF(2^4) coset, both directions | 264,576 plans |
| GF8 dyadic boundary-prefix plans with shifted cosets through N=256 | 1,674 plans |
| GF8 irregular sparse-mask plans with shifted cosets through N=256 | 378 plans |
| Unit input vectors compared with the direct matrix | 1,107,672 |
| Direct requested-symbol comparisons | 8,315,600 |
| Secondary padded vectors | 4,814 |
| Secondary padded requested-symbol comparisons | 26,070 |
| Zero vectors | 268,002 |
| Compact API rejection checks | 3 |
| Partial-subtree dense materializations | 0 |

“Exhaustive GF(2^4)” here has a precise scope.  Every prefix geometry is covered
for every legal parent through N=16 and every aligned coset.  Every pair of
sparse masks is covered through N=8.  Enumerating all pairs of 16-bit sparse
masks would require more than four billion plan geometries, so N=16 sparse
masks are covered by all prefixes while broader irregular behavior is covered
in the GF8 sweep.  Unit-vector equality is exhaustive for each tested linear
map; it is not a claim that all `16^N` messages were enumerated.

## Operation and memory accounting

Counts below are exact idealized transform-kernel counts emitted by the
instrumented scalar schedule and exclude plan setup.  A fixed multiplication
excludes constants zero and one.  Loads and stores are only logical operand and
result accesses in butterflies and compact boundary combinations.  Compact
input/output packing, dense gather/scatter, leaf accesses, dictionary merges and
lookups, final output scatter, and Python object overhead are not counted.  The
padded oracle uses a complete radix-2 parent.  `Dense` is the largest local
complete-subtree vector, not a claim about final C++ scratch allocation or total
memory traffic.

| Direction | N | Active / requested | Fixed mul candidate / padded | XOR candidate / padded | Loads candidate / padded | Stores candidate / padded | Complete blocks | Largest dense candidate / padded |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Forward | 16 | 9 / 7 | 5 / 17 | 16 / 49 | 32 / 64 | 31 / 64 | 2 | 4 / 16 |
| Forward | 64 | 33 / 17 | 32 / 129 | 80 / 321 | 158 / 384 | 143 / 384 | 1 | 16 / 64 |
| Forward | 256 | 129 / 65 | 192 / 769 | 448 / 1,793 | 766 / 2,048 | 703 / 2,048 | 1 | 64 / 256 |
| Inverse | 16 | 9 / 7 | 5 / 17 | 16 / 49 | 30 / 64 | 29 / 64 | 3 | 4 / 16 |
| Inverse | 64 | 33 / 17 | 32 / 129 | 80 / 321 | 128 / 384 | 113 / 384 | 1 | 16 / 64 |
| Inverse | 256 | 129 / 65 | 192 / 769 | 448 / 1,793 | 640 / 2,048 | 577 / 2,048 | 1 | 64 / 256 |

In the representative prefix cases, candidate logical loads fall to 31-50% of
the padded count, stores to 28-48%, and the largest dense subtree to one quarter
of the parent.  These are algebraic/schedule results, not runtime measurements.
The Python implementation uses dictionaries at ragged boundaries and does not
establish the peak scratch, cache traffic, SIMD regularity, or throughput of a
production implementation.

## Standalone C++ production-relevance checkpoint

`experiments/leopard2/c2_truncated_cpp/c2_truncated.cpp` translates the compact
recursion to a reusable, immutable flat C++ schedule without changing the core
library or default CMake graph.  It is compiled manually against a test-enabled
Leopard library.  Complete full/full subtrees call the actual production GF8 or
GF16 transform kernels.  Ragged pairs use an independent scalar polynomial
field implementation and lane-varying 64-byte shards; GF16 reads and writes all
32 low/high ALTMAP pairs independently.

The C++ checkpoint independently reconstructs Leopard's skew factors in the
declared Cantor coordinates and compares all 65,790 GF8/GF16 entries with the
production test hooks before executing a plan.  Each candidate result is then
compared byte-for-byte with a complete, explicitly zero-padded production parent
transform.  This is a parent-wire-identity test, not an exact-profile test.

The bounded campaign contains 22 cases per backend and compares 684,416 output
bytes per backend.  It covers both directions for:

- GF8 prefixes and first/last aligned cosets through N=256;
- GF16 first/last aligned cosets at N=256, 1,024, and 4,096, plus a deep
  last-coset prefix at N=8,192;
- GF16 medium and large cells at N=1,024 and N=4,096; and
- irregular GF16 input/output masks at N=1,024.

Auto, forced scalar, forced SSSE3, and forced AVX2 builds produce the same
standard FNV-1a-64 correctness digest, `0xeebedd5febef7e07`.  The standalone C++ source and
linked test-enabled library also pass ASan plus UBSan.  These values are from
the retained post-freeze run bound to core-matrix fingerprint
`076d16f953d643485dcdd1ab084f887c572ff2558f715380cee7236c511c5d58`.
The raw and merged results named in the reproduction section are authoritative
for this bounded checkpoint.

Each raw result embeds the SHA-256 of this C++ source, the frozen core Git SHA,
the SHA-256 of the exact linked static archive, compiler identity, and sanitizer
mode.  The merger recomputes the source/archive hashes and rejects unbound or
cross-wired evidence.  It also requires the exact 22-case geometry, validates
scratch formulas and numeric ranges, and requires operation/scratch/schedule
parity across all four backends, the sanitizer run, and the benchmark run.

### Bounded schedule and scratch accounting

Plan construction first emits SSA-like compact values and then performs a
deterministic liveness allocation.  A boundary operation has four simultaneously
live 64-byte shard objects (two inputs and two outputs) before reusing dead input
slots.  A complete-subtree operation is
in-place and needs only a pointer vector for its largest full block.  Therefore
the reported logical execution-scratch bound is:

    arena_slots * 64 + 4 * 64 + maximum_complete_block * sizeof(void*)

Caller-owned compact inputs and outputs are excluded on both candidate and
padded sides.  The padded comparison counts `N * 64 + N * sizeof(void*)`.
Serialized schedule bytes use a fixed documented representation (32-byte header,
32-bit maps/slots, 16-byte operation headers); C++ allocator headers are not
included.  Final `std::vector` capacity is reported separately as resident plan
bytes excluding allocator headers.  Ordinary C++ call frames, exception runtime
state, and compiler spills are not part of `execution_scratch_bytes`.  GCC 13.3
with `-O3 -fstack-usage` reports a 384-byte static frame for the inlined
`Executor::execute` path: 256 bytes are the four explicitly counted shards and
128 bytes are other frame/spill space.  It reports 112 bytes for the padded
wrapper.  Applications budgeting total per-thread transient memory must add the
compiler-reported stack bound; the result field describes explicit algorithmic
scratch, not the entire thread stack.

For prefix cases, bounded execution scratch is 47.28-51.39% of the padded
workspace.  The irregular cases retain 77.17-77.43%, and their schedules contain
no complete subtree at all.  The deepest retained GF16 prefix uses 278,912 bytes
instead of 589,824 bytes at N=8,192.  Its inverse serialized schedule is 368,712
bytes, demonstrating that naive per-pair flat schedules can cost more metadata
than they save in scratch.

### Setup versus execution

The retained single-thread run was pinned to one allowed logical CPU with
`OMP_NUM_THREADS=1`.  Inputs, outputs, schedules, and scratch were allocated
before execution timing.  Candidate timing includes compact input copy, ragged
work, complete-subtree kernels, and output scatter.  Padded timing includes zero
fill, active-input copy, full production transform, and requested-output scatter.
Plan setup is reported separately.  Medians and median absolute deviations are
in microseconds; the host was not exclusively reserved, so these are diagnostic
crossovers rather than release benchmark claims.

Pinning is verified from the raw result rather than inferred from the command:
the benchmark JSON records `sched_getaffinity(0) = [15]`,
`OMP_NUM_THREADS = 1`, and `omp_get_max_threads() = 1`, together with hostname,
kernel/architecture, and CPU model.  The merger requires a one-CPU raw affinity,
OpenMP values of one, and identical host/machine identity across release,
sanitizer, and benchmark results.

| Field / N | Direction | Setup median / MAD | Candidate median / MAD | Padded median / MAD | Padded / candidate |
| --- | --- | ---: | ---: | ---: | ---: |
| GF8 / 256 | Forward | 14.650 / 0.270 | 48.171 / 0.331 | 1.500 / 0.010 | 0.031x |
| GF8 / 256 | Inverse | 27.730 / 0.159 | 116.731 / 2.200 | 1.620 / 0.010 | 0.014x |
| GF16 / 1,024 | Forward | 46.180 / 0.310 | 745.499 / 4.890 | 12.640 / 0.069 | 0.017x |
| GF16 / 1,024 | Inverse | 116.641 / 0.910 | 765.569 / 2.999 | 12.630 / 0.090 | 0.016x |
| GF16 / 4,096 | Forward | 171.362 / 3.310 | 3,069.914 / 10.411 | 59.911 / 0.530 | 0.020x |
| GF16 / 4,096 | Inverse | 461.335 / 1.580 | 3,185.055 / 4.970 | 61.001 / 0.430 | 0.019x |

The scalar boundary path is approximately 32.1-72.1 times slower than the padded
production transform in these cells.  This is a negative result for the current
executor, even though scratch improves.  A fused SIMD boundary kernel and a
compressed complete-block-plus-boundary schedule were not implemented, so their
crossover remains inconclusive rather than rejected.

## Disposition and limitations

The bounded C++ prototype validates:

- the compact mask-ordered API contract;
- forward/backward mask propagation using exact butterfly coefficients;
- descriptors for full dyadic subtrees plus ragged boundary pairs; and
- zero/one multiplier specialization with the existing active-coset skew table.

Do not promote either experiment executor or infer a dispatcher threshold.  The
following work remains before production consideration:

- replace scalar boundary multiplication with fused GF8/GF16 SIMD kernels;
- compress per-pair metadata into complete-block descriptors and compact
  boundary runs with an explicit code-size budget;
- add arbitrary byte tails, aliasing, guard-page, sanitizer, fuzz, and immutable
  concurrent-plan tests beyond the standalone 64-byte checkpoint;
- compare against C1 pruning and C3 basis-conversion alternatives; and
- run end-to-end encoder/decoder cells across neighboring geometries before any
  dispatch promotion.

The inverse recurrence is complete for the declared zero-padded parent-transform
semantics and all tested masks.  It is not an exact-K interpolation algorithm and
does not replace the separate C3 inverse/basis-conversion study.

## Reproduction

From the repository root:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c2/truncated_transform.py \
        --output /tmp/c2-self-test.json

    sha256sum \
        experiments/leopard2/non_power_of_two/c2/truncated_transform.py \
        /tmp/c2-self-test.json

Compare `/tmp/c2-self-test.json` byte-for-byte with the retained result.  The
program uses only the Python standard library.  It has no benchmark mode and no
default CMake target imports it.

For the retained source and result:

    359a9a3f14c57de690cc85e1e2cb1a33e6bab734a3ff3712b215ed1cfeeaffa1  truncated_transform.py
    cfd20eaecb97e143fa1b9d3a4bd180b056616942a9dd31a038e27ebea9a99474  self_test.json

The C++ checkpoint is intentionally outside CMake.  Build test-enabled backend
libraries, then compile the standalone source against each archive:

The retained run linked the frozen archives under
`build/resume-backends/{auto,scalar,ssse3,avx2}` and
`build/resume-asan`; their hashes are listed below.  The commands here rebuild
equivalent archives in `build/c2-checkpoint`, embed their newly computed hashes
in fresh raw results, and validate those bindings rather than assuming binary
archive identity across build directories.

    JOBS="$(nproc)"; if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    SOURCE=experiments/leopard2/c2_truncated_cpp/c2_truncated.cpp
    SOURCE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
    CORE_SHA="$(git rev-parse HEAD)"
    CORE_MATRIX_SHA=076d16f953d643485dcdd1ab084f887c572ff2558f715380cee7236c511c5d58
    for variant in auto scalar ssse3 avx2; do
        cmake -S . -B "build/c2-checkpoint/$variant" -G Ninja \
            -DCMAKE_BUILD_TYPE=Release \
            -DLEO2_BUILD_TESTS=ON \
            -DLEO2_BUILD_BENCHMARKS=OFF \
            -DLEO2_BACKEND_VARIANT="$variant"
        cmake --build "build/c2-checkpoint/$variant" \
            -j "$JOBS" --target libleopard
        LIBRARY="build/c2-checkpoint/$variant/liblibleopard.a"
        LIBRARY_SHA="$(sha256sum "$LIBRARY" | awk '{print $1}')"
        g++ -std=c++11 -O3 -DNDEBUG -DLEO2_ENABLE_TEST_HOOKS=1 \
            -DLEO2_C2_SOURCE_SHA256="\"$SOURCE_SHA\"" \
            -DLEO2_C2_CORE_GIT_SHA="\"$CORE_SHA\"" \
            -DLEO2_C2_CORE_MATRIX_SHA256="\"$CORE_MATRIX_SHA\"" \
            -DLEO2_C2_LIBRARY_SHA256="\"$LIBRARY_SHA\"" \
            -DLEO2_C2_SANITIZER_MODE="\"none\"" \
            -I. "$SOURCE" "$LIBRARY" \
            -fopenmp -pthread -o "build/c2-checkpoint/c2-$variant"
    done

Run correctness-only backend checks concurrently; these are deterministic and
are not timing samples:

    mkdir -p experiments/leopard2/c2_truncated_cpp/results
    for variant in auto scalar ssse3 avx2; do
        env OMP_NUM_THREADS=1 "build/c2-checkpoint/c2-$variant" \
            --mode correctness --backend-label "$variant" \
            --output "experiments/leopard2/c2_truncated_cpp/results/$variant.json" &
    done
    wait

Run the diagnostic timing alone on one allowed CPU.  CPU 15 was allowed on the
retained host; select another member of `os.sched_getaffinity(0)` when needed:

    env OMP_NUM_THREADS=1 taskset -c 15 \
        build/c2-checkpoint/c2-auto \
        --mode all --backend-label auto \
        --output experiments/leopard2/c2_truncated_cpp/results/benchmark.json

Configure the sanitizer library and compile the standalone source with matching
instrumentation:

    cmake -S . -B build/c2-checkpoint/asan -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
        -DLEO2_BUILD_TESTS=ON \
        -DLEO2_BUILD_BENCHMARKS=OFF \
        -DLEO2_BACKEND_VARIANT=auto
    cmake --build build/c2-checkpoint/asan \
        -j "$JOBS" --target libleopard
    ASAN_LIBRARY=build/c2-checkpoint/asan/liblibleopard.a
    ASAN_LIBRARY_SHA="$(sha256sum "$ASAN_LIBRARY" | awk '{print $1}')"
    g++ -std=c++11 -O1 -g -DLEO2_ENABLE_TEST_HOOKS=1 -I. \
        -DLEO2_C2_SOURCE_SHA256="\"$SOURCE_SHA\"" \
        -DLEO2_C2_CORE_GIT_SHA="\"$CORE_SHA\"" \
        -DLEO2_C2_CORE_MATRIX_SHA256="\"$CORE_MATRIX_SHA\"" \
        -DLEO2_C2_LIBRARY_SHA256="\"$ASAN_LIBRARY_SHA\"" \
        -DLEO2_C2_SANITIZER_MODE="\"asan-ubsan\"" \
        -fsanitize=address,undefined -fno-omit-frame-pointer \
        -fno-sanitize-recover=all \
        "$SOURCE" "$ASAN_LIBRARY" \
        -fopenmp -pthread -fsanitize=address,undefined \
        -o build/c2-checkpoint/c2-asan-ubsan
    env OMP_NUM_THREADS=1 \
        ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
        UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
        build/c2-checkpoint/c2-asan-ubsan \
        --mode correctness --backend-label asan-ubsan \
        --output experiments/leopard2/c2_truncated_cpp/results/asan-ubsan.json

Validate the fail-closed merger and produce the retained machine summary:

    PYTHONDONTWRITEBYTECODE=1 python3 -X dev \
        experiments/leopard2/c2_truncated_cpp/test_checkpoint.py -v
    PYTHONDONTWRITEBYTECODE=1 python3 -X dev \
        tools/leopard2_c2_checkpoint.py \
        --backend experiments/leopard2/c2_truncated_cpp/results/auto.json \
        --backend experiments/leopard2/c2_truncated_cpp/results/scalar.json \
        --backend experiments/leopard2/c2_truncated_cpp/results/ssse3.json \
        --backend experiments/leopard2/c2_truncated_cpp/results/avx2.json \
        --library auto=build/c2-checkpoint/auto/liblibleopard.a \
        --library scalar=build/c2-checkpoint/scalar/liblibleopard.a \
        --library ssse3=build/c2-checkpoint/ssse3/liblibleopard.a \
        --library avx2=build/c2-checkpoint/avx2/liblibleopard.a \
        --sanitizer experiments/leopard2/c2_truncated_cpp/results/asan-ubsan.json \
        --sanitizer-library build/c2-checkpoint/asan/liblibleopard.a \
        --benchmark experiments/leopard2/c2_truncated_cpp/results/benchmark.json \
        --source "$SOURCE" \
        --core-git-sha "$CORE_SHA" \
        --core-matrix-sha256 "$CORE_MATRIX_SHA" \
        --repository . \
        --output experiments/leopard2/c2_truncated_cpp/results/checkpoint.json

Reproduce the compiler stack accounting separately.  This does not alter the
explicit scratch result:

    g++ -std=c++11 -O3 -DNDEBUG -DLEO2_ENABLE_TEST_HOOKS=1 -I. \
        -fstack-usage -c "$SOURCE" \
        -o build/c2-checkpoint/c2-stack.o
    rg 'Executor::execute|full_padded_execute' \
        build/c2-checkpoint/c2-stack.su

The retained core revision is
`506d3272804efa47e21d6ffbd954b33debcfaf7e`.  Exact source, checker, and result
hashes are:

    932eaf7eb594d53a0aad86ca426785a7258212ec90bba64fe8aa14f005fffe74  c2_truncated.cpp
    2cb5e2de7023b26ebc17155f6980c76659915754f256fe9601a62f82e30fd2e0  test_checkpoint.py
    67a0a26d71feb82d1f953b2c2367f86f1ad3e9054592337fec4385dc61eeaf09  leopard2_c2_checkpoint.py
    c3aa6e1baee766dca717b9b82d85c8f840522da46105dc4a895c5c8cdda9800c  results/auto.json
    4963cec441db221ac081421a312848efff581b880508d90cc889857a924ce317  results/scalar.json
    d09547ce600dff935b37d0ac34b83c6b3b0dd9ded7a2bac23e407cc4b905965e  results/ssse3.json
    2a8ca8b05db7ba1478c471d1bfa5a264cfe75a5119c3e6038b054235fb44ecee  results/avx2.json
    467245e1f2db0358745892193ba0be5ec6b3158c5c9defa20525f324bae656ed  results/asan-ubsan.json
    3023abc7b5242811e2317b235e15b3b8e5dc436fd11a8a24ff5f1b15b2e2d36a  results/benchmark.json
    761466245a079770fb17001467c9efb3186a08c3bca552271db1f3351b916ec6  results/c2-stack.su
    132b2237dffce14bdba8fb621b145935fa8fc811a3849a16ccb61598df45b959  results/checkpoint.json

The linked frozen archive hashes embedded in those raw results are:

    auto    2585cd5059404b7154e11d738545a45f12f24c12a987ed94fa84b8cdf8f3d349
    scalar  f013f8df2dc483bcea0ea416dc4eecd9f0b66d9a78672a2e8fa240b1194fb66d
    ssse3   29a2ce3fe03f1beab1a44eb3644d9f814ef0022477986e85765f5cb9727f4249
    avx2    85cfc5bc2ce21ee0a2874927104e3118403a9d8821ceb9df2746e448aa6122cf
    asan    d0de519fed38d3a0a9e7e91dbd994ab81766bb7ca4b45272cddd7957093981e3
