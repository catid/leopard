# Leopard2 bounded direct systematic encoding

Status: implemented as a bounded, allocation-free execution path. Production
AUTO dispatch promotes only the conservative low-profile region established by
the crossover evidence below; every other cell retains the transform encoder,
which also remains the wire oracle.

## Scope

The direct implementation is bounded to public counts satisfying

    1 <= K <= 16
    1 <= R <= 16.

The bound is an implementation and measurement boundary, not a mathematical
limit. Test-enabled builds prepare rows for every bounded high/low and GF8/GF16
combination so differential coverage can force either encoder. A production
build allocates the table only for a codec that can enter the AUTO region:
Low V1 with scalar `K >= 3`, or qualified SSSE3/AVX2/AVX512 `K >= 2`. Thus legacy-high,
`K=1`, scalar `K=2`, and unmeasured-backend codecs pay no encoder-table setup or
footprint. A codec outside the bound retains the normal transform encoder
without a partial direct table or a changed wire profile.

For a requested parity coordinate `q`, direct execution evaluates the
systematic generator row

    parity(q) = sum_(0 <= i < K) G(q,i) message(i).

The immutable codec stores all `R*K` nonzero coefficients in row-major order and
in the selected Leopard field's logarithmic representation. Coefficient one is
represented by log zero and specializes to copy for the first term or XOR for a
later term. Other terms use the existing fixed multiply or multiply-add byte
kernels. The result is therefore the same codeword as interpolation followed by
evaluation; this path does not define a new code.

## Full-parent barycentric formula

Coordinates are elements in Leopard's declared Cantor representation. Field
subtraction is addition, and the stored coordinate difference between points
`a` and `b` is `a xor b`.

Let `S` contain every systematic coordinate of the parent code, including the
coordinates shortened to known zero. For each public systematic coordinate
`x_i`, codec setup computes

    w_i = inverse(product over s in S, s != x_i of (x_i xor s)).

For a transmitted parity coordinate `q`, define

    Z_S(q) = product over s in S of (q xor s),

    G(q,i) = Z_S(q) * inverse(q xor x_i) * w_i.

This is the Lagrange cardinal polynomial for `x_i`, evaluated at `q`. Every
factor is nonzero because the systematic and parity coordinate sets are
disjoint. The implementation nevertheless rejects a zero denominator or
coefficient during setup and leaves the transform fallback available.

The profile-specific sets are:

| Profile | Parent systematic set `S` | Public `x_i` | Parity `q_j` |
| --- | --- | --- | --- |
| Legacy high V1 | `[T,N)` | `T+i` | `j` |
| Low V1 | `[0,P)` | `i` | `P+j` |

Here `0 <= i < K` and `0 <= j < R`. In the high profile, the shortened
coordinates are `[T+K,N)`; in the low profile, they are `[K,P)`. They remain in
`S` even though their message values are zero, so they contribute to both
`Z_S(q)` and each derivative denominator represented by `w_i`. In particular,
the high systematic set must not be replaced by an additive-subspace derivative
identity merely because its endpoints are powers of two.

Puncturing only omits parity rows from transmission. High coordinates `[R,T)`
and low coordinates `[P+R,N)` are not members of the public recovery array and
do not alter a transmitted row `G(q_j,i)`. Thus shortening changes the
interpolation problem, while puncturing changes only which evaluations are
visible.

The formula is the same full-parent derivation used by bounded direct repair,
but encoding precomputes every public parity row rather than constructing a
loss-specific square system.

## Requested outputs: count versus prefix

The public recovery pointer array is a mask: a null entry skips that output. The
encoder derives two different quantities from it:

- `Q`, the number of non-null recovery pointers; and
- the transform prefix, one plus the highest requested recovery index.

Direct execution scans the stored rows, skips null outputs, and performs exactly
`K*Q` fixed-coefficient shard operations. It neither computes nor materializes
holes below the highest selected row. The high and low transform encoders retain
their prefix because the truncated transform schedule may need every dependency
up to the highest requested coordinate. Consequently, two masks with the same
`Q` can have the same direct arithmetic cost but different transform costs.

Automatic direct dispatch currently receives the actual `Q`, not the transform
prefix. Any promoted cost rule must be justified for sparse as well as dense
masks; if mask position materially changes the crossover, the dispatcher input
must be extended rather than hiding that dependence in an optimistic threshold.

An all-null recovery mask is a successful no-output call after the ordinary
argument, scratch, padding, and overlap checks. No parity row is executed.

## Setup, execution, and buffer contract

Generator values and logarithms are computed once during codec creation. Setup
may allocate the bounded coefficient vector. Encoding reads it without mutation
and performs no allocation. A constructed codec may therefore be shared across
concurrent encode calls when each call supplies distinct outputs and scratch.

Direct selection occurs only after the existing encode layout and buffer checks.
The public contract is intentionally identical to the transform path:

- `shard_bytes` must be positive and representable by `size_t`;
- GF8 accepts arbitrary positive physical byte counts;
- GF16 requires a positive even physical byte count;
- the explicitly padded-odd GF16 layout is passed as its even physical wire
  size, and every systematic shard must retain its required final zero byte;
- every original pointer is required, while each recovery pointer may be null;
- original inputs may alias one another;
- non-null outputs must be mutually disjoint and disjoint from every input and
  from scratch; and
- scratch must satisfy `leo2_encode_scratch_size()` and
  `leo2_scratch_alignment()` even when direct execution is forced.

The final requirement keeps path selection internal and deterministic: callers
do not need a different allocation strategy for a crossover decision. Direct
execution uses only the range-validation portion of that scratch and writes each
requested output in place; it does not stage transform shard slots.

The byte kernels receive the exact physical `shard_bytes`, not the 64-byte
transform rounding. Unaligned shard buffers and GF8 vector tails are handled by
the existing fixed-multiplier kernels. GF16 complete tiles retain the legacy
ALTMAP representation and even compact tails retain its defined low-byte/high-
byte layout. No internal basis or representation conversion changes stored
parity bytes.

The first source term initializes an output by copy or multiplication. Remaining
terms XOR or multiply-add into it. Validation is complete before this first
write, so an invalid alias, pad byte, size, or scratch argument cannot leave a
partially encoded output.

## Diagnostic hooks and production dispatch

`Leopard2Direct.h` declares private hooks only when
`LEO2_ENABLE_TEST_HOOKS` is defined:

- set a codec to AUTO, forced direct, or forced transform mode;
- query whether the codec owns a complete direct table; and
- query the path for a physical shard size and requested-output count.

These names are not part of `leopard2.h`, the stable ABI, or a production build.
The mode setter is diagnostic mutation: tests set it before sharing a codec with
worker threads and do not change it during concurrent execution. Forcing direct
on an incapable codec returns `LEO2_UNSUPPORTED`; normal builds have neither the
mode field nor the hook symbols.

The production AUTO predicate is a deterministic built-in threshold. It selects
direct encoding only when all of these conditions hold:

- the profile is Low V1;
- exactly one parity output is requested;
- the immutable codec has a complete bounded table, so `1 <= K,R <= 16`;
- the physical shard size is at least 1,024 bytes and a multiple of the
  64-byte SIMD tile;
- `K >= 3` on the scalar backend; or
- `K >= 2` on the qualified SSSE3 or AVX2 backend.

The high profile, `K=1`, two or more outputs, smaller or ragged shards, and
unmeasured backends retain the transform encoder. In particular, scalar `K=2`
and GF8 remainder-63 tails showed real regressions. The threshold does not
benchmark on the caller's first encode and cannot change the fixed profile,
field, coordinate map, or parity bytes.

## Correctness gates

The direct encoder must be checked by a mathematically independent systematic
generator-matrix or interpolation oracle. The transform encoder is an important
differential comparator but cannot be the only oracle. Acceptance coverage must
include:

- both profiles and both fields for every `1 <= K,R <= 16`;
- every systematic basis message plus deterministic random messages;
- dense, prefix, suffix, singleton, and holey recovery masks, including no
  requested output;
- the `K=1,R>1` and `R=1,K>1` fringes and the accepted 16/fallback 17 bounds;
- byte-for-byte comparison with the transform encoder, and with legacy Leopard
  for every high-profile compatibility claim;
- GF8 arbitrary tails and GF16 even compact tails around SIMD and page
  boundaries;
- unaligned buffers, padded-odd wires, rejected odd native GF16 sizes, scratch
  errors, and every allowed or rejected alias class;
- repeated and concurrent execution of one immutable codec;
- an allocation counter around execution; and
- identical parity under every qualified scalar and SIMD backend.

Re-encoding after recovery must reproduce every requested parity byte. A test
that derives expected values through the same barycentric helper or optimized
field byte kernel is not independent and does not satisfy the oracle gate.

## Crossover and promotion gate

Pinned comparisons force direct and transform execution for the exact same
profile, field, `K`, `R`, shard bytes, output mask, backend, compiler, and input.
They report codec setup separately from execution and include at least dense and
holey masks at equal `Q`. Runs use a single isolated physical core for
cache-sensitive latency, controlled warmups, AB/BA order, several repetitions,
and median/MAD or confidence intervals. Batch and multicore measurements are a
separate scaling question and do not replace the isolated crossover data.

The benchmark also reports the direct row-term, initialization, accumulation,
field-symbol, and modeled source/output byte counts. These are an analytic
streaming-kernel model before cache effects, and unit coefficients may specialize
to copy/XOR. Transform operation counts and optional hardware counters are
explicitly null here and remain tracked by the general benchmark-harness bead
`leopard-79h.16`; no cache, cycle, or memory-bandwidth claim is inferred from
the logical GB/s denominator.

The reproducible runner's `--r all --full` screen program covers the bounded
`K,R <= 16` grid, GF8 sizes from one byte through large shards, legal GF16
sizes from two bytes through large shards, every available backend, and
representative reuse and batch counts. The recorded concurrent discovery used
broad and targeted subsets to locate boundaries, followed by the compact pinned
confirmation described below. The forced-direct benchmark deliberately rejects
`K=17` or `R=17`; the accepted-16/fallback-17 boundary is instead covered by
correctness and AUTO dispatch tests. Measurements account for codec setup, the
public scratch contract, bytes read and written, and the difference between `Q`
and transform prefix.

AUTO may be extended only after all correctness gates pass and a parameter
region shows a statistically credible execution improvement of at least five
percent, with no unexplained regression greater than two percent in neighboring
cells. The current dispatcher retains transform fallback outside the promoted
region and at every measured boundary. More complex mask- or backend-specific
rules require correspondingly stronger evidence and a compact maintainable
table. Negative or inconclusive cells remain documented; they are not converted
into a favorable threshold by averaging unrelated workloads.

The discovery sweep executed 3,744 broad cells and 64,662 targeted cells, with
separate forced-direct and forced-transform invocations and no parity or backend
resolution failure. Those concurrent runs were used only to locate boundaries;
their throughput is explicitly non-authoritative. The targeted result and
analysis SHA-256 values are respectively
`9da8a9823ab0defd4388ee9aedec8b5269137c96f9ef33fd3ab9c3e8716623ed`
and
`4e89e3ae141fa6b6f28a09ca05807b18492f0538ece8be46a058840ba48b8783`.

Authoritative confirmation used the checked-in runner, pinned every benchmark
process to CPU 15, and left sibling CPU 31 idle. Each cell reused one stable
seed for direct and transform, ran three ABBA cycles, 15 timed samples per
process, four warmups, and reuse 64. All 75 jobs passed raw-artifact hash,
parameter, parity, untouched-output, backend, and source-stability checks. All
45 candidate cells exceeded the five-percent promotion gate: the minimum
speedup was 19.42%, the median was 113.85%, and no ABBA round regressed. Here
speedup is `(transform_time / direct_time - 1) * 100`, so it may exceed 100%.
The 30 deliberately excluded neighbors were mixed: 16 regressed, with a worst
speedup of -46.48%. In particular, all six scalar `K=2` cells regressed by
20.78% to 31.10%, GF8 ragged tails included severe losses, and `Q=2` was mixed.
This supports the narrow AUTO predicate while retaining transform fallback.

The authoritative source fingerprint is
`61dc00e7d517e63ac50f2c18413817393e7fe0473cf1993e0bd84ee4e63de1b2`.
The manifest, complete matrix, and deterministic analysis SHA-256 values are
respectively
`c9077ba547a36715638f1e8589aff35049a93d854e9e976b3ff1ff71efea1fad`,
`d48df87800b42f1cdd0d72bd66c838212f36f5571a20855bc0708e503f72f5cb`,
and
`e07ead8f010b8e1ff562417560ece852b6b19670b77e0eb27ac3eec3cb18e139`.

An earlier manually scripted pinned run is retained only as historical evidence:
it changed the seed between ABBA processes and its summary mislabeled a
dispersion field. Its timings are not used to justify AUTO. That audit failure
motivated the runner's stable per-cell seed, semantic raw-JSON revalidation,
and artifact hashing; negative results remain visible in the committed
checkpoint rather than being silently discarded.

A historical, non-authoritative tests-off production ABBA comparison against
commit `4b4aa3a` measured codec setup, rather than comparing two forced modes
that both own the test table. It is retained only as setup-cost context, not as
dispatch authority. Ineligible legacy-high and scalar `K=2` codecs were
unchanged. The eligible scalar Low `K=3,R=7` setup median rose from 0.13 to 0.25
microseconds; `K=16,R=16` rose from 0.42 to 1.765 microseconds. At the latter
cell, one 1-KiB direct encode saved 3.257 microseconds, so the 1.345-microsecond
incremental setup cost was recovered within the first execution. The setup and
execution summary SHA-256 values are
`bfb50391486338870a20d759ef284b82ebfd4b70b4dd64f81f2bf1ae986f085e`
and
`4772307616db23d97a018e8cbd5210df4ac09c01f76f934e6651b51a7cb63dd3`.
This older setup study is retained as a historical estimate, not as execution
authority for the dispatcher; the same-seed run above is the promotion gate.

### Reproducing the crossover program

Configure one Release build per qualified backend. `JOBS` should be the
allowed CPU count, capped at 128. This host exposed 32 logical CPUs, so builds
and parallel tests used at most 32 rather than the target-machine maximum of
128. The historical discovery sweep intentionally used 16 workers, one per
physical core, to reduce SMT interference; the new non-authoritative screen
runner defaults to all 32 allowed logical CPUs unless `--workers` overrides it.

```sh
JOBS="$(nproc)"
if [ "$JOBS" -gt 128 ]; then JOBS=128; fi

for backend in scalar ssse3 avx2; do
    cmake -S . -B "build/direct-encode-$backend" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BACKEND_VARIANT="$backend" \
        -DLEO2_BUILD_TESTS=ON \
        -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF \
        -DLEO2_ENABLE_CUDA=OFF
    cmake --build "build/direct-encode-$backend" \
        --target bench_leopard2_direct_encode -j "$JOBS"
done

python3 tools/leopard2_direct_encode_crossover.py self-test
python3 tools/leopard2_direct_encode_crossover.py screen \
    --source . --build-root build \
    --result-dir results/leopard2/direct-encode-crossover/screen \
    --backends scalar,ssse3,avx2 --r all --full --workers "$JOBS"
python3 tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir results/leopard2/direct-encode-crossover/screen
```

Screening is intentionally concurrent and non-authoritative. Before the pinned
phase, stop the screen and every other compute- or memory-intensive job. Choose
one allowed physical CPU, keep its SMT sibling idle, and substitute that CPU for
`ISO_CPU`; the runner records the allowed affinity set and topology and uses one
worker. Do not assume CPU numbers are contiguous.

```sh
ISO_CPU=15
python3 tools/leopard2_direct_encode_crossover.py pinned \
    --source . --build-root build \
    --result-dir results/leopard2/direct-encode-crossover/pinned \
    --backends scalar,ssse3,avx2 --r 1,16 \
    --cpu "$ISO_CPU" --abba-rounds 3 \
    --iterations 15 --warmups 4 --reuse 64
python3 tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir results/leopard2/direct-encode-crossover/pinned
```

This compact pinned boundary program is 75 cells and 900 benchmark process
launches. Adding `--full` expands it to 1,444 cells and 17,328 launches; treat
that as an optional long isolated campaign, not the ordinary reproduction step.

Each result directory is configuration-specific and resumable. Reusing a
directory with changed source, flags, executables, grid settings, or machine
identity is rejected instead of mixing evidence. Passed jobs retain hashed raw
benchmark JSON and stdout/stderr; analysis revalidates those hashes and semantic
parameters. The compact, committed checkpoint at
`experiments/leopard2/direct_encode/results/checkpoint.json` preserves the
promotion evidence needed to resume on another machine, while the much larger
raw result trees remain generated artifacts.

## Legacy-high GF8/AVX2 full-output experiment

The legacy-high full-output candidate remains disabled by default and does not
change the ordinary production dispatcher. It is compiled with:

```sh
-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=ON
```

`LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO` defaults to `ON` only inside that
default-off experiment. Setting it to `OFF` retains the same generator tables
and AVX2 object code while forcing production AUTO back to the transform
encoder. This supplies a same-layout timing control without changing the
stable API or wire profile.

The candidate is limited to legacy-high GF8, an explicit AVX2 backend,
`5 <= K <= 16`, `5 <= R <= 8`, and all `R` output pointers present. For
shards through 64 bytes, an AVX2 output-group kernel keeps one through five
destinations in registers and evaluates exact-width 32-, 16-, 8-, 4-, 2-, and
1-byte pieces. A three-byte remainder deliberately uses scalar products because
the XMM form lost in the whole-call screen. The `K=5,R=5` kernel fixes both
counts at compile time; other shapes use groups of four outputs. Larger
eligible shards use the existing source-major multi-output multiply-add
backend. Partial-output calls retain the mature row-major or transform paths.

AUTO contains only the measured lengths that cleared the same-source
five-percent discovery gate. `K=5,R=5` has a separate broader region; the
other 41 shapes use a compact exact-length switch through 65 bytes. The
`K=16,R=8` three- and 64-byte cells remain transform fallbacks. This irregular
table is one reason the complete experiment remains default-off pending a
stronger isolated campaign and a solution to the residual exact-main gap.

Correctness does not depend on timing evidence. The Release direct-generator
suite passed 1,024 profiles, 8,704 basis messages, 48,749,869 parity symbols,
8,192 output masks, unaligned buffers, concurrency, allocation auditing, and
AUTO boundaries. A production-only matrix then checked all 42 eligible `K,R`
shapes at 16 tail and SIMD boundaries: AUTO-on executed 462 direct and 210
transform cells, while the same-layout AUTO-off build executed all 672 through
the transform path. Both matched the independent systematic generator oracle.
The same 672-cell AUTO-on matrix passed Clang 18 ASan+UBSan at about 36 MiB
RSS. The much broader sanitizer executable is not required for this candidate;
it reached the session's 384 MiB safety ceiling and was replaced by this bounded
gate rather than rerun with an unsafe memory allowance.

The final frozen production executables had these SHA-256 values:

- candidate: `1c61bbd0c5c39081769c31de3d9579850fbbfc758ae51be031068c84f6be40a2`
- same-layout AUTO-off control:
  `3db952c1149bf49894207ba561d3979e12835edff580836aca34a78ca0f05d71`
- exact Leopard main:
  `be4be156bf873d02ab6b11c95fcc805070c947501f6567a37181450ea7008d9e`

Candidate and control used byte-identical AVX2 backend objects. The compact
CPU-4 run held the canonical benchmark lock, pinned every process, verified
all three executable hashes before and after, and required identical original,
parity, and recovered-output digests. It is directional rather than promotion
authority while the coordinator-exclusion issue remains open. Ratios greater
than one mean the candidate is faster:

| K,R | bytes | exact main / candidate | AUTO-off / candidate |
| --- | ---: | ---: | ---: |
| 5,5 | 64 | 0.897 | 1.819 |
| 5,5 | 256 | 0.832 | 1.244 |
| 5,5 | 1,024 | 1.045 | 1.204 |
| 5,5 | 4,096 | 1.288 | 1.060 |
| 8,8 | 64 | 0.513 | 1.233 |
| 8,8 | 256 | 0.661 | 1.020 |
| 8,8 | 1,024 | 0.919 | 0.996 |
| 8,8 | 4,096 | 1.115 | 1.020 |
| 12,8 | 64 | 0.593 | 1.291 |
| 16,8 | 64 | 0.545 | 1.047 |

Thus the arithmetic kernel substantially improves Leopard2's own tiny
full-output path, but fixed safe-API validation/setup and `K*R` direct
arithmetic still leave 64-byte `K>=8` cells behind Leopard main. It is not
promoted merely because it wins against Leopard2's transform fallback.

Rejected variants are retained as negative results rather than production
code: a whole-T=8 callback, two- and three-output grouping, an XMM three-byte
tail, forced source-loop unrolling, and a generated coefficient-major
`K=5,R=5` circuit. Each either failed the five-percent whole-call gate or
regressed through extra source loads, register pressure, or code size.
