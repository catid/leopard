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
Authoritative crossover acquisition and replay accept only single-configuration
CMake generators, matching the runner-owned Ninja Release build below.
Multi-config generators remain valid for non-authoritative screening, but their
configuration-qualified artifact layouts are not part of the authoritative
provenance-closure contract.

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

Screening is intentionally concurrent and non-authoritative. Before an
authoritative phase, stop the screen and every other compute- or
memory-intensive job. Choose one allowed physical CPU and its SMT sibling,
substitute them for `ISO_CPU` and `ISO_SIBLING`, and keep the sibling idle. The
runner records the allowed affinity set and topology, holds the canonical
campaign and CPU-pair locks, creates a clean runner-owned Release/explicit-AVX2
build, freezes its executable, and uses one worker. Do not assume CPU numbers
are contiguous. Controlled-build schema v9 records `/usr/bin/cmake` as its
lexical `argv[0]` so CMake can locate its installed modules, while a separate
`/proc/self/fd` executable field proves that the immutable sealed CMake bytes
were executed. Historical v7 replay retains its original procfd-as-`argv[0]`
contract. Across replayable schemas, the audited build umask produces the
controlled executable with exact owner-only mode `0700`; the runner validates
that mode before freezing a read-only execution copy.

```sh
ISO_CPU=15
ISO_SIBLING=47
python3 tools/leopard2_direct_encode_crossover.py historical-avx2 \
    --source . --cpu "$ISO_CPU" --sibling "$ISO_SIBLING"
python3 tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir results/leopard2/direct-encode-crossover/historical-avx2
```

The separate sparse-high discovery campaign is intentionally frozen: its 91
cells include exact K/R/byte/mask boundaries and batches 1, 4, and 16. The
runner enables only `LEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE=ON` and leaves
`LEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE_AUTO=OFF`; changing its backend,
grid, batch, iterations, warmups, reuse, threshold, or worker count is rejected
before topology or artifact handling.

```sh
python3 tools/leopard2_direct_encode_crossover.py sparse-high-avx2 \
    --source . --cpu "$ISO_CPU" --sibling "$ISO_SIBLING"
python3 tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir results/leopard2/direct-encode-crossover/sparse-high-avx2
```

This campaign is a forced-direct versus forced-transform discovery study over
`metrics.encode_execution.median_us_per_batch_call`. Its 95% Student-t
intervals are marginal per-cell intervals with two degrees of freedom, not a
simultaneous guarantee. The frozen campaign itself does not measure production
AUTO and its results do not authorize a dispatcher or production promotion.
Current v5 analysis therefore uses neutral decision-threshold field names for
sparse, historical, and screen results; historical v4 artifacts retain their
original field names for replay. Exit status 2 means that an authoritative
decision threshold was not met. For `sparse-high-avx2`, that is a valid
discovery outcome and is not a statement that a production promotion failed or
was even evaluated.

The separately versioned production-AUTO qualification mode consumes the
ordinary no-hook benchmark rather than either forced path. It creates one clean
Release/Ninja explicit-AVX2 build with tests off, GF8 on, sparse tables on,
sparse AUTO compiled-default off, full-output high direct off, OpenMP on, and
GF16 off.
The runner proves that `bench_leopard2_high_sparse_auto` links the ordinary
`libleopard.a` plus its independent oracle object, rejects every test-hook
object/archive/definition, and freezes the executable and archive together.

```sh
python3 tools/leopard2_direct_encode_crossover.py \
    sparse-high-production-auto-avx2 \
    --source . --cpu "$ISO_CPU" --sibling "$ISO_SIBLING"
python3 tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir \
    results/leopard2/direct-encode-crossover/production-auto-avx2
```

The v10 campaign is fixed at 88 logical cells. Its 74 candidates cover all 36
qualified tuples through one-shot row zero, eight middle/last-row samples, and
batch plus binding APIs at batches 1, 4, and 16 over five anchor shapes. Its 14
route-negative controls cross all seven API/batch lanes at `(2,16,4096,row0)`
for caller-explicit AVX2/thread one and caller-AUTO/thread four. Every omitted
row, API/tuple cross, reuse count, backend, thread count, hardware target, and
out-of-table neighbor is explicitly outside the authorization scope. The
ordered grid has SHA-256
`ff6ab52d98a915afd1a86003cd6820dae1a63c66faa3276a8d97b056155db6d2`;
changed timing, grid, build, threshold, timeout, or worker options are rejected
before topology, locking, source capture, or building.

Each cell executes three rounds containing two ABBA blocks. The first compares
tables-on/AUTO-on with tables-on/AUTO-off to measure the actual production
route; the second compares tables-on/AUTO-off with tables-off/AUTO-off to expose
table preparation and setup cost. The primary observation is
`metrics.amortized.derived_median_us_per_api_call`; execution, codec setup, and
binding setup remain in every raw record. Route, preparation, and their
round-wise summed net contrast use marginal two-sided 95% Student-t intervals
at two degrees of freedom. A candidate qualifies only if both route and net
lower bounds reach five percent. All controls remain visible, and an unexplained
net regression worse than two percent prevents a passing decision. Reanalysis
cannot lower the preregistered five-percent gate, and every policy invocation
must retain the same non-vacuous input and parity checksums. This mode does not
flip a compiled default or authorize a production promotion; exit status 2
records that the evidence gate was not met.

Each result directory is configuration-specific and resumable. Reusing a
directory with changed source, flags, executables, grid settings, or machine
identity is rejected instead of mixing evidence. Passed jobs retain hashed raw
benchmark JSON and stdout/stderr; analysis revalidates those hashes and semantic
parameters. Replay validates the recorded taskset digest and lock UID
structurally instead of requiring the current host to have the same taskset
bytes or user ID; acquisition-time execution still verifies the live sealed
taskset snapshot and held lock. Replay also keeps the Release/tests/canonical-
Git and required object/archive/link-graph closure checks, while deliberately
omitting only the original source-versus-executable mtime comparison. The
authoritative reader rejects multi-config metadata before consulting that
artifact inventory.

Campaign-specific compact checkpoints preserve authenticated projections, but
do not replace the retained raw trees required for full replay. The LOW AUTO
evidence remains in
`experiments/leopard2/direct_encode/results/checkpoint.json`, the forced sparse
discovery is in
`experiments/leopard2/direct_encode/results/sparse_high_avx2_checkpoint_20260831.json`,
and the production-AUTO decision below is in
`experiments/leopard2/direct_encode/results/production_auto_avx2_checkpoint_20260831.json`.

### Production-AUTO qualification outcome

The 2026-08-31 v10 campaign completed with all 88 jobs structurally valid, but
the preregistered performance gate was not met. Only 56 of 74 candidates had
both production-route and summed-net marginal 95-percent lower bounds at or
above five percent. All 14 route-negative controls were retained, and one had
a net point estimate below the minus-two-percent guard. The final strict resume
and both deterministic replays returned exit status 2; the initial pass returned
1 solely for 11 archived isolation rejections. This is a valid performance
no-promotion result, not a correctness or acquisition failure: sparse-high AUTO
remains compiled-default `OFF`, and no production route or installed API changed.

The decisive cells include:

| Cell | Route gain (95% interval) | Net gain (95% interval) | Decision |
| --- | ---: | ---: | --- |
| one-shot `K=3,R=2`, 4096 B | -27.06% (-31.88% to -21.90%) | -25.13% (-36.02% to -12.39%) | weakest candidate route lower bound; fail |
| binding batch 1 `K=16,R=2`, 4096 B | -23.29% (-27.83% to -18.48%) | -24.54% (-37.25% to -9.25%) | weakest candidate net lower bound; fail |
| one-shot `K=16,R=4`, 4096 B | 21.74% (14.73% to 29.19%) | 21.40% (12.34% to 31.18%) | weakest passing candidate lower bound |
| explicit-AVX2 binding batch 16 `K=2,R=16`, 4096 B | -0.38% (-22.40% to 27.88%) | -2.19% (-24.74% to 27.10%) | control point guard triggered |

The 36 one-shot tuple candidates split 24 passing and 12 failing. The 30
batch/binding candidates split 24 passing and six failing; those six are all
the `K=16,R=2,4096` anchor across batch and binding at batches 1, 4, and 16.
All eight sampled parity-row candidates passed. Passing cells are not averaged
with these failures and do not authorize a post-hoc narrower selector.

The accepted source was commit
`b455108e61017c711e238ccf159d45da50e77ca2`, tree
`4c2caeb67fa22379f9ac9de74e7a4039207fae19`, source fingerprint
`ca8f40ff251c035ac14b222225ac15e12828b827189b092cad77ad607b2a3d01`,
and initialized `sse2neon` gitlink
`cad518a93b326f0f644b7972d488d04eaa2b0475`. Timing used CPU109 on
`foureyes.lan` with exact SMT sibling CPU45 reserved. The first full pass
accepted 77 jobs and rejected 11 solely for sibling work; a quiet strict resume
reran all 24 invocations for those 11 jobs, and every accepted sibling delta
was zero. The 77 already-valid jobs were revalidated and retained byte for byte;
only the 11 rejected job attempts were replaced, and those failed measurements
were neither pooled nor used.

Fresh exact-source GCC 13 gates passed 5/5 focused Release tests and 4/4 focused
ASan+UBSan tests. Every campaign job passed independent-oracle parity,
route/table witnesses, stable non-vacuous input and parity identities,
input/output canaries, process containment, raw/log hashes, frozen-artifact
checks, source stability, and isolation. Normal and `python -O` reanalysis both
returned 2 and reproduced `analysis.json` byte for byte at SHA-256
`a6da9452f530e3ac3af183a268822094c0715da84b06b19a01ad26f8eff6a853`.

The accepted manifest, matrix, analysis, and controlled-build SHA-256 values
are respectively
`db5d472d73ea0646e782665e6960dc9c7c0c736b9987ef4280a80e3a21789d82`,
`346326eff8f0bdc4b402bf641f0a49463cd7d998b272c004e1736317d7306433`,
`a6da9452f530e3ac3af183a268822094c0715da84b06b19a01ad26f8eff6a853`,
and
`7eab103780f3dd29262f8e4d8a75d992b17a1561312c92cbb59fd5c8a7926916`.
The 174,689,363-byte raw working tree remains on `foureyes.lan`; its canonical
6,490-file inventory digest is
`003cbb457cfa057957b98794ea185595d9c1d4a0fbdad7b41363a50fc024b4cd`.
Its 4,401,675-byte accepted archive, acquisition logs, discarded-attempt
archive, and replay records are retained in owner-only durable storage under
`foureyes.lan:/home/catid/leopard-evidence/production-auto-b455108-20260831`;
that 20-file payload is bound by manifest SHA-256
`5e8adb4eba3fe99c71091af4c11afd3793c8ee001badc7450f57c902faadc2e8`.
The compact checkpoint retains every cell and discarded-attempt record, while
[optimization report 40](../experiments/leopard2/optimization_log/40-sparse-high-production-auto-qualification.md)
records the full disposition.

This evidence is bounded to the frozen rows, K/R/byte/API/batch crosses, reuse,
backend, thread counts, host, and OpenMP runtime. Every omitted row or cross,
other reuse count, backend, thread count, out-of-table neighbor, and hardware
target remains unauthorized. The intervals are marginal Student-t intervals at
two degrees of freedom, not a simultaneous campaign guarantee. Any redesigned
predicate or implementation requires a fresh preregistered campaign.

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
When the separately gated sparse-Q1 candidate below is enabled, qualifying
large aligned Q=1 calls use that same mature row-major direct executor; they
never enter this source-major full-output kernel.

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

- candidate: `aadeca79b0df0d6207a9205102dc468e147c131d3ca8f2ab2d2dfa2fd559f150`
- same-layout AUTO-off control:
  `04178d8a1a04e6f040572a73d149b1fe53cca9bbe3e9c87f45a7586c9bc92f1c`
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
| 5,5 | 64 | 0.874 | 1.803 |
| 5,5 | 256 | 0.825 | 1.159 |
| 5,5 | 1,024 | 1.036 | 1.172 |
| 5,5 | 4,096 | 1.289 | 1.056 |
| 8,8 | 64 | 0.503 | 1.288 |
| 8,8 | 256 | 0.659 | 1.015 |
| 8,8 | 1,024 | 0.915 | 1.019 |
| 8,8 | 4,096 | 1.087 | 0.992 |
| 12,8 | 64 | 0.624 | 1.338 |
| 16,8 | 64 | 0.520 | 1.009 |

Thus the arithmetic kernel substantially improves Leopard2's own tiny
full-output path, but fixed safe-API validation/setup and `K*R` direct
arithmetic still leave 64-byte `K>=8` cells behind Leopard main. It is not
promoted merely because it wins against Leopard2's transform fallback.

Rejected variants are retained as negative results rather than production
code: a whole-T=8 callback, two- and three-output grouping, an XMM three-byte
tail, forced source-loop unrolling, and a generated coefficient-major
`K=5,R=5` circuit. Each either failed the five-percent whole-call gate or
regressed through extra source loads, register pressure, or code size.

## Legacy-high GF8/AVX2 sparse-Q1 AUTO experiment

Sparse legacy-high direct encoding has two independent, default-off controls:

```sh
-DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE=ON
-DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE_AUTO=ON
```

The first control prepares the existing bounded generator table at codec
creation; by itself it leaves production AUTO routing unchanged and supports
forced-path measurement. The second control authorizes the experimental AUTO
dispatcher, and is rejected unless table preparation is also enabled.

AUTO dispatch is deliberately narrower than table preparation. It requires a
caller-requested AUTO context that resolves to AVX2, one context thread, GF8
legacy-high native layout, flags zero, and exactly one non-null parity output.
The attested tuple set is exactly:

- 4,096-byte shards, `K` in `{2,3,4,8,12,16}`, and `R` in `{2,4,8,16}`;
- `(K,R)` equal to `(2,16)` or `(16,2)`, with shard bytes in
  `{1024,1088,2048,4032,4160,65536}`.

Explicit backend requests, pooled contexts, extra parity outputs, flags,
ragged sizes, and every unlisted K/R/byte tuple retain the transform path. The
candidate executes the existing allocation-free row-major direct encoder; it
adds no arithmetic kernel and does not enter the full-output source-major
experiment.

This option extracts sparse behavior that the older full-output option had
admitted incidentally through the generic Q=1 selector. It is measurement
machinery, not promotion authority: no performance gain is claimed here, and
the committed direct-encode checkpoint and LOW-profile production AUTO region
remain unchanged. Raw direct-versus-transform measurements can be made with:

```sh
bench_leopard2_direct_encode --profile high --k 2 --r 16 --q 1 \
    --bytes 4096 --force-direct
bench_leopard2_direct_encode --profile high --k 2 --r 16 --q 1 \
    --bytes 4096 --force-transform
```

An AUTO-path experiment must use an attested build with both controls enabled;
parent ON/AUTO OFF is the same-layout route control. Production defaults remain
OFF because the 2026-08-31 production-AUTO qualification did not meet its
evidence gate; see [Production-AUTO qualification outcome](#production-auto-qualification-outcome).

The standalone `bench_leopard2_high_sparse_auto` target supplies the no-hook
telemetry input for that qualification. It links the ordinary production
archive and the independent direct systematic generator oracle, rejects cells
outside the exact 36-tuple sparse-high table, and emits
`leopard2-high-sparse-auto-benchmark-v1`. Its three context-local policies are
`tables-off-auto-off`, `tables-on-auto-off`, and `tables-on-auto-on`; the first
two are production-route controls and the third exercises ordinary AUTO. The
one-shot, batch, and reusable-binding APIs report codec setup, binding setup
where applicable, execution, modeled reuse, requested and effective backend,
thread count, prepared rows, and an untimed route witness that is disabled
before timing.

```sh
bench_leopard2_high_sparse_auto --api one-shot --batch 1 \
    --backend auto --threads 1 --policy tables-off-auto-off
bench_leopard2_high_sparse_auto --api batch --batch 4 \
    --backend auto --threads 1 --policy tables-on-auto-off
bench_leopard2_high_sparse_auto --api binding --batch 4 \
    --backend auto --threads 1 --policy tables-on-auto-on
```

Every raw document says `"authoritative": false`: a pinned paired runner must
still freeze and hash the executable and archive, establish CPU isolation, and
apply the qualification decision rule. That runner is
`sparse-high-production-auto-avx2`, described above. The strict smoke parser
runs normally
and with Python assertions disabled, independently re-derives schema keys,
layout, memory, route counts, setup/amortization formulas, and the full oracle
result. Neither the target nor its tests change the installed API or production
defaults.

## Promoted GF8/AVX2 T=4 batch-table amortization

Dense reusable legacy-high bindings now amortize the fixed AVX2 preparation
cost of the register-resident T=4 encoder. The selector requires GF8, AVX2, no
context worker pool, `R` equal to 3 or 4, `K` in
`{3,4,5,6,7,9,10,11}`, a homogeneous shard length divisible by 32 from 32
through the measured per-shape ceiling below, and every transmitted parity
output present. Any sparse output, heterogeneous byte count, unsupported
shape, shard above its ceiling, or pooled context retains the ordinary
prevalidated executor.

| K | R=3 maximum bytes | R=4 maximum bytes |
| ---: | ---: | ---: |
| 3 | 16,384 | 2,048 |
| 4 | 12,288 | 6,144 |
| 5 | 6,144 | 4,096 |
| 6 | 16,384 | 8,192 |
| 7 | 16,384 | 4,096 |
| 9 | 4,096 | 3,072 |
| 10 | 8,192 | 4,096 |
| 11 | 6,144 | 2,048 |

The arithmetic kernel was already capable of keeping the four inverse
accumulators and final forward transform in registers. The previous binding
loop nevertheless entered that complete encoder once per stripe, rebuilding
the same code-dependent nibble tables each time. The new callback prepares the
three possible inverse-block table groups and the forward table group once per
binding execution, then walks flat item-major pointer arrays. An `R=3`
instantiation stores only the transmitted three-coordinate parity prefix.
Neither variant changes the field, coordinate order, skew constants,
shortening, puncturing, or parity bytes.

Complete AVX2 passes in ordinary encoding below 2 KiB also enter the existing
register-resident kernel. A ragged final public-byte tail is deliberately
excluded from that coarse route: it is padded and executed exactly once by the
tail path, not mistaken for another independent 64-byte encode. Backend startup
tests cover every promoted `K,R` combination with two unaligned guarded
stripes. The production test contributes 201 independent systematic-generator
checks, including 2,049- and 2,080-byte fallbacks, changed source bytes,
punctured `R=3`, and sparse null outputs. A separate allocation counter proves
that repeated binding execution allocates no memory.

The exact candidate source passed four focused GNU 13 Release tests and the
same four tests under Clang 18 ASan, UBSan, and leak detection. GCC 13 and
Clang 18 `-Werror` builds passed, as did GF8-only and GF16-only `-Werror`
libraries. Candidate and diagnostic control binaries reported the same clean
commit and tree and had byte-identical executable ELF sections; only a nonzero
initialized-data selector word differed.

Authoritative timing used CPU 4 with SMT sibling 20 reserved and idle. The
checkpoint ran 2,058 processes over 104 selected cells and twelve 2,112-byte
route-off boundaries. Three noisy cells whose point gains exceeded five
percent received a fresh independent-seed, nine-round holdout of 162 more
processes. After replacing only those predeclared estimates, every selected
same-source cell cleared the five-percent lower-confidence gate: speedups
ranged from 1.062x to 6.834x, the median was 1.270x, and the smallest lower
95-percent bound was 1.051x. Exact Leopard main was available for 99 cells
(`K=3,R=4` is rejected by its `R<=K` API); Leopard2 won all 99, from 1.131x
to 4.474x with a 1.729x median. The smallest exact-main lower confidence bound
was 1.121x.

| K,R | bytes | batch | control / candidate | Leopard main / candidate |
| --- | ---: | ---: | ---: | ---: |
| 4,4 | 64 | 1 | 3.112x | 1.706x |
| 4,4 | 64 | 64 | 5.704x | 3.725x |
| 4,4 | 1,024 | 1 | 1.311x | 1.217x |
| 4,4 | 1,024 | 64 | 1.350x | 1.619x |
| 11,4 | 2,048 | 1 | 1.062x | 1.196x |

These are execution comparisons after reusable setup. Leopard2 executes one
prevalidated binding call containing the stated batch, while exact Leopard
main executes the same number of legacy encode calls. Ordinary one-shot calls
still pay validation and table setup, so this result is not generalized to
every tiny non-binding invocation. At the excluded 2,112-byte boundary,
same-source point ratios ranged from 0.986x to 1.012x and no cell had a
credible regression greater than two percent.

The larger per-shape ceilings were evaluated separately after the common
2,048-byte promotion. A 102-cell screen retained every result, including
thirteen cells that initially missed the five-percent lower-confidence gate.
An independent nine-round holdout supported eight of those cells; unsupported
ceilings were reduced instead of weakening the gate. The final-source
nine-round holdout then accepted all four changed maxima and all four
immediately unselected neighbors. To remove an observed ELF page-placement
bias, candidate and control were hard links to the same immutable executable;
a context-local setup-only diagnostic flag selected the old path before codec
creation.

Across the final 49 selected target cells above 2 KiB, same-source speedups
ranged from 1.059x to 1.265x with a 1.101x median. The smallest lower
95-percent confidence bound was 1.050x. Leopard2 beat exact Leopard main in
all 49 cells by 1.355x to 1.897x with a 1.632x median; its smallest lower
bound was 1.334x. The final-source maxima were:

| K,R | bytes | control / candidate | Leopard main / candidate |
| --- | ---: | ---: | ---: |
| 9,3 | 4,096 | 1.094x `[1.079,1.110]` | 1.568x `[1.541,1.596]` |
| 9,4 | 3,072 | 1.062x `[1.054,1.069]` | 1.675x `[1.659,1.690]` |
| 10,4 | 4,096 | 1.080x `[1.074,1.086]` | 1.426x `[1.416,1.436]` |
| 11,4 | 2,048 | 1.096x `[1.090,1.102]` | 1.465x `[1.457,1.472]` |

The final selector passed the four focused Release tests, three focused Clang
18 ASan/UBSan/LSan tests, strict GCC and Clang builds, and GF8-only and
GF16-only field builds. The reduced-field gate found and fixed an unconditional
reference to the GF8 diagnostic state in a GF16-only context. Relative to the
common 2-KiB implementation, archive text grew by 608 bytes (0.079 percent)
and benchmark executable text by 1,264 bytes (0.133 percent).

The compact result is
`experiments/leopard2/gf8_high_encode/results/`
`t4_batch_checkpoint_20260731.json`. The reproducible runner is
`experiments/leopard2/gf8_high_encode/run_t4_batch_abba.py`; it rejects
changed binary hashes, executable-section differences, dirty or mismatched
source attestation, wrong dispatch, parity-digest differences, incorrect CPU
topology, and any non-idle work observed on the reserved sibling.

The per-shape extension result is
`experiments/leopard2/gf8_high_encode/results/`
`t4_extended_checkpoint_20260731.json`; its runner is
`experiments/leopard2/gf8_high_encode/run_t4_extended_abba.py`.

## Promoted GF8/AVX2 T=8 one-block binding

The reusable prevalidated encode-batch binding now enters the existing
one-block AVX2 transform directly for dense legacy-high GF8 codes with
`5 <= K <= 8`, `5 <= R <= 8`, and shard lengths 64, 128, 192, 256, 320, 384,
448, or 512 bytes. Shortened inputs use an immutable zero shard and punctured
outputs use three bounded stack shards. The shortcut is all-or-nothing across
a heterogeneous binding: one unsupported byte count, backend, shape, or null
parity output retains the mature per-item path for the whole binding. Ordinary
one-shot encoding is unchanged.

Correctness was checked against the independent direct systematic generator,
not against the transform itself. The Release candidate and same-source
control each passed seven focused suites. The candidate covered 1,088 exact
`K=8,R=8` parity checks and 912 checks across the complete
`K=5..8,R=5..8` shape grid, plus changed input bytes, sparse-output
nonselection, unaligned guards, mixed-size bindings, allocation-free repeated
execution, thread-pool execution, and legacy golden vectors. Three focused
Clang 18 ASan+UBSan suites also passed under a 256 MiB process-group limit.

The authoritative CPU-4 campaign used three mirrored rounds, 2,532 fresh
processes, one stable seed per cell, and an idle SMT sibling. Candidate and
control had byte-identical executable sections; only an initialized-data
selector marker differed. All 112 target cells and 64 neighboring cells passed
their gates. Same-source control/candidate speedups ranged from 1.121x to
1.999x. Among the 70 shapes supported by exact Leopard main, Leopard
main/candidate speedups ranged from 1.040x to 1.558x, and the smallest lower
95-percent confidence bound was 1.023x.

| K,R | bytes | same-source control / candidate | Leopard main / candidate |
| --- | ---: | ---: | ---: |
| 5,5 | 128 | 1.774x | 1.383x |
| 5,5 | 320 | 1.245x | 1.090x |
| 5,5 | 512 | 1.143x | 1.113x |
| 8,8 | 128 | 1.999x | 1.558x |
| 8,8 | 320 | 1.328x | 1.187x |
| 8,8 | 512 | 1.177x | 1.212x |

The Leopard-main numbers are an end-to-end API comparison: Leopard2 executes
one prevalidated 64-item binding, while the exact main harness invokes the
legacy API 64 times. They are not presented as a pure inner-kernel comparison.
The promoted transform body is one 3,811-byte VEX-only AVX2 function with no
hot-path calls. The partial-profile wrapper has a fixed 1,696-byte frame, of
which 1,536 bytes are the three discarded-parity shards.

The compact machine-readable checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_one_block_extended_checkpoint_20260730.json`. The full generated raw tree
is intentionally not committed. The checked-in authoritative runner is
`experiments/leopard2/gf8_high_encode/`
`run_t8_one_block_extended_abba.py`; it verifies clean source identity,
immutable binary hashes, all executable ELF sections, per-cell dispatch,
workload digests, CPU/sibling isolation, and failure-output persistence before
accepting a result.

## Promoted K=5,R=5, 1088-byte AVX2 T=8 tile

The only remaining post-1-KiB T=8 Leopard1 loss was the legacy-high GF8
`K=5,R=5`, 1088-byte cell. The padded one-block callback still loaded three
shortened zero rows, stored three punctured parity rows, and required caller
storage for those rows. The promoted fixed-profile kernel keeps the three
zeros in registers, reads exactly five systematic rows, and stores exactly
five parity rows. It performs the same complete parent-code arithmetic and
therefore emits the existing legacy parity bytes.

The backend KAT compares the callback with the generic transform using
unaligned 65-byte rows, covering two AVX2 vectors and the scalar tail. The
production suite compares every parity row with the independent direct
systematic generator at 1088 bytes and covers changed-source reexecution,
sparse-output fallback, unaligned guards, a heterogeneous `{64,1088}`
binding, an all-or-nothing `{64,513}` binding, and four-worker execution.
Focused GNU 13 Release and Clang 18 ASan+UBSan+LSan gates each passed all
three backend, production, and AUTO-route tests.

Two clean-source, independent-seed CPU-4 campaigns retained 5,976 raw
process records in total. CPU 20, the SMT sibling, had zero non-idle jiffies
in every accepted round. Ratios greater than one mean the new kernel is
faster:

| campaign | same-source control / candidate | exact Leopard main / candidate |
| --- | ---: | ---: |
| discovery | 1.385x `[1.379, 1.391]` | 1.430x `[1.406, 1.454]` |
| holdout | 1.381x `[1.334, 1.431]` | 1.434x `[1.407, 1.461]` |

Candidate and control had identical executable sections; a nonzero data
selector disabled the bound path in the control. All workload digests matched,
and neither seed family had a regression-neighbor failure. The broad runner
also requalified every older 576-through-1024 selector, so its aggregate status
was rejected when several unchanged historical targets missed the five-percent
gate during these repeats. That aggregate result is retained rather than
hidden, but it does not invalidate the independently gated new 1088-byte cell.

The GCC 13 kernel is 1,954 bytes and uses a 256-byte stack frame. Assembly
review found no out-of-bounds, aliasing, or ISA defect, but it did identify the
next bounded avenue: three dead forward updates, two hot-loop YMM spills, and
an avoidable 3,232-byte padded-fallback caller frame. The compact checkpoint,
including frozen hashes, confidence intervals, test commands, and raw-result
hashes, is
`experiments/leopard2/gf8_high_encode/results/`
`t8_k5r5_1088_checkpoint_20260730.json`.

## Spill-free fused K=5,R=5 T=8 circuit

The first fixed-profile kernel above retained two pairs of YMM spills in its
32-byte loop. A direct dependency derivation removes those spills without
changing the parent transform. Write `Mx(v)` for bytewise multiplication by
the legacy GF8 element at log-table index `x`. After the lower inverse
radix-4, the shortened upper input has the exact form

    [E, 0, 0, 0] -> [M119(E), M136(E), M238(E), E].

For each distance-four inverse/forward pair, the intermediate values needed by
the requested outputs reduce to

    t = u XOR s
    z = u XOR M85(t)
    b = t XOR z.

The four lower outputs then use

    o0  = u0
    o1  = u1 XOR u0
    l02 = u2 XOR u0
    l13 = u3 XOR u1
    o2  = l02 XOR M85(l13)
    o3  = l13 XOR o2,

and the sole requested upper output uses

    c4 = b0 XOR M85(b2)
    c5 = b1 XOR M85(b3)
    o4 = c4 XOR M17(c5).

The implementation scatters `o0..o3` while the log-85 tables are live, ending
those vector lifetimes before computing `o4`. GCC 13 therefore emits an
835-byte multiple-of-32 core with no stack frame, stack canary, or vector
spill. The predecessor was 1,954 bytes with a 256-byte frame and 128 bytes of
hot stack traffic per 32 input bytes. A separate scalar-tail function retains
arbitrary positive byte lengths and the backend KAT covers an unaligned
65-byte row. The 3,232-byte padded top-level caller frame remains: splitting
it measured approximately neutral and was reverted.

A clean nine-round, pinned holdout at 1088 bytes measured `1.0509x`
`[1.0455,1.0563]` over the prior Leopard2 kernel and `1.5237x`
`[1.5107,1.5369]` over exact Leopard main. The latter compares one
prevalidated Leopard2 call for 64 stripes with 64 legacy API calls, so it is an
end-to-end API result rather than an inner-kernel ratio. Neighbor campaigns at
1024, 1087, 1089, and 1152 bytes found no credible regression greater than two
percent. The mean gain clears five percent, while the lower confidence bound
misses it narrowly at 4.55 percent; promotion uses the task's simple-circuit
exception because the assembly defect is eliminated, the exact-main win is
large, and dispatch outside the exact target cell is unchanged.

The final source passes the three focused GNU 13 Release tests, the same tests
under Clang 18 ASan+UBSan+LSan, and strict GCC/Clang warning builds. The
machine-readable derivation, frozen binary identities, raw-manifest hashes,
round results, neighbor screen, validation commands, and rejected variants are
in
`experiments/leopard2/gf8_high_encode/results/`
`t8_k5r5_spillfree_checkpoint_20260731.json`.

## Promoted T=8 one-kibibyte selector extension

The reusable prevalidated legacy-high GF8/AVX2 binding now also selects the
existing T=8 direct-input transform at exactly 1024 bytes for four additional
profiles:

    (K,R) = (6,5), (6,6), (10,8), (16,5).

This is a dispatch extension, not a new arithmetic kernel or wire profile.
Candidate and diagnostic-control executables had byte-identical instruction
sections; a read-only initialized-data marker disabled only these four routes
in the control. All other legal `K=5..16`, `R=5..min(K,8)` shapes were
neighbor controls, including profiles already selected by earlier production
tables.

The final current-source CPU-4 campaign used nine counterbalanced target
rounds, three neighbor rounds, 672 fresh processes, batch/reuse 64, and a
reserved idle SMT sibling. Every original, parity, and recovered-output digest
matched. All four targets cleared the five-percent lower-confidence
same-source gate, all four beat exact Leopard main, and none of the 38
neighbors failed the two-percent regression gate:

| K,R | control / candidate, 95% CI | exact main / candidate, 95% CI |
| --- | ---: | ---: |
| 6,5 | 1.0615x `[1.0552,1.0679]` | 1.1396x `[1.1300,1.1492]` |
| 6,6 | 1.0631x `[1.0573,1.0689]` | 1.1404x `[1.1302,1.1508]` |
| 10,8 | 1.0765x `[1.0687,1.0843]` | 1.1582x `[1.1503,1.1661]` |
| 16,5 | 1.0932x `[1.0828,1.1038]` | 1.3028x `[1.2934,1.3123]` |

Two attractive point estimates were deliberately excluded. `K=11,R=5`
measured 1.0601x but had lower bound 1.0475x, and `K=11,R=6` measured
1.0552x with lower bound 1.0477x. Their rejected summary and raw hashes remain
in the compact result rather than being hidden.

The final selector passed four focused GNU 13 Release tests, the feature-off
production test, the GF8-only production test, and the production test under
Clang 18 ASan+UBSan+LSan. The compact current-source checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_one_kib_checkpoint_20260731.json`; the checked-in runner is
`experiments/leopard2/gf8_high_encode/run_t8_two_block_abba.py`.

## Promoted T=8 ragged 65--928-byte selector

The same prepared legacy-high GF8/AVX2 binding now covers the following
non-power-of-two byte tiers for dense `K=5..16`,
`R=5..min(K,8)` profiles:

    65..191, 193..224, 257..352, 416, 449..480,
    513..544, 577..608, 641..672, 736, 769..800,
    864, 897..928.

This changes only the deterministic selector. It reuses the already-audited
one- and two-block callbacks and does not change arithmetic, wire bytes, or
scratch. Five measured shape/byte combinations remain on the mature transform
path:

    (K,R,bytes) = (5,5,191), (6,5,319), (6,6,319),
                  (7,5,319), (7,6,319).

The streamed discovery and predeclared holdout covered 2,058 cells and 32,004
accepted benchmark processes. Applying the five exclusions yields 1,213
qualified target cells and 845 route neighbors. Across that population the
same-source selector-off/control speedup has a `2.0296x` geometric mean and
`1.0504x` minimum lower 95-percent confidence bound. The application-equivalent
padded exact-main speedup has a `1.2851x` geometric mean and `1.000026x`
minimum lower bound. Neighbor control/candidate results retain the required
two-percent safety floor.

A separate final-source confirmation at commit
`a808eac051469330583b7241d2591e30d7b6a354` exercised 14 targets and five
exclusions in 816 fresh processes. All digests and isolation gates passed. Its
same-source geometric mean is `1.7715x`, and its padded exact-main geometric
mean is `1.0598x`; their minimum lower bounds are `1.0834x` and `1.00091x`.
The fifth exclusion, `(7,6,319)`, was added because its earlier exact-main
lower bound was `0.99897x`, even though it was substantially faster than the
selector-off Leopard2 path.

Fresh selector-on and selector-off GNU 13 Release gates and selector-on
Clang 18 ASan+UBSan+LSan gates passed. Binary hashes were verified before and
after the campaign, and the reserved SMT sibling remained idle. Full
provenance, confidence summaries, canonical population digests, rejected
preliminary runs, and raw-result hashes are in
`experiments/leopard2/gf8_high_encode/results/`
`t8_ragged_checkpoint_20260731.json`.

## Promoted T=8 arbitrary-tail selector

The reusable legacy-high GF8/AVX2 binding now covers every legal
`K=5..16`, `R=5..min(K,8)` profile for shard byte counts
`1,2,3,7,8,15,16,17,31,32,33,63`. The one-block callback executes one
shifted inverse transform and one forward transform. The two-block callback
fuses both shifted inverse transforms, their coefficient XOR, and the final
forward transform. For 33 through 63 bytes, an overlapping final in-range
32-byte vector avoids temporary tail storage; shorter inputs use bounded
staging.

The three-round frozen campaign covered 504 target cells and 84 route
neighbors in 10,080 accepted processes. Two exact-main comparisons were
statistically inconclusive, so only those two predeclared cells were repeated
for nine rounds. The combined disposition has no target or neighbor failure:

| comparison | geometric-mean speedup | minimum point | minimum 95% lower bound |
| --- | ---: | ---: | ---: |
| Same-source selector-off control / Leopard2 | 2.9609x | 1.8349x | 1.6586x |
| Padded exact Leopard main / Leopard2 | 1.5669x | 1.1947x | 1.0263x |

Exact Leopard main cannot process a sub-64-byte shard. The second row
therefore times a legal, zero-padded 64-byte legacy call while using the same
logical input stream and comparing original, parity, and recovered-output
digests over the requested prefix. It is an application-equivalent comparison,
not a claim that both codecs physically processed the same number of bytes.
The same-source comparison does use the exact requested byte count in both
executables.

The production-default binary from commit
`59a745da8c6ecbc4342c64da428f762c3e16df36` has the same executable-section
digest as the frozen candidate. Fresh default Release, feature-off, GF8-only,
GF16-only, and Clang 18 ASan+UBSan+LSan gates passed. Full provenance,
confidence intervals, rejected preliminary runs, and raw-result hashes are in
`experiments/leopard2/gf8_high_encode/results/`
`t8_tiny_checkpoint_20260731.json`.

## Promoted T=32 two-message-block 256-byte encoder

The legacy-high GF8/AVX2 encoder now has a bounded complete transform for
exactly `K=64`, `R=T=32`, and 256-byte shards.  It inverse-transforms the first
32-message block directly into the recovery accumulator, folds the shifted
inverse transform of the second block through one 32-row temporary slab, and
finishes with one forward transform.  Packed base pointers replace the general
64-entry work-pointer traversal.  The wire profile, skew constants, arithmetic,
public scratch size, and parity bytes are unchanged.

The ordinary encode and multi-item batch entry points retain their full
address and alias validation.  A reusable batch binding proves the packed
layout once during setup and stores the work offset, so repeated execution is
byte-heavy arithmetic only.  This setup/execution separation matters for a
fair Leopard1 comparison: Leopard1's encode call has no equivalent batch-wide
alias preflight.  Setup remains explicitly measured rather than treated as
free.

The frozen CPU-13 campaign used exact Leopard main at commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, nine counterbalanced target
rounds, three neighbor rounds, and 216 fresh processes.  No process or round
was discarded, every workload/parity/recovery digest matched, the reserved
SMT sibling remained idle, and no byte/count neighbor had a credible
two-percent regression:

| execution path | median candidate | setup | control / candidate, 95% CI | exact main / candidate, 95% CI |
| --- | ---: | ---: | ---: | ---: |
| reusable binding, batch 1 | 0.832 us | 0.500 us | 1.1230x `[1.1164,1.1296]` | 1.1046x `[1.0969,1.1124]` |
| reusable binding, batch 8 | 6.756 us | 4.455 us | 1.1486x `[1.1418,1.1554]` | 1.1878x `[1.1855,1.1901]` |
| ordinary one-shot, batch 1 | 0.862 us | included | 1.1032x `[1.0975,1.1089]` | 1.0680x `[1.0609,1.0751]` |

The generated function is 9,671 bytes and its object audit found no EVEX,
ZMM, opmask, or YMM stack-spill instruction.  Arbitrary shard layouts, sparse
parity requests, adjacent byte/count shapes, non-AVX2 backends, and bindings
with a worker pool retain the mature transform or scheduler.  Full provenance
and confidence intervals are in
`experiments/leopard2/gf8_high_encode/results/`
`t32_two_block_b256_checkpoint_20260803.json`.

The fresh production-default Release archive is byte-identical to the frozen
candidate archive.  Default and selector-disabled direct-oracle tests,
two-worker batch execution, general batch-alias and context/backend tests,
GF8-only and GF16-only builds, GCC 13 and Clang 18 `-Werror` builds, and the
Clang 18 ASan+UBSan focused test pass.  The feature-off and GF16-only archives
contain neither the generated object nor its symbol.
