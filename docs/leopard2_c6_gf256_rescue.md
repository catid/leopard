# Leopard2 C6 exact-GF256 field-boundary rescue

## Outcome

C6 identifies a compelling exact-GF8 region, but promotes no production code,
wire profile, or dispatcher rule.

For transmitted codes with `K+R <= 256`, both frozen dyadic profiles have
10,795 parameter cells whose parent is 512 and therefore selects GF16.  A
separately defined exact prefix-evaluation code over Leopard's legacy GF8
representation is much faster when the smaller side is two through four shards
on this host.  Every one of 24 target cells with smaller side three clears the
10% credible-gain rule, and every one of 14 side-two/four neighboring cells
avoids the 2% regression limit.  The wide side-17 comparison is mixed: only
two of twelve encoder cells clear 10%, with credible gains from -22.07% to
+16.85%.  This is a region result, not a universal exact-size win.

The original-repair executor clears the same target and neighbor gates against
the forced specialized padded-GF16 decoder.  All 56 decoder cells are positive,
including every maximum-loss cell.  This is deliberately an algorithm/profile
comparison, not a claim against the production AUTO dispatcher: a later C10
run must also compare GF16 direct repair and AUTO selection.

The candidate remains experimental because it changes field, coordinate set,
generator matrix, and parity bytes.  C6 now includes an immutable direct
original-repair plan, but it has no frozen serialized code identity, mixed
original/parity-erasure planner, second-host timing, production API review, or
calibrated dispatcher.  The default build does not compile the C6 sources;
`CMakeLists.txt`, the public ABI, and production AUTO selection are unchanged.

## The exact code studied

The candidate is a systematic prefix-evaluation Reed-Solomon code in Leopard's
legacy GF8 representation.  Integers `0..255` below are the symbol bit patterns
of Leopard's legacy Cantor-coordinate representation: `i xor s` denotes
`omega_i + omega_s`, not XOR of polynomial-basis integer encodings.

    source evaluation coordinates = 0 .. K-1
    parity evaluation coordinates = K .. K+R-1
    polynomial degree              < K

For source point `x_i=i`, define

    w_i = inverse(product over 0 <= s < K, s != i of (i xor s))

and for parity point `q=K+j`, define

    Z(q)   = product over 0 <= s < K of (q xor s)
    G(j,i) = Z(q) * inverse(q xor i) * w_i.

Encoding is `p_j = sum_i G(j,i) d_i`.  Setup stores all `K*R` nonzero
coefficients in Leopard logarithmic form.  Execution uses the production
fixed-multiplier backend, performs no allocation, and needs no transform
scratch.  The table is deliberately dense: C6 is a direct crossover experiment,
not an assertion that dense evaluation scales to balanced codes.

For `L` missing originals, decode setup deterministically selects parity rows
`0..L-1`, inverts their `L`-column repair minor, and folds each recovered output
into fixed coefficients over those parity shards and every surviving original.
The immutable plan stores only the nonzero terms and output offsets.  Execution
performs fixed multiply-adds, restores missing originals only, performs no
allocation, and requires no scratch.  A no-loss plan returns before inspecting
any shard pointer.  The current checkpoint assumes the selected parity prefix
is present; mixed original/parity loss selection remains future profile work.

All `K+R` evaluation points are distinct and lie in GF8.  Hence any `K`
evaluation rows form a nonsingular generalized Vandermonde matrix, proving the
effective `[K+R,K]` code is MDS.  The experiment does not rely on that theorem
alone: GF(2^4) rank-tests all 65,534 subsets of `K` coordinates in the full
16-coordinate code across `1 <= K < 16`, and legacy-GF8 checks 80 independently
selected repair minors over 24 boundary geometries.  The Python decoder oracle
also constructs 16 loss-pattern plans independently, inverts their minors, and
checks all 4,933 folded nonzero execution terms.  The byte-heavy C++ plan is a
separate implementation rather than a relabelled rank test.

The coordinate-set SHA-256 values in `results/algebra.json` are research
fingerprints prefixed with `LEO2-EXACT-PREFIX-GF8-UNFROZEN`.  They are not
serialized production identifiers or a compatibility promise.  Experiment W
explicitly excludes the unfrozen exact profile.  C7/C8 must assign a new
profile family/version and coordinate-map identity before any promotion.

## Why this cannot be legacy compatible

In every target cell the declared high or low parent has 512 coordinates.  A
GF8 field contains only 256 elements, so no GF8 truncated transform can evaluate
that same parent coordinate set.  Truncating computation while retaining the
legacy wire must remain in GF16 and cannot provide a GF8 field-boundary rescue.
Changing to the 256-point GF8 prefix code changes both the alphabet and the
interpolation set.

The independent algebra artifact records an explicit unequal generator
coefficient for every one of its 24 GF8/GF16 boundary cases.  Every measured
stripe also produces unequal exact-GF8 and padded-GF16 output digests.  These
are incompatibility witnesses, not merely different implementation paths.

## What was compared

| Requested C6 method | Executable disposition |
| --- | --- |
| Legacy padded GF16 | Measured through the public Leopard2 API with the requested frozen high/low profile and explicit GF16 field. |
| Exact/direct systematic GF8 | Implemented as the standalone allocation-free direct encoder described above; independently checked byte by byte. |
| Wire-compatible truncated GF8 | Impossible for the target 512-coordinate parents. C2 remains a same-field parent-preserving experiment and was not relabeled. |
| Tang--Han epsilon GF8 | C3b implements Appendix-A Algorithm 8 as exact inverse algebra. Its best scalar inverse gain was 0.862 and it has no end-to-end exact codec/profile identity, so it was not relabeled as a C6 encoder. |
| Binary-block GF8 | C5 rejected the materialized join; a fully fused exact join yields the dense direct evaluator measured here. No duplicate method name is used for the same map. |
| Generic GF8 interpolation/evaluation | Used as the independent exact polynomial oracle. Once setup is compiled into `G`, it has the same code definition as direct evaluation; per-stripe generic interpolation is not presented as a separate optimized plan. C3b's generic inverse best gain was 1.066, below the 10% rule. |

Decoder execution compares that exact GF8 plan with the public Leopard2
specialized decoder, forced to the declared frozen profile and explicit GF16
field.  Each side first encodes its own profile, then receives identical source
payloads and loss locations; equality is required only for restored originals,
not for the intentionally different parity bytes.  Plan setup and byte-heavy
execution are timed separately.

This separation matters: fewer field operations in a prior scalar transform
stage are not an end-to-end GF8 codec measurement, and successful decoding does
not make a new generator matrix legacy compatible.

## Correctness evidence

The independent Python algebra uses carryless polynomial multiplication with
polynomial `0x11d` behind the eight declared Cantor basis values.  The C++
candidate separately builds rows through production GF8 field helpers and then
compares them with another carryless polynomial-basis implementation.

| Gate | Evidence |
| --- | ---: |
| Affected public pairs enumerated | 10,795 high + 10,795 low |
| Exhaustive GF(2^4) full-length generators | 15 |
| Exhaustive GF(2^4) `K`-coordinate subsets | 65,534 |
| Implied GF(2^4) public `(K,R)` geometries | 120 |
| GF8 boundary geometries | 24 |
| GF8 nonzero parity coefficients | 123,596 |
| GF8 repair-minor rank checks | 80 |
| GF8 direct evaluation comparisons | 68 |
| Explicit GF8/GF16 coefficient witnesses | 24 |
| Independent exact decode plans | 16 |
| Independently folded nonzero decode terms | 4,933 |
| Independently recovered scalar values | 70 |
| C++ boundary cases per backend | 12 |
| C++ coefficients per backend | 66,920 |
| C++ byte comparisons per backend | 278,074 |
| Deterministic C++ digest | `0xc58c8188e359fdf2` |
| C++ decode executions per backend | 320 |
| C++ recovered-byte comparisons per backend | 236,500 |
| Deterministic C++ decode digest | `0x29158fba4df259c1` |
| No-loss NULL-pointer no-op calls per backend | 96 |
| Maximum-original-loss cases per backend | 96 |
| Tracked hot-path allocations | 0 |
| Backends with identical evidence | scalar, SSSE3, AVX2, AUTO |
| ASan+UBSan | pass, leak detection enabled |

The encoder checks use 257-byte unaligned GF8 buffers, one-byte guards on both
sides, and an independently derived field product for every expected output
byte.  Decoder lengths are `1,2,3,7,31,64,65,257` bytes.  They cover both
`R>K` and `K>R`, losses `1,2,4` where legal, and the maximum number of missing
originals.  Guards and untouched surviving-output buffers are checked after
every call.  Normal scalar/SIMD artifacts interpose global `new` while executing
the plan; the ASan+UBSan artifact disables that interposition explicitly so the
sanitizer retains its allocator-mismatch detector.  Strict GCC 13.3 compilation
uses `-Wall -Wextra -Wpedantic -Werror`.

The merger independently derives the exact unique 50-cell encode and 56-cell
decode benchmark geometries, parent sizes, coefficient/plan term counts,
payload/input/output bytes, selected parity prefixes, maximum-loss coverage,
and public GF16 scratch layouts.  Eighteen test methods (one happy path and 17
adversarial groups) exercise 67 evidence mutations covering coordinated
derived-field changes, missing/extra/duplicate/relabelled cells, nonfinite
uncertainty, accounting forgeries, encode incompatibility and decode equality,
backend/affinity/threading changes, allocation-instrumentation labels,
source/core/library changes, sanitizer labels, and exact executable/result/log
manifest bindings.

## Pinned AVX2 result

This checkpoint is deliberately bound to baseline core commit
`48803c06fbd7a6802b4438af60e3104895938c9d`, a historical pre-SIMD-follow-up
snapshot.  The executable and linked archives record that SHA.  These numbers
must not be presented as performance of the later/current root or used to
override its dispatcher without a fresh comparison after integration.

The authoritative run used GCC 13.3 and forced AVX2 on an AMD Ryzen 9 9950X3D.
The runner pinned the benchmark child to allowed CPU 15; sibling CPU 31 was
reserved idle.  `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, and OpenMP max threads
was one.  Compilation and other memory-intensive work were stopped.  Setup is
outside execution.  Each cell uses alternating order, two warmups, nine samples,
and median/MAD.  Credible gain is

    100 * ((padded_median - padded_MAD) /
           (exact_median + exact_MAD) - 1).

The 50 cells cover 64 B and 1 KiB with batch eight, 64 KiB with batch one,
reuse counts one through eight, and 1 MiB for the two/three-shard side.  Counts
129, 130, and 131 exercise adjacent just-over-power-of-two cases.  Side two and
four are the target's direct neighbors; side 17 records the crossover failure.

| Profile/cell | Bytes; batch/reuse | Exact us | Padded GF16 us | Credible gain |
| --- | ---: | ---: | ---: | ---: |
| high `(3,129)` | 64; 8/8 | 12.439 | 60.338 | +379.32% |
| high `(3,129)` | 1 KiB; 8/4 | 47.956 | 237.048 | +377.43% |
| high `(3,129)` | 64 KiB; 1/1 | 347.184 | 1,731.779 | +377.42% |
| high `(3,129)` | 1 MiB; 1/1 | 9,322.334 | 67,216.130 | +614.32% |
| low `(129,3)` | 64; 8/8 | 12.301 | 66.271 | +434.14% |
| low `(129,3)` | 1 KiB; 8/4 | 45.251 | 261.201 | +469.98% |
| low `(129,3)` | 64 KiB; 1/1 | 611.287 | 2,687.880 | +240.85% |
| low `(129,3)` | 1 MiB; 1/1 | 10,534.627 | 82,541.852 | +671.02% |
| high `(17,129)` | 64; 8/8 | 82.465 | 65.411 | -22.07% |
| high `(17,129)` | 1 KiB; 8/4 | 234.515 | 233.853 | -0.81% |
| low `(129,17)` | 64; 8/8 | 71.727 | 70.628 | -2.75% |
| low `(129,17)` | 64 KiB; 1/1 | 2,133.144 | 2,720.790 | +16.85% |

Across target side-three encoder cells, credible gain ranges from +225.58% to
+671.02% with median +406.26%.  Side-two/four neighbors range from +237.82% to
+1,011.70%.  The wide-side losses prevent a universal direct threshold and
validate the region-based premise of C10.

The decoder matrix covers losses `1,2,3,4,17`, including 36 maximum-loss cells.
Representative execution-only results are:

| Profile/cell; losses | Bytes; batch/reuse | Exact GF8 us | Specialized GF16 us | Credible gain |
| --- | ---: | ---: | ---: | ---: |
| high `(3,129)`; 1 | 64; 8/8 | 0.095 | 175.017 | +179,126.54% |
| high `(3,129)`; 3 | 64; 8/8 | 0.301 | 150.019 | +49,039.63% |
| high `(3,129)`; 3 | 64 KiB; 1/1 | 12.450 | 5,287.899 | +34,609.65% |
| high `(3,129)`; 3 | 1 MiB; 1/1 | 331.934 | 182,395.337 | +52,868.34% |
| low `(129,3)`; 1 | 64; 8/8 | 4.330 | 151.134 | +3,364.51% |
| low `(129,3)`; 3 | 64; 8/8 | 13.255 | 163.893 | +1,126.59% |
| low `(129,3)`; 3 | 64 KiB; 1/1 | 735.288 | 6,551.563 | +660.94% |
| low `(129,3)`; 3 | 1 MiB; 1/1 | 10,555.788 | 190,857.082 | +1,681.66% |
| high `(17,129)`; 17 | 64 KiB; 1/1 | 300.423 | 5,783.454 | +1,763.01% |
| low `(129,17)`; 17 | 64; 8/8 | 77.963 | 198.022 | +153.34% |

Across the 24 side-three decoder cells, execution credible gain ranges from
+657.51% to +179,126.54%, with median +19,516.05%.  Side-two/four neighbors
have a minimum +529.37%; side-17 cells range from +153.34% to +29,778.58%.
Every maximum-loss cell is positive, from +153.34% to +107,772.10%.  The very
large percentages reflect comparison with a forced full specialized GF16
decoder, not with AUTO or the production direct-repair path.

Exact decode-plan setup medians range from 0.08 to 74.39 microseconds; padded
specialized-plan setup ranges from 1.40 to 51.67 microseconds.  Dense matrix
setup is the exact candidate's weak point at 17 losses: for low `(129,17)`, 17
losses, 64-byte shards, exact setup is 74.07 microseconds versus 7.45, while
execution is 77.96 versus 198.02.  Including setup once, every measured cell
still favors exact GF8 by median time: the minimum one-shot ratio is 1.351x and
the matrix median is 33.37x.  Amortizing setup by each cell's declared reuse
count (one, four, or eight) raises those figures to 2.281x and 36.12x.  These
are ratios of medians, not additional confidence bounds.

Exact encoder-plan setup medians range from 1.09 to 51.01 microseconds; padded
codec setup ranges from 38.04 to 75.89 microseconds.  Exact coefficient tables range
from 258 to 2,210 bytes and execution scratch is zero.  The public padded GF16
encoder scratch ranges from 38,976 bytes to 536,878,272 bytes in this matrix.
Caller-owned source and parity buffers are excluded from both scratch figures.
Exact decode execution also uses zero scratch; its folded plan payload ranges
from 60 to 26,596 bytes, versus 51,520 to 676,343,936 bytes of public padded
GF16 decode scratch.  The large scratch reductions are specific to the dense
direct execution plans.

## Evidence binding and reproduction

The retained benchmark manifest binds SHA-256 hashes of the exact executable,
linked `liblibleopard.a`, C++ source, runner, raw result, stdout, and stderr.
The checkpoint also freezes the exact executable and library hashes and rejects
artifact relabeling.  The root default build remains unaffected.

From the repository root, use at most eight build jobs.  Build forced variants
as in the C5/C2 experiment docs, then compile `c6_gf256.cpp` against each static
archive with the four `LEO2_C6_*` provenance macros shown in the source.  Run
correctness-only artifacts with:

    OMP_NUM_THREADS=1 build/c6/c6-avx2 --backend avx2 \
      experiments/leopard2/non_power_of_two/c6/results/avx2.json \
      --correctness-only

Run the sanitizer candidate with both the library and standalone source built
using ASan+UBSan, `LEO2_C6_SANITIZER_MODE="asan-ubsan"`, and
`LEO2_C6_DISABLE_GLOBAL_NEW_TRACKING=1`.  Run it with
`ASAN_OPTIONS=detect_leaks=1:halt_on_error=1` and
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`.  Normal backend builds keep
the global-new allocation tracker enabled.

For authoritative timing, first reserve one physical core and its SMT sibling;
then run only:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c6/run_benchmark.py \
      --cpu 15 --executable build/c6/c6-avx2 \
      --library build/c6/avx2/liblibleopard.a \
      --result experiments/leopard2/non_power_of_two/c6/results/benchmark.json \
      --manifest experiments/leopard2/non_power_of_two/c6/results/benchmark-manifest.json \
      --stdout experiments/leopard2/non_power_of_two/c6/results/benchmark.stdout.txt \
      --stderr experiments/leopard2/non_power_of_two/c6/results/benchmark.stderr.txt

Regenerate algebra and validate the merger:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c6/algebra.py \
      --output experiments/leopard2/non_power_of_two/c6/results/algebra.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c6/test_checkpoint.py -v

The full checkpoint command is the `checkpoint.py` invocation documented by
its `--help`; it consumes all four backend artifacts, ASan+UBSan, algebra, the
benchmark result, manifest, and both retained logs.  Finally run
`sha256sum -c experiments/leopard2/non_power_of_two/c6/SHA256SUMS`.

## Disposition and next work

Recommend closing Bead `leopard-79h.18.1.7` as a completed experiment, without
promoting production code.  Its acceptance bullets resolve as follows:

- The algebra enumerates all 10,795 affected `(K,R)` cells per frozen profile;
  the hash-bound runner measures representative power-boundary, byte, batch,
  reuse, and loss cells against padded GF16.
- Direct systematic GF256 (GF8) encode and original repair are independently
  implemented and measured end to end.
- Wire-compatible truncated GF8 is impossible for a 512-coordinate parent;
  same-wire pruning remains the same-field GF16 C2 result and is not relabelled.
- Tang--Han epsilon GF8 is retained as C3b's negative inverse-stage result; it
  has no frozen end-to-end code identity and is not relabelled as C6.
- C5 found the materialized binary-block join dense.  Its fully fused exact join
  is the direct generator evaluated here, so C6 does not double-count it as a
  distinct algorithm.
- Generic interpolation/evaluation supplies the independent oracle.  C3b's
  generic inverse best gain of 1.066 is below the 10% promotion threshold.
- Exhaustive small-field MDS, independent GF8 algebra, coordinate fingerprints,
  and explicit GF8/GF16 inequality witnesses establish both correctness and new
  wire identity.
- The direct candidate clears the 10% encode and decode target-region gates,
  but no changed generator is exposed as a public or serialized profile.  Its
  identity remains explicitly unfrozen, so production promotion waits for
  C7/C8, C10 dispatch, and W serialization.

The previous reason to keep C6 open was the missing byte-heavy decoder.  That
gap is now closed; remaining work belongs to the named dependent experiments
and should not keep the C6 comparison bead open.

- Preserve the exact GF8 direct encoder and repair plan as default-off C7-C9
  input.
- Do not describe it as legacy compatible or route production AUTO to it.
- Freeze and serialize a new profile identity before producing persistent data.
- Extend plan selection to mixed original/parity erasures and compare a
  transform alternative before any codec promotion.
- Add a second CPU family and the full C10 dispatcher matrix.  The current
  side-17 losses must remain visible.
- Revisit Tang--Han or a factored exact transform only as a distinct executable;
  prior inverse-stage operation counts are not end-to-end codec evidence.
