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
four of twelve cells clear 10%, with credible gains from -26.12% to +34.47%.
This is a region result, not a universal exact-size win.

The candidate remains experimental because it changes field, coordinate set,
generator matrix, and parity bytes.  It has no frozen serialized code identity,
production erasure decoder, second-host timing, or calibrated dispatcher.  The
default build does not compile the C6 sources; `CMakeLists.txt`, the public ABI,
and production AUTO selection are unchanged.

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

All `K+R` evaluation points are distinct and lie in GF8.  Hence any `K`
evaluation rows form a nonsingular generalized Vandermonde matrix, proving the
effective `[K+R,K]` code is MDS.  The experiment does not rely on that theorem
alone: GF(2^4) rank-tests all 65,534 subsets of `K` coordinates in the full
16-coordinate code across `1 <= K < 16`, and legacy-GF8 checks 80 independently
selected repair minors over 24 boundary geometries.  A production decoder would
still need an immutable byte-execution plan and end-to-end loss benchmarks; C6
does not disguise the rank oracle as that missing implementation.

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
| C++ boundary cases per backend | 12 |
| C++ coefficients per backend | 66,920 |
| C++ byte comparisons per backend | 278,074 |
| Deterministic C++ digest | `0xc58c8188e359fdf2` |
| Backends with identical evidence | scalar, SSSE3, AVX2, AUTO |
| ASan+UBSan | pass, leak detection enabled |

The C++ checks use 257-byte unaligned GF8 buffers, one-byte guards on both
sides, and an independently derived field product for every expected output
byte.  Strict GCC 13.3 compilation uses `-Wall -Wextra -Wpedantic -Werror`.

The merger independently derives the exact unique 50-cell benchmark geometry,
parent size, coefficient table, term count, payload/input/output bytes, and
public GF16 scratch layout.  Twelve test methods exercise 35 adversarial
evidence mutations covering coordinated
derived-field changes, missing/extra/duplicate/relabelled cells, nonfinite
uncertainty, accounting forgeries, output-identity forgeries, backend/affinity/
threading changes, source/core/library changes, sanitizer labels, and exact
executable/result/log manifest bindings.

## Pinned AVX2 result

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
| high `(3,129)` | 64; 8/8 | 12.468 | 60.957 | +381.30% |
| high `(3,129)` | 1 KiB; 8/4 | 47.611 | 236.773 | +379.58% |
| high `(3,129)` | 64 KiB; 1/1 | 352.814 | 1,752.150 | +373.86% |
| high `(3,129)` | 1 MiB; 1/1 | 9,388.161 | 68,502.033 | +623.01% |
| low `(129,3)` | 64; 8/8 | 12.921 | 64.706 | +399.17% |
| low `(129,3)` | 1 KiB; 8/4 | 46.996 | 275.746 | +461.05% |
| low `(129,3)` | 64 KiB; 1/1 | 759.739 | 3,222.148 | +237.53% |
| low `(129,3)` | 1 MiB; 1/1 | 12,959.893 | 93,967.875 | +601.62% |
| high `(17,129)` | 64; 8/8 | 87.617 | 65.628 | -25.28% |
| high `(17,129)` | 1 KiB; 8/4 | 238.940 | 233.625 | -2.57% |
| low `(129,17)` | 64; 8/8 | 74.710 | 73.248 | -2.51% |
| low `(129,17)` | 64 KiB; 1/1 | 2,386.058 | 3,550.742 | +34.47% |

Across target side-three cells, credible gain ranges from +229.23% to +623.01%
with median +392.25%.  Side-two/four neighbors range from +189.34% to +954.14%.
The wide-side losses prevent a universal direct threshold and validate the
region-based premise of C10.

Exact plan setup medians range from 1.11 to 51.86 microseconds; padded codec
setup ranges from 37.91 to 76.06 microseconds.  Exact coefficient tables range
from 258 to 2,210 bytes and execution scratch is zero.  The public padded GF16
encoder scratch ranges from 38,976 bytes to 536,878,272 bytes in this matrix.
Caller-owned source and parity buffers are excluded from both scratch figures.
The large scratch reduction is specific to the dense direct execution plan and
does not include the still-missing exact decode plan.

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
using ASan+UBSan and `LEO2_C6_SANITIZER_MODE="asan-ubsan"`.

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

- Preserve the exact GF8 direct encoder as default-off C7/C8 input.
- Do not describe it as legacy compatible or route production AUTO to it.
- Freeze and serialize a new profile identity before producing persistent data.
- Implement and benchmark an immutable exact erasure plan, including direct
  small-loss repair and a transform alternative, before any codec promotion.
- Add a second CPU family and the full C10 dispatcher matrix.  The current
  side-17 losses must remain visible.
- Revisit Tang--Han or a factored exact transform only as a distinct executable;
  prior inverse-stage operation counts are not end-to-end codec evidence.
