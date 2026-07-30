# Leopard2 bounded direct repair

## Scope and dispatch

Leopard2 retains the specialized low/high LCH decoders and the generic
`O(N log N)` decoder.  An immutable direct-repair plan is selected by default
in three measured regions:

- the original bounded region, `2 <= K <= 16`, for one through four missing
  originals and a parent systematic dimension no larger than 256;
- a GF8/AVX2 legacy-high region for exactly one missing original when
  `17 <= K <= 128`, `P=T`, and the parent is `N=2P`; and
- the separately measured GF8/AVX2 `K=65`, `P=T=128` region for one through
  eight missing originals.

Two broader one-loss candidates are available only when configured with
`-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=ON`:

- unequal legacy-high GF8/AVX2 parents with `K > 16`, `R > 1`, and parent
  length no larger than 256; and
- low-profile GF8/AVX2 parents with `K > 16`, `R > K`, and parent length no
  larger than 256.

Forced generic, specialized, tiled, or materialized decoding disables direct
repair.  The existing `R=1` XOR and `K=1` copy paths remain earlier dispatch
choices.  Other cases use the locator and transform decoder.  The added
experimental one-loss regions are deliberately tied to GF8 and the effective
AVX2 operations table.  Diagnostic screens are not promotion evidence; the
ordinary build keeps its established selector until the frozen same-source,
exact-main, reuse, large-shard, and neighboring-shape gates pass.  They are not
a claim that direct repair wins for multiple losses, GF16, GFNI, AVX-512,
SSSE3, or scalar execution.  The historical global `LEO2_GFNI_VARIANT`
diagnostic also leaves these selectors disabled even if the experiment option
is set.

## Generator-row derivation

This is a new derivation from the declared wire-profile coordinate sets.  It
does not call or copy the test-only direct oracle.

Let `S` be every parent systematic coordinate, including shortened zero
coordinates.  Thus `S=[T,N)` for the legacy high profile and `S=[0,P)` for the
low profile.  Coordinates are legacy Leopard Cantor-coordinate field elements,
so subtraction is bitwise XOR.  Define

    Z_S(q) = product over s in S of (q xor s)

and precompute one codec-time barycentric weight per public original coordinate
`x_j`:

    w_j = inverse(product over s in S, s != x_j of (x_j xor s)).

For a transmitted parity coordinate `q`, the systematic generator coefficient
is the Lagrange cardinal polynomial evaluated at `q`:

    G(q,j) = Z_S(q) * inverse(q xor x_j) * w_j.

This construction includes every shortened parent coordinate in `S`, while
only public originals need stored weights.  Puncturing changes which parity
rows are available; it does not change `G`.

For `L` missing originals, plan setup selects the lowest-index `L` available
parity shards and forms only the `L x L` matrix `A` containing their generator
coefficients at the missing columns.  Gaussian elimination in the legacy
field produces `A^-1`.  For each missing output `o`, plan setup folds the byte
equation into fixed terms:

    x_o = sum_e A^-1[o,e] parity_e
          + sum over surviving j of
            (sum_e A^-1[o,e] G(e,j)) original_j.

The plus sign is also subtraction in characteristic two.  Zero coefficients
are removed.  Nonzero coefficients are stored in Leopard logarithmic form so
execution does no field setup, matrix work, locator work, or allocation.
Coefficient one is ordered first when present and is specialized to copy/XOR.

When the complete parent systematic set is an aligned additive coset
`a + V_d`, XOR by any `x` in the coset permutes the other coordinates onto
`V_d` without zero.  Therefore every public systematic coordinate has the
same denominator

    c_d = product over nonzero u in V_d of u.

Codec setup computes this product and its inverse once, then fills the public
weight vector.  The general per-coordinate calculation remains the fallback
for a systematic interval that is not an aligned power-of-two coset.  This
reduces setup from `O(K D)` to `O(D + K)` for every low-profile bounded direct
encoder/repair codec and for equal-rounded `P=T` legacy-high direct repair;
it does not change generator coefficients or wire bytes.

The unequal legacy-high parent has another exact structure.  Its complete
systematic coordinate set is

    S = [T,N) = V_n \ V_t.

Let `s_n`, `s_t`, and `Z_S` be the vanishing polynomials of `V_n`, `V_t`, and
`S`.  Since `s_n = s_t Z_S`, differentiating at any `x` in `S` gives

    c_n = s_n'(x) = s_t(x) Z_S'(x),

because `Z_S(x)=0` and the derivative `c_n` of a complete additive subspace is
constant.  The needed barycentric weight is therefore

    w_x = inverse(Z_S'(x)) = s_t(x) / c_n.

The additive polynomial `s_t` is constant on each aligned `T`-coordinate
coset.  Production computes `c_n` once, computes one `s_t` value per occupied
message block, and fills the public prefix of that block.  This is `O(N+D)`
field work rather than `O(K D)`.  It applies equally to shortened public
messages because the fixed shortened suffix remains part of `S`; puncturing
does not alter `S`.  The independent direct-product test enumerates every
power-of-two `N<=256`, every `T<N`, and every `x` in `S` in GF(2^4), legacy
GF8, and legacy GF16, comparing this quotient identity with the original
denominator product.

The MDS property makes any valid selected `L x L` parity submatrix invertible.
The implementation nevertheless treats a zero pivot as a plan-preparation
failure and safely falls back to the transform decoder.

## Execution and memory

Immutable codecs store the bounded barycentric weights.  Immutable decode plans
store the sorted missing originals, per-output term offsets, tagged source
indices, and fixed-multiplier logarithms.  Plan execution writes each missing
original directly to its caller-provided disjoint output and uses scratch only
for the existing overlap/range validation array.

GF8 fixed multiply and multiply-add accept any positive byte count.  Complete
64-byte tiles use the current SSSE3/AVX2/AVX512 backend.  The measured extended
GF8/AVX2 regions zero-pad a trailing partial tile on the stack and run the same
backend kernel; the original small region retains its scalar field-element
tail.

With the generalized one-loss experiment enabled, pure AVX2 additionally
groups four
fixed-coefficient sources per call for shard sizes 1 through 63 except 7 bytes.
Each group loads the output once, uses XMM or scalar tails without padding or
out-of-range access, and either initializes or XOR-accumulates explicitly.
The remaining zero through three sources retain their original term order.
Seven-byte and 64-byte-or-larger shards use the mature output-major executor.
The callback is absent from scalar, SSSE3, AVX-512, GFNI, GF8-disabled, and
global-GFNI experiment tables, so an AVX2-looking diagnostic table cannot
accidentally interpret affine GFNI tables as nibble tables.

GF16 uses optimized complete-tile kernels and handles an even compact tail as
`q` low bytes followed by `q` high bytes.  Odd GF16 byte counts remain
unsupported.  No-loss, `R=1`, `K=1`, specialized LCH, and generic fallback
semantics are unchanged.

## Correctness evidence

The production API test now checks:

- every GF8 element pair and every nonzero GF8 inverse against the independent
  polynomial-basis/Cantor-coordinate oracle;
- deterministic GF16 multiply/inverse samples against that oracle;
- fixed multiply and multiply-add over GF8 sizes `1,63,64,65,257`;
- fixed multiply and multiply-add over compact GF16 sizes
  `2,34,64,66,1026`;
- high and low direct repair in GF8 and GF16, including skipped low-index parity
  rows, repeated immutable plan execution, and comparison with the forced
  generic decoder;
- an independent direct-algebra check that enumerates aligned cosets through
  dimension 256 in GF8 and GF16 and verifies every barycentric denominator
  against `product(V_d \\ {0})`; and
- an independent complement-subspace check for every supported power-of-two
  `N,T` through the field order—`N<=16` in GF(2^4), and `N<=256` in GF8 and
  GF16—verifying `s_t(x)/c_n` against direct products for every legacy-high
  parent systematic coordinate; and
- exact dispatch boundaries for the bounded region, the equal-rounded one-loss
  region, unequal legacy-high and redundancy-dominant low one-loss regions,
  and the `K=65` eight-loss exception; and
- in an experiment-enabled build, targeted generalized one-loss execution
  through high-profile
  `(K,R)=(240,16)` and low-profile `(31,200)`, including first and last
  systematic coordinates, first and last transmitted parity rows, arbitrary
  tails, and shards through 64 KiB plus one byte, with every result compared
  to the independent/generic oracle; and
- in that build, the four-source callback at every selected boundary from 1
  through 63 bytes,
  the 7- and 64-byte fallbacks, all four possible term-count remainders,
  unaligned inputs and outputs, guard bytes, repeated source pointers, every
  multiplier and input byte, and explicit callback-count attribution.

Validation commands completed on 2026-07-16:

    OMP_NUM_THREADS=1 ./build/release/leopard2_api_test
    ./build/release/leopard2_random_test --seed 903176249 --cases 2048 --threads 32
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      ctest --test-dir build/direct-repair-asan -j 32 --output-on-failure

The release random run completed 2,048 patterns, 4,096 repeated plan
executions, 24,962 recovered shards, and 4,452 generic comparisons.  The
ASan+UBSan build completed all 11 configured tests without a finding.  GCC 13.3
compiled `leopard2.cpp` with `-Wall -Wextra -Wpedantic -Wconversion
-Wsign-conversion -Wshadow -Werror`; both FF translation units passed
`-Wall -Wextra -Wpedantic -Werror` after suppressing the repository's existing
`FFTSkew - 1` array-bounds diagnostic.

A 2026-07-21 follow-up rebuilt and passed `leopard2_api`,
`leopard2_direct_oracle`, and `leopard2_direct_repair` in both Release and
Clang 18 ASan+UBSan configurations after extending the selector.  This focused
optimization pass intentionally did not run a new broad fuzz campaign.

## Pinned crossover evidence

The benchmark process was pinned to CPU 0, one codec thread, AVX2, GCC 13.3,
with median/MAD reporting.  Only CPUs 0-31 were available in the worker's
affinity set, so these are intentionally isolated single-core measurements,
not 128-core scaling claims.  The comparator is the retained forced-generic
transform decoder for the exact same profile, field, counts, losses, and shard
bytes.

| Cell | Direct median (us) | Generic median (us) | Speedup | Direct / generic scratch |
| --- | ---: | ---: | ---: | ---: |
| high GF8, K=16 R=8 L=1, 4 KiB | 1.241 | 9.042 | 7.28x | 640 / 230528 |
| high GF8, K=16 R=8 L=4, 1 MiB | 997.048 | 3060.098 | 3.07x | 640 / 58721408 |
| low GF8, K=8 R=248 L=1, 4 KiB | 1.564 | 70.260 | 44.93x | 4224 / 2105472 |
| low GF8, K=8 R=248 L=4, 1 MiB | 510.456 | 58523.025 | 114.65x | 4224 / 536879232 |
| low GF16, K=16 R=240 L=4, 66 B | 1.455 | 10.274 | 7.06x | 4352 / 73984 |
| low GF16, K=16 R=240 L=4, 64 KiB | 72.606 | 1311.237 | 18.06x | 4352 / 33562880 |

Machine-readable results are in the ignored directory
`.research/leopard2/direct_repair_bench/`.  Direct plan setup medians ranged
from 0.28 to 1.45 microseconds in these cells.  The forced-generic setup ranged
from 0.19 to 36.05 microseconds; setup and execution are reported separately.

## Equal-rounded one-loss promotion

A 30-worker directional matrix compared direct repair with the same-binary
forced translated-low transform over 1,950 cells: the requested 13 `K` values
from 17 through 128, three equal-rounded `R` boundaries per `K`, loss counts
1, 4, 8, 12, and 16, ten shard sizes from 1 byte through 1 MiB, and reuse
counts 1, 8, 64, and 1024.  Every data, parity, and repaired-output digest
matched.  One loss won all 390 cells, with minimum execution and reuse-one
speedups of 1.369x and 1.577x respectively.  Multi-loss results were
fragmented, so they were not generalized beyond the existing `K=65`
exception.

The promotion gate then compared against the exact Leopard1 `main` codec at
commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, not merely against another
Leopard2 path.  Thirty representative cells covered ten `(K,R)` pairs and
64-byte, 4-KiB, and 1-MiB shards.  Each cell used three pinned A/B/B/A rounds,
15 measured iterations per process, four warmups, identical workload digests,
and candidate time equal to execution plus plan setup divided by reuse.  A
cell was discarded and rerun if the reserved SMT sibling accumulated even one
non-idle jiffy.  All 30 accepted cells had zero such activity, selected the
direct path, and matched all digests.  Mean paired Leopard1-over-Leopard2
speedups ranged from 6.128x to 58.769x; even the weakest cell's 95-percent
interval was `[5.858x, 6.399x]`.

The measured executable and static archive were then rebuilt from clean commit
`b665ed0fb3c8d2479eec6129fb43672a1a328630` and matched byte for byte.  The
full per-cell medians, confidence intervals, source and binary hashes,
raw-bundle manifests, exclusions, directional summary, and selector definition
are in
`experiments/leopard2/direct_repair/results/one_loss_equal_rounded_exact_main.json`.
