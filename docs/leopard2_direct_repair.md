# Leopard2 bounded direct repair

## Scope and dispatch

Leopard2 retains the specialized low/high LCH decoders and the generic
`O(N log N)` decoder.  An immutable direct-repair plan is selected only when
all of the following deterministic conditions hold:

- generic decoding was not explicitly forced;
- `2 <= K <= 16`;
- one through four original shards are missing;
- the parent systematic dimension is at most 256; and
- neither the existing `R=1` XOR path nor the existing `K=1` copy path applies.

Every other case falls back to the existing locator and transform decoder.
The deliberately small production region was chosen after measuring the cells
below and keeps the initial derivation, codec setup, and review surface bounded.
It is not a claim that `K=17` is a mathematical or performance crossover.
Expanding the region requires additional pinned measurements and randomized
coverage, not a wire-format change.

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
64-byte tiles use the current SSSE3/AVX2/AVX512 backend and a trailing partial tile uses
scalar field elements.  GF16 uses the same optimized complete-tile kernels and
handles an even compact tail as `q` low bytes followed by `q` high bytes.  Odd
GF16 byte counts remain unsupported.  No-loss, `R=1`, `K=1`, specialized LCH,
and generic fallback semantics are unchanged.

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
- the accepted `K=16` and fallback `K=17` scratch/dispatch boundary.

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
