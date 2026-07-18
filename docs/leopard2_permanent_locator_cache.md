# Leopard2 permanent-locator cache disposition

This checkpoint resolves the apparent dense permanent-locator cache gap without
adding a dominated production cache.  The production sparse cache is complete
for the current deterministic plan construction; neither a coordinate-domain
nor a transform-domain dense cache removes a Walsh transform.

## Linearity and the exact dynamic count

For permanent erasures `P`, pattern-specific real or virtual erasures `D`, and
active parent `V_n`, the logarithmic locator is

    lambda_E[x] = sum over e in E of ell[x xor e] modulo q,
    ell[0] = 0,
    q = 2^m - 1.

`P` and `D` are disjoint in plan storage.  The zero logarithm omits the self
factor at a locator root and therefore preserves the derivative value used by
the reveal factor.  Directly from the sum,

    lambda_(P union D) = lambda_P + lambda_D modulo q.

The current transmitted code has `K+R` coordinates.  Plan construction retains
every surviving original and then the lowest-index surviving parity coordinates
until exactly `K` public coordinates remain.  Every other received parity is a
virtual erasure.  Therefore exactly `R` transmitted coordinates are missing or
virtual, independent of the real loss pattern.  Permanent shortened or
punctured parent coordinates are outside that transmitted set, so

    |D| = (K + R) - K = R.

The codec's `IsDirectLocatorPreferred(N, recovery_count)` decision consequently
predicts the plan's dynamic direct/Walsh decision exactly; it is not a loss-count
heuristic.  A debug-only production assertion and an exhaustive replicated
selection test guard this invariant.  The test enumerates 9,647 valid presence
patterns for every `1 <= K,R <= 6`, including no loss, mixed loss, maximum loss,
and surplus received parity.

## Dense candidates

Let `H` be the unnormalized `N`-point Walsh transform and let

    C = N^-1 H(ell),
    N^-1 = 2^m / N modulo q.

The production union-mask path computes

    H(C pointwise-multiplied-by H(1_(P union D))).

Two strongest reusable alternatives were implemented as executable test
prototypes:

1. Coordinate cache: cache `lambda_P`, compute the complete active-Walsh
   `lambda_D`, then add the two `N`-entry outputs.
2. Transform cache: cache `C * H(1_P)`, compute `H(1_D)`, fuse an addition of
   the cached vector into the pointwise multiplication, then run the second
   Walsh transform.

For `b=log2(N)`, one Walsh transform produces `N*b` modular add/subtract
results.  Setup of the permanent cache is excluded from all three reusable-plan
rows, just as codec setup is excluded from plan setup.

| Reusable plan construction | Walsh modular results | Pointwise products | Extra modular additions | Persistent codec cache |
| --- | ---: | ---: | ---: | ---: |
| Union mask (production) | `2*N*b` | `N` | `0` | `0` |
| Coordinate cache | `2*N*b` | `N` | `N` | `N` field symbols |
| Transform cache | `2*N*b` | `N` | `N` | `N` field symbols |

Both dense caches retain both transforms and add an `O(N)` modular pass.  The
transform form merely moves that pass between the transforms.  The union-mask
initialization already writes the caller-owned `N`-entry locator output and
needs no separate temporary.  Thus the dense candidates increase arithmetic
and persistent storage without reducing asymptotic or pass count.  No dense
cache is promoted.  Timing was intentionally not used to override this strict
operation/pass dominance while unrelated sanitizer/fuzz jobs shared the host.

The sparse direct path is different.  With `p=|P|` and `d=|D|`, a fresh direct
union costs `(p+d)*(N-1)` modular additions.  The existing codec cache copies
`N` field symbols and performs `d*(N-1)` additions, saving exactly
`p*(N-1)` additions.  Because current plans always have `d=R`, codec setup
materializes this cache precisely when the measured direct cutoff can use it.
Unconditionally storing another `N` symbols would not expand the current
production reuse region and would penalize every punctured dense codec.

## Independent correctness evidence

`leopard2_locator_test` now validates the decomposition independently of the
production Walsh implementation:

- every disjoint permanent/dynamic assignment through `N=8` in test-only
  `GF(2^4)` (`3^2 + 3^4 + 3^8` assignments), plus 4,096 deterministic `N=16`
  decompositions;
- every `GF(2^4)` union erasure subset through `N=16` against direct products;
- add, multi-add, remove, swap, dynamic-only, union-mask, puncture-like, and
  dense transitions at GF8 `N=16,64,256` and GF16
  `N=32,256,2,048,65,536`;
- production direct, active-Walsh, dispatched permanent-base, coordinate-cache,
  and transform-cache results against the direct product;
- the active-parent `2^(m-n)` normalization in both cache prototypes;
- 256 concurrent preparations sharing one immutable permanent base; and
- structural counts at every power of two from 2 through 65,536.

The focused Release result is:

    leopard2 locator tests passed: locator_entries_compared=9894532
      permanent_cache_structural_cases=16
      sparse_additions_saved=715811501
      selection_patterns=9647 concurrent_preparations=256

`sparse_additions_saved` is a deterministic aggregate over the structural test's
synthetic `(N,p,d)` cells, not a throughput claim.  It makes regressions in the
derived operation model visible.

The focused gate also passed:

- GCC 13.3 Release with `-Wall -Wextra -Werror -Wpedantic`: locator, random
  differential, and public API contract, 3/3;
- GCC 13.3 Debug with `DEBUG` enabled, exercising the production exact-`R`
  assertion: the same 3/3;
- Clang 18.1 ASan+UBSan with leak detection and non-recovering undefined
  behavior: the same 3/3; and
- Clang 18.1 TSan with `ignore_noninstrumented_modules=1`: the focused locator
  test, including 256 immutable-base concurrent preparations.

Representative reproduction:

    cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$(nproc)" --target leopard2_locator_test
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
      build/release/leopard2_locator_test

If a future decoder retains more than `K` public coordinates, changes virtual
selection, or finds a kernel that eliminates a transform, it must revisit the
exact-`R` proof and this negative dense result rather than silently reusing the
present cutoff.
