# Leopard2 product-tree and epsilon-assisted locator algebra checkpoint

This document records the non-timing checkpoint for `leopard-79h.29`.  It does
not change the production locator dispatcher or its thresholds.  The isolated
prototype is
`experiments/leopard2/locator_construction/locator_product_tree.cpp`; its
deterministic result is
`experiments/leopard2/locator_construction/results/algebra.json`.

## Product-tree construction

For an erasure-coordinate set `E`, the prototype constructs the monic locator

    Lambda(x) = product over e in E of (x + omega_e).

Subtraction is addition in characteristic two.  Leopard's Cantor coordinate
order has `omega_i + omega_j = omega_(i xor j)`, so an independent direct factor
at coordinate `i` is

    product over e in E, e != i of omega_(i xor e).

At a surviving coordinate this is `Lambda(omega_i)`.  At an erased coordinate,
`Lambda(omega_i)` is zero and the same product is `Lambda'(omega_i)`.  The
prototype obtains the latter independently by differentiating the completed
dense polynomial: only odd-degree coefficients survive in characteristic two.
This is a new, elementary derivation implemented directly from the product
definition; it does not reuse Leopard's logarithmic convolution.

The factors `(x + omega_e)` are combined in a balanced binary tree in the
monomial basis.  This checkpoint deliberately uses schoolbook dense polynomial
multiplication and Horner evaluation.  That makes the algebra obvious and gives
a useful lower-risk oracle, but it is not a claim that this form should be used
in production.

Five independently reached representations are cross-checked at every active
coordinate:

1. the product-tree polynomial, using its formal derivative at roots;
2. a direct field-element product over all other erasures;
3. the production direct logarithmic construction;
4. the production active-parent Walsh construction; and
5. the epsilon-assisted LCH construction described below.

The test-only GF(2^4) implementation has its own polynomial (`0x13`), log and
exponent tables, multiplication, direct-product evaluator, and normalized
active Walsh convolution.  GF8 and GF16 use only Leopard's scalar field
primitive for product-tree multiplication and compare the result with the two
production locator APIs.

## Epsilon-assisted locator candidate

R13 does not itself construct a polynomial from a known root set.  Appendix A,
Algorithm 8 accepts `epsilon` contiguous evaluations of an already-defined
degree-`< epsilon` polynomial and returns its normalized-LCH coefficients plus
the completed enclosing `q = ceil_pow2(epsilon)` evaluation block.  Algorithm 9
accepts a complementary suffix under the coefficient-zero constraint used by
that paper's arbitrary-parameter encoder.  Neither algorithm accepts roots.

That source boundary does not rule out a composition with direct root products.
The prototype now implements the following clearly labeled new derivation.  Let
`d = |E|` and let `Lambda` be the monic degree-`d` locator.

When `d` is not a power of two, set `epsilon = d + 1` and
`q = ceil_pow2(d)`.  Directly evaluate `Lambda` at the first `epsilon` Cantor
coordinates, then apply R13 Algorithm 8 literally.  The result is the exact
normalized-LCH coefficient vector for `Lambda` and its first completed `q`
evaluations.  Evaluate the same coefficients on every remaining aligned
`q`-coset with a shifted `q`-point LCH transform.

The naive construction would double `q` when `d` is itself a power of two.
Instead, let `d = q` and interpolate only the first `q` locator values.  The
unique degree-`<q` interpolant is

    g(x) = Lambda(x) + s_log2(q)(x),

because `Lambda` and the active-subspace polynomial are both monic of degree
`q`, and `s_log2(q)` vanishes on the first `q` points.  Additivity gives

    s_log2(q)(B + v) = s_log2(q)(B), for v in V_log2(q),

so a shifted transform of `g` is converted back to `Lambda` by XORing one
precomputed constant into every output of that aligned block.  This preserves
`q = d` at the power-of-two boundary.

At erased coordinates, the zero locator evaluation is replaced by the direct
root product `Lambda'(omega_e)`.  Every such value is independently compared
with the formal derivative of the dense product-tree polynomial.  `d=0` is the
constant-one locator.  For `d=N`, `d+1` is outside the active parent and may
exceed the field; the locator is the active subspace polynomial, whose derivative
is one constant.  The prototype computes that constant once as the product of
all nonzero active coordinates and scatters it.  The full-parent special case is
exercised in GF(2^4), full-field GF8, and active-parent GF16 cases.

The scalar LCH butterflies do not call Leopard's transform tables.  Their skew
factors are derived from

    s_0(x) = x
    s_(j+1)(x) = s_j(x) * (s_j(x) + s_j(omega_(2^j)))

and normalized by `s_j(omega_(2^j))`.  Thus the candidate provides an
independent algebraic differential while retaining Leopard's field and Cantor
coordinate representation.

R15 still matters for a future fast monomial/LCH product-tree implementation,
but it does not directly provide a known-root builder.  The current candidate
is a new composition of the root definition, R13 Algorithm 8, and the additive
subspace identity above; it is not attributed to either paper as a published
locator algorithm.

Sources consulted for this boundary:

- R13, Appendix A, Algorithm 8, and its locator/Forney discussion:
  https://arxiv.org/pdf/2207.11079
- R15, Sections 3.3--3.6 on truncated Lagrange/LCH and monomial/LCH basis
  conversion: https://arxiv.org/pdf/1807.07785

Local text extracts are cached as
`.research/leopard2/text/R13-2207.11079.txt` and
`.research/leopard2/text/R15-1807.07785.txt`.

## Deterministic correctness result

The frozen result contains no wall-clock measurements:

| Coverage | Cases |
| --- | ---: |
| GF(2^4), every subset at `N=2,4,8,16` | 65,812 |
| GF8 differential, every active parent through 256 | 57 |
| GF16 differential, every active parent through 65,536 | 113 |
| Total cases | 65,982 |
| Locator representations cross-checked | 9,868,260 entries |
| Derivative-at-erasure checks | 528,588 entries |
| Epsilon candidate cases cross-checked | 65,982 |
| Full-parent epsilon special cases | 20 |

The GF8 sweep includes empty, sparse, quarter, half, and full erasure sets.  The
GF16 sweep uses the same dense sets through `N=256`; larger parents include
deterministic counts through 64 erasures so that the `N=65,536` parent is covered
without turning a quadratic oracle into a timing experiment.

The evidence file SHA-256 at this checkpoint is
`9b7b337aad145cd639ef3001b08b39e1f34b6d7a5ba90e1c9784d9f054145830`.
The same matrix passed a Clang 18 ASan+UBSan build with leak detection and
halt-on-error enabled; that run reproduced the identical result hash.

## Deterministic operation and storage accounting

The accounting distinguishes operation domains rather than pretending they
have equal cost:

- product-tree and product-evaluation counts are field multiplications and
  field additions (XORs);
- direct counts are the equivalent field-product oracle operations and the
  production log-table/modular-add operations;
- Walsh counts are modular additions in the exponent ring plus pointwise
  integer multiply/reductions; and
- the Walsh kernel precomputation is excluded, matching production reuse.

The epsilon accounting separates reusable fixed setup from one erasure-pattern
execution.  `fixed_setup_*` counts the independently generated subspace
normalizers, inverses, and unique skew/subspace constants used by a precompiled
schedule.  Direct prefix products, Algorithm 8, shifted-coset transforms, and
root derivatives are the per-pattern arithmetic used to construct a decode
plan.  They are not shard-byte decode execution.  Field additions are XORs;
fixed multiplications by zero or one and products whose dynamic operand is zero
are specialized away by the scalar oracle.

`logical_scratch_elements` is the two `q`-element blocks needed by a reusable
in-place realization: one coefficient block and one value/work block.  As with
the product-tree storage count, the common `N`-entry locator output is excluded.
The isolated recursive prototype deliberately uses `std::vector` temporaries,
so this is a deterministic target-space model rather than allocator-RSS
evidence.

`peak_tree_elements` is the maximum combined coefficient storage for one
complete tree level and the partially constructed next level.  The common
`N`-entry locator output is not included.  `locator_polynomial_elements` and
`derivative_polynomial_elements` are reported separately.  This is deterministic
storage accounting, not allocator-RSS evidence.

Representative cells are:

| Field | N | E | Product build + evaluation field multiplies | Direct factors/log adds | Walsh modular adds | Walsh pointwise | Peak tree elements |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GF8 | 256 | 7 | 1,881 | 1,785 | 4,096 | 256 | 25 |
| GF8 | 256 | 128 | 58,175 | 32,640 | 4,096 | 256 | 448 |
| GF8 | 256 | 256 | 165,759 | 65,280 | 4,096 | 256 | 896 |
| GF16 | 2,048 | 14 | 29,012 | 28,658 | 45,056 | 2,048 | 49 |
| GF16 | 65,536 | 14 | 917,844 | 917,490 | 2,097,152 | 65,536 | 49 |
| GF16 | 65,536 | 64 | 4,200,799 | 4,194,240 | 2,097,152 | 65,536 | 224 |

The dense-coefficient product tree plus Horner evaluation does not reduce
algebraic work relative to the existing choices.  At sparse counts it closely
tracks direct `N*E` evaluation while using more expensive field multiplication;
at dense GF8 counts it is decisively worse than both direct log accumulation and
the active Walsh convolution.  It is retained as an independent oracle and as a
baseline for a future fast multipoint evaluator, not promoted.

The epsilon candidate has a materially different symbolic region.  Selected
cells below report its per-pattern field arithmetic; fixed setup is shown
separately and is reusable for the same `(field,N,q)` schedule.

| Field | N | E | epsilon | q | Pattern field multiplies | Pattern XORs | Fixed-setup multiplies / XORs / inverses | Scratch elements |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GF8 | 256 | 7 | 8 | 8 | 479 | 785 | 352 / 131 / 3 | 16 |
| GF8 | 256 | 128 | 128 | 128 | 25,373 | 2,689 | 524 / 275 / 8 | 256 |
| GF8 | 256 | 256 | special | 0 | 255 | 0 | 0 / 0 / 0 | 1 |
| GF16 | 2,048 | 14 | 15 | 16 | 4,360 | 8,235 | 3,330 / 1,414 / 4 | 32 |
| GF16 | 65,536 | 14 | 15 | 16 | 127,366 | 262,187 | 106,498 / 45,062 / 4 | 32 |
| GF16 | 65,536 | 64 | 64 | 64 | 204,798 | 459,009 | 131,079 / 64,527 / 7 | 128 |

These counts make the candidate promising enough to benchmark: for example,
the `N=65,536,E=14` path replaces roughly 917,000 direct factors with about
127,000 field multiplications plus 262,000 XORs.  They do not establish a win.
GF16 multiplication is much more expensive than a production log-table add,
and the fixed schedule/table footprint is substantial.  No dispatcher or
production code selects this path at this checkpoint.

## Reproduction and remaining gates

From a configured Release tree:

    cmake --build build/release -j "$(nproc)" --target leopard
    c++ -std=c++11 -O2 -Wall -Wextra -Wpedantic -I. \
      experiments/leopard2/locator_construction/locator_product_tree.cpp \
      build/release/libleopard.a -fopenmp -pthread \
      -o build/release/leopard2_locator_product_tree
    build/release/leopard2_locator_product_tree --output \
      experiments/leopard2/locator_construction/results/algebra.json
    jq -e '.schema == "leopard2.locator-construction.algebra.v2" and
      .timing_evidence == false and
      .epsilon_candidate == "implemented_nonproduction_correctness_checkpoint_pending_timing" and
      .summary.gf2_4_exhaustive_cases == 65812 and
      (.cases | length) == .summary.accounting_cases' \
      experiments/leopard2/locator_construction/results/algebra.json

The sanitizer replay used an independently instrumented library and prototype:

    cmake -S . -B build/locator-asan -G Ninja \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_COMPILER=clang++-18 \
      -DCMAKE_CXX_FLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined" \
      -DLEO2_BUILD_BENCHMARKS=OFF \
      -DLEO2_ENABLE_CUDA=OFF
    cmake --build build/locator-asan -j "$(nproc)" \
      --target leopard leopard2_locator_test
    clang++-18 -std=c++11 -O1 -g -fno-omit-frame-pointer \
      -fsanitize=address,undefined -Wall -Wextra -Wpedantic -I. \
      experiments/leopard2/locator_construction/locator_product_tree.cpp \
      build/locator-asan/libleopard.a -fopenmp -pthread \
      -o build/locator-asan/leopard2_locator_product_tree
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      build/locator-asan/leopard2_locator_product_tree --output \
      /tmp/leopard2-locator-algebra-asan.json
    cmp experiments/leopard2/locator_construction/results/algebra.json \
      /tmp/leopard2-locator-algebra-asan.json
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      ctest --test-dir build/locator-asan -R '^leopard2_locator$' \
      --output-on-failure

This checkpoint intentionally leaves the Bead open.  Completion still requires
isolated setup/reuse timing against production direct and active Walsh paths,
memory measurements beyond deterministic element counts, and deterministic
promotion thresholds.  The epsilon-assisted construction is now an executable
correctness candidate, but it remains isolated and default-off until those
measurements establish a useful region and a dispatcher can avoid regressions.
