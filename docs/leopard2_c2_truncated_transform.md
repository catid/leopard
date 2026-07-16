# Leopard2 C2 parent-preserving truncated transforms

Status: scalar correctness checkpoint complete; retain as an experiment and
advance only the compact recursive schedule to a bounded C++ prototype.  This
does not define an exact-length code, does not change a coordinate set, and is
not enabled by the default build.

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

## Disposition and limitations

Promote only a bounded C++ prototype of:

- the compact mask-ordered API contract;
- forward/backward mask propagation using exact butterfly coefficients;
- descriptors for full dyadic subtrees plus ragged boundary pairs; and
- zero/one multiplier specialization with the existing active-coset skew table.

Do not promote this Python executor or infer a dispatcher threshold.  No timing
was run, intentionally, because this checkpoint was developed while independent
repository validation was using the host.  The following work remains before
production consideration:

- translate to GF8/GF16 C++ and compare every available scalar/SIMD backend;
- define a reusable flat schedule and a rigorously bounded scratch allocator;
- add arbitrary byte tails, aliasing, guard-page, sanitizer, fuzz, and immutable
  concurrent-plan tests;
- measure setup separately from execution and compare against C1 pruning,
  padded transforms, and C3 basis-conversion alternatives; and
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
