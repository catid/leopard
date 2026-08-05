# Leopard2 bounded direct repair

## Scope and dispatch

Leopard2 retains the specialized low/high LCH decoders and the generic
`O(N log N)` decoder.  Direct repair is available by default in five measured
regions:

- the original bounded region, `2 <= K <= 16`, for one through four missing
  originals and a parent systematic dimension no larger than 256;
- a dual-route GF8/AVX2 legacy-high region with `5 <= K <= 16`,
  `5 <= R <= 8`, and five through eight missing originals, except full loss
  at `K=7` or `K=8` where the translated transform remains faster;
- a GF8/AVX2 legacy-high region for exactly one missing original when
  `17 <= K <= 128`, `P=T`, and the parent is `N=2P`;
- the separately measured GF8/AVX2 `K=65`, `P=T=128` region for one through
  eight missing originals;
- generalized one-loss GF8/AVX2 repair for unequal legacy-high parents with
  `K > 16`, `R > 1`, and for low-profile parents with `K > 16`, `R > K`,
  provided the parent length is no larger than 256.

The small dual-route region retains both direct coefficients and the required
transform metadata in one immutable reusable plan.  The exact execution byte
count deterministically selects transform below 1024 bytes and source-major
direct repair at 1024 bytes and above.  Because those routes have different
workspace layouts, callers must pass the intended byte count to
`leo2_decode_plan_scratch_size` and use that result for the matching execution;
batch items resolve their routes independently.  The one-shot `leo2_decode`
path accounts for its non-amortized setup cost and uses transform below 4096
bytes, preparing a direct-only plan at 4096 bytes and above.  Its codec-level
scratch query remains conservative and pattern-independent.

`-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=OFF` restores transform-only production
behavior for this loss-five-through-eight region.  The diagnostic
`LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1` and `=2` settings instead force the
output-major and source-major direct controls, respectively.  If direct-term
preparation cannot complete, production safely retains the transform route.

The generalized selector is enabled in normal CMake and checked-in Visual
Studio builds.  `-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF` retains the
same-source transform control for reproducible comparisons.

Forced generic, specialized, tiled, or materialized decoding disables direct
repair.  The existing `R=1` XOR and `K=1` copy paths remain earlier dispatch
choices.  Other cases use the locator and transform decoder.  The generalized
one-loss regions are deliberately tied to GF8 and the effective AVX2 operations
table.  Outside the qualified small dual-route and `K=65` regions, they are not
a claim that direct repair wins for multiple losses, GF16, GFNI, AVX-512,
SSSE3, or scalar execution.  The historical global
`LEO2_GFNI_VARIANT` diagnostic also leaves these selectors disabled even when
the CMake option is enabled.

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

For the promoted `L=1` region, the generic matrix expression reduces exactly
to one row normalization.  Codec setup evaluates `log Z_S(q)` once per aligned
recovery block and stores its public `R`-entry prefix.  Plan setup then forms
the selected row entirely in logarithmic form:

    log G(q,j) = log Z_S(q) + log w_j - log(q xor x_j)
    log inverse(G(q,m)) = -log G(q,m)

for missing original `m`.  The recovery coefficient is the inverse above and
each surviving-original coefficient is
`log G(q,j) - log G(q,m)`.  This removes the temporary field-element row,
generic 1-by-1 elimination, exponent-table conversion, and repeated vanishing
product.  The generic solver remains the fallback and correctness oracle.

## Execution and memory

Immutable codecs store the bounded barycentric weights.  Direct-capable decode
plans store the sorted missing originals, per-output term offsets, tagged source
indices, and fixed-multiplier logarithms.  A qualified reusable small dual plan
also retains its locator and transform schedules; both representations remain
immutable and safe for concurrent execution.  When direct repair is selected,
execution writes each missing original directly to its caller-provided disjoint
output and uses scratch only for the existing overlap/range validation array.
When the byte-aware resolver selects transform, the ordinary transform
workspace applies.

GF8 fixed multiply and multiply-add accept any positive byte count.  Complete
64-byte tiles use the current SSSE3/AVX2/AVX512 backend.  The measured extended
GF8/AVX2 regions zero-pad a trailing partial tile on the stack and run the same
backend kernel; the original small region retains its scalar field-element
tail.

For generalized one-loss repair, pure AVX2 additionally groups four
fixed-coefficient sources per call for shard sizes 1 through 63 except 7 bytes.
Each group loads the output once, uses XMM or scalar tails without padding or
out-of-range access, and either initializes or XOR-accumulates explicitly.
The remaining zero through three sources retain their original term order.
Seven-byte shards use the mature output-major executor.

The bounded direct-repair dispatcher itself is kept out of line.  Without this
boundary, GCC inlined the dispatcher into the general decode executor and made
unselected paths inherit the layout cost of the four-source callback.  In the
matched control and candidate, `DecodePlanExecuteInternal` is now exactly
4,281 text bytes in both binaries; the separate direct-repair functions are
1,634 and 1,613 bytes.  This is an execution-layout constraint only and does
not change the plan, arithmetic, or wire profile.

For larger generalized one-loss repairs, a separate pure-AVX2 callback groups
two fixed-coefficient sources per output pass.  It uses YMM operations on each
32-byte prefix, integrated XMM operations on 16- and 8-byte tails, and scalar
code only for the final seven or fewer bytes.  Pairing preserves source order,
halves output loads/stores and indirect dispatches, and adds no shard-data
scratch.  Production selects it only for the exact measured profile/count
shapes and byte thresholds:

- high `(K,R)=(192,64),(200,30),(224,32)` at 64 bytes or larger;
- high `(240,16)` at 64 KiB or larger;
- low `(31,200)` at 64 bytes or larger;
- low `(17,31)` and `(127,128)` at 1 MiB or larger.

Low `(32,224)` was rejected rather than assigned a production threshold.  Its
current-source results were non-monotonic: it missed the confidence gate at
8 KiB and 256 KiB and regressed at 16 KiB, despite wins at several larger
sizes.  This shape retains the mature output-major direct executor.

The exact profile/count threshold is derived once while constructing the
immutable codec and cached as one 64-bit value.  Execution therefore performs
only a threshold load and branch; it does not rescan the measured-shape table.
This matters for tiny unselected neighbors: before caching, low `(31,199)` at
64 bytes had a credible 2.8-percent layout/dispatch regression; afterward its
execution geomean was 1.016x, while selected `(31,200)` remained 1.298x
faster.

The low `(127,128)` threshold is intentionally 1 MiB: current-source screens
found only about one to six percent at 64 KiB and one to four percent at
256 KiB, but at least 1.32x in the measured 1-MiB rounds.  Immediate count and
byte neighbors retain the mature executor.  The pair callback is absent from
scalar, SSSE3, AVX-512, GFNI, GF8-disabled, and global-GFNI experiment tables.
It changes neither repair coefficients nor the wire profile.

The four-source callback has the same backend exclusions, so an AVX2-looking
diagnostic table cannot accidentally interpret affine GFNI tables as nibble
tables.

An optional plan-time source-major transpose was also evaluated for the
measured `K=R=65`, two-, four-, and eight-loss AVX2 region.  It removes the
per-execution stack transpose, but reserves 65 20-byte source-row records,
adds 3,012 bytes of object text, and makes plan setup about 1.44x to 1.83x
slower.  In a 254-pair confirmation its best execution result was 1.0443x
with a 95-percent interval of `[1.0398x, 1.0488x]`; the upper bound therefore
misses the five-percent promotion threshold.  Plan-amortized gains were at
most 2.81 percent in the confirmed cells.  The option remains default-off.
Machine-readable negative evidence is
`experiments/leopard2/direct_repair/results/avx2_source_plan_negative_20260730.json`.

GF16 uses optimized complete-tile kernels and handles an even compact tail as
`q` low bytes followed by `q` high bytes.  Odd GF16 byte counts remain
unsupported.  No-loss, `R=1`, `K=1`, specialized LCH, and generic fallback
semantics are unchanged.

## Correctness evidence

The production test suite now checks:

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
  the `K=65` eight-loss exception, and the small dual-route count, loss, backend,
  1024-byte reusable-plan, and 4096-byte one-shot boundaries; and
- targeted generalized one-loss execution
  through high-profile
  `(K,R)=(240,16)` and low-profile `(31,200)`, including first and last
  systematic coordinates, first and last transmitted parity rows, arbitrary
  tails, and shards through 64 KiB plus one byte, with every result compared
  to the independent/generic oracle; and
- the four-source callback at every selected boundary from 1
  through 63 bytes,
  the 7- and 64-byte fallbacks, all four possible term-count remainders,
  unaligned inputs and outputs, guard bytes, repeated source pointers, every
  multiplier and input byte, and explicit callback-count attribution; and
- the two-source callback at zero through 257 bytes around every vector
  boundary, all 65,025 live multiplier-log pairs, coefficient suppression and
  identity, replace and accumulate modes, unaligned guard bytes, repeated
  source pointers, null pointers for a zero-byte operation, and exact
  integration callback counts at every production threshold.

The 2026-07-30 pure-AVX2 pair-fusion screen used frozen, hashed binaries pinned
to CPU 4 with its SMT sibling idle.  All seven selected threshold cells won:
execution geomean speedups ranged from 1.0831x to 1.5564x, while
plan-amortized geomean speedups ranged from 1.0831x to 1.4647x.  The weakest
execution and amortized 95-percent lower bounds were 1.0790x.  Eight unselected
count/byte neighbors had no geomean regression larger than 1.0 percent.  Every
workload digest matched.
The AVX2 object text grew by 3,464 bytes (3.57 percent), and disassembly found
no EVEX prefix.  Machine-readable evidence is
`experiments/leopard2/direct_repair/results/avx2_pair_production_20260730.json`.
Focused Clang 18 ASan+UBSan+LSan backend/API tests passed serially under a
3-GiB cgroup ceiling; their largest recorded RSS was 1,559,312 KiB.  A
GF16-only production compile also passed, confirming that the optional Ops
member does not break reduced-field builds.

The four-source tiny-shard promotion used the same pinned CPU and frozen-binary
discipline.  Across eight high/low shapes, five representative selected byte
sizes, and plan reuse 1, 8, and 64, all 120 selected execution cells won.
Execution geomeans ranged from 1.0723x to 2.5058x and the weakest 95-percent
lower bound was 1.0595x.  Every plan-amortized geomean also won, from 1.0215x
to 2.3906x.  Six reuse-one cells did not establish a five-percent
plan-amortized lower bound because unchanged plan construction dominates a
one- to six-byte execution; none regressed, and every execution interval
cleared the five-percent gate.

After two high-confidence confirmations, 102 unselected byte/count neighbors
had no raw or credible regression greater than two percent.  The prior
two-source milestone also retained at least a 1.0898x lower bound in all seven
selected cells.  Release, option-off, Clang 18 ASan+UBSan+LSan, context-backend,
and GF16-only gates passed under one- to three-GiB cgroup ceilings; the largest
observed RSS was 1,559,272 KiB.  Machine-readable evidence is
`experiments/leopard2/direct_repair/results/avx2_four_tiny_production_20260730.json`.

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

## Equal-rounded multi-loss promotion

The AVX2 source-major executor is also selected for two through eight missing
originals when `17 <= K <= 128`, `ceil_pow2(K) == ceil_pow2(R)`, and the
separately calibrated `K=65` family does not apply.  The final candidate is
commit `00191aff90d8b20b547fa28b4693e3d7b6b4ebcf`, tree
`2db2b9eafaebd24578002788bd2f12bbc2e5bc6e`.  Candidate and initialized-data
control binaries have identical executable sections; only the immutable
selector datum differs.

The authoritative CPU-13 campaign compared that candidate with both its
same-source transform control and exact Leopard1 `main` commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  It covered 47 target cells and
five selector neighbors in 906 accepted processes.  One three-round
33-byte cell was inconclusive against Leopard1, so a predeclared nine-round
holdout added 54 accepted processes and established a 1.1302x point speedup
with a 1.1041x lower confidence bound.  All original, parity, and recovery
digests matched.

After replacing that cell with its holdout, the same-source transform over
candidate execution geomean is 3.9574x and its weakest lower confidence bound
is 1.3724x.  Exact Leopard1 over candidate execution has a 5.0310x geomean
and a 1.1041x weakest lower bound.  At 64 KiB, the exact-main geomean is
7.5633x and the weakest lower bound is 2.8748x.  At 1 KiB they are 6.4154x
and 2.8550x.  The five unselected neighbors have no credible regression over
two percent.

For logical byte counts 1 through 65, exact Leopard1 is invoked on its required
64-byte-rounded storage while logical digests remain equal.  Candidate
execution wins all 14 tail cells, with a 2.0651x exact-main geomean and a
1.1041x weakest lower bound.  A brand-new plan does not yet win those tiny
cells: its first-use geomean is 0.6228x relative to Leopard1 because about four
microseconds of coefficient/scatter setup dominates.  The conservative
observed reuse crossover is two through fifteen executions.  This limitation
is explicit in `leopard-79h.38.5.10.44`; it does not negate the reusable-plan
execution promotion.

The benchmark retains target CPU and SMT-sibling scheduler accounting,
interrupt deltas, and a `wait4` cross-check for every process.  The effective
bounds are `max(20 us, 50 ppm)` for unexplained target runtime and
`max(50 us, 200 ppm)` for sibling runtime.  These cap possible contamination
at 0.66 percent of the shortest retained invocation and 0.02 percent of long
invocations, below the two-percent neighbor gate, while rejecting the
retained scheduler outliers.  Both raw bundles independently verify.
Machine-readable identities, hashes, confidence intervals, correctness gates,
and reproduction commands are in
`experiments/leopard2/direct_repair/results/equal_rounded_avx2_authoritative_20260731.json`.

## Small-code dual-plan setup and fallback schedule

Reusable legacy-high GF8/AVX2 plans in the qualified `K=5..16`, `R=5..8`,
loss-five-through-eight region own both a transform fallback for tiny shards
and direct coefficients for larger shards.  The initial dual implementation
constructed those coefficients a second time from generator rows and an
`L`-by-`L` inverse even though its transform locator already contains the
needed interpolation factors.  For missing execution coordinate `x` and
selected survivor `s`, Leopard2 now emits the direct multiplier

    Lambda(s) / ((x + s) Lambda'(x)).

The complete derivation and locator-log convention are in
`docs/leopard2_math_and_sources.md`.  A 214-plan pattern sweep compares the
new rows with the retained matrix solver; seven directed native/translated,
punctured, shortened, missing-parity, virtual-erasure, and explicit high-profile
`R > K` cases additionally expose 444 coefficients through source-coordinate
unit inputs against the independent direct systematic generator oracle.
All four locator/fallback selector combinations pass the focused Release gate;
the promoted locator-on/fallback-off configuration also passes focused
ASan+UBSan and TSan execution, including concurrent immutable-plan reuse.

Pinned directional measurements reduced plan setup by 1.58x at K5/R5 and by
3.35--3.47x at K16/R8.  The formerly losing 1024-byte/reuse-one cells moved to
1.259x and 1.278x versus exact Leopard1.  The 1023-byte/reuse-eight K16/R8
point moved from a regression to a 1.048x point win.  These measurements use
immutable binaries and the canonical campaign lock, but still await the
external affinity-supervisor seal required for final authoritative status.

An experimental selector can omit pruned transform graphs where the regular
Algorithm 4/5 fallback appeared faster.  A 400-cell boundary screen covered
bytes `1,63,64,65,255,256,257,1023`; all 352 newly affected cells improved,
with geometric setup, execution, and amortized speedups of 2.778x, 1.198x, and
1.348x.  However, the 48 intended-inert controls showed 18 execution point
regressions, seven exceeding 2% and a worst point of 0.949x.  The two binaries
had different hot-function text layouts and the sequential 11-sample runs had
no confidence intervals, so those controls expose an unresolved measurement
confound.  The selector is therefore disabled by default pending an isolated
same-binary or section-stable ABBA replay with confidence intervals, including
explicit high-profile `R > K` shapes.

Machine-readable directional results are in
`experiments/leopard2/direct_repair/results/small_dual_locator_terms_directional_20260805.json`
and
`experiments/leopard2/direct_repair/results/small_dual_regular_fallback_directional_20260805.json`.
An attempted specialized-only metadata reduction was rejected because its
amortized result was neutral to negative; that result is retained in
`experiments/leopard2/direct_repair/results/small_dual_specialized_metadata_negative_20260805.json`.

## Scratch-resident tiny native-high one-shot path

For GF8/AVX2 legacy-high `N=32`, `T=8` codecs with `K=9..16`, `R=5..8`,
losses `5..R`, and at most 256 bytes per shard, the public one-shot decoder now
constructs its deterministic selection, locator, block prefixes, and output
descriptors in caller scratch and invokes the tiled Algorithm 5 executor
without allocating a heap plan.  Full 64-byte prefixes reveal directly into
the requested outputs; a ragged tail is staged in existing scratch and
gathered.  The reusable plan API and every ineligible one-shot call retain the
ordinary immutable-plan path.

Across 768 allocation-audited public cases and 168 independent-oracle cases,
the new path reproduced reusable Algorithm 5 and direct systematic recovery,
with zero eligible hot-path allocations.  Focused Release, ASan+UBSan, TSan,
and strict-warning gates pass.  A zero-sibling-jiffy same-source ABBA measured
1.342x at `K=16,R=8,L=8,B=64`; its 257-byte inert control was neutral.  A
separate pure-AVX2 comparison against exact Leopard1 measured 1.188x at the
original 64-byte target, with a 95% interval of 1.175x--1.202x.

The same exact-main replay found a 0.829x ratio at 65 bytes.  That is an
adjacent absolute crossover gap rather than a candidate regression: the raw
path still improves 65-byte one-shot time by 1.251x over its heap-plan control,
but its extra ragged transform tile costs more than Leopard1.  Direct repair of
only the ragged tail remains the leading follow-up.

See `experiments/leopard2/optimization_log/27-raw-native-high-one-shot.md` and
`experiments/leopard2/direct_repair/results/raw_native_high_one_shot_20260805.json`.
