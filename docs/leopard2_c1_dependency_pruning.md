# Leopard2 C1 parent-preserving dependency pruning

Status: scalar experiment plus bounded C++/SIMD flat-schedule and all-level
fused implementation complete. Immutable exact-mask schedules are consumed by
the production high Algorithm 5 and low Algorithm 4 decoders. The Algorithm 5
inverse schedules also have a promoted GF8 SSSE3 and GF8/GF16 AVX2 accumulating
boundary; rejected or unmeasured backends retain the materialized oracle.
Encoder integration and the broader end-to-end crossover map remain.

The original scalar experiment is
`experiments/leopard2/non_power_of_two/c1/dependency_pruning.py`. It answers the
C1 question left open by the symbolic C0 model: can Leopard omit work caused by
shortened inputs and punctured outputs while computing the exact existing
dyadic parent transform? The answer is yes for the tested scalar model. This
experiment does not define an exact-length code and does not change a wire
profile.

## Preserved mathematics and wire identity

For a transform of length `N = 2^n` in `GF(2^m)`, `n <= m`, the evaluator keeps:

- Leopard's public Cantor-coordinate representation;
- the normalized LCH basis `Xbar_i`;
- the enclosing power-of-two additive subspace or aligned coset;
- the existing coordinate order and skew factors;
- every padded input as an explicit mathematical zero; and
- every unrequested output as present in the parent code, merely not computed.

The scalar skew generator reconstructs the field multipliers represented by
`FFTSkew` in `LeopardFF8.cpp`. Leopard stores their logarithms for table
multiplication; C1 stores the equivalent field elements so multipliers zero and
one are visible to schedule construction. The implementation then uses the
same forward butterfly

    x' = x + m*y
    y' = y + x'

and inverse butterfly

    y' = y + x
    x' = x + m*y'

in Leopard's existing stage order. A shift is an aligned additive-coset
coordinate, not a renumbered block. Child shifts are `shift` and
`shift + N/2`, where addition is XOR in the stored Cantor coordinates.

An independent oracle constructs every subspace polynomial directly,
normalizes every LCH basis polynomial, converts LCH coefficients to a monomial
polynomial, and evaluates it by Horner's rule. It shares neither the skew
generator nor butterfly graph with the candidate. This follows the active
subspace definitions in `docs/leopard2_math_and_sources.md` and the LCH basis in
R16/R17. Existing Leopard source is used only to freeze its schedule and wire
conventions.

## Exact pruning analysis

Plan setup performs two passes over the actual scalar butterfly DAG.

1. A forward structural pass marks whether each intermediate can be nonzero
   from the supplied active-input mask. Coefficient zero is honored exactly.
2. A backward pass starts at the requested-output mask and marks which incoming
   intermediates are required. It uses each butterfly's four exact matrix
   coefficients rather than assuming that both outputs depend on both inputs.

The intersection is a proof-carrying execution schedule: an operation is kept
only when a structurally live input contributes to a requested result. Identity
writes, all-zero outputs, coefficient-zero products, and unused single outputs
are removed during setup. The byte loop contains no mask branches. Plan setup
and byte execution remain separate in all forms.

The model supports arbitrary sparse input and output masks, although Leopard's
current encoding profiles primarily produce prefixes. A plan rejects input
data that violates its declared known-zero mask.

## Execution forms compared

| Form | Representation | Result |
| --- | --- | --- |
| Recursive | Pruned transform tree with a branch per visited node | Correct; small scalar branch overhead at tiny payloads. |
| Flat | Precompiled required-butterfly list | Correct; simplest C++ prototype candidate and usually the fastest or tied form. |
| Hybrid | Complete dyadic subtransforms plus scalar boundary operations | Correct; no clear Python advantage, but exposes regular regions to existing radix-4/SIMD kernels. |
| Generated | Experiment-only Python source with indices and constants embedded | Correct for selected common schedules; 6.9-8.2x faster than the generic scalar-symbol interpreter when called on prevalidated input, but source grows from 902 bytes at N=16 to 25,086 bytes at N=256. |

The generated callable validates the known-zero input mask by default, like the
other execution forms. Its isolated micro-row explicitly disables duplicate
validation after the benchmark has constructed valid input. No generated code
is used by the library or default build. A production study
should use templates, constexpr data, or offline generation with an explicit
code-size budget; it must not add runtime executable memory.

## Correctness evidence

The retained deterministic result is
`experiments/leopard2/non_power_of_two/c1/results/self_test.json`. Its SHA-256
is:

    6632eeef83cbb308ab5cc206b1b319aa334ed396ef6c125c1d929b687caa0aea

The 2026-07-16 run passed:

| Check | Count |
| --- | ---: |
| Direct-polynomial versus full LCH transform vectors | 134,404 |
| Full forward/inverse round trips | 134,404 |
| Every GF4 prefix geometry, forward and inverse | 816 plans |
| Every GF4 sparse output mask through N=8 for every input prefix | 2,400 plans |
| Every sparse input/output mask through GF4 N=4, plus all N=8 input masks against ten irregular output masks | 5,672 plans |
| GF8 dyadic boundary/prefix sweeps through N=256 | 856 plans |
| GF8 arbitrary sparse-input sweeps through N=256 | 144 plans |
| Recursive/flat/hybrid/generated requested-output comparisons | 29,687 executions and 187,666 symbols |
| Arbitrary-tail byte comparisons at 1, 7, 31, 65, 257, and 1,025 bytes | 253,686 bytes |
| Generated-kernel known-zero-mask rejection regression | 1 |

GF4 covers every aligned coset for each supported parent. GF8 covers the first
and last aligned cosets for sub-field-size parents and the full N=256 field.
For N=1 and N=2, the oracle sweep includes all field values; consequently the
GF8 N=2 portion alone covers all 65,536 two-coefficient messages. Inputs include
all-zero plans, empty output masks, prefixes on both sides of powers of two,
alternating sparse masks, edge-only masks, and deterministic irregular masks.

The deterministic self-test was run twice after the final execution changes
and reproduced the same JSON hash. Python development-mode checks and bytecode
compilation also passed. This isolated experiment allocates managed Python
objects, so C/C++ ASan and UBSan are not applicable; production translation
still requires the repository sanitizer matrix.

## Operation and memory results

Counts below are exact for the scalar schedules, per payload byte. Padded counts
already specialize Leopard's zero multiplier and identity destination rather
than charging a deliberately naive butterfly.

| N | Active inputs | Requested outputs | Kept / padded butterflies | Fixed multiplies pruned / padded | XORs pruned / padded | Loads+stores pruned / padded | Flat schedule |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 9 | 7 | 12 / 32 | 5 / 17 | 16 / 49 | 40 / 113 | 192 B |
| 64 | 33 | 17 | 63 / 192 | 32 / 129 | 80 / 321 | 206 / 705 | 1,008 B |
| 256 | 129 | 65 | 319 / 1,024 | 192 / 769 | 448 / 1,793 | 1,086 / 3,841 | 5,104 B |

These cases remove 62.5-68.8% of butterflies and 64.6-71.7% of modeled shard
loads and stores. The hybrid representation finds one complete regular subtree
in each case: 4, 32, and 192 butterflies respectively, leaving 8, 31, and 127
boundary butterflies. This supports calling existing fused kernels for complete
subtrees instead of interpreting every operation individually.

## Scalar timing evidence

The retained host-dependent result is
`experiments/leopard2/non_power_of_two/c1/results/benchmark_amd_9950x3d.json`.
Its SHA-256 is:

    d92695c50efe9a5a5f835066f0b815bfeea1f3c2a9ffdecb766a860a9b12ac2d

The run used Python 3.12.3 on an AMD Ryzen 9 9950X3D. CPU 0 was isolated with
`taskset`; the process affinity contained one CPU while the host exposed 32.
Other agents paused CPU- and memory-heavy work for this short phase. Every
form received two warmups, nine interleaved repetitions through 1 KiB and seven
above it. Execution order rotated each repetition. The maximum median absolute
deviation was 1.33% of the median; every other cell was below 0.7%. Setup is
excluded. Timed byte transforms are forward transforms; inverse transforms are
covered by the exhaustive correctness program.

| N | Active / requested | Shard bytes | Padded | Recursive | Flat | Hybrid | Best speedup |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 9 / 7 | 64 | 101.1 us | 36.9 us | 34.6 us | 36.1 us | 2.92x |
| 16 | 9 / 7 | 1,024 | 1.533 ms | 0.546 ms | 0.545 ms | 0.546 ms | 2.81x |
| 16 | 9 / 7 | 16,384 | 25.717 ms | 9.199 ms | 9.197 ms | 9.174 ms | 2.80x |
| 16 | 9 / 7 | 65,536 | 103.437 ms | 37.319 ms | 36.980 ms | 37.081 ms | 2.80x |
| 64 | 33 / 17 | 64 | 624.4 us | 177.2 us | 167.5 us | 181.1 us | 3.73x |
| 64 | 33 / 17 | 1,024 | 9.816 ms | 2.856 ms | 2.848 ms | 2.858 ms | 3.45x |
| 64 | 33 / 17 | 16,384 | 165.357 ms | 48.630 ms | 48.599 ms | 48.465 ms | 3.41x |
| 256 | 129 / 65 | 64 | 3.406 ms | 0.921 ms | 0.882 ms | 0.956 ms | 3.86x |
| 256 | 129 / 65 | 1,024 | 54.481 ms | 15.347 ms | 15.300 ms | 15.375 ms | 3.56x |

These are deliberately scalar Python timings. They confirm that the reduced
memory and arithmetic work survives schedule traversal overhead, but they do
not predict a C++ SIMD speedup, cache behavior, or production crossover. The
large gap between generated and generic single-symbol interpreters is mostly
Python object/property overhead; it justifies specialization research, not a
7-8x production claim.

## Bounded C++ translation

The parent-preserving flat form is now translated into the internal C++11
implementation in `Leopard2Plan.cpp` and `Leopard2Plan.h`. GF8 and GF16 bind
the shared compiler to their existing `FFTSkewStorage` tables through
`PreparePrunedTransformPlan`; the experiment therefore consumes the exact
legacy logarithms rather than regenerating constants in a second production
path. It is not exposed by `leopard2.h`; decode plans own and dispatch the
immutable schedules without changing profile identity or the public ABI.

Setup expands the existing radix-2 graph once, records input liveness before
each operation, propagates exact output dependencies backward, removes
identity writes, and publishes the candidate plan only after complete
validation. The immutable plan owns masks, flat operation descriptors, fused
start indices, requested structurally-zero output indices, and, for inverse
plans, the compact root dependency flags described below. Byte execution
allocates nothing. It calls the selected scalar, SSSE3, or AVX2 fixed-product and
butterfly table explicitly, so a process-global default cannot leak into a
lower-backend test.

The final operation vector reserves only the retained operation count, rather
than retaining capacity for the full padded graph after pruning. The temporary
raw and dependency vectors remain setup-only storage and are released before
the immutable plan is published.

The follow-up hybrid compiler recognizes every complete live and requested
four-coordinate two-layer group. It records a start index on the unchanged four
radix-2 descriptors and executes them with the selected backend's mature
`Butterfly4` kernel. Forward scheduling groups the two cross-half operations
with their two child-pair operations before descending four grandchildren;
inverse scheduling visits the grandchildren first and emits the exact reverse
four-operation group. Operations at distinct offsets and in distinct additive
subspaces are disjoint, so this reordering preserves the same radix-2 DAG while
making complete groups contiguous at every pair of transform layers. Setup
validates all four coordinates, multiplier equality, liveness, and output
dependencies before emitting a descriptor. Execution validates each selected
descriptor against the retained radix-2 entries. Inverse plans additionally
retain one liveness/dependency byte for each root pair so an accumulating sink
can reproduce identity rows that the compact in-place operation vector omits.
Ragged boundaries and an odd
unpaired transform layer retain their specialized two-way descriptors. The
sorted start-index list costs four bytes per fused group rather than one marker
byte per radix-2 operation.

Leopard stores `m` as a logarithm, with 255 and 65,535 as the GF8/GF16 zero
sentinels and log zero representing field element one. Setup needs only the
following coefficient predicates; it does not convert the whole skew table
back to field elements:

- `m == 0` iff the logarithm is the field sentinel;
- `m == 1` iff the logarithm is zero;
- `m + 1 == 0` iff `m == 1`; and
- `m + 1 == 1` iff `m == 0`.

When both inputs are live, execution uses the mature two-way backend butterfly.
Zero and one multipliers reduce to XORs. A one-live-input boundary uses the
dead peer as bounded temporary storage for copy, fixed multiply, or
multiply-add. Requested outputs proven structurally zero are cleared after the
flat schedule because a dead slot may have served as that temporary. The
caller must provide disjoint shard slots and keep masked-off inputs at
mathematical zero; dead unrequested outputs are not preserved.

`tests/leopard2/test_pruned_transform.cpp` compares requested outputs both with
the full padded production transform and with the independent direct-polynomial
oracle in `tests/leopard2/direct_oracle.cpp`. The latter constructs normalized
LCH basis polynomials, converts to/from the monomial basis, evaluates by Horner,
and interpolates by Lagrange; it shares neither the skew table nor butterfly
graph with the candidate. The deterministic gate covers:

| C++ check | Result |
| --- | ---: |
| Qualified backends | scalar, SSSE3, AVX2 |
| Compiled/executed plans | 19,728 |
| Requested bytes compared | 67,095,177 |
| Independent direct-polynomial symbols | 13,512 |
| Fused four-way descriptors exercised | 215,205 |
| Effective execution descriptors | 787,494 (four radix-2 entries count as one fused step) |
| Largest owned plan in the gate | 6,946,992 bytes on this libstdc++ host for a complete GF16 N=65,536 plan |
| Exhaustive sparse masks | every input/output mask at N=2 and N=4, forward and inverse, first and last aligned cosets, GF8 and GF16 |
| Larger parents | GF8 through N=256; GF16 through N=1,024 |
| Real profile masks | high message-tail IFFT and transmitted/holey parity FFT; low shortened-message IFFT and final/holey parity-block FFT for GF8 (100,30), (17,100) and GF16 (1000,200), (257,700) |
| Shard tails | GF8 1, 7, 17, 63, 64, 65, 129 bytes; GF16 2, 18, 62, 64, 66, 130 bytes |
| Shared-plan concurrency | eight threads, sixteen executions each, including read-only accumulating execution for GF8 and GF16 per available backend |
| Fused descriptor integrity | complete GF8/GF16 forward and inverse plans emit `(N/4)*floor(log2(N)/2)` sorted, non-overlapping descriptors across all paired layers |

The complete Debug CTest graph passed 49/49 after initializing the checkout's
`sse2neon` submodule. GCC 13.3 and Clang 18.1 strict Release builds passed with
`-Werror -Wpedantic`. Clang 18 ASan+UBSan and TSan builds each repeated all
19,728 prototype cases and the concurrency gate without a diagnostic.
A tests-off Release archive also built cleanly, exported no test-hook symbol,
and passed the repository's fail-closed SSE2/SSSE3/AVX2 member-isolation check.

The initial fused-leaf follow-up repeated the focused matrix after its compact
start-index representation under GCC 13.3 strict Release, Clang 18.1
ASan+UBSan, and Clang 18.1 TSan (OpenMP disabled in sanitizer builds); all
passed. The complete Release CTest graph was effectively 49/49: the first run
passed 48 tests and failed only because the fresh worktree lacked the tracked
`sse2neon` contents, then the exact Visual Studio project test passed after
`git submodule update --init sse2neon`. Its tests-off archive again contained
no test hooks and passed portable ISA isolation. Independent review is still
required before that follow-up can be integrated or benchmarked for promotion.
The subsequent all-level grouping repeats these gates separately before it can
supersede the leaf-only checkpoint.

The all-level candidate repeated the 19,728-plan focused matrix and the maximum
GF16 plan under GCC 13.3 strict Release, Clang 18.1 ASan+UBSan, and Clang 18.1
TSan; all passed. Its fresh Release graph passed 49/49, including CUDA
optionality and the Visual Studio project scanner after initializing the tracked
submodule. A tests-off archive contained no test hooks and passed portable ISA
isolation. This is candidate validation, not promotion: independent algebra and
implementation review plus isolated end-to-end crossover evidence are still
required.

These are correctness results, not timing evidence. Other workers were active
on the host, so no cache-sensitive or authoritative crossover number was
collected for this checkpoint.

## Production decoder integration

`leo2_decode_plan_create` records the exact parent coordinates selected for a
specialized decode. For each incomplete nonempty P- or T-sized input block it
compiles an inverse schedule with the exact selected-input mask and all block
outputs requested. Complete and empty blocks retain their simpler mature
paths. Algorithm 4 additionally compiles one P-point forward schedule with all
accumulator inputs live and exactly the missing systematic coordinates
requested. Algorithm 5 compiles one shifted T-point output schedule for each
block containing a missing original. A schedule is retained only when its
operation count or zero/one specialization removes byte-heavy work; otherwise
execution falls back to the mature prefix/dependency transform.

Both the materialized and side-sized tiled Algorithm 4/5 kernels consume the
same immutable schedules for GF8 and GF16. Execution does not allocate, and
the selected context backend is passed explicitly through every retained
butterfly. Algorithm 5 begins each exact sparse output schedule with an
immutable-source out-of-place root layer, then resumes the retained schedule;
this removes its former whole-T evaluator copy without replacing the exact
output mask with a prefix. The decode-plan schedule gate compares the pruned
kernels with the unpruned prepared kernels and reports, for its deterministic
GF8/GF16 sparse patterns:

| Decoder | Retained / padded-equivalent operations |
|---|---:|
| Low Algorithm 4 | 3,330 / 3,696 |
| High Algorithm 5 | 3,608 / 4,752 |

The low acceptance gate independently compares the public production path to
direct interpolation across materialized/tiled modes, byte tails, GF8/GF16,
one-shot and reused plans, batch execution, parity re-encoding, legal aliases,
guards, and concurrent use. Context-backend tracing exercises the integrated
low path under every qualified scalar/SSSE3/AVX2 table. These are correctness
and structural operation-count results; the separate exact inverse-sink
crossover below applies only to Algorithm 5 syndrome construction.

### Exact Algorithm 5 inverse accumulation boundary

Later nonempty Algorithm 5 input blocks are dead immediately after their
inverse transform is XOR-reduced into the T-shard syndrome accumulator. The
original exact-mask integration nevertheless materialized every requested root
output and then made a separate T-shard XOR pass. The inverse executor now
stops before that root layer and applies the exact rows directly to the
accumulator. For one root pair the normalized inverse butterfly is

    y' = x + y
    x' = x + m * (x + y)

The executor specializes structural zero inputs and multipliers zero and one.
When both inputs and outputs are live it calls the selected backend's mature
read-only `IFFTButterfly2Xor`; sparse rows use its fixed multiply-add and XOR
primitives. A complete fused-four descriptor that crosses the root is split
into its two independent child butterflies and two accumulating root rows;
complete groups below the root remain fused. This also covers T=2 and an odd
number of transform layers without a special algebra.

Plan setup records the full root's flags separately because an identity row
may be correct but absent from the compact write schedule. Execution validates
the masks, every buffer pointer, operation descriptor, fused-four descriptor,
root flags, field, byte count, and required backend entries before its first
write. Rejection therefore leaves both work and accumulator untouched and the
existing materialize-then-XOR path remains an atomic fallback. The plan and
executor allocate nothing and are immutable under concurrent use.

The production crossover is deliberately narrower than executor correctness.
Pinned measurements promote the sink for GF8 AVX2 at 1 KiB and above, GF8
SSSE3 at 64 KiB and above, and GF16 AVX2 at 1 KiB and above. Smaller vector
shards keep the materialized schedule because the fixed plan preflight
otherwise dominates. GF16 SSSE3 improved by 4.86% at its 64-KiB target, but
that does not clear the project's default 5% complexity threshold, so it also
keeps materialization. Scalar initially regressed by 2.44% for GF8 and 1.78%
for GF16, so scalar plans discard the unused root metadata and keep the prior
materialized schedule. Rejected GF16 SSSE3 and unmeasured NEON plans discard it
as well. The final GF8 scalar confirmation was 0.41% faster than the matched
control, which is effectively neutral. NEON is left on the materialized path
until native hardware can establish an end-to-end crossover; cross-compilation
alone is not a performance claim.

The deterministic C++ gate now executes both materialized and accumulating
forms for every inverse case. It covers every input/output mask at T=2 and T=4,
prefixes and deterministic holes through GF8 T=256 and GF16 T=1,024, shifted
cosets, arbitrary requested outputs, GF8 byte tails, GF16 compact tails,
scalar/SSSE3/AVX2 tables, malformed-boundary atomic rejection, and shared-plan
concurrency. Atomic rejection covers corrupt root flags, out-of-range retained
operations, and invalid fused descriptors without changing either caller
workspace. The final candidate passed the complete Release graph (73/73),
focused Clang 18 ASan+UBSan and TSan builds, and strict Clang 18 and GCC 13
builds. Tests-off dual-field, GF8-only, and GF16-only archives built cleanly and
exported no test-hook symbols. A fresh Clang 18 ASan+UBSan libFuzzer campaign
completed 2,000 inputs with 6,119 covered edges and no crash. After the final
dispatcher cleanup, another 2,000-run continuation rebuilt the exact source,
reused that corpus, reached 6,132 edges, and also found no crash.

The authoritative crossover used physical core 14 on the 2026-07-18 host with
SMT sibling 30 reserved and idle. Both CPUs were 100% idle in the before/after
`mpstat` samples. The governor was `powersave` under `amd-pstate-epp`; perf
counters were unavailable because `perf_event_paranoid` was 4. Candidate and
control ran in ABBA order with one context thread, OpenMP dynamic teams
disabled, two or three warmups, nine, eleven, or fifteen retained samples per
outer run, and four outer runs per side. The control was commit
`354274b2765d79ad`.
Pooled medians below are decode execution only; setup remained separately
reported by the benchmark JSON. The crossover-candidate and control executable
SHA-256 values were `8b5a2a85ef85895988def2b44682537c79ee76cae776cb6cff103501981b8df0`
and `64bba105f60b7968daa5d783a502dfbaf8b377212c5f017e6b685d353c3b4a99`.
The two R=32 AVX2 target cells and scalar fallback were repeated after the
final scalar metadata policy; the broader vector matrix preceded that
metadata-only correction and exercises the same vector execution path.

| Field/backend | K,R,loss | Bytes | Kernel | Control | Candidate | Delta |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| GF8 AVX2 | 224,32,16 | 64 KiB | tiled | 1,670.6 us | 1,577.0 us | -5.60% |
| GF8 AVX2 | 192,64,32 | 64 KiB | tiled | 2,014.7 us | 1,936.3 us | -3.89% |
| GF8 AVX2 | 240,16,8 | 64 KiB | tiled | 1,108.8 us | 1,060.7 us | -4.33% |
| GF8 AVX2 | 240,16,8 | 64 KiB | materialized | 1,120.6 us | 1,084.5 us | -3.22% |
| GF8 AVX2 | 240,16,8 | 4 KiB | tiled | 60.35 us | 59.27 us | -1.79% |
| GF8 AVX2 | 129,2,1 | 64 KiB | tiled | 198.11 us | 197.24 us | -0.44% |
| GF16 AVX2 | 257,32,16 | 64 KiB | tiled | 2,369.8 us | 2,191.9 us | -7.51% |
| GF16 AVX2 | 1000,200,100 | 64 KiB | tiled | 19,926.4 us | 19,336.6 us | -2.96% |
| GF8 SSSE3 | 224,32,16 | 64 KiB | tiled | 2,719.3 us | 2,597.0 us | -4.50% |
| GF16 SSSE3 | 257,32,16 | 64 KiB | tiled | 4,433.3 us | 4,211.7 us | -5.00% |
| GF8 scalar, final policy | 224,32,16 | 64 KiB | tiled | 11,964.2 us | 11,915.3 us | -0.41% |

Negative deltas are improvements. The qualified vector cells range from
neutral at T=2 to 7.51% faster, with no decode regression greater than 2% in
the measured neighbors. Encode is not changed by this dispatch; its pooled
medians stayed within 1.73% of control in the same alternating runs. The
strongest meaningful region clears the project's 5% promotion rule, while the
deterministic backend gate avoids the measured scalar loss.

The independent review then challenged the unmeasured tiny-shard region. The
same ABBA runner repeated 64 B through 64 KiB with longer reuse. These rows show
the raw sink candidate before its byte threshold; positive deltas are
regressions and therefore select the materialized fallback in production.

| Field/backend | Bytes | Control | Raw sink | Delta |
| --- | ---: | ---: | ---: | ---: |
| GF8 AVX2 | 64 B | 6.450 us | 7.659 us | +18.73% |
| GF8 AVX2 | 256 B | 9.658 us | 10.749 us | +11.31% |
| GF8 AVX2 | 512 B | 13.980 us | 14.905 us | +6.61% |
| GF8 AVX2 | 1 KiB | 25.994 us | 25.748 us | -0.95% |
| GF8 SSSE3 | 64 B | 7.469 us | 8.701 us | +16.49% |
| GF8 SSSE3 | 4 KiB | 165.603 us | 170.280 us | +2.82% |
| GF8 SSSE3 | 32 KiB | 1,324.258 us | 1,327.051 us | +0.21% |
| GF8 SSSE3 | 64 KiB | 2,719.269 us | 2,596.995 us | -4.50% |
| GF16 AVX2 | 64 B | 12.964 us | 14.112 us | +8.86% |
| GF16 AVX2 | 256 B | 18.769 us | 19.764 us | +5.30% |
| GF16 AVX2 | 512 B | 27.198 us | 27.497 us | +1.10% |
| GF16 AVX2 | 1 KiB | 44.089 us | 42.782 us | -2.96% |
| GF16 SSSE3 | 64 B | 14.200 us | 15.120 us | +6.48% |
| GF16 SSSE3 | 256 B | 25.496 us | 25.844 us | +1.36% |
| GF16 SSSE3 | 512 B | 42.880 us | 42.588 us | -0.68% |
| GF16 SSSE3 | 1 KiB | 76.767 us | 74.774 us | -2.60% |

The deterministic thresholds above exclude every credible regression. The
T=2 64-KiB edge was repeated separately and ranged from neutral to 3.01%
faster across GF8/GF16 and SSSE3/AVX2, so the small transform shape does not
require another dispatcher dimension.

A final ABBA confirmation rebuilt both sides with tests disabled, so the
measured archives were the production library rather than the test-instrumented
counter build. The candidate and control executable SHA-256 values were
`a7606bfe7410363ae715245798679c3c7eaca5e6d9de8fdd8975248d39d5f67d`
and
`1776039cfcce03b25e9a62be49e2a81d33262124fe80d9089602d453ef51044f`.
The crossover manifest SHA-256 was
`4d9c8ea08bf15feec4ea4d2cb0e3409833abd5e7c5e8dc70d5227c2f8640bc67`;
the ordered raw-result-set SHA-256 was
`b912a6c09600144e8aea78d4c7c4e2d102df03b7c33cd422bfed6e4db8163761`.
All candidate/control original, parity, and recovered-output digests matched.

| Production crossover cell | Measured path | Control | Candidate | Delta |
| --- | --- | ---: | ---: | ---: |
| GF8 AVX2, 256 B | materialized | 10.116 us | 10.084 us | -0.31% |
| GF8 AVX2, 1 KiB | sink | 26.993 us | 26.730 us | -0.97% |
| GF8 SSSE3, 4 KiB | materialized | 175.867 us | 169.519 us | -3.61% |
| GF8 SSSE3, 64 KiB | sink | 2,875.195 us | 2,650.424 us | -7.82% |
| GF16 AVX2, 256 B | materialized | 19.035 us | 19.040 us | +0.03% |
| GF16 AVX2, 1 KiB | sink | 44.063 us | 42.219 us | -4.19% |
| GF16 SSSE3, 64 B | materialized | 14.413 us | 14.103 us | -2.15% |
| GF16 SSSE3, 1 KiB | raw sink, rejected | 76.854 us | 74.580 us | -2.96% |

The boundary sink rows establish a safe crossover only; they are not promotion
evidence by themselves. Separate tests-off 64-KiB target cells applied the
project's default 5% rule per backend. Their manifest SHA-256 was
`0d5d192fda0257a9cd980c3958918d22d5ffb4b9e90180c3da72d0545e182f80`
and ordered raw-result-set SHA-256 was
`ac26c5592fbeec51e58c52b774c4d0223221e837128f704bea37d00674a566da`.

| Production target cell | Control | Sink | Delta | Decision |
| --- | ---: | ---: | ---: | --- |
| GF8 AVX2, 64 KiB | 1,759.390 us | 1,638.989 us | -6.84% | promote |
| GF8 SSSE3, 64 KiB | 2,875.195 us | 2,650.424 us | -7.82% | promote |
| GF16 AVX2, 64 KiB | 2,256.830 us | 2,098.504 us | -7.02% | promote |
| GF16 SSSE3, 64 KiB | 4,404.303 us | 4,190.257 us | -4.86% | reject |

The retained target wins are all above 5%; deterministic field/backend/byte
gates avoid every rejected region. Tests-off fallback cells had no decode
regression greater than 0.03%. The unchanged encode path stayed within 1% of
control in every corresponding production cell, confirming that larger
variations seen in test-hook executables were instrumentation-layout effects
rather than a production regression.

## C++ transform benchmark checkpoint

`bench/leopard2/pruned_transform_benchmark.cpp` now provides the missing
primitive-level A/B instrument needed before codec integration.  It is built
only when both tests and benchmarks are enabled because it deliberately uses
test-only access to the full padded LCH kernels.  It is not installed, is not
linked by production callers, and cannot change dispatch.

For one exact prefix or deterministic-holey mask cell it constructs the
immutable plan once, removes the fused descriptors from a copy to obtain the
flat executor, and verifies both copies against a separately executed full
padded transform before taking any timing sample.  For a holey mask the full
baseline executes through the last live input or requested output, matching
the existing prefix-truncation contract.  It reports:

- plan setup samples independently of byte execution;
- full, flat, and all-level-fused raw samples in rotating order;
- padded butterflies, retained operations, fused descriptors, effective
  execution steps, and exact vector-capacity footprint;
- the requested and actually qualified backend;
- a requested-output digest and source SHA/dirty state; and
- the number of preallocated workspaces admitted by an explicit memory cap.

Resets occur outside each timed interval.  Multiple preallocated workspaces
avoid applying a sparse-input plan to its own dead temporary outputs on the
next iteration.  The JSON labels itself non-authoritative: host topology,
singleton affinity, idle SMT sibling, reservation ownership, frequency state,
cache/TLB counters, executable identity, timeout capture, and deterministic
replay remain the responsibility of the external evidence runner.  In
particular, development smoke numbers collected while other workers are active
must not be used for the 10 percent promotion threshold.

Build and validate the JSON contract with:

    cmake -S . -B build/c1-benchmark -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/c1-benchmark \
        --target bench_leopard2_pruned_transform -j 128
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
        python3 tools/leopard2_pruned_transform_benchmark_json_test.py \
        build/c1-benchmark/bench_leopard2_pruned_transform

Replace 128 with the allowed CPU count when fewer CPUs are available.  A
representative high/low boundary cell can then be emitted with explicit
parameters, for example:

    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
        build/c1-benchmark/bench_leopard2_pruned_transform \
        --field gf8 --size 64 --shift 0 --direction forward \
        --input-active 33 --output-requested 17 --bytes 65536 \
        --backend avx2 --iterations 31 --samples 9 \
        --setup-iterations 31 --memory-mib 512

This checkpoint closes only the raw transform-instrumentation gap.  The
external authoritative runner, real encode/decode integration, alias/scatter
and public-codec integration fuzz gates, end-to-end neighboring cells, and a
deterministic dispatcher remain open.

## Plan/executor fuzz boundary

`tests/leopard2/fuzz_pruned_transform.cpp` exercises the internal plan compiler
and executor without exposing them through the public ABI.  Each input selects
GF8 or GF16, a dyadic parent through 256 coordinates, an aligned coset,
forward or inverse direction, arbitrary sparse masks, byte tails, and one
qualified scalar/SSSE3/AVX2 backend.  It compares the flat and fused schedules
with a full padded transform on every requested coordinate.  It also checks
that malformed mask values leave an existing plan byte-for-byte equivalent,
and that invalid stored indices, multipliers where representable, zero-output
indices, and fused descriptors fail before an out-of-range access.

The ordinary CTest graph builds a deterministic replay driver.  A separate
Clang libFuzzer target compiles the library and isolated SIMD members directly
with coverage, ASan, and UBSan, so the fuzz boundary is not an uninstrumented
archive call.  Both targets require the internal test-hook ABI for the full
transform oracle and are never installed.

Run the deterministic gate with:

    cmake --build build/c1-benchmark \
        --target leopard2_pruned_fuzz_smoke -j 128
    OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
        build/c1-benchmark/leopard2_pruned_fuzz_smoke \
        0xc1f022d17 1024

Build the coverage-guided target with Clang and run one bounded campaign with:

    cmake -S . -B build/c1-fuzz -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF \
        -DLEO2_BUILD_FUZZERS=ON -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++
    cmake --build build/c1-fuzz \
        --target leopard2_pruned_plan_fuzzer -j 128
    mkdir -p build/c1-fuzz/corpus
    ASAN_OPTIONS=detect_leaks=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
        build/c1-fuzz/leopard2_pruned_plan_fuzzer \
        build/c1-fuzz/corpus -runs=2000 -max_len=256 -timeout=5 \
        -rss_limit_mb=8192

Run the fuzzer under an external per-job memory cap; the libFuzzer RSS option
is not a substitute for that containment.  Clang 18's internal RSS accounting
on the 2026-07-17 host reported about 4,173--4,180 MiB from process start and
incorrectly tripped its default 2,048 MiB limit.  A fixed-input control measured
45,312 KiB maximum resident set with `/usr/bin/time -v`.  The explicit 8 GiB
override above avoids that false internal limit while retaining an external
memory cap.

Replace 128 with the allowed CPU count.  With the RSS override, the 2026-07-17
development campaign completed 2,000 inputs under Clang 18 ASan+UBSan with no
crash, and the strict GCC plus Clang sanitizer replay each completed 1,024
deterministic inputs.  The process affinity exposed 28 CPUs rather than the
requested 128; all 28 ran distinct ASan+UBSan seeds concurrently for another
7,168 inputs.  The 45,312 KiB fixed-input measurement is a control, not a claim
about peak campaign memory. These earlier campaigns are fuzz/correctness
evidence; the later Algorithm 5 sink crossover and production policy are
reported separately above.

## Disposition

The following are now promoted behind tests and explicit production dispatch:

- exact forward-live/backward-needed plan construction;
- flat boundary operation lists;
- complete-subtransform descriptors that call existing fused kernels at every
  paired transform layer; and
- zero/one multiplier plus identity-write specialization;
- exact input/output schedules in the Algorithm 4 and Algorithm 5 decoders;
  and
- the GF8 SSSE3 and GF8/GF16 AVX2 exact Algorithm 5 inverse accumulation
  boundary.

Do not promote the Python executor, measured timing thresholds, or generated
source. Keep recursive execution as an oracle/debug form. The hybrid form is
algebraically useful but not independently faster in this scalar environment;
its production case depends on replacing complete groups with real radix-4 or
SIMD kernels. Generated common-pair kernels remain a separate code-size-limited
experiment.

This clears the C1 scalar correctness threshold, C++ translation, decoder
integration, backend determinism, arbitrary tails, immutable-plan concurrency,
sanitizers, fuzzing, and the bounded Algorithm 5 sink crossover. Remaining C1
work is:

- encoder integration using real profile masks and shifted blocks;
- a broader measured choice between fused and boundary descriptors by backend,
  shard size, and plan sparsity;
- native NEON measurement before enabling its exact inverse sink;
- code/table footprint and instruction-cache measurement; and
- the full region dispatcher study covering neighboring parameters where
  little work prunes.

## Reproduction

Run the deterministic correctness program:

    PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c1/dependency_pruning.py \
        self-test --output /tmp/c1-self-test.json

Compare it with the retained artifact:

    sha256sum /tmp/c1-self-test.json

Run the scalar timing program on one otherwise idle allowed CPU:

    PYTHONHASHSEED=0 taskset -c 0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c1/dependency_pruning.py \
        benchmark --output /tmp/c1-benchmark.json

CPU 0 is valid only on the recorded host. On another machine, select one CPU
from the process's actual affinity mask. Timing JSON is expected to vary by
host; mathematical counts and the self-test JSON are deterministic.

No default CMake target imports this experiment, and it requires only the
Python standard library.

Build and run the default-off C++ prototype gate:

    cmake -S . -B build/c1b -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=OFF
    cmake --build build/c1b --target leopard2_pruned_transform_test -j 128
    OMP_NUM_THREADS=1 build/c1b/leopard2_pruned_transform_test

Replace 128 with the allowed CPU count when fewer CPUs are available.
