# Leopard2 coarse transform backend boundary

Status: implemented and correctness-qualified.  The exact-main comparison is
complete; the smaller matched-parent performance effect of this boundary is
not yet qualified because system activity repeatedly reached the measured
core's SMT sibling.  This change does not select a new codec path, profile,
field, or wire format.

## Problem and scope

The portable Leopard2 build keeps SSSE3 and AVX2 instructions in dedicated
translation units so that runtime dispatch cannot execute an unsupported ISA.
Before this change, the field scheduler invoked a context-selected function
pointer for every radix-four butterfly.  The exact main-branch implementation
at `6e5725e` instead emits transform arithmetic directly from its field
translation units.  Disassembly of matched production objects confirmed that
the portable path repeatedly loaded and called the private `Ops` entry inside
the transform loop, while main used local direct calls.  This dispatch
granularity is one concrete part of the remaining main/current transform gap.

This checkpoint moves the loop over one complete contiguous radix-four group
behind the private backend boundary.  A field transform now makes one
context-routed call for the group; the scalar, SSSE3, or AVX2 translation unit
then walks its independent shard lanes with direct backend-kernel calls.  The
following operations have coarse entries:

- GF8 forward radix-four groups;
- GF8 inverse radix-four groups;
- GF8 inverse radix-four groups fused with XOR accumulation; and
- GF16 forward and inverse radix-four groups, including the existing
  backend-and-byte-size fused/split policy.

The distance-one stage stays on the established leaf entry because its group
contains only one butterfly.  Pruned transforms, error-bit transforms,
radix-two remainder stages, and the out-of-place first coefficient layer also
retain their existing leaf entries.  The change is therefore a deliberately
small production boundary, not a claim that all exact-main code-generation
differences have been removed.

## Structural effect

For a complete transform of `m = 256` shards, four radix-four stages previously
made `64 * 4 = 256` indirect four-way calls.  The coarse form makes one range
call for each group in the first three stages (`1 + 4 + 16 = 21`) and 64 leaf
calls for the distance-one stage, for 85 context-routed calls.  This is a 66.8
percent reduction in this class of indirect dispatch.  Arithmetic, memory
traffic within each butterfly, and transform order are unchanged.

Test-only tracing counts both quantities independently.  It preserves the
logical butterfly count used by the existing transform-callsite assertions,
and separately requires representative GF8/GF16 high/low encodes to execute a
nonzero number of range entries and fewer range calls than logical leaves.

## Portability and context selection

The range callbacks are private members of the immutable `Ops` table; no
public C or C++ ABI changes.  Scalar, SSSE3, and AVX2 tables all implement the
entries, and builds omitting GF8 or GF16 publish null entries for the omitted
field.  Startup known-answer tests reject a backend before publication when a
required entry is absent or disagrees with its already-qualified leaf kernels.

On portable x86 builds, the field translation units contain no SSSE3 or AVX2
instructions and dispatch into the isolated backend object.  Targets that
intentionally retain legacy in-field SIMD, including SSE2NEON-style builds,
keep the established default-context leaf route rather than accidentally
calling the scalar range entry.  Explicit non-default contexts still use their
selected private backend.

The startup range tests cover modulus sentinels, nonzero multipliers, forward,
inverse, inverse-plus-XOR, fused and split GF16 forms, unaligned data, vector
tails, and guard bytes.  Codec tests then cover both profiles and fields,
legacy parity, high/low decode, arbitrary counts, concurrent execution, and
transform differential behavior.

The range KAT deliberately caps its GF8 cases at 129 bytes.  The independently
qualified leaf KAT retains the 1,025-byte cutoff case, while the range case
still spans multiple vectors and a tail without placing roughly 50 KiB of
temporary arrays on the initializing thread's stack.  GF16 range cases include
64, 128, and 130 bytes to cover both exact fused thresholds and the first
larger split-policy tail.

## Qualification checkpoint

The implementation passed these source gates before integration:

- GCC 13 Release: 70/70 tests, including the complete legacy golden-vector,
  field-option, CUDA-optionality, and project-graph gates.
- GCC 13 with `-Wall -Wextra -Wpedantic -Werror`: build clean and 9/9 focused
  tests.
- Clang 18 with `-Wall -Wextra -Wpedantic -Werror`, ASan, and UBSan: build
  clean and 9/9 focused tests.
- GF8-only and GF16-only Release builds.
- Production tests-off build and object-code inspection confirming a single
  field-side range dispatch followed by direct calls to the AVX2 range's local
  kernel.

The focused set included backend startup tests, mixed runtime contexts,
high/low acceptance, legacy-high pruning and parity, GF16 legacy encoding,
arbitrary-count profiles, transform differential testing, and concurrent
encoding.

The final code commit is `07816757` on coordinator parent `049331e8`.  A fresh
Release build passed 70/70 tests with `OMP_NUM_THREADS=1` and 30-way CTest
scheduling.  This includes batch-aliasing coverage.  Limiting each test's
implicit OpenMP team avoided nested oversubscription while independent test
processes used the available CPUs.

The backend matrix then rebuilt and tested `auto`, `scalar`, `ssse3`, and
`avx2` variants.  All four passed, with 39, 38, 38, and 38 recorded tests
respectively and no source mismatch.  The generated matrix is
`.research/leopard2/coarse_transform/backend-matrix-5fd34d0/matrix.json`, has
SHA-256
`fc818df4479970121895e6395962ee40c23b0b1ff87d56c14fdcf75633dd0ab6`,
and records source fingerprint
`0deecf762e4e634d14e844b98c67ac51c61b851aa0c26e5baa9ea8b8e3486822`.
The associated evidence-contract tests passed their canonical and historical
replays plus 89 adversarial manifest mutations.

## Exact-main comparison

An isolated ABBA campaign compared exact main `6e5725eb` with Leopard2 code
`07816757` on CPU 14.  CPU 30, its SMT sibling, recorded zero non-idle jiffies
over the measurement interval.  At capture, the independent verifier accepted
the build, source, executable, affinity, and sample closure.  The generated
manifest is
`.research/leopard2/coarse_transform/exact-main-0781675-final-cpu14-v4/manifest.json`,
has file SHA-256
`e09c484c79659f1fe2691b52153e2efea995add63d7dc2f745c7b9ca4174c2f9`,
and records internal digest
`2bbd414819cce46efa64d6e3c93a317ebf77416ff0e62a0b2c7e66615c2fc8b1`.

The following values are exact-main time divided by Leopard2 time, so a value
above one means Leopard2 was faster.  Parentheses give the 95 percent
bootstrap confidence interval.

| Cell | Encode | First decode | Reused-plan decode |
| --- | ---: | ---: | ---: |
| GF8 XOR | 1.006 (0.972-1.042) | 1.016 (0.982-1.051) | 1.018 (0.984-1.053) |
| GF8 high, one loss | 0.798 (0.787-0.808) | 2.040 (2.009-2.071) | 2.045 (2.014-2.076) |
| GF8 high, full loss | 0.795 (0.780-0.810) | 1.683 (1.658-1.710) | 1.696 (1.670-1.723) |
| GF8 balanced, full loss | 0.940 (0.937-0.943) | 0.767 (0.754-0.781) | 0.769 (0.755-0.782) |
| GF16 inflation, eight losses | 0.927 (0.918-0.936) | 2.770 (2.760-2.781) | 2.790 (2.780-2.801) |
| GF16 high, one loss | 0.899 (0.888-0.910) | 3.814 (3.806-3.822) | 3.826 (3.819-3.834) |
| GF16 high, full loss | 0.899 (0.896-0.901) | 2.741 (2.728-2.755) | 2.753 (2.739-2.767) |
| GF16 large, eight losses | 0.852 (0.846-0.858) | 2.818 (2.798-2.838) | 3.079 (3.058-3.100) |

This supports the specialized high-rate decode advantage over exact main, but
not a broad speedup claim.  Transform encoding was 6.0 to 20.5 percent slower
in the listed non-XOR cells, and balanced GF8 full-loss decoding was about 23
percent slower.  Those gaps remain production optimization targets.

## Matched-parent qualification

The coarse-callback change itself must be measured against its immediate
parent, not inferred from the exact-main comparison.  Three matched-parent
attempts were deliberately rejected:

- v1 completed but CPU 30 recorded two non-idle jiffies, and it exposed a
  benchmark-provenance bug involving resolved external-library symlinks;
- v2 was terminated after later-disclosed overlapping compile and smoke-test
  work invalidated isolation; and
- v3 completed all 1,664 invocations and passed the corrected provenance
  checks, but CPU 30 recorded 36 non-idle jiffies out of 36,262 while an
  account-wide affinity watchdog continually moved visible user processes away
  from the benchmark pair.

The v3 diagnostic is retained at
`.research/leopard2/coarse_transform/matched-parent-a30ec74-cpu14-v3/EXTERNAL_ISOLATION_INVALIDATED.json`.
No v3 ratio is authoritative or used for promotion.  A valid retry requires
system-wide or privileged isolation that can also exclude activity outside the
current account.  Consequently, an authoritative matched-parent residual and
the requested balanced-GF8 forced-generic versus forced-specialized residual
remain open.

Remaining implementation opportunities include radix-two stage grouping,
out-of-place coefficient-stage grouping through a different execution layout,
pruned-schedule batching, and fully inlined or generated per-size backend
loops.

## Rejected out-of-place range follow-up

A July 2026 follow-up tested a private range callback around the first
out-of-place forward radix-four group used by low-profile evaluation and the
mature high-profile evaluator.  The first version moved the lane loop into
each backend and called its qualified leaf directly.  Two subsequent variants
forced the leaf inline and then duplicated the SSSE3/AVX2 GF8 vector body so
the fixed multiplication tables were loaded once per range.  The initial
implementation passed the full Release suite, strict GCC and Clang sanitizer
builds, GF8-only and GF16-only builds, and the AUTO/scalar/SSSE3/AVX2 backend
matrix.  Each follow-up rebuilt and passed the focused backend qualification
and explicit-context tests before its timing screen.

The variants were nevertheless rejected.  CPU-0-pinned, single-thread,
reversed-order diagnostics compared the candidate with its immediate parent
`d9ef32e` using 51 samples, 1,024 reuses per sample, and eight warmups.  The
final table-hoisted AVX2 candidate measured as follows; ranges show the two
counterbalanced process runs in microseconds per encode call:

| Cell | Parent | Candidate | Result |
| --- | ---: | ---: | --- |
| low GF8 (17,33), 64 B | 0.520-0.523 | 0.531-0.543 | 2-4% slower |
| low GF8 (64,192), 64 B | 2.450-2.476 | 2.647-2.657 | 7-8% slower |
| high GF8 (192,64), 64 B | 2.431-2.438 | 2.508-2.512 | about 3% slower |
| low GF8 (64,192), 1 KiB | 14.556-14.607 | 14.648-14.688 | neutral/slower |

GF16 screens were likewise neutral or slightly slower.  The callback reduced
the number of indirect calls, but that structural count did not translate to
a whole-codec win and violated the neighboring-regression gate.  The code was
therefore reverted rather than hidden behind production dispatch.  A future
attempt should change the data traversal or generate a complete first-stage
kernel, and must repeat immediate-parent whole-codec measurements instead of
assuming dispatch-count reduction is sufficient.
