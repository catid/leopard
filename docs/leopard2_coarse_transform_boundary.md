# Leopard2 coarse transform backend boundary

Status: implemented and correctness-qualified.  Performance promotion remains
pending an isolated, exact-source comparison with the main-branch codec.  This
change does not select a new codec path, profile, field, or wire format.

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

- GCC 13 Release: 66/66 tests, including the complete legacy golden-vector,
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

After rebasing onto coordinator commit `866dd3a`, a fresh Release build passed
70/70 tests with `OMP_NUM_THREADS=1` and 30-way CTest scheduling.  This includes
the newly integrated batch-aliasing coverage.  Limiting each test's implicit
OpenMP team avoided nested oversubscription while the independent test
processes used every allowed CPU.

No authoritative timing was collected while other workers occupied the host.
The candidate must be compared against exact main and a matched parent control
on an isolated physical core before any throughput claim or wider promotion.
Remaining optimization opportunities include radix-two stage grouping,
out-of-place coefficient-stage grouping, pruned-schedule batching, and fully
inlined/generated per-size backend loops.
