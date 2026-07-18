# High-rate evaluator copy elimination

High-rate Algorithm 5 constructs one `T`-coefficient evaluator polynomial and
reuses it for every message coset containing a requested missing original.
Previously, both the materialized and tiled production kernels copied all `T`
coefficient shards into the destination workspace before every such forward
transform.  For `Q` requested output blocks and `B` physical bytes per shard,
that pass alone read and wrote `Q*T*B` bytes, or `2*Q*T*B` bytes of logical
memory traffic.

The evaluator block is now immutable.  A mature prefix transform begins with
the existing backend out-of-place butterfly machinery:

- `T=2` uses one two-way butterfly;
- GF8 `T>=4` uses `T/4` out-of-place fused two-layer butterflies;
- GF16 retains its measured fusion policy: a 64-byte kernel tile and a
  128-byte AVX2 tile use `T/4` four-way operations, while other sizes use
  `T/2` out-of-place two-way operations followed by the established in-place
  second layer; and
- every later transform layer resumes in place in the output workspace.

An exact sparse C1 output schedule uses the same immutable-source contract but
executes only its root radix-2 layer out of place.  The executor then skips the
equivalent root operations in the immutable plan and resumes its retained
operations, including qualified child four-way groups.  This preserves the
exact requested-coordinate mask instead of replacing it with a prefix
transform.  The first root multiplier is stored in plan setup; byte execution
does not allocate, construct masks, or inspect shard values.

The materialized kernel writes each evaluated block into its existing parent
workspace slots.  The tiled kernel writes into its reusable `T`-slot tile.
Both keep the evaluator accumulator unchanged across output blocks.  A trusted
plan-validation failure retains the previous copy-plus-plan executor as a
compatibility fallback, but test-only counters require that fallback count to
remain zero in production-path tests.

This is an execution-only change.  Field representation, active-parent basis,
coordinate order, locator values, reveal factors, parity bytes, requested
originals, and the legacy-high-v1 wire profile do not change.

## Receive-boundary fusion

Algorithm 5's first block reduction consumes received coordinate rows before
any locator multiplication.  On qualified SSSE3 and AVX2 backends, complete
four-row GF8 groups now enter the first two inverse LCH layers directly through
the established out-of-place backend operation.  For a group beginning at
`r`, the operation uses exactly the existing distance-one skew entries at
`r+1`, `r+3`, and `r+2`; the ordinary in-place decoder schedule resumes at
distance four.  Input rows remain immutable and the destination is the
existing caller-provided decode scratch.  No value-dependent plan setup or
allocation was added.

Groups containing a missing row retain explicit copy/zero staging before the
same two inverse layers.  The backend operation deliberately does not accept
null rows, and sharing one zero scratch row with four disjoint outputs would
weaken its alias contract.  Exact-mask pruned blocks also retain their existing
staging until a pruned immutable-source inverse executor is independently
qualified.

This boundary is enabled only for GF8 SSSE3 and AVX2.  An isolated scalar
screen was positive but did not clear the 5% promotion threshold, so scalar and
unmeasured backends such as NEON retain copy-first staging.  The equivalent
GF16 candidate was implemented and tested, but isolated same-source
measurements improved AVX2
execution by only 1.0--1.8% through 16 KiB and were neutral at 32 KiB.  Those
results did not meet the 5% production promotion rule, so every GF16 backend
uses the deterministic copy-first path.  Test hooks and the operation model
make that negative disposition explicit rather than leaving a dormant
size-dependent branch.  The materialized GF16 kernel performs that staging in
one parent-wide OpenMP pass, matching the control schedule and covering empty
and pruned blocks with zero rows before block transforms begin.  The tiled
kernel retains its established per-tile staging because it has no N-row
workspace.

Algorithm 4 is not changed by this fusion.  Its receive values are multiplied
by plan-owned locator factors before the first inverse transform.  Algorithm
5's later locator-weighted inverse has the same dependency.  Correctly fusing
either boundary requires a new backend operation that applies four independent
fixed multipliers before the inverse butterfly; the current unweighted
operation cannot be reused by reordering the diagonal multiplication.  That is
a separate measured experiment, not an assumed algebraic optimization.

For qualified GF8 K=240, R=16 with original zero missing, the incomplete
parity block and incomplete first message block both own exact-pruned input
plans, so their complete T=16 workspaces remain staged.  The other 14 full
message blocks contribute 56 four-row source operations.  Sixteen live rows
are copied and 16 absent rows are zeroed.  Relative to copy-first staging this
removes 224 shard reads and 224 shard writes per execution.  This is a boundary delta; the broader
operation model's pre-existing absolute decode copy total is not used as a
correctness oracle.  The scoped source guard checks the decoder helper itself
and rejects loss of the call, widening beyond SSSE3/AVX2, or accidental
promotion of the rejected GF16 variant.

## Reveal and scatter

The tiled Algorithm 5 path already fuses the final reveal-factor multiply with
its scatter into the per-request kernel-layout retention slot by calling the
out-of-place fixed multiplier once.  The materialized path keeps an in-place
multiply because its parent-coordinate work slot is the source consumed by the
common public gather.

The public gather remains separate.  It converts compact GF16 tails and
implements the documented public alias behavior; bypassing it only for some
aligned destinations would add an alias/layout branch without removing the
retention requirement for other calls.  No additional scatter fusion was
promoted without isolated end-to-end evidence.

## Correctness and structural evidence

The pruned-transform gate compares immutable-source execution to a separately
executed full transform for GF8 and GF16, scalar/SSSE3/AVX2, transform sizes 2
through 1,024 (plus GF8 size 256), aligned shifts, sparse output masks, and
non-vector-aligned byte counts.  It verifies that the source coefficients stay
byte-identical and reuses one plan concurrently in 128 executions.

The high decoder acceptance and plan-schedule gates retain the direct
polynomial/interpolation oracle, generic decoder, prepared decoder, and
materialized/tiled differential comparisons.  Their exercised matrix includes
GF8/GF16, compact tails, exact sparse output plans, repeated and concurrent
plans, parity re-encoding, guards, and supported public aliases.  Context
backend tracing requires every requested output block to enter an out-of-place
backend operation and requires zero compatibility-copy fallbacks.  The
operation-count self-test removes the former per-block copy charge and rejects
a source mutation that restores an unconditional whole-`T` copy.

Focused Release validation on the implementation branch completed:

- high acceptance: 4 profiles, 8 specialized and 4 generic executions, 32
  concurrent executions, 38 restored shards, 272 direct parity symbols, 187
  direct recovered symbols, 4 parity rebuilds, 2,234 guard checks, and 4
  overlap checks;
- decode-plan schedules: 1,726,644 dependency queries, 128 concurrent plan
  executions, no execution allocation, and high pruning retained 3,608 of
  4,752 padded-equivalent operations; and
- context tracing: both profiles, both fields, and every qualified local
  scalar/SSSE3/AVX2 backend passed with no copy fallback.

The final policy also passed the full 73-test Release suite, a strict Clang 18
`-Wall -Wextra -Wpedantic -Werror` build, and all 68 applicable Clang
ASan/UBSan tests.  The sanitizer matrix excludes only the portable-ISA binary
scanner because sanitizer instrumentation itself inserts `lahf`; its normal
Release instance passed.

## Isolated timing evidence

The source-boundary candidate at `a09d705` was compared with its immediate
control `3a00fd8` using three independent ABBA rounds on one pinned physical
core while reserving its SMT sibling.  Each sample retained source, binary,
archive, workload, and `/proc/stat` identities; all sibling non-idle deltas were
zero.  Ratios below are control time divided by candidate time, with 95%
intervals:

| Backend, field, and workload | Shard bytes | Ratio |
| --- | ---: | ---: |
| AVX2 GF8 K=240, R=16, one lost original | 64 | 1.027 [1.016, 1.038] |
| AVX2 GF8 K=240, R=16, one lost original | 64 KiB | 1.200 [1.196, 1.205] |
| AVX2 GF8 K=240, R=16, all 16 losses | 64 KiB | 1.039 [1.023, 1.055] |
| AVX2 GF16 K=1000, R=200, one lost original | 4 KiB | 1.018 [1.013, 1.023] |
| AVX2 GF16 K=1000, R=200, one lost original | 16 KiB | 1.017 [1.015, 1.020] |
| AVX2 GF16 K=1000, R=200, one lost original | 32 KiB | 1.002 [0.998, 1.006] |

The final backend-qualification screen used the same three-round ABBA design.
CPU 4 ran the benchmark while SMT sibling 20 accumulated zero non-idle
jiffies.  SSSE3 improved 4 KiB one-loss decode by 1.103x [1.100, 1.106] and
64 KiB by 1.096x [1.090, 1.103], clearing the target threshold.  Scalar
improved only 1.018x [1.015, 1.022] at 4 KiB and 1.015x [0.997, 1.034] at
64 KiB; it therefore remains copy-first despite the absence of a measured
regression.  At 64 KiB with all 16 losses, SSSE3 was 1.011x [1.000, 1.022]
and scalar was 1.006x [0.999, 1.014].

An initial GF16 rejection patch preserved copy-first arithmetic but
accidentally split materialized staging into one OpenMP region per nonempty T
block.  An eight-core forced-materialized screen found regressions of 3.3% at
64-byte shards, 1.4% at 1 KiB for K=1000, R=200, and 1.9% at 64 bytes for
K=4096, R=512.  Restoring the single parent-wide staging region removed that
schedule regression; the structural test requires exactly K selected-row
copies and N-K zero rows even when later input blocks are empty or pruned.

The GF8 K=240, R=16 candidate was also compared against an exact detached
Leopard main build at `6e5725e` under the same fail-closed isolation protocol.
One-loss decode was 2.680x [2.651, 2.709] faster on first use and 2.689x
[2.660, 2.719] faster with plan reuse.  At all 16 losses it was 1.725x
[1.713, 1.737] faster on first use and 1.739x [1.727, 1.752] with reuse.
These main-branch ratios describe the complete Leopard2 decoder, not the
isolated contribution of this one fusion; the immediate-control table is the
promotion evidence for the change itself.
