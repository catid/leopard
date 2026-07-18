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

These are correctness and structural results, not authoritative timing.  The
host was running other project workers during implementation, so an isolated
same-source A/B crossover sweep remains required before making an end-to-end
throughput claim.
