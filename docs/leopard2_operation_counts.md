# Leopard2 operation-count model

`tools/leopard2_operation_counts.py` is a deterministic model of Leopard2's
byte-heavy execution paths.  Its transform and arithmetic schedule is
ISA-independent; a separate layer applies the selected production backend's
promoted decode-fusion predicates to logical traffic.  It reports:

- legacy-compatible high-rate encoding;
- the specialized high-rate original-recovery decoder;
- low-rate encoding and specialized low-rate decoding;
- the generic parent-transform decoder and a standalone generic transform; and
- the bounded direct-repair dispatcher path.

The model is deliberately separate from the production kernels.  It can reveal
whether a change altered the algorithmic schedule without confusing that change
with table layout, cache, or compiler effects.  Backend-qualified traffic is
identified explicitly and guarded against production-source policy drift.  It
is not a replacement for hardware counters or elapsed-time benchmarks.

## Classification of results

Every JSON/CSV metric carries a classification.

- `exact_schedule` means the value is derived from the actual radix-4 loop
  topology, prefix/mask pruning, coordinate layout, copies, zero fills, or
  public scratch-layout formula.  Fused transforms are expanded into equivalent
  radix-2 butterfly counts.  Decode scratch bytes are exact for the pointer
  width recorded in the report.
- `exact_schedule` also labels the reveal/scatter and Algorithm 5 syndrome
  payload split after applying a named backend's exact production crossover
  predicates.  It describes bytes crossing those software boundaries, not
  cache-line transfers or instruction counts.
- `derived_lower_bound` and `derived_upper_bound` bracket arithmetic affected
  by transform skew constants.  A nonzero-skew butterfly performs two field
  additions and one fixed multiplication; Leopard specializes a zero skew to
  one addition and no multiplication.  This standalone tool intentionally does
  not read a backend's private skew table.
- `estimated_lower_bound` and `estimated_upper_bound` describe logical shard
  traffic.  They exclude cache-line amplification, metadata, allocation,
  prefetch, write allocation, multiplication-table reads, and application-side
  traffic.  The lower transform bound assumes a fused kernel loads and stores
  both butterfly operands once.  The upper bound represents separate logical
  primitives that reread operands.
- `estimated_upper_bound` on direct repair reflects coefficient-dependent term
  cancellation.  The dense bound has K terms for each of L repaired originals;
  the generated plan may remove zero coefficients or specialize coefficient 1.

`estimated_full_workspace_passes_*` divides logical shard accesses by the
execution working-slot count.  It is a normalization for comparing algorithms,
not a cache-pass measurement.

Plan setup is excluded.  In particular, locator construction, direct matrix
inversion, coefficient generation, and erasure-plan allocation are not charged
to execution.  The benchmark suite reports setup separately.

## Profiles and deterministic decode pattern

The geometry follows the wire definitions:

- high: T = ceil_pow2(R), N = ceil_pow2(K + T), with parity at 0..T-1
  and public originals at T..T+K-1;
- low: P = ceil_pow2(K), N = ceil_pow2(P + R), with public originals at
  0..K-1 and parity at P..P+R-1.

For a loss count L, the default deterministic pattern loses originals 0..L-1,
keeps every other original, and selects public parity 0..L-1.  This mirrors the
plan's lowest-index deterministic parity selection and leaves exactly K input
coordinates.  `--loss-mask` selects a different original pattern.  It still
uses the lowest L parity coordinates.

Encoding accepts `--requested-parity all`, `none`, or comma-separated indices
and inclusive ranges such as `0-3,7`.  High encoding trims only the unused
suffix, matching the production encoder.  Low encoding additionally skips
empty parity blocks.  The generic forward transform accepts an arbitrary parent
coordinate mask through `--transform-output-mask`.

Transform kernels execute complete 64-byte tiles, so reports include both
requested shard bytes and `kernel_shard_bytes`; direct repair/XOR/copy consumes
the exact public length.  GF16 reports reject odd
physical byte counts because the native wire profile requires complete two-byte
symbols.  A padded-odd GF16 application payload must be translated to its even
physical wire size before it is passed to this model.

## Decode selection, fusion traffic, and scratch accounting

Schema version 3 retains schema v2's corrected pointer-versus-copy accounting
and adds the actual selected decode rule plus backend-qualified reveal and
syndrome boundaries.  A
transform decode owns an N-entry coordinate pointer map, but complete aligned
passes point those entries directly at the selected K caller shards.  Mapping a
pointer is metadata work; it is not a K-shard data copy.  The transform kernel
may subsequently copy or multiply those selected rows into its work area, and
the model charges that path-specific work separately.

`logical_copy_vectors` now covers only complete kernel-work copies.  Public
output and tail-boundary copies are reported in payload bytes so a one-byte
tail is not misrepresented as another whole rounded shard.  Algorithm 5 gathers
L complete requested payloads.  Algorithm 4 reveals every aligned prefix
directly and gathers only `L * (shard_bytes mod 64)` tail bytes.  Generic decode
retains the in-scratch reveal and scatter except where qualified GF8
SSSE3/AVX2 execution removes both for the aligned prefix.

Reports separate:

- `decode_reveal_inplace_temporary_payload_bytes`, the payload copied through
  the 64-byte alias-safe temporary used by in-place fixed multiplication;
- `decode_reveal_scatter_payload_bytes`, the public payload later copied out of
  transform scratch;
- `decode_reveal_direct_payload_bytes`, the aligned payload multiplied directly
  into caller output;
- `decode_syndrome_fused_accumulation_payload_bytes`, Algorithm 5 rows whose
  final inverse layer accumulates into the syndrome workspace; and
- `decode_syndrome_materialized_xor_payload_bytes`, rows retaining a separate
  materialize-then-XOR boundary.

The `estimated_backend_bytes_*` metrics apply those exact boundary deltas to the
broader logical traffic bounds.  They remain estimates because they do not
model caches, write allocation, tables, or backend instruction sequences.
Schema v3 applies only reveal and syndrome deltas to those totals; the existing
receive-source and weighted-locator fusion effects remain explicit detail
fields and are not silently folded into the aggregate.

The exact public scratch metrics mirror `leo2_decode_plan_scratch_size()` and
`leo2_decode_scratch_size()`:

- validation reserves 2K+R `AddressRange` records;
- transform plans reserve N coordinate pointers plus one pointer per work slot;
- materialized and generic decoders use N work slots;
- a forced tiled low decoder uses 2P work slots;
- a forced tiled high plan uses 2T+L slots, while the pattern-independent
  one-shot codec query conservatively uses 2T+R;
- ordinary specialized AUTO retains N slots when N is no larger than its tiled
  layout;
- direct XOR, copy, and bounded-repair plans retain only range-validation
  metadata and no shard-data or pointer-map storage; and
- a no-loss plan reports zero scratch, although its codec's pattern-independent
  one-shot query remains nonzero.

For a ragged final 64-byte tile, transform scratch reserves K+R fixed 64-byte
public-coordinate slots.  Exactly K selected coordinates are populated for a
valid plan; absent, surplus, shortened, and punctured coordinates do not cause
additional copies.  Reports therefore distinguish reserved tail bytes from
selected payload bytes and zero-padding bytes.  Work slots are sized to
`max(aligned_prefix_bytes, 64)` for a ragged shard, not to the rounded total
shard length.  `--pointer-bytes` declares the 4- or 8-byte ABI used to account
for `AddressRange`, pointer maps, alignment, and the corresponding `size_t`
limit.  The model mirrors production's checked add, multiply, and alignment
operations: it rejects a zero shard length, a value outside `uint64_t`, a value
that cannot round up to 64 bytes in the declared `size_t`, or a complete scratch
layout that overflows.  It never labels an unrepresentable Python integer as an
exact public-query byte count.

`--backend scalar|ssse3|avx2|neon` chooses the backend predicates used for
traffic accounting.  It never changes the structural butterfly count.
`--decode-selection path` forces the path named by the report and
`--decode-workspace materialized|tiled` selects that specialized workspace.
`--decode-selection auto` mirrors production AUTO, including bounded direct
repair and measured crossover rules; `specialized` mirrors
`LEO2_CODEC_FORCE_SPECIALIZED_DECODE`.  `--multi-item-batch` applies the
multi-stripe selector exception.  Every decode report records the selected
path, rule, matching automatic-rule mask, and required work slots.

The standalone selector mirror returns an immutable no-loss plan before byte
geometry validation, matching private path introspection and execution.  Public
plan scratch rejects zero bytes, while the pattern-independent codec scratch
query still validates and rounds its transform-capable byte geometry.  A full
report includes that codec query and therefore rejects an unrepresentable byte
length even when its selected plan path is no-op.

The production-linked
`leopard2_decode_scratch_probe` additionally emits AUTO rows.  Those rows record
the context's selected backend and the actual deterministic policy path at each
shard size.  The independent Python selector is cross-checked against those
production rows at every configured boundary; AUTO backend choice itself is a
host observation, never a universal crossover claim.

## Usage

Run internal invariants and the independent test suite:

    python3 tools/leopard2_operation_counts.py self-test
    python3 experiments/leopard2/operation_counts/test_operation_counts.py

The configured CTest suite also builds a production-linked public-query probe
and compares it to the Python formulas, including adversarial source mutations:

    ctest --test-dir build/release -R leopard2_decode_scratch_crosscheck \
      --output-on-failure

Report a full high-rate encode cell as deterministic JSON:

    python3 tools/leopard2_operation_counts.py report \
      --path legacy_high_encode --k 240 --r 16 --field gf8 \
      --shard-bytes 1024 --requested-parity all

Report a sparse low-rate parity request as CSV:

    python3 tools/leopard2_operation_counts.py report \
      --path low_encode --k 8 --r 248 --field gf8 \
      --shard-bytes 65536 --requested-parity 0,247 --format csv

Compare specialized and generic decoding for the same loss pattern:

    python3 tools/leopard2_operation_counts.py report \
      --path legacy_high_decode --k 240 --r 16 --field gf8 \
      --backend avx2 --shard-bytes 1024 --loss-mask 0,15 \
      --decode-workspace tiled

    python3 tools/leopard2_operation_counts.py report \
      --path generic_decode --profile high --k 240 --r 16 --field gf8 \
      --backend avx2 --shard-bytes 1024 --loss-mask 0,15

Report the actual production AUTO rule for a batch-shaped cell:

    python3 tools/leopard2_operation_counts.py report \
      --path legacy_high_decode --k 224 --r 32 --field gf8 \
      --backend avx2 --shard-bytes 32768 --loss-count 8 \
      --decode-selection auto --multi-item-batch

Model the generic parent FFT with a non-prefix output mask:

    python3 tools/leopard2_operation_counts.py report \
      --path generic_transform --profile low --k 100 --r 156 --field gf8 \
      --shard-bytes 4096 --direction forward \
      --transform-output-mask 0,7,63,127

The output is stable: JSON object keys and CSV metric rows are sorted, mask
ranges are canonical, floating pass estimates are rounded, and the schema has
an explicit version.  Schema v3 includes schema v2's exact scratch
byte/component metrics and adds selected-rule/backend-traffic fields.  Use
`--output FILE` to write either format directly.

## Representative validated counts

The self-test fixes two useful schedule checks:

- high K=240, R=16 encode: 15 full 16-point inverse transforms plus one
  16-point forward transform = 512 radix-2 butterfly equivalents, with 224
  fused accumulation XOR vectors;
- low K=8, R=248 encode: one 8-point inverse transform plus 31 full 8-point
  forward transforms = 384 radix-2 butterfly equivalents.  Its immutable
  coefficient block is read directly by each parity block's out-of-place first
  FFT layer, so the model charges no separate 8-shard coefficient copy per
  block.  A source-level mutation guard rejects the former whole-P copy loop
  before that zero-copy model is accepted.  The remaining 256 logical copies
  are the 8 systematic inputs staged for interpolation and the 248 requested
  parity outputs scattered to caller buffers.
- GF8 high K=240, R=16 decode with original 0 missing: on qualified SSSE3 and
  AVX2 backends, the incomplete parity and first-message blocks retain their
  two exact-pruned T=16 staging workspaces.  The remaining full blocks provide
  56 complete four-row receive groups that enter the first two inverse layers
  directly.  The receive boundary copies 16 live rows and zeroes 16 absent
  rows, removing 224 logical copy vectors (224 shard reads and 224 shard
  writes) from the former copy-first schedule.  These are exact
  boundary deltas.  Schema v3's absolute accounting treats the K received
  coordinate entries as pointer mappings, and charges Algorithm 5's separate
  copy-first work staging before applying this qualified SIMD delta.  Scalar,
  NEON, exact-mask plans, and other unqualified backends retain copy-first
  staging.
- GF16 reports retain the deterministic copy-first receive boundary.  The
  otherwise equivalent source-staging candidate was measured but did not meet
  the production promotion threshold, so GF16 reports charge every selected
  receive-row copy and report zero removed copy vectors.  Materialized decode
  stages the full parent in one regular pass; tiled decode stages only live
  tiles and skips completely empty later blocks.

It also checks every full transform size from 1 through 65,536 against the
closed form N log2(N) / 2, compares prefix schedules with independently built
prefix masks, checks monotonic pruning, exercises direct-repair bounds, verifies
stable JSON/CSV serialization, and covers production GF16 byte-length rejection.

## Limitations

The operation model follows the current algorithms, not every backend micro-op.
Its production-source guards intentionally fail when a guarded reveal or
syndrome crossover changes, requiring a schema/model review.  It does not count
locator-plan setup, field-table initialization,
thread synchronization, validation instruction cost, or NUMA effects.  Exact
decode scratch size does include the storage occupied by validation range
metadata; logical byte-traffic estimates do not treat metadata scans as shard
traffic.
It does not infer exact zero/one multiplier counts from private skew and direct
coefficient tables; those are explicitly bounded.  The model also assumes the
current deterministic received-subset policy rather than optimizing the subset.

When generated transform schedules become production inputs, their serialized
operation lists should expose exact zero/one multiplier counts.  At that point
this model can consume those schedules and tighten arithmetic bounds without
changing the output schema.
