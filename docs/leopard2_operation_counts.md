# Leopard2 operation-count model

`tools/leopard2_operation_counts.py` is a deterministic, ISA-independent model
of Leopard2's byte-heavy execution paths.  It reports the algorithmic work for:

- legacy-compatible high-rate encoding;
- the specialized high-rate original-recovery decoder;
- low-rate encoding and specialized low-rate decoding;
- the generic parent-transform decoder and a standalone generic transform; and
- the bounded direct-repair dispatcher path.

The model is deliberately separate from the production kernels.  It can reveal
whether a change altered the algorithmic schedule without confusing that change
with AVX2, SSSE3, NEON, table layout, cache, or compiler effects.  It is not a
replacement for hardware counters or elapsed-time benchmarks.

## Classification of results

Every JSON/CSV metric carries a classification.

- `exact_schedule` means the value is derived from the actual radix-4 loop
  topology, prefix/mask pruning, coordinate layout, copies, zero fills, or
  scratch-slot formula.  Fused transforms are expanded into equivalent radix-2
  butterfly counts.
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

The current byte kernels execute complete 64-byte tiles, so reports include
both requested shard bytes and `kernel_shard_bytes`.  GF16 reports reject odd
byte counts because the current wire profile requires complete two-byte symbols.
Scratch-slot counts exclude pointer/range metadata and alignment padding.

## Usage

Run internal invariants and the independent test suite:

    python3 tools/leopard2_operation_counts.py self-test
    python3 experiments/leopard2/operation_counts/test_operation_counts.py

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
      --shard-bytes 1024 --loss-mask 0,15

    python3 tools/leopard2_operation_counts.py report \
      --path generic_decode --profile high --k 240 --r 16 --field gf8 \
      --shard-bytes 1024 --loss-mask 0,15

Model the generic parent FFT with a non-prefix output mask:

    python3 tools/leopard2_operation_counts.py report \
      --path generic_transform --profile low --k 100 --r 156 --field gf8 \
      --shard-bytes 4096 --direction forward \
      --transform-output-mask 0,7,63,127

The output is stable: JSON object keys and CSV metric rows are sorted, mask
ranges are canonical, floating pass estimates are rounded, and the schema has
an explicit version.  Use `--output FILE` to write either format directly.

## Representative validated counts

The self-test fixes two useful schedule checks:

- high K=240, R=16 encode: 15 full 16-point inverse transforms plus one
  16-point forward transform = 512 radix-2 butterfly equivalents, with 224
  fused accumulation XOR vectors;
- low K=8, R=248 encode: one 8-point inverse transform plus 31 full 8-point
  forward transforms = 384 radix-2 butterfly equivalents.

It also checks every full transform size from 1 through 65,536 against the
closed form N log2(N) / 2, compares prefix schedules with independently built
prefix masks, checks monotonic pruning, exercises direct-repair bounds, verifies
stable JSON/CSV serialization, and covers production GF16 byte-length rejection.

## Limitations

The operation model follows the current algorithms, not every possible backend
micro-op.  It does not count locator-plan setup, field-table initialization,
thread synchronization, pointer validation, range metadata, or NUMA effects.
It does not infer exact zero/one multiplier counts from private skew and direct
coefficient tables; those are explicitly bounded.  The model also assumes the
current deterministic received-subset policy rather than optimizing the subset.

When generated transform schedules become production inputs, their serialized
operation lists should expose exact zero/one multiplier counts.  At that point
this model can consume those schedules and tighten arithmetic bounds without
changing the output schema.
