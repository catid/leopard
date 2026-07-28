# 03 — GF8 out-of-place radix-eight staging

**Disposition: LANDED** (GFNI variant, gated `m == 32`)

## Idea
The staging kernel `AVX2FF8Butterfly4Out`, which reads caller source shards into
the work slots, owns ~40% of a legacy-high GF8 encode call. An isolation screen
showed it runs **within 4% of a kernel that only moves the bytes** at 64 KiB per
shard — it is bound by load/store throughput, not arithmetic. No instruction
change helps a movement floor; moving less does.

A radix-eight group carries three transform layers per load/store round instead
of two: 8 loads + 8 stores for 24 shard-layers against 16 for two radix-four
rounds, a third less traffic per unit of transform work.

## Versus leopard1 (exact main), encode
| Cell | before | after |
| --- | ---: | ---: |
| GF8 K=224 R=32, 64 KiB | ~1.14 | **~1.52** |

Measured directly against the same build without radix-eight: **1.332x encode,
1.129x decode**. Other cells neutral.

## Why register pressure was not the obstacle
Out-of-place means the eight inputs are consumed into registers before any
output is written, so eight data vectors are live rather than sixteen, and every
multiplier matrix is a memory operand — free on a memory-bound kernel. The
earlier rejection of radix-eight (correct: it does not change vector-ALU ops per
byte-layer) simply did not apply to a kernel that is not ALU bound.

## Why the gate is exactly `m == 32`
- `m == 8`: **correctness hazard.** The group consumes every layer, so the
  accumulating final layer that folds a message block into `xor_result` never
  runs and every block after the first is silently dropped. The randomized
  cross-build differential caught this before any timing run.
- `m == 16`: three layers leave one, displacing the fused `dist4 == m`
  accumulating radix-four with a radix-two sweep plus a separate XOR.
  Measured **0.88x-0.94x**.
- `m == 128`: marginal, 1.05x-1.07x with one neutral cell.

## Validation
950 randomized cross-build shapes across 4 seeds, 0 mismatches. Stock suite
100/100. GFNI suite 99/100 (only the expected `leopard2_portable_isa` ISA audit).
