# 06 — Cache blocking the GF16 radix-four schedules

**Disposition: REJECTED — measured neutral**

## Idea
Two restructurings, both aimed at keeping a radix-four group cache-resident:
1. Byte-block `AVX2FF16Butterfly4Split` into 8 KiB chunks so the four-shard
   group stays in L1 across all four radix-two layers.
2. Chunk `AVX2FF16Butterfly4Range`'s coordinate range so all four layers run
   over a group of coordinates before advancing, targeting a 256 KiB live set.

## Result
Encode ratios versus the same build without the change, twelve GF8/GF16 cells,
clustered ABBA:
- form 1: **0.980-1.015** (neutral to slightly negative)
- form 2: **0.985-1.005**

No leopard1 column: the change measured neutral against Leopard2 itself, so it
cannot move the leopard1 ratio.

## What it taught
Form 1 also exposed a targeting error worth remembering: `UseFusedButterfly4`
returns false above 128 bytes, so `AVX2FF16Butterfly4Split` **is not on the hot
path** for 4 KiB or 64 KiB shards. Form 2 targeted the path that is, and was
still neutral.

At these shard sizes the GF16 AVX2 transform is ALU/port bound, not
L1/L2-miss bound; prefetchers and the out-of-order window already absorb the
layer-major reuse distance. **Do not re-propose blocking here without a profile
showing stalls.**
