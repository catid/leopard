# 13 — GF16 decode workspace tiling

**Disposition: SHIPPED, then SUPERSEDED IN PART by report 15**

> **Read report 15 alongside this.** The `padded_side` 256/512 key promoted here
> was calibrated at 64 KiB shards only, and it excluded sides 64 and 128 at
> *every* shard size. Those sides genuinely fit L3 at 64 KiB but overshoot it
> 4-8x at 1 MiB, where tiling them is worth up to **3.226x** — the best decode
> cell of the whole campaign. Report 15 keeps the two entries below (they catch
> below-threshold re-sweep wins a size-independent rule cannot see) and adds a
> live-set rule covering the region this table never reached. The measurements
> in this report stand; the *policy* is now a union.

## Idea
GF16 legacy-high decode retains `2T` work slots at the *full* aligned shard
prefix. At 64 KiB shards that is a live set of

| `padded_side` | slots | live set at 64 KiB | vs 32 MB L3 |
| --- | ---: | ---: | --- |
| 256 | 512 | 33.6 MB (46.7 MB at max loss) | 1.05x - 1.46x over |
| 512 | 1024 | 67.7 MB (100.0 MB at max loss) | 2.1x - 3.1x over |

Every transform layer therefore streams the entire workspace from DRAM. The
decode executor was *already* field-agnostic — it passes `pass_bytes` straight
into `ExecuteTransformDecodePass` — so the GF8-only restriction in the tile
policy was a calibration boundary, not a correctness one. Splitting the payload
into passes that keep the retained slot set cache-resident should recover that
traffic without touching the field arithmetic or the coordinate schedule.

## The false start that mattered
The first screen came back dead neutral, which looked like another static model
that failed to measure. It was not: **the tile never activated.** The reported
decode scratch was unchanged (33.7 MB in, 33.7 MB out), which is impossible if a
tile had been applied. The global guard rejected the shape before reaching the
side switch —

```c
aligned_prefix_bytes < 128U * 1024U   /* GF8-derived floor */
```

— and a 64 KiB shard has a 64 KiB prefix. **Scratch is the activation witness
for any tiling experiment; a neutral timing result with unchanged scratch is a
non-result, not a negative result.** Every subsequent screen in this campaign
asserted on scratch before reading a ratio.

## Result
Base is the current GFNI + radix-eight build, untiled. Decode execution median
us, best-of-three, one pinned logical CPU (14) of an AMD Ryzen 9 9950X3D,
64 KiB shards, reuse 8, one thread. **All recovered-original digests matched
the untiled build in every cell.**

| Cell | side | blks | base us | scratch | 4 KiB | 8 KiB | 16 KiB | 32 KiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K=255 R=129 l=1 | 256 | 1 | 3490.9 | 33.6 MB | **1.398** | 1.389 | 1.325 | 1.276 |
| K=255 R=129 l=129 | 256 | 1 | 4930.9 | 33.6 MB | **1.349** | 1.329 | 1.296 | 1.301 |
| K=200 R=200 l=1 | 256 | 1 | 3180.1 | 33.6 MB | 1.396 | **1.463** | — | 1.326 |
| K=200 R=200 l=200 | 256 | 1 | 4185.6 | 33.6 MB | **1.416** | 1.384 | — | 1.357 |
| K=250 R=250 l=1 | 256 | 1 | 3354.6 | 33.6 MB | **1.392** | 1.388 | — | 1.264 |
| K=400 R=200 l=1 | 256 | 2 | 4855.7 | 33.6 MB | 1.195 | 1.136 | — | **1.225** |
| K=400 R=200 l=200 | 256 | 2 | 8961.8 | 46.7 MB | 1.272 | 1.280 | — | **1.322** |
| K=1000 R=200 l=1 | 256 | 4 | 8321.7 | 33.7 MB | 1.089 | 1.099 | 1.154 | **1.182** |
| K=1000 R=200 l=200 | 256 | 4 | 14544.9 | 46.7 MB | — | 1.180 | 1.219 | **1.254** |
| K=2000 R=500 l=8 | 512 | 4 | 36791.9 | 67.7 MB | — | 1.913 | **1.913** | 1.538 |
| K=2000 R=500 l=500 | 512 | 4 | 49014.9 | 100.0 MB | — | 1.746 | **1.766** | 1.458 |
| K=4096 R=512 l=8 | 512 | 8 | 58140.7 | 67.8 MB | — | — | **1.824** | 1.521 |
| K=200 R=50 l=8 *(control)* | 64 | — | 1550.3 | 8.9 MB | — | — | 0.998 | 1.016 |

`padded_side == 64` is an untouched control: its scratch is unchanged and it
measures neutral, confirming the guard did not widen beyond the two calibrated
sides.

## Two regimes, and why the shipped constant is neither optimum
The sweep splits cleanly on **message-block count** (`original_count` vs
`padded_side`) — the same structural predicate this file already uses for GF8
encode tiling:

- **1 block** (`K <= side`): small tiles win, monotonically. 32 KiB is the worst
  tested tile in all five cells.
- **2+ blocks**: 32 KiB wins at side 256, and the preference *strengthens* with
  block count (2 blocks: 1.225 vs 1.195; 4 blocks: 1.182 vs 1.089). A multi-block
  traversal re-walks the retained workspace once per message block, so a larger
  tile amortises the re-walk; a single-block traversal has no re-walk to amortise.
- **side 512** is flat between 8 and 16 KiB and falls off hard at 32 KiB, where
  the live set (34 MB) re-crosses L3.

The shipped policy is a **uniform 16 KiB** for GF16 legacy-high AVX2 tiled
decode at sides 256 and 512, floored at 64 KiB shards. It is deliberately
neither per-regime optimum:

- it is within about 5% of the best measured tile in every one of the 13 cells;
- it takes the entire side-512 win (1.77x - 1.91x), which is the bulk of the value;
- it needs no extrapolation into (side, block-count) combinations that were
  never measured — notably single-block at side 512.

Shipping the per-regime table would be precisely the overfit that report 12
declined to commit off a four-cell screen. The per-cell optima are recorded
above so a future offline sweep across `padded_side` x shard size x block count
can promote them on real evidence.

Scratch falls with the tile: 67.7 MB -> 17 MB at K=2000 R=500, 33.6 MB -> 8 MB
at side 256.

## Versus leopard1 — the shipped build
Measured on the promoted 16 KiB build, not the sweep. `vs L1` is leopard1
decode time (which includes its locator setup, since leopard1 has no reusable
erasure plan) divided by Leopard2's decode amortised at reuse 8; above 1.0
favours Leopard2. **`correctness.legacy_comparison` reported `matched` for both
builds in every cell** — that is a direct wire comparison against leopard1's
output, stronger than a self-consistent digest.

| Cell | side | untiled us | tiled us | speedup | scratch | vs L1 before | vs L1 after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K=2000 R=500 l=8 | 512 | 36409.7 | 19088.5 | **1.907** | 67.7 -> 17.0 MB | 3.36 | **6.42** |
| K=4096 R=512 l=8 | 512 | 57425.3 | 31600.5 | **1.817** | 67.8 -> 17.1 MB | 4.70 | **8.59** |
| K=2000 R=500 l=500 | 512 | 48329.5 | 27353.0 | **1.767** | 100.0 -> 25.1 MB | 2.85 | **5.04** |
| K=200 R=200 l=1 | 256 | 3199.7 | 2212.3 | 1.446 | 33.6 -> 8.4 MB | 2.19 | **3.17** |
| K=200 R=200 l=200 | 256 | 4219.0 | 2949.0 | 1.431 | 33.6 -> 8.4 MB | 1.78 | **2.54** |
| K=250 R=250 l=1 | 256 | 3379.1 | 2431.2 | 1.390 | 33.6 -> 8.4 MB | 2.20 | **3.16** |
| K=255 R=129 l=1 | 256 | 3452.0 | 2541.8 | 1.358 | 33.6 -> 8.4 MB | 2.07 | **2.88** |
| K=400 R=200 l=200 | 256 | 9030.0 | 6772.8 | 1.333 | 46.7 -> 11.7 MB | 2.52 | **3.34** |
| K=255 R=129 l=129 | 256 | 4874.6 | 3677.7 | 1.325 | 33.6 -> 8.4 MB | 1.66 | **2.12** |
| K=400 R=200 l=1 | 256 | 4858.7 | 3975.9 | 1.222 | 33.6 -> 8.4 MB | 4.18 | **5.07** |
| K=1000 R=200 l=200 | 256 | 14279.0 | 11849.1 | 1.205 | 46.7 -> 11.7 MB | 4.38 | **5.30** |
| K=1000 R=200 l=1 | 256 | 8375.0 | 7237.2 | 1.157 | 33.7 -> 8.5 MB | 6.25 | **7.11** |
| *K=200 R=50 l=8* | *64* | *1501.5* | *1496.6* | *1.003* | *8.9 -> 8.9 MB* | *4.40* | *4.40* |
| *K=300 R=100 l=8* | *128* | *2710.7* | *2705.6* | *1.002* | *17.3 -> 17.3 MB* | *2.75* | *2.75* |
| *K=224 R=32 l=8 (GF8)* | *32* | *1012.6* | *1010.6* | *1.002* | *16.8 -> 16.8 MB* | *2.03* | *2.03* |

**Speedup over the twelve targeted cells: min 1.157, median 1.333, max 1.907.**

The three italic rows are controls at `padded_side` 32, 64 and 128 — two GF16
and one GF8. All three are neutral to within 0.3% **with scratch byte-identical**,
which is the direct evidence that the new guard admitted exactly the two
calibrated sides and widened nothing.

Headline: **GF16 `K=4096 R=512` decode goes from 4.70x to 8.59x versus
leopard1**, and `K=2000 R=500` from 3.36x to 6.42x. These are the cells where
the untiled workspace was 2-3x larger than L3, so they had the most DRAM traffic
to recover.

## Validation

**Randomized differential, untiled vs tiled build: 260 shapes, 0 mismatches.**
Coverage by `padded_side`: `{1:25, 2:15, 4:13, 16:9, 32:14, 64:13, 128:26,
256:75, 512:25}`, of which **67 shapes activated the new tile**. Compared
original/parity/recovered digests, the round-trip flag, and shape acceptance
across GF8/GF16, scalar/SSSE3/AVX2/auto, 1-2 threads, batch 1-2, and shard sizes
straddling the 64 KiB floor (4 KiB, 32 KiB, 65472, 65536, 65600, 98304, 128 KiB).
Tiling changes only how the payload is split into passes, so any digest
difference would be a break; there were none.

**Suites** (stock and GFNI, run serially on a quiet machine):

- `leopard2_field_options_matrix` — load-bearing here, because the change
  restructured the tile helper's `#ifdef` nesting from one whole-body
  `LEO_HAS_FF8` guard into separate per-field blocks. This is the only gate that
  compiles the `NO_LEO_HAS_FF8` / `NO_LEO_HAS_FF16` combinations.
- `leopard2_public_api_contract` — decode scratch-size slope assertions.
- `leopard2_operation_counts_self_test` and `leopard2_decode_scratch_crosscheck`
  **failed and were fixed properly.** `tools/leopard2_operation_counts.py`
  name-pins the tile helper in its schema-v2 source model and documents its
  scope as GF8-only; both needed updating for the rename to
  `AVX2DecodeExecutionTileBytes` and the widened field coverage. The contract
  did its job — it is designed to fail when the tile policy's scope changes.

**Two failures remain, both verified as non-regressions:**

- `leopard2_portable_isa` in the GFNI build — fails *identically in the
  pre-change untiled GFNI build*, because the experimental variant passes
  `-mgfni` globally and `vgf2p8affineqb` lands in `Leopard2BackendAVX2.cpp.o`.
  This is the standing GFNI shipping blocker recorded in
  `docs/leopard2_gfni_codec.md`; it needs a separate translation unit with
  runtime CPUID dispatch. The stock build passes.
- `leopard2_lab_self_test` — "oversubscription evidence did not record the
  team". Passes in isolation and passed in the first run; a load-sensitive flake
  under `ctest -j12`.
