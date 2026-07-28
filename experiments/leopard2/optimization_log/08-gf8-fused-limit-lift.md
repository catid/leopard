# 08 — GF8 fused-limit lift

**Disposition: REJECTED — encode neutral**

## Idea
The direct analogue of report 02. GF16 carried a size cap on radix-four
butterfly fusion that turned out to be stale once GFNI made the multiply cheap,
and lifting it paid a median +13.8% across nine GF16 cells. GF8 carries its own
fused-limit constant, so the obvious next move was to lift it the same way.

The shard-touch model predicts the same shape of win: fusing two layers halves
the number of times the payload is touched.

## Result
**Encode neutral across the GF8 cells.** The GF8 column of the same 24-cell
screen that showed the GF16 win moved by a median of 1.7% with four small
regressions (`results/main_compare_screen_gfni_fused.json`), and lifting the
GF8-specific limit on top of that did not move it further.

## Why it failed where GF16 succeeded — the temperature discount
This is the clearest case in the campaign of the touch model being right about
the *count* and wrong about the *price*.

A GF8 element is one byte where a GF16 element is two, so at equal shard count
and equal shard size a GF8 transform's working set is half the size. The GF8
cells in the screen sit at `padded_side` 16-128 with 64 KiB shards: their live
transform rows fit in L2. Saving a touch on **warm, cache-resident** data is
worth ~0.25 us in the profiled cell; saving one on **cold, DRAM-streamed** data
is worth ~0.55 us. GF16 at side 256-512 is squarely in the cold regime — which
is also exactly why report 13's decode tiling later paid 1.15-1.91x on the same
geometry.

**Rule extracted:** before promising a win from fewer payload touches, compute
the live working set and compare it to L2/L3. Touch reduction pays roughly 2x
more on cold data, and near nothing when the data was already resident. This one
fact explains the rejections in reports 06, 08, 09 and 11.

**To revisit:** GF8 at large `padded_side` with large shards — the regime where
its working set does leave cache — was not swept. The cap is unmeasured there,
not shown neutral.
