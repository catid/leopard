# 02 — GF16 radix-four fusion at all sizes

**Disposition: LANDED**

## Idea
The GF16 transform capped radix-four butterfly fusion below a size threshold,
falling back to pairs of radix-two layers. That cap was calibrated when the
per-butterfly arithmetic was nibble-table lookups. GFNI (report 01) made the
multiply dramatically cheaper, which shifts the balance: with cheap arithmetic
the win from fusing two layers into one pass — halving the number of times the
shard payload is touched — is no longer offset by the fused kernel's higher
register pressure.

This is the shard-touch model applied directly: a transform stage costs `2m`
touches regardless of radix, so `L` layers cost `2m*ceil(L/2)` fused versus
`2m*L` unfused. The cap was a good candidate for being **stale** rather than
wrong, which is the same reasoning that later paid off for the GF16 decode
workspace (report 13).

## Result
Gain is the fused build's leopard1 encode ratio divided by the unfused GFNI
build's, over the 24-cell clustered-ABBA screen. Raw data:
`results/main_compare_screen_gfni.json` and
`results/main_compare_screen_gfni_fused.json`.

| Cell | field | GFNI | +fused | gain |
| --- | --- | ---: | ---: | ---: |
| K=255 R=129, 64 KiB | gf16 | 1.240 | 1.464 | **+18.1%** |
| K=255 R=129, 4 KiB | gf16 | 1.385 | 1.616 | **+16.6%** |
| K=4096 R=512, 4 KiB | gf16 | 1.462 | 1.695 | **+16.0%** |
| K=1000 R=200, 4 KiB | gf16 | 1.555 | 1.781 | **+14.5%** |
| K=2000 R=500, 4 KiB | gf16 | 1.520 | 1.729 | **+13.8%** |
| K=300 R=100, 64 KiB | gf16 | 1.374 | 1.554 | **+13.1%** |
| K=1000 R=200, 64 KiB | gf16 | 1.382 | 1.547 | **+12.0%** |
| K=200 R=50, 64 KiB | gf16 | 1.577 | 1.764 | **+11.8%** |
| K=1000 R=999, 4 KiB | gf16 | 1.560 | 1.689 | **+8.3%** |
| K=100 R=100, 64 KiB | gf8 | 1.507 | 1.610 | +6.8% |
| K=192 R=64, 64 KiB | gf8 | 1.515 | 1.609 | +6.2% |
| K=129 R=1, 64 KiB | gf8 | 1.169 | 1.141 | -2.4% |
| K=8 R=8, 1 KiB | gf8 | 1.146 | 1.121 | -2.2% |

**GF16 median +13.8% over 9 cells. GF8 median +1.7% over 15 cells.**

## Reading the result
The effect is cleanly field-scoped, which is what makes it credible rather than
noise: the change targets the GF16 fusion cap, and every one of the nine GF16
cells improved by 8-18% while the fifteen GF8 cells moved by a median 1.7%
with four small regressions. The GF8 column is the control, and it behaves like
one.

The four GF8 regressions (-0.6% to -2.4%) are all cells where the transform is
not the dominant cost — `K=129 R=1` is a one-parity shape and `K=8 R=8` at
1 KiB is dominated by fixed per-call cost (see report 07). They sit inside the
round-to-round spread of the screen.

## Caveat
These are ratios against leopard1 from the clustered-ABBA **screen**, not the
isolated promotion campaign. See the README's standing caveats.
