# 07 — Fixed per-call cost at tiny K (the only cells still losing)

**Disposition: DIAGNOSED, NOT FIXABLE BY KERNEL WORK**

## Why this report exists
After GFNI, GF16 fusion and radix-eight, exactly two cells in the 24-cell
headline screen still lose to leopard1: `K=2 R=2` (0.606) and `K=4 R=4` (0.851)
at 1 KiB. Every remaining kernel idea was aimed at making the transform faster.
This report establishes that **no kernel change can reach those cells**, so that
future effort is not spent there.

## Method
Sweep each tiny shape across byte sizes and fit `time = fixed + per_KiB * KiB`
separately for exact main and Leopard2. If the loss is transform cost, the
per-KiB slope is worse. If the loss is API cost, the intercept is worse and the
slope is fine. Raw data: `results/fixed_cost_attribution.txt`.

## Result — it is entirely the intercept

| Shape | main fixed | leo2 fixed | delta | main per-KiB | leo2 per-KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| K=2 R=2 | 0.0069 us | 0.0394 us | +0.0325 | 0.0282 us | **0.0137 us** |
| K=4 R=4 | 0.0070 us | 0.0484 us | +0.0414 | 0.0852 us | **0.0779 us** |
| K=8 R=8 | 0.0048 us | 0.0478 us | +0.0430 | 0.2944 us | **0.2902 us** |
| K=16 R=4 | -0.0018 us | 0.0725 us | +0.0743 | 0.3604 us | **0.3015 us** |
| K=16 R=16 | -0.0343 us | 0.0878 us | +0.1220 | 0.8248 us | **0.7821 us** |
| K=32 R=32 | 0.0123 us | 0.2697 us | +0.2574 | 2.2331 us | **1.8934 us** |

**Leopard2's byte-proportional encode cost already equals or beats exact main in
every shape measured.** The entire residual is a 0.03-0.26 us fixed cost that
main does not pay.

The ratio column in the raw data makes the crossover explicit: `K=2 R=2` climbs
0.255 -> 1.287 as bytes go 64 -> 4096, crossing 1.0 near 2 KiB. The cell is only
"losing" because the headline table measures it at 1 KiB.

## What the fixed cost is
It is the price of an API leopard1 does not have: buffer validation walking `K`
and `R` pointers, plan lookup, and dispatch/backend selection. An attribution
build measured encode buffer validation alone at **1.1-1.8 ns per shard
pointer**, and removing validation entirely moved GF8 `K=32 R=32` at 64 bytes
from 0.67 to 0.94 — a large relative move that still does not reach 1.0, and one
that is not available anyway because the validation is a safety contract.

## Conclusion
These two cells are a **measurement-framing artifact plus an API-surface cost**,
not a codec deficiency. Do not spend kernel effort on them. If they matter for a
real workload, the lever is batching (amortising validation and dispatch across
items), not arithmetic.
