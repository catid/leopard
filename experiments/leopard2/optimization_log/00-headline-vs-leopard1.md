# 00 — Headline: current build versus leopard1

Refreshed after the inverse radix-eight landed. Ratios are exact-main
(leopard1, commit `6e5725eb`) time divided by Leopard2 time, so **above 1.0
favours Leopard2**. Clustered ABBA, one pinned logical CPU of an AMD Ryzen 9
9950X3D, three rounds per cell, five samples per child, reuse 8, one thread.
All three workload digests matched exact main in every cell.

Screen, not promotion evidence: no CPU-pair lease, no build-closure
verification, no sibling-idleness gate. The isolated
`experiments/leopard2/main_compare/run_abba.py` campaign is still required
before any of this enters `docs/leopard2_vs_main_benchmark.md`.

## Encode, by cumulative idea

| Cell | stock | +GFNI | +GF16 fusion | +radix-8 | dec stock | dec now |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GF16 K=200 R=50, 64 KiB | 1.194 | 1.577 | 1.764 | **1.815** | 3.10 | **4.62** |
| GF8 K=224 R=32, 64 KiB | 1.134 | 1.338 | 1.361 | **1.788** | 1.46 | **1.72** |
| GF16 K=1000 R=200, 4 KiB | 1.091 | 1.555 | 1.781 | **1.766** | 2.28 | **3.66** |
| GF16 K=2000 R=500, 4 KiB | 1.097 | 1.520 | 1.729 | **1.735** | 1.97 | **2.81** |
| GF16 K=1000 R=999, 4 KiB | 1.034 | 1.560 | 1.689 | **1.685** | 1.67 | **2.41** |
| GF16 K=4096 R=512, 4 KiB | 1.119 | 1.462 | 1.695 | **1.679** | 3.42 | **4.87** |
| GF8 K=192 R=64, 64 KiB | 1.093 | 1.515 | 1.609 | **1.662** | 1.39 | **1.82** |
| GF16 K=255 R=129, 4 KiB | 1.003 | 1.385 | 1.616 | **1.633** | 2.04 | **2.82** |
| GF8 K=100 R=100, 64 KiB | 1.070 | 1.507 | 1.610 | **1.596** | 1.27 | **1.46** |
| GF8 K=128 R=128, 64 KiB | 1.070 | 1.598 | 1.625 | **1.576** | 1.25 | **1.45** |
| GF16 K=300 R=100, 64 KiB | 1.077 | 1.374 | 1.554 | **1.567** | 2.01 | **2.94** |
| GF16 K=1000 R=200, 64 KiB | 1.128 | 1.382 | 1.547 | **1.559** | 4.16 | **6.13** |
| GF8 K=240 R=16, 4 KiB | 1.168 | 1.540 | 1.529 | **1.512** | 3.16 | **4.02** |
| GF8 K=64 R=64, 64 KiB | 1.141 | 1.506 | 1.542 | **1.505** | 1.35 | **1.50** |
| GF8 K=128 R=128, 4 KiB | 0.972 | 1.549 | 1.551 | **1.505** | 1.24 | **1.55** |
| GF16 K=255 R=129, 64 KiB | 0.946 | 1.240 | 1.464 | **1.470** | 1.67 | **2.17** |
| GF8 K=240 R=16, 64 KiB | 1.226 | 1.423 | 1.414 | **1.402** | 3.03 | **3.48** |
| GF8 K=16 R=16, 4 KiB | 1.019 | 1.307 | 1.318 | **1.314** | 1.27 | **1.59** |
| GF8 K=129 R=1, 64 KiB | 1.152 | 1.169 | 1.141 | **1.165** | 1.15 | **1.14** |
| GF8 K=129 R=1, 4 KiB | 1.008 | 1.045 | 1.071 | **1.147** | 1.02 | **1.16** |
| GF8 K=8 R=8, 1 KiB | 0.885 | 1.146 | 1.121 | **1.143** | 4.42 | **5.76** |
| GF8 K=16 R=4, 1 KiB | 0.823 | 1.028 | 1.077 | **1.062** | 4.16 | **5.46** |
| GF8 K=4 R=4, 1 KiB | 0.682 | 0.853 | 0.859 | **0.851** | 5.27 | **6.44** |
| GF8 K=2 R=2, 1 KiB | 0.574 | 0.577 | 0.603 | **0.606** | 5.98 | **7.38** |

## Summary

- **Median encode versus leopard1: 1.074 -> 1.536.** The codec as a whole is a
  median **1.46x faster than it was** against the same baseline.
- **Median decode versus leopard1: 1.99 -> 2.81.**
- Best encode cell: GF16 K=200 R=50 at 64 KiB, **1.815x**.
- Best decode cell: GF8 K=2 R=2 at 1 KiB, **7.38x** (leopard1 has no reusable
  erasure plan, so its decode always re-runs locator setup).

## The only two cells still losing

`K=2 R=2` (0.606) and `K=4 R=4` (0.851) at 1 KiB. These are **not**
transform-bound: the whole encode call is ~70 ns against exact main's ~40 ns.
A size sweep separating fixed from byte-proportional cost showed Leopard2's
byte-proportional encode cost already equals or beats exact main in every shape
measured (K=2 R=2: 0.0137 vs 0.0282 us/KiB). The entire residual is fixed
per-call cost — buffer validation walks over K and R, plan and dispatch
selection — which is the price of an API leopard1 does not have. An attribution
build measured encode buffer validation at 1.1-1.8 ns per shard pointer, and
removing it entirely moves GF8 K=32 R=32 at 64 bytes from 0.67 to 0.94.

No kernel change reaches these cells. See report 07.
