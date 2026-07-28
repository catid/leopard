# 18 — There are two different leopard1 baselines, and they differ by ~50% on encode

**Disposition: METHODOLOGY CORRECTION — no code change; affects how every
vs-leopard1 number in this log should be read**

## The discrepancy
Refreshing report 00's 24-cell headline against the final build produced a
**median encode of 1.011x**, where report 00 records **1.536x**. Decode roughly
agreed (2.737x vs 2.81x). A 50% gap on encode across the same cells is not
noise, and it is not a regression — the two tables use **different leopard1
baselines**.

## Root cause
`CMakeLists.txt` builds `leopard.h` (the original leopard1 API) and
`LeopardFF8.cpp` / `LeopardFF16.cpp` into the **same library target** as
Leopard2. The benchmark's in-process `legacy` section therefore runs

> leopard1's *algorithm* (Algorithm 4/5 scheduling, no reusable erasure plan)
> on **Leopard2's compiled kernels**.

Measured directly — the same leopard1 cells timed in the GFNI build and the
stock build:

| Cell | leopard1 in GFNI build | leopard1 in stock build | leopard1 gain from `-mgfni` |
| --- | ---: | ---: | ---: |
| K=200 R=50, 64 KiB | 776.5 us | 1333.5 us | **1.72x** |
| K=224 R=32, 64 KiB | 517.3 us | 771.4 us | **1.49x** |
| K=1000 R=200, 64 KiB | 6945.1 us | 9836.4 us | **1.42x** |
| K=192 R=64, 64 KiB | 539.9 us | 707.4 us | **1.31x** |

Building the tree with `-mgfni` makes the *baseline* 1.3x-1.7x faster. The
in-process comparison is therefore **same kernels, different orchestration**.

## The two baselines answer different questions

| | Report 00 (`main_compare_screen.py`) | In-process `legacy` |
| --- | --- | --- |
| leopard1 binary | separate build, pinned commit `6e5725eb` | same tree, same flags |
| leopard1 kernels | leopard1-era: no GFNI, no radix-8, no GF16 fusion | **identical to Leopard2's** |
| Measures | total codec advantage, kernels included | **architecture only** |
| Median encode, 24 cells | 1.536x | 1.011x |
| Median decode, 24 cells | 2.81x | 2.737x |

Both are legitimate. Report 00 answers *"how much faster is this codec than the
reference leopard1 implementation"* — the user-facing question. The in-process
run answers *"how much of that is orchestration rather than shared kernel
work"*.

The gap between them is the honest measure of how much the **shared kernel
improvements** (GFNI, GF16 radix-four fusion, radix-eight) contribute: they live
in the shared library, so leopard1's API benefits from them too when built in
the same tree. On encode at 64 KiB that is nearly the entire advantage.

## Consequence for reports 13-16 — they are conservative, not inflated
Every vs-leopard1 figure in reports 13 and 14 was taken with the **in-process**
legacy in the GFNI build, i.e. against a leopard1 that already had GFNI and the
other shared kernel work. Those ratios therefore isolate the layout change and
are measured against the **tougher** of the two baselines:

- GF16 `K=4096 R=512` decode **8.59x** — versus a GFNI-equipped leopard1.
- GF16 encode going from parity (0.99-1.05) to **2.44x-2.54x** at 256 KiB —
  that parity figure is precisely the point: with identical kernels the two
  architectures encode at the same speed at 64 KiB, and the entire separation at
  large shards is the tiling work.

## Refreshed headline, in-process baseline, final build
All 33 cells reported `legacy_comparison: matched`.

| Group | median encode | median decode | n |
| --- | ---: | ---: | ---: |
| Original 24 cells (<= 64 KiB) | 1.011x | 2.737x | 24 |
| Large-shard cells (256 KiB - 1 MiB) | **2.244x** | **5.866x** | 9 |

| Large-shard cell | encode | decode |
| --- | ---: | ---: |
| GF16 K=1000 R=200, 1 MiB | **2.960x** | **9.306x** |
| GF16 K=4096 R=512, 256 KiB | **2.570x** | **9.244x** |
| GF16 K=2000 R=500, 256 KiB | **2.536x** | **7.400x** |
| GF16 K=1000 R=200, 256 KiB | **2.533x** | **8.158x** |
| GF16 K=255 R=129, 256 KiB | **2.244x** | **4.429x** |
| GF16 K=300 R=100, 1 MiB | **2.229x** | **5.866x** |
| GF8 K=224 R=32, 1 MiB | **2.088x** | **5.504x** |
| GF16 K=300 R=100, 256 KiB | **1.965x** | **4.808x** |
| GF16 K=200 R=50, 256 KiB | 1.046x | **5.351x** |

**The original 24-cell table stops at 64 KiB, which is near the worst shard size
for this campaign's work.** At 256 KiB - 1 MiB the same codec is 2-3x on encode
where the 64 KiB table shows parity. Any headline drawn only from the original
cells understates the current codec at realistic shard sizes.

## Rule
**Check what your baseline is compiled from.** A baseline that shares a build
tree with the thing under test will absorb the optimisations you are trying to
measure. Here it moved a headline by 50% and, had it gone unnoticed, would have
looked like a large encode regression.
