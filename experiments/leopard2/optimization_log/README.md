# Leopard2 AVX2 optimization log

One markdown file per idea, kept so a rejected idea can be revisited with its
original numbers instead of being re-derived or re-tried blindly.

Every report states the effect **versus exact Leopard main (leopard1)** where
that comparison is possible, and versus stock Leopard2 otherwise.

**Three** standing caveats on the leopard1 column:

- **There are two different leopard1 baselines and they differ by ~50% on
  encode.** `main_compare_screen.py` uses a separately-built leopard1 pinned at
  commit `6e5725eb`; the benchmark's in-process `legacy` section runs leopard1's
  algorithm on **Leopard2's own compiled kernels**, because `leopard.h` and
  `LeopardFF8/FF16.cpp` share a library target. Building with `-mgfni` therefore
  makes the in-process baseline 1.3x-1.7x faster. Report 00 uses the external
  baseline (measures total codec advantage); reports 13-16 use the in-process
  one (measures orchestration only, and is the tougher baseline). See report 18
  before comparing numbers across reports.

- leopard1 rejects `R > K`, so **no low-rate (`LOW_V1`) cell has a leopard1
  ratio**. Low-rate ideas are reported against stock Leopard2 only.
- All ratios here come from the clustered-ABBA screen
  (`experiments/leopard2/gfni_codec/main_compare_screen.py`), pinned to one
  logical CPU with the machine quiet. That is a **screen, not promotion
  evidence**: it does not hold the CPU-pair lease, verify build closure, or gate
  on sibling idleness. The isolated
  `experiments/leopard2/main_compare/run_abba.py` campaign is still required
  before any number here enters `docs/leopard2_vs_main_benchmark.md`.

Ratios above 1.0 favour Leopard2.

## Current position

**Median encode versus leopard1: 1.074 -> 1.536** (a 1.46x improvement of the
codec against a fixed baseline). Median decode 1.99 -> 2.81. Best encode cell
1.815x, best decode cell 7.38x. Two cells still lose, both 1 KiB tiny-K, both
fixed per-call cost rather than transform cost. Full table: `00-headline-vs-leopard1.md`.

Reports 13, 14, 15 and 16 all landed after that table was taken. They share one
root cause — a full-payload transform workspace held against a 32 MB L3 — and
together they move the GF16 legacy-high cells substantially. All ratios below
were verified against leopard1's actual output (`legacy_comparison: matched`).

**Decode versus leopard1, 64 KiB shards** (report 13):

| Cell | before | after |
| --- | ---: | ---: |
| GF16 K=4096 R=512 | 4.70 | **8.59** |
| GF16 K=1000 R=200 | 6.25 | **7.11** |
| GF16 K=2000 R=500 | 3.36 | **6.42** |
| GF16 K=2000 R=500, max loss | 2.85 | **5.04** |
| GF16 K=200 R=200 | 2.19 | **3.17** |

**Encode versus leopard1, explicit AVX2 both sides** (report 14). These cells
were at *parity* before tiling:

| Cell | shard | before | after |
| --- | --- | ---: | ---: |
| GF16 K=4096 R=512 | 256 KiB | 0.999 | **2.541** |
| GF16 K=1000 R=200 | 256 KiB | 0.990 | **2.522** |
| GF16 K=2000 R=500 | 256 KiB | 1.007 | **2.470** |
| GF16 K=2000 R=500 | 64 KiB | 0.998 | **1.731** |

**Live-set rule** (report 15) then covered the sides and shard sizes both tables
excluded, against the untiled build: decode up to **3.226x** and encode up to
**2.735x**, with the best legacy-high cell — GF16 `K=300 R=100` decode at
1 MiB — at a `padded_side` of 128 that every earlier table skipped.

**Low-rate profile** (report 16). Both guards required `LEGACY_HIGH_V1`, so
`LOW_V1` never tiled at any shard size — an exclusion with no measurement behind
it, inherited from a campaign scoped to legacy-high. Low-rate has **no leopard1
comparison at all** (leopard1 rejects `R > K`), so this is versus stock
Leopard2:

| Cell | side | encode | decode |
| --- | ---: | ---: | ---: |
| GF16 K=256 R=768, 1 MiB | 256 | 2.515x | **3.099x** |
| GF16 K=128 R=896, 1 MiB | 128 | 2.365x | **2.954x** |
| GF16 K=200 R=800, 256 KiB | 256 | 2.293x | **2.589x** |
| GF16 K=1000 R=2000, 64 KiB | 1024 | 2.328x | **2.159x** |

Scratch falls with all of it: encode 1073.8 MB -> 8.5 MB (126x) at legacy-high
`K=2000 R=500` with 1 MiB shards; low-rate 537 MB -> 17 MB; decode 100 MB ->
25 MB at maximum loss.

`00-headline-vs-leopard1.md` predates all of these and is not restated.

## Index

| # | Idea | Disposition | Best effect |
| --- | --- | --- | --- |
| 01 | GFNI affine multiplication | **landed** | 1.15-1.67x whole-codec |
| 02 | GF16 radix-four fusion at all sizes | **landed** | +13.8% median, 9 GF16 cells |
| 03 | GF8 out-of-place radix-eight staging | **landed** | 1.33x at K=224 R=32 |
| 04 | GF8 forward radix-eight (low-rate) | **landed** | 1.14-1.23x low-rate |
| 05 | Batch preflight sort | **landed** | 1.36-1.76x batch |
| 06 | Cache blocking, two forms | rejected | neutral |
| 07 | Fixed per-call cost at tiny K | diagnosed, not kernel-fixable | slope already wins |
| 08 | GF8 fused-limit lift | rejected | encode neutral |
| 09 | Dispatch-bound unrolling | rejected | 1.01x |
| 10 | Non-temporal stores | rejected | 2-3x worse |
| 11 | In-place radix-eight | rejected | neutral |
| 12 | GF8 encode byte tiling at 64 KiB | recorded, not shipped | +4-5%, sharp optimum |
| 13 | GF16 decode workspace tiling | **landed** | 1.15-1.91x decode |
| 14 | GF16 encode workspace tiling (AVX2) | **landed** | up to 2.92x encode, 126x less scratch |
| 15 | Live-set tile rule replacing side tables | **landed** | up to 3.38x decode, 2.74x encode |
| 16 | Low-rate profile tiling (GF16 + GF8, both paths) | **landed** | up to 3.15x decode, 2.86x encode |
| 17 | Thread-scaled live-set target | refuted, not shipped | no effect; saved a knob |
| 18 | Two different leopard1 baselines | methodology correction | ~50% encode gap explained |
| 19 | GFNI production backend member | **landed** | vs leopard1: enc 1.10->1.71, dec 3.42->5.00 median |
| 20 | High-profile AUTO direct encode | **invalid screen; rerun open** | setup median was mislabeled as encode execution |
| 21 | AUTO-GFNI calibrated-host candidate | mechanism shipped, default-off | 44-cell screen: 0 regressions |
| 22 | Dead block-0 pruned schedule (finding 1) | **landed** | plan setup 1.06-1.15x |
| 23 | Writable-partitioned batch preflight | **landed** | 1.01-1.06x batch; model gap recorded |
| 24 | GF16 affine table packing | **landed** | 8 MiB -> 2 MiB, perf neutral, requirement 2 done |
| 25 | Affinity-aware GF16 L3 tiling | **landed** | avoids 20.7-26.1% large-cache regressions; retains up to 1.11x |
| 26 | Stack-owned transient one-shot plan | rejected and reverted | target aggregate 0.988x; two regressions |
| 27 | Scratch-resident native-high one-shot Algorithm 5 | **landed** | 1.188x vs main at 64 B; exposed 65-B crossover |
| 28 | Whole-shard native-high locator-direct repair | **landed** | 1.289-4.513x vs exact main; 65 B now 2.228x |
| 29 | Native-high whole-direct crossover extension | **landed through 7168 B** | selected 1.124-4.972x vs prior; 1.886-5.133x vs main |
| 30 | Native-high source-major raw one-shot | **landed at 12-16 KiB** | 1.128x overall vs plan fallback; 3.484x overall vs exact main |

## Method notes worth keeping

- The **shard-touch model** predicts these results to ~1.4%: a transform stage
  costs `2m` touches regardless of radix, so `L = log2(m)` layers cost
  `2m*ceil(L/2)` with radix-four rounds and `2m*ceil(L/3)` with radix-eight.
  Use it to size a candidate.
- Discount it by **temperature**: cold streaming touches cost ~0.55 us each in
  the profiled cell, warm in-place touches ~0.25 us. Saving touches only pays on
  cold data. This single fact explains every rejection in 06, 11 and 12.
- Five separate static models predicted large wins that measured flat. Every
  idea that paid was validated on a standalone kernel screen **before**
  integration. Measure first.
- **Never screen a replacement policy only at the parameter value the old
  policy was pinned to.** Report 14's predecessor tiled GF16 encode at exactly
  64 KiB, so 64 KiB was the natural place to measure its replacement — and it is
  the *worst* size for the change: 1.14x at 64 KiB versus 2.92x at 1 MiB. The
  pinned value is the least informative point in the sweep.
- **When a guard is keyed on a shape parameter, ask what physical quantity it
  is a proxy for, and whether that proxy survives the range you did not sweep.**
  This happened three times in one campaign (report 15): an exact-shard-size
  bound, two `padded_side` tables, and a message-block exclusion were all
  proxies for `2 * padded_side * shard_bytes` versus L3, and all three broke
  outside the swept range. Keying on the live set directly predicted an
  untested region correctly before it was measured.
- **Verify the change actually engaged before believing a neutral result.**
  Report 13's first screen measured dead neutral and looked like another failed
  model. It was not — the tile never activated, because a global 128 KiB
  aligned-prefix floor rejected the 64 KiB shards before the new policy was
  reached. The tell was that reported scratch was *unchanged*, which is
  impossible if a tile had applied. Once fixed, the same idea measured up to
  1.91x. Find a side-channel that must move if the code path ran — scratch
  bytes, a counter, a size — and assert on it. **A neutral result with an
  unchanged witness is a non-result, not a negative result.**
- **Gate a timing run on the test log printing its summary line, not on a
  process poll with a retry bound.** Two measurements were discarded this
  session because a bounded wait loop timed out and proceeded while a ctest
  suite was still running, inflating base times 2-4x. Check the machine is quiet
  *before* timing, not after the numbers look strange.
- **A guard relaxation can silently widen a policy you were not editing.**
  Report 16 relaxed one profile test at the top of the decode tile helper to
  admit low-rate to the GF16 branch. That also handed low-rate to the *GF8*
  table further down the same function, producing a real but uncalibrated 2.540x
  that looked like a free win. It surfaced only because the screen carried a GF8
  control cell whose decode moved under an ostensibly GF16-only edit. Keep a
  control in every screen for the shapes the change is *not* supposed to touch,
  and check them before reading the wins.
- **Check what your baseline is compiled from.** A baseline sharing a build
  tree with the thing under test absorbs the optimisations you are measuring.
  Report 18: the in-process leopard1 runs on Leopard2's kernels, so `-mgfni`
  sped the *baseline* up 1.3x-1.7x and moved a headline by 50%. Unnoticed, it
  would have read as a large encode regression.
- **Extract benchmark fields by an exact, versioned schema path.** Report 20's
  recursive “first median” lookup selected codec setup instead of encode
  execution, producing a plausible near-1.0 table and a false refutation.
  Require exact metric paths, nonempty digest identities, immutable binary
  provenance, and adversarial parser tests before interpreting a re-screen.
- **A mechanism that predicts an effect is not evidence the effect exists.**
  Report 17's thread-scaled tile target had a clean mechanical story (T threads
  share L3, so the aggregate live set is T * target), a knob the context could
  already drive, and a plausible case. A three-point screen showed shrinking the
  target eightfold changes nothing — the high-thread fall-off is DRAM bandwidth,
  not L3 contention. Screening the knob before building it turned a feature into
  a comment. A refuted hypothesis is maintenance you never take on.
- **A `SKIP` in a sweep is missing coverage, not a neutral result.** Report 16
  had two cells fail to construct, which would have left an entire `padded_side`
  untested while the results table looked complete. Reproducing them by hand
  showed both builds rejected them identically — a pre-existing API constraint
  (`LOW_V1` needs `R` as a multiple of `P`), not the change under test — and
  picking a conforming shape covered the side properly.
- When a contract test blocks a change, derive what its counters actually count
  before restating it. Report 04 records four attempts; the first three were
  guesses and the fourth, which worked, came from reading `fetch_add(dist)` at
  the increment site.
