# 20 — High-profile AUTO direct encode

**Disposition: INVALID MEASUREMENT — no performance conclusion. The experiment
remains default-off while a corrected explicit-AVX2 ABBA campaign is run.**

## 2026-07-28 audit correction

The screen below did not measure encode execution.  Its ad hoc recursive JSON
extractor selected the first key containing `median`; in the benchmark schema
that is `metrics.codec_setup.median_us`, which precedes
`metrics.encode_execution.median_us_per_batch_call`.  Test-hook builds prepare
the bounded direct rows regardless of the forced execution mode, so the
near-1.0 ratios are expected setup ratios and say nothing about the direct
encoder's execution ceiling.

The screen also omitted the earlier `K=2,R=16,Q=1,4 KiB` cell that had measured
17.88x, used mutable AUTO-variant build outputs rather than lane-owned frozen
binaries, had no authoritative ABBA confidence interval, and compared empty
outer digest objects.  The benchmark itself checked parity internally, but the
screen's claimed independent digest comparison was vacuous.  Therefore neither
the 0.875x-1.017x range nor the explanation that later transform work consumed
the direct-path margin is valid evidence.  Bead `leopard-79h.42.1` tracks the
corrected campaign.

## The candidate
Open finding 5 and the 2026-07-16 direct-encode checkpoint held a striking
number: the `excluded_high_profile` region — high-profile cells inside the
existing K,R <= 16 direct-encode bounds — measured **median +541%, max +1688%**
over six cells, blocked from production AUTO only by
`CanAutoDirectEncodeCodec`'s `profile != LEO2_PROFILE_LOW_V1` test. The one
measured regression (scalar K=3 R=1, -19.6%) was cleanly excludable.

Two scope facts sharpened it. First, the high-profile direct schedule is
already prepared and differentially tested in every test-hook build — only
production AUTO caps it. Second, `AutoDirectEncodePreferred` requires
`requested_recovery_count == 1`: AUTO direct encode is a **single-parity-
request (Q=1)** optimization. A first screen at full parity measured dead
neutral for exactly this reason — the checkpoint's `excluded_q2` cell had
already shown Q=2 regressing -30%, and full-parity requests never enter the
path at all.

## The invalid screen (retained for auditability)
`bench_leopard2_direct_encode` (the purpose-built harness with per-path
introspection), Q=1, high profile, GF8, AVX2, nine cells over K,R <= 16 and
1 KiB-64 KiB. The decisive column is the **ceiling**: force-direct versus
force-transform on the same experiment build, which measures the direct
encoder's intrinsic margin with no dispatch policy in the way.

The following values are codec-setup timings mislabeled as execution timings.
They must not be used for dispatch or optimization decisions.

| Cell | base setup | exp setup | setup ratio | setup ratio (forced modes) |
| --- | ---: | ---: | ---: | ---: |
| K=2 R=1, 1 KiB | 0.070 us | 0.070 us | 1.000x | 1.000x |
| K=3 R=1, 4 KiB | 0.080 | 0.080 | 1.000x | 0.875x |
| K=4 R=4, 1 KiB | 0.160 | 0.180 | 0.889x | 0.944x |
| K=8 R=8, 1 KiB | 0.380 | 0.420 | 0.905x | 1.000x |
| K=16 R=16, 1 KiB | 1.210 | 1.300 | 0.931x | 0.992x |
| K=16 R=16, 4 KiB | 1.210 | 1.300 | 0.931x | 1.000x |
| K=8 R=8, 64 KiB | 0.390 | 0.410 | 0.951x | 1.000x |
| K=16 R=4, 4 KiB | 1.160 | 1.200 | 0.967x | 1.017x |
| K=16 R=16, 64 KiB | 1.920 | 1.300 | 1.477x* | 0.985x |

\* The setup sample is noisy.  No execution-path inference is valid for this
row or any other row in the table.

No encode-execution comparison follows from this table.  It neither confirms
nor refutes whether later transform changes reduced the older direct-path
advantage.

## What shipped anyway
- The experiment gate (`LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE`) remains
  default-off pending valid current-tree execution evidence.
- A real latent fix in `bench_leopard2_direct_encode`'s policy mirror
  (`ExpectedAutoDirectPath`): it did not know `LEO2_BACKEND_GFNI`, so AUTO
  direct-encode benchmarking on an explicit GFNI context would have failed its
  own introspection check — a promotion gap the GFNI test sweep missed because
  this mirror lives in the bench, not the test matrix.
- Verification that the harness's self-check works as designed: it failed the
  experiment build precisely because the library's actual path selection
  diverged from the bench's policy model, i.e. it caught a real policy change
  the moment one occurred.

## Conditions to revisit
Re-establish the force-direct execution ceiling on exact historical cells and
neighbors using explicit AVX2, exact schema paths, frozen binaries, non-vacuous
parity comparison, isolated ABBA rounds, and confidence intervals.  AUTO
admission remains out of scope until that corrected ceiling is comfortably
above 1.0 across a deterministic region.
