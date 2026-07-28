# 20 — High-profile AUTO direct encode

**Disposition: REFUTED ON THE CURRENT TREE — the campaign's transform work ate
the margin. Default-off experiment retained as the record.**

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

## The screen that refuted it
`bench_leopard2_direct_encode` (the purpose-built harness with per-path
introspection), Q=1, high profile, GF8, AVX2, nine cells over K,R <= 16 and
1 KiB-64 KiB. The decisive column is the **ceiling**: force-direct versus
force-transform on the same experiment build, which measures the direct
encoder's intrinsic margin with no dispatch policy in the way.

| Cell | base AUTO | exp AUTO | AUTO effect | ceiling (f-trans/f-direct) |
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

\* base-side noise: the experiment's AUTO time equals its own force-transform
time, so the base's 1.92 us is an outlier, not a direct-path win.

**Force-direct never beats force-transform beyond noise (0.875x-1.017x).** The
+541% opportunity was measured against the 2026-07-16 transform — before GFNI
affine multiplication, radix-eight staging, GF16 fusion at all sizes, and this
campaign's tiling. Those changes made the transform 1.5-2x faster on exactly
these shapes while the direct generator-row path did not move. The margin the
checkpoint recorded no longer exists to harvest.

This is the temperature/proxy lesson in a third form: **stored evidence ages.
A measured advantage over a component is only as durable as that component's
stagnation** — and this campaign spent two days making the comparison baseline
faster.

## What shipped anyway
- The experiment gate (`LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE`, default-off) with
  the refutation recorded at the gate itself, so the +541% checkpoint number
  cannot re-seduce a future reader without the re-screen result in view.
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
Only if a future change makes the direct path cheaper (e.g. GFNI-accelerated
generator-row multiply-add with fewer per-call fixed costs) or the transform
substantially slower on tiny-K: re-establish a force-direct ceiling win first;
AUTO admission is not worth discussing until the ceiling is comfortably above
1.0 across the region.
