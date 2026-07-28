# 17 — Does the live-set target need to scale with thread count?

**Disposition: HYPOTHESIS REFUTED — no change shipped, and that is the result**

## The concern
Reports 13-16 promote a tile sized so the retained transform rows land near a
16 MB live set. Every one of those measurements was **single-threaded on one
pinned core**, so the target implicitly assumes one thread owning the whole L3.

Concurrent stripes share L3. At `T` threads the aggregate live set is
`T * target`, so 8 threads at 16 MB nominally wants 128 MB against this die's
32 MB — which predicts the tiled build should itself start thrashing, and that a
thread-scaled target (`16 MB / T`) should recover the loss.

That is a plausible, mechanically-reasoned argument for adding a knob. It is
also wrong.

## Host detail that matters
This is a 9950X3D with **asymmetric L3 across its two dies**:

| CCD | CPUs | L3 |
| --- | --- | ---: |
| CCD0 (3D V-cache) | 0-7, 16-23 | 96 MB |
| CCD1 | 8-15, 24-31 | 32 MB |

All single-core work in this campaign pinned CPU 14, i.e. **CCD1, the 32 MB
die**. The threaded runs below pin `8-15` so all eight threads share that one
32 MB L3 — the worst case for the contention hypothesis, deliberately.

## Step 1: does the benefit survive concurrency at all?
Tiled versus untiled, CCD1, batch scaled with threads:

| Cell | 1 thr | 2 thr | 4 thr | 8 thr |
| --- | ---: | ---: | ---: | ---: |
| K=1000 R=200, 256 KiB — encode | 2.551x | **2.932x** | 1.979x | 1.379x |
| K=1000 R=200, 256 KiB — decode | 2.575x | 4.756x | **5.874x** | 3.186x |
| K=2000 R=500, 64 KiB — encode | 1.656x | 3.336x | **3.602x** | 2.527x |
| K=2000 R=500, 64 KiB — decode | 1.974x | **2.764x** | 1.934x | 1.422x |

Two things stand out. Tiling **never becomes a loss** — 1.38x to 5.87x across
the board. And the benefit is *larger* at 2-4 threads than at 1, peaking at
**5.874x**, because the untiled layout degrades faster under contention than the
tiled one does. But it falls off at 8 threads, which is exactly what the
contention hypothesis predicts.

## Step 2: the test that refutes it
If the 8-thread fall-off were L3 contention, shrinking the target so that
`T * target` fits L3 should recover it. At 8 threads a 2 MB target gives a
16 MB aggregate — comfortably resident. Ratios versus the untiled base:

| Cell | thr | 16 MB (shipped) | 4 MB | 2 MB |
| --- | ---: | --- | --- | --- |
| K=1000 R=200, 256 KiB | 8 | e1.39x d2.74x | e1.38x d2.75x | e1.38x d2.75x |
| K=2000 R=500, 64 KiB | 8 | e2.52x d1.41x | e2.52x d1.42x | e2.53x d1.42x |
| K=2000 R=500, 64 KiB | 4 | e3.57x d1.96x | e3.56x d1.93x | e3.55x d1.94x |

**Shrinking the target eightfold changes nothing.** The fall-off is DRAM
bandwidth saturation, not L3 contention: at 8 cores the memory path is already
the binding constraint and no tile size moves it.

## Why this is worth a report
The thread-scaled target was a well-motivated feature with a clean mechanical
story, a knob the context already had the information to drive, and a plausible
several-percent case. Building it would have been reasonable. It would also
have been **pure cost** — added policy surface, another calibration axis, and a
constant that later readers would assume was measured.

One three-point screen killed it. The shipped constant is correct as-is,
including under concurrency, and now says so in the source with the numbers
attached so the question does not get re-opened from first principles.

**Rule: a mechanism that predicts an effect is not evidence the effect exists.
Screen the knob before adding it — a refuted hypothesis is a feature you do not
have to maintain.**

## Caveat recorded, not fixed
The constants remain calibrated for a **32 MB L3**, and this CPU has a 96 MB die
as well. On CCD0 the 64 MB threshold would tile a 70 MB live set that would have
stayed resident — a small pessimisation, not a correctness issue. The
live-set formulation makes the fix obvious (query L3 at context init and scale
both constants), which would also handle the CCD asymmetry on this machine and
cross-host portability at once. Left open; the existing GF8 tables have the same
host-pinned property and predate this work.
