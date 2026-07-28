# 16 — Low-rate profile tiling

**Disposition: ENCODE SHIPPED (1.64x-2.55x). Decode measured separately.**

## Idea
After report 15 unified the GF16 tile policy on the live set, one region was
still excluded by construction: **both tiling guards required
`LEO2_PROFILE_LEGACY_HIGH_V1`**, so the low-rate profile never tiled at any
shard size.

That exclusion was not a measurement. `EncodeLayout` sets
`work_count = padded_side * 2` for *every* profile, so `LOW_V1` holds its rows
on the identical formula and reaches the identical live set. It had simply been
outside the scope of the campaign that produced the guards.

Low-rate matters disproportionately here: **leopard1 rejects `R > K` outright**,
so this profile has no leopard1 comparison at all and is measured against stock
Leopard2 only. It is also the regime the `LOW_V1` profile exists to serve.

## Result — encode
Against the shipped build (which tiles legacy-high but not low), explicit AVX2,
all parity digests matching. Only the live-set branch applies: the two
`padded_side` entries from report 15 stay gated on `LEGACY_HIGH_V1` because they
were calibrated on that traversal.

| Cell | side | live set | speedup | encode scratch |
| --- | ---: | ---: | ---: | --- |
| K=256 R=768, 1 MiB | 256 | 537 MB | **2.553x** | 537 MB -> 17 MB |
| K=128 R=896, 1 MiB | 128 | 268 MB | **2.365x** | 268 MB -> 17 MB |
| K=1000 R=2000, 64 KiB | 1024 | 134 MB | **2.306x** | 134 MB -> 17 MB |
| K=200 R=800, 256 KiB | 256 | 134 MB | **2.304x** | 134 MB -> 17 MB |
| K=128 R=896, 256 KiB | 128 | 67 MB | **1.643x** | 67 MB -> 17 MB |
| *K=32 R=224, 1 MiB* | *32* | *67 MB* | *1.001x* | *67 MB -> 67 MB* |

Encode scratch collapses to a flat 17 MB regardless of shard size, because the
tile is chosen to hold the retained rows at the target rather than to divide the
payload.

**The rule generalised past every table that preceded it.** `K=1000 R=2000` has
`padded_side` 1024 — a side no tile table ever mentioned — and it tiled at 8 KiB
for 2.306x with no new entry.

## The control row is the useful one
`K=32 R=224` did **not** tile, despite a 67 MB live set that clears the
threshold. Its scratch is byte-identical (67 MB -> 67 MB), which is the
activation witness from report 13 doing its job.

The reason: `K + R = 256`, so that shape resolves to **GF8**, and the live-set
branch is GF16-only. So the guard behaved correctly — and the row exposes the
next gap:

> **GF8 low-rate never tiles either.** The GF8 tile table is also gated on
> `LEGACY_HIGH_V1`. GF8's own thresholds target a ~4 MB retained set rather than
> GF16's 16 MB (see report 15's note on GF8 already being live-set calibrated),
> so the GF16 constants must not simply be extended to it. This needs its own
> sweep and is left open.

## Result — decode
Decode was measured separately, with the base build already carrying low-rate
*encode* tiling so the comparison isolates the decode change. Encode read
0.984x-1.011x across that run, which is the control confirming the isolation.

`kDecodeRuleWorkspaceTiled` is reachable for `LOW_V1` — `Leopard2Dispatch.h`
excludes only *other* profiles — so the path condition was already satisfied and
the profile test was the sole blocker.

## Result — promoted policy, both paths, versus the untiled base

| Cell | side | live set | encode | decode | scratch |
| --- | ---: | ---: | ---: | ---: | --- |
| K=256 R=768, 1 MiB | 256 | 537 MB | 2.515x | **3.099x** | 537 MB -> 17 MB |
| K=128 R=896, 1 MiB | 128 | 268 MB | 2.365x | **2.954x** | 268 MB -> 17 MB |
| K=200 R=800, 256 KiB | 256 | 134 MB | 2.293x | **2.589x** | 134 MB -> 17 MB |
| K=1000 R=2000, 64 KiB | 1024 | 134 MB | 2.328x | **2.159x** | 134 MB -> 17 MB |
| K=128 R=896, 256 KiB | 128 | 67 MB | 1.637x | **2.080x** | 67 MB -> 17 MB |
| *K=32 R=224, 1 MiB (GF8)* | *32* | *67 MB* | *1.001x* | *0.996x* | *unchanged* |

## A silent widening that the control row caught
Relaxing the profile test at the top of `AVX2DecodeExecutionTileBytes` does not
only admit low-rate to the GF16 live-set branch — it also admits it to the
**GF8 decode table further down**, which is calibrated on the legacy-high
traversal. That happened, and it was not visible in the headline numbers: the
GF8 low-rate cell jumped to **2.540x** in the experiment run.

It surfaced only because the screen carried a GF8 cell and its decode moved
despite an edit that looked GF16-only. The fix was an explicit
`profile == LEO2_PROFILE_LEGACY_HIGH_V1` test on the GF8 branch, and the control
row above (0.996x, scratch unchanged) is the evidence it held.

**GF8 low-rate is therefore left open, not taken**, despite a measured 2.540x.
One favourable cell is not evidence that the other GF8 sides and shard sizes
behave, and GF8's floors target a ~4 MB retained set against GF16's 16 MB (see
report 15). Extending it needs its own sweep. This is the same judgement that
kept side-256 single-block excluded in report 14 — which, notably, later turned
out to be a 1.86x *win* at a different shard size, so the follow-up here may
well pay.

## Validation
- **Stock suite: 100% passed, 0 failures out of 100.**
- GFNI suite: 99%, sole failure `leopard2_portable_isa`, which fails identically
  in the unmodified GFNI build (`-mgfni` applied globally puts
  `vgf2p8affineqb` in `Leopard2BackendAVX2.cpp.o`). Pre-existing and documented
  in `docs/leopard2_gfni_codec.md`.
- **Randomized differential: 300 shapes, 0 mismatches**, sides
  `{1:17, 2:16, 4:16, 16:19, 32:14, 64:12, 128:24, 256:81, 512:33}`, of which
  **69 activated a tile**. The generator draws `R` independently of `K`, so
  low-rate shapes are covered.
- Every timing row above independently confirmed matching parity and recovered
  digests against the untiled build.
- The GF8 low-rate control (`K=32 R=224`, 1 MiB) read 1.001x encode / 0.996x
  decode with byte-identical scratch at the point GF8 was still excluded,
  confirming the `LEGACY_HIGH_V1` gate added to the GF8 branch held.

Re-validated after GF8 low-rate decode was subsequently promoted:

- **Stock suite: 100% passed, 0 failures out of 100.**
- GFNI suite: 99%, sole failure the pre-existing `leopard2_portable_isa`.
- **Randomized differential: 320 shapes, 0 mismatches**, sides
  `{1:15, 2:13, 4:20, 16:19, 32:24, 64:17, 128:29, 256:99, 512:29}`, of which
  **83 activated a tile**.

## GF8 low-rate — the follow-up, swept and shipped
The accidental 2.540x above was not promoted on its own. The sweep it called
for was then run: **GF8's own table** applied to the low profile — not GF16's
constants, whose 16 MB target does not apply to a field calibrated at ~4 MB.

Against the shipped build (GF8 low still excluded), all recovered digests
matching:

| Cell | side | live set | decode | scratch |
| --- | ---: | ---: | ---: | ---: |
| K=64 R=192, 1 MiB | 64 | 134 MB | **3.145x** | 4 MB |
| K=32 R=224, 1 MiB | 32 | 67 MB | **2.577x** | 4 MB |
| K=16 R=240, 1 MiB | 16 | 34 MB | **1.696x** | 4 MB |
| K=32 R=224, 512 KiB | 32 | 34 MB | **1.597x** | 4 MB |
| K=64 R=192, 256 KiB | 64 | 34 MB | **1.488x** | 4 MB |
| *K=32 R=224, 256 KiB* | *32* | *17 MB* | *0.998x* | *below GF8's own floor* |
| *K=128 R=128, 1 MiB* | *128* | *268 MB* | *1.001x* | *`R == K` -> legacy-high* |

Five gains across four `padded_side` values, **no regression**. Both neutral
rows are neutral for identifiable reasons rather than by dismissal: side 32
requires a 384 KiB prefix in GF8's existing table, so 256 KiB is that
calibration correctly declining the shape; and `K=128 R=128` has `R == K`, so
it resolves to legacy-high and serves as an untouched control proving the change
did not widen into the high profile.

**3.145x is the largest decode gain of the campaign.**

## GF8 low-rate encode — the last gap, also swept and shipped
The GF8 encode table is a separate guard in `EncodeLayout`, gated on
`LEGACY_HIGH_V1` the same way. Its floors are conservative for low-rate:
`original_count > padded_side` is always false when `P = ceil_pow2(K) >= K`, so
low-rate always draws the larger single-block minimums (4 MiB at side 16, 2 MiB
at side 32, 1 MiB at side 64).

| Cell | side | shard | speedup | scratch |
| --- | ---: | --- | ---: | ---: |
| K=65 R=128 | 128 | 1 MiB | **2.864x** | 17 MB |
| K=100 R=128 | 128 | 1 MiB | **2.752x** | 17 MB |
| K=100 R=128 | 128 | 512 KiB | **2.699x** | 17 MB |
| K=64 R=192 | 64 | 2 MiB | **2.477x** | 4 MB |
| K=64 R=192 | 64 | 1 MiB | **2.461x** | 4 MB |
| K=32 R=224 | 32 | 2 MiB | **2.105x** | 4 MB |
| K=16 R=240 | 16 | 4 MiB | **1.838x** | 4 MB |
| *K=64 R=192* | *64* | *512 KiB* | *0.999x* | *67 MB, below floor* |
| *K=128 R=128* | *128* | *1 MiB* | *0.988x* | *`R == K` -> legacy-high* |

Gains at every `padded_side`, **no regression**.

### A sweep gap worth chasing
Two cells first came back as `SKIP` (`K=100 R=156` and `K=120 R=136`, both side
128) — which would have left an entire `padded_side` untested while the table
looked complete. Running them by hand showed **both builds reject them
identically** with `Invalid or unsupported shard counts (-3)`, so it was a
pre-existing API constraint rather than anything this change did: `LOW_V1`
requires `R` to be a multiple of `P`. The working cells confirm it —
`K=64 R=192` is `3*P`, `K=32 R=224` is `7*P`, `K=16 R=240` is `15*P`. Choosing
`R = 128 = 1*P` covered side 128 properly.

**A `SKIP` in a sweep is missing coverage, not a neutral result.** Reproduce it
in both builds before accepting it.

## Disposition summary
- **GF16 low-rate encode and decode: shipped.**
- **GF8 low-rate encode and decode: shipped**, by extending GF8's own calibrated
  table to the low profile — not by importing GF16's constants, whose 16 MB
  retained-set target does not apply to a field calibrated at ~4 MB.
- GF8 *legacy-high* is unchanged: its table was already live-set calibrated (see
  report 15), which is independent corroboration of the model rather than a
  missed opportunity.
