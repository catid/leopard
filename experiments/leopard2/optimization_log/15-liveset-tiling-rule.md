# 15 — Replacing the `padded_side` tile tables with a live-set rule

**Disposition: SHIPPED — closes a blind spot in reports 13 and 14, worth up to 3.2x**

## The bug in my own two findings
Reports 13 and 14 both promoted a tile keyed on `padded_side`:

```
decode:  side 256 or 512, prefix >= 64 KiB  -> 16 KiB
encode:  side 512 -> 8 KiB;  side 256 && K > side -> 32 KiB
```

Both tables were calibrated at **64 KiB shards only**. Sides 64 and 128
measured neutral there and were excluded. That exclusion was correct at 64 KiB
and wrong everywhere else, because the quantity that actually governs the win is
not the side — it is the **live set**:

```
live_set = 2 * padded_side * aligned_prefix_bytes
```

Against a 32 MB L3:

| side | 64 KiB | 256 KiB | 1 MiB |
| ---: | ---: | ---: | ---: |
| 64 | 8.4 MB | 33.6 MB | 134 MB |
| 128 | 16.8 MB | 67 MB | 268 MB |
| 256 | 33.6 MB | 134 MB | 537 MB |
| 512 | 67 MB | 268 MB | 1074 MB |

Sides 64 and 128 genuinely fit L3 at 64 KiB. At 1 MiB they overshoot it by 4x
and 8x — and a side-keyed table excluded them at *every* size.

## The rule, and the prediction that tested it
Every measured optimum in both reports lands near the same live set, not the
same tile:

| Path | side | best tile | 2 * side * tile |
| --- | ---: | ---: | ---: |
| decode, multi-block | 256 | 32 KiB | 16.8 MB |
| decode | 512 | 16 KiB | 16.8 MB |
| encode | 512 | 8 KiB | 8.4 MB |
| encode, multi-block | 256 | 32 KiB | 16.8 MB |

So: **tile so the retained rows land near 16 MB, once the untiled set exceeds
64 MB (2x L3).** This was promoted only after a falsifiable prediction:

> side 128 at 256 KiB is a 67.1 MB live set, essentially identical to side 512
> at 64 KiB (67.2 MB), which measured 1.73x. If the rule is right, tiling side
> 128 there should gain about the same.

Measured: **1.895x** (encode), **2.281x** (decode). The rule predicted an
untested region before it was measured.

## Result
Against the untiled build, explicit AVX2 backend, all digests matching.
Rows in *italics* are below the 64 MB threshold and correctly stay untiled.

### Encode
| Cell | side | live set | before | after |
| --- | ---: | ---: | ---: | ---: |
| K=2000 R=500, 1 MiB | 512 | 1074 MB | 1.00 | **2.735** |
| K=1000 R=200, 256 KiB | 256 | 134 MB | 1.00 | **2.504** |
| K=300 R=100, 1 MiB | 128 | 268 MB | 0.999 | **2.497** |
| K=200 R=50, 1 MiB | 64 | 134 MB | 1.002 | **2.491** |
| K=300 R=100, 256 KiB | 128 | 67 MB | 0.997 | **1.903** |
| K=250 R=250, 256 KiB | 256 | 134 MB | 1.002 | **1.864** |
| K=2000 R=500, 64 KiB | 512 | 67 MB | 1.00 | **1.718** |
| K=1000 R=200, 64 KiB | 256 | 34 MB | 1.00 | **1.154** |
| *K=250 R=250, 64 KiB* | *256* | *34 MB* | *1.015* | *1.003* |
| *K=200 R=50, 256 KiB* | *64* | *34 MB* | *1.000* | *1.005* |

### Decode
| Cell | side | live set | before | after |
| --- | ---: | ---: | ---: | ---: |
| K=300 R=100, 1 MiB | 128 | 268 MB | ~1.0 | **3.226** |
| K=200 R=50, 1 MiB | 64 | 134 MB | ~1.0 | **3.044** |
| K=2000 R=500, 256 KiB | 512 | 268 MB | ~1.0 | **2.683** |
| K=1000 R=200, 256 KiB | 256 | 134 MB | ~1.0 | **2.553** |
| K=250 R=250, 256 KiB | 256 | 134 MB | ~1.0 | **2.444** |
| K=300 R=100, 256 KiB | 128 | 67 MB | ~1.0 | **2.281** |
| K=2000 R=500, 64 KiB | 512 | 67 MB | 1.0 | **1.945** |
| K=250 R=250, 64 KiB | 256 | 34 MB | 1.0 | **1.367** |
| K=1000 R=200, 64 KiB | 256 | 34 MB | 1.0 | **1.124** |
| *K=200 R=50, 256 KiB* | *64* | *34 MB* | *1.0* | *1.006* |

**Best cell of the entire campaign: GF16 `K=300 R=100` decode at 1 MiB,
3.226x** — a `padded_side` of 128, which both shipped tables excluded outright.

## Post-ship confirmation at small sides
The shipped guard has no side restriction, so it tiles combinations the
promotion sweep never covered — small GF16 sides at large shards. That is
exactly the "guard applied outside its measured range" failure mode this
campaign hit three times, so it was checked rather than assumed. Sides 8-64,
forced GF16, all digests matching:

| Cell | side | live set | encode | decode |
| --- | ---: | ---: | ---: | ---: |
| K=100 R=32, 4 MiB | 32 | 268 MB | 2.395x | **3.376x** |
| K=100 R=16, 4 MiB | 16 | 134 MB | 2.503x | **3.189x** |
| K=200 R=64, 1 MiB | 64 | 134 MB | 2.458x | **3.140x** |
| K=100 R=32, 1 MiB | 32 | 67 MB | 1.985x | 2.515x |
| K=100 R=16, 2 MiB | 16 | 67 MB | 2.023x | 2.504x |
| K=60 R=8, 4 MiB | 8 | 67 MB | 1.934x | 2.401x |

No regression anywhere in the extrapolated region; the rule is confirmed across
`padded_side` 8 through 1024. The 3.376x decode at side 32 supersedes
`K=300 R=100` (3.226x) as the best cell of the campaign.

## Why it is a union, not a replacement
The rule does not subsume the calibrated entries. Side 256 with multiple message
blocks at 64 KiB is only a 33.6 MB live set — below the threshold — yet tiling
still pays 1.12x-1.16x, because each message block re-sweeps the workspace and
the reuse is worth capturing even when the set nearly fits. So:

```
side 512                     -> 8 KiB (encode) / 16 KiB (decode)   [calibrated]
side 256 && K > side         -> 32 KiB (encode) / 16 KiB (decode)  [calibrated]
otherwise, live set >= 64 MB -> 16 MB / (2 * side)                 [derived]
```

The calibrated entries win where they were measured; the live set governs the
region they never covered. Where both apply the rule reproduces the table within
noise (2.502 vs 2.511, 1.687 vs 1.712), which is the evidence that they are the
same phenomenon.

## The threshold is load-bearing
64 MB is not a tuning knob — it is what keeps the one measured regression out.
Side-256 single-block at 64 KiB is a 33.6 MB set swept about twice; tiling it
costs **0.80x-0.98x**. It sits below the threshold and stays untiled, confirmed
at 1.003x.

Note that this same shape is a **1.86x win at 256 KiB**. My original single-block
exclusion in report 14 was therefore *also* size-specific — the identical error
as the exact-size bound I had just criticised, committed one policy dimension
over. The live-set formulation removes the dimension instead of patching it.

## Validation
- **Stock suite: 100% passed, 0 failures out of 100**, with both the encode
  union policy and the decode live-set rule in place.
- GFNI suite: 99%, the sole failure being `leopard2_portable_isa`, which fails
  identically in the unmodified GFNI build (`-mgfni` applied globally puts
  `vgf2p8affineqb` in `Leopard2BackendAVX2.cpp.o`). Pre-existing, documented in
  `docs/leopard2_gfni_codec.md`.
- **Randomized differential: 300 shapes, 0 mismatches**, sides
  `{1:17, 2:19, 4:14, 16:21, 32:16, 64:13, 128:27, 256:86, 512:28}`, of which
  **71 activated a tile**. Plus a separate 200-shape run against the stock
  (non-GFNI) build.
- Every timing row above independently confirmed matching parity/recovered
  digests against the untiled build.

A note on how the intermediate state was caught: while the constants were still
local to `EncodeLayout`, `leopard2_field_options_matrix` and
`leopard2_test_hook_isolation` failed to compile. Those two gates rebuild
`leopard2.cpp` from the live source tree rather than from a build directory,
so they are the only ones that detect a mid-edit inconsistency the main build
has already cached past. Worth remembering when a change spans two functions.

## Method note worth keeping
Three times in this campaign a policy was keyed on a proxy that held only at the
sizes swept: the exact-64-KiB encode bound, the side tables, and the
single-block exclusion. Each time the fix was to find the physical quantity the
proxy stood in for and key on that instead.

**When a guard is keyed on a shape parameter, ask what physical quantity it is a
proxy for, and whether that proxy survives the range you did not sweep.**
