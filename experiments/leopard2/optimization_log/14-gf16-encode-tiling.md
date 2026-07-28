# 14 — GF16 encode workspace tiling on the AVX2 path

**Disposition: SHIPPED — up to 2.9x encode, and the largest memory win of the
campaign. The side key and the single-block exclusion below were EXTENDED by
report 15.**

> The `padded_side` table promoted here excluded sides 64 and 128 at every
> shard size, and excluded side-256 single-block counts at every shard size.
> Both exclusions were calibrated at 64 KiB and are wrong above it: side 128 at
> 256 KiB gains **1.90x**, and side-256 single-block — a 0.93x *regression* at
> 64 KiB — is a **1.86x win** at 256 KiB. Report 15 keeps the entries below and
> adds a live-set rule for the region they never covered.

## Idea
Report 13 tiled the GF16 *decode* workspace and recovered 1.16x-1.91x. Encode
has the identical geometry — `work_count = 2 * padded_side` slots held at the
full aligned prefix — so the same argument applies. Reading the guard showed the
existing GF16 encode tile was pinned to a single cell:

```c
codec->original_count == 1000 && codec->recovery_count == 200 &&
codec->padded_side == 256 &&
CodecMayUseAutoAVX512Encode(codec) && ... &&
geometry.aligned_prefix_bytes == kHighGF16EncodeCellBytes   /* exactly 64 KiB */
```

One `(K, R)`, one side, one exact shard size. Everything else — including all of
AVX2 at `padded_side` 512, where the live set is 67 MB against a 32 MB L3 — ran
a single full-payload pass. The code comment said the narrowness was deliberate
*"until neighboring sizes and code shapes have equivalent same-source
evidence"*, i.e. it was an evidence gap, not a finding.

## The first predicate was wrong, and the boundary cells caught it
Decode split on message-block count, so I assumed encode would too and expected
single-block cells to regress. Measured at 64 KiB, tile sizes 4/8/32 KiB:

| Cell | side | blks | 4 KiB | 8 KiB | 32 KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| K=1500 R=300 | 512 | 3 | 1.604 | **1.782** | 1.404 |
| K=2000 R=500 | 512 | 4 | 1.587 | **1.724** | 1.382 |
| K=4096 R=512 | 512 | 8 | 1.492 | **1.647** | 1.351 |
| K=500 R=400 | 512 | **1** | 1.097 | 1.207 | **1.255** |
| K=512 R=512 | 512 | **1** | 1.036 | 1.130 | **1.184** |
| K=600 R=200 | 256 | 3 | 1.065 | 1.066 | **1.164** |
| K=400 R=200 | 256 | 2 | 1.095 | 1.090 | **1.156** |
| K=1000 R=200 | 256 | 4 | 0.970 | 1.040 | **1.132** |
| K=250 R=250 | 256 | **1** | 0.800 | 0.819 | 0.933 |
| K=255 R=129 | 256 | **1** | 0.897 | 0.936 | 0.978 |
| K=200 R=200 | 256 | **1** | 0.875 | 0.909 | 0.963 |
| *K=200 R=50* | *64* | — | *1.002* | *1.004* | *1.003* |
| *K=300 R=100* | *128* | — | *1.013* | *1.012* | *1.008* |

**Single-block at side 512 gains 1.18x-1.26x — it does not regress.** The
regression is confined to *side 256 + single block* (0.80x-0.98x at every tile
size). Had I shipped the block-count predicate I would have left 1.18x-1.26x on
the table at side 512 while believing I was avoiding a regression there.

The real split is **how far the untiled live set overshoots L3**:

- side 512 at 64 KiB = 67.2 MB, about 2.1x L3 — pays at every block count.
- side 256 at 64 KiB = 33.6 MB, about 1.05x L3 — already nearly resident, so it
  only pays when multiple message blocks force repeated re-sweeps of the
  workspace. A single-block count sweeps it roughly twice: there is no reuse to
  recover and the extra passes are pure cost.

## The 64 KiB screen was badly understating the win
Because the *shipped* policy was pinned to exactly 64 KiB, the natural screen
was at 64 KiB. That is the worst size for this change. Sweeping shard size on
the promoted build (explicit AVX2 backend, all parity digests matching):

| Cell | side | 64 KiB | 128 KiB | 256 KiB | 1 MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| K=1000 R=200 | 256 | 1.138 | 1.923 | 2.521 | **2.921** |
| K=600 R=200 | 256 | 1.156 | 1.943 | 2.423 | **2.831** |
| K=4096 R=512 | 512 | 1.662 | 2.243 | 2.532 | **2.748** |
| K=2000 R=500 | 512 | 1.710 | 2.240 | 2.487 | **2.719** |

The gain rises monotonically with shard size in all four cells, because the
untiled live set grows linearly with the shard while the tiled one is pinned to
the requested pass. The ragged neighbour (65600 B) also tiles and gains
1.12x-1.67x.

**Memory is the other half of the result.** Encode scratch at 1 MiB shards:

| Cell | untiled | tiled | reduction |
| --- | ---: | ---: | ---: |
| K=2000 R=500 | 1073.8 MB | 8.5 MB | **126x** |
| K=4096 R=512 | 1073.9 MB | 8.5 MB | **126x** |
| K=1000 R=200 | 536.9 MB | 16.8 MB | **32x** |

A gigabyte of encode scratch for a 2 GB payload was the real cost of the
single-pass layout.

## Shipped policy
```
GF16 + LEGACY_HIGH_V1 + ops->kind == AVX2 + aligned prefix >= 64 KiB:
    padded_side 512                        -> 8 KiB   (any block count)
    padded_side 256 && K > padded_side      -> 32 KiB
    otherwise                               -> untiled
```
Side-256 single-block is excluded because it *regresses*, not because it is
untuned. Sides 32/64/128 were measured neutral with byte-identical scratch and
keep the single-pass layout. The pre-existing pinned AVX-512 cell runs first and
is left untouched; the new guard is skipped when a tile is already selected.

## The contract test that had to change
`leopard2_direct_encode` failed with *"GF16 byte tiling crossed its exact
measured-size boundary"*. It was correct to fail: it asserts (a) tiling is
active only at exactly 64 KiB, and (b) no explicit backend ever tiles. This
change alters both.

It was resolved with measurement, not by relaxing the assertion:

- The exact-size bound was replaced by a `balanced_pass_bytes` helper mirroring
  `ComputeBalancedExecutionTiles`, so neighbouring and larger sizes assert their
  own exact balanced-pass geometry. A **new** bound was added that did not exist
  before — `neighboring_scratch <= one_pass_scratch`, i.e. work rows may never
  exceed the requested tile.
- The fixed-backend loop now expects AVX2 to tile and continues to require that
  scalar, SSSE3, AVX-512 and NEON each retain one complete pass. That second
  half is the load-bearing part: it is what stops the guard widening to
  backends with no evidence.
- The stale comment claiming promotion was limited to the AUTO path was
  corrected. It never was, mechanically — the AUTO context's baseline ops *are*
  AVX2 on this host, which is exactly why the new guard fired there.

The `ragged_scratch` assertion, the AVX-512-dense/AVX2-partial dispatch check,
and the test's own ~300 MiB byte-for-byte differential against the old API are
all unchanged.

## Validation

**Randomized differential, untiled vs tiled build: 240 shapes, 0 mismatches.**
Coverage by `padded_side`: `{1:6, 2:12, 4:17, 16:13, 32:9, 64:11, 128:31,
256:86, 512:16}`, of which **77 shapes activated the new tile**. Compares
original/parity/recovered digests, the round-trip flag and shape acceptance
across GF8/GF16, scalar/SSSE3/AVX2/auto, 1-2 threads, batch 1-2, and sizes
straddling the 64 KiB floor. Tiling only changes how the payload is split into
passes, so any digest difference would be a defect.

**Suites, run serially on a quiet machine:**

- **Stock: 100% passed, 0 failures out of 100.**
- GFNI: 99%, the single failure being `leopard2_portable_isa`, which fails
  identically in the unmodified GFNI build (`-mgfni` applied globally puts
  `vgf2p8affineqb` in `Leopard2BackendAVX2.cpp.o`). This is the standing GFNI
  blocker in `docs/leopard2_gfni_codec.md`, not a regression.
- `leopard2_direct_encode` — the contract discussed above; passes with the
  updated assertions, including its own ~300 MiB byte-for-byte differential
  against the old API on the tiled path.
- Every shard-size row in the sweep above independently confirmed matching
  parity digests against the untiled build.

## Versus leopard1
Explicit AVX2 backend on both sides — leopard1 is an AVX2-era codec, so this is
the apples-to-apples comparison for an AVX2 change. **It is a different baseline
from report 00**, whose GF16 encode ratios use AUTO and therefore include the
AVX-512 encode path plus the pre-existing pinned tile; do not read the two
tables as contradicting each other.

`legacy_comparison` reported `matched` for both builds in every cell.

| Cell | side | shard | vs L1 before | vs L1 after |
| --- | ---: | --- | ---: | ---: |
| K=4096 R=512 | 512 | 256 KiB | 0.999 | **2.541** |
| K=1000 R=200 | 256 | 256 KiB | 0.990 | **2.522** |
| K=2000 R=500 | 512 | 256 KiB | 1.007 | **2.470** |
| K=600 R=200 | 256 | 256 KiB | 1.001 | **2.441** |
| K=2000 R=500 | 512 | 64 KiB | 0.998 | **1.731** |
| K=4096 R=512 | 512 | 64 KiB | 1.003 | **1.675** |
| K=1000 R=200 | 256 | 64 KiB | 1.052 | **1.210** |
| K=600 R=200 | 256 | 64 KiB | 0.998 | **1.173** |
| *K=200 R=50* | *64* | *256 KiB* | *1.036* | *1.027* |
| *K=300 R=100* | *128* | *256 KiB* | *1.011* | *1.011* |

On the AVX2 path these GF16 legacy-high cells were sitting at **parity** with
leopard1 (0.99-1.05). Tiling turns that into **2.44x-2.54x** at 256 KiB shards.

The two italic rows were the untiled controls for *this* report; report 15
subsequently tiles both of them at large shards (side 128 at 256 KiB reaches
1.90x), so they are no longer controls in the shipped policy.
