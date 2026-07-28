# 12 — GF8 encode byte tiling at 64 KiB shards

**Disposition: RECORDED, NOT SHIPPED — real but a sharp optimum**

## Idea
GF8 encode tiling does not engage at 64 KiB: for `padded_side == 32` it requires
>= 512 KiB shards. The work array is then 32 x 64 KiB = 2.1 MB, exceeding the
1 MB L2, so every in-place round goes to L3. That threshold table was generated
offline with **nibble** kernels; GFNI and radix-eight have since made the
arithmetic much cheaper, so the calibration was a good candidate for being
stale — the same reasoning that correctly predicted the GF16 fusion cap was
stale.

## Result — it is only slightly stale
Against the current GFNI + radix-eight build, all digests matching:

| Cell | untiled | 8 KiB tile | 16 KiB tile | 32 KiB tile |
| --- | ---: | ---: | ---: | ---: |
| K=224 R=32, 64 KiB | 446.4 us | **1.043x** | 0.900x | 0.970x |
| K=224 R=32, 16 KiB | 114.3 us | **1.040x** | 0.989x | 0.990x |
| K=200 R=32, 64 KiB | 426.3 us | **1.053x** | 0.919x | 0.838x |
| K=160 R=32, 64 KiB | 339.1 us | **1.044x** | 0.910x | 0.840x |

## Why it is not shipped
8 KiB is consistently worth 4-5%, but 16 KiB and 32 KiB lose 10-16%. That is a
sharp optimum measured on four cells against a policy the project derived from a
pinned offline sweep. Shipping it would be exactly the overfit the existing
side-specific table exists to prevent.

**To revisit:** run the same offline sweep across `padded_side` x shard size x
block count that produced the current table. The 4-5% is real and repeatable;
only the evidence base is too narrow.
