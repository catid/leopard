# 10 — Non-temporal stores in the staging kernel

**Disposition: REJECTED — 2-3x worse**

## Idea
The out-of-place staging kernel writes work slots it does not read back within
the same pass, so `_mm256_stream_si256` should avoid polluting cache.

## Result
Combined read+write throughput, isolated screen:
| bytes/shard | ordinary stores | non-temporal |
| ---: | ---: | ---: |
| 1024 | 153.4 GB/s | 39.0 GB/s |
| 4096 | 157.6 | 54.7 |
| 16384 | 156.7 | 60.6 |
| 65536 | 130.4 | 61.2 |

## Why it failed
The premise is wrong. The work slots **are** re-read — by the following
transform layers. Bypassing the cache forces those layers to re-fetch from
memory. Non-temporal stores only pay for data that is genuinely write-once, and
nothing in this transform is.
