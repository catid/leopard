# 05 — Batch preflight: hand-rolled heapsort to std::sort

**Disposition: LANDED** (applies to all backends, not GFNI-gated)

## Idea
`SortBatchRangesWithoutAllocation` was a hand-rolled in-place heapsort. Its
documented justification, "without relying on an allocator", does not hold:
`<algorithm>` is already included and `PrepareEncodeBatchRanges` already calls
`std::sort` on the same allegedly allocation-free path. `BatchRangeLess` is a
strict total order, so the sorted sequence is identical and the max-end sweep
cannot observe a different tie order.

The sort matters because the scalable preflight runs as a **serial section
inside the thread pool**: measured 8.5-9.5 us/item, 91.3% of a 64-byte-shard
batch, an Amdahl ceiling of 1.09x at any thread count.

## Result (whole batch, not just the sort)
| Cell | encode | decode |
| --- | ---: | ---: |
| K=3 R=2, B=64 | **1.76x** | **2.58x** |
| K=3 R=2, B=1024 | 1.73x | 1.35x |
| K=100 R=28, B=64 | 1.36x | 1.30x |
| K=100 R=28, B=1024 | **1.72x** | **1.52x** |

`B=1` and `B=8` use the compatibility path and are unchanged.

No leopard1 column: leopard1 has no batch API, so there is nothing to compare
against. This is a Leopard2-versus-Leopard2 improvement.

## Left on the table
A three-pass LSD radix on `range.begin` measured 16.6x against the heapsort
(0.52 us/item vs 8.75), but needs a second equal-size buffer and therefore an
additive change to `leo2_encode_batch_preflight_scratch_size`. Open.
