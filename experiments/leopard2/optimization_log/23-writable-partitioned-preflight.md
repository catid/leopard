# 23 — Writable-partitioned batch alias preflight (findings 2+3, second step)

**Disposition: LANDED — 1.01x-1.06x end-to-end batch, 53 lines removed.
The findings' Amdahl model did not reproduce; recorded at measured value.**

## Change
`SortBatchRangesWithoutAllocation` + `ValidateSortedBatchRanges` (introsort over
every record, then a two-accumulator sweep) are replaced by
`ValidateBatchRangesPartitioned`: `std::partition` writable records to the
front, sort only those, sweep them for writable/writable overlap, then give
each readable record one predecessor probe via `std::upper_bound` — valid
because the writable set must be pairwise disjoint anyway, so sorted it is a
disjoint ascending interval set and only the predecessor can intersect.

Readable/readable overlap is permitted by the contract, so ordering the
readable majority (K per item versus R+1 writables) bought nothing.
Complexity: O(N log N) -> O(W log W + (N-W) log W) with W = B*(R+1).
Allocation-free as before.  Zero-length-range semantics preserved exactly
(a writable beginning at a readable's end intersects under neither
formulation).  `leopard2_batch_aliasing` and the public-API contract pass
unchanged — the behavioral gate, since no source-text contract pins the
implementation.

## Measured (K=100 R=28, B=1024, AVX2, digests matching)

| bytes / threads | base | new | speedup |
| --- | ---: | ---: | ---: |
| 8192 / 1 | 230.7 ms | 218.4 ms | **1.056x** |
| 64 / 8 | 136.2 ms | 130.2 ms | **1.045x** |
| 64 / 1 | 140.0 ms | 135.0 ms | **1.037x** |
| 1024 / 8 | 143.7 ms | 140.5 ms | 1.023x |
| 1024 / 1 | 151.2 ms | 150.4 ms | 1.006x |

## The model-versus-measurement gap, stated plainly
Findings 2+3 projected much more: "heapsort 98.6% of the preflight",
"preflight 91.3% of the batch at 64-byte shards, Amdahl ceiling 1.09x", and a
measured 3.1x sort-level win for this exact scheme over `std::sort`.  The
end-to-end batch entry point moves only 1.01x-1.06x.  Two reasons, both
instructive:

1. **Report 05 already banked the big step.**  The findings' decomposition was
   taken against the heapsort era; the heapsort -> introsort landing halved
   the sort before this change arrived.
2. **The isolated-component measurement overstated the share.**  The 8.75
   us/item sort figure came from timing the preflight components alone; the
   production batch path around it (encode work, per-item bookkeeping, the
   thread-pool run) is evidently much larger than the model's 9.6 us/item
   total, so even a 3x sort-level win moves the end-to-end number by percent,
   not halves.  **A component-level measurement is not an end-to-end claim
   until the component's share is re-measured in place.**

Kept because the change is strictly-less-work with a measured gain everywhere
and a negative nowhere, plus a real simplification.  Findings 2 and 3 should
be treated as partially resolved with their headline shares corrected; the
16.6x LSD-radix option is NOT worth its ABI change at the true share and is
hereby declined.
