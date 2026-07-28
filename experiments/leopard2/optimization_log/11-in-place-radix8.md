# 11 — In-place radix-eight

**Disposition: REJECTED — neutral**

## Idea
Report 03's out-of-place radix-eight butterfly landed 1.33x at `K=224 R=32`. Its
design keeps only eight data vectors live: all eight inputs are consumed into
registers before any output is written. The natural follow-up was an **in-place**
variant, eliminating the separate output buffer and its store traffic.

## Result
**Neutral.** No promotion.

## Why
Out-of-place was never paying a penalty that in-place could remove.

1. The kernel is **memory-bound**, and the byte count it moves is the same
   either way: eight vectors in, eight vectors out. In-place changes the
   destination address, not the traffic.
2. Out-of-place writes are streaming stores to a destination that is about to be
   read by the next layer, so they stay in cache. There is no write-allocate
   penalty to recover.
3. In-place costs register pressure. Holding eight inputs while producing eight
   outputs into the same registers forces either a longer dependency chain or a
   spill, and on this kernel the matrices are already consumed as **memory
   operands** — which is free precisely because the kernel is memory-bound and
   has spare issue slots. In-place erodes that headroom.

This is the same temperature argument as report 08 from the other side: the
out-of-place kernel's extra buffer is warm, so its cost was already near zero,
and removing a near-zero cost yields a near-zero win.

**Rule extracted:** "eliminate the temporary buffer" is only a win if the
temporary is cold or the traffic actually shrinks. Count bytes moved, not
buffers allocated.
