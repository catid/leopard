# 39 — Sparse-high forced-path discovery

## Result

The bounded legacy-high GF8 sparse-Q1 forced-path discovery passed on the
explicit AVX2 build. All 54 candidate cells cleared the preregistered
five-percent decision threshold using the lower endpoint of their marginal
paired-log 95-percent Student-t intervals. The weakest candidate point gain was
25.23 percent, the median gain was 123.50 percent, and the weakest lower
endpoint was 17.64 percent.

This result makes **no production change**. The campaign forced direct and
transform execution through test hooks, did not invoke production AUTO, and
did not measure production table-preparation cost. Its intervals are per-cell
marginal intervals at `df=2`, not a simultaneous campaign guarantee. The
current default-off predicate also spans more K/R/byte/API cases than this
finite grid. Production-AUTO qualification remains separately tracked in
Bead `leopard-79h.42.9`.

Leopard main has no equivalent partial-output Q=1 API, so this campaign has no
Leopard1 ratio. It compares the two Leopard2 execution paths on identical
inputs and requested parity rows.

## Accepted campaign

The accepted source was commit
`c462f0ddef75979749c8cdd8e7199326ff1be892`, tree
`7626c9036c698c8d6278d47d1dbe11c7386a212f`, with source fingerprint
`7d7579e27905758e7301b9849b2969597072293640abd85f8f8093a8d7fcc9c7`.
The controlled executable SHA-256 was
`3705ca83195da93457f4afbfbebc260a3323c68fadd3c47e7fc15c2754ee7b86`
and its exact owner-only mode was `0700`.

The host was an AMD Ryzen Threadripper PRO 9985WX. A read-only 30-second
`/proc/stat` preflight found 20 SMT pairs with zero non-idle jiffies in all 120
windows. Timing used CPU 14 with sibling 78 reserved, three ABBA rounds, 15
iterations after four warmups, reuse 64, one worker and one codec thread.
Batches 1, 4, and 16 were represented by the frozen 91-cell grid.

Representative weakest candidate results are:

| Cell | Direct us | Transform us | Gain | Marginal 95% interval |
| --- | ---: | ---: | ---: | ---: |
| K=4 R=2 Q=1, 4096 B, batch 1 | 0.249360 | 0.312266 | 25.23% | 19.67%–29.10% |
| K=16 R=2 Q=1, 1088 B, batch 1 | 0.271031 | 0.348180 | 28.46% | 17.64%–36.57% |
| K=16 R=2 Q=1, 4096 B, batch 16 | 19.268532 | 27.118289 | 40.74% | 32.03%–52.29% |

Every accepted job passed selected/reference parity, direct/transform parity,
non-vacuous parity-identity, unrequested-output, process-containment, frozen
binary, source, and isolation checks. Every accepted sibling delta was zero.
Normal and `python -O` deterministic reanalysis reproduced `analysis.json`
byte for byte.

## Excluded neighbors

All 37 excluded neighbors passed correctness and evidence validation. Twelve
were slower under forced direct execution; all material losses lie on already
excluded selector boundaries:

- batch-1 `R=1`: 10 of 12 cells regressed, with a minimum point gain of
  -29.04 percent. The two batch 4/16 controls improved slightly but their
  intervals included zero. The predicate requires `R>1`.
- prefix `Q=2`: two of five cells regressed, reaching -22.80 percent. The
  predicate requires exactly one requested output.
- all three holey `Q=2` controls improved under forced direct, but they remain
  outside the singleton-output predicate and do not erase the prefix losses.
- the remaining profile, field, K=1, and ragged/threshold controls stay excluded
  even when forced direct happened to win; discovery results do not authorize
  broadening a guard.

The weakest excluded result was K=16/R=1 at 1024 bytes: -29.04 percent with a
marginal interval of -30.98 to -26.88 percent. The K=16/R=2 prefix-Q2 control
was -22.80 percent with an interval of -29.90 to -17.82 percent. These negative
cells are retained in the compact checkpoint rather than omitted from the
headline.

## Isolation retries

No contaminated measurement entered inference. An initial CPU63/127 diagnostic
rejected 15 jobs. The first quiet-pair pass rejected 12 jobs for one or two
sibling jiffies; a strict resume revalidated and retained 79 clean jobs and
reran only those 12. The second pass rejected one job for one sibling jiffy;
the final resume reran that job cleanly. Thus the accepted corpus contains 79,
11, and one jobs from the three CPU14 passes respectively.

The two pre-resume CPU14 bundles were archived before replacement. Their
SHA-256 values are
`59e4b9c50e6edb1ac2558f51a9098430d6bff5534038165c51cebf183fee7252`
and
`2871fe376d7694d424d23bbea9dc6ed85a441f829e13c52b9c05b3c86c321786`.
The compact checkpoint lists every discarded job and marks all three rejected
states as neither used nor pooled.

## Evidence and next boundary

The accepted top-level hashes are:

- manifest: `61a2d0e8fcd7d6d225e260a935a5a37a95cae5976960165bd29bdfac821603e8`;
- matrix: `639fdc0fb717c7d180010a713192920f1f6b45f0f785c37482a13a741f8f511a`;
- analysis: `3f1c60ebc9f9bc780daccef8d22dcc558e7d150db749f22a26fff2cd0bb513fc`;
- controlled build: `21663f5e202a90eb73b739506731091cbb63e79d20d15867761fc345e98595c8`.

The 160 MiB raw bundle remains generated evidence on `ripper.lan`. The
committed projection, including every cell and discarded isolation record, is
`experiments/leopard2/direct_encode/results/sparse_high_avx2_checkpoint_20260831.json`.

The next gate must freeze a finite production predicate, resolve whether AVX2
means caller-requested or effective backend, add a tables-on/AUTO-off control,
and measure actual production AUTO plus setup/amortization across one-shot,
batch, and binding APIs. It remains default-off until that independent gate.
