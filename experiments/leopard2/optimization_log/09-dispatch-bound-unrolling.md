# 09 — Dispatch-bound unrolling

**Disposition: REJECTED — 1.01x**

## Idea
The transform executors call through a backend `Ops` function-pointer table once
per butterfly group. At small `padded_side` the group is short, so the indirect
call and its argument marshalling were hypothesised to be a meaningful fraction
of the work. Unrolling the dispatch — issuing several groups per indirect call —
would amortise it.

## Result
**1.01x.** Inside run-to-run noise; not promoted.

## Why
The indirect call is correctly predicted by the branch predictor. The table is
written once at context init and never mutated afterwards, so every call site
sees the same target for the lifetime of the process. An indirect branch with a
single observed target costs essentially the same as a direct one after warmup,
and the benchmark's warmup iterations guarantee that state.

The measurement that should have preceded the idea: GF8 encode was already
established at **~99.9% vector-kernel time**, which caps any dispatch-overhead
idea at 0.1% before a line is written. That attribution existed before this
experiment was run — it simply was not consulted.

**Rule extracted:** check the existing time attribution before proposing an
overhead-reduction idea. The ceiling is often already measured.
