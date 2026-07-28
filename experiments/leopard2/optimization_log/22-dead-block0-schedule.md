# 22 — Dead block-0 pruned input schedule (open finding 1)

**Disposition: SHIPPED — decode plan setup 1.06x-1.15x faster on native-high,
execution and digests untouched**

## The finding
Open finding 1, from the source audit: `PreparePlanExecutionMetadata` compiled
a pruned inverse schedule for block 0 on essentially every native-high
transform decode plan, and **all four Algorithm 5 executors discard that entry
unexecuted** — block 0's inverse is elided outright by the
`FFT_0(IFFT_0(F_0) + H_later) = F_0 + FFT_0(H_later)` cancellation.  Commit
`bd13ef7` moved the executors to block 1 and added the discard guard but never
updated the planner.  Cost: a `Theta((T/2) log2 T)` sparsity-independent
compile (2,304 raw butterflies at T=512) plus up to ~27 KB of retained plan
storage per plan lifetime, with no consumer.

## Safety argument, verified in source before editing
The executors' discard is **conditional** — `input_plans[index].block == 0`
(`LeopardFF8.cpp:4554`, `LeopardFF16.cpp:3878/:4264/:4538`) — so a plan that
simply never contains a block-0 entry flows through unchanged; the guard
becomes dead-but-harmless.  Had the discard been an unconditional
`++input_plan_index`, this change would have skipped a live block-1 schedule —
which is why the executor sites were read before the one-line planner edit.
The low-profile branch keeps block 0: its executors genuinely run it.

## Result
A/B, identical trees except the one hunk, AVX2 backend, best-of-3 x 7
iterations, all recovered digests matching:

| Cell (1 KiB unless noted) | plan setup before | after | speedup | decode exec |
| --- | ---: | ---: | ---: | --- |
| K=224 R=32 | 9.71 us | 8.46 us | **1.148x** | unchanged |
| K=240 R=16 | 5.62 us | 4.99 us | **1.126x** | unchanged |
| K=4096 R=512 | 283.90 us | 257.74 us | **1.102x** | unchanged |
| K=300 R=100 | 32.26 us | 29.71 us | **1.086x** | unchanged |
| K=2000 R=500 | 210.25 us | 194.00 us | **1.084x** | unchanged |
| K=1000 R=200, 64 KiB | 98.78 us | 91.42 us | **1.081x** | unchanged |
| K=1000 R=200, 4 KiB | 98.07 us | 91.16 us | **1.076x** | unchanged |
| K=1000 R=200 | 96.37 us | 90.87 us | **1.061x** | unchanged |

Setup-only, as the finding predicted: decode execution deltas are within noise
in both directions.  The win is largest at small sides, where the dead compile
was the biggest share of setup, and it also removes the retained block-0
schedule bytes from every native-high plan's lifetime footprint.

This closes open finding 1 in
`docs/leopard2_open_optimization_findings.md`.
