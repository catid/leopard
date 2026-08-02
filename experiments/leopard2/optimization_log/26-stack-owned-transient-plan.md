# Stack-owned transient one-shot decode plan

## Hypothesis

`leo2_decode()` builds and destroys a temporary immutable decode plan around
one synchronous execution.  Supplying a fresh, default-constructed plan object
from the wrapper's stack should remove one `new`/`delete` pair without changing
the plan's vector allocations, execution, scratch, wire format, or reusable
public-plan lifetime.

## Candidate and correctness

Candidate commit `395fc3fb61fecd7fa3de53cc72fd2e3d13a95db3` added
conditional ownership inside `DecodePlanCreateInternal`; public and diagnostic
factories remained heap-owned, while only the synchronous one-shot wrapper
provided stack storage.  An independent failure-path audit found no escaping
pointer or ownership defect.  Focused GNU Release and Clang ASan+UBSan gates
each passed 6/6: one-shot no-loss, forced-hook no-loss, public API, random
differential, dense-plan policy, and high/low duality.

## Isolated same-source result

The clean pre-change control was commit `f11ce43e245886af8ec087cb1bf70de694dbb5b4`.
Both revisions were built serially with GCC 13.3, Release, explicit AVX2
selection, benchmarks enabled, tests/fuzzers/CUDA disabled, then copied with
their archives into a read-only artifact directory before timing.

- control executable SHA-256:
  `b1af83bf64b3eea63c78e34fab8415ec24536229e0649aae150a4d450ce25547`
- candidate executable SHA-256:
  `e8f61bb14942adb6a88be2020259796ddcf82785e0abd2c3ea065ee9d1b4d428`
- control archive SHA-256:
  `1281783ac1bb39505314026b3ef02e68128066e9dfc6b578e20c92168fe1397b`
- candidate archive SHA-256:
  `a27b225f421b80bca5049552414a49e6b384b32325b55940d1308fa5ca5eb64a`

The campaign ran three ABBA/BAAB rounds on CPU 14, with sibling 30 required to
gain zero non-idle jiffies for every accepted process.  Each invocation used
41 retained samples, 16 warmups, reuse 64, a fixed seed, GF8, legacy-high,
explicit AVX2, 64-byte shards, and the public one-shot timing metric.  All 84
processes passed without retry; schema, resolved identity, scratch, round trip,
and all workload digests matched.  Ratios are control/candidate, so values
above one favor the stack candidate.

| Cell | one-shot ratio | 95% CI | Result |
| --- | ---: | ---: | --- |
| K=24 R=20, loss 20 | 0.9807x | 0.9709--0.9907 | regression |
| K=40 R=12, loss 12 | 0.9711x | 0.9649--0.9774 | regression |
| K=95 R=95, loss 95 | 1.0128x | 1.0067--1.0190 | small win |
| K=12 R=10, loss 10 control | 0.9864x | 0.9795--0.9934 | regression |
| K=48 R=40, loss 40 control | 0.9866x | 0.9728--1.0005 | inconclusive/regressive |
| K=106 R=53, loss 53 control | 1.0203x | 1.0188--1.0217 | small win |
| K=24 R=20, one-loss bypass | 1.0178x | 0.9948--1.0412 | neutral |

The aggregate of the three target cells is **0.9881x**, 95% CI
**0.9860--0.9901**.  This is below the predeclared 1.02x kill threshold and
contains two clear regressions, including one beyond the neighboring-regime
2% limit.  The candidate was therefore reverted by commit `529bfb0`; no
exact-main follow-up was run because the same-source experiment had already
disqualified the mechanism.

Removing a visible allocation was not sufficient evidence of a win.  The
cause of the mixed result was not isolated; plausible contributors include the
larger wrapper frame and code-layout changes, while the allocator's thread
cache makes the removed object allocation cheap.  Do not retry this exact
stack-object form without a distinct mechanism and a kernel-level witness.

Machine-readable summary:
`experiments/leopard2/decoder_dispatch/results/stack_owned_transient_plan_negative_20260802.json`.
