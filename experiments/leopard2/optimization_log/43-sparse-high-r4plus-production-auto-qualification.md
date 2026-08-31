# 43 — Sparse-high R4+ production-AUTO qualification

## Result

The fresh mechanism-bounded production-AUTO campaign did **not** meet its
preregistered gate. All 88 cells passed acquisition, correctness, route,
provenance, and isolation validation. All 56 R in `{4,8,16}` candidates cleared
both five-percent marginal lower-confidence-bound requirements, and all 18 R=2
structural controls retained the required transform route. One of the 14
performance controls nevertheless crossed the fixed minus-two-percent net
point-regression guard.

This is a valid negative promotion decision, not an acquisition or candidate
performance failure. The compiled sparse-high AUTO default remains `OFF`; no
production route, public API, installed behavior, or frozen decision rule was
changed after observing the data.

The decisive control was caller-AUTO/thread-four, binding API, batch 16 at
`K=2,R=16,4096`. Its two route arms both retained transform execution, but its
setup-inclusive net point estimate was -2.982%, beyond the -2% veto. The
control is intentionally sensitive to noise or uncontrolled differences
between supposedly equivalent routes. Its result cannot be waived after the
fact merely because every actual candidate passed.

## Accepted campaign

The source was detached and clean at commit
`09ef2c6191671e72874828384eb73a5e7fcdc5d7`, tree
`1fe65a93b9e79ed11795d9f156bbca8435749440`, with initialized `sse2neon`
gitlink `cad518a93b326f0f644b7972d488d04eaa2b0475`. The 1,201-file source
fingerprint was
`e16688fb7d692f633e839092cdbdbb76327bbdf31dfd3191b5021bb18bc36296`;
the worktree and before/after source identities remained unchanged. The runner
SHA-256 was
`93e4c2037e70e538281a8c848418053e03a3f8374391ae46d20df09493b5610e`.

The runner created one tests-off Release/Ninja explicit-AVX2 ordinary archive
with GF8 and OpenMP enabled, GF16 and full-output high direct disabled, sparse
table preparation enabled, and sparse AUTO compiled-default off. It rejected
test-hook linkage and froze the executable and ordinary archive together:

- controlled-build document:
  `6a0838999f63b6cb5c70b3ddddc10d2c25897f952d4aec1fc90b09a2e9f2cffb`;
- benchmark executable:
  `8760337e91224d6f946a0b0a5f6847ed8c4346387e86b0c1e1f409b966a160ae`;
- ordinary `libleopard.a`:
  `86bf8861b8acf398f61e1343f2364564f8eae58492cddf66349b6aba40a74782`.

Timing ran on `foureyes.lan`, an AMD Ryzen Threadripper PRO 9985WX, with every
benchmark process pinned to CPU 109 and exact SMT sibling CPU 45 reserved.
Both are core 45 on socket 0. A housekeeping-pinned 20-second preflight
observed zero non-idle jiffies on both CPUs. The runner held the canonical
campaign lock and exclusive pair-wide lease. Every one of the 88 accepted jobs
recorded exactly zero reserved-sibling non-idle jiffies; benchmark-CPU non-idle
work ranged from 25 to 110 jiffies. No attempt was discarded or resumed.

The frozen 88-cell v11 grid has SHA-256
`5cb8acb9cdb2fc00de22753c9789f59bed4b515d0f62a663d45c2934c3e8b627`.
It contains 56 candidates, 18 R=2 structural controls, and 14 performance
controls. Each cell used three rounds of two ABBA contrasts, 15 timing
iterations after four warmups, 15 setup iterations, four API calls per sample,
modeled reuse 64, one worker, and a 120-second per-process timeout. The exact
inventory contains 88 job documents and 2,112 raw benchmark documents, with 24
serialized invocations per cell.

## Decision evidence

A candidate qualified only when both the production-route and summed-net
marginal 95-percent lower bounds reached five percent. Only the 14
requested-backend/thread controls received the minus-two-percent net point
guard. R=2 timings were retained but could not create a performance veto.

| Cell | Route gain (95% interval) | Table gain (95% interval) | Net gain (95% interval) | Decision |
| --- | ---: | ---: | ---: | --- |
| one-shot `K=8,R=4`, 4096 B | +31.68% (+12.67% to +53.89%) | — | +31.18% (+10.92% to +55.13%) | weakest candidate route lower bound; pass |
| one-shot `K=16,R=4`, 4096 B | +21.60% (+15.19% to +28.36%) | — | +20.23% (+10.90% to +30.34%) | weakest candidate net lower bound; pass |
| explicit-AVX2 batch 16 `K=2,R=16`, 4096 B | +0.17% (-8.81% to +10.02%) | -2.10% (-9.30% to +5.68%) | -1.94% (-4.91% to +1.13%) | within control point guard |
| AUTO/thread-four binding batch 16 `K=2,R=16`, 4096 B | -1.93% (-6.61% to +2.99%) | -1.07% (-4.14% to +2.09%) | **-2.982%** (-6.70% to +0.88%) | control point guard triggered |

Aggregate candidate point gains were at least 21.41% for the route contrast and
20.23% net. Their weakest lower bounds were 12.67% and 10.90%, respectively.
All 24 qualified-tuple cells, 24 public-API cells, and eight parity-row cells
qualified separately; no aggregation rescued a failing candidate.

The AUTO/thread-four batch-16 control had the widest net interval and worst
lower endpoint (-21.55%), but its point gain was +0.199% and it correctly did
not trigger the guard. The frozen control rule uses the net point estimate, not
its interval; only the -2.982% binding control crossed that rule.

All 18 padded-side-two controls passed the independent parity oracle, source
immutability, canary, route witness, and raw-schema checks. Under every policy
they reported zero prepared direct rows, no direct witness calls, and transform
execution. Their timing point estimates ranged below and above zero, but were
excluded from the performance veto exactly as preregistered.

## Correctness and sanitizer gates

A separate clean exact-commit lane on `ripper.lan` used GCC/G++ 13.3.0 and held
that host's canonical lock through configure, build, and test. Release passed
6/6 focused tests in 91.54 seconds. ASan+UBSan passed 5/5 executable tests in
211.31 seconds with leak detection, strict-string checks, and sanitizer
halt-on-error enabled; no ASan, LeakSanitizer, or UBSan diagnostic appeared.
Source and submodule identities and all pre/post artifact hashes matched.

The owner-only frozen correctness root is
`ripper.lan:/home/catid/leopard-evidence/correctness-09ef2c-ripper.LvcwzG`.
Its 58-file critical manifest verified 58/58 and has SHA-256
`2e9fd975997edff50a30363c884568fdf54ed9ec7daaba08841c97a1fcd5b4b9`.
The Release and sanitizer test-log hashes are, respectively,
`4ba9bf0e184beaef32c6713b696b066c9241d6e401d682b60b62a6944c42ee91`
and `0c1c1c414021e0b845ff84c768569e32c93a11292df96a8d8d6ecaea9f4192a8`.

## Replay and evidence lineage

Normal and assertion-disabled replay both returned decision status 2 and wrote
outside the result tree. Their analysis documents are byte-identical to each
other and to the acquisition-time retained analysis, with SHA-256
`23b7fd64830125b08acbfd6dafe1665cfd83d5d8558212ae95a82d5b7f9c7294`.
The accepted top-level result hashes are:

- manifest: `6fa2a837d03e180701242cc65d5b862d1aaefff79b955d7cc85b4220f3b4a51e`;
- matrix: `85a22b51a6c461c53d75f367deceafb79682643a26767b64b6e5fdacbdba32c0`;
- analysis: `23b7fd64830125b08acbfd6dafe1665cfd83d5d8558212ae95a82d5b7f9c7294`.

The canonical 6,490-file result-tree inventory was byte-identical before and
after both replays, with SHA-256
`e9c8cf93288c23f0fb4b68a29c040810849e7e88d2f9d8027eb64d048db16b18`.
The live accepted root is
`foureyes.lan:/tmp/leopard-v11-09ef2c6.I0N59X/result`.

The repository checkpoint is
`experiments/leopard2/direct_encode/results/production_auto_r4plus_avx2_checkpoint_20260831.json`.
Its SHA-256 is
`7c6055eb6e7b9701728acc4ead5795794517c561090bff3847bc138e5c3acac7`.

An owner-only immutable archive of the result, replay validation, and original
coordinator logs is at
`foureyes.lan:/home/catid/leopard-evidence/production-auto-r4plus-09ef2c-20260831.SPOW4C`.
The 5,000,725-byte archive SHA-256 is
`8498573ef52bd80ebd7fa5a1cee024c8ef60ab7b00faf68b484691655b2fd195`;
its three-entry critical-manifest SHA-256 is
`0f34e7ec3732e78b7ef60550ddef9e68d9aff3a02f030c9400f4b0631480076e`.

The tmux wrapper's separate, non-evidence exit marker contains literal `2n`
rather than `2` plus newline because of coordinator quoting. The original bytes
are archived unchanged and disclosed in the observation record. This defect is
outside the result tree and decision contract: acquisition stdout, retained
matrix/analysis, and both independently captured replay exit statuses all agree
on decision status 2.

## Independent review

Two independent read-only Codex reviewers converged without a blocker. The
decision reviewer reconstructed every three-round paired-log contrast and
Student-t interval from all 2,112 raw records, reproduced the 56/18/14 role
split and every region aggregate, and confirmed that exactly one eligible
performance control crossed the fixed veto. The integrity reviewer separately
validated the exact source and controlled-build provenance, the complete job,
raw, and log inventories, all route/correctness/parity witnesses, byte-identical
normal and assertion-disabled replay, and every member and manifest entry in
the frozen archive.

Claude was not invoked: the user required a Claude Max/OAuth path with no API
billing, and no unverified billing path was used.

## Disposition and scope

The conclusion applies only to this frozen explicit-AVX2 host, OpenMP runtime,
grid, and caller set. The intervals are per-cell marginal two-sided 95-percent
Student-t intervals with two degrees of freedom, not a simultaneous guarantee.
No measurements are pooled with the v10 corpus, and no control is deleted or
reclassified after seeing its timing.

The R>=4 mechanism remains compelling candidate evidence, but the complete
promotion gate failed. Sparse-high AUTO therefore remains compiled-default
`OFF`. Any future attempt to change the control rule, increase replication, or
qualify another caller/runtime must use a newly versioned preregistered campaign
rather than reinterpret this negative result.
