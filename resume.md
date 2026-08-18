# Session resume state — 2026-07-30, branch `codex/leopard2`

Beads is the durable task source. Use the Beads 1.x binary explicitly:

    /home/catid/.nvm/versions/node/v24.1.0/bin/bd

The legacy `~/.local/bin/bd` is 0.47 and must not touch this checkout.

## 2026-08-18 GF8/AVX2 T16 Q3/Q4 B64 checkpoint

Commit `81c8fef` lands the legacy-wire-compatible packed 64-byte encoder for
`K=33..65`, `R=9..16`.  A generated AVX2 Q3/Q4 circuit fuses the outer inverse
layer, XOR accumulation, and first parity-transform layer.  `K=65` adds the
independently verified coordinate-80 Lagrange column for the same `[128,112]`
parent.  The public scratch contract remains 32 T16 shard rows and no public
API, field representation, coordinate map, or parity ordering changed.

The frozen 25-round confirmation reports exact-Leopard1/candidate ratios of
1.598x at K62/R16, 1.859x at K65/R9, and 1.760x at K65/R16.  Same-binary mature
controls are 1.875x, 2.022x, and 1.918x respectively.  The adjacent T64
K62/R33 repair is 1.089x versus exact main and 1.156x versus its mature control.
Every predeclared inactive predecessor comparison retained a lower CI95 above
0.98.  Evidence is in
`build/regression-high-b64-v23-k65-layout-neutral-confirm25-raw/summary.json`
and optimization-log report 36.

Final gates: Release focused 7/7; CMake/Visual Studio graph audit 165/165 in
normal and Python `-O`; strict Clang 18; feature-disabled Release; and Clang 18
ASan+UBSan+LSan.  The expanded sanitizer matrix passed in 305.12 seconds at
71,680 KiB peak RSS with no swap, and the production-archive sanitizer test
passed separately.  Direct-oracle coverage includes every K33..65 at R16,
every K33..64 at R32 and R33, R9 boundaries, packed/unaligned/batch/fallback,
alias and scratch errors, and a K65/R9 tail-basis regression against Leopard1.

Beads remains the task source.  The next largest measured exact-Leopard1 GF8
encode deficit is K62/R8/B64 (about 0.833x), followed by K66/R16/B64 (about
0.91x).  Keep builds/tests/timings serialized through
`/tmp/leopard-gf8-authoritative.lock`; use Release `-j2` under a 512 MiB VA cap,
runtime tests under 256 MiB, and sanitizer builds `-j1`.  Do not apply an
address-space cap to ASan runtime.

## 2026-08-02 GF8/AVX2 transient-plan and active Walsh checkpoint

Commit `949afe3` makes public one-shot transform plans byte-aware, retains only
the selected execution-route metadata, skips dead tiny/full-loss pruned
schedules, removes redundant presence/selection scans, and routes dense GF8
locator setup through a parent-sized AVX2 Walsh kernel.  Reusable public plans
remain a byte-independent superset; the byte identity is diagnostic/transient
metadata and is enforced by scratch queries, ordinary/scalable execution, and
batch bindings.  The same milestone restores the documented historical global
`LEO2_GFNI_VARIANT` build, including its high-direct fallback, without changing
production GFNI ownership or AUTO dispatch.

Focused GNU Release gates pass 10/10: portable ISA and registration,
backend ops, public API, random differential, locator, dense-plan policy,
low-decoder acceptance, high/low duality, and benchmark JSON.  The latest
failure-atomicity audit found that the diagnostic factory could retain a
caller sentinel on zero bytes or an invalid presence value; the commit clears
the handle after canonical overlap validation and adds both regressions.  That
test also passes under Clang ASan+UBSan.  Expanded known-count locator coverage
performs 162 comparisons over 816,480 coordinates, including every GF8
power-of-two parent and GF16 up to 65,536.  The documented global-GFNI
`bench_leopard2` build and a GF8 K=16/R=8 round trip pass; the default,
in-place-GFNI, high-direct-GFNI, GF8-only, and GF16-only focused variants were
also rebuilt serially.

Pinned same-binary CPU14 ABBA at
`/tmp/leopard2-oneshot-walsh-abba-v1` matched all workload digests in 21/21
cells across N=32/64/128/256 and 64/256/1024-byte shards.  Enabling the active
Walsh kernel improves full one-shot time by 1.095x--1.157x at 64 bytes,
1.062x--1.097x at 256 bytes, and 1.025x--1.054x at 1024 bytes; isolated setup
improves 1.217x--1.939x.  A preliminary exact-main screen leaves only three
sampled 64-byte losses (K/R 24/20, 40/12, and 95/95); because that candidate
was frozen before the final commit, retain it as directional only.

Bead `leopard-79h.38.2.5` remains open until a final clean-source exact-main
ABBA is complete.  Stack ownership for the synchronous ephemeral one-shot plan
was implemented in `395fc3f`, passed focused Release and ASan+UBSan 6/6, then
failed its isolated same-source gate: the three-target aggregate was 0.9881x
(95% CI 0.9860--0.9901), with K/R 24/20 and 40/12 regressing.  Commit `529bfb0`
reverts it; report 26 and
`decoder_dispatch/results/stack_owned_transient_plan_negative_20260802.json`
retain the negative result.  The next step is the final clean-source exact-main
ABBA for the promoted transient-plan/Walsh code.  Builds and tests must stay
`-j1`, memory-capped, and serialized through
`/tmp/leopard-gf8-authoritative.lock`; prior compilation parallelism repeatedly
OOMed the session.

## 2026-08-02 GF8/AVX2 R=1 small-reduction closure

Commit `78a38547b9b11fafa5b8143324066af6d9fac509`, tree
`46a3b65dde6073ec6625af28fbe28af3e839eb2a`, promotes only the four cells
that cleared the predeclared all-lane 95% confidence floor: K=2 at exactly
64, 256, and 1024 bytes uses the terminal reduction, and K=4 at exactly 1024
bytes uses the dense reduction.  K=3, K=4 at 64/256 bytes, and the broad
small-byte coarse thresholds remain deliberately unpromoted.

The frozen pure-AVX2 campaign at
`/tmp/leopard2-r1-small-8ff4ed9-abba-r3-v1` retained 111 cells and 2,664
accepted processes.  All input, parity, and recovered-data digests match;
CPU 14 was pinned with sibling 30 recording zero non-idle jiffies.  Candidate
arithmetic improved the same-source public lanes by 1.16x--2.13x in the four
promoted cells, although exact Leopard main still wins these tiny calls by
0.494x--0.954x for one-shot encode and 0.571x--0.892x for one-shot decode.
The next R=1 work is therefore validation/dispatch/setup, not another broad
reduction threshold.

GNU Release/API gates and the final Clang 18 ASan+UBSan+LSan R1 test pass;
the latter was rerun after the multi-item/backend/GF16 coverage expansion and
completed 1/1 in 13.14 seconds.  Builds remain serialized at `-j1` under
`/tmp/leopard-gf8-authoritative.lock` because earlier parallel compilation
repeatedly exhausted memory.  Bead `leopard-79h.38.5.10.23.1` is closed.

The generic exact-main ABBA runner correctly rejected a final-source attempt
before timing because its v9 contract predates current benchmark-attestation
v4, the production selector defaults, an AVX-512-disabled effective-AVX2
candidate, and a public pure-AVX2 main selector.  No result from that failed
preflight was retained.  Follow-up `leopard-79h.16.2.2` owns the evidence
schema migration.  The next ready codec optimization is P0
`leopard-79h.38.2.5`: reduce tiny GF8 maximum-loss one-shot plan overhead,
starting with a byte-aware one-shot setup that builds only the selected path.

## 2026-07-31 GF8/AVX2 Cauchy-setup checkpoint

Commit `568877cb36264d95aa980be58b6265d56fee704c`, tree
`ccedce4d9f6331f973560c2ac16571b3dd42b7eb`, reuses every direct-repair
Cauchy cross-difference log across its row product, column product, and final
denominator and visits each within-set derivative pair once.  At eight losses
this reduces the relevant difference-log calls from 304 to 120 without
changing coefficients, plan layout, scratch, the ABI, or the wire format.

The frozen identical-text candidate/control campaign retained 2,490 accepted
processes plus 225 fail-closed isolation retries.  Across the broad 60-target
one-shot matrix the Cauchy refactor's same-source geometric mean is 1.0510x,
its worst point ratio is 0.9953x, and every target beats exact Leopard1 with a
minimum lower confidence bound of 1.00056.  The required K=17/K=128, L=8,
1--65-byte cells at reuse 1,2,4,8,16,64 improve by 1.1035x--1.1097x
geometrically and have a minimum exact-main lower bound of 1.31287.

The pre-existing sole failure `(K,R,L,bytes)=(32,17,8,33)` passed a separate
nine-round holdout: control/candidate 1.10550x with lower bound 1.09314, and
exact-main/candidate 1.02774x with lower bound 1.01856.  Release ON/OFF,
strict GCC/Clang reduced-field, ASan+UBSan, runner, benchmark-JSON, random,
direct-oracle, arbitrary-tail, guard, and concurrent one-shot gates pass.
The broader TSan context test exceeded its deliberately imposed 384 MiB
cgroup after the production high-direct TSan test passed; no result from the
terminated test was accepted or retried.

Compact evidence is
`experiments/leopard2/direct_repair/results/`
`cauchy_log_reuse_checkpoint_20260731.json`.  Close
`leopard-79h.38.5.10.44` after the evidence commit.  The next active CPU/GF8
step is `leopard-79h.38.5.10.27`, the final exact-Leopard1 all-K map and
deterministic crossover table.

## 2026-07-31 GF8/AVX2 equal-rounded multi-loss checkpoint

The reusable-plan direct-repair selector is source commit
`00191aff90d8b20b547fa28b4693e3d7b6b4ebcf`, tree
`2db2b9eafaebd24578002788bd2f12bbc2e5bc6e`.  For legacy-high GF8/AVX2 it
selects source-major direct repair when `17 <= K <= 128`,
`ceil_pow2(K) == ceil_pow2(R)`, `K != 65`, and two through eight originals
are missing.  Candidate and initialized-data control executable sections are
byte-identical.

The final 47-target/five-neighbor CPU-13 campaign retained 906 accepted
processes; a predeclared nine-round holdout retained another 54 and resolved
the only inconclusive exact-main cell.  The combined same-source and exact
Leopard1 execution geomeans are 3.9574x and 5.0310x.  Their weakest lower
confidence bounds are 1.3724x and 1.1041x.  All logical original, parity, and
recovery digests match, and all neighbors pass.

At one through 65 bytes, reusable execution wins all 14 cells but first use
does not: exact-main-over-Leopard2 first-use geomean is 0.6228x because plan
setup dominates.  The measured reuse crossover is two through fifteen calls.
Follow-up `leopard-79h.38.5.10.44` owns that explicit one-shot gap.  Compact
evidence is
`experiments/leopard2/direct_repair/results/`
`equal_rounded_avx2_authoritative_20260731.json`.

Close `leopard-79h.38.5.10.30` after the evidence commit.  The next CPU/GF8
step is the final all-K map and deterministic crossover table; do not mark the
overall AVX2 goal complete.

## 2026-07-31 GF8/AVX2 T=8 ragged 65--928-byte checkpoint

The final selector is commit
`a808eac051469330583b7241d2591e30d7b6a354`, tree
`36db13839c8a6d70441cff7a8737f924d55b7125`. Dense legacy-high GF8/AVX2
T=8 profiles now use the prepared binding over qualified ragged tiers from
65 through 928 bytes. Five cells remain explicitly excluded:
`(5,5,191)`, `(6,5,319)`, `(6,6,319)`, `(7,5,319)`, and `(7,6,319)`.

The discovery plus predeclared holdout covered 2,058 cells and 32,004 accepted
processes, yielding 1,213 qualified targets and 845 neighbors. The selector's
same-source geometric mean is 2.0296x with a 1.0504x minimum lower confidence
bound; application-equivalent padded exact-main geometric mean is 1.2851x with
a 1.000026x minimum lower bound.

The final-source confirmation at
`/tmp/leopard2-t8-ragged-a808eac/final-selector-v4-pinned-r1` covered 14
targets and five neighbors in 816 accepted processes. All digests, frozen
binary rehashes, and isolation checks passed, with zero sibling work. Its
same-source/exact-main geometric means are 1.7715x/1.0598x and minimum lower
bounds are 1.0834x/1.00091x. GNU Release selector-on passed 3/3,
selector-off passed 2/2, and Clang 18 ASan+UBSan+LSan passed 2/2.

Compact evidence is
`experiments/leopard2/gf8_high_encode/results/`
`t8_ragged_checkpoint_20260731.json`. Attach the result to
`leopard-79h.38.5.10.43.16` and close that child after the evidence commit.
Keep parent `.43` open for the remaining sub-4-KiB and final all-K work. The
next highest-confidence AVX2/GF8 avenue is equal-rounded multi-loss direct
repair under `leopard-79h.38.5.10.30`.

## 2026-07-31 GF8/AVX2 T=8 arbitrary-tail checkpoint

The production selector is commit
`59a745da8c6ecbc4342c64da428f762c3e16df36`, tree
`4201a2afba6648830bdbeba0a11ae02fed896af4`. Every dense legacy-high GF8/AVX2
T=8 profile with `K=5..16`, `R=5..min(K,8)` now uses the fused prepared
binding at byte counts `1,2,3,7,8,15,16,17,31,32,33,63`.

The frozen three-round campaign covered 504 targets and 84 64/65-byte
neighbors in 10,080 accepted processes. A predeclared nine-round holdout
resolved the only two initially inconclusive exact-main cells. The combined
same-source geomean is 2.9609x with a 1.6586x minimum lower confidence bound;
the padded-64 exact-main geomean is 1.5669x with a 1.0263x minimum lower
bound. Exact Leopard main physically processes 64 bytes for the sub-64-byte
comparison, while all logical input/parity/recovery digests match.

Fresh Release default, feature-off, GF8-only, GF16-only, and Clang 18
ASan+UBSan+LSan gates passed. The final default binary's executable sections
are byte-identical to the timed candidate. Compact evidence is
`experiments/leopard2/gf8_high_encode/results/`
`t8_tiny_checkpoint_20260731.json`.

Attach the evidence to `leopard-79h.38.5.10.43.15`, close that child, and keep
the parent `leopard-79h.38.5.10.43` open for the remaining sub-4-KiB high-rate
GF8/AVX2 gap. Continue to use one substantial process at a time, Release
`MemoryMax=512M` with `-j2`, sanitizer `-j1`, and narrow-test/timing
`MemoryMax=256M`.

## 2026-07-31 GF8/AVX2 T=8 one-kibibyte checkpoint

The final selector source is commit
`cf86e4d2c7d9b2b906a491d5b3312f31c662e57f`, tree
`572676891889dbf6fe4e7a439e2d67cbf218eb2b`. At exactly 1024 bytes the
reusable legacy-high GF8/AVX2 binding now extends the existing direct-input
T=8 route to `(K,R)=(6,5),(6,6),(10,8),(16,5)`.

The final frozen CPU-4 campaign is
`/tmp/leopard2-t8-1k-final-cf86e4d/evidence-final-r9-cpu4`: 42 K/R cells,
672 processes, four targets, 38 neighbors, all workload/parity/recovery
digests equal, zero sibling activity, no target failure, and no neighbor
failure. Same-source gains are 1.0615x--1.0932x and exact-main gains are
1.1396x--1.3028x. The rejected `K=11,R=5` and `K=11,R=6` candidates missed
the five-percent lower-confidence floor and remain on the prior path.

Final focused gates passed: GNU 13 Release 4/4, explicit feature-off 1/1,
GF8-only 1/1, and Clang 18 ASan+UBSan+LSan 1/1. The compact evidence is
`experiments/leopard2/gf8_high_encode/results/`
`t8_one_kib_checkpoint_20260731.json`.

Attach this result to `leopard-79h.38.5.10.43`, commit and push the evidence
milestone, then continue the same bead's still-open 1--256-byte exact-Leopard1
gap or the final all-K map under `leopard-79h.38.5.10.27`. Do not close the
parent bead: the full sub-4-KiB acceptance matrix is not yet complete.

## 2026-07-31 GF8/AVX2 T=4 extended-binding checkpoint

The final production selector is commit
`789d25da08a71f8479df78b2d1cdccbc8d9ef6c8`. Reusable dense legacy-high
bindings now keep the prepared T=4 AVX2 tables through per-`K,R` ceilings up
to 16 KiB. The exact table, merged 49-cell evidence, raw hashes, code-size
delta, and reproduction arguments are in
`experiments/leopard2/gf8_high_encode/results/`
`t4_extended_checkpoint_20260731.json`.

The final-source nine-round same-executable holdout is
`/tmp/leopard2-t4-extended-789d25d-final-holdout-r9`: 8/8 cells accepted,
360 processes, all digests equal, no sibling activity, no route-off
regression, and no Leopard1 loss. The four target same-source gains are
1.062x to 1.096x; exact-main gains are 1.426x to 1.675x. The final focused
Release gate passed 4/4, ASan+UBSan+LSan passed 3/3 at 110,840 KiB RSS, strict
GCC/Clang passed, and reduced-field builds passed. The reduced-field gate also
found and fixed a GF16-only compile failure caused by an unconditional
reference to the GF8 diagnostic selector.

Before taking new work, attach this evidence to
`leopard-79h.38.5.10.35.2`, close it, sync Beads, and commit/push the evidence
files. Then resume the largest remaining exact-Leopard1 GF8 gap, currently
`leopard-79h.38.5.10.43` (sub-4-KiB high encode), using the same-inode
candidate/control method.

## 2026-07-30 OOM-safe source-plan negative result

The default-off `LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN` candidate was screened
from frozen same-source AVX2 binaries after Release and Clang 18
ASan+UBSan+LSan gates passed 4/4.  Preformatting the K=65 multi-loss
source-major schedule improved confirmed execution by only 1.05 to 4.43
percent while making plan setup 1.58 to 1.82 times slower.  Even the best
execution interval had an upper bound below the five-percent promotion gate.
The candidate remains default-off and is rejected; evidence is in
`experiments/leopard2/direct_repair/results/avx2_source_plan_negative_20260730.json`.

The diagnostic used one process at a time, a one-GiB benchmark ceiling, CPU 4,
and at most 42,496 KiB RSS.  Sanitizer tests were serial under a three-GiB
ceiling and peaked at 1,612,700 KiB.  No OOM event occurred.  Continue the
active AVX2 goal through another bounded avenue rather than running the
deferred 16-MiB/reuse-64 campaign on this controller.

## 2026-07-30 OOM-safe AVX2 four-source checkpoint

The generalized GF8/AVX2 one-loss path now isolates its bounded direct-repair
dispatcher out of line, preventing the four-source tiny-shard callback from
changing unrelated decoder layout.  Frozen same-source screens cover 120
selected cells and 102 inactive neighbors.  All selected execution cells won
with a minimum 95-percent lower bound of 1.0595x; every plan-amortized geomean
also won.  No inactive neighbor regressed by more than two percent after
high-confidence confirmation, and the earlier two-source milestone retained
its gains.  Evidence is in
`experiments/leopard2/direct_repair/results/avx2_four_tiny_production_20260730.json`.

Final bounded gates passed: Release 3/3, option-off Release 2/2, Clang 18
ASan+UBSan+LSan 3/3, and GF16-only compilation.  Builds used `-j2` or sanitizer
`-j1`; tests were serialized under one- to three-GiB cgroup ceilings.  Peak RSS
was 1,559,272 KiB and the enclosing cgroup still reported zero OOM events.
Keep the same policy: no subworkers, one native process at a time, and no broad
16-MiB/reuse-64 or complete provenance campaign from the controller.

## 2026-07-30 OOM-safe AVX2 pair-fusion checkpoint

The controlling session reset twice during allocation-heavy work.  Continue
with no subworkers, `-j2` at most (`-j1` for sanitizers), one native process at
a time, and an explicit cgroup memory ceiling for every test or benchmark.
Do not rerun the full provenance suite or the 16-MiB/reuse-64 campaign from the
controlling session.  The kernel still reports zero cgroup OOM/OOM-kill events,
but that is not permission to run unbounded children.

The current dirty tree adds a pure-AVX2 two-source/one-output callback for
generalized GF8 one-loss repair.  It preserves term order and repair
coefficients while halving output passes.  The selector is restricted to seven
measured K/R/profile shapes and conservative byte thresholds.  Frozen
same-source ABBA screens found 1.0831x to 1.5564x execution geomean speedups in
selected cells; the weakest 95-percent lower bound is 1.0790x.  No unselected
neighbor regressed by more than 1.0 percent.  All workload digests matched.
The measured threshold is cached once in each immutable codec, removing a
hot-path table scan that had caused a rejected 2.8-percent neighbor regression.
LOW `(127,128)` was moved from 64 KiB to 1 MiB after a current-host threshold
sweep.  LOW `(32,224)` was rejected after non-monotonic current-source results
missed the confidence gate at 8 and 256 KiB and regressed at 16 KiB.  The
pure-AVX2 object text delta is 3,464 bytes (3.57 percent) with no EVEX prefix.

Release backend/API tests pass both with the generalized path enabled and
disabled.  Focused Clang 18 ASan+UBSan+LSan backend/API tests also pass under
a 3-GiB cgroup ceiling; the largest observed RSS was 1,559,312 KiB.  A
GF16-only production build passed with 178,832 KiB peak RSS.  Evidence is
staged in
`experiments/leopard2/direct_repair/results/avx2_pair_production_20260730.json`.
The immediate safe next step is review, Beads update, commit, and push.

## 2026-07-30 low-memory AVX2 promotion checkpoint

After two further controller restarts, keep this campaign in a strict
low-memory lane even though the host currently reports about 88 GiB available
and the enclosing cgroup still reports zero `oom` and `oom_kill` events:

- do not launch subworkers from the controlling session;
- build with `-j2` at most, or `-j1` for sanitizer builds;
- run only one test or benchmark process at a time;
- put native sanitizer/test children in a systemd scope with a 4 GiB
  `MemoryMax`; and
- do not run the complete provenance suite until
  `leopard-79h.40.36.12` resolves its controller-restart behavior.

The current working tree promotes generalized GF8/AVX2 one-loss repair and
specializes plan construction in logarithmic form.  An immutable codec cache
stores one vanishing-polynomial logarithm per transmitted parity coordinate;
plan setup no longer materializes a generator row, performs a generic 1x1
elimination, or converts that row from field elements back to logarithms.
The generic solver remains the fallback.

Frozen same-source evidence covers 864 broad cells and 110 counterbalanced
neighbor cells.  Plan-amortized direct repair won every broad cell (minimum
1.0311x, geomean 5.1139x); the 110-cell ABBA screen had a 1.2168x minimum.
Against the exact Leopard `master` adapter, all 48 reusable-codec/plan cells
won with a 1.6951x minimum.  Completely cold context+codec+plan+one-execution
timing retains two slight K=254/R=2/64-byte losses, documented in
`experiments/leopard2/direct_repair/results/general_one_loss_promotion_20260730.json`.

Fresh default Release AVX2 tests passed the direct oracle, backend operations,
and API suites.  A serial Clang 18.1.3 ASan+UBSan+LSan API build/test was
bounded by a 4 GiB cgroup and passed; build peak RSS was 633,716 KiB and test
peak RSS was 1,362,656 KiB.  The Visual Studio graph and benchmark-attestation
refresh tests also passed.  Commit/push and Beads closure remain the immediate
next steps.

## 2026-07-30 restart-safe checkpoint

Do **not** run the complete `tools/leopard2_build_provenance_test.py` suite
from the controlling Codex session.  Two attempts immediately terminated and
restarted that session.  The enclosing cgroup subsequently reported
`memory.events` values `oom=0` and `oom_kill=0` with `memory.max=max`, so the
failure is not proven to be host-memory exhaustion; process-containment fault
injection may instead be affecting the controller.  The isolation work is
tracked by `leopard-79h.40.36.12`.  Until that bead is complete, use a separate
session/process boundary plus strict address-space, RSS, and time limits for
individual provenance tests only.

The current bug-first working tree contains:

- a default-OFF generalized GF8/AVX2 one-loss direct-repair candidate and its
  four-source tiny-shard callback;
- exact v3 CMake selector attestation propagated through the benchmark/evidence
  runners;
- retained-descriptor and terminal-precedence fixes in the lab, fuzz, and
  provenance tooling;
- parallel fuzz campaign/audit v5 with a signed, closed-world probe execution
  policy; and
- regression coverage for recycled numeric file descriptors, process-tree
  teardown, historical evidence contracts, and current schema closure.

Safe validation completed after the restart:

- `tools/leopard2_lab.py self-test`: PASS under normal and optimized Python,
  peak RSS 58,368 and 58,112 KiB;
- `tools/leopard2_fuzz_campaign.py self-test`: PASS under normal and optimized
  Python, peak RSS 53,760 and 53,248 KiB;
- Python compilation and `git diff --check`: PASS.

The older statement below that the pre-optimization bug hunt is complete is
historical.  The current gate remains open until the saved working tree is
committed/pushed, memory-safe compiled correctness checks pass, and
`leopard-79h.40.36.12` explains the unsafe provenance-suite launch.

## 2026-07-29 bug-first resume gate

The pre-optimization bug hunt is complete.  Nine reviewable commits capture
the previously uncommitted AVX2 experiment and evidence work:

- `4de69ba` — default-off small-direct AVX2 dispatch/executor introspection;
- `9cf3f67` — compiler-replay transport and semantic compile closure;
- `96c67cf` — exact-main immutable Git/executable/evidence capture;
- `fca72cc` — all-K identity-bound descriptor and atomic evidence ownership;
- `4e25515` — bounded small-direct ABBA/exhaustive experiment runners;
- `8b91f0a` — CMake experiment, sanitizer, architecture, and attestation model;
- `30e61f2` — direct-encode crossover evidence hardening;
- `f9a22c1` — balanced-promotion evidence validation; and
- `046158f` — remaining benchmark helper fail-closed boundaries.

Concrete bugs fixed in the final pass included:

- PairLease and Git-symlink descriptor transfer gaps;
- all-K evidence-directory acquisition/close interruption gaps;
- benchmark-attestation mutation through a symlinked output ancestor;
- quoted/custom/multi-config sanitizer misclassification;
- Apple target-architecture selection;
- partial PID-marker publication and short-lived lab observation fixtures;
- stale CMake attestation identity; and
- a structural-verifier error that collapsed exact `STREQUAL` comparisons
  against empty/zero into truthiness, incorrectly equating `OFF` with an empty
  build type.

Fresh validation at the committed source content:

- GNU 13 Release: 113/113 CTests passed in 423.22 seconds using `-j30`;
- exact-main runner: 125/125 under normal Python and 125/125 under `-O`;
- all-K runner: 78/78 under normal Python and 78/78 under `-O`;
- provenance: 71/71 on Python 3.12 and 3.13, normal and `-O`;
- Visual Studio/CMake graph verifier: 145/145 normal and `-O`;
- Clang 18 ASan+UBSan mode 2: API, dense-policy, and high/low-duality tests
  passed; the exhaustive GF8 direct-repair executable checked 1,982,812
  erasure patterns / 162,888,934 basis symbols without a sanitizer finding;
- benchmark schema v6 reported the expected `source_major` executor for
  K=8, R=8, L=5, 65-byte shards; and
- both Python 3.12/3.13 compilation and `git diff --check` passed.

The independent C++ audit found no correctness, UB, bounds, aliasing,
concurrency, or backend-dispatch defect in the current direct-repair changes.
The experimental small-direct selector remains disabled by default.

One independent audit found a distinct, still-open authoritative replay gap:
`leopard-79h.38.5.10.38.2.3.15.4.1`.  The selected Make target can still reach
an added literal recipe, `cmake -E env` can launch an unretained helper, and a
mutated reachable `DependInfo.cmake` can execute CMake language.  Do not treat
clean-replay timings as authoritative until the complete selected
target/prerequisite/recipe and executable CMake-input inventory is bound.
This does not invalidate codec correctness tests or frozen-binary diagnostic
screens.

Resume the active AVX2 GF(256) goal with
`leopard-79h.38.5.10.39`, `.39.2`, and `.39.2.2`: prototype and screen
one-loss direct-L1 fusion against Leopard1 using lane-owned frozen binaries.
Keep the default production build unchanged until a candidate meets the
correctness and neighboring-regression gates.

## P0 pre-resume bug gate completed

`leopard-79h.49` and all fifteen children are closed. The reviewable code
milestones are:

- `a9e0e3e` — defined `NextPow2` edge behavior, overflow-safe GF16 live-set
  policy, warning-clean reduced-field builds, and backend-correct scalar
  validation;
- `a8881a2` — fail-closed lab process/TID/affinity containment and real runtime
  observation evidence;
- `cf96061` — build-time benchmark source attestation, exact current/historical
  compile provenance, and shadow/nested/superproject defenses.

Final validation at those commits:

- GNU 13 Release: 104/104 CTests passed twice;
- Clang 18 ASan+UBSan: 36/36 focused production tests passed;
- portable strict dual-field build: 97/97;
- strict GF8-only and GF16-only: 7/7 each, warning-free;
- forced-scalar high-pruned regression: 1/1;
- run_abba: 56/56; Visual Studio/CMake graph: 134/134;
- lab, benchmark-matrix, balanced-promotion, attestation-refresh, CUDA-optional,
  ISA-isolation, and field-option self-tests passed;
- three independent review partitions and a final repeat found no unresolved
  actionable defect.

Disposable validation builds are under:

- `/tmp/leopard-bughunt-final-release-1785263016`
- `/tmp/leopard2-bughunt-asan-ubsan-1871043`
- `/tmp/leopard2-bughunt-matrix-1876619`

## Affinity-aware L3 tiling milestone

`leopard-79h.46` has reached its evidence-backed implementation point.  The
fail-closed evidence runner is committed as `f674e57`.  The proportional
candidate at `767161c` was validly measured and rejected: on the 96 MiB cache
domain it slowed the 96 MiB live-set cell by 26.1% encode and 20.7% decode and
slowed 256 MiB encode by 5.9%.

The promoted hybrid candidate is commit
`cb2ef79cc180166326acd80a3d798392ce5157fe`, tree
`9d48810b25b227db8a27e3081be6ccd48e3276a6`.  It:

- snapshots the minimum L3 visible to an explicitly requested AVX2 context;
- clamps the policy to the measured 32--96 MiB range;
- retains the calibrated 16 MiB execution target;
- uses `max(effective L3, 64 MiB)` as the generic crossover; and
- keeps measured side overrides cache-gated and scales their pass widths.

The reviewed production/model integration is `1a2e7c2`.  Supporting bug-fix
commits are `1c83b9d` (single-command benchmark recipes), `4090bee` (stable lab
runtime fixtures), and `45fb2ab` (defined FF16 reference-helper tail).

An authoritative eight-cell F-to-H campaign on CPU 4/sibling 20 accepted the
hybrid with no credible regression over 2%.  It retained 1.1066x high
side-512/64-KiB encode and 1.0946x maximum-loss encode gains.  An independent
CPU 14/sibling 30 campaign measured a 1.1855x decode improvement for LOW
K=200, R=800, loss count 9, 64 KiB, explicit AVX2, and the 32 MiB policy.  The
final bug audit found that the campaign exercised only one erasure mask, while
other nine-erasure masks produce materially different pruned schedules.  The
LOW diagnostic is therefore retained but not promoted.  Production keeps the
generic LOW rule pending a counterbalanced multi-pattern campaign.

Evidence and disposition are tracked in:

- `experiments/leopard2/optimization_log/25-affinity-aware-l3-tiling.md`;
- `experiments/leopard2/l3_tiling/results/checkpoint.json`; and
- the generated `/tmp/leopard-l3-evidence-*` bundles named by the checkpoint.

The operation-count model is now schema v4.  It mirrors cache-derived GF16
passes and the GF8 side table, reports independent plan/codec pass geometry,
uses maximum balanced-pass bytes for scratch slots, preserves actual batch
execution rows separately from conservative scratch rows, and assigns zero
data-slot width to no-op/direct/K=1 terminal paths.  Private build-tree
introspection obtains the same geometry from production codec/plan objects
without changing the public C ABI or installed header.

The final integrated correctness gate passed:

- schema-v4 pure tests: 23/23 under normal Python and 23/23 under `python -O`;
- operation-model self-test: 397 checks under each Python mode;
- production-linked scratch/selector/cache-policy cross-check passed;
- two independent 100,000-case randomized model/production-invariant sweeps
  passed;
- current GNU Release: 104/104 CTests;
- current focused Clang ASan+UBSan: 5/5, with the immediately preceding full
  integrated sanitizer build at 97/97;
- strict GCC, strict Clang, GF8-only, and GF16-only builds passed, with 7/7
  reduced-field tests in each field build; and
- the lab containment fixture passed in the full suite and 40/40 times while
  unrelated workloads saturated every allowed CPU.  The fake performance-
  counter fixture independently passed 5/5 under the same full-CPU contention.

The lab fixture had one real scheduling-dependent test bug: an escaped child
lived for only 0.4 seconds, so a loaded supervisor could miss it.  The child
and active peer now remain observable until quarantine, while a ten-second
manifest timeout bounds a failed run.  A second 80-ms fake-counter workload
could likewise finish before its child was observed; it now remains visible
for two seconds and failure output preserves the full runtime/counter detail.

The final bug gate also:

- rejected impossible GF8 parent lengths in the direct scratch-model helper;
- withdrew the LOW side-256 override whose single-mask evidence did not support
  all nine-loss schedules;
- rejected signed multi-command executable recipes in producer, offline
  replay, and normalized benchmark-promotion scope (57/57 normal and optimized
  Python runner tests); and
- initialized the complete FF16 butterfly reference buffer, keeping the test
  helper defined for future odd-byte adversarial probes.

The final static/model/API/diff reviews found no unresolved actionable defect.

## Next production optimization priorities

The active goal is AVX2 GF(256), so after `.46` continue in this order:

1. finish the `.38.2.4.1` evidence-validator closure, then use `.38.2.4` to
   complete the current balanced full-loss generic-dispatch matrix;
2. `leopard-79h.38.5.10.39` and `.39.2` — restage and evaluate the broader
   one-loss direct-repair/pair-fusion region against transforms and Leopard1;
3. close `leopard-79h.42` as a measured negative unless new at-scale sparse
   rows establish a distinct useful region;
4. `leopard-79h.43` — exact-main campaign, requiring a bare terminal with no
   Codex/Claude session alive because the fail-closed affinity supervisor sees
   session process churn;
5. `leopard-79h.44` — explicit-GFNI campaign and possible AUTO promotion;
6. `leopard-79h.47` — second-microarchitecture validation;
7. `leopard-79h.48` — final verified Leopard1 comparison documentation; and
8. return to GF16-only `leopard-79h.50` and multi-pattern follow-up `.52` after
   the active GF(256) goal.

`leopard-79h.45`, the GF16 table-footprint milestone, is closed.  Never use
screen measurements in the published Leopard1 comparison.  Only verified,
immutable `run_abba.py` evidence bundles are authoritative.

## Host and coordination

The allowed CPU set currently contains 30 logical CPUs:
`0-14,16-30` (15 SMT pairs), one NUMA node. Cap parallel build/test work at
30. Serialize builds, correctness runs, and timings with:

    /tmp/leopard-gf8-authoritative.lock

Authoritative timing additionally requires frozen lane-owned binaries and the
provenance checks in `AGENTS.md`. The exact-main campaign remains a bare-
terminal-only task on this host.

## 2026-08-18 GF8 R=32 regression closure

The current-source regression line `leopard-79h.38.22` reached its acceptance
point after the K62/R33, T16 Q3/Q4, K62/R8, K66/R16, and T32 Q3 boundary
changes.  The authoritative decision is the standalone 75-round campaign in
`build/regression-current-atlas-final97-970107e-confirm75-v1`; the earlier
25-round v2 bundle is a rejected pilot and must never be pooled with it.

The accepted candidate is commit `970107e9c295ff3b4cec601836ba1df17dee7564`,
tree `0ea848d4a9eb5c56f4b85ec6ecf7c7bd784a0a39`, executable SHA-256
`a8d9fe801cfe1662e8548894d37b7c5fdabf9686abe0a8e2e4bca4ec68dfec4d`.
The exact-main executable SHA-256 is
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.

The campaign covered 97 predeclared GF8 legacy-high AVX2 `R=32` workloads,
64- and 1024-byte shards, and both encode and setup-inclusive decode.  All 194
exact-main comparisons passed and all 194 identical-binary confidence
intervals stayed inside `[1/1.02,1.02]`.  The weakest exact-main lower bounds
were 1.031403x for encode and 1.054421x for decode.  Five attempts with one
reserved-sibling non-idle jiffy were retained as discarded evidence and
retried; all 7,275 accepted rounds recorded zero sibling activity.

Evidence identities:

- raw SHA-256 `bcebaf2a3f9a60e0d6d39ba27cc8a1ba0c6841390db0469cca676e1c55634dd7`;
- summary SHA-256 `27c5cdc7927300ecd3e6da94fbbd7b21c3da0df3cc11fe9b5323059bdaac259d`;
- ordered artifact SHA-256 `39e6afd3f8a6cb0ce30d8ee339ff71b61179d07664dcb7377ea6a5bdfdeaaa25`;
- matrix SHA-256 `fe21e525ff084564af37a4e7b43f17df8e971f96a252199125e26d50e445da8b`.

The compact checkpoint and reproduction note are
`experiments/leopard2/gf8_high_encode/results/current_atlas_final97_confirmation_checkpoint_20260818.json`
and
`experiments/leopard2/optimization_log/38-current-atlas-final97-regression-closure.md`.
The 541 MiB raw bundle remains ignored build evidence; validate it by the
hashes above when present.  Follow-up resumable journal work is tracked in
`leopard-27x`; broader open work remains under the Beads graph rather than this
file.
