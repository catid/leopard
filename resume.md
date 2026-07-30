# Session resume state — 2026-07-30, branch `codex/leopard2`

Beads is the durable task source. Use the Beads 1.x binary explicitly:

    /home/catid/.nvm/versions/node/v24.1.0/bin/bd

The legacy `~/.local/bin/bd` is 0.47 and must not touch this checkout.

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
