# Actual AUTO GFNI boundary qualification

Bead: `leopard-79h.38.5.4.17.1`. Preregistered 2026-09-09, before any
new candidate timing. This qualifies the actual default-off implementation
at `35f53fb5c12f3f336f5a0e6235c5c243993a604a`, not the already consumed
explicit-GFNI screen. No old timings are reused.

## Candidate and correctness

The boundary-only control adds K=1000/R=200/32-KiB and
K=1000/R=199/64-KiB to the existing R=200/64-KiB AUTO GFNI route. The old
target remains enabled in both modes. R=199/32-KiB remains inactive. All
existing host/model, thread, field, profile, layout, flags, full-output and
API guards remain. Explicit backend requests retain their meaning.

The candidate's 60 focused Release and both-field ASan/UBSan/LSan checks
are retained in `auto_gfni_boundary_candidate.md`. The new dual-linked driver
adds 40 untimed checks across all eight workloads below: all 134,938,624
comparison bytes match independently linked Leopard1, and all 16 off/on
sanitizer records match Release. The actual public-encode route probe is
checked and normalized to production mode before clocks; no timed route
counter updates are permitted.

Only small drivers were newly compiled. Release archive SHA-256 is
`80bffc9e873585d9a18fcf6a294413b8c2a76d84fa6e7f9a437fb96ca8458533`;
Leopard1 remains exact commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.
The JSON plan pins both executables, archives, expected records and candidate
source/header. Native build peak was 177,090,560 bytes / 512 MiB; preflight
peak was 141,082,624 bytes / 256 MiB. All six memory-event counters and swap
were zero. Seven pure collector tests pass normally and under Python `-O`.

## Fixed protocol

The plan, driver and collector are committed and pushed before measurement.
Fresh lane-owned readonly inputs bind to that preregistration commit and are
rehash-checked before and after every child. One attempt only; never resume,
overwrite, retry, change CPU, trim or pool partial attempts.

| Cell | K / R / shard bytes | Public API | Role; off/on routes |
| --- | --- | --- | --- |
| 0 | 1000 / 200 / 32768 | ordinary encode | target; AVX2/GFNI |
| 1 | 1000 / 199 / 65536 | ordinary encode | target; AVX2/GFNI |
| 2 | 1000 / 200 / 32768 | ordinary one-item batch | target; AVX2/GFNI |
| 3 | 1000 / 199 / 65536 | ordinary one-item batch | target; AVX2/GFNI |
| 4 | 1000 / 200 / 65536 | ordinary encode | unchanged old target; GFNI/GFNI |
| 5 | 1000 / 199 / 32768 | ordinary encode | unchanged hole; AVX2/AVX2 |
| 6 | 1000 / 200 / 32768 | explicit AVX2 encode | unchanged; AVX2/AVX2 |
| 7 | 4096 / 512 / 4096 | ordinary encode | unchanged; AVX2/AVX2 |

- Local `work`, Threadripper 9980X model 08h, kernel 6.8.0-137-generic.
  Child CPU 26, sibling 90, controller 0. Hold the canonical build/test/timing
  lock and CPU-pair lease. No other process's affinity or host setting changes.
- Slipgate's three services must remain inactive/disabled and its two
  containers stopped with restart `no`, checked before and after.
- Run 24 untimed main/off/on identity checks, then a ten-second passive
  zero-sibling gate. Every timed child's sibling non-idle delta must be zero.
- In cell order, run three rounds each. Each round has off/on/on/off and
  same-on A/B/B/A (identical executable, inode and mode argument). Cells 0/1
  additionally have Leopard1/on/on/Leopard1. Total: **216 timed children**.
- Each child uses aligned unique inputs, seed 20260906, one initial encode,
  four warmups and 21 samples of one full-output public encode. Setup,
  allocation, parity dumps and hashing are outside clocks. Leopard1 exposes
  its first R work buffers without an artificial final copy. Its native build
  is not ISA-matched to an explicitly restricted Leopard2 backend.
- Compute process medians, geometric ABBA ratios, then geometric means of
  three rounds. All eight same-on aggregate controls and all four unchanged
  off/on neighbors must be within `[1/1.02, 1.02]`. All four target off/on
  aggregates must be at least 1.05, with every target round above one.
- Decision precedence: invalid controls; rejected neighbor gate; rejected
  target gate; otherwise continue to production integration. No directional
  conclusion if controls fail. Leopard1 comparisons are separately identified
  and only apply to ordinary encode, not a nonexistent Leopard1 batch API.
- One 256-MiB/no-swap scope; codec CPU cap 30 seconds; all six memory events
  must be zero. Any process/provenance/isolation/resource failure invalidates
  the attempt. Partial data receives no ratios or performance conclusion.

A positive qualification is not default-on production integration, a
confidence interval, an independent-model review, another-host result, or
broad v19 closure. Codex self-review and deterministic/adversarial replay
replace Claude review under the user's explicit opt-out. No subagents or
SSH worker hosts are used. The unavailable unattended-research skill is not
replaced with a broad controller; this is one bounded local attempt.

Raw preparation: `/tmp/leopard-auto-boundary-screen.2TfJyu`. The existing
candidate correctness reference is the local readonly 218-entry bundle
`.research/leopard-79h/auto-gfni-boundary-checks.bJ6wiH`, outer manifest
`ae09b1b76a0dd877dcf0acff648df60c90cfdbe4442fbfbb0dee2f2bb71da512`.

## Qualification result

Preregistration `6bff9ecb8658c1313e88d680b47413d3f2531065` was pushed
before the sole launch. All 24 preflights and 216 timed children passed;
every timed sibling delta was zero. Sibling 90 stayed at 570178 non-idle
jiffies through the 10.000532241-second passive window. Shutdown snapshots
matched. Scope exit was zero after 89.87 seconds, peak 128,065,536 bytes /
256 MiB, all six memory events zero, swap zero.

Ratios are candidate throughput divided by comparator throughput:

| Workload | Candidate / boundary-off | Candidate / Leopard1 | Same-on control |
| --- | ---: | ---: | ---: |
| R200 / 32 KiB, ordinary | 1.536713 | 1.478279 | 0.996733 |
| R199 / 64 KiB, ordinary | 1.484029 | 1.428760 | 1.001166 |
| R200 / 32 KiB, one-item batch | 1.542486 | not compared | 1.000267 |
| R199 / 64 KiB, one-item batch | 1.460807 | not compared | 0.998436 |
| Existing R200 / 64 KiB | 1.000317 | not compared | 1.002194 |
| Inactive R199 / 32 KiB | 0.995878 | not compared | 1.002430 |
| Explicit AVX2 R200 / 32 KiB | 1.002416 | not compared | 0.999166 |
| K4096 / R512 / 4 KiB | 0.997726 | not compared | 0.999869 |

All four targets clear 5%, every target round is positive, and all eight
controls plus four unchanged neighbors pass the 2% equivalence bound.
Decision: **continue to production integration**. The code remains
default-off at this evidence checkpoint; qualification is not yet a claim
that the enabled default has passed its clean integration build.

The independent retained-only replayer imports no collector and executes no
codec. Normal, Python `-O`, and readonly-copy runs agree on all 54 round
ratios, 18 aggregates and the decision. It verifies the 14 frozen inputs,
seven build artifacts, 40 preparation records, 24 preflights, 216 raw timing
records, full Leopard1 parity, resource limits and shutdown state. Fifteen
malformed-row/claim mutations are rejected; four synthetic cases check the
positive gate and all three rejection decisions.

Evidence is local readonly `.research/leopard-79h/auto-gfni-boundary-screen.fvwFjq`:
671 manifest entries, 258,806,789 bytes including the manifest; outer SHA-256
`bf66d846bf2907a222d5ee87fcfba60abe3b7fc99966498c26982564af296268`.
Retention peak was 36,311,040 bytes / 256 MiB and readonly replay peak was
90,300,416 bytes / 256 MiB, both with all events/swap zero. The exact journal,
scope and pin hashes and unrounded results are in
`results/auto_gfni_boundary_screen_20260909.json`.

The attempt is consumed. Next integrate only the qualified default routes,
update existing default assumptions/docs, and perform clean Release plus
focused sanitizer checks. The broader goal and `.17.1` remain open.

## Production integration

The two qualified boundaries are now **default-on**. A clean GCC 13.3
Release build with both fields and the production/test-hook archives passed.
The entire production archive differs from the qualified default-off archive
by exactly one byte, at offset 285580: `2` becomes `1`. Independent ELF
inspection identifies it as the boundary-mode initializer in `.data`, not
executable code. All other archive bytes, including all kernel code, are
identical. The archive is
`d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a`.
This preserves the measured same-binary on-state; no new timing is asserted.

Final checks passed:

- 60 focused native checks, 30 Release and 30 both-field ASan/UBSan/LSan,
  including all three default routes, new boundaries, exclusions, guarded
  tails/partials, API paths, fallback returns and shared-codec concurrency.
- 13 CTests: hookless production parity for all three cells, the focused
  existing backend suite, existing inert/fallback cases, and six new real
  backend KAT/FF8-allocation/FF16-allocation failure cases. Each new case
  checks fault consumption, successful AVX2 fallback and cached-failure reuse.
- The separate production portable-ISA audit. Final native scope peaks were
  175,161,344 bytes for the 60 checks, 230,985,728 for the 13 CTests, and
  70,586,368 for ISA, each below 256 MiB with all six events and swap zero.

One initial combined CTest scope hit the 256-MiB cap and is preserved as a
failure, not a passing resource result. The new focused entry initially
created unrelated SSSE3/AVX-512 contexts before branching; those table
initializations are now avoided. Comparator parity buffers also leave scope
before batch/concurrency/decode fixtures. The same checks pass without
raising the cap or reducing their coverage. The broader monolithic tests
were not run. The final comment-only source refresh reproduced exactly the
same production archive and all three CMake test executable hashes.

The retained-only production verifier independently parses archive members
and the ELF data symbol, checks all 60 raw records and the CTest XML, rejects
skipped integration cases, and verifies the initial failure plus final
resource envelopes. Normal, Python `-O`, and readonly-copy runs pass; the
sealed copy's manifest and final source snapshots match the worktree.
This remains Codex self-review and deterministic evidence under the user's
Claude opt-out, not independent-model `CONVERGED`.

Production evidence is local readonly
`.research/leopard-79h/auto-gfni-boundary-production.B5NV9R`: 1179 entries,
85,466,116 bytes including the manifest; outer SHA-256
`b6c584e2d705e9104d306b57d8cb9fc15bd24b1da7ed3c8655d8a9cdec90d16c`.
The sanitizer checkpoint changes only the core object; its other 23 objects
remain byte-identical to the previously qualified both-field archive. It
predates only the final two-line comment correction, not a behavior change.

This completes `.17.1` and its parent `.17`, not the overall performance
goal. Next is `.38.5.4.18`: attribute and reduce the remaining explicitly
requested AVX2 deficit. Explicit AVX2 is still explicit AVX2; routing it to
GFNI would not solve that task. Broader CPU/size coverage and v19 exact-main
release qualification remain separate, open work.
