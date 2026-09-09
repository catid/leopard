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
