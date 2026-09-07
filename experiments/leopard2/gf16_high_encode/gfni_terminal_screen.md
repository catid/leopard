# GFNI terminal accumulation timing screen

Bead: `leopard-79h.38.5.4.14`. Date: 2026-09-07.
Preregistration: one terminal-only OFF/ON screen; no production promotion.

This contrasts the default-off field overlay with the qualified terminal
kernel from correctness milestone `726a202`. First-stage forwarding remains
unchanged/OFF. Source staging, tiling, forward arithmetic and all other
backends remain unchanged. It is not a retry of `.10` or `.13` and cannot
establish a pristine-production or Leopard1 speedup by itself.

## Measurement boundary and fixed controls

`gfni_terminal_timing.cpp` delegates the existing public measurement loop:
one initial encode, four warmups and 21 measured encodes. There are no
individual-callback, public-encode, libc or linker wrappers. The bounded
once-per-source-policy-pass hook remains in both modes; trace validation and
printing occur outside timing. OFF/ON share one path, inode, code layout and
archive. Only `--terminal=0` versus `--terminal=1` differs. OFF itself contains
the experimental branch and hook; it is not unmodified production.

The plan fixes foureyes CPU22/sibling86, controller CPU0, ten seconds of
passive sibling observation, one attempt, six cells and three ABBA rounds for
both OFF/ON and same-OFF comparisons: 144 timed invocations. All six same-OFF
and five mechanically unchanged OFF/ON aggregate controls must be within
`[1/1.02,1.02]`. Only with valid controls does a target ratio at least 1.05,
positive in all three rounds, select future qualification. Otherwise the
decision is rejection or inconclusive controls. No per-control-round
equivalence requirement or confidence interval is implied.

Launching the collector consumes the sole attempt, even on an early failure.
No retries, CPU substitution, partial analysis, pooling or threshold changes.
The canonical lock and physical-pair lease cover the attempt. All sixteen
pinned plan/artifact inputs are rehashed before and after each child;
executable device/inode/size/mode must remain stable. Every timed child must
observe zero nonidle sibling jiffies. No unrelated workload or affinity changes.

## Pre-timing qualification

- All 24 Release/sanitizer OFF/ON records match the six original workloads.
  Twelve full parity files match exact independently linked Leopard1,
  totaling 122,028,032 bytes.
- Four clock-free target exercises run 26 separate check workloads each and
  retain exactly 52 passes. They test trace capacity, not timed-loop allocation
  behavior. Both 16- and 64-pass capacities pass both modes, exact records,
  overflow refusal, reset and neighbor tests in both builds.
- The rebuilt default-16 Release callback diagnostic is byte-identical to
  `726a202` (`9219486c...`). Sixteen malformed/invalid-timing requests refuse
  before workload entry. No timed request was run during qualification.
- Seven pure collector protocol tests and retained-only front-end replay pass
  normally and with Python `-O`. Post-build hashes remain unchanged; symbol
  inspection confirms no callback/linker wrappers in the timing executable.
- The 33 native qualification scopes peak at 141,201,408 bytes under 256 MiB;
  front-end builds peak at 115,245,056 under 512 MiB. All six memory-event
  counters and swap are zero; builds/checks are serialized.

The preceding sanitizer-archive scan was resolved by verified **existing**
project policy, not by weakening any check. Sanitizer correctness and
unsanitized Release ISA audits are separate; see
`gfni_terminal_sanitizer_policy.md`. Original raw failures remain retained.
No codec/archive rebuild or sanitizer/compiler flag change occurred here.

The codec archive is
`bf189e85b309eda5dea55c8a47e5b78e38779822d3556f86f85c457e081be127`.
The timing executable is
`683d3cef87436cae6236ba83227dc1866bb9944d19d7598f02cd34de81b1788f`.
The plan is `18ab975e9d7d00ca66565c2a73bd5055edf9254220120fff2a5be572b87661ae`;
its frozen inventory is `e77e07384ea55259428c76d15e8fc2e8425a540acc3963fa458a651d7850f7ae`.
The separate stdlib result replayer pins both, imports no collector and
executes no codec. Full source/correctness provenance remains in
`gfni-terminal.nEDGNa`, outer manifest
`c669fff28621c67577535f2924cff758ece2708dc2e217279a6c1a0cb03e3e96`.

Review is Codex self-review plus deterministic/adversarial/sanitizer checks,
not independent-model `CONVERGED`. The parent objective and actual broader
Release metadata/v19/exact-Leopard1 qualifications remain open.
