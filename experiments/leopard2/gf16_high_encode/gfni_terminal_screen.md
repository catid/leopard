# GFNI terminal accumulation timing screen

Bead: `leopard-79h.38.5.4.14`. Date: 2026-09-07.
Status: completed clean screen, below the fixed 5% threshold; no promotion.

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

## Completed result

Preregistration `ef963c6` was pushed before the sole server launch. All twelve
untimed checks passed. The passive sibling counter stayed at 194060 over
10.000063972 seconds, and all 144 timed invocations observed zero sibling
work. Frozen hashes and executable identity remained unchanged throughout.
The decision is `reject_for_this_screen`.

Ratios are overlay-OFF time / terminal-ON time; greater than one favors fusion.
They are geometric means of three round contrasts, not confidence intervals.

| Cell / route | OFF / terminal | Same-OFF control |
| --- | ---: | ---: |
| K1000/R200/64 KiB AUTO/GFNI target | 1.042655 | 1.003660 |
| Same shape, explicit AVX2 | 1.001940 | 1.007697 |
| Same shape, explicit AVX-512 | 0.997797 | 0.999345 |
| K1000/R200/32 KiB AUTO | 0.999727 | 0.999303 |
| K1000/R199/64 KiB AUTO | 0.993877 | 1.007691 |
| K4096/R512/4 KiB AUTO | 0.999601 | 1.002030 |

All eleven aggregate controls meet the fixed 2% band. Individual control
rounds need not: for example, the target's first same-OFF round is 1.026034
and the R199 neighbor's second OFF/ON round is 0.977368. Their aggregate
controls pass as preregistered. The target's three ratios are 1.047415,
1.039629 and 1.040937. The 1.042654726 aggregate is directionally positive
in all three rounds but below 1.05, so no production policy is promoted.
This is not evidence of no benefit, nor a pristine-production or Leopard1 win.

The independent stdlib replay rehashed all sixteen inputs, checked every raw
record and trace in exact order, and independently recomputed log-median
round/aggregate contrasts. Normal and optimized Python agree with the
collector. The server scope exited zero and peaked at 133,791,744 bytes under
256 MiB; the larger initial result replay peak is 12,275,712 bytes. All six
memory-event counters and swap are zero. Production code remains unchanged.

The complete read-only bundle `.research/leopard-79h/gfni-terminal-screen.33r6nh`
is retained locally and on ripper, including front-end qualification, frozen
inputs, all raw server output, independent replay and resource logs. The
preceding `gfni-terminal.nEDGNa` correctness bundle and
`sanitizer-isa-policy.DVv0Pp` correction bundle remain separately retained.

- Attempt journal: `b44dc82dee0abb25f588b797586175aa353a5a62e0f5b4b935140368640d5f00`.
- Server scope: `c0f5eaa5f55615bb027c59b630e3e8602e525a43be30fb411a49243ef0668ecf`.

The terminal-only `.14` evaluation closes as a completed below-threshold
experiment, not completion of the performance objective. Follow-up `.16`
evaluates first-stage and terminal fusion together using a new four-state
contrast. The separate 1.033262x and 1.042655x observations motivate that
experiment but must not be multiplied, added, pooled or used to bypass fresh
correctness, interaction and timing evidence. The original gates stay fixed.
