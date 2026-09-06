# GFNI first-inverse distance-one timing filter

Bead: `leopard-79h.38.5.4.13`. Date: 2026-09-06.
Status: completed clean screen; below the preregistered 5% threshold.

Correctness milestone `c393c5c` proves that forwarding the target's first
inverse groups replaces 2,000 split two-way calls with 500 existing in-place
GFNI range calls. This screen asks whether that exact change improves public
encode time. It does not change source staging, tiling, final accumulation,
forward arithmetic or backend selection. It is separate from the rejected
source-staging experiment and the exhausted `.10` Leopard1 comparison.

## Timing boundary and controls

The new executable delegates the unchanged `current_route_screen.cpp` workload
and its initial one encode, four additional warmups and 21 timed public
encodes. It has no linker, libc, public-encode or individual-backend-callback
wrappers. A bounded once-per-source-policy-pass hook remains in the frozen
experimental FF16 object. Both modes retain that hook and its pass trace;
validation and printing happen after the measurement loop. This is a
diagnostic with small common instrumentation, not zero-instrumentation
production performance.

OFF and ON use the same executable path, inode, code layout and archive. Only
`--fuse=0` versus `--fuse=1` changes. OFF itself contains the experimental
field branch and hook, so its time must not be described as pristine current
production time. A successful screen only selects future production-code
qualification and a separate independently linked Leopard1 comparison.

The committed plan fixes foureyes CPU22/sibling86, controller CPU0, a ten-second
passive gate, one attempt, six cells, three OFF/ON ABBA rounds and three
same-OFF ABBA rounds per cell: 144 timed invocations. All six same-OFF and five
unchanged OFF/ON aggregate controls must lie in `[1/1.02,1.02]`. With valid
controls, the target needs at least 5% aggregate gain and positive gain in all
three rounds. Otherwise the result is rejection or inconclusive controls as
specified in the plan. No per-control-round equivalence gate is implied.

Starting the collector consumes the sole attempt, including an early failure.
No retry, CPU substitution, partial inference, pooling or threshold relaxation.
No unrelated workload or affinity changes. The canonical campaign lock and
CPU-pair lease cover the attempt; frozen inputs are rehashed before and after
every child and the executable's device/inode/size/mode must remain stable.
Every timed child must observe zero nonidle sibling jiffies.

## Pre-timing validation

- All 24 Release/sanitizer × OFF/ON × six-cell records match the existing
  public workload. Twelve full parity files match standalone exact Leopard1,
  totaling 122,028,032 bytes.
- Four clock-free target exercises each run 26 check workloads and retain
  exactly 52 source-policy passes. They test capacity, not timed-loop
  allocation behavior. Both 16- and 64-pass capacities pass selection, exact
  records, atomic overflow refusal, reset and neighbor tests in both builds.
- The recompiled default-16 Release callback diagnostic is byte-identical to
  the prior frozen binary (`64d0767a...`). Sixteen malformed/invalid-timing
  requests refuse before entering the workload. Seven pure collector tests
  and retained-only qualification replay pass normally and with Python `-O`.
- Build peak is 114,929,664 bytes under 512 MiB; the 33 successful native
  validation scopes peak at 141,471,744 bytes under 256 MiB. All six memory
  event counters and swap are zero. No field or sanitizer check was removed.

The experimental archive remains `39cfa6d4d5b06aae3837f036254b94db7730685271eeae82f68858265a26881f`.
The timing executable is `3e45d325a2291759e49fcb02fc301e5444e6e4e9bdaf171efe5bbcbf0d0cabd9`.
The preregistered plan is `caaac7545c1ce8655689ee965702fdc3a8c8a640dedf2ea736ce7c10432f0389`;
its frozen inventory is `b9faeccb48e9261f13829766a025605820363798e8c5bb93b14aa0506f11e54f`.
The separate stdlib result replayer pins both, imports no collector and
executes no codec.

Review is Codex self-review and deterministic/adversarial/sanitizer checks
under the user's Claude opt-out, not an independent-model `CONVERGED` claim.
The parent performance gap and existing v19 requirements remain open.

## Completed result

Preregistration `b6a22c1` was pushed before the sole server launch. All twelve
untimed checks passed. The sibling counter stayed at 192335 over
10.000071536 seconds, and all 144 timed invocations observed zero sibling
work. Every frozen input and executable identity remained unchanged, and all
workload/trace checks passed. The decision is `reject_for_this_screen`.

Ratios are overlay-OFF time / fused-ON time; larger than one favors fusion.
They are geometric means of three round contrasts, not confidence intervals.

| Cell / route | OFF / fused | Same-OFF control |
| --- | ---: | ---: |
| K1000/R200/64 KiB AUTO/GFNI target | 1.033262 | 1.003129 |
| Same shape, explicit AVX2 | 0.998708 | 1.000139 |
| Same shape, explicit AVX-512 | 1.004082 | 0.997000 |
| K1000/R200/32 KiB AUTO | 0.992013 | 0.998582 |
| K1000/R199/64 KiB AUTO | 1.000750 | 0.998494 |
| K4096/R512/4 KiB AUTO | 1.008114 | 1.000048 |

All eleven aggregate controls meet the fixed 2% equivalence band. The rule
does not require each individual control round to meet that band: the 32-KiB
neighbor's first OFF/ON round is 0.968212 and the K4096 second round is
1.026739. The target's three ratios are 1.039564, 1.018694 and 1.041686.
All are positive, but the 1.033262 aggregate falls short of the required
1.05. This is a directional improvement in the experimental contrast, not
evidence of no benefit, a pristine-production speedup, or a Leopard1 win.
The candidate remains experiment-only; no policy is promoted.

The independent stdlib replay rehashed all thirteen plan/artifact inputs,
checked every raw check/trace and ordered timed record, and independently
recomputed all round and aggregate contrasts using log medians. Normal and
optimized Python agree with the collector's decision. The server scope
exited zero and peaked at 133,742,592 bytes under 256 MiB; the larger replay
peak is 11,964,416 bytes. All six memory-event counters and swap are zero.
Production source/header/CMake differences from base `36dc0c8` remain empty.

The complete read-only bundle is retained locally and on ripper at
`.research/leopard-79h/gfni-inverse-distance1-screen.MSr1fM`, including the
new front-end source/binaries, native qualification, every raw server output,
frozen inputs, independent replay and resource logs. The preceding complete
correctness/source bundle `gfni-inverse-distance1.EN8c4J` remains separately
retained with outer hash `7bdbe9d17c8d9f39436fd5ee049ce8dac40d96e8f17543b2eca643154763437b`.

- Attempt journal: `50e5b846a3bf6dac24631cd8fad268d482f8cca6c5f095324d20dfc6352a373e`.
- Server scope log: `1a773929a9612b4bf6eb5d84749bad3491ebfe81980cdc000318eddab364c4a3`.

The `.13` evaluation is closed as a completed below-threshold experiment,
not as completion of the performance objective. Follow-up `.38.5.4.14`
investigates a GFNI-only final inverse accumulation kernel: the current split
boundary still materializes its first layer before its accumulating second
layer. That is distinct from the existing in-place range kernel and from the
old rejected AVX-512 nibble-table accumulation experiment `.5`. The GFNI
version is not assumed faster; it needs its own correctness and isolated
timing evidence. Any later combination with first-stage forwarding needs a
new contrast, not multiplication of historical ratios or relaxation of gates.
