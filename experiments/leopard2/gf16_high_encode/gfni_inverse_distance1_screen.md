# GFNI first-inverse distance-one timing filter

Bead: `leopard-79h.38.5.4.13`. Date: 2026-09-06.
Status: preregistered candidate filter, not a production or Leopard1 result.

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
