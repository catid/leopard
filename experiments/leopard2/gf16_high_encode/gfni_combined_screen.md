# Combined GFNI four-mode timing screen

Bead: `leopard-79h.38.5.4.16`. Date: 2026-09-07.
Status: sole attempt exhausted at the passive gate; **zero timed encodes**.

Correctness milestone `b46721b` validates the combination without changing
production. This new screen measures neither fusion, first-stage only,
terminal only, and both, as modes 0/1/2/3. Historical `.13`/`.14` ratios are
not combined or reused. The exhausted `.10` Leopard1 comparison is not retried.

## Measurement contract

All modes use the same executable path, inode, code layout and full dual-field
archive. The front end delegates the unchanged public workload: one initial
encode, four additional warmups and 21 measured encodes. It has no libc,
public-encode, per-callback or linker wrappers. The bounded once-per-pass
experimental hook remains in every mode; trace validation and output occur
outside timing. Mode 0 is an experimental OFF control, not pristine production.

The fixed host is foureyes CPU22/sibling86, controller CPU0. Six fixed cells
and three rounds each use the following eight mode labels, followed by an
identically labeled/positioned eight-slot control that always executes mode 0:

| Round | Mirrored mode order |
| --- | --- |
| 0 | 0, 1, 2, 3, 3, 2, 1, 0 |
| 1 | 1, 2, 3, 0, 0, 3, 2, 1 |
| 2 | 2, 3, 0, 1, 1, 0, 3, 2 |

There are 24 untimed server checks and 288 timed invocations. Each mode has
two observations per block with equal mean slot position. The complete
sequence is cell, round, factorial block, matched all-OFF block.

The protocol takes each process's median of 21 samples, then the geometric
mean of the two mirrored observations for each mode in each round. Ratios
are mode-0 time divided by mode-1/2/3 time, aggregated geometrically across
three rounds. **Mode 0 / mode 3 is the primary comparison.**

The interaction factor is `t1*t2/(t0*t3)` within each fresh round, aggregated
geometrically. Above one means the observed combined benefit exceeds the
product of the two individual benefits in that same round. This is diagnostic
attribution, not a selection gate, historical-ratio multiplication or a
confidence interval.

All 18 matched all-OFF aggregate contrasts and all 15 mechanically unchanged
neighbor contrasts must remain within `[1/1.02, 1.02]`. With those 33 controls
valid, primary gain at least 5% and positive gain in all three rounds select
future qualification. Otherwise the result is rejection or inconclusive
controls. Individual control rounds have no separate equivalence requirement.

The sole attempt is consumed at collector launch, including early failure.
No retry, CPU substitution, partial analysis, pooling or threshold change.
The canonical lock and physical-pair lease cover the attempt. All sixteen
plan/artifact pins and executable identity are checked before and after every
child. A ten-second passive sibling gate immediately precedes timing; every
timed child must observe zero sibling work. No unrelated affinity or workload
changes are authorized or performed.

## Pre-timing qualification

- All 48 Release/sanitizer fixed records match the original six workloads;
  24 full parity files match independently linked exact Leopard1, totaling
  244,056,064 bytes.
- Eight clock-free target exercises execute 26 independent checks each and
  retain exactly 52 passes. These test capacity, not timed-loop allocation or
  warmup behavior. Both 16- and 64-pass bounds pass all four modes, exact
  records, overflow refusal, reset and neighbor checks in both builds.
- The default-16 Release callback diagnostic is byte-identical to `b46721b`
  (`46c9b86c...`). Sixteen malformed/invalid-timing requests refuse before
  entering the workload. No timed request was executed during qualification.
- Eleven pure protocol tests include an independent log-space derivation,
  interaction attribution, every aggregate control mode, incomplete/order/
  isolation failures, typed metadata and fixed thresholds. Normal and Python
  `-O` runs agree. Retained-only front-end replay agrees in both modes.
- The 61 native qualification scopes peak at 141,209,600 bytes under 256 MiB.
  Front-end builds peak at 114,200,576 bytes under 512 MiB. All six memory-event
  counters and swap are zero. Two pre-compiler scratch setup failures are
  retained; the copied directory needed write permission. No codec or kernel
  rebuild, compiler-flag workaround or production change occurred.

The full Release archive remains
`aa3621ca58da46f22bba976f89f258ddc1317f72e74cc2a90602e909856d2c8c`.
The timing executable is
`959cb7b3bbce599e9d489e5b8b1d3e25dcf1185775f5a56289733084a59392e1`.
The plan is `029b62d7370664f56f950846b6061302d6e6f8cb46127d5f5fc72ee43f0b6f81`;
its frozen inventory is `de90ed1c119b0fd58613aa61930ce374b961e76044e6990f365989c9c964add8`.
The independent result replayer pins both and imports only the standard
library, with no collector import or codec execution.

Initial kernel/public correctness and Release ISA evidence remain in
`gfni-combined.7FNBT8`; the separate copy-resource clarification remains in
`gfni-combined-retention.Vl85Q7`. Both are retained locally and on ripper.
Review is Codex self-review plus deterministic/adversarial/sanitizer checks
under the Claude opt-out, not independent-model `CONVERGED`.
Production integration and broader v19/exact-Leopard1 qualification remain open.

## Completed attempt, no performance conclusion

Preregistration `14c12fa` was pushed before launch. All 24 untimed server
checks passed with the expected mode traces. The passive sibling counter
then increased from 194551 to 194560 over 10.000070824 seconds. The collector
stopped as specified, with `complete:false`, no analysis and zero timed
invocations. Neither the combined gain nor its interaction has been measured.
This does not reject the candidate or establish a Leopard1 comparison.

Independent standard-library replay in normal and optimized Python verifies
all sixteen frozen pins, raw preflight records, exact failure and absence of
timings/analysis. The server scope peaked at 132,927,488 bytes under 256 MiB;
result replay peaked at 12,242,944 bytes. All six memory-event counters and
swap are zero. The sole attempt is consumed; no retries or host changes were made.

- Attempt journal: `571fcbeb1a29da2bfcdec57e52de454670fbda8e8cea93bff7c59fb318b3d9d9`.
- Server scope: `416be8bebcc1c8881f806f43a09e1d94eaa8a8094d9af25b6f24176332f6d9a3`.

The evidence bundle is `.research/leopard-79h/gfni-combined-screen-failed.ZYQvWW`.
The read-only snapshot contains complete byte copies; a hard-link attempt
was refused across the temporary/workspace filesystem boundary and retained
as a setup failure. The second-host copy is separately hash-verified.
Build setup failures, qualification, frozen artifacts, raw server output and
independent replay are all retained. The original correctness bundles remain.

Post-copy clarification: copying the evidence reached the separate 256-MiB
artifact-copy cap, recording 2,832 `memory.max` events with zero OOM/kill/swap.
The copied bytes and manifest verified; this is not an all-zero resource
result. It does not alter the native build/check/attempt results above. The
bundled report predates this clarification; the copy log is retained separately.

Follow-up `leopard-79h.38.5.4.16.1` requires explicit user approval and a new
preregistration for any controlled-core successor. The user has been asked
whether unrelated user-space threads may temporarily be excluded from the
fixed CPU pair and restored afterward; no approval is assumed. The combined
task and overall performance objective remain open.
