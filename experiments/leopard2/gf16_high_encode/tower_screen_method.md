# Direct tower encoder performance screen

Tracker `leopard-79h.38.5.4.18.4.3`. This is a new, local, single-attempt
warm-throughput experiment. No performance clocks have run under this method.
The protocol, tools and exact artifact pins must be committed **and pushed**
before its attempt starts. Correctness qualification is `6917d49` (encoder)
and `7d17cda` (public frontend); neither is a speed claim.

## Question and comparators

Does the default-OFF internal-tower encoder improve the actual public encode
operation by at least 10% versus original Leopard2 and native Leopard1 at the
remaining slow AVX2 cases, without disturbing excluded paths?

The three executable identities come from the read-only qualification bundle
`tower-public-qualified.KPkPef`. The freezer makes private executable/archive
copies into `/tmp/leopard-tower-screen.0mjTlC/frozen`, never links to another
lane's executable. SHA-256, source commits, shared driver objects, original
link recipes, qualification records and manifest are pinned and checked before
and after every invocation. No codec or driver rebuild is needed.

Original L2 and tower Release consume the identical compiled public driver.
Native L1 uses its native API-specific driver and unchanged native archive at
`6e5725eb`; it is not replaced with an AVX2-restricted baseline. These are
product comparisons, not an ISA-matched causal ablation.

Tower OFF and ON use one executable. OFF and excluded paths still pass through
the tower overlay's CopySource wrapper: **OFF is not pristine production**.
Original/ON and native/ON must be measured directly. Original/OFF is a separate
diagnostic, not a gain that can be added to OFF/ON. Previous scheduling and
GFNI attempts remain consumed and are never reused, pooled, trimmed or retried.

## Fixed workloads and roles

| Cell | K / R / bytes | L2 public request | Tower ON | Role |
| --- | --- | --- | --- | --- |
| 0 | 1000 / 200 / 32768 | AVX2 ordinary | selected | target |
| 1 | 1000 / 199 / 65536 | AVX2 ordinary | selected, two tiles | target |
| 2 | 1000 / 200 / 65536 | AVX2 ordinary | selected, two tiles | target |
| 3 | 4096 / 512 / 4096 | AVX2 ordinary | excluded | unchanged neighbor |
| 4 | 1000 / 199 / 32768 | AUTO, actual AVX2 | selected | target |
| 5 | 1000 / 200 / 32768 | explicit GFNI | excluded | unchanged neighbor |
| 6 | 1000 / 200 / 32768 | AUTO, actual GFNI | excluded | unchanged neighbor |
| 7 | 17 / 7 / 64 | GF8 AVX2 | excluded | unchanged neighbor |
| 8 | same as cell0 | actual one-item L2 batch | selected | target |

Cell4 is a measured remaining AUTO deficit and is therefore a full target,
including a direct native comparison, not an affected neighbor. Cell3 is
excluded, not affected. There are no affected-only neighbors in this screen.
Excluded means no tower arithmetic; it does not assert unchanged wrapper/code
layout. Their original/ON comparison still has a symmetric 2% gate.
Native cell8 calls its ordinary encoder, not a nonexistent native batch API.

## Warm cost boundary

Each process makes four preflight calls, four warmup passes and 21 measured
passes with four slots each. The four preflight tower snapshots must have the
exact lifetime-initialization sequence: e.g. eligible 0110 has 0,1,1,1; excluded,
original and all-OFF calls stay at zero. All untraced work counters remain zero.

Each clock span encloses complete public encode calls, including source and
final output conversions, public result checks, per-call counters and loop
overhead. Selection/inspection, buffer setup, normalization/storage, guards and
full parity comparisons are outside. Parity is checked after **every group**;
those memory reads influence later cache state. No overhead subtraction.

The group size is 1 for GF16 and 256 for tiny GF8 cell7, fixed from qualified
workloads without new timing-based tuning. Each process retains all 84 spans,
with 104 public calls for group1 or 25,604 for group256. Grouped averages are
not single-call latency samples. Batching remains repeated one-item calls.

The additional **6 MiB cache** remains a real cost; canonical tables stay.
Its cold initialization occurs in preflight, outside warm spans. A pass does
not establish cold-start latency, memory savings or context-creation speed.
Any cold-start experiment needs a separately specified method.

## Fixed order, estimates and gates

Cells run 0 through 8, each with three consecutive rounds. In each round:

| Comparison | Process schedules | Cells |
| --- | --- | --- |
| Paired forward; reverse | 0110; 1001 | all |
| Same OFF | 0000 four times | all |
| Same ON | 1111 four times | all |
| Original / ON | PPPP,1111,1111,PPPP | all |
| Original / OFF diagnostic | PPPP,0000,0000,PPPP | all |
| Same original | PPPP four times | all |
| Native / ON | NNNN,1111,1111,NNNN | targets 0,1,2,4,8 |
| Same native | NNNN four times | targets 0,1,2,4,8 |

This is **714 timed processes / 59,976 spans**, preceded by 36 clock-free
preflights (native/original/OFF/ON at each cell). Mixed-process estimates are
the median of 21 per-pass sqrt(OFF-product/ON-product) ratios. Keep the two
paired directions separate. Homogeneous process cost is the median of all 84
normalized spans, averaging the two middle values. Four-process ABBA uses
sqrt(cost0*cost3/(cost1*cost2)). Each same-path process also supplies a median
of 21 internal outer/middle ratios. Each distinct comparison aggregates its
three rounds geometrically, retaining every individual round.

Fixed decision precedence:

1. All **32 cross-process and 128 within-process** same-path aggregates must
   lie in [1/1.02,1.02], otherwise `inconclusive_controls`.
2. Every excluded neighbor 3,5,6,7 must lie in that symmetric interval for
   each paired direction and direct original/ON, otherwise reject neighbor.
3. Every target 0,1,2,4,8 must reach **1.10** for each paired direction.
4. Every target must reach **1.10** for direct original/ON.
5. Every target must reach **1.10** for direct native/ON.

For all target gates every corresponding individual round must also exceed 1.
All 201 aggregates / 603 round ratios remain in the result, including outlying
individual control rounds. A passing screen permits production-candidate
qualification only, not automatic promotion, merging, deployment or release.
No confidence intervals, identified noise cause or broad host claim is made.
Grouped GF8 has failed a prior cross-process control; grouping here is not
claimed to solve that problem, and its unchanged gate remains mandatory.

## Isolation, resources and review

Only local `work`, Threadripper9980X, CPU26/sibling90/controller0. Hold the
canonical `/tmp/leopard-gf8-authoritative.lock` and existing CPU-pair lease.
Verify Slipgate/OBS services disabled/inactive and the two containers stopped
with restart=no, before and after. Require one initial ten-second passive
zero-sibling window after preflights and before the first timing. Separately,
require zero sibling non-idle tick deltas around every timed process. Do not
move, stop or reconfigure unrelated processes or visit SSH worker hosts.

Native preparation and timing are serialized under 256 MiB/no swap. Each child
has 30-second CPU, 60-second wall, zero core-dump and 1 MiB output-file limits.
All six memory-event counters and swap must be zero; no resource cap is raised.
Read-only freezes, replays and pure tests hold the same canonical lock.

Exclusive creation of `/tmp/leopard-tower-screen.0mjTlC/attempt1` consumes the
sole attempt, including failure in preflight/isolation. There is no resume,
rerun, CPU search, adaptive performance stopping, trimming or pooling. Partial
attempts get no performance analysis. Successful timing requires terminal
exit0 and valid post-run resource/identity evidence, not just a passing ratio.

Review uses Codex self-review, the user-authorized local read-only bug reviewer
and deterministic/adversarial tests under the explicit Claude opt-out. It is
not an independent-model fixed-point CONVERGED claim. The independent replayer
imports no collector or qualification validator and executes no codec.

Before publication, test every gate, control, target role, record field,
initialization snapshot, zero diagnostic count, exact threshold, raw terminal
inventory and resource failure; compare collector and independent replay,
including a complete maximum-duration synthetic record under the existing
4 MiB JSON bound. Then freeze actual binaries and validate all 36 untimed
checks through both readers. This document is not evidence those checks passed
until their actual results are recorded.

### Pre-freeze review results

All 21 pure tests pass in normal and optimized Python. This includes every
target's four gain comparisons, exact 1.10 and adjacent integer cases, rejection
of the old 1.05 gate, all 160 stability controls, all four excluded neighbors,
and adversarial lifetime-cache snapshots/totals for all nine cells, six orders
and both check/measure record forms. The full 714-process maximum-interval
fixture fits the unchanged 4 MiB parser bound; raw-file replay agrees with the
collector. Normal/optimized scope peaks are 56,061,952 / 55,730,176 bytes under
256 MiB, exit0 with all six memory events and swap zero.

The local read-only reviewer found an isolation wording ambiguity and missing
adversarial initialization-counter coverage. Both were corrected; a narrow
follow-up review confirmed resolution with no new defect. No native job or
clock was executed by the reviewer. The initial 20-test passing log is retained
separately from the expanded 21-test gates. Actual frozen-byte identity and the
36 native preflights will be recorded separately before preregistration is
published; neither this paragraph nor synthetic fixtures substitute for them.
