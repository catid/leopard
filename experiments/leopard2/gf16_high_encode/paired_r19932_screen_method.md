# Paired AUTO R199/32 KiB successor: method review and preregistration

Bead: `leopard-79h.38.5.4.19.1.2`. Review provenance: Codex self-review and
deterministic/adversarial tests. The user opted out of Claude and remote workers;
this is not independent-model `CONVERGED`.

The previous `aa2b034` attempt is consumed and inconclusive. Its native,
explicit-AVX2 and tiny-GF8 stability failures remain recorded; the cause of the
millisecond whole-process shifts is unknown. Do not repeat that attempt, pool
its samples with this experiment, omit outliers, or relax its gates. Candidate
`45e2eff` remains default-OFF. This experiment does not change codec source.

## What changes, and what does not

Use the exact plain steady-clock frontends qualified in `553d0df`, with original
native Leopard1 and both-field Leopard2 Release archives. Do not rebuild or
substitute pure-AVX2 Leopard1. Copy executable bytes into this lane's immutable
directory, pin their SHA-256 and source identities, and verify before and after
every invocation. Both OFF and ON use the same Leopard2 executable.

OFF/ON comparisons now alternate quiescent diagnostic state in a single process
using one codec and fixed buffers. Both orders `0110` and `1001` run separately.
Every process has four route-probed preflight calls, four warmup schedule passes,
and 21 sampled passes of four slots. The tiny GF8 cell groups **256** complete
public calls per span; every other cell groups **one**. These choices come from
untimed-qualified options, not new timing-based tuning. Every process therefore
has 84 spans/168 clock endpoints and 104 total calls (25,604 for GF8).

Inside each duration: start clock, N complete public API calls, end clock. Public
result checks, per-call counters, and repetition-loop overhead are included.
State selection/inspection, normalization, sample storage, slot accounting,
full parity comparison and allocation guards are outside. Parity comparisons
after **every group** touch output memory and influence subsequent cache state.
No overhead subtraction. Grouped per-call averages are not single-call latency.
The batch target still performs one-item batch calls, never a 256-item batch.
OFF includes diagnostic-path overhead; it is not pristine production OFF.

## Fixed inventory and order

Use all nine original cells, in numerical order, with three consecutive rounds
per cell. In each round run these processes in this exact order:

| Comparison | Four-slot schedule in each successive process | Cells |
| --- | --- | --- |
| paired forward | `0110` | all nine |
| paired reverse | `1001` | all nine |
| same OFF | `0000`, `0000`, `0000`, `0000` | all nine |
| same ON | `1111`, `1111`, `1111`, `1111` | all nine |
| native / ON | `NNNN`, `1111`, `1111`, `NNNN` | targets 0 and 1 |
| same native | `NNNN`, `NNNN`, `NNNN`, `NNNN` | targets 0 and 1 |

Thus **318 timed processes**, 26,712 measured spans, plus **27 untimed
preflights** (`NNNN`, `0000`, `1111` for every cell). No adaptive schedule,
stopping on observed speed, reruns, CPU search, or sample removal. An execution,
identity or isolation failure aborts the attempt and prevents analysis.

## Estimators, controls and decision gates

All 84 spans are retained and participate. For each paired process, compute
`sqrt(OFF-slot-product / ON-slot-product)` separately for each of 21 passes,
then take the median of those 21 ratios. Reverse order uses OFF slots 1 and 2.
Keep the two paired directions separate throughout analysis and gating.

For every homogeneous process, its scalar cost is the median of **all 84**
normalized span durations (arithmetic mean of the two middle values). Each
cross-process four-process comparison is `sqrt(cost0*cost3/(cost1*cost2))`.
The 20 original cross-process control aggregates remain: same OFF and same ON
for nine cells, plus same native for two targets. Within-process pairing is
not evidence that the separately linked native comparison is stable.

Additionally, for each of the four process slots in every same-path comparison,
calculate the median of 21 within-process outer/middle ABBA ratios. Keep each
process slot separate across rounds: **80 additional within-process control
aggregates**. These are derived only from same-path controls, not mixed-state
target processes; they do not replace or mask the 20 cross-process controls.

Each comparison/control aggregates its three round ratios by geometric mean.
Retain every individual round ratio, including failures; there is no sample
trimming, pooling across comparisons, or reuse of old experiment data. Gates:

- All 100 same-path control aggregates must lie in `[1/1.02, 1.02]`.
- All seven unchanged neighbors must lie in that same interval for **each**
  paired direction separately.
- Both targets must reach `OFF/ON >= 1.05` for each paired direction, with
  `OFF/ON > 1` in every round for each direction.
- Both targets must reach `native/ON >= 1.05`, with `native/ON > 1` every round.

Decision priority: invalid/partial data have no analysis; otherwise failed
controls mean `inconclusive_controls`, followed by neighbor, paired-target,
and native rejection gates, respectively. Only all-pass permits the next
**actual default-ON artifact qualification**, not automatic production promotion.
No confidence interval, noise-cause identification, pristine-OFF, or broad-v19
authoritative claim. The additional controls are deliberately conservative;
they can reject a promising candidate rather than establish an improvement.

## Resource, isolation and publication gates

Local `work` only, CPU26/sibling90/controller0, canonical
`/tmp/leopard-gf8-authoritative.lock` plus the pair lease. No host settings,
unrelated processes or affinity changes. Serial, one 256 MiB scope, swap zero,
30-second CPU limit/60-second wall limit per child. The wrapper records exit,
peak, all memory events and swap. Memory-limit events are never described as
all-zero. The terminal raw replay requires exit zero and all-zero events for
this timing attempt. A complete collector record alone cannot bypass that gate.

Slipgate/OBS/Forge stopped-and-disabled state must match before and after.
Require ten seconds of zero sibling non-idle ticks before timing and zero
sibling ticks around every measured invocation. No unrelated CPU displacement.

`paired_r19932_screen_plan.json` fixes the attempt root
`/tmp/leopard-paired-integration.la42Cz/attempt1`, one exclusive creation and no
resume. `freeze_paired_r19932_screen.py` only copies pinned qualified artifacts
and executes 27 four-call `--check` preflights, without benchmark clocks. Pure
tests and the independent raw replayer must pass normally and under `-O`.
Commit **and push** this method, plan, collector, replayer and tests before any
`--measure` call. The collector verifies their bytes against that commit and
that the commit is an ancestor of the fetched origin topic branch.

The fixed-plan and raw-record replayer has no collector import or native
execution. It independently reconstructs public identity, call counts,
sample-major ratios, per-process medians and every gate. Frozen inputs bind to
the 751-entry qualified timer manifest (`6c86e08b...`) and the previously
validated build/check/replay records. This is reuse of completed qualification,
not a claim to repeat its 1.03 GB full-parity validation during a timing run.

## Pre-clock review evidence

Sixteen pure tests pass normally and with `-O`. They exercise all 20 cross-
process and 80 within-process gates, each neighbor in both paired orders and
both regression directions, both targets/native thresholds and individual
round failures, exact 5% threshold behavior, sample-major pairing, median84,
typed identity/call/route/clock-source checks, fractional GF8 normalization,
malformed intervals, missing/reordered records, every fixed-plan field and
input pin, duplicate JSON keys, and altered resource/derived claims.

A full synthetic raw-file fixture exercises 318 process outputs and 27
preflights through the independent replayer, including terminal-record/raw-file
disagreements, incomplete data, sibling contamination, stderr, disabled-service
state and unexpected-file failures. Native execution is forbidden in these
tests; input provenance is checked separately against actual frozen files.
Both modes passed in one 256 MiB/no-swap scope: peak 42,913,792 bytes, exit zero,
all six memory events and swap zero. No real performance samples were used.
