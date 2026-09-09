# AVX2 forward/accumulating scheduling: direct performance method

Bead `leopard-79h.38.5.4.18.3.2`. Review is Codex self-review plus
deterministic/adversarial tests, under the user's Claude opt-out. It is not
independent-model `CONVERGED`. This document precedes all benchmark clocks for
this candidate; registration requires this method, plan and tools to be both
committed and pushed.

The prior inverse-only experiment and both inconclusive GFNI attempts remain
consumed. Their timings are not reused, added, pooled, trimmed or retried here.
The distinct forward/accumulating runtime candidate passed focused qualification
in `bd4c175` and public frontend qualification in `a6634c7`. Do not rebuild or
repeat those matrices. No production source or default changes in this screen.

## Comparators and cost boundary

Freeze actual private copies of the exact native Leopard1, original current
Leopard2, and runtime Leopard2 plain executables from the `a6634c7` evidence.
Retain their archives, build/check/replay records, qualification manifest, and
native/shared-Leopard2 driver objects. SHA-256 and source/commit identities are
fixed in `avx2_adjacent_screen_plan.json`. Verify all frozen inputs before and
after every invocation. Only the untraced Release links can supply timings.

OFF and ON run in the same runtime executable. Original Leopard2 and runtime
Leopard2 consume the same compiled public driver object, but are separate
processes with different codec archives and link-specific control shims.
**Runtime OFF is not original production**: the prepared-range loop changed
from 39 instructions/3 stack references to 40/4. Original/ON is a direct
comparison; do not infer it by multiplying OFF/ON with historical results.
Original/OFF is separately observed as diagnostic overhead, never a promotion
gate or substitute for Original/ON. Native Leopard1 remains the product target,
not the restricted-ISA variant. This is not an ISA-matched causal ablation.

Use the qualified sample-major schedule: four preflight calls (first may
initialize caches), four warmup passes, and 21 sampled passes, with four slots
per pass. Each span surrounds N complete public calls with two clock endpoints.
Public result checks, per-call counters and repetition-loop overhead are inside.
State selection/inspection, normalization/storage, slot bookkeeping, full parity
comparison and guards are outside. Full parity checks after **every group**
touch output memory and influence subsequent cache state. Initialization and
state-selection cost are excluded; this is not context-creation throughput.

Choose N=1 for all GF16 cells, and N=256 only for tiny explicit-AVX2 GF8 cell7.
This preserves the qualified power-of-two normalization and amortizes endpoint
cost for the small control. It is not a claimed cure for the previous GF8
cross-process shift: grouped-256 controls failed in that different experiment,
and the cause remains unknown. No choice is tuned using new timings. Retain
all 84 spans per process, no overhead subtraction. Grouped per-call averages
are not single-call latency samples. Totals are 104 public calls for group1,
25,604 for group256. Cell8 uses repeated **one-item** batch calls, never a
multi-item batch. Native's cell8 counterpart is its ordinary public encoder;
do not label it a native batch API.

## Fixed cells, order and estimators

Preserve the eight AVX2-screen cells exactly and add the qualified one-item
batch as cell8. Targets are 0/1/2/8; affected neighbors 3/4; unchanged neighbors
5/6/7. This retains the earlier AVX2 affected/unchanged distinction, not the
different GFNI experiment's role assignment. GF8 explicitly requests AVX2;
AUTO R19932 remains AVX2 with the unpromoted GFNI extension OFF.

Run cells 0 through 8 in order, three consecutive rounds per cell. Within each
round run these processes in the following fixed order:

| Comparison | Four-slot schedule in each successive process | Cells |
| --- | --- | --- |
| paired forward | `0110` | all |
| paired reverse | `1001` | all |
| same OFF | `0000`, `0000`, `0000`, `0000` | all |
| same ON | `1111`, `1111`, `1111`, `1111` | all |
| original / ON | `PPPP`, `1111`, `1111`, `PPPP` | all |
| original / OFF, diagnostic | `PPPP`, `0000`, `0000`, `PPPP` | all |
| same original | `PPPP`, `PPPP`, `PPPP`, `PPPP` | all |
| native / ON | `NNNN`, `1111`, `1111`, `NNNN` | targets 0/1/2/8 |
| same native | `NNNN`, `NNNN`, `NNNN`, `NNNN` | targets 0/1/2/8 |

There are **690 timed processes / 57,960 spans**, preceded by **36 untimed
preflights**: `NNNN`, `PPPP`, `0000`, `1111` for every cell. Every scheduled
process is retained. No adaptive stopping on performance, CPU search, resuming,
reruns or sample exclusion. Execution, identity or isolation failure ends the
single attempt without performance analysis.

For each mixed paired process, compute `sqrt(OFF-product / ON-product)` per
four-slot pass, then median across its 21 passes. Keep `0110` and `1001`
separate through all gates. For homogeneous processes, scalar cost is the
median of **all 84** normalized spans (mean of the two middle values). A
four-process ABBA ratio is `sqrt(cost0*cost3/(cost1*cost2))`.

For every process slot in each same-path control, also retain its median of 21
within-process outer/middle ABBA ratios. This yields **31 cross-process and
124 within-process control aggregates**, kept separate. Within-process pairing
does not demonstrate stability of separately linked original/native processes.
Each comparison aggregates its three round ratios by geometric mean. All
individual round ratios remain reported; they are not separately discarded.

## Decision gates and limitations

Gates apply in this fixed priority:

1. All 155 same-path control aggregates must lie in `[1/1.02, 1.02]`; otherwise
   `inconclusive_controls` takes priority over apparently passing gains.
2. Unchanged neighbors 5/6/7 must lie within the same interval versus OFF in
   **each** paired direction and versus original Leopard2 directly.
3. Affected neighbors 3/4 must be at least `1/1.02` on those same three
   comparisons. Improvements are allowed, as in the earlier AVX2 screen;
   this is not relaxation of a symmetric unchanged-path gate.
4. All four targets must reach OFF/ON >=1.05 in each paired direction, with
   every corresponding round >1.
5. All four targets must reach original/ON >=1.05, with every round >1.
6. All four targets must reach native/ON >=1.05, with every round >1.

Native's 5% gate strengthens the earlier inverse-only screen, whose native
comparison was descriptive. Original/ON adds the missing product-baseline
gate. A pass permits only **actual production-candidate qualification**, not
automatic promotion, merging or deployment. Any failed prerequisite bars using
passing subsets as promotion evidence. No confidence intervals, identified
noise cause, pure-production-OFF claim, broad host generalization or v19
authoritative-campaign claim is made.

## Resource, isolation and pre-clock review

Local host `work` only: CPU26/sibling90/controller0, canonical
`/tmp/leopard-gf8-authoritative.lock` and the existing CPU-pair lease. No
unrelated process, affinity, service, host-setting, or SSH-worker changes.
Require the stopped/disabled Slipgate/OBS/Forge condition before and after;
ten seconds of zero sibling non-idle ticks before timing and zero sibling
ticks around every measured invocation. Only this controller/child's own
affinity is set. No automatic repair or retry of a failed isolation condition.

Run serially in one 256 MiB, no-swap scope; each native child has 30-second CPU,
60-second wall and 1 MiB output-file limits. The terminal timing replay requires
exit0 and all-zero memory events/swap. The previous qualification retention's
6,738 memory-pressure events remain disclosed, not relabeled as all-zero or
carried into this timing gate. No builds are needed; any required build would
remain serial under 512 MiB with the prescribed GCC memory-scheduling flags.

The exact attempt root is `/tmp/leopard-adjacent-screen.04vrZl/attempt1`.
Its exclusive creation consumes the sole attempt even if a subsequent
preflight/isolation check fails; no resume. The freezer only copies qualified
inputs and runs the 36 clock-free four-call checks. Freeze, code review and
pure tests must precede a **committed and pushed** preregistration. The
collector checks source bytes against that commit and ancestry in the fetched
origin topic branch before any `--measure` call.

The independent replayer has no collector/qualification-validator imports or
native execution. It reconstructs fixed public identity, call and route/count
metadata, grouped normalization, all ratios and gates from raw output files,
and checks the qualified manifest and actual shared-driver link recipes.

Pre-clock tests cover every cross/within control, every neighbor comparison,
all four target comparisons and per-round positivity, original/OFF being
diagnostic only, every protocol/pin field, typed identities and zero diagnostic
counts, sample-major order, median84, raw-file/terminal-record disagreement,
missing/reordered outputs, changed conditions and resource counters.

The initial synthetic raw fixture exposed pretty-printed 690-process records
exceeding the unchanged 4 MiB parser bound. Compact JSON preserves **every**
sample and metadata field; a complete maximum-interval fixture verifies it
fits and parses in both readers. No reader/resource/performance limit is
raised. Initial failure logs/source are retained; no timing attempt existed.

All 20 pure tests pass normally and under Python `-O`, including a complete
690-process/36-preflight raw-file fixture and the maximum-interval size check.
The normal/optimized scopes peaked at 58,589,184 / 57,290,752 bytes under
256 MiB, each exit0, all six memory events0 and swap0. The initial oversized
fixture failure peaked at 61,612,032 bytes, exit1, events0/swap0; it remains
retained, not a consumed performance attempt. Source self-review also corrected
the freezer's inherited pin-schema name before its first execution. Actual
frozen-byte identity and 36 native preflights are required before publication.
