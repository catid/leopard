# AUTO R=199 / 32 KiB integration screen

Bead: `leopard-79h.38.5.4.19.1`. **Completed, inconclusive controls; default
remains off.** Preregistration `aa2b03438503338ab622af29ebbe2f41234fc37d` was
committed and pushed before clocks. The single attempt is consumed and may not
be retried. The original protocol below remains unchanged.

The [default-off candidate](auto_r19932_candidate.md), source commit
`45e2effd869859c9b3aa48190eff6f4738817c61`, has passed focused Release and full
ASan/UBSan/LSan correctness, full native Leopard1 parity, real backend-failure
tests and Release ISA checks. That evidence is retained unchanged, including
the native qualification's 612 memory.max events and its no-OOM/no-swap result.

## Fixed comparison

The [machine-readable plan](auto_r19932_screen_plan.json) fixes all nine cells,
identities, orders, resource/host requirements and the single attempt location.
The existing qualified frontend and archives are copied, not rebuilt:

| Input | SHA-256 |
| --- | --- |
| Same-binary candidate frontend | `6caf97218b3240116bc7e6bf873cc82e4236a2bedff6690e07ec3332c59966a7` |
| Native Leopard1 frontend | `4f132c94c568913740b9266e7998357cc214ae6861b97f4301b943180c2bc695` |
| Candidate archive | `89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334` |
| Original native Leopard1 archive | `3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1` |

Targets are ordinary encode and one-item batch at K=1000/R=199/32768 bytes.
The seven unchanged neighbors are the three established GFNI cells,
R=198/32768 bytes, explicit AVX2 at the target shape, K=4096/R=512/4096 bytes,
and GF8 K=17/R=7/64 bytes. OFF is the changed candidate's disabled state,
not a claim of pristine production binary identity. Native Leopard1 remains
the original product comparator, not the ISA-restricted Leopard1 variant.
For one-item batch, the native comparator is one ordinary `leo_encode` call.

Every cell uses OFF/ON/ON/OFF, OFF/OFF/OFF/OFF and ON/ON/ON/ON groups.
Both targets additionally use native/ON/ON/native and a four-native control.
There are three rounds and 21 measured samples per process after one untimed
route check and four warmups: **372 timed invocations and 27 untimed preflights**.
The operation, source generation, output semantics, scratch size, encode count
and route identity must match the clock-free qualification.

Each process contributes its median. A four-process group produces
`sqrt((median0 / median1) * (median3 / median2))`; three round ratios combine
geometrically. Ratios greater than one favor ON. All values and individual
round excursions are retained; no trimming, pooling or retry is permitted.

## Pass conditions, fixed before clocks

- Each of the two targets must improve at least 5% over OFF, with every round
  positive. Each must also improve at least 5% over native Leopard1, with every
  round positive. The native criterion is explicitly part of this fresh plan.
- All 20 same-path aggregate controls must lie in `[1/1.02, 1.02]`.
- All seven unchanged OFF/ON neighbor aggregates must lie in that same interval.
- Ten-second passive observation and every timed invocation must have zero
  non-idle SMT sibling jiffies. No unrelated process affinity may be changed.
- Local `work`, CPU26/sibling90, controller CPU0, canonical campaign lock plus
  pair lease; one 256-MiB/no-swap scope, serial children and CPU/file-size caps.
- Original source and binary SHA-256 identities are verified before and after
  execution. Slipgate/OBS remain stopped, disabled and non-restarting.

Decision precedence is control failure, neighbor failure, AUTO target failure,
native target failure, then `qualify_default_on_artifact`. Passing permits
preparing and verifying the actual default-on artifact; it does not itself
flip the default. No confidence intervals or broad v19 qualification are claimed.

Frozen files: `/tmp/leopard-auto-r19932-screen.oSPt31/frozen`.
Only allowed attempt: `/tmp/leopard-auto-r19932-screen.oSPt31/attempt1`.
The collector verifies that the exact plan and sources belong to a commit
already pushed to the topic branch before beginning. Changing any gate after
samples are observed would require a distinct authorized experiment.

Pre-timing checks: 20 pure protocol/independent-replay/qualification tests pass
in normal and optimized Python modes. All frozen inputs and the complete
earlier qualification were replayed without codec execution. Preparation peak
80,752,640 bytes; freeze/replay/tests peak 90,071,040 bytes, both below 256 MiB
with all six memory events zero and swap disabled. Review is Codex self-review
and deterministic/adversarial tests; Claude and subagents remain opted out.

## Completed result: do not promote

All 372 timed invocations and 27 preflights completed, with exact recorded
workload/route identities and zero sibling non-idle jiffies. The passive
observation was 10,000,475,534 ns, with the sibling counter unchanged at 570683.
The fixed decision is **`inconclusive_controls`**, not an accepted optimization
and not evidence that the candidate regresses.

| Target API | OFF / ON | Native / ON | Qualification |
| --- | ---: | ---: | --- |
| Ordinary encode | 1.5388704673 | 1.4549241629 | Inconclusive controls |
| One-item batch | 1.5263170544 | 1.4709829783 | Inconclusive controls |

The numerical target ratios are positive in all three rounds and exceed both
5% thresholds. All seven unchanged-neighbor aggregates are inside the fixed
2% interval. However, three of the 20 same-path aggregate controls fail:

| Unchanged control | Aggregate ratio | Required interval |
| --- | ---: | --- |
| Native Leopard1, ordinary target (cell 0) | 0.9746447598 | [0.9803921569, 1.02] |
| Explicit AVX2, same ON (cell 6) | 0.9758334085 | [0.9803921569, 1.02] |
| GF8, same OFF (cell 8) | 0.9516373686 | [0.9803921569, 1.02] |

No outlier, round or control is removed. The apparent 52.6–53.9% improvement
over OFF and 45.5–47.1% improvement over native are **not qualified speedup
claims**. The previous independent direct screen is not pooled with this run.

Read-only inspection identifies different symptoms, not a proven host cause:

- Native cell 0, round 2 same-path process medians are 4.132661, 4.405398,
  4.363876 and 4.158098 ms. The slower middle processes contain sustained
  slower samples, not merely one maximum that could explain their medians.
- Explicit AVX2 cell 6, round 0 medians are 4.476764, 4.808288, 4.780066 and
  4.498366 ms. This millisecond-scale shift cannot be explained by the tiny
  GF8 call's nanosecond-scale granularity.
- GF8 same-OFF process medians range from 150 to 190 ns. Round 1 medians are
  150, 190, 170 and 161 ns; its control ratio is 0.864683768. Timing individual
  calls at this scale is too coarse and variable here for a stable 2% claim.
  This observation does not justify weakening the bound.

No system-frequency, scheduler, allocation-placement, cache or thermal cause
was isolated by these records. The result does not authorize unrelated process
affinity changes, host configuration changes, CPU movement or a timing retry.

## Retained evidence and next work

[Machine-readable raw replay](results/auto_r19932_screen_20260909.json).
Collector-free replay passes in normal and optimized Python and agrees with
the collector's fixed decision and all 93 round / 31 aggregate ratios. It also
revalidates all 22 frozen inputs and the earlier 223-positive qualification,
including 297,765,056 full parity comparison bytes and its resource caveat.
Read-only retained-copy replays and its full manifest check pass.

Native scope peak: 128,798,720 / 268,435,456 bytes. Raw replay peak:
98,013,184 bytes. Retention plus sealed replays/manifest: 119,652,352 bytes.
These scopes exited zero, with all six memory events zero and no swap.

Read-only bundle:
`.research/leopard-79h/auto-r19932-screen-inconclusive.rexj13hd`.
842 files, 7,127,981 bytes; manifest SHA-256:
`fca53b9ad4da9da6190878e64f01be1c16111f008230459a89d69388c13f54ee`.
Raw attempt journal SHA-256:
`00a9b20bb3f1d077a30154c14d8b6d08896b6f967b9f01b8e3487fb96afe39e0`.
Scope log SHA-256:
`5baf809bebf1416b353d770c9148c3eadbb9f0da0baa44c0d8da15ed424c2293`.
Delivery logs: `/tmp/leopard-auto-r19932-screen-delivery.9nv5Aj`.

Next Bead `leopard-79h.38.5.4.19.1.1` is **untimed methodology qualification**:
preserve the qualified codec and all nine cases, check a paired OFF/ON frontend
using the same codec/buffers within one process, and amortize clock overhead
across grouped repeated public calls for tiny workloads. This is a hypothesis
for improving measurement stability, not a demonstrated fix for the native or
AVX2 process shifts. Group averages must not be mislabeled single-call latency.

That child may perform clock-free correctness, call-count, route, bounds,
parity and clock-abort checks only. It may not execute a benchmark, retry this
attempt, pool data or relax existing gates. A timed successor requires its own
explicit review and committed/pushed preregistration. Actual default-on
artifact verification and production enablement remain after that gate. The
integration, parent performance task and full user goal remain open.
