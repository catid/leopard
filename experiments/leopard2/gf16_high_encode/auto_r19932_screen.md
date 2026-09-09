# AUTO R=199 / 32 KiB integration screen

Bead: `leopard-79h.38.5.4.19.1`. Fresh preregistration; no integration samples
collected when this plan was written. This is not a retry of the exhausted
direct R199 screen or an automatic production promotion.

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
