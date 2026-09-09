# Paired AUTO R199/32 KiB result: inconclusive controls

Bead: `leopard-79h.38.5.4.19.1.3`.

The single preregistered attempt completed, but **does not qualify the candidate
for promotion**. One of 100 same-path stability controls failed. Candidate
`45e2eff` remains default-OFF; all codec sources and archives are unchanged.
Do not retry, trim samples, pool attempts, or relax the 2% gate because the
aggregate failure is close to the threshold.

Preregistration `28bc6a98a0db00cbf5ddb15484a26134d67606d5` was committed and
pushed before any timing. Its [method review](paired_r19932_screen_method.md)
fixes both paired orders, grouped GF8 timing, all seven neighbors, the original
native comparator, 20 cross-process and 80 within-process control aggregates.
This is a distinct successor to consumed/inconclusive `aa2b034`, not its retry.

## Result

All 318 timed processes, 27 preflights and 26,712 measured spans completed with
exact record identities, unchanged frozen binaries and zero sibling activity.
Both paired directions passed the target and neighbor gates. Both native
comparisons passed their target gates. All target rounds were positive.

These are descriptive ratios from an **inconclusive** experiment, not qualified
speedup claims:

| Target | OFF/ON, `0110` | OFF/ON, `1001` | Native/ON |
| --- | ---: | ---: | ---: |
| K1000/R199/32768 ordinary encode | 1.522865 | 1.522156 | 1.419257 |
| Same shape, one-item batch | 1.512145 | 1.518331 | 1.415007 |

All 80 within-process controls and 19 of 20 cross-process controls passed.
The sole failed aggregate was cell 8 (K17/R7/64 GF8), cross-process same-OFF:
**1.0203606459901169**, exceeding the preregistered upper bound **1.02**.
Its three round ratios were `0.9894841931`, `1.0536901825`, `1.0189181673`.
Every raw sample and round remains retained.

Round 1's four process medians, each using all 84 grouped per-call averages,
were `167.046875`, `150.45703125`, `151.0859375`, `151.0859375` ns/call.
Those durations are grouped averages over 256 calls, not single-call latency.
All eight GF8 within-process control aggregates were close to parity. This
shows why within-process stability cannot substitute for the cross-process
control: the retained medians changed between processes even though their
internal ABBA controls passed. The cause remains unknown. Grouping did not
establish sufficient cross-process stability under the fixed gate.

The previous native/explicit-AVX2 aggregate failures did not recur in this
attempt. That is not proof of a noise fix or a causal explanation. Neither
attempt can supply promotion evidence by combining only its passing parts.

## Validation and evidence

Sixteen synthetic/adversarial tests passed normally and under `-O`, both from
source and frozen copies, before clocks. The independent raw replayer imports
no collector and executes no codec. Raw and read-only-bundle replays each pass
normally and under `-O`; all four outputs are byte-identical. They reconstruct
all 120 comparison/control aggregates and 360 round ratios, checking raw files,
sample normalization, call/route/workload identity, conditions, pins, source
commit, complete inventory and resource logs.

| Scope | Peak bytes | Limit | Memory events / swap |
| --- | ---: | --- | --- |
| Timed attempt | 143,904,768 | 256 MiB | all zero |
| Raw replay, both modes | 25,313,280 | 256 MiB | all zero |
| Retention and sealed replay, both modes | 45,113,344 | 256 MiB | all zero |

Every scope exited zero. The run was local to `work`, CPU26/sibling90/controller0,
under the canonical lock and pair lease, serial, with per-child CPU limits and
no swap. Slipgate/OBS/Forge disabled state matched before and after. No host
settings, unrelated process affinities, remote workers or production code changed.
Codex self-review and deterministic checks follow the user's Claude opt-out;
this is not independent-model `CONVERGED`.

Raw: `/tmp/leopard-paired-integration.la42Cz`.
Read-only: `.research/leopard-79h/paired-r19932-inconclusive.hGfIWZ`.
Its 784-file, 9,559,347-byte snapshot has manifest SHA-256
`2c1dd4857bc6dc35de7d6ce33421f6ef7de50e448f438d7743835b842df2d5c4`.
Sealed-replay/retention logs: `/tmp/leopard-paired-delivery.O6A0UM`.
The [machine-readable result](results/paired_r19932_screen_20260909.json)
contains every aggregate/round ratio and artifact/resource pins.

## Roadmap

Keep this R199/32 KiB AUTO extension off and retain the unexplained GF8
cross-process shift as a separate follow-up. No new timing attempt is authorized
by this failed result; any materially different future protocol needs its own
explicit review and committed/pushed preregistration, preserving the gates.

Resume the separate AVX2 forward/accumulating scheduling candidate (`b81f0a5`,
`leopard-79h.38.5.4.18.3`), whose untimed codegen/correctness qualification is
already complete. Its next step is a same-binary runtime-control qualification,
then a separate direct performance gate—not adding the earlier inverse-only
gain or inheriting these GFNI ratios. The three previously shipped AUTO GFNI
cells remain unchanged. The broader Leopard1/Leopard2 improvement goal stays open.
