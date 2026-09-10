# Three-epoch paired control-shift diagnostic: reviewed design

Tracker `leopard-79h.38.5.4.19.1.4.2`, 2026-09-10.
Status: design review only, **not an executable preregistration or permission
to read benchmark clocks**. No codec/default change. R199/32KiB stays OFF.

## Evidence and question

The consumed paired screen failed a GF8 same-OFF cross-process aggregate of
1.0203606459901169 against its fixed1.02 bound. The differing process medians
persisted through all84 grouped spans; within-process controls passed.
The [retained audit](gf8_shift_audit.md) cannot recover those processes' maps,
buffers or frequency state. Tower/adjacent GF8 requested explicit AVX2, whereas
paired GF8 requested AUTO: their absolute durations are not interchangeable.

The [metadata qualification](paired_runtime_metadata.md), pushed as9057974,
now provides observations without real clocks. An exhaustive projection of
its24 positive Release GF8 records (all four schedules, both groups, and all
three check/exercise/clock-exercise modes) finds:

| Recorded quantity | All24 new Release GF8 processes |
| --- | --- |
| Distinct executable load biases |24 |
| Source/reference data offsets modulo4096 |896 /2176 |
| Scratch/output data offsets modulo4096 |1280 /3200 |
| Input/output pointer-array offsets modulo4096 |2704 /3728 |
| Scratch minus source data address |12672 bytes |
| Output minus source data address |14592 bytes |

The16 abort records share main/encode ELF offsets19392/186176; the8 synthetic
records share19456/186752. These are two differently linked binaries, not
within-binary code drift. Every recorded endpoint in each process agrees.
These observations concern the new instrumented binaries, not historical maps.
They do not rule out page-level placement, physical cache mapping, predictor,
frequency or scheduling effects, and are not performance observations.

Projection inputs are the sealed `checks/release-8-*.stdout` files and
`checks/checks.json` in `.research/leopard-79h/paired-runtime-metadata.oQioRv`,
whose full manifest/replay verification is recorded in the qualification.
Select only returncode0 records with args profile=`release`, cell=`8`;
group by the explicit clock-binary kind, never by an observed duration.
Compute each address modulo4096 and each difference with exact integers.
The accompanying JSON checkpoint retains the two complete projection groups.

The next question is narrower and testable: **does a between-process cost
offset decay, or persist, across three fixed complete observation epochs while
that process retains its codec, buffers and recorded mappings?** This changes
the observation horizon. It is not extra warmup followed by another attempt
to pass the old integration screen.

## Fixed proposed workload

Keep the original318-process inventory, order, three rounds, nine workloads,
both paired directions, original native-L1 comparator and all cross-process
same-path controls from `paired_r19932_screen_method.md`. Do not omit a cell,
move CPU, substitute restricted L1, change GF8 AUTO to AVX2, or retry old data.

Allocate/initialize once per fresh process, then execute exactly three epochs
indexed0,1,2 with no adaptive extension, sleeping, printing, new driver buffer
allocation or intentional eviction between epochs. Every epoch contains:

1. Before-epoch endpoint snapshot.
2. The original four single-call preflights, one per schedule slot, including
   actual GFNI route probes and normalization back to production mode.
3. Four warmup schedule passes, followed by21 sampled schedule passes.
4. After-epoch snapshot. Defer all JSON/file output until all epochs finish.

Use group256 only for GF8, group1 elsewhere. A group still means repeated
complete public calls; the batch cell uses one-item batch calls. Preserve
public-call lambdas, grouping, result checks and per-call accounting, and all
post-group parity/guard checks. State selection and metadata remain outside
each duration. There is no overhead subtraction or host-wide cache operation.

Reference parity is established by epoch0's first preflight. Later epochs must
verify their first output against it before any reference refresh, so repeated
preflights cannot hide a persistent wrong result. Require all input hashes and
reference/output parity semantics to remain unchanged across epochs.

Each process has312 selection records, six complete snapshots,252 measured
spans and504 clock endpoints. Public encode calls total312 with group1, or
76812 with group256: three times `(4 + 100*group)`. Per-epoch/per-slot calls
are `1 + 25*group`. Total318-process measured spans:80136. Repeated epochs are
correlated observations; there are still318 launches, not954 independent ones.
The27 untimed campaign preflights remain separate from timed processes.

## Analysis and decisions fixed before any timing

Reconstruct the original estimators independently for each epoch index:
median21 within-process paired pass ratios; median84 homogeneous process cost;
the same cross-process ABBA ratios and three-round geometric aggregates.
Keep paired direction, cell, process slot and round separate. Never average
epochs together, drop epoch0, select the fastest epoch or pool old attempts.

Apply every inherited performance requirement separately to every epoch:
100 same-path controls per epoch (300 total),2% symmetric unchanged-neighbor
equivalence in each paired direction,5% paired-target/native gates with every
target round positive. Retain individual control excursions even when their
predeclared aggregate passes. Failure precedence stays unchanged. These gates
provide diagnostic context; **even all-pass cannot enable the candidate**.

For every homogeneous-control process, also retain epoch1/epoch0 and
epoch2/epoch0 cost ratios, plus each epoch's median cost, raw84 spans, phase
identity and six endpoint snapshots. Report every process, not a selected
slow subset. Do not apply those ratios to mixed OFF/ON process medians or
interpret three correlated epochs as three independent trials.

No causal automatic classifier is proposed. Complete data can show an offset
that changes over this window or persists through epoch2. It cannot identify
ASLR, frequency, caches, predictor state or a launch defect as the cause.
Stable endpoints exclude only observed endpoint relocation, not transient
relocation or unobserved state changes. Snapshot/probe/parity work itself can
condition later execution even though it is outside the timed spans.
This instrumented three-epoch window may still be too short to discriminate.
An inconclusive outcome remains useful retained evidence, never retry authority.

## Required next qualification and publication boundary

The existing metadata frontend supports only two snapshots/104 selections,
permits only abort/synthetic clocks, and deliberately refuses `--measure`.
It is not a timing-capable three-epoch frontend. Do not weaken or overwrite it.
Implement a separate pinned overlay and explicit epoch-indexed schema.

Before timing, qualify the expanded frontend with unchanged native-L1 and
both-field L2 Release/fullASan+UBSan+LSan archives, all nine cells, both GF8
groups, all schedules and check/exercise/synthetic/abort modes. Verify:

- All312 selections, six snapshots and epoch/pass/slot/public-call witnesses;
  actual observed probes remain distinct from queried selection.
- Equality across all six snapshots, not just each endpoint pair; native
  parity aliases and up-to1024 native output pointers remain intact.
- Exact252/504 span/clock counts and strict epoch-major accounting, including
    malformed-clock, missing/duplicated/reordered epoch and capacity failures.
- Full native-reference parity, input preservation, guard/tail/API semantics,
  and real object/ELF/source/build provenance. No machine-code identity claim
  from source-identical lambdas or unchanged library objects alone.
- No driver print/allocation/metadata/probe work inside grouped public-call spans;
  all six snapshots retained until final output, with precise resource limits.

Those complete-path totals apply to complete exercise/synthetic/timing modes.
Check-only and expected early-abort modes need their own explicit inventories;
do not demand252 spans or312 selections from a path that intentionally stops
earlier. The final reviewer specifically calls out this implementation boundary.

Keep serial local canonical-lock execution, builds512MiB/checks256MiB/no swap,
compilerGC10/4096, CPU/wall bounds, and no unrelated affinity or host changes.
Then separately qualify the steady-clock frontend and collector/replayer,
freeze immutable artifacts, independently review the exact method and commit
AND push the one-attempt preregistration before clocks. Preserve local CPU26/
sibling90, zero-sibling and stopped-service gates. All failed evidence remains.

The user-authorized local read-only reviewer found the three-epoch contrast
scientifically distinct **as a diagnostic**, with the interpretation and
qualification restrictions above. Main adopts its six-snapshot, complete-
epoch and homogeneous-control recommendations. Review is bounded Codex
inspection, not Claude or independent-model fixed-point convergence.
The overall optimization goal and R19932 integration remain incomplete.

## Design verification

The final bounded local reviewer found no material method/count/interpretation
error; it independently checked the complete-path arithmetic, not the recorded
address projection. A separate stdlib projection verifier checked every selected
stdout hash against both the qualified checks inventory and pinned outer
manifest, reconstructed both clock-kind groups and all values above, verified
actual AUTO/AVX2/non-GFNI records and checked the proposed count arithmetic.
Normal/optimized outputs agree; peaks10,366,976 /13,238,272 bytes under256MiB,
exit0/all six memory events0/swap0, under the canonical lock. No codec executed.
These checks do not qualify an unimplemented three-epoch frontend.

Projection source/logs: `/tmp/leopard-paired-epoch-method.mgjVSL`.
Source SHA256:
`1a9204a5904517ea081f2d72c40310e8743c8dc7401ac0845ca0026f6dec0cc5`.
Normal log SHA256:
`d7f861c63d81e2aca5450d7202ec65f452515847fb9653907f73745afe0e88a2`.
Next separate implementation/qualification child:
`leopard-79h.38.5.4.19.1.4.3`. This design's publication is not a timing
preregistration and does not authorize launching the old or proposed screen.
