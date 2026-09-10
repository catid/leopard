# Direct tower encoder screen: negative result

Tracker `leopard-79h.38.5.4.18.4.4`, 2026-09-10.
Preregistration `076c6df289a375b37956dacac6366b2ceec346c9` was pushed
before the sole attempt. Production sources and defaults are unchanged.

## Decision

All 714 timed processes, 36 preflights and 59,976 spans completed, with exact
frozen identities and zero sibling activity. All 32 cross-process and 128
within-process aggregate stability controls pass. The fixed decision is
**reject_unchanged_neighbor**: excluded GF8 cell7 original/ON is 1.079044,
outside the symmetric 2% interval. This positive shift is not a slowdown, but
it is not an allowed unchanged-path equivalence result.

Independently, all five targets fail all 10% gain gates. Direct original/ON
is essentially flat at 32KiB and about 7.8% lower throughput at 64KiB. The
candidate stays **experimental/default-OFF**. This rejects this implementation
for promotion, not the algebra's correctness or every possible tower design.

Ratios are candidate throughput divided by comparator throughput; larger is
better. OFF is not pristine production. No ratios are multiplied across runs.

| Target K1000 | Paired0110 ON/OFF | Paired1001 ON/OFF | ON/original L2 | ON/native L1 |
| --- | ---: | ---: | ---: | ---: |
| R200/32KiB AVX2 | 0.985183 | 0.985459 | 0.998354 | 0.915587 |
| R199/64KiB AVX2 | 0.936308 | 0.938683 | 0.922334 | 0.887791 |
| R200/64KiB AVX2 | 0.932571 | 0.935021 | 0.921905 | 0.879972 |
| R199/32KiB AUTO | 0.981550 | 0.987725 | 0.996105 | 0.917097 |
| R200/32KiB one-item batch | 0.979768 | 0.978033 | 0.996999 | 0.918795 |

Native cell8 is ordinary encode, not a native batch API. Every native/ON target
round is below one. The actual runtime OFF/ON comparison is negative too;
fewer shuffles and no loop spills did not produce an end-to-end benefit.
Extra XOR/dependency work, conversions and cache behavior are possible costs,
not separately timed causal attributions. The added6MiB table cache is still
present and cold initialization remains unmeasured.

## Excluded neighbors and retained control outliers

| Cell | ON/original L2 | Paired0110 | Paired1001 | Gate |
| --- | ---: | ---: | ---: | --- |
| 3 | 0.998377 | 0.999867 | 1.000197 | Pass |
| 5 | 0.998740 | 0.999758 | 0.999940 | Pass |
| 6 | 0.996921 | 0.998597 | 1.001618 | Pass |
| 7 | 1.079044 | 1.000134 | 1.000053 | Fail: direct original comparison |

GF8 original/OFF is also shifted positively, at1.073332. Mixed paired ratios
are almost exactly one. This is a useful observation about original versus
experimental executable comparisons, not proof of a code-layout, ASLR,
frequency, allocator or cache cause. Old GFNI/adjacent failures are not pooled.

The preregistered stability gate applies to aggregates. Four individual
control rounds fall outside2% and are retained: cell4 same-native round2
1.0291638264; cell7 same-original rounds0/1/2 at1.0205538357,
1.0551423800 and0.9484256168. No sample or round was excluded. There are no
confidence intervals or a claim that every individual control round passed.

## Verification and scope

Both collector-free raw replays pass in normal and optimized Python. They
independently reconstruct201 aggregates and603 round ratios; verify51 frozen
inputs, full raw records, preflight identities, link/source provenance,
shutdown conditions and resource evidence; and agree with the collector.

| Scope | Peak bytes | Cap | Outcome |
| --- | ---: | --- | --- |
| Timing, 542.54seconds wall | 149626880 | 256MiB | exit0/events0/swap0 |
| Raw replay normal / optimized | 32477184 /34328576 | 256MiB | exit0/events0/swap0 |
| Six retention tests normal / optimized | 15589376 /21250048 | 256MiB | exit0/events0/swap0 |

The prior21 protocol tests in each Python mode and36 frozen clock-free
preflights remain valid, not rerun. No codec or driver was rebuilt.
Correctness remains the separately retained public/encoder qualification.

Raw root: `/tmp/leopard-tower-screen.0mjTlC`.
Read-only bundle: `.research/leopard-79h/tower-screen-negative.V1230H`,
1650 files /15311540 bytes, outer manifest
`d79c3ff53e51fc25df18d3867f4738e2971bafb0f0f104895954afd600bfca75`.
All1649 manifest members equal raw files, and the entire namespace is verified
read-only. Executables and other files are private copies, not raw hard links.
The machine-readable result preserves every ratio.

Raw and retained normal/optimized replays produce four identical canonical
results, SHA256 `f44596d137ea03b8f980528f08f5f71ee9f09eff3af20918f929c5315b4cf688`.
Retention peak17047552, sealed replays46575616/33275904, final audit28979200
bytes: all under256MiB, exit0, all six memory events0 and swap0. Delivery logs
and the independent manifest/namespace/raw-equality audit are at
`/tmp/leopard-tower-screen-delivery.LlssYu`.

The read-only retention reviewer found a reserved-root-manifest collision;
the retainer now rejects it before copying. All seven final retention tests
pass in normal/optimized Python (15822848/20664320-byte peaks, all events0,
swap0). The initial six-test sources/logs and one wrong-working-directory test
invocation are preserved separately; that invocation failed before importing
the test module and is not relabeled passing. No timing was repeated.

The sealed report/result capture the pre-retention checkpoint. This current
repository report/result append verified final retention metadata rather than
mutating that sealed snapshot. No qualification or benchmark source changed.

## Next action

Do not retry this tower attempt, trim it, loosen its10%/2% gates, or enable the
candidate. Resume the existing untimed tiny-GF8 cross-process investigation
`leopard-79h.38.5.4.19.1.4`, now informed by a shifted original/OFF comparison
as well as original/ON. Diagnose the actual frozen code/comparison boundaries
before proposing any distinct future measurement method. The promising
R199/32KiB AUTO-GFNI extension stays OFF until it has valid integration evidence.

Local work only. The read-only reviewer and deterministic checks are not a
Claude fixed-point CONVERGED claim. The overall Leopard1 performance goal and
remaining explicit-AVX2/AUTO deficits stay open.
