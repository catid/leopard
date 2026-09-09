# AVX2 forward and accumulating pair codegen prototype

Tracking: `leopard-79h.38.5.4.18.3`, still in progress. Date: 2026-09-09.
Compile-time mode3 now passes the focused Release/full-sanitizer qualification
below. It remains **untimed and experiment-only**; no production source is
modified. Modes1/2 retain codegen-only evidence, not separate native qualification.

The [inverse-only experiment](avx2_pair_screen.md) completed below its fixed5%
target gate. Its candidate is not included here. This follow-up changes the
untouched forward pair and accumulating inverse pair, whose actual production
loops still reload spilled table vectors. The accumulating loop also reloads
source-x vectors. Prior results must not be added to predict a combined gain.

## Candidate and actual compiler output

[avx2_adjacent_schedule.patch](avx2_adjacent_schedule.patch) applies only to a
private exact3a2f064 source copy. Bits1/2 of the compile-time experimental
control select the forward/accumulating schedules; zero is unchanged.
GFNI and AVX-512 variants are excluded. The existing inverse in-place body,
tables, field representation, scalar tails, transform order, source policy
and tiling are unchanged.

The forward helper directly accumulates each nibble product into x. Its
empty GNU-asm dependency boundaries emit no hardware instruction. Data vectors
are passed by reference so the forward caller can reuse the physical y
vectors without retaining duplicate pre-boundary values. The accumulating
inverse form completes the disjoint y-output early and accumulates directly
into x. The four-buffer disjointness contract is in `Leopard2Backend.h`;
the execution checks below cover input preservation and scalar behavior.

The same original GCC Release recipe compiled all four modes, serially with
both fields retained. Mode0's entire AVX2 object is byte-identical to production
`bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d`.
All four actual objects contain no EVEX, high YMM, ZMM, ternary logic or GFNI.

| Actual 64-byte loop | Mode0 instructions / stack references | Enabled instructions / stack references |
| --- | ---: | ---: |
| In-place inverse, unchanged | 40 / 4 | 40 / 4 |
| In-place forward | 40 / 4 | 36 / 0 |
| Prepared-range forward | 39 / 3 | 36 / 0 |
| Accumulating inverse | 45 / 4 | 39 / 0 |

All listed loops retain eight byte shuffles. The forward-only and
accumulating-only modes isolate these changes; mode3 includes both.
Its object SHA is
`514f4de1d07b5325673d2a828ed447854db5c69427be7acf5edfe080c2fb1c1a`.
Full addresses, all mode hashes and raw static counts are in the
[codegen result](results/avx2_adjacent_codegen_20260909.json).
Instruction reduction is not a time share, speedup estimate or correctness
proof. A longer dependency chain or changed code layout may offset it.

## Checks and remaining gates

Normal and optimized read-only audits re-disassemble the actual objects,
check the saved disassemblies, exact source/object hashes and full-object ISA
ceiling. The zero-context tracked patch applies without fuzz to production
and reproduces the compiled source exactly:
`be5a07759cd0968299ebf243be77b299afda78065de704dcb4a319ca551ba205`.
Eight parser/adversarial tests pass in normal and optimized Python, including
missing/duplicate functions, nested-range versus disjoint-loop ambiguity,
malformed instruction bodies and each excluded ISA family.

The build peaked at182,116,352 bytes under512MiB, with all six memory-event
counters and swap zero. Parser tests and audits ran under256MiB/no-swap.
An initial private-source identity failure stopped before compilation; the
retained source copy had two older comment lines and was read-only. Correcting
only that owned copy restored exact3a2f064 input identity. The initial audit
also correctly rejected a range's nested loops as ambiguous under the older
single-pair parser. The new parser identifies the sole innermost eight-shuffle
loop and continues to reject distinct competing loops; its regression tests
cover both cases. Both failed logs remain retained.

Current execution counts and focused semantic qualification are now complete
as described below. A fresh same-binary control and committed-and-pushed
preregistration are still required before any timing; retain5% target/2%
control-neighbor and zero-sibling gates. This milestone is not promotion.

Raw workspace: `/tmp/leopard-avx2-adjacent.BBKeav`. Codex self-review and
deterministic/adversarial checks are the review provenance; the user's Claude
opt-out remains in force, with no independent-model `CONVERGED` claim.
All work is local, with no subagents, worker SSH or unrelated affinity changes.

## Retained milestone

The local read-only bundle `.research/leopard-79h/avx2-adjacent.gibc2g30`
contains72 files totaling20,771,674 bytes. Its outer manifest SHA-256 is
`ee1ace97d5464e44583ede141d24a30b6921de749350a5d698a36d93fbedafe1`.
It includes exact private sources, all four objects/disassemblies/recipes,
both initial failures, final audits/tests and their parser dependencies.
The sealed-copy audit reproduces the numerical result exactly. Retention
peaked at12,632,064 bytes under256MiB, with all six memory-event counters
and swap zero. That initial bundle predates native qualification and contains
no benchmark or promotion evidence.

## Focused semantic qualification completed

The mode3 Release archive is
`935cbcac7fe6ca0ec6ca7255a32b6991e60d72a843ba4321ee67070999bbc175`;
the full ASan/UBSan archive is
`88d4cb0e8ddd4d8bd040346af06527782102d32bb586b7628f1105322dd15b0f`.
Only the AVX2 member is replaced; each archive's other23 members remain
byte-identical to the original both-field production archive. The Release
member remains the previously inspected mode3 object. Neither archive contains
the driver-only compatibility control or callback observer.

Both profiles passed62 positive native records in total. Per profile:

- Every65,535 ordinary log for forward pairs and accumulating inverse pairs,
  compared with the unchanged scalar backend; input preservation and guards.
- 918 boundary accumulations covering zero length, scalar/vector boundaries,
  long payloads, misalignment and three repeated XOR applications, including
  exact cancellation after the second application.
- 384 forward-range cases, all eight zero-skew masks, both fused hints and
  distances1/4/16; plus the retained66,147 inverse/pair and256 inverse-range cases.
- Eight public full/subset/bounds shapes, three decode round trips and
  four-thread both-field round trips.
- Original-production observer delegation tests and eight candidate plus
  eight observed-production single-encode parity checks.

All32 full comparisons, totaling217,712,384 bytes, match separately linked
original native Leopard1 for all eight cells. Two deliberate clock guards
abort before reading a benchmark clock, and20 malformed/timing CLI requests
are rejected. Total raw records:84. The reused parity frontend's `off` argument
sets only its driver-compatibility stub; **the candidate library is fixed at
compile-time mode3**, and its archive/member hashes bind that active body.

The unchanged production archives supply callback counts. An optional workload
include reuses the existing16-entry exact-delegation observer with the new
eight-cell frontend; the default old workload remains available. Both original
FF16 compile recipes disable legacy in-field SIMD, so the observer's private
Ops table does not change the selected algorithm. No shared Ops table changes.
Release and sanitizer observations agree exactly and match the independent
structural traversal model. Range-internal pair counts below are derived from
observed distances and zero-skew masks, not independently instrumented again.

| Cell | K/R/bytes and route | Forward pairs / 64-byte blocks | Accumulating pairs / 64-byte blocks |
| --- | --- | ---: | ---: |
| 0 | 1000/200/32768 AVX2 | 665 / 340480 | 384 / 196608 |
| 1 | 1000/199/65536 AVX2 | 1330 / 680960 | 768 / 393216 |
| 2 | 1000/200/65536 AVX2 | 1330 / 680960 | 768 / 393216 |
| 3 | 4096/512/4096 AVX2 | 1793 / 114752 | 1792 / 114688 |
| 4 | 1000/199/32768 AUTO→AVX2 | 665 / 340480 | 384 / 196608 |
| 5–6 | explicit GFNI / AUTO→GFNI | no affected AVX2 pair | no affected AVX2 pair |
| 7 | GF8 | GF16 observer bypassed | GF16 observer bypassed |

These counts are not runtime percentages, memory-traffic measurements or a
bound on speedup. No effects are pooled with the earlier inverse-only screen.
The [qualification result](results/avx2_adjacent_qualification_20260909.json)
contains all operation counts and archive identities.

The serial build peaked at236,920,832/512MiB; native qualification peaked at
245,288,960/256MiB. All six memory-event counters and swap were zero. The
collector-free replay rechecks raw records, exact cell/parity bindings,
production and candidate members, source pins, full parity and the independently
derived counts without executing a codec. Normal and optimized replay pass.
Twenty pure model/parser/adversarial tests also pass in both Python modes;
the old six-cell observer/model tests remain included. Self-review tightened
the new replayer's label-to-cell, parity-order and production-archive bindings;
mutated and swapped public-record tests cover these guards.

The separate native-qualification bundle is
`.research/leopard-79h/avx2-adjacent-qualified.wxkx9v2b`:331 files,
341,920,323 bytes, outer manifest SHA-256
`6c798d28d47f5d69834711402cc4f1969cba621b16ea67e1e5acadb1999a5bc8`.
It retains the exact linked drivers, archives, recipes, raw native outputs,
parity bytes, tests and replay dependencies, alongside the earlier static
evidence and failures. Its read-only full manifest and sealed-copy semantic
replay both pass and reproduce the published numerical result exactly.
Normal/optimized delivery replays peaked at61,423,616/64,155,648 bytes;
retention at35,930,112; full-manifest plus sealed replay at178,704,384.
Each stayed below256MiB with all six memory-event counters and swap zero.
Separate retention and verification logs are in
`/tmp/leopard-adjacent-delivery.xLB13y`. The original native-L1 reference and
production archives remain separately pinned local dependencies of replay.

## Next performance priority

Cell4 is a distinct remaining AUTO boundary, not a measured current deficit.
Follow-up tracking: `leopard-79h.38.5.4.19`.
Production `UseAutoGF16GFNIEncode` and the new observer both confirm R199/32KiB
still uses AVX2. Earlier qualification left it unchanged and did not time it
against GFNI or native Leopard1. Measure those three unchanged production
paths under a fresh preregistration before considering any default-policy
extension. This is motivated by nearby GFNI wins, not proof they transfer.
The qualified AVX2 prototype stays available for a later direct family-wide
experiment; its timing and production gates remain open.
