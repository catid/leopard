# AVX2 forward and accumulating pair codegen prototype

Tracking: `leopard-79h.38.5.4.18.3`, still in progress. Date: 2026-09-09.
This is a distinct, **untimed and not yet correctness-qualified** prototype.
No codec has been executed using it and no production source is modified.

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
preserving input bytes and scalar behavior still needs execution tests.

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

Next, establish current execution counts for the explicit-AVX2 targets and
real AUTO-AVX2 neighbors; qualify the new forward and accumulating schedules
against ordinary-log scalar oracles, boundaries/guards, partials, both fields,
concurrency and exact independently linked Leopard1 parity in Release and
full ASan/UBSan/leak builds. Only then consider a fresh same-binary experiment
with committed-and-pushed preregistration,5% target/2% control-neighbor and
zero-sibling gates. This stage authorizes no timing attempt or promotion.

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
and swap zero. Nothing from this milestone has been benchmarked or promoted.
