# GFNI inverse distance-one forwarding experiment

Bead: `leopard-79h.38.5.4.13`. Date: 2026-09-06.
Status: initial correctness/structural milestone passed; timing not performed.

The preceding callback observer identified 2,000 first-stage two-way calls
which bypass the existing fused GFNI range kernel. This experiment changes
only that dispatch, retaining the copy-first policy, two 32-KiB passes, final
accumulation, forward transform and all existing backend kernels. It is not a
retry of the rejected out-of-place source-staging candidate.

## Source and activation boundary

Production files are unchanged. `gfni_inverse_distance1.patch` is an overlay
for a lane-owned copy of `LeopardFF16.cpp` from base codec
`36dc0c8f66604b8d974468e687c51e6183ecb61d`. Its SHA-256 is
`050bba9f88c69aeb681a6f1f912adc781b78eabc71f415249798cbfda2e823e4`.
The full dual-field Release/sanitizer archives are copied from the prior
current-route bundle, then only their FF16 object is replaced. All other 23
members and their ordering are byte-identical. The recorded original FF16
compile commands are reused with only source/output/include paths relocated
to the frozen lane. Both portable in-field SIMD exclusions remain enabled.

A default-off driver hook records and selects each high-profile source-policy
pass. It matches GFNI, K1000/R200, all 200 outputs requested, side256,
32-KiB execution bytes, 64-KiB source-policy bytes, and a present sparse-plan
descriptor with zero blocks. The boolean is passed through private encoder
helpers; only its nonaccumulating distance-one branch forwards directly to
the original `ff16_ifft_butterfly4_range` with the unchanged work pointers,
skews, bytes and false hint. Other branches call the old helper. The separate
low-profile call uses the default false argument. No global fused selector,
Ops identity, arithmetic kernel or public API is changed.

The match describes the aligned prefix, not every public shard length:
explicit GFNI at 65,538 bytes also matches its first two passes and leaves its
two-byte tail unchanged. AUTO at that neighboring length remains AVX2 and
does not match. Directed tests cover both distinctions.

The check driver includes the qualified callback observer and accepts only
`--check cell --fuse=0|--fuse=1 [parity_file]`. The source-policy hook has a
checked sixteen-pass capacity. It neither reads a clock nor permits the
existing workload's timing mode. The unchanged workload JSON's codec commit
identifies the **base**; the overlay and candidate archive/binary hashes below
identify the actual experimental source and execution, not an unmodified
production build.

## Verified structural result

At the AUTO/GFNI target, OFF reproduces the baseline's 4,132 callbacks. ON
has 2,632: exactly 2,000 first inverse two-way calls become 500 existing
four-way range calls at distance one. Every other histogram bucket and every
source-policy field is unchanged. All five other fixed routes retain the
baseline histogram in both modes. Counts are not CPU time, traffic estimates
or a demonstrated speedup.

- All 24 Release/sanitizer × OFF/ON × six-cell workload records match the
  original. Twelve full Release parity files match standalone exact Leopard1
  (`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`), totaling 122,028,032 bytes.
  Each sanitizer trace is byte-identical to its Release counterpart.
- Each build passes 480 in-place range comparisons against an explicit scalar
  two-layer schedule: fifteen byte lengths from zero to 65,598, exact allocation
  ends, offsets zero/17, all eight zero-skew masks, and both fusion hints.
- Each build passes fourteen public directed cases with unaligned exact-end
  buffers, repeated immutable inputs, two-byte tails, partial/omitted parity,
  shape/route neighbors, native odd-byte rejection, short scratch and overlap
  rejection. Omitted outputs, canaries, inputs and OFF/ON parity stay exact.
- The hook passes seventeen negative predicate tests, both mode selections,
  null versus empty descriptor handling, reset, and atomic refusal at sixteen
  passes. Four structural-model adversarial tests and the retained-only replay
  pass normally and with Python `-O`. Fourteen malformed or timing requests
  fail before entering the workload.
- All 54 successful native matrix/kernel/public scopes stay below 256 MiB,
  maximum 142,548,992 bytes, with all six memory-event counters and swap zero.
  Archive/driver builds peak at 127,016,960 bytes; unit builds peak at
  89,796,608, with the final corrected fixture at 88,375,296 bytes, all under
  512 MiB. No field, validation, tail or sanitizer check was removed.

## Retained diagnostic failures

The first OFF workload completed parity generation but failed the experiment's
selection-count assertion: my initial matcher incorrectly required a null
sparse pointer. Inspection of `leopard2.cpp` confirmed that the public encoder
passes a descriptor by address even when it has zero blocks. The corrected
matcher explicitly requires that observed descriptor state; all final hook
records prove it. The original source, binaries and failed scope are retained.

The first scalar-reference kernel test failed because my reference directly
passed the 65535 zero-skew sentinel to `Butterfly2`, whose documented contract
accepts ordinary fixed multipliers only. The corrected independent schedule
specializes zero skew into its retained XOR edge, as the field caller does.
The failure and old fixture/binaries remain retained. No production arithmetic
was changed to make either diagnostic pass.

## Frozen identities and retention

| Artifact | SHA-256 |
| --- | --- |
| Release archive | `39cfa6d4d5b06aae3837f036254b94db7730685271eeae82f68858265a26881f` |
| ASan+UBSan archive | `4ff63bbf6347fb528dde177385fa5104a83c0cd223a61ef7fb01749e1d29c84f` |
| Release check driver | `64d0767a3e611c6cf3b7fa20856ac44d10800449ca7272a67a9833cf73e438d9` |
| Sanitizer check driver | `c1f0778afc056ecb5106b332dcbfea60c83f86d907e083f7443baa686dae70d0` |
| Release directed test | `c954108c570112b4f81c21fb963b90be281ad4dbe091c5196ddfed258ca28cbe` |
| Sanitizer directed test | `64156ff11bd225a914bd1712e625884871b7f150823d468b64c6ab9da2a13740` |

The read-only bundle `.research/leopard-79h/gfni-inverse-distance1.EN8c4J`
is retained locally and on ripper, with all sources, archives, executables,
recipes, failed and final checks, raw parity/traces, model tests and replay
dependencies authenticated by `SHA256SUMS`. The referenced baseline bundle
remains `gf16-current-route-failed.STc10h`, outer manifest
`e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`.
Run its frozen `replay/verify_gfni_inverse_distance1.py` with the experiment
bundle path and baseline bundle's `preflight` directory, inside the established
256-MiB/no-swap scope. Normal and optimized replay execute no codec or collector.

## Next gate

The task remains open. A separate driver without per-callback instrumentation,
qualified trace capacity for the full warmup/measurement loop, and a newly
preregistered frozen-binary candidate/control screen are still required.
OFF/ON layout controls must pass before any speed interpretation. Production
promotion requires the existing target/neighbor gates and broader integration
checks; current-versus-Leopard1 and v19 qualification are not replaced by this
milestone. The old exhausted Leopard1 timing attempt has not been rerun.

Review uses Codex self-review and deterministic/adversarial validation under
the user's Claude opt-out, not independent-model `CONVERGED`. No Claude/API,
kernel setting, unrelated workload or affinity was used or changed.
