# Bounded AUTO GFNI candidate: implementation checkpoint

Bead: `leopard-79h.38.5.4.17.1`. Date: 2026-09-09.
Status: implemented and correctness-checked, **default off**; performance
qualification and production enablement remain open.

The preceding [explicit-backend screen](gfni_boundary_screen.md) found GFNI
54.5% and 47.1% faster than current AUTO in the two measured deficit cells.
This checkpoint implements the corresponding AUTO selection, without enabling
it by default or treating explicit-backend measurements as AUTO qualification.

## Change

`g_auto_gf16_gfni_boundary_mode` is a separate atomic initialized-data control:
1 means enabled and 2 means disabled. The new private diagnostic setter must
be used before codec creation while operations and route inspection are
quiescent. The original global AUTO/GFNI switch, host qualification, AUTO
AVX2 baseline, thread count 1, K=1000/T=256, legacy-high GF16, native layout,
flags zero, full-output checks and API exclusions all remain in force.

| R / shard bytes | Boundary mode off | Boundary mode on |
| --- | --- | --- |
| 200 / 65536 | Existing GFNI route | Existing GFNI route |
| 200 / 32768 | AVX2 | New GFNI route |
| 199 / 65536 | AVX2 | New GFNI route |
| 199 / 32768 | AVX2 | AVX2 |

Other shapes and explicit backend requests retain their previous behavior.
R=199 codec setup only qualifies/caches the optional table when the boundary
and global controls are enabled. A codec created without that table cannot
bypass qualification by enabling the control later. The two new private
functions return false in GF16-disabled builds; their definitions and state
remain field-guarded.

## Build and checks

Only `leopard2.cpp.o` was rebuilt using the retained per-profile compiler
commands and fresh immutable source copies. An independent GNU-ar parser
verified that the other 23 of 24 object payloads are byte-identical to the
original Release and both-field ASan/UBSan archives. Adding private declarations
to `Leopard2Direct.h` does not change another translation unit's definitions.

- Core source SHA-256:
  `0e4cdc4485e96e6d1ca711ad0b5192925df01e7a4951bca3d40b35e580c061b7`.
- Private header SHA-256:
  `0c6ff3efdfbb8754cd5abebaf520f05241f16e4f5487f27b093aae4dbd2ca065`.
- Release archive:
  `80bffc9e873585d9a18fcf6a294413b8c2a76d84fa6e7f9a437fb96ca8458533`.
- ASan/UBSan archive:
  `9a87a6baa349195e02c810fab8131fe5d4edfaf86de5a0efd1b7cafedb207c38`.

All 60 focused native commands passed: 30 Release and 30 with ASan, UBSan and
leak detection. They cover the old default-on route; boundary-off/on and
global-off behavior; late-enable safety; negative context, codec, byte and
processor identities; actual public route-call accounting; partial outputs;
ordinary one-item batches; excluded scalable aliases, reusable bindings and
multi-item batches; invalid-batch atomicity; decode round trips; and two
threads using one immutable codec with independent buffers. Eight guarded
shapes run in both boundary modes, covering aligned/unaligned buffers, even
tails, GF8/GF16, short scratch and odd GF16 rejection. The prior guard fixture
was parameterized for AUTO without changing its original explicit-GFNI default.

Fallback cases use link wrappers around the host predicate and optional
`GetQualifiedOps` return boundary, injecting unavailable/OOM/KAT failure
results. They verify that the new selector retains AVX2 and correct bytes.
They are not a fresh replay of the underlying backend allocation fault,
whole-backend KAT fault, or cached-failure consumption tests. The wrappers are
only in check executables, not in the candidate archives or future timings.

Archive build: 15.59 seconds, peak 181,243,904 bytes / 512 MiB. Check-driver
build: peak 94,679,040 bytes / 512 MiB. All native checks: 33.29 seconds, peak
175,464,448 bytes / 256 MiB. Every successful scope exited zero with all six
memory-event counters and swap zero. The canonical campaign lock serialized
all builds and checks. No SSH worker, other-process affinity change, Claude
call, or subagent was used.

`verify_auto_gfni_boundary_checks.py` passes normally and under `python -O`.
It independently hashes 33 source inputs, archive payloads, compiled checks,
all raw output/error records and resource counters; it also rejects four
typed-record mutations. Codex self-review and deterministic checks are used
under the user's Claude opt-out; there is no independent-model `CONVERGED`
claim.

Two preparation errors are retained: the initial staging audit incorrectly
included a CMake file whose only historical difference added unrelated Python
CTest registrations (no compiler was launched); the first offline verifier
assumed aligned scratch geometry for ragged buffers. The corrected audit uses
the actual recorded compilation input list, and the verifier retains the
pre-existing 64,000-byte tail-staging requirement, confirmed against the prior
guarded checkpoint. No native failure or performance sample was discarded.

## Remaining qualification

There have been **zero new timed encode invocations** for this candidate.
The next step in the same Bead is a new candidate-specific driver and frozen
preregistration, with actual AUTO route probes normalized before timing.
The boundary-only control preserves the old GFNI target as an unchanged
control. Both new ordinary-encode targets and their admitted one-item batch
paths need measurement, together with representative inactive neighbors and
same-binary controls. Keep the 5-percent target, 2-percent control/neighbor,
zero-sibling and resource gates. Do not rerun the consumed explicit-backend
diagnostic or pool its samples.

Fresh full Leopard1 parity checks for the actual AUTO candidate and a clearly
identified independently linked comparison remain required. Production
enablement also needs appropriate existing-test/default-assumption and API
documentation updates, plus clean build/integration checks. This object-level
checkpoint is not a clean full-CMake build, broad v19 closure, performance
qualification, or a claim that the production default is faster already.

Raw preparation is retained at `/tmp/leopard-auto-gfni-boundary.0JZpcL` and
references the unchanged original bundle
`.research/leopard-79h/gf16-current-route-failed.STc10h`.

The local immutable checkpoint is
`.research/leopard-79h/auto-gfni-boundary-checks.bJ6wiH`; it includes both
archives, object manifests, source snapshots, check executables, all 60 raw
checks, failed/successful preparation logs and the independent audit. No
second-host copy is made under the user's local-only instruction.

The sealed checkpoint has 218 manifest entries and 74,147,297 bytes. Its
outer `SHA256SUMS` hash is
`ae09b1b76a0dd877dcf0acff648df60c90cfdbe4442fbfbb0dee2f2bb71da512`.
The sealed-copy audit passed. Retention peaked at 81,117,184 bytes / 256 MiB
with zero memory events/swap. This manifest metadata supplements the sealed
report snapshot; the immutable bundle is not rewritten.
