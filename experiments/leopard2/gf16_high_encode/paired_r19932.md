# Untimed paired AUTO R199/32 KiB frontend

Bead: `leopard-79h.38.5.4.19.1.1` (still open).

The new **untimed schedule prototype** passes Release and full Leopard2
ASan/UBSan/LSan checks. It does not establish a speedup or fix the cause of
the failed controls in `auto_r19932_screen.md`. That attempt remains consumed,
and the codec candidate remains default-OFF. OFF denotes the changed binary
with its diagnostic selector disabled, not pristine production.

## Implemented and checked

`paired_r19932.cpp` creates one codec and one fixed set of input, output and
scratch buffers per process. It executes `0110`, `1001`, `0000` or `1111`
schedules, changing only the R199/32 KiB diagnostic flag while execution and
inspection are quiescent. The separately linked, unchanged original native
Leopard1 comparator uses `NNNN`; it is not the pure-AVX2 attribution variant.

Each slot first performs one route-probed encode. The probe is normalized
before exercise. Four warmup passes and 21 exercise passes then execute the
entire four-slot schedule in sample-major order. Groups contain either one
complete public API call or, for the tiny GF8 cell only, 256 complete public
calls. Thus exercise performs 104 or 25,604 public operations per process,
including four preflight calls. The one-item batch case still calls the public
batch API with exactly **one** item on every operation. Grouping is not a
multi-item batch, and any future group duration divided by its operation count
would be a grouped per-call average, not single-call latency. Neither grouping
factor has been chosen for a new timing protocol.

The external test-only `paired_public_witness.cpp` checks actual public-entry
call order, API and state counts, and unchanged codec, pointer arrays, shard
pointers and scratch identity. The plain and witnessed executables share the
same compiled driver object. Allocation-edge canaries and ASan poisoning guard
the contiguous buffers; these are not per-shard redzones. Every exercised group
is compared byte-for-byte with the process's preflight parity. Source integrity
is checked at the end. Final parity for every witnessed schedule is retained
and compared in full with the native comparator.

Coverage preserves all nine cells: the ordinary and one-item-batch R199/32 KiB
targets, the three already-enabled GFNI neighbors, R198/32 KiB, explicit AVX2,
K4096/R512/4 KiB, and GF8 K17/R7/64 B. Both mixed-state orders and both same-path
controls pass. No codec, header, ISA flags, routing baseline or previous
single-call frontend was changed.

## Evidence

- 180 positive frontend invocations; 18 clock aborts; 39 CLI refusals.
- 80 complete Leopard2/native parity comparisons: **486,808,576 bytes**.
- Actual canary corruption rejected in Release and sanitizer builds; actual
  underflow/overflow reads rejected by ASan as `use-after-poison`.
- Nine pure protocol/adversarial tests pass normally and with Python `-O`.
- Collector-free raw and read-only-copy replays pass in both modes, producing
  identical results. They validate the frozen archives, driver artifacts,
  records, schedule/order accounting, full parity and resource logs.
- Both executable variants link a clock-abort wrapper. There is no `--measure`
  interface; `--clock-guard` aborts at the first exercise clock, after exactly
  four probes and four warmup schedule passes. No benchmark clock was read.

The unchanged archives retain SHA-256 identities:

| Profile | SHA-256 |
| --- | --- |
| Native Leopard1 | `3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1` |
| Leopard2 Release | `89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334` |
| Leopard2 sanitizer | `c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9` |

Resource accounting is **not all-zero**:

| Scope | Peak bytes | Limit | `memory.events max` |
| --- | ---: | --- | ---: |
| Driver build | 165,138,432 | 512 MiB | 0 |
| Frontend checks | 268,435,456 | 256 MiB | 1,519 |
| Guard build | 86,040,576 | 512 MiB | 0 |
| Guard checks | 21,446,656 | 256 MiB | 0 |
| Raw replay/tests | 268,435,456 | 256 MiB | 1,761 |
| Retention/sealed replay | 268,435,456 | 256 MiB | 10,966 |

Every scope exited zero, with no OOM, OOM kill or swap and no other memory
events. The cap was not increased. These observations do not isolate why the
limit events occurred. All substantial jobs were local, serial, and protected
by `/tmp/leopard-gf8-authoritative.lock`.

Raw evidence: `/tmp/leopard-paired-r19932.AqUmvR`.
Read-only copy: `.research/leopard-79h/paired-r19932-qualified.jAfhHN`.
Its 609-entry `SHA256SUMS` has digest
`590f01c47dc27dfe3d656a17a422e1a46bcebeea084a1480bf13a2873871bd08`.
The terminal retention-scope log is in the raw directory; the copy contains
only its initial snapshot. Replays also rehash the original pinned build
inputs at their recorded local paths; this is not a standalone/hermetic bundle.
Machine-readable outcome: `results/paired_r19932_checks_20260909.json`.

## Remaining gate

This is not yet a qualified timing executable. A real grouped start/end timer
adapter, its exact cost boundary and normalization still need clock-free
qualification (including deterministic synthetic-clock tests). The prototype
currently performs full parity/guard checks between groups; their placement
must be explicit in any timing adapter rather than silently changing the
workload. Plain grouped execution must also be exercised without the external
witness. The active Bead remains open for that work.

Any subsequent performance experiment requires a separate explicit review and
a committed **and pushed** preregistration. Preserve both native and same-path
controls, all nine cells, the 5% target / 2% control-neighbor / zero-sibling
gates and a fixed attempt budget. Do not rerun or pool the consumed attempt,
trim its samples, infer a speedup from its target ratios, or enable the
candidate without a valid performance gate and actual default-on artifact
verification. Native/AVX2 whole-process slowdowns remain causally unexplained.

Review provenance: Codex self-review and deterministic/adversarial checks.
Claude was explicitly opted out; no independent-model `CONVERGED` claim, remote
workers, unrelated process changes, host-setting changes or new timings.
