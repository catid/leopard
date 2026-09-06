# Clock-free GF16 callback attribution

Bead: `leopard-79h.38.5.4.12`. Date: 2026-09-06.
Result: work attribution selects a new dispatch hypothesis; no timing claim.

The source-staging screen was flat despite removing 65,536,000 explicit copy
bytes. This observer instead counts the current encoder's actual backend
callbacks. It links the unchanged portable dual-field production archive at
`36dc0c8f66604b8d974468e687c51e6183ecb61d` and delegates the unchanged
`current_route_screen.cpp` check workload. It accepts no timing mode.

## Observation boundary and restrictions

The linker wraps only `ReedSolomonEncodeWithSourcePolicy`. For that call, a
stack-private copy of `Ops` substitutes sixteen counting callbacks, each of
which delegates the exact arguments to the original function pointer. Null
callbacks stay null; kind, name and all unrelated table bytes remain unchanged.
The original published table is checked for mutation after each call. Source
policy, input/output pointers, tile sizes and all arithmetic are preserved.
The driver forces one OpenMP thread locally, rejects reentry, and bounds the
record to sixteen passes, 256 histogram buckets and one million callbacks.

This technique is qualified only for these frozen portable archives. Their
actual Release and sanitizer FF16 compile commands both define
`LEO2_DISABLE_AVX2_CODEGEN=1` and `LEO2_DISABLE_SSSE3_CODEGEN=1`. Consequently
the address-based legacy in-field SIMD bypasses are compiled out. Substituting
an Ops address in a legacy/native in-field SIMD build would not establish the
same execution path and is not supported by this evidence.

Counts stop at the Ops boundary: range-internal helper calls are not counted
again. Libc copies, API validation, initialization, allocation and setup are
outside the observation scope. Range distance denotes lane groups, not bytes
transferred. A zero-skew mask suppresses multiplication, not XOR. Neither call
counts nor logical butterfly edges are time shares or measured memory traffic.

## Findings

For K1000/R200/64 KiB AUTO/GFNI, the public encode has two 32-KiB passes:

| Callback | Calls | Lane groups |
| --- | ---: | ---: |
| Inverse two-way | 2,768 | 2,768 |
| Inverse two-way with accumulation | 768 | 768 |
| Inverse four-way range | 160 | 1,144 |
| Forward four-way range | 36 | 360 |
| Forward two-way | 394 | 394 |
| XOR | 6 | 6 |
| Total callbacks | 4,132 | — |

The first inverse stage accounts for 2,000 ordinary two-way calls. The final
accumulation accounts for another 768 ordinary and 768 accumulating two-way
calls. Intermediate inverse and forward ranges cover 1,504 four-way lane
groups. These plus the first/final inverse and forward leaves total 9,952
logical two-way edges; that is a structural count, not a speedup bound.

The field's `IFFT_DIT4_Range` immediately diverts distance-one groups through
`IFFT_DIT4`, which selects the split two-way path at these byte sizes. Other
groups reach `AVX2FF16Butterfly4Range`: its already-shipped GFNI implementation
uses the fused four-way kernel at every size, ignoring the field's
`prefer_fused=false` hint. The first-stage early return predates the GFNI
backend. No new broad-fusion kernel or table-packing implementation is needed
to test this dispatch distinction.

Selected follow-up: let only the target GFNI encoder's first inverse
distance-one groups reach that existing in-place range kernel. Preserve the
copy-first schedule, tiling, terminal accumulation, forward transform, and
other routes. This differs from the rejected out-of-place source-staging
candidate; nevertheless, it is only a hypothesis until separately measured.

## Validation

- Twelve Release/sanitizer records for six fixed routes match the original
  public workload records; all six full Release parity files match standalone
  exact Leopard1, totaling 61,014,016 bytes. Release and sanitizer callback
  records are byte-identical for every route.
- An independent structural traversal model derives every bucket, distance,
  zero-skew mask and hint for the six shapes, including the odd-layer R512
  transform. The retained-only replay executes no collector or codec.
- Native unit tests check all sixteen exact delegations, unrelated metadata,
  null preservation, skew masks, both hints, aggregation, reentry, reset and
  atomic counter/bucket overflow refusal. Release and ASan+UBSan+LSan pass.
- Four adversarial model tests and the full replay pass with normal Python
  and `-O`. Twelve malformed/timing requests fail before workload execution.
- Build peaks are 90,324,992 bytes for the probe and 88,694,784 for the unit
  binaries under 512 MiB. Fourteen successful native scopes peak at
  141,656,064 bytes under 256 MiB; all six memory-event counters and swap are
  zero. No field, sanitizer or workload check was disabled.

Review is Codex self-review plus deterministic/adversarial checks under the
user's Claude opt-out, not an independent-model `CONVERGED` claim. Hardware
perf remains unavailable; no kernel settings or other workloads were changed.
The exhausted current-versus-Leopard1 timing attempt is not rerun or relabeled.
The parent performance gap and v19 qualification remain open.

## Retention and replay

The read-only bundle is `.research/leopard-79h/gf16-callback-probe.U9S0pT`,
retained locally and on ripper. It includes frozen observer/unit binaries and
sources, model/replay dependencies, build/check recipes, raw records/parity,
resource logs, and the original archive pins. `SHA256SUMS` authenticates all
bundle files. It references the separately retained complete source/build
bundle `.research/leopard-79h/gf16-current-route-failed.STc10h` (outer manifest
`e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`).

From the repository root, replay using the frozen scripts:

```bash
python3 .research/leopard-79h/gf16-callback-probe.U9S0pT/replay/verify_gf16_callback_probe.py \
  .research/leopard-79h/gf16-callback-probe.U9S0pT \
  .research/leopard-79h/gf16-current-route-failed.STc10h/preflight
```

Use the campaign's bounded 256-MiB/no-swap scope for execution and inspect its
resource result. The old absolute scratch paths in the retained build recipes
document actual commands; rebuilding elsewhere requires fresh output paths,
the canonical build lock, the same frozen codec inputs and new artifact hashes.
