# Untimed GF16 public-API copy attribution

Bead: `leopard-79h.38.5.4.10.1`. Date: 2026-09-06.
Result: copy-path attribution completed; **no speed or regression claim**.

The current-route timing plan stopped at its passive gate and remains
exhausted. This follow-up executes only `--check`, using the same six complete
public-API workloads and unchanged codec archives. It does not retry that
plan, relax a timing rule, or change any production route or workload affinity.

## Method and validation

`current_route_copy_probe.cpp` includes the original workload driver and is
compiled independently against current Leopard2 and standalone Leopard1.
GNU ld wraps each public encode entry plus `memcpy`, `memset`, `__memcpy_chk`
and `__memset_chk`. Counters are active only during one public encode. The
probe forces its own OpenMP thread count to one and dynamic teams off; it does
not alter other processes. A linker self-check covers all four libc symbols
before each workload. `--measure` and out-of-range cells are rejected.

Only external symbol calls from the linked objects are observed. Inlined
copies, vector arithmetic loads/stores, memory operations internal to shared
libraries, cache traffic and latency are outside this measurement. Source
inspection is required to interpret an absence of copy calls. The probe
binaries are instrumented diagnostics, not timing artifacts.

Codec source identities and archive hashes are those in
`current_route_screen.md`: Leopard2 `36dc0c8`, Leopard1 `6e5725eb`. Driver-only
`-fno-builtin-memcpy/-memset` flags do not recompile either codec archive.
Final source and artifact hashes were verified again after all checks.

All12 Release output files compare byte-for-byte against the retained exact
Leopard1 parity, totaling122,028,032 compared bytes. All18 workload JSON records
match their original expectations. The six ASan+UBSan+LSan records and copy
counts match current Release exactly. Normal and optimized Python replays
validate all18 records, parity bytes and zero-event resource envelopes, and
reject16 altered count records each. Nine invalid/timing invocations are
rejected; an inherited64-thread/dynamic OpenMP configuration is safely
overridden and produces the same cell5 result and counts.

## Observed explicit input copies

| Cell | K/R/shard bytes | Current route | Leopard1 calls / bytes | Leopard2 calls / bytes |
| --- | --- | --- | ---: | ---: |
| 0 | 1000/200/65536 | AUTO GFNI | 1000 / 65,536,000 | 2000 / 65,536,000 |
| 1 | 1000/200/65536 | AVX2 | 1000 / 65,536,000 | 2000 / 65,536,000 |
| 2 | 1000/200/65536 | AVX-512 | 1000 / 65,536,000 | 1000 / 65,536,000 |
| 3 | 1000/200/32768 | AUTO AVX2 | 1000 / 32,768,000 | 1000 / 32,768,000 |
| 4 | 1000/199/65536 | AUTO AVX2 | 1000 / 65,536,000 | 2000 / 65,536,000 |
| 5 | 4096/512/4096 | AUTO AVX2 | 4096 / 16,777,216 | 0 / 0 |

All observed copy sources lie in the original input slab; there are no
other-source copy calls or fortified calls during these encodes. Calls whose
destination lies in the parity slab are a subset of the initial input copies,
**not an additional final scatter**. In `leopard2.cpp`, the aligned single-pass
and tiled high-profile paths bind `work[i]` directly to the caller's recovery
buffer or its tile offset. There is no separate final parity-copy pass here.

The K1000 partial final message block zeroes24 unused rows in both codecs.
Tiled current paths split those24 clears into48 calls, preserving total zeroed
bytes:1,572,864 at64KiB, or786,432 at32KiB. Cell5 has no zeroing calls.

`LeopardFF16.cpp::IFFT_DIT_Encoder_Impl` explains the input-copy difference:
for large high-profile AVX2/AVX-512/GFNI transforms it retains copy-first when
the public source-policy size exceeds16KiB. The two32-KiB execution passes
therefore preserve the64-KiB copy-first policy. Cell5 instead consumes sources
through the fused out-of-place first inverse stage, avoiding the explicit
16-MiB copy pass. This does not mean the input is never read.

## Consequence for the next experiment

Do not optimize a nonexistent final output copy, and do not remove byte
tiling merely because it doubles call count: its prior cache/scratch benefit
does not disappear, and these counts measure no call overhead.

The GFNI copy-first gate is a specific candidate for investigation. The
source comment attributes its crossover to AVX2, while GFNI owns a distinct
affine first-stage kernel; commit2942f35 added GFNI to the shared gate. This
does not prove the gate is wrong. Follow-up `leopard-79h.38.5.4.11` first audits
prior evidence, then tests a GFNI-only, isolated source-staging contrast with
the same byte tiling. A valid controlled-host comparison is still required
before any performance conclusion or production promotion. The existing v19
and K65 gates, failed-attempt budget and host-approval boundary are unchanged.

## Resource and failure retention

The initial self-check failed before codec execution because fortified
fixed-size libc operations were inlined to vector stores, even with the
driver's no-builtin flags. Retained disassembly proves this. Volatile indirect
symbol calls fixed the self-check; fortified symbol wrappers were also added.
No codec source changed.

The second, combined validation scope passed at268,374,016/268,435,456 bytes,
all events and swap zero, but had insufficient headroom for reuse. Final
validation runs each shape in its own256-MiB/no-swap scope under the canonical
lock: maximum141,275,136 bytes. The final serial build peaked at146,530,304
under512MiB; CLI/thread-policy guards peaked at35,786,752 under256MiB. All final
resource counters and swap are zero. No field, shape or sanitizer was removed.

The171-file read-only bundle, including the failed self-check, combined run,
final source/binaries, all output/parity/resource logs, recipes and replay, is
retained locally and on ripper at
`.research/leopard-79h/gf16-copy-probe.zE4w1M`. Its outer `SHA256SUMS` is
`b9c8751c14e138c7dd3645e478ccfe2cb1149ff134665d7be548445164fd737d`.
It references the previously sealed158-file bundle
`gf16-current-route-failed.STc10h` for unchanged codec archives/source and the
original parity oracle; that bundle is retained alongside it on both hosts.

Review is Codex self-review plus deterministic/sanitizer checks under the
user's Claude opt-out, not independent-model `CONVERGED`. The unavailable
`deli-auto-research` skill is not represented as an automatic watchdog.
