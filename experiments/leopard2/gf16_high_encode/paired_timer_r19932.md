# Paired timer adapter qualification

Bead: `leopard-79h.38.5.4.19.1.1`.

This is untimed measurement-method qualification, not a performance result.
The R199/32 KiB AUTO candidate stays default-OFF; the consumed `aa2b034`
integration attempt stays inconclusive and must not be repeated or pooled.
Its failed native/AVX2 whole-process controls remain causally unexplained.

## Exact cost boundary

The separate `paired_timer_r19932.cpp` preserves the preceding
`paired_r19932.cpp` prototype, its fixed codec/buffers, allocation-edge guards,
four-state schedules, route probes, full parity comparisons and call counts.
`PairedGroupTiming.h` supplies the actual shared grouped-call loop:

1. State selection and route inspection happen outside the duration.
2. Start clock, N complete public calls, end clock.
3. Duration validation, normalization, sample storage, slot accounting and
   full parity/guard checks happen outside the duration.

The public result checks, per-call counter increment and repetition-loop
overhead remain **inside** the duration. There is no timer-overhead subtraction.
The one-item batch API still receives one item on each public call. The group
count is 1, or 256 for the tiny GF8 cell only. A duration divided by 256 is a
grouped per-call average, not single-call latency.

Four route-probed preflight calls precede four warmup schedule passes and
21 sample schedule passes. A sampling invocation therefore contains 84 spans
and 168 clock calls, with exactly 104 or 25,604 total public encode operations.
Full output comparisons remain after **every** group, as in the preceding
prototype. Although excluded from the measured span, these checks touch output
memory and can affect subsequent cache state. A future timing protocol must
disclose this placement and use the same frontend; this is not an attribution
experiment isolating timer granularity from all other effects.

All three clock variants consume the same compiled driver object. `plain`
links the real steady-clock symbol, `synthetic` replaces it with deterministic
ticks and a public-entry witness, and `abort` exits before reading that clock.
Qualification executes `plain --exercise`, `synthetic --clock-exercise` and
`abort --clock-guard`. It never invokes `--measure`, even as a negative CLI
test. Clock-kind checks prevent synthetic/abort output from being presented as
real timing. The unwrapped executable refuses `--clock-guard` before allocation
or any clock call.

The synthetic clock records actual public-entry counts at each boundary.
Adjacent endpoints prove that precisely N operations lie inside every span
and none occur between the end of one sampled group and the start of the next.
Expected durations are deliberately nonmultiples of 256, checking fractional
normalization without integer truncation. Nonpositive, reversed, negative and
overlarge intervals are rejected before unsafe signed subtraction. Durations
are bounded by `2^53-1`; conversion to binary64 and division by 1 or 256 are
exact over this accepted range.

The synthetic/public-entry wrappers add test overhead and are not benchmark
binaries. Native Leopard1 remains the original separately linked product
comparator, not the pure-AVX2 attribution build. Within-process pairing cannot
by itself establish cross-process stability for that product comparison.

## Validation status

Untimed qualification passed:

- 182 positive invocations: 180 all-nine-cell native/Release/full Leopard2
  ASan+UBSan+LSan frontend exercises and two grouped-timer unit executables.
- 37 C++ unit cases per Release/sanitizer profile; ten Python protocol and
  adversarial tests pass normally and with `-O`.
- 7,560 synthetic spans, each bracketed around the exact expected number of
  actual public-entry calls; 18 abort guards, 12 malformed-clock rejections
  and 42 CLI refusals, including the new plain `--clock-guard` refusal.
- 170 complete native parity comparisons totaling **1,034,468,224 bytes**.
- Collector-free raw and read-only-copy replays pass in both Python modes,
  checking archive/source/driver identities, shared-object link recipes,
  record inventory, call/order/clock boundaries, normalization, full parity
  and resource logs. The four replay results are byte-identical.
- Allocation-guard implementation is byte-identical to the earlier qualified
  prototype (tested explicitly), whose actual canary/ASan negative checks
  remain preserved. The new frontend's per-group guards pass for all cells.

All three codec archives, `leopard2.cpp` and `Leopard2Direct.h` are unchanged
from the earlier candidate qualification. No benchmark clock was read.

| Scope | Peak bytes | Limit | `memory.events max` |
| --- | ---: | --- | ---: |
| Corrected driver build | 158,175,232 | 512 MiB | 0 |
| Corrected native checks | 183,197,696 | 256 MiB | 0 |
| Raw replay/tests | 268,435,456 | 256 MiB | 7,385 |
| Retention/sealed replay | 268,435,456 | 256 MiB | 24,599 |

These four scopes exited zero with no OOM, OOM kill, swap or other memory
events. Replay/retention reached their limit: this is **not** an all-zero
resource claim, and the limit was not increased. The qualification collector
uses `fsync` and `POSIX_FADV_DONTNEED` on its own closed-process parity evidence
files; it does not change host settings, other processes or codec buffers.

Read-only evidence: `.research/leopard-79h/paired-timer-qualified.D601n0`.
Its 751-entry `SHA256SUMS` has digest
`6c86e08b8fcef3a7d3a9c345b5820bf9e8af1310133b6c8835cd4adcc2fab37b`.
Raw evidence is `/tmp/leopard-paired-timer-fixed.nU5rLs`; its terminal
`retention.log` supersedes the initial snapshot copied into the bundle.
Replays still rehash recorded original local input paths; this is not a
standalone/hermetic bundle. Result and record pins are in
`results/paired_timer_r19932_checks_20260909.json`.

Self-review caught an initial CLI provenance bug: the unwrapped executable
accepted `--clock-guard`, which could have produced real durations labeled
untimed. No such combination was executed. The initial qualification was
deliberately stopped by terminating its verified owned child; its collector
exited 1, peak 155,901,952 bytes under 256 MiB, all memory events/swap zero.
The initial build and partial logs remain at
`/tmp/leopard-paired-timer.eRs5C1`. They are superseded, not qualification
evidence. The corrected run is separate at
`/tmp/leopard-paired-timer-fixed.nU5rLs`; no old files or consumed timing
attempts were overwritten.

Review provenance is Codex self-review plus deterministic/adversarial checks.
Claude is explicitly opted out; this is not independent-model `CONVERGED`.
Work is local only, serial under the canonical lock, with 512 MiB builds,
256 MiB checks, per-child CPU limits, and no swap. No host settings, unrelated
processes, remote workers, production codec or previous frontends changed.

## Next performance gate

`leopard-79h.38.5.4.19.1.2` covers the separate explicit method review and new
preregistration. It depends on this untimed qualification. The new protocol
must retain both target APIs, all seven neighbors, the original native
comparator, and controls appropriate to every within- and cross-process
comparison. Preserve 5% target, 2% control/neighbor and zero-sibling gates,
the fixed attempt budget and immutable executable provenance. Commit and push
the preregistration before any benchmark clock reads. Promotion additionally
requires a valid performance gate and validation of the actual default-on
artifact. The broader Leopard1/Leopard2 improvement goal remains open.
