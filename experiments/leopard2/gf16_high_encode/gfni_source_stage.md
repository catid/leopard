# GFNI-only source-staging contrast

Bead: `leopard-79h.38.5.4.11`. Date: 2026-09-06.
Status: **correctness validated; subsequent timing filter rejected promotion**.

The untimed evidence below remains valid. The later separately preregistered
filter completed 144 clean timed invocations and found an essentially flat
0.998812x current/staged ratio, with passing aggregate controls. No production
change was promoted. See `gfni_source_stage_screen.md` for that distinct result.

## Why this is a distinct experiment

The previous copy probe showed that current AUTO/GFNI copies 65,536,000 input
bytes before its first inverse stage at K1000/R200/64 KiB. Output already binds
directly to the caller's buffers, so this is not a final-scatter optimization.
The 2026-07-22 rejected source-staging experiment used AVX-512 nibble-table
arithmetic, not GFNI. Its retained ratio was 0.998799, directional rather than
authoritative. That negative result remains rejected and is not pooled here.

Commit `2942f35` extended AVX2-family policy predicates to the new GFNI identity,
including the greater-than-16-KiB copy-first crossover. Report 19 explains why:
the earlier experimental GFNI member had run under the AVX2 identity. It does
not provide a GFNI-specific staging crossover comparison. Current GFNI owns
its affine multiplication tables and out-of-place first-stage kernel. This
motivates a new contrast; it does not prove the inherited threshold is wrong.

## Isolated mechanism

No production source, archive, route, or CMake configuration changes. GNU ld
wraps the existing cross-object `ReedSolomonEncodeWithSourcePolicy` call. The
exact mangled name is checked against the pinned archive before building.
The wrapper changes only the source-policy argument from 65,536 to 16,384 when
all of these internal pass properties match:

- GFNI operation table;
- K1000, recovery prefix 200, requested output count 200, padded side 256;
- no sparse schedule blocks;
- execution pass 32,768 bytes and original source policy 65,536 bytes.

The two 32-KiB execution passes, work pointers, byte count, caller validation,
arithmetic table and wire format are unchanged. The existing encoder then
loads sources through `ff16_ifft_butterfly4_out` instead of copying them and
running the first stage in place. Mode 0 forwards the original policy; mode 1
activates this diagnostic predicate. The wrapper is default-off and exists
only in separately linked experiment executables.

This is an internal pass predicate, **not a proposed production AUTO selector**.
For example, explicit GFNI with 65,538 bytes has the same two aligned passes
plus an unchanged compact tail, so those aligned passes also change in the
directed test. Current AUTO at that size stays AVX2 and does not change.

## Observed result

At the exact AUTO/GFNI target:

| Untimed observation | Control | Candidate |
| --- | ---: | ---: |
| Execution passes / bytes each | 2 / 32,768 | 2 / 32,768 |
| Explicit input-copy calls | 2,000 | 0 |
| Explicit input-copy bytes | 65,536,000 | 0 |
| Zero-fill calls / total bytes | 48 / 1,572,864 | 48 / 1,572,864 |
| Full parity output | Exact Leopard1 match | Exact Leopard1 match |

All five other fixed cells retain their original copy counts and policy
arguments, including explicit AVX2/AVX-512, AUTO 32-KiB, R199, and K4096/R512/4-KiB.
Both Release modes and both sanitizer modes match all six original workload
records: 24 records total. The 12 Release parity files compare every byte with
the independently linked Leopard1 oracle, totaling 122,028,032 compared bytes.
The source-policy trace proves that only the two intended target passes change.

These are explicit libc-call counts, not total memory traffic or timings.
The fused kernel still reads the inputs and performs arithmetic. The copy
instrumentation adds per-call overhead and **must not be used for timing**.

## Correctness and failure retention

`test_gfni_source_stage.cpp` passed in Release and ASan+UBSan+LSan:

- 224 first-stage kernel comparisons against scalar arithmetic, covering all
  eight zero-skew combinations, fourteen lengths from 0 through 65,598 bytes,
  zero/17-byte offsets and exact allocation ends;
- 16 one-field predicate boundary negatives plus the exact positive;
- 14 public cases: exact AUTO and explicit GFNI, partial-output masks, K/R
  neighbors, explicit AVX-512, even compact tails around 64 KiB, and native-odd
  rejection. Twelve cases encode in both modes; two odd cases correctly reject
  query and execution before any transform. Outputs, source hashes, prefix
  canaries and omitted-output contents match;
- Too-small scratch and input/output aliasing reject before the wrapper is
  reached and leave source/output bytes unchanged.

An initial fixture incorrectly expected native GF16 to accept 65,535 bytes.
Its scratch query returned `LEO2_UNSUPPORTED`, consistent with the documented
complete-symbol native layout. That failed test/source/binary is retained.
The final fixture preserves both odd cases as rejection checks and adds the
even 65,534/65,538 AUTO neighbors. No codec or candidate change was needed.

The four fixed-workload binaries also reject 12 malformed/timing requests.
The separate stdlib replay checks all 24 workload/copy/policy records, full
parity files, matching Release/sanitizer directed results and 54 successful
resource envelopes. Normal and optimized Python both pass. Every native job
uses the canonical lock, one substantial process, and no swap. Maximum test
scope peak was 141,651,968 bytes under 256 MiB; the public-probe build peaked at
111,902,720 under 512 MiB, and final directed-test build at 88,297,472. All six
memory-event counters and swap were zero. No coverage or memory cap was removed.

Review is Codex self-review plus deterministic and sanitizer checks under the
user's Claude opt-out, not independent-model `CONVERGED`. The unavailable
`deli-auto-research` skill is not represented as an automatic watchdog.

## Evidence and next gate

Current codec source remains `36dc0c8`; exact Leopard1 remains `6e5725eb`.
All initial artifact/source hashes match again after validation. The 174-file
read-only bundle is retained locally and on ripper at
`.research/leopard-79h/gfni-source-stage.NcBQms`, outer `SHA256SUMS`:
`b4dd94d114c439f5a1b597fc3f8f7413c8a39ea8de198d5b84be6424b198b542`.
It retains initial failures, final sources/binaries, recipes and every raw
record/log. Its two unchanged reference bundles are stored alongside it:
`gf16-current-route-failed.STc10h` and `gf16-copy-probe.zE4w1M`.

At this correctness milestone, the next gate was a separately preregistered
candidate-versus-current timing filter,
using fresh immutable binaries without the libc-copy wrappers. This would
test this new implementation, not rerun the exhausted current-versus-Leopard1
plan. Retain a fixed attempt budget, same-binary controls and passive sibling
gate; stop on contamination without changing other workloads. Any affinity
intervention still requires user approval. Existing v19 qualification and
authoritative exact-Leopard1 gap-closure gates remain open and unchanged.
That subsequent filter is now complete and negative. The candidate remains
experiment-only; the next investigation is callback-cost attribution, not a
retry of this rejected policy replacement.
