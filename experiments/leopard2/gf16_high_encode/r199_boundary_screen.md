# Remaining AUTO R199 / 32 KiB boundary

Tracking: `leopard-79h.38.5.4.19`, in progress. Date: 2026-09-09.
The new comparison frontend is correctness-qualified, **not timed**.
No production source, route policy, compiler flags or codec archive changed.

Current `UseAutoGF16GFNIEncode` admits K=1000/R=200 at 32 KiB and
K=1000/R=199 or 200 at 64 KiB on this model-08 host. It deliberately leaves
R=199 at 32 KiB on AVX2. Nearby GFNI wins motivate measuring this remaining
case; they do not prove a current Leopard1 deficit or transferable speedup.
The [qualified AVX2 prototype](avx2_adjacent_schedule.md) remains separate.

## Actual comparison inputs

[r199_boundary_screen.cpp](r199_boundary_screen.cpp) uses one frontend,
separately linked to unchanged original native Leopard1 and current Leopard2.
K=1000, R=199, bytes=32768, GF16 legacy-high full parity, one thread.
Original native Leopard1 remains the product comparator, including its native
ISA/tuning; this is not a pure-AVX2 or single-instruction comparison.

| Archive | SHA-256 |
| --- | --- |
| Native Leopard1, source `6e5725e` | `3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1` |
| Production Leopard2 Release, source `3a2f064` | `d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a` |
| Production Leopard2 full sanitizer | `e2d1bacb1142f71f53ab8aaf927dc0254c37fa9abfdea89727633faf634724a6` |

The normal Release frontend is `d2edecdccc9d366fdd03d0412bb51a8439fe1aef15ed2ace9f34c0513d623490`;
native frontend is `ea73abe01e1903ae0bf2fcd756ed978d0697095c798be07d282e894ae949747c`.
Both are built but have only executed their single-check, untimed mode.
Their separately linked clock-guard variants also exercise the complete
1 check + 4 warmup + 21 sample-call schedule without reading a benchmark clock.
Five deliberate clock attempts terminate at the guard with status 86.

Only the public encode call is inside the future sample interval. Initialization,
allocation, input generation, scratch query, hashes and parity dumps are outside.
Native Leopard1 returns parity in its first R work buffers; Leopard2 uses
separate parity outputs. Both conventions remain intact. Scratch sizes are
16,777,216 and 16,808,512 bytes respectively; these are not runtime costs.

## Correctness evidence

49 positive records, 82 total; 28 malformed/timing CLI refusals and five clock
guards. All nine full parity comparisons, totaling 58,687,488 bytes, pass:
the fresh native result matches the prior pinned native reference, and both
AUTO and GFNI match it in Release and full ASan/UBSan/LSan builds, including
the separate driver-only callback observers.

Each Leopard2 profile passes eight new guarded shapes: the exact target,
unaligned target, 32770/32766-byte tails, small 64/66-byte GF16 payloads,
and GF8/GF16 small-shape controls. Each covers six full/subset/empty output
masks, input preservation, poisoned or checked guards, scratch geometry,
short-scratch rejection and odd-GF16 rejection. The shared experimental test
was factored into `CheckShape`; all eight original shapes are rerun in each
profile and their records agree between Release and sanitizer builds.

The actual production observers confirm one pass with backend kind 3 for
AUTO and 6 for explicit GFNI. Both have 2,066 callbacks and match the
independent structural traversal model. The original FF16 archives compile
out the legacy address-based SIMD paths, so the private observing Ops table
does not select another algorithm. Counting callbacks is not timing attribution.
The timing-capable executables contain neither observer nor clock wrapper.

The [collector-free replay](verify_r199_boundary_checks.py) checks the exact
archive/source/header/driver pins, raw record inventory, CLI outcomes, actual
routes, all parity bytes, native reference manifest and resource envelopes.
Normal and optimized replays agree. Eight pure replay/model tests pass in
each Python mode, including strict-type, route, scratch, source-identity,
call-count and accidental-timing mutations.

The serial driver build peaked at 155,054,080/512 MiB; native qualification
at 180,563,968/256 MiB. All six memory-event counters and swap were zero.
[Numerical qualification result](results/r199_boundary_checks_20260909.json).
Raw workspace: `/tmp/leopard-r199-boundary.v13Wsv`.
All work is local and serial under the canonical campaign lock. Claude and
subworkers remain opted out; this is Codex self-review plus deterministic
and adversarial evidence, not independent-model `CONVERGED`.

The local read-only bundle
`.research/leopard-79h/r199-boundary-qualified.s3qjbngr` retains 230 files,
142,110,648 bytes, outer manifest SHA-256
`a509866ab0dbf2967a8a25258816134927f4b29ac8ce50c7bb3c2ce54de2fad3`.
All manifest entries and read-only permissions pass, and the sealed-copy
semantic replay reproduces the tracked result and both earlier replay outputs.
Retention peaked at 36,192,256/256 MiB with all memory events and swap zero.
Separate delivery checks and retention logs:
`/tmp/leopard-r199-delivery.awlZwh`. Replay retains explicitly pinned local
dependencies on the original native reference and codec archives.

Delivery removed one extra trailing blank line from the tracked frontend.
Recompiling that normalized source at the original logical source path
reproduced **all three complete frontend objects byte-for-byte**, including
the sanitizer object. Its source SHA is now
`cce006d06c014de738ad497b4cad6c3f12c258a5dda822c03e02af2d0f84ebe5`;
the retained original source SHA is
`dff959176e2cc705468bfa59e0669c8c4d537a5e72409512d1ec3ec1eaf49a9b`.
Two earlier cross-directory reproduction checks matched the unsanitized
objects but not sanitizer source-path metadata; both logs are retained.
The accepted same-path check peaked at 136,261,632/512 MiB, all events/swap
zero. The owned temporary source copy was restored to its exact original
bytes; no codec archive or previously qualified executable was modified.

## Remaining performance gate

Before clocks, qualify the collector/protocol and commit **and push** a fresh
immutable one-attempt preregistration. Measure AUTO/native, GFNI/AUTO and
GFNI/native directly, plus all three same-path controls; three rounds of
21 samples with balanced ABBA ordering. Local CPU 26/sibling 90, zero sibling
work, 256 MiB/no swap, unchanged 2% control bound and 5% candidate threshold.
Do not retry or pool any older attempt, alter unrelated thread affinity, or
infer a result from the neighboring cases. A positive initial screen requires
separate API/boundary/neighbor and same-binary integration qualification before
any AUTO policy extension. The full performance goal remains open.
