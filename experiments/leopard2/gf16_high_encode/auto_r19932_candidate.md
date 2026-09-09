# AUTO GFNI R=199 / 32 KiB qualification

Bead: `leopard-79h.38.5.4.19.1`. Status: correctness checkpoint,
**default off and not performance-qualified**.

The earlier [direct screen](r199_boundary_screen.md) measured native Leopard1
8.48% faster than production AUTO at K=1000/R=199/32 KiB. Explicit GFNI was
53.07% faster than AUTO and 42.28% faster than native Leopard1. Those results
motivate this integration; they do not establish its speed or authorize default
promotion. That experiment is exhausted and will not be retried or pooled.

## Implementation

A separate internal atomic control adds only K=1000/R=199/32 KiB to the
existing GFNI selector. It starts disabled. The three previously qualified
cells (R=200/32 KiB, R=199/64 KiB, R=200/64 KiB) retain their old controls.
Existing cached-table, actual host family/model, one-thread, AUTO,
GF16, legacy-high profile, native-layout, flags, and full-output guards remain.
Codec setup already considered R=199 for its existing 64-KiB route, so no
additional setup eligibility or qualification request was introduced.

The new control must be changed only while codec execution and inspection are
quiescent. A codec without a qualified cached table cannot acquire one through
late enabling. Ordinary encode and one-item batch are the only added APIs;
partial output, scalable/reusable/multi-item batch, decode, explicit backends,
other field/layout/profile/flags/thread/shape/size cases remain excluded.

Default OFF is a state of the changed binary, **not pristine production binary
identity**. Integration timing must compare OFF and ON within this same binary.

## Deterministic evidence

The retained clock-free qualification has 223 positive native invocations,
259 total records, 117 public records, five clock aborts before any benchmark
clock read, and 31 malformed/unsupported CLI refusals. It compares 44 complete
parity buffers, totaling 297,765,056 bytes, against independently linked native
Leopard1. Eight fresh native cells also match earlier retained native parity;
R=198/32 KiB has a newly produced native comparator.

Both Release and full ASan/UBSan/LSan exercise:

- Added/old route eligibility, cached-table and late-enable behavior, host and
  qualification wrappers, full output and one-item batch routing.
- Partial output, scalable and reusable bindings, two-item batch, invalid batch
  atomicity, decode roundtrip, shared immutable-codec concurrency across GF8/GF16.
- Eight new and eight established guarded shapes in both control states:
  six output masks, alignment offsets, even tails, short scratch and odd GF16
  rejection. Guard regions are ASan-poisoned in the sanitizer build.
- Nine public frontend cells in both candidate states, including both target
  APIs, the three old GFNI cells, R=198, explicit AVX2, K=4096/R=512/4 KiB, and
  GF8. Every clock-free exercise executes the intended 1+4+21 encode schedule.

The public timing frontend uses aligned allocations, not the separate focused
test's poisoned outer guards. Timing-capable native/Release executables ran
only `--check`; no integration benchmark samples have been collected.

Twelve additional real Release backend-hook checks pass. Nine cover KAT,
GF8-table allocation and GF16-table allocation faults at the two old boundary
cells and the new cell. They verify cached qualification status, no partial
table publication, one fault consumption, full fallback parity against explicit
AVX2, and subsequent codec setup. Three cover disabled/ineligible inert setup
and the established production path. These are actual backend fault hooks,
not the focused test's qualification wrappers. They are Release-only; the
hookless candidate's focused tests have full sanitizer coverage.

Three new CTest registrations were configured and executed successfully with
the pinned manually linked hook executable. The fresh CMake core recipe matches
the actual hook build recipe. This is not a complete fresh CMake library build.

The strict existing Release archive ISA scan passes. Only `leopard2.cpp.o`
changes in each 24-member production archive; the other 23 are byte-identical.
The separate hook archive has 23 members, with 22 unchanged. Its initial build
harness wrongly assumed 24; the failed record is retained, and the successful
core compilation was reused after checking the actual original inventory.

Seven pure verifier tests pass in normal and optimized Python modes, including
record/route/CLI/type/resource rejection cases. Collector-free raw and read-only
copy replays agree in both modes. They validate exact inventories, source/build/
artifact identities, real archive payloads, raw statuses and records, complete
parity bytes, route/call schedules, CTest results and resource records.
Review provenance is Codex self-review and deterministic/adversarial checks;
Claude is explicitly opted out. No independent-model `CONVERGED` claim.

## Artifacts and resource caveats

| Artifact | SHA-256 |
| --- | --- |
| Candidate Release archive | `89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334` |
| Candidate sanitizer archive | `c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9` |
| Original native Leopard1 archive | `3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1` |
| Real-hook candidate archive | `c876e7d394fe765515485a35f20d0628ae9d261ec277c7e445ac1482431ac34f` |

Raw: `/tmp/leopard-auto-r19932.9EGR0p`.
Read-only local bundle:
`.research/leopard-79h/auto-r19932-qualified.515m59j2`.
Its manifest covers 1,572 files / 422,122,956 bytes; manifest SHA-256:
`4799b38994df2804ea77df521200fb28e66ba4fb4f1a9b1bf600d4399b0a9f18`.
Retention/replay logs:
`/tmp/leopard-auto-r19932-delivery.ZBIs9I`.
[Machine-readable result](results/auto_r19932_checks_20260909.json).

All substantial work was serial on local `work`, under the canonical lock and
512-MiB build / 256-MiB check caps, with swap disabled. The initial native
qualification **reached 268,435,456 bytes and recorded 612 memory.max events**.
It exited successfully with no OOM, OOM kill or swap. This is not an all-zero
resource run and is retained as such. No cap was raised or test removed.

Candidate build peak: 274,698,240 bytes. Real-hook checks plus ISA:
78,254,080 bytes. Final raw replay/tests: 94,797,824 bytes. Retention plus
read-only replays: 227,409,920 bytes. Those scopes have all six memory events zero
and no swap. The failed hook-inventory scope also stayed below its cap.

Replay locally, without executing any codec:

```sh
python3 experiments/leopard2/gf16_high_encode/verify_auto_r19932_checks.py \
  .research/leopard-79h/auto-r19932-qualified.515m59j2
```

The bundle retains exact build/preparation scripts, commands, source copies,
archives/executables, raw output/parity, the initial hook-harness failure,
corrected continuation, CMake registration and validation logs. The replayer
also requires its pinned earlier native/archive references under `.research`.
Large raw artifacts remain local; source, result and replay code are in Git.

## Next gate and roadmap

The integration remains open. The fresh [same-binary integration
screen](auto_r19932_screen.md), preregistration `aa2b034`, completed all 372
timed invocations but failed three unchanged-path control aggregates. It is
inconclusive, consumed and may not be retried. No default promotion follows.
Next is untimed qualification of a paired, clock-amortized measurement frontend
under `leopard-79h.38.5.4.19.1.1`, with unchanged codecs and safety gates. Any
timed successor requires a separate reviewed, committed and pushed plan.

Only after that gate passes may the default change. The actual default-on
artifact must then be inspected and rerun through focused safety/route tests;
the existing production test intentionally still requires R=199/32 KiB to be
disabled at this checkpoint. Field-reduced configurations and broad project
qualification are not claimed by this checkpoint.

After closing the default-route gap, the already-qualified, untimed AVX2
forward/accumulating prototype remains available for the explicit-AVX2 residual.
No earlier below-threshold result is being added, multiplied or retried.
The parent issue and full Leopard1-versus-Leopard2 goal remain open.
