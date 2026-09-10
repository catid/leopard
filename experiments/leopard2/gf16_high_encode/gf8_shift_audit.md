# Retained GF8 shift investigation — 2026-09-10

Tracker: `leopard-79h.38.5.4.19.1.4`. **Untimed evidence, not a resolved
performance cause or qualification for promotion.** No codec was executed,
rebuilt or modified; no consumed timing attempt was repeated.

## Findings that change the next action

1. **The paired and tower GF8 cells have different requested routes.** Both
   perform K17/R7/64-byte ordinary encoding with AVX2 arithmetic, but the paired
   driver requests AUTO for cell8; tower requests explicit AVX2 for cell7.
   `auto_requested` therefore differs, and AUTO traverses more eligibility
   predicates before rejecting the GF16/AVX-512 routes. This is a cross-driver
   distinction, not an OFF/ON difference within either screen. Do not call the
   two frontends identical or use their absolute durations interchangeably.
2. **Tower preserves the GF8 implementation objects.** All23 original members
   other than `LeopardFF16.cpp.o` are byte-identical. Tower adds two private
   experiment objects. The paired screen's complete24-member codec archive is
   byte-identical to the tower screen's original-Leopard2 archive.
3. **Linked placement changes despite unchanged objects.** All35 inspected
   non-main function instances have identical instruction offsets and
   symbol-normalized disassembly in original, tower and paired executables.
   Their loaded-relative addresses differ. Original/tower share the same
   34,088-byte driver object; their grouped public-call loop has the same13
   normalized instructions but different64-byte alignment.
4. **The paired failure is a sustained between-process shift.** Its failed
   round's first process stays near167ns per grouped call; the other three stay
   near150–151ns. It is not explained by deleting one maximum or one initial
   sample. All12 same-OFF processes and all1,008 spans are retained in the
   diagnostic projection; the original control decision remains unchanged.
5. **A causal explanation is still unavailable.** Those old processes did
   not retain their mappings, allocation addresses or frequency/counter state.
   Static inter-binary placement alone cannot explain variation between runs
   of one fixed paired executable. ASLR, caches, allocation placement,
   predictor state and frequency are hypotheses, not established causes.

## Exact call path and selection boundary

The local read-only reviewer traced the source path for this shape:

`leo2_encode → EncodeInternal → validation → SelectTransformEncodeOps →
ExecuteTransformEncodePass → ff8::ReedSolomonEncode`

Contiguous source/output slabs satisfy packed-buffer validation. There is one
64-byte tile and16 work slots; the inverse processes blocks8+8+1, followed by
the seven-output forward transform. The special two-block T8 terminal requires
K9–16, so it does **not** handle K17. No GF16 tower conversion or butterfly is
on this path. The audit includes the relevant validation, field driver,
weighted/output inverse, forward-range, two-way and XOR callback symbols,
plus nearby GF8 functions. It retains every duplicate local-symbol occurrence
from different backend variants rather than selecting an ambiguous name.

There is nevertheless a GF16-related global read: `UseAutoGF16GFNIEncode`
loads `g_auto_gf16_gfni_encode_mode` before testing codec eligibility.
Both frontends normalize this mode to1 before measured groups. The GF8 shape
rejects eligibility before reaching the boundary/R19932 mode reads. Thus
neither “the GF8 path reads no GF16-related state” nor “toggling R19932 changes
GF8 arithmetic” is supported.

Sources: frozen `leopard2.cpp`, paired `paired_timer_r19932.cpp`, tower's
`avx2_adjacent_public.cpp` base and pinned build record. The tower build record
binds the same driver object into both original and experimental links. The
qualified tower overlay changes diagnostics/schema, not backend requests or
the public-call loop. The shared `PairedGroupTiming.h` encloses complete public
calls, result checks, counters and repetition; selection, route inspection,
parity checking and output are outside those spans. This source review is
not a dynamic execution trace.

## Binary layout evidence

Addresses below are ELF-relative values, **not captured runtime addresses**.

| Location | Original L2 | Tower executable | Observation |
| --- | ---: | ---: | --- |
| `main` | `0x48d0` | `0x4970` | Both7,128 bytes/1,500 instructions; address changes160 bytes |
| Grouped public-call loop | `0x57b0` | `0x5850` | Same13 normalized instructions; offset within64-byte line48→16 |
| Grouped `leo2_encode` call | `0x57da` | `0x587a` | Offset within64-byte line26→58 |
| `leo2_encode` | `0x2c000` | `0x2c200` | Both7,102 bytes/1,707 instructions; still64-byte aligned, page offset0→512 |
| `SelectTransformEncodeOps` | `0x13390` | `0x13590` | Both1,259 bytes/272 instructions; normalized disassembly equal |

The paired executable's `leo2_encode` is at`0x2bbc0`, size7,102, also
normalized-equal. Its frontend is genuinely different: `main` is6,442 bytes.
The tower/original main comparison retains eight differing symbolic string
offset annotations; it is not mislabeled fully normalized-equal. The shared
driver-object identity and the separately matched group loop are the narrower
evidence used here.

Normalization removes only symbol-annotated RIP-relative displacement spelling
and direct branch/call address spelling. It preserves constants, registers,
instruction offsets, symbol names and target offsets. All raw disassembly and
encoding digests are retained. This is **not** a formal semantic-equivalence
proof, runtime indirect-target validation, or attribution of cycles to layout.

## Paired failed round, without trimming

Each entry summarizes84 spans, each averaging256 complete public calls.
These are not individual-call latency samples. Round numbering is zero-based.

| Round1 same-OFF process slot | All-span median ns/call | First four-span pass median | Last four-span pass median |
| --- | ---: | ---: | ---: |
| 0 |167.046875|167.494141|168.062500|
| 1 |150.457031|150.341797|150.458984|
| 2 |151.085938|151.259766|151.044922|
| 3 |151.085938|151.261719|151.476562|

Within each process, all four schedule-slot medians are close. The full result
also includes rounds0/2, minima/maxima, all process and schedule-slot medians,
and hashes binding each raw output to the pinned retained manifest and journal.
No new performance ratios, exclusion rule, pooled estimate or promotion
decision are computed.

## Validation and next experiment boundary

Ten pure/adversarial tests pass normally and under Python`-O`. Complete
retained-only audits agree byte-for-byte in both modes, inspecting108 function
instances across three binaries. Final audit peaks46,125,056/50,708,480 bytes;
test peaks16,863,232/21,389,312, all under256MiB with all six memory events and
swap zero. Jobs are serial under the canonical lock. There are **zero codec
executions and zero new timings**.

Reviewer findings fixed: duplicate/type-coercing JSON comparisons, and missing
trailing instruction differences. Tests additionally caught overbroad address
normalization; the failed test/source snapshot is retained, not relabeled
passing. Final bounded local read-only review found no further actionable
defects. No Claude or independent-model`CONVERGED` claim.

Next qualify **untimed** diagnostic metadata in a separate frontend version:
explicit requested/context/operation routes, actual executable load placement,
relevant function addresses and source/output/scratch allocation relationships.
Preserve the paired campaign's actual AUTO GF8 request and all nine workloads;
do not silently substitute explicit AVX2 to obtain a passing control. Verify
metadata outside grouped calls with synthetic/abort clocks and original-codec
parity. Such metadata may support a later discriminating method, but cannot
recover missing old process state or itself establish stability.

Any subsequent real timing requires separate method review and a fresh
committed-and-pushed preregistration. Preserve target/native5% and controls/
neighbors2% requirements; tower's rejected experiment retains its separate10%
gate. R19932 remains default-OFF, rejected tower stays experiment-only, and
the parent investigation and full performance objective remain open.

Raw audit workspace: `/tmp/leopard-gf8-shift-audit.SuEN62`. Final full result
SHA-256:`4812ff5c64cc115303936c7376c85da6c5336952413d502ef4f89da9c282eb8b`.
The compact repository result records local read-only retention metadata and
references the unchanged paired/tower evidence bundles; it does not duplicate
or mutate their historical snapshots.

Final retention: `.research/leopard-79h/gf8-shift-audit.w9805W`,123 files /
4,751,088 bytes; manifest
`33a353fda2dc44995b6c06482ccd1bccfffed9f9cc86845db64ea34e1b84bd0d`.
The first manifest mistakenly included itself and failed verification; that
manifest and log are preserved. The corrected manifest excludes its own file,
passes verification and covers the read-only tree. The final replay from that
read-only copy reproduces the exact result above at46,047,232-byte peak, all
events/swap zero. An earlier replay ran before sealing and is labeled preseal,
not counted as read-only evidence. Retention peaked4,882,432 bytes with all
events/swap zero. The sealed report is the pre-retention checkpoint; this
appendix and the [compact result](results/gf8_shift_audit_20260910.json) add final
delivery metadata without rewriting it. Untimed follow-up:
`leopard-79h.38.5.4.19.1.4.1`.
