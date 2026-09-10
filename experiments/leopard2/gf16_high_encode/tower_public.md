# Tower encoder public frontend qualification

Tracker `leopard-79h.38.5.4.18.4.2`, 2026-09-10. This is **untimed
qualification**, not a performance result or timing preregistration. The tower
encoder remains default OFF; production codec files and all six reused archives
are unchanged. The broader Leopard1 performance objective remains open.

## Actual comparison boundary

`tower_public_overlay.py` adapts the exact previously qualified public frontend
by SHA-256. Its generated source changes the link interface and diagnostic
records, not workloads, public calls, buffers, schedules or clock boundaries.
Original L2, tower Release and traced tower links consume **one shared Release
driver object**. Original-sanitizer and tower-sanitizer links share a separately
fully instrumented driver. Native Leopard1 has its own API-specific driver and
unchanged native archive, not an AVX2-restricted substitute.

All nine previous cells are retained:

| Cell | K / R / bytes | Requested path | Tower when ON |
| --- | --- | --- | --- |
| 0 | 1000 / 200 / 32768 | AVX2 ordinary | Selected |
| 1 | 1000 / 199 / 65536 | AVX2 ordinary | Selected, two byte tiles |
| 2 | 1000 / 200 / 65536 | AVX2 ordinary | Selected, two byte tiles |
| 3 | 4096 / 512 / 4096 | AVX2 ordinary | Excluded by source-policy size |
| 4 | 1000 / 199 / 32768 | AUTO, actual AVX2 | Selected; GFNI extension remains OFF |
| 5 | 1000 / 200 / 32768 | Explicit GFNI | Excluded |
| 6 | 1000 / 200 / 32768 | AUTO, actual GFNI | Excluded |
| 7 | 17 / 7 / 64 | GF8 AVX2 | Excluded |
| 8 | Cell0 | Actual one-item L2 batch | Selected |

Native cell8 uses ordinary `leo_encode`, not a nonexistent native batch API.
Native parity is read directly from its first R work buffers, with no extra
output-copy penalty. Runtime OFF still calls the tower overlay's CopySource
wrapper: it is not an untouched original execution control.

## What the frontend measures and checks

Four preflight calls record the actual AUTO route and tower work separately.
All four tower snapshots include process-lifetime cache initialization. For
example, an eligible `0110` schedule reports initialization counts `0,1,1,1`;
an all-OFF or excluded workload remains zero. Reset clears work counts only.
Untraced Release still exposes cache lifetime, while its work counts are zero.

The independent structural oracle derives source/output conversion bytes and
all inverse, forward and accumulating pair calls. Tower range wrappers execute
all four edges even at zero skew; the oracle does not reuse the canonical
range's skip counts. At cell0 one selected public encode has one pass, 1000
source conversions /32,768,000 bytes, 200 output conversions /6,553,600 bytes,
and 3672/917/384 inverse/forward/accumulating pair callbacks.

The subsequent execution loop has four warmup and21 sampled-position passes,
each containing four slots. Selection, route inspection, result normalization,
storage, full parity comparison and guard checks are outside each span. Public
result checks, per-call bookkeeping and repetition are inside. A span encloses
one complete encode, except GF8 also qualifies groups of256 complete calls.
That means104 or25,604 calls per exercised process and84 potential spans.
The synthetic link verifies168 endpoints against independently witnessed actual
public-call counts and stable buffers. No real benchmark-clock request is run.

Full parity is checked after **every group**. This touches output memory and
influences subsequent cache state. Any later timing using this frontend is warm
throughput under that policy, not single-call latency, cold initialization or
context creation. Source conversion and final output conversion are inside
public encode. The added6MiB table cache and cold setup remain real costs;
initialization is outside these warm spans and must not be claimed free.

## Validation and scope

The fixed inventory covers native, original Release, original sanitizer,
tower Release, trace and full ASan/UBSan/LSan profiles. It includes all runtime
schedules `0110`, `1001`, `0000`, `1111`; original `PPPP`; native `NNNN`;
plain clock-free execution, synthetic clocks, intentional abort clocks,
malformed requests, grouped-clock arithmetic failures and actual guard tests.
Every process has its own256MiB/no-swap scope,60-second CPU and120-second wall
limit, core dumps disabled, and the canonical build/test lock. The controller
does not alter unrelated processes or access remote worker hosts.

`verify_tower_public.py` executes no codec or collector. It checks actual child
exit markers, stderr, all resource counters, source/archive/executable hashes,
shared-driver link recipes, records and full parity bytes. New native outputs
are also compared with the pinned earlier native evidence. Reference
self-comparisons are excluded from reported byte totals. Retained replay uses
source snapshots, not mutable live source files. Retention shares identical
parity files only within its private bundle; byte totals count logically
distinct collected outputs, not physically distinct retained inodes.

The user-requested local read-only Codex reviewer found no actionable static
bug in the frontend/link/build, oracle/runner, or final retained-input/resource
checks. This is static
review plus deterministic/native evidence, not a Claude fixed-point
`CONVERGED` claim. The user's Claude opt-out remains in force.

The preceding whole-encoder boundary, masks, concurrency, decode, ISA and
kernel qualification in `tower_encoder.md` is carried forward through exact
archive identity, not rerun or broadened by this workload-only frontend.
No TSAN, cross-host or full release/API matrix is claimed.

## Completed native result

All501 processes completed with the expected outcomes:334 successes,
45 intentional clock aborts,24 malformed-clock refusals,96 malformed-CLI
refusals and2 expected sanitizer poisoned-read failures. The334 include330
public workloads,2 grouped-timer unit processes and2 canary checks.
All321 non-reference-self parity comparisons pass, totaling1,886,452,800 bytes.
Fresh native outputs also match60,981,696 bytes of independently retained
native-L1 reference output across all nine cells.

Ten pure adversarial tests pass normally and under Python `-O`. Raw normal
and optimized collector-free replays agree on all records, counts, full-byte
comparisons and resource evidence. Native checks peak at199,360,512 bytes /
256MiB; their controller at40,976,384 /256MiB; the serial frontend build at
326,787,072 /512MiB. All six memory-event counters and swap are zero.

Raw evidence: `/tmp/leopard-tower-public.30Gv2i`. Build metadata SHA-256:
`d304b72038eb460efd69ccd2526bd26c187fec083be7ba2f2724c7085a479818`.
Native record journal SHA-256:
`71c9d4d37e96322f5381558a452c3fd7e51da24a7b4f9ac9c3066dc7f0e7320b`.

| Exact future plain frontend | SHA-256 |
| --- | --- |
| Native Leopard1 | `9f777340c054041870f15c916374ae5728ea3c4cb6a5291c12ffe1cd0e942115` |
| Original Leopard2 | `dd0c82b408e00744716f809c1c87a7f972eb5d7be010079b5154c44cbd290f68` |
| Tower Release, runtime OFF/ON | `45ba62393647dcc64873eb592402e551b33de37719245540310ca0bf988920b9` |

The shared L2 Release driver is `99223acc12f2f616c0e4e5c0cb39d7a1a5207aa27570eedf3f08a1ec9290aa75`.
These executables have been exercised without reading real benchmark clocks.
Their future timing comparison is not established by native validation.

The final read-only bundle is
`.research/leopard-79h/tower-public-qualified.KPkPef`:1775 files,
2,198,711,456 logical bytes and296,603,648 allocated bytes. Private parity-only
deduplication saves1,906,113,600 bytes; executables are separate copies.
Manifest SHA-256:
`c52a32f869cdd5c8aa4886fcc9f4bc2ee77fa7f1f9a2bc27cd295600b29a1c20`.
All1774 manifest members and the complete read-only namespace verify.

Both read-only replays pass and produce exactly the same result as both raw
replays. Ten tests also pass in both Python modes using only copied retained
Python modules and the generated-build directory's pinned original driver
fixture. The asset is staged beside the test module; no source is rewritten.
Largest replay peak46,923,776 bytes, retention67,280,896 bytes, retained tests
17,432,576 bytes, each under256MiB with all events/swap zero. Final audit and
delivery logs are in `/tmp/leopard-tower-public-delivery.bwlqN2`.
This paragraph and the structured result's retention metadata supplement the
sealed report snapshot, which records the earlier retention-pending checkpoint.

## Next gate

Follow-up `leopard-79h.38.5.4.18.4.3` requires separate method review and committed **and pushed**
preregistration must fix target roles and original/native comparisons before
any performance run. Cell3 is now excluded; cell4 is selected. Do not inherit
the old scheduling experiment's role assignment or5% threshold. Tower retains
its10% end-to-end gain and2% control/neighbor gates. No old timing attempts are
retried, pooled, trimmed or added to this candidate.
