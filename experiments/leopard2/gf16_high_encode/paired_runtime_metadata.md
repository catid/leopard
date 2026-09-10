# Untimed paired runtime-metadata qualification

Tracker `leopard-79h.38.5.4.19.1.4.1`, 2026-09-10. This supplies missing
observations for the [GF8 control-shift investigation](gf8_shift_audit.md).
It does not establish a timing cause, speedup, stable future timing method,
or promotion. R199/32KiB AUTO remains default-OFF.

## Implementation and observation boundary

`paired_metadata_overlay.py` accepts only the SHA-pinned original paired
timer source. It generates a separate frontend without editing that source
or any codec. All nine cells, seed, allocations, GF8's actual AUTO request,
one-item batch API, four schedules and GF8 groups 1/256 remain intact.
The two public-call lambdas and grouped execution region are source-identical.
This is not a claim that adding instrumentation preserves machine-code layout.

Only abort and synthetic clock variants are linked. `--measure` is rejected
before workload allocation or metadata collection, including when otherwise
well-formed. Independent ELF inspection rejects a real steady-clock import.
Intentional clock-abort cases and synthetic public-call witnesses check the
existing grouping boundary; no real benchmark clock is read.

Fixed-capacity storage records all 4/104 selections with phase/pass/slot
identity: requested backend, resolved context baseline, candidate state and
queried GFNI operation selection. The four observed preflight GFNI call counts
are separate from those queries. Native L1 uses explicit unavailable sentinels
and the label `original_native_compiler_policy`, not a fictional AVX2 request.

Two snapshots, after allocation and after the checked/exercised workload,
record source/reference/scratch/output allocation geometry, parity and pointer
array spans, every input/output row pointer, and real public-function, main
and clearly named metadata-anchor addresses. No metadata capture is inserted
inside the encode lambda. The original test-only per-call witness is retained.

`dladdr` image base and `dl_iterate_phdr` load bias are recorded separately.
The collector-free replayer reads each exact executable's ELF program and
symbol tables and verifies complete function spans against executable load
segments. It checks typed schema completeness, every pointer/length/alignment,
overflow, disjoint allocation ranges and intentional native aliases. Snapshot
equality compares fields, never C++ padding.

## Results

| Check | Verified result |
| --- | --- |
| Profiles | Original native L1, unchanged dual-field L2 Release, full L2 ASan/UBSan with leak detection |
| Positive processes | 270: all nine cells, every permitted schedule/group, check/exercise/synthetic |
| Positive public calls | 478,080, with exact independent state/API/order witnesses |
| Metadata | 19,080 selection records; 540 endpoint snapshots |
| Parity | 261 distinct-output comparisons against native references; 1,582,128,320 bytes |
| Synthetic spans | 7,560; exact clock/public-call ordering, no actual time samples |
| Expected failures | 18 clock aborts, 12 malformed-clock refusals, 42 CLI refusals |
| Supplemental C++ units | 22 capacity/overflow/equality cases per profile; zero public encodes |
| Python tests | 24 per normal/optimized mode, including ten real-record replay mutations |

All positive endpoints match. In the new GF8 frontend, AUTO request 0 resolves
to context AVX2 3, the operation query is non-GFNI, and all four actual GFNI
probe counts are zero. Native parity aliases scratch as expected; the K4096
case has 1024 output pointers, not just its R=512 parity pointers.

The new sanitizer executable demonstrates why image-base terminology matters:
its non-PIE load bias is zero while `dladdr` reports image base `0x400000`.
Release/native PIE placements resolve correctly too. These are observations
of these newly linked frontends. They do not recover historical process maps.

Endpoint equality does not prove uninterrupted allocation stability between
snapshots. The query records and preflight counts do not add per-call route
tracing. There are no frequency/counter observations or causal attribution to
ASLR, alignment, caches or allocation placement. Old timing attempts remain
consumed and untouched.

## Resources and review corrections

All work ran locally; builds/checks were serialized under
`/tmp/leopard-gf8-authoritative.lock`. Builds used 512MiB/no swap and compiler
GC10/4096. Native processes used 256MiB/no swap, CPU limits and wall watchdogs.

| Scope | Peak bytes | Final memory events / swap |
| --- | ---: | --- |
| Frontend build | 153,325,568 | all zero |
| 342-process matrix, maximum | 171,573,248 | all zero |
| Supplemental unit build | 160,346,112 | all zero |
| Supplemental units, maximum | 7,831,552 | all zero |
| Final normal / optimized replay | 55,205,888 / 56,311,808 | all zero |

The initial result replay reached the 256MiB cap with 5,075 `memory.max`
events, zero OOM/kill and no swap. Its source snapshot and log remain retained;
it is not an all-zero resource result. Streaming comparisons now request
eviction only for this lane's own read files, avoiding accumulated parity-file
cache without increasing the cap. Only the checker was rerun, not codecs.

The user-authorized local read-only reviewer found four replay gaps: unchecked
recorded commands/peaks, ignored build-input mappings, incomplete digest-map
coverage and nine reference self-comparisons counted as parity comparisons.
These are fixed and mutation-tested. The native runner already used the
recorded guards. A final review also required exact supplemental unit compile/
link recipes, including source, object, include and sanitizer flags; these are
now checked with three additional real-record mutations. The C++ snapshot-padding
issue was fixed before compilation.
This is Codex review and deterministic evidence under the user's Claude opt-out,
not independent-model `CONVERGED`.

## Evidence and next boundary

Raw root: `/tmp/leopard-paired-metadata.MnmrEa`. Original collection tools,
the initial replay version/failure, final tools, source/header copies,
executables/archives, all outputs/parity, resource logs and final replays are
retained separately. Build/check manifests:

- Build: `f8f2529034e462ccbcfcc94ec3951182a3497a34e1258e01b32b49e3354ae543`.
- Checks: `6d2c4d39169c7159ace7ac36c4962099050e7f42edc124bd98abaaffb961a355`.
- Final tools: `ec9368409da9785c2ef400d92fc66eb4ed9e53ca6c87793d156053297698d94d`.

The next task is a separately reviewed control-shift measurement method using
this qualified metadata, not a retry of either consumed R199 screen. Any real
timing still needs committed AND pushed preregistration, all nine workloads
including AUTO GF8, 5% target/native gates and 2% controls/neighbors. Metadata
qualification alone does not authorize timing or default promotion. The broader
native-Leopard1 performance gap remains open.

## Verified delivery

Read-only bundle: `.research/leopard-79h/paired-runtime-metadata.oQioRv`.
It contains 1,114 files / 1,743,975,663 bytes, including `SHA256SUMS` with hash
`368ff74830ea3ed52abbbff9bece0a14c405f5f8c393c1ee420261d8168f71d6`.
Files are private byte copies, mode0444; all directories are mode0555.
The independent delivery audit verifies the exact namespace, every manifest
entry against both raw and retained bytes, sizes, permissions and non-aliasing
inodes. No original evidence was deleted or rewritten.

Both sealed-copy replays match both raw replays byte-for-byte, canonical SHA
`2450608f2051aef3aad2b87a59f88014bae2cc384c3a3d176cd30cf74d6d3567`.
Their peaks are 55,549,952 / 56,635,392 bytes; the independent full delivery
audit peaks at37,687,296 bytes. All exit0 with all six memory events0 and swap0,
under the same256MiB cap and canonical lock. Retention previously peaked at
67,096,576 bytes with the same zero-event/no-swap result.
All24 tests also pass from the sealed tools in normal and optimized Python,
including the ten real-record mutations, with no skips. Peaks54,005,760 /
55,504,896 bytes; exit0, all memory events0 and swap0. Tests write only their
own disposable fixtures, not the sealed evidence.

The final local read-only reviewer confirmed the complete supplemental unit
compile/link recipes and the three new mutation cases cover the last reported
provenance gap, with no remaining material defect in that bounded static scope.
This is not a full independent-model convergence claim.

Delivery logs and the independent audit source are retained at
`/tmp/leopard-paired-metadata-delivery.TJ0r4O`. Audit source SHA
`9cc5dc2843e4df75d985670081ace042c54a04eeb004455e460a7ed3eed10599`;
audit log SHA
`b653da6260000f36c5071c178ddb09230c51f8c8c8503a9e5fe14d197b84563b`.
The sealed report/result are the pre-delivery snapshots; this appendix and
the repository result append delivery evidence without modifying the bundle.
Next method-review child: `leopard-79h.38.5.4.19.1.4.2`.
