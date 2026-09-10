# Internal-tower GF16 encoder: first untimed public qualification

Tracker `leopard-79h.38.5.4.18.4`, 2026-09-10. **Default-off experiment;
no production worktree source changes, timing, or speedup claim.**
This advances the [qualified butterflies](tower_butterflies.md) into an
isolated full high-rate encoder. The broader Leopard1 performance goal stays
open. This is scoped correctness evidence, not a complete release/API matrix.

## Implementation and costs

The strict overlay copies the pinned original `LeopardFF16.cpp` into a private
build directory. It selects tower operations only for AVX2, transform side at
least256, source-policy bytes above16384, pass bytes at least256 and divisible
by64, and no active sparse plan. This is the existing large copy-first,
non-fused high-profile schedule. Small/tail passes, other operation backends,
low-profile execution and decode stay canonical. Padded-odd framing can enter
the selected aligned prefix and is explicitly tested.

Both actual source-copy sites convert canonical inputs into tower coordinates.
The representation remains internal across inverse transforms, accumulation
and forward FFT. The entire recovery **prefix**, including scratch-backed
holes, converts back before returning. Existing public wire bytes and scratch
geometry are unchanged; there is no K-shard preconversion allocation.

One immutable coefficient cache is built through `std::call_once` before Ops
publication. Actual Release BSS reserves `0x600000` bytes (6MiB) for its tables,
plus a32-byte conversion map and416-byte Ops object. Existing canonical
backend tables remain needed: **this adds6MiB; it does not save2MiB of total
resident memory**. Setup uses a64KiB local subfield array. Cold initialization
latency and cache effects have not been timed.

Source conversion replaces an existing copy, but final output conversion is
new work. The actual AVX2 butterfly object is byte-identical to the preceding
qualified object `fc363f48…`; six-shuffle multiplication and two-shuffle
conversion therefore retain that code-generation evidence. New field/control
objects were independently checked against actual disassembly: no AVX in
baseline-ISA selection/control, and no excluded wider ISA in the AVX2 kernel.

Runtime OFF and excluded paths still use the new `CopySource` wrapper. They
are semantic controls, **not unchanged execution controls**. Separate original
archives and executables are retained for comparisons and any future timing.

## Checks actually run

| Check | Evidence |
| --- | --- |
| Public candidate |21shapes +6special cases in Release, traced Release and full ASan/UBSan/LSan:81 successful processes |
| Original Leopard2 |21shapes in original Release and original full-sanitizer archives:42 successful processes |
| Cross-profile parity |411,831,404 bytes compared across the four non-reference profiles against original Release; no reference-self comparison counted |
| Native Leopard1 |Five matching retained native workloads, compared against all five profiles:206,602,240 bytes equal |
| Actual cached coefficients |Per profile:65,536 log/sentinel cases ×32 input-pair basis lanes =2,097,152 pairs, checked against independent polynomial and original scalar arithmetic |
| Selection |9,000 directed predicate cases per profile; selection alone never initializes the cache |
| Unsupported callbacks |Nine deliberately aborting callbacks/fused-hint violations per profile; core dumps disabled |
| CLI refusal |Three malformed/timing requests per kernel profile reject before initialization |
| Replay/overlay tests |Seven tests pass normally and under Python `-O`, from source and retained copies |

Public shapes cover all four K1000/R199-or200/32-or64KiB AVX2 combinations,
misaligned compact tails,16384/16448 policy boundary,128/130-byte exclusions,
R128/R129 side boundary, single/partial first message blocks, a selected
side512 transform, the previously enabled AUTO/GFNI cell, explicit GFNI,
scalar, and GF8. All21 compare full parity and scratch size with original L2.

Each candidate shape tests six output masks (full, prefix1, prefixR-1,
middle/last, alternating, none), plus full/subset one-item batches. Traced
profiles check per-mask conversion bytes and output **prefix**, not selected
output count. AUTO full output uses GFNI while partial output falls back to
AVX2 and can initialize tower; lifetime cache assertions account for that.
Source hashes, output/source/scratch guards, unrequested outputs, short-scratch
and odd-GF16 rejection, and three payload-overlap rejections are checked.

Special cases cover padded-odd encode/decode and bad systematic-pad rejection,
low-profile and GF8 canonical exclusion, ordinary GF16 decode, cold concurrent
initialization from two shared-codec callers with disjoint buffers, warm
encode/decode overlap, and a two-item pool batch with different payloads,
ragged lengths and a sparse output mask. Ordinary, preflight-scratch and bound
batch execution each start with poisoned outputs; binding creation must not
encode. An invalid second item must preserve both items' fresh output/scratch
poison. Worker-local diagnostic counts are checked for the explicit caller
threads, not represented as totals across internal pool/OpenMP workers.

## Limits and retained failures

The exact native-L1 comparisons cover K1000/R200/32768, K1000/R199/65536,
K1000/R200/65536, K1000/R199/32768 and K4096/R512/4096, using the prior
`avx2-adjacent-public.NFkTeG` native archive/seed evidence. No new native L1
execution or timing occurred. Other shapes use original-L2 controls.
Decode checks restore one or three originals with every parity shard
available; they are not an exhaustive erasure-pattern suite. No TSAN,
cross-host, large-thread-count or broad monolithic API suite was run.
Exhaustive cached-log testing exercises the inverse pair; other forms rely
on preceding isolated-kernel qualification plus public-transform coverage.

The user-requested local read-only Codex reviewer found test defects, not a
new arithmetic defect: AUTO lifetime-cache accounting and batch checks that
could miss a no-op or idempotent premature execution. Those were corrected
before final qualification. An initializer-list pointer-type issue was also
fixed before compilation. This is static review plus deterministic tests,
not a Claude fixed-point `CONVERGED` claim; the user's no-Claude instruction
remains in force.

Two failed native harness generations remain retained: the initial AUTO cache
assertion, and an overstrong scratch-preservation assertion after payload
overlap. General validation can use scratch range tables; the public API only
promises untouched scratch for **metadata** overlap rejection. Corrected
tests retain size-rejection scratch preservation, public data/guard checks
and zero tower entry. A replay syntax error was fixed with its original
source/log retained. No codec arithmetic was changed to accommodate these.

## Evidence and next gate

Final raw root: `/tmp/leopard-tower-encoder-correctness.cvYHeW`.
Read-only bundle: `.research/leopard-79h/tower-encoder-qualified.kvqrsu`.
653files /677,242,732bytes, manifest SHA-256
`56e78fef5472278effb208ba6e68609838b8e2b56e254dd95b0924c2d0d09abb`.
It includes original archives, actual objects/disassembly, sources, raw check
records, parity files, selected native-L1 oracle files, replay tools and failed
harness source/log snapshots. Earlier large failed-generation binaries remain
at their original temporary roots; none were deleted.

Final build peak175,407,104/536,870,912bytes. Maximum native scope
257,683,456/268,435,456bytes; those build/native scopes had zero memory events and swap.
This is close to the cap, so checks must remain separate processes. Retention
peak37,126,144bytes. The initial retained-only replays completed with matching
results but reached256MiB through file cache (2,685/939 `memory.max` events,
zero OOM/kill/swap); these are not all-zero resource successes. The current
replayer drops verified file-cache pages with `posix_fadvise` and passes fresh
normal/optimized retained-only replays at118,079,488/117,260,288bytes with all
events/swap zero. The corrected script and old/new replay/resource logs are
retained separately in `tower-encoder-streaming-replay.FMS613` under the same
research directory; script SHA-256 is `0b5cefac455ce1d24d0630819f98ee5bd58ef80ab1e8f7b31a4a224007b9b8b1`.
The original
bundle's older replayer and capped logs remain retained; use the current
repository replayer for resource-bounded reproduction. Both raw and all
retained-only outputs agree; the full retained manifest verifies.
Replay rechecks actual archive members and
disassembly, raw records/resources, exit evidence and full parity bytes,
without executing a codec. Native39 kernel cases include deliberate nonzero
exits; only successful/expected outcomes are claimed.

No timing is authorized by this correctness result alone. The next tracked
work is public measurement-front-end qualification and method review, then a
distinct committed **and pushed** preregistration using original L2/native L1
comparators. Actual conversion/setup/cache costs must remain in scope.
ExperimentT's eventual10% end-to-end gain and2% control/neighbor gates remain
unchanged. No previous scheduling screen will be retried, pooled or added to
this result. The experiment stays default-off.
