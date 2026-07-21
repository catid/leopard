# Leopard2 plan-to-production gap audit

This audit treats the tracked `plan.md` as the acceptance contract rather than
using historical issue closure as proof. The plan checkpoint has SHA-256
`670ce14053ec528bea1cf72ba70b8fcd90eec7349bec03c80332adbb318ef753`.
The first complete source pass was made at `cdf2211` on 2026-07-18, and this
revision was independently rechecked against integrated source through
`1962218`, then reconciled with the Algorithm 4 pruning integration at
`122d7de`, the Algorithm 5 no-copy integration at `4dcd887`, and the batch
alias-safety checkpoint at `866dd3a`. The repeated bug campaign continued
changing source after the historical test snapshots cited below, so a final
current-HEAD test repeat is still required. The durable work items cited below
are Beads issues; this document is an evidence map, not a second task list.

## Production algorithm check

| Plan requirement | Current production evidence | Disposition |
|---|---|---|
| Legacy high wire profile | `leo2_codec_create` derives `T=ceil_pow2(R)`, `N=ceil_pow2(K+T)`, and the parity/message/shortened coordinate map. `ReedSolomonEncode` retains the block-IFFT accumulator and truncated final FFT. Legacy golden, API, arbitrary-count, GF16 matrix, and exact-main diagnostic tests compare parity bytes. | The profile and encoder are implemented. Final current-HEAD compatibility repetition remains part of `leopard-79h.15`/`.40`. |
| Low wire profile | Codec setup derives `P=ceil_pow2(K)`, `N=ceil_pow2(P+R)`, with the systematic prefix and punctured parity suffix. `ReedSolomonEncodeLow` performs one padded-P inverse transform and shifted evaluations, including `R>K`. | The profile and arithmetic are implemented. Commit `bd681a7` removes the former whole-P coefficient copy, but its default-production performance promotion remains provisional pending `leopard-79h.26.1.1`. |
| R10 Algorithm 4 | `ReedSolomonDecodeLowPrunedPlanned` and its tiled form perform exact-mask P-point block inverse transforms, the active-parent derivative/weighted reduction, an exact requested-output P-point transform, and original-only recovery from plan-owned locator factors. The immutable GF8/GF16 schedules retain the mature prefix/dependency transforms as fallbacks when pruning does not save byte-heavy work. | Arithmetic and active-parent behavior are implemented and independently tested against direct interpolation, materialized/tiled fallbacks, byte tails, repeated/concurrent plans, and qualified context backends. Production performance promotion remains provisional pending isolated crossover evidence under `leopard-79h.18.1.12.1`. |
| R10 Algorithm 5 message-only form | `ReedSolomonDecodeHighPrunedPlanned` and its tiled form construct the T-wide accumulator, apply locator factors, invert to the evaluator, and evaluate only output blocks containing missing originals. Commit `e2ce390` adds immutable exact-mask C1 schedules for GF8/GF16 input and output transforms. Commit `4dcd887` removes the former whole-T copy by evaluating the immutable coefficient accumulator through an out-of-place root butterfly in prepared, materialized, tiled, and exact-sparse variants. Complete GF8 receive groups fuse source staging with the first two inverse layers only on qualified SSSE3/AVX2/AVX512; scalar, NEON, and GF16 retain copy-first after their candidates were below threshold or unmeasured. The later locator-weighted inverse boundary fuses four row weights with its first two inverse layers for measured GF8 AVX2 and the explicit AVX512 recompilation at T>=64, at least half-live rows, and 16--256 KiB kernel passes; every neighboring backend, field, loss, and size case retains scale-then-IFFT. Commit `1000e66` makes the aligned materialized reveal multiply write the validated public destination directly, removing its alias temporary and final gather; compact tails retain scratch reveal and gather. | Algorithm 5 arithmetic and all three promoted no-copy boundaries are implemented, isolated, and correctness-tested. The direct-reveal gate found 5.5--9.6% gains in maximum-loss materialized cells, 2.9--3.4% in AUTO-selected production cells, and neutral forced-tiled controls. Algorithm 4 remains independent and unchanged by the weighted Algorithm 5 boundary. |
| Active-parent derivation | `docs/leopard2_math_and_sources.md` defines active N, shifts, normalization, coordinate maps, shortening/puncturing, locator and reveal factors. Direct GF(2^4), GF8/GF16 transform, locator, arbitrary-count, high/low acceptance, and XDRS differential evidence exercise those identities. | This pass found no arithmetic or coordinate-map contradiction. That is not a final proof or source audit: formula traceability and the final literature refresh remain `leopard-79h.27`. |
| Independent fallback/oracle | The retained generic active-parent decoder remains selectable. Direct field, generator/interpolation, repair, transform-differential, GF(2^4), GF8 and GF16 tests do not use the optimized path as their only oracle. | Implemented. |

The specialized decoders take the codec's active parent `N`, not the field
order. Codec setup precomputes permanent puncture/shortening erasure vectors
and caches their locator contribution in the sparse-direct region. Dense
locator setup calls the active-N Walsh construction, so the former unconditional
256/65,536-entry behavior is gone.  The `leopard-79h.11.1` checkpoint proves
that current plan selection always contributes exactly `R` non-permanent real
plus virtual erasures, so the codec's sparse-cache decision is exact.  Both
coordinate- and transform-domain dense cache prototypes retain both Walsh
transforms, add an `N`-entry modular pass, and consume `N` persistent field
symbols.  `docs/leopard2_permanent_locator_cache.md` records the exhaustive
decomposition tests and negative dense promotion result; rebuilding the explicit
union is the non-dominated dense production path.

The default-off locator experiment under `leopard-79h.29` reached a negative
algebraic result for its dense product-tree form and a promising epsilon
correctness checkpoint. It still lacks isolated setup/reuse timing and a final
disposition. It is a comparison experiment, not evidence that either candidate
is in production, and should not block the CPU release once active-N Walsh and
sparse setup have passed their production gates.

## API, memory, and execution check

| Plan area | Verified behavior | Remaining boundary |
|---|---|---|
| Stable API and identity | `leopard2.h` exposes an opaque, thread-safe context and immutable codecs/decode plans; profile, field, layout, parent and padded-side introspection; scratch queries; one-shot/reusable decode; batch calls; result strings; and per-context scalar/SSSE3/AVX2/AVX512 selection. Legacy `leo_*` signatures and contracts remain preserved, although `leopard.h` itself has received additive documentation changes. AUTO is deterministic and never selects a CPU-dependent wire profile. | The context is not immutable: its pool starts and grows lazily. Endian-stable serialized identity remains experimental under `leopard-79h.18.23`; without it, profile/field identity lives in the in-memory codec and caller metadata. |
| Allocation and aliases | Scalar encode and reusable plan execution use caller scratch and perform no allocation. Inputs may alias inputs; outputs and the full supplied scratch range are rejected on overlap. No-loss decode is a true zero-scratch no-op. R=1 and K=1 direct paths use range-only scratch after `8f24962`. | A batch call can lazily create/grow context workers, and no caller-executor adapter exists. The strict setup-only allocation guarantee remains `leopard-79h.14.2`. The one-shot decode wrapper intentionally allocates plan setup state. |
| Arbitrary byte lengths | GF8 accepts all positive physical lengths. Native GF16 accepts complete two-byte symbols. The separately versioned padded-odd layout plus pack/unpack helpers maps every positive application payload to an even physical shard without silently changing native wire identity. SIMD and ragged 64-byte tails are tested. | Native GF16 `leo2_encode`/decode rejects an odd physical `shard_bytes` by design. Decision `leopard-79h.39.1` accepts the versioned padded-odd wire layout as the arbitrary-positive-application-payload contract while retaining complete-symbol native GF16 semantics. |
| Scratch target | Specialized decode uses `min(N,2P)` low slots or `min(N,2T+L)` high slots where profitable; aligned caller inputs are referenced directly and only a ragged 64-byte tail stages public coordinates. Direct repair has no shard-data scratch. Encode reuses 2P/2T work and bounds tail staging. Algorithm 5 now evaluates its immutable T-wide coefficient accumulator without copying it for each requested output block. | Further encoder/output pruning is `leopard-79h.26.4`. Ragged staging retains fixed `O(K+R)` one-tile terms. The broader transform fallback still has configurations whose staging exceeds the plan's aspirational tiled `O(min(P,T))` target. |
| Parity subsets and rebuild | Null recovery outputs are never written. Missing parity is not rebuilt during decode; tests explicitly call `leo2_encode` after combining surviving and restored originals and compare rebuilt parity byte-for-byte. | Decision `leopard-79h.39.1` accepts ordinary `leo2_encode` as the explicit parity-rebuild operation instead of duplicating codec logic in a second entry point. Sparse masks still compute unused prefixes inside encoder transform blocks (`leopard-79h.26.4`). |
| Runtime backends | Baseline objects are portable; SSSE3, AVX2, and explicit AVX512 live in qualified members behind immutable context ops tables and startup KATs. Different contexts can select lower qualified backends concurrently. AUTO reports AVX2 as its baseline; the current model-scoped Zen 5 candidate may widen a complete-output legacy-high GF16 transform encode after optional AVX512 qualification, while unknown models and explicit requests retain their exact table. Compile-time GF8/GF16 options and default CUDA absence are tested. | The Zen 5 selector still requires independent review and fresh exact-main evidence under `leopard-79h.38.5.4`. Required native NEON remains production work in `leopard-79h.13.7`; Windows compiler isolation is `leopard-79h.13.5`. PMULL/SVE/SVE2 in `leopard-79h.18.25` are separate optional research and must not be presented as the native-NEON production gate. |
| Batch and multicore | Batch encode/decode use a context pool after lazy start, deterministic error selection, per-item caller scratch and immutable shared codecs/plans. Commit `1e69cd5` makes default sizing honor the allowed CPU set and avoids starting workers for unused/single-item contexts. Commits `cb821e3` through `9dbb2bd` add whole-batch atomic alias preflight: immutable inputs may be shared, but every writable shard and scratch range is rejected on any cross-item read/write/metadata overlap before execution. API version 4 adds a size-queryable caller-owned `O(M log M)` interval preflight for large batches while retaining the correct allocation-free compatibility entry points and their exact small-batch path. The lab runner records affinity/topology, stable seeds, resumable jobs and per-job output. | The pool still uses a shared atomic `fetch_add` task queue rather than the planned static schedule. A caller executor/prewarm path, worker pinning, first-touch/node partitioning, per-NUMA pools or accumulators, and final scaling evidence remain under `leopard-79h.14`, `.14.2`, `.14.4`, `.23`, and `.24`. |

The production dispatcher is deterministic, but it is not yet the plan's full
region model. Current decisions use selected combinations of K, R, N, loss
count, shard bytes, field and backend, with one measured single-versus-multiple
batch exception. Decision `leopard-79h.38.4.2` intentionally does not add a V1
expected-reuse or persistent exact-batch-count hint: the batch call already
knows its exact count, every item must fit the ordinary plan scratch query, and
an immutable plan may be shared across callers with different actual reuse.
Reuse and exact batch count remain benchmark/reporting axes; a future measured
setup policy needs an additive versioned plan-options API. Thread count and
most batch shapes do not yet inform direct/generic/specialized selection.
`leopard-79h.38.4` owns the production existing-path dispatcher independently
of exact-profile C10 research (`leopard-79h.18.1.11`).

## Correctness and evidence check

The production graph contains separate direct-field/interpolation, active-LCH,
locator, generator-row, generic-decode and optimized-codec tests. The test-only
GF(2^4) programs exhaust field identities, active-parent normalization, locator
subsets, MDS/rank properties and bounded direct-repair profiles. GF8/GF16
differential suites cover coordinate boundaries, tails, requested parity
subsets, virtual erasures, repeated plans and concurrent execution.
`docs/leopard2_xdrs_differential.md` records the pinned XDRS algorithms and the
coordinate/code-definition difference instead of treating research defaults as
an oracle.

This does not yet satisfy the plan's request to compare an incremental locator
plan against fresh generic/low/high/direct/pruned paths: incremental add/remove/
swap plans remain an experiment under `leopard-79h.18.5`.

Compact golden vectors and exact-main adapter digests protect legacy parity and
recovery. Compile-time GF8-only/GF16-only matrices verify that omission does not
silently select a different canonical field. Running
`python3 tools/leopard2_operation_counts.py self-test` executes 375 checks that
distinguish butterflies, fixed
multiplications, XORs and byte traffic from ISA changes; the low-copy source
guard prevents the former whole-P copy from returning unnoticed.

Beads notes for the combined low-copy/high-C1 checkpoint report Release 65/65
and focused ASan+UBSan 17/17. Those are historical milestone results, not a
committed result bundle bound to `1962218`. The independent strict/static pass
under `leopard-79h.40.2` separately reports GCC and Clang 65/65 plus focused
sanitizer evidence on its audited candidate. Subsequent initialization, option,
size-validation and other bug fixes changed production source, so the final
current-HEAD Release/sanitizer/fuzz repeat remains pending under
`leopard-79h.15`, `.40.3`, and `.40.4`. MSan, native ARM runtime, external
128-CPU/NUMA, and the final 128-seed campaign must not be inferred from local
x86 evidence.

## Exact Leopard1 performance check

The bounded diagnostic in
`experiments/leopard2/main_compare/results/2fce390-allk-diagnostic/` compares
exact Leopard main `6e5725e` with Leopard2 candidate `2fce390`, not current
integrated source. It completed 2,760 cells: every GF8 K plus 238 representative
GF16 cells, for 11,040 child invocations. All round trips and comparable
original/parity/recovery digests match. Thirty workers intentionally saturated
the allowed CPUs, so these timings are attribution and prioritization evidence,
not authoritative promotion measurements.

Within that exact 2,760-cell population, median exact-main/current-tree-legacy
ratios are 0.891 for encode and 0.926 for decode. The older 0.919 decode value in
`leopard-79h.38.1` describes a partial 1,853-cell population and must not be
mixed with the committed full-population value without labeling the population.
The diagnostic separates exact main, the retained legacy entry points in the
`2fce390` tree, and Leopard2 so portable-backend loss is not mislabeled as an
Algorithm 4/5 cost.

The evidence does not support a universal-speedup claim:

- the retained legacy path in `2fce390` is already slower than exact main at
  the aggregate medians, isolating a portable/out-of-line backend throughput
  gap (`leopard-79h.38.1`);
- Leopard2 improves the retained legacy decode aggregate median, but the forced
  balanced full-loss diagnostic only establishes that specialized loses to
  generic at T=16 and T=32 for K=15,16,30,31,32. A rule for every T at or above
  eight is a hypothesis requiring the pinned matrix in `leopard-79h.38.2`, not
  a measured result;
- the historical R=1 change at `8f24962` reduced K=129, 4 KiB time by about
  13% encode and 28% decode versus `2fce390`, while exact main remained about
  7–8% faster; and
- exact main rejects R>K, so target low-rate Algorithm 4/low-encoder cells have
  no valid Leopard1 timing counterpart. Forced low at the balanced boundary is
  not evidence for AUTO.

The high-C1 Bead comment reports a 24.1% scheduled-operation reduction and
directional 1.7–19.3% prior-Leopard2 gains. The low-copy Bead comment reports
directional gains of roughly 3–23% and a -0.17% high-rate control. No committed
raw/summary bundle supporting those exact ranges was found in this audit, and
the isolated low-copy run has not completed its strict host-isolation gate.
These numbers must therefore remain labeled non-authoritative and provisional;
they do not complete the respective promotion criteria. Final current-source
exact-main and dispatch evidence remains `leopard-79h.16.2`, `.38`, `.38.1`,
and `.38.2`.

The full plan benchmark matrix is also open. The saturated all-K diagnostic is
not a replacement for isolated 64 B through 16 MiB, loss/reuse/batch/thread
crossovers, setup separation, PMU/cache/TLB or memory-bandwidth counters where
available, scratch and code/table footprint, and 1-to-128 scaling. Those gates
remain in `leopard-79h.16`, `.23`, and `.24`.

## Production gates and dispositions represented in Beads

| Requirement family | Durable Beads coverage |
|---|---|
| Dense active-locator permanent-contribution checkpoint | `leopard-79h.11.1`; sparse reuse is complete and dense caching has a tested negative disposition |
| Sparse high/low encoder output pruning | `leopard-79h.26.4` |
| Algorithm 5 current-source isolated crossover and aligned direct-reveal disposition | `leopard-79h.26.5`; compact tails retain the gather |
| Algorithm 4 C1 isolated crossover and promotion evidence | `leopard-79h.18.1.12`, `.18.1.12.1` |
| Provisional low-copy isolated final A/B | `leopard-79h.26.1`, `.26.1.1` |
| Production all-K existing-path dispatch and exact-main gaps | `leopard-79h.38`, `.38.1`, `.38.2`, `.38.4` |
| GF16 padded-odd and explicit encode-as-parity-rebuild API decision | `leopard-79h.39.1` |
| Required native NEON and Windows validation | `leopard-79h.13`, `.13.5` |
| Optional PMULL/SVE/SVE2 research | `leopard-79h.18.25` |
| Allocation-free executor setup, scalable batch alias metadata, static scheduling, affinity, 128-core and NUMA | `leopard-79h.14`, `.14.1`, `.14.2`, `.14.4`, `.14.5`, `.23`, `.24` |
| Full Leopard1-centric benchmark matrix, counters and footprint | `leopard-79h.16`, `.16.2`, `.16.4`, `.23`, `.24` |
| Full sanitizer/fuzz/compatibility and release gates | `leopard-79h.15`, `.17`, `.40` |
| Product-tree/epsilon locator comparison | `leopard-79h.29`, as an unpromoted experiment rather than a CPU-release blocker |
| Incremental locator correctness comparison | `leopard-79h.18.5` |
| Exact-size/profile and deterministic crossover research | `leopard-79h.18.1.8`, `.18.1.10`, `.18.1.11` |
| Remaining A-W experiment disposition, including negative results | `leopard-79h.18.26` and its open experiment children |
| Final bibliography/formula refresh and evidence report | `leopard-79h.25`, `.27` |
| Optional CUDA work after the CPU codec | `leopard-79h.19.3` through `.19.9`; disabled by default and documented in `README_CUDA.md` |

## Dependency hygiene and priority

The Beads graph was reconciled after `1962218` with the user's core-codec and
Leopard1-first priority. `bd dep cycles` reports no cycles:

- final exact-main refresh `leopard-79h.16.2` now waits for the production
  dispatcher and pruned/tiled work (`.38` and `.26`), not optional Jerasure or
  exact-profile research;
- CPU release `leopard-79h.17` now waits for the focused exact-main refresh
  `.16.2` and the full Leopard1-centric matrix `.16.4`, not the broad benchmark
  parent or locator experiment;
- production `.26` and `.18.1.12` use the narrow closed runner/backend
  prerequisites `.16.1`, `.36`, and `.13.2` rather than broad unfinished
  parents; and
- scheduler `.14`, scaling `.23`, formula gate `.27`, product-tree/epsilon
  `.29`, and exact C10 crossover work use the narrow artifacts they actually
  consume.

Jerasure is end-stage optional comparison work: `leopard-79h.16.3` is open at
P4 and no longer blocks exact-main or CPU release work. Product-tree/epsilon and
the exact low profile are open P3 experiments. Unknown-error correction remains
correctly isolated under `leopard-79h.18.14` and must not enter the production
erasure hot path.

## Audit conclusion

The high and low mathematical codecs, active-parent normalization, systematic
wire profiles, reusable locators/plans, direct fallback paths and side-sized
decode execution are present. The definition of done is not yet met. Remaining
production source work includes sparse encoder pruning, a complete existing-path
dispatcher, restoring exact-main backend throughput, native NEON, and
allocation-free/static/NUMA-aware batch execution. High output evaluation no
longer copies the T-wide
coefficient block; the retained common gather is a deliberate compatibility
boundary rather than an unimplemented fusion. Decision `leopard-79h.39.1` freezes odd application
payloads as the explicit padded-GF16 wire layout and parity rebuild as an
ordinary encode after original recovery; neither requires a duplicate hot path.

Production evidence also remains incomplete: current-HEAD compatibility and
sanitizer/fuzz repeats, authoritative exact-main and full benchmark matrices,
128-core/multi-NUMA execution, final formula/literature traceability, and release
documentation. Speculative exact-code, Jerasure, CUDA and other A-W experiments
must not block the core codec except where the plan explicitly places them after
CPU release.
