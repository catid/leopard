# Leopard2 plan-to-production gap audit

This audit treats the tracked `plan.md` as the acceptance contract rather than
using historical issue closure as proof.  The plan checkpoint has SHA-256
`670ce14053ec528bea1cf72ba70b8fcd90eec7349bec03c80332adbb318ef753`.
The first complete source pass was made at `cdf2211` on 2026-07-18.  The durable
work items cited below are Beads issues; this document is an evidence map, not a
second task list.

## Production algorithm check

| Plan requirement | Current production evidence | Disposition |
|---|---|---|
| Legacy high wire profile | `leo2_codec_create` derives `T=ceil_pow2(R)`, `N=ceil_pow2(K+T)`, and the parity/message/shortened coordinate map. `ReedSolomonEncode` retains the block-IFFT accumulator and truncated final FFT. Legacy golden, API, arbitrary-count, GF16 matrix, and exact-main diagnostic tests compare parity bytes. | Implemented. |
| Low wire profile | Codec setup derives `P=ceil_pow2(K)`, `N=ceil_pow2(P+R)`, with the systematic prefix and punctured parity suffix. `ReedSolomonEncodeLow` performs one padded-P inverse transform and shifted evaluations, including `R>K`. | Implemented. Commit `bd681a7` removes the former whole-P coefficient copy for each parity block. |
| R10 Algorithm 4 | `ReedSolomonDecodeLowPlanned` and its tiled form perform P-point block inverse transforms, the active-parent derivative/weighted reduction, one P-point output transform, and original-only recovery from plan-owned locator factors. | Arithmetic and active-parent behavior are implemented and independently tested. Exact sparse C1 input/output schedules are not yet consumed; see `leopard-79h.18.1.12.1`. |
| R10 Algorithm 5 message-only form | `ReedSolomonDecodeHighPrunedPlanned` and its tiled form construct the T-wide accumulator, apply locator factors, invert to the evaluator, and evaluate only output blocks containing missing originals. Commit `e2ce390` adds immutable exact-mask C1 schedules for GF8/GF16 input and output transforms. | Implemented. A remaining whole-T copy before each output-block evaluation and unfused reveal/scatter are tracked by `leopard-79h.26.5`. |
| Active-parent derivation | `docs/leopard2_math_and_sources.md` defines active N, shifts, normalization, coordinate maps, shortening/puncturing, locator and reveal factors. Direct GF(2^4), GF8/GF16 transform, locator, arbitrary-count, high/low acceptance, and XDRS differential evidence exercise those identities. | No arithmetic or coordinate-map gap found in this pass. Final formula/source and literature refresh remains `leopard-79h.27`. |
| Independent fallback/oracle | The retained generic active-parent decoder remains selectable. Direct field, generator/interpolation, repair, transform-differential, GF(2^4), GF8 and GF16 tests do not use the optimized path as their only oracle. | Implemented. |

The specialized decoders take the codec's active parent `N`, not the field
order.  Codec setup precomputes permanent puncture/shortening erasure vectors
and caches their locator contribution in the sparse-direct region.  Dense
locator setup calls the active-N Walsh construction, so the former unconditional
256/65,536-entry behavior is gone, but it still rebuilds the permanent component
with each dense pattern; `leopard-79h.11.1` owns that remaining separation.
The independently derived product-tree candidate remains an unpromoted
comparison.

## API, memory, and execution check

| Plan area | Verified behavior | Remaining boundary |
|---|---|---|
| Stable API and identity | `leopard2.h` exposes immutable contexts, codecs and decode plans; profile, field, layout, parent and padded-side introspection; scratch queries; one-shot/reusable decode; batch calls; result strings; and per-context scalar/SSSE3/AVX2 selection. Old `leopard.h` remains intact. AUTO is deterministic and never selects a CPU-dependent wire profile. | Serialized identity remains correctly experimental under `leopard-79h.18.23`. |
| Allocation and aliases | Scalar encode and plan execution use caller scratch and perform no allocation. Inputs may alias inputs; outputs and the full supplied scratch range are rejected on overlap. No-loss decode is a true zero-scratch no-op. R=1 and K=1 direct paths use range-only scratch after `8f24962`. | A context pool can lazily create/grow workers from the first larger batch and no caller-executor adapter exists. The strict setup-only allocation guarantee is tracked by `leopard-79h.14.2`. |
| Arbitrary byte lengths | GF8 accepts all positive lengths. Native GF16 accepts complete two-byte symbols; the versioned padded-odd layout plus pack/unpack helpers supports every positive application payload without silently changing native wire identity. SIMD and ragged 64-byte tails are tested. | This is a documented framing refinement of the plan's “arbitrary byte length” requirement, not silent native-GF16 reinterpretation. |
| Scratch target | Specialized decode uses `min(N,2P)` low slots or `min(N,2T+L)` high slots where profitable; aligned caller inputs are referenced directly and only a ragged 64-byte tail stages public coordinates. Direct repair has no shard-data scratch. Encode reuses 2P/2T work and bounds tail staging. | High Algorithm 5 still copies T shards for each requested output block (`leopard-79h.26.5`). Further encoder/output pruning is `leopard-79h.26.4`. |
| Parity subsets and rebuild | Null recovery outputs are never written. Missing parity is not rebuilt during decode; callers explicitly invoke encode after combining surviving and restored originals, and high/low acceptance tests compare rebuilt parity byte-for-byte. | Sparse masks can still compute an unused prefix inside an encoder transform block. Exact allocation-free pruning is `leopard-79h.26.4`. |
| Runtime backends | Baseline objects are portable; SSSE3 and AVX2 live in qualified members behind immutable context ops tables and startup KATs. Different contexts can select lower qualified backends concurrently. Compile-time GF8/GF16 options and default CUDA absence are tested. | Native NEON remains open in `leopard-79h.13`/`.18.25`; Windows compiler isolation is `leopard-79h.13.5`; exact-main throughput lost by out-of-line portable granularity is `leopard-79h.38.1`. |
| Batch and multicore | Batch encode/decode use a persistent context pool after lazy start, deterministic error selection, per-item scratch and immutable shared objects. The lab runner records affinity/topology, stable seeds, resumable jobs and per-job output. | Affinity-aware pools, caller execution, final 1..128 scaling and multi-NUMA evidence are `leopard-79h.14.1`, `.14.2`, `.23`, and `.24`. This host currently exposes 30 allowed logical CPUs and one NUMA node, so it cannot prove the external 128-core/NUMA gates. |

The public “optional parity rebuild” requirement is satisfied as an explicit
encode operation rather than a second arithmetic API.  This avoids putting
parity-only derivative work into decode and is exercised by the high, low,
arbitrary-count, and GF16 padded-layout tests.

## Correctness and evidence check

The production graph contains separate direct-field/interpolation, active-LCH,
locator, generator-row, generic-decode and optimized-codec tests.  The
test-only GF(2^4) programs exhaust field identities, active-parent
normalization, locator subsets, MDS/rank properties and small repair profiles;
GF8/GF16 differential suites cover coordinate boundaries, tails, requested
parity subsets, virtual erasures, repeated plans and concurrent execution.
`docs/leopard2_xdrs_differential.md` records the pinned XDRS algorithm and
coordinate conversion instead of treating research defaults as an oracle.

Compact golden vectors and exact-main adapter digests protect legacy parity and
recovery.  Compile-time GF8-only/GF16-only matrices verify that omission does
not silently select a different canonical field.  The operation-count model's
331 self-checks distinguish butterflies, fixed multiplications, XORs and byte
traffic from ISA changes; the new low-copy source-immutability KAT prevents a
whole-P copy from returning unnoticed.

The integrated Release graph passed 65/65 tests after the low-copy and high-C1
changes, and focused ASan+UBSan passed 17/17.  These are milestone gates rather
than the final campaign: the repeated manual/static/sanitizer/fuzz passes are
owned by `leopard-79h.40`, while the full release correctness and compatibility
gate remains `leopard-79h.15`.  MSan, native ARM runtime, external 128-CPU/NUMA,
and the final 128-seed campaign must not be inferred from local x86 evidence.

## Exact Leopard1 performance check

The reproducible runner and bounded results in
`experiments/leopard2/main_compare/results/2fce390-allk-diagnostic/` compare
against the exact main-branch codec at commit `6e5725e`.  Across 2,760 total
cells (every GF8 K plus 238 representative GF16 cells) and 11,040 invocations,
all round trips and comparable digests match.
The diagnostic separates exact main, the current-tree retained legacy entry
points, and Leopard2 so build/portable-backend loss is not mislabeled as an
Algorithm 4/5 cost.

The current evidence does not support a claim of universal speedup:

- current-tree legacy itself has median exact-main/current ratios of 0.891 for
  encode and 0.926 for decode, which isolates a portable/out-of-line backend
  throughput gap (`leopard-79h.38.1`);
- Leopard2 improves current-tree legacy decode at the aggregate median, but
  balanced `K=R`, full-loss cells select a specialized path that loses for
  every measured T band at or above eight (`leopard-79h.38.2`);
- the R=1 range/direct milestone reduces K=129, 4 KiB time by about 13% encode
  and 28% decode versus pre-change Leopard2, while exact main remains about
  7–8% faster because the backend gap is still present;
- main rejects `R>K`, so target low-rate Algorithm 4/low-encoder cells have no
  valid Leopard1 timing counterpart.  Forced low at the balanced boundary is
  not used as evidence for AUTO.

Commit `e2ce390` reduces padded-equivalent high decode schedule operations by
24.1%; audited directional cells improved 1.7–19.3% versus the prior Leopard2
without a measured neighbor regression.  Commit `bd681a7` improved audited
low encode cells by about 3–23% versus the prior Leopard2 and left a high-rate
control unchanged.  Final isolated exact-main refresh and promotion tables
remain `leopard-79h.16.2`, `.38.1`, `.38.2`, and `.18.1.11`.

## Open production gates represented in Beads

| Requirement family | Durable Beads coverage |
|---|---|
| Exact-mask Algorithm 4 production pruning | `leopard-79h.18.1.12.1` |
| Dense active-locator permanent-contribution reuse | `leopard-79h.11.1` |
| Sparse high/low encoder output pruning | `leopard-79h.26.4` |
| Algorithm 5 no-copy output evaluation and reveal/scatter fusion | `leopard-79h.26.5` |
| Low-copy isolated final A/B | `leopard-79h.26.1.1` |
| All-K exact-main backend and balanced/full-loss crossovers | `leopard-79h.38.1`, `.38.2` |
| Native NEON and remaining platform validation | `leopard-79h.13`, `.13.5`, `.18.25` |
| Allocation-free executor setup, affinity, 128-core and NUMA | `leopard-79h.14`, `.14.1`, `.14.2`, `.23`, `.24` |
| Full sanitizer/fuzz/compatibility and release gates | `leopard-79h.15`, `.17`, `.40` |
| Exact-size/profile and deterministic crossover research | `leopard-79h.18.1.8`, `.18.1.10`, `.18.1.11` |
| Remaining A-W experiment disposition, including negative results | `leopard-79h.18.26` and its open experiment children |
| Final bibliography/formula refresh and evidence report | `leopard-79h.25`, `.27` |
| Optional CUDA work after the CPU codec | `leopard-79h.19.3` through `.19.9`; disabled by default and documented in `README_CUDA.md` |

Jerasure work remains deliberately end-stage (`leopard-79h.16.3.2`, P4).  It
is not a dependency for implementing or validating Algorithm 4/5.  Unknown
error correction likewise stays isolated under `leopard-79h.18.14` and does
not enter the production erasure path.

## Audit conclusion

The high and low mathematical codecs, active-parent normalization, systematic
wire profiles, reusable locators/plans, direct fallback paths and side-sized
decode execution are present.  The principal gaps are now implementation
quality and selection gaps rather than missing Algorithm 4/5 identities:
exact-mask low decoding, no-copy/fused high output evaluation, sparse encoder
pruning, restoring exact-main backend throughput, broadening balanced dispatch,
and completing executor/platform/release evidence.  Repeated bug passes are
tracked separately by `leopard-79h.40`; their findings must be resolved before
this audit and the CPU release gate can close.
