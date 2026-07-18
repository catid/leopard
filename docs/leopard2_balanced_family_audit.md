# Leopard2 balanced full-recovery family audit

Status: source and correctness audit at `d9ef32e`.  No production dispatch
change is proposed from the retained diagnostic timings.  The current
singleton policy remains intact pending a fresh, isolated comparison of the
current generic, materialized Algorithm 5, and tiled Algorithm 5 kernels. The
authenticated external-matrix protocol for that comparison is implemented in
`leopard2_balanced_forced_evidence.md`; no broad authoritative timing has yet
been run with it.

## Outcome

The existing `ShouldUseBalancedGenericDecode` policy is an evidence boundary,
not an algebraic boundary.  It selects the generic active-parent decoder only
for legacy-high GF8 with `K=R=T=128`, `N=256`, all originals missing, rounded
shard bytes from 256 through 1 MiB, and a scalar, SSSE3, AVX2, or explicit
AVX512 backend.
Explicit force flags still override it.

For every balanced legacy-high GF8 code with `5 <= K=R <= 128`, let
`T=ceil_pow2(K)`.  The profile construction gives `N=2T`, with `T-K`
shortened systematic coordinates and the same number of punctured parity
coordinates.  If every original is missing, all `K` transmitted parity shards
must be present.  Both decoders therefore use the same selected coordinates,
locator, requested systematic prefix, and output scatter.  Choosing a decoder
cannot change parity, field representation, coordinates, or wire identity.

AUTO plan setup already prepares both execution paths.  The generic plan owns
its parent output dependencies and input prefix; Algorithm 5 owns block input
masks, output-block descriptors, exact pruned schedules, and fixed factors.
`UseGenericDecode` is deterministic in the immutable codec/plan and rounded
byte count, so the scratch query and execution make the same choice.

Algorithm 4 is outside this decision.  It remains the low-profile decoder and
must not be inferred from balanced legacy-high results.  Likewise, partial
recovery retains Algorithm 5's message-only output advantage and is a control,
not part of the candidate region.

## Source-level candidate matrix

The table uses the repository's ISA-independent radix-2 topology model.  `G`
is generic and `H` is high-rate Algorithm 5; `H/G` is scheduled butterfly
equivalents, not a latency prediction.  Non-power-of-two rows may also use the
new exact input/output pruning plans, so current-source timing is required even
where the structural ratio is favorable.

| Public balanced K range | T | N | G butterflies | H butterflies | H/G range | Evidence and disposition |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 1 | 2 | n/a | n/a | n/a | R=1 XOR path; exclude. |
| 2-4 | 2-4 | 4-8 | n/a | n/a | n/a | AUTO bounded direct repair handles full loss; exclude. |
| 5-8 | 8 | 16 | 48 | 51-52 | 1.063-1.083 | Structurally eligible; no current forced-path timing. |
| 9-14 | 16 | 32 | 97-103 | 120-128 | 1.225-1.243 | Structurally eligible; no current forced-path timing. |
| 15-16 | 16 | 32 | 104 | 128 | 1.231 | Old saturated AVX2 diagnostic favored generic; remeasure current code. |
| 17-29 | 32 | 64 | 232-256 | 309-335 | 1.309-1.336 | Structurally eligible; no current forced-path timing. |
| 30-32 | 32 | 64 | 256 | 335-336 | 1.309-1.313 | Old saturated AVX2 diagnostic favored generic; remeasure current code. |
| 33-64 | 64 | 128 | 461-544 | 680-768 | 1.412-1.475 | Structurally eligible; no current forced-path timing. |
| 65-127 | 128 | 256 | 1,064-1,280 | 1,629-1,856 | 1.449-1.532 | Structurally eligible; no current forced-path timing. |
| 128 | 128 | 256 | 1,280 | 1,856 | 1.450 | Historical isolated crossover produced the current singleton; requalify after kernel changes. |

The retained all-K artifact at
`experiments/leopard2/main_compare/results/2fce390-allk-diagnostic/` is useful
for prioritization but not promotion.  Its forced cells used AVX2 at 4 KiB and
64 KiB.  Generic was 13.3 to 41.7 percent faster than materialized specialized
decode for `K=15,16,30,31,32`, but the host was intentionally saturated.
More importantly, its current candidate was `2fce390`.  Since then the high
decoder gained exact input/output pruning and immutable evaluator execution,
the generic decoder gained fused reveal/scatter, and both paths gained coarse
backend traversal.  These changes affect the two sides differently.

The historical `K=128` isolated artifact also predates those changes.  It
remains the provenance for the existing conservative policy, but it cannot
justify widening the policy on current source without a new comparison.

## Scratch and correctness constraints

Let `B` be an aligned shard byte count.  For a ragged shard, let
`A=floor(B/64)*64` and `W=max(A,64)`.  For balanced full recovery `L=K` and
`N=2T`, so the current layouts reduce to:

| Path | Aligned data bytes | Ragged data bytes | Work slots |
| --- | ---: | ---: | ---: |
| Generic | `B*N` | `W*N + 128*K` | `N` |
| Normal/materialized Algorithm 5 | `B*N` | `W*N + 128*K` | `N` |
| Forced tiled Algorithm 5 | `B*(N+K)` | `W*(N+K) + 128*K` | `N+K` |

Generic and normal materialized execution also have the same `2N` pointer
entries and `3K` range entries.  Their complete scratch queries are therefore
identical, not merely asymptotically equal.  The pattern-independent one-shot
query substitutes `R` for `L`; because `R=K`, it also reserves exactly the
full-loss AUTO layout.  Widening AUTO inside this balanced/full-loss family
cannot under-allocate scratch.

Forced tiled execution is different.  It deliberately retains `K` output
slots in addition to two T-sized tiles and therefore uses more scratch than
either generic or materialized execution.  A performance campaign must compare
all three kernels.  The old forced diagnostic compared generic only with
materialized specialized decode; it cannot establish that generic beats the
current tiled/pruned Algorithm 5 implementation.

The new API test enumerates all 124 GF8 geometries `K=R=5..128`.  For every K
it checks AUTO, forced generic, forced materialized, and forced tiled recovery
at the ragged 193-byte boundary.  It also checks scratch at 193, 256, and 4,097
bytes, proves exact AUTO/generic/materialized equality, proves the extra tiled
capacity, and proves that the one-shot query covers the full-loss AUTO plan.
The focused API target passed in three independent configurations: the
Release AUTO backend, a GCC 13 scalar build with
`-Wall -Wextra -Wpedantic -Werror`, and a Clang 18 SSSE3 build with ASan,
UBSan, the same strict warnings, and fail-fast sanitizer settings.  The full
Release graph then passed all 69 CTest targets.  These are correctness and
build gates, not performance measurements; the existing backend matrix will
still exercise scalar, SSSE3, and AVX2 variants after integration.

No new independent algebra is needed to select between already-qualified
kernels, but promotion still requires generic/materialized/tiled byte equality
through every target and neighbor cell, legacy parity golden vectors, alias and
ragged-tail checks, full Release and sanitizer graphs, and concurrent immutable
plan execution.

## Exact measurement boundary

The smallest family worth measuring is:

    profile == legacy_high_v1
    field == GF8
    K == R
    missing_original_count == K
    N == 2*T
    T >= 8
    backend in { scalar, SSSE3, AVX2 }

This is a measurement domain, not a production predicate.  The deterministic
campaign must contain the following axes:

| Axis | Required values |
| --- | --- |
| K=R candidates | `5,7,8,9,14,15,16,17,29,30,31,32,33,62,63,64,65,126,127,128` |
| Actual shard bytes | `64,192,193,256,4032,4033,4096,4097,65536,1048576,1048577` |
| Missing originals | unique values from `1,4,floor(K/2),K-1,K` |
| Rate controls | `R=K-1`, `R=K`, and `R=max(1,floor(K/2))` |
| Backends | scalar, SSSE3, AVX2; NEON as an unpromoted control when available |
| Execution modes | forced generic, forced materialized, forced tiled, AUTO, current-tree legacy, exact main |
| Reuse | `1,8,64` with setup and execution reported separately |
| Batch | `1,8,64`, because the current generic decision does not inspect batch size |

The 192/193 pair is required because dispatch receives a 64-byte-rounded size:
193 bytes currently enters the 256-byte policy bucket even though only a
192-byte aligned prefix and one tail byte execute.  The 4,032/4,033 and
4,096/4,097 controls cover the separate generic reveal/scatter threshold and a
ragged tail.  The 1 MiB+1 control proves the upper rounded-size exclusion.

Forced codec flags omit metadata for the path they cannot execute, whereas an
AUTO plan owns both metadata sets and fixed factors.  Forced comparisons are
valid for byte-heavy execution but their setup times are not a substitute for
AUTO plan setup.  After deriving a proposed table, a candidate build must run
the actual AUTO path and verify its selected kernel, output, scratch, and
setup-amortized result.

Promotion should be by compact measured regions, not by the broad structural
hypothesis alone: the chosen path needs a 95 percent confidence lower bound of
at least 1.05 over both alternatives, no unexplained neighboring regression
over two percent, and no loss against exact Leopard main that is incorrectly
attributed to dispatch.  Shared-host or saturated timings remain diagnostic.

## Tooling finding

Schema v2 of `tools/leopard2_operation_counts.py` removes the historical
`api_scratch_data_slots=K+R+N` metric.  It now mirrors `DecodeLayout`,
`DirectDecodeLayout`, and `ComputeSplitScratchLayout` byte-for-byte for a
declared pointer width: aligned caller inputs remain pointer mappings, ragged
transform input reserves K+R fixed 64-byte slots while populating exactly K,
and work is reported separately as N, 2P, or 2T+L slots.  A production-linked
probe cross-checks plan and one-shot codec queries, including AUTO rows labeled
with the host-selected backend.  Those exact scratch metrics may be used in
this crossover; AUTO timing thresholds remain host/backend observations.
