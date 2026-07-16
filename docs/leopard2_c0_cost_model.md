# Leopard2 C0 symbolic cost and dependency model

Status: completed experiment scaffold and prioritization evidence. This model
does not promote an exact-size wire profile or select a production dispatcher.

The C0 simulator is
`experiments/leopard2/non_power_of_two/c0/simulate.py`. It answers a narrow
question: which non-power-of-two experiments deserve implementation and
measurement first? It is an architecture-neutral systematic-encoding model,
not a prediction of cycles or throughput. All candidates still require direct
algebra tests, exhaustive MDS tests where applicable, wire-compatibility tests,
and end-to-end benchmarks.

## Code geometry modeled

For the legacy high profile, the simulator uses the frozen profile definition:

- `T = ceil_pow2(R)`;
- `N = ceil_pow2(K + T)`;
- `D = N - T`;
- `D - K` shortened message coordinates;
- `T - R` punctured parity coordinates.

Encoding is represented as `ceil(K/T)` inverse `T`-point transforms followed
by one forward `T`-point transform requesting only `R` outputs. A partial final
message block is an active input prefix.

For the low profile, it uses:

- `P = ceil_pow2(K)`;
- `N = ceil_pow2(P + R)`;
- `P - K` shortened message coordinates;
- `N - P - R` punctured parity coordinates.

Encoding is represented as one inverse `P`-point transform with `K` active
inputs followed by `ceil(R/P)` shifted forward transforms. The last parity
block requests only its transmitted output prefix.

The exact candidates instead use the smaller public side, `R` for high and `K`
for low. These are marked `new_exact_profile`. They are not legacy compatible
unless a future generator-matrix equivalence proof says otherwise.

## Metrics and score

Every count is per payload byte, so shard byte length is deliberately absent.
The matrix records:

- radix-2 butterfly equivalents;
- XOR byte operations;
- fixed-coefficient multiplications;
- byte loads and stores;
- temporary shard slots;
- irregular boundary operations;
- transform passes;
- copied and zero-filled shard bytes;
- the profile transform side and the candidate transform side (these differ for
  exact methods when the public smaller side is not a power of two);
- parent and transmitted lengths, parent inflation, and the field selected by
  each geometry;
- whether dyadic padding forces GF16 even though the exact transmitted code
  fits GF8.

The stable symbolic ranking score is:

    XORs + 2 * fixed multiplications + loads + stores
        + 4 * irregular boundary operations

Butterfly equivalents are not added to this score because their component
arithmetic and memory operations are already counted. The weights make result
ordering reproducible; they do not encode the throughput of any CPU. Memory
traffic, coefficient specialization, cache behavior, SIMD regularity, and
setup/execution amortization must be measured separately.

## Candidate definitions and confidence

| Method | Wire scope | Confidence and intended next stage |
| --- | --- | --- |
| `padded_dyadic_parent` | parent-preserving | Complete regular radix-2 network; exact structural count and baseline. |
| `prefix_pruned_parent` | parent-preserving | Exact structural DAG count after known-zero input pruning; C1 must apply actual skew-factor zero/one specialization. |
| `dependency_pruned_parent` | parent-preserving | Exact forward-live/backward-needed DAG intersection for input and output prefixes; C1 execution candidate. |
| `recursive_truncated_parent` | parent-preserving | Same exact arithmetic DAG as dependency pruning, with a lower-bound packed-frontier scratch estimate and explicit ragged-path penalty; C2 candidate. |
| `binary_dyadic_exact` | new exact profile | Full transforms on the binary decomposition of the exact side plus an optimistic linear cross-block join. This is a lower bound until C5 derives the joining map. |
| `full_block_tail_exact` | new exact profile | Largest full dyadic block plus direct tail and a conservative dense head/tail join; C4 candidate. |
| `direct_matrix_exact` | new exact profile | Exact dense systematic generator arithmetic, `K*R` products and `R*(K-1)` XORs; useful as a tiny-code and GF8-rescue comparison. |
| `subproduct_tree_exact` | new exact profile | Conservative `3*q*ceil(log2(q))^2` generic polynomial proxy; comparison/fallback, not an LCH count. |

The structural simulator treats a generic radix-2 butterfly as making both
outputs depend on both inputs. It intentionally does not guess which real LCH
skew factors are zero or one. That coefficient-aware refinement belongs in C1.

## Reproduction

Run the invariant suite:

    python3 experiments/leopard2/non_power_of_two/c0/simulate.py self-test

Run the retained compact 256-by-256 sweep on every allowed CPU, capped at 128
when necessary:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    python3 experiments/leopard2/non_power_of_two/c0/simulate.py compact \
        --workers "$JOBS" \
        --output-dir experiments/leopard2/non_power_of_two/c0/results

Generate full heat maps and the complete 1,048,576-row long-form matrix without
adding the large outputs to the repository:

    python3 experiments/leopard2/non_power_of_two/c0/simulate.py compact \
        --workers "$JOBS" --heatmaps --full-matrix \
        --output-dir /tmp/leopard2-c0-full

The compressed matrix is ordered by `K`, `R`, profile, then method. Gzip uses a
zero timestamp and an empty embedded filename. Worker scheduling and worker
count therefore do not change any mathematical artifact.

Retained compact files are:

- `results/summary.json`: aggregate winners, method means, gain thresholds, and
  top-gain cells for all 65,536 GF8-range public parameter pairs;
- `results/focus.csv`: all method rows for the cross product of
  `3,5,7,9,15,17,31,33,63,65,127,129,255`;
- `results/gf16_samples.json`: 21 representative larger parameter pairs,
  including high, low, balanced, odd, and just-over-power-of-two cases;
- `results/manifest.json`: sizes and SHA-256 hashes.

## Validation evidence

On 2026-07-16, the allowed CPU set contained 32 logical CPUs. The complete run
used all 32 workers. Phases that did not use 32 CPUs were deterministic output
validation and the independent 7-worker rerun; neither was a benchmark.

The invariant suite passed:

- 1,440 input/output-prefix dependency networks checked;
- every network through length 32 compared with an independently represented
  backward graph trace;
- boundedness and monotonicity spot checks through length 256;
- 4,896 complete cost rows checked across boundary and representative GF16
  counts;
- geometry, field-boundary, direct-matrix monotonicity, serialization, and
  repeated-worker-payload invariants passed.

The full sweep produced:

- 65,536 `(K,R)` pairs;
- 131,072 profile cells;
- 1,048,576 method rows;
- 16 gain heat maps, each 256 by 256;
- a 24,312,909-byte compressed long-form matrix containing 1,048,577 CSV lines
  including its header;
- elapsed 7.24 seconds, 100.15 user CPU seconds, and 86,108 KiB maximum
  coordinator RSS for the 32-worker run.

A second full run with 7 workers reproduced byte-identical summary, focus,
GF16 sample, all 16 heat-map, and compressed-matrix files. The full matrix hash
was:

    938d12d58a563c2a3843da54c384b269dd679881f188469099bdecf362a73e39

The retained compact hashes are in `results/manifest.json`; the principal
hashes are:

| File | SHA-256 |
| --- | --- |
| `summary.json` | `56176ecfe27473da175230cc2182c0b44f4de41f6e9330eddc430eedf9714e6d` |
| `focus.csv` | `538b7c539c5bb22dba4b976d0686256b4cf5c0d6e33ab98acabd7d064ebe9d4d` |
| `gf16_samples.json` | `7bd13398a297202b33c9499a992661adfef9c6e755e482e05f8209cefcd3ae98` |

## Sweep results

The winner counts are symbolic-score winners, not promoted dispatcher cells:

| Profile | Method | Winning cells |
| --- | --- | ---: |
| legacy high | padded dyadic parent | 1,435 |
| legacy high | prefix-pruned parent | 52,668 |
| legacy high | binary dyadic exact | 7,166 |
| legacy high | direct matrix exact | 4,267 |
| low | padded dyadic parent | 1,831 |
| low | prefix-pruned parent | 28,338 |
| low | dependency-pruned parent | 20,835 |
| low | binary dyadic exact | 9,651 |
| low | direct matrix exact | 4,881 |

The methods with the broadest modeled gains were the wire-compatible pruning
variants:

| Profile and method | Mean gain | Cells at least 10% better |
| --- | ---: | ---: |
| legacy high prefix-pruned | 1.190 | 38,940 |
| legacy high dependency-pruned | 1.190 | 38,833 |
| legacy high recursive truncated | 1.183 | 37,831 |
| low prefix-pruned | 1.066 | 20,540 |
| low dependency-pruned | 1.099 | 26,894 |
| low recursive truncated | 1.093 | 25,397 |
| legacy high binary dyadic exact | 0.902 | 11,832 |
| low binary dyadic exact | 0.901 | 11,704 |
| legacy high direct matrix exact | 0.573 | 5,559 |
| low direct matrix exact | 0.560 | 4,790 |

The exact candidates have poor global means because they target regions, not
the whole square. The optimistic binary candidate wins many cells just above a
power of two. For balanced counts, its modeled gain is 1.73 at `(17,17)`, 1.79
at `(33,33)`, 1.82 at `(65,65)`, and 1.85 at `(129,129)`. Those values are a
hypothesis for C5, not evidence that the required cross-block algebra is sparse.

Direct matrices dominate the smallest dense cells in this score: modeled gain
is 2.36 at `(3,3)`, 2.31 at `(5,5)`, and 1.81 at `(9,9)`. This supports a small
bounded direct encoder/repair crossover study rather than an unbounded dense
fallback.

Across the 256-by-256 square, 21,590 profile cells and 16,256 distinct public
parameter pairs fit GF8 by `K+R` but are pushed into GF16 by at least one dyadic
profile. There are 10,795 such cells in each profile. This is the strongest
non-throughput reason to prioritize the C6/D GF8 field-boundary work. Maximum
parent inflation was 3.938 at legacy-high `(K,R)=(1,129)`, where 130
transmitted coordinates use a 512-coordinate parent.

The generic subproduct proxy never reached a 10% modeled gain. Full-block plus
dense-tail reached that threshold in only 1,763 high and 1,819 low cells and
never beat the other candidates under this score. Neither result is a proof of
poor measured performance, but both lower implementation priority.

## Implementation priorities from C0

1. Implement C1 prefix and bidirectional dependency schedules against the
   existing padded transform. They preserve the parent code, have exact DAG
   counts, and cover the broadest modeled region. Validate every output
   byte-for-byte before measuring.
2. Implement C2 packed recursive/truncated execution only after C1 establishes
   real operation lists. Its arithmetic savings are similar, while its scratch,
   recursion, branch, and schedule-size estimates remain uncertain.
3. Run C6/D GF8 boundary rescue in parallel with C1/C2. The 16,256 affected
   public pairs justify direct GF8 and exact/truncated GF8 measurements even
   when their field-operation count is not minimal.
4. Retain bounded direct generator/matrix paths for tiny `K`, tiny `R`, and
   small byte counts. Measure setup separately and cap coefficient tables.
5. Before implementing optimized C5 kernels, derive and scalar-test the
   cross-block join for binary additive cosets. Kill the route if the map is
   dense; the current winning regions use an explicitly optimistic lower bound.
6. Restrict C4 full-block-plus-tail work to prioritized small-tail and
   field-boundary cells. Try direct, padded, and recursive tails with measured
   cross-block factors rather than one universal tail method.
7. Keep generic subproduct trees as an independent exact-size comparison and
   possible setup-heavy fallback, not a production hot-path priority.

## Limitations

- The current model covers systematic encoding. Decoder locator construction,
  erasure count, pattern shape, received-subset selection, plan reuse, and
  repaired-output scatter need separate models and measurements.
- A field multiplication is one symbolic unit in both GF8 and GF16. The field
  boundary flag is therefore more informative than the score when padding
  changes the field.
- Memory counts model logical shard-byte traffic. They do not predict cache
  lines, TLBs, write allocation, SIMD lane utilization, fusion, or NUMA traffic.
- Prefix masks model the frozen profile layouts. Sparse general masks may have
  different DAG shapes.
- `recursive_truncated_parent` has exact arithmetic counts but only a scratch
  lower bound. `binary_dyadic_exact`, `full_block_tail_exact`, and
  `subproduct_tree_exact` have labeled heuristic algebra costs.
- An exact candidate may recover data correctly and still be wire-incompatible.
  It requires a new serialized profile identifier unless parity equivalence is
  proved for every claimed parameter set.
- No symbolic winner should enter AUTO. Promotion still requires correctness,
  memory-traffic measurement, at least a credible 10% target-region gain, and
  a deterministic dispatcher that avoids neighboring regressions.
