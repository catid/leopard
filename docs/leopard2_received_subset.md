# Leopard2 received-subset experiment

Status: experiment-F correctness and symbolic-cost checkpoint. Nothing in this
checkpoint changes the production decoder, default dispatcher, public ABI, or
wire format. No wall-clock performance measurement was performed.

## Question and coordinate model

When more than `K` public shards survive, a transform decoder needs only `K` of
them. The production baseline keeps every surviving systematic shard, adds the
lowest recovery indices needed to reach `K`, and treats the rest as virtual
erasures. Experiment F asks whether selecting a different received subset can
reduce active transform blocks and memory traffic without changing decoded
data.

The checkpoint uses the frozen parent maps rather than treating public shard
indices as parent coordinates:

- legacy high: `T=ceil_pow2(R)`, `N=ceil_pow2(K+T)`, recovery coordinates
  `0..R-1`, and original coordinates `T..T+K-1`;
- low: `P=ceil_pow2(K)`, `N=ceil_pow2(P+R)`, original coordinates `0..K-1`,
  and recovery coordinates `P..P+R-1`.

Blocks have width `T` or `P`, respectively. Shortened and punctured parent
coordinates are fixed profile state and are never selectable.

## Policies and objective

The program compares four deterministic policies:

1. `lowest_parent_coordinate` selects the lowest available parent coordinates.
2. `prefer_systematic` is the production rule: all surviving originals, then
   the lowest recovery indices required to reach `K`.
3. `block_aligned_greedy` repeatedly chooses the block contributing the most
   still-needed coordinates, breaking ties by the fused inverse-transform cost,
   prefix slots, and block index.
4. `exact_block_dp` performs dynamic programming over blocks. For a fixed
   number selected from one block, the lowest local coordinates minimize both
   the staged prefix and the deterministic coordinate tie-break. The DP is
   therefore exact for this objective.

The objective is lexicographic:

1. number of active transform blocks;
2. radix-2 butterfly equivalents in the current fused DIT inverse-transform
   prefix schedule;
3. sum of staged input-prefix slots;
4. lexicographically lowest coordinate tuple.

This ordering deliberately prioritizes eliminating whole block transforms.
It is a symbolic schedule model, not a CPU cost model. It excludes locator
setup, cache effects, SIMD occupancy, output-block pruning, shard byte length,
batch size, plan reuse, and the current full-workspace staging cost.

## Independent algebra proof

The test-only field is `GF(2^4)` with irreducible polynomial `x^4+x+1` and
coordinate basis `[1,2,4,8]`. Every valid positive `(K,R)` geometry whose
active parent fits 16 coordinates is included: 85 legacy-high cells and 85
low-profile cells.

For each cell, the program builds the shortened-parent generator twice:

- directly with Lagrange basis products over every parent source coordinate,
  retaining the first `K` columns and fixing the shortened suffix to zero;
- independently by inverting the complete source-point Vandermonde matrix in
  the monomial basis.

All 9,758 resulting generator entries agree. The program then enumerates every
`K`-row public subset, inverts it, and verifies the inverse times the selected
generator is the identity. This checks 110,868 subsets and 789,908 basis-vector
decodes. For `K<=2`, it additionally decodes every possible field message from
every subset, for 218,672 exhaustive message decodes. Because the maps are
linear, the basis proof covers every message for larger `K` as well.

Finally, it enumerates every public receive pattern with at least `K` survivors
and verifies that each policy returns exactly `K` available coordinates from a
subset already covered by the exhaustive MDS/decode proof. The totals are
758,926 receive patterns and 3,035,704 policy decode cases. This proves the
selection decision cannot change the decoded message in the tested field;
production GF8/GF16 differential tests remain required before integration.

## Symbolic result

Across all GF(2^4) receive patterns, exact block DP was strictly better than the
prefer-systematic baseline in 309,047 patterns and equal in 449,879. It never
lost the declared objective. Aggregate counts were:

| Policy | Active blocks | IFFT butterflies | Prefix slots |
| --- | ---: | ---: | ---: |
| prefer systematic | 1,382,377 | 8,268,647 | 4,892,850 |
| lowest parent coordinate | 1,368,453 | 8,188,829 | 4,900,615 |
| block-aligned greedy | 1,203,851 | 7,636,433 | 4,552,245 |
| exact block DP | 1,203,851 | 7,530,777 | 4,504,504 |

Relative to prefer-systematic, exact DP reduced the aggregate modeled counts by
12.91% active blocks, 8.92% butterflies, and 7.94% prefix slots. The split was:

| Profile | Patterns improved | Active blocks | Butterflies | Prefix slots |
| --- | ---: | ---: | ---: | ---: |
| legacy high | 52,520 / 163,796 | -4.70% | -5.60% | -3.13% |
| low | 256,527 / 595,130 | -15.65% | -11.21% | -10.64% |

These aggregate small-field ratios are prioritization evidence only. They are
not a throughput claim and do not meet the experiment's promotion threshold.
The greedy policy matched the exact objective in 702,601 of 758,926 patterns;
the remaining 56,325 patterns quantify the potential value of the DP setup.

## Reproduction

Run the focused invariants and DP-vs-brute-force tests:

    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v \
        experiments/leopard2/received_subset/test_received_subset.py

Validate retained coverage, combinatorial counts, aggregate consistency, and
the source-content binding:

    python3 experiments/leopard2/received_subset/received_subset.py validate

Run the exhaustive checkpoint using every allowed CPU, capped at 128:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    python3 experiments/leopard2/received_subset/received_subset.py run \
        --workers "$JOBS" \
        --output experiments/leopard2/received_subset/results/checkpoint.json

The retained result is
`experiments/leopard2/received_subset/results/checkpoint.json`, SHA-256
`aed08c1e646c574163c38e0a550ebbbcb1f1b6adce2124d09e2520b5dbe09f66`.
The allowed set on the evidence host contained 28 CPUs, so the primary Python
3.12 run used 28 workers. A second Python 3.13 run with seven workers produced a
byte-identical file. Neither run was timed as a benchmark.

## Promotion gaps

Before any production use, the experiment still needs:

- integration with immutable locator and transform schedules, without changing
  profile identity or decode bytes;
- large random GF8/GF16 differential tests against direct interpolation and the
  existing deterministic plan;
- setup-versus-execution measurements at multiple plan reuse counts, shard byte
  sizes, patterns, batches, backends, and thread counts;
- an authoritative isolated timing campaign demonstrating at least a credible
  10% improvement in a meaningful region, with no neighboring regression;
- a compact deterministic dispatcher rule and a decision whether DP setup cost
  is justified over the much simpler greedy policy.

Until those gates pass, `prefer_systematic` remains the production selection
rule.
