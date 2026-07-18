# All-K exact-main gap diagnostic

This is a diagnostic checkpoint for `leopard-79h.38.2`, not promotion
evidence. It compares exact Leopard main
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` with current Leopard2
`2fce390c9855b6c86b7e20fa86625db500757859`. Thirty independent workers
intentionally saturated every allowed logical CPU, so individual timing
outliers must not set production thresholds.

The independently linked binaries were:

- exact-main archive `8dc99f69...d5d4e` and executable `93d8c6b2...1eee30`;
- current archive `28b8e96c...d5d4e` and executable `fe569e1a...91526`.

The matrix completed 2,760/2,760 cells and 11,040 child invocations. It covers
every GF8 `K=1..255`, three legal redundancy bands, 4 KiB and 64 KiB shards,
and one/max loss, plus 238 representative GF16 cells at 512 B and 4 KiB.
Every main/current invocation produced the same original, parity, and recovered
digests; every round trip and current-tree legacy comparison passed.

Speedup below is left-hand time divided by right-hand time, so values above one
favor the right-hand implementation.

| Ratio | Encode | Decode first use | Decode at reuse |
| --- | ---: | ---: | ---: |
| exact main / Leopard2, median | 0.852 | 1.287 | 1.303 |
| exact main / current-tree legacy, median | 0.891 | 0.926 | — |
| current-tree legacy / Leopard2, median | 0.970 | 1.330 | 1.337 |

Leopard2 failed the 1.05 exact-main threshold in 2,258 encode cells (81.8%),
1,089 first-use decode cells (39.5%), and 1,078 reused-plan decode cells
(39.1%). The attribution is important: most broad encode loss already exists
between exact main's whole-build native code and the current portable-core
legacy API. Leopard2 adds a smaller median encode/API cost, while its planned
decoder wins substantially in the median.

## Ranked gaps

1. Restore coarse native-kernel performance without raising the archive ISA
   floor. Exact main/current legacy is 0.891 encode and 0.926 decode. A causal
   build should test moving complete FFT/IFFT/encode/decode traversals into
   ISA-specific translation units, reducing hot callback granularity.
2. Generalize balanced full-loss dispatch. AUTO's generic fallback is restricted
   to exact `K=R=128`, but same-seed, same-core ABBA diagnostics at
   `K=15,16,30,31,32` find generic 13.3% to 41.7% faster than specialized at
   4 KiB and 64 KiB. Generic still trails current legacy by 17% to 27%, so
   dispatch and residual generic/backend work are separate gaps.
3. Make the `R=1` XOR path genuinely lean. At GF8/4 KiB, current legacy / L2 is
   0.821 encode and 0.613 first-use decode; at GF16/512 B it is 0.367 and 0.128.
   The direct path still pays general range sorting, pointer setup, and
   transform-shaped scratch geometry. It needs direct scratch and linear safe
   validation before a kernel rewrite.
4. Reduce common encode API overhead. Current legacy / L2 encode is 0.970 in
   the overall median and 0.920 for GF16/512 B. Merge range construction with
   pointer population and test a monotonic-range fast path with a safe fallback.
5. Treat high-`T` full recovery separately from Algorithm 5 partial recovery.
   Full output removes message-only pruning; at GF8 `T=128`, current legacy / L2
   is 0.751 at 4 KiB and 0.827 at 64 KiB for max loss.

## Independent dispatch audit

`ShouldUseBalancedGenericDecode` is intentionally evidence-scoped to GF8
`K=R=T=128`, `N=256`, all originals missing, 256 B through 1 MiB, and the
scalar/SSSE3/AVX2 backends. `docs/leopard2_balanced_dispatch.md` explicitly
says neighbors stayed specialized because they were not measured. This is not
a mathematical restriction.

AUTO plans already retain both execution metadata sets. In balanced full loss,
specialized scratch clamps `min(N, 2T+L)` to `N`, while generic also uses `N`,
so selecting generic does not increase full-shard scratch. Both paths consume
the same locator and coordinate selection and have identical wire behavior.

The smallest structural hypothesis worth validating, but not merging, is:

    legacy-high && GF8 && K == R && missing == K && N == 2*T && T >= 8
        && 256 <= rounded_bytes <= 1 MiB
        && backend in { scalar, SSSE3, AVX2 }

This excludes direct-repair `T<=4`, partial-loss Algorithm 5 wins, GF16, and
unmeasured backends. A pinned campaign must cover `K` values
`5,7,8,9,14,15,16,17,29,30,31,32,33,62,63,64,65,126,127,128`; byte sizes
64, 256, 4 KiB, 64 KiB, 1 MiB, and 1 MiB+64; all three x86 backends; loss
neighbors `1,4,K/2,K-1,K`; rate neighbors; reuse `1,8,64`; and GF16 boundary
controls. Compare exact main, current legacy, forced generic, forced specialized
materialized/tiled, and AUTO. Promotion requires a credible 5% target gain, no
unexplained 2% neighbor regression, selected-path threshold tests, byte equality,
scratch/wire gates, full CTest, ASan/UBSan, and TSan.

The existing `.38.2` acceptance criteria therefore remain open. This sweep is
not pinned, is AVX2-only, and covers only two sizes and one/max loss. It does
satisfy the requirement not to promote from saturated diagnostics and provides
the concrete matrix and root causes for the authoritative campaign.

The forced diagnostic used `K` order `15,16,30,31,32`, then byte order
`4096,65536`; cell index `i` used CPU `i` and seed `989855744+i`. Each cell ran
`--force-generic, --force-specialized, --force-specialized, --force-generic`
with otherwise identical arguments:

```sh
taskset -c CPU build/gap-current/bench_leopard2 \
  --k K --r K --bytes BYTES --loss K --batch 1 --reuse 16 \
  --iterations 9 --warmup 2 --threads 1 --seed SEED \
  --profile high --field auto --backend auto FORCE_FLAG \
  --retain-samples --json OUTPUT
```

## Reproduction

```sh
CURRENT_WT="${TMPDIR:-/tmp}/leopard-wt-allk-gap-map"
MAIN_WT="${TMPDIR:-/tmp}/leopard-wt-allk-main"

git worktree add --detach "$CURRENT_WT" \
  2fce390c9855b6c86b7e20fa86625db500757859
git worktree add --detach "$MAIN_WT" \
  6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198

cmake -S "$CURRENT_WT" -B "$CURRENT_WT/build/gap-current" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF \
  -DLEO2_BUILD_BENCHMARKS=ON -DLEO2_BUILD_FUZZERS=OFF \
  -DLEO2_ENABLE_CUDA=OFF -DLEO2_BACKEND_VARIANT=auto
cmake --build "$CURRENT_WT/build/gap-current" \
  -j 30 --target bench_leopard2

cmake -S "$CURRENT_WT/experiments/leopard2/main_compare" \
  -B "$CURRENT_WT/build/gap-main" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DLEOPARD_MAIN_SOURCE_DIR="$MAIN_WT"
cmake --build "$CURRENT_WT/build/gap-main" \
  -j 30 --target leopard_main_benchmark

cd "$CURRENT_WT"
python3 experiments/leopard2/main_compare/run_allk_gap.py \
  --main build/gap-main/leopard_main_benchmark \
  --current build/gap-current/bench_leopard2 \
  --output .research/leopard2/allk-gap-attribution-2fce390-vs-6e5725e \
  --workers 30 --timeout 120 --with-current-legacy
```

The concise machine-readable record is `summary.json`. Its ignored raw bundle
is `.research/leopard2/allk-gap-attribution-2fce390-vs-6e5725e`; the raw merged
cell SHA-256 is `69a2c3ee...afd8a` and the raw summary SHA-256 is
`c44a3003...3ddc6`.
