# Low-P16 partial-output campaign

`run_abba.py` measures the low-profile GF8/AVX2 P=16 partial-output
optimization from one frozen Leopard2 executable. Runtime mode `0` retains the
mature scratch/scatter path and mode `1` enables direct final-radix output, so
the comparison does not confound different builds. The exact Leopard main
adapter is included for the three target cells as descriptive context.

The low-v1 and Leopard-main legacy-high-v1 profiles have different parity wire
formats. The runner therefore requires parity equality between Leopard2 modes,
and original/recovered-data equality with main, but does not require low-v1
parity to equal legacy-high-v1 parity.

The default campaign uses CPU 13, reserves SMT sibling 29, retains 21 samples
after five warmups (4096 calls per sample), runs nine mirrored target rounds
and three mirrored neighbor rounds, and holds
`/tmp/leopard-gf8-authoritative.lock`. Pin the
coordinator away from the measured pair, for example:

    python3 experiments/leopard2/low_p16_partial_output/run_abba.py self-test

    taskset -c 12 python3 \
      experiments/leopard2/low_p16_partial_output/run_abba.py run \
      --leopard2 /absolute/path/to/bench_leopard2 \
      --main /absolute/path/to/leopard_main_benchmark \
      --source-commit "$(git rev-parse HEAD)" \
      --source-tree "$(git rev-parse 'HEAD^{tree}')" \
      --output /absolute/lane-owned/result-directory

Resume an interrupted run with the same arguments plus `--resume`; the
original executable paths may be omitted because the runner uses its frozen
artifacts. Replay all retained JSON, route, hash, source, timing-sample, digest,
and analysis checks with:

    python3 experiments/leopard2/low_p16_partial_output/run_abba.py verify \
      --output /absolute/lane-owned/result-directory

Exit status `0` means all promotion gates passed, `2` means the evidence was
valid but a performance gate failed, and `1` means evidence validation failed.
Exact-main results never affect promotion.

## Qualified result

The 2026-08-03 clean-source campaign passed all promotion gates: 270 accepted
processes, six objectively rejected and retried isolation samples, no digest
failure, no target failure, and no neighboring regression failure. Across the
three target cells, direct partial output improved reusable execution by
1.087x to 1.270x and public one-shot encoding by 1.091x to 1.238x. See
`results/promotion_20260803.json` for exact confidence intervals, artifact
hashes, correctness coverage, host policy, and the Leopard-main comparison.

The exact-main comparison remains deliberately descriptive. Explicit LOW_V1
was still slower than Leopard main on the three overlapping small `R <= K`
cells; the new API's AUTO policy uses LEGACY_HIGH_V1 in that regime. LOW_V1's
distinct capability is `R > K`, which the Leopard-main API rejects.
