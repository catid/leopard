# C1 inverse-prefix encoder experiment

## Disposition

Production promotion is rejected. The parent-code-preserving inverse-prefix
candidate did not meet C1's 10% improvement rule and regressed neighboring
AVX2 cells by more than the permitted 2%. Leopard2 `AUTO` remains on the
mature encoder.

The experimental implementation is intentionally absent from the production
tree. It is preserved on the archival branch
`agent/sparse-input-prefix-inverse`, whose history contains measured commit
`67eb32957824fe2b9a2a51945a6e1980aaef73a1` and final evidence commit
`f010d2671da71df9c931409476029f66dd12f5fb`. This avoids permanent plan
storage, internal field plumbing, and an unreachable executor for a rejected
path while keeping the exact source available for a future experiment.

## Hypothesis and wire identity

The mature high and low encoders materialize the shortened suffix of one
inverse-transform block as zero. The candidate compiled an immutable C1
schedule that consumed the active source prefix directly, fused complete
four-source groups, copied at most three ragged sources, and treated the suffix
as structural zero.

The experiment did not change the field, coordinate map, transform
normalization, parity order, or wire profile. Low rate used the strict
`K < P` input prefix. Legacy high rate used only the partial final `T`-shard
message block, with shift `T + floor(K / T) * T`; complete blocks retained the
mature path.

## Correctness evidence

The candidate remained reachable only through the diagnostic forced-transform
path and the internal benchmark; production `AUTO` was unchanged. GCC 13.3
Release validation covered 25,407 pruned plans, 33,399 executions, and
99,581,751 compared bytes across scalar, SSSE3, and AVX2. It included every
strict prefix through size 64, larger GF8/GF16 boundaries, ragged prefixes,
arbitrary GF8 byte tails, staged GF16 tails, shared immutable-plan concurrency,
and direct systematic parity comparisons. Focused Clang 18 ASan, UBSan,
LeakSanitizer, and compatible TSan cases passed.

The first performance screen was discarded because its archive contained
atomic test-hook instrumentation. Those counters affected the mature and
candidate loops differently. The corrected run linked a tests-off production
archive and reported `library_test_hooks=false`; no number from the discarded
screen is used below.

## Pinned production-linked result

The corrected runner pinned one worker to CPU 4 on an AMD Ryzen 9 9950X3D,
with SMT sibling CPU 20 and competing work attested idle. Each cell used 15
rotated samples of 32 iterations, three warmups, and separately measured plan
setup. Candidate/control is elapsed time, so values below one favor the
candidate. Gain is `control / candidate - 1`.

| Field/backend | Profile `(K,R)` | Bytes | Candidate/control | Gain |
|---|---:|---:|---:|---:|
| GF8 scalar | low `(17,18)` | 64 KiB | 0.9158 | +9.19% |
| GF8 AVX2 | low `(33,34)` | 64 KiB | 1.0805 | -7.45% |
| GF8 AVX2 | low `(65,66)` | 1 MiB | 0.9665 | +3.47% |
| GF8 AVX2 | high `(33,32)` | 64 KiB | 1.0938 | -8.57% |
| GF16 AVX2 | low `(129,130)` | 64 KiB | 0.9968 | +0.32% |

The median execution gain was 0.316%. The best observed gain was 9.19%, below
the required 10%; the median at reuse one became a 0.119% regression after
setup amortization. The two material AVX2 regressions independently disqualify
promotion. This is a negative result rather than an inconclusive one.

## Retained provenance

The committed matrix embeds every timing sample, median, MAD, parity digest,
command, and job-stream hash. Its companion manifest supplies the source
fingerprint, executable hash, machine identity, cell-manifest hash, settings,
and isolation attestation that the matrix itself references by configuration
ID.

- measured source: `67eb32957824fe2b9a2a51945a6e1980aaef73a1`, clean;
- compiler: GCC 13.3.0, C++11, test hooks disabled;
- source fingerprint: `daaf3e46fb2e17e5e02c2a3baa8421bdfb7fbb74aad55879317cdc78960b980c`;
- executable SHA-256: `26ffabb38864b59c3328a2b6f213dda22e76e967cea573f19f7a78fc6f007221`;
- configuration ID: `5b9833df6ea20f61b2bc18529306b3eec530e1070d65c5390edf3810391eae77`;
- cell manifest SHA-256: `3d6f7ed92f243a2e765f1c2dbc4d8f8b3d8bd71f760238def4b33b0d98e4da67`;
- retained run manifest SHA-256: `7c695ffdd3f4b77cdcb92f179888c169aec30e26e9860580292ce8d787a86b33`;
- retained matrix SHA-256: `004c10e8dd19fa225e75a3011dd68c65f3dddfa688f96c9b839ff6f9dd61673b`.

Artifacts:

- `experiments/leopard2/non_power_of_two/c1/inverse_prefix_cells.json`;
- `experiments/leopard2/non_power_of_two/c1/results/inverse_prefix_pinned_amd_9950x3d_manifest.json`;
- `experiments/leopard2/non_power_of_two/c1/results/inverse_prefix_pinned_amd_9950x3d.json`.

The raw manifest contains historical absolute worktree paths as observed by
the runner. They are provenance, not portable build instructions. No secret or
credential is present.

## Reproduction

After the archival branch is published, reproduce the measured implementation
without adding it to production:

    git fetch origin agent/sparse-input-prefix-inverse
    git switch --detach 67eb32957824fe2b9a2a51945a6e1980aaef73a1
    cmake -S . -B build/inverse-prefix -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/inverse-prefix -j 128 \
      --target bench_leopard2_sparse_encode

Prepare a current isolation attestation for an allowed physical CPU and its
SMT sibling, then run:

    python3 tools/leopard2_sparse_encode_crossover.py pinned \
      --source . \
      --executable build/inverse-prefix/bench_leopard2_sparse_encode \
      --result-dir build/inverse-prefix-results \
      --cell-manifest \
        experiments/leopard2/non_power_of_two/c1/inverse_prefix_cells.json \
      --iterations 32 --samples 15 --warmups 3 --setup-iterations 8 \
      --reuse 1,8,64 --memory-mib 1024 --workers 1 --cpu CPU \
      --isolation-attestation ATTESTATION.json --no-resume

The next technically distinct hypothesis is a compact regular-subtransform
plus boundary schedule, and for high rate a direct sink into the root
accumulator. Either requires a new correctness and neighboring-cell gate; this
rejected flat executor is not a production starting point.
