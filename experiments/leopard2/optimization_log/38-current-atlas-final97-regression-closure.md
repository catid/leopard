# 38 — Current-atlas final97 regression closure

## Question and scope

After the bounded GF8 encoder and decoder fixes in reports 27–37, does current
Leopard2 still lose to the standalone Leopard main codec in any workload that
the final current-atlas regression matrix identified?

The answer is **no for this matrix**.  This is deliberately narrower than a
whole-codec or cross-machine claim.  The accepted campaign covers:

- GF8, `legacy_high_v1`, explicit AVX2;
- `R=32`, the exact 97 predeclared target and neighbor cells;
- `K=32..209` at the specific K values recorded in the checkpoint;
- 64- and 1024-byte shards, one-item batch, one thread;
- the committed deterministic one-, two-, approximately ten-percent, and
  maximum-loss erasure patterns for those cells;
- full-output encode execution and one-shot decode with pattern setup charged.

It does **not** establish a result for GF16, `LOW_V1`, other redundancy
counts, other byte sizes, larger batches, multicore execution, another backend,
or another processor.

## Accepted confirmation

The accepted source is commit
`970107e9c295ff3b4cec601836ba1df17dee7564`, tree
`0ea848d4a9eb5c56f4b85ec6ecf7c7bd784a0a39`.  Its frozen Leopard2
executable SHA-256 is
`a8d9fe801cfe1662e8548894d37b7c5fdabf9686abe0a8e2e4bca4ec68dfec4d`;
the linked archive is
`60ef83eeb4fbba7bafea09c5256acf577abefa3a2889c541af1802d8eb3e89d5`.
The comparator is the standalone exact Leopard main adapter at commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, executable SHA-256
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.
It is not the in-process legacy API that shares Leopard2's newly compiled
kernels.

The design used 75 round-level contrasts per cell: 25 repetitions of three
balanced six-process orders.  Each process used 31 timed iterations after 64
warmups.  CPU 13 was pinned, SMT sibling 29 was reserved, the candidate and
same-binary control shared one frozen path and inode, and the runner held both
the canonical campaign lock and the CPU-pair lease.  Confidence intervals are
two-sided Student intervals over the 75 round log contrasts
(`df=74`, `t=1.992543495180924`).  There was no trimming, pooling, cell
selection after observing outcomes, or statistical retry.

All 194 exact-main intervals—encode and decode for every cell—had a lower 95%
confidence bound at or above 1.0.  All 194 identical-binary control intervals
were wholly inside the predeclared `[1/1.02,1.02]` equivalence band.

| Result | exact main / current Leopard2, geometric mean [CI95] |
| --- | ---: |
| Weakest encode, all cells: K34/R32/B1024/loss1 neighbor | 1.033204 [1.031403,1.035008] |
| Weakest encode target: K35/R32/B1024/loss1 | 1.047110 [1.045990,1.048231] |
| Weakest decode: K95/R32/B64/loss32 | 1.056045 [1.054421,1.057671] |
| Median across 97 encode comparisons | 1.120843 |
| Median across 97 decode comparisons | 1.486131 |

The tightest same-binary interval was K47/R32/B64/loss1 encode:
1.005302 [0.996307,1.014377], still wholly inside the equivalence band.

The journal contains 7,275 accepted round artifacts and 43,650 accepted child
processes.  Five attempts (30 processes) were discarded before inference:
K77/B1024/loss8 round 70, K92/B64/loss1 round 73, K79/B1024/loss2 round 10,
K35/B64/loss32 round 71, and K53/B64/loss32 round 37.  Each observed exactly
one non-idle jiffy on reserved sibling 29.  Every accepted round observed zero
sibling activity, and every workload digest matched.  Independent full-journal,
artifact, and statistical replays were clean.

## Rejected evidence retained

The earlier 25-round v2 pilot is not pooled into the result.  It had no
exact-main failure, but its stricter same-binary negative-control gate rejected
three intervals, each containing 1.0 and each associated with a transient
candidate-labelled process plateau:

| Pilot failure | pilot control/candidate CI95 | fresh 75-round control/candidate CI95 | fresh exact-main/current CI95 |
| --- | ---: | ---: | ---: |
| K41/B64 encode | [0.973393,1.005860] | [0.996283,1.004088] | [1.114557,1.126111] |
| K36/B64 encode | [0.980113,1.002963] | [0.992262,1.002403] | [1.121348,1.136340] |
| K66/B64/loss32 companion encode | [0.979513,1.001706] | [0.998357,1.004342] | [1.143648,1.153306] |

The rejected pilot raw and summary hashes are
`5823b85c41b2f3cff081dfc1ee776b139a61354427afdfa294703081be3d4d1c`
and
`a4df8880ca7c70d892f6b0da0647792daed163e0316612beeb6e9f9bb529c29c`.
Its ordered-artifact hash is
`713b35a81cefcf59964699bbc5eed6f5bc364a8a4dec37fe09293db7b1bea2ad`.
It remains useful diagnostic context but contributes zero observations to the
accepted confirmation.

One still earlier launch failed before timing because
`CMAKE_CXX_COMPILER` had cache type `STRING`, while the production
provenance contract requires `FILEPATH`.  Its failure journal SHA-256 is
`ac5a3ee6856c743e0726185fad7df2fa23e73e1a4f4d13ffb962a6938fbeebed`;
zero cells were timed and no result was retained.

## Evidence

The complete ignored raw and summary inputs have SHA-256:

- raw:
  `bcebaf2a3f9a60e0d6d39ba27cc8a1ba0c6841390db0469cca676e1c55634dd7`;
- summary:
  `27c5cdc7927300ecd3e6da94fbbd7b21c3da0df3cc11fe9b5323059bdaac259d`;
- ordered accepted/discarded artifact sequence:
  `39e6afd3f8a6cb0ce30d8ee339ff71b61179d07664dcb7377ea6a5bdfdeaaa25`;
- fixed matrix:
  `fe21e525ff084564af37a4e7b43f17df8e971f96a252199125e26d50e445da8b`;
- runner:
  `7cfed68b4cf47dc88c75e8f62a846550e3b55c96a5bc2fbf206d72987903077d`;
- effective build configuration:
  `4172a7bb940fccaf70994df229f0da0d813ef8d18e2b11977aa72f419872426b`.

The committed compact projection, including all 97 loss patterns and all 388
exact-main/control point estimates and intervals, is
`experiments/leopard2/gf8_high_encode/results/current_atlas_final97_confirmation_checkpoint_20260818.json`.
It is a checkpoint of the authenticated full bundle, not a substitute for its
round journals.

Both hardened collector suites passed under normal Python and `python -O`.
The checkpoint was independently compared cell-for-cell with the full summary,
and all 97 workload records and 388 intervals matched exactly.  The final
confirmation commits only evidence machinery; the focused C++ correctness,
strict-compiler, and sanitizer gates for the selected routes are recorded in
reports 35–37.

## Reproduction

Run from a clean checkout at the measured commit.  The explicit cache types are
part of the evidence contract.  The serial build is intentional: the runner
performs its own reproducible-build replay, and this avoids unnecessary peak
memory.

    git switch --detach 970107e9c295ff3b4cec601836ba1df17dee7564
    rm -rf build/current-atlas-final97-confirm-970107e-v1
    flock /tmp/leopard-gf8-authoritative.lock bash -c '
      cmake -S . -B build/current-atlas-final97-confirm-970107e-v1 \
        -G "Unix Makefiles" \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/cc \
        -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++ \
        -DCMAKE_MAKE_PROGRAM:FILEPATH=/usr/bin/gmake \
        -DCMAKE_EXPORT_COMPILE_COMMANDS:UNINITIALIZED=ON \
        -DLEO2_BACKEND_VARIANT:STRING=avx2 \
        -DLEO2_BUILD_BENCHMARKS:BOOL=ON \
        -DLEO2_BUILD_TESTS:BOOL=ON \
        -DLEO2_ENABLE_CUDA:BOOL=OFF \
        -DLEOPARD_ENABLE_GF8:BOOL=ON \
        -DLEOPARD_ENABLE_GF16:BOOL=ON \
        -DLEO2_FLAG_MAVX512F:UNINITIALIZED=FALSE \
        -DLEO2_FLAG_MAVX512BW:UNINITIALIZED=FALSE \
        -DLEO2_FLAG_MAVX512VL:UNINITIALIZED=FALSE &&
      cmake --build build/current-atlas-final97-confirm-970107e-v1 \
        --target bench_leopard2 -- -j1
    '
    sha256sum \
      build/current-atlas-final97-confirm-970107e-v1/bench_leopard2

The expected candidate digest is
`a8d9fe801cfe1662e8548894d37b7c5fdabf9686abe0a8e2e4bca4ec68dfec4d`.
Build the standalone pure-AVX2 Leopard main adapter as described in
`experiments/leopard2/main_compare/README.md`; do not substitute the
in-process legacy path.  The timing command below intentionally refuses a main
or candidate executable whose digest differs from the accepted identities.
The confirmation runner acquires the canonical lock itself, so do not wrap this
command in a second `flock`.

    rm -rf build/regression-current-atlas-final97-970107e-confirm75-v1
    /usr/bin/taskset -c 13 \
      experiments/leopard2/gf8_high_encode/run_current_atlas_final97_confirmation_abba.py \
      --candidate build/current-atlas-final97-confirm-970107e-v1/bench_leopard2 \
      --candidate-sha256 a8d9fe801cfe1662e8548894d37b7c5fdabf9686abe0a8e2e4bca4ec68dfec4d \
      --control build/current-atlas-final97-confirm-970107e-v1/bench_leopard2 \
      --control-sha256 a8d9fe801cfe1662e8548894d37b7c5fdabf9686abe0a8e2e4bca4ec68dfec4d \
      --main build/regression-k62r8-b64-fused-8306863-final25-v4/frozen/main \
      --main-sha256 78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93 \
      --source-commit 970107e9c295ff3b4cec601836ba1df17dee7564 \
      --source-tree 0ea848d4a9eb5c56f4b85ec6ecf7c7bd784a0a39 \
      --output build/regression-current-atlas-final97-970107e-confirm75-v1 \
      --cpu 13 --sibling 29 --iterations 31 --warmup 64 --rounds 75

On a different processor, use a fresh output and treat the result as a new
campaign rather than relabeling it as reproduction of the accepted 9950X3D
timings.

## Result

The finite current-atlas GF8/R32 regression set is closed on the measured AVX2
host: the weakest confirmed encode and setup-inclusive decode lower confidence
bounds are 1.0314x and 1.0544x over standalone exact Leopard main.  The result
does not erase the broader release matrix; it establishes that the specific
historical GF8 regressions targeted by reports 27–37 are no longer present.
