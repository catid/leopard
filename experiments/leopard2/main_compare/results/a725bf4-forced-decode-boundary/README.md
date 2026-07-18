# Exact-main forced decode boundary checkpoint

This checkpoint validates schema-v4 forced decode attribution in the exact
Leopard-main ABBA runner.  It is deliberately narrow: aligned 64 KiB shards,
balanced GF8 legacy-high codes at `K=R=127` and `K=R=128`, every original
missing, batch one, reuse eight, and the production AUTO AVX2 backend.  It is a
harness qualification and a dispatcher lead, not the broad promotion matrix.

The independently linked baseline is exact Leopard main commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  The candidate is Leopard2 commit
`a725bf46846714c4f1e438e5a6a43c2e29b37613`.  Each mode ran three independent
`baseline,candidate,candidate,baseline` rounds.  The runner authenticated the
source, object, archive, executable, runtime-library, exact argv, force flags,
workload digests, and CPU-isolation closure before accepting a result.

## Result

Ratios are exact-main time divided by Leopard2 time, so values above one favor
Leopard2.

| Candidate mode | K | First-use ratio (95% CI) | Reuse-8 ratio (95% CI) |
| --- | ---: | ---: | ---: |
| AUTO | 127 | 0.83276 (0.83173--0.83379) | 0.83590 (0.83490--0.83691) |
| AUTO | 128 | 1.04909 (1.04727--1.05091) | 1.05198 (1.05015--1.05382) |
| generic | 127 | 1.05153 (1.04855--1.05452) | 1.05235 (1.04941--1.05530) |
| generic | 128 | 1.05109 (1.05053--1.05166) | 1.05189 (1.05134--1.05244) |
| materialized | 127 | 0.83512 (0.82916--0.84113) | 0.83808 (0.83208--0.84412) |
| materialized | 128 | 0.90479 (0.89723--0.91241) | 0.90675 (0.89917--0.91440) |
| tiled | 127 | 0.93231 (0.92688--0.93777) | 0.93595 (0.93054--0.94139) |
| tiled | 128 | 0.96709 (0.96036--0.97388) | 0.96932 (0.96253--0.97615) |

The source-linked introspection pass reports AUTO as
`materialized/workspace_materialized` at K=127 and
`generic/balanced_generic` at K=128.  Forced generic, materialized, and tiled
report `forced_generic`, `forced_materialized`, and `forced_tiled`
respectively.  All eight introspection runs recovered the exact original data
and produced identical per-cell data, parity, and recovery digests.

Generic is the only current kernel that beats exact main in both cells.  Its
K=127 geometric improvement is about 5.2 percent, but the 95-percent lower
bounds of 1.04855 and 1.04941 narrowly miss the default 1.05 promotion gate.
K=128 clears the gate.  The K=127 AUTO/materialized result identifies a real
dispatcher gap for the broader balanced crossover task; this checkpoint alone
does not define a new threshold.

AUTO, generic, and materialized used physical-core pair 4/20.  Tiled used
pair 14/30 after two pair-4/20 attempts were rejected.  Every accepted bundle
records zero non-idle jiffies on its reserved sibling.  Each mode's ABBA ratio
against exact main is accepted independently.  Direct cross-mode latency
claims remain directional because tiled ran on a different physical core.

## Fail-closed evidence

Four complete 24-invocation attempts were rejected before analysis:

- unshielded AUTO on 14/30: CPU 30 accumulated five user and one system jiffy;
- unshielded AUTO on 4/20: CPU 20 accumulated two user and two system jiffies;
- shielded tiled on 4/20: CPU 20 accumulated two user and two system jiffies;
- shielded tiled retry on 4/20: CPU 20 accumulated four user jiffies.

All four failure bundles pass `verify-failure`; none contributes a performance
number above.  During accepted runs the coordinator temporarily restricted
only current user-owned Codex-descendant agent/tool processes away from the
reserved pair.  It excluded kernel and SSH/session infrastructure, keyed every
saved affinity by PID plus process start time, restored the exact original mask
in an exit/signal trap, and retained the runner's zero-sibling-jiffy rule.

The complete accepted and rejected bundles are retained under
`/home/catid/leopard-results/`.  `summary.json` binds every bundle by file
count, byte count, deterministic file-tree digest, and manifest or failure
digest.  The compact checkpoint is committed instead of roughly 12 MiB of
per-invocation generated streams.  Copying a bundle to another path preserves
structural replay; current-input replay additionally requires the authenticated
build and source paths recorded by that bundle.

## Reproduction

Build the independently linked baseline and current candidate:

    cmake -S experiments/leopard2/main_compare \
        -B /home/catid/leopard-builds/forced-decode-main-a725bf4 \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEOPARD_MAIN_SOURCE_DIR=/home/catid/leopard-wt-main-baseline \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /home/catid/leopard-builds/forced-decode-main-a725bf4 -j4
    ctest --test-dir /home/catid/leopard-builds/forced-decode-main-a725bf4 \
        --output-on-failure

    cmake -S . \
        -B /home/catid/leopard-builds/forced-decode-candidate-a725bf4 \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /home/catid/leopard-builds/forced-decode-candidate-a725bf4 \
        -j26 --target bench_leopard2

Validate the runner:

    python3 -m py_compile \
        experiments/leopard2/main_compare/run_abba.py \
        experiments/leopard2/main_compare/test_run_abba.py
    python3 experiments/leopard2/main_compare/test_run_abba.py

After reserving both SMT threads of one physical core, run one mode by replacing
`MODE`, `CPU`, `SIBLING`, `RESERVATION`, and `OUTPUT`:

    taskset -c 0-14,16-30 python3 \
        experiments/leopard2/main_compare/run_abba.py run \
        --baseline /home/catid/leopard-builds/forced-decode-main-a725bf4/leopard_main_benchmark \
        --candidate /home/catid/leopard-builds/forced-decode-candidate-a725bf4/bench_leopard2 \
        --baseline-archive /home/catid/leopard-builds/forced-decode-main-a725bf4/libleopard_main_exact.a \
        --candidate-archive /home/catid/leopard-builds/forced-decode-candidate-a725bf4/libleopard.a \
        --baseline-build-dir /home/catid/leopard-builds/forced-decode-main-a725bf4 \
        --candidate-build-dir /home/catid/leopard-builds/forced-decode-candidate-a725bf4 \
        --baseline-source-root /home/catid/leopard-wt-main-baseline \
        --candidate-source-root /home/catid/leopard-worktrees/forced_decode_campaign \
        --candidate-commit a725bf46846714c4f1e438e5a6a43c2e29b37613 \
        --candidate-mode MODE --reservation-file RESERVATION \
        --output OUTPUT --cpu CPU --reserved-sibling SIBLING \
        --cell k127-full-64k:127:127:65536:127:107 \
        --cell k128-full-64k:128:128:65536:128:109 \
        --reuse 8 --iterations 9 --warmup 2 --timeout 300

Replay accepted or rejected bundles while their build inputs remain present:

    python3 experiments/leopard2/main_compare/run_abba.py verify \
        --manifest OUTPUT/manifest.json
    python3 experiments/leopard2/main_compare/run_abba.py verify-failure \
        --failure OUTPUT/failure.json
