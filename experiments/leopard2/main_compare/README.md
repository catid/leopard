# Exact Leopard master benchmark adapter

This standalone project builds the legacy codec only from detached Leopard
revision `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. It does not link the
Leopard2 checkout's library and it is not included by the repository's normal
CMake build.

The adapter mirrors `bench/leopard2/benchmark.cpp` where the APIs overlap:
deterministic XorShift64 input generation and missing-shard selection, aligned
storage allocated before timing, every transmitted parity shard supplied to
decode, batch/reuse/warmup/iteration semantics, and median/MAD/minimum/maximum
summaries. Correctness and workload digests are computed before warmup.
Legacy decode necessarily includes its locator and other pattern-dependent
setup on every call; the JSON records that fact explicitly.

Workload fingerprints use FNV-1a-64 with offset
`14695981039346656037` and prime `1099511628211`, rendered as 16 lowercase
hexadecimal digits. Traversal is stripe index, shard index, then byte. The
`recovered_originals` shard order is the sorted
`missing_original_indices` list. The fields are:

    workload_digests.algorithm
    workload_digests.original_data
    workload_digests.transmitted_parity
    workload_digests.recovered_originals

Every encode and decode iteration is retained chronologically under
`samples_us_per_batch_call`; summary values are derived from those arrays.

## Build and smoke test

Create a detached baseline outside both the Leopard2 source and build trees:

    git worktree add --detach /tmp/leopard-main-exact \
        6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198

Configure the standalone project. Its configure step rejects the wrong commit,
a symbolic branch checkout, or tracked modifications. The exact master sources
form a private static archive; only that archive is linked into the adapter.
The build retains `compile_commands.json`, and its smoke test requires every
codec and adapter translation unit to carry the canonical master's effective
Release flags: `-march=native -Wall -Wextra -fopenmp -g -O0 -O3`, without
`-DNDEBUG`.

    cmake -S experiments/leopard2/main_compare \
        -B build/leopard-main-compare -DCMAKE_BUILD_TYPE=Release \
        -DLEOPARD_MAIN_SOURCE_DIR=/tmp/leopard-main-exact
    cmake --build build/leopard-main-compare -j "$(nproc)"
    ctest --test-dir build/leopard-main-compare --output-on-failure

For a single JSON cell:

    build/leopard-main-compare/leopard_main_benchmark \
        --k 240 --r 16 --bytes 64KiB --loss 3 --batch 1 \
        --reuse 8 --iterations 9 --warmup 2 --threads 1 \
        --seed 1 --json -

## Counterbalanced comparison

`run_abba.py` supplies the other half of the comparison. It accepts only fresh
Release artifacts with the expected compile, object, archive, link, source,
runtime-library, and field/profile identities. Each cell runs three independent
`baseline,candidate,candidate,baseline` rounds on one pinned logical CPU while
the sibling is reserved. Version-2 evidence additionally holds a pair-wide,
per-user lock at
`/run/user/UID/leopard2-cpu-leases/leopard2-cpu-pair-UID-A-B.lock`, derived
from both sorted logical CPU numbers. It validates the owned runtime and lease
directories, binds the open file's device and inode to its path, and rechecks
that identity after every child. The runner samples both CPUs from `/proc/stat`
immediately around the child sequence, and fails closed unless the benchmark
CPU did work and the reserved sibling accumulated zero non-idle jiffies. The
pair-wide lock prevents cooperating runners that implement this protocol from
using the physical core through different coordinator reservation files. A
coordinator must still exclude older or unrelated tools.
Both implementations must report the same deterministic loss set and all three
workload digests. A strict child environment prevents ambient OpenMP, loader,
allocator, or profiling settings from silently changing one executable.

Build Leopard2 separately with production test hooks disabled:

    cmake -S . -B /tmp/leopard2-production \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /tmp/leopard2-production -j "$(nproc)" \
        --target bench_leopard2

Choose a physical core from the allowed CPU set and reserve both of its SMT
threads. The reservation must be canonical JSON without a trailing newline;
replace the CPU numbers, owner, and nonce below. Keep unrelated user work off
both logical CPUs for the entire run. The runner is unprivileged: it cannot
prevent a root-owned task or kernel thread from being scheduled there, and
Linux scheduler counters have jiffy resolution. It therefore records and
checks the observed sibling counters rather than claiming OS-exclusive CPU
ownership.
The collector also holds the same stable
`/run/user/UID` directory lock as the current C7 and backend-butterfly runners
for the whole campaign, so replacing this reservation file does not split their
exclusive authority.

    python3 -c 'import json,sys; sys.stdout.write(json.dumps({"benchmark_cpu":15,"nonce":"replace-with-unique-value","owner":"benchmark coordinator","reserved_sibling":31,"schema":"leopard2-cpu-reservation/v1","status":"held"},sort_keys=True,separators=(",",":")))' \
        > build/leopard2-main-reservation.json

Run from an initial affinity that contains both reserved CPUs. The
`--preset representative` matrix covers exact-wire GF8 and GF16 high-rate, balanced, XOR,
field-inflation, and larger-parent cases. Custom cells use
`ID:K:R:BYTES:LOSSES:SEED` and must remain within the exact-main API's
`R <= K` and 64-byte-size restrictions.

    taskset -c 0-31 python3 \
        experiments/leopard2/main_compare/run_abba.py run \
        --baseline /tmp/leopard-main-compare/leopard_main_benchmark \
        --candidate /tmp/leopard2-production/bench_leopard2 \
        --baseline-archive /tmp/leopard-main-compare/libleopard_main_exact.a \
        --candidate-archive /tmp/leopard2-production/liblibleopard.a \
        --baseline-build-dir /tmp/leopard-main-compare \
        --candidate-build-dir /tmp/leopard2-production \
        --baseline-source-root /tmp/leopard-main-exact \
        --candidate-source-root "$PWD" \
        --candidate-commit "$(git rev-parse HEAD)" \
        --reservation-file build/leopard2-main-reservation.json \
        --output /tmp/leopard2-vs-main --cpu 15 --reserved-sibling 31 \
        --preset representative --reuse 8 --iterations 9 --warmup 2

Verify while the exact build inputs still exist:

    python3 experiments/leopard2/main_compare/run_abba.py verify \
        --manifest /tmp/leopard2-vs-main/manifest.json

The version-2 verifier recomputes the pair-lock identity, every per-field CPU
counter delta, the zero-non-idle-sibling decision, workload identities, and all
statistics. Version-1 retained bundles remain structurally replayable for
historical comparison, but they do not acquire version-2 isolation semantics
retroactively. A sibling-activity or child failure produces signed
`failure.json` diagnostics that bind every retained invocation stream. Verify
those diagnostics with:

    python3 experiments/leopard2/main_compare/run_abba.py verify-failure \
        --failure /tmp/leopard2-vs-main/failure.json

The verifier labels them explicitly as failed diagnostics, not performance
evidence. A failed campaign must not be converted into a passing run by
discarding samples.

The reported ratio is baseline time divided by Leopard2 time, so values above
one favor Leopard2. Encode excludes codec construction for both implementations.
Legacy decode inherently rebuilds its erasure setup on every call. Leopard2
`decode_first_use` is explicitly a derived median plan-create plus one execution
with an already-created codec; `decode_reuse_amortized` spreads that same plan
setup over the requested reuse count. Those components are measured in separate
loops and are not represented as a jointly timed call. Confidence intervals use
one clustered mean log contrast per ABBA round (three observations, Student-t
with two degrees of freedom), not six falsely independent adjacent pairs.
