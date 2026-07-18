# Exact sparse-output encoder crossover experiment

`bench_leopard2_sparse_encode` compares exact requested-output schedules with
the mature prefix forward transform. Both paths are linked into the same
binary, use the same field tables and selected backend, consume the same source
shards, and are checked byte-for-byte before and after timing. The benchmark is
diagnostic, test-hook-only, and is not installed.

The cell benchmark reports four measurements separately:

- exact schedule setup, including mask scanning and compilation;
- mature prefix transform execution;
- exact execution with the schedule already prepared;
- the current call-local exact total, which recompiles the schedule before
  each execution.

It also models setup amortization at explicit reuse counts. The model is
`prepared_execution + setup / reuse`; it is labeled as a model and does not
pretend that the current public encoder persists an output-mask plan. The
measured call-local total is the relevant comparison for the implementation as
it exists today.

The harness calls the existing GF8/GF16 high- and low-profile encoder kernels
directly. It therefore isolates schedule setup and byte-heavy transform work,
but does not include public-API validation, arbitrary-tail staging, batch
scheduling, or application copies. A promotion decision still requires the
ordinary end-to-end Leopard2 benchmark matrix after a kernel threshold is
selected.

## Matrix

`tools/leopard2_sparse_encode_crossover.py` supplies two sparse mask shapes for
each of:

- high GF8: `K=192, R=64`;
- low GF8: `K=32, R=192`;
- high GF16: `K=1000, R=200`;
- low GF16: `K=128, R=896`.

The screen defaults to 64-byte and 1-KiB shards. The pinned matrix defaults to
64 B, 1 KiB, 64 KiB, and 256 KiB. Reuse points default to 1, 8, and 64. The
backend list is selectable; every raw result records both the requested and
resolved backend. Each cell uses rotating three-way order for mature prefix,
prepared exact, and call-local exact samples.

Structural accounting includes full, mature-prefix, retained, one-output, and
fused-four butterfly counts plus schedule and dependency-workspace bytes. A
matrix cell fails if its sparse mask does not reduce the mature prefix graph.

## Evidence identity and authority

The runner fails closed:

- a deterministic manifest hashes every source file that can affect the
  experiment, the executable, settings, machine identity, and cells;
- compiler identity is reported by the cell and the runner retains the CMake
  cache hash plus selected build/field/backend flags when available;
- job IDs and merged ordering are deterministic;
- stdout and stderr are retained and hashed for every job;
- resume rejects stale configurations or modified artifacts;
- the source fingerprint, Git state, and executable hash are rechecked after
  the run;
- every benchmark result must report the Git SHA and dirty state embedded at
  CMake configure time;
- pinned mode rejects a dirty source tree and a binary whose embedded SHA does
  not equal the current clean source SHA;
- pinned mode requires one worker, `taskset`, an allowed CPU, and an explicit
  isolation attestation naming that CPU and affirming that its SMT sibling and
  competing work were idle.

The standalone C++ benchmark always says `authoritative: false`. Only the
runner may elevate a complete pinned matrix, and only after every identity and
isolation gate passes. Screen results are never authoritative.

An isolation attestation has this shape:

    {
      "schema": "leopard2-benchmark-isolation-attestation/v1",
      "cpu": 16,
      "smt_sibling_idle": true,
      "competing_work_idle": true,
      "operator": "operator name or automation identity",
      "timestamp_utc": "2026-07-18T12:00:00Z"
    }

The attestation is evidence, not a mechanism that makes the host idle. The
operator remains responsible for checking topology, sibling placement,
governor/turbo state, and concurrent workloads.

## Reproduction

Configure and build on the exact committed source that will be measured. A
post-commit reconfigure is required so the benchmark embeds the clean commit
identity.

    cmake -S . -B build/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release --target bench_leopard2_sparse_encode -j 128

Run the correctness/schema smoke and runner self-test:

    python3 tools/leopard2_sparse_encode_benchmark_json_test.py \
      build/release/bench_leopard2_sparse_encode
    python3 tools/leopard2_sparse_encode_crossover.py self-test

Run a parallel, non-authoritative screen:

    OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 \
    python3 tools/leopard2_sparse_encode_crossover.py screen \
      --source . \
      --executable build/release/bench_leopard2_sparse_encode \
      --result-dir results/leopard2/sparse-encode-screen \
      --backends scalar,ssse3,avx2 \
      --workers 128

Run an isolated pinned matrix after preparing the attestation:

    OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 \
    python3 tools/leopard2_sparse_encode_crossover.py pinned \
      --source . \
      --executable build/release/bench_leopard2_sparse_encode \
      --result-dir results/leopard2/sparse-encode-pinned \
      --backends avx2 \
      --cpu 16 \
      --workers 1 \
      --isolation-attestation /absolute/path/isolation.json

Revalidate and summarize retained artifacts without rerunning jobs:

    python3 tools/leopard2_sparse_encode_crossover.py analyze \
      --result-dir results/leopard2/sparse-encode-pinned

Generated result directories are ignored. Evidence selected for a project
checkpoint should be copied intentionally, accompanied by its manifest and
hashes, rather than committing an incidental screen.

## Promotion rule

Exact schedules remain outside `AUTO` until a clean pinned matrix shows at
least a credible 5% gain in a useful region, no unexplained regression greater
than 2% in neighboring cells, and the end-to-end public encoder confirms the
gain. Prepared-only wins do not justify promotion when call-local setup erases
them; they instead motivate a separate reusable encode-plan design.
