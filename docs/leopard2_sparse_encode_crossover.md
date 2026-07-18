# Exact sparse-output encoder crossover experiment

`bench_leopard2_sparse_encode` compares exact requested-output schedules with
the mature prefix forward transform. Both internal kernel paths are linked into
the same binary, use the same field tables and selected backend, consume the
same source shards, and are checked byte-for-byte before and after timing. The
benchmark is
diagnostic, links the ordinary production `leopard` archive, and is not
installed. It is available when benchmarks are enabled even if the test suite
is disabled. This does not mean exact scheduling is wired into the production
public encoder: public scratch reservation, schedule compilation, and selection
remain `LEO2_ENABLE_TEST_HOOKS`-only. Production `AUTO` still selects the mature
prefix path.

The benchmark-v4 schema reports `build.library_test_hooks=false`. The CMake
target sets that marker only while linking production `leopard`; manual builds
default to `true`, and the pinned runner rejects a true or missing marker. This
matters because atomic test counters perturb the mature and candidate paths by
different amounts and can bias a crossover result.

The benchmark-v4 primary comparison is pairwise and setup-inclusive. Each cell
retains exactly three independent rounds in ABBA, BAAB, ABBA order. For each
round it computes the paired log-time contrast
`mean(log(prefix)) - mean(log(call_local))`. The runner estimates the mean
contrast across rounds with a two-sided 95% Student-t interval (`df=2`) and
exponentiates the estimate and interval to speedup/gain space. Endpoint medians,
gain medians, and MAD are descriptive only; pooled rotating three-way samples
and ratios of global medians are not promotion evidence.

Prepared execution remains a separate diagnostic block. The cell reports:

- exact schedule setup, including mask scanning and compilation;
- mature prefix transform execution;
- exact execution with the schedule already prepared;
- the current call-local exact total, which recompiles the schedule before
  each execution.

It also models setup amortization at explicit reuse counts. The model is
`prepared_execution + setup / reuse`; it is labeled as a model and does not
pretend that the current public encoder persists an output-mask plan. The
measured call-local total is the only candidate in the primary ABBA pair for the
implementation as it exists today. Prepared-only wins can motivate a separate
immutable encode-plan API but cannot promote call-local plumbing.

The harness calls the existing GF8/GF16 high- and low-profile encoder kernels
directly. Before and after timing it checks deterministic symbols from every
requested output against an independent direct-generator construction,
performs a production public-decode recovery using candidate parity, records
SHA-256 digests, verifies immutable inputs, and checks allocation guards and
applicable unrequested-output canaries.
Timing still isolates schedule setup and byte-heavy transform work and does not
include public-API validation, arbitrary-tail staging, batch scheduling, or
application copies. A promotion decision still requires the ordinary end-to-end
Leopard2 benchmark matrix after a kernel threshold is selected.

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
resolved backend. Sparse candidate masks are accompanied by dense/prefix
negative controls. Candidate cells must reduce the mature prefix graph;
controls are retained precisely so the runner can detect an implausible blanket
advantage.

Structural accounting includes full, mature-prefix, retained, one-output, and
fused-four butterfly counts plus schedule and dependency-workspace bytes. A
matrix cell fails if its sparse mask does not reduce the mature prefix graph.

## Evidence identity and authority

The runner fails closed:

- a deterministic manifest hashes every source file that can affect the
  experiment, the executable, settings, machine identity, and cells;
- compiler identity is reported by the cell and the runner retains and hashes
  the CMake cache, production archive/object inputs, and executable;
- job IDs and merged ordering are deterministic;
- stdout and stderr are retained and hashed for every job;
- resume rejects stale configurations, path escapes, modified artifacts, and
  artifacts from a different current root;
- source, archive/object, executable, and runtime identity are checked before
  every retained invocation and rechecked after the run;
- every benchmark result must report the Git SHA and dirty state embedded at
  CMake configure time;
- benchmark-v4 must report that its linked library has no test hooks;
- pinned mode rejects a dirty source tree and a binary whose embedded SHA does
  not equal the current clean source SHA;
- pinned mode requires one worker, an allowed CPU, reservation of that CPU and
  every topology-reported SMT sibling, runtime affinity proof from the cell,
  and an explicit isolation attestation naming the whole reserved sibling set
  and affirming that competing work was idle.

The standalone C++ benchmark always says `authoritative: false`. Only the
runner may elevate a complete pinned matrix, and only after every identity and
isolation gate passes. Screen results are never authoritative.

An isolation attestation has this shape:

    {
      "schema": "leopard2-benchmark-isolation-attestation/v1",
      "cpu": 16,
      "reserved_cpus": [16, 80],
      "smt_sibling_idle": true,
      "competing_work_idle": true,
      "operator": "operator name or automation identity",
      "timestamp_utc": "2026-07-18T12:00:00Z"
    }

The attestation is evidence, not a mechanism that makes the host idle. The
runner holds the project-wide abstract-socket plus `flock` CPU-pair lease,
constrains its own orchestration to housekeeping CPUs outside that pair, checks
sibling scheduler counters, requires exactly zero non-idle sibling jiffies, and
fails if the benchmark's runtime affinity or observed topology differs, but the
operator remains responsible for governor/turbo state and concurrent workloads.

## Reproduction

Configure and build on the exact committed source that will be measured. A
post-commit reconfigure is required so the benchmark embeds the clean commit
identity.

    cmake -S . -B build/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=OFF \
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

Do not run an authoritative sparse matrix while the weighted matched-layout
benchmark owns the isolation lease. After that campaign releases the host and
the experiment owner explicitly authorizes timing, run an isolated pinned
matrix with an attestation covering the complete sibling set:

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

Exact schedules remain outside `AUTO` until the lower endpoint of the
exponentiated 95% Student-t interval clears a 5% gain in a useful sparse region,
dense/prefix negative controls show no implausible blanket speedup, no
unexplained neighboring regression exceeds 2%, and the end-to-end public encoder
confirms the gain. Prepared-only wins do not justify promotion when call-local
setup erases them; they instead motivate a separate reusable encode-plan design.
