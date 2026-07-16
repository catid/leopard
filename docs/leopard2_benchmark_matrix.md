# Leopard2 benchmark matrix runner

Status: production experiment infrastructure. Benchmark results remain
machine-specific evidence and do not change a wire profile or dispatcher.

`tools/leopard2_benchmark_matrix.py` generates deterministic job specifications
for `tools/leopard2_lab.py` and collects each benchmark's JSON from the lab
runner's durable per-job stdout. This keeps the benchmark executable focused on
one code cell while the runner owns affinity, memory limits, timeouts, resume,
failure logs, and deterministic merging.

## Presets

- `smoke` has four specialized/generic jobs for one low- and one high-rate
  cell. It validates the pipeline, not performance.
- `checkpoint` has 45 isolated jobs: low `(16,240)`, balanced `(128,128)`, and
  high `(240,16)`; 4 KiB, 64 KiB, and 1 MiB shards; zero, one, and eight losses;
  plus paired generic-decoder measurements for nonzero loss. It is the bounded
  authoritative local comparison.
- `balanced-crossover` has 144 paired jobs around the only checkpoint
  regression: `(128,128)` high profile, 256 B through 64 KiB, and one through
  128 missing originals. It supplies evidence for a deterministic dispatch
  threshold rather than turning one noisy cell into policy.
- `required` has 2,483 jobs spanning every count group in the mission, ten
  shard sizes from 64 B through 16 MiB, all required loss classes, paired
  generic decoding, reuse/batch samples, and 1-to-128-thread scaling samples.
  It is designed for resumable execution across appropriate hosts. Cells whose
  memory footprint is unsuitable for a specific machine should remain recorded
  as unavailable or be run on a larger host; they must not be silently replaced
  by a different code.

Specialized and generic members of a pair use the same deterministic benchmark
seed. The collector refuses stale lab result digests, failed round trips,
unknown benchmark schemas, or missing result files. A reported speedup is
`generic median / specialized median` for the same `(K,R,profile,field,backend,
bytes,loss,batch,reuse,threads,seed)` tuple.

## Reproduction

Configure and build a Release benchmark first:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    cmake -S . -B build/release \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=ON \
        -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$JOBS"

Generate the authoritative local checkpoint with one timing process at a time.
This intentionally uses one worker: simultaneous cache- and bandwidth-sensitive
jobs would invalidate the comparison. Each benchmark process is still pinned by
the lab runner, and full-machine parallel work should resume after this phase.
Select one allowed physical CPU from the process affinity mask; CPU 0 below is
an example, not a portable assumption.

    mkdir -p results/leopard2/benchmark-checkpoint
    python3 tools/leopard2_benchmark_matrix.py generate \
        --benchmark "$PWD/build/release/bench_leopard2" \
        --preset checkpoint --workers 1 --pinned-cpu 0 \
        --output results/leopard2/benchmark-checkpoint/spec.json
    python3 tools/leopard2_lab.py manifest \
        --spec results/leopard2/benchmark-checkpoint/spec.json \
        --output results/leopard2/benchmark-checkpoint/manifest.json \
        --root "$PWD" --workers 1
    python3 tools/leopard2_lab.py run \
        --manifest results/leopard2/benchmark-checkpoint/manifest.json \
        --output-dir results/leopard2/benchmark-checkpoint/run \
        --workers 1
    python3 tools/leopard2_benchmark_matrix.py collect \
        --manifest results/leopard2/benchmark-checkpoint/manifest.json \
        --results-dir results/leopard2/benchmark-checkpoint/run \
        --output results/leopard2/benchmark-checkpoint/collected.json

Rerunning the lab command resumes completed jobs whose job digest still matches.
Use `--rerun-failed` only after retaining and understanding the original stderr.

Generate the larger portable specification with:

    python3 tools/leopard2_benchmark_matrix.py generate \
        --benchmark "$PWD/build/release/bench_leopard2" \
        --preset required --workers 1 \
        --output results/leopard2/benchmark-required/spec.json

The target host builds its own manifest so its allowed CPU set and topology are
recorded. Thread counts 64 and 128 intentionally require a host exposing those
CPUs. This 32-CPU, single-NUMA-node host cannot provide that evidence.

## Validation

The permanent self-test checks deterministic generation, the 45-, 144-, and
2,483-job preset cardinalities, count/loss boundaries, memory-limit floors,
collector schema and digest gates, specialized/generic seed pairing, and the
speedup calculation. It runs through CTest as
`leopard2_benchmark_matrix_self_test`.

The benchmark executable separately records codec and decode-plan setup,
execution, setup amortized at the selected reuse count, input and generated or
repaired output throughput, median/MAD/minimum/maximum timing, selected
profile/field/backend, scratch, legacy availability, and round-trip status.
Logical operation counts are provided independently by
`tools/leopard2_operation_counts.py`; neither tool presents estimates as
hardware counters.
