# Leopard2 benchmark matrix runner

Status: production experiment infrastructure. Benchmark results remain
machine-specific evidence and do not change a wire profile or dispatcher.

`tools/leopard2_benchmark_matrix.py` generates deterministic job specifications
for `tools/leopard2_lab.py` and collects each benchmark's JSON from the lab
runner's durable per-job stdout. This keeps the benchmark executable focused on
one code cell while the runner owns affinity, memory limits, timeouts, resume,
failure logs, and deterministic merging.

## Presets

- `smoke` has six jobs for one low- and one high-rate cell: automatic,
  forced-specialized, and forced-generic for each. It validates the pipeline,
  not performance.
- `checkpoint` has 63 jobs intended for isolated execution: low `(16,240)`,
  balanced `(128,128)`, and
  high `(240,16)`; 4 KiB, 64 KiB, and 1 MiB shards; zero, one, and eight losses;
  automatic no-loss rows and all three request modes for nonzero loss. It is
  the intended bounded local comparison when it is run with the isolation rules
  below. Merely completing this preset is not evidence that the host was
  isolated.
- `balanced-crossover` has 216 jobs for `(128,128)` high profile, 256 B through
  64 KiB, and one through 128 missing originals, with all three request modes.
  It is a focused dispatcher diagnostic; its existence is not evidence for a
  threshold.
- `required` currently has 2,483 jobs spanning the listed CPU base count groups, ten
  shard sizes from 64 B through 16 MiB, all required loss classes,
  forced-specialized/forced-generic pairs for nonzero loss, automatic no-loss
  rows, reuse/batch samples, and automatic 1-to-128-thread scaling samples.
  Every thread count for a rate case performs the same fixed work: batch 128,
  reuse 8, 4 KiB shards, and eight losses.
  It is designed for resumable execution across appropriate hosts. Cells whose
  memory footprint is unsuitable for a specific machine are recorded as
  `unavailable` before launch or can be run on a larger host; they must not be
  silently replaced by a different code. Per-job address-space limits use the
  full conservative estimate and are not capped at 64 GiB. The lab runner also
  limits the sum of concurrently active job estimates to 80% of the readable
  physical/cgroup capacity. Preflight-unavailable cells are counted separately,
  not as executed processes.

Automatic, forced-specialized, and forced-generic members use the same
deterministic benchmark seed and a shared CPU-assignment group, so non-pinned
manifest generation does not accidentally rotate them onto different cores.
The collector's primary comparison is between the two forced, therefore known,
paths. A reported `forced_specialized_speedup_vs_forced_generic` is
`forced-generic median / forced-specialized median`.

When all three modes exist, a separate dispatcher-check record reports the
automatic median, the best forced median in that run, and their ratio. It does
not infer which internal path automatic policy selected from timing alone.

Pairing is bound to the pair group scheduled in the manifest and requires
equality of the requested parameters and exact missing-index list; resolved
profile, field, backend, thread count, parent and padded side; benchmark build
identity; benchmark executable SHA-256; assigned CPU set; and seed. If two
successful members of a scheduled group emit different identities, collection
fails instead of silently treating them as unrelated cells. Duplicate members
are errors, as are non-finite or non-positive timing, stale result or manifest
digests, failed round trips, unknown schemas, and missing result files.

The lab manifest v2 copies the generator's source schema, digest, and metadata.
Each job records the resolved executable path, size, and content SHA-256 inside
its job digest. The runner verifies that identity at campaign start,
immediately before each launch, and after execution. Terminal results hash
stdout and stderr and carry their own result digest; resume, merge, and collect
all reject mutated evidence. Each matrix job also carries an expected-cell
object in its job digest, and collection compares every emitted request
parameter with it before pairing; the collector requires the complete expected
field set. Regenerate older v1 manifests; they cannot establish these identities
and are intentionally not accepted by this collector.

## Known matrix follow-up

The exact 2,483-row `required` preset is a reproducible CPU checkpoint, not yet
the final Definition-of-Done matrix. Bead `leopard-79h.16.1` owns the follow-up:
add larger GF16 low-rate analogues; explicit reuse 1/one-stripe-per-plan plus
8/64/1024 reuse with batch varied independently; automatic nonzero-loss
dispatcher rows; and counterbalanced forced-path order (or repeated AB/BA
runs). Its new exact cardinality and independently enumerated dimension set
must replace 2,483; the current count is not a compatibility constraint.

## Reproduction

Configure and build a Release benchmark first:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    cmake -S . -B build/release \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEO2_BUILD_TESTS=ON \
        -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$JOBS"

Generate a candidate local checkpoint with one timing process at a time. This
intentionally uses one worker: simultaneous cache- and bandwidth-sensitive jobs
would invalidate the comparison. Each benchmark process is still pinned by the
lab runner, and full-machine parallel work should resume after this phase.
Authority additionally requires externally confirming host isolation and
recording frequency/governor state; the generator records serial execution but
does not claim isolation. Select one allowed physical CPU from the process
affinity mask; CPU 0 below is an example, not a portable assumption.

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

Rerunning the lab command resumes completed jobs whose job digest still matches,
including the executable content identity.
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

The permanent self-test checks deterministic generation, the exact 63-, 216-,
and 2,483-job preset cardinalities, fixed-work thread scaling, uncapped large
memory estimates and unavailable preflight, executable content identity,
source-metadata preservation, current cgroup v1/v2 limiting-ancestor discovery,
three-mode seed/CPU grouping, forced-path and dispatcher-check calculations.
Negative cases exercise bad benchmark schemas,
stale manifest, job, result, and output identities, delayed executable
replacement, request-parameter mismatches, failed round trips, duplicate pair
members, scheduled-pair identity drift, incomplete expected cells, and zero and
non-finite timing. It runs through CTest as
`leopard2_benchmark_matrix_self_test`.

The benchmark executable separately records codec and decode-plan setup,
execution, setup amortized at the selected reuse count, input and generated or
repaired output throughput, median/MAD/minimum/maximum timing, selected
profile/field/backend, scratch, legacy availability, and round-trip status.
Logical operation counts are provided independently by
`tools/leopard2_operation_counts.py`; neither tool presents estimates as
hardware counters.
