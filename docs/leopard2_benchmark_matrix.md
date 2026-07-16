# Leopard2 benchmark matrix runner

Status: production experiment infrastructure. Benchmark results remain
machine-specific evidence and do not change a wire profile or dispatcher.

`tools/leopard2_benchmark_matrix.py` generates deterministic job specifications
for `tools/leopard2_lab.py` and collects each benchmark's JSON from the lab
runner's durable per-job stdout. This keeps the benchmark executable focused on
one code cell while the runner owns affinity, memory limits, timeouts, resume,
failure logs, and deterministic merging.

Linux `perf stat` counters are an opt-in runner feature. The generator signs
the exact `perf` executable and requested event list into every job. Before
timing, the lab runner probes the same events on the job's assigned CPU set.
When an optional probe is denied or an event is unsupported, the benchmark
still runs and its result records `performance_counters.status=unavailable`
plus the probe command, exit code, and diagnostic; it never substitutes zero
counts. A successful probe wraps the benchmark, retains the delimiter-format
`perf-stat.txt`, hashes that file into the terminal result, and parses each
requested event as `counted`, `not-counted`, or `missing`. Collection rejects
mutated raw counter output, request/evidence mismatches, or paired jobs with
different counter availability. Use `--require-perf-counters` only for a
campaign where a missing PMU measurement should make a cell explicitly
`unavailable` rather than execute it without counters.

## Presets

- `smoke` has ten jobs for one low- and one high-rate cell. Each cell runs the
  counterbalanced forced-path repetitions plus one automatic row described
  below. It validates the pipeline, not performance.
- `checkpoint` has 99 jobs intended for isolated execution: low `(16,240)`,
  balanced `(128,128)`, and
  high `(240,16)`; 4 KiB, 64 KiB, and 1 MiB shards; zero, one, and eight losses;
  automatic no-loss rows and counterbalanced nonzero-loss comparisons. It is
  the intended bounded local comparison when it is run with the isolation rules
  below. Merely completing this preset is not evidence that the host was
  isolated.
- `balanced-crossover` has 360 jobs for `(128,128)` high profile, 256 B through
  64 KiB, and one through 128 missing originals, with counterbalanced forced
  comparisons and one automatic row per logical cell.
  It is a focused dispatcher diagnostic; its existence is not evidence for a
  threshold.
- `required` has exactly 7,134 jobs: 6,870 main count/size/loss/mode rows,
  240 reuse/batch rows, and 24 thread-scaling rows. It spans the listed CPU
  base count groups, ten shard sizes from 64 B through 16 MiB, all required
  loss classes, counterbalanced forced paths plus one automatic row for every
  nonzero loss, automatic no-loss rows, a full reuse-by-batch grid, and
  automatic 1-to-128-thread scaling samples.
  The larger GF16 low-rate analogues are `(512,1536)`, `(1600,2496)`, and
  `(2032,2064)`. Their low-profile padded sides are 512, 2048, and 2048; their
  parents are 2048, 8192, and 8192 coordinates respectively.
  Reuse and batch vary independently over `{1,8,64,1024}` for each checkpoint
  low/balanced/high case. Batch 1 and 8 use 64 KiB shards, batch 64 uses 4 KiB,
  and batch 1024 uses 256 B, bounding allocation while retaining a full 4x4
  grid. All five scheduled rows are present for each of those 48 logical
  cells, so this portion can take a long time even on a fast host.
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

For a nonzero-loss logical cell, `A` means forced-specialized and `B` means
forced-generic. Lexicographically sorted job IDs encode the exact slots:

    order-ab.slot00-forced-specialized
    order-ab.slot01-forced-generic
    order-ab.slot02-automatic
    order-ba.slot00-forced-generic
    order-ba.slot01-forced-specialized

Thus a serial runner executes adjacent `S/G/automatic` followed by `G/S`; the
automatic observation exists once and belongs to the AB comparison. AB and BA
use the same deterministic benchmark seed and shared CPU-assignment group.
Collection separates the two order trials while retaining that shared affinity.
The five rows also share a signed `resume_group`. Every terminal result records
an opaque `run_epoch` identifying one lab-runner invocation. A group resumes
only when all five results are complete and carry the same epoch; a missing row
or mixed epochs reschedule all five rows. This preserves the temporal meaning
of AB/BA after interruption while retaining resumability at logical-cell
boundaries. Jobs without `resume_group` retain the lab runner's ordinary
job-granular resume behavior.
For authoritative timing, generate with `--workers 1 --pinned-cpu <allowed-cpu>`
and also run the lab with `--workers 1`. Multiple workers preserve the repeated
trials but do not preserve their temporal ordering. The pinned CPU must be an
allowed physical CPU selected from the process affinity mask, and host isolation
must still be verified externally.

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
The larger low-rate analogue rows additionally require resolved GF16 with
`(padded_side,parent_count)` equal to `(512,2048)`, `(2048,8192)`, and
`(2048,8192)` respectively; collection rejects a different resolution.

The lab manifest v2 copies the generator's v3 source schema, digest, and metadata.
Each job records the resolved executable path, size, and content SHA-256 inside
its job digest. The runner verifies that identity at campaign start,
immediately before each launch, and after execution. Terminal results hash
stdout and stderr and carry their own result digest; resume, merge, and collect
all reject mutated evidence. Each matrix job also carries an expected-cell
object in its job digest, and collection compares every emitted request
parameter with it before pairing; the collector requires the complete expected
field set. Regenerate older v1/v2 benchmark specifications; they do not define
trial/order semantics and are intentionally not accepted by this v3 collector.

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
        --perf-stat /usr/bin/perf \
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

Rerunning the lab command resumes completed logical cells whose job digests and
executable content identities still match and whose terminal results share one
run epoch. If any counterbalanced row is absent or comes from another runner
invocation, the entire five-row cell is rerun. Existing result files are
overwritten only for that incomplete or mixed-epoch cell.
Use `--rerun-failed` only after retaining and understanding the original stderr.
The default counter request covers cycles, instructions, generic cache and
branch events, page/context events, and data-TLB loads/misses. Site-specific
uncore memory-bandwidth events are not portable; supply a comma-separated list
with `--perf-events` when the host PMU and event naming have been verified.
The runner does not require privileged installation or attempt to change
`perf_event_paranoid`.

The 2026-07-16 checkpoint on this host exercised the complete ten-job smoke
pipeline with `/usr/bin/perf` pinned to CPU 15. All ten benchmarks succeeded;
all ten counter records were explicitly unavailable because the kernel reports
`perf_event_paranoid=4` and the probe exited 255. The ignored, machine-specific
evidence is under
`results/leopard2/benchmark-perf-smoke-20260716/`. This is pipeline evidence,
not an authoritative performance result; other implementation workers were
active during the smoke run.

Generate the larger required specification with an allowed pinned CPU:

    python3 tools/leopard2_benchmark_matrix.py generate \
        --benchmark "$PWD/build/release/bench_leopard2" \
        --preset required --workers 1 --pinned-cpu 0 \
        --output results/leopard2/benchmark-required/spec.json

The target host builds its own manifest so its allowed CPU set and topology are
recorded. `--pinned-cpu` pins every single-thread comparison, including both
AB/BA repetitions, while the multi-thread scaling rows retain their requested
topology-aware CPU counts. Thread counts 64 and 128 intentionally require a
host exposing those CPUs. A smaller host cannot provide that evidence.

## Validation

The permanent self-test checks deterministic generation; the exact 99-, 360-,
and 7,134-job preset cardinalities; an independently enumerated 7,134-signature
dimension/mode `Counter`; exact lexicographic AB/BA slots; one automatic row;
shared seeds, CPU groups, and pinned CPUs; the three GF16 parent expectations;
the full independent 4x4 reuse/batch grid; fixed-work thread scaling; uncapped
large memory estimates and unavailable preflight; executable content identity,
source-metadata preservation, current cgroup v1/v2 limiting-ancestor discovery,
group-atomic partial/mixed-epoch resume with ungrouped job-granular fallback,
and forced-path and dispatcher-check calculations over the complete five-row
AB/BA fixture.
Negative cases exercise bad benchmark schemas,
stale manifest, job, result, and output identities, delayed executable
replacement, request-parameter mismatches, failed round trips, duplicate pair
members, scheduled-pair identity drift, incomplete expected cells, mixed run
epochs, and zero and non-finite timing. It runs through CTest as
`leopard2_benchmark_matrix_self_test`.
The lab self-test additionally uses a content-addressed fake `perf` provider to
verify parsed and hashed counter evidence, atomic resume, post-run corruption
rejection, optional-denied bare execution, required-denied preflight, and
counter-executable replacement rejection. It also rejects positional relabeling
of a different reported event and rejects available/partial evidence without an
available probe, a complete measurement list, and retained hashed raw output.
The exact preflight command is bound into each manifest request. Resume and
matrix collection reparse `perf-stat.txt` through one shared canonical parser
and require the recorded measurements, status, and detail to equal that
raw-derived result; coordinated JSON/digest edits cannot substitute a different
value while leaving the retained raw bytes unchanged. The matrix self-test
applies the same evidence invariants during collection.
Lab manifests therefore use `leopard2-lab-manifest/v3`; older manifests must be
regenerated rather than being resumed under weaker evidence semantics.
CTest also runs the independent operation-count model's schedule invariants;
those counts remain modeled bounds rather than PMU observations.

The benchmark executable separately records codec and decode-plan setup,
execution, setup amortized at the selected reuse count, input and generated or
repaired output throughput, median/MAD/minimum/maximum timing, selected
profile/field/backend, scratch, legacy availability, and round-trip status.
Logical operation counts are provided independently by
`tools/leopard2_operation_counts.py`. Those estimates remain labeled as a
model; only successful, signed `perf stat` observations are presented as
hardware counters.

## External-library comparison policy

ISA-L, Jerasure, FastECC, and other libraries may be added as separate job
executables, but a result is comparable only when the adapter records the same
public `(K,R)`, shard bytes, loss pattern, batch/reuse semantics, thread count,
and amount of generated or repaired output. Field or code-family limitations,
setup included by one side but not the other, and unsupported `R>K` or parent
sizes must be explicit exclusions rather than silently adjusted cells. The
default-off, single-thread ISA-L checkpoint adapter is documented in
`docs/leopard2_isal_comparison.md`. It does not create a production dependency
and does not complete the external-library matrix. Jerasure and the remaining
ISA-L multicore cells still have no reviewed adapter, and a bounded result
does not become a cross-library claim for unmeasured cells. That remaining work
keeps the benchmark-harness Bead open.

The checked-in ISA-L JSON artifacts are the coordinated V2 bounded checkpoint.
They contain raw timing samples from both providers, immutable/direct scalar
correctness oracles, detached recorded-commit builds with before/after
Git-object replay, the actual file-input and transitive runtime-library
closures, a controlled child environment, symmetric Leopard `--skip-legacy`
work, a separately supplied full-correctness artifact binding, advisory CPU
leases, and a post-timing integrity recheck. Strict trusted-cache validation
and the independent external audit pass. The six single-thread cells are
accepted bounded evidence; they do not complete the required matrix or support
claims about multicore execution, other machines, or wire compatibility.

`tools/leopard2_external_comparison.py` makes that boundary machine-readable.
It deterministically regenerates the selected matrix, classifies every job, and
reports aggregate `adapter-available-unmeasured`/`adapter-required`/`excluded`,
reason, and qualification counts. The internal per-cell classifier records the resolved Leopard2
profile/field/parent, the provider field, whether wire compatibility exists
(currently false for every external provider), and the reasons or
qualifications. Its `isa-l-checkpoint` command requires the bootstrap cache,
reconstructs and rehashes the complete trusted build provenance, and invokes
the fail-closed bounded-artifact validator with the exact 128-case correctness
artifact. It will not label a merely portable/self-consistent tool, link, or
runtime relabel as verified. The matrix audit itself still contains no timing
measurements. On the 7,134-job required matrix, the audit classifies
5,003 ISA-L rows as supported by the default-off adapter but explicitly
unmeasured, retains 21 multicore adapter gaps,
and excludes 2,110 whose `K+R` exceeds GF(256); all 7,134 required GF8/GF16
rows are Jerasure adapter candidates because their shard sizes are multiples of the deterministic
eight-byte adapter region contract. ISA-L remains eligible for `K+R<=256`
when Leopard2's dyadic parent inflation selects GF16, but the report requires that field
advantage to be disclosed. Jerasure likewise retains GF8 for those public
lengths rather than silently following Leopard2 into GF16. The ISA-L adapter
uses the MDS-safe `gf_gen_cauchy1_matrix`; ISA-L documents that
`gf_gen_rs_matrix` does not guarantee every submatrix is invertible for many
larger `(K,R)` pairs. It is standalone and default-off, so ordinary Leopard
builds have no ISA-L dependency. FastECC is excluded from all current rows
because its prime-field NTT stores wider parity sectors than source sectors, and
ECC-Benchmark is excluded as an independent provider because it is a historical
harness with different timing/loss semantics. Adapter availability is
explicitly not a measurement.

Jerasure's public region API requires longword-multiple sizes and longword-
aligned pointers. To make offline classification independent of audit-host
`sizeof(long)`, the prospective comparison adapter uses a conservative
eight-byte region granularity and eight-byte alignment contract. Arbitrary
unaligned shard lengths are excluded from direct-region comparison; any future
staging or padding experiment must report copied and padded bytes explicitly.

Reproduce the policy and optional source-cache checks with:

    python3 tools/leopard2_external_comparison.py self-test
    python3 tools/leopard2_external_comparison.py matrix --preset required
    python3 tools/leopard2_external_comparison.py cache \
        --path .research/leopard2

The cache check performed on 2026-07-16 verified all four recorded commits.
This host lacks system NASM, so the ISA-L checkpoint uses the pinned official
NASM 2.16.03 archive and locally built binary described in the dedicated
comparison document. The host still lacks the GF-Complete package required by
Jerasure 2.0. These are recorded host/toolchain facts, not favorable fairness
exclusions.
