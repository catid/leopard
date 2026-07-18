# Nested-thread-safe sanitizer fuzz campaigns

`tools/leopard2_fuzz_campaign.py` creates a content-addressed manifest for the
deterministic API and pruned-transform replay fuzzers. It delegates execution
to `tools/leopard2_lab.py`, so every seed retains its stdout, stderr, timeout,
CPU set, memory policy, executable hash, and terminal JSON independently. A
rerun resumes only valid completed results, and merge order is deterministic.

The default campaign creates one API and one pruned-transform job per allowed
CPU, up to 128 CPUs. Every job receives one CPU and declares one aggregate
workload thread. The lab runner replaces inherited nested-runtime settings with
a signed environment including:

- `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`, `OMP_DYNAMIC=FALSE`,
  `OMP_NESTED=FALSE`, and `OMP_MAX_ACTIVE_LEVELS=1`;
- one-thread caps for OpenBLAS, GotoBLAS, MKL, BLIS, Accelerate/vecLib,
  NumExpr, and Rayon; and
- `MKL_DYNAMIC=FALSE`.

The effective values are recorded in both the job manifest and terminal result.
On Linux, the runner samples the complete process session through `/proc`,
records peak process count, aggregate thread count, RSS, and process affinity,
and rejects a successful result if the workload exceeds its declared thread or
CPU allocation. A sampled RSS cap is available separately from `RLIMIT_AS`;
this matters for ASan because its shadow memory reserves a very large virtual
address range.

The generic runner retains a launch-floor observation after its affinity hook
successfully reaches `exec`, which keeps sub-millisecond utility jobs resumable.
The fuzz campaign auditor is stricter: it requires at least one later `/proc`
sample with observed affinity. RSS is sampled at the same time (and may read as
zero if a very short process has already become a zombie). A campaign cannot
pass merely on the launch floor.

## Reproduce a 30-by-2 campaign

First build `leopard2_fuzz_smoke` and `leopard2_pruned_fuzz_smoke` with Clang
ASan and UBSan applied to both the codec and executables. Then run:

    python3 tools/leopard2_fuzz_campaign.py create \
      --api-executable build/asan-ubsan/leopard2_fuzz_smoke \
      --pruned-executable build/asan-ubsan/leopard2_pruned_fuzz_smoke \
      --seeds-per-target 30 \
      --iterations 8192 \
      --minimum-memory-mb 512 \
      --rss-limit-mb 2048 \
      --output build/fuzz-campaign/manifest.json

    python3 tools/leopard2_lab.py run \
      --manifest build/fuzz-campaign/manifest.json \
      --output-dir build/fuzz-campaign/results

    python3 tools/leopard2_fuzz_campaign.py audit \
      --manifest build/fuzz-campaign/manifest.json \
      --output-dir build/fuzz-campaign/results \
      --output build/fuzz-campaign/audited-results.json

The 60 stable seeds are distinct. CPU assignment may repeat between jobs, but
the scheduler never runs jobs with overlapping CPU sets simultaneously. Even
the generic `--allow-cpu-overlap` mode will serialize work when aggregate
declared thread demand would exceed the union of allocated CPUs.

## Intentional internal teams

Generic lab specifications remain able to exercise a deliberate internal
thread team. The opt-in is explicit and must fit the job CPU set:

    {
      "cpu_count": 4,
      "thread_runtime": {
        "max_threads": 4,
        "allow_internal_team": true
      }
    }

`max_threads` counts the main workload thread, helper threads, and descendant
process main threads in aggregate. Requesting more than one thread without the
opt-in, requesting more threads than CPUs, or placing a conflicting controlled
value in `env` is a manifest error. Runtime excess makes the terminal outcome
`evidence_invalid`; RSS excess makes it `memory_limit`. Neither can be audited
as successful campaign evidence.

Lab manifest, result, and merge schemas are respectively v4, v2, and v2.
Earlier artifacts must be regenerated because they do not bind the effective
nested-runtime environment or runtime observations.
