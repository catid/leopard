# Nested-thread-safe sanitizer fuzz campaigns

`tools/leopard2_fuzz_campaign.py` creates a content-addressed manifest for the
deterministic API and pruned-transform replay fuzzers. It delegates execution
to `tools/leopard2_lab.py`, so every seed retains its stdout, stderr, timeout,
CPU set, memory policy, executable hash, and terminal JSON independently. A
rerun resumes only valid completed results, and merge order is deterministic.

Campaign creation fails closed unless each target identifies itself as the
expected deterministic replay role and reports both compile-time and linked
runtime ASan/UBSan support. A function compiled in the linked Leopard2 core
translation unit must independently report the same compile-time sanitizer
coverage, preventing an instrumented replay wrapper from laundering an
uninstrumented codec library. A separate `nm` probe requires the corresponding
ASan and UBSan global runtime symbols. The executable hash, role attestation,
core attestation, symbol-table digest, and probing-tool identity are signed into the manifest;
audit repeats the live probes and rejects changed binaries or tools. Merely
setting sanitizer environment variables is not instrumentation evidence.
The v5 producer opens each target and `nm` exactly once, executes it through a
retained `/proc/self/fd` name, and revalidates the same descriptor after all
descendants are contained. Probe stdout and stderr are retained regular files,
not pipes: the manifest signs their per-role limits, hard and soft
`RLIMIT_FSIZE`, default `SIGXFSZ`, retained-name checks, and bounded descriptor
reads. A subreaper plus an unpredictable containment token binds detached
descendants to the probe before capture evidence is read.
The destination is cleared before those checks, so a failed create cannot
leave a stale manifest. An output path that resolves to either fuzz executable
is rejected without modifying that executable.

The parallel campaign covers AddressSanitizer and UndefinedBehaviorSanitizer,
but deliberately sets `detect_leaks=0`. On Linux, LeakSanitizer's
stop-the-world teardown can briefly clone a same-executable helper process. It
is not a codec worker, but it is indistinguishable from an application fork to
the generic `/proc` session monitor. Letting that helper through would weaken
the campaign's one-thread evidence: a real extra application process could use
the same allowance. LeakSanitizer therefore runs as the separate, serial,
full-seed companion campaign documented below. Both audited phases are
required for the complete sanitizer gate.

The default campaign creates one API and one pruned-transform job per allowed
CPU, up to 128 CPUs. Every job receives one CPU and declares one aggregate
workload thread. The lab runner replaces inherited nested-runtime settings with
a signed environment including:

- `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`, `OMP_DYNAMIC=FALSE`,
  `OMP_NESTED=FALSE`, and `OMP_MAX_ACTIVE_LEVELS=1`;
- one-thread caps for OpenBLAS, GotoBLAS, MKL, BLIS, Accelerate/vecLib,
  NumExpr, and Rayon; and
- `MKL_DYNAMIC=FALSE`.

The signed sanitizer environment is
`ASAN_OPTIONS=abort_on_error=1:detect_leaks=0:halt_on_error=1:suppressions=`,
`LSAN_OPTIONS=detect_leaks=0:suppressions=`, and
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=`. In
compiler-rt, explicit sanitizer environment settings take precedence over
linked runtime defaults, so both the failure policy and empty suppression-file
paths are bound. Inherited settings cannot silently re-enable leak checking or
add a suppression file inside the audited parallel phase.

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

The auditor also reconstructs the complete versioned campaign contract before
reading results. It binds the exact API and pruned job IDs, stable seeds,
content-addressed executable assigned to each role, `{seed}` and iteration
arguments, sanitizer environment, timeout, memory policy, one-thread runtime
policy, and deterministic CPU assignment. A generic lab manifest cannot become
sanitizer evidence merely by copying campaign-looking metadata. Historical
artifacts used `leopard2-fuzz-campaign/v4` and
`leopard2-fuzz-campaign-audit/v4`. New artifacts use
`leopard2-fuzz-campaign/v5` and
`leopard2-fuzz-campaign-audit/v5`. The v5 manifest signs the exact hardened
probe-execution policy; the audit repeats that policy, records whether it was
bound by the source manifest or supplied by a secure live replay of a
historical v4 source, and signs the terminal artifact with `audit_digest`.
An old v4 manifest is accepted only by the audit path that immediately repeats
the current live probes. Offline consumers and new leak companions reject it,
and new code never emits a v4 audit. The metadata and audit output explicitly
record that ASan and UBSan are active,
LSan is inactive, and leak checking is a separate phase.

The strict campaign is currently Linux/ELF evidence: independent
instrumentation checks use `nm`, and live process, affinity, thread, and RSS
validation requires `/proc`. Other hosts fail closed instead of emitting a
weaker artifact.

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

The final output uses `leopard2-fuzz-campaign-audit/v5`, embeds the ordinary
deterministic lab merge, the exact probe policy and its binding origin, and
contains an `audit_digest` over the complete terminal object. It is written
atomically only after every stricter
campaign gate passes. The destination is cleared before manifest validation,
so a failed audit cannot leave a stale artifact that could be mistaken for
accepted fuzz evidence. The destination may not overwrite the manifest or any
per-job result. Create and audit also refuse to remove an existing executable;
audit explicitly protects both fuzz targets and the recorded instrumentation
tool, including resolved symlink aliases.

The 60 stable seeds are distinct. CPU assignment may repeat between jobs, but
the scheduler never runs jobs with overlapping CPU sets simultaneously. Even
the generic `--allow-cpu-overlap` mode will serialize work when aggregate
declared thread demand would exceed the union of allocated CPUs.

## Full-coverage leak-check companion

Derive the LeakSanitizer manifest from the already validated v5 parallel
manifest. The companion reproduces every source job ID, stable seed, iteration
count, executable identity, command, timeout, working directory, and assigned
CPU. For the 30-by-2 example above this is the same 60 jobs and 491,520 replay
iterations, not a two-test smoke sample:

    python3 tools/leopard2_fuzz_campaign.py leak-create \
      --manifest build/fuzz-campaign/manifest.json \
      --output build/fuzz-campaign/leak-manifest.json

    python3 tools/leopard2_fuzz_campaign.py leak-run \
      --manifest build/fuzz-campaign/leak-manifest.json \
      --output-dir build/fuzz-campaign/leak-results

    python3 tools/leopard2_fuzz_campaign.py leak-audit \
      --manifest build/fuzz-campaign/leak-manifest.json \
      --output-dir build/fuzz-campaign/leak-results \
      --output build/fuzz-campaign/leak-audited-results.json

### Integrated evidence checkpoint (2026-07-18)

The tree at commit `72597bf` was rebuilt with Clang 18.1.3 and
ASan+UBSan, then exercised on the 30-CPU allowed set. The parallel v4 phase
passed and audited 60 of 60 distinct jobs (30 API and 30 pruned-transform),
8,192 iterations per job and 491,520 iterations total. A second run resumed
all 60 results without execution. Its manifest digest is
`33dc9e4450f57dca3c3c18f8beed31cfc02d5cd1c5ec17a1df285edfd7e592da`;
the manifest and audit file SHA-256 values are respectively
`582af8d51b28e95f5c8555a3279617a175dd4a52d1a88284c50ade25f44308d1`
and
`343a034fad04a4df87906961b824cb8797c58f41f127b456a8a27ea144070a43`.

The then-current serial LSan v3 companion at that historical checkpoint passed
and audited the same 60 jobs and 491,520 iterations, followed by a 60-of-60
zero-execution resume. It signed
both LSan control hooks as undefined. Its manifest digest is
`2c7a5a11cb1e13eebe7e4184b798a5ff7d4105801469ad18e07757126fb8ded0`;
the leak manifest and audit file SHA-256 values are respectively
`59e9688312a4b1cdf85add16d407d3de895382469fc83e76990a53fd91b9a390`
and
`47e6e82cb76feb50507417657e391505777ddb08472a56fecc5eeaccb6a208c5`.
The ignored local artifacts are under
`build/fuzz-campaign-integrated-72597bf/`; regenerate them after any codec or
instrumentation change rather than treating these hashes as current forever.

`leak-create` repeats the live role/core ASan+UBSan attestation and independent
ELF symbol probe before signing the companion. It then invokes a separate
linked-core canary that deliberately leaves one unreachable 12,345-byte
allocation and requires compiler-rt's own LeakSanitizer diagnostic, the linked
core canary frame, and exact exit status 86 from a recoverable check.
`leak-run` and `leak-audit` repeat both the canary and the symbol check before
accepting retained work, so `detect_leaks=1` is evidence-backed rather than
inferred from an environment string. Defined
`__lsan_is_turned_off` and `__lsan_default_suppressions` control hooks are
rejected; undefined weak runtime references are accepted and signed.
`leak-run` is deliberately serial and sets
`ASAN_OPTIONS=abort_on_error=1:detect_leaks=1:halt_on_error=1:suppressions=`,
`LSAN_OPTIONS=detect_leaks=1:suppressions=`, and
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=`, together with the same
one-thread native/OpenMP environment used by the parallel campaign. Ambient
sanitizer or thread settings cannot weaken those values.

The intentional canary uses its own separately signed probe environment:
`ASAN_OPTIONS=abort_on_error=0:detect_leaks=1:halt_on_error=1:suppressions=` and
`LSAN_OPTIONS=detect_leaks=1:exitcode=86:suppressions=`, together with
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=`. Disabling
abort is deliberate: an expected leak must not create a core file or invoke a
host crash reporter. The replay flushes the diagnostic and exits directly with
86 after the recoverable check, avoiding a second teardown report. The manifest stores the canonical
exit, stdout hash/size, allocation size, required diagnostic markers, and a
normalized diagnostic-evidence hash. Actual replay jobs retain
`abort_on_error=1` and fail on any sanitizer diagnostic.

The companion intentionally permits transient LeakSanitizer helpers and does
not make a process-count claim. It still pins the replay to its signed CPU and
requires native/OpenMP thread limits of one. The v4 parallel audit remains the
historical checkpoint for one-process runtime behavior; current campaigns use
the v5 parallel audit. Never substitute the
companion for that audit or remove `detect_leaks=0` from the parallel manifest.

Each companion job retains content-addressed stdout, stderr, executable
identity before and after execution, and a terminal result. A rerun resumes
only an exact successful result. Any stderr (including a LeakSanitizer
diagnostic), nonzero exit, timeout, changed executable, malformed success
marker, stale result, symlink/hard-link alias, missing job, or extra job fails
closed. Audit clears a stale destination before validation, checks the complete
matrix, merges in signed manifest order, and emits
`leopard2-fuzz-leak-campaign-audit/v4`. The manifest, per-job result, and merge
schemas are respectively `leopard2-fuzz-leak-campaign/v4`,
`leopard2-fuzz-leak-result/v2`, and `leopard2-fuzz-leak-merge/v2`.
Earlier companions must be regenerated. V1 did not independently prove that
the target had not disabled LSan; v2 did not attest the selective
default-suppression hook; v3 did not bind or enforce finite stdout and stderr
captures. The v4 contract retains the v3 statically linked suppression-hook
and inherited-suppression-file checks and additionally signs a 2 MiB limit for
each captured stream. The child receives the corresponding `RLIMIT_FSIZE` and
default `SIGXFSZ` disposition. Cleanup validates the retained file descriptors
before reading at most the signed limit plus one byte; a signal, `EFBIG`,
exact-limit saturation, oversized retained inode, or concurrent change makes
the result invalid. Detached writers are contained before any capture is read,
so neither timeout nor resume can wait for pipe EOF or allocate from unbounded
output. A privileged loader configuration or a runtime-injected library
remains outside this artifact's trust boundary and must be controlled by the
test host.

For a quick local leak smoke only, the two deterministic CTests may still be
run serially with the explicit ASan/LSan/UBSan environment. They are useful
during development but are not the full-coverage leak artifact and do not
replace the three companion commands above:

    ASAN_OPTIONS=abort_on_error=1:detect_leaks=1:halt_on_error=1:suppressions= \
    LSAN_OPTIONS=detect_leaks=1:suppressions= \
    UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions= \
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE \
    ctest --test-dir build/asan-ubsan -j1 --output-on-failure \
      -R '^leopard2_(fuzz_smoke|pruned_fuzz_smoke)$'

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
