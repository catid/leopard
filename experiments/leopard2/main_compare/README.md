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

For ISA-matched AVX2 attribution on hosts that also expose AVX-512, configure
the separate, default-off exact-main profile with
`-DLEO_MAIN_PURE_AVX2=ON`.  It keeps the same exact main sources and API while
using the explicit compiler ceiling
`-march=x86-64 -mtune=generic -mavx2 -mno-avx512f`.  This option is for
controlled comparison builds; the default continues to reproduce main's
historical `-march=native` policy. Current version-16 ABBA evidence requires
this pure-AVX2 profile and the matching `--baseline-pure-avx2` runner selector;
the native build above is useful for historical reproduction but is not a
version-16 baseline. A suitable baseline build is:

    cmake -S experiments/leopard2/main_compare \
        -B /tmp/leopard-main-pure-avx2 -G "Unix Makefiles" \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEOPARD_MAIN_SOURCE_DIR=/tmp/leopard-main-exact \
        -DLEO_MAIN_PURE_AVX2=ON
    cmake --build /tmp/leopard-main-pure-avx2 -j "$(nproc)"
    ctest --test-dir /tmp/leopard-main-pure-avx2 --output-on-failure

The broad, CPU-saturating `run_allk_gap.py` diagnostic is deliberately separate
from the isolated promotion runner below. It requires the Leopard2 source root
and requested commit explicitly. The source must be the clean Git top level at
both ends of the run with no tracked modifications, every Leopard2 result must
embed that exact commit and
tree, every exact-main result must embed the requested baseline commit, and both
executable byte/metadata identities must remain unchanged. Every child executes
an immutable sealed snapshot captured into the run contract, so replacing an
executable pathname during a campaign cannot mix binaries. Existing output can
be resumed only when its complete run contract matches; stale cell files are
rejected rather than mixed into a new result.
The runner also rejects Git replace refs and stale or mixed production objects:
it binds the exact Release compile/source/object/archive/link closure, extracts
and compares every archive member with its external object, and performs a
second runner-owned configure/build in an empty directory.  All linked object,
archive, and executable bytes must reproduce before timing begins.  This
fail-closed parser currently requires the Unix Makefiles generator.

Configure this default-off diagnostic explicitly and build its attested target:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi

    cmake -S . -B /tmp/leopard2-all-k -G "Unix Makefiles" \
        -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DLEO2_BUILD_BENCHMARKS=ON -DLEO2_BUILD_ALLK_DIAGNOSTIC=ON \
        -DLEO2_BACKEND_VARIANT=avx2 -DLEO2_ENABLE_CUDA=OFF \
        -DLEO2_FLAG_MAVX512F=FALSE -DLEO2_FLAG_MAVX512BW=FALSE \
        -DLEO2_FLAG_MAVX512VL=FALSE
    cmake --build /tmp/leopard2-all-k -j "$JOBS" \
        --target bench_leopard2_allk

    python3 experiments/leopard2/main_compare/run_allk_gap.py \
        --main /tmp/leopard-gf8-final-map-9b1d439-20260722T030328Z/main-build/leopard_main_benchmark \
        --main-commit 6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198 \
        --main-sha256 a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910 \
        --current /tmp/leopard2-all-k/bench_leopard2_allk \
        --current-source-root "$PWD" \
        --current-commit "$(git rev-parse HEAD)" \
        --output /tmp/leopard2-all-k-gap --workers "$JOBS" \
        --reproducible-build-jobs 1 \
        --gf8-only --with-current-legacy \
        --lock /tmp/leopard-gf8-authoritative.lock

The final GF8 run contains 2,522 cells spanning every `K=1..255`, three
redundancy bands where distinct, 4 KiB and 64 KiB shards, and one/maximum
losses.  The runner requires the frozen pure-AVX2 Leopard-main executable SHA,
rejects EVEX in either executable, forces the Leopard2 GF8/AVX2 context, and
holds the canonical benchmark lock for the complete campaign.  It intentionally
saturates the allowed CPUs and identifies regions for follow-up; its output is
not authoritative single-core performance evidence.  Every near or losing
region still requires the isolated ABBA runner before a promotion or final
disposition.  One untimed schema-v5 preflight binds the sealed Leopard2
executable to the exact clean source commit/tree.  Timed calls then use the
distinct schema-v3 decode-path report, so summaries record the implementation
actually selected instead of inferring a potentially stale route from K/R.

The current all-K run-contract and manifest generation is v19. It requires the
production-default exact T32/B256 terminal and T16/Q2 fused selectors to be
`ON`, the legacy-high sparse direct encoder to be `OFF`, and legacy-high direct
AUTO routing to be `ON`; sparse direct AUTO remains `OFF`. Historical v18 stays
pinned to benchmark configuration v15 and has no sparse-AUTO field. Historical
v16 remains pinned to configuration v13, while v17 binds configuration v14 and
the default-off sparse selector; relabeling any body as another generation is
rejected.
The 2,522-cell matrix uses 4 KiB and 64 KiB shards, so the exact 256-byte
T32 terminal and at-most-64-byte T16/Q2 kernel cannot affect a timed cell, but
their values still belong to the authenticated production selector identity.

The final-source GF8 checkpoint is
`results/b1f334a-final-allk-gf8/summary.json`.  Its broad 2,522-cell run is a
discovery sweep; every one of its 63 apparent losing cells was subsequently
replayed on one isolated physical core.  No encode or reused-decode regression
survived outside the degenerate K=1,R=1 copy case.  The K=2,3,4,5,7 R=1 cells
that lose only when reusable-plan construction is charged to a single call all
beat exact main through the comparable allocation-free one-shot API.  K=1,R=1
retains Leopard2's stronger overlap, range, scratch, and failure-atomicity
validation around an otherwise identical copy; the checkpoint records both
that residual and the previously exhausted safe specialization candidates.

## Counterbalanced comparison

`run_abba.py` supplies the other half of the comparison. Its current evidence
contract is raw/manifest/failure version 18; version 17 remains the frozen
native-main AUTO GF16-GFNI predecessor. It accepts only fresh Release
artifacts with the expected compile, object, archive, link, source,
runtime-library, field/profile, and requested decoder-path identities. Each cell
runs three independent
`baseline,candidate,candidate,baseline` rounds on one pinned logical CPU while
the sibling is reserved. Current evidence additionally holds a pair-wide,
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

Versions 12 through 16 are frozen effective-AVX2 comparisons for source
revisions where GFNI remained explicit-only, not merely requests for the AUTO
backend. They must not be used with a source revision whose operation-specific
AUTO GF16/GFNI selector is production-enabled. Leopard main must use the
pure-AVX2 build above. Leopard2 must use
`LEO2_BACKEND_VARIANT=auto`, and the benchmark must resolve that AUTO request to
`avx2`. The AVX-512F/BW/VL and preferred-vector-width probes are forced false,
so `Leopard2BackendAVX512.cpp` cannot enter the build closure. The
VEX-encoded GFNI translation unit may remain in the archive with
`-mavx2 -mgfni -mno-avx512f`. On family 1Ah/model 08h, current production
source can borrow that table for the exact K=1000/R=200/64-KiB full-output
cell while the context still reports `avx2`; schema v16 cannot observe that
operation-specific substitution. A successor schema must either pin the
selector off and attest the zero-route control, or attest the production
default GFNI route under a separately named comparison. Until then, a current
AUTO-GFNI candidate cannot be labeled version-16 AVX2 evidence.

### Version 17 production GFNI and passive shared-host observations

Version 17 is the exact production-route comparison for the single GF16
legacy-high cell K=1000, R=200, B=65536, full parity, one thread. Its Leopard1
arm is detached commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`
built under Leopard1's native `-march=native` Release policy. Its Leopard2 arm
requests AUTO: the context reports AVX2, while every candidate child attests
that the exact encode operation used the startup-qualified VEX-256 GFNI table.
The two reported ratios use the same Leopard1 `leo_encode` observations as the
numerator and separately measure Leopard2 ordinary one-item
`leo2_encode_batch` and one-shot `leo2_encode`. They are correlated estimates
and must not be multiplied or stacked.

The original wrapper uses affinity supervisor v14 to confine every same-UID
process away from CPU52 and its SMT sibling CPU116. On the current OBS
streaming host that supervisor cannot operate non-disruptively: its v4 attempt
rejected a Chromium GPU process with nonuniform per-thread masks before the
runner launched, then its fail-closed restoration path terminated that
unrelated process. The sealed v4 envelope contains no campaign directory and
no timing samples. Retrying the active supervisor while OBS is running is
therefore neither safe nor useful.

The same wrapper historically provided the following weaker v17 acquisition
generation. This command is retained here only to identify the sealed contract;
fresh use is now rejected because the live producer has advanced to v18:

    experiments/leopard2/main_compare/run_authoritative_v17_gfni_main_compare.sh \
        --passive-shared-host \
        /home/catid/leopard/.research/leopard-79h/<commit>-v17-passive-main-v1

This mode changes affinity only for the wrapper and its owned descendants. It
pins the wrapper to housekeeping CPUs, starts the owned runner with the
original allowed set, and lets the unchanged v17 runner pin every benchmark
child to CPU52. It never executes the active supervisor, never changes or
signals an unrelated process, requires JSON-null campaign supervision, and
retains the runner's pair lease, reservation, source/build/digest, route,
20-ms timer-window, three-round ABBA, and zero-nonidle-jiffy CPU116 gates.
Read-only pre/post endpoint censuses disclose same-UID affinity eligibility and
an enclosing `/proc/stat` interval. They are non-atomic and explicitly do not
claim complete process history.

Passive evidence uses the unchanged
`leopard2-main-compare-{raw,manifest}/v17` schemas because null supervision is
an existing, tested v17 state. Its surrounding assurance is separately
versioned as
`leopard2-main-compare-v17-passive-independent-audit/v1`,
`leopard2-v17-gfni-main-passive-core-manifest/v4`, and
`leopard2-v17-gfni-main-passive-not-promoted-envelope/v2`. A successful
passive acquisition always publishes `NOT_PROMOTED.json`, even when both
observed lower confidence bounds exceed 1.05. The audit, result summary, core,
and terminal envelope all bind `promotion_eligible:false`,
`cpu_pair_exclusive:false`, and
`causal_performance_claim_eligible:false`.

The defensible interpretation is limited: CPU116 accumulated zero observed
nonidle `/proc/stat` jiffies, and CPU52's descriptive nonidle-jiffy excess over
the retained child-wall-time ceiling stayed within its rejection threshold in
both the runner interval and the enclosing census interval. This is a
one-sided rejection screen, not process attribution or an interference upper
bound, and neither condition proves CPU52 exclusivity. Foreign work can fit
inside child wall time or below the 100-Hz jiffy resolution, and active work
elsewhere can change package boost/thermal residency, shared LLC state, and
memory bandwidth. The policy discloses CPU52 nice, iowait, IRQ, softirq, and
steal deltas: nice/IRQ/softirq/steal are already included in aggregate
nonidle, while iowait is diagnostic idle accounting; none has a separate
zero-tolerance gate. The
three-round Student-t interval is consequently reported as a nominal clustered
interval under shared-host load, not a quiet-host or causal performance
guarantee. Any number from this path must be described as a host-, compiler-,
API-, and workload-specific passive observation; promotion still requires a
dedicated or safely isolated environment.

The retained passive-v17 attempt at
`.research/leopard-79h/2cc900f-v17-passive-main-v1` is a sealed, verified
failure and supplies no ratio or confidence interval. All twelve children
completed, but CPU116 accumulated 388 nonidle jiffies during the 2293-second
runner interval. The retained child wall durations total only 50.755 seconds;
even if the sibling gate had passed, CPU52 accumulated 155 more nonidle jiffies
than the old aggregate child-wall ceiling, above the v17 limit of 16. Those
samples are diagnostic only and must never be pooled, re-summarized, or reused
by a successor campaign.

The earlier active-v17 diagnostic envelopes are also retained as rejections:
`.research/leopard-79h/b9167ad-v17-gfni-main-v1`,
`.research/leopard-79h/4db018f-v17-gfni-main-v2`,
`.research/leopard-79h/76b0cb4-v17-gfni-main-v3`, and
`.research/leopard-79h/0319e39-v17-gfni-main-v4`. Together with the passive-v17
envelope above, they must be cited but never imported into a v18 lane.

### Version 18 prospective per-invocation passive screen

Version 18 replaces the aggregate passive-v17 contamination screen with a
prospectively fixed launch-through-reap screen around each of the twelve ABBA
children. Each endpoint is one fresh read-to-EOF of `/proc/stat` containing
both CPU52 and CPU116, with its read-start and read-finish monotonic boundaries
retained. For invocation `i`, with child wall duration `d_i` nanoseconds, the
contract is exact integer arithmetic:

    cap_i = (d_i * 100) // 1000000000 + 1
    excess_i = max(0, cpu52_nonidle_i - cap_i)
    accept_i = (excess_i == 0 and cpu116_nonidle_i == 0)

There is no tuned allowance in those expressions. For wall duration `d`, at
most `floor(d * f) + 1` points of a frequency-`f` accounting lattice can fall
inside the interval, independent of its phase. At Linux `SC_CLK_TCK=100`, that
forces the CPU52 cap above; precise virtual CPU accounting is also bounded by
wall time plus the same one-jiffy conversion allowance. The benchmark owns no
work on CPU116, and one jiffy is the instrument's smallest positive value, so
zero measurable sibling activity is the only non-arbitrary finite screen.
Allowing one or more sibling jiffies would introduce a free policy parameter.

Against the retained v17 failure, this is a tighter observed domain, not a
post-hoc relaxation:

| | gated window | aggregate allowance |
|---|---:|---:|
| v17 | 2,293.2 s | `ceil(sum(d_i) * 100 / 1e9) + 16` = 5,076 + 16 = **5,092** |
| v18 | 50.755 s (**45.2 times smaller**) | `sum(floor(d_i * 100 / 1e9) + 1)` = **5,081** |

Version 18 removes 11 jiffies of aggregate slack and 2,242 seconds of
unmeasured exposure. The removed interval was not retained measurement time:
each child path recomputed full Git identities for both source roots and build
provenance for both build trees before and after the short benchmark, producing
about 186 seconds of input-snapshot/hash dead time between child measurements.
That diagnosed interval error, rather than either observed failure count, is why
v18 changes the measured domain.

A valid campaign requires all twelve `accept_i` values to be true. CPU52 and
CPU116 whole-run deltas and their residuals outside the retained child windows
are nonnegative structural/disclosure records, not contamination gates. The
pre/post same-UID census remains a fail-closed endpoint safety check: no thread
may be confined to the reserved pair and a persistent retained thread may not
change affinity between endpoints. It is still non-atomic, does not record
execution history, and does not establish exclusivity. The result remains a
host-, compiler-, API-, and workload-specific passive nominal observation;
`promotion_eligible`, `promotion_passed`, `cpu_pair_exclusive`, and
`causal_performance_claim_eligible` are always false.

The campaign procedure is preregistered in the same commit as the implementation:

1. **Thresholds are frozen at commit time.**
   `cap_i = floor(d_i * 100 / 1e9) + 1`, per-window CPU52 excess `== 0`, and
   per-window CPU116 nonidle `== 0`. They may not be changed after any attempt.
   If they are ever changed, every prior attempt is void and the campaign
   restarts from a new preregistration commit.
2. **Attempt budget `N_max = 3`.** The wrapper enforces
   `--attempt n --attempt-budget 3`; the attempt and budget are recorded in the
   core manifest and both terminal envelope kinds. Lanes are named
   `.research/leopard-79h/<commit7>-v18-passive-main-a<N>`.
3. **All attempts are retained and cited.** Any report of the observation must
   state “attempt N of at most 3” and list every sealed envelope, including
   `FAILED` ones. This keeps fail-fast from degenerating into retry-until-pass.
4. **No host mutation.** No affinity change to any foreign process, no signals,
   no renice, and no cgroup edits. The wrapper mutates only its own affinity and
   its owned descendants (`affinity_mutation_scope` is unchanged).
5. **If all three attempts reject, stop.** Do not retune. Either escalate to a
   newly preregistered pre-flight pair-qualification screen, committed before it
   is run, or to a dedicated/quiet host. Both require a fresh preregistration
   commit and, where applicable, external authority.
6. **Reporting language is fixed in advance.** This is a host-, compiler-, API-,
   and workload-specific *passive nominal observation*; its clustered
   three-round Student-t interval is reported as nominal under shared-host load;
   the two ratios are correlated and must not be multiplied or stacked; there is
   no exclusivity, causal, or promotion claim.
7. The pre-existing v17 sealed envelopes, including
   `.research/leopard-79h/2cc900f-v17-passive-main-v1`, are cited as prior
   rejections. Their timing samples are not reused, re-summarized, or pooled with
   any v18 result.

For a committed source revision, attempt one is launched as:

    experiments/leopard2/main_compare/run_authoritative_v17_gfni_main_compare.sh \
        --passive-shared-host --attempt 1 --attempt-budget 3 \
        /home/catid/leopard/.research/leopard-79h/<commit7>-v18-passive-main-a1

The historical wrapper filename and `v17-gfni-encode` preset string are frozen
workload identifiers; they do not imply a v17 evidence schema. Fresh active-v17
acquisition is deliberately disabled because the live producer now emits the
passive-only v18 contract. Existing active and passive v17 envelopes remain
verifiable by the versioned replay branches.

Version 18 uses `leopard2-main-compare-{raw,manifest,failure}/v18`, isolation
schema `leopard2-main-compare-isolation/v2`, audit schema
`leopard2-main-compare-v18-passive-independent-audit/v1`, census policy schema
`leopard2-passive-shared-host-policy/v2`, result schema
`leopard2-v18-gfni-main-passive-result-summary/v1`, core schema
`leopard2-v18-gfni-main-passive-core-manifest/v1`, success-terminal schema
`leopard2-v18-gfni-main-passive-not-promoted-envelope/v1`, and failure-terminal
schema `leopard2-v18-gfni-main-failed-envelope/v1`. The four scalar counters in
the core and success terminal bind the zero windowed CPU52 excess, zero windowed
CPU116 nonidle count, and the two disclosed out-of-window nonidle residuals.
Every v18 outcome also binds
`leopard2-v18-gfni-main-attempt-lineage/v1`; early wrapper failures additionally
retain `leopard2-v18-gfni-main-passive-wrapper-failure/v1`. The lineage commits
to each earlier attempt's canonical path and outer `SHA256SUMS`, and durable
verification reopens the complete failed-attempt chain.

### Version 19 conditioned-passive pair-qualified screen (dormant, NOT ARMED)

The standalone `v19_host_preflight.py` now implements the read-only pre-lane
host/resource gate. It consumes the exact NOT ARMED preregistration bytes,
requires ripper's full 128-thread / 64-core topology, one NUMA node and eight
L3 domains, and checks the process's actual unified cgroup membership. The leaf
must have a 512 MiB memory limit, zero swap limit and usage, and clean memory
events; visible ancestors cannot tighten that memory cap. Host swap must be
disabled and the soft open-file limit must be at least 65,536.

The collector makes two matching normalized observations, including namespace,
scope-directory and kernel-file identity checks. Its reads are bounded and
reject symlink traversal and observed file/parent replacement. These checks
follow the [kernel's cgroup v2 memory semantics](https://docs.kernel.org/admin-guide/cgroup-v2.html#memory-interface-files).
They are **not** an atomic snapshot, continuous resource authority, or proof of
host exclusivity. Replay checks internal consistency, not authenticity of a
claimed host. The file-identity digest is only a capture fingerprint.

This gate does not create a campaign lane, execute another program, build,
qualify a pair, change affinity, or run a benchmark. The normal and optimized
`leopard2_v19_host_preflight_*self_test` CTests use synthetic observations and
temporary filesystem fixtures only. Integration into the still-dormant host /
build / artifact orchestrator and its sealed controller closure remains
required; later stages must revalidate the relevant live identities and limits.

`v19_owned_artifacts.py` supplies the next dormant ownership primitives. Its
context captures the host before acquiring the canonical lock, verifies the
actual descriptor's exclusive lock record, and refreshes host/resource facts
at boundaries. Lock path mutation history is retained; a harmless competing
open/close does not invalidate it. Unlocks, downgrades, changed identities, and
observed pathname replacement do. This is cooperative serialization, not a
defense against arbitrary same-process unlock/relock between observations.

The artifact owner copies the four fixed build outputs into a fresh
`v19-artifacts` subtree of an existing private lane, requires the pinned
preflight hashes, seals the file/directory modes, and retains guarded files
plus sealed executable memfds. Each boundary rehashes the held files, not
cached copies: [inotify does not report mmap-based writes](https://man7.org/linux/man-pages/man7/inotify.7.html).
It never lends the lock descriptor or executes
anything. Each artifact is bounded to 16 MiB and the total to 32 MiB. Partial
copies remain for diagnosis and cannot be reused by the failed owner. The
normal/optimized `leopard2_v19_owned_artifacts_*self_test` CTests exercise real
temporary files, locks and kernel seals with synthetic host/build inputs.

`v19_retained_preflight.py` now supplies a read-only owner for the historical
build preflight. It verifies the preregistered outer checksum manifest, exact
bounded inventory, source/build identities, successful terminal/resource
records, candidate provenance, and all first/replay artifact hashes. It holds
guarded descriptors and rehashes current bytes at every boundary, including
empty diagnostic files. Limits are 64 manifest-listed files, 16 MiB per file,
32 MiB total, and 1 MiB per JSON document. Files must be unaliased and mode
0444; directories must be owner-only. The normal/optimized
`leopard2_v19_retained_preflight_*self_test` CTests use synthetic filesystem
fixtures, including coherently resealed contradictions, pathname replacement
and writable-mapping mutations. No retained script or executable is run.

The returned source pins and canonical build path describe historical inputs,
not a new build or a live execution authority. The actual preflight remains on
ripper under `.research/leopard-79h/cf7a705-v19-build-preflight-ripper-a3`.
Its canonical build tree must not be overwritten: the baseline output is
path-sensitive. A [bounded relocation experiment](v19_build_relocation.md)
reproduced all four pinned outputs from fresh scratch sources using an explicit
baseline-only compiler path map and task-owned source-cache eviction. That
changed recipe is not yet authorized by the frozen acquisition contract and
does not complete source/runtime ownership. A caller must hold and revalidate
the host/lock owner around
this context and preserve the historical evidence while staging fresh builds.

`v19_streamed_sources.py` adds a staged-file lifetime owner for the fixed
`candidate-source` and `leopard1-source` children of a private workspace. It
retains no-follow directory/file descriptors for their complete inventories,
including ordinary `.git` trees, and checks mutation history, metadata,
pathnames and current content at boundaries. Hash reads are limited to 64 KiB;
file contents are not retained. Limits are 4,096 files, 1,024 directories,
64 MiB per file, 256 MiB total, and 32 path components. Links, shared file
inodes, special files, world-writable entries and filesystem crossings fail
closed. Failed owners cannot be reused even if the changed bytes are restored.

Optional cache eviction flushes and advises only the exact held nlink=1 file
descriptors, after hashing; it never reopens a pathname for eviction or uses
global cache controls. It avoids rehashing file contents after eviction,
which would repopulate the cache before compilation. These are sequential
observations, not an atomic snapshot or a cache-residency guarantee. The
normal/optimized `leopard2_v19_streamed_sources_*self_test` CTests include
mutation, mmap, fd-cleanup, eviction and a 32 MiB bounded-buffer fixture.

This owner proves neither initial Git commit/tree authenticity nor independent
source creation, and its records say so. Those gates must run while the owner
is held. It does not invoke Git or compilers, authorize the mapped recipe,
verify runtime closure, or arm acquisition. Its preflight/host/ownership
dependencies share one module chain to avoid duplicating the large provenance
module inside the narrow build memory budget.

The 2026-09-05 native ripper diagnostic held this owner across fresh candidate
and baseline builds: all four complete output hashes matched the original
preflight, and the 1,341-file source records and owner exit checks passed.
Peak cgroup memory was 433,385,472 bytes within 512 MiB, with every memory-event
counter and swap usage zero. Focused CTest passed 8/8 (including the 16-case
source suite in both Python modes); the project-graph suite passed 172/172
normally and with assertions disabled. The sealed 106-entry diagnostic bundle
is on ripper at
`.research/leopard-79h/v19-streamed-sources.vCqsXe`, with outer `SHA256SUMS` hash
`0f0489cc57bbf84303a99f78e3ea29b1119439dd16295fad0f5acdc39c2c53c5`.
This is build/lifetime evidence only, not a codec timing or authorization of
the relocated recipe. The initial test-harness allocation failure is retained
with the corrected test results; no memory limit was relaxed.

`v19_source_identity.py` authenticates initial content while borrowing the live
streamed-source and retained-preflight owners. It independently reconstructs
Git blob/tree IDs from held file bytes using 64 KiB buffers, checks the complete
candidate SHA-256/mode/size manifest against the pinned preflight, and binds the
candidate, baseline and initialized `sse2neon` content trees to authenticated
raw commit objects. Read-only commit/index queries use one retained, sealed
`/usr/bin/git` executable and held root/gitdir descriptors; index entries and
default flags must exactly match the independently reconstructed worktrees.
Executable modes use Git's owner-execute bit, not group/other execute bits.
Git status or an index flag cannot excuse different source bytes or Git modes.

This is the exact ordinary detached-clone dialect used by the native diagnostic:
known local origin paths and config bytes, no active hooks, linked metadata,
alternate object stores, grafts, replacement refs or additional Git roots.
Unknown nested `.git` files are not silently excluded. The baseline's committed
`sse2neon` gitlink is bound in its exact uninitialized (empty-directory) state;
the candidate's matching gitlink is initialized and its content authenticated.
Other untracked files and empty worktree directories are rejected. Commit/query
output and individual tree objects are bounded to 1 MiB, as is the sum of
reconstructed path bytes.
The record retains raw commit payloads and reconstructed entries for audit;
later boundaries revalidate source, preflight and Git-tool bytes, and failures
latch. The normal/optimized `leopard2_v19_source_identity_*self_test` CTests use
real miniature Git repositories with a synthetic retained-preflight input.

Authentication does not establish independent source creation, executable
runtime-library closure, atomic observations, mapped-recipe authorization or
permission to run a workload. Those remain separate gates; the borrowed owners
must remain alive, and no historical version's validators are relaxed.

The exact final authentication module passed 21 real-Git cases in both Python
modes, all 10 focused ownership CTests, and a native fresh-build diagnostic
reproducing all four pinned artifacts. The unchanged project-graph files passed
172/172 tests in both modes. Final native peak memory was 518,365,184 bytes,
with all memory-event counters and swap usage zero. That diagnostic released
unused controller heap with process-local `gc.collect()`/`malloc_trim(0)` before
each build; compiler flags and artifact hashes were unchanged. An earlier
untrimmed run reached the limit (`max=3`, no OOM) and was rejected, not reused.
This is one bounded build observation, not a universal memory guarantee.

The 189-entry sealed bundle on ripper at
`.research/leopard-79h/v19-source-identity.uTtcjw` retains the final proof,
both rejected runs, the superseded positive, the heap diagnostic and the
executable-mode regression reproduction. Its outer `SHA256SUMS` hash is
`876bf2ed47f1fadf27f83797c2a9b98b26f3d1cc3489da1acd948329fe2fe688`.
Normal and assertion-disabled local replay verified the retained commit/tree,
index, candidate-manifest, artifact, source and diagnostic-disposition records.
No benchmark executable was run.

`v19_fresh_build.py` now composes the held preflight, host/lock, authenticated
sources and frozen outputs into a dormant fresh-build owner. It creates an
exclusive private container with sibling work and artifact lanes, runs the
fixed independent clone/detached-checkout sequence, and retains source owners
through both one-job builds and artifact freeze. Each cloned directory is held
before checkout or a later nested clone can use it as a writable parent; a
replacement is rejected before that later child can write elsewhere.
Failure directories are kept;
neither the historical build tree nor the caller's parent is overwritten or
removed. This module has no CLI or acquisition dispatch.

Its explicit `leopard2-v19-relocated-build-recipe/v1` binds every stage, argv,
environment override and umask. Only the baseline receives the exact physical
workspace-to-original-root file-prefix map. The candidate recipe is unchanged.
All 32 candidate compile entries are compared against the held, hash-pinned
original compile database after physical-path relocation; all five baseline
entries have explicit expected arguments. Effective cache values and archive/
executable link arguments are checked too. The baseline adapter's Boolean
compile-export cache entry has a separate parser: no candidate or historical
validator is relaxed. Generated make/flag files and metadata remain guarded
and are rehashed after the contained child returns. Metadata is bounded to
128 files, 1 MiB each and 4 MiB total.

Clone children use umask `0002` to preserve the pinned tracked-file modes;
configure/build children use `0022` to produce safe artifact permissions. The
parent's umask is restored on failure as well as success. Each build performs
the previously validated process-local free-heap preparation. Launcher files
use retained descriptors and bounded streaming hashes; the Git executable must
match its pinned hash. Root-owned packaged launchers may have hard links, whose
count remains checked; independently staged source files still require one
link. This is observed launcher-byte identity, **not immutable tool execution
or complete compiler/loader/runtime ownership**.

The normal/optimized `leopard2_v19_fresh_build_*self_test` CTests use synthetic
host/compiler responses with real filesystem descriptors and mutation guards.
They cover stage ordering, exact mapping and cache dialects, link changes,
source/makefile mmap drift, write/restore history, resource loss, redirected
clone destinations, retained failures, cleanup, and false acquisition claims.

The final owner passed 18 cases in both Python modes, the existing ownership
suites, and the unchanged 172-case project-graph checks in both modes. A native
fresh build on ripper completed all ten stages, retained 59 metadata files and
froze four outputs matching the original full-file hashes. Its 512 MiB/no-swap
cgroup peaked at 427,016,192 bytes, with every memory-event counter zero.
The sealed 135-entry evidence bundle is
`.research/leopard-79h/v19-fresh-build.LdVxE3` on ripper; its outer `SHA256SUMS`
hash is `d82634a6b77b9bf08364a3dabbeb15be88ca68d4f19b727bc7d69767cdabdb0e`.
It retains the failed descriptor-limit, launcher and cache-dialect attempts,
the superseded pre-directory-guard positive, and the reproduced source-parent
redirection regression. This establishes a bounded fresh-build path, not a
codec speedup or permission to acquire timings.

These primitives are **not yet connected to acquisition**. Physical
v18-lineage verification, runtime closure, lifetime-owned builder handoff,
and wrapper/controller-closure integration remain required.
Their records explicitly do not assert source/build history, continuous
resource authority, or permission to run a workload. The owning context must
remain alive through later consumers and terminal sealing.

Version 19 is preregistered but deliberately incapable of live acquisition at
this revision. Its canonical preregistration record has
`live_acquisition_armed:false` and SHA-256
`27c1b7d76a0ecdbe194d6e6b62c01e48b1c7d10fc8ef99ebad4d76238669f0c1`.
An invocation of `--conditioned-passive-v19 --attempt N --attempt-budget 2`
validates that fixed point and then exits 2 with
`conditioned v19 authoritative acquisition is NOT ARMED at this preregistration
commit; no scan, lane, or workload was created`. Arming requires a new
preregistration commit. The exhausted `passive-v2` producer also refuses a fresh
`--passive-shared-host` run; the three rejected v18 attempts are frozen by
`v18_failure_lineage_sha256`
`0743214bb7b3a37dff7c00a4f1d302ec300c903c153477c3d57c5a0d10a2780c`.

If a later preregistration arms acquisition, its ordering is fixed as follows:

1. acquire `/tmp/leopard-gf8-authoritative.lock`;
2. freeze identities, build, test, and seal the lane-owned executables;
3. qualify a pair for exactly twelve nominal seven-second windows;
4. acquire a bridge of exactly two nominal one-second windows, with the scan
   tail equal to the bridge head and the bridge tail serving as the campaign
   presample;
5. independently verify the selected pair;
6. acquire its dynamic pair lease and reservation;
7. narrow the controller away from that pair;
8. capture the pre-campaign census; and
9. start the first campaign window no more than 5,000,000,000 ns after the
   bridge tail.

The canonical lock begins before qualification and remains held through
terminal sealing. No pair lease may be held during qualification or the bridge;
it is acquired only after an accepted bridge. Executables are immutable,
lane-owned artifacts whose source identity and SHA-256 are verified both before
and after any future timing.

The attempt budget is exactly two. Attempt one selects the lowest primary from
a fresh qualifying screen. Attempt two either freezes the selected pair with no
fallback, or performs one fresh retry only when attempt one selected no pair.
Reports must say `attempt N of 2`, cite the v17 and v18 rejections, and retain
lanes as `.research/leopard-79h/<commit7>-v19-conditioned-main-a<N>`. Thresholds
and policy may not be retuned. After two attempts the campaign stops; any
successor requires both a new preregistration and a dedicated or
OS-exclusive environment.

The false-claim ceiling is immutable:
`promotion_eligible:false`, `host_exclusivity_proved:false`,
`whole_campaign_interval_observed:false`, and
`causal_performance_claim_allowed:false`. Its SHA-256 is
`c23366107b618bc6666d44675a3b301f49cc5bebd0e0f0a36b8a4bc1d4858257`.
A success terminal is therefore always `NOT_PROMOTED`. The two ratios remain
correlated and must not be multiplied; the strongest permitted classification
is a host/compiler/API/workload-specific conditioned-passive nominal
observation on a shared host.

The frozen resource envelope permits one substantial process, 536,870,912
bytes of memory, zero swap, one Release build job (within an authorized cap of
two), and one sanitizer job. The retained preflight is
`.research/leopard-79h/cf7a705-v19-build-preflight-ripper-a3`, binding candidate
commit `cf7a7056e0bd7f54b8da436a39cae857beab10c1` and exact-main baseline commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` on `ripper`. CPU 0 and sibling 64
are housekeeping; primaries 1 through 63 use sibling offset 64 with no host
fallback or foreign-process mutation.

Version 19 uses `leopard2-main-compare-{raw,manifest,failure}/v19`, independent
audit schema
`leopard2-main-compare-v19-conditioned-passive-independent-audit/v1`, census
policy `leopard2-passive-shared-host-policy/v3`, and generation `passive-v3`.
Its complete and failed cores, success and failure terminals, attempt lineage,
exit status, wrapper preregistration, and controller closure all have separately
versioned v19 schemas. The CTest entries
`leopard2_v19_end_to_end_self_test`,
`leopard2_v19_end_to_end_optimized_self_test`,
`leopard2_v19_wrapper_replay_self_test`, and
`leopard2_v19_wrapper_replay_optimized_self_test` exercise synthetic evidence
only. Wrapper fixtures live under `/tmp/leopard-v19-wrapper-replay.*`; these
tests acquire no campaign lock or lease and perform no scan, build, benchmark,
timing, or workload.

The fresh v16 producer enforces that boundary before capturing executables or
starting timing: the Git-bound `leopard2.cpp` must contain exactly one selector
initializer in the historical default-disabled `2U` form.  Missing, duplicated,
reformatted, stale, or default-on selectors fail closed.  This producer-only
guard does not change structural replay of retained v12-v16 evidence.

For reproduction from a pre-default-on source revision, build Leopard2
separately with production test hooks disabled and the complete version-16
selector set fixed explicitly:

    cmake -S . -B /tmp/leopard2-production -G "Unix Makefiles" \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_OPENMP=ON \
        -DLEOPARD_ENABLE_GF8=ON -DLEOPARD_ENABLE_GF16=ON \
        -DLEO2_BACKEND_VARIANT=auto \
        -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
        -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
        -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF \
        -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF \
        -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF \
        -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=ON \
        -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=ON \
        -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=ON \
        -DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON \
        -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF \
        -DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=ON \
        -DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=ON \
        -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=OFF \
        -DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=ON \
        -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=ON \
        -DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=ON \
        -DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=ON \
        -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0 \
        -DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON \
        -DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON \
        -DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF \
        -DLEO2_FLAG_MAVX512F=FALSE \
        -DLEO2_FLAG_MAVX512BW=FALSE \
        -DLEO2_FLAG_MAVX512VL=FALSE \
        -DLEO2_FLAG_MPREFER_VECTOR_WIDTH_256=FALSE \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /tmp/leopard2-production -j "$(nproc)" \
        --target bench_leopard2

The version-12 compile-command record uses schema v8 and fixes exactly 26
candidate compile actions: 17 production-library sources and nine
non-library actions.  The candidate build-configuration record and generated
configuration file both use schema v6.  Four of the 17 library sources are
separate GF8/AVX2 objects added after the version-11 closure:
`Leopard2BackendAVX2T32B256.cpp`, `Leopard2BackendAVX2T16B64.cpp`,
`Leopard2LowP32B64AVX2.cpp`, and `Leopard2BackendAVX2T2K4.cpp`.  Each has the
exact isolated `-mavx2 -mno-avx512f -falign-functions=64` profile.  The first
three are selected respectively by
`LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=ON`,
`LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=ON`, and
`LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=ON`; the T=2/K=4 object is present
whenever the GF8 AVX2 backend is present.  The T=32 two-block disable selector
must be `OFF`, while the older T=32 generated selector and its disable selector
remain `OFF`.  Changing one of these selectors, merging one of these objects
into another translation unit, or omitting an action changes the version-12
build identity.

The version-13 record advances the compile-command schema to v9 and
the build-configuration record/file schema to v7. It fixes 27 candidate
compile actions: 18 production-library sources and the same nine non-library
actions. The added source is the isolated fixed-shape
`Leopard2BackendAVX2T8K8B1024.cpp` object, compiled with the same pure-AVX2
flags and with `LEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1`. The linked portable core
receives that capability definition, while exact version-12 replay neither
contains the object nor advertises the capability.

The version-14 record advances the compile-command schema to v10 and
the build-configuration record/file schema to v8. It keeps the same 27-action
source closure, requires `LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON` in both the
CMake cache and canonical generated configuration, and requires
`LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1` on the portable `leopard2.cpp`
translation unit. Exact version-13 replay contains none of that new selector
or compile definition.

Version 15 advances the compile-command schema to v11 and compile profile to
v7. It keeps the version-14 source closure and enables the small-dual locator
terms while requiring the regular fallback selector to remain `OFF`. The
T32/B256 generated terminal remains `OFF` in this historical contract.

Current version 16 advances the compile-command schema to v12 and compile
profile to v8. It changes only the production selector identity needed here:
`LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON` in the CMake cache and generated
configuration, with the matching definition on `leopard2.cpp` and the isolated
T32/B256 translation unit. Historical version-15 evidence remains replayable
only with that selector and definition absent.

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

Run from an initial affinity that contains both reserved CPUs.  Create the
report's parent directory first; the supervisor refuses to begin if it cannot
durably fsync both the report and its directory entry.  The
`--preset representative` matrix covers exact-wire GF8 and GF16 high-rate, balanced, XOR,
field-inflation, and larger-parent cases. Custom cells use
`ID:K:R:BYTES:LOSSES:SEED` and must remain within the exact-main API's
`R <= K` restriction. Versions 12 through 16 accept any positive logical shard byte
count.  Leopard2 processes exactly that logical count.  Because exact Leopard
main requires 64-byte shards, its adapter rounds the physical count up to the
next multiple of 64, fills the suffix with zeroes, and fingerprints and checks
only the application-visible logical prefix.  A non-multiple-of-64 cell is
therefore labeled **application-equivalent**: the codeword and recovered data
match over the requested logical bytes, but exact main deliberately processes
additional zero-padded physical bytes.  It is not a same-physical-byte-work
comparison.  Historical version-11 cells retain their original requirement
that shard bytes be multiples of 64.

Versions 12 through 16 bind one explicit `--candidate-mode`: `auto`, `generic`,
`materialized`, or `tiled` (the default is `auto`). The two forced workspace
modes also force the specialized decoder. The runner verifies the exact child
arguments and all four emitted force booleans, so evidence for one path cannot
be relabeled as another.

    taskset -c 0-31 python3 -I -S \
        tools/leopard2_affinity_supervisor.py run \
        --report /tmp/leopard2-vs-main-affinity.json \
        --reserved-cpus 15,31 -- \
        python3 \
        experiments/leopard2/main_compare/run_abba.py run \
        --baseline /tmp/leopard-main-pure-avx2/leopard_main_benchmark \
        --candidate /tmp/leopard2-production/bench_leopard2 \
        --baseline-archive /tmp/leopard-main-pure-avx2/libleopard_main_exact.a \
        --candidate-archive /tmp/leopard2-production/libleopard.a \
        --baseline-build-dir /tmp/leopard-main-pure-avx2 \
        --candidate-build-dir /tmp/leopard2-production \
        --baseline-source-root /tmp/leopard-main-exact \
        --candidate-source-root "$PWD" \
        --candidate-commit "$(git rev-parse HEAD)" \
        --baseline-pure-avx2 \
        --candidate-mode auto \
        --reservation-file build/leopard2-main-reservation.json \
        --output /tmp/leopard2-vs-main --cpu 15 --reserved-sibling 31 \
        --preset representative --reuse 8 --iterations 9 --warmup 2

The supervisor is an unprivileged, same-user aid, not an OS-exclusive CPU
shield. Before its first mutation it requires every existing process to have a
uniform original mask across that process's threads. It journals every mask
before changing it. A new thread may be restored only when its exact original
mask follows from that uniform process or a still-identifiable uniform creator.
A newcomer whose already-restricted mask has
no such provenance invalidates the run and is boundedly terminated, rather
than being left with a guessed or widened mask.
Every run and recovery also holds one persistent same-UID lock under
`/run/user/UID`; its exact inode is part of the journal. This prevents two
cooperating supervisors from racing each other's original masks.

Launch the supervisor itself with `python3 -I -S` as shown above. Every evidence
operation fails closed without those flags, preventing ambient Python startup
hooks from running before the supervisor can establish its own invariants.

The runner is launched through a fork/pipe gate. The child creates its new
session but cannot execute benchmark code until its PID,
start-time, procfs directory inodes, session, command, boot ID, and
PID-namespace identity are durable in
the report. The executable and optional Python script are copied into
write-sealed memfds whose identities, seals, sizes, and hashes are bound in the
v10 through v13 reports; the working directory is held by descriptor. Execution
uses those
immutable content snapshots, so neither pathname replacement nor same-inode
overwrite can select different code. Before copying a source that the current
UID could write, the supervisor must acquire a Linux read lease. Capture fails
closed while any writable descriptor or writable `MAP_SHARED` mapping exists;
a later write-open lease-break request is observed through blocked `SIGIO` and
also invalidates capture. Sources independently proven inaccessible to same-UID
writers record that distinct guard policy. This closes dirty-page mutations for
which inode timestamps need not change. A small fixed Python bootstrap preserves
the original script's `__file__`, `sys.argv`, and `sys.path[0]`. Python starts
with `-I -S`, so ambient `PYTHONPATH`, `PYTHONHOME`, user-site, and
`sitecustomize` code cannot run before that bootstrap. The exec boundary does
not inherit the supervisor's environment. It constructs the exact environment
recorded by the main-comparison runner and adds only the report-bound execution
nonce. The current version-16 child environment is exactly `LANG=C`, `LC_ALL=C`,
`OMP_DYNAMIC=FALSE`, `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`,
`OMP_PROC_BIND=TRUE`, `PATH=/usr/bin:/bin`, and `TZ=UTC`; it deliberately does
not set `OMP_PLACES`.  Retained version-11 evidence continues to replay with
its historical `OMP_PLACES=cores` entry and without `OMP_THREAD_LIMIT`.
Consequently ambient Python hooks and dynamic loader variables such as
`LD_PRELOAD`, `LD_LIBRARY_PATH`, and `LD_AUDIT` cannot inject unrecorded code.
The v10 through v13 reports bind this strict base environment. The child also
verifies each live snapshot identity and seal set immediately before
`execve`. The supervisor requests Linux `F_SEAL_EXEC` and tests it when the
kernel supports that seal; an atomic `EINVAL` fallback retains all mandatory
write/grow/shrink/seal protections, with executable mode rechecked at the
launch boundary. The supported command grammar is either one captured native
ELF executable or a captured ELF interpreter with a `.py` script directly in
`argv[1]`. A shebang and an option-prefixed or later Python script are rejected,
so neither can redirect execution through an uncaptured interpreter or source.
Arguments after a captured direct script remain ordinary application arguments.
The gated child
establishes and verifies launch affinity directly rather than executing an
unbound `taskset` pathname. Parent EOF aborts the gate, which closes the SIGKILL
window before that journal entry. Before writing the gate byte, the supervisor
durably changes the child flag to mean "release may have occurred"; crash
recovery must therefore reject the scope even if the supervisor dies before it
observes release completion. Inherited Python SIGINT, SIGTERM, and SIGHUP
handlers are replaced by default child dispositions while those signals are
blocked, before the gate can wait or release. On every ordinary exit the
supervisor boundedly signals,
then SIGKILLs if necessary, every surviving process in the recorded child
session and every retained descendant that created another session
before restoring any same-user mask. This lets the main runner use its own
bounded per-invocation sessions without those timed children being mistaken for
unrelated work. SIGINT, SIGTERM, and SIGHUP observed before the explicit
post-cleanup signal boundary invalidate evidence and use the same bounded
cleanup. Signals after that boundary retain the caller's original disposition;
they cannot retroactively alter an already-finished benchmark transaction.
Failure to inspect pending signals at the start or either final durability
boundary is itself a terminal evidence failure, but it cannot bypass child
cleanup, affinity restoration, terminal journaling, immutable-descriptor close,
or exact caller-mask restoration.
The supervisor is a Linux child subreaper and retains exact descendant
identities, so a double-forked child remains a cleanup target after reparenting.
Restoration performs the full bounded rescan budget after any mutation and
repairs late threads it observes. This is deliberately not certified as globally
complete: `clone(2)` can copy a restricted mask and publish its task after any
finite last `/proc` scan. Consequently a v10 or v11 report is authoritative only when
all surrounding same-UID work was already outside the reserved CPUs and no
non-supervisor mask was changed. The bounded exception is the supervisor's own
single main thread: the implementation controls its sole fork point and creates
no other thread. A run that changes any outside mask finishes as failed evidence
after best-effort restoration. Likewise, crash recovery of a
released child remains failed even after known identities are cleaned, because
the replacement process is no longer the original subreaper and cannot exclude
an unjournaled cross-session orphan. A future cgroup-backed containment mode can
lift these restrictions when it can prove persistent membership.
The 25 ms monitor is read-only on an unchanged scan. It does not rewrite or
fsync the journal for each poll, for ordinary child-process churn, or merely to
retain audit-only provenance; it writes only when recovery-relevant state or an
actual anomaly changes.

The wrapper must finish with an accepted v13 report and its separate
`.accepted.json` commit seal. An `accepted: true` body without this
hash-bound, directory-fsynced seal is not evidence. The
wrapper restores the caller's signal handlers and exact entry mask before
installing the seal, then reads the exact report/seal pair back. Once that
read-back commits acceptance, a later coordinator-lock teardown error resolves
to the committed result rather than returning failure beside usable evidence.
The ordinary ABBA
manifest must independently pass its unchanged zero-sibling-jiffy verifier
before timings are usable. If the supervisor itself is killed with SIGKILL,
restore its journal before continuing:

    python3 -I -S tools/leopard2_affinity_supervisor.py restore \
        --report /tmp/leopard2-vs-main-affinity.json

Recovery is deliberately limited to the same kernel boot and PID namespace.
It first verifies and empties the recorded child session, then retries both
pending and previously failed restores. A reboot or namespace transition is a
hard refusal: numeric PIDs/TIDs no longer name the journaled objects. The tool
never changes another UID's threads and deliberately does not weaken the
runner's rejection of kernel, other-user, or otherwise unexcluded sibling
activity. Linux has no pidfd affinity operation for non-leader TIDs; a detected
identity change across a numeric-TID affinity syscall therefore leaves an
explicit unrecoverable uncertainty and can never produce accepted evidence.

After both files verify, create a canonical joint binding. It authenticates the
accepted supervisor report hash and schema, the exact command identity, the
main manifest and raw-bundle hashes, the runner source identity, CPU pair, all
numeric campaign parameters, build/source paths, reservation, and output path.
Versions 5 through 16 additionally carry a fresh 256-bit nonce from the gated supervisor
into the runner and bind the same runner PID, enclosing monotonic intervals,
launch and reserved CPU sets, campaign hash, and held reservation payload. A
compatible campaign from a different supervisor execution cannot be rebound.
These same-execution guarantees are retained by the current version-7
joint-binding format. Version 2 binds reports for main-comparison versions
5 through 11, version 3 binds version 12, version 4 binds version 13,
version 5 binds version 14, version 6 binds version 15, and version 7 binds
version 16;
historical bindings are not silently upgraded.

    python3 -I -S tools/leopard2_affinity_supervisor.py verify-report \
        --report /tmp/leopard2-vs-main-affinity.json
    python3 -I -S tools/leopard2_affinity_supervisor.py bind \
        --report /tmp/leopard2-vs-main-affinity.json \
        --manifest /tmp/leopard2-vs-main/manifest.json \
        --output /tmp/leopard2-vs-main/affinity-binding.json
    python3 -I -S tools/leopard2_affinity_supervisor.py verify-binding \
        --binding /tmp/leopard2-vs-main/affinity-binding.json

Verify while the exact build inputs still exist:

    python3 experiments/leopard2/main_compare/run_abba.py verify \
        --manifest /tmp/leopard2-vs-main/manifest.json \
        --affinity-binding /tmp/leopard2-vs-main/affinity-binding.json

The current version-16 verifier recomputes the pair-lock identity, every per-field CPU
counter delta, the zero-non-idle-sibling decision, workload identities, and all
statistics. Versions 3 through 16 also bind the canonical CMake target `leopard`, archive
`libleopard.a`, and `leopard.dir` dependency closure. It retains the exact,
bounded UTF-8 archive link-recipe content, binds its byte length and SHA-256 to
the recipe-file identity, and parses those bytes to require the declared
archive, ordinary target directory, and matching `ranlib` command. Retained
version-2 bundles replay with their original `libleopard`/`liblibleopard.a`
identity, record shape, and isolation semantics. Version-3 bundles replay as
AUTO-only evidence without retroactively acquiring the decoder-mode field;
version-4 bundles retain their original hardened archive closure. Current
version 16 requires exact source, compiler, CMake, compile/object, archive,
member, executable-link, output, tool, and runtime-dependency record shapes in
offline verification. The executable link is consumed by a fail-closed grammar:
only the canonical compiler, value-free build flags, one benchmark object, one
project archive, one output, and the explicit GCC OpenMP and pthread operands
are accepted. The resolved libgomp and libpthread files retain complete byte
identities and the complete lexical symlink chain. Verification reads each
bounded file into a private immutable byte snapshot, consumes its digest and
file-format checks from that snapshot, then re-resolves and directly compares
the symlink, descriptor, mode, and terminal-path identities before closing the
external-input closure. Baseline and candidate records must be identical, and
the external inputs may not be newer than either executable.

The compile-command proof has its own versioned schema and explicit baseline or
candidate profile. Every configured compile-database entry must match its exact
producer-known ordered argv, working directory, source, and CMake output
operand. For every required translation unit the proof additionally retains
that contract with the strict source/object byte identities. A target-only
build may leave unrelated configured objects absent; only their existence is
conditional, not their compile-command contract.
There must be exactly one positional source, exactly one `-c`, and exactly one
`-o`; the source and output must match both the compile-database entry and the
producer-known object path. Definitions, includes, C++ language mode,
optimization sequence, OpenMP flags, and per-TU ISA isolation flags must equal
the profile in order. There is no catch-all `-D` or `-I` acceptance. The
retained CMake Release flags and executable recipe remain context-specific
closed grammars. GCC/Clang long aliases, unknown optimization controls,
sanitizer, profile, LTO, instrumentation, inline/vector disables, coverage,
and last-option overrides fail closed.
Response files, linker forwarding, scripts, search/library
controls, specs, alternate tool roots/loaders, plugins, and undeclared inputs
are rejected. The producer-known translation-unit sets are exact: all five
baseline and all 27 current candidate CMake actions are fixed. Those candidate
actions comprise 18 production-library sources and nine non-library actions;
the required linked artifact closure has all 18 archive members plus the
benchmark object;
a coherent source/object/member/link truncation is rejected. Each source retains
the bounded raw Git commit object, whose Git SHA-1 and named tree are recomputed
offline. Version 5 retains the bounded raw `ldd` output for historical replay.
Versions 12 through 16 strictly parse every complete `ldd` line and retain a bounded,
versioned canonical transcript that replaces only each terminal ASLR load
address with the literal `<ASLR_LOAD_ADDRESS>` token. SONAMEs, resolved paths,
line kinds, and line order remain intact, so immediate snapshots of unchanged
binaries are stable while a missing address, path corruption, line removal, or
raw/summary relabel fails closed. The transcript and dependency summary must
agree, include at least one shared library, and identify exactly one dynamic
loader. Every hashed dependency file retains the
exact canonical loader-declared path through which it was read, so swapping two
otherwise valid shared-library records is rejected offline. This retained pair
does not independently prove that an ordinary `ldd` line was never removed by a
writer able to coherently rewrite and re-sign both the canonical output and its
summary; current-input verification detects such a change against the live
closure.

Its host identity includes the exact sysfs online CPU/node list bytes, every
retained cache index and cache record, a separately retained canonical listing
of each CPU's sysfs cache directory, NUMA-node directory inventory and placement,
shared-CPU domain, and heterogeneous-core class. The cache hierarchy, summary,
and listing must agree three ways, but remain one internally signed evidence
bundle rather than an independently authenticated snapshot of sysfs. Campaign,
cell, CPU, reservation, and host count fields are bounded exact integers; JSON
booleans are never accepted as integer values. Removing any schema-required
field, or changing only one of the retained paired representations, fails even
with `--no-current-input-check`. Version-1 bundles remain
structurally replayable without acquiring later isolation semantics
retroactively. A schema/path relabel that leaves the historical recipe bytes
unchanged is rejected.

Bundle digests are unkeyed integrity checks, not an independent authenticity
anchor. They detect edits relative to the retained evidence, and the semantic
recipe binding prevents relabeling old recipe bytes under versions 3 through 16. They
cannot prevent a hostile writer from replacing every evidence byte and
recomputing every internally consistent digest. Preserve evidence in an
immutable or independently authenticated store, or add an external digital
signature or transparency log, when adversarial provenance is in scope. A
sibling-activity or child failure produces signed
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
versions 12 through 16 pass `--measure-one-shot-decode` to the candidate and define
`decode_first_use` as the median directly timed public `leo2_decode` one-shot
call, including its plan setup with the codec already created. It is not the
sum of separately measured plan and execution medians. The separate
`decode_reuse_amortized` metric remains the median plan execution plus median
plan-create divided by the requested reuse count; codec setup is excluded from
both metrics. Confidence intervals use one clustered mean log contrast per ABBA
round (three observations, Student-t with two degrees of freedom), not six
falsely independent adjacent pairs.

Retained version-11 manifest/raw/failure bundles remain replayable by `verify`
and `verify-failure`.  They keep compile-command schema v7, build-configuration
record/file schema v5, their 13-library-source/22-compile-action candidate
closure, 64-byte-multiple cell policy, and historical child environment.  They
are not retroactively assigned the four new isolated AVX2 objects, arbitrary
logical tails, build-configuration v6, or compile-command v8 semantics.

Retained version-12 bundles likewise remain exact historical evidence. They
keep compile-command schema v8, build-configuration record/file schema v6, and
their 17-library-source/26-action closure; verification never relabels them as
version 13 or inserts the fixed-shape K8 object or capability.

Retained version-13 bundles keep compile-command schema v9,
build-configuration record/file schema v7, and their 18-library-source/27-action
closure. Verification never relabels them as version 14 or inserts the
small-dual-direct selector or portable-core compile definition.

Retained version-14 bundles keep compile-command schema v10, compile profile
v6, build-configuration record/file schema v8, and the small-dual-direct
selector. Verification never inserts the version-15 locator-term selectors.

Retained version-15 bundles keep compile-command schema v11, compile profile
v7, build-configuration record/file schema v9, the locator-term selectors, and
the historical `LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF` identity.
Verification never relabels them as version 16 or inserts the generated
T32/B256 definitions.

Retained version-9 bundles also keep their original statistics contract:
`decode_first_use` is the derived sum of separately timed plan-create and one
execution medians. Verification never upgrades a version-9 bundle to
version-12, version-13, version-14, version-15, or version-16 pure-AVX2,
effective-AVX2, or directly timed one-shot semantics.
New `run` campaigns always emit version 16.
