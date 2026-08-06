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
historical `-march=native` policy. Current version-14 ABBA evidence requires
this pure-AVX2 profile and the matching `--baseline-pure-avx2` runner selector;
the native build above is useful for historical reproduction but is not a
version-14 baseline. A suitable baseline build is:

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

The current all-K run-contract and manifest generation is v12. It requires the
production-default exact T32/B256 terminal selector to be `ON`. Historical v11
remains replayable only with that selector `OFF`; relabeling either body as the
other generation is rejected. The 2,522-cell matrix uses 4 KiB and 64 KiB
shards, so the exact 256-byte terminal cannot affect a timed cell, but its value
still belongs to the authenticated production selector identity.

## Counterbalanced comparison

`run_abba.py` supplies the other half of the comparison. Its current evidence
contract is raw/manifest/failure version 14. It accepts only fresh Release
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

Versions 12 through 14 are effective-AVX2 comparisons, not merely requests for the AUTO
backend. Leopard main must use the pure-AVX2 build above. Leopard2 must use
`LEO2_BACKEND_VARIANT=auto`, and the benchmark must resolve that AUTO request to
`avx2`. The AVX-512F/BW/VL and preferred-vector-width probes are forced false,
so `Leopard2BackendAVX512.cpp` cannot enter the build closure. The
VEX-encoded GFNI translation unit may remain in the archive with
`-mavx2 -mgfni -mno-avx512f`; GFNI remains explicit-only in this profile and
AUTO may not select it. The runner attests the CMake cache, complete compile
graph, per-TU
ISA flags, archive members, linked executable, and reported runtime backend,
so a native, AVX-512, or AUTO-GFNI candidate cannot be relabeled as version-14
AVX2 evidence.

Build Leopard2 separately with production test hooks disabled and the complete
version-14 selector set fixed explicitly:

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
        -DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF \
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

The current version-14 record advances the compile-command schema to v10 and
the build-configuration record/file schema to v8. It keeps the same 27-action
source closure, requires `LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON` in both the
CMake cache and canonical generated configuration, and requires
`LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1` on the portable `leopard2.cpp`
translation unit. Exact version-13 replay contains none of that new selector
or compile definition.

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
`R <= K` restriction. Versions 12 through 14 accept any positive logical shard byte
count.  Leopard2 processes exactly that logical count.  Because exact Leopard
main requires 64-byte shards, its adapter rounds the physical count up to the
next multiple of 64, fills the suffix with zeroes, and fingerprints and checks
only the application-visible logical prefix.  A non-multiple-of-64 cell is
therefore labeled **application-equivalent**: the codeword and recovered data
match over the requested logical bytes, but exact main deliberately processes
additional zero-padded physical bytes.  It is not a same-physical-byte-work
comparison.  Historical version-11 cells retain their original requirement
that shard bytes be multiples of 64.

Versions 12 through 14 bind one explicit `--candidate-mode`: `auto`, `generic`,
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
v10 and v11 reports; the working directory is held by descriptor. Execution
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
nonce. The current version-14 child environment is exactly `LANG=C`, `LC_ALL=C`,
`OMP_DYNAMIC=FALSE`, `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`,
`OMP_PROC_BIND=TRUE`, `PATH=/usr/bin:/bin`, and `TZ=UTC`; it deliberately does
not set `OMP_PLACES`.  Retained version-11 evidence continues to replay with
its historical `OMP_PLACES=cores` entry and without `OMP_THREAD_LIMIT`.
Consequently ambient Python hooks and dynamic loader variables such as
`LD_PRELOAD`, `LD_LIBRARY_PATH`, and `LD_AUDIT` cannot inject unrecorded code.
The v10 and v11 reports bind this strict base environment. The child also
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

The wrapper must finish with an accepted v11 report and its separate
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
Versions 5 through 14 additionally carry a fresh 256-bit nonce from the gated supervisor
into the runner and binds the same runner PID, enclosing monotonic intervals,
launch and reserved CPU sets, campaign hash, and held reservation payload. A
compatible campaign from a different supervisor execution cannot be rebound.
These same-execution guarantees are retained by the current version-5
joint-binding format. Version 2 binds reports for main-comparison versions
5 through 11, version 3 binds version 12, version 4 binds version 13, and
version 5 binds version 14;
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

The current version-14 verifier recomputes the pair-lock identity, every per-field CPU
counter delta, the zero-non-idle-sibling decision, workload identities, and all
statistics. Versions 3 through 14 also bind the canonical CMake target `leopard`, archive
`libleopard.a`, and `leopard.dir` dependency closure. It retains the exact,
bounded UTF-8 archive link-recipe content, binds its byte length and SHA-256 to
the recipe-file identity, and parses those bytes to require the declared
archive, ordinary target directory, and matching `ranlib` command. Retained
version-2 bundles replay with their original `libleopard`/`liblibleopard.a`
identity, record shape, and isolation semantics. Version-3 bundles replay as
AUTO-only evidence without retroactively acquiring the decoder-mode field;
version-4 bundles retain their original hardened archive closure. Current
version 14 requires exact source, compiler, CMake, compile/object, archive,
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
Versions 12 through 14 strictly parse every complete `ldd` line and retain a bounded,
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
recipe binding prevents relabeling old recipe bytes under versions 3 through 12. They
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
versions 12 through 14 pass `--measure-one-shot-decode` to the candidate and define
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

Retained version-9 bundles also keep their original statistics contract:
`decode_first_use` is the derived sum of separately timed plan-create and one
execution medians. Verification never upgrades a version-9 bundle to
version-12, version-13, or version-14 pure-AVX2, effective-AVX2, or directly timed one-shot
semantics.
New `run` campaigns always emit version 14.
