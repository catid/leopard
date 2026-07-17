# Reproducing the LOW_V1 coefficient-copy comparison

## Claim and compared sources

This campaign measures the encoder-only effect of replacing each LOW_V1
parity block's full `P`-shard coefficient copy with an out-of-place first
butterfly layer.  The fixed comparison is:

- control: `4070e4e527935026fb87593567587558f0a08d51`;
- candidate: `6d3afee213b94d486cf5f8145ac18078883ebc20`.

Both sides are Leopard2.  This is not the separate Leopard-main comparison.
The collector rejects other commits, dirty source trees, non-Release builds,
enabled tests/fuzzers/CUDA, mismatched compilers or effective flags, missing
compile commands, instrumentation, non-production backend flags, incomplete
object/archive/link closure, unreproducible linked executables, or changed
runtime dependencies. Provenance validation binds the Python coordinator, Git,
compiler, archiver, ranlib, every required source/object pair, each archive
member's bytes, the exact normalized link recipe and external link inputs, and
a byte-identical clean recompilation of every required object followed by a
byte-identical clean-room relink of each benchmark. The clean recompilation
executes the exact retained compiler command and working directory with only
its output redirected to a private temporary file; a digest or size difference
invalidates the run before measurement.

The authoritative bundles use `leopard2-low-encode-copy-raw/v3` and
`leopard2-low-encode-copy-manifest/v3`. Failed-run diagnostics use
`leopard2-low-encode-copy-failure/v3`, whose signed lifecycle is an exact prefix
of the declared phases and records the failed phase plus teardown result.

The retained benchmark output is schema `leopard2-benchmark-v2`.  Every child
must retain its raw timing samples, resolve the explicitly requested LOW_V1
field/backend, pass its round trip, and emit identical original/parity/recovery
digests on both sides.  Codec setup and encode execution are analyzed
separately; setup is reported but does not decide promotion.

## Fixed matrix and policy

Every one of scalar, SSSE3, and AVX2 runs these eight aligned shard sizes:

| Role | Field | K | R | Bytes |
|---|---:|---:|---:|---:|
| target | GF8 | 8 | 248 | 64 |
| neighbor | GF8 | 16 | 240 | 256 |
| target | GF8 | 32 | 224 | 1,024 |
| target | GF8 | 64 | 192 | 4,096 |
| target | GF16 | 100 | 156 | 16,384 |
| neighbor | GF16 | 127 | 129 | 65,536 |
| neighbor | GF16 | 128 | 128 | 262,144 |
| neighbor | GF16 | 129 | 127 | 1,048,576 |

Thus the fixed matrix has 24 cells.  Each cell runs five independent
`A1/B1/B2/A2` rounds.  An invocation contributes its median retained sample;
the two paired log contrasts in a round are averaged, and the five independent
round contrasts use a two-sided 95 percent Student-t interval (`df=4`).  A
target passes when the lower confidence bound on control-time/candidate-time is
at least 1.05.  A neighbor passes when that lower bound is at least `1/1.02`,
the predeclared no-regression threshold.  One failed cell makes the manifest a
valid `policy_failed` result and makes verification exit 2.  Performance never
changes the structural validity of the evidence.

Loss coordinates use the benchmark's xorshift64 sequence (shifts 13, 7, and
17). The coordinator self-test checks hard-coded vectors produced by an
independently compiled C++ implementation before it invokes the Python mock,
so a shared coordinator/mock RNG defect cannot authenticate itself.

Benchmark-v2's `original_data` digest covers all K source shards, while
`recovered_originals` covers only the L missing originals. They are deliberately
different digest domains. The harness checks every repaired byte against its
source internally; the collector requires that checked round trip and compares
all three digests independently across every control and candidate invocation.

## Isolation and publication contract

Launch affinity must contain exactly the chosen measured SMT pair plus at
least one housekeeping CPU.  The collector verifies Linux topology reports
exactly two siblings, moves itself to housekeeping CPUs, and pins each child to
one measured logical CPU.  It acquires the common pair-wide lease at
`/run/user/UID/leopard2-cpu-leases/leopard2-cpu-pair-UID-A-B.lock`, so it cannot
overlap exact-main, butterfly, Jerasure, or other audited Leopard2 collectors
on that physical core. A Linux abstract Unix-domain socket is held for the same
CPU-pair/root identity for the full lease lifetime. It remains exclusive even
if the diagnostic filesystem lock or its containing directory is renamed or
replaced; the filesystem inode lock is retained for reviewable identity and
continuous tamper checks.

The coordinator reservation is a canonical JSON file without a trailing
newline.  Its parent must be owned mode 0700 and the file must be an owned,
single-link mode-0600 regular file.  The collector takes an exclusive lock and
binds the parent and file device/inode, mode, ownership, link count, payload,
and digest. Both reservation and pair lease are revalidated before and after
every measured child. The coordinator's affinity is also required to equal the
housekeeping set at those boundaries. The reservation also has a stable Linux
abstract-socket lease, so replacing its path cannot create an overlapping
reservation.

After the initial attestation, the collector copies each benchmark into a
private staged file, revalidates its bytes, opens it read-only, unlinks its
name, and executes the inherited descriptor through `/proc/self/fd`. A source
path replacement after attestation therefore cannot change the program that a
child executes. Every invocation retains and revalidates that immutable
execution identity.

Pre/post `/proc/stat` evidence covers the first eight non-double-counted Linux
CPU counters.  The measured CPU must accrue non-idle work, the reserved sibling
must accrue time, and the sibling must accrue exactly zero non-idle jiffies.
These counters cannot attribute every jiffy to the child, so the external host
reservation and an otherwise idle machine remain required.

Children and build/provenance helpers have timeouts and separately bounded
stdout and stderr collection. On Linux the runner temporarily becomes a child
subreaper, follows process identities and ancestry through procfs, signals by
pidfd, reaps adopted children, and requires two empty scans before restoring the
prior subreaper state. This contains descendants that call `setsid()` or
double-fork, not just the initial process group. Platforms without those Linux
facilities fail closed before an authoritative child is spawned; portable
evidence replay remains available. Every cleanup deadline is finite and no
unbounded `wait` is used. Retained JSON has byte, depth, node, string, integer, floating-point,
path, and collection limits. Retained-artifact hashing is bounded, rejects
symlinks, hardlinks, FIFOs and other special files before a nonblocking open,
binds the open descriptor to the named inode, and rejects concurrent
metadata or size changes. The retained-file inventory rejects symlinks,
unexpected names, empty directories, excess entries, and files outside the
declared caps.

Output is first written under a private mode-0700 staging directory and fully
validated. Both lock contexts have exited, immutable executable descriptors
are closed, and the coordinator's exact launch affinity is restored before
success can become visible. The stage tree and output parent are fsynced before
Linux `renameat2(RENAME_NOREPLACE)`, which is the final publication operation.
There is no failable validation or teardown after that rename, and a run
cannot replace an existing result. A pre-publication or teardown failure
instead publishes signed lifecycle diagnostics; it can never leave an
acceptable success manifest visible.

## Build the exact inputs

The paths below match the runner defaults.  Do not rebuild while a campaign is
running.

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi

    git -C /home/catid/leopard worktree add --detach \
      /home/catid/leopard-wt-low-copy-baseline \
      4070e4e527935026fb87593567587558f0a08d51
    git -C /home/catid/leopard worktree add \
      -b reproduce/leopard2-low-copy-candidate \
      /home/catid/leopard-wt-low-copy \
      6d3afee213b94d486cf5f8145ac18078883ebc20

    cmake -S /home/catid/leopard-wt-low-copy-baseline \
      -B /home/catid/leopard-wt-low-copy-baseline/build/low-copy-authoritative \
      -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DENABLE_OPENMP=ON \
      -DLEO2_BACKEND_VARIANT=auto -DLEO2_BUILD_BENCHMARKS=ON \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_FUZZERS=OFF \
      -DLEO2_ENABLE_CUDA=OFF
    cmake --build \
      /home/catid/leopard-wt-low-copy-baseline/build/low-copy-authoritative \
      -j "$JOBS"

    cmake -S /home/catid/leopard-wt-low-copy \
      -B /home/catid/leopard-wt-low-copy/build/low-copy-authoritative \
      -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DENABLE_OPENMP=ON \
      -DLEO2_BACKEND_VARIANT=auto -DLEO2_BUILD_BENCHMARKS=ON \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_FUZZERS=OFF \
      -DLEO2_ENABLE_CUDA=OFF
    cmake --build \
      /home/catid/leopard-wt-low-copy/build/low-copy-authoritative \
      -j "$JOBS"

The candidate worktree must be clean at exactly
`6d3afee213b94d486cf5f8145ac18078883ebc20`; the control is required to be
clean and detached at its declared commit.

## Fast non-timing validation

Run this while other jobs are active if desired.  It never invokes either real
benchmark:

    python3 experiments/leopard2/low_encode_copy/run_abba.py self-test

Then close the mock/production semantic gap with the tiny real-harness smoke.
It executes benchmark timing machinery only as a correctness/schema probe and
does not retain or report timings, compare speed, reserve a core, or constitute
authoritative evidence:

    python3 experiments/leopard2/low_encode_copy/run_abba.py production-smoke \
      --backend scalar --backend ssse3 --backend avx2

## Run the authoritative campaign

Do this only after all other workers and memory-intensive processes stop.
Select a physical SMT pair from the process's allowed affinity, and retain at
least one other allowed CPU for housekeeping.  The example uses CPUs 0 and 16;
replace them with a pair actually reported by this host.

    CPU=0
    SIBLING=16
    RESERVATION_DIR="$PWD/.research/leopard2/low-copy-reservation"
    RESERVATION="$RESERVATION_DIR/reservation.json"
    OUTPUT="$PWD/.research/leopard2/low-copy-final-$(date -u +%Y%m%dT%H%M%SZ)"
    install -d -m 700 "$RESERVATION_DIR"
    python3 - "$CPU" "$SIBLING" "$RESERVATION" <<'PY'
    import json, secrets, sys
    from pathlib import Path
    cpu, sibling, path = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
    payload = {
        "benchmark_cpu": cpu,
        "nonce": secrets.token_hex(32),
        "owner": "leopard2-low-copy-coordinator",
        "reserved_sibling": sibling,
        "schema": "leopard2-cpu-reservation/v1",
        "status": "held",
    }
    path.write_bytes(json.dumps(payload, sort_keys=True,
                                separators=(",", ":")).encode())
    path.chmod(0o600)
    PY

    HOUSEKEEPING="$(python3 - "$CPU" "$SIBLING" <<'PY'
    import os, sys
    pair = {int(sys.argv[1]), int(sys.argv[2])}
    print(min(set(os.sched_getaffinity(0)) - pair))
    PY
    )"
    taskset -c "$CPU,$SIBLING,$HOUSEKEEPING" \
      python3 experiments/leopard2/low_encode_copy/run_abba.py run \
      --cpu "$CPU" --reserved-sibling "$SIBLING" \
      --reservation-file "$RESERVATION" --output "$OUTPUT" \
      --iterations 5 --warmup 1 --reuse 1 --timeout 180

The campaign deliberately uses fewer than all available cores: it is a
cache-sensitive single-core measurement whose sibling must remain idle.  Resume
full-machine parallel work immediately afterward.

## Replay

Replay with live build/source closure checking:

    python3 experiments/leopard2/low_encode_copy/run_abba.py verify \
      --manifest "$OUTPUT/manifest.json"

Portable structural inspection after moving the evidence to another machine:

    python3 experiments/leopard2/low_encode_copy/run_abba.py verify \
      --manifest "$OUTPUT/manifest.json" --no-current-input-check

Live replay exits 0 for a policy pass and 2 for valid negative evidence.
Portable inspection deliberately exits 1 even when internally consistent,
because it cannot revalidate the declared source/build/tool closure and is not
authoritative evidence. A malformed or inconsistent bundle also exits 1. A
failed run publishes `failure.json`; verify it with:

    python3 experiments/leopard2/low_encode_copy/run_abba.py verify-failure \
      --failure "$OUTPUT/failure.json"

That last command exits 1 by design because diagnostics are never valid timing
evidence.
