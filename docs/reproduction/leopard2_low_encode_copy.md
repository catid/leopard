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
a byte-identical clean-room relink of each benchmark.

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

## Isolation and publication contract

Launch affinity must contain exactly the chosen measured SMT pair plus at
least one housekeeping CPU.  The collector verifies Linux topology reports
exactly two siblings, moves itself to housekeeping CPUs, and pins each child to
one measured logical CPU.  It acquires the common pair-wide lease at
`/run/user/UID/leopard2-cpu-leases/leopard2-cpu-pair-UID-A-B.lock`, so it cannot
overlap exact-main, butterfly, Jerasure, or other audited Leopard2 collectors
on that physical core.

The coordinator reservation is a canonical JSON file without a trailing
newline.  Its parent must be owned mode 0700 and the file must be an owned,
single-link mode-0600 regular file.  The collector takes an exclusive lock and
binds the parent and file device/inode, mode, ownership, link count, payload,
and digest. Both reservation and pair lease are revalidated before and after
every measured child. The coordinator's affinity is also required to equal the
housekeeping set at those boundaries and around evidence publication.

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

Children have separately bounded stdout and stderr collection. Retained JSON
has byte, depth, node, string, integer, floating-point, and collection limits.
The retained-file inventory rejects symlinks, unexpected names, empty
directories, excess entries, and files outside the declared caps. Output is
first written under a private staging directory, fully validated and fsynced,
and then published with Linux `renameat2(RENAME_NOREPLACE)`. A run cannot
replace an existing result.

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
