# Leopard2 / Intel ISA-L comparison checkpoint

This is a bounded external comparison, not a wire-compatibility result and not
a replacement for the required Leopard2 benchmark matrix. Intel ISA-L and
Leopard2 both receive the same public `K`, `R`, shard byte count, deterministic
source bytes, missing-original indices, batch size, plan/table reuse count, and
single-thread execution constraint. Their parity bytes and generator matrices
differ. Only public input, generated-output, and repaired-output throughput is
comparable.

## Isolation and exact source identity

The adapter is disabled by construction: the production `CMakeLists.txt` never
descends into `experiments/leopard2/isal_compare`. That standalone project must
be configured explicitly and links a private static ISA-L archive. A normal
Leopard build does not find, compile, link, install, or load ISA-L. The CTest
`leopard2_isal_comparison_self_test` validates policy and artifact parsers using
synthetic data only, so it also runs on machines with no ISA-L or NASM.

The checkpoint uses:

| Component | Exact identity | License / hash |
|---|---|---|
| Intel ISA-L | https://github.com/intel/isa-l at `e8cc5e87fc64b4da434f32bc1fa18184622a4998` (reported version 2.32.1) | BSD-3-Clause; the full notice is retained beside the adapter |
| NASM | https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz | archive SHA-256 `1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148` |

When system NASM is absent, `bootstrap` downloads that exact official archive,
checks the pinned digest before extraction, builds it under the ignored
`.research/leopard2/toolchains` cache, then verifies `nasm -v`. ISA-L must be a
clean checkout at the exact commit. Its static library, license, NASM binary,
and both benchmark executables are hashed into ignored provenance and the
retained checkpoint. No third-party source, archive, or library is committed.
Bootstrap refuses any tracked or untracked repository change. It records the
clean Leopard commit and a content-addressed bundle covering the production and
standalone CMake files, both benchmark sources, runner, audit, and notice.
Checkpoint and correctness validation require those local source hashes to
still match; executable hashes alone are not accepted as reproducibility
evidence.

## Codec and timing semantics

The adapter intentionally calls `gf_gen_cauchy1_matrix`, not
`gf_gen_rs_matrix`. ISA-L's own documentation warns that the latter does not
guarantee invertibility of every submatrix for many larger parameter pairs.
Encoding setup includes construction of the systematic Cauchy generator matrix
and `ec_init_tables`. Decode-plan setup deterministically chooses surviving
systematic rows followed by the lowest parity rows until exactly `K` rows are
available, inverts that `K` by `K` matrix, and prepares only the requested
missing-original rows. Execution reuses those tables for all stripes and reuse
iterations.

Both providers operate directly on 64-byte-aligned application shard buffers;
there is no format conversion, padding, transfer, or staging copy. Encoding
reads `K * shard_bytes` and generates `R * shard_bytes` per stripe. Decode is
offered all `(K - loss_count + R)` available shards, selects exactly `K`, and
generates `loss_count * shard_bytes`. Reports retain offered and selected byte
counts separately. Setup and byte-heavy execution are always separate, and
amortized decode is derived at the declared reuse count.

Each authoritative checkpoint cell runs four independent provider pairs in
ABBA order. Each child performs two warmups and nine timed samples. The runner
requires one explicit allowed CPU and one explicitly reserved SMT sibling. It
sets its own affinity to the singleton CPU, verifies a child inherits exactly
that set, runs every provider child by inheritance, restores its original
affinity, sets OpenMP to one thread, and records topology and governor data when
readable. ISA-L retains every raw timing sample;
validation recomputes median, MAD, extrema, rates, setup amortization, cell
cardinality, ABBA order, identities, and aggregate values. Leopard2 results use
the already validated production benchmark schema. Both implementations
restore and compare every missing byte to independently retained deterministic
source data, and the paired missing-index lists must be identical.

The bounded cells are high `(240,16)`, balanced `(128,128)`, low `(64,192)`,
tiny-shard batch `(240,16)`, and padding-boundary `(129,100)` and `(225,30)`.
For the last two, public `K+R` is at most 256 and ISA-L remains GF(256), while
the requested Leopard legacy-high V1 dyadic parent is 512 and therefore uses
GF16. Results label that field advantage explicitly; it must not be presented
as an ISA-only kernel advantage.

## Reproduction

From a Leopard topic-branch checkout:

    python3 tools/leopard2_isal_compare.py self-test
    python3 tools/leopard2_isal_compare.py bootstrap --jobs 8
    python3 tools/leopard2_isal_compare.py correctness --cases 128 \
        --output experiments/leopard2/isal_compare/correctness_result.json
    python3 tools/leopard2_isal_compare.py validate-correctness \
        experiments/leopard2/isal_compare/correctness_result.json
    python3 tools/leopard2_isal_compare.py run \
        --cpu <isolated-allowed-cpu> \
        --reserved-idle-cpu <reserved-SMT-sibling> \
        --output experiments/leopard2/isal_compare/checkpoint_result.json
    python3 tools/leopard2_isal_compare.py validate \
        experiments/leopard2/isal_compare/checkpoint_result.json
    python3 tools/leopard2_external_comparison.py isa-l-checkpoint \
        experiments/leopard2/isal_compare/checkpoint_result.json

`bootstrap` is capped at eight build jobs. The single-core timing phase uses
fewer cores intentionally because cache-sensitive provider comparisons are not
valid under concurrent memory-intensive load.

## Scope limits

This adapter currently freezes only one-thread execution and positive-loss
decode. The machine-readable required-matrix audit therefore marks supported
single-thread rows `adapter-available-unmeasured`, keeps multicore rows
`adapter-required`, and excludes public lengths beyond ISA-L's GF(256)
bound. The full signed 7,134-cell matrix, persistent-pool multicore scheduling,
64/128-core batch scaling, other x86 microarchitectures, and available PMU and
memory-bandwidth evidence remain release-level work. The retained artifact also
lists these gates so a bounded checkpoint cannot be mistaken for their
completion.

ISA-L's public `ec_encode_data` call takes a byte length, source/table counts,
tables, and buffers; it has no thread-count or executor argument. The current
adapter likewise contains no persistent pool for distributing independent
stripes. Consequently the 21 required-matrix rows above one thread are real
implementation gaps and remain `adapter-required` and unmeasured. They are not
excluded from mathematical comparison, and no single-thread result is
substituted for them.
