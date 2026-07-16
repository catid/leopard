# Leopard2 / Intel ISA-L comparison protocol

> **Evidence status:** the checked-in `checkpoint_result.json` and
> `correctness_result.json` use the superseded V1 schemas. They predate the
> source-replay, raw-Leopard-sample, immutable-oracle, and post-run integrity
> gates described below. The V2 validators intentionally reject them. Their
> numbers are retained only as a provisional historical observation until a
> clean V2 rebuild and coordinated pinned run replaces both files; they are not
> accepted performance or release evidence.

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
Bootstrap refuses every source/build-input or unrelated tracked/untracked
change (the two generated evidence paths are the only exception), deletes the
ISA-L build/install tree plus both benchmark build trees, and materializes the
recorded Leopard commit as a new detached Git worktree under the ignored cache.
Both benchmarks are configured and compiled from that detached materialization,
not from the invoking checkout. Before and after compilation the runner
reconstructs and compares the complete Git-object input bundle; either a dirty
materialization or a changed object/hash aborts bootstrap. It records the clean
Leopard commit and a content-addressed bundle covering the
production codec/backends, production and standalone CMake files, both
benchmark sources, runner, audit, and notice. Every bundle entry carries its
Git mode, blob ID (or gitlink commit), and SHA-256. Validation reconstructs the
bundle with `git ls-tree` and `git cat-file` from the recorded commit and checks
its recorded tree, so replay remains valid after later documentation commits
without trusting the current checkout contents. That portable validation mode
checks artifact shape, self-derived hashes, and committed source bytes; by
itself it does not independently attest the recorded compiler, actual link
commands, or runtime-library paths. The optional
`--require-local-build-match` gate loads the ignored bootstrap cache, rehashes
and reconstructs its complete provenance, exact-compares the artifact's tools,
normalized compile and actual link commands, static inputs, runtime dependency
paths/hashes, and executable hashes, and also requires the current checkout to
be the clean benchmark-source commit. Executable hashes alone are not accepted
as reproducibility evidence.

The build identity retains normalized `compile_commands.json` entries and the
actual CMake-generated link commands for the ISA-L archive, standalone adapter,
Leopard library, and Leopard benchmark. It requires the adapter link to name the
private `${ISA_L_INSTALL}/lib/libisal.a` exactly once and binds that archive's
SHA-256 to the provider identity. `ldd` plus `readelf` are themselves identified
by resolved path, version output, and executable hash; every `DT_NEEDED`
dependency records its resolved path, real path, SONAME, and file hash. A
post-run rebuild of this provenance must match byte-for-byte.

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

Correctness does not trust that ISA-L-generated matrix as its parity oracle.
The adapter independently derives every parity coefficient as
`inverse((K + parity_index) XOR source_index)` in scalar GF(256), using the
`0x11d` polynomial and a separately implemented Fermat inverse, and checks all
`K*R` generated coefficients. The correctness campaign recomputes every parity
byte from those independently derived coefficients and immutable source bytes.
Performance children instead check a deterministic projection of at most 64
boundary/random byte positions per parity shard and stripe before and after
timing; their JSON records projection mode and exact checked/total byte counts.
Projected performance checks are never presented as the full correctness gate.

Both providers operate directly on 64-byte-aligned application shard buffers;
there is no format conversion, padding, transfer, or staging copy. Encoding
reads `K * shard_bytes` and generates `R * shard_bytes` per stripe. Decode is
offered all `(K - loss_count + R)` available shards, selects exactly `K`, and
generates `loss_count * shard_bytes`. Reports retain offered and selected byte
counts separately. Setup and byte-heavy execution are always separate, and
amortized decode is derived at the declared reuse count.

For this external comparison the Leopard child is always invoked with
`--skip-legacy --retain-samples`. The first flag is recorded in its JSON and
prevents all correctness, warmup, allocation, and timed work in the old Leopard
benchmark path; ISA-L has no corresponding second codec. Normal invocations of
`bench_leopard2` do not enable this mode and retain the prior default JSON
structure and legacy-oracle behavior. CTest includes a regression for both
shapes. Validation also recomputes the missing-index permutation directly from
`K`, loss count, and seed; verifies the high-profile padded side is
`ceil_pow2(R)` and the low-profile side is `ceil_pow2(K)`; accepts only a real
resolved backend; and requires that backend to remain identical across every
repetition.

Each V2 checkpoint cell runs four independent provider pairs in
ABBA order. Each child performs exactly two untimed warmups and exactly nine
timed samples, all of which must be retained; other counts are rejected. The runner
requires one explicit allowed CPU and one explicitly reserved SMT sibling. It
sets its own affinity to the singleton CPU, verifies a child inherits exactly
that set, runs every provider child by inheritance, restores its original
affinity, sets OpenMP to one thread, and records topology and governor data when
readable. Evidence children inherit no ambient environment. They receive only
an exact recorded allowlist (`LANG=C`, `LC_ALL=C`, one-thread OpenMP settings,
`MALLOC_ARENA_MAX=1`, and the pinned GNU affinity). Collection rejects ambient
`LD_PRELOAD`, `LD_LIBRARY_PATH`, `GLIBC_TUNABLES`, loader controls, OpenMP
controls, allocator controls, sanitizer settings, and common math-runtime
tuning variables before launching anything. It holds advisory `fcntl` locks for both logical siblings throughout
measurement and the post-timing source/tool/library/executable recheck. Those
locks serialize cooperating Leopard2 lab jobs; they are explicitly not an
OS-exclusive CPU reservation. The `fcntl` import is conditional: policy
self-tests remain portable, while pinned collection explicitly requires Linux.
ISA-L and Leopard2 retain every raw timing sample;
validation recomputes median, MAD, extrema, rates, setup amortization, cell
cardinality, ABBA order, identities, and aggregate values. Leopard2 results use
the already validated production benchmark schema. Both implementations
restore and compare every missing byte to independently retained deterministic
source data, and the paired missing-index lists must be identical.

A performance run requires `--correctness-artifact`. Its checkpoint embeds the
correctness artifact's canonical SHA-256 and exact source-bundle, build,
library, NASM, and ISA-L executable identities. Checkpoint validation likewise
requires that same correctness artifact and rejects either an identity mismatch
or a missing gate; a projected timing oracle cannot stand in for the full-byte
campaign. Authoritative checkpoint collection and validation require exactly
the documented deterministic 128-case campaign. Smaller 16-to-512-case
artifacts remain useful for standalone correctness development, but cannot
authorize or validate timing evidence; only the synthetic internal self-test
has a reduced-fixture override.

Host metadata records the scaling driver, governor, energy-performance
preference, cpuinfo and policy min/max frequencies in kHz, AMD P-state status,
and explicit nullable `boost` and Intel `no_turbo` controls. Labeled pre/post
snapshots retain both readable current-frequency sources in kHz. The post
snapshot is taken immediately when the final child returns, before parsing that
child's result or rehashing provenance. A current
frequency is a point-in-time observation, not a claim that frequency remained
fixed during a run; validation checks field shape, units, positivity, and range
ordering without requiring the two snapshots to match. On AMD P-state hosts,
`amd-pstate-epp`, EPP, and `/sys/devices/system/cpu/amd_pstate/status` therefore
remain visible even when the generic boost path is absent.

The bounded cells are high `(240,16)`, balanced `(128,128)`, low `(64,192)`,
tiny-shard batch `(240,16)`, and padding-boundary `(129,100)` and `(225,30)`.
For the last two, public `K+R` is at most 256 and ISA-L remains GF(256), while
the requested Leopard legacy-high V1 dyadic parent is 512 and therefore uses
GF16. Results label that field advantage explicitly; it must not be presented
as an ISA-only kernel advantage.

## Provisional V1 result (not accepted evidence)

The superseded V1 checkpoint was measured on an AMD Ryzen 9 9950X3D under Linux
6.8.0-134. The runner and every child had singleton affinity to CPU 15; its SMT
sibling CPU 31 was explicitly reserved idle. The readable scaling governor was
`powersave`. The build used clean Leopard commit
`1f41ddd41e2b3c040ff06e15f3dbaa58a8b05863` and source-bundle SHA-256
`b284300601d18ef6edfb6e786067611bb3d59764c64de3b5a76ef660d283ec5c`.

The following provisional rates are decimal GB/s. `Enc out` counts generated parity bytes.
The three execution-only decode columns count, respectively, every offered
received byte, the deterministic `K`-row subset actually consumed, and repaired
original bytes. `Plan us` is not included in those execution rates. `Amort out`
includes plan setup amortized over the cell's declared reuse count of eight.
Each number is the median of four independent run medians; each run contains
nine timed samples after two warmups. The checkpoint retains the corresponding
MADs and all raw ISA-L samples.

| Cell | Provider | Enc out | Dec offered | Dec selected | Dec repaired | Plan us | Amort out |
|---|---:|---:|---:|---:|---:|---:|---:|
| high 240/16, 64 KiB x1, L=1 | ISA-L | 0.759 | 38.511 | 36.246 | 0.151 | 8225.907 | 0.045 |
| high 240/16, 64 KiB x1, L=1 | Leopard2 | 0.895 | 8.512 | 8.012 | 0.033 | 1.231 | 0.033 |
| high 240/16, 4 KiB x8, L=4 | ISA-L | 0.869 | 33.428 | 31.836 | 0.531 | 5831.250 | 0.134 |
| high 240/16, 4 KiB x8, L=4 | Leopard2 | 1.054 | 9.708 | 9.245 | 0.154 | 1.270 | 0.154 |
| balanced 128/128, 64 KiB x1, L=8 | ISA-L | 1.716 | 37.143 | 19.171 | 1.198 | 835.975 | 0.966 |
| balanced 128/128, 64 KiB x1, L=8 | Leopard2 | 6.542 | 5.007 | 2.584 | 0.162 | 1.290 | 0.162 |
| low 64/192, 64 KiB x1, L=8 | ISA-L | 3.607 | 75.422 | 19.464 | 2.433 | 129.931 | 2.262 |
| low 64/192, 64 KiB x1, L=8 | Leopard2 | 9.108 | 12.062 | 3.113 | 0.389 | 1.270 | 0.389 |
| boundary 129/100, 64 KiB x1, L=4 | ISA-L GF8 | 1.709 | 62.860 | 36.040 | 1.118 | 856.779 | 0.769 |
| boundary 129/100, 64 KiB x1, L=4 | Leopard2 GF16 | 3.249 | 2.752 | 1.578 | 0.049 | 21.446 | 0.049 |
| boundary 225/30, 64 KiB x1, L=2 | ISA-L GF8 | 0.971 | 40.322 | 35.859 | 0.319 | 4884.030 | 0.128 |
| boundary 225/30, 64 KiB x1, L=2 | Leopard2 GF16 | 1.206 | 4.864 | 4.325 | 0.038 | 7.130 | 0.038 |

The provisional V1 data reported that Leopard2 generated parity faster in all
six cells: ISA-L's encode throughput was
0.262x to 0.849x Leopard2. ISA-L's byte-heavy decode execution was 3.443x to
22.839x faster by repaired-output rate. ISA-L's much larger matrix-inversion and
table setup cost matters at low reuse: after amortization it was 0.872x
Leopard2 in the 4 KiB batch cell, while remaining 1.345x to 15.724x in the
other five cells. The 22.839x and 15.724x boundary figures also include GF8
versus padded-parent GF16 and are not kernel-only comparisons. These six cells
are a bounded checkpoint, not evidence for the unmeasured matrix or other
machines.

The provisional V1 deterministic adapter campaign reported 128/128 cases, including
eight no-loss cases, sixteen maximum-loss cases, and 34 cases where dyadic
padding gives ISA-L GF8 versus Leopard2 GF16. Every case verifies the systematic
generator prefix and restores each requested source byte against independently
retained input. It does not compare parity bytes because the wire formats and
generator matrices intentionally differ.

The superseded V1 artifacts use a canonical digest that excludes their self-identifying
`artifact_sha256` value. The fail-closed validators recompute it. For independent
transport checks, the complete serialized-file digest is also listed:

| Provisional V1 artifact | Historical canonical SHA-256 | Serialized file SHA-256 |
|---|---|---|
| `checkpoint_result.json` | `eb3b8083425321a922ded8eaca9ae1144658b5615f596219f9245b07d5c32834` | `e67efd48af1d7333144d958adfc0b850c08d2048b87974eec41bb33147b220a8` |
| `correctness_result.json` | `3349b66454a15c8121170f743ad759c843bd157654290a591fb88f9497055ee4` | `4392d9d5320434ab8a944194911fac0d6e01a75ab5d8f36d84d202b5f8d723d4` |

## Reproduction

From a Leopard topic-branch checkout:

    python3 tools/leopard2_isal_compare.py self-test
    python3 tools/leopard2_isal_compare.py bootstrap --jobs 8
    python3 tools/leopard2_isal_compare.py correctness --cases 128 \
        --output experiments/leopard2/isal_compare/correctness_result.json
    python3 tools/leopard2_isal_compare.py validate-correctness \
        experiments/leopard2/isal_compare/correctness_result.json \
        --cache .research/leopard2 --require-local-build-match
    python3 tools/leopard2_isal_compare.py run \
        --cpu <isolated-allowed-cpu> \
        --reserved-idle-cpu <reserved-SMT-sibling> \
        --correctness-artifact \
            experiments/leopard2/isal_compare/correctness_result.json \
        --output experiments/leopard2/isal_compare/checkpoint_result.json
    python3 tools/leopard2_isal_compare.py validate \
        experiments/leopard2/isal_compare/checkpoint_result.json \
        --correctness-artifact \
            experiments/leopard2/isal_compare/correctness_result.json \
        --cache .research/leopard2 --require-local-build-match
    python3 tools/leopard2_external_comparison.py isa-l-checkpoint \
        experiments/leopard2/isal_compare/checkpoint_result.json \
        --correctness-artifact \
            experiments/leopard2/isal_compare/correctness_result.json \
        --cache .research/leopard2

`bootstrap` is capped at eight build jobs. The pinned `run` command is Linux-only;
the parser/self-test and artifact policy remain usable elsewhere. The
single-core timing phase uses
fewer cores intentionally because cache-sensitive provider comparisons are not
valid under concurrent memory-intensive load.

Without `--require-local-build-match`, the two validation commands replay the
recorded Git source commit and verify internal artifact consistency, but they
do not treat recorded tool/link/runtime labels as independently trusted. The
strict flag is the local evidence-audit mode: it requires the bootstrap cache,
reconstructs that cache's live build provenance, exact-compares the complete
build identity, and requires current `HEAD` and its build inputs to be clean and
equal to the benchmark-source commit. The external `isa-l-checkpoint` audit
also requires `--cache`, reconstructs the same trusted provenance, enforces the
exact 128-case gate, and reports a trusted-cache verification status; it cannot
upgrade a portable-only or self-consistently relabeled artifact to verified.
Do not run the timing command
concurrently with any other cache- or memory-sensitive job.

## Scope limits

This adapter currently freezes only one-thread execution. Positive-loss and
no-loss decode are supported. The machine-readable required-matrix audit
therefore marks supported single-thread rows `adapter-available-unmeasured`,
keeps multicore rows `adapter-required`, and excludes public lengths beyond
ISA-L's GF(256) bound. The full signed 7,134-cell matrix, persistent-pool
multicore scheduling, 64/128-core batch scaling, other x86 microarchitectures,
and available PMU and memory-bandwidth evidence remain release-level work. The
retained artifact also lists these gates so a bounded checkpoint cannot be
mistaken for their completion.

ISA-L's public `ec_encode_data` call takes a byte length, source/table counts,
tables, and buffers; it has no thread-count or executor argument. The current
adapter likewise contains no persistent pool for distributing independent
stripes. Consequently the 21 required-matrix rows above one thread are real
implementation gaps and remain `adapter-required` and unmeasured. They are not
excluded from mathematical comparison, and no single-thread result is
substituted for them.
