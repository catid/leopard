# Leopard2 balanced full-recovery dispatch

Status: production crossover enabled only in the measured GF8 region described
below. This choice changes execution strategy, not field mathematics,
coordinates, parity bytes, or the wire profile.

## Decision

For the legacy-high profile, the retained generic active-parent decoder is used
automatically only when all of the following are true:

- `K = R = T = 128` and `N = 256`;
- the field is legacy GF8;
- every original is missing;
- the rounded shard size is 256 bytes through 1 MiB; and
- the selected backend is scalar, SSSE3, or AVX2.

All neighboring counts, partial-loss patterns, GF16 codes, NEON, and sizes
outside that interval retain the profile-specific high decoder. Explicit
`LEO2_CODEC_FORCE_GENERIC_DECODE` and
`LEO2_CODEC_FORCE_SPECIALIZED_DECODE` flags override automatic selection for
differential tests and benchmarks. They are mutually exclusive and do not alter
the code identity.

The reason is aggregate work, not the absence of a smaller transform side.
The operation model counts 1,280 generic versus 1,856 specialized butterflies
at this point, with estimated upper workspace traffic of 46.5 versus 51.5
passes. Algorithm 5 still saves one transform level because `T = N/2`, but full
original recovery removes its message-only output pruning.

## Historical isolated evidence

The retained summary is
`experiments/leopard2/decoder_dispatch/results/balanced_amd_9950x3d.json`; its
file SHA-256 is:

    76abe8d88f6c2e8d7220b6b6f784e8a0f3f2f502953491a13f280efa36bf5c0c

It contains decode-execution and decode-plan-setup medians and MADs separately,
SHA-256 identities of the raw benchmark JSON inputs, build identities, exact
run ordering and parameters, and a canonical content digest. It predates the
current authenticated three-mode protocol. The former analyzer inferred
specialized execution from `force_generic_decode=false`; that is insufficient
after AUTO can select generic and after materialized and tiled Algorithm 5
became separate paths. The retained artifact is historical provenance for this
already-narrow rule only. It cannot be replayed through the current analyzer,
relabeled as current evidence, or used to widen dispatch.

Three separate Release binaries were built from pre-dispatch commit `50e7858`,
so `force_generic=false` unambiguously selected the specialized path. Their
SHA-256 identities are retained in the result. Each cell ran on CPU 0 alone,
with `OMP_NUM_THREADS=1`, three warmups, 15 timing samples, and four executions
per sample. A second complete run reversed both size and decoder order. The
machine was an AMD Ryzen 9 9950X3D with 32 allowed logical CPUs, one socket and
one NUMA node. CPU 0 used `amd-pstate-epp`, the `powersave` governor, and a
reported 600 MHz to 5.752 GHz range. `perf_event_paranoid=4` prevented
unprivileged hardware-counter collection.

Mean generic improvement over specialized across the two runs:

| Backend | 64 B | 256 B | 4 KiB | 64 KiB | 1 MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| scalar | 28.74% | 30.94% | 31.58% | 32.13% | 26.11% |
| SSSE3 | 6.20% | 14.17% | 16.69% | 16.73% | 5.18% |
| AVX2 | 0.38% | 8.83% | 14.31% | 19.00% | 7.75% |

Maximum MAD was 0.89% of the median. The 64-byte AVX2 row was effectively tied
and is deliberately excluded from automatic dispatch. Every promoted measured
anchor clears the 5% experiment threshold in both-run aggregate. Neighboring
counts, loss patterns, fields, and backends were not part of this isolated
artifact and therefore retain the previous specialized selection. Sizes between
the measured anchors use the same linear byte kernels; the upper and lower
bounds remain conservative and deterministic.

An earlier 45-job pipeline run and a 144-job crossover run overlapped other
agents' correctness work. They remain useful pipeline and round-trip evidence,
but their timing was discarded and is not a basis for this decision. In
particular, the original Bead text naming a 4 KiB, one-loss regression came from
that discarded run; the isolated sweep instead showed that the robust crossover
is the full-original-loss case.

## Correctness gates

After the dispatch change:

- the full Release CTest suite passed 13/13;
- the ASan+UBSan CTest suite passed 13/13;
- forced `auto`, scalar, SSSE3, and AVX2 builds passed deterministic golden,
  API, random, and transform-differential tests, with no cross-backend mismatch;
  the locally retained, policy-ignored
  `results/leopard2/backend-matrix/matrix.json` has SHA-256
  `5bfa8137f2c223b061e9eca4b19623c33129e1572ed10a6aef2b9b6ec7b52037`
  and source fingerprint
  `7061c28ac8747c48cbe08a796674e12c42d1ceee65105670440d992743efcc53`;
- the API regression recovers all 128 originals at an arbitrary 257-byte size
  and compares automatic, forced-generic, and forced-specialized outputs; and
- mutually exclusive and unknown codec flags are rejected; and
- a pure policy test covers both byte boundaries and rejects neighboring
  counts, partial recovery, GF16, the low profile, and unmeasured backends.

The forced-backend matrix was regenerated after the final evidence-scope
tightening and diagnostic flag addition; its source fingerprint remained stable
for the entire run.

Regenerate that ignored correctness artifact with:

    python3 tools/leopard2_backend_matrix.py run --jobs "$(nproc)" --no-resume

## Current reproduction contract

Fresh promotion evidence uses the same-binary external-matrix runner described
in `leopard2_balanced_forced_evidence.md`. It authenticates exact generic,
materialized, and tiled force tuples in benchmark argv and output, performs
three clustered ABBA/BAAB/ABBA rounds, binds source/object/executable/runtime
identity, and rejects an active SMT sibling. The committed smoke matrix is a
protocol check, not a replacement for the full matrix.

To reconstruct the historical input only, create a detached worktree at the
recorded pre-dispatch commit, then configure one benchmark per backend:

    git worktree add --detach /tmp/leopard2-pre-dispatch 50e7858
    for backend in scalar ssse3 avx2; do
        cmake -S /tmp/leopard2-pre-dispatch \
            -B "/tmp/leopard2-pre-dispatch/build/dispatch-$backend" \
            -G Ninja -DCMAKE_BUILD_TYPE=Release \
            -DLEO2_BACKEND_VARIANT="$backend" \
            -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON \
            -DLEO2_ENABLE_CUDA=OFF
        cmake --build "/tmp/leopard2-pre-dispatch/build/dispatch-$backend" \
            -j "$(nproc)" --target bench_leopard2
    done

For every backend and size `64 256 4096 65536 1048576`, run both the default
pre-dispatch specialized codec and `--force-generic`, using the parameters
retained under `method` in the JSON summary. Reverse mode and size order for the
second run. The historical raw format has intentionally not been grandfathered
into the new analyzer because doing so would permit the ambiguous
AUTO-as-specialized classification. Host-specific medians are expected to
vary; exact forced path, clean build closure, raw observations, and isolation
are now the reproducible contract.
