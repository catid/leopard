# 40 — Sparse-high production-AUTO qualification

## Result

The bounded legacy-high GF8 sparse-Q1 production-AUTO qualification did **not**
meet its preregistered performance gate. All 88 cells passed acquisition and
integrity validation, but only 56 of 74 candidates cleared both five-percent
lower-bound requirements, and one of 14 route-negative controls crossed the
prespecified minus-two-percent net point-regression guard. The final strict
resume acquisition and both deterministic replays therefore returned exit
status 2; the initial pass returned 1 solely for its 11 archived isolation
rejections.

This is a valid no-promotion performance decision, not a benchmark acquisition
or correctness failure. Production AUTO was measured, but the compiled
sparse-high AUTO default remains `OFF`; no production route, public API, or
installed behavior changed.

The failures are not isolated statistical edge cases. The 36 qualified-tuple
one-shot cells split 24 passing and 12 failing, while the 30 public batch and
binding cells split 24 passing and six failing. All six public-API failures are
the `K=16,R=2,4096` anchor, one in each batch/binding lane at batches 1, 4, and
16. All eight parity-row samples passed. These cells remain separate in the
checkpoint and are not averaged into a favorable aggregate.

## Accepted campaign

The accepted source was detached at commit
`b455108e61017c711e238ccf159d45da50e77ca2`, tree
`4c2caeb67fa22379f9ac9de74e7a4039207fae19`, with initialized `sse2neon`
gitlink `cad518a93b326f0f644b7972d488d04eaa2b0475` and source fingerprint
`ca8f40ff251c035ac14b222225ac15e12828b827189b092cad77ad607b2a3d01`.
The worktree was clean and its before/after fingerprints matched. The runner
SHA-256 was
`28edc37f2c51de483e1d602a8066fee1fee9243ff140a5412cc017cc7fa024b0`.

The runner created one tests-off Release/Ninja explicit-AVX2 build with GF8 and
OpenMP enabled, GF16 and full-output high direct disabled, sparse table
preparation enabled, and sparse AUTO compiled-default off. It rejected test-hook
objects and froze the ordinary production executable and archive together:

- executable: `8432e5a2cc58f31df340130d329d3351cc57a565a443cff56f0bd02ae94ceb39`;
- ordinary `libleopard.a`:
  `1b43f5ee4d07e0ca64fac5ebcb25aee8a67e9d684a4e21eb71fd174911d93067`;
- frozen-bundle provenance:
  `6847a7e494321c28aba152c984a212e5f1857b15eed50329404ee6f812d8e430`.

Timing ran on `foureyes.lan`, an AMD Ryzen Threadripper PRO 9985WX, with every
benchmark process pinned to CPU 109 and its exact SMT sibling CPU 45 reserved.
Both are core 45 on socket 0. A pinned 20-second preflight observed zero
non-idle jiffies on both CPUs. The runner held the canonical campaign lock and
an exclusive pair-wide lease, used one worker, and accepted only jobs whose
reserved sibling accumulated zero non-idle jiffies.

The fixed 88-cell grid has SHA-256
`ff6ab52d98a915afd1a86003cd6820dae1a63c66faa3276a8d97b056155db6d2`.
Each cell used three rounds of two ABBA contrasts, 15 timing iterations after
four warmups, 15 setup iterations, four API calls per sample, modeled reuse 64,
and a 120-second per-process timeout. That produced 24 serialized processes per
cell and 2,112 raw benchmark records overall. The primary metric was amortized
derived median microseconds per API call.

## Correctness and replay gates

A fresh clean clone of the exact source passed the focused GCC 13 Release gate
on `ripper.lan`, 5/5 tests: the production sparse-AUTO and sparse-encode tests,
benchmark registration, and normal plus `python -O` smoke parsers. It used
backend AUTO, OpenMP, both fields, sparse tables on, and sparse AUTO
compiled-default off. Its log, CMake cache, ordinary archive, and benchmark
SHA-256 values were respectively:

- `d357eeabf08c19c98bb3934894de8fd51b4483099fcf5db9dc00ed208d9336f8`;
- `d54ca2b45250b8ad79257b51628ecb1e36d190a933352c63e03dbc573c1dad74`;
- `fe8a9fc3adc5b2650807a596e897090385982808e2a83dcd03403e4be84d9b44`;
- `02b38273a5a7aac4b1832314a729763f031e8180c36b16e9014cad8698d6be38`.

The fresh GCC 13 ASan+UBSan lane passed 4/4 instrumented tests: both production
tests and both smoke parsers. Leak detection, strict string checks, ASan/UBSan
halt-on-error, UBSan stack traces, and one OpenMP thread were enabled. The
recursive registration test intentionally remained in the Release lane because
it creates an independent unsanitized Release build. The sanitizer log, CMake
cache, ordinary archive, and benchmark SHA-256 values were respectively:

- `80ffdd26075ea3db0cfb0348ea4e506296fe7ef155248605effb37ced57d7a6d`;
- `89de5c8af934f0f78fdb9a919dd1c624b13173b74fff23e2b42851d42356dd79`;
- `6805f34d6a3a8efbe1f9e3589741a55895ffbb88c43363a320951c5742ea948c`;
- `fa13a305689a32ba5dc715b39438a616f1dfb76dde9634b82f664fecb49f68cc`.

The sanitizer configure command supplied global C/C++ flags
`-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer`, while CMake's
`RelWithDebInfo` flags appended `-O2 -g -DNDEBUG`; the retained `build.ninja`
therefore used effective `-O2` as its last optimization option. This was a
correctness sanitizer gate and had no effective-O1 requirement.

The full unrelated CTest suite was not substituted for these focused gates;
clean-clone failures outside this lane remain tracked in Bead `leopard-28m`.

The focused-gate reproduction commands reconstructed from the retained CMake
caches and CTest logs are below. `SOURCE`, `BUILD`, and the canonical lock
holder must refer to fresh, lane-owned locations; the coordinator serialized
builds and tests under `/tmp/leopard-gf8-authoritative.lock`.

```sh
/usr/bin/cmake -S "$SOURCE" -B "$BUILD" -G Ninja \
    -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 \
    -DCMAKE_BUILD_TYPE=Release \
    -DLEO2_BACKEND_VARIANT=auto \
    -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON \
    -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
    -DENABLE_OPENMP=ON -DLEOPARD_ENABLE_GF8=ON \
    -DLEOPARD_ENABLE_GF16=ON \
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF \
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE=ON \
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE_AUTO=OFF
/usr/bin/cmake --build "$BUILD" --parallel 128 --target \
    leopard2_high_sparse_auto_production_test \
    leopard2_sparse_encode_production_test \
    bench_leopard2_high_sparse_auto
/usr/bin/ctest --test-dir "$BUILD" --output-on-failure -R \
    '^leopard2_(high_sparse_auto_(production|benchmark_(registration|smoke|smoke_optimized))|sparse_encode_production)$'
```

The sanitizer lane uses the same options and targets except for these configure
and test differences:

```sh
/usr/bin/cmake -S "$SOURCE" -B "$BUILD" -G Ninja \
    -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    '-DCMAKE_C_FLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer' \
    '-DCMAKE_CXX_FLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer' \
    -DLEO2_BACKEND_VARIANT=auto \
    -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON \
    -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
    -DENABLE_OPENMP=ON -DLEOPARD_ENABLE_GF8=ON \
    -DLEOPARD_ENABLE_GF16=ON \
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF \
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE=ON \
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE_AUTO=OFF
/usr/bin/cmake --build "$BUILD" --parallel 128 --target \
    leopard2_high_sparse_auto_production_test \
    leopard2_sparse_encode_production_test \
    bench_leopard2_high_sparse_auto
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 OMP_NUM_THREADS=1 \
    /usr/bin/ctest --test-dir "$BUILD" --output-on-failure -R \
    '^leopard2_(high_sparse_auto_(production|benchmark_(smoke|smoke_optimized))|sparse_encode_production)$'
```

The selected correctness files now live under owner-only durable storage at
`ripper.lan:/home/catid/leopard-evidence/production-auto-b455108-20260831-correctness`.
Its manifest covers ten files totaling 125,326,031 bytes and has SHA-256
`b2598c741a8acd1046dcca1c10999183f3b198c14910a8eb0359b857213d34de`.
That manifest binds both exit markers, logs, caches, archives, and benchmark
binaries. The original shell argv, sanitizer environment block, lock
acquisition, and before/after Git commands were not serialized in the focused
gate logs; the command, environment, and lock details above are coordinator
records and reproducibility metadata, not acquisition-time cryptographic
attestations.

Every accepted campaign job passed the independent systematic-generator and
post-timing parity oracles, route/table witnesses, stable non-vacuous input and
parity identities, input-immutability and output-canary checks, process
containment, raw/log hashing, frozen-artifact validation, source stability, and
isolation. The exact inventory contains 88 job documents, 2,112 raw JSON
documents, 2,112 stdout logs, and 2,112 stderr logs, with no extras.

Normal and assertion-disabled replay wrote outside the evidence tree:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
    tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir /tmp/leopard-production-auto-qualification-b455108-foureyes-v2 \
    --promotion-percent 5 \
    --output /tmp/leopard-production-auto-qualification-b455108-foureyes-v2-validation/analysis-normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O \
    tools/leopard2_direct_encode_crossover.py analyze \
    --result-dir /tmp/leopard-production-auto-qualification-b455108-foureyes-v2 \
    --promotion-percent 5 \
    --output /tmp/leopard-production-auto-qualification-b455108-foureyes-v2-validation/analysis-optimized.json
```

Both returned 2. The retained, normal, and optimized analysis files were
byte-identical with SHA-256
`a6da9452f530e3ac3af183a268822094c0715da84b06b19a01ad26f8eff6a853`,
and the retained document canonically matched `matrix.json.analysis`. A full
pre/post inventory proved that replay did not change the result tree.

## Candidate and route-negative control results

A candidate qualified only when both the production-route and summed-net
marginal 95-percent lower bounds reached five percent. The control guard used
the net point estimate, not its interval. Representative decisive cells are:

| Cell | Route gain (95% interval) | Table gain (95% interval) | Net gain (95% interval) | Decision |
| --- | ---: | ---: | ---: | --- |
| one-shot `K=3,R=2`, 4096 B | -27.06% (-31.88% to -21.90%) | 2.65% (-9.39% to 16.30%) | -25.13% (-36.02% to -12.39%) | weakest candidate route lower bound; fail |
| binding batch 1 `K=16,R=2`, 4096 B | -23.29% (-27.83% to -18.48%) | -1.62% (-13.15% to 11.44%) | -24.54% (-37.25% to -9.25%) | weakest candidate net lower bound; fail |
| one-shot `K=16,R=4`, 4096 B | 21.74% (14.73% to 29.19%) | -0.29% (-2.77% to 2.27%) | 21.40% (12.34% to 31.18%) | weakest passing candidate lower bound |
| explicit-AVX2 binding batch 16 `K=2,R=16`, 4096 B | -0.38% (-22.40% to 27.88%) | -1.82% (-8.33% to 5.15%) | -2.19% (-24.74% to 27.10%) | control point guard triggered |

The forced-path discovery in report 39 measured encode execution after forcing
either arithmetic path. This campaign instead measured the ordinary production
AUTO route through public APIs and separately accounted for table preparation.
The negative route contrasts above show that the forced-path win does not carry
through this bounded production route in every preregistered cell. Passing
neighbors do not authorize deleting the failing cells or narrowing the policy
after seeing these data.

## Scope fences

The conclusion applies only to the frozen 74 candidates and 14 controls on this
host and OpenMP runtime. Unmeasured parity rows, K/R/byte/API/batch crosses,
reuse counts, backends, thread counts, out-of-table neighbors, hardware, and
OpenMP runtimes remain outside the evidence. The intervals are per-cell
marginal two-sided 95-percent Student-t intervals with two degrees of freedom,
not a simultaneous campaign guarantee. This negative decision supplies no
authority for a narrower post-hoc selector; any redesigned predicate or
implementation requires a fresh preregistered campaign.

## Isolation retries

No contaminated measurement entered inference. A first launch failed before
build or timing because the reused clone's `sse2neon` working directory was
uninitialized and did not resolve to its indexed gitlink. Its empty result root,
log, and exit marker are retained as a pre-evidence rejection and were never
pooled.

The first complete CPU109 pass accepted 77 jobs and rejected 11 solely because
reserved sibling CPU45 accumulated 1–18 non-idle jiffies. The entire pre-resume
state was archived before replacement. Its durable copy is
`foureyes.lan:/home/catid/leopard-evidence/production-auto-b455108-20260831/discarded-attempts/attempt1-result-tree.tar.zst`,
size 4,398,717 bytes, SHA-256
`d6993f4267026935b2ee796608c0106b8e995beeaa54832fdfb2fb697897ddba`.
Its failure projection SHA-256 is
`1f1ca506c7e04130b1e07b4d7cf20aca9b5ed23dec2af03a55861311c0e641af`.

After a pinned five-second preflight again observed zero non-idle jiffies on
both CPUs, an observer-free strict resume reran all 24 invocations for only
those 11 jobs. All passed with zero sibling work, yielding an accepted
distribution of 77 initial jobs plus 11 complete resume jobs. The 77 passed job
documents, 1,848 raw records, and 3,696 logs from the initial pass are
byte-identical in the final corpus and were revalidated for inference. Only the
11 failed job attempts were discarded; their 264 raw records were replaced by
complete reruns and were neither used nor pooled. The pre-resume archive and
external acquisition logs provide attribution; the runner's retained `resumed`
field remains false and is not used as proof.

## Evidence lineage and disposition

The accepted top-level SHA-256 values are:

- manifest: `db5d472d73ea0646e782665e6960dc9c7c0c736b9987ef4280a80e3a21789d82`;
- matrix: `346326eff8f0bdc4b402bf641f0a49463cd7d998b272c004e1736317d7306433`;
- analysis: `a6da9452f530e3ac3af183a268822094c0715da84b06b19a01ad26f8eff6a853`;
- controlled build:
  `7eab103780f3dd29262f8e4d8a75d992b17a1561312c92cbb59fd5c8a7926916`.

The accepted raw working tree remains at
`foureyes.lan:/tmp/leopard-production-auto-qualification-b455108-foureyes-v2`.
Its 6,490 regular files total 174,689,363 bytes; the SHA-256 of the canonical
sorted relative-path/mode/size/file-hash inventory is
`003cbb457cfa057957b98794ea185595d9c1d4a0fbdad7b41363a50fc024b4cd`.
A compressed accepted archive is retained durably at
`foureyes.lan:/home/catid/leopard-evidence/production-auto-b455108-20260831/accepted/result-tree.tar.zst`.
It is 4,401,675 bytes with SHA-256
`1aedb13aa49017a73e16c2bd847ffb3a21a137dbe2fd3a9c7360dddb36c65e54`.
The same owner-only durable root contains both acquisition attempts, their exit
markers, the discarded-attempt projection, both external replay outputs, and
the full pre/post file inventories. Its manifest covers 20 files totaling
11,601,776 bytes and has SHA-256
`5e8adb4eba3fe99c71091af4c11afd3793c8ee001badc7450f57c902faadc2e8`.

The committed projection, including all 88 cell estimates plus the pre-evidence
rejection and isolation-retry lineages, is
`experiments/leopard2/direct_encode/results/production_auto_avx2_checkpoint_20260831.json`.
The raw tree restored from the durable archive is still required for full
replay. The final disposition is no promotion and no production change;
sparse-high AUTO remains default-off.
