# C1 inverse-prefix encoder evaluation

## Scope and hypothesis

The high and low encoders already consume complete four-source groups through
the out-of-place inverse radix-four backend, but the mature schedule still
materializes every shortened suffix shard as zero before later layers read it.
The C1 hypothesis was that an immutable dependency schedule could preserve the
source-staged prefix, treat the suffix as structural zero, and remove that
memory pass profitably.

This experiment is parent-code preserving.  It does not change the field,
coordinate map, transform normalization, parity order, or wire profile.

## Implementation and oracle

`PrunedTransformPlan` now records the active inverse prefix and the three fixed
multipliers for each complete four-source group.  The allocation-free source
executor:

1. reads only the active source pointer array;
2. sends complete groups directly through the selected scalar, SSSE3, or AVX2
   out-of-place inverse radix-four operation;
3. copies at most three ragged active shards;
4. skips the consumed leaf operations and executes the remaining immutable C1
   schedule without reading or zeroing the inactive suffix; and
5. leaves caller sources unchanged.

The low profile owns one plan for `K < P`.  Legacy-high owns one plan for a
partial final `T`-shard message block, with shift
`T + floor(K / T) * T`.  Complete blocks retain the mature regular schedule.
The candidate is selected only by the diagnostic `FORCE_TRANSFORM` path and by
the internal benchmark.  Production `AUTO` remains on the mature encoder.

The independent oracle is the padded full inverse transform followed by the
existing direct systematic generator-matrix parity oracle.  Tests publish a
source pointer array of exactly the active-prefix length, so ASan detects even
validation of a nonexistent shortened pointer.

## Correctness evidence

Release validation on GCC 13.3.0 passed:

- `leopard2_pruned_transform_test`: 3 backends, 25,407 plans, 33,399
  executions, and 99,581,751 compared bytes;
- every strict prefix through transform size 64, plus ragged/four-way boundary
  classes through GF8 size 256 and GF16 size 1024;
- arbitrary GF8 byte tails and complete-symbol GF16 tails;
- concurrent execution of one shared immutable inverse plan;
- `leopard2_direct_encode_test`: 1,024 profiles, 8,704 basis messages, 1,024
  random messages, 433,513 parity symbols, 8,192 mask executions, aliases,
  unaligned inputs, batches, and concurrent calls;
- deterministic pruned-transform fuzz replay: 16,384 iterations, seed
  `14627333968358193854`; and
- all 74 ordinary CTest entries passed after initializing the clean worktree's
  `sse2neon` submodule; the initial missing-submodule project-test failure was
  a worktree setup issue, not a codec mismatch;
- focused Clang 18 ASan, UBSan, and leak-sanitizer runs passed the pruned-plan,
  direct-encode, 16,384-case fuzz replay, and concurrent-encode tests; and
- Clang 18 TSan passed the immutable-plan and concurrent-encode tests.  The
  direct-encode test deliberately overrides global allocation operators and
  therefore cannot link with TSan's own replacements; its concurrent calls
  remain covered by the compatible tests.

Instrumentation proves that the candidate writes each active source shard once,
copies at most three of those shards, zero-fills no shortened shard, and skips
exactly `side - active_prefix` zero-fill shard writes.  Parity remained
byte-identical in every exercised high/low, GF8/GF16, scalar/SSSE3/AVX2 cell.

## Initial instrumented screen (not quantitative evidence)

The first pinned screen linked the benchmark against the test-hook archive.
That archive performs atomic instrumentation updates in encoder loops, and the
mature and candidate inverse forms execute different numbers of those updates.
The bias favors the candidate, so the screen remains useful as a conservative
reason not to promote it, but its numerical ratios are invalid as production
throughput claims.  The table is retained to make the discarded evidence and
the bug-fix trail explicit; a production-linked rerun follows in the evidence
commit after the benchmark's fail-closed identity gate.

Measurements used one logical CPU, CPU 4, with its SMT sibling CPU 20 idle.
Both are physical core 4 on socket/NUMA node 0.  The process affinity was
`{4}`, `OMP_NUM_THREADS=1`, the reported governor was `powersave`, and the host
advertised 600--5752 MHz.  No other worker used CPU 4/20 during these cells.
The benchmark rotates mature prefix, mature exact-output, call-local exact,
pruned-inverse prefix, and pruned-inverse exact forms.  It reports medians and
MAD over 15 samples and keeps inverse-plan setup separate.  Ratio is mature
execution time divided by pruned-inverse execution time, so values above one
favor C1.

| Field/backend | Profile `(K,R)` | Shard bytes | Prefix ratio | Exact-output ratio | Result |
|---|---:|---:|---:|---:|---|
| GF8 AVX2 | low `(33,34)` | 64 KiB | 0.922 | 0.924 | 7.8% slower |
| GF8 AVX2 | low `(33,100)` | 64 KiB | 0.966 | 0.979 | 3.4% slower |
| GF8 AVX2 | low `(33,34)` | 1 MiB | 0.993 | 1.002 | neutral |
| GF8 AVX2 | low `(65,66)` | 1 MiB | 1.036 | 1.022 | 3.6% faster |
| GF16 AVX2 | low `(129,130)` | 64 KiB | 1.003 | 1.005 | neutral |
| GF8 scalar | low `(9,10)` | 64 KiB | 1.069 | 1.075 | 6.9% faster |
| GF8 scalar | low `(17,18)` | 64 KiB | 1.092 | 1.101 | 9.2% / 10.1% faster |
| GF8 scalar | low `(18,19)` | 64 KiB | 1.088 | 1.094 | 8.8% faster |
| GF8 scalar | low `(65,66)` | 64 KiB | 1.071 | 1.075 | 7.1% faster |
| GF8 scalar | low `(65,66)` | 1 MiB | 1.065 | 1.073 | 6.5% faster |
| GF16 scalar | low `(65,66)` | 64 KiB | 1.053 | 1.051 | 5.3% faster |
| GF8 AVX2 | high `(33,32)` | 64 KiB | 0.914 | 0.917 | 8.6% slower |
| GF8 AVX2 | high `(65,64)` | 1 MiB | 0.948 | 0.913 | 5.2% slower |

For low GF8 `(33,100)` the immutable inverse plan contained 143 retained
operations, eight source groups, and skipped 31 zero-shard writes; median plan
setup was 2.083 microseconds.  Setup is amortizable, but execution still lost
3.4%.  The flat boundary operations and indirect dispatch offset the removed
zero pass on AVX2.  High rate also loses the mature fused final-layer
accumulation because the current C1 candidate materializes the partial block
before XORing it into the accumulator.

## Production-linked authoritative evidence

The corrected benchmark is buildable with `LEO2_BUILD_TESTS=OFF`, reports
`library_test_hooks=false`, and is rejected by the pinned runner if that marker
is true or missing.  Commit `67eb32957824fe2b9a2a51945a6e1980aaef73a1`
was clean before and after the run.  The runner pinned one worker to CPU 4,
verified CPU 20 was its idle SMT sibling, used 15 rotated samples of 32
iterations with three warmups, and produced `authoritative=true` for all five
cells.  Candidate/control is elapsed time, so less than one favors C1;
control/candidate and the amortized columns are speed ratios, so greater than
one favors C1.

| Field/backend | Profile `(K,R)` | Bytes | Candidate/control | Control/candidate | Amortized control/candidate at reuse 1 / 8 / 64 | Setup median |
|---|---:|---:|---:|---:|---:|---:|
| GF8 scalar | low `(17,18)` | 64 KiB | 0.9158 | 1.0919 | 1.0912 / 1.0918 / 1.0919 | 0.793 us |
| GF8 AVX2 | low `(33,34)` | 64 KiB | 1.0805 | 0.9255 | 0.9195 / 0.9247 / 0.9254 | 2.003 us |
| GF8 AVX2 | low `(65,66)` | 1 MiB | 0.9665 | 1.0347 | 1.0346 / 1.0347 / 1.0347 | 4.015 us |
| GF8 AVX2 | high `(33,32)` | 64 KiB | 1.0938 | 0.9143 | 0.9115 / 0.9139 / 0.9142 | 0.685 us |
| GF16 AVX2 | low `(129,130)` | 64 KiB | 0.9968 | 1.0032 | 0.9988 / 1.0026 / 1.0031 | 9.346 us |

Across the five cells, median execution gain was 0.316%, ranging from an 8.57%
regression to a 9.19% gain.  Reuse-one setup amortization changes the median to
a 0.119% regression.  The exact-output combination ranged from an 8.28%
regression to a 9.91% gain, also below the 10% C1 threshold and not attributable
to inverse pruning alone.

The artifact retains every raw sample plus median/MAD.  Dispersion is small in
the decisive cells: for scalar `(17,18)`, mature and candidate medians/MADs are
1,289,759.844/307.812 ns and 1,181,192.344/333.156 ns; for AVX2 high `(33,32)`
they are 206,631.062/95.000 ns and 226,005.344/49.969 ns.  The harness does not
claim a confidence interval, so no stronger statistical claim is made.  That
limitation cannot rescue promotion: the best observed mature-prefix gain is
already below 10%, while two neighboring AVX2 cells regress by more than 2%.

Machine-readable results, including the isolation attestation, complete build
identity, parity digests, samples, medians, MADs, and merged summary, are in
`experiments/leopard2/non_power_of_two/c1/results/inverse_prefix_pinned_amd_9950x3d.json`.

## Final disposition

Reject production promotion.  The authoritative run did not clear the 10%
rule in any mature-prefix cell, and the AVX2 low `(33,34)` and high `(33,32)`
neighbors regress by 7.45% and 8.57%, respectively.  `AUTO` remains unchanged
on the mature encoder.  This is a negative promotion result, not an
inconclusive result.

Retain the correct executor, tests, counters, and benchmark form as an
experimental building block.  The next worthwhile implementation is a compact
regular-subtransform-plus-boundary schedule, followed by high-rate integration
with the inverse root accumulation sink.  Either must be requalified across
neighboring counts and backends before production dispatch changes.

## Reproduction

For correctness and raw diagnostic output from a Release build with tests and
benchmarks enabled:

    OMP_NUM_THREADS=1 ./build/leopard2_pruned_transform_test
    OMP_NUM_THREADS=1 ./build/leopard2_direct_encode_test
    ./build/leopard2_pruned_fuzz_smoke 14627333968358193854 16384
    taskset -c 4 env OMP_NUM_THREADS=1 ./build/bench_leopard2_sparse_encode \
      --profile low --field gf8 --k 17 --r 18 --bytes 65536 \
      --requested-parity 0-17 --backend scalar --iterations 32 \
      --samples 15 --warmups 3 --setup-iterations 8 --reuse 1,8,64

For production-linked pinned evidence, configure the benchmark without test
hooks, then provide a current isolation attestation to the fail-closed runner:

    cmake -S . -B build/benchprod -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/benchprod -j 128 \
      --target bench_leopard2_sparse_encode
    python3 tools/leopard2_sparse_encode_crossover.py pinned \
      --source . \
      --executable build/benchprod/bench_leopard2_sparse_encode \
      --result-dir build/inverse-prefix-pinned \
      --cell-manifest \
        experiments/leopard2/non_power_of_two/c1/inverse_prefix_cells.json \
      --iterations 32 --samples 15 --warmups 3 --setup-iterations 8 \
      --reuse 1,8,64 --memory-mib 1024 --workers 1 --cpu 4 \
      --isolation-attestation PATH_TO_ATTESTATION.json --no-resume
