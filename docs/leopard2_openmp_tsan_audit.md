# GF16 OpenMP ThreadSanitizer audit

This audit resolves `leopard-79h.40.11`.  GCC 13 ThreadSanitizer reports data
races when the GF16 encoder executes consecutive OpenMP worksharing loops.  The
same report occurs in the standalone control
`tests/leopard2/openmp_tsan_barrier_probe.cpp`, which contains only two
standards-valid consecutive `omp parallel for` loops over the same array.

## Disposition

No Leopard data race was confirmed.  The GCC report is a libgomp
happens-before false positive:

- `LEO_OPENMP_PARALLEL_FOR` expands to `omp parallel for` without `nowait`.
  The implicit barrier completes before the calling thread starts the next
  region or advances the transform-stage loop.
- The first source-staging region writes suffix slots
  `[m_truncated, m)`.  Complete source groups write distinct four-slot ranges
  below that suffix.  A partial group runs serially after both regions.
- Each radix-four stage iteration advances by `dist4` and touches exactly its
  own `[r, r + dist4)` shard-pointer block.  Stages depend on one another, but
  their separate worksharing regions are ordered by the implicit barriers.
- The final radix-two region assigns each `i` the distinct pair
  `{i, i + dist}`.  The vector XOR helper likewise partitions by distinct
  output shards.

The first GCC report is a write by the caller to the next outlined-loop stack
descriptor racing a worker's alleged read from the prior descriptor.  A
non-halting diagnostic then reports the expected downstream transform
dependencies as races because the same missing libgomp happens-before edges
are propagated.  Within a single stage, the ranges above are disjoint.

The minimal control is deliberately not part of the default CMake graph.  It
exists to distinguish this tool/runtime limitation from Leopard code and can be
compiled directly with either toolchain.

## Reproduction

The stable host invocation disables address randomization because ordinary
GCC TSan startup intermittently fails before `main` with `unexpected memory
mapping`.

GCC 13 codec reproduction:

    cmake -S . -B /home/catid/leopard-builds/root-openmp-tsan -G Ninja \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_C_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
      -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
      -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=OFF \
      -DENABLE_OPENMP=ON
    cmake --build /home/catid/leopard-builds/root-openmp-tsan \
      --target leopard2_gf16_legacy_encoder_matrix_test -j 28
    env TSAN_OPTIONS=halt_on_error=1:history_size=7 \
      OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 \
      taskset -c 0-3,16-19 setarch x86_64 -R \
      /home/catid/leopard-builds/root-openmp-tsan/\
leopard2_gf16_legacy_encoder_matrix_test

This exits 66 at the first report.  The ignored raw report has SHA-256
`0a772d9280d54ad7b1d9d3676b7f6e8099ed7236ca6ec07b7980164985736600`.

Minimal GCC control:

    g++-13 -std=c++11 -O1 -g -fopenmp -fsanitize=thread \
      -fno-omit-frame-pointer \
      tests/leopard2/openmp_tsan_barrier_probe.cpp \
      -o /home/catid/leopard-builds/root-openmp-tsan/\
openmp_tsan_barrier_probe
    env TSAN_OPTIONS=halt_on_error=1:history_size=7 \
      OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 \
      taskset -c 0-3,16-19 setarch x86_64 -R \
      /home/catid/leopard-builds/root-openmp-tsan/\
openmp_tsan_barrier_probe

It exits 66 with the same prior-region descriptor pattern.  Its ignored raw
report has SHA-256
`eee0a32bf0f017d9dfee3a919181cb1cde3f6fcc79b6ce7b65b799f379d317ef`.

Clang 18 warns that the dynamically linked, non-instrumented OpenMP runtime
requires `ignore_noninstrumented_modules=1` to avoid runtime-internal false
positives.  With that setting and the same `setarch`/affinity/OMP environment,
all of these completed without a TSan report:

- the standalone two-loop control;
- the 20-profile, 360-case GF16 legacy compatibility matrix;
- 528 concurrent encode executions across four profiles; and
- the main API test with 710 recovered shards and 72 plan executions.

This setting weakens inspection inside libomp, so it is not the sole basis for
the disposition.  The independent minimal reproduction, OpenMP barrier
semantics, explicit disjoint-range audit, multi-thread exact-parity runs, and
no-OpenMP tests are required corroborating evidence.

## Correctness controls

The Release GF16 compatibility matrix passed all 360 cases independently with
`OMP_NUM_THREADS=1,2,4,8`, producing the same 2,607,300 parity bytes and 52,200
subset-parity bytes each time.  A separate `ENABLE_OPENMP=OFF` Release build
passed the same matrix plus 528 concurrent encode executions.  The complete
current-source Release suite passed 84 of 84 tests before this documentation
change.

No synchronization or hot-path instruction was added.  Adding another barrier
would duplicate an existing language guarantee and would not make the GCC
runtime visible to TSan; suppressing application frames would instead risk
hiding real defects.  Future reports that name a within-stage overlapping
range, rather than a dependency separated by an OpenMP barrier, must be treated
as new evidence and investigated independently.
