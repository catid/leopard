# R10 RS(256,K) AVX2 reproduction

This directory records the 2026-08-11 reproduction of the decoding workload
from R10. The compact summary is human- and machine-readable; the compressed
bundle retains every raw JSON record used by it.

- Candidate source: `a0d781c1ad6f68c49d35d7fe165433aa8b7fc981`
- Exact Leopard main: `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`
- R10 harness SHA-256: `369944dc1d84c8e933a251b61e074fea5ece5fa38c2e39683018cc09233038ee`
- Raw bundle SHA-256: `8548a5f7df1908dbddbdc2cf78d3e6ecbfb3ef883dd448ab8e4d75c24b5aeb3e`

The candidate build explicitly disables GFNI and AVX-512 code generation and
forces the AVX2 backend. It was built from a clean detached worktree. Both
frozen candidate executables contain zero GFNI, ZMM, opmask, and `vpternlog`
instructions.

## Candidate build

From a clean checkout at the candidate commit:

    cmake -S . -B build/r10-avx2 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=ON \
      -DLEO2_BUILD_FUZZERS=OFF \
      -DLEO2_ENABLE_CUDA=OFF \
      -DLEO2_BACKEND_VARIANT=avx2 \
      -DLEO2_FLAG_MAVX512F=FALSE \
      -DLEO2_FLAG_MAVX512BW=FALSE \
      -DLEO2_FLAG_MAVX512VL=FALSE \
      -DLEO2_FLAG_MGFNI=FALSE
    cmake --build build/r10-avx2 --target bench_leopard2 -j1
    c++ -std=c++11 -O3 -DNDEBUG \
      -Wall -Wextra -Wpedantic -Werror \
      -DLEO2_R10_SOURCE_COMMIT=\"a0d781c1ad6f68c49d35d7fe165433aa8b7fc981\" \
      -I. experiments/leopard2/r10_rs256_benchmark.cpp \
      build/r10-avx2/libleopard.a -fopenmp -lpthread \
      -o build/r10-avx2/r10_rs256_benchmark

## One paper-matrix cell

The authoritative campaign pinned this command to one allowed physical CPU and
serialized it with the repository benchmark lock:

    taskset -c 13 build/r10-avx2/r10_rs256_benchmark \
      --k 224 --bytes 1024 --profile high --mode auto \
      --pattern-seed 2 --data-seed 17 \
      --reuse 512 --iterations 25 --warmup 16 \
      --json k224-auto.json

Use `--mode generic` for the same-binary generic control. Low-rate paper cells
use `--profile low`. Five pattern seeds were used: 2, 3, 5, 7, and 11.

## Validate retained evidence

    python3 -m json.tool a0d781c-avx2-summary.json >/dev/null
    sha256sum a0d781c-avx2-raw.tar.gz
    tar -tzf a0d781c-avx2-raw.tar.gz

The summary explicitly separates the exact detached-main comparison from the
in-tree legacy timing used for the four fixed missing-original counts.
