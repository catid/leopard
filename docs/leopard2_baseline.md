# Leopard2 legacy baseline

This baseline was captured on 2026-07-16 before changing Leopard's algorithms or
public API.  The source revision was
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` from `origin/master`.

## Host and toolchain

- Ubuntu 24.04.4 LTS, Linux 6.8.0-134-generic, x86-64.
- AMD Ryzen 9 9950X3D, one socket, 16 physical cores and 32 logical CPUs.
- The process affinity mask exposes CPUs 0-31 in one NUMA node.  The requested
  128-core and multi-node measurements cannot be performed on this host.
- GCC/G++ 13.3.0, CMake 3.28.3, Ninja 1.11.1.
- Runtime features include AVX2, AVX-512BW/VBMI, GFNI, and VPCLMULQDQ.
- The governor reported `powersave` with `amd-pstate-epp`.
- `perf_event_paranoid=4`, so unprivileged hardware counters are unavailable.

Both Release and Debug builds used the repository's existing CMake configuration,
OpenMP, `-march=native`, and 32 parallel build jobs.  `ctest` reported that the
legacy project defines no CTest tests.  The existing random test executable and
experiment executable are intentionally unbounded, so they were run under a
timeout:

| Configuration | Limit | Successful decodes | Reported errors |
| --- | ---: | ---: | ---: |
| Release random test | 30 s | 331 | 0 |
| Debug random test | 15 s | 91 | 0 |
| Release experiment | 15 s | 7,828 | 0 |

The compiler reported existing warnings in the tests and a pointer-before-array
idiom around `FFTSkew - 1` in the FF8 and FF16 transform code.  These are baseline
observations, not new Leopard2 warnings.

## Pinned legacy throughput sample

The following compatibility baseline used CPU 0, `OMP_NUM_THREADS=1`, a 64 KiB
shard, seven repetitions, and the old benchmark harness.  Values are median input
throughput in decimal MB/s.  The legacy harness reports its best inner trial, so
these numbers are useful for regression orientation but are not confidence-interval
quality release results.

| Field/profile-shaped case | K | R | Encode MB/s | Decode MB/s |
| --- | ---: | ---: | ---: | ---: |
| GF8 balanced | 128 | 128 | 5,691.05 | 2,430.07 |
| GF8 high rate | 240 | 16 | 15,839.5 | 4,228.13 |
| GF16 high rate | 1,000 | 200 | 5,613.84 | 1,404.18 |
| XOR special case | 129 | 1 | 79,756.1 | 82,883.8 |

## Golden vectors

`tests/leopard2/legacy_golden.csv` records deterministic FNV-1a hashes of all
encoded recovery shards and recovered original shards for a fixed 64-byte input
generator.  It spans power-of-two boundaries, the GF8/GF16 selection boundary,
balanced codes, and the one-parity special case.  The generator and full local
logs are retained under the ignored `.research/leopard2/baseline` directory.

## Reproduction

From the repository root:

    cmake -S . -B build/baseline-release -G Ninja \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-march=native
    cmake --build build/baseline-release -j 32
    ctest --test-dir build/baseline-release -j 32 --output-on-failure
    timeout 30s build/baseline-release/leopard_test

The CPU count of 32 is the detected affinity limit on this machine, rather than a
hard-coded project limit.
