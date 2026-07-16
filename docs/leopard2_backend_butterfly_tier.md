# Leopard2 isolated two-way butterfly tier

Status: promoted as the second bounded production x86 backend tier. Fused-four
butterflies, native NEON, and wider experimental ISAs remain open work.

## Scope and call-site audit

The default portable build deliberately disables whole-translation-unit SSSE3
and AVX2 code generation in `LeopardFF8.cpp` and `LeopardFF16.cpp`. Before this
tier, fixed multiply, multiply-add, and XOR already used the immutable runtime
backend, but the two-way LCH butterflies fell back to separate baseline memory
passes. Because the radix-four wrappers decompose through these functions when
their old whole-TU intrinsic branches are unavailable, the audited two-way
call sites cover every nonzero-skew fallback butterfly in both fields, not only
the final odd transform layer.

This tier adds five private operations to the existing immutable table:

- forward and inverse GF8 two-way butterflies;
- the GF8 inverse-butterfly XOR accumulator used by the encoder; and
- forward and inverse GF16 two-way butterflies over legacy ALTMAP tiles.

Scalar, SSSE3, and AVX2 implementations preserve the equations and legacy
symbol representation. They fuse the multiply and both XOR dependencies in a
single pass over each pair of shards. GF8 accepts every byte count internally;
GF16 accepts complete symbols and retains the existing 64-byte ALTMAP tile plus
compact-tail layout. The accumulating GF8 operation treats all inputs as
read-only and requires four pairwise-disjoint buffers, matching its production
encoder call sites. No public API, coordinate map, field, or wire identifier
changed.

Fused-four extraction was intentionally excluded. Its six-table live set and
more complicated zero-skew specialization deserve a separate evidence gate.

## Correctness and isolation gates

Startup known-answer tests now cover forward, inverse, and accumulating
butterflies with guard bytes and unaligned pointers. GF8 covers all 256
logarithm values and lengths around 16/32/64-byte vector boundaries through
521 bytes. GF16 covers boundary logs, full ALTMAP tiles, and compact tails
through 194 bytes. The accumulating input buffers are checked for immutability.

Evidence collected on the 2026-07-16 x86 host:

- strict GCC 13.3 `-Wall -Wextra -Wpedantic -Werror`: 32/32 CTests passed;
- ASan plus UBSan, halt-on-error with leak detection: 32/32 passed;
- TSan under disabled ASLR: the 16-thread immutable-backend test completed
  1,024 executions and the four-profile encoder test completed 528 executions
  with no report;
- four forced variants (`auto,scalar,ssse3,avx2`) passed all 19 deterministic
  comparison tests, startup KATs, archive classification, and output matching;
  source fingerprint
  `ba4f96d2f38b79f76dc2ce8e4a81a5ab4670904c2945ba529132f42d0b6e2678`,
  merged matrix SHA-256
  `9594c8d4921ca472df1fec219753d4c2d5d5fd9717aec6647550bc0c8de18937`;
- the portable archive checker passed for Release and sanitizer builds. The
  sanitizer build exposed VEX `vmovd`; the exact move is now allowed because
  the existing runtime contract already proves AVX and OS-managed XMM/YMM
  state. A positive fixture was added while the fail-closed VEX allowlist and
  FMA/F16C/EVEX negative controls remain in force; and
- AArch64/SSE2NEON cross-compilation completed at submodule commit
  `cad518a93b326f0f644b7972d488d04eaa2b0475`. This is compile-only evidence,
  not a native ARM correctness or performance claim.

## Pinned end-to-end result

The authoritative machine-readable record is
`experiments/leopard2/backend_butterfly/results/abba_manifest.json`; its SHA-256
is `bcf67fcc07910b1d46fbac22239be10f3c54f2406533d89eea8ca3cbb6577643`.
The compact adjacent summary points to that manifest. The fail-closed runner
set and verified its own affinity as exactly CPU 15, exported
`OMP_NUM_THREADS=1` and `OMP_DYNAMIC=FALSE`, verified that CPUs 15 and 31 are
thread siblings, and recorded the coordinator's reservation of sibling 31 as
idle. It hashed the exact baseline and candidate binaries, every argv plus its
environment and affinity, and every raw stdout and stderr. The stable raw
evidence digest is
`668c534950fb59292161ba6d0a7425497d7af70346d3f36851adbe8ab11598b6`.

The runner rejects missing, duplicate, or relabelled ABBA slots and replays
all promotion statistics from raw benchmark JSON. Its self-test ran a complete
48-job mock campaign, then proved rejection of a modified stdout, missing raw
file, duplicate/relabelled manifest entry, and missing manifest entry. Runner
SHA-256 is
`40e2432ebd5361fa40dbb0a7f1e91e9967c88e77d7472c5935ea9536541a5935`;
the exact candidate source fingerprint it bound is
`97f994005bc68cda349e7d98425f3e3500cbedd6a09622ed8aacebbf85dc348f`.

Measurements used three A-B-B-A rounds, six invocation medians per build,
seven samples per invocation, three warmups, reuse eight, and one thread. All
48 invocations reported AVX2 and passed the round trip.

| Cell | Encode baseline to tier | Encode gain | Decode baseline to tier | Decode gain |
| --- | ---: | ---: | ---: | ---: |
| high GF8, K=240 R=16, 64 KiB, L=4 | 1171.602 to 935.422 us | 25.25% | 2151.791 to 1892.822 us | 13.68% |
| low GF8, K=32 R=224, 64 KiB, L=16 | 1247.569 to 994.001 us | 25.51% | 765.713 to 702.801 us | 8.95% |
| high GF16, K=1000 R=200, 16 KiB, L=8 | 2714.694 to 2218.318 us | 22.38% | 6781.538 to 5968.193 us | 13.63% |
| low GF16, K=128 R=1024, 16 KiB, L=16 | 2210.335 to 1912.097 us | 15.60% | 1283.556 to 1215.780 us | 5.57% |

The measured region has no regression; all eight encode/decode metrics improve
by at least 5%. This clears the normal production promotion threshold without
changing dispatcher behavior.

## Reproduction

Build and run the normal correctness gate:

    cmake -S . -B build/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build build/release -j8
    ctest --test-dir build/release -j8 --output-on-failure

Run the forced-backend differential matrix:

    python3 tools/leopard2_backend_matrix.py run \
      --source . \
      --build-root build/backend-matrix \
      --result-dir build/backend-matrix-results \
      --variants auto,scalar,ssse3,avx2 \
      --jobs 8 --variant-workers 4 --timeout 900 --no-resume

Exercise the evidence runner and its mutation tests before timing:

    python3 experiments/leopard2/backend_butterfly/run_abba.py self-test

The full runner command requires explicit baseline/candidate paths and hashes,
the baseline commit, backend-matrix fingerprint and digest, output directory,
benchmark CPU, reserved sibling, and coordinator reservation evidence. Run
`--help` for the complete fail-closed interface. Recheck retained raw evidence
and exact binaries with:

    python3 experiments/leopard2/backend_butterfly/run_abba.py verify \
      --manifest results/leopard2/backend-butterfly-abba-20260716-v2/campaign/manifest.json \
      --baseline results/leopard2/backend-butterfly-abba-20260716-v2/bin/baseline-bench_leopard2 \
      --candidate results/leopard2/backend-butterfly-abba-20260716-v2/bin/candidate-bench_leopard2
