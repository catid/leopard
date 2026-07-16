# Leopard2 backend isolation checkpoint

Status: bounded production-safety checkpoint; the SIMD/backend Bead remains
open.  This checkpoint deliberately chooses a portable default over silently
shipping a host-specific archive.

## Default contract

The CMake build no longer probes for or appends `-march=native`.  On x86-64,
the default `LEO2_BACKEND_VARIANT=auto` archive is compiled at the platform
SSE2 baseline.  Its fixed-multiplier path is scalar and its XOR path may use
baseline SSE2.  SSSE3 and AVX2 are not compiled into that archive unless the
caller explicitly supplies an ISA compilation boundary.

The existing `ssse3` and `avx2` CMake variants remain diagnostic, opt-in
whole-archive builds.  They are useful for correctness and performance
comparison, but are not portable binaries and must not be redistributed as if
they were runtime-dispatched.  The production follow-up is to move these
kernels into ISA-specific translation units or target functions, leaving all
public API, initialization, and dispatch control flow at the platform baseline.

The default archive is checked by
`tools/check_leopard2_portable_isa.sh`.  The check rejects project-supplied
ISA-raising compiler flags found in Make flags, `compile_commands.json`, or
Ninja metadata, and disassembles every object in the archive. Its mnemonic
denylist covers SSE3, SSSE3, SSE4, AVX-family, BMI, LZCNT, POPCNT, CX16, ADX,
AES, SHA, carryless multiplication, XGETBV/XSAVE, and other common post-SSE2
compiler targets. Real synthetic archives exercise SSSE3, SSE4.1, and AVX2
negative controls; separate fixtures exercise all three supported metadata
formats and `-march=native`.

The denylist is deliberately conservative, not a proof about every future x86
instruction spelling. It matches only disassembled mnemonic fields, avoiding
symbol/operand false positives, and rejects every `-march` or `-mcpu` value
rather than guessing whether a toolchain treats it as baseline-compatible.
CTest registers the gate for default and scalar x86-64 builds only when both
`objdump`/`llvm-objdump` and a POSIX `sh` are available. The backend matrix
treats a missing or unexecuted gate as a failure rather than accepting CTest's
successful "No tests were found" status.

## Runtime probing correction

AVX2 diagnostic builds now require all of the following before setting
`CpuHasAVX2`:

1. CPUID maximum basic leaf is at least 7.
2. CPUID leaf 1 reports AVX and OSXSAVE.
3. XGETBV reports both XMM and YMM state enabled in XCR0.
4. CPUID leaf 7 subleaf 0 reports AVX2.

XGETBV is therefore never executed on a CPU or OS that did not advertise the
required contract. The portable default does not compile XGETBV at all. In
particular, MSVC's general intrinsic availability no longer enables SSSE3 or
AVX2 in `auto`: only explicitly forced MSVC diagnostic variants enable those
paths. A preprocessing-only `_MSC_VER=1930` simulation produced no optional
`LEO_TRY_*` macro for `auto` or `scalar`, `LEO_TRY_SSSE3` for forced SSSE3, and
both optional macros for forced AVX2. No MSVC or clang-cl toolchain was
available on the evidence host, so this simulation and source review still
need native Windows build and runtime evidence.

## Correctness and safety evidence

Evidence collected on 2026-07-16 with GCC 13.3.0 on an AMD Ryzen 9 9950X3D
(32 allowed logical CPUs, 16 physical cores, one NUMA node):

- default Release build and CTest: 25/25 passed, current-source final pass
  completed in 4.43 seconds;
- strict `-Wall -Wextra -Wpedantic -Werror` build and CTest: 25/25 passed,
  current-source final pass completed in 53.45 seconds with OpenMP enabled;
- ASan plus UBSan build and CTest: 25/25 passed, current-source final pass
  completed in 9.97 seconds;
- explicit SSSE3 build and CTest: 24/24 passed on CPUs 0-15;
- explicit AVX2 build and CTest: 24/24 passed on CPUs 16-31;
- compile-only AArch64/SSE2NEON check: pinned submodule commit
  `cad518a93b326f0f644b7972d488d04eaa2b0475` and GCC/G++ 13.3 cross-build of
  `libleopard` completed 6/6 build steps; no runtime or performance claim;
- four-variant backend matrix (`auto,scalar,ssse3,avx2`): passed with no
  mismatches using 4 concurrent workers and 8 build jobs per variant; both
  `auto` and `scalar` explicitly executed and passed the archive gate.  The
  source fingerprint was
  `ae94466b71d24bec50521d223f8196c19e0f0bba6a86bab7f03fdf6775f13bd7`;
  SHA-256 of the merged `matrix.json` was
  `01a96f8d86b3f19a57af8904e71d4f81eb9264dfab8fef666cc873354e14f4b7`;
- default archive disassembly: no forbidden mnemonic and no `-march`,
  `-mssse3`, `-mavx*`, `-mbmi*`, `-mlzcnt`, or `-mpopcnt` compile flag.

Exercising the scalar path under UBSan exposed pre-existing alignment-unsafe
GF8 reference loads and stores.  They now use alignment-safe copies; the same
unaligned random test that found the issue passes under ASan/UBSan.

All profile, direct-oracle, legacy-golden, encoder, decoder, arbitrary-count,
and transform-differential tests passed in the default and both opt-in SIMD
builds.  The ISA change does not alter field arithmetic, coordinates, or wire
bytes.

The strict build uses the repository's normal OpenMP-enabled configuration.
An additional strict build with OpenMP disabled stopped on the existing
unguarded `#pragma omp` directives under `-Werror=unknown-pragmas`; this is a
configuration-warning limitation, not a failed codec test.

## Measured performance cost

The safety fallback is not a performance endpoint.  The table below reports
single-thread medians from CPU 16, with its sibling idle, `OMP_NUM_THREADS=1`,
9 measured samples, 3 warmups, and reuse 8.  GB/s is aggregate input for encode
and offered received input for decode.

| Case | Backend | Encode GB/s | Decode GB/s |
| --- | ---: | ---: | ---: |
| high GF8 K=240 R=16, 64 KiB, L=4 | scalar default | 3.005 | 2.395 |
| high GF8 K=240 R=16, 64 KiB, L=4 | SSSE3 opt-in | 11.955 | 8.151 |
| high GF8 K=240 R=16, 64 KiB, L=4 | AVX2 opt-in | 20.657 | 13.047 |
| low GF8 K=32 R=224, 64 KiB, L=16 | scalar default | 0.336 | 2.780 |
| low GF8 K=32 R=224, 64 KiB, L=16 | SSSE3 opt-in | 1.217 | 9.727 |
| low GF8 K=32 R=224, 64 KiB, L=16 | AVX2 opt-in | 1.924 | 15.476 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | scalar default | 0.530 | 0.329 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | SSSE3 opt-in | 3.536 | 2.075 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | AVX2 opt-in | 6.591 | 3.861 |

Relative to the AVX2 diagnostic build, the portable scalar default is about
82-92% slower in these cells.  That regression is the reason this is only a
safety checkpoint and the SIMD/backend Bead remains open; silently retaining a
host-only default would instead risk illegal instructions before dispatch.

## Remaining production gates

- Isolate SSSE3 and AVX2 kernels in ISA-specific translation units or target
  functions and connect them to baseline runtime dispatch.
- Keep backend self-tests and deterministic wire-output comparisons at each
  isolated boundary.
- Re-run strict, sanitizer, fuzz, old/new compatibility, and end-to-end
  performance matrices after isolation; recover the measured SIMD throughput
  without raising the archive-wide ISA floor.
- Implement and test native NEON separately.  The existing SSE2NEON path is a
  translation layer and this checkpoint makes no native-NEON performance claim.
- Treat AVX-512, GFNI, SVE/SVE2, and other experimental kernels as later
  evidence-gated work rather than extending the default ISA contract.

## Reproduction

Default build, static gate, and tests:

    cmake -S . -B build/portable-isa-release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build build/portable-isa-release -j32
    ctest --test-dir build/portable-isa-release -j32 --output-on-failure
    sh tools/check_leopard2_portable_isa.sh \
      "$(command -v objdump)" \
      build/portable-isa-release/liblibleopard.a \
      build/portable-isa-release \
      "$(command -v cc)" \
      "$(command -v ar)"

Explicit diagnostic builds replace the build directory and add either
`-DLEO2_BACKEND_VARIANT=ssse3` or `-DLEO2_BACKEND_VARIANT=avx2`.

Fresh four-variant matrix used for the fingerprint above:

    python3 tools/leopard2_backend_matrix.py run \
      --source . \
      --build-root build/portable-isa-matrix-v3 \
      --result-dir build/portable-isa-matrix-results-v3 \
      --variants auto,scalar,ssse3,avx2 \
      --jobs 32 \
      --variant-workers 4 \
      --timeout 900 \
      --no-resume

Compile-only AArch64/SSE2NEON preservation check:

    git submodule update --init --depth 1 sse2neon
    cmake -S . -B build/review-portable-aarch64 -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=aarch64-linux-gnu-gcc \
      -DCMAKE_CXX_COMPILER=aarch64-linux-gnu-g++ \
      -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=OFF \
      -DENABLE_OPENMP=OFF
    cmake --build build/review-portable-aarch64 -j32 --target libleopard
