# Leopard2 backend isolation checkpoint

Status: first production runtime-dispatch checkpoint complete; the SIMD/backend
Bead remains open for fused butterflies, native NEON, and platform gates.

## Default contract

The CMake build does not append `-march=native`. On x86-64, baseline control,
field, API, and scheduler translation units remain at the platform SSE2 floor.
Three private fixed-kernel implementations are separate translation units:
scalar, SSSE3 (`-mssse3 -mno-avx`), and AVX2 (`-mavx2 -mno-avx512f`). The
feature probe is a fourth named member and is the only member allowed to
contain XGETBV.

`LEO2_BACKEND_VARIANT=auto` is now a production runtime-dispatched archive.
It selects AVX2, then SSSE3, then scalar according to the guarded feature probe.
The diagnostic `scalar`, `ssse3`, and `avx2` variants force selection of the
same isolated ops tables and fail initialization rather than substituting a
different backend. Backend choice changes kernels only, never wire identity.

Every variant archive is checked by `tools/check_leopard2_portable_isa.sh`.
The checker extracts and classifies each member. Baseline members reject
post-SSE2 instructions; the feature member permits XGETBV only; the SSSE3
member permits SSE3/SSSE3 but not AVX/SSE4; and the AVX2 member permits VEX
instructions while rejecting AVX-512 registers/masks and feature-probe
instructions. Metadata rejects target-raising flags outside the named objects,
unrelated flags inside them, all `-march`/`-mcpu`, and LTO leakage.

The denylist is deliberately conservative, not a proof about every future x86
instruction spelling. It matches only disassembled mnemonic fields, avoiding
symbol/operand false positives, and rejects every `-march` or `-mcpu` value
rather than guessing whether a toolchain treats it as baseline-compatible.
CTest registers the gate for every x86-64 variant when both
`objdump`/`llvm-objdump` and a POSIX `sh` are available. The backend matrix
treats a missing or unexecuted gate as a failure rather than accepting CTest's
successful "No tests were found" status.

## Runtime probing correction

AVX2 selection requires all of the following:

1. CPUID maximum basic leaf is at least 7.
2. CPUID leaf 1 reports AVX and OSXSAVE.
3. XGETBV reports both XMM and YMM state enabled in XCR0.
4. CPUID leaf 7 subleaf 0 reports AVX2.

XGETBV is therefore never executed on a CPU or OS that did not advertise the
required contract, and it is compiled only in the feature-probe member. Pure
classifier tests cover absent leaves, missing AVX/OSXSAVE, missing XMM/YMM
state, and the complete contract. No MSVC or clang-cl toolchain was available
on the evidence host, so native Windows build and runtime evidence remains open.

## Correctness and safety evidence

Runtime-dispatch checkpoint evidence on 2026-07-16 with GCC 13.3.0 on the
32-CPU allowed set:

- strict Release `-Wall -Wextra -Wpedantic -Werror`: 29/29 CTests passed;
- ASan plus UBSan: 29/29 CTests passed with leak and halt-on-error enabled;
- TSan, run with ASLR disabled for this host's runtime limitation: the
  16-thread fixed-kernel gate (1,024 executions) and the existing four-profile
  concurrent encoder gate (528 executions) passed without a report;
- four concurrent forced variants (`auto,scalar,ssse3,avx2`) passed 19 common
  deterministic tests each with no output mismatch; all four passed archive
  classification, and `auto` also passed optional-CUDA isolation. Source
  fingerprint: `42c4a1700dde44e25748e9f1135e3b27e608285b9c41c396b46a9c05acdb75b5`;
  merged matrix SHA-256:
  `003af6f7e03cb8c958e7d280c11350b408a59975d4e54891ed5b0b00e9a8ecaa`;
- startup KATs cover all 256 GF8 multiplier logs and all 256 byte values,
  multiply and multiply-add, unaligned buffers and tails; GF16 covers complete
  ALTMAP tiles plus compact tails and boundary logs; Common XOR covers zero
  through 257 bytes;
- the pure CPUID/XCR0 classifier and public ops-derived backend introspection
  passed; `auto` selected AVX2 on the evidence host;
- Clang 18 built sanitizer-instrumented libFuzzer copies for all four backend
  variants; 256 deterministic smoke executions per variant completed under
  ASan/UBSan, and coverage confirmed execution inside the selected isolated
  scalar, SSSE3, or AVX2 boundary;
- AArch64/SSE2NEON compile-only preservation passed at submodule commit
  `cad518a93b326f0f644b7972d488d04eaa2b0475`. This is not a native-NEON
  runtime or performance claim.

The preceding portable-safety checkpoint evidence is retained below for the
before/after record:

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

## Measured performance checkpoint

The table below reports the first isolated fixed-kernel tier, not the eventual
SIMD endpoint. Measurements were pinned to allowed CPU 15 with no simultaneous
benchmark, one thread, 7 measured samples, 3 warmups, and reuse 8. GB/s is
aggregate input for encode and offered received input for decode. The forced
variants select the same ops tables used by `auto`; `auto` selected AVX2 on the
evidence host. The scalar and SSSE3 rows are forced-variant measurements; the
AVX2 rows are measurements of the production `auto` build after it reported
AVX2, not of the forced AVX2 diagnostic build.

| Case | Backend | Encode GB/s | Decode GB/s |
| --- | ---: | ---: | ---: |
| high GF8 K=240 R=16, 64 KiB, L=4 | scalar | 3.036 | 2.365 |
| high GF8 K=240 R=16, 64 KiB, L=4 | SSSE3 | 8.252 | 5.468 |
| high GF8 K=240 R=16, 64 KiB, L=4 | AVX2 (`auto`) | 13.534 | 7.784 |
| low GF8 K=32 R=224, 64 KiB, L=16 | scalar | 0.338 | 6.338 |
| low GF8 K=32 R=224, 64 KiB, L=16 | SSSE3 | 0.964 | 16.396 |
| low GF8 K=32 R=224, 64 KiB, L=16 | AVX2 (`auto`) | 1.685 | 20.623 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | scalar | 0.423 | 0.256 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | SSSE3 | 3.223 | 1.746 |
| high GF16 K=1000 R=200, 16 KiB, L=8 | AVX2 (`auto`) | 6.047 | 2.841 |

The isolated AVX2 fixed-multiply/multiply-add and XOR tier materially recovers
the portable fallback cost without raising the archive ISA floor. It does not
yet recover the earlier whole-translation-unit diagnostic ceiling in every
cell because two-way butterflies, IFFT-XOR butterflies, and fused-four
butterflies still execute their baseline implementations. Those kernels are
the next extraction target. The scalar GF8 packing correction made during this
checkpoint restored the previous scalar encode range instead of accepting a
table-dispatch regression.

## Remaining production gates

- Extract the two-way, IFFT-XOR, and fused-four/radix butterfly families into
  the same private ops boundary and remeasure the whole-codec ceiling.
- Implement and test native NEON separately. The AArch64 evidence above is a
  compile-only preservation gate for the existing SSE2NEON translation path,
  not a native-NEON runtime or performance claim.
- Add native MSVC and clang-cl build/runtime evidence; neither toolchain was
  available on the checkpoint host.
- Add allocation-failure injection for backend table construction and remove
  the transitional duplicate legacy scalar multiplication tables once every
  legacy pre-publication caller has moved behind the ops boundary.
- Keep startup self-tests, static member classification, deterministic wire
  comparisons, sanitizers, and pinned performance cells as gates for each new
  isolated kernel family.
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
      -DCMAKE_SYSTEM_NAME=Linux \
      -DCMAKE_SYSTEM_PROCESSOR=aarch64 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=aarch64-linux-gnu-gcc \
      -DCMAKE_CXX_COMPILER=aarch64-linux-gnu-g++ \
      -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=OFF \
      -DENABLE_OPENMP=OFF
    cmake --build build/review-portable-aarch64 -j32 --target libleopard
