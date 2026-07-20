# Leopard2 backend determinism matrix

`LEO2_BACKEND_VARIANT` is a diagnostic CMake cache setting with five values:
`auto`, `scalar`, `ssse3`, `avx2`, and `avx512`. The default is `auto`. In that mode the
project builds baseline control flow at the x86-64 SSE2 floor and adds named
SSSE3, AVX2, and AVX-512VL object members compiled with target-local ISA flags
when the compiler supports them. Startup CPUID/XCR0 probing selects the
immutable context-baseline table only after its GF8/GF16/XOR known-answer tests
pass. A forced variant caps both `AUTO` and explicit context selection at that
backend; `avx512` may still qualify lower tables. In an ordinary `auto` build,
an eligible codec may additionally qualify one immutable operation table as
described below. Neither choice changes the field, coordinate profile, or wire
bytes.

The forced variants apply to the x86 library and independently built fuzzer
copy. `scalar` selects table-backed scalar fixed multiplication and portable
XOR. `ssse3` requires runtime SSSE3 and selects isolated 128-bit nibble-table
kernels. `avx2` additionally requires AVX, OSXSAVE, XMM/YMM XCR0 state, and the
AVX2 CPUID bit before selecting 256-bit fixed multiplication and Common XOR.
`avx512` requires AVX2 plus AVX-512F/BW/VL and OS-enabled opmask/ZMM state. It
selects the explicit AVX-512VL table, whose current data path remains 256 bits
wide while using the expanded register file. `AUTO` reports AVX2 as the
context baseline on a qualifying x86 host. On AMD family 1Ah/model 44h only,
a legacy-high GF16 codec may select the qualified AVX-512VL operation table for
a full-output encode when `K >= 8`, `N >= 16`, `2 <= R <= 4096`, and the shard
length is an exact 64-byte multiple from 64 bytes through 4 MiB inclusive.
Unknown CPU models and all cells outside those fixed bounds retain AVX2.
All baseline FF/control objects have legacy whole-TU optional code generation
disabled.

The current CMake diagnostic variants are x86-only. CMake checks that the
compiler accepts each isolated object's target flags; it does not require the
build host to implement that ISA. The Python matrix separately skips execution
of a SIMD diagnostic when the host lacks the selected feature, and records
`avx512` as unavailable when any required compiler flag or CPU feature is
missing. A direct CMake request fails at configure time when its isolated
backend cannot be built. `scalar`
performs neither optional-feature nor compiler-ISA check. Variants do not
silently substitute a different implementation. All five archives have the
same available member isolation; the matrix derives its exact compile-source
contract from the CMake probe results so an unsupported optional object is not
mistaken for graph drift. The variants are deterministic selection diagnostics,
while `auto` is the production runtime-dispatched binary.

Run the standard-library-only matrix from the repository root:

    python3 tools/leopard2_backend_matrix.py self-test
    python3 tools/leopard2_backend_matrix.py run

The runner detects the process affinity instead of assuming CPU numbers,
checks host and compiler support for the SIMD variants, and creates one
isolated build per variant.
It builds and runs the startup-KAT/synthetic-feature/concurrency gate, frozen
legacy golden vectors, the public API suite, a
fixed-seed random smoke suite, the independent production-constant and bare-LCH
differential, the direct-generator transform differential, and the AUTO encode
selector's positive/negative calibration boundaries. Every variant
runs the static per-member portable-ISA archive gate. The matrix also
compares the reusable decode-plan schedule differential, allocation, and
concurrency test across all executable backends. The `auto` build
runs `leopard2_cuda_optional`, proving that a normal build does not need a CUDA
compiler or toolkit. Each forced build also makes the public API suite assert
that runtime backend introspection reports the requested scalar, SSSE3, or
AVX2 backend; the AVX-512 variant is checked the same way when available.
Lower forced builds also assert that an explicit AVX-512 context is rejected,
so matching output alone cannot conceal a failed force or qualification cap.
The failure matrix verifies that an optional AVX-512 KAT failure during
eligible AUTO codec setup falls back to AVX2 without changing codec creation.
If CMake did not register or execute the portable-ISA gate (for example because
`objdump` or a POSIX `sh` is unavailable), the matrix fails with an actionable
reason; it does not accept CTest's zero exit status for an empty selection.
Test processes are pinned to distinct allowed CPUs and normally use one OpenMP
worker.  The context-backend isolation test deliberately uses four workers to
exercise explicit ops-table propagation inside transforms.  The isolated
builds run concurrently and divide up to 128 build jobs between them. The heavier boundary differential remains a
normal CTest target (`leopard2_boundaries`) and can use the full test affinity.

Results are written under `results/leopard2/backend-matrix/`, with full command
logs, compiler identity, sorted CPU flags, selected CMake value, executable
hashes, and normalized stdout/stderr hashes. `matrix.json` fails if any test
fails, source inputs change during the run, or normalized output from a forced
variant differs from `auto`. Unsupported forced variants are recorded as
`unavailable`, not emulated. Completed results with the same source, compiler,
CPU, generator, and job configuration are resumed; pass `--no-resume` for a
fresh execution.

Useful restricted runs include:

    python3 tools/leopard2_backend_matrix.py run --variants auto,scalar
    python3 tools/leopard2_backend_matrix.py run --jobs 128 --no-resume

This matrix is a correctness and wire-determinism gate. It is not an
authoritative performance benchmark: timings are deliberately omitted, and
test pinning is intended to isolate correctness output rather than measure
throughput.
