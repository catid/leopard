# Leopard2 backend determinism matrix

`LEO2_BACKEND_VARIANT` is a diagnostic CMake cache setting with four values:
`auto`, `scalar`, `ssse3`, and `avx2`. The default is `auto`. In that mode the
project adds no backend definitions or target-local ISA flags. Since the
project no longer appends `-march=native`, a normal x86-64 build is limited to
the SSE2 platform baseline and uses the scalar fixed-multiplier path. A forced
variant changes implementation kernels only; it never changes the field,
coordinate profile, or wire bytes.

The forced variants currently apply to the x86 library and independently built
fuzzer copy. `scalar` disables the legacy SSSE3 arithmetic selection and AVX2
code generation and requires no optional CPU feature or compiler ISA flag. Its
x86-64 archive contains only baseline SSE2 memory kernels plus scalar field
arithmetic. `ssse3` enables the legacy SSSE3 selection but disables AVX code
generation. `avx2` enables AVX2 code generation while disabling AVX-512
generation. Runtime CPUID checks remain active.

The current CMake diagnostic variants are x86-only. CMake checks that the
compiler accepts each SIMD variant's target flags; it does not require the
build host to implement that ISA. The Python matrix separately skips execution
of a SIMD diagnostic when the host lacks the selected feature. `scalar`
performs neither optional-feature nor compiler-ISA check. Variants do not
silently substitute a different implementation. SSSE3 and AVX2 are still
opt-in whole-archive diagnostics, not production runtime-dispatched binaries.
The default and scalar variants do not broaden the binary's ISA requirements.

Run the standard-library-only matrix from the repository root:

    python3 tools/leopard2_backend_matrix.py self-test
    python3 tools/leopard2_backend_matrix.py run

The runner detects the process affinity instead of assuming CPU numbers,
checks host and compiler support for the SIMD variants, and creates one
isolated build per variant.
It builds and runs the frozen legacy golden vectors, the public API suite, a
fixed-seed random smoke suite, the independent production-constant and bare-LCH
differential, and the direct-generator transform differential. The `auto` and
`scalar` builds also run the static portable-ISA archive gate. The matrix also
compares the reusable decode-plan schedule differential, allocation, and
concurrency test across all executable backends. The `auto` build
runs `leopard2_cuda_optional`, proving that a normal build does not need a CUDA
compiler or toolkit. Each forced build also makes the public API suite assert
that runtime backend introspection reports the requested scalar, SSSE3, or
AVX2 backend; matching output alone cannot conceal a failed force control.
If CMake did not register or execute the portable-ISA gate (for example because
`objdump` or a POSIX `sh` is unavailable), the matrix fails with an actionable
reason; it does not accept CTest's zero exit status for an empty selection.
Test processes are pinned to distinct allowed CPUs and
use one OpenMP worker; the isolated builds run concurrently and divide up to
128 build jobs between them. The heavier boundary differential remains a
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
