# Leopard2 backend determinism matrix

`LEO2_BACKEND_VARIANT` is a diagnostic CMake cache setting with four values:
`auto`, `scalar`, `ssse3`, and `avx2`. The default is `auto`. In that mode the
project adds no backend definitions or target-local ISA flags, preserving the
historical build behavior exactly. A forced variant changes implementation
kernels only; it never changes the field, coordinate profile, or wire bytes.

The forced variants currently apply to the x86 library and independently built
fuzzer copy. `scalar` disables the legacy SSSE3 arithmetic selection and AVX2
code generation. The library still contains its baseline SSE2 memory kernels
and dormant SSSE3 translation units, so this is an algorithmic oracle rather
than a distributable no-SIMD binary. `ssse3` enables the legacy SSSE3 selection
but disables AVX code generation. `avx2` enables AVX2 code generation while
disabling AVX-512 generation. Runtime CPUID checks remain active.

Forced variants fail configuration on a non-x86 target or when the compiler
does not accept the required target flags. They do not silently substitute a
different implementation. The target-local flags follow the repository's
historical global `-march=native`; this diagnostic option does not broaden the
default binary's ISA requirements.

Run the standard-library-only matrix from the repository root:

    python3 tools/leopard2_backend_matrix.py self-test
    python3 tools/leopard2_backend_matrix.py run

The runner detects the process affinity instead of assuming CPU numbers,
checks host and compiler support, and creates one isolated build per variant.
It builds and runs the frozen legacy golden vectors, the public API suite, a
fixed-seed random smoke suite, and the independent direct-generator transform
differential. The `auto` build also runs
`leopard2_cuda_optional`, proving that a normal build does not need a CUDA
compiler or toolkit. Each forced build also makes the public API suite assert
that runtime backend introspection reports the requested scalar, SSSE3, or
AVX2 backend; matching output alone cannot conceal a failed force control.
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
