# Optional CUDA backend

The Leopard2 CUDA backend is an end-of-project experimental track for large,
batched encode and erasure-decode workloads.  CPU correctness, compatibility,
and release readiness do not depend on CUDA, and the ordinary library has no
CUDA runtime or toolkit dependency.

## Build gate

CUDA is disabled by default.  A normal configuration does not probe for a CUDA
compiler:

    cmake -S . -B build/release

CMake project options use `-D`; opt in explicitly with:

    cmake -S . -B build/cuda -DLEO2_ENABLE_CUDA=ON

The spelling `cmake --enable-cuda` is not accepted by CMake for a project-defined
option; the command above is its canonical CMake equivalent.  Enabling the gate
requires a CMake release with first-class CUDA-language support and a usable CUDA
toolchain.  `LEO2_ENABLE_CUDA=ON` currently validates and enables only that
language gate.  CUDA codec sources, runtime dispatch, tests, and benchmarks will
be added behind it by the end-of-project CUDA Beads epic; until those tasks
close, no CUDA backend is advertised by the public API.

## Required CUDA completion gates

The future backend must remain isolated from the CPU wire profiles and scalar
oracle.  It must:

- reproduce CPU parity and recovered originals byte-for-byte;
- include host-to-device and device-to-host transfer time in end-to-end results;
- target large shards and realistic batches rather than tiny transfer-free
  microbenchmarks;
- handle unavailable drivers/devices as an optional-backend condition;
- add no CUDA headers, libraries, compiler probes, or runtime initialization to
  the default build; and
- keep CUDA-specific mathematics and storage choices from changing serialized
  profile identity.

The `leopard2_cuda_optional` test sets an invalid `CUDACXX`, configures without
specifying the CUDA option, verifies that CUDA remains disabled and unprobed,
and builds the ordinary CPU library to catch accidental CUDA compile or link
dependencies.  It also verifies that explicitly enabling CUDA with that invalid
compiler fails with a CUDA-specific diagnostic instead of silently omitting the
requested backend.
