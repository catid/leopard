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

## CPU-only install contract

The default install remains a standalone CPU package:

    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
    cmake --build build/release --target install

It installs `leopard.h`, `leopard2.h`, the static CPU library, and a relocatable
CMake config package.  Consumers can use `find_package(leopard CONFIG REQUIRED)`
and link `leopard::libleopard`.  The exported target carries only the CPU
threading/OpenMP dependencies selected when the library was built; it does not
contain CUDA targets, headers, compiler settings, runtime libraries, or device
initialization.  A future CUDA backend must remain in a separate opt-in export
set so enabling or installing it cannot change this default target contract.

Both public headers provide a C ABI, so application sources may be C or C++.
The static library itself is implemented in C++, however, and its CMake package
requires the consuming project to enable the `CXX` language so CMake can select
the portable C++ link driver and runtime.  A C application should therefore use
`project(my_app LANGUAGES C CXX)` and may still compile its source as `.c`.  A
C-only CMake project is rejected during `find_package` with that instruction
instead of failing later with unresolved C++ or OpenMP runtime symbols.

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

The `leopard2_cuda_optional` test configures once with an invalid `CUDACXX` and
once with that environment variable absent, without specifying the CUDA option,
and verifies that CUDA remains disabled and unprobed.  It builds and
stage-installs the ordinary CPU library, rejects CUDA-named
install artifacts and CUDA references in installed headers or CMake metadata,
relocates the complete install tree, then configures, builds, links, and runs a
separate `find_package` consumer with CUDA still unavailable.  This catches
accidental CUDA compile, link, export, public-header, and absolute-install-path
dependencies.  The test also verifies that explicitly enabling CUDA with that
invalid compiler fails with a CUDA-specific diagnostic instead of silently
omitting the requested backend.
