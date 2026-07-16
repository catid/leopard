# Windows builds

The CMake build is the authoritative Windows build for Leopard and Leopard2.
It discovers the current MSVC or clang-cl toolchain, builds the ordinary CPU
library without a CUDA dependency, and keeps optional x86 instruction sets in
dedicated translation units so runtime dispatch does not raise the library's
baseline ISA.

For example, from an x64 Native Tools prompt for Visual Studio 2022:

    cmake -S . -B build\vs2022 -G "Visual Studio 17 2022" -A x64 -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build\vs2022 --config Release --parallel
    ctest --test-dir build\vs2022 -C Release --output-on-failure

`LEO2_ENABLE_CUDA` remains `OFF` unless it is explicitly enabled.  See
`README_CUDA.md` for the separate experimental CUDA gate.

## Checked-in Visual Studio 2015-format solution

`proj/Leopard.sln` retains metadata intended for Visual Studio 2015 legacy
consumers.
The checked-in legacy metadata is pinned consistently to Visual Studio 14,
MSBuild ToolsVersion 14.0, and the v140 platform toolset; this follows the
solution's own `VisualStudioVersion = 14.0.25420.1` declaration rather than the
stale ToolsVersion 15.0 value inherited when the projects were downgraded from
their original Visual Studio 2017 metadata.
Its `Leopard` static-library project includes both public APIs and the same
production C++ translation units as CMake.  It deliberately keeps the baseline
files free of optional SSSE3/AVX2 code generation, compiles the SSSE3 backend in
its own file, and applies `/arch:AVX2` only to
`Leopard2BackendAVX2.cpp`.  Whole-program optimization is disabled for the
static library so link-time optimization cannot dissolve that ISA boundary.

The legacy solution does not contain the full Leopard2 correctness suite,
fuzzers, benchmarks, package-install tests, or optional CUDA work.  Use CMake
for release validation.  The host-runnable structural gate can be run on any
platform with:

    python3 tests/proj/test_leopard_vcxproj.py

It compares the hand-maintained project with CMake's production source graph
and rejects missing or duplicate sources, filter drift, a project-wide ISA
increase, unproved package or object-library link edges, or a CUDA dependency.
Package discovery and link commands inside unmodeled CMake function or macro
scopes fail closed.  This structural check is not native Windows compiler or
runtime evidence; final MSVC and clang-cl build, test, and dispatch validation
still requires a Windows host.

Open limitation: this repository has no captured native Visual Studio 2015
load/build result.  The structural gate proves the checked-in metadata and
source/option graph, but the legacy-IDE compatibility claim remains pending a
real VS2015 build.  The approved conventional per-user property-sheet import
is also an external trust boundary whose contents are not present for this
host-side check.
