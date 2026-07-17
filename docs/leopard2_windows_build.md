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
The CMake proof follows compile, include, and link mutations onto every object
library attached to `libleopard`; only the checked-in baseline definitions,
include paths, runtime links, and the AVX2 object's `/arch:AVX2` option are
accepted.  Directory-wide compile/link options, source or target property
bypasses, compiler-flag and IPO rewrites, generated/custom build commands, and
package discovery or build mutations inside unmodeled function/macro scopes
fail closed.  Required options, compiler probes, package/install commands, and
protected cache/control assignments are exact manifests rather than presence
checks, and their evaluation order is pinned.  Project options and non-forced
cache defaults remain symbolic so every supported command-line setting is
checked, as are all x86 processor spellings recognized by the build.  When
CUDA is off, tool discovery is limited to the checked objdump and shell probes
plus the benchmark-only Git lookup used to stamp locator measurements.  That
Git path admits exactly `rev-parse HEAD` and
`status --porcelain --untracked-files=normal`, with fixed result/output
variables, source working directory, guards, and evaluation order;
CUDA source properties, toolkit commands in tests, and CUDA/NVIDIA markers
anywhere in the reachable CPU source/header text are rejected.  The sole
marker-bearing test is the exact checked CUDA-optional self-test manifest.  The
Visual Studio proof likewise pins the solution projects and configuration
mappings, XML namespace, complete root evaluation phase order, item types,
configuration tools, properties, and Windows case-insensitive source aliases.
The constrained lexer understands arbitrary-delimiter CMake
bracket comments and balances bracket arguments; bracket arguments themselves
fail closed until the graph model has an explicit semantic use for them.

This is a proof about mutations authored in the checked-in `CMakeLists.txt`,
not a claim that an arbitrary invocation is trustworthy.  Environment flags,
compiler launchers selected by a toolchain file, compiler or linker wrappers,
and imported-package implementations are external trust boundaries and must be
controlled by the release environment.  The CMake graph interpreter currently
models the MSVC compiler branch only; clang-cl and the GNU-compatible compiler
branch still require their own exact mutation manifests.  The structural check
is also not native Windows compiler or runtime evidence; final MSVC and
clang-cl build, test, and dispatch validation still requires a Windows host.

Open limitation: this repository has no captured native Visual Studio 2015
load/build result.  The structural gate proves the checked-in metadata and
source/option graph, but the legacy-IDE compatibility claim remains pending a
real VS2015 build.  The approved conventional per-user property-sheet import
is also an external trust boundary whose contents are not present for this
host-side check.
