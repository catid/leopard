cmake_minimum_required(VERSION 3.7)
project(leopard)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(cmake/Leopard2BenchmarkAttestation.cmake)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 11)

set(LIB_SOURCE_FILES
        leopard.cpp
        leopard.h
        leopard2.cpp
        leopard2.h
        Leopard2Backend.cpp
        Leopard2Backend.h
        Leopard2BackendScalar.cpp
        Leopard2CpuFeatures.cpp
        Leopard2Direct.h
        Leopard2Plan.cpp
        Leopard2Plan.h
        LeopardCommon.cpp
        LeopardCommon.h)

set(BENCH_SOURCE_FILES
        tests/benchmark.cpp)

set(EXPERIMENT_SOURCE_FILES
        tests/experiments.cpp)

option(LEO2_BUILD_TESTS "Build Leopard2 correctness tests" ON)
option(LEO2_BUILD_BENCHMARKS "Build Leopard benchmark programs" ON)
option(LEO2_BUILD_ALLK_DIAGNOSTIC
    "Build the source-attested all-K Leopard1 comparison benchmark" OFF)
option(LEO2_BUILD_FUZZERS "Build Leopard2 libFuzzer targets" OFF)
option(LEO2_ENABLE_CUDA "Build the optional Leopard2 CUDA backend" OFF)
option(LEOPARD_ENABLE_GF8 "Include the GF(2^8) codec" ON)
option(LEOPARD_ENABLE_GF16 "Include the GF(2^16) codec" ON)
option(LEO2_PORTABLE_ISA_RELEASE_AUDIT
    "Require the strict x86-64 portable-ISA release audit" OFF)
option(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN
    "Preformat AVX2 direct-repair source schedules in decode plans" OFF)
option(LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT
    "Fuse tiny equal-rounded GF8 one-shot repair without a heap plan" ON)
if(LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT)
    add_definitions(-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1)
else()
    add_definitions(-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=0)
endif()
option(LEO2_EXPERIMENT_CAUCHY_LOG_REUSE
    "Reuse Cauchy cross-difference logs during direct-repair setup" ON)
if(LEO2_EXPERIMENT_CAUCHY_LOG_REUSE)
    add_definitions(-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1)
else()
    add_definitions(-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=0)
endif()
option(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE
    "Enable default-off candidate legacy-high GF8/AVX2 direct diagnostics"
    OFF)
option(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
    "Allow the legacy-high direct-encode experiment to change AUTO dispatch"
    ON)
option(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR
    "Disable the promoted AVX2 T=8 encoder for same-source diagnostics"
    OFF)
option(LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING
    "Enable the promoted shortened/punctured AVX2 T=8 binding fast path"
    ON)
option(LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
    "Enable the promoted two-message-block AVX2 T=8 binding fast path"
    ON)
option(LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING
    "Enable the qualified 65..1024-byte AVX2 T=8 ragged binding selector"
    ON)
option(LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT
    "Enable promoted generalized GF8/AVX2 one-loss direct repair"
    ON)
option(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL
    "Build the default-off GF8/AVX2 Algorithm 4 P32/N64/B64 terminal"
    OFF)
set(LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE "0" CACHE STRING
    "Default-off small GF8 direct-repair experiment: 0=transform, 1=output-major, 2=source-major")
set_property(CACHE LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE PROPERTY STRINGS
    0 1 2)
if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE MATCHES "^[012]$")
    message(FATAL_ERROR
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE must be 0, 1, or 2")
endif()

if(LEO2_BUILD_ALLK_DIAGNOSTIC AND NOT LEO2_BUILD_BENCHMARKS)
    message(FATAL_ERROR
        "LEO2_BUILD_ALLK_DIAGNOSTIC requires LEO2_BUILD_BENCHMARKS=ON")
endif()

if(NOT LEOPARD_ENABLE_GF8 AND NOT LEOPARD_ENABLE_GF16)
    message(FATAL_ERROR
        "At least one of LEOPARD_ENABLE_GF8 or LEOPARD_ENABLE_GF16 must be ON")
endif()
if(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL AND NOT LEOPARD_ENABLE_GF8)
    message(FATAL_ERROR
        "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL requires GF8")
endif()

# Diagnostic-only backend selection.  auto is the portable default and does
# not add definitions or target-local ISA options.  Forced variants hold the
# wire format fixed while making backend determinism testable in isolated
# builds.
set(LEO2_BACKEND_VARIANT "auto" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512")
set_property(CACHE LEO2_BACKEND_VARIANT PROPERTY STRINGS
    auto scalar ssse3 avx2 avx512)
string(TOLOWER "${LEO2_BACKEND_VARIANT}" LEO2_BACKEND_VARIANT_NORMALIZED)
if(NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "auto" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "scalar" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "ssse3" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "avx2" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "avx512")
    message(FATAL_ERROR
        "Invalid LEO2_BACKEND_VARIANT='${LEO2_BACKEND_VARIANT}'. "
        "Expected auto, scalar, ssse3, avx2, or avx512")
endif()
set(LEO2_BACKEND_VARIANT "${LEO2_BACKEND_VARIANT_NORMALIZED}" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512" FORCE)
set_property(CACHE LEO2_BACKEND_VARIANT PROPERTY STRINGS
    auto scalar ssse3 avx2 avx512)

# CUDA is deliberately opt-in.  In particular, the default configuration must
# not probe for a CUDA compiler or make CUDA a transitive build dependency.
# Backend sources will be added behind this gate by the end-of-project CUDA
# work tracked in Beads.
if(LEO2_ENABLE_CUDA)
    enable_language(CUDA)
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()

# Do not raise the ISA floor of the whole library.  The default archive must be
# safe to move between machines that implement the platform baseline.  Faster
# diagnostic variants are explicit below; production runtime-dispatched ISA
# backends must live behind translation-unit or function target boundaries.

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    check_cxx_compiler_flag("/O2" CXX_FLAG_O2)
    check_cxx_compiler_flag("/Oy" CXX_FLAG_Oy)
    check_cxx_compiler_flag("/Zi" CXX_FLAG_Zi)
    check_cxx_compiler_flag("/W4" CXX_FLAG_W4)

    if(CXX_FLAG_O2)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    endif()

    if(CXX_FLAG_Oy)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Oy")
    endif()

    if(CXX_FLAG_Zi)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
    endif()

    if(CXX_FLAG_W4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    check_cxx_compiler_flag("-G" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Minform=inform" CXX_FLAG_MINFORM_INFORM)

    if(CXX_FLAG_g)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif(CXX_FLAG_g)

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif(CXX_FLAG_O3)

    if(CXX_FLAG_MINFORM_INFORM)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
    endif(CXX_FLAG_MINFORM_INFORM)
else()
    check_cxx_compiler_flag("-G" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Wall" CXX_FLAG_Wall)
    check_cxx_compiler_flag("-Wextra" CXX_FLAG_Wextra)

    if(CXX_FLAG_g)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif(CXX_FLAG_g)

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif(CXX_FLAG_O3)

    if(CXX_FLAG_Wall)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    endif(CXX_FLAG_Wall)

    if(CXX_FLAG_Wextra)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    endif(CXX_FLAG_Wextra)
endif()

find_package(OpenMP)
find_package(Threads REQUIRED)
cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" ON "OPENMP_FOUND" OFF)
if(ENABLE_OPENMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP support requested, but could not be found")
    endif(NOT OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(ENABLE_OPENMP)

# This data-only classifier is configuration-aware.  It is loaded after all
# project-wide C++ flags are finalized so a later CTest registration can select
# the archive audit independently for each active configuration.
include(cmake/Leopard2SanitizerClassification.cmake)

# BEGIN LEO2_TARGET_ARCHITECTURE_MODEL
set(LEO2_TARGET_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}")
set(LEO2_APPLE_UNIVERSAL_TARGET OFF)
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    # CMAKE_SYSTEM_PROCESSOR describes the host in common Apple cross-build
    # workflows.  CMAKE_OSX_ARCHITECTURES is the actual target contract.  A
    # multi-slice archive cannot safely consume an OBJECT library whose
    # -mssse3/-mavx2 option is valid for only its x86 slice, so universal builds
    # deliberately retain the portable scalar/native-ARM graph.
    if("${CMAKE_OSX_ARCHITECTURES}" MATCHES ";")
        set(LEO2_APPLE_UNIVERSAL_TARGET ON)
        set(LEO2_TARGET_PROCESSOR "")
        message(STATUS
            "Leopard2 Apple universal build: optional x86 ISA objects disabled")
    else()
        set(LEO2_TARGET_PROCESSOR "${CMAKE_OSX_ARCHITECTURES}")
    endif()
endif()

set(LEO2_X86_TARGET OFF)
set(LEO2_X86_64_TARGET OFF)
if(LEO2_TARGET_PROCESSOR STREQUAL "x86" OR
   LEO2_TARGET_PROCESSOR STREQUAL "X86" OR
   LEO2_TARGET_PROCESSOR STREQUAL "x86_64" OR
   LEO2_TARGET_PROCESSOR STREQUAL "x86_64h" OR
   LEO2_TARGET_PROCESSOR STREQUAL "amd64" OR
   LEO2_TARGET_PROCESSOR STREQUAL "AMD64" OR
   LEO2_TARGET_PROCESSOR STREQUAL "i386" OR
   LEO2_TARGET_PROCESSOR STREQUAL "i486" OR
   LEO2_TARGET_PROCESSOR STREQUAL "i586" OR
   LEO2_TARGET_PROCESSOR STREQUAL "i686")
    set(LEO2_X86_TARGET ON)
endif()
if(LEO2_TARGET_PROCESSOR STREQUAL "x86_64" OR
   LEO2_TARGET_PROCESSOR STREQUAL "x86_64h" OR
   LEO2_TARGET_PROCESSOR STREQUAL "amd64" OR
   LEO2_TARGET_PROCESSOR STREQUAL "AMD64")
    set(LEO2_X86_64_TARGET ON)
endif()
# END LEO2_TARGET_ARCHITECTURE_MODEL

# This opt-in gate is intentionally narrower than ordinary portable builds.
# It is a release/CI evidence mode for single-config, native x86-64 builds with
# GNU-compatible archive/disassembly metadata.  Normal Windows, cross, and
# developer builds do not acquire any of these tool requirements.
if(LEO2_PORTABLE_ISA_RELEASE_AUDIT)
    if(MSVC OR WIN32)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT is not implemented for "
            "MSVC/Windows COFF archives; use the default portable build and "
            "the Visual Studio structural audit")
    endif()
    if(NOT UNIX OR APPLE OR CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT currently requires a native "
            "non-Apple Unix build")
    endif()
    if(NOT LEO2_X86_64_TARGET)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT currently supports x86-64 only")
    endif()
    if(NOT ((CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
             CMAKE_C_COMPILER_ID STREQUAL "Clang") AND
            (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
             CMAKE_CXX_COMPILER_ID STREQUAL "Clang")))
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires supported "
            "GCC or Clang C/C++ toolchains")
    endif()
    if(CMAKE_GENERATOR MATCHES "Multi-Config")
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires a single-configuration "
            "generator")
    endif()
    if(NOT LEO2_BUILD_TESTS)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires LEO2_BUILD_TESTS=ON")
    endif()
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires "
            "CMAKE_BUILD_TYPE=Release")
    endif()
    if(NOT CMAKE_EXPORT_COMPILE_COMMANDS)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires "
            "CMAKE_EXPORT_COMPILE_COMMANDS=ON")
    endif()
    if(NOT CMAKE_AR OR NOT EXISTS "${CMAKE_AR}")
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires a usable CMAKE_AR")
    endif()
endif()

add_library(leopard STATIC ${LIB_SOURCE_FILES})
# Keep source-tree consumers that linked the historical target working while
# making `leopard` the canonical CMake target.  Installed consumers use the
# conventional `leopard::leopard` export.  ALIAS targets do not create a
# second archive.
add_library(libleopard ALIAS leopard)
if(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN)
    target_compile_definitions(leopard PRIVATE
        LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=1)
endif()
if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
    target_compile_definitions(leopard PRIVATE
        LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1
        LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO=$<BOOL:${LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO}>)
endif()
target_compile_definitions(leopard PRIVATE
    LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING}>
    LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING}>
    LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING}>)
if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE STREQUAL "0")
    # This diagnostic selector is consumed only by leopard2.cpp.  Keep its
    # trusted compile-command delta source-scoped for both the production and
    # test-hook archives; unrelated library objects must remain byte-for-byte
    # reproducible across modes.
    set_property(SOURCE leopard2.cpp APPEND PROPERTY COMPILE_DEFINITIONS
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=${LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE}")
endif()
if(LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT)
    # The candidate changes production plan selection and publishes a
    # pure-AVX2 callback.  Keep every consumer of the selector/KAT macro in
    # lockstep while leaving the ordinary archive free of experimental code.
    set_property(SOURCE
        leopard2.cpp
        Leopard2Backend.cpp
        Leopard2BackendAVX2.cpp
        tests/leopard2/test_api.cpp
        APPEND PROPERTY COMPILE_DEFINITIONS
        LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1)
endif()
if(LEOPARD_ENABLE_GF8)
    target_sources(leopard PRIVATE LeopardFF8.cpp LeopardFF8.h)
else()
    target_compile_definitions(leopard PRIVATE NO_LEO_HAS_FF8=1)
endif()
if(LEOPARD_ENABLE_GF16)
    target_sources(leopard PRIVATE LeopardFF16.cpp LeopardFF16.h)
else()
    target_compile_definitions(leopard PRIVATE NO_LEO_HAS_FF16=1)
endif()
if(LEO2_X86_TARGET)
    target_compile_definitions(leopard PRIVATE
        LEO2_DISABLE_SSSE3_CODEGEN=1
        LEO2_DISABLE_AVX2_CODEGEN=1)
endif()
target_include_directories(leopard PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(leopard PUBLIC Threads::Threads)

# Optional x86 kernels are compiled in dedicated object libraries.  None of
# leopard's baseline translation units is built above the x86-64 SSE2 floor;
# runtime feature checks select one immutable ops table after startup KATs.
set(LEO2_HAVE_SSSE3_BACKEND OFF)
set(LEO2_HAVE_AVX2_BACKEND OFF)
set(LEO2_HAVE_AVX512_BACKEND OFF)
set(LEO2_HAVE_GFNI_BACKEND OFF)
if(LEO2_X86_TARGET)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        check_cxx_compiler_flag("-mssse3" LEO2_FLAG_MSSSE3)
        check_cxx_compiler_flag("-mno-avx" LEO2_FLAG_MNO_AVX)
        if(LEO2_FLAG_MSSSE3 AND LEO2_FLAG_MNO_AVX)
            add_library(leopard2_backend_ssse3 OBJECT
                Leopard2BackendSSSE3.cpp)
            target_include_directories(leopard2_backend_ssse3 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_options(leopard2_backend_ssse3 PRIVATE
                -mssse3 -mno-avx)
            set(LEO2_HAVE_SSSE3_BACKEND ON)
        endif()

        check_cxx_compiler_flag("-mavx2" LEO2_FLAG_MAVX2)
        check_cxx_compiler_flag("-mno-avx512f" LEO2_FLAG_MNO_AVX512F)
        check_cxx_compiler_flag("-falign-functions=64"
            LEO2_FLAG_FALIGN_FUNCTIONS_64)
        if(LEO2_FLAG_MAVX2 AND LEO2_FLAG_MNO_AVX512F)
            add_library(leopard2_backend_avx2 OBJECT
                Leopard2BackendAVX2.cpp)
            if(LEOPARD_ENABLE_GF8)
                target_sources(leopard2_backend_avx2 PRIVATE
                    Leopard2BackendAVX2Xor.cpp)
            endif()
            target_include_directories(leopard2_backend_avx2 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_options(leopard2_backend_avx2 PRIVATE
                -mavx2 -mno-avx512f)
            if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
                target_compile_options(leopard2_backend_avx2 PRIVATE
                    -falign-functions=64)
            endif()
            set(LEO2_HAVE_AVX2_BACKEND ON)
        endif()

        check_cxx_compiler_flag("-mavx512f" LEO2_FLAG_MAVX512F)
        check_cxx_compiler_flag("-mavx512bw" LEO2_FLAG_MAVX512BW)
        check_cxx_compiler_flag("-mavx512vl" LEO2_FLAG_MAVX512VL)
        check_cxx_compiler_flag("-mprefer-vector-width=256"
            LEO2_FLAG_MPREFER_VECTOR_WIDTH_256)
        if(LEO2_HAVE_AVX2_BACKEND AND LEO2_FLAG_MAVX512F AND
           LEO2_FLAG_MAVX512BW AND LEO2_FLAG_MAVX512VL AND
           LEO2_FLAG_MPREFER_VECTOR_WIDTH_256)
            # Compile the reviewed AVX2 algorithms a second time with the
            # AVX-512VL register file.  The data path remains 256 bits wide;
            # F/BW/VL are required so the compiler may address ymm16..ymm31.
            add_library(leopard2_backend_avx512 OBJECT
                Leopard2BackendAVX512.cpp)
            target_include_directories(leopard2_backend_avx512 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_definitions(leopard2_backend_avx512 PRIVATE
                LEO2_HAVE_AVX2_BACKEND=1)
            target_compile_options(leopard2_backend_avx512 PRIVATE
                -mavx2 -mavx512f -mavx512bw -mavx512vl
                -mprefer-vector-width=256)
            if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
                target_compile_options(leopard2_backend_avx512 PRIVATE
                    -falign-functions=64)
            endif()
            set(LEO2_HAVE_AVX512_BACKEND ON)
        endif()

        check_cxx_compiler_flag("-mgfni" LEO2_FLAG_MGFNI)
        if(LEO2_HAVE_AVX2_BACKEND AND LEO2_FLAG_MGFNI)
            # Compile the reviewed AVX2 algorithms a third time with the GFNI
            # affine multiplication kernels and its own backend identity.  The
            # data path stays VEX 256-bit; -mno-avx512f keeps the member
            # EVEX-free so the portable-ISA audit can hold it to the AVX2
            # operand contract plus the GF2P8 family.
            add_library(leopard2_backend_gfni OBJECT
                Leopard2BackendGFNI.cpp)
            target_include_directories(leopard2_backend_gfni PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_definitions(leopard2_backend_gfni PRIVATE
                LEO2_HAVE_AVX2_BACKEND=1
                LEO2_HAVE_GFNI_BACKEND=1)
            target_compile_options(leopard2_backend_gfni PRIVATE
                -mavx2 -mgfni -mno-avx512f)
            if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
                target_compile_options(leopard2_backend_gfni PRIVATE
                    -falign-functions=64)
            endif()
            set(LEO2_HAVE_GFNI_BACKEND ON)
        endif()
    elseif(MSVC)
        # x64 MSVC accepts SSSE3 intrinsics without raising other translation
        # units.  AVX2 remains isolated behind /arch:AVX2.
        add_library(leopard2_backend_ssse3 OBJECT
            Leopard2BackendSSSE3.cpp)
        target_include_directories(leopard2_backend_ssse3 PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        set(LEO2_HAVE_SSSE3_BACKEND ON)
        check_cxx_compiler_flag("/arch:AVX2" LEO2_FLAG_ARCH_AVX2)
        if(LEO2_FLAG_ARCH_AVX2)
            add_library(leopard2_backend_avx2 OBJECT
                Leopard2BackendAVX2.cpp)
            if(LEOPARD_ENABLE_GF8)
                target_sources(leopard2_backend_avx2 PRIVATE
                    Leopard2BackendAVX2Xor.cpp)
            endif()
            target_include_directories(leopard2_backend_avx2 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_options(leopard2_backend_avx2 PRIVATE /arch:AVX2)
            set(LEO2_HAVE_AVX2_BACKEND ON)
        endif()
    endif()
endif()

if(LEO2_HAVE_SSSE3_BACKEND)
    if(NOT LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_backend_ssse3 PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_backend_ssse3 PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_ssse3>)
    target_compile_definitions(leopard PRIVATE LEO2_HAVE_SSSE3_BACKEND=1)
endif()
if(LEO2_HAVE_AVX2_BACKEND)
    target_compile_definitions(leopard2_backend_avx2 PRIVATE
        LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING}>)
    if(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL)
        # Keep the default-off exact terminal in its own object.  Enabling the
        # experiment must not perturb GCC's call graph, constprop clones, or
        # function placement in the mature AVX2 backend translation unit.
        add_library(leopard2_low_p32_b64_avx2 OBJECT
            Leopard2LowP32B64AVX2.cpp)
        target_include_directories(leopard2_low_p32_b64_avx2 PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_low_p32_b64_avx2 PRIVATE
            LEO2_HAVE_AVX2_BACKEND=1
            LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1)
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_low_p32_b64_avx2 PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
            target_compile_options(leopard2_low_p32_b64_avx2 PRIVATE
                -mavx2 -mno-avx512f)
            if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
                target_compile_options(leopard2_low_p32_b64_avx2 PRIVATE
                    -falign-functions=64)
            endif()
        elseif(MSVC)
            target_compile_options(leopard2_low_p32_b64_avx2 PRIVATE
                /arch:AVX2)
        endif()
        target_sources(leopard PRIVATE
            $<TARGET_OBJECTS:leopard2_low_p32_b64_avx2>)
    endif()
    if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
        target_compile_definitions(leopard2_backend_avx2 PRIVATE
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1)
    endif()
    if(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
        target_compile_definitions(leopard2_backend_avx2 PRIVATE
            LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_backend_avx2 PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_backend_avx2 PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_avx2>)
    target_compile_definitions(leopard PRIVATE
        LEO2_HAVE_AVX2_BACKEND=1
        LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=$<BOOL:${LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL}>)
elseif(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL)
    message(FATAL_ERROR
        "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL requires the AVX2 backend")
endif()
if(LEO2_HAVE_AVX512_BACKEND)
    if(NOT LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_backend_avx512 PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_backend_avx512 PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_avx512>)
    target_compile_definitions(leopard PRIVATE LEO2_HAVE_AVX512_BACKEND=1)
endif()
if(LEO2_HAVE_GFNI_BACKEND)
    if(NOT LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_backend_gfni PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_backend_gfni PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_gfni>)
    target_compile_definitions(leopard PRIVATE LEO2_HAVE_GFNI_BACKEND=1)
endif()
if(ENABLE_OPENMP)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(leopard PUBLIC OpenMP::OpenMP_CXX)
    else()
        # FindOpenMP did not provide imported targets before CMake 3.9.  Keep
        # the 3.7 project floor usable by propagating its compiler-driver flag;
        # for static archives this is also the link flag that selects the
        # OpenMP runtime on GCC and Clang.
        target_link_libraries(leopard PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Keep the installed CPU package independent of optional backends.  CUDA
# targets, once implemented, must use a separate opt-in export set rather than
# becoming transitive dependencies of leopard::leopard.
set(LEOPARD_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"
    CACHE STRING "Install directory for Leopard CMake package files")
configure_package_config_file(
    cmake/leopardConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake"
    INSTALL_DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")
install(TARGETS leopard
    EXPORT leopardTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES leopard.h leopard2.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake"
    DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")
install(EXPORT leopardTargets
    FILE leopardTargets.cmake
    NAMESPACE leopard::
    DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")

if(NOT LEO2_BACKEND_VARIANT STREQUAL "auto")
    if(NOT LEO2_X86_TARGET)
        message(FATAL_ERROR
            "LEO2_BACKEND_VARIANT=${LEO2_BACKEND_VARIANT} is unavailable on "
            "target processor '${LEO2_TARGET_PROCESSOR}'; use auto")
    endif()

    if(LEO2_BACKEND_VARIANT STREQUAL "scalar")
        target_compile_definitions(leopard PRIVATE
            LEO2_BACKEND_FORCE_SCALAR=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "ssse3")
        if(NOT LEO2_HAVE_SSSE3_BACKEND)
            message(FATAL_ERROR
                "LEO2_BACKEND_VARIANT=ssse3 is unavailable for this compiler")
        endif()
        target_compile_definitions(leopard PRIVATE
            LEO2_BACKEND_FORCE_SSSE3=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "avx2")
        if(NOT LEO2_HAVE_AVX2_BACKEND)
            message(FATAL_ERROR
                "LEO2_BACKEND_VARIANT=avx2 is unavailable for this compiler")
        endif()
        target_compile_definitions(leopard PRIVATE
            LEO2_BACKEND_FORCE_AVX2=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "avx512")
        if(NOT LEO2_HAVE_AVX512_BACKEND)
            message(FATAL_ERROR
                "LEO2_BACKEND_VARIANT=avx512 is unavailable for this compiler")
        endif()
        target_compile_definitions(leopard PRIVATE
            LEO2_BACKEND_FORCE_AVX512=1)
    endif()
endif()

if(LEO2_BUILD_BENCHMARKS)
    add_executable(bench_leopard ${BENCH_SOURCE_FILES})
    target_link_libraries(bench_leopard leopard)

    add_executable(bench_leopard2 bench/leopard2/benchmark.cpp)
    target_link_libraries(bench_leopard2 leopard)
    leopard2_enable_benchmark_source_attestation(bench_leopard2)

    add_executable(bench_leopard2_prevalidated_batch
        bench/leopard2/benchmark.cpp)
    target_compile_definitions(bench_leopard2_prevalidated_batch PRIVATE
        LEO2_BENCHMARK_PREVALIDATED_BATCH=1)
    target_link_libraries(bench_leopard2_prevalidated_batch leopard)
    leopard2_enable_benchmark_source_attestation(
        bench_leopard2_prevalidated_batch)

    if(LEOPARD_ENABLE_GF8)
        add_executable(bench_leopard2_batch_preflight
            bench/leopard2/batch_preflight_benchmark.cpp)
        target_link_libraries(bench_leopard2_batch_preflight leopard)
    endif()

    add_executable(bench_leopard2_locator
        bench/leopard2/locator_benchmark.cpp)
    target_link_libraries(bench_leopard2_locator leopard)
    if(NOT LEOPARD_ENABLE_GF8)
        target_compile_definitions(bench_leopard2_locator PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(NOT LEOPARD_ENABLE_GF16)
        target_compile_definitions(bench_leopard2_locator PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    set(LEO2_LOCATOR_SOURCE_GIT_SHA "unknown")
    set(LEO2_LOCATOR_SOURCE_GIT_TREE "unknown")
    set(LEO2_LOCATOR_SOURCE_DIRTY 1)
    find_program(LEO2_LOCATOR_GIT_EXECUTABLE NAMES git)
    if(LEO2_LOCATOR_GIT_EXECUTABLE)
        execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_GIT_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_GIT_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(LEO2_LOCATOR_GIT_RESULT EQUAL 0)
            set(LEO2_LOCATOR_SOURCE_GIT_SHA ${LEO2_LOCATOR_GIT_OUTPUT})
        endif()
        execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE} rev-parse HEAD^{tree}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_TREE_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_TREE_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(LEO2_LOCATOR_TREE_RESULT EQUAL 0)
            set(LEO2_LOCATOR_SOURCE_GIT_TREE ${LEO2_LOCATOR_TREE_OUTPUT})
        endif()
        execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE}
                status --porcelain --untracked-files=normal
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_STATUS_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_STATUS_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(LEO2_LOCATOR_STATUS_RESULT EQUAL 0 AND
           LEO2_LOCATOR_STATUS_OUTPUT STREQUAL "")
            set(LEO2_LOCATOR_SOURCE_DIRTY 0)
        endif()
    endif()
    # The standard benchmark supports an untimed schema-v5 identity preflight
    # as well as its ordinary timed schemas.  Its generated header is refreshed
    # at build time, rather than freezing configure-time Git state.
    target_compile_definitions(bench_leopard2_locator PRIVATE
        LEO2_LOCATOR_SOURCE_GIT_SHA="${LEO2_LOCATOR_SOURCE_GIT_SHA}"
        LEO2_LOCATOR_SOURCE_DIRTY=${LEO2_LOCATOR_SOURCE_DIRTY})

    if(LEO2_BUILD_ALLK_DIAGNOSTIC)
        add_executable(bench_leopard2_allk bench/leopard2/benchmark.cpp)
        target_link_libraries(bench_leopard2_allk leopard)
        leopard2_enable_benchmark_source_attestation(bench_leopard2_allk)
    endif()

    # The timed paths use ordinary internal transform entry points; untimed
    # closure also links the independent direct-generator oracle and exercises
    # production public decode. The target deliberately links the production
    # archive even when the test suite is enabled: test-hook counters perturb
    # the compared encoder paths by different amounts and invalidate crossover
    # measurements.
    if(LEOPARD_ENABLE_GF8 AND LEOPARD_ENABLE_GF16)
        add_library(leopard2_sparse_encode_benchmark_object OBJECT
            bench/leopard2/sparse_encode_benchmark.cpp)
        add_library(leopard2_sparse_encode_oracle_object OBJECT
            tests/leopard2/direct_oracle.cpp)
        target_include_directories(leopard2_sparse_encode_benchmark_object PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            tests/leopard2)
        target_include_directories(leopard2_sparse_encode_oracle_object PRIVATE
            tests/leopard2)
        target_compile_definitions(leopard2_sparse_encode_benchmark_object PRIVATE
            LEO2_SPARSE_ENCODE_SOURCE_GIT_SHA="${LEO2_LOCATOR_SOURCE_GIT_SHA}"
            LEO2_SPARSE_ENCODE_SOURCE_DIRTY=${LEO2_LOCATOR_SOURCE_DIRTY}
            LEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS=0)
        add_executable(bench_leopard2_sparse_encode
            $<TARGET_OBJECTS:leopard2_sparse_encode_benchmark_object>
            $<TARGET_OBJECTS:leopard2_sparse_encode_oracle_object>)
        target_link_libraries(bench_leopard2_sparse_encode leopard)
        set_property(TARGET bench_leopard2_sparse_encode APPEND PROPERTY
            LINK_DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteSparseEncodeEvidenceSidecar.cmake)
        add_custom_command(TARGET bench_leopard2_sparse_encode POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                -DOUTPUT=$<TARGET_FILE:bench_leopard2_sparse_encode>.leopard2-evidence
                -DEXECUTABLE=$<TARGET_FILE:bench_leopard2_sparse_encode>
                -DPRODUCTION_ARCHIVE=$<TARGET_FILE:leopard>
                -DBENCHMARK_OBJECT=$<TARGET_OBJECTS:leopard2_sparse_encode_benchmark_object>
                -DORACLE_OBJECT=$<TARGET_OBJECTS:leopard2_sparse_encode_oracle_object>
                -DBENCHMARK_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/bench_leopard2_sparse_encode.dir/link.txt
                -DPRODUCTION_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/leopard.dir/link.txt
                -DBUILD_PROGRAM=${CMAKE_MAKE_PROGRAM}
                -DBUILD_ROOT=${CMAKE_CURRENT_BINARY_DIR}
                -DBUILD_GENERATOR=${CMAKE_GENERATOR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteSparseEncodeEvidenceSidecar.cmake
            VERBATIM)
    endif()

    if(LEOPARD_ENABLE_GF16)
        add_executable(bench_leopard2_gf16_high_attribution
            experiments/leopard2/gf16_high_encode/attribution.cpp)
        target_link_libraries(
            bench_leopard2_gf16_high_attribution leopard)
    endif()

    add_executable(experiment_leopard ${EXPERIMENT_SOURCE_FILES})
endif()

if(LEO2_BUILD_TESTS)
    enable_testing()

    # Production-library public-query probe.  The Python operation model reads
    # its JSONL output and independently reconstructs every byte of range,
    # pointer-map, ragged-input, and path-specific work scratch.
    add_executable(leopard2_decode_scratch_probe
        tests/leopard2/test_decode_scratch_probe.cpp)
    target_link_libraries(leopard2_decode_scratch_probe leopard)

    # Test hooks change codec construction, scratch sizing, backend setup, and
    # a few diagnostic dispatch decisions in addition to exporting symbols.
    # Keep the production/install target invariant when tests are enabled by
    # compiling a separate, non-installed archive for the tests that need
    # those hooks.  Its x86 backends are separate objects as well: sharing the
    # production objects would either omit their hooks or contaminate the
    # runtime-dispatched production archive.
    add_library(leopard_test_hooks STATIC EXCLUDE_FROM_ALL
        ${LIB_SOURCE_FILES})
    target_compile_definitions(leopard_test_hooks PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=1)
    endif()
    if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO=$<BOOL:${LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO}>)
    endif()
    target_compile_definitions(leopard_test_hooks PRIVATE
        LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING}>
        LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING}>
        LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING}>)
    if(LEOPARD_ENABLE_GF8)
        target_sources(leopard_test_hooks PRIVATE LeopardFF8.cpp LeopardFF8.h)
    else()
        target_compile_definitions(leopard_test_hooks PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(LEOPARD_ENABLE_GF16)
        target_sources(leopard_test_hooks PRIVATE
            LeopardFF16.cpp LeopardFF16.h)
    else()
        target_compile_definitions(leopard_test_hooks PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()
    if(LEO2_X86_TARGET)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_DISABLE_SSSE3_CODEGEN=1
            LEO2_DISABLE_AVX2_CODEGEN=1)
    endif()
    target_include_directories(leopard_test_hooks PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(leopard_test_hooks PUBLIC Threads::Threads)

    if(LEO2_HAVE_SSSE3_BACKEND)
        add_library(leopard2_backend_ssse3_test_hooks OBJECT EXCLUDE_FROM_ALL
            Leopard2BackendSSSE3.cpp)
        target_include_directories(leopard2_backend_ssse3_test_hooks PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_ssse3_test_hooks PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_ssse3_test_hooks PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_ssse3_test_hooks PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
            target_compile_options(leopard2_backend_ssse3_test_hooks PRIVATE
                -mssse3 -mno-avx)
        endif()
        target_sources(leopard_test_hooks PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_ssse3_test_hooks>)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_HAVE_SSSE3_BACKEND=1)
    endif()
    if(LEO2_HAVE_AVX2_BACKEND)
        add_library(leopard2_backend_avx2_test_hooks OBJECT EXCLUDE_FROM_ALL
            Leopard2BackendAVX2.cpp)
        if(LEOPARD_ENABLE_GF8)
            target_sources(leopard2_backend_avx2_test_hooks PRIVATE
                Leopard2BackendAVX2Xor.cpp)
        endif()
        target_include_directories(leopard2_backend_avx2_test_hooks PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_avx2_test_hooks PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING}>)
        if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
            target_compile_definitions(
                leopard2_backend_avx2_test_hooks PRIVATE
                LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1)
        endif()
        if(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            target_compile_definitions(
                leopard2_backend_avx2_test_hooks PRIVATE
                LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_avx2_test_hooks PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_avx2_test_hooks PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
           CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
            target_compile_options(leopard2_backend_avx2_test_hooks PRIVATE
                -mavx2 -mno-avx512f)
            if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
                target_compile_options(leopard2_backend_avx2_test_hooks PRIVATE
                    -falign-functions=64)
            endif()
        elseif(MSVC)
            target_compile_options(leopard2_backend_avx2_test_hooks PRIVATE
                /arch:AVX2)
        endif()
        target_sources(leopard_test_hooks PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx2_test_hooks>)
        if(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL)
            target_sources(leopard_test_hooks PRIVATE
                $<TARGET_OBJECTS:leopard2_low_p32_b64_avx2>)
        endif()
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_HAVE_AVX2_BACKEND=1
            LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=$<BOOL:${LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL}>)
    endif()
    if(LEO2_HAVE_AVX512_BACKEND)
        add_library(leopard2_backend_avx512_test_hooks OBJECT EXCLUDE_FROM_ALL
            Leopard2BackendAVX512.cpp)
        target_include_directories(leopard2_backend_avx512_test_hooks PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_avx512_test_hooks PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HAVE_AVX2_BACKEND=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_avx512_test_hooks
                PRIVATE NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_avx512_test_hooks
                PRIVATE NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_avx512_test_hooks PRIVATE
            -mavx2 -mavx512f -mavx512bw -mavx512vl
            -mprefer-vector-width=256)
        if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
            target_compile_options(leopard2_backend_avx512_test_hooks PRIVATE
                -falign-functions=64)
        endif()
        target_sources(leopard_test_hooks PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx512_test_hooks>)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_HAVE_AVX512_BACKEND=1)
    endif()
    if(LEO2_HAVE_GFNI_BACKEND)
        add_library(leopard2_backend_gfni_test_hooks OBJECT EXCLUDE_FROM_ALL
            Leopard2BackendGFNI.cpp)
        target_include_directories(leopard2_backend_gfni_test_hooks PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_gfni_test_hooks PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HAVE_AVX2_BACKEND=1
            LEO2_HAVE_GFNI_BACKEND=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_gfni_test_hooks
                PRIVATE NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_gfni_test_hooks
                PRIVATE NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_gfni_test_hooks PRIVATE
            -mavx2 -mgfni -mno-avx512f)
        if(LEO2_FLAG_FALIGN_FUNCTIONS_64)
            target_compile_options(leopard2_backend_gfni_test_hooks PRIVATE
                -falign-functions=64)
        endif()
        target_sources(leopard_test_hooks PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_gfni_test_hooks>)
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_HAVE_GFNI_BACKEND=1)
    endif()
    if(ENABLE_OPENMP)
        if(TARGET OpenMP::OpenMP_CXX)
            target_link_libraries(leopard_test_hooks PUBLIC OpenMP::OpenMP_CXX)
        else()
            target_link_libraries(leopard_test_hooks PUBLIC
                "${OpenMP_CXX_FLAGS}")
        endif()
    endif()
    if(LEO2_BACKEND_VARIANT STREQUAL "scalar")
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_BACKEND_FORCE_SCALAR=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "ssse3")
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_BACKEND_FORCE_SSSE3=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "avx2")
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_BACKEND_FORCE_AVX2=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "avx512")
        target_compile_definitions(leopard_test_hooks PRIVATE
            LEO2_BACKEND_FORCE_AVX512=1)
    endif()

    # Structural reachability checks link the exact production archive and
    # the deliberately richer hook archive separately.  The former must elide
    # native Algorithm-5 metadata in translated P=T,N=2P codecs; the latter
    # must retain it because its codec mode can switch after construction.
    add_executable(leopard2_translated_metadata_test
        tests/leopard2/test_translated_metadata.cpp)
    target_link_libraries(leopard2_translated_metadata_test leopard)
    add_executable(leopard2_translated_metadata_hooks_test
        tests/leopard2/test_translated_metadata.cpp)
    target_compile_definitions(leopard2_translated_metadata_hooks_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1
        LEO2_EXPECT_HOOK_NATIVE_METADATA=1)
    target_link_libraries(leopard2_translated_metadata_hooks_test
        leopard_test_hooks)
    if(LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_translated_metadata_test PRIVATE
            LEO2_EXPECT_GF8=1)
        target_compile_definitions(
            leopard2_translated_metadata_hooks_test PRIVATE
            LEO2_EXPECT_GF8=1)
    else()
        target_compile_definitions(leopard2_translated_metadata_test PRIVATE
            LEO2_EXPECT_GF8=0)
        target_compile_definitions(
            leopard2_translated_metadata_hooks_test PRIVATE
            LEO2_EXPECT_GF8=0)
    endif()
    if(LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_translated_metadata_test PRIVATE
            LEO2_EXPECT_GF16=1)
        target_compile_definitions(
            leopard2_translated_metadata_hooks_test PRIVATE
            LEO2_EXPECT_GF16=1)
    else()
        target_compile_definitions(leopard2_translated_metadata_test PRIVATE
            LEO2_EXPECT_GF16=0)
        target_compile_definitions(
            leopard2_translated_metadata_hooks_test PRIVATE
            LEO2_EXPECT_GF16=0)
    endif()
    add_test(NAME leopard2_translated_metadata
        COMMAND leopard2_translated_metadata_test)
    add_test(NAME leopard2_translated_metadata_hooks
        COMMAND leopard2_translated_metadata_hooks_test)
    set_tests_properties(
        leopard2_translated_metadata leopard2_translated_metadata_hooks
        PROPERTIES ENVIRONMENT
            "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    add_executable(leopard2_field_options_test
        tests/leopard2/test_field_options.cpp)
    if(LEOPARD_ENABLE_GF8)
        target_compile_definitions(leopard2_field_options_test PRIVATE
            LEO2_EXPECT_GF8=1)
    else()
        target_compile_definitions(leopard2_field_options_test PRIVATE
            LEO2_EXPECT_GF8=0)
    endif()
    if(LEOPARD_ENABLE_GF16)
        target_compile_definitions(leopard2_field_options_test PRIVATE
            LEO2_EXPECT_GF16=1)
    else()
        target_compile_definitions(leopard2_field_options_test PRIVATE
            LEO2_EXPECT_GF16=0)
    endif()
    target_link_libraries(leopard2_field_options_test leopard)
    add_test(NAME leopard2_field_options COMMAND leopard2_field_options_test)
    set_tests_properties(leopard2_field_options PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_executable(leopard2_concurrent_initialization_test
        tests/leopard2/test_concurrent_initialization.cpp)
    target_link_libraries(leopard2_concurrent_initialization_test leopard)
    add_test(NAME leopard2_concurrent_initialization
        COMMAND leopard2_concurrent_initialization_test)
    set_tests_properties(leopard2_concurrent_initialization PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_test(NAME leopard2_release_flag_isolation
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/release_flag_isolation
            -DLEO2_C_COMPILER=${CMAKE_C_COMPILER}
            -DLEO2_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DLEO2_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
            -DLEO2_GENERATOR=${CMAKE_GENERATOR}
            -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_release_flags.cmake)
    add_test(NAME leopard2_sanitizer_classification
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_sanitizer_classification.cmake)
    add_test(NAME leopard2_target_architecture_model
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/target-architecture-model
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_target_architecture_model.cmake)

    # The installed archive keeps optional ISAs in isolated object files, but
    # historical/manual whole-translation-unit SSSE3 builds remain supported.
    # Build that path explicitly so aligned-table OOM behavior cannot regress
    # merely because the production target no longer executes it.
    if(LEO2_X86_TARGET AND LEOPARD_ENABLE_GF8 AND LEOPARD_ENABLE_GF16 AND
       (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
        CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
        CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") AND
       LEO2_FLAG_MSSSE3 AND LEO2_FLAG_MNO_AVX)
        add_library(leopard_legacy_simd_init_test_lib STATIC
            ${LIB_SOURCE_FILES}
            LeopardFF8.cpp LeopardFF8.h
            LeopardFF16.cpp LeopardFF16.h)
        target_include_directories(leopard_legacy_simd_init_test_lib PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard_legacy_simd_init_test_lib PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_PORTABLE_ISA_PRIVILEGED_FIXTURE=1)
        target_compile_options(leopard_legacy_simd_init_test_lib PRIVATE
            -mssse3 -mno-avx)
        target_link_libraries(leopard_legacy_simd_init_test_lib PUBLIC
            Threads::Threads)

        add_executable(leopard2_legacy_simd_init_failure_test
            tests/leopard2/test_legacy_simd_init_failure.cpp)
        target_compile_definitions(leopard2_legacy_simd_init_failure_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_PORTABLE_ISA_PRIVILEGED_FIXTURE=1)
        target_compile_options(leopard2_legacy_simd_init_failure_test PRIVATE
            -mssse3 -mno-avx)
        target_link_libraries(leopard2_legacy_simd_init_failure_test
            leopard_legacy_simd_init_test_lib)
        add_test(NAME leopard2_legacy_simd_init_failure_ff8
            COMMAND leopard2_legacy_simd_init_failure_test ff8)
        add_test(NAME leopard2_legacy_simd_init_failure_ff16
            COMMAND leopard2_legacy_simd_init_failure_test ff16)
    endif()

    # Static ISA auditing is meaningful only for production code generation.
    # Sanitizers may insert target-specific helper sequences (Clang ASan/UBSan
    # has emitted LAHF in an otherwise SSSE3-only object), so scanning an
    # instrumented archive can report the sanitizer runtime as a codec ISA
    # leak.  Keep this classification in the test-registration scope so it
    # cannot affect the production target/source graph audited by the Visual
    # Studio project verifier.
    if(LEO2_PORTABLE_ISA_RELEASE_AUDIT AND
       LEO2_PRODUCTION_ARCHIVE_SANITIZED)
        message(FATAL_ERROR
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires an unsanitized Release "
            "archive; run sanitizer correctness and portable-ISA release "
            "audits in separate build directories")
    endif()

    # The portable archive audit is useful for reduced-field archives too, so
    # keep it outside the dual-field-only exhaustive-test graph below.
    if(LEO2_X86_64_TARGET)
        find_program(LEO2_OBJDUMP_EXECUTABLE NAMES objdump llvm-objdump)
        find_program(LEO2_SH_EXECUTABLE NAMES sh)
        if(LEO2_PORTABLE_ISA_RELEASE_AUDIT AND
           (NOT LEO2_OBJDUMP_EXECUTABLE OR
            NOT EXISTS "${LEO2_OBJDUMP_EXECUTABLE}"))
            message(FATAL_ERROR
                "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires objdump or "
                "llvm-objdump")
        endif()
        if(LEO2_PORTABLE_ISA_RELEASE_AUDIT AND
           (NOT LEO2_SH_EXECUTABLE OR NOT EXISTS "${LEO2_SH_EXECUTABLE}"))
            message(FATAL_ERROR
                "LEO2_PORTABLE_ISA_RELEASE_AUDIT requires a POSIX sh")
        endif()
        if(LEO2_PORTABLE_ISA_RELEASE_AUDIT AND
           NOT EXISTS
               "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh")
            message(FATAL_ERROR
                "LEO2_PORTABLE_ISA_RELEASE_AUDIT checker script is missing")
        endif()
        if(LEO2_OBJDUMP_EXECUTABLE AND
           EXISTS "${LEO2_OBJDUMP_EXECUTABLE}" AND
           LEO2_SH_EXECUTABLE AND EXISTS "${LEO2_SH_EXECUTABLE}")
            set(LEO2_PORTABLE_ISA_SELF_TEST_ARGS)
            if((CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
                CMAKE_C_COMPILER_ID STREQUAL "Clang") AND NOT MSVC)
                list(APPEND LEO2_PORTABLE_ISA_SELF_TEST_ARGS
                    "${CMAKE_C_COMPILER}" "${CMAKE_AR}")
            endif()
            set(LEO2_PORTABLE_ISA_EXPECTED_CLASSES cpu_features)
            if(LEO2_HAVE_SSSE3_BACKEND)
                list(APPEND LEO2_PORTABLE_ISA_EXPECTED_CLASSES ssse3)
            endif()
            if(LEO2_HAVE_AVX2_BACKEND)
                list(APPEND LEO2_PORTABLE_ISA_EXPECTED_CLASSES avx2)
                if(LEOPARD_ENABLE_GF8)
                    list(APPEND LEO2_PORTABLE_ISA_EXPECTED_CLASSES avx2_xor)
                endif()
            endif()
            if(LEO2_HAVE_AVX512_BACKEND)
                list(APPEND LEO2_PORTABLE_ISA_EXPECTED_CLASSES avx512)
            endif()
            if(LEO2_HAVE_GFNI_BACKEND)
                list(APPEND LEO2_PORTABLE_ISA_EXPECTED_CLASSES gfni)
            endif()
            string(REPLACE ";" "," LEO2_PORTABLE_ISA_EXPECTED_CLASSES_ARG
                "${LEO2_PORTABLE_ISA_EXPECTED_CLASSES}")
            set(LEO2_PORTABLE_ISA_METADATA_ARGS)
            if(LEO2_PORTABLE_ISA_RELEASE_AUDIT)
                list(APPEND LEO2_PORTABLE_ISA_METADATA_ARGS required)
            endif()
            if(LEO2_SANITIZED_CONFIGURATION_COUNT GREATER 0)
                # Keep the checker executable and its adversarial fixtures in
                # sanitizer CTest matrices, but do not mislabel sanitizer-
                # generated instructions as production archive evidence.
                if(LEO2_MULTI_CONFIG)
                    add_test(NAME leopard2_portable_isa_checker_self_test
                        COMMAND "${LEO2_SH_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh"
                            --self-test-only
                            "${LEO2_OBJDUMP_EXECUTABLE}"
                            "${CMAKE_C_COMPILER}" "${CMAKE_AR}"
                        CONFIGURATIONS ${LEO2_SANITIZED_CONFIGURATIONS})
                else()
                    add_test(NAME leopard2_portable_isa_checker_self_test
                        COMMAND "${LEO2_SH_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh"
                            --self-test-only
                            "${LEO2_OBJDUMP_EXECUTABLE}"
                            "${CMAKE_C_COMPILER}" "${CMAKE_AR}")
                endif()
            endif()
            if(LEO2_UNSANITIZED_CONFIGURATION_COUNT GREATER 0)
                if(LEO2_MULTI_CONFIG)
                    add_test(NAME leopard2_portable_isa
                        COMMAND "${LEO2_SH_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh"
                            "${LEO2_OBJDUMP_EXECUTABLE}"
                            "$<TARGET_FILE:leopard>"
                            "${CMAKE_CURRENT_BINARY_DIR}"
                            ${LEO2_PORTABLE_ISA_SELF_TEST_ARGS}
                            "${LEO2_PORTABLE_ISA_EXPECTED_CLASSES_ARG}"
                            ${LEO2_PORTABLE_ISA_METADATA_ARGS}
                        CONFIGURATIONS ${LEO2_UNSANITIZED_CONFIGURATIONS})
                else()
                    add_test(NAME leopard2_portable_isa
                        COMMAND "${LEO2_SH_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh"
                            "${LEO2_OBJDUMP_EXECUTABLE}"
                            "$<TARGET_FILE:leopard>"
                            "${CMAKE_CURRENT_BINARY_DIR}"
                            ${LEO2_PORTABLE_ISA_SELF_TEST_ARGS}
                            "${LEO2_PORTABLE_ISA_EXPECTED_CLASSES_ARG}"
                            ${LEO2_PORTABLE_ISA_METADATA_ARGS})
                endif()
            endif()
            string(REPLACE ";" ","
                LEO2_SANITIZED_CONFIGURATIONS_ARG
                "${LEO2_SANITIZED_CONFIGURATIONS}")
            string(REPLACE ";" ","
                LEO2_UNSANITIZED_CONFIGURATIONS_ARG
                "${LEO2_UNSANITIZED_CONFIGURATIONS}")
            add_test(NAME leopard2_portable_isa_registration
                COMMAND ${CMAKE_COMMAND}
                    -DLEO2_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}
                    -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                    -DLEO2_MULTI_CONFIG=${LEO2_MULTI_CONFIG}
                    -DLEO2_SANITIZED_CONFIGURATIONS=${LEO2_SANITIZED_CONFIGURATIONS_ARG}
                    -DLEO2_UNSANITIZED_CONFIGURATIONS=${LEO2_UNSANITIZED_CONFIGURATIONS_ARG}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_portable_isa_registration.cmake)
        endif()
    endif()

    if(LEOPARD_ENABLE_GF8)
        # This target deliberately links the ordinary production archive, not
        # leopard_test_hooks.  It proves that the default-off historical
        # candidate controls real table preparation and AUTO dispatch in the
        # shipped graph, while both selector states retain oracle parity.
        add_executable(leopard2_high_direct_production_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_high_direct_production.cpp)
        target_include_directories(leopard2_high_direct_production_test PRIVATE
            tests/leopard2)
        if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_DIRECT_PRODUCTION=1)
        else()
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_DIRECT_PRODUCTION=0)
        endif()
        if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE AND
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO)
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_DIRECT_AUTO=1)
        else()
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_DIRECT_AUTO=0)
        endif()
        if(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_VECTOR=0)
        else()
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_VECTOR=1)
        endif()
        if(LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING AND
           NOT LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING=1)
        else()
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING=0)
        endif()
        if(LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING AND
           NOT LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING=1)
        else()
            target_compile_definitions(
                leopard2_high_direct_production_test PRIVATE
                LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING=0)
        endif()
        target_compile_definitions(
            leopard2_high_direct_production_test PRIVATE
            LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=$<BOOL:${LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING}>)
        target_link_libraries(leopard2_high_direct_production_test leopard)
        add_test(NAME leopard2_high_direct_production
            COMMAND leopard2_high_direct_production_test)
        set_tests_properties(leopard2_high_direct_production PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        if(LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL)
            add_executable(leopard2_low_p32_b64_terminal_test
                tests/leopard2/direct_oracle.cpp
                tests/leopard2/test_low_p32_b64_terminal.cpp)
            target_include_directories(
                leopard2_low_p32_b64_terminal_test PRIVATE tests/leopard2)
            target_compile_definitions(
                leopard2_low_p32_b64_terminal_test PRIVATE
                LEO2_ENABLE_TEST_HOOKS=1
                LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1)
            target_link_libraries(
                leopard2_low_p32_b64_terminal_test leopard_test_hooks)
            add_test(NAME leopard2_low_p32_b64_terminal
                COMMAND leopard2_low_p32_b64_terminal_test)
            set_tests_properties(leopard2_low_p32_b64_terminal PROPERTIES
                ENVIRONMENT
                    "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")
        endif()

        add_executable(leopard2_k5r5_b64_terminal_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_k5r5_b64_terminal.cpp)
        target_include_directories(leopard2_k5r5_b64_terminal_test PRIVATE
            tests/leopard2)
        target_compile_definitions(leopard2_k5r5_b64_terminal_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_k5r5_b64_terminal_test
            leopard_test_hooks)
        add_test(NAME leopard2_k5r5_b64_terminal
            COMMAND leopard2_k5r5_b64_terminal_test)
        set_tests_properties(leopard2_k5r5_b64_terminal PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        if(NOT LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            add_executable(leopard2_k8r8_b64_terminal_test
                tests/leopard2/direct_oracle.cpp
                tests/leopard2/test_k8r8_b64_terminal.cpp)
            target_include_directories(
                leopard2_k8r8_b64_terminal_test PRIVATE tests/leopard2)
            target_compile_definitions(
                leopard2_k8r8_b64_terminal_test PRIVATE
                LEO2_ENABLE_TEST_HOOKS=1)
            target_link_libraries(leopard2_k8r8_b64_terminal_test
                leopard_test_hooks)
            add_test(NAME leopard2_k8r8_b64_terminal
                COMMAND leopard2_k8r8_b64_terminal_test)
            set_tests_properties(leopard2_k8r8_b64_terminal PROPERTIES
                ENVIRONMENT
                    "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

            # Link the ordinary archive as a separate gate for the fixed,
            # portable production scratch geometry and public/batch ABI.
            add_executable(leopard2_k8r8_b64_terminal_production_test
                tests/leopard2/direct_oracle.cpp
                tests/leopard2/test_k8r8_b64_terminal_production.cpp)
            target_include_directories(
                leopard2_k8r8_b64_terminal_production_test PRIVATE
                tests/leopard2)
            target_link_libraries(
                leopard2_k8r8_b64_terminal_production_test leopard)
            add_test(NAME leopard2_k8r8_b64_terminal_production
                COMMAND leopard2_k8r8_b64_terminal_production_test)
            set_tests_properties(
                leopard2_k8r8_b64_terminal_production PROPERTIES
                ENVIRONMENT
                    "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")
        endif()

        add_executable(leopard2_k5r4_b64_terminal_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_k5r4_b64_terminal.cpp)
        target_include_directories(leopard2_k5r4_b64_terminal_test PRIVATE
            tests/leopard2)
        target_compile_definitions(leopard2_k5r4_b64_terminal_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_k5r4_b64_terminal_test
            leopard_test_hooks)
        add_test(NAME leopard2_k5r4_b64_terminal
            COMMAND leopard2_k5r4_b64_terminal_test)
        set_tests_properties(leopard2_k5r4_b64_terminal PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        add_executable(leopard2_t4_packed_terminal_family_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_t4_packed_terminal_family.cpp)
        target_include_directories(
            leopard2_t4_packed_terminal_family_test PRIVATE tests/leopard2)
        target_compile_definitions(
            leopard2_t4_packed_terminal_family_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_t4_packed_terminal_family_test
            leopard_test_hooks)
        add_test(NAME leopard2_t4_packed_terminal_family
            COMMAND leopard2_t4_packed_terminal_family_test)
        set_tests_properties(leopard2_t4_packed_terminal_family PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        # This target intentionally links the regular production archive.  In
        # particular it exercises the non-hook T=4 terminal's fixed scratch
        # geometry instead of the diagnostic EncodeLayout branch above.
        add_executable(leopard2_t4_packed_terminal_family_production_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_t4_packed_terminal_family_production.cpp)
        target_include_directories(
            leopard2_t4_packed_terminal_family_production_test PRIVATE
            tests/leopard2)
        target_link_libraries(
            leopard2_t4_packed_terminal_family_production_test leopard)
        add_test(NAME leopard2_t4_packed_terminal_family_production
            COMMAND leopard2_t4_packed_terminal_family_production_test)
        set_tests_properties(
            leopard2_t4_packed_terminal_family_production PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        add_executable(leopard2_t32_final_ifft_range_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_t32_final_ifft_range.cpp)
        target_include_directories(
            leopard2_t32_final_ifft_range_test PRIVATE tests/leopard2)
        target_compile_definitions(
            leopard2_t32_final_ifft_range_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_t32_final_ifft_range_test
            leopard_test_hooks)
        add_test(NAME leopard2_t32_final_ifft_range
            COMMAND leopard2_t32_final_ifft_range_test)
        set_tests_properties(leopard2_t32_final_ifft_range PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        add_executable(leopard2_balanced_b64_terminal_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_balanced_b64_terminal.cpp)
        target_include_directories(
            leopard2_balanced_b64_terminal_test PRIVATE tests/leopard2)
        target_compile_definitions(
            leopard2_balanced_b64_terminal_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_balanced_b64_terminal_test
            leopard_test_hooks)
        add_test(NAME leopard2_balanced_b64_terminal
            COMMAND leopard2_balanced_b64_terminal_test)
        set_tests_properties(leopard2_balanced_b64_terminal PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        # Exercise the fixed non-hook scratch geometry and ordinary archive;
        # the hook test above intentionally derives its larger diagnostic
        # workspace through EncodeLayout.
        add_executable(leopard2_balanced_b64_terminal_production_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_balanced_b64_terminal_production.cpp)
        target_include_directories(
            leopard2_balanced_b64_terminal_production_test PRIVATE
            tests/leopard2)
        target_link_libraries(
            leopard2_balanced_b64_terminal_production_test leopard)
        add_test(NAME leopard2_balanced_b64_terminal_production
            COMMAND leopard2_balanced_b64_terminal_production_test)
        set_tests_properties(
            leopard2_balanced_b64_terminal_production PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        add_executable(leopard2_t2_packed_terminal_test
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_t2_packed_terminal.cpp)
        target_include_directories(leopard2_t2_packed_terminal_test PRIVATE
            tests/leopard2)
        target_compile_definitions(leopard2_t2_packed_terminal_test PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(leopard2_t2_packed_terminal_test
            leopard_test_hooks)
        add_test(NAME leopard2_t2_packed_terminal
            COMMAND leopard2_t2_packed_terminal_test)
        set_tests_properties(leopard2_t2_packed_terminal PROPERTIES
            ENVIRONMENT
                "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

        if(LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING AND
           NOT LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
            add_executable(leopard2_k16r8_b256_terminal_test
                tests/leopard2/direct_oracle.cpp
                tests/leopard2/test_k16r8_b256_terminal.cpp)
            target_include_directories(
                leopard2_k16r8_b256_terminal_test PRIVATE tests/leopard2)
            target_compile_definitions(
                leopard2_k16r8_b256_terminal_test PRIVATE
                LEO2_ENABLE_TEST_HOOKS=1)
            target_link_libraries(leopard2_k16r8_b256_terminal_test
                leopard_test_hooks)
            add_test(NAME leopard2_k16r8_b256_terminal
                COMMAND leopard2_k16r8_b256_terminal_test)
            set_tests_properties(leopard2_k16r8_b256_terminal PROPERTIES
                ENVIRONMENT
                    "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

            add_executable(leopard2_k9r5_b256_terminal_test
                tests/leopard2/direct_oracle.cpp
                tests/leopard2/test_k9r5_b256_terminal.cpp)
            target_include_directories(
                leopard2_k9r5_b256_terminal_test PRIVATE tests/leopard2)
            target_compile_definitions(
                leopard2_k9r5_b256_terminal_test PRIVATE
                LEO2_ENABLE_TEST_HOOKS=1)
            target_link_libraries(leopard2_k9r5_b256_terminal_test
                leopard_test_hooks)
            add_test(NAME leopard2_k9r5_b256_terminal
                COMMAND leopard2_k9r5_b256_terminal_test)
            set_tests_properties(leopard2_k9r5_b256_terminal PROPERTIES
                ENVIRONMENT
                    "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")
        endif()
    endif()

    # Bounded experiment-only verifier for the AVX2 loss-5-through-8 direct
    # repair screen. It is intentionally excluded from ordinary builds and
    # CTest: the campaign launches deterministic shards and verifies every
    # selected coefficient matrix against the independent direct-RS oracle.
    # This target only exercises GF8 and must remain available to reduced
    # GF8-only configurations.
    if(LEOPARD_ENABLE_GF8 AND LEO2_HAVE_AVX2_BACKEND)
        add_executable(leopard2_small_direct_exhaustive_test EXCLUDE_FROM_ALL
            tests/leopard2/direct_oracle.cpp
            tests/leopard2/test_small_direct_exhaustive.cpp)
        target_include_directories(leopard2_small_direct_exhaustive_test
            PRIVATE tests/leopard2)
        if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE STREQUAL "0")
            set_property(
                SOURCE tests/leopard2/test_small_direct_exhaustive.cpp
                APPEND PROPERTY COMPILE_DEFINITIONS
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=${LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE}")
        endif()
        target_link_libraries(leopard2_small_direct_exhaustive_test leopard)
    endif()

    # The exhaustive internal suite intentionally cross-checks both concrete
    # field implementations in the same process.  Reduced archives instead
    # build the focused public/legacy contract test above; the default dual
    # build also runs a recursive install-and-consume matrix for both reduced
    # configurations.
    if(LEOPARD_ENABLE_GF8 AND LEOPARD_ENABLE_GF16)
    if(LEO2_BUILD_BENCHMARKS)
        add_executable(bench_leopard2_direct_encode
            bench/leopard2/direct_encode_benchmark.cpp)
        target_compile_definitions(bench_leopard2_direct_encode PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_BENCHMARK_BUILD_VARIANT="${LEO2_BACKEND_VARIANT}")
        if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
            target_compile_definitions(bench_leopard2_direct_encode PRIVATE
                LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1
                LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO=$<BOOL:${LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO}>)
        endif()
        target_link_libraries(bench_leopard2_direct_encode leopard_test_hooks)
        leopard2_enable_benchmark_source_attestation(
            bench_leopard2_direct_encode)

        # This diagnostic benchmark is built only with the test hooks that
        # expose the full padded LCH oracle; it is never installed.  Production
        # Algorithm 5 owns immutable schedules through its ordinary plan API.
        add_executable(bench_leopard2_pruned_transform
            bench/leopard2/pruned_transform_benchmark.cpp)
        target_compile_definitions(bench_leopard2_pruned_transform PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_PRUNED_SOURCE_GIT_SHA="${LEO2_LOCATOR_SOURCE_GIT_SHA}"
            LEO2_PRUNED_SOURCE_DIRTY=${LEO2_LOCATOR_SOURCE_DIRTY})
        target_link_libraries(bench_leopard2_pruned_transform
            leopard_test_hooks)

        # Same-source Algorithm 5 copy/no-copy attribution.  This deliberately
        # recompiles the ordinary JSON benchmark against the private hooks
        # archive; neither the selector nor this executable is installed or
        # linked into the production archive.
        add_executable(bench_leopard2_high_decode_copy_attribution
            bench/leopard2/benchmark.cpp)
        target_compile_definitions(
            bench_leopard2_high_decode_copy_attribution PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HIGH_DECODE_COPY_ATTRIBUTION=1)
        target_link_libraries(
            bench_leopard2_high_decode_copy_attribution leopard_test_hooks)

        # Same-source P=T,N=2P comparison between ordinary Algorithm 5 and a
        # translated Algorithm 4 plan over the unchanged legacy-high wire.
        add_executable(bench_leopard2_high_low_duality
            bench/leopard2/benchmark.cpp)
        target_compile_definitions(
            bench_leopard2_high_low_duality PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HIGH_LOW_DUALITY_ATTRIBUTION=1)
        target_link_libraries(
            bench_leopard2_high_low_duality leopard_test_hooks)

    endif()

    add_executable(leopard2_direct_oracle_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_direct_oracle.cpp)
    target_include_directories(leopard2_direct_oracle_test PRIVATE tests/leopard2)
    add_test(NAME leopard2_direct_oracle COMMAND leopard2_direct_oracle_test)

    add_executable(leopard2_backend_ops_test
        tests/leopard2/test_backend_ops.cpp)
    target_link_libraries(leopard2_backend_ops_test leopard)
    add_test(NAME leopard2_backend_ops COMMAND leopard2_backend_ops_test)
    set_tests_properties(leopard2_backend_ops PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_backend_failures_test
        tests/leopard2/test_backend_failures.cpp)
    target_compile_definitions(leopard2_backend_failures_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_backend_failures_test leopard_test_hooks)
    add_test(NAME leopard2_backend_failure_scalar_ff8_allocation
        COMMAND leopard2_backend_failures_test scalar ff8-allocation)
    add_test(NAME leopard2_backend_failure_scalar_ff16_allocation
        COMMAND leopard2_backend_failures_test scalar ff16-allocation)
    add_test(NAME leopard2_backend_failure_scalar_kat
        COMMAND leopard2_backend_failures_test scalar kat)
    add_test(NAME leopard2_backend_failure_ssse3_ff8_allocation
        COMMAND leopard2_backend_failures_test ssse3 ff8-allocation)
    add_test(NAME leopard2_backend_failure_ssse3_ff16_allocation
        COMMAND leopard2_backend_failures_test ssse3 ff16-allocation)
    add_test(NAME leopard2_backend_failure_ssse3_kat
        COMMAND leopard2_backend_failures_test ssse3 kat)
    add_test(NAME leopard2_backend_failure_avx2_ff8_allocation
        COMMAND leopard2_backend_failures_test avx2 ff8-allocation)
    add_test(NAME leopard2_backend_failure_avx2_ff16_allocation
        COMMAND leopard2_backend_failures_test avx2 ff16-allocation)
    add_test(NAME leopard2_backend_failure_avx2_kat
        COMMAND leopard2_backend_failures_test avx2 kat)
    if(LEO2_HAVE_AVX512_BACKEND)
        add_test(NAME leopard2_backend_failure_avx512_ff8_allocation
            COMMAND leopard2_backend_failures_test avx512 ff8-allocation)
        add_test(NAME leopard2_backend_failure_avx512_ff16_allocation
            COMMAND leopard2_backend_failures_test avx512 ff16-allocation)
        add_test(NAME leopard2_backend_failure_avx512_kat
            COMMAND leopard2_backend_failures_test avx512 kat)
        add_test(NAME leopard2_backend_auto_avx512_kat_fallback
            COMMAND leopard2_backend_failures_test
                auto-avx512-kat-fallback)
    endif()
    if(LEO2_HAVE_GFNI_BACKEND)
        add_test(NAME leopard2_backend_failure_gfni_ff8_allocation
            COMMAND leopard2_backend_failures_test gfni ff8-allocation)
        add_test(NAME leopard2_backend_failure_gfni_ff16_allocation
            COMMAND leopard2_backend_failures_test gfni ff16-allocation)
        add_test(NAME leopard2_backend_failure_gfni_kat
            COMMAND leopard2_backend_failures_test gfni kat)
    endif()
    set_tests_properties(
        leopard2_backend_failure_scalar_ff8_allocation
        leopard2_backend_failure_scalar_ff16_allocation
        leopard2_backend_failure_scalar_kat
        leopard2_backend_failure_ssse3_ff8_allocation
        leopard2_backend_failure_ssse3_ff16_allocation
        leopard2_backend_failure_ssse3_kat
        leopard2_backend_failure_avx2_ff8_allocation
        leopard2_backend_failure_avx2_ff16_allocation
        leopard2_backend_failure_avx2_kat
        PROPERTIES ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")
    if(LEO2_HAVE_AVX512_BACKEND)
        set_tests_properties(
            leopard2_backend_failure_avx512_ff8_allocation
            leopard2_backend_failure_avx512_ff16_allocation
            leopard2_backend_failure_avx512_kat
            leopard2_backend_auto_avx512_kat_fallback
            PROPERTIES ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")
    endif()
    if(LEO2_HAVE_GFNI_BACKEND)
        set_tests_properties(
            leopard2_backend_failure_gfni_ff8_allocation
            leopard2_backend_failure_gfni_ff16_allocation
            leopard2_backend_failure_gfni_kat
            PROPERTIES ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")
    endif()

    add_executable(leopard2_context_backends_test
        tests/leopard2/test_context_backends.cpp)
    target_compile_definitions(leopard2_context_backends_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_context_backends_test leopard_test_hooks)
    add_test(NAME leopard2_context_backends
        COMMAND leopard2_context_backends_test)
    set_tests_properties(leopard2_context_backends PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=4")

    add_executable(leopard2_auto_encode_backend_test
        tests/leopard2/test_auto_encode_backend.cpp)
    target_compile_definitions(leopard2_auto_encode_backend_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_auto_encode_backend_test leopard_test_hooks)
    add_test(NAME leopard2_auto_encode_backend
        COMMAND leopard2_auto_encode_backend_test)
    set_tests_properties(leopard2_auto_encode_backend PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    if(LEOPARD_ENABLE_GF16)
        # This target deliberately links the ordinary production archive, not
        # leopard_test_hooks.  Sparse AUTO scratch and dispatch must therefore
        # remain correct when diagnostic hook builds do not conservatively
        # reserve schedules for every codec/backend.
        add_executable(leopard2_sparse_encode_production_test
            tests/leopard2/test_sparse_encode_production.cpp)
        target_link_libraries(leopard2_sparse_encode_production_test leopard)
        add_test(NAME leopard2_sparse_encode_production
            COMMAND leopard2_sparse_encode_production_test)
        set_tests_properties(leopard2_sparse_encode_production PROPERTIES
            ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")
    endif()

    add_executable(leopard2_r1_xor_test
        tests/leopard2/test_r1_xor.cpp)
    target_link_libraries(leopard2_r1_xor_test leopard)
    add_test(NAME leopard2_r1_xor COMMAND leopard2_r1_xor_test)
    set_tests_properties(leopard2_r1_xor PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_executable(leopard2_direct_encode_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_direct_encode.cpp)
    target_compile_definitions(leopard2_direct_encode_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
        target_compile_definitions(leopard2_direct_encode_test PRIVATE
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=1
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO=$<BOOL:${LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO}>)
    endif()
    if(LEO2_BACKEND_VARIANT STREQUAL "auto" AND
       LEO2_HAVE_AVX512_BACKEND)
        target_compile_definitions(leopard2_direct_encode_test PRIVATE
            LEO2_TEST_AUTO_MAY_USE_AVX512_ENCODE=1)
    endif()
    target_include_directories(leopard2_direct_encode_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_direct_encode_test leopard_test_hooks)
    add_test(NAME leopard2_direct_encode COMMAND leopard2_direct_encode_test)
    set_tests_properties(leopard2_direct_encode PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_active_lch_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_active_lch.cpp)
    target_compile_definitions(leopard2_active_lch_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_include_directories(leopard2_active_lch_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_active_lch_test leopard_test_hooks)
    add_test(NAME leopard2_active_lch COMMAND leopard2_active_lch_test)
    set_tests_properties(leopard2_active_lch PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_api_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_api.cpp)
    target_compile_definitions(leopard2_api_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN)
        target_compile_definitions(leopard2_api_test PRIVATE
            LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=1)
    endif()
    if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE STREQUAL "0")
        set_property(SOURCE tests/leopard2/test_api.cpp
            APPEND PROPERTY COMPILE_DEFINITIONS
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=${LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE}")
    endif()
    target_include_directories(leopard2_api_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_api_test leopard_test_hooks)
    add_test(NAME leopard2_api COMMAND leopard2_api_test)

    add_executable(leopard2_public_api_contract_test
        tests/leopard2/test_public_api_contract.cpp)
    target_compile_definitions(leopard2_public_api_contract_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_public_api_contract_test leopard_test_hooks)
    add_test(NAME leopard2_public_api_contract
        COMMAND leopard2_public_api_contract_test)

    # The production-linked executable observes the public one-shot
    # allocation path.  The hook-linked companion also proves that forced
    # diagnostic decoder modes retain canonical plan construction.
    add_executable(leopard2_one_shot_no_loss_test
        tests/leopard2/test_one_shot_no_loss.cpp)
    target_link_libraries(leopard2_one_shot_no_loss_test leopard)
    add_executable(leopard2_one_shot_no_loss_hooks_test
        tests/leopard2/test_one_shot_no_loss.cpp)
    target_compile_definitions(
        leopard2_one_shot_no_loss_hooks_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1
        LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD=1)
    target_link_libraries(
        leopard2_one_shot_no_loss_hooks_test leopard_test_hooks)
    if(LEOPARD_ENABLE_GF8)
        target_compile_definitions(
            leopard2_one_shot_no_loss_test PRIVATE
            LEO2_TEST_EXPECT_GF8=1)
        target_compile_definitions(
            leopard2_one_shot_no_loss_hooks_test PRIVATE
            LEO2_TEST_EXPECT_GF8=1)
    endif()
    if(LEOPARD_ENABLE_GF16)
        target_compile_definitions(
            leopard2_one_shot_no_loss_test PRIVATE
            LEO2_TEST_EXPECT_GF16=1)
        target_compile_definitions(
            leopard2_one_shot_no_loss_hooks_test PRIVATE
            LEO2_TEST_EXPECT_GF16=1)
    endif()
    add_test(NAME leopard2_one_shot_no_loss
        COMMAND leopard2_one_shot_no_loss_test)
    add_test(NAME leopard2_one_shot_no_loss_hooks
        COMMAND leopard2_one_shot_no_loss_hooks_test)
    set_tests_properties(
        leopard2_one_shot_no_loss
        leopard2_one_shot_no_loss_hooks PROPERTIES
        ENVIRONMENT
            "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    add_executable(leopard2_batch_aliasing_test
        tests/leopard2/test_batch_aliasing.cpp)
    if(LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING AND
       NOT LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
        target_compile_definitions(leopard2_batch_aliasing_test PRIVATE
            LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING=1)
    else()
        target_compile_definitions(leopard2_batch_aliasing_test PRIVATE
            LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING=0)
    endif()
    target_link_libraries(leopard2_batch_aliasing_test leopard)
    add_test(NAME leopard2_batch_aliasing
        COMMAND leopard2_batch_aliasing_test)
    set_tests_properties(leopard2_batch_aliasing PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND ENABLE_OPENMP)
        add_executable(leopard2_initialization_threads_test
            tests/leopard2/test_initialization_threads.cpp)
        target_link_libraries(leopard2_initialization_threads_test leopard)
        add_test(NAME leopard2_initialization_threads_legacy
            COMMAND leopard2_initialization_threads_test legacy)
        add_test(NAME leopard2_initialization_threads_explicit
            COMMAND leopard2_initialization_threads_test explicit)
        add_test(NAME leopard2_initialization_threads_default
            COMMAND leopard2_initialization_threads_test default)
        add_test(NAME leopard2_initialization_threads_gf16_high_codec
            COMMAND leopard2_initialization_threads_test gf16-high-codec)
        add_test(NAME leopard2_initialization_threads_gf16_high_plan
            COMMAND leopard2_initialization_threads_test gf16-high-plan)
        add_test(NAME leopard2_initialization_threads_gf16_low_plan
            COMMAND leopard2_initialization_threads_test gf16-low-plan)
    endif()

    add_executable(leopard2_random_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_random.cpp)
    target_include_directories(leopard2_random_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_random_test leopard)
    add_test(NAME leopard2_random COMMAND leopard2_random_test)

    add_executable(leopard2_locator_test tests/leopard2/test_locator.cpp)
    target_link_libraries(leopard2_locator_test leopard)
    add_test(NAME leopard2_locator COMMAND leopard2_locator_test)

    add_executable(leopard2_decode_plan_schedule_test
        tests/leopard2/test_decode_plan_schedule.cpp)
    target_link_libraries(leopard2_decode_plan_schedule_test leopard)
    add_test(NAME leopard2_decode_plan_schedule
        COMMAND leopard2_decode_plan_schedule_test)
    set_tests_properties(leopard2_decode_plan_schedule PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_dense_plan_policy_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_dense_plan_policy.cpp)
    target_include_directories(leopard2_dense_plan_policy_test
        PRIVATE tests/leopard2)
    if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE STREQUAL "0")
        set_property(SOURCE tests/leopard2/test_dense_plan_policy.cpp
            APPEND PROPERTY COMPILE_DEFINITIONS
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=${LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE}")
    endif()
    target_link_libraries(leopard2_dense_plan_policy_test leopard)
    add_test(NAME leopard2_dense_plan_policy
        COMMAND leopard2_dense_plan_policy_test)
    set_tests_properties(leopard2_dense_plan_policy PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    add_executable(leopard2_gf16_tails_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_gf16_tails.cpp)
    target_include_directories(leopard2_gf16_tails_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_gf16_tails_test leopard)
    add_test(NAME leopard2_gf16_tails COMMAND leopard2_gf16_tails_test)

    add_executable(leopard2_gf16_padded_odd_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_gf16_padded_odd.cpp)
    target_include_directories(leopard2_gf16_padded_odd_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_gf16_padded_odd_test leopard)
    add_test(NAME leopard2_gf16_padded_odd COMMAND leopard2_gf16_padded_odd_test)
    set_tests_properties(leopard2_gf16_padded_odd PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_gf16_legacy_encoder_matrix_test
        tests/leopard2/test_encoder_gf16_legacy_matrix.cpp)
    target_link_libraries(leopard2_gf16_legacy_encoder_matrix_test leopard)
    add_test(NAME leopard2_gf16_legacy_encoder_matrix
        COMMAND leopard2_gf16_legacy_encoder_matrix_test)
    set_tests_properties(leopard2_gf16_legacy_encoder_matrix PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_low_gf16_direct_rows_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_low_gf16_direct_rows.cpp)
    target_include_directories(leopard2_low_gf16_direct_rows_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_low_gf16_direct_rows_test leopard)
    add_test(NAME leopard2_low_gf16_direct_rows
        COMMAND leopard2_low_gf16_direct_rows_test)
    set_tests_properties(leopard2_low_gf16_direct_rows PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_decode_high_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_decode_high_acceptance.cpp)
    target_include_directories(leopard2_decode_high_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_decode_high_acceptance_test leopard)
    add_test(NAME leopard2_decode_high_acceptance
        COMMAND leopard2_decode_high_acceptance_test)
    set_tests_properties(leopard2_decode_high_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_high_pruned_legacy_test
        tests/leopard2/test_high_pruned_legacy.cpp)
    target_link_libraries(leopard2_high_pruned_legacy_test leopard)
    add_test(NAME leopard2_high_pruned_legacy
        COMMAND leopard2_high_pruned_legacy_test)
    set_tests_properties(leopard2_high_pruned_legacy PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    # Preserve the production-library compatibility gate above, and execute
    # the same matrix once against the private hook archive so the
    # parent-wide GF16 receive-stage counters are not compiled out.  The hook
    # archive is never installed and cannot contaminate the production target.
    add_executable(leopard2_high_pruned_stage_test
        tests/leopard2/test_high_pruned_legacy.cpp)
    target_compile_definitions(leopard2_high_pruned_stage_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1
        LEO2_REQUIRE_HIGH_PRUNED_STAGE_HOOKS=1)
    target_link_libraries(leopard2_high_pruned_stage_test leopard_test_hooks)
    add_test(NAME leopard2_high_pruned_stage
        COMMAND leopard2_high_pruned_stage_test)
    set_tests_properties(leopard2_high_pruned_stage PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    add_executable(leopard2_decode_low_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_decode_low_acceptance.cpp)
    target_include_directories(leopard2_decode_low_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_decode_low_acceptance_test leopard)
    add_test(NAME leopard2_decode_low_acceptance
        COMMAND leopard2_decode_low_acceptance_test)
    set_tests_properties(leopard2_decode_low_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    # Correctness and control-path proof for the production IT2026 Algorithm 4
    # coset translation and its retained test-only Algorithm 5 comparison.
    add_executable(leopard2_high_low_duality_test
        tests/leopard2/test_high_low_duality.cpp)
    target_compile_definitions(leopard2_high_low_duality_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE STREQUAL "0")
        set_property(SOURCE tests/leopard2/test_high_low_duality.cpp
            APPEND PROPERTY COMPILE_DEFINITIONS
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=${LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE}")
    endif()
    target_link_libraries(leopard2_high_low_duality_test leopard_test_hooks)
    add_test(NAME leopard2_high_low_duality
        COMMAND leopard2_high_low_duality_test)
    set_tests_properties(leopard2_high_low_duality PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")

    add_executable(leopard2_arbitrary_counts_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_arbitrary_counts_acceptance.cpp)
    target_include_directories(leopard2_arbitrary_counts_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_arbitrary_counts_acceptance_test leopard)
    add_test(NAME leopard2_arbitrary_counts_acceptance
        COMMAND leopard2_arbitrary_counts_acceptance_test)
    set_tests_properties(leopard2_arbitrary_counts_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_max_counts_test tests/leopard2/test_max_counts.cpp)
    target_link_libraries(leopard2_max_counts_test leopard)
    add_test(NAME leopard2_max_counts COMMAND leopard2_max_counts_test)
    set_tests_properties(leopard2_max_counts PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_encode_concurrency_test
        tests/leopard2/test_encode_concurrency.cpp)
    target_link_libraries(leopard2_encode_concurrency_test leopard)
    add_test(NAME leopard2_encode_concurrency
        COMMAND leopard2_encode_concurrency_test)
    set_tests_properties(leopard2_encode_concurrency PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_codec_options_abi_test
        tests/leopard2/test_codec_options_abi.c)
    set_target_properties(leopard2_codec_options_abi_test PROPERTIES
        LINKER_LANGUAGE CXX)
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
       CMAKE_C_COMPILER_ID STREQUAL "Clang" OR
       CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
        target_compile_options(leopard2_codec_options_abi_test PRIVATE
            -fshort-enums)
    endif()
    target_link_libraries(leopard2_codec_options_abi_test leopard)
    add_test(NAME leopard2_codec_options_abi
        COMMAND leopard2_codec_options_abi_test)

    add_executable(leopard2_direct_repair_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/direct_repair.cpp
        tests/leopard2/test_direct_repair.cpp)
    target_include_directories(leopard2_direct_repair_test PRIVATE tests/leopard2)
    add_test(NAME leopard2_direct_repair COMMAND leopard2_direct_repair_test)

    add_executable(leopard2_legacy_golden_test
        tests/leopard2/test_legacy_golden.cpp)
    target_link_libraries(leopard2_legacy_golden_test leopard)
    add_test(NAME leopard2_legacy_golden COMMAND leopard2_legacy_golden_test)

    add_executable(leopard2_boundaries_test
        tests/leopard2/test_boundaries.cpp)
    target_link_libraries(leopard2_boundaries_test leopard)
    add_test(NAME leopard2_boundaries COMMAND leopard2_boundaries_test)

    add_executable(leopard2_transform_differential_test
        tests/leopard2/test_transform_differential.cpp
        tests/leopard2/direct_oracle.cpp)
    target_compile_definitions(leopard2_transform_differential_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_transform_differential_test
        leopard_test_hooks)
    add_test(NAME leopard2_transform_differential
        COMMAND leopard2_transform_differential_test)

    add_executable(leopard2_pruned_transform_test
        tests/leopard2/test_pruned_transform.cpp
        tests/leopard2/direct_oracle.cpp)
    target_compile_definitions(leopard2_pruned_transform_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_include_directories(leopard2_pruned_transform_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_pruned_transform_test leopard_test_hooks)
    add_test(NAME leopard2_pruned_transform
        COMMAND leopard2_pruned_transform_test)
    set_tests_properties(leopard2_pruned_transform PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_executable(leopard2_fuzz_smoke
        tests/leopard2/fuzz_api.cpp
        tests/leopard2/fuzz_replay.cpp)
    target_compile_definitions(leopard2_fuzz_smoke PRIVATE
        LEO2_FUZZ_REPLAY_ROLE="api")
    # The deterministic campaign driver attests sanitizer state from the
    # linked core (and may invoke a deliberate leak canary).  Keep those
    # diagnostics in the private hook archive, never in production leopard.
    target_link_libraries(leopard2_fuzz_smoke leopard_test_hooks)
    add_test(NAME leopard2_fuzz_smoke COMMAND leopard2_fuzz_smoke)

    add_executable(leopard2_pruned_fuzz_smoke
        tests/leopard2/fuzz_pruned_transform.cpp
        tests/leopard2/fuzz_replay.cpp)
    target_compile_definitions(leopard2_pruned_fuzz_smoke PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1
        LEO2_FUZZ_REPLAY_ROLE="pruned")
    target_link_libraries(leopard2_pruned_fuzz_smoke leopard_test_hooks)
    add_test(NAME leopard2_pruned_fuzz_smoke
        COMMAND leopard2_pruned_fuzz_smoke 0xc1f022d17 256)
    set_tests_properties(leopard2_pruned_fuzz_smoke PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")

    add_test(
        NAME leopard2_cuda_optional
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/cuda-optional-test
            -DLEO2_GENERATOR=${CMAKE_GENERATOR}
            -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
            -DLEO2_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}
            -DLEO2_STATIC_LIBRARY_PREFIX=${CMAKE_STATIC_LIBRARY_PREFIX}
            -DLEO2_STATIC_LIBRARY_SUFFIX=${CMAKE_STATIC_LIBRARY_SUFFIX}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_cuda_optional.cmake)

    # GNU ar/nm command-line syntax is used only for a Linux-hosted proof that
    # disabled field objects and symbols are absent.  Other tool families still
    # run the complete reduced build/install/consumer matrix below.
    set(LEO2_ARCHIVE_SYMBOL_CHECKS OFF)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT MSVC AND
       CMAKE_AR AND CMAKE_NM)
        set(LEO2_ARCHIVE_SYMBOL_CHECKS ON)
    endif()
    add_test(
        NAME leopard2_field_options_matrix
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/field-options-test
            -DLEO2_GENERATOR=${CMAKE_GENERATOR}
            -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
            -DLEO2_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}
            -DLEO2_ARCHIVE_SYMBOL_CHECKS=${LEO2_ARCHIVE_SYMBOL_CHECKS}
            -DLEO2_EXPECT_SMALL_DIRECT_EXHAUSTIVE_TARGET=${LEO2_HAVE_AVX2_BACKEND}
            -DLEO2_AR_COMMAND=${CMAKE_AR}
            -DLEO2_NM_COMMAND=${CMAKE_NM}
            -DLEO2_STATIC_LIBRARY_PREFIX=${CMAKE_STATIC_LIBRARY_PREFIX}
            -DLEO2_STATIC_LIBRARY_SUFFIX=${CMAKE_STATIC_LIBRARY_SUFFIX}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_field_options.cmake)

    add_test(
        NAME leopard2_test_hook_isolation
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/test-hook-isolation
            -DLEO2_GENERATOR=${CMAKE_GENERATOR}
            -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
            -DLEO2_ENABLE_OPENMP=${ENABLE_OPENMP}
            -DLEO2_ARCHIVE_SYMBOL_CHECKS=${LEO2_ARCHIVE_SYMBOL_CHECKS}
            -DLEO2_NM_COMMAND=${CMAKE_NM}
            -DLEO2_STATIC_LIBRARY_PREFIX=${CMAKE_STATIC_LIBRARY_PREFIX}
            -DLEO2_STATIC_LIBRARY_SUFFIX=${CMAKE_STATIC_LIBRARY_SUFFIX}
            -DLEO2_EXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_hook_isolation.cmake)

    if(LEO2_BUILD_BENCHMARKS)
        add_test(
            NAME leopard2_benchmark_attestation_refresh
            COMMAND ${CMAKE_COMMAND}
                -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/benchmark-attestation-refresh
                -DLEO2_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DLEO2_GENERATOR=${CMAKE_GENERATOR}
                -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
                -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
                -DLEO2_EXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_benchmark_attestation_refresh.cmake)
        add_test(
            NAME leopard2_sparse_encode_benchmark_registration
            COMMAND ${CMAKE_COMMAND}
                -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/sparse-encode-registration
                -DLEO2_C_COMPILER=${CMAKE_C_COMPILER}
                -DLEO2_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DLEO2_GENERATOR=${CMAKE_GENERATOR}
                -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
                -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
                -DLEO2_ARCHIVE_SYMBOL_CHECKS=${LEO2_ARCHIVE_SYMBOL_CHECKS}
                -DLEO2_NM_COMMAND=${CMAKE_NM}
                -DLEO2_STATIC_LIBRARY_PREFIX=${CMAKE_STATIC_LIBRARY_PREFIX}
                -DLEO2_STATIC_LIBRARY_SUFFIX=${CMAKE_STATIC_LIBRARY_SUFFIX}
                -DLEO2_EXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_sparse_encode_benchmark_registration.cmake)
    endif()

    # The registered Python tools use postponed annotations, PEP 585
    # collections, PEP 604 unions, and Path.is_relative_to().  Do not
    # accidentally register them with Python 2 (possible with old
    # FindPythonInterp defaults) or an older Python 3 that cannot parse them.
    set(LEO2_PYTHON_EXECUTABLE "")
    if(CMAKE_VERSION VERSION_LESS 3.12)
        find_package(PythonInterp 3.10 QUIET)
        if(PYTHONINTERP_FOUND)
            set(LEO2_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
        endif()
    else()
        find_package(Python3 3.10 COMPONENTS Interpreter QUIET)
        if(Python3_Interpreter_FOUND)
            set(LEO2_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
        endif()
    endif()
    if(LEO2_PYTHON_EXECUTABLE)
        # This is a regression for the Linux-only authoritative benchmark
        # evidence helper, which deliberately relies on POSIX descriptor
        # inheritance and canonical /usr/bin tools.  Keep it out of ordinary
        # codec-only builds and non-FHS/non-Linux CTest inventories.
        if(LEO2_BUILD_BENCHMARKS AND
           CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux" AND
           EXISTS "/usr/bin/cmake" AND EXISTS "/usr/bin/cc" AND
           EXISTS "/usr/bin/c++")
            add_test(
                NAME leopard2_build_provenance_compiler_replay
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_build_provenance_test.py
                    ${CMAKE_CURRENT_SOURCE_DIR})
            # This test deliberately retains an inotify history for the
            # checkout's Git metadata.  Other CTests exercise Git and may
            # legitimately create index.lock, so concurrent execution would
            # look exactly like the mutation this test must reject.
            set_tests_properties(
                leopard2_build_provenance_compiler_replay PROPERTIES
                RUN_SERIAL TRUE
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                TIMEOUT 300)
        endif()
        if(LEO2_BUILD_BENCHMARKS)
            add_test(
                NAME leopard2_locator_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_locator_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2_locator>)
            set_tests_properties(leopard2_locator_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
            if(TARGET bench_leopard2_allk)
                add_test(
                    NAME leopard2_benchmark_json_regression
                    COMMAND ${LEO2_PYTHON_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_benchmark_json_test.py
                        $<TARGET_FILE:bench_leopard2>
                        $<TARGET_FILE:bench_leopard2_allk>)
            else()
                add_test(
                    NAME leopard2_benchmark_json_regression
                    COMMAND ${LEO2_PYTHON_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_benchmark_json_test.py
                        $<TARGET_FILE:bench_leopard2>)
            endif()
            set_tests_properties(leopard2_benchmark_json_regression PROPERTIES
                ENVIRONMENT
                    "LEO2_EXPECTED_SOURCE_ATTESTATION_HEADER=${LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER};OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
            add_test(
                NAME leopard2_pruned_transform_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_pruned_transform_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2_pruned_transform>)
            set_tests_properties(
                leopard2_pruned_transform_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
            add_test(
                NAME leopard2_sparse_encode_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_sparse_encode_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2_sparse_encode>)
            set_tests_properties(
                leopard2_sparse_encode_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
        endif()
        add_test(
            NAME leopard2_cost_model_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/cost_model.py
                self-test)
        # The lab runner's self-test exercises Linux process affinity and
        # /proc-based runtime attestation.  The strict sanitizer campaign also
        # attests ELF instrumentation with nm and POSIX alias semantics.  Keep
        # both host-tool tests out of non-Linux CTest inventories; portable
        # codec and Python tests remain registered independently.
        if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
            add_test(
                NAME leopard2_allk_gap_identity_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/main_compare/test_run_allk_gap.py)
            add_test(
                NAME leopard2_allk_gap_identity_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/main_compare/test_run_allk_gap.py)
            add_test(
                NAME leopard2_t8_ragged_runner_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/gf8_high_encode/test_run_t8_two_block_abba.py)
            add_test(
                NAME leopard2_t8_ragged_runner_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/gf8_high_encode/test_run_t8_two_block_abba.py)
            add_test(
                NAME leopard2_affinity_supervisor_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_affinity_supervisor.py
                    self-test)
            # The supervisor's fake transaction/recovery tests deliberately use
            # explicit exceptions rather than assert statements.  Keep an
            # optimized-Python registration so that invariant cannot regress.
            add_test(
                NAME leopard2_affinity_supervisor_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_affinity_supervisor.py
                    self-test)
            add_test(
                NAME leopard2_lab_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_lab.py
                    self-test)
            add_test(
                NAME leopard2_fuzz_campaign_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_fuzz_campaign.py
                    self-test)
            # Both harnesses deliberately create short-lived and detached
            # process trees while scanning global /proc state.  Running them
            # beside the other process-heavy Python CTests can let a fixture
            # complete before it receives one real observation.  Serialize
            # the two tests without weakening their fail-closed runtime gates.
            set_tests_properties(
                leopard2_lab_self_test
                PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            set_tests_properties(
                leopard2_fuzz_campaign_self_test
                PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            # The retained GF16 neighbor-evidence replay imports the Linux
            # affinity/lease support used by its fail-closed runner.
            add_test(
                NAME leopard2_gf16_neighbor_evidence_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/gf16_high_encode/test_run_neighbor_abba.py)
            # The direct-encode crossover runner uses fcntl leases and Linux
            # affinity/procfs checks even in self-test mode.
            add_test(
                NAME leopard2_direct_encode_crossover_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_direct_encode_crossover.py
                    self-test)
            add_test(
                NAME leopard2_direct_encode_crossover_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_direct_encode_crossover.py
                    self-test)
            add_test(
                NAME leopard2_small_direct_abba_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/test_run_small_direct_abba.py)
            set_tests_properties(
                leopard2_small_direct_abba_self_test PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            add_test(
                NAME leopard2_small_direct_abba_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/test_run_small_direct_abba.py)
            set_tests_properties(
                leopard2_small_direct_abba_optimized_self_test PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            add_test(
                NAME leopard2_equal_rounded_abba_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/run_equal_rounded_abba.py
                    self-test)
            add_test(
                NAME leopard2_equal_rounded_abba_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/run_equal_rounded_abba.py
                    self-test)
            set_tests_properties(
                leopard2_equal_rounded_abba_self_test
                leopard2_equal_rounded_abba_optimized_self_test PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            add_test(
                NAME leopard2_small_direct_exhaustive_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/test_run_small_direct_exhaustive.py)
            set_tests_properties(
                leopard2_small_direct_exhaustive_self_test PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
            add_test(
                NAME leopard2_small_direct_exhaustive_optimized_self_test
                COMMAND ${LEO2_PYTHON_EXECUTABLE} -O
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/direct_repair/test_run_small_direct_exhaustive.py)
            set_tests_properties(
                leopard2_small_direct_exhaustive_optimized_self_test
                PROPERTIES
                ENVIRONMENT
                    "PYTHONDONTWRITEBYTECODE=1;PYTHONWARNINGS=error::ResourceWarning"
                RUN_SERIAL TRUE
                TIMEOUT 300)
        endif()
        add_test(
            NAME leopard2_benchmark_matrix_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_benchmark_matrix.py
                self-test)
        add_test(
            NAME leopard2_sparse_encode_crossover_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_sparse_encode_crossover.py
                self-test)

        add_test(
            NAME leopard2_r1_xor_crossover_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_r1_xor_crossover.py
                self-test)
        add_test(
            NAME leopard2_external_comparison_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_external_comparison.py
                self-test)
        # Pure policy/evidence validation only: this test neither finds nor
        # links ISA-L.  The actual adapter lives in a standalone default-off
        # experiment CMake project.
        add_test(
            NAME leopard2_isal_comparison_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_isal_compare.py
                self-test)
        # Pure policy/evidence validation only: this test never discovers,
        # compiles, links, or loads Jerasure/GF-Complete.  The real adapter is
        # an explicitly configured standalone experiment project.
        add_test(
            NAME leopard2_jerasure_comparison_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_jerasure_compare.py
                self-test)
        add_test(
            NAME leopard2_jerasure_default_optionality_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_jerasure_compare.py
                optionality-test)
        add_test(
            NAME leopard2_operation_counts_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_operation_counts.py
                self-test)
        if(LEOPARD_ENABLE_GF8)
            add_test(
                NAME leopard2_decode_scratch_crosscheck
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/operation_counts/test_decode_scratch_crosscheck.py
                    --probe $<TARGET_FILE:leopard2_decode_scratch_probe>)
        endif()
        add_test(
            NAME leopard2_tiled_high_runner_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/decoder_dispatch/run_tiled_high_abba.py
                --self-test)
        add_test(
            NAME leopard2_high_decode_copy_contract_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/high_decode_copy/benchmark_contract.py
                self-test)
        add_test(
            NAME leopard2_high_decode_copy_runner_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/high_decode_copy/run_abba.py
                self-test)
        add_test(
            NAME leopard2_high_decode_copy_composite_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/high_decode_copy/compose.py
                self-test)
        if(TARGET bench_leopard2_high_decode_copy_attribution)
            add_test(
                NAME leopard2_high_decode_copy_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/high_decode_copy/benchmark_contract.py
                    smoke
                    --benchmark
                    $<TARGET_FILE:bench_leopard2_high_decode_copy_attribution>)
            set_tests_properties(
                leopard2_high_decode_copy_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")
            # The build-identity checker authenticates CMake's Unix Makefiles
            # link.txt recipes.  Ninja expresses the same commands only through
            # its generated rule graph, so do not register a test that cannot
            # establish the recipe provenance it promises.
            if(CMAKE_EXPORT_COMPILE_COMMANDS AND
               CMAKE_BUILD_TYPE STREQUAL "Release" AND
               CMAKE_GENERATOR STREQUAL "Unix Makefiles")
                add_test(
                    NAME leopard2_high_decode_copy_build_identity
                    COMMAND ${LEO2_PYTHON_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/high_decode_copy/run_abba.py
                        build-smoke
                        --source-root ${CMAKE_CURRENT_SOURCE_DIR}
                        --binary
                        $<TARGET_FILE:bench_leopard2_high_decode_copy_attribution>
                        --hook-archive $<TARGET_FILE:leopard_test_hooks>)
            endif()
        endif()
        add_test(
            NAME leopard2_tiled_high_analyzer_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/decoder_dispatch/analyze_tiled_high.py
                --self-test)
        add_test(
            NAME leopard2_balanced_forced_runner_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/decoder_dispatch/run_balanced_abba.py
                --self-test)
        add_test(
            NAME leopard2_balanced_forced_analyzer_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/decoder_dispatch/analyze_balanced.py
                --self-test)
        add_test(
            NAME leopard2_balanced_promotion_plan_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py
                self-test)
        add_test(
            NAME leopard2_visual_studio_project_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/proj/test_leopard_vcxproj.py)
        add_test(
            NAME leopard2_canonical_library_docs_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/proj/test_canonical_library_docs.py)
    endif()
    endif()
endif()

if(LEO2_BUILD_FUZZERS)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
       NOT CMAKE_C_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR
            "LEO2_BUILD_FUZZERS requires upstream Clang C/C++ libFuzzer support")
    endif()

    # Build an independently instrumented copy of the library into the fuzz
    # executable.  Applying sanitizer flags only to fuzz_api.cpp would leave
    # the codec itself invisible to ASan, UBSan, and coverage-guided mutation;
    # applying them to leopard would in turn break ordinary targets that
    # link the uninstrumented sanitizer runtime.
    set(LEO2_FUZZ_SANITIZER_FLAGS
        -fsanitize=fuzzer,address,undefined
        -fno-omit-frame-pointer)
    add_executable(leopard2_api_fuzzer
        tests/leopard2/fuzz_api.cpp
        ${LIB_SOURCE_FILES})
    target_include_directories(leopard2_api_fuzzer PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_options(leopard2_api_fuzzer PRIVATE
        ${LEO2_FUZZ_SANITIZER_FLAGS})
    # The backend fuzz object members are shared with the pruned-plan fuzzer,
    # whose full-transform oracle requires the test-hook ABI.  Compile the API
    # fuzzer's private core copy with the same hooks so those instrumented
    # members have one consistent, fully sanitized definition.
    target_compile_definitions(leopard2_api_fuzzer PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_X86_TARGET)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_DISABLE_SSSE3_CODEGEN=1
            LEO2_DISABLE_AVX2_CODEGEN=1)
    endif()

    add_executable(leopard2_pruned_plan_fuzzer
        tests/leopard2/fuzz_pruned_transform.cpp
        ${LIB_SOURCE_FILES})
    target_include_directories(leopard2_pruned_plan_fuzzer PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_options(leopard2_pruned_plan_fuzzer PRIVATE
        ${LEO2_FUZZ_SANITIZER_FLAGS})
    target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_X86_TARGET)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            LEO2_DISABLE_SSSE3_CODEGEN=1
            LEO2_DISABLE_AVX2_CODEGEN=1)
    endif()

    # The fuzzer executables compile private sanitizer-instrumented copies of
    # the core rather than linking leopard.  Mirror the library's field source
    # and omission definitions so reduced-field builds exercise the same graph
    # and never retain a disabled field by accident.
    if(LEOPARD_ENABLE_GF8)
        target_sources(leopard2_api_fuzzer PRIVATE
            LeopardFF8.cpp LeopardFF8.h)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            LeopardFF8.cpp LeopardFF8.h)
    else()
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            NO_LEO_HAS_FF8=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            NO_LEO_HAS_FF8=1)
    endif()
    if(LEOPARD_ENABLE_GF16)
        target_sources(leopard2_api_fuzzer PRIVATE
            LeopardFF16.cpp LeopardFF16.h)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            LeopardFF16.cpp LeopardFF16.h)
    else()
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            NO_LEO_HAS_FF16=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            NO_LEO_HAS_FF16=1)
    endif()

    # Experiment W's parser is deliberately absent from leopard and the
    # installed API.  This opt-in target nevertheless gives its complete C99
    # deserialize/validate/serialize boundary native libFuzzer, ASan, and
    # UBSan coverage.
    add_executable(leopard2_code_identity_fuzzer
        experiments/leopard2/code_identity/fuzz_code_identity.c
        experiments/leopard2/code_identity/code_identity.c)
    target_include_directories(leopard2_code_identity_fuzzer PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/code_identity)
    target_compile_options(leopard2_code_identity_fuzzer PRIVATE
        ${LEO2_FUZZ_SANITIZER_FLAGS}
        -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
        -Wshadow -Werror)

    # Do not silently turn a forced fuzzer into a scalar-only build.  Compile
    # sanitizer-instrumented copies of the same isolated backend members so
    # coverage, ASan, and UBSan include the fixed-kernel boundary too.
    if(LEO2_HAVE_SSSE3_BACKEND)
        add_library(leopard2_backend_ssse3_fuzz OBJECT
            Leopard2BackendSSSE3.cpp)
        target_include_directories(leopard2_backend_ssse3_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_ssse3_fuzz PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_ssse3_fuzz PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_ssse3_fuzz PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_ssse3_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS} -mssse3 -mno-avx)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_ssse3_fuzz>)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_ssse3_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_SSSE3_BACKEND=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            LEO2_HAVE_SSSE3_BACKEND=1)
    endif()
    if(LEO2_HAVE_AVX2_BACKEND)
        add_library(leopard2_backend_avx2_fuzz OBJECT
            Leopard2BackendAVX2.cpp)
        if(LEOPARD_ENABLE_GF8)
            target_sources(leopard2_backend_avx2_fuzz PRIVATE
                Leopard2BackendAVX2Xor.cpp)
        endif()
        target_include_directories(leopard2_backend_avx2_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_avx2_fuzz PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_avx2_fuzz PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_avx2_fuzz PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_avx2_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS} -mavx2 -mno-avx512f)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx2_fuzz>)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx2_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_AVX2_BACKEND=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            LEO2_HAVE_AVX2_BACKEND=1)
    endif()
    if(LEO2_HAVE_AVX512_BACKEND)
        add_library(leopard2_backend_avx512_fuzz OBJECT
            Leopard2BackendAVX512.cpp)
        target_include_directories(leopard2_backend_avx512_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_avx512_fuzz PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HAVE_AVX2_BACKEND=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_avx512_fuzz PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_avx512_fuzz PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_avx512_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS}
            -mavx2 -mavx512f -mavx512bw -mavx512vl
            -mprefer-vector-width=256)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx512_fuzz>)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx512_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_AVX512_BACKEND=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            LEO2_HAVE_AVX512_BACKEND=1)
    endif()
    if(LEO2_HAVE_GFNI_BACKEND)
        add_library(leopard2_backend_gfni_fuzz OBJECT
            Leopard2BackendGFNI.cpp)
        target_include_directories(leopard2_backend_gfni_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_definitions(leopard2_backend_gfni_fuzz PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_HAVE_AVX2_BACKEND=1
            LEO2_HAVE_GFNI_BACKEND=1)
        if(NOT LEOPARD_ENABLE_GF8)
            target_compile_definitions(leopard2_backend_gfni_fuzz PRIVATE
                NO_LEO_HAS_FF8=1)
        endif()
        if(NOT LEOPARD_ENABLE_GF16)
            target_compile_definitions(leopard2_backend_gfni_fuzz PRIVATE
                NO_LEO_HAS_FF16=1)
        endif()
        target_compile_options(leopard2_backend_gfni_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS}
            -mavx2 -mgfni -mno-avx512f)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_gfni_fuzz>)
        target_sources(leopard2_pruned_plan_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_gfni_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_GFNI_BACKEND=1)
        target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
            LEO2_HAVE_GFNI_BACKEND=1)
    endif()
    set_property(TARGET leopard2_api_fuzzer APPEND_STRING PROPERTY
        LINK_FLAGS " -fsanitize=fuzzer,address,undefined")
    set_property(TARGET leopard2_pruned_plan_fuzzer APPEND_STRING PROPERTY
        LINK_FLAGS " -fsanitize=fuzzer,address,undefined")
    set_property(TARGET leopard2_code_identity_fuzzer APPEND_STRING PROPERTY
        LINK_FLAGS " -fsanitize=fuzzer,address,undefined")
    target_link_libraries(leopard2_api_fuzzer Threads::Threads)
    target_link_libraries(leopard2_pruned_plan_fuzzer Threads::Threads)
    if(NOT LEO2_BACKEND_VARIANT STREQUAL "auto")
        if(LEO2_BACKEND_VARIANT STREQUAL "scalar")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SCALAR=1)
            target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SCALAR=1)
        elseif(LEO2_BACKEND_VARIANT STREQUAL "ssse3")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SSSE3=1)
            target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SSSE3=1)
        elseif(LEO2_BACKEND_VARIANT STREQUAL "avx2")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_AVX2=1)
            target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_AVX2=1)
        elseif(LEO2_BACKEND_VARIANT STREQUAL "avx512")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_AVX512=1)
            target_compile_definitions(leopard2_pruned_plan_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_AVX512=1)
        endif()
    endif()
endif()
