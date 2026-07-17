cmake_minimum_required(VERSION 3.7)
project(leopard)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

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
        LeopardCommon.h
        LeopardFF16.cpp
        LeopardFF16.h
        LeopardFF8.cpp
        LeopardFF8.h)

set(BENCH_SOURCE_FILES
        tests/benchmark.cpp)

set(EXPERIMENT_SOURCE_FILES
        tests/experiments.cpp)

option(LEO2_BUILD_TESTS "Build Leopard2 correctness tests" ON)
option(LEO2_BUILD_BENCHMARKS "Build Leopard benchmark programs" ON)
option(LEO2_BUILD_FUZZERS "Build Leopard2 libFuzzer targets" OFF)
option(LEO2_ENABLE_CUDA "Build the optional Leopard2 CUDA backend" OFF)

# Diagnostic-only backend selection.  auto is the portable default and does
# not add definitions or target-local ISA options.  Forced variants hold the
# wire format fixed while making backend determinism testable in isolated
# builds.
set(LEO2_BACKEND_VARIANT "auto" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, or avx2")
set_property(CACHE LEO2_BACKEND_VARIANT PROPERTY STRINGS auto scalar ssse3 avx2)
string(TOLOWER "${LEO2_BACKEND_VARIANT}" LEO2_BACKEND_VARIANT_NORMALIZED)
if(NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "auto" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "scalar" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "ssse3" AND
   NOT LEO2_BACKEND_VARIANT_NORMALIZED STREQUAL "avx2")
    message(FATAL_ERROR
        "Invalid LEO2_BACKEND_VARIANT='${LEO2_BACKEND_VARIANT}'. "
        "Expected auto, scalar, ssse3, or avx2")
endif()
set(LEO2_BACKEND_VARIANT "${LEO2_BACKEND_VARIANT_NORMALIZED}" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, or avx2" FORCE)
set_property(CACHE LEO2_BACKEND_VARIANT PROPERTY STRINGS auto scalar ssse3 avx2)

# CUDA is deliberately opt-in.  In particular, the default configuration must
# not probe for a CUDA compiler or make CUDA a transitive build dependency.
# Backend sources will be added behind this gate by the end-of-project CUDA
# work tracked in Beads.
if(LEO2_ENABLE_CUDA)
    enable_language(CUDA)
endif()

if(NOT CMAKE_BUILD_TYPE)
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
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3")
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
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3")
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

set(LEO2_X86_TARGET OFF)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "X86" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "AMD64" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "i386" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "i486" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "i586" OR
   CMAKE_SYSTEM_PROCESSOR STREQUAL "i686")
    set(LEO2_X86_TARGET ON)
endif()

add_library(libleopard STATIC ${LIB_SOURCE_FILES})
if(LEO2_X86_TARGET)
    target_compile_definitions(libleopard PRIVATE
        LEO2_DISABLE_SSSE3_CODEGEN=1
        LEO2_DISABLE_AVX2_CODEGEN=1)
endif()
target_include_directories(libleopard PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(libleopard PUBLIC Threads::Threads)

# Optional x86 kernels are compiled in dedicated object libraries.  Nothing in
# libleopard's baseline translation units is built above the x86-64 SSE2 floor;
# runtime feature checks select one immutable ops table after startup KATs.
set(LEO2_HAVE_SSSE3_BACKEND OFF)
set(LEO2_HAVE_AVX2_BACKEND OFF)
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
            target_include_directories(leopard2_backend_avx2 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})
            target_compile_options(leopard2_backend_avx2 PRIVATE /arch:AVX2)
            set(LEO2_HAVE_AVX2_BACKEND ON)
        endif()
    endif()
endif()

if(LEO2_HAVE_SSSE3_BACKEND)
    target_sources(libleopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_ssse3>)
    target_compile_definitions(libleopard PRIVATE LEO2_HAVE_SSSE3_BACKEND=1)
endif()
if(LEO2_HAVE_AVX2_BACKEND)
    target_sources(libleopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_avx2>)
    target_compile_definitions(libleopard PRIVATE LEO2_HAVE_AVX2_BACKEND=1)
endif()
if(ENABLE_OPENMP)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(libleopard PUBLIC OpenMP::OpenMP_CXX)
    else()
        # FindOpenMP did not provide imported targets before CMake 3.9.  Keep
        # the 3.7 project floor usable by propagating its compiler-driver flag;
        # for static archives this is also the link flag that selects the
        # OpenMP runtime on GCC and Clang.
        target_link_libraries(libleopard PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Keep the installed CPU package independent of optional backends.  CUDA
# targets, once implemented, must use a separate opt-in export set rather than
# becoming transitive dependencies of leopard::libleopard.
set(LEOPARD_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"
    CACHE STRING "Install directory for Leopard CMake package files")
configure_package_config_file(
    cmake/leopardConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake"
    INSTALL_DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")
install(TARGETS libleopard
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
            "processor '${CMAKE_SYSTEM_PROCESSOR}'; use auto")
    endif()

    if(LEO2_BACKEND_VARIANT STREQUAL "scalar")
        target_compile_definitions(libleopard PRIVATE
            LEO2_BACKEND_FORCE_SCALAR=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "ssse3")
        if(NOT LEO2_HAVE_SSSE3_BACKEND)
            message(FATAL_ERROR
                "LEO2_BACKEND_VARIANT=ssse3 is unavailable for this compiler")
        endif()
        target_compile_definitions(libleopard PRIVATE
            LEO2_BACKEND_FORCE_SSSE3=1)
    elseif(LEO2_BACKEND_VARIANT STREQUAL "avx2")
        if(NOT LEO2_HAVE_AVX2_BACKEND)
            message(FATAL_ERROR
                "LEO2_BACKEND_VARIANT=avx2 is unavailable for this compiler")
        endif()
        target_compile_definitions(libleopard PRIVATE
            LEO2_BACKEND_FORCE_AVX2=1)
    endif()
endif()

if(LEO2_BUILD_BENCHMARKS)
    add_executable(bench_leopard ${BENCH_SOURCE_FILES})
    target_link_libraries(bench_leopard libleopard)

    add_executable(bench_leopard2 bench/leopard2/benchmark.cpp)
    target_link_libraries(bench_leopard2 libleopard)

    add_executable(bench_leopard2_locator
        bench/leopard2/locator_benchmark.cpp)
    target_link_libraries(bench_leopard2_locator libleopard)
    set(LEO2_LOCATOR_SOURCE_GIT_SHA "unknown")
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
    target_compile_definitions(bench_leopard2_locator PRIVATE
        LEO2_LOCATOR_SOURCE_GIT_SHA="${LEO2_LOCATOR_SOURCE_GIT_SHA}"
        LEO2_LOCATOR_SOURCE_DIRTY=${LEO2_LOCATOR_SOURCE_DIRTY})

    add_executable(experiment_leopard ${EXPERIMENT_SOURCE_FILES})
endif()

if(LEO2_BUILD_TESTS)
    enable_testing()
    target_compile_definitions(libleopard PRIVATE LEO2_ENABLE_TEST_HOOKS=1)
    if(LEO2_HAVE_SSSE3_BACKEND)
        target_compile_definitions(leopard2_backend_ssse3 PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
    endif()
    if(LEO2_HAVE_AVX2_BACKEND)
        target_compile_definitions(leopard2_backend_avx2 PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
    endif()

    if(LEO2_BUILD_BENCHMARKS)
        add_executable(bench_leopard2_direct_encode
            bench/leopard2/direct_encode_benchmark.cpp)
        target_compile_definitions(bench_leopard2_direct_encode PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1)
        target_link_libraries(bench_leopard2_direct_encode libleopard)

        # C1 is internal and absent from production dispatch.  This benchmark
        # is built only with the test hooks that expose the full padded LCH
        # oracle; it is never installed and cannot affect production behavior.
        add_executable(bench_leopard2_pruned_transform
            bench/leopard2/pruned_transform_benchmark.cpp)
        target_compile_definitions(bench_leopard2_pruned_transform PRIVATE
            LEO2_ENABLE_TEST_HOOKS=1
            LEO2_PRUNED_SOURCE_GIT_SHA="${LEO2_LOCATOR_SOURCE_GIT_SHA}"
            LEO2_PRUNED_SOURCE_DIRTY=${LEO2_LOCATOR_SOURCE_DIRTY})
        target_link_libraries(bench_leopard2_pruned_transform libleopard)

    endif()

    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR
       CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64" OR
       CMAKE_SYSTEM_PROCESSOR STREQUAL "AMD64")
        find_program(LEO2_OBJDUMP_EXECUTABLE NAMES objdump llvm-objdump)
        find_program(LEO2_SH_EXECUTABLE NAMES sh)
        if(LEO2_OBJDUMP_EXECUTABLE AND LEO2_SH_EXECUTABLE)
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
            endif()
            string(REPLACE ";" "," LEO2_PORTABLE_ISA_EXPECTED_CLASSES_ARG
                "${LEO2_PORTABLE_ISA_EXPECTED_CLASSES}")
            add_test(NAME leopard2_portable_isa
                COMMAND "${LEO2_SH_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_leopard2_portable_isa.sh"
                    "${LEO2_OBJDUMP_EXECUTABLE}"
                    "$<TARGET_FILE:libleopard>"
                    "${CMAKE_CURRENT_BINARY_DIR}"
                    ${LEO2_PORTABLE_ISA_SELF_TEST_ARGS}
                    "${LEO2_PORTABLE_ISA_EXPECTED_CLASSES_ARG}")
        endif()
    endif()

    add_executable(leopard2_direct_oracle_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_direct_oracle.cpp)
    target_include_directories(leopard2_direct_oracle_test PRIVATE tests/leopard2)
    add_test(NAME leopard2_direct_oracle COMMAND leopard2_direct_oracle_test)

    add_executable(leopard2_backend_ops_test
        tests/leopard2/test_backend_ops.cpp)
    target_link_libraries(leopard2_backend_ops_test libleopard)
    add_test(NAME leopard2_backend_ops COMMAND leopard2_backend_ops_test)
    set_tests_properties(leopard2_backend_ops PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_backend_failures_test
        tests/leopard2/test_backend_failures.cpp)
    target_compile_definitions(leopard2_backend_failures_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_backend_failures_test libleopard)
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

    add_executable(leopard2_context_backends_test
        tests/leopard2/test_context_backends.cpp)
    target_compile_definitions(leopard2_context_backends_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_context_backends_test libleopard)
    add_test(NAME leopard2_context_backends
        COMMAND leopard2_context_backends_test)
    set_tests_properties(leopard2_context_backends PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=4")

    add_executable(leopard2_r1_xor_test
        tests/leopard2/test_r1_xor.cpp)
    target_compile_definitions(leopard2_r1_xor_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_link_libraries(leopard2_r1_xor_test libleopard)
    add_test(NAME leopard2_r1_xor COMMAND leopard2_r1_xor_test)
    set_tests_properties(leopard2_r1_xor PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_executable(leopard2_direct_encode_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_direct_encode.cpp)
    target_compile_definitions(leopard2_direct_encode_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_include_directories(leopard2_direct_encode_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_direct_encode_test libleopard)
    add_test(NAME leopard2_direct_encode COMMAND leopard2_direct_encode_test)
    set_tests_properties(leopard2_direct_encode PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_active_lch_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_active_lch.cpp)
    target_compile_definitions(leopard2_active_lch_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_include_directories(leopard2_active_lch_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_active_lch_test libleopard)
    add_test(NAME leopard2_active_lch COMMAND leopard2_active_lch_test)
    set_tests_properties(leopard2_active_lch PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_api_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_api.cpp)
    target_include_directories(leopard2_api_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_api_test libleopard)
    add_test(NAME leopard2_api COMMAND leopard2_api_test)

    add_executable(leopard2_public_api_contract_test
        tests/leopard2/test_public_api_contract.cpp)
    target_link_libraries(leopard2_public_api_contract_test libleopard)
    add_test(NAME leopard2_public_api_contract
        COMMAND leopard2_public_api_contract_test)

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND ENABLE_OPENMP)
        add_executable(leopard2_initialization_threads_test
            tests/leopard2/test_initialization_threads.cpp)
        target_link_libraries(leopard2_initialization_threads_test libleopard)
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
    target_link_libraries(leopard2_random_test libleopard)
    add_test(NAME leopard2_random COMMAND leopard2_random_test)

    add_executable(leopard2_locator_test tests/leopard2/test_locator.cpp)
    target_link_libraries(leopard2_locator_test libleopard)
    add_test(NAME leopard2_locator COMMAND leopard2_locator_test)

    add_executable(leopard2_decode_plan_schedule_test
        tests/leopard2/test_decode_plan_schedule.cpp)
    target_link_libraries(leopard2_decode_plan_schedule_test libleopard)
    add_test(NAME leopard2_decode_plan_schedule
        COMMAND leopard2_decode_plan_schedule_test)
    set_tests_properties(leopard2_decode_plan_schedule PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_gf16_tails_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_gf16_tails.cpp)
    target_include_directories(leopard2_gf16_tails_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_gf16_tails_test libleopard)
    add_test(NAME leopard2_gf16_tails COMMAND leopard2_gf16_tails_test)

    add_executable(leopard2_gf16_padded_odd_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_gf16_padded_odd.cpp)
    target_include_directories(leopard2_gf16_padded_odd_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_gf16_padded_odd_test libleopard)
    add_test(NAME leopard2_gf16_padded_odd COMMAND leopard2_gf16_padded_odd_test)
    set_tests_properties(leopard2_gf16_padded_odd PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_gf16_legacy_encoder_matrix_test
        tests/leopard2/test_encoder_gf16_legacy_matrix.cpp)
    target_link_libraries(leopard2_gf16_legacy_encoder_matrix_test libleopard)
    add_test(NAME leopard2_gf16_legacy_encoder_matrix
        COMMAND leopard2_gf16_legacy_encoder_matrix_test)
    set_tests_properties(leopard2_gf16_legacy_encoder_matrix PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_low_gf16_direct_rows_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_low_gf16_direct_rows.cpp)
    target_include_directories(leopard2_low_gf16_direct_rows_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_low_gf16_direct_rows_test libleopard)
    add_test(NAME leopard2_low_gf16_direct_rows
        COMMAND leopard2_low_gf16_direct_rows_test)
    set_tests_properties(leopard2_low_gf16_direct_rows PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_decode_high_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_decode_high_acceptance.cpp)
    target_include_directories(leopard2_decode_high_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_decode_high_acceptance_test libleopard)
    add_test(NAME leopard2_decode_high_acceptance
        COMMAND leopard2_decode_high_acceptance_test)
    set_tests_properties(leopard2_decode_high_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_decode_low_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_decode_low_acceptance.cpp)
    target_include_directories(leopard2_decode_low_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_decode_low_acceptance_test libleopard)
    add_test(NAME leopard2_decode_low_acceptance
        COMMAND leopard2_decode_low_acceptance_test)
    set_tests_properties(leopard2_decode_low_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_arbitrary_counts_acceptance_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_arbitrary_counts_acceptance.cpp)
    target_include_directories(leopard2_arbitrary_counts_acceptance_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_arbitrary_counts_acceptance_test libleopard)
    add_test(NAME leopard2_arbitrary_counts_acceptance
        COMMAND leopard2_arbitrary_counts_acceptance_test)
    set_tests_properties(leopard2_arbitrary_counts_acceptance PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_max_counts_test tests/leopard2/test_max_counts.cpp)
    target_link_libraries(leopard2_max_counts_test libleopard)
    add_test(NAME leopard2_max_counts COMMAND leopard2_max_counts_test)
    set_tests_properties(leopard2_max_counts PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(leopard2_encode_concurrency_test
        tests/leopard2/test_encode_concurrency.cpp)
    target_link_libraries(leopard2_encode_concurrency_test libleopard)
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
    target_link_libraries(leopard2_codec_options_abi_test libleopard)
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
    target_link_libraries(leopard2_legacy_golden_test libleopard)
    add_test(NAME leopard2_legacy_golden COMMAND leopard2_legacy_golden_test)

    add_executable(leopard2_boundaries_test
        tests/leopard2/test_boundaries.cpp)
    target_link_libraries(leopard2_boundaries_test libleopard)
    add_test(NAME leopard2_boundaries COMMAND leopard2_boundaries_test)

    add_executable(leopard2_transform_differential_test
        tests/leopard2/test_transform_differential.cpp
        tests/leopard2/direct_oracle.cpp)
    target_link_libraries(leopard2_transform_differential_test libleopard)
    add_test(NAME leopard2_transform_differential
        COMMAND leopard2_transform_differential_test)

    add_executable(leopard2_pruned_transform_test
        tests/leopard2/test_pruned_transform.cpp
        tests/leopard2/direct_oracle.cpp)
    target_compile_definitions(leopard2_pruned_transform_test PRIVATE
        LEO2_ENABLE_TEST_HOOKS=1)
    target_include_directories(leopard2_pruned_transform_test
        PRIVATE tests/leopard2)
    target_link_libraries(leopard2_pruned_transform_test libleopard)
    add_test(NAME leopard2_pruned_transform
        COMMAND leopard2_pruned_transform_test)
    set_tests_properties(leopard2_pruned_transform PROPERTIES
        ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=1")

    add_executable(leopard2_fuzz_smoke
        tests/leopard2/fuzz_api.cpp
        tests/leopard2/fuzz_replay.cpp)
    target_link_libraries(leopard2_fuzz_smoke libleopard)
    add_test(NAME leopard2_fuzz_smoke COMMAND leopard2_fuzz_smoke)

    add_test(
        NAME leopard2_cuda_optional
        COMMAND ${CMAKE_COMMAND}
            -DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/cuda-optional-test
            -DLEO2_GENERATOR=${CMAKE_GENERATOR}
            -DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
            -DLEO2_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_cuda_optional.cmake)

    if(CMAKE_VERSION VERSION_LESS 3.12)
        find_package(PythonInterp QUIET)
        set(LEO2_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
    else()
        find_package(Python3 COMPONENTS Interpreter QUIET)
        set(LEO2_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
    endif()
    if(LEO2_PYTHON_EXECUTABLE)
        if(LEO2_BUILD_BENCHMARKS)
            add_test(
                NAME leopard2_locator_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_locator_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2_locator>)
            set_tests_properties(leopard2_locator_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
            add_test(
                NAME leopard2_benchmark_json_regression
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2>)
            add_test(
                NAME leopard2_pruned_transform_benchmark_smoke
                COMMAND ${LEO2_PYTHON_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_pruned_transform_benchmark_json_test.py
                    $<TARGET_FILE:bench_leopard2_pruned_transform>)
            set_tests_properties(
                leopard2_pruned_transform_benchmark_smoke PROPERTIES
                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=FALSE")
        endif()
        add_test(
            NAME leopard2_cost_model_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/cost_model.py
                self-test)
        add_test(
            NAME leopard2_lab_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_lab.py
                self-test)
        add_test(
            NAME leopard2_benchmark_matrix_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_benchmark_matrix.py
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
        add_test(
            NAME leopard2_direct_encode_crossover_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/leopard2_direct_encode_crossover.py
                self-test)
        add_test(
            NAME leopard2_visual_studio_project_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/proj/test_leopard_vcxproj.py)
    endif()
endif()

if(LEO2_BUILD_FUZZERS)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR
            "LEO2_BUILD_FUZZERS requires upstream Clang libFuzzer support")
    endif()

    # Build an independently instrumented copy of the library into the fuzz
    # executable.  Applying sanitizer flags only to fuzz_api.cpp would leave
    # the codec itself invisible to ASan, UBSan, and coverage-guided mutation;
    # applying them to libleopard would in turn break ordinary targets that
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
    if(LEO2_X86_TARGET)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_DISABLE_SSSE3_CODEGEN=1
            LEO2_DISABLE_AVX2_CODEGEN=1)
    endif()

    # Do not silently turn a forced fuzzer into a scalar-only build.  Compile
    # sanitizer-instrumented copies of the same isolated backend members so
    # coverage, ASan, and UBSan include the fixed-kernel boundary too.
    if(LEO2_HAVE_SSSE3_BACKEND)
        add_library(leopard2_backend_ssse3_fuzz OBJECT
            Leopard2BackendSSSE3.cpp)
        target_include_directories(leopard2_backend_ssse3_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_options(leopard2_backend_ssse3_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS} -mssse3 -mno-avx)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_ssse3_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_SSSE3_BACKEND=1)
    endif()
    if(LEO2_HAVE_AVX2_BACKEND)
        add_library(leopard2_backend_avx2_fuzz OBJECT
            Leopard2BackendAVX2.cpp)
        target_include_directories(leopard2_backend_avx2_fuzz PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
        target_compile_options(leopard2_backend_avx2_fuzz PRIVATE
            ${LEO2_FUZZ_SANITIZER_FLAGS} -mavx2 -mno-avx512f)
        target_sources(leopard2_api_fuzzer PRIVATE
            $<TARGET_OBJECTS:leopard2_backend_avx2_fuzz>)
        target_compile_definitions(leopard2_api_fuzzer PRIVATE
            LEO2_HAVE_AVX2_BACKEND=1)
    endif()
    set_property(TARGET leopard2_api_fuzzer APPEND_STRING PROPERTY
        LINK_FLAGS " -fsanitize=fuzzer,address,undefined")
    target_link_libraries(leopard2_api_fuzzer Threads::Threads)
    if(NOT LEO2_BACKEND_VARIANT STREQUAL "auto")
        if(LEO2_BACKEND_VARIANT STREQUAL "scalar")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SCALAR=1)
        elseif(LEO2_BACKEND_VARIANT STREQUAL "ssse3")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_SSSE3=1)
        elseif(LEO2_BACKEND_VARIANT STREQUAL "avx2")
            target_compile_definitions(leopard2_api_fuzzer PRIVATE
                LEO2_BACKEND_FORCE_AVX2=1)
        endif()
    endif()
endif()
