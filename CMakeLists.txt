cmake_minimum_required(VERSION 3.7)
project(leopard)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 11)

set(LIB_SOURCE_FILES
        leopard.cpp
        leopard.h
        leopard2.cpp
        leopard2.h
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

check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

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

add_library(libleopard STATIC ${LIB_SOURCE_FILES})
target_include_directories(libleopard PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libleopard PUBLIC Threads::Threads)

if(LEO2_BUILD_BENCHMARKS)
    add_executable(bench_leopard ${BENCH_SOURCE_FILES})
    target_link_libraries(bench_leopard libleopard)

    add_executable(bench_leopard2 bench/leopard2/benchmark.cpp)
    target_link_libraries(bench_leopard2 libleopard)

    add_executable(experiment_leopard ${EXPERIMENT_SOURCE_FILES})
endif()

if(LEO2_BUILD_TESTS)
    enable_testing()
    add_executable(leopard2_direct_oracle_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_direct_oracle.cpp)
    target_include_directories(leopard2_direct_oracle_test PRIVATE tests/leopard2)
    add_test(NAME leopard2_direct_oracle COMMAND leopard2_direct_oracle_test)

    add_executable(leopard2_api_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_api.cpp)
    target_include_directories(leopard2_api_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_api_test libleopard)
    add_test(NAME leopard2_api COMMAND leopard2_api_test)

    add_executable(leopard2_random_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_random.cpp)
    target_include_directories(leopard2_random_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_random_test libleopard)
    add_test(NAME leopard2_random COMMAND leopard2_random_test)

    add_executable(leopard2_locator_test tests/leopard2/test_locator.cpp)
    target_link_libraries(leopard2_locator_test libleopard)
    add_test(NAME leopard2_locator COMMAND leopard2_locator_test)

    add_executable(leopard2_gf16_tails_test
        tests/leopard2/direct_oracle.cpp
        tests/leopard2/test_gf16_tails.cpp)
    target_include_directories(leopard2_gf16_tails_test PRIVATE tests/leopard2)
    target_link_libraries(leopard2_gf16_tails_test libleopard)
    add_test(NAME leopard2_gf16_tails COMMAND leopard2_gf16_tails_test)

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
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/test_cuda_optional.cmake)

    if(CMAKE_VERSION VERSION_LESS 3.12)
        find_package(PythonInterp QUIET)
        set(LEO2_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
    else()
        find_package(Python3 COMPONENTS Interpreter QUIET)
        set(LEO2_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
    endif()
    if(LEO2_PYTHON_EXECUTABLE)
        add_test(
            NAME leopard2_cost_model_self_test
            COMMAND ${LEO2_PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/experiments/leopard2/cost_model.py
                self-test)
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
    add_executable(leopard2_api_fuzzer
        tests/leopard2/fuzz_api.cpp
        ${LIB_SOURCE_FILES})
    target_include_directories(leopard2_api_fuzzer PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_options(leopard2_api_fuzzer PRIVATE
        -fsanitize=fuzzer,address,undefined
        -fno-omit-frame-pointer)
    set_property(TARGET leopard2_api_fuzzer APPEND_STRING PROPERTY
        LINK_FLAGS " -fsanitize=fuzzer,address,undefined")
    target_link_libraries(leopard2_api_fuzzer Threads::Threads)
endif()
