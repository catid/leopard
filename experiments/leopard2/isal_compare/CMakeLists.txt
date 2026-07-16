cmake_minimum_required(VERSION 3.12)
project(leopard2_isal_compare LANGUAGES CXX)

include(CTest)

# This is a standalone, default-off research build.  The production Leopard
# CMake project never descends into this directory and therefore never probes
# for, links, or installs ISA-L.
set(LEO2_ISAL_SOURCE_ROOT "" CACHE PATH
    "Pinned Intel ISA-L source checkout used to prove source identity")
set(LEO2_ISAL_INSTALL_ROOT "" CACHE PATH
    "Private ISA-L installation containing include/ and lib/")
set(LEO2_ISAL_EXPECTED_COMMIT
    "e8cc5e87fc64b4da434f32bc1fa18184622a4998" CACHE STRING
    "Exact ISA-L commit accepted by this experiment")

if(NOT LEO2_ISAL_SOURCE_ROOT OR NOT LEO2_ISAL_INSTALL_ROOT)
    message(FATAL_ERROR
        "Set LEO2_ISAL_SOURCE_ROOT and LEO2_ISAL_INSTALL_ROOT.  "
        "See tools/leopard2_isal_compare.py bootstrap.")
endif()

execute_process(
    COMMAND git -C "${LEO2_ISAL_SOURCE_ROOT}" rev-parse HEAD
    RESULT_VARIABLE LEO2_ISAL_GIT_RESULT
    OUTPUT_VARIABLE LEO2_ISAL_ACTUAL_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE LEO2_ISAL_GIT_ERROR)
if(NOT LEO2_ISAL_GIT_RESULT EQUAL 0 OR
   NOT LEO2_ISAL_ACTUAL_COMMIT STREQUAL LEO2_ISAL_EXPECTED_COMMIT)
    message(FATAL_ERROR
        "ISA-L source identity mismatch: expected ${LEO2_ISAL_EXPECTED_COMMIT}, "
        "got '${LEO2_ISAL_ACTUAL_COMMIT}' (${LEO2_ISAL_GIT_ERROR})")
endif()

set(LEO2_ISAL_HEADER
    "${LEO2_ISAL_INSTALL_ROOT}/include/isa-l/erasure_code.h")
set(LEO2_ISAL_LIBRARY "${LEO2_ISAL_INSTALL_ROOT}/lib/libisal.a")
if(NOT EXISTS "${LEO2_ISAL_HEADER}" OR NOT EXISTS "${LEO2_ISAL_LIBRARY}")
    message(FATAL_ERROR
        "The private ISA-L installation is incomplete: "
        "${LEO2_ISAL_INSTALL_ROOT}")
endif()
file(SHA256 "${LEO2_ISAL_LIBRARY}" LEO2_ISAL_LIBRARY_SHA256)

add_library(leopard2_isal_external STATIC IMPORTED GLOBAL)
set_target_properties(leopard2_isal_external PROPERTIES
    IMPORTED_LOCATION "${LEO2_ISAL_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LEO2_ISAL_INSTALL_ROOT}/include")

add_executable(leopard2_isal_benchmark isal_benchmark.cpp)
target_compile_features(leopard2_isal_benchmark PRIVATE cxx_std_11)
target_compile_options(leopard2_isal_benchmark PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic>)
target_compile_definitions(leopard2_isal_benchmark PRIVATE
    LEO2_ISAL_COMMIT="${LEO2_ISAL_ACTUAL_COMMIT}"
    LEO2_ISAL_LIBRARY_SHA256="${LEO2_ISAL_LIBRARY_SHA256}")
target_link_libraries(leopard2_isal_benchmark PRIVATE leopard2_isal_external)

if(BUILD_TESTING)
    add_test(NAME isal_compare_roundtrip_boundary
        COMMAND leopard2_isal_benchmark
            --k 129 --r 100 --profile high --bytes 65 --loss 100
            --batch 2 --reuse 2 --iterations 3 --warmup 1 --threads 1
            --seed 501 --json
            "${CMAKE_CURRENT_BINARY_DIR}/roundtrip-boundary.json")
    add_test(NAME isal_compare_roundtrip_no_loss
        COMMAND leopard2_isal_benchmark
            --k 64 --r 192 --profile low --bytes 17 --loss 0
            --batch 3 --reuse 2 --iterations 3 --warmup 1 --threads 1
            --seed 502 --json
            "${CMAKE_CURRENT_BINARY_DIR}/roundtrip-no-loss.json")
    add_test(NAME isal_compare_rejects_unimplemented_threads
        COMMAND leopard2_isal_benchmark
            --k 8 --r 4 --profile high --bytes 64 --loss 1
            --batch 1 --reuse 1 --iterations 3 --warmup 1 --threads 2
            --seed 503 --json
            "${CMAKE_CURRENT_BINARY_DIR}/must-not-exist.json")
    set_tests_properties(isal_compare_rejects_unimplemented_threads
        PROPERTIES WILL_FAIL TRUE)
endif()
