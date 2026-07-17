cmake_minimum_required(VERSION 3.12)
project(leopard2_jerasure_compare LANGUAGES C CXX)

include(CTest)
include(CheckCCompilerFlag)

# This standalone project is intentionally absent from Leopard's production
# CMake graph.  It builds the pinned third-party sources into private static
# archives; a normal Leopard build never probes for, links, or loads them.
set(LEO2_JERASURE_SOURCE_ROOT "" CACHE PATH
    "Clean pinned Jerasure source checkout")
set(LEO2_GF_COMPLETE_SOURCE_ROOT "" CACHE PATH
    "Clean pinned GF-Complete source checkout")
set(LEO2_JERASURE_EXPECTED_COMMIT
    "de1739cc8483696506829b52e7fda4f6bb195e6a" CACHE STRING
    "Exact Jerasure commit accepted by this experiment")
set(LEO2_GF_COMPLETE_EXPECTED_COMMIT
    "a6862d10c9db467148f20eef2c6445ac9afd94d8" CACHE STRING
    "Exact GF-Complete commit accepted by this experiment")
set(LEO2_JERASURE_EXPECTED_LICENSE_SHA256
    "83b6b3ff237848fbccfa889bb52cfb13c331c1f83544b907617c7f8f31eb1769"
    CACHE STRING "Exact Jerasure BSD-3-Clause license digest")
set(LEO2_GF_COMPLETE_EXPECTED_LICENSE_SHA256
    "cb9790699b4a3d56a43bba1dd859f7f41361cd224e8745a24eef933ea134a280"
    CACHE STRING "Exact GF-Complete BSD-3-Clause license digest")

foreach(root IN ITEMS LEO2_JERASURE_SOURCE_ROOT LEO2_GF_COMPLETE_SOURCE_ROOT)
    if(NOT ${root})
        message(FATAL_ERROR "Set ${root}; see tools/leopard2_jerasure_compare.py")
    endif()
endforeach()

function(leo2_require_clean_source root expected_commit prefix)
    execute_process(COMMAND git -C "${root}" rev-parse HEAD
        RESULT_VARIABLE result OUTPUT_VARIABLE actual
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE error)
    if(NOT result EQUAL 0 OR NOT actual STREQUAL expected_commit)
        message(FATAL_ERROR
            "${prefix} source identity mismatch: expected ${expected_commit}, "
            "got '${actual}' (${error})")
    endif()
    execute_process(COMMAND git -C "${root}" status --porcelain
            --untracked-files=all
        RESULT_VARIABLE status_result OUTPUT_VARIABLE status
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE status_error)
    if(NOT status_result EQUAL 0 OR NOT status STREQUAL "")
        message(FATAL_ERROR "${prefix} checkout must be clean: ${status_error}${status}")
    endif()
    execute_process(COMMAND git -C "${root}" rev-parse "HEAD^{tree}"
        RESULT_VARIABLE tree_result OUTPUT_VARIABLE tree
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE tree_error)
    if(NOT tree_result EQUAL 0)
        message(FATAL_ERROR "Cannot identify ${prefix} tree: ${tree_error}")
    endif()
    set(${prefix}_COMMIT "${actual}" PARENT_SCOPE)
    set(${prefix}_TREE "${tree}" PARENT_SCOPE)
endfunction()

leo2_require_clean_source("${LEO2_JERASURE_SOURCE_ROOT}"
    "${LEO2_JERASURE_EXPECTED_COMMIT}" JERASURE)
leo2_require_clean_source("${LEO2_GF_COMPLETE_SOURCE_ROOT}"
    "${LEO2_GF_COMPLETE_EXPECTED_COMMIT}" GF_COMPLETE)

file(SHA256 "${LEO2_JERASURE_SOURCE_ROOT}/License.txt" JERASURE_LICENSE_SHA256)
file(SHA256 "${LEO2_GF_COMPLETE_SOURCE_ROOT}/License.txt" GF_COMPLETE_LICENSE_SHA256)
file(SHA256 "${LEO2_JERASURE_SOURCE_ROOT}/include/jerasure.h" JERASURE_HEADER_SHA256)
file(SHA256 "${LEO2_JERASURE_SOURCE_ROOT}/include/reed_sol.h" REED_SOL_HEADER_SHA256)
file(SHA256 "${LEO2_GF_COMPLETE_SOURCE_ROOT}/include/gf_complete.h"
    GF_COMPLETE_HEADER_SHA256)
if(NOT JERASURE_LICENSE_SHA256 STREQUAL LEO2_JERASURE_EXPECTED_LICENSE_SHA256)
    message(FATAL_ERROR "Jerasure license identity mismatch")
endif()
if(NOT GF_COMPLETE_LICENSE_SHA256 STREQUAL LEO2_GF_COMPLETE_EXPECTED_LICENSE_SHA256)
    message(FATAL_ERROR "GF-Complete license identity mismatch")
endif()

set(GF_COMPLETE_INCLUDE "${LEO2_GF_COMPLETE_SOURCE_ROOT}/include")
set(JERASURE_INCLUDE "${LEO2_JERASURE_SOURCE_ROOT}/include")

# Match the upstream split: gf_method deliberately avoids broad SIMD flags,
# while the arithmetic archive may use the host-supported x86 flags.
add_library(leopard2_gf_complete_util STATIC
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_method.c")
target_include_directories(leopard2_gf_complete_util PRIVATE "${GF_COMPLETE_INCLUDE}")
target_compile_options(leopard2_gf_complete_util PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-O3;-fPIC;-Wsign-compare>)

add_library(leopard2_gf_complete_external STATIC
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_wgen.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w4.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w8.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w16.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w32.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w64.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_w128.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_rand.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_general.c"
    "${LEO2_GF_COMPLETE_SOURCE_ROOT}/src/gf_cpu.c")
target_include_directories(leopard2_gf_complete_external PUBLIC "${GF_COMPLETE_INCLUDE}")
target_link_libraries(leopard2_gf_complete_external PUBLIC leopard2_gf_complete_util)
target_compile_options(leopard2_gf_complete_external PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-O3;-fPIC;-Wsign-compare>)

set(LEO2_GF_COMPLETE_SIMD_FLAGS "")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i.86)$")
    foreach(flag IN ITEMS -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mpclmul)
        string(MAKE_C_IDENTIFIER "${flag}" flag_id)
        check_c_compiler_flag("${flag}" "LEO2_HAVE_${flag_id}")
        if(LEO2_HAVE_${flag_id})
            target_compile_options(leopard2_gf_complete_external PRIVATE "${flag}")
            list(APPEND LEO2_GF_COMPLETE_SIMD_FLAGS "${flag}")
        endif()
    endforeach()
endif()
string(JOIN "," LEO2_GF_COMPLETE_SIMD_FLAGS_TEXT ${LEO2_GF_COMPLETE_SIMD_FLAGS})

add_library(leopard2_jerasure_external STATIC
    "${LEO2_JERASURE_SOURCE_ROOT}/src/galois.c"
    "${LEO2_JERASURE_SOURCE_ROOT}/src/jerasure.c"
    "${LEO2_JERASURE_SOURCE_ROOT}/src/reed_sol.c"
    "${LEO2_JERASURE_SOURCE_ROOT}/src/cauchy.c"
    "${LEO2_JERASURE_SOURCE_ROOT}/src/liberation.c")
target_include_directories(leopard2_jerasure_external PUBLIC
    "${JERASURE_INCLUDE}" "${GF_COMPLETE_INCLUDE}")
target_link_libraries(leopard2_jerasure_external PUBLIC leopard2_gf_complete_external)
target_compile_options(leopard2_jerasure_external PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-O3;-fPIC;-Wsign-compare>)

add_executable(leopard2_jerasure_benchmark jerasure_benchmark.cpp)
target_compile_features(leopard2_jerasure_benchmark PRIVATE cxx_std_11)
target_compile_options(leopard2_jerasure_benchmark PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic>)
target_compile_definitions(leopard2_jerasure_benchmark PRIVATE
    LEO2_JERASURE_COMMIT="${JERASURE_COMMIT}"
    LEO2_JERASURE_TREE="${JERASURE_TREE}"
    LEO2_JERASURE_HEADER_SHA256="${JERASURE_HEADER_SHA256}"
    LEO2_REED_SOL_HEADER_SHA256="${REED_SOL_HEADER_SHA256}"
    LEO2_JERASURE_LICENSE_SHA256="${JERASURE_LICENSE_SHA256}"
    LEO2_GF_COMPLETE_COMMIT="${GF_COMPLETE_COMMIT}"
    LEO2_GF_COMPLETE_TREE="${GF_COMPLETE_TREE}"
    LEO2_GF_COMPLETE_HEADER_SHA256="${GF_COMPLETE_HEADER_SHA256}"
    LEO2_GF_COMPLETE_LICENSE_SHA256="${GF_COMPLETE_LICENSE_SHA256}"
    LEO2_GF_COMPLETE_SIMD_FLAGS="${LEO2_GF_COMPLETE_SIMD_FLAGS_TEXT}")
target_link_libraries(leopard2_jerasure_benchmark PRIVATE leopard2_jerasure_external)

if(BUILD_TESTING)
    add_test(NAME jerasure_compare_roundtrip_gf8
        COMMAND leopard2_jerasure_benchmark
            --k 129 --r 100 --profile high --bytes 64 --loss 8
            --batch 2 --reuse 2 --iterations 3 --warmup 0 --threads 1
            --seed 501 --oracle full --json
            "${CMAKE_CURRENT_BINARY_DIR}/roundtrip-gf8.json")
    add_test(NAME jerasure_compare_roundtrip_gf16
        COMMAND leopard2_jerasure_benchmark
            --k 200 --r 100 --profile high --bytes 72 --loss 4
            --batch 1 --reuse 2 --iterations 3 --warmup 0 --threads 1
            --seed 502 --oracle full --json
            "${CMAKE_CURRENT_BINARY_DIR}/roundtrip-gf16.json")
    add_test(NAME jerasure_compare_rejects_unaligned_region
        COMMAND leopard2_jerasure_benchmark
            --k 8 --r 4 --profile high --bytes 65 --loss 1
            --batch 1 --reuse 1 --iterations 3 --warmup 0 --threads 1
            --seed 503 --oracle full --json
            "${CMAKE_CURRENT_BINARY_DIR}/invalid.json")
    set_tests_properties(jerasure_compare_rejects_unaligned_region
        PROPERTIES WILL_FAIL TRUE)
    add_test(NAME jerasure_compare_rejects_threads
        COMMAND leopard2_jerasure_benchmark
            --k 8 --r 4 --profile high --bytes 64 --loss 1
            --batch 1 --reuse 1 --iterations 3 --warmup 0 --threads 2
            --seed 504 --oracle full --json
            "${CMAKE_CURRENT_BINARY_DIR}/invalid-thread.json")
    set_tests_properties(jerasure_compare_rejects_threads PROPERTIES WILL_FAIL TRUE)
endif()
