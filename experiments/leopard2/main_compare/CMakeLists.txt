cmake_minimum_required(VERSION 3.7)
project(leopard2_main_compare CXX)

include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
    "Retain exact-main comparison compile evidence" FORCE)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR
        "Exact-main comparison requires -DCMAKE_BUILD_TYPE=Release")
endif()

# Reproduce the effective Release flags in the canonical master CMake file.
# That file constructs Release flags from its Debug string, yielding this
# intentional -g/-O0/-O3 sequence and no -DNDEBUG.
set(CMAKE_CXX_FLAGS_RELEASE "-g -O0 -O3" CACHE STRING
    "Exact effective Leopard master Release optimization flags" FORCE)

set(LEOPARD_MAIN_SOURCE_DIR "" CACHE PATH
    "Detached Leopard master source at the required comparison revision")
# This is deliberately not a cache option: callers may choose the source path,
# but cannot redefine what "main baseline" means for this adapter.
set(LEOPARD_MAIN_EXPECTED_COMMIT
    "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198")

if(NOT LEOPARD_MAIN_SOURCE_DIR)
    message(FATAL_ERROR "LEOPARD_MAIN_SOURCE_DIR is required")
endif()
get_filename_component(LEOPARD_MAIN_SOURCE_DIR
    "${LEOPARD_MAIN_SOURCE_DIR}" REALPATH)

find_package(Git REQUIRED)
execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${LEOPARD_MAIN_SOURCE_DIR}" rev-parse HEAD
    RESULT_VARIABLE LEO_MAIN_HEAD_RESULT
    OUTPUT_VARIABLE LEO_MAIN_HEAD
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE LEO_MAIN_HEAD_ERROR)
if(NOT LEO_MAIN_HEAD_RESULT EQUAL 0 OR
   NOT LEO_MAIN_HEAD STREQUAL LEOPARD_MAIN_EXPECTED_COMMIT)
    message(FATAL_ERROR
        "Legacy source must be exact detached Leopard master commit "
        "${LEOPARD_MAIN_EXPECTED_COMMIT}; found '${LEO_MAIN_HEAD}': "
        "${LEO_MAIN_HEAD_ERROR}")
endif()

execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${LEOPARD_MAIN_SOURCE_DIR}"
        symbolic-ref -q HEAD
    RESULT_VARIABLE LEO_MAIN_SYMBOLIC_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(LEO_MAIN_SYMBOLIC_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Legacy source must be detached at ${LEOPARD_MAIN_EXPECTED_COMMIT}")
endif()

execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${LEOPARD_MAIN_SOURCE_DIR}"
        status --porcelain=v1 --untracked-files=no
    RESULT_VARIABLE LEO_MAIN_STATUS_RESULT
    OUTPUT_VARIABLE LEO_MAIN_STATUS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE LEO_MAIN_STATUS_ERROR)
if(NOT LEO_MAIN_STATUS_RESULT EQUAL 0 OR NOT LEO_MAIN_STATUS STREQUAL "")
    message(FATAL_ERROR
        "Legacy source has tracked changes and is not an exact baseline: "
        "${LEO_MAIN_STATUS}${LEO_MAIN_STATUS_ERROR}")
endif()

set(LEOPARD_MAIN_FILES
    leopard.cpp leopard.h
    LeopardCommon.cpp LeopardCommon.h
    LeopardFF8.cpp LeopardFF8.h
    LeopardFF16.cpp LeopardFF16.h)
foreach(LEO_MAIN_FILE IN LISTS LEOPARD_MAIN_FILES)
    if(NOT EXISTS "${LEOPARD_MAIN_SOURCE_DIR}/${LEO_MAIN_FILE}")
        message(FATAL_ERROR
            "Required main source is missing: ${LEO_MAIN_FILE}")
    endif()
endforeach()

# Build the baseline sources as their own archive.  No source or object from the
# Leopard2 checkout is part of this target.
add_library(leopard_main_exact STATIC
    "${LEOPARD_MAIN_SOURCE_DIR}/leopard.cpp"
    "${LEOPARD_MAIN_SOURCE_DIR}/leopard.h"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardCommon.cpp"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardCommon.h"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardFF8.cpp"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardFF8.h"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardFF16.cpp"
    "${LEOPARD_MAIN_SOURCE_DIR}/LeopardFF16.h")
target_include_directories(leopard_main_exact PUBLIC
    "${LEOPARD_MAIN_SOURCE_DIR}")

# Match the master build's native-code policy rather than quietly benchmarking
# a differently compiled baseline.
check_cxx_compiler_flag("-march=native" LEO_MAIN_HAS_MARCH_NATIVE)
check_cxx_compiler_flag("-Wall" LEO_MAIN_HAS_WALL)
check_cxx_compiler_flag("-Wextra" LEO_MAIN_HAS_WEXTRA)
if(NOT LEO_MAIN_HAS_MARCH_NATIVE OR NOT LEO_MAIN_HAS_WALL OR
   NOT LEO_MAIN_HAS_WEXTRA)
    message(FATAL_ERROR
        "Compiler cannot reproduce canonical master -march=native/-Wall/-Wextra flags")
endif()
target_compile_options(leopard_main_exact PRIVATE
    -march=native -Wall -Wextra)

find_package(OpenMP REQUIRED)
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(leopard_main_exact PUBLIC OpenMP::OpenMP_CXX)
else()
    target_compile_options(leopard_main_exact PUBLIC "${OpenMP_CXX_FLAGS}")
    target_link_libraries(leopard_main_exact PUBLIC "${OpenMP_CXX_FLAGS}")
endif()

add_executable(leopard_main_benchmark legacy_main_benchmark.cpp)
target_compile_definitions(leopard_main_benchmark PRIVATE
    LEOPARD_MAIN_SOURCE_COMMIT="${LEOPARD_MAIN_EXPECTED_COMMIT}")
target_compile_options(leopard_main_benchmark PRIVATE
    -march=native -Wall -Wextra)
target_link_libraries(leopard_main_benchmark PRIVATE leopard_main_exact)

find_program(LEO_MAIN_PYTHON_EXECUTABLE NAMES python3 python)
if(LEO_MAIN_PYTHON_EXECUTABLE)
    enable_testing()
    add_test(NAME leopard_main_benchmark_smoke
        COMMAND "${LEO_MAIN_PYTHON_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_legacy_main_benchmark.py"
            --benchmark $<TARGET_FILE:leopard_main_benchmark>
            --compile-commands "${CMAKE_BINARY_DIR}/compile_commands.json")
endif()
