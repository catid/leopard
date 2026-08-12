cmake_minimum_required(VERSION 3.15)
project(leopard2_wirehair_atlas_adapter C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
    "Retain Wirehair atlas compile evidence" FORCE)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "Wirehair atlas adapter requires a Release build")
endif()

set(WIREHAIR_SOURCE_DIR "" CACHE PATH
    "Clean detached Wirehair source at the atlas revision")
set(WIREHAIR_EXPECTED_COMMIT
    "067ca7cdb66aed424ec23f97557429bf791c6f0c")
if(NOT WIREHAIR_SOURCE_DIR)
    message(FATAL_ERROR "WIREHAIR_SOURCE_DIR is required")
endif()
get_filename_component(WIREHAIR_SOURCE_DIR "${WIREHAIR_SOURCE_DIR}" REALPATH)

find_package(Git REQUIRED)
execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${WIREHAIR_SOURCE_DIR}" rev-parse HEAD
    RESULT_VARIABLE WH_HEAD_RESULT OUTPUT_VARIABLE WH_HEAD
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE WH_HEAD_ERROR)
if(NOT WH_HEAD_RESULT EQUAL 0 OR NOT WH_HEAD STREQUAL WIREHAIR_EXPECTED_COMMIT)
    message(FATAL_ERROR
        "Wirehair source must be exact ${WIREHAIR_EXPECTED_COMMIT}; "
        "found '${WH_HEAD}': ${WH_HEAD_ERROR}")
endif()
execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${WIREHAIR_SOURCE_DIR}"
        status --porcelain=v1 --untracked-files=all
    RESULT_VARIABLE WH_STATUS_RESULT OUTPUT_VARIABLE WH_STATUS
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE WH_STATUS_ERROR)
if(NOT WH_STATUS_RESULT EQUAL 0 OR NOT WH_STATUS STREQUAL "")
    message(FATAL_ERROR
        "Wirehair source is not clean: ${WH_STATUS}${WH_STATUS_ERROR}")
endif()

include(CheckCXXCompilerFlag)
foreach(flag IN ITEMS -march=x86-64 -mtune=generic -mavx2
        -mno-avx512f -mno-avx512bw -mno-avx512vl)
    string(MAKE_C_IDENTIFIER "WH_ATLAS_${flag}" probe)
    check_cxx_compiler_flag("${flag}" "${probe}")
    if(NOT ${probe})
        message(FATAL_ERROR "Compiler does not support ${flag}")
    endif()
endforeach()
set(WH_ATLAS_ISA_FLAGS
    -march=x86-64 -mtune=generic -mavx2
    -mno-avx512f -mno-avx512bw -mno-avx512vl)
add_compile_options(${WH_ATLAS_ISA_FLAGS})

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_CODEC_V2 OFF CACHE BOOL "" FORCE)
set(MARCH_NATIVE OFF CACHE BOOL "" FORCE)
set(WIREHAIR_BUILD_BOTH OFF CACHE BOOL "" FORCE)
set(WIREHAIR_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(WIREHAIR_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(WIREHAIR_ENABLE_SCHEDULED_TESTS OFF CACHE BOOL "" FORCE)
set(WIREHAIR_ENABLE_LIBFUZZER OFF CACHE BOOL "" FORCE)
add_subdirectory("${WIREHAIR_SOURCE_DIR}" wirehair-source EXCLUDE_FROM_ALL)

add_executable(wirehair_v1_benchmark wirehair_v1_benchmark.cpp)
target_compile_definitions(wirehair_v1_benchmark PRIVATE
    LEO2_ATLAS_WIREHAIR_COMMIT="${WIREHAIR_EXPECTED_COMMIT}"
    LEO2_ATLAS_PURE_AVX2=0
    LEO2_ATLAS_WIREHAIR_V1_COMPACT_PATH=1)
target_include_directories(wirehair_v1_benchmark PRIVATE
    "${WIREHAIR_SOURCE_DIR}")
target_compile_options(wirehair_v1_benchmark PRIVATE
    -Wall -Wextra -Wpedantic -Werror)
target_link_libraries(wirehair_v1_benchmark PRIVATE wirehair::wirehair)

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(NAME leopard2_atlas_python_unit
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_generate_atlas.py" -v)
    add_test(NAME leopard2_atlas_generator_self_test
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/generate_atlas.py" self-test)
    add_test(NAME leopard2_atlas_wirehair_adapter_smoke
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_wirehair_adapter.py"
            --binary "$<TARGET_FILE:wirehair_v1_benchmark>"
            --compile-commands
                "${CMAKE_BINARY_DIR}/compile_commands.json"
            --wirehair-source "${WIREHAIR_SOURCE_DIR}"
            --expected-commit "${WIREHAIR_EXPECTED_COMMIT}")

    set_tests_properties(
        leopard2_atlas_python_unit
        leopard2_atlas_generator_self_test
        leopard2_atlas_wirehair_adapter_smoke
        PROPERTIES
            PROCESSORS 1
            RUN_SERIAL TRUE
            TIMEOUT 120)
    set_tests_properties(leopard2_atlas_wirehair_adapter_smoke
        PROPERTIES SKIP_RETURN_CODE 77)
endif()
