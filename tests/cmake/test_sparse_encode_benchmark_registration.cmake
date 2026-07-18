# Prove that the sparse encoder benchmark remains buildable without the test
# suite, links an uninstrumented Leopard archive, and reports that identity in
# its fail-closed JSON schema.

if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_C_COMPILER OR NOT DEFINED LEO2_CXX_COMPILER OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_EXECUTABLE_SUFFIX)
    message(FATAL_ERROR
        "sparse benchmark registration test arguments are incomplete")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(configure_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(DEFINED LEO2_GENERATOR_PLATFORM AND NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND configure_command
    "-DCMAKE_C_COMPILER=${LEO2_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    -DLEO2_BUILD_TESTS=OFF
    -DLEO2_BUILD_BENCHMARKS=ON
    -DLEO2_ENABLE_CUDA=OFF
    "${LEO2_SOURCE_DIR}")

execute_process(
    COMMAND ${configure_command}
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "benchmark-only configure failed (${configure_result})\n"
        "stdout:\n${configure_stdout}\n"
        "stderr:\n${configure_stderr}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build .
        --target bench_leopard2_sparse_encode --config Release
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE build_result
    OUTPUT_VARIABLE build_stdout
    ERROR_VARIABLE build_stderr)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR
        "benchmark-only target build failed (${build_result})\n"
        "stdout:\n${build_stdout}\n"
        "stderr:\n${build_stderr}")
endif()

set(compile_database "${LEO2_BINARY_DIR}/compile_commands.json")
if(NOT EXISTS "${compile_database}")
    message(FATAL_ERROR "benchmark-only build omitted compile_commands.json")
endif()
file(READ "${compile_database}" compile_commands)
string(FIND "${compile_commands}" "LEO2_ENABLE_TEST_HOOKS" hook_marker)
string(FIND "${compile_commands}"
    "LEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS=0" identity_marker)
if(NOT hook_marker EQUAL -1 OR identity_marker EQUAL -1)
    message(FATAL_ERROR
        "benchmark-only compile identity is instrumented or unmarked")
endif()

file(GLOB_RECURSE benchmark_candidates
    "${LEO2_BINARY_DIR}/bench_leopard2_sparse_encode${LEO2_EXECUTABLE_SUFFIX}")
list(LENGTH benchmark_candidates benchmark_count)
if(NOT benchmark_count EQUAL 1)
    message(FATAL_ERROR
        "expected one sparse benchmark executable, found ${benchmark_count}")
endif()
list(GET benchmark_candidates 0 benchmark)
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1
        "${benchmark}"
        --profile low --field gf8 --k 17 --r 18 --bytes 64
        --requested-parity 0-17 --backend scalar --iterations 1
        --samples 3 --warmups 1 --setup-iterations 1 --reuse 1,8,64
        --memory-mib 16
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE benchmark_json
    ERROR_VARIABLE run_stderr)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR
        "benchmark-only executable failed (${run_result}): ${run_stderr}")
endif()
foreach(required_text
    "\"schema\": \"leopard2-sparse-encode-benchmark-v2\""
    "\"library_test_hooks\": false"
    "\"inverse_pruned_parity_match\": true")
    string(FIND "${benchmark_json}" "${required_text}" required_position)
    if(required_position EQUAL -1)
        message(FATAL_ERROR
            "benchmark-only JSON omitted '${required_text}'")
    endif()
endforeach()
