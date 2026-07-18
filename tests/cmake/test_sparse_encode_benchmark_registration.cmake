# Prove that the sparse encoder benchmark exists when the test suite is
# disabled, links the production archive, and reports that identity through
# its fail-closed JSON schema.

if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_C_COMPILER OR NOT DEFINED LEO2_CXX_COMPILER OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_EXECUTABLE_SUFFIX OR
   NOT DEFINED LEO2_ARCHIVE_SYMBOL_CHECKS OR
   NOT DEFINED LEO2_STATIC_LIBRARY_PREFIX OR
   NOT DEFINED LEO2_STATIC_LIBRARY_SUFFIX)
    message(FATAL_ERROR
        "sparse benchmark registration test arguments are incomplete")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(configure_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND configure_command
    "-DCMAKE_C_COMPILER=${LEO2_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}"
    -DCMAKE_BUILD_TYPE=Release
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

# The JSON marker is deliberately fail-closed, but it is still compiled into
# the executable.  Where the host offers GNU-compatible nm, independently
# prove that the archive linked by this tests-off graph exports no hook ABI.
# On other generators the tests-off configure remains the portable structural
# proof: leopard_test_hooks is not defined and therefore cannot satisfy the
# benchmark's link dependency.
file(GLOB_RECURSE hook_archive_candidates
    "${LEO2_BINARY_DIR}/*leopard_test_hooks*")
list(LENGTH hook_archive_candidates hook_archive_count)
if(NOT hook_archive_count EQUAL 0)
    message(FATAL_ERROR
        "tests-off benchmark build unexpectedly produced a hook archive")
endif()
if(LEO2_ARCHIVE_SYMBOL_CHECKS)
    if(NOT DEFINED LEO2_NM_COMMAND OR LEO2_NM_COMMAND STREQUAL "")
        message(FATAL_ERROR "archive symbol checks require nm")
    endif()
    file(GLOB_RECURSE production_archive_candidates
        "${LEO2_BINARY_DIR}/${LEO2_STATIC_LIBRARY_PREFIX}leopard${LEO2_STATIC_LIBRARY_SUFFIX}")
    list(LENGTH production_archive_candidates production_archive_count)
    if(NOT production_archive_count EQUAL 1)
        message(FATAL_ERROR
            "expected one production Leopard archive, found "
            "${production_archive_count}")
    endif()
    list(GET production_archive_candidates 0 production_archive)
    execute_process(
        COMMAND "${LEO2_NM_COMMAND}" -g --defined-only "${production_archive}"
        RESULT_VARIABLE nm_result
        OUTPUT_VARIABLE archive_symbols
        ERROR_VARIABLE nm_stderr)
    if(NOT nm_result EQUAL 0)
        message(FATAL_ERROR
            "cannot inspect production archive (${nm_result}): ${nm_stderr}")
    endif()
    string(REGEX MATCH "TestOnly|leo2_test_" hook_symbol "${archive_symbols}")
    if(NOT hook_symbol STREQUAL "")
        message(FATAL_ERROR
            "production archive retained a test-hook symbol: ${hook_symbol}")
    endif()
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
        --profile high --field gf8 --k 48 --r 16 --bytes 64
        --requested-parity 0,7,15 --backend scalar --iterations 1
        --rounds 3 --warmups 1 --setup-iterations 1 --reuse 1,8,64
        --memory-mib 16
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE benchmark_json
    ERROR_VARIABLE run_stderr)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR
        "benchmark-only executable failed (${run_result}): ${run_stderr}")
endif()
foreach(required_text
    "\"schema\": \"leopard2-sparse-encode-benchmark-v4\""
    "\"library_test_hooks\": false"
    "\"exact_prefix_parity_match\": true"
    "\"direct_generator_parity_match\": true"
    "\"encode_decode_recovery_match\": true")
    string(FIND "${benchmark_json}" "${required_text}" required_position)
    if(required_position EQUAL -1)
        message(FATAL_ERROR
            "benchmark-only JSON omitted '${required_text}'")
    endif()
endforeach()
