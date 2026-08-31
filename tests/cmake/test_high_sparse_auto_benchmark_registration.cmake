# Prove that the sparse-high production-AUTO benchmark is available in a
# benchmark-only build, links the ordinary production archive, and reports
# no-hook route evidence for each supported runtime policy.

if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_C_COMPILER OR NOT DEFINED LEO2_CXX_COMPILER OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_EXECUTABLE_SUFFIX OR
   NOT DEFINED LEO2_ARCHIVE_SYMBOL_CHECKS OR
   NOT DEFINED LEO2_STATIC_LIBRARY_PREFIX OR
   NOT DEFINED LEO2_STATIC_LIBRARY_SUFFIX)
    message(FATAL_ERROR
        "sparse-high AUTO benchmark registration arguments are incomplete")
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
    -DLEOPARD_ENABLE_GF8=ON
    -DLEO2_BACKEND_VARIANT=auto
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE=ON
    -DLEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE_AUTO=OFF
    "${LEO2_SOURCE_DIR}")

execute_process(
    COMMAND ${configure_command}
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "sparse-high benchmark-only configure failed (${configure_result})\n"
        "stdout:\n${configure_stdout}\n"
        "stderr:\n${configure_stderr}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build .
        --target bench_leopard2_high_sparse_auto --config Release
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE build_result
    OUTPUT_VARIABLE build_stdout
    ERROR_VARIABLE build_stderr)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR
        "sparse-high benchmark-only target build failed (${build_result})\n"
        "stdout:\n${build_stdout}\n"
        "stderr:\n${build_stderr}")
endif()

file(GLOB_RECURSE hook_archive_candidates
    "${LEO2_BINARY_DIR}/*leopard_test_hooks*")
list(LENGTH hook_archive_candidates hook_archive_count)
if(NOT hook_archive_count EQUAL 0)
    message(FATAL_ERROR
        "tests-off sparse-high benchmark produced a hook archive")
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
    "${LEO2_BINARY_DIR}/bench_leopard2_high_sparse_auto${LEO2_EXECUTABLE_SUFFIX}")
list(LENGTH benchmark_candidates benchmark_count)
if(NOT benchmark_count EQUAL 1)
    message(FATAL_ERROR
        "expected one sparse-high AUTO benchmark, found ${benchmark_count}")
endif()
list(GET benchmark_candidates 0 benchmark)

function(assert_benchmark_json benchmark_json expected_policy expected_api
         expected_batch expected_rows expected_auto expected_route
         expected_public_executions expected_item_calls
         expected_direct_calls expected_transform_calls
         expect_binding_setup)
    if(expect_binding_setup)
        set(expected_binding_text "\"binding_setup\": {")
    else()
        set(expected_binding_text "\"binding_setup\": null")
    endif()
    foreach(required_text
        "\"schema\": \"leopard2-high-sparse-auto-benchmark-v1\""
        "\"authoritative\": false"
        "\"backend_variant\": \"auto\""
        "\"build_type\": \"Release\""
        "\"build_configuration_schema\": \"leopard2-benchmark-build-configuration/v16\""
        "\"library_test_hooks\": false"
        "\"high_sparse_tables_compiled\": true"
        "\"high_sparse_auto_compiled_default\": false"
        "\"requested_backend\": \"auto\""
        "\"policy\": \"${expected_policy}\""
        "\"api\": \"${expected_api}\""
        "\"batch\": ${expected_batch}"
        "\"requested_thread_count\": 1"
        "\"effective_backend\": \"avx2\""
        "\"thread_count\": 1"
        "\"direct_generator_rows\": ${expected_rows}"
        "\"auto_direct_selected\": ${expected_auto}"
        "\"selected_route\": \"${expected_route}\""
        "\"route_witness_armed\": true"
        "\"witness_public_executions\": ${expected_public_executions}"
        "\"expected_item_calls\": ${expected_item_calls}"
        "\"direct_calls\": ${expected_direct_calls}"
        "\"transform_calls\": ${expected_transform_calls}"
        "\"witness_disabled_before_timing\": true"
        "\"independent_oracle_match\": true"
        "\"input_immutable\": true"
        "\"unrequested_outputs_untouched\": true"
        "\"post_timing_recheck_match\": true"
        "\"codec_setup\": {"
        "${expected_binding_text}"
        "\"execution\": {"
        "\"amortized\": {")
        string(FIND "${benchmark_json}" "${required_text}" required_position)
        if(required_position EQUAL -1)
            message(FATAL_ERROR
                "sparse-high benchmark JSON omitted '${required_text}'")
        endif()
    endforeach()
endfunction()

set(common_arguments
    --k 2 --r 16 --bytes 1024 --parity-index 7
    --backend auto --threads 1
    --iterations 1 --setup-iterations 1 --calls-per-sample 1
    --warmups 0 --reuse 1 --memory-mib 16)

# Run the tables-off lane first so a non-AVX2 host can retain the structural
# registration proof without pretending to have executed an unsupported ISA.
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1
        "${benchmark}" ${common_arguments}
        --api one-shot --batch 1 --policy tables-off-auto-off
    RESULT_VARIABLE first_result
    OUTPUT_VARIABLE first_json
    ERROR_VARIABLE first_stderr)
if(NOT first_result EQUAL 0)
    string(STRIP "${first_stderr}" first_error)
    if(first_error STREQUAL
       "leopard2 sparse-high AUTO benchmark: effective backend must be avx2 for sparse-high telemetry")
        string(STRIP "${first_json}" first_output)
        if(NOT first_output STREQUAL "")
            message(FATAL_ERROR
                "unsupported AVX2 host emitted benchmark JSON")
        endif()
        message(STATUS
            "sparse-high AUTO runtime lanes skipped: host AVX2 unavailable")
        return()
    endif()
    message(FATAL_ERROR
        "tables-off sparse-high benchmark failed (${first_result})\n"
        "stdout:\n${first_json}\n"
        "stderr:\n${first_stderr}")
endif()
assert_benchmark_json(
    "${first_json}"
    tables_off_auto_off one_shot 1
    0 false transform
    1 1 0 1 false)

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1
        "${benchmark}" ${common_arguments}
        --api batch --batch 4 --policy tables-on-auto-off
    RESULT_VARIABLE second_result
    OUTPUT_VARIABLE second_json
    ERROR_VARIABLE second_stderr)
if(NOT second_result EQUAL 0)
    message(FATAL_ERROR
        "AUTO-off sparse-high batch failed (${second_result}): "
        "${second_stderr}")
endif()
assert_benchmark_json(
    "${second_json}"
    tables_on_auto_off batch 4
    16 false transform
    1 4 0 4 false)

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1
        "${benchmark}" ${common_arguments}
        --api binding --batch 4 --policy tables-on-auto-on
    RESULT_VARIABLE third_result
    OUTPUT_VARIABLE third_json
    ERROR_VARIABLE third_stderr)
if(NOT third_result EQUAL 0)
    message(FATAL_ERROR
        "AUTO-on sparse-high binding failed (${third_result}): "
        "${third_stderr}")
endif()
assert_benchmark_json(
    "${third_json}"
    tables_on_auto_on binding 4
    16 true direct
    2 8 8 0 true)
