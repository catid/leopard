if(NOT DEFINED LEO2_SOURCE_DIR OR LEO2_SOURCE_DIR STREQUAL "")
    message(FATAL_ERROR "LEO2_SOURCE_DIR is required")
endif()
if(NOT DEFINED LEO2_BINARY_DIR OR LEO2_BINARY_DIR STREQUAL "")
    message(FATAL_ERROR "LEO2_BINARY_DIR is required")
endif()
if(NOT DEFINED LEO2_GENERATOR OR LEO2_GENERATOR STREQUAL "")
    message(FATAL_ERROR "LEO2_GENERATOR is required")
endif()

find_program(LEO2_TEST_GIT_EXECUTABLE NAMES git)
if(NOT LEO2_TEST_GIT_EXECUTABLE)
    message(STATUS
        "Skipping benchmark attestation refresh regression: Git unavailable")
    return()
endif()

function(leo2_run)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "Command failed (${result}): ${ARGN}\n"
            "stdout:\n${output}\n"
            "stderr:\n${error}")
    endif()
    set(LEO2_LAST_STDOUT "${output}" PARENT_SCOPE)
    set(LEO2_LAST_STDERR "${error}" PARENT_SCOPE)
endfunction()

function(leo2_git)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${fixture_source}" ${ARGN})
    set(LEO2_LAST_STDOUT "${LEO2_LAST_STDOUT}" PARENT_SCOPE)
endfunction()

function(leo2_expected_identity commit_out tree_out dirty_out)
    leo2_git(rev-parse --verify HEAD)
    string(STRIP "${LEO2_LAST_STDOUT}" commit)
    leo2_git(rev-parse --verify HEAD^{tree})
    string(STRIP "${LEO2_LAST_STDOUT}" tree)
    leo2_git(status --porcelain=v1 --untracked-files=normal
        --ignore-submodules=none)
    if(LEO2_LAST_STDOUT STREQUAL "")
        set(dirty 0)
    else()
        set(dirty 1)
    endif()
    set(${commit_out} "${commit}" PARENT_SCOPE)
    set(${tree_out} "${tree}" PARENT_SCOPE)
    set(${dirty_out} "${dirty}" PARENT_SCOPE)
endfunction()

function(leo2_read_effective_configuration
        build_directory digest_out build_type_out material_out)
    set(configuration_file
        "${build_directory}/generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt")
    if(NOT EXISTS "${configuration_file}")
        message(FATAL_ERROR
            "Missing effective build configuration: ${configuration_file}")
    endif()
    file(READ "${configuration_file}" configuration)
    string(REGEX MATCH
        "^schema=leopard2-benchmark-build-configuration/v10\nsha256=([0-9a-f]+)\n"
        header "${configuration}")
    set(declared_digest "${CMAKE_MATCH_1}")
    string(LENGTH "${declared_digest}" declared_digest_length)
    if(NOT header OR NOT declared_digest MATCHES "^[0-9a-f]+$" OR
       NOT declared_digest_length EQUAL 64)
        message(FATAL_ERROR
            "Invalid effective build-configuration header:\n${configuration}")
    endif()
    string(LENGTH "${header}" header_length)
    string(SUBSTRING "${configuration}" "${header_length}" -1 material)
    string(SHA256 actual_digest "${material}")
    if(NOT actual_digest STREQUAL declared_digest)
        message(FATAL_ERROR
            "Effective build-configuration digest mismatch: "
            "${actual_digest} != ${declared_digest}\n"
            "material:\n${material}")
    endif()
    file(STRINGS "${configuration_file}" build_type_lines
        REGEX "^CMAKE_BUILD_TYPE=")
    file(STRINGS "${configuration_file}" direct_source_plan_lines
        REGEX "^LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=")
    file(STRINGS "${configuration_file}" high_t16_generated_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=")
    file(STRINGS "${configuration_file}" high_t16_q2_fused_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=")
    file(STRINGS "${configuration_file}" high_t32_two_block_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=")
    file(STRINGS "${configuration_file}" high_t32_disable_two_block_lines
        REGEX "^LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=")
    file(STRINGS "${configuration_file}" low_p32_terminal_lines
        REGEX "^LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=")
    file(STRINGS "${configuration_file}" one_shot_equal_rounded_lines
        REGEX "^LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=")
    file(STRINGS "${configuration_file}" cauchy_log_reuse_lines
        REGEX "^LEO2_EXPERIMENT_CAUCHY_LOG_REUSE=")
    file(STRINGS "${configuration_file}" high_direct_encode_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=")
    file(STRINGS "${configuration_file}" high_t8_disable_vector_lines
        REGEX "^LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=")
    file(STRINGS "${configuration_file}" high_t8_partial_binding_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=")
    file(STRINGS "${configuration_file}" high_t8_two_block_binding_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=")
    file(STRINGS "${configuration_file}" high_t8_ragged_binding_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=")
    file(STRINGS "${configuration_file}" high_t32_generated_lines
        REGEX "^LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=")
    file(STRINGS "${configuration_file}" high_t32_disable_generated_lines
        REGEX "^LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=")
    file(STRINGS "${configuration_file}" general_one_loss_direct_lines
        REGEX "^LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=")
    file(STRINGS "${configuration_file}" small_direct_mode_lines
        REGEX "^LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=")
    file(STRINGS "${configuration_file}" small_dual_direct_lines
        REGEX "^LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=")
    file(STRINGS "${configuration_file}" small_dual_locator_lines
        REGEX "^LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=")
    file(STRINGS "${configuration_file}" small_dual_fallback_lines
        REGEX "^LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=")
    list(LENGTH build_type_lines build_type_line_count)
    list(LENGTH direct_source_plan_lines direct_source_plan_line_count)
    list(LENGTH high_t16_generated_lines high_t16_generated_line_count)
    list(LENGTH high_t16_q2_fused_lines high_t16_q2_fused_line_count)
    list(LENGTH high_t32_two_block_lines high_t32_two_block_line_count)
    list(LENGTH high_t32_disable_two_block_lines
        high_t32_disable_two_block_line_count)
    list(LENGTH low_p32_terminal_lines low_p32_terminal_line_count)
    list(LENGTH one_shot_equal_rounded_lines
        one_shot_equal_rounded_line_count)
    list(LENGTH cauchy_log_reuse_lines cauchy_log_reuse_line_count)
    list(LENGTH high_direct_encode_lines high_direct_encode_line_count)
    list(LENGTH high_t8_disable_vector_lines high_t8_disable_vector_line_count)
    list(LENGTH high_t8_partial_binding_lines
        high_t8_partial_binding_line_count)
    list(LENGTH high_t8_two_block_binding_lines
        high_t8_two_block_binding_line_count)
    list(LENGTH high_t8_ragged_binding_lines
        high_t8_ragged_binding_line_count)
    list(LENGTH high_t32_generated_lines high_t32_generated_line_count)
    list(LENGTH high_t32_disable_generated_lines
        high_t32_disable_generated_line_count)
    list(LENGTH general_one_loss_direct_lines
        general_one_loss_direct_line_count)
    list(LENGTH small_direct_mode_lines small_direct_mode_line_count)
    list(LENGTH small_dual_direct_lines small_dual_direct_line_count)
    list(LENGTH small_dual_locator_lines small_dual_locator_line_count)
    list(LENGTH small_dual_fallback_lines small_dual_fallback_line_count)
    if(NOT build_type_line_count EQUAL 1 OR
       NOT direct_source_plan_line_count EQUAL 1 OR
       NOT high_t16_generated_line_count EQUAL 1 OR
       NOT high_t16_q2_fused_line_count EQUAL 1 OR
       NOT high_t32_two_block_line_count EQUAL 1 OR
       NOT high_t32_disable_two_block_line_count EQUAL 1 OR
       NOT low_p32_terminal_line_count EQUAL 1 OR
       NOT one_shot_equal_rounded_line_count EQUAL 1 OR
       NOT cauchy_log_reuse_line_count EQUAL 1 OR
       NOT high_direct_encode_line_count EQUAL 1 OR
       NOT high_t8_disable_vector_line_count EQUAL 1 OR
       NOT high_t8_partial_binding_line_count EQUAL 1 OR
       NOT high_t8_two_block_binding_line_count EQUAL 1 OR
       NOT high_t8_ragged_binding_line_count EQUAL 1 OR
       NOT high_t32_generated_line_count EQUAL 1 OR
       NOT high_t32_disable_generated_line_count EQUAL 1 OR
       NOT general_one_loss_direct_line_count EQUAL 1 OR
       NOT small_direct_mode_line_count EQUAL 1 OR
       NOT small_dual_direct_line_count EQUAL 1 OR
       NOT small_dual_locator_line_count EQUAL 1 OR
       NOT small_dual_fallback_line_count EQUAL 1)
        message(FATAL_ERROR
            "Effective build configuration omits CMAKE_BUILD_TYPE or has "
            "an invalid direct-selector record count")
    endif()
    list(GET build_type_lines 0 effective_build_type)
    string(REGEX REPLACE "^CMAKE_BUILD_TYPE=" ""
        effective_build_type "${effective_build_type}")
    list(GET small_direct_mode_lines 0 small_direct_mode)
    string(REGEX REPLACE
        "^LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=" ""
        small_direct_mode "${small_direct_mode}")
    if(NOT small_direct_mode MATCHES "^[012]$")
        message(FATAL_ERROR
            "Effective build configuration has an invalid small-direct mode")
    endif()
    list(GET small_dual_direct_lines 0 small_dual_direct)
    string(REGEX REPLACE
        "^LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=" ""
        small_dual_direct "${small_dual_direct}")
    if(NOT small_dual_direct MATCHES "^(ON|OFF)$")
        message(FATAL_ERROR
            "Effective build configuration has an invalid small "
            "dual-direct selector")
    endif()
    foreach(selector_lines IN ITEMS
            small_dual_locator_lines small_dual_fallback_lines)
        list(GET ${selector_lines} 0 selector_value)
        string(REGEX REPLACE "^[^=]+=" "" selector_value "${selector_value}")
        if(NOT selector_value MATCHES "^(ON|OFF)$")
            message(FATAL_ERROR
                "Effective build configuration has an invalid small-dual "
                "experiment selector")
        endif()
    endforeach()

    file(STRINGS "${build_directory}/CMakeCache.txt" cached_digest_lines
        REGEX
        "^LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:INTERNAL=")
    file(STRINGS "${build_directory}/CMakeCache.txt" cached_schema_lines
        REGEX
        "^LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:INTERNAL=")
    list(LENGTH cached_digest_lines digest_line_count)
    list(LENGTH cached_schema_lines schema_line_count)
    if(NOT digest_line_count EQUAL 1 OR NOT schema_line_count EQUAL 1)
        message(FATAL_ERROR
            "Effective build-configuration cache binding is incomplete")
    endif()
    list(GET cached_digest_lines 0 cached_digest)
    list(GET cached_schema_lines 0 cached_schema)
    string(REGEX REPLACE "^[^=]+=" "" cached_digest "${cached_digest}")
    string(REGEX REPLACE "^[^=]+=" "" cached_schema "${cached_schema}")
    if(NOT cached_digest STREQUAL declared_digest OR
       NOT cached_schema STREQUAL
           "leopard2-benchmark-build-configuration/v10")
        message(FATAL_ERROR
            "Effective build configuration differs from its cache binding")
    endif()
    set(${digest_out} "${declared_digest}" PARENT_SCOPE)
    set(${build_type_out} "${effective_build_type}" PARENT_SCOPE)
    set(${material_out} "${material}" PARENT_SCOPE)
endfunction()

function(leo2_check_probe
        executable build_directory
        expected_commit expected_tree expected_dirty)
    leo2_read_effective_configuration(
        "${build_directory}" expected_configuration expected_build_type
        ignored_material)
    if(ARGC GREATER 5)
        set(expected_build_type "${ARGV5}")
    endif()
    leo2_run("${executable}")
    string(REPLACE "\r\n" "\n" probe "${LEO2_LAST_STDOUT}")
    set(expected
        "commit=${expected_commit}\ntree=${expected_tree}\ndirty=${expected_dirty}\n")
    string(APPEND expected
        "build_type=${expected_build_type}\n"
        "build_configuration=${expected_configuration}\n")
    if(NOT probe STREQUAL expected)
        message(FATAL_ERROR
            "Unexpected attestation from ${executable}\n"
            "expected:\n${expected}\nactual:\n${probe}")
    endif()
endfunction()

function(leo2_build_and_check)
    leo2_run(${fixture_build_command})
    leo2_expected_identity(expected_commit expected_tree expected_dirty)
    leo2_check_probe("${standard_executable}"
        "${fixture_build}"
        "${expected_commit}" "${expected_tree}" "${expected_dirty}")
    leo2_check_probe("${allk_executable}"
        "${fixture_build}"
        "${expected_commit}" "${expected_tree}" "${expected_dirty}")
endfunction()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(fixture_source "${LEO2_BINARY_DIR}/source")
set(fixture_build "${LEO2_BINARY_DIR}/build")
file(MAKE_DIRECTORY "${fixture_source}")
file(MAKE_DIRECTORY "${fixture_build}")

file(TO_CMAKE_PATH
    "${LEO2_SOURCE_DIR}/cmake/Leopard2BenchmarkAttestation.cmake"
    attestation_module)
file(WRITE "${fixture_source}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.7)\n"
"project(leopard2_attestation_fixture CXX)\n"
"include(\"${attestation_module}\")\n"
"set(LEO2_TEST_EFFECTIVE_SUFFIX \"ONE\" CACHE STRING \"\")\n"
"if(MSVC)\n"
"  set(CMAKE_CXX_FLAGS \"\${CMAKE_CXX_FLAGS} /DLEO2_EFFECTIVE_\${LEO2_TEST_EFFECTIVE_SUFFIX}=1\")\n"
"else()\n"
"  set(CMAKE_CXX_FLAGS \"\${CMAKE_CXX_FLAGS} -DLEO2_EFFECTIVE_\${LEO2_TEST_EFFECTIVE_SUFFIX}=1\")\n"
"endif()\n"
"add_executable(bench_leopard2 probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2)\n"
"add_executable(bench_leopard2_allk probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2_allk)\n")
file(WRITE "${fixture_source}/probe.cpp"
    "#include <iostream>\n"
    "#if !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER)\n"
    "#error \"missing exact generated attestation header path\"\n"
    "#endif\n"
    "#if !defined(LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256) || !defined(LEO2_BENCHMARK_BUILD_TYPE)\n"
    "#error \"missing effective build-configuration attestation\"\n"
    "#endif\n"
    "#include LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER\n"
    "int main() {\n"
"  std::cout << \"commit=\" << LEO2_BENCHMARK_SOURCE_COMMIT << \"\\n\"\n"
"            << \"tree=\" << LEO2_BENCHMARK_SOURCE_TREE << \"\\n\"\n"
"            << \"dirty=\" << LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY << \"\\n\"\n"
"            << \"build_type=\" << LEO2_BENCHMARK_BUILD_TYPE << \"\\n\"\n"
"            << \"build_configuration=\" << LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256 << \"\\n\";\n"
    "  return 0;\n"
    "}\n")
file(WRITE "${fixture_source}/identity_input.txt" "initial\n")
file(WRITE "${fixture_source}/.gitignore"
    "/leopard2_benchmark_source_attestation.h\n")

leo2_git(init)
leo2_git(config user.name "Leopard2 attestation regression")
leo2_git(config user.email "leopard2-attestation@example.invalid")
leo2_git(config commit.gpgsign false)
leo2_git(add CMakeLists.txt probe.cpp identity_input.txt .gitignore)
leo2_git(commit -m initial)

# A source-local ignored header must not shadow the absolute generated compile
# input.  The fixture remains Git-clean, and compilation fails immediately if
# quoted-include lookup is accidentally restored.
file(WRITE
    "${fixture_source}/leopard2_benchmark_source_attestation.h"
    "#error \"source-local ignored attestation header was consumed\"\n")

set(configure_command
    "${CMAKE_COMMAND}" -E chdir "${fixture_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}" -DCMAKE_BUILD_TYPE=Release)
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND configure_command "${fixture_source}")
leo2_run(${configure_command})

set(fixture_build_command "${CMAKE_COMMAND}" --build "${fixture_build}")
set(fixture_output_directory "${fixture_build}")
if(LEO2_GENERATOR MATCHES "Visual Studio|Xcode|Multi-Config")
    list(APPEND fixture_build_command --config Release)
    set(fixture_output_directory "${fixture_build}/Release")
endif()
set(standard_executable
    "${fixture_output_directory}/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}")
set(allk_executable
    "${fixture_output_directory}/bench_leopard2_allk${LEO2_EXECUTABLE_SUFFIX}")
leo2_build_and_check()
leo2_read_effective_configuration(
    "${fixture_build}" initial_configuration initial_build_type
    initial_configuration_material)
if(NOT initial_build_type STREQUAL "Release" OR
   NOT initial_configuration_material MATCHES
       "CMAKE_CXX_FLAGS=[^\n]*LEO2_EFFECTIVE_ONE=1" OR
   NOT initial_configuration_material MATCHES
       "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=\nLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=\nLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=\nLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=\nLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=" OR
   NOT initial_configuration_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "Initial sidecar omitted the effective non-cache CXX flag:\n"
        "${initial_configuration_material}")
endif()

# These selectors change distinct production translation units or routing.
# Each one must independently perturb the effective-configuration digest; an
# omitted selector could otherwise let two different benchmark binaries share
# one attested build identity.
foreach(selector
        LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED
        LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
        LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK
        LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK
        LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL)
    leo2_run(${configure_command} "-D${selector}=ON")
    leo2_read_effective_configuration(
        "${fixture_build}" selector_configuration selector_build_type
        selector_material)
    if(selector_configuration STREQUAL initial_configuration OR
       NOT selector_build_type STREQUAL "Release" OR
       NOT selector_material MATCHES "${selector}=ON\n")
        message(FATAL_ERROR
            "${selector} did not change the attested effective "
            "configuration")
    endif()
    leo2_run(${configure_command} "-D${selector}=")
    leo2_read_effective_configuration(
        "${fixture_build}" restored_configuration restored_build_type
        restored_material)
    if(NOT restored_configuration STREQUAL initial_configuration OR
       NOT restored_build_type STREQUAL initial_build_type OR
       NOT restored_material STREQUAL initial_configuration_material)
        message(FATAL_ERROR
            "Clearing ${selector} did not restore the original effective "
            "configuration")
    endif()
endforeach()

# The default-on GF8 dual-direct route is a production code selector.  Its
# spelling is deliberately exact so CMake truth aliases cannot produce a
# different binary under the same attested material.
leo2_run(${configure_command} -DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=OFF)
leo2_read_effective_configuration(
    "${fixture_build}" dual_direct_off_configuration
    dual_direct_off_build_type dual_direct_off_material)
if(dual_direct_off_configuration STREQUAL initial_configuration OR
   NOT dual_direct_off_build_type STREQUAL initial_build_type OR
   NOT dual_direct_off_material MATCHES
       "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=OFF\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "GF8 small dual-direct OFF state was not exactly attested")
endif()
leo2_run(${configure_command} -DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON)
leo2_read_effective_configuration(
    "${fixture_build}" dual_direct_restored_configuration
    dual_direct_restored_build_type dual_direct_restored_material)
if(NOT dual_direct_restored_configuration STREQUAL initial_configuration OR
   NOT dual_direct_restored_build_type STREQUAL initial_build_type OR
   NOT dual_direct_restored_material STREQUAL initial_configuration_material)
    message(FATAL_ERROR
        "Restoring GF8 small dual-direct ON did not restore the original "
        "effective configuration")
endif()

# Both small-dual subselectors change leopard2.cpp and must independently alter
# the effective configuration.  Restore each canonical default before the next
# fixture phase so later digest comparisons remain stable.
foreach(selector IN ITEMS
        LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS
        LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK)
    if(selector STREQUAL "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS")
        set(changed_value OFF)
        set(default_value ON)
    else()
        set(changed_value ON)
        set(default_value OFF)
    endif()
    leo2_run(${configure_command} "-D${selector}=${changed_value}")
    leo2_read_effective_configuration(
        "${fixture_build}" changed_configuration changed_build_type
        changed_material)
    if(changed_configuration STREQUAL initial_configuration OR
       NOT changed_build_type STREQUAL initial_build_type OR
       NOT changed_material MATCHES "${selector}=${changed_value}\n")
        message(FATAL_ERROR
            "${selector} did not change the attested effective configuration")
    endif()
    leo2_run(${configure_command} "-D${selector}=${default_value}")
    leo2_read_effective_configuration(
        "${fixture_build}" restored_configuration restored_build_type
        restored_material)
    if(NOT restored_configuration STREQUAL initial_configuration OR
       NOT restored_build_type STREQUAL initial_build_type OR
       NOT restored_material STREQUAL initial_configuration_material)
        message(FATAL_ERROR
            "Restoring ${selector}=${default_value} did not restore the "
            "initial effective configuration")
    endif()
endforeach()

set(attestation_header
    "${fixture_build}/generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h")
set(build_configuration_sidecar
    "${fixture_build}/generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt")
if(NOT EXISTS "${attestation_header}")
    message(FATAL_ERROR "Fixture did not generate its attestation header")
endif()
if(NOT EXISTS "${build_configuration_sidecar}")
    message(FATAL_ERROR
        "Fixture did not generate its effective-configuration sidecar")
endif()
file(GLOB build_configuration_temporaries
    "${build_configuration_sidecar}.*.tmp")
if(build_configuration_temporaries)
    message(FATAL_ERROR
        "Effective-configuration publication left temporary files: "
        "${build_configuration_temporaries}")
endif()

# The selector probes above intentionally rewrite target compile definitions
# without rebuilding between every configure.  Settle the restored graph once
# before recording the no-op baseline; otherwise a correctly stale object can
# be mistaken for needless work by the following unchanged refresh.
leo2_build_and_check()

# A no-op reconfigure and unchanged refresh must not touch either generated
# identity input or relink either benchmark.  Exact-main freshness requires
# the benchmark object to remain at least as new as this sidecar.
file(TIMESTAMP "${attestation_header}" stable_header_time "%s" UTC)
file(TIMESTAMP "${build_configuration_sidecar}"
    stable_configuration_time "%s" UTC)
file(TIMESTAMP "${standard_executable}" stable_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" stable_allk_time "%s" UTC)
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command})
leo2_build_and_check()
file(TIMESTAMP "${attestation_header}" repeat_header_time "%s" UTC)
file(TIMESTAMP "${build_configuration_sidecar}"
    repeat_configuration_time "%s" UTC)
file(TIMESTAMP "${standard_executable}" repeat_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" repeat_allk_time "%s" UTC)
if(NOT stable_header_time STREQUAL repeat_header_time OR
   NOT stable_configuration_time STREQUAL repeat_configuration_time OR
   NOT stable_standard_time STREQUAL repeat_standard_time OR
   NOT stable_allk_time STREQUAL repeat_allk_time)
    message(FATAL_ERROR
        "Unchanged source identity caused needless benchmark recompilation")
endif()

# Reconfiguration must update the sidecar, embedded digest, cache binding, and
# both benchmark objects when an effective non-cache flag changes.
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command} -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO)
leo2_read_effective_configuration(
    "${fixture_build}" changed_configuration changed_build_type
    changed_configuration_material)
if(changed_configuration STREQUAL initial_configuration OR
   NOT changed_build_type STREQUAL "Release" OR
   NOT changed_configuration_material MATCHES
       "CMAKE_CXX_FLAGS=[^\n]*LEO2_EFFECTIVE_TWO=1" OR
   changed_configuration_material MATCHES "LEO2_EFFECTIVE_ONE=1")
    message(FATAL_ERROR
        "Reconfiguration did not update the exact effective configuration")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}" changed_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" changed_allk_time "%s" UTC)
if(changed_standard_time STREQUAL repeat_standard_time OR
   changed_allk_time STREQUAL repeat_allk_time)
    message(FATAL_ERROR
        "Effective build-configuration change did not relink benchmarks")
endif()

# Each Boolean codec experiment changes production code and therefore
# must independently change the attested digest and relink both benchmarks.
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=ON
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" direct_source_configuration
    direct_source_build_type direct_source_material)
if(direct_source_configuration STREQUAL changed_configuration OR
   NOT direct_source_build_type STREQUAL "Release" OR
   NOT direct_source_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=ON\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "Direct-source-plan selector was not uniquely and positionally "
        "attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    direct_source_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    direct_source_allk_time "%s" UTC)
if(direct_source_standard_time STREQUAL changed_standard_time OR
   direct_source_allk_time STREQUAL changed_allk_time)
    message(FATAL_ERROR
        "Direct-source-plan change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=ON
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" high_direct_configuration
    high_direct_build_type high_direct_material)
if(high_direct_configuration STREQUAL changed_configuration OR
   high_direct_configuration STREQUAL direct_source_configuration OR
   NOT high_direct_build_type STREQUAL "Release" OR
   NOT high_direct_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=ON\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "High-direct-encode selector was not uniquely and positionally "
        "attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_direct_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_direct_allk_time "%s" UTC)
if(high_direct_standard_time STREQUAL direct_source_standard_time OR
   high_direct_allk_time STREQUAL direct_source_allk_time)
    message(FATAL_ERROR
        "High-direct-encode change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=ON
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" high_t8_disable_configuration
    high_t8_disable_build_type high_t8_disable_material)
if(high_t8_disable_configuration STREQUAL changed_configuration OR
   high_t8_disable_configuration STREQUAL direct_source_configuration OR
   high_t8_disable_configuration STREQUAL high_direct_configuration OR
   NOT high_t8_disable_build_type STREQUAL "Release" OR
   NOT high_t8_disable_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=ON\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "High-T8-vector diagnostic-disable selector was not uniquely and "
        "positionally attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_t8_disable_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_t8_disable_allk_time "%s" UTC)
if(high_t8_disable_standard_time STREQUAL high_direct_standard_time OR
   high_t8_disable_allk_time STREQUAL high_direct_allk_time)
    message(FATAL_ERROR
        "High-T8-vector diagnostic-disable change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=ON
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" high_t8_partial_configuration
    high_t8_partial_build_type high_t8_partial_material)
if(high_t8_partial_configuration STREQUAL changed_configuration OR
   high_t8_partial_configuration STREQUAL direct_source_configuration OR
   high_t8_partial_configuration STREQUAL high_direct_configuration OR
   high_t8_partial_configuration STREQUAL high_t8_disable_configuration OR
   NOT high_t8_partial_build_type STREQUAL "Release" OR
   NOT high_t8_partial_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=ON\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "High-T8 partial-binding selector was not uniquely and positionally "
        "attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_t8_partial_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_t8_partial_allk_time "%s" UTC)
if(high_t8_partial_standard_time STREQUAL high_t8_disable_standard_time OR
   high_t8_partial_allk_time STREQUAL high_t8_disable_allk_time)
    message(FATAL_ERROR
        "High-T8 partial-binding change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=ON
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=ON
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" high_t8_two_block_configuration
    high_t8_two_block_build_type high_t8_two_block_material)
if(high_t8_two_block_configuration STREQUAL changed_configuration OR
   high_t8_two_block_configuration STREQUAL
       direct_source_configuration OR
   high_t8_two_block_configuration STREQUAL high_direct_configuration OR
   high_t8_two_block_configuration STREQUAL
       high_t8_disable_configuration OR
   high_t8_two_block_configuration STREQUAL
       high_t8_partial_configuration OR
   NOT high_t8_two_block_build_type STREQUAL "Release" OR
   NOT high_t8_two_block_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=ON\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=ON\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "High-T8 two-block-binding selector was not uniquely and "
        "positionally attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_t8_two_block_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_t8_two_block_allk_time "%s" UTC)
if(high_t8_two_block_standard_time STREQUAL
       high_t8_partial_standard_time OR
   high_t8_two_block_allk_time STREQUAL high_t8_partial_allk_time)
    message(FATAL_ERROR
        "High-T8 two-block-binding change did not relink benchmarks")
endif()

# The generated T32 arithmetic selector and its layout-stable diagnostic
# control both affect benchmark interpretation.  Attest them independently so
# an ON candidate cannot be mistaken for its same-text OFF control.
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" high_t32_configuration high_t32_build_type
    high_t32_material)
if(high_t32_configuration STREQUAL high_t8_two_block_configuration OR
   NOT high_t32_build_type STREQUAL "Release" OR
   NOT high_t32_material MATCHES
       "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\n")
    message(FATAL_ERROR
        "High-T32 generated selector was not independently attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_t32_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_t32_allk_time "%s" UTC)
if(high_t32_standard_time STREQUAL high_t8_two_block_standard_time OR
   high_t32_allk_time STREQUAL high_t8_two_block_allk_time)
    message(FATAL_ERROR
        "High-T32 generated selector change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=ON)
leo2_read_effective_configuration(
    "${fixture_build}" high_t32_disabled_configuration
    high_t32_disabled_build_type high_t32_disabled_material)
if(high_t32_disabled_configuration STREQUAL high_t32_configuration OR
   NOT high_t32_disabled_build_type STREQUAL "Release" OR
   NOT high_t32_disabled_material MATCHES
       "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=ON\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=ON\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\n")
    message(FATAL_ERROR
        "High-T32 diagnostic-disable selector was not independently attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    high_t32_disabled_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    high_t32_disabled_allk_time "%s" UTC)
if(high_t32_disabled_standard_time STREQUAL high_t32_standard_time OR
   high_t32_disabled_allk_time STREQUAL high_t32_allk_time)
    message(FATAL_ERROR
        "High-T32 diagnostic-disable change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=ON
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0)
leo2_read_effective_configuration(
    "${fixture_build}" general_one_loss_configuration
    general_one_loss_build_type general_one_loss_material)
if(general_one_loss_configuration STREQUAL changed_configuration OR
   general_one_loss_configuration STREQUAL direct_source_configuration OR
   general_one_loss_configuration STREQUAL high_direct_configuration OR
   general_one_loss_configuration STREQUAL high_t8_disable_configuration OR
   general_one_loss_configuration STREQUAL high_t8_partial_configuration OR
   general_one_loss_configuration STREQUAL
       high_t8_two_block_configuration OR
   NOT general_one_loss_build_type STREQUAL "Release" OR
   NOT general_one_loss_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=ON\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "General one-loss selector was not uniquely and positionally "
        "attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}"
    general_one_loss_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}"
    general_one_loss_allk_time "%s" UTC)
if(general_one_loss_standard_time STREQUAL
       high_t8_two_block_standard_time OR
   general_one_loss_allk_time STREQUAL high_t8_two_block_allk_time)
    message(FATAL_ERROR
        "General one-loss selector change did not relink benchmarks")
endif()

# The default-off direct-repair selectors are part of the code identity.  Keep
# every prior material record in order and append the versioned selector records so
# small-direct modes 0/1/2 differ only where intended.
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1)
leo2_read_effective_configuration(
    "${fixture_build}" mode_one_configuration mode_one_build_type
    mode_one_configuration_material)
if(mode_one_configuration STREQUAL changed_configuration OR
   NOT mode_one_build_type STREQUAL "Release" OR
   NOT mode_one_configuration_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "Small-direct mode 1 was not uniquely and positionally attested")
endif()
leo2_build_and_check()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF
    -DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF
    -DLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF
    -DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF
    -DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF
    -DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2)
leo2_read_effective_configuration(
    "${fixture_build}" mode_two_configuration mode_two_build_type
    mode_two_configuration_material)
if(mode_two_configuration STREQUAL changed_configuration OR
   mode_two_configuration STREQUAL mode_one_configuration OR
   NOT mode_two_build_type STREQUAL "Release" OR
   NOT mode_two_configuration_material MATCHES
       "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=OFF\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=OFF\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=OFF\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=OFF\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "Small-direct mode 2 was not uniquely and positionally attested")
endif()
leo2_build_and_check()

# The two production direct-repair setup selectors are independent pieces of
# the wire-identical code identity and must each perturb the attested digest.
file(TIMESTAMP "${standard_executable}" mode_two_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" mode_two_allk_time "%s" UTC)
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=ON
    -DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=OFF)
leo2_read_effective_configuration(
    "${fixture_build}" one_shot_configuration one_shot_build_type
    one_shot_material)
if(one_shot_configuration STREQUAL mode_two_configuration OR
   NOT one_shot_build_type STREQUAL "Release" OR
   NOT one_shot_material MATCHES
       "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=ON\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=OFF\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "One-shot equal-rounded selector was not independently attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}" one_shot_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" one_shot_allk_time "%s" UTC)
if(one_shot_standard_time STREQUAL mode_two_standard_time OR
   one_shot_allk_time STREQUAL mode_two_allk_time)
    message(FATAL_ERROR
        "One-shot equal-rounded selector change did not relink benchmarks")
endif()

leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_run(${configure_command}
    -DLEO2_TEST_EFFECTIVE_SUFFIX=TWO
    -DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=OFF
    -DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=ON)
leo2_read_effective_configuration(
    "${fixture_build}" cauchy_configuration cauchy_build_type
    cauchy_material)
if(cauchy_configuration STREQUAL one_shot_configuration OR
   NOT cauchy_build_type STREQUAL "Release" OR
   NOT cauchy_material MATCHES
       "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=OFF\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=ON\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
    message(FATAL_ERROR
        "Cauchy-log-reuse selector was not independently attested")
endif()
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}" cauchy_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" cauchy_allk_time "%s" UTC)
if(cauchy_standard_time STREQUAL one_shot_standard_time OR
   cauchy_allk_time STREQUAL one_shot_allk_time)
    message(FATAL_ERROR
        "Cauchy-log-reuse selector change did not relink benchmarks")
endif()

# Unstaged and staged tracked changes are dirty.
file(APPEND "${fixture_source}/identity_input.txt" "unstaged\n")
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}" dirty_standard_time "%s" UTC)
leo2_git(add identity_input.txt)
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_build_and_check()
file(TIMESTAMP "${standard_executable}" staged_standard_time "%s" UTC)
if(NOT dirty_standard_time STREQUAL staged_standard_time)
    message(FATAL_ERROR
        "Equivalent dirty identity caused needless recompilation after staging")
endif()

# A commit changes both committed identities and returns to clean.
leo2_git(commit -m staged-change)
leo2_build_and_check()

# A tracked edit and a plain Git revert must both refresh without CMake
# reconfiguration.
file(APPEND "${fixture_source}/identity_input.txt" "reverted\n")
leo2_build_and_check()
leo2_git(checkout -- identity_input.txt)
leo2_build_and_check()

# Relevant untracked source participates in the established dirty contract.
file(WRITE "${fixture_source}/new_source.cpp" "int untracked_source;\n")
leo2_build_and_check()
file(REMOVE "${fixture_source}/new_source.cpp")
leo2_build_and_check()

# A source archive may be unpacked below an unrelated repository.  Git walks
# parent directories, but that enclosing identity must not be attributed to
# the archive source root.
set(outer_repository "${LEO2_BINARY_DIR}/outer-repository")
set(archive_source "${outer_repository}/archive-source")
set(archive_build "${LEO2_BINARY_DIR}/archive-build")
file(MAKE_DIRECTORY "${archive_source}")
file(MAKE_DIRECTORY "${archive_build}")
file(WRITE "${outer_repository}/.gitignore" "/archive-source/\n")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}" init)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}"
    config user.name "Leopard2 outer-repository regression")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}"
    config user.email "leopard2-outer@example.invalid")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}"
    config commit.gpgsign false)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}"
    add .gitignore)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${outer_repository}"
    commit -m outer-root)
file(WRITE "${archive_source}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.7)\n"
"project(leopard2_nested_archive_fixture CXX)\n"
"include(\"${attestation_module}\")\n"
"add_executable(bench_leopard2 probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
file(COPY "${fixture_source}/probe.cpp" DESTINATION "${archive_source}")

set(archive_configure_command
    "${CMAKE_COMMAND}" -E chdir "${archive_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}" -DCMAKE_BUILD_TYPE=Release)
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND archive_configure_command
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND archive_configure_command
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND archive_configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND archive_configure_command "${archive_source}")
leo2_run(${archive_configure_command})

set(archive_build_command "${CMAKE_COMMAND}" --build "${archive_build}")
set(archive_output_directory "${archive_build}")
if(LEO2_GENERATOR MATCHES "Visual Studio|Xcode|Multi-Config")
    list(APPEND archive_build_command --config Release)
    set(archive_output_directory "${archive_build}/Release")
endif()
leo2_run(${archive_build_command})
leo2_check_probe(
    "${archive_output_directory}/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}"
    "${archive_build}"
    "unknown" "unknown" 1)

# When Leopard is consumed with add_subdirectory(), the attested source is the
# inner Leopard checkout, not the enclosing superproject.
set(super_source "${LEO2_BINARY_DIR}/superproject-source")
set(inner_source "${super_source}/inner-leopard")
set(super_build "${LEO2_BINARY_DIR}/superproject-build")
file(MAKE_DIRECTORY "${inner_source}")
file(MAKE_DIRECTORY "${super_build}")
file(WRITE "${super_source}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.7)\n"
"project(leopard2_attestation_superproject CXX)\n"
"add_subdirectory(inner-leopard)\n")
file(WRITE "${super_source}/.gitignore" "/inner-leopard/\n")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}" init)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}"
    config user.name "Leopard2 superproject regression")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}"
    config user.email "leopard2-superproject@example.invalid")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}"
    config commit.gpgsign false)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}"
    add CMakeLists.txt .gitignore)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${super_source}"
    commit -m superproject-root)

file(WRITE "${inner_source}/CMakeLists.txt"
"include(\"${attestation_module}\")\n"
"add_executable(bench_leopard2 probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
file(COPY "${fixture_source}/probe.cpp" DESTINATION "${inner_source}")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}" init)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    config user.name "Leopard2 inner-source regression")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    config user.email "leopard2-inner@example.invalid")
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    config commit.gpgsign false)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    add CMakeLists.txt probe.cpp)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    commit -m inner-leopard)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    rev-parse --verify HEAD)
string(STRIP "${LEO2_LAST_STDOUT}" inner_commit)
leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${inner_source}"
    rev-parse --verify HEAD^{tree})
string(STRIP "${LEO2_LAST_STDOUT}" inner_tree)

set(super_configure_command
    "${CMAKE_COMMAND}" -E chdir "${super_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}" -DCMAKE_BUILD_TYPE=Release)
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND super_configure_command
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND super_configure_command
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND super_configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND super_configure_command "${super_source}")
leo2_run(${super_configure_command})

set(super_build_command "${CMAKE_COMMAND}" --build "${super_build}")
set(super_output_directory "${super_build}/inner-leopard")
if(LEO2_GENERATOR MATCHES "Visual Studio|Xcode|Multi-Config")
    list(APPEND super_build_command --config Release)
    set(super_output_directory "${super_build}/inner-leopard/Release")
endif()
leo2_run(${super_build_command})
leo2_check_probe(
    "${super_output_directory}/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}"
    "${super_build}"
    "${inner_commit}" "${inner_tree}" 0)

# The line-oriented v2 sidecar cannot represent embedded line delimiters.
# Reject both controls before writing rather than allowing an ambiguous digest
# material stream.
foreach(control_escape n r)
    set(invalid_source
        "${LEO2_BINARY_DIR}/invalid-${control_escape}-source")
    set(invalid_build
        "${LEO2_BINARY_DIR}/invalid-${control_escape}-build")
    file(MAKE_DIRECTORY "${invalid_source}" "${invalid_build}")
    file(WRITE "${invalid_source}/CMakeLists.txt"
        "cmake_minimum_required(VERSION 3.7)\n"
        "project(leopard2_invalid_attestation_fixture CXX)\n"
        "include(\"${attestation_module}\")\n"
        "set(CMAKE_CXX_FLAGS \"prefix\\${control_escape}suffix\")\n"
        "add_executable(bench_leopard2 probe.cpp)\n"
        "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
    file(WRITE "${invalid_source}/probe.cpp" "int main() { return 0; }\n")
    set(invalid_configure_command
        "${CMAKE_COMMAND}" -E chdir "${invalid_build}" "${CMAKE_COMMAND}"
        -G "${LEO2_GENERATOR}")
    if(DEFINED LEO2_GENERATOR_PLATFORM AND
       NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
        list(APPEND invalid_configure_command
            -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(DEFINED LEO2_GENERATOR_TOOLSET AND
       NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
        list(APPEND invalid_configure_command
            -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    if(DEFINED LEO2_CXX_COMPILER AND
       NOT LEO2_CXX_COMPILER STREQUAL "")
        list(APPEND invalid_configure_command
            "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
    endif()
    list(APPEND invalid_configure_command "${invalid_source}")
    execute_process(
        COMMAND ${invalid_configure_command}
        RESULT_VARIABLE invalid_result
        OUTPUT_VARIABLE invalid_output
        ERROR_VARIABLE invalid_error)
    if(invalid_result EQUAL 0)
        message(FATAL_ERROR
            "Attestation accepted an embedded \\${control_escape} value")
    endif()
    set(invalid_log "${invalid_output}\n${invalid_error}")
    if(NOT invalid_log MATCHES "Cannot attest CMAKE_CXX_FLAGS:" OR
       NOT invalid_log MATCHES "not contain CR or LF characters")
        message(FATAL_ERROR
            "Embedded \\${control_escape} failed for an unexpected reason:\n"
            "${invalid_log}")
    endif()
endforeach()

# Reject a syntactically well-framed but unsupported mode before it can be
# hashed into a seemingly valid sidecar.
set(invalid_mode_source "${LEO2_BINARY_DIR}/invalid-mode-source")
set(invalid_mode_build "${LEO2_BINARY_DIR}/invalid-mode-build")
file(MAKE_DIRECTORY "${invalid_mode_source}" "${invalid_mode_build}")
file(WRITE "${invalid_mode_source}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.7)\n"
    "project(leopard2_invalid_mode_attestation_fixture CXX)\n"
    "include(\"${attestation_module}\")\n"
    "set(LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE \"3\")\n"
    "add_executable(bench_leopard2 probe.cpp)\n"
    "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
file(WRITE "${invalid_mode_source}/probe.cpp"
    "int main() { return 0; }\n")
set(invalid_mode_configure_command
    "${CMAKE_COMMAND}" -E chdir "${invalid_mode_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND invalid_mode_configure_command
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND invalid_mode_configure_command
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND
   NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND invalid_mode_configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND invalid_mode_configure_command "${invalid_mode_source}")
execute_process(
    COMMAND ${invalid_mode_configure_command}
    RESULT_VARIABLE invalid_mode_result
    OUTPUT_VARIABLE invalid_mode_output
    ERROR_VARIABLE invalid_mode_error)
if(invalid_mode_result EQUAL 0)
    message(FATAL_ERROR
        "Attestation accepted unsupported small-direct mode 3")
endif()
set(invalid_mode_log "${invalid_mode_output}\n${invalid_mode_error}")
if(NOT invalid_mode_log MATCHES
       "Cannot attest LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE:" OR
   NOT invalid_mode_log MATCHES "expected exactly 0, 1,")
    message(FATAL_ERROR
        "Unsupported small-direct mode failed for an unexpected reason:\n"
        "${invalid_mode_log}")
endif()

# CMake accepts several truthy BOOL spellings, but build attestation accepts
# only the canonical wire spellings ON and OFF.
set(invalid_dual_source "${LEO2_BINARY_DIR}/invalid-dual-source")
set(invalid_dual_build "${LEO2_BINARY_DIR}/invalid-dual-build")
file(MAKE_DIRECTORY "${invalid_dual_source}" "${invalid_dual_build}")
file(WRITE "${invalid_dual_source}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.7)\n"
    "project(leopard2_invalid_dual_attestation_fixture CXX)\n"
    "include(\"${attestation_module}\")\n"
    "set(LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT \"TRUE\")\n"
    "add_executable(bench_leopard2 probe.cpp)\n"
    "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
file(WRITE "${invalid_dual_source}/probe.cpp"
    "int main() { return 0; }\n")
set(invalid_dual_configure_command
    "${CMAKE_COMMAND}" -E chdir "${invalid_dual_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND invalid_dual_configure_command
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND invalid_dual_configure_command
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND
   NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND invalid_dual_configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND invalid_dual_configure_command "${invalid_dual_source}")
execute_process(
    COMMAND ${invalid_dual_configure_command}
    RESULT_VARIABLE invalid_dual_result
    OUTPUT_VARIABLE invalid_dual_output
    ERROR_VARIABLE invalid_dual_error)
if(invalid_dual_result EQUAL 0)
    message(FATAL_ERROR
        "Attestation accepted noncanonical dual-direct value TRUE")
endif()
set(invalid_dual_log "${invalid_dual_output}\n${invalid_dual_error}")
if(NOT invalid_dual_log MATCHES
       "Cannot attest LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT:" OR
   NOT invalid_dual_log MATCHES "expected exactly ON or OFF")
    message(FATAL_ERROR
        "Noncanonical dual-direct value failed for an unexpected reason:\n"
        "${invalid_dual_log}")
endif()

# One build tree cannot silently reuse the first target's source attestation
# for a target that requests a different source root.
set(conflicting_source "${LEO2_BINARY_DIR}/conflicting-source-root-source")
set(conflicting_build "${LEO2_BINARY_DIR}/conflicting-source-root-build")
file(MAKE_DIRECTORY
    "${conflicting_source}/first" "${conflicting_source}/second"
    "${conflicting_build}")
file(WRITE "${conflicting_source}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.7)\n"
    "project(leopard2_conflicting_attestation_fixture CXX)\n"
    "include(\"${attestation_module}\")\n"
    "file(WRITE \"\${CMAKE_CURRENT_BINARY_DIR}/probe.cpp\" \"int main() { return 0; }\\n\")\n"
    "add_executable(first_benchmark \"\${CMAKE_CURRENT_BINARY_DIR}/probe.cpp\")\n"
    "leopard2_enable_benchmark_source_attestation(first_benchmark \"\${CMAKE_CURRENT_SOURCE_DIR}/first\")\n"
    "add_executable(second_benchmark \"\${CMAKE_CURRENT_BINARY_DIR}/probe.cpp\")\n"
    "leopard2_enable_benchmark_source_attestation(second_benchmark \"\${CMAKE_CURRENT_SOURCE_DIR}/second\")\n")
set(conflicting_configure_command
    "${CMAKE_COMMAND}" -E chdir "${conflicting_build}" "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(DEFINED LEO2_GENERATOR_PLATFORM AND
   NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND conflicting_configure_command
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND
   NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND conflicting_configure_command
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
if(DEFINED LEO2_CXX_COMPILER AND
   NOT LEO2_CXX_COMPILER STREQUAL "")
    list(APPEND conflicting_configure_command
        "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
endif()
list(APPEND conflicting_configure_command "${conflicting_source}")
execute_process(
    COMMAND ${conflicting_configure_command}
    RESULT_VARIABLE conflicting_result
    OUTPUT_VARIABLE conflicting_output
    ERROR_VARIABLE conflicting_error)
if(conflicting_result EQUAL 0)
    message(FATAL_ERROR
        "Attestation accepted inconsistent target source roots")
endif()
set(conflicting_log "${conflicting_output}\n${conflicting_error}")
if(NOT conflicting_log MATCHES
       "requested inconsistent attested source roots")
    message(FATAL_ERROR
        "Inconsistent source roots failed for an unexpected reason:\n"
        "${conflicting_log}")
endif()

if(UNIX)
    # Reject an ancestor symlink before MAKE_DIRECTORY can follow it and
    # create the attestation directory under an unrelated victim.
    set(ancestor_symlink_source
        "${LEO2_BINARY_DIR}/ancestor-symlink-source")
    set(ancestor_symlink_build
        "${LEO2_BINARY_DIR}/ancestor-symlink-build")
    set(ancestor_symlink_victim
        "${LEO2_BINARY_DIR}/ancestor-symlink-victim")
    file(MAKE_DIRECTORY
        "${ancestor_symlink_source}" "${ancestor_symlink_build}"
        "${ancestor_symlink_victim}")
    file(WRITE "${ancestor_symlink_source}/CMakeLists.txt"
        "cmake_minimum_required(VERSION 3.7)\n"
        "project(leopard2_ancestor_symlink_fixture CXX)\n"
        "include(\"${attestation_module}\")\n"
        "add_executable(bench_leopard2 probe.cpp)\n"
        "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
    file(WRITE "${ancestor_symlink_source}/probe.cpp"
        "int main() { return 0; }\n")
    leo2_run("${CMAKE_COMMAND}" -E create_symlink
        "${ancestor_symlink_victim}"
        "${ancestor_symlink_build}/generated")
    set(ancestor_symlink_configure_command
        "${CMAKE_COMMAND}" -E chdir "${ancestor_symlink_build}"
        "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
    if(DEFINED LEO2_GENERATOR_PLATFORM AND
       NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
        list(APPEND ancestor_symlink_configure_command
            -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(DEFINED LEO2_GENERATOR_TOOLSET AND
       NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
        list(APPEND ancestor_symlink_configure_command
            -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    if(DEFINED LEO2_CXX_COMPILER AND
       NOT LEO2_CXX_COMPILER STREQUAL "")
        list(APPEND ancestor_symlink_configure_command
            "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
    endif()
    list(APPEND ancestor_symlink_configure_command
        "${ancestor_symlink_source}")
    execute_process(
        COMMAND ${ancestor_symlink_configure_command}
        RESULT_VARIABLE ancestor_symlink_result
        OUTPUT_VARIABLE ancestor_symlink_output
        ERROR_VARIABLE ancestor_symlink_error)
    if(ancestor_symlink_result EQUAL 0)
        message(FATAL_ERROR
            "Attestation accepted a generated-directory symbolic link")
    endif()
    set(ancestor_symlink_log
        "${ancestor_symlink_output}\n${ancestor_symlink_error}")
    if(NOT ancestor_symlink_log MATCHES
           "attestation directory must not traverse a symbolic link" OR
       EXISTS
           "${ancestor_symlink_victim}/leopard2-benchmark-attestation")
        message(FATAL_ERROR
            "Ancestor symlink rejection had an external side effect or "
            "failed unexpectedly:\n${ancestor_symlink_log}")
    endif()

    # A pre-existing sidecar symlink must neither be followed nor replaced.
    set(symlink_source "${LEO2_BINARY_DIR}/sidecar-symlink-source")
    set(symlink_build "${LEO2_BINARY_DIR}/sidecar-symlink-build")
    set(symlink_generated
        "${symlink_build}/generated/leopard2-benchmark-attestation")
    set(symlink_victim "${LEO2_BINARY_DIR}/sidecar-symlink-victim.txt")
    file(MAKE_DIRECTORY "${symlink_source}" "${symlink_generated}")
    file(WRITE "${symlink_victim}" "unchanged\n")
    file(WRITE "${symlink_source}/CMakeLists.txt"
        "cmake_minimum_required(VERSION 3.7)\n"
        "project(leopard2_sidecar_symlink_fixture CXX)\n"
        "include(\"${attestation_module}\")\n"
        "add_executable(bench_leopard2 probe.cpp)\n"
        "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
    file(WRITE "${symlink_source}/probe.cpp" "int main() { return 0; }\n")
    leo2_run("${CMAKE_COMMAND}" -E create_symlink
        "${symlink_victim}"
        "${symlink_generated}/leopard2_benchmark_build_configuration.txt")
    set(symlink_configure_command
        "${CMAKE_COMMAND}" -E chdir "${symlink_build}" "${CMAKE_COMMAND}"
        -G "${LEO2_GENERATOR}")
    if(DEFINED LEO2_GENERATOR_PLATFORM AND
       NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
        list(APPEND symlink_configure_command
            -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(DEFINED LEO2_GENERATOR_TOOLSET AND
       NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
        list(APPEND symlink_configure_command
            -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    if(DEFINED LEO2_CXX_COMPILER AND
       NOT LEO2_CXX_COMPILER STREQUAL "")
        list(APPEND symlink_configure_command
            "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
    endif()
    list(APPEND symlink_configure_command "${symlink_source}")
    execute_process(
        COMMAND ${symlink_configure_command}
        RESULT_VARIABLE symlink_result
        OUTPUT_VARIABLE symlink_output
        ERROR_VARIABLE symlink_error)
    if(symlink_result EQUAL 0)
        message(FATAL_ERROR "Attestation accepted a sidecar symbolic link")
    endif()
    set(symlink_log "${symlink_output}\n${symlink_error}")
    file(READ "${symlink_victim}" symlink_victim_content)
    if(NOT symlink_log MATCHES
           "effective-configuration sidecar must not be a symbolic" OR
       NOT symlink_victim_content STREQUAL "unchanged\n")
        message(FATAL_ERROR
            "Sidecar symlink rejection was unsafe or failed unexpectedly:\n"
            "${symlink_log}")
    endif()
endif()

# Exercise the generator classification itself with real Ninja single- and
# multi-config builds.  CMAKE_CONFIGURATION_TYPES may be populated manually in
# a single-config project, so it is not a reliable classifier.
leo2_run("${CMAKE_COMMAND}" --help)
string(FIND "${LEO2_LAST_STDOUT}" "Ninja Multi-Config"
    ninja_multi_config_index)
find_program(LEO2_TEST_NINJA_EXECUTABLE NAMES ninja ninja-build)
if(LEO2_TEST_NINJA_EXECUTABLE)
    set(ninja_source "${LEO2_BINARY_DIR}/ninja-source")
    file(MAKE_DIRECTORY "${ninja_source}")
    file(WRITE "${ninja_source}/CMakeLists.txt"
        "cmake_minimum_required(VERSION 3.7)\n"
        "project(leopard2_ninja_attestation_fixture CXX)\n"
        "set(CMAKE_CONFIGURATION_TYPES \"Debug;Release\" CACHE STRING \"\" FORCE)\n"
        "include(\"${attestation_module}\")\n"
        "add_executable(bench_leopard2 probe.cpp)\n"
        "leopard2_enable_benchmark_source_attestation(bench_leopard2)\n")
    file(COPY "${fixture_source}/probe.cpp" DESTINATION "${ninja_source}")
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}" init)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        config user.name "Leopard2 Ninja regression")
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        config user.email "leopard2-ninja@example.invalid")
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        config commit.gpgsign false)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        add CMakeLists.txt probe.cpp)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        commit -m ninja-configurations)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        rev-parse --verify HEAD)
    string(STRIP "${LEO2_LAST_STDOUT}" ninja_commit)
    leo2_run("${LEO2_TEST_GIT_EXECUTABLE}" -C "${ninja_source}"
        rev-parse --verify HEAD^{tree})
    string(STRIP "${LEO2_LAST_STDOUT}" ninja_tree)

    set(single_build "${LEO2_BINARY_DIR}/ninja-single-build")
    file(MAKE_DIRECTORY "${single_build}")
    set(single_configure_command
        "${CMAKE_COMMAND}" -E chdir "${single_build}" "${CMAKE_COMMAND}"
        -G Ninja -DCMAKE_BUILD_TYPE=Release
        "-DCMAKE_MAKE_PROGRAM=${LEO2_TEST_NINJA_EXECUTABLE}")
    if(DEFINED LEO2_CXX_COMPILER AND
       NOT LEO2_CXX_COMPILER STREQUAL "")
        list(APPEND single_configure_command
            "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
    endif()
    list(APPEND single_configure_command "${ninja_source}")
    leo2_run(${single_configure_command})
    leo2_read_effective_configuration(
        "${single_build}" single_configuration single_build_type
        single_configuration_material)
    if(NOT single_build_type STREQUAL "Release" OR
       NOT single_configuration_material MATCHES
           "CMAKE_GENERATOR=Ninja\nCMAKE_CONFIGURATION_TYPES=Debug;Release\n" OR
       NOT single_configuration_material MATCHES
           "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$")
        message(FATAL_ERROR
            "Ninja single-config sidecar was misclassified:\n"
            "${single_configuration_material}")
    endif()
    leo2_run("${CMAKE_COMMAND}" --build "${single_build}")
    leo2_check_probe(
        "${single_build}/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}"
        "${single_build}" "${ninja_commit}" "${ninja_tree}" 0)

    if(NOT ninja_multi_config_index EQUAL -1)
        set(multi_build "${LEO2_BINARY_DIR}/ninja-multi-build")
        file(MAKE_DIRECTORY "${multi_build}")
        set(multi_configure_command
            "${CMAKE_COMMAND}" -E chdir "${multi_build}" "${CMAKE_COMMAND}"
            -G "Ninja Multi-Config"
            "-DCMAKE_MAKE_PROGRAM=${LEO2_TEST_NINJA_EXECUTABLE}")
        if(DEFINED LEO2_CXX_COMPILER AND
           NOT LEO2_CXX_COMPILER STREQUAL "")
            list(APPEND multi_configure_command
                "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}")
        endif()
        list(APPEND multi_configure_command "${ninja_source}")
        leo2_run(${multi_configure_command})
        leo2_read_effective_configuration(
            "${multi_build}" multi_configuration multi_build_type
            multi_configuration_material)
        string(REGEX MATCH
            "CMAKE_CONFIGURATION_TYPES=([^\n]*)\n"
            ignored_configuration_types "${multi_configuration_material}")
        set(multi_configuration_types "${CMAKE_MATCH_1}")
        list(FIND multi_configuration_types Debug multi_debug_index)
        list(FIND multi_configuration_types Release multi_release_index)
        if(NOT multi_build_type STREQUAL "" OR
           NOT multi_configuration_material MATCHES
               "CMAKE_GENERATOR=Ninja Multi-Config\n" OR
           NOT multi_configuration_material MATCHES
               "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=\nLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR=\nLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=\nLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=\nLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=\nLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=\nLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=\nLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=\nLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=\nLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0\nLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=ON\nLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=ON\nLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=OFF\n$" OR
           multi_debug_index EQUAL -1 OR multi_release_index EQUAL -1)
            message(FATAL_ERROR
                "Multi-config sidecar did not declare an empty "
                "CMAKE_BUILD_TYPE and Debug/Release configuration set:\n"
                "${multi_configuration_material}")
        endif()

        leo2_run("${CMAKE_COMMAND}" --build "${multi_build}" --config Debug)
        leo2_run("${CMAKE_COMMAND}" --build "${multi_build}" --config Release)
        leo2_check_probe(
            "${multi_build}/Debug/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}"
            "${multi_build}" "${ninja_commit}" "${ninja_tree}" 0 Debug)
        leo2_check_probe(
            "${multi_build}/Release/bench_leopard2${LEO2_EXECUTABLE_SUFFIX}"
            "${multi_build}" "${ninja_commit}" "${ninja_tree}" 0 Release)
    else()
        message(STATUS
            "Skipping Ninja Multi-Config attestation regression: "
            "generator unavailable")
    endif()
else()
    message(STATUS
        "Skipping Ninja attestation regressions: executable unavailable")
endif()
