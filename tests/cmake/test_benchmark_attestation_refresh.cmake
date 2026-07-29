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
        "^schema=leopard2-benchmark-build-configuration/v1\nsha256=([0-9a-f]+)\n"
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
    list(LENGTH build_type_lines build_type_line_count)
    if(NOT build_type_line_count EQUAL 1)
        message(FATAL_ERROR
            "Effective build configuration omits CMAKE_BUILD_TYPE")
    endif()
    list(GET build_type_lines 0 effective_build_type)
    string(REGEX REPLACE "^CMAKE_BUILD_TYPE=" ""
        effective_build_type "${effective_build_type}")

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
           "leopard2-benchmark-build-configuration/v1")
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
       "CMAKE_CXX_FLAGS=[^\n]*LEO2_EFFECTIVE_ONE=1")
    message(FATAL_ERROR
        "Initial sidecar omitted the effective non-cache CXX flag:\n"
        "${initial_configuration_material}")
endif()

set(attestation_header
    "${fixture_build}/generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h")
if(NOT EXISTS "${attestation_header}")
    message(FATAL_ERROR "Fixture did not generate its attestation header")
endif()

# An unchanged refresh must not touch the generated header or relink either
# benchmark.
file(TIMESTAMP "${attestation_header}" stable_header_time "%s" UTC)
file(TIMESTAMP "${standard_executable}" stable_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" stable_allk_time "%s" UTC)
leo2_run("${CMAKE_COMMAND}" -E sleep 1)
leo2_build_and_check()
file(TIMESTAMP "${attestation_header}" repeat_header_time "%s" UTC)
file(TIMESTAMP "${standard_executable}" repeat_standard_time "%s" UTC)
file(TIMESTAMP "${allk_executable}" repeat_allk_time "%s" UTC)
if(NOT stable_header_time STREQUAL repeat_header_time OR
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

# The line-oriented v1 sidecar cannot represent embedded line delimiters.
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
           "CMAKE_GENERATOR=Ninja\nCMAKE_CONFIGURATION_TYPES=Debug;Release\n")
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
