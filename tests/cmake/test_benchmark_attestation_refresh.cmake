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

function(leo2_check_probe executable expected_commit expected_tree expected_dirty)
    leo2_run("${executable}")
    string(REPLACE "\r\n" "\n" probe "${LEO2_LAST_STDOUT}")
    set(expected
        "commit=${expected_commit}\ntree=${expected_tree}\ndirty=${expected_dirty}\n")
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
        "${expected_commit}" "${expected_tree}" "${expected_dirty}")
    leo2_check_probe("${allk_executable}"
        "${expected_commit}" "${expected_tree}" "${expected_dirty}")
endfunction()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(fixture_source "${LEO2_BINARY_DIR}/source")
set(fixture_build "${LEO2_BINARY_DIR}/build")
file(MAKE_DIRECTORY "${fixture_source}")

file(TO_CMAKE_PATH
    "${LEO2_SOURCE_DIR}/cmake/Leopard2BenchmarkAttestation.cmake"
    attestation_module)
file(WRITE "${fixture_source}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.7)\n"
"project(leopard2_attestation_fixture CXX)\n"
"include(\"${attestation_module}\")\n"
"add_executable(bench_leopard2 probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2)\n"
"add_executable(bench_leopard2_allk probe.cpp)\n"
"leopard2_enable_benchmark_source_attestation(bench_leopard2_allk)\n")
file(WRITE "${fixture_source}/probe.cpp"
    "#include <iostream>\n"
    "#if !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER)\n"
    "#error \"missing exact generated attestation header path\"\n"
    "#endif\n"
    "#include LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER\n"
    "int main() {\n"
"  std::cout << \"commit=\" << LEO2_BENCHMARK_SOURCE_COMMIT << \"\\n\"\n"
"            << \"tree=\" << LEO2_BENCHMARK_SOURCE_TREE << \"\\n\"\n"
"            << \"dirty=\" << LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY << \"\\n\";\n"
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
    "${CMAKE_COMMAND}" -S "${fixture_source}" -B "${fixture_build}"
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
    "${CMAKE_COMMAND}" -S "${archive_source}" -B "${archive_build}"
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
    "unknown" "unknown" 1)

# When Leopard is consumed with add_subdirectory(), the attested source is the
# inner Leopard checkout, not the enclosing superproject.
set(super_source "${LEO2_BINARY_DIR}/superproject-source")
set(inner_source "${super_source}/inner-leopard")
set(super_build "${LEO2_BINARY_DIR}/superproject-build")
file(MAKE_DIRECTORY "${inner_source}")
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
    "${CMAKE_COMMAND}" -S "${super_source}" -B "${super_build}"
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
    "${inner_commit}" "${inner_tree}" 0)
