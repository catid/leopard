if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_CTEST_COMMAND OR
   NOT DEFINED LEO2_STATIC_LIBRARY_PREFIX OR
   NOT DEFINED LEO2_STATIC_LIBRARY_SUFFIX)
    message(FATAL_ERROR
        "CUDA optional test requires source, binary, generator, CTest, and "
        "static-library naming inputs")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(install_prefix "${LEO2_BINARY_DIR}/stage")

# Use the traditional source-directory configure form rather than -S/-B so
# this regression test remains runnable at the project's CMake 3.7 floor.
# LEO2_ENABLE_CUDA is intentionally omitted: the test verifies the declared
# default, rather than forcing the option OFF itself.
set(configure_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND configure_command
    -DLEO2_BUILD_TESTS=OFF
    -DLEO2_BUILD_BENCHMARKS=OFF
    -DCMAKE_BUILD_TYPE=Release
    "-DCMAKE_INSTALL_PREFIX=${install_prefix}"
    "${LEO2_SOURCE_DIR}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        "CUDACXX=${LEO2_BINARY_DIR}/definitely-not-a-cuda-compiler"
        ${configure_command}
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)

if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "Default non-CUDA configure failed (${configure_result})\n"
        "stdout:\n${configure_stdout}\n"
        "stderr:\n${configure_stderr}")
endif()

set(cpu_cache_file "${LEO2_BINARY_DIR}/CMakeCache.txt")
file(READ "${cpu_cache_file}" cache_contents)
string(FIND "${cache_contents}" "LEO2_ENABLE_CUDA:BOOL=OFF" option_position)
if(option_position EQUAL -1)
    message(FATAL_ERROR "Nested configure did not retain LEO2_ENABLE_CUDA=OFF")
endif()
file(STRINGS "${cpu_cache_file}" cuda_cache_entries REGEX "^CMAKE_CUDA")
if(cuda_cache_entries)
    message(FATAL_ERROR
        "CUDA-disabled configure unexpectedly created CUDA cache entries: "
        "${cuda_cache_entries}")
endif()

# Also cover the ordinary environment where CUDACXX is wholly absent.  The
# invalid-value run above is the stronger probe-resistance case; this second
# configure protects against accidentally making the variable mandatory.
set(cuda_unset_dir "${LEO2_BINARY_DIR}-without-cudacxx")
file(REMOVE_RECURSE "${cuda_unset_dir}")
file(MAKE_DIRECTORY "${cuda_unset_dir}")
unset(ENV{CUDACXX})
execute_process(
    COMMAND ${configure_command}
    WORKING_DIRECTORY "${cuda_unset_dir}"
    RESULT_VARIABLE cuda_unset_result
    OUTPUT_VARIABLE cuda_unset_stdout
    ERROR_VARIABLE cuda_unset_stderr)
if(NOT cuda_unset_result EQUAL 0)
    message(FATAL_ERROR
        "Default configure without CUDACXX failed (${cuda_unset_result})\n"
        "stdout:\n${cuda_unset_stdout}\n"
        "stderr:\n${cuda_unset_stderr}")
endif()
set(cuda_unset_cache_file "${cuda_unset_dir}/CMakeCache.txt")
file(READ "${cuda_unset_cache_file}" cuda_unset_cache)
string(FIND "${cuda_unset_cache}" "LEO2_ENABLE_CUDA:BOOL=OFF"
    cuda_unset_option_position)
file(STRINGS "${cuda_unset_cache_file}" cuda_unset_cache_entries
    REGEX "^CMAKE_CUDA")
if(cuda_unset_option_position EQUAL -1 OR cuda_unset_cache_entries)
    message(FATAL_ERROR
        "Default configure without CUDACXX did not remain CUDA-free")
endif()

# Configuration alone would not catch an accidental CUDA header or link
# dependency added to the ordinary CPU library.  Build that target as part of
# the optionality contract as well.
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build . --target leopard --config Release
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE build_result
    OUTPUT_VARIABLE build_stdout
    ERROR_VARIABLE build_stderr)

if(NOT build_result EQUAL 0)
    message(FATAL_ERROR
        "Default non-CUDA library build failed (${build_result})\n"
        "stdout:\n${build_stdout}\n"
        "stderr:\n${build_stderr}")
endif()

# Inspect the fresh build tree before install so a stale staged artifact cannot
# satisfy the naming assertion.
set(canonical_library_name
    "${LEO2_STATIC_LIBRARY_PREFIX}leopard${LEO2_STATIC_LIBRARY_SUFFIX}")
set(legacy_library_name
    "${LEO2_STATIC_LIBRARY_PREFIX}libleopard${LEO2_STATIC_LIBRARY_SUFFIX}")
set(saw_built_canonical_library FALSE)
set(saw_built_legacy_library FALSE)
file(GLOB_RECURSE built_files
    RELATIVE "${LEO2_BINARY_DIR}"
    "${LEO2_BINARY_DIR}/*")
foreach(built_file ${built_files})
    get_filename_component(built_name "${built_file}" NAME)
    if(built_name STREQUAL "${canonical_library_name}")
        set(saw_built_canonical_library TRUE)
    elseif(built_name STREQUAL "${legacy_library_name}")
        set(saw_built_legacy_library TRUE)
    endif()
endforeach()
if(NOT saw_built_canonical_library OR saw_built_legacy_library)
    message(FATAL_ERROR
        "Fresh CPU build naming contract failed: expected "
        "'${canonical_library_name}' and no '${legacy_library_name}'")
endif()

# A successful install with no install rules would not prove package
# optionality.  Stage the normal CPU package and inspect its actual headers and
# CMake export before consuming it from a separate project.
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build . --target install --config Release
    WORKING_DIRECTORY "${LEO2_BINARY_DIR}"
    RESULT_VARIABLE install_result
    OUTPUT_VARIABLE install_stdout
    ERROR_VARIABLE install_stderr)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR
        "Default non-CUDA install failed (${install_result})\n"
        "stdout:\n${install_stdout}\n"
        "stderr:\n${install_stderr}")
endif()

file(GLOB_RECURSE installed_files
    RELATIVE "${install_prefix}"
    "${install_prefix}/*")
if(NOT installed_files)
    message(FATAL_ERROR
        "Default non-CUDA install produced no artifacts")
endif()

set(saw_leopard_header FALSE)
set(saw_leopard2_header FALSE)
set(saw_package_config FALSE)
set(saw_targets_export FALSE)
set(saw_canonical_library FALSE)
set(saw_legacy_library FALSE)
foreach(installed_file ${installed_files})
    string(TOLOWER "${installed_file}" installed_file_lower)
    string(FIND "${installed_file_lower}" "cuda" cuda_name_position)
    if(NOT cuda_name_position EQUAL -1)
        message(FATAL_ERROR
            "CUDA-disabled install unexpectedly produced '${installed_file}'")
    endif()

    get_filename_component(installed_name "${installed_file}" NAME)
    if(installed_name STREQUAL "leopard.h")
        set(saw_leopard_header TRUE)
    elseif(installed_name STREQUAL "leopard2.h")
        set(saw_leopard2_header TRUE)
    elseif(installed_name STREQUAL "leopardConfig.cmake")
        set(saw_package_config TRUE)
    elseif(installed_name STREQUAL "leopardTargets.cmake")
        set(saw_targets_export TRUE)
    elseif(installed_name STREQUAL "${canonical_library_name}")
        set(saw_canonical_library TRUE)
    elseif(installed_name STREQUAL "${legacy_library_name}")
        set(saw_legacy_library TRUE)
    endif()

    if(installed_file MATCHES "\\.(cmake|h)$")
        file(READ "${install_prefix}/${installed_file}" installed_text)
        string(TOLOWER "${installed_text}" installed_text_lower)
        if(installed_text_lower MATCHES "cuda|cudart|nvidia")
            message(FATAL_ERROR
                "CUDA-disabled installed metadata/header '${installed_file}' "
                "contains a CUDA dependency")
        endif()
    endif()
endforeach()
if(NOT saw_canonical_library OR saw_legacy_library)
    message(FATAL_ERROR
        "CPU archive naming contract failed: expected "
        "'${canonical_library_name}' and no '${legacy_library_name}'")
endif()
if(NOT saw_leopard_header OR NOT saw_leopard2_header OR
   NOT saw_package_config OR NOT saw_targets_export)
    message(FATAL_ERROR
        "CPU package is incomplete: leopard.h=${saw_leopard_header}, "
        "leopard2.h=${saw_leopard2_header}, "
        "config=${saw_package_config}, targets=${saw_targets_export}")
endif()

# Consume the package only after moving its complete install tree.  This turns
# the documented relocatability property into a regression test and catches
# absolute source/build/install paths embedded in either config or export.
set(relocated_prefix "${LEO2_BINARY_DIR}/relocated-stage")
file(REMOVE_RECURSE "${relocated_prefix}")
file(RENAME "${install_prefix}" "${relocated_prefix}")
if(EXISTS "${install_prefix}" OR NOT EXISTS "${relocated_prefix}")
    message(FATAL_ERROR "Could not relocate the staged CPU package")
endif()

# Validate the exported dependency surface instead of relying only on textual
# inspection.  This fresh consumer has CUDA explicitly unavailable, imports the
# installed target, checks its transitive links, compiles against both public
# headers, links the default OpenMP-enabled static archive, and runs it.
set(consumer_source_dir "${LEO2_BINARY_DIR}/consumer-source")
set(consumer_binary_dir "${LEO2_BINARY_DIR}/consumer-build")
file(MAKE_DIRECTORY "${consumer_source_dir}" "${consumer_binary_dir}")
file(WRITE "${consumer_source_dir}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.7)
project(leopard_cpu_package_consumer LANGUAGES C CXX)

find_package(leopard CONFIG REQUIRED)
find_package(leopard CONFIG REQUIRED)
if(NOT TARGET leopard::leopard)
    message(FATAL_ERROR "Installed package did not export leopard::leopard")
endif()
get_target_property(leopard_type leopard::leopard TYPE)
if(NOT leopard_type STREQUAL "STATIC_LIBRARY")
    message(FATAL_ERROR
        "Canonical target has unexpected type '${leopard_type}'")
endif()
get_target_property(leopard_links
    leopard::leopard INTERFACE_LINK_LIBRARIES)
string(TOLOWER "${leopard_links}" leopard_links_lower)
if(leopard_links_lower MATCHES "cuda|cudart|nvidia")
    message(FATAL_ERROR
        "CPU target exposes a CUDA link dependency: ${leopard_links}")
endif()

add_executable(leopard_cpu_consumer main.c)
target_link_libraries(leopard_cpu_consumer PRIVATE leopard::leopard)

# The old imported name remains a forwarding compatibility target, not a
# second archive or canonical export.
if(NOT TARGET leopard::libleopard)
    message(FATAL_ERROR
        "Installed package lost leopard::libleopard compatibility target")
endif()
get_target_property(leopard_compat_type leopard::libleopard TYPE)
if(NOT leopard_compat_type STREQUAL "INTERFACE_LIBRARY")
    message(FATAL_ERROR
        "Compatibility target has unexpected type '${leopard_compat_type}'")
endif()
get_target_property(leopard_compat_links
    leopard::libleopard INTERFACE_LINK_LIBRARIES)
if(NOT leopard_compat_links STREQUAL "leopard::leopard")
    message(FATAL_ERROR
        "Compatibility target does not forward to leopard::leopard: "
        "${leopard_compat_links}")
endif()
add_executable(leopard_legacy_name_consumer main.c)
target_link_libraries(
    leopard_legacy_name_consumer PRIVATE leopard::libleopard)
enable_testing()
add_test(NAME leopard_cpu_consumer COMMAND leopard_cpu_consumer)
add_test(NAME leopard_legacy_name_consumer COMMAND leopard_legacy_name_consumer)
]=])
file(WRITE "${consumer_source_dir}/main.c" [=[
#include <leopard.h>
#include <leopard2.h>

int main()
{
    return leo_result_string(Leopard_Success) != 0 &&
        leo2_result_string(LEO2_SUCCESS) != 0 ? 0 : 1;
}
]=])

set(consumer_configure_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND consumer_configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND consumer_configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND consumer_configure_command
    -DCMAKE_BUILD_TYPE=Release
    "-DCMAKE_PREFIX_PATH=${relocated_prefix}"
    "${consumer_source_dir}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        "CUDACXX=${consumer_binary_dir}/definitely-not-a-cuda-compiler"
        ${consumer_configure_command}
    WORKING_DIRECTORY "${consumer_binary_dir}"
    RESULT_VARIABLE consumer_configure_result
    OUTPUT_VARIABLE consumer_configure_stdout
    ERROR_VARIABLE consumer_configure_stderr)
if(NOT consumer_configure_result EQUAL 0)
    message(FATAL_ERROR
        "Installed CPU package consumer configure failed "
        "(${consumer_configure_result})\n"
        "stdout:\n${consumer_configure_stdout}\n"
        "stderr:\n${consumer_configure_stderr}")
endif()
file(STRINGS "${consumer_binary_dir}/CMakeCache.txt" consumer_cuda_cache_entries
    REGEX "^CMAKE_CUDA")
if(consumer_cuda_cache_entries)
    message(FATAL_ERROR
        "Installed CPU package consumer unexpectedly created CUDA cache entries: "
        "${consumer_cuda_cache_entries}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build . --config Release
    WORKING_DIRECTORY "${consumer_binary_dir}"
    RESULT_VARIABLE consumer_build_result
    OUTPUT_VARIABLE consumer_build_stdout
    ERROR_VARIABLE consumer_build_stderr)
if(NOT consumer_build_result EQUAL 0)
    message(FATAL_ERROR
        "Installed CPU package consumer build failed (${consumer_build_result})\n"
        "stdout:\n${consumer_build_stdout}\n"
        "stderr:\n${consumer_build_stderr}")
endif()
execute_process(
    COMMAND "${LEO2_CTEST_COMMAND}" -C Release --output-on-failure
    WORKING_DIRECTORY "${consumer_binary_dir}"
    RESULT_VARIABLE consumer_test_result
    OUTPUT_VARIABLE consumer_test_stdout
    ERROR_VARIABLE consumer_test_stderr)
if(NOT consumer_test_result EQUAL 0)
    message(FATAL_ERROR
        "Installed CPU package consumer test failed (${consumer_test_result})\n"
        "stdout:\n${consumer_test_stdout}\n"
        "stderr:\n${consumer_test_stderr}")
endif()

# A C source is supported, but importing a static C++ implementation from a
# project that has not enabled CXX cannot be linked portably.  Require the
# package to reject that configuration immediately with a useful diagnostic,
# rather than allowing an eventual unresolved-runtime failure from the C link
# driver.
set(c_only_source_dir "${LEO2_BINARY_DIR}/c-only-consumer-source")
set(c_only_binary_dir "${LEO2_BINARY_DIR}/c-only-consumer-build")
file(MAKE_DIRECTORY "${c_only_source_dir}" "${c_only_binary_dir}")
file(WRITE "${c_only_source_dir}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.7)
project(leopard_c_only_consumer LANGUAGES C)
find_package(leopard CONFIG REQUIRED)
]=])
set(c_only_configure_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND c_only_configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND c_only_configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND c_only_configure_command
    "-DCMAKE_PREFIX_PATH=${relocated_prefix}"
    "${c_only_source_dir}")
execute_process(
    COMMAND ${c_only_configure_command}
    WORKING_DIRECTORY "${c_only_binary_dir}"
    RESULT_VARIABLE c_only_configure_result
    OUTPUT_VARIABLE c_only_configure_stdout
    ERROR_VARIABLE c_only_configure_stderr)
if(c_only_configure_result EQUAL 0)
    message(FATAL_ERROR
        "C-only package consumer unexpectedly configured without CXX")
endif()
set(c_only_diagnostic "${c_only_configure_stdout}\n${c_only_configure_stderr}")
if(NOT c_only_diagnostic MATCHES "CXX enabled|LANGUAGES C CXX")
    message(FATAL_ERROR
        "C-only package rejection did not explain the CXX requirement:\n"
        "${c_only_diagnostic}")
endif()

# Opting in is the only point at which a CUDA toolchain becomes mandatory.
# Force a missing compiler even on CUDA-equipped test hosts and require a
# clean configure failure that names CUDA rather than falling through to a CPU
# build with the requested backend silently absent.
set(cuda_enabled_dir "${LEO2_BINARY_DIR}-enabled-without-toolkit")
file(REMOVE_RECURSE "${cuda_enabled_dir}")
file(MAKE_DIRECTORY "${cuda_enabled_dir}")
set(cuda_enabled_command
    "${CMAKE_COMMAND}"
    -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND cuda_enabled_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND cuda_enabled_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND cuda_enabled_command
    -DLEO2_ENABLE_CUDA=ON
    -DLEO2_BUILD_TESTS=OFF
    -DLEO2_BUILD_BENCHMARKS=OFF
    "${LEO2_SOURCE_DIR}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        "CUDACXX=${cuda_enabled_dir}/definitely-not-a-cuda-compiler"
        ${cuda_enabled_command}
    WORKING_DIRECTORY "${cuda_enabled_dir}"
    RESULT_VARIABLE cuda_enabled_result
    OUTPUT_VARIABLE cuda_enabled_stdout
    ERROR_VARIABLE cuda_enabled_stderr)

if(cuda_enabled_result EQUAL 0)
    message(FATAL_ERROR
        "CUDA-enabled configure unexpectedly succeeded without a CUDA compiler")
endif()
set(cuda_enabled_diagnostic "${cuda_enabled_stdout}\n${cuda_enabled_stderr}")
string(FIND "${cuda_enabled_diagnostic}" "CUDA" cuda_diagnostic_position)
if(cuda_diagnostic_position EQUAL -1)
    message(FATAL_ERROR
        "CUDA-enabled failure did not provide a CUDA-specific diagnostic:\n"
        "${cuda_enabled_diagnostic}")
endif()
