if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_CTEST_COMMAND OR
   NOT DEFINED LEO2_ARCHIVE_SYMBOL_CHECKS OR
   NOT DEFINED LEO2_STATIC_LIBRARY_PREFIX OR
   NOT DEFINED LEO2_STATIC_LIBRARY_SUFFIX)
    message(FATAL_ERROR "Field-option test is missing required tool inputs")
endif()
if(LEO2_ARCHIVE_SYMBOL_CHECKS AND
   (NOT DEFINED LEO2_AR_COMMAND OR NOT DEFINED LEO2_NM_COMMAND))
    message(FATAL_ERROR "Archive inspection was requested without ar/nm tools")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")

foreach(field_variant gf8 gf16)
    if(field_variant STREQUAL "gf8")
        set(enable_gf8 ON)
        set(enable_gf16 OFF)
        set(expected_mask 1)
        set(disabled_object "LeopardFF16.cpp")
        set(disabled_symbol "leopard::ff16::")
    else()
        set(enable_gf8 OFF)
        set(enable_gf16 ON)
        set(expected_mask 2)
        set(disabled_object "LeopardFF8.cpp")
        set(disabled_symbol "leopard::ff8::")
    endif()

    set(build_dir "${LEO2_BINARY_DIR}/${field_variant}-build")
    set(install_dir "${LEO2_BINARY_DIR}/${field_variant}-install")
    file(MAKE_DIRECTORY "${build_dir}")
    set(configure_command "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
    if(LEO2_GENERATOR_PLATFORM)
        list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(LEO2_GENERATOR_TOOLSET)
        list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    list(APPEND configure_command
        -DLEO2_BUILD_TESTS=ON
        -DLEO2_BUILD_BENCHMARKS=OFF
        -DLEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1=OFF
        -DLEO2_EXPERIMENTAL_GF8_AVX2_DIRECT_L1_PAIR_FUSION=OFF
        -DLEOPARD_ENABLE_GF8=${enable_gf8}
        -DLEOPARD_ENABLE_GF16=${enable_gf16}
        -DCMAKE_BUILD_TYPE=Release
        "-DCMAKE_INSTALL_PREFIX=${install_dir}"
        "${LEO2_SOURCE_DIR}")
    execute_process(
        COMMAND ${configure_command}
        WORKING_DIRECTORY "${build_dir}"
        RESULT_VARIABLE configure_result
        OUTPUT_VARIABLE configure_stdout
        ERROR_VARIABLE configure_stderr)
    if(NOT configure_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} configure failed (${configure_result})\n"
            "${configure_stdout}\n${configure_stderr}")
    endif()

    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build . --config Release
        WORKING_DIRECTORY "${build_dir}"
        RESULT_VARIABLE build_result
        OUTPUT_VARIABLE build_stdout
        ERROR_VARIABLE build_stderr)
    if(NOT build_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} build failed (${build_result})\n"
            "${build_stdout}\n${build_stderr}")
    endif()
    execute_process(
        COMMAND "${LEO2_CTEST_COMMAND}" -C Release -R
            "^leopard2_field_options$" --output-on-failure
        WORKING_DIRECTORY "${build_dir}"
        RESULT_VARIABLE test_result
        OUTPUT_VARIABLE test_stdout
        ERROR_VARIABLE test_stderr)
    if(NOT test_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} contract test failed (${test_result})\n"
            "${test_stdout}\n${test_stderr}")
    endif()

    if(LEO2_ARCHIVE_SYMBOL_CHECKS)
        set(library_name
            "${LEO2_STATIC_LIBRARY_PREFIX}leopard${LEO2_STATIC_LIBRARY_SUFFIX}")
        file(GLOB_RECURSE built_libraries "${build_dir}/${library_name}")
        list(LENGTH built_libraries built_library_count)
        if(NOT built_library_count EQUAL 1)
            message(FATAL_ERROR
                "${field_variant} expected one ${library_name}, found "
                "${built_library_count}: ${built_libraries}")
        endif()
        list(GET built_libraries 0 built_library)
        execute_process(
            COMMAND "${LEO2_AR_COMMAND}" t "${built_library}"
            RESULT_VARIABLE ar_result
            OUTPUT_VARIABLE archive_members
            ERROR_VARIABLE ar_stderr)
        if(NOT ar_result EQUAL 0)
            message(FATAL_ERROR
                "Could not inspect ${built_library}: ${ar_stderr}")
        endif()
        string(FIND "${archive_members}" "${disabled_object}" disabled_member)
        if(NOT disabled_member EQUAL -1)
            message(FATAL_ERROR
                "${field_variant} archive retained disabled object "
                "${disabled_object}:\n${archive_members}")
        endif()
        execute_process(
            COMMAND "${LEO2_NM_COMMAND}" -C --defined-only "${built_library}"
            RESULT_VARIABLE nm_result
            OUTPUT_VARIABLE archive_symbols
            ERROR_VARIABLE nm_stderr)
        if(NOT nm_result EQUAL 0)
            message(FATAL_ERROR "Could not inspect symbols: ${nm_stderr}")
        endif()
        string(FIND "${archive_symbols}" "${disabled_symbol}"
            disabled_symbol_pos)
        if(NOT disabled_symbol_pos EQUAL -1)
            message(FATAL_ERROR
                "${field_variant} archive retained disabled field symbols "
                "matching ${disabled_symbol}")
        endif()
    endif()

    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build . --target install --config Release
        WORKING_DIRECTORY "${build_dir}"
        RESULT_VARIABLE install_result
        OUTPUT_VARIABLE install_stdout
        ERROR_VARIABLE install_stderr)
    if(NOT install_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} install failed (${install_result})\n"
            "${install_stdout}\n${install_stderr}")
    endif()

    set(consumer_source "${LEO2_BINARY_DIR}/${field_variant}-consumer-source")
    set(consumer_build "${LEO2_BINARY_DIR}/${field_variant}-consumer-build")
    file(MAKE_DIRECTORY "${consumer_source}" "${consumer_build}")
    file(WRITE "${consumer_source}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.7)
project(leopard_field_consumer LANGUAGES CXX)
find_package(leopard CONFIG REQUIRED)
if(NOT DEFINED EXPECT_GF8 OR NOT DEFINED EXPECT_GF16 OR
   NOT DEFINED EXPECT_MASK)
    message(FATAL_ERROR "consumer expectations missing")
endif()
if(NOT "${leopard_GF8_ENABLED}" STREQUAL "${EXPECT_GF8}" OR
   NOT "${leopard_GF16_ENABLED}" STREQUAL "${EXPECT_GF16}")
    message(FATAL_ERROR
        "package capability mismatch: gf8=${leopard_GF8_ENABLED}, "
        "gf16=${leopard_GF16_ENABLED}")
endif()
add_executable(field_consumer main.cpp)
target_compile_definitions(field_consumer PRIVATE EXPECT_MASK=${EXPECT_MASK})
target_link_libraries(field_consumer PRIVATE leopard::leopard)
enable_testing()
add_test(NAME field_consumer COMMAND field_consumer)
]=])
    file(WRITE "${consumer_source}/main.cpp" [=[
#include <leopard2.h>
int main()
{
    leo2_context* context = 0;
    if (leo2_context_create(0, &context) != LEO2_SUCCESS)
        return 1;
    const unsigned mask = leo2_context_field_mask(context);
    leo2_context_destroy(context);
    return mask == EXPECT_MASK ? 0 : 2;
}
]=])
    set(consumer_configure "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
    if(LEO2_GENERATOR_PLATFORM)
        list(APPEND consumer_configure -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(LEO2_GENERATOR_TOOLSET)
        list(APPEND consumer_configure -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    list(APPEND consumer_configure
        -DCMAKE_BUILD_TYPE=Release
        "-DCMAKE_PREFIX_PATH=${install_dir}"
        -DEXPECT_GF8=${enable_gf8}
        -DEXPECT_GF16=${enable_gf16}
        -DEXPECT_MASK=${expected_mask}
        "${consumer_source}")
    execute_process(
        COMMAND ${consumer_configure}
        WORKING_DIRECTORY "${consumer_build}"
        RESULT_VARIABLE consumer_configure_result
        OUTPUT_VARIABLE consumer_configure_stdout
        ERROR_VARIABLE consumer_configure_stderr)
    if(NOT consumer_configure_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} consumer configure failed "
            "(${consumer_configure_result})\n${consumer_configure_stdout}\n"
            "${consumer_configure_stderr}")
    endif()
    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build . --config Release
        WORKING_DIRECTORY "${consumer_build}"
        RESULT_VARIABLE consumer_build_result
        OUTPUT_VARIABLE consumer_build_stdout
        ERROR_VARIABLE consumer_build_stderr)
    if(NOT consumer_build_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} consumer build failed "
            "(${consumer_build_result})\n${consumer_build_stdout}\n"
            "${consumer_build_stderr}")
    endif()
    execute_process(
        COMMAND "${LEO2_CTEST_COMMAND}" -C Release --output-on-failure
        WORKING_DIRECTORY "${consumer_build}"
        RESULT_VARIABLE consumer_test_result
        OUTPUT_VARIABLE consumer_test_stdout
        ERROR_VARIABLE consumer_test_stderr)
    if(NOT consumer_test_result EQUAL 0)
        message(FATAL_ERROR
            "${field_variant} consumer failed (${consumer_test_result})\n"
            "${consumer_test_stdout}\n${consumer_test_stderr}")
    endif()
endforeach()

# A fieldless archive has no useful mathematical codec and is rejected at
# configure time rather than producing a misleading linkable target.
set(fieldless_dir "${LEO2_BINARY_DIR}/fieldless-build")
file(MAKE_DIRECTORY "${fieldless_dir}")
set(fieldless_configure "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND fieldless_configure -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND fieldless_configure -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND fieldless_configure
    -DLEOPARD_ENABLE_GF8=OFF -DLEOPARD_ENABLE_GF16=OFF
    -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF
    "${LEO2_SOURCE_DIR}")
execute_process(
    COMMAND ${fieldless_configure}
    WORKING_DIRECTORY "${fieldless_dir}"
    RESULT_VARIABLE fieldless_result
    OUTPUT_VARIABLE fieldless_stdout
    ERROR_VARIABLE fieldless_stderr)
if(fieldless_result EQUAL 0)
    message(FATAL_ERROR "Fieldless configuration unexpectedly succeeded")
endif()
set(fieldless_diagnostic "${fieldless_stdout}\n${fieldless_stderr}")
if(NOT fieldless_diagnostic MATCHES "At least one.*GF8.*GF16")
    message(FATAL_ERROR
        "Fieldless rejection did not explain the requirement:\n"
        "${fieldless_diagnostic}")
endif()

# A GF8-specific diagnostic must not silently become an inert or misleading
# option in a field-disabled build.  Reject the contradiction at configure
# time, before its focused API test could request an unavailable codec.
set(experimental_without_gf8_dir
    "${LEO2_BINARY_DIR}/experimental-without-gf8-build")
file(MAKE_DIRECTORY "${experimental_without_gf8_dir}")
set(experimental_without_gf8_configure
    "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND experimental_without_gf8_configure
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND experimental_without_gf8_configure
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND experimental_without_gf8_configure
    -DLEOPARD_ENABLE_GF8=OFF -DLEOPARD_ENABLE_GF16=ON
    -DLEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1=ON
    -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF
    "${LEO2_SOURCE_DIR}")
execute_process(
    COMMAND ${experimental_without_gf8_configure}
    WORKING_DIRECTORY "${experimental_without_gf8_dir}"
    RESULT_VARIABLE experimental_without_gf8_result
    OUTPUT_VARIABLE experimental_without_gf8_stdout
    ERROR_VARIABLE experimental_without_gf8_stderr)
if(experimental_without_gf8_result EQUAL 0)
    message(FATAL_ERROR
        "GF8 direct-L1 experiment unexpectedly accepted GF8=OFF")
endif()
set(experimental_without_gf8_diagnostic
    "${experimental_without_gf8_stdout}\n${experimental_without_gf8_stderr}")
if(NOT experimental_without_gf8_diagnostic MATCHES
       "GF8_AVX2_GENERAL_DIRECT_L1 requires.*GF8=ON")
    message(FATAL_ERROR
        "GF8 direct-L1 rejection did not explain the requirement:\n"
        "${experimental_without_gf8_diagnostic}")
endif()

# Pair fusion is an implementation variant of the general direct-L1
# experiment, never an independently selectable production behavior.
set(pair_without_general_dir
    "${LEO2_BINARY_DIR}/pair-without-general-build")
file(MAKE_DIRECTORY "${pair_without_general_dir}")
set(pair_without_general_configure
    "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
if(LEO2_GENERATOR_PLATFORM)
    list(APPEND pair_without_general_configure
        -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(LEO2_GENERATOR_TOOLSET)
    list(APPEND pair_without_general_configure
        -T "${LEO2_GENERATOR_TOOLSET}")
endif()
list(APPEND pair_without_general_configure
    -DLEOPARD_ENABLE_GF8=ON -DLEOPARD_ENABLE_GF16=ON
    -DLEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1=OFF
    -DLEO2_EXPERIMENTAL_GF8_AVX2_DIRECT_L1_PAIR_FUSION=ON
    -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF
    "${LEO2_SOURCE_DIR}")
execute_process(
    COMMAND ${pair_without_general_configure}
    WORKING_DIRECTORY "${pair_without_general_dir}"
    RESULT_VARIABLE pair_without_general_result
    OUTPUT_VARIABLE pair_without_general_stdout
    ERROR_VARIABLE pair_without_general_stderr)
if(pair_without_general_result EQUAL 0)
    message(FATAL_ERROR
        "Direct-L1 pair fusion unexpectedly accepted general-L1=OFF")
endif()
set(pair_without_general_diagnostic
    "${pair_without_general_stdout}\n${pair_without_general_stderr}")
if(NOT pair_without_general_diagnostic MATCHES
       "DIRECT_L1_PAIR_FUSION requires.*GENERAL_DIRECT_L1=ON")
    message(FATAL_ERROR
        "Direct-L1 pair-fusion rejection did not explain the requirement:\n"
        "${pair_without_general_diagnostic}")
endif()
