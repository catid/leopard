if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_STATIC_LIBRARY_PREFIX OR
   NOT DEFINED LEO2_STATIC_LIBRARY_SUFFIX OR
   NOT DEFINED LEO2_EXECUTABLE_SUFFIX OR NOT DEFINED LEO2_ENABLE_OPENMP OR
   NOT DEFINED LEO2_ARCHIVE_SYMBOL_CHECKS)
    message(FATAL_ERROR "Test-hook isolation test is missing required inputs")
endif()
if(LEO2_ARCHIVE_SYMBOL_CHECKS AND
   (NOT DEFINED LEO2_NM_COMMAND OR NOT EXISTS "${LEO2_NM_COMMAND}"))
    message(FATAL_ERROR "Test-hook archive inspection requires a usable nm")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")
set(production_library_name
    "${LEO2_STATIC_LIBRARY_PREFIX}leopard${LEO2_STATIC_LIBRARY_SUFFIX}")
set(hook_library_name
    "${LEO2_STATIC_LIBRARY_PREFIX}leopard_test_hooks${LEO2_STATIC_LIBRARY_SUFFIX}")
set(production_fingerprints)

foreach(test_mode tests-on tests-off)
    if(test_mode STREQUAL "tests-on")
        set(build_tests ON)
    else()
        set(build_tests OFF)
    endif()
    set(build_dir "${LEO2_BINARY_DIR}/${test_mode}-build")
    set(install_dir "${LEO2_BINARY_DIR}/${test_mode}-install")
    file(MAKE_DIRECTORY "${build_dir}")
    set(configure_command "${CMAKE_COMMAND}" -G "${LEO2_GENERATOR}")
    if(LEO2_GENERATOR_PLATFORM)
        list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
    endif()
    if(LEO2_GENERATOR_TOOLSET)
        list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
    endif()
    list(APPEND configure_command
        -DLEO2_BUILD_TESTS=${build_tests}
        -DLEO2_BUILD_BENCHMARKS=OFF
        -DLEO2_BUILD_FUZZERS=OFF
        -DLEO2_ENABLE_CUDA=OFF
        -DENABLE_OPENMP=${LEO2_ENABLE_OPENMP}
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
            "${test_mode} configure failed (${configure_result})\n"
            "${configure_stdout}\n${configure_stderr}")
    endif()

    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build . --target leopard --config Release
        WORKING_DIRECTORY "${build_dir}"
        RESULT_VARIABLE build_result
        OUTPUT_VARIABLE build_stdout
        ERROR_VARIABLE build_stderr)
    if(NOT build_result EQUAL 0)
        message(FATAL_ERROR
            "${test_mode} production build failed (${build_result})\n"
            "${build_stdout}\n${build_stderr}")
    endif()
    file(GLOB_RECURSE production_libraries
        "${build_dir}/${production_library_name}")
    list(LENGTH production_libraries production_library_count)
    if(NOT production_library_count EQUAL 1)
        message(FATAL_ERROR
            "${test_mode} expected one ${production_library_name}, found "
            "${production_library_count}: ${production_libraries}")
    endif()
    list(GET production_libraries 0 production_library)

    if(test_mode STREQUAL "tests-on")
        execute_process(
            COMMAND "${CMAKE_COMMAND}" --build . --target leopard_test_hooks
                --config Release
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE hook_build_result
            OUTPUT_VARIABLE hook_build_stdout
            ERROR_VARIABLE hook_build_stderr)
        if(NOT hook_build_result EQUAL 0)
            message(FATAL_ERROR
                "Hook archive build failed (${hook_build_result})\n"
                "${hook_build_stdout}\n${hook_build_stderr}")
        endif()
        file(GLOB_RECURSE hook_libraries "${build_dir}/${hook_library_name}")
        list(LENGTH hook_libraries hook_library_count)
        if(NOT hook_library_count EQUAL 1)
            message(FATAL_ERROR
                "Expected one ${hook_library_name}, found "
                "${hook_library_count}: ${hook_libraries}")
        endif()
        list(GET hook_libraries 0 hook_library)
    endif()

    if(LEO2_ARCHIVE_SYMBOL_CHECKS)
        execute_process(
            COMMAND "${LEO2_NM_COMMAND}" -C --defined-only
                "${production_library}"
            RESULT_VARIABLE nm_result
            OUTPUT_VARIABLE production_symbols
            ERROR_VARIABLE nm_stderr)
        if(NOT nm_result EQUAL 0)
            message(FATAL_ERROR
                "Could not inspect ${production_library}: ${nm_stderr}")
        endif()
        if(production_symbols MATCHES
           "TestOnly|leo2_test_|TestSetSetupFault|TestShouldFailAllocation|TestGet(Scalar|SSSE3|AVX2)TableState")
            message(FATAL_ERROR
                "${test_mode} production archive exports test-hook symbols")
        endif()
        if(test_mode STREQUAL "tests-on")
            execute_process(
                COMMAND "${LEO2_NM_COMMAND}" -C --defined-only
                    "${hook_library}"
                RESULT_VARIABLE hook_nm_result
                OUTPUT_VARIABLE hook_symbols
                ERROR_VARIABLE hook_nm_stderr)
            if(NOT hook_nm_result EQUAL 0)
                message(FATAL_ERROR
                    "Could not inspect ${hook_library}: ${hook_nm_stderr}")
            endif()
            if(NOT hook_symbols MATCHES "TestOnly" OR
               NOT hook_symbols MATCHES "leo2_test_" OR
               NOT hook_symbols MATCHES
                   "TestOnlySetHighDecodeCopyFallback" OR
               NOT hook_symbols MATCHES
                   "TestOnlyHighDecodeCopyFallbackEnabled")
                message(FATAL_ERROR
                    "Hook archive does not contain its expected diagnostics")
            endif()
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
            "${test_mode} install failed (${install_result})\n"
            "${install_stdout}\n${install_stderr}")
    endif()
    file(GLOB_RECURSE installed_files RELATIVE "${install_dir}"
        "${install_dir}/*")
    set(installed_production_libraries)
    foreach(installed_file ${installed_files})
        get_filename_component(installed_name "${installed_file}" NAME)
        if(installed_name STREQUAL "${hook_library_name}" OR
           installed_file MATCHES "leopard_test_hooks" OR
           installed_file MATCHES
               "bench_leopard2_high_decode_copy_attribution")
            message(FATAL_ERROR
                "${test_mode} installed private hook artifact ${installed_file}")
        endif()
        if(installed_name STREQUAL "${production_library_name}")
            list(APPEND installed_production_libraries
                "${install_dir}/${installed_file}")
        endif()
    endforeach()
    list(LENGTH installed_production_libraries installed_library_count)
    if(NOT installed_library_count EQUAL 1)
        message(FATAL_ERROR
            "${test_mode} install expected one production archive, found "
            "${installed_library_count}: ${installed_production_libraries}")
    endif()

    set(consumer_build "${LEO2_BINARY_DIR}/${test_mode}-consumer")
    file(MAKE_DIRECTORY "${consumer_build}")
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
        "${LEO2_SOURCE_DIR}/tests/cmake/test_hook_consumer")
    execute_process(
        COMMAND ${consumer_configure}
        WORKING_DIRECTORY "${consumer_build}"
        RESULT_VARIABLE consumer_configure_result
        OUTPUT_VARIABLE consumer_configure_stdout
        ERROR_VARIABLE consumer_configure_stderr)
    if(NOT consumer_configure_result EQUAL 0)
        message(FATAL_ERROR
            "${test_mode} consumer configure failed "
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
            "${test_mode} consumer build failed (${consumer_build_result})\n"
            "${consumer_build_stdout}\n${consumer_build_stderr}")
    endif()
    set(consumer_name
        "leopard_test_hook_consumer${LEO2_EXECUTABLE_SUFFIX}")
    file(GLOB_RECURSE consumer_executables "${consumer_build}/${consumer_name}")
    list(LENGTH consumer_executables consumer_executable_count)
    if(NOT consumer_executable_count EQUAL 1)
        message(FATAL_ERROR
            "${test_mode} expected one ${consumer_name}, found "
            "${consumer_executable_count}: ${consumer_executables}")
    endif()
    list(GET consumer_executables 0 consumer_executable)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
            OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 "${consumer_executable}"
        RESULT_VARIABLE consumer_result
        OUTPUT_VARIABLE consumer_fingerprint
        ERROR_VARIABLE consumer_stderr)
    if(NOT consumer_result EQUAL 0)
        message(FATAL_ERROR
            "${test_mode} consumer failed (${consumer_result})\n"
            "${consumer_fingerprint}\n${consumer_stderr}")
    endif()
    string(STRIP "${consumer_fingerprint}" consumer_fingerprint)
    list(APPEND production_fingerprints "${consumer_fingerprint}")
endforeach()

list(GET production_fingerprints 0 tests_on_fingerprint)
list(GET production_fingerprints 1 tests_off_fingerprint)
if(NOT tests_on_fingerprint STREQUAL tests_off_fingerprint)
    message(FATAL_ERROR
        "Production API changed with LEO2_BUILD_TESTS:\n"
        "tests ON:\n${tests_on_fingerprint}\n"
        "tests OFF:\n${tests_off_fingerprint}")
endif()
