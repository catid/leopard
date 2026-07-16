if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_GENERATOR)
    message(FATAL_ERROR "CUDA optional test requires source, binary, and generator paths")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")

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

file(READ "${LEO2_BINARY_DIR}/CMakeCache.txt" cache_contents)
string(FIND "${cache_contents}" "LEO2_ENABLE_CUDA:BOOL=OFF" option_position)
if(option_position EQUAL -1)
    message(FATAL_ERROR "Nested configure did not retain LEO2_ENABLE_CUDA=OFF")
endif()
string(FIND "${cache_contents}" "CMAKE_CUDA_COMPILER" compiler_position)
if(NOT compiler_position EQUAL -1)
    message(FATAL_ERROR "CUDA-disabled configure unexpectedly probed a CUDA compiler")
endif()

# Configuration alone would not catch an accidental CUDA header or link
# dependency added to the ordinary CPU library.  Build that target as part of
# the optionality contract as well.
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build . --target libleopard --config Release
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
