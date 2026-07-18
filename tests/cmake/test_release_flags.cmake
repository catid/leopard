# Verify that Leopard's compiler branches never derive Release flags from the
# Debug cache entry, then exercise the active compiler branch through a nested
# OpenMP-disabled Release build.  GNU-family builds make unknown pragmas fatal
# so an unguarded OpenMP directive cannot silently re-enter the serial build.

if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_C_COMPILER OR NOT DEFINED LEO2_CXX_COMPILER OR
   NOT DEFINED LEO2_GENERATOR OR NOT DEFINED LEO2_CXX_COMPILER_ID)
    message(FATAL_ERROR "release flag isolation test arguments are incomplete")
endif()

file(READ "${LEO2_SOURCE_DIR}/CMakeLists.txt" leopard_cmake)
string(REGEX MATCHALL
    "set\\(CMAKE_CXX_FLAGS_RELEASE[ \t\r\n]+\"\\$\\{CMAKE_CXX_FLAGS_DEBUG\\}"
    contaminated_assignments "${leopard_cmake}")
if(contaminated_assignments)
    message(FATAL_ERROR
        "a compiler branch still derives Release flags from Debug flags: "
        "${contaminated_assignments}")
endif()

file(REMOVE_RECURSE "${LEO2_BINARY_DIR}")
set(configure_command
    "${CMAKE_COMMAND}"
    -S "${LEO2_SOURCE_DIR}"
    -B "${LEO2_BINARY_DIR}"
    -G "${LEO2_GENERATOR}"
    "-DCMAKE_C_COMPILER=${LEO2_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${LEO2_CXX_COMPILER}"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    -DCMAKE_CXX_FLAGS_DEBUG=-DLEO2_DEBUG_FLAG_MUST_NOT_LEAK
    -DCMAKE_CXX_FLAGS_RELEASE=-DLEO2_RELEASE_FLAG_MUST_SURVIVE
    -DLEO2_BUILD_TESTS=OFF
    -DLEO2_BUILD_BENCHMARKS=OFF
    -DLEO2_ENABLE_CUDA=OFF
    -DENABLE_OPENMP=OFF)
if(LEO2_CXX_COMPILER_ID STREQUAL "GNU" OR
   LEO2_CXX_COMPILER_ID STREQUAL "Clang" OR
   LEO2_CXX_COMPILER_ID STREQUAL "AppleClang")
    list(APPEND configure_command
        "-DCMAKE_CXX_FLAGS=-Werror=unknown-pragmas")
endif()
if(DEFINED LEO2_GENERATOR_PLATFORM AND NOT LEO2_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND configure_command -A "${LEO2_GENERATOR_PLATFORM}")
endif()
if(DEFINED LEO2_GENERATOR_TOOLSET AND NOT LEO2_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND configure_command -T "${LEO2_GENERATOR_TOOLSET}")
endif()

execute_process(
    COMMAND ${configure_command}
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "nested Release configure failed (${configure_result})\n"
        "stdout:\n${configure_stdout}\n"
        "stderr:\n${configure_stderr}")
endif()

set(compile_database "${LEO2_BINARY_DIR}/compile_commands.json")
if(NOT EXISTS "${compile_database}")
    message(FATAL_ERROR "nested configure did not generate compile_commands.json")
endif()
file(READ "${compile_database}" compile_commands)
string(FIND "${compile_commands}" "LEO2_RELEASE_FLAG_MUST_SURVIVE"
    release_marker)
string(FIND "${compile_commands}" "LEO2_DEBUG_FLAG_MUST_NOT_LEAK"
    debug_marker)
if(release_marker EQUAL -1)
    message(FATAL_ERROR "caller-supplied Release flags were discarded")
endif()
if(NOT debug_marker EQUAL -1)
    message(FATAL_ERROR "caller-supplied Debug flags leaked into Release")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${LEO2_BINARY_DIR}"
        --target leopard --parallel 2
    RESULT_VARIABLE build_result
    OUTPUT_VARIABLE build_stdout
    ERROR_VARIABLE build_stderr)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR
        "OpenMP-disabled strict Release build failed (${build_result})\n"
        "stdout:\n${build_stdout}\n"
        "stderr:\n${build_stderr}")
endif()
