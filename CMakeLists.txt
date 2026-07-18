cmake_minimum_required(VERSION 3.7)
project(leopard)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 11)

set(LIB_SOURCE_FILES
        leopard.cpp
        leopard.h
        LeopardCommon.cpp
        LeopardCommon.h
        LeopardFF16.cpp
        LeopardFF16.h
        LeopardFF8.cpp
        LeopardFF8.h
        LeopardFF8Xor.cpp
        LeopardFF8Xor.h
        leopard_ff8xor.cpp
        leopard_ff8xor.h)

set(BENCH_SOURCE_FILES
        tests/benchmark.cpp)

set(EXPERIMENT_SOURCE_FILES
        tests/experiments.cpp)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Keep the experiment honest about its AVX2 register budget even when the host
# compiler enables AVX-512 through -march=native.  GCC's reassociation pass also
# lengthens live ranges across the generated CNOT sequence and spills the 16
# plane registers, so preserve the generator's verified operation ordering.
set(FF8XOR_COMPILE_FLAGS "")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
    check_cxx_compiler_flag("-mno-avx512f" COMPILER_SUPPORTS_NO_AVX512F)
    if(COMPILER_SUPPORTS_NO_AVX512F)
        set(FF8XOR_COMPILE_FLAGS "${FF8XOR_COMPILE_FLAGS} -mno-avx512f")
    endif()
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    check_cxx_compiler_flag("-fno-tree-reassoc" COMPILER_SUPPORTS_NO_TREE_REASSOC)
    if(COMPILER_SUPPORTS_NO_TREE_REASSOC)
        set(FF8XOR_COMPILE_FLAGS "${FF8XOR_COMPILE_FLAGS} -fno-tree-reassoc")
    endif()
endif()
if(FF8XOR_COMPILE_FLAGS)
    set_source_files_properties(LeopardFF8Xor.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_COMPILE_FLAGS}")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    check_cxx_compiler_flag("/O2" CXX_FLAG_O2)
    check_cxx_compiler_flag("/Oy" CXX_FLAG_Oy)
    check_cxx_compiler_flag("/Zi" CXX_FLAG_Zi)
    check_cxx_compiler_flag("/W4" CXX_FLAG_W4)

    if(CXX_FLAG_O2)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    endif()

    if(CXX_FLAG_Oy)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Oy")
    endif()

    if(CXX_FLAG_Zi)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
    endif()

    if(CXX_FLAG_W4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    check_cxx_compiler_flag("-G" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Minform=inform" CXX_FLAG_MINFORM_INFORM)

    if(CXX_FLAG_g)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif(CXX_FLAG_g)

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3")
    endif(CXX_FLAG_O3)

    if(CXX_FLAG_MINFORM_INFORM)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
    endif(CXX_FLAG_MINFORM_INFORM)
else()
    check_cxx_compiler_flag("-G" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Wall" CXX_FLAG_Wall)
    check_cxx_compiler_flag("-Wextra" CXX_FLAG_Wextra)

    if(CXX_FLAG_g)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif(CXX_FLAG_g)

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3")
    endif(CXX_FLAG_O3)

    if(CXX_FLAG_Wall)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    endif(CXX_FLAG_Wall)

    if(CXX_FLAG_Wextra)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    endif(CXX_FLAG_Wextra)
endif()

find_package(OpenMP)
cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" ON "OPENMP_FOUND" OFF)
if(ENABLE_OPENMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP support requested, but could not be found")
    endif(NOT OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(ENABLE_OPENMP)

add_library(libleopard STATIC ${LIB_SOURCE_FILES})

add_executable(bench_leopard ${BENCH_SOURCE_FILES})
target_link_libraries(bench_leopard libleopard)

add_executable(experiment_leopard ${EXPERIMENT_SOURCE_FILES})

enable_testing()

add_executable(test_leopard_ff8xor tests/test_ff8xor.cpp)
target_link_libraries(test_leopard_ff8xor libleopard)
add_test(NAME leopard_ff8xor COMMAND test_leopard_ff8xor)

add_executable(test_leopard_ff8xor_statistics tests/test_ff8xor_statistics.cpp)
target_link_libraries(test_leopard_ff8xor_statistics libleopard)
add_test(NAME leopard_ff8xor_statistics COMMAND test_leopard_ff8xor_statistics)

add_executable(bench_leopard_ff8xor tests/benchmark_ff8xor.cpp)
target_link_libraries(bench_leopard_ff8xor libleopard)

# The normal build consumes the checked-in generated source and does not need
# Python.  These optional targets only regenerate or verify it for developers.
# FindPython3 was added after this project's CMake 3.7 minimum.
if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3 QUIET)
    set(FF8XOR_PYTHON_EXECUTABLE "${PYTHON_EXECUTABLE}")
else()
    find_package(Python3 COMPONENTS Interpreter QUIET)
    set(FF8XOR_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
endif()

if(FF8XOR_PYTHON_EXECUTABLE)
    add_custom_target(generate_ff8_xor_circuits
        COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8_xor_circuits.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

    add_custom_target(check_ff8_xor_circuits
        COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8_xor_circuits.py"
                --check
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()
