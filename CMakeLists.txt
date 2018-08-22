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
        LeopardFF8.h)

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
