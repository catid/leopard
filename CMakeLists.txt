cmake_minimum_required(VERSION 3.7)
project(leopard)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

# CMAKE_CONFIGURATION_TYPES is intentionally populated above even for some
# single-configuration generators, so use the generator identity where later
# logic must distinguish configuration models.
set(FF8XOR_MULTI_CONFIG OFF)
if(CMAKE_GENERATOR MATCHES
   "^(Visual Studio|Xcode|Ninja Multi-Config|Green Hills MULTI)")
    set(FF8XOR_MULTI_CONFIG ON)
endif()
set(FF8XOR_RELEASE_TEST_CONFIG_ARGS "")
set(FF8XOR_BUILD_CONFIG_ARGS "")
if(FF8XOR_MULTI_CONFIG)
    list(APPEND FF8XOR_RELEASE_TEST_CONFIG_ARGS CONFIGURATIONS Release)
    list(APPEND FF8XOR_BUILD_CONFIG_ARGS --config "$<CONFIG>")
endif()

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
        LeopardFF8XorAVX2.cpp
        LeopardFF8XorAVX2.h
        LeopardFF8XorAVX512.cpp
        LeopardFF8XorAVX512.h
        LeopardFF8XorTranspose.cpp
        LeopardFF8XorTranspose.h
        LeopardFF8XorTransposeAVX2.cpp
        LeopardFF8XorTransposeAVX2.h
        LeopardFF8XorTransposeAVX512.h
        LeopardFF8XorTransposeBITALG.cpp
        LeopardFF8XorTransposeVBMI.cpp
        leopard_ff8xor.cpp
        leopard_ff8xor.h)

set(BENCH_SOURCE_FILES
        tests/benchmark.cpp)

set(EXPERIMENT_SOURCE_FILES
        tests/experiments.cpp)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" FF8XOR_BUILD_TYPE_UPPER)

check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Compile the experimental dispatcher, portable uint64_t path, and 128-bit
# tail for the baseline x86-64 ISA.  This trailing source flag resets the
# inherited project-wide -march=native; AVX2 and AVX-512 are confined to the
# separately targeted translation units below.
set(FF8XOR_BASELINE_COMPILE_FLAGS "")
set(FF8XOR_AVX2_COMPILE_FLAGS "")
set(FF8XOR_BUILD_AVX2 OFF)
option(LEO_FF8XOR_ENABLE_AVX2
    "Build isolated experimental AVX2 FF8 XOR circuit and boundary-transpose kernels" ON)
if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
       CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(FF8XOR_BASELINE_COMPILE_FLAGS "-march=x86-64")
        check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2_FLAG)
        if(LEO_FF8XOR_ENABLE_AVX2 AND COMPILER_SUPPORTS_AVX2_FLAG)
            set(FF8XOR_BUILD_AVX2 ON)
            set(FF8XOR_AVX2_COMPILE_FLAGS "-march=x86-64 -mavx2")
        endif()
    elseif(MSVC)
        # MSVC x64 already has an SSE2 baseline.  Keep AVX2 in the isolated
        # source so the public dispatcher remains safe on pre-AVX2 CPUs.
        check_cxx_compiler_flag("/arch:AVX2" COMPILER_SUPPORTS_MSVC_AVX2_FLAG)
        if(LEO_FF8XOR_ENABLE_AVX2 AND COMPILER_SUPPORTS_MSVC_AVX2_FLAG)
            set(FF8XOR_BUILD_AVX2 ON)
            set(FF8XOR_AVX2_COMPILE_FLAGS "/arch:AVX2")
        endif()
    endif()
endif()

# A caller may add an extension flag such as -mssse3 globally after choosing a
# generic -march.  Explicitly negate every extension this baseline object must
# not inherit; probe each spelling so older GCC/Clang versions are not rejected
# merely because they predate a newer extension.
if(FF8XOR_BASELINE_COMPILE_FLAGS AND
   (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(FF8XOR_BASELINE_DISABLE_FLAGS
        -mno-sse3
        -mno-ssse3
        -mno-sse4.1
        -mno-sse4.2
        -mno-avx
        -mno-avx2
        -mno-avx512f
        -mno-avx512bw
        -mno-avx512cd
        -mno-avx512dq
        -mno-avx512er
        -mno-avx512ifma
        -mno-avx512pf
        -mno-avx512vl
        -mno-avx512bitalg
        -mno-avx512vbmi
        -mno-avx512vbmi2
        -mno-avx512vpopcntdq
        -mno-avx512vnni
        -mno-avx512bf16
        -mno-avx512fp16
        -mno-avx512vp2intersect
        -mno-avx5124fmaps
        -mno-avx5124vnniw
        -mno-avxvnni
        -mno-avxvnniint8
        -mno-avxvnniint16
        -mno-bmi
        -mno-bmi2
        -mno-lzcnt
        -mno-popcnt
        -mno-aes
        -mno-pclmul
        -mno-fma
        -mno-gfni
        -mno-vaes
        -mno-vpclmulqdq)
    foreach(FF8XOR_BASELINE_DISABLE_FLAG ${FF8XOR_BASELINE_DISABLE_FLAGS})
        string(MAKE_C_IDENTIFIER
            "FF8XOR_SUPPORTS_${FF8XOR_BASELINE_DISABLE_FLAG}"
            FF8XOR_BASELINE_DISABLE_VARIABLE)
        check_cxx_compiler_flag(
            "${FF8XOR_BASELINE_DISABLE_FLAG}"
            ${FF8XOR_BASELINE_DISABLE_VARIABLE})
        if(${FF8XOR_BASELINE_DISABLE_VARIABLE})
            set(FF8XOR_BASELINE_COMPILE_FLAGS
                "${FF8XOR_BASELINE_COMPILE_FLAGS} ${FF8XOR_BASELINE_DISABLE_FLAG}")
        endif()
    endforeach()
endif()

# GCC's reassociation pass lengthens live ranges across the generated CNOT
# sequence and can spill the 16 plane registers.  Preserve the generator's
# verified operation ordering in every circuit-bearing object.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    check_cxx_compiler_flag("-fno-tree-reassoc" COMPILER_SUPPORTS_NO_TREE_REASSOC)
    if(COMPILER_SUPPORTS_NO_TREE_REASSOC)
        set(FF8XOR_BASELINE_COMPILE_FLAGS
            "${FF8XOR_BASELINE_COMPILE_FLAGS} -fno-tree-reassoc")
        if(FF8XOR_BUILD_AVX2)
            set(FF8XOR_AVX2_COMPILE_FLAGS
                "${FF8XOR_AVX2_COMPILE_FLAGS} -fno-tree-reassoc")
        endif()
    endif()
endif()
if(FF8XOR_BASELINE_COMPILE_FLAGS)
    # Keep both public API gateways and the experimental dispatcher safe to
    # enter before any runtime feature check has succeeded.
    set_source_files_properties(
        leopard.cpp
        leopard_ff8xor.cpp
        LeopardFF8Xor.cpp
        LeopardFF8XorTranspose.cpp
        PROPERTIES
        COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
endif()

if(FF8XOR_BUILD_AVX2)
    set_source_files_properties(LeopardFF8XorAVX2.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX2_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "LEO_FF8XOR_BUILD_AVX2=1")
    set_source_files_properties(LeopardFF8XorTransposeAVX2.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX2_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "LEO_FF8XOR_BUILD_AVX2_TRANSPOSE=1")
    set_property(SOURCE LeopardFF8Xor.cpp APPEND PROPERTY
        COMPILE_DEFINITIONS "LEO_FF8XOR_HAS_AVX2_KERNELS=1")
    message(STATUS "FF8 XOR AVX2 kernels: enabled (isolated object)")
else()
    if(FF8XOR_BASELINE_COMPILE_FLAGS)
        set_source_files_properties(LeopardFF8XorAVX2.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
        set_source_files_properties(LeopardFF8XorTransposeAVX2.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
    endif()
    message(STATUS "FF8 XOR AVX2 kernels: unavailable")
endif()

# Preserve the old metadata variable name for benchmark output and downstream
# build scripts while making its baseline meaning explicit.
set(FF8XOR_COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")

# AVX-512 is kept in its own translation unit so the portable dispatcher and
# the AVX2 object cannot acquire instructions that are unsafe on older CPUs.
# Unsupported compilers build small stubs from this source; runtime
# CPUID/XCR0 checks still gate every call into a real kernel.
option(LEO_FF8XOR_ENABLE_AVX512
    "Build experimental AVX-512 FF8 XOR circuit and boundary-transpose kernels" ON)
set(FF8XOR_BUILD_AVX512 OFF)
if(LEO_FF8XOR_ENABLE_AVX512 AND FF8XOR_BUILD_AVX2 AND
   CMAKE_SIZEOF_VOID_P EQUAL 8 AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND
    (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)
    check_cxx_compiler_flag("-mavx512vl" COMPILER_SUPPORTS_AVX512VL)
    if(COMPILER_SUPPORTS_AVX512F AND
       COMPILER_SUPPORTS_AVX512VL)
        set(FF8XOR_BUILD_AVX512 ON)
    endif()
endif()

if(FF8XOR_BUILD_AVX512)
    # Reset the global -march=native before enabling the exact features checked
    # by the runtime dispatcher.  Otherwise an AVX-512-capable build host could
    # leak BW/DQ/VBMI/BITALG instructions into this isolated object.
    set(FF8XOR_AVX512_COMPILE_FLAGS
        "-march=x86-64 -mavx2 -mavx512f -mavx512vl")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
       COMPILER_SUPPORTS_NO_TREE_REASSOC)
        set(FF8XOR_AVX512_COMPILE_FLAGS
            "${FF8XOR_AVX512_COMPILE_FLAGS} -fno-tree-reassoc")
    endif()
    set_source_files_properties(LeopardFF8XorAVX512.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX512_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "LEO_FF8XOR_BUILD_AVX512=1")
    set_property(SOURCE LeopardFF8Xor.cpp APPEND PROPERTY
        COMPILE_DEFINITIONS "LEO_FF8XOR_HAS_AVX512_KERNELS=1")
    message(STATUS "FF8 XOR AVX-512VL/ZMM kernels: enabled")
else()
    if(FF8XOR_BASELINE_COMPILE_FLAGS)
        set_source_files_properties(LeopardFF8XorAVX512.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
    endif()
    message(STATUS "FF8 XOR AVX-512VL/ZMM kernels: unavailable")
endif()

# Packed-boundary transposes use two narrower AVX-512 feature contracts than
# the generated circuit kernels above.  BITALG and VBMI are compiled in
# separate objects so runtime dispatch never requires one feature merely
# because the other direction uses it.  AVX-512VL is not required.
set(FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE OFF)
set(FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE OFF)
set(FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS "")
set(FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS "")
if(LEO_FF8XOR_ENABLE_AVX512 AND
   CMAKE_SIZEOF_VOID_P EQUAL 8 AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND
   (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)
    check_cxx_compiler_flag("-mavx512bw"
        COMPILER_SUPPORTS_AVX512BW_TRANSPOSE)
    check_cxx_compiler_flag("-mavx512bitalg"
        COMPILER_SUPPORTS_AVX512BITALG_TRANSPOSE)
    check_cxx_compiler_flag("-mavx512vbmi"
        COMPILER_SUPPORTS_AVX512VBMI_TRANSPOSE)
    if(COMPILER_SUPPORTS_AVX512F AND
       COMPILER_SUPPORTS_AVX512BW_TRANSPOSE AND
       COMPILER_SUPPORTS_AVX512BITALG_TRANSPOSE)
        set(FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE ON)
        set(FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS
            "-march=x86-64")
    endif()
    if(COMPILER_SUPPORTS_AVX512F AND
       COMPILER_SUPPORTS_AVX512BW_TRANSPOSE AND
       COMPILER_SUPPORTS_AVX512VBMI_TRANSPOSE)
        set(FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE ON)
        set(FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS
            "-march=x86-64")
    endif()

    # A later -march does not necessarily cancel extension flags explicitly
    # supplied in global CXX flags.  Negate every probed non-contract feature
    # so adversarial build flags cannot broaden either runtime gate.  Append
    # the required feature set last: Clang otherwise treats some negative
    # options (notably -mno-fma) as revoking AVX-512F prerequisites.  Do not
    # negate AVX2 here: GCC models AVX-512F as depending on its target option
    # even though these inspected kernels emit no AVX2-width instructions.
    set(FF8XOR_AVX512_TRANSPOSE_COMMON_DISABLE_FLAGS
        -mno-avx512cd
        -mno-avx512dq
        -mno-avx512er
        -mno-avx512ifma
        -mno-avx512pf
        -mno-avx512vl
        -mno-avx512vpopcntdq
        -mno-avx512vnni
        -mno-avx512bf16
        -mno-avx512fp16
        -mno-avx512vp2intersect
        -mno-avx5124fmaps
        -mno-avx5124vnniw
        -mno-avxvnni
        -mno-avxvnniint8
        -mno-avxvnniint16
        -mno-gfni
        -mno-vaes
        -mno-vpclmulqdq
        -mno-bmi
        -mno-bmi2
        -mno-lzcnt
        -mno-popcnt
        -mno-aes
        -mno-pclmul
        -mno-fma)
    set(FF8XOR_AVX512_BITALG_TRANSPOSE_DISABLE_FLAGS
        ${FF8XOR_AVX512_TRANSPOSE_COMMON_DISABLE_FLAGS}
        -mno-avx512vbmi
        -mno-avx512vbmi2)
    set(FF8XOR_AVX512_VBMI_TRANSPOSE_DISABLE_FLAGS
        ${FF8XOR_AVX512_TRANSPOSE_COMMON_DISABLE_FLAGS}
        -mno-avx512bitalg
        -mno-avx512vbmi2)

    if(FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE)
        foreach(FF8XOR_TRANSPOSE_DISABLE_FLAG
                ${FF8XOR_AVX512_BITALG_TRANSPOSE_DISABLE_FLAGS})
            string(MAKE_C_IDENTIFIER
                "FF8XOR_TRANSPOSE_SUPPORTS_${FF8XOR_TRANSPOSE_DISABLE_FLAG}"
                FF8XOR_TRANSPOSE_DISABLE_VARIABLE)
            check_cxx_compiler_flag(
                "${FF8XOR_TRANSPOSE_DISABLE_FLAG}"
                ${FF8XOR_TRANSPOSE_DISABLE_VARIABLE})
            if(${FF8XOR_TRANSPOSE_DISABLE_VARIABLE})
                set(FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS
                    "${FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS} ${FF8XOR_TRANSPOSE_DISABLE_FLAG}")
            endif()
        endforeach()
        set(FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS
            "${FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS} -mavx512f -mavx512bw -mavx512bitalg")
    endif()
    if(FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE)
        foreach(FF8XOR_TRANSPOSE_DISABLE_FLAG
                ${FF8XOR_AVX512_VBMI_TRANSPOSE_DISABLE_FLAGS})
            string(MAKE_C_IDENTIFIER
                "FF8XOR_TRANSPOSE_SUPPORTS_${FF8XOR_TRANSPOSE_DISABLE_FLAG}"
                FF8XOR_TRANSPOSE_DISABLE_VARIABLE)
            check_cxx_compiler_flag(
                "${FF8XOR_TRANSPOSE_DISABLE_FLAG}"
                ${FF8XOR_TRANSPOSE_DISABLE_VARIABLE})
            if(${FF8XOR_TRANSPOSE_DISABLE_VARIABLE})
                set(FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS
                    "${FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS} ${FF8XOR_TRANSPOSE_DISABLE_FLAG}")
            endif()
        endforeach()
        set(FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS
            "${FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS} -mavx512f -mavx512bw -mavx512vbmi")
    endif()
endif()

if(FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE)
    set_source_files_properties(LeopardFF8XorTransposeBITALG.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS}"
        COMPILE_DEFINITIONS
            "LEO_FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE=1")
    set_property(SOURCE LeopardFF8XorTranspose.cpp APPEND PROPERTY
        COMPILE_DEFINITIONS
            "LEO_FF8XOR_HAS_AVX512_BITALG_TRANSPOSE=1")
    message(STATUS
        "FF8 XOR AVX-512BITALG packed-to-plane transpose: enabled")
else()
    if(FF8XOR_BASELINE_COMPILE_FLAGS)
        set_source_files_properties(LeopardFF8XorTransposeBITALG.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
    endif()
    message(STATUS
        "FF8 XOR AVX-512BITALG packed-to-plane transpose: unavailable")
endif()

if(FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE)
    set_source_files_properties(LeopardFF8XorTransposeVBMI.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS}"
        COMPILE_DEFINITIONS
            "LEO_FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE=1")
    set_property(SOURCE LeopardFF8XorTranspose.cpp APPEND PROPERTY
        COMPILE_DEFINITIONS
            "LEO_FF8XOR_HAS_AVX512_VBMI_TRANSPOSE=1")
    message(STATUS
        "FF8 XOR AVX-512VBMI plane-to-packed transpose: enabled")
else()
    if(FF8XOR_BASELINE_COMPILE_FLAGS)
        set_source_files_properties(LeopardFF8XorTransposeVBMI.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_BASELINE_COMPILE_FLAGS}")
    endif()
    message(STATUS
        "FF8 XOR AVX-512VBMI plane-to-packed transpose: unavailable")
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
    check_cxx_compiler_flag("-g" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Minform=inform" CXX_FLAG_MINFORM_INFORM)

    if(CXX_FLAG_g AND
       NOT " ${CMAKE_CXX_FLAGS_DEBUG} " MATCHES " -g ")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif()

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif(CXX_FLAG_O3)

    if(CXX_FLAG_MINFORM_INFORM)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
    endif(CXX_FLAG_MINFORM_INFORM)
else()
    check_cxx_compiler_flag("-g" CXX_FLAG_g)
    check_cxx_compiler_flag("-O0" CXX_FLAG_O0)
    check_cxx_compiler_flag("-O3" CXX_FLAG_O3)
    check_cxx_compiler_flag("-Wall" CXX_FLAG_Wall)
    check_cxx_compiler_flag("-Wextra" CXX_FLAG_Wextra)

    if(CXX_FLAG_g AND
       NOT " ${CMAKE_CXX_FLAGS_DEBUG} " MATCHES " -g ")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif()

    if(CXX_FLAG_O0)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    endif(CXX_FLAG_O0)

    if(CXX_FLAG_O3)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
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

add_test(NAME leopard_experiment_init
    COMMAND experiment_leopard --init-only)
set_tests_properties(leopard_experiment_init PROPERTIES TIMEOUT 30)

add_executable(test_leopard_ff8xor tests/test_ff8xor.cpp)
target_link_libraries(test_leopard_ff8xor libleopard)
add_test(NAME leopard_ff8xor COMMAND test_leopard_ff8xor)

add_executable(test_leopard_ff8xor_statistics tests/test_ff8xor_statistics.cpp)
target_link_libraries(test_leopard_ff8xor_statistics libleopard)
add_test(NAME leopard_ff8xor_statistics COMMAND test_leopard_ff8xor_statistics)

add_executable(test_leopard_api_counts tests/test_api_counts.cpp)
target_link_libraries(test_leopard_api_counts libleopard)
add_test(NAME leopard_api_counts COMMAND test_leopard_api_counts)

add_executable(bench_leopard_ff8xor tests/benchmark_ff8xor.cpp)
target_link_libraries(bench_leopard_ff8xor libleopard)

# These remain OFF/empty on builds without the isolated vector sources.
if(NOT DEFINED FF8XOR_BUILD_AVX2)
    set(FF8XOR_BUILD_AVX2 OFF)
endif()
if(NOT DEFINED FF8XOR_AVX2_COMPILE_FLAGS)
    set(FF8XOR_AVX2_COMPILE_FLAGS "")
endif()
if(NOT DEFINED FF8XOR_BUILD_AVX512)
    set(FF8XOR_BUILD_AVX512 OFF)
endif()
if(NOT DEFINED FF8XOR_AVX512_COMPILE_FLAGS)
    set(FF8XOR_AVX512_COMPILE_FLAGS "")
endif()

# Isolated AVX-512 XOR3 experiment.  The scalar driver is forced back to the
# baseline x86-64 ISA before runtime feature detection; only the second
# translation unit may contain VPTERNLOG or ZMM instructions.
add_executable(bench_leopard_ff8xor_xor3
    tests/benchmark_ff8xor_xor3.cpp
    tests/benchmark_ff8xor_xor3.h
    tests/benchmark_ff8xor_xor3_avx512.cpp)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND
   (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    # The trailing -march replaces the project-wide -march=native and is also
    # understood by older compilers that may not recognize -mno-avx512vl.
    set(FF8XOR_XOR3_BASELINE_FLAGS "-march=x86-64")
    set_source_files_properties(tests/benchmark_ff8xor_xor3.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_XOR3_BASELINE_FLAGS}")
    if(NOT FF8XOR_BUILD_AVX512)
        set_source_files_properties(
            tests/benchmark_ff8xor_xor3_avx512.cpp PROPERTIES
            COMPILE_FLAGS "${FF8XOR_XOR3_BASELINE_FLAGS}")
    endif()
endif()
if(FF8XOR_BUILD_AVX512)
    set_source_files_properties(tests/benchmark_ff8xor_xor3_avx512.cpp PROPERTIES
        COMPILE_FLAGS "${FF8XOR_AVX512_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "LEO_FF8XOR_XOR3_AVX512=1")
    # This definition only tells the baseline driver that the guarded entry
    # point exists.  It does not enable any vector ISA in that translation unit.
    set_source_files_properties(tests/benchmark_ff8xor_xor3.cpp PROPERTIES
        COMPILE_DEFINITIONS "LEO_FF8XOR_XOR3_AVX512_BUILT=1")
endif()
add_test(NAME leopard_ff8xor_xor3
    COMMAND bench_leopard_ff8xor_xor3 --verify-only)
set_tests_properties(leopard_ff8xor_xor3 PROPERTIES SKIP_RETURN_CODE 77)
add_test(NAME leopard_ff8xor_xor3_bad_cpu_empty
    COMMAND "${CMAKE_COMMAND}"
            -DPROGRAM=$<TARGET_FILE:bench_leopard_ff8xor_xor3>
            -DCASE=empty
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_xor3_bad_cpu.cmake")
add_test(NAME leopard_ff8xor_xor3_bad_cpu_whitespace
    COMMAND "${CMAKE_COMMAND}"
            -DPROGRAM=$<TARGET_FILE:bench_leopard_ff8xor_xor3>
            -DCASE=whitespace
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_xor3_bad_cpu.cmake")
add_test(NAME leopard_ff8xor_xor3_bad_cpu_overflow
    COMMAND "${CMAKE_COMMAND}"
            -DPROGRAM=$<TARGET_FILE:bench_leopard_ff8xor_xor3>
            -DCASE=overflow
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_xor3_bad_cpu.cmake")
add_test(NAME leopard_ff8xor_xor3_bad_cpu_erange
    COMMAND "${CMAKE_COMMAND}"
            -DPROGRAM=$<TARGET_FILE:bench_leopard_ff8xor_xor3>
            -DCASE=erange
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_xor3_bad_cpu.cmake")

# Embed a reproducible summary of the configured CMake flags in benchmark
# metadata.  The source-specific FF8 XOR flags are included explicitly; a
# verbose build remains authoritative for compiler-driver-added options.
# Select the actual configuration flags at generate/build time.  Multi-config
# generators use the project's fixed Debug/Release set.  Single-config builds
# retain dynamic lookup so custom configurations such as Coverage report their
# real CMAKE_CXX_FLAGS_<CONFIG> value rather than an empty string.
if(FF8XOR_MULTI_CONFIG)
    set(FF8XOR_BENCH_CONFIG_FLAGS
        "$<$<CONFIG:Debug>:${CMAKE_CXX_FLAGS_DEBUG}>$<$<CONFIG:Release>:${CMAKE_CXX_FLAGS_RELEASE}>")
else()
    set(FF8XOR_BENCH_CONFIG_FLAGS
        "${CMAKE_CXX_FLAGS_${FF8XOR_BUILD_TYPE_UPPER}}")
endif()
set(FF8XOR_BENCH_BUILD_FLAGS
    "global=[${CMAKE_CXX_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} config=[${FF8XOR_BENCH_CONFIG_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} cxx_standard=${CMAKE_CXX_STANDARD}")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} ff8xor=[${FF8XOR_COMPILE_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} avx2_built=${FF8XOR_BUILD_AVX2}")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} avx2=[${FF8XOR_AVX2_COMPILE_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} avx512_built=${FF8XOR_BUILD_AVX512}")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} avx512=[${FF8XOR_AVX512_COMPILE_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} transpose_bitalg_built=${FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE}")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} transpose_bitalg=[${FF8XOR_AVX512_BITALG_TRANSPOSE_FLAGS}]")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} transpose_vbmi_built=${FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE}")
set(FF8XOR_BENCH_BUILD_FLAGS
    "${FF8XOR_BENCH_BUILD_FLAGS} transpose_vbmi=[${FF8XOR_AVX512_VBMI_TRANSPOSE_FLAGS}]")
string(REPLACE "\\" "\\\\" FF8XOR_BENCH_BUILD_FLAGS_ESCAPED
    "${FF8XOR_BENCH_BUILD_FLAGS}")
string(REPLACE "\"" "\\\"" FF8XOR_BENCH_BUILD_FLAGS_ESCAPED
    "${FF8XOR_BENCH_BUILD_FLAGS_ESCAPED}")
target_compile_definitions(bench_leopard_ff8xor PRIVATE
    LEO_BENCH_CMAKE_BUILD_TYPE="$<CONFIG>"
    LEO_BENCH_BUILD_FLAGS="${FF8XOR_BENCH_BUILD_FLAGS_ESCAPED}")

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

set(FF8XOR_ASSEMBLY_REQUIRE_ARGS "")
if(FF8XOR_BUILD_AVX2)
    list(APPEND FF8XOR_ASSEMBLY_REQUIRE_ARGS --require-avx2)
endif()
if(FF8XOR_BUILD_AVX512)
    list(APPEND FF8XOR_ASSEMBLY_REQUIRE_ARGS --require-avx512)
endif()

# Sanitizers deliberately inject helper instructions and runtime calls that do
# not describe the production ISA contract (Clang ASan, for example, inserts
# LAHF into these objects).  Functional sanitizer coverage remains registered;
# run the opcode census on a corresponding unsanitized build.  Multi-config
# generators can instrument one configuration without changing the configure-
# time CMAKE_BUILD_TYPE, so inspect every configured flag set.  Conservatively
# omit artifact-shape targets when any configuration is instrumented: this
# avoids a Release target silently inspecting a Debug artifact, and parser unit
# tests remain available in every configuration.
set(FF8XOR_SANITIZER_INSTRUMENTED OFF)
set(FF8XOR_CONFIGURED_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if(FF8XOR_MULTI_CONFIG)
    foreach(FF8XOR_CONFIGURATION ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER "${FF8XOR_CONFIGURATION}"
            FF8XOR_CONFIGURATION_UPPER)
        set(FF8XOR_CONFIGURED_CXX_FLAGS
            "${FF8XOR_CONFIGURED_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${FF8XOR_CONFIGURATION_UPPER}}")
    endforeach()
else()
    set(FF8XOR_CONFIGURED_CXX_FLAGS
        "${FF8XOR_CONFIGURED_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${FF8XOR_BUILD_TYPE_UPPER}}")
endif()
if(FF8XOR_CONFIGURED_CXX_FLAGS MATCHES "fsanitize")
    set(FF8XOR_SANITIZER_INSTRUMENTED ON)
endif()

find_program(FF8XOR_NM_EXECUTABLE NAMES llvm-nm nm)
find_program(FF8XOR_OBJDUMP_EXECUTABLE NAMES llvm-objdump objdump)
find_program(FF8XOR_AR_EXECUTABLE NAMES llvm-ar ar)
set(FF8XOR_CAN_INSPECT_X86_ASSEMBLY OFF)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND
   (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang") AND
   FF8XOR_NM_EXECUTABLE AND FF8XOR_OBJDUMP_EXECUTABLE AND
   FF8XOR_AR_EXECUTABLE)
    set(FF8XOR_CAN_INSPECT_X86_ASSEMBLY ON)
endif()

if(FF8XOR_PYTHON_EXECUTABLE)
    add_custom_target(generate_ff8_xor_circuits
        COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8_xor_circuits.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_custom_target(check_ff8_xor_circuits
        COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8_xor_circuits.py"
                --check
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_test(NAME ff8xor_circuits
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8_xor_circuits.py"
                --check)

    add_custom_target(generate_ff8xor_schedule_corpus
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8xor_schedule_corpus.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_custom_target(check_ff8xor_schedule_corpus
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8xor_schedule_corpus.py"
                --check
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_test(NAME ff8xor_schedule_corpus
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8xor_schedule_corpus.py"
                --check)

    add_custom_target(generate_ff8xor_cost_profiles
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8xor_cost_profiles.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_custom_target(check_ff8xor_cost_profiles
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ff8xor_cost_profiles.py"
                --check
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        VERBATIM)

    add_test(NAME ff8xor_cost_model
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8_xor_cost_model.py")

    add_test(NAME ff8xor_cost_profiles
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_cost_profiles.py")

    add_test(NAME ff8xor_cost_generator_selection
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_cost_generator_selection.py")

    add_test(NAME ff8xor_cost_evaluator
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_cost_evaluator.py")

    add_test(NAME ff8xor_cost_evidence
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_cost_evidence.py")

    add_test(NAME ff8xor_benchmark_metadata
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_benchmark_metadata.py"
                "$<TARGET_FILE:bench_leopard_ff8xor>")

    add_test(NAME ff8xor_assembly_inspector
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_assembly_inspector.py")

    add_test(NAME ff8xor_xor_batch_inspector
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_xor_batch_inspector.py")

    add_test(NAME ff8xor_synthesis
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8_xor_synthesis.py")

    add_test(NAME ff8xor_xor3_schedule
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8_xor3_schedule.py")

    add_test(NAME ff8xor_baseline_isa_inspector
        COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                "${FF8XOR_PYTHON_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_baseline_isa_inspector.py")

    if(FF8XOR_CAN_INSPECT_X86_ASSEMBLY)
        add_test(NAME ff8xor_cmake_sanitizer_gating
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ff8xor_cmake_sanitizer_gating.py")
    endif()

    if(FF8XOR_CAN_INSPECT_X86_ASSEMBLY AND
       NOT FF8XOR_SANITIZER_INSTRUMENTED)
        # Debug and less aggressive compilers may leave baseline-attributed
        # initialization helpers out of line.  Keep the complete helper set in
        # the census instead of relying on Release inlining to hide omissions.
        set(FF8XOR_BASELINE_COMMON_SYMBOLS
            "leopard::(SIMDSafe(Allocate|Free)|InitializeCPUArch|IsAVX2Supported|_cpuid|_xgetbv0)")
        set(FF8XOR_BASELINE_FF8_SYMBOLS
            "leopard::(SIMDSafe(Allocate|Free)|ff8::(AddMod|SubMod|MultiplyLog|FWHT|FFTInitialize|Initialize))")
        set(FF8XOR_BASELINE_FF16_SYMBOLS
            "leopard::(SIMDSafe(Allocate|Free)|ff16::(AddMod|SubMod|MultiplyLog|FWHT|FFTInitialize|Initialize))")

        add_test(NAME ff8xor_baseline_isa_census
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_baseline_isa.py"
                    "$<TARGET_FILE:libleopard>"
                    --member leopard.cpp.o
                    --member LeopardCommon.cpp.o
                    --member LeopardFF8.cpp.o
                    --member LeopardFF16.cpp.o
                    --member LeopardFF8Xor.cpp.o
                    --member LeopardFF8XorTranspose.cpp.o
                    --member leopard_ff8xor.cpp.o
                    --member-symbol
                        "LeopardCommon.cpp.o=${FF8XOR_BASELINE_COMMON_SYMBOLS}"
                    --member-symbol
                        "LeopardFF8.cpp.o=${FF8XOR_BASELINE_FF8_SYMBOLS}"
                    --member-symbol
                        "LeopardFF16.cpp.o=${FF8XOR_BASELINE_FF16_SYMBOLS}"
                    --ar "${FF8XOR_AR_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict)

        add_custom_target(check_ff8xor_baseline_isa
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_baseline_isa.py"
                    "$<TARGET_FILE:libleopard>"
                    --member leopard.cpp.o
                    --member LeopardCommon.cpp.o
                    --member LeopardFF8.cpp.o
                    --member LeopardFF16.cpp.o
                    --member LeopardFF8Xor.cpp.o
                    --member LeopardFF8XorTranspose.cpp.o
                    --member leopard_ff8xor.cpp.o
                    --member-symbol
                        "LeopardCommon.cpp.o=${FF8XOR_BASELINE_COMMON_SYMBOLS}"
                    --member-symbol
                        "LeopardFF8.cpp.o=${FF8XOR_BASELINE_FF8_SYMBOLS}"
                    --member-symbol
                        "LeopardFF16.cpp.o=${FF8XOR_BASELINE_FF16_SYMBOLS}"
                    --ar "${FF8XOR_AR_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
                    --json
            DEPENDS libleopard
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            VERBATIM)

        add_test(NAME ff8xor_baseline_isa_custom_target
            COMMAND "${CMAKE_COMMAND}" --build
                    "${CMAKE_CURRENT_BINARY_DIR}"
                    ${FF8XOR_BUILD_CONFIG_ARGS}
                    --target check_ff8xor_baseline_isa)
    elseif(FF8XOR_CAN_INSPECT_X86_ASSEMBLY)
        message(STATUS
            "FF8 XOR baseline ISA census: skipped for sanitizer-instrumented objects")
    endif()

    if(FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE AND
       FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE AND
       FF8XOR_CAN_INSPECT_X86_ASSEMBLY AND
       NOT FF8XOR_SANITIZER_INSTRUMENTED AND
       (FF8XOR_MULTI_CONFIG OR
        FF8XOR_BUILD_TYPE_UPPER STREQUAL "RELEASE"))
        add_test(NAME ff8xor_transpose_avx512_assembly
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_transpose_avx512.py"
                    "$<TARGET_FILE:libleopard>"
                    --ar "${FF8XOR_AR_EXECUTABLE}"
                    --nm "${FF8XOR_NM_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
            ${FF8XOR_RELEASE_TEST_CONFIG_ARGS})

        add_custom_target(check_ff8xor_transpose_avx512_assembly
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_transpose_avx512.py"
                    "$<TARGET_FILE:libleopard>"
                    --ar "${FF8XOR_AR_EXECUTABLE}"
                    --nm "${FF8XOR_NM_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
                    --json
            DEPENDS libleopard
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            VERBATIM)
    endif()

    if(FF8XOR_BUILD_AVX512)
        if(FF8XOR_CAN_INSPECT_X86_ASSEMBLY AND
           NOT FF8XOR_SANITIZER_INSTRUMENTED)
            # Register the spill/code-shape gate only for optimized builds.
            # Debug deliberately materializes vector arguments on the stack.
            if(FF8XOR_MULTI_CONFIG OR
               FF8XOR_BUILD_TYPE_UPPER STREQUAL "RELEASE")
                add_test(NAME ff8xor_xor3_assembly
                    COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                            "${FF8XOR_PYTHON_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_xor3.py"
                            "$<TARGET_FILE:bench_leopard_ff8xor_xor3>"
                            --nm "${FF8XOR_NM_EXECUTABLE}"
                            --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                            --strict
                    ${FF8XOR_RELEASE_TEST_CONFIG_ARGS})

                add_custom_target(check_ff8xor_xor3_assembly
                    COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                            "${FF8XOR_PYTHON_EXECUTABLE}"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_xor3.py"
                            "$<TARGET_FILE:bench_leopard_ff8xor_xor3>"
                            --nm "${FF8XOR_NM_EXECUTABLE}"
                            --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                            --strict
                            --json
                    DEPENDS bench_leopard_ff8xor_xor3
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    VERBATIM)
            endif()
        elseif(NOT FF8XOR_CAN_INSPECT_X86_ASSEMBLY)
            message(STATUS
                "FF8 XOR3 assembly inspection: nm/objdump unavailable")
        endif()
    endif()

    if(NOT FF8XOR_SANITIZER_INSTRUMENTED)
        add_custom_target(inspect_ff8xor_assembly
            COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_assembly.py"
                    "$<TARGET_FILE:libleopard>"
                    --json
                    --output "${CMAKE_CURRENT_BINARY_DIR}/ff8xor_assembly_census.json"
                    ${FF8XOR_ASSEMBLY_REQUIRE_ARGS}
            DEPENDS libleopard
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            VERBATIM)

        add_custom_target(check_ff8xor_assembly
            COMMAND "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_assembly.py"
                    "$<TARGET_FILE:libleopard>"
                    --strict
                    --fail-on-spills
                    ${FF8XOR_ASSEMBLY_REQUIRE_ARGS}
            DEPENDS libleopard
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            VERBATIM)
    endif()

    if(FF8XOR_CAN_INSPECT_X86_ASSEMBLY AND FF8XOR_BUILD_AVX2 AND
       NOT FF8XOR_SANITIZER_INSTRUMENTED AND
       (FF8XOR_MULTI_CONFIG OR
        FF8XOR_BUILD_TYPE_UPPER STREQUAL "RELEASE"))
        add_test(NAME ff8xor_xor_batch_assembly
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_xor_batch.py"
                    "$<TARGET_FILE:libleopard>"
                    --nm "${FF8XOR_NM_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
                    ${FF8XOR_ASSEMBLY_REQUIRE_ARGS}
                    ${FF8XOR_RELEASE_TEST_CONFIG_ARGS})

        add_custom_target(check_ff8xor_xor_batch_assembly
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_xor_batch.py"
                    "$<TARGET_FILE:libleopard>"
                    --nm "${FF8XOR_NM_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
                    ${FF8XOR_ASSEMBLY_REQUIRE_ARGS}
            DEPENDS libleopard
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            VERBATIM)
    endif()

    if(FF8XOR_CAN_INSPECT_X86_ASSEMBLY AND
       NOT FF8XOR_SANITIZER_INSTRUMENTED AND
       (FF8XOR_MULTI_CONFIG OR
        FF8XOR_BUILD_TYPE_UPPER STREQUAL "RELEASE"))
        add_test(NAME ff8xor_assembly_census
            COMMAND "${CMAKE_COMMAND}" -E env "PYTHONDONTWRITEBYTECODE=1"
                    "${FF8XOR_PYTHON_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/inspect_ff8xor_assembly.py"
                    "$<TARGET_FILE:libleopard>"
                    --nm "${FF8XOR_NM_EXECUTABLE}"
                    --objdump "${FF8XOR_OBJDUMP_EXECUTABLE}"
                    --strict
                    ${FF8XOR_ASSEMBLY_REQUIRE_ARGS}
                    ${FF8XOR_RELEASE_TEST_CONFIG_ARGS})
    endif()

    if(FF8XOR_SANITIZER_INSTRUMENTED)
        message(STATUS
            "FF8 XOR production assembly checks: skipped for sanitizer-instrumented objects")
    endif()
endif()
