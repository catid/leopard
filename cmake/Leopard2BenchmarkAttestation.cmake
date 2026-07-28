if(COMMAND leopard2_enable_benchmark_source_attestation)
    return()
endif()

set(_LEO2_BENCHMARK_ATTESTATION_MODULE_DIR
    "${CMAKE_CURRENT_LIST_DIR}")

function(leopard2_enable_benchmark_source_attestation target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR
            "Cannot attest missing Leopard2 benchmark target '${target}'")
    endif()

    if(ARGC GREATER 1)
        set(source_dir "${ARGV1}")
    else()
        set(source_dir "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
    get_filename_component(source_dir "${source_dir}" ABSOLUTE)

    find_program(LEO2_BENCHMARK_GIT_EXECUTABLE NAMES git)

    set(build_configuration_material "")
    foreach(variable
            CMAKE_BUILD_TYPE
            CMAKE_CXX_COMPILER
            CMAKE_CXX_FLAGS
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS_MINSIZEREL
            ENABLE_OPENMP
            LEOPARD_ENABLE_GF8
            LEOPARD_ENABLE_GF16
            LEO2_BACKEND_VARIANT
            LEO2_BENCHMARK_GIT_EXECUTABLE
            LEO2_BUILD_BENCHMARKS
            LEO2_BUILD_TESTS
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE)
        string(APPEND build_configuration_material
            "${variable}=${${variable}}\n")
    endforeach()
    string(SHA256 build_configuration_sha256
        "${build_configuration_material}")

    set(output_dir
        "${CMAKE_BINARY_DIR}/generated/leopard2-benchmark-attestation")
    set(output_file
        "${output_dir}/leopard2_benchmark_source_attestation.h")
    set(refresh_target
        leopard2_benchmark_source_attestation_refresh)

    if(NOT TARGET "${refresh_target}")
        add_custom_target("${refresh_target}"
            COMMAND "${CMAKE_COMMAND}" -E env
                GIT_CONFIG_GLOBAL=/dev/null
                GIT_CONFIG_NOSYSTEM=1
                GIT_CONFIG_SYSTEM=/dev/null
                GIT_NO_REPLACE_OBJECTS=1
                GIT_OPTIONAL_LOCKS=0
                "${CMAKE_COMMAND}"
                "-DLEO2_SOURCE_DIR=${source_dir}"
                "-DLEO2_OUTPUT_FILE=${output_file}"
                "-DLEO2_GIT_EXECUTABLE=${LEO2_BENCHMARK_GIT_EXECUTABLE}"
                -P
                "${_LEO2_BENCHMARK_ATTESTATION_MODULE_DIR}/GenerateBenchmarkSourceAttestation.cmake"
            BYPRODUCTS "${output_file}"
            COMMENT "Refreshing Leopard2 benchmark source attestation"
            VERBATIM)
        set_source_files_properties("${output_file}" PROPERTIES
            GENERATED TRUE
            HEADER_FILE_ONLY TRUE)
    endif()

    # A top-level target dependency is intentional: the always-run refresh
    # must complete before the generator decides whether benchmark.cpp and its
    # generated-header dependency are up to date.
    add_dependencies("${target}" "${refresh_target}")
    target_sources("${target}" PRIVATE "${output_file}")
    target_compile_definitions("${target}" PRIVATE
        LEO2_BENCHMARK_SOURCE_ATTESTATION=1
        "LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\"${build_configuration_sha256}\""
        "LEO2_BENCHMARK_BUILD_TYPE=\"${CMAKE_BUILD_TYPE}\""
        "LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER=\"${output_file}\"")

    set(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER
        "${output_file}" PARENT_SCOPE)
endfunction()
