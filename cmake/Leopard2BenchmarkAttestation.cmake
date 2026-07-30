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
    get_filename_component(source_dir "${source_dir}" REALPATH)
    get_property(previous_source_dir GLOBAL PROPERTY
        LEO2_BENCHMARK_ATTESTED_SOURCE_DIR)
    if(previous_source_dir AND NOT previous_source_dir STREQUAL source_dir)
        message(FATAL_ERROR
            "Leopard2 benchmark targets requested inconsistent attested "
            "source roots: '${previous_source_dir}' and '${source_dir}'")
    endif()
    set_property(GLOBAL PROPERTY LEO2_BENCHMARK_ATTESTED_SOURCE_DIR
        "${source_dir}")

    find_program(LEO2_BENCHMARK_GIT_EXECUTABLE NAMES git)

    set(build_configuration_material "")
    foreach(variable
            CMAKE_BUILD_TYPE
            CMAKE_GENERATOR
            CMAKE_CONFIGURATION_TYPES
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
            LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN
            LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE
            LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT
            LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE)
        set(variable_value "${${variable}}")
        if(variable STREQUAL "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE")
            if(NOT DEFINED LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE)
                set(variable_value "0")
            elseif(NOT LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
                    MATCHES "^[012]$")
                message(FATAL_ERROR
                    "Cannot attest "
                    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE: expected "
                    "exactly 0, 1, or 2")
            endif()
        endif()
        string(FIND "${variable_value}" "\n" newline_index)
        string(FIND "${variable_value}" "\r" carriage_return_index)
        if(NOT newline_index EQUAL -1 OR
           NOT carriage_return_index EQUAL -1)
            message(FATAL_ERROR
                "Cannot attest ${variable}: effective build-configuration "
                "values must not contain CR or LF characters")
        endif()
        string(APPEND build_configuration_material
            "${variable}=${variable_value}\n")
    endforeach()
    string(SHA256 build_configuration_sha256
        "${build_configuration_material}")
    get_property(previous_build_configuration_sha256 GLOBAL PROPERTY
        LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256)
    if(previous_build_configuration_sha256 AND
       NOT previous_build_configuration_sha256 STREQUAL
           build_configuration_sha256)
        message(FATAL_ERROR
            "Leopard2 benchmark targets observed different effective "
            "build configurations")
    endif()
    set_property(GLOBAL PROPERTY
        LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256
        "${build_configuration_sha256}")
    set(LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA
        "leopard2-benchmark-build-configuration/v3"
        CACHE INTERNAL
        "Leopard2 benchmark effective-configuration schema"
        FORCE)
    set(LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256
        "${build_configuration_sha256}"
        CACHE INTERNAL
        "Leopard2 benchmark effective-configuration digest"
        FORCE)

    set(output_dir
        "${CMAKE_BINARY_DIR}/generated/leopard2-benchmark-attestation")
    set(output_parent "${CMAKE_BINARY_DIR}/generated")
    foreach(attestation_directory
            "${CMAKE_BINARY_DIR}" "${output_parent}" "${output_dir}")
        # Check each existing ancestor before creating its child.  Checking
        # only after MAKE_DIRECTORY can already have followed a generated/
        # symlink and mutated an unrelated directory.
        if(IS_SYMLINK "${attestation_directory}")
            message(FATAL_ERROR
                "Leopard2 benchmark attestation directory must not traverse "
                "a symbolic link: '${attestation_directory}'")
        endif()
        if(EXISTS "${attestation_directory}")
            if(NOT IS_DIRECTORY "${attestation_directory}")
                message(FATAL_ERROR
                    "Leopard2 benchmark attestation path is not a directory: "
                    "'${attestation_directory}'")
            endif()
        else()
            file(MAKE_DIRECTORY "${attestation_directory}")
        endif()
        get_filename_component(attestation_directory_absolute
            "${attestation_directory}" ABSOLUTE)
        get_filename_component(attestation_directory_realpath
            "${attestation_directory}" REALPATH)
        if(IS_SYMLINK "${attestation_directory}" OR
           NOT attestation_directory_absolute STREQUAL
               attestation_directory_realpath)
            message(FATAL_ERROR
                "Leopard2 benchmark attestation directory must not traverse "
                "a symbolic link: '${attestation_directory}'")
        endif()
    endforeach()
    set(build_configuration_file
        "${output_dir}/leopard2_benchmark_build_configuration.txt")
    if(IS_SYMLINK "${build_configuration_file}")
        message(FATAL_ERROR
            "Leopard2 benchmark effective-configuration sidecar must not be "
            "a symbolic link")
    endif()
    set(build_configuration_lock_file
        "${build_configuration_file}.lock")
    if(IS_SYMLINK "${build_configuration_lock_file}")
        message(FATAL_ERROR
            "Leopard2 benchmark effective-configuration lock must not be a "
            "symbolic link")
    endif()
    file(LOCK "${build_configuration_lock_file}"
        GUARD FUNCTION
        TIMEOUT 120
        RESULT_VARIABLE build_configuration_lock_result)
    if(NOT build_configuration_lock_result STREQUAL "0")
        message(FATAL_ERROR
            "Cannot lock Leopard2 benchmark effective-configuration "
            "sidecar: ${build_configuration_lock_result}")
    endif()
    string(RANDOM LENGTH 32 ALPHABET 0123456789abcdef
        build_configuration_temporary_suffix)
    set(build_configuration_temporary
        "${build_configuration_file}.${build_configuration_temporary_suffix}.tmp")
    if(EXISTS "${build_configuration_temporary}" OR
       IS_SYMLINK "${build_configuration_temporary}")
        message(FATAL_ERROR
            "Refusing to use an existing Leopard2 benchmark "
            "effective-configuration temporary")
    endif()
    file(WRITE "${build_configuration_temporary}"
        "schema=${LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA}\n"
        "sha256=${build_configuration_sha256}\n"
        "${build_configuration_material}")
    file(SHA256 "${build_configuration_temporary}"
        new_build_configuration_file_sha256)
    set(build_configuration_changed TRUE)
    if(EXISTS "${build_configuration_file}")
        file(SHA256 "${build_configuration_file}"
            old_build_configuration_file_sha256)
        if(new_build_configuration_file_sha256 STREQUAL
           old_build_configuration_file_sha256)
            set(build_configuration_changed FALSE)
        endif()
    endif()
    if(build_configuration_changed)
        file(RENAME "${build_configuration_temporary}"
            "${build_configuration_file}")
    else()
        file(REMOVE "${build_configuration_temporary}")
    endif()
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
    if(CMAKE_GENERATOR MATCHES "Visual Studio|Xcode|Multi-Config")
        set(benchmark_build_type "$<CONFIG>")
    else()
        set(benchmark_build_type "${CMAKE_BUILD_TYPE}")
    endif()
    target_compile_definitions("${target}" PRIVATE
        LEO2_BENCHMARK_SOURCE_ATTESTATION=1
        "LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\"${build_configuration_sha256}\""
        "LEO2_BENCHMARK_BUILD_TYPE=\"${benchmark_build_type}\""
        "LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER=\"${output_file}\"")

    set(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER
        "${output_file}" PARENT_SCOPE)
endfunction()
