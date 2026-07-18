if(NOT DEFINED OUTPUT OR NOT DEFINED EXECUTABLE OR
   NOT DEFINED PRODUCTION_ARCHIVE OR NOT DEFINED BENCHMARK_OBJECT OR
   NOT DEFINED ORACLE_OBJECT)
    message(FATAL_ERROR "sparse evidence sidecar arguments are incomplete")
endif()

foreach(required_path "${EXECUTABLE}" "${PRODUCTION_ARCHIVE}"
                      "${BENCHMARK_OBJECT}" "${ORACLE_OBJECT}")
    if(NOT EXISTS "${required_path}" OR IS_DIRECTORY "${required_path}")
        message(FATAL_ERROR
            "sparse evidence sidecar input is not a file: ${required_path}")
    endif()
endforeach()

file(SHA256 "${EXECUTABLE}" executable_sha256)
file(SIZE "${EXECUTABLE}" executable_size)
file(SHA256 "${PRODUCTION_ARCHIVE}" production_archive_sha256)
file(SIZE "${PRODUCTION_ARCHIVE}" production_archive_size)
file(SHA256 "${BENCHMARK_OBJECT}" benchmark_object_sha256)
file(SIZE "${BENCHMARK_OBJECT}" benchmark_object_size)
file(SHA256 "${ORACLE_OBJECT}" oracle_object_sha256)
file(SIZE "${ORACLE_OBJECT}" oracle_object_size)

function(hash_optional_file path prefix)
    if(EXISTS "${path}" AND NOT IS_DIRECTORY "${path}")
        file(SHA256 "${path}" digest)
        file(SIZE "${path}" size)
        set(${prefix}_sha256 "${digest}" PARENT_SCOPE)
        set(${prefix}_size "${size}" PARENT_SCOPE)
    else()
        # Some CMake generators do not materialize link.txt.  The benchmark
        # remains usable for non-authoritative screens, while the pinned
        # runner rejects a sidecar without an attested generator graph.
        set(${prefix}_sha256 "missing" PARENT_SCOPE)
        set(${prefix}_size "missing" PARENT_SCOPE)
    endif()
endfunction()

if(EXISTS "${BENCHMARK_LINK_RECIPE}" AND
   EXISTS "${PRODUCTION_LINK_RECIPE}")
    set(link_recipe_kind "cmake-link-txt-v1")
    hash_optional_file("${BENCHMARK_LINK_RECIPE}" benchmark_link_recipe)
    hash_optional_file("${PRODUCTION_LINK_RECIPE}" production_link_recipe)
elseif(BUILD_GENERATOR MATCHES "Ninja" AND DEFINED BUILD_PROGRAM AND
       DEFINED BUILD_ROOT)
    set(link_recipe_kind "ninja-tool-commands-v1")
    set(benchmark_commands "${OUTPUT}.benchmark-commands.tmp")
    set(production_commands "${OUTPUT}.production-commands.tmp")
    execute_process(
        COMMAND "${BUILD_PROGRAM}" -C "${BUILD_ROOT}" -t commands
                bench_leopard2_sparse_encode
        RESULT_VARIABLE benchmark_result
        OUTPUT_FILE "${benchmark_commands}")
    execute_process(
        COMMAND "${BUILD_PROGRAM}" -C "${BUILD_ROOT}" -t commands leopard
        RESULT_VARIABLE production_result
        OUTPUT_FILE "${production_commands}")
    if(NOT benchmark_result EQUAL 0 OR NOT production_result EQUAL 0)
        file(REMOVE "${benchmark_commands}" "${production_commands}")
        message(FATAL_ERROR "cannot capture Ninja link graph commands")
    endif()
    hash_optional_file("${benchmark_commands}" benchmark_link_recipe)
    hash_optional_file("${production_commands}" production_link_recipe)
    file(REMOVE "${benchmark_commands}" "${production_commands}")
else()
    set(link_recipe_kind "missing")
    set(benchmark_link_recipe_sha256 "missing")
    set(benchmark_link_recipe_size "missing")
    set(production_link_recipe_sha256 "missing")
    set(production_link_recipe_size "missing")
endif()

file(WRITE "${OUTPUT}"
    "schema=leopard2-sparse-encode-link-sidecar/v1\n"
    "executable_sha256=${executable_sha256}\n"
    "executable_size=${executable_size}\n"
    "production_archive_sha256=${production_archive_sha256}\n"
    "production_archive_size=${production_archive_size}\n"
    "benchmark_object_sha256=${benchmark_object_sha256}\n"
    "benchmark_object_size=${benchmark_object_size}\n"
    "oracle_object_sha256=${oracle_object_sha256}\n"
    "oracle_object_size=${oracle_object_size}\n"
    "link_recipe_kind=${link_recipe_kind}\n"
    "benchmark_link_recipe_sha256=${benchmark_link_recipe_sha256}\n"
    "benchmark_link_recipe_size=${benchmark_link_recipe_size}\n"
    "production_link_recipe_sha256=${production_link_recipe_sha256}\n"
    "production_link_recipe_size=${production_link_recipe_size}\n")
