if(NOT DEFINED LEO2_SOURCE_DIR)
    message(FATAL_ERROR
        "Sanitizer-classification test is missing LEO2_SOURCE_DIR")
endif()

set(classifier
    "${LEO2_SOURCE_DIR}/cmake/Leopard2SanitizerClassification.cmake")
if(NOT EXISTS "${classifier}")
    message(FATAL_ERROR "Sanitizer classifier is missing: ${classifier}")
endif()

# CMake treats values such as OFF and 0 as false constants.  The project must
# default only an actually empty single-config build type, or those legitimate
# custom configuration names are rewritten to Release before this classifier
# sees them.
file(READ "${LEO2_SOURCE_DIR}/CMakeLists.txt" project_cmake)
string(FIND "${project_cmake}"
    "if(\"\${CMAKE_BUILD_TYPE}\" STREQUAL \"\")\n    set(CMAKE_BUILD_TYPE Release)\nendif()"
    explicit_empty_build_type_guard)
if(explicit_empty_build_type_guard EQUAL -1)
    message(FATAL_ERROR
        "Project build-type default does not preserve false-like custom "
        "configuration names")
endif()

function(leopard2_expect_lists case_name expected_sanitized expected_clean)
    if(NOT "${LEO2_SANITIZED_CONFIGURATIONS}" STREQUAL
           "${expected_sanitized}")
        message(FATAL_ERROR
            "${case_name}: sanitized configurations differ: "
            "'${LEO2_SANITIZED_CONFIGURATIONS}' != '${expected_sanitized}'")
    endif()
    if(NOT "${LEO2_UNSANITIZED_CONFIGURATIONS}" STREQUAL
           "${expected_clean}")
        message(FATAL_ERROR
            "${case_name}: clean configurations differ: "
            "'${LEO2_UNSANITIZED_CONFIGURATIONS}' != '${expected_clean}'")
    endif()
endfunction()

function(leopard2_case_clean_release)
    set(CMAKE_GENERATOR "Ninja")
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    include("${classifier}")
    leopard2_expect_lists(clean-release "" "Release")
    if(LEO2_PRODUCTION_ARCHIVE_SANITIZED OR LEO2_MULTI_CONFIG)
        message(FATAL_ERROR "clean-release: scalar classification differs")
    endif()
endfunction()

function(leopard2_case_custom_asan)
    set(CMAKE_GENERATOR "Ninja")
    set(CMAKE_BUILD_TYPE "ASan")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_ASAN "-O1 -fsanitize=address,undefined")
    # A different known configuration must not influence the active custom one.
    set(CMAKE_CXX_FLAGS_RELEASE "-fsanitize=thread")
    include("${classifier}")
    leopard2_expect_lists(custom-asan "ASan" "")
    if(NOT LEO2_PRODUCTION_ARCHIVE_SANITIZED OR LEO2_MULTI_CONFIG)
        message(FATAL_ERROR "custom-asan: scalar classification differs")
    endif()
endfunction()

function(leopard2_case_multi_mixed)
    set(CMAKE_GENERATOR "Ninja Multi-Config")
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release;ASan")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_ASAN "-O1 -fsanitize=address")
    include("${classifier}")
    leopard2_expect_lists(multi-mixed "ASan" "Debug;Release")
    if(LEO2_PRODUCTION_ARCHIVE_SANITIZED OR NOT LEO2_MULTI_CONFIG)
        message(FATAL_ERROR "multi-mixed: scalar classification differs")
    endif()
endfunction()

function(leopard2_case_multi_global_sanitizer)
    set(CMAKE_GENERATOR "Xcode")
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
    set(CMAKE_CXX_FLAGS "/fsanitize=address")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    include("${classifier}")
    leopard2_expect_lists(multi-global "Debug;Release" "")
    if(NOT LEO2_PRODUCTION_ARCHIVE_SANITIZED OR NOT LEO2_MULTI_CONFIG)
        message(FATAL_ERROR "multi-global: scalar classification differs")
    endif()
endfunction()

function(leopard2_case_non_sanitizer_spelling)
    set(CMAKE_GENERATOR "Unix Makefiles")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS "-fno-sanitize-recover=all")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2")
    include("${classifier}")
    leopard2_expect_lists(non-sanitizer-spelling "" "RelWithDebInfo")
endfunction()

function(leopard2_case_quoted_sanitizers)
    set(CMAKE_GENERATOR "Ninja Multi-Config")
    set(CMAKE_CONFIGURATION_TYPES "DoubleQuoted;SingleQuoted;Clean")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DOUBLEQUOTED
        "-O1 \"-fsanitize=address,undefined\"")
    set(CMAKE_CXX_FLAGS_SINGLEQUOTED
        "-O1 '/fsanitize=address'")
    set(CMAKE_CXX_FLAGS_CLEAN "-O3")
    include("${classifier}")
    leopard2_expect_lists(
        quoted-sanitizers "DoubleQuoted;SingleQuoted" "Clean")
    if(LEO2_PRODUCTION_ARCHIVE_SANITIZED OR NOT LEO2_MULTI_CONFIG)
        message(FATAL_ERROR
            "quoted-sanitizers: mixed classification differs")
    endif()
endfunction()

function(leopard2_case_empty_build_type)
    set(CMAKE_GENERATOR "Ninja")
    set(CMAKE_BUILD_TYPE "")
    set(CMAKE_CXX_FLAGS "-O1 -fsanitize=undefined")
    # No configuration-specific flag is active when the build type is empty.
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    include("${classifier}")
    leopard2_expect_lists(
        empty-build-type "__LEO2_DEFAULT__" "")
    if(NOT LEO2_PRODUCTION_ARCHIVE_SANITIZED OR LEO2_MULTI_CONFIG)
        message(FATAL_ERROR "empty-build-type: scalar classification differs")
    endif()
endfunction()

function(leopard2_case_default_sentinel_name_collision)
    set(CMAKE_GENERATOR "Ninja")
    set(CMAKE_BUILD_TYPE "__LEO2_DEFAULT__")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS___LEO2_DEFAULT__
        "-O1 -fsanitize=address")
    include("${classifier}")
    leopard2_expect_lists(
        default-sentinel-name-collision "__LEO2_DEFAULT__" "")
    if(NOT LEO2_PRODUCTION_ARCHIVE_SANITIZED OR LEO2_MULTI_CONFIG OR
       LEO2_AUDIT_USES_DEFAULT_FLAGS)
        message(FATAL_ERROR
            "default-sentinel-name-collision: custom configuration was "
            "treated as the empty-build-type sentinel")
    endif()
endfunction()

function(leopard2_case_false_named_custom_configuration)
    set(CMAKE_GENERATOR "Ninja")
    set(CMAKE_BUILD_TYPE "OFF")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_OFF "-O1 -fsanitize=address")
    include("${classifier}")
    leopard2_expect_lists(false-named-custom-configuration "OFF" "")
    if(NOT LEO2_PRODUCTION_ARCHIVE_SANITIZED OR LEO2_MULTI_CONFIG OR
       LEO2_AUDIT_USES_DEFAULT_FLAGS OR
       NOT LEO2_SANITIZED_CONFIGURATION_COUNT EQUAL 1 OR
       NOT LEO2_UNSANITIZED_CONFIGURATION_COUNT EQUAL 0)
        message(FATAL_ERROR
            "false-named-custom-configuration: CMake false constant was "
            "treated as empty or omitted from the classification")
    endif()
endfunction()

leopard2_case_clean_release()
leopard2_case_custom_asan()
leopard2_case_multi_mixed()
leopard2_case_multi_global_sanitizer()
leopard2_case_non_sanitizer_spelling()
leopard2_case_quoted_sanitizers()
leopard2_case_empty_build_type()
leopard2_case_default_sentinel_name_collision()
leopard2_case_false_named_custom_configuration()
