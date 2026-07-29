if(NOT DEFINED LEO2_SOURCE_DIR OR NOT DEFINED LEO2_BINARY_DIR)
    message(FATAL_ERROR
        "Target-architecture model test is missing required inputs")
endif()

file(READ "${LEO2_SOURCE_DIR}/CMakeLists.txt" cmake_source)
set(begin_marker "# BEGIN LEO2_TARGET_ARCHITECTURE_MODEL")
set(end_marker "# END LEO2_TARGET_ARCHITECTURE_MODEL")
string(FIND "${cmake_source}" "${begin_marker}" model_begin)
string(FIND "${cmake_source}" "${end_marker}" model_end)
if(model_begin EQUAL -1 OR model_end EQUAL -1 OR
   NOT model_begin LESS model_end)
    message(FATAL_ERROR
        "Target-architecture model markers are missing or out of order")
endif()
string(LENGTH "${end_marker}" end_marker_length)
math(EXPR model_length
    "${model_end} - ${model_begin} + ${end_marker_length}")
string(SUBSTRING "${cmake_source}" ${model_begin} ${model_length}
    target_architecture_model)

file(MAKE_DIRECTORY "${LEO2_BINARY_DIR}")

function(leopard2_check_target_architecture
        case_name system_processor apple osx_architectures
        expected_processor expected_x86 expected_x86_64 expected_universal)
    set(case_script "${LEO2_BINARY_DIR}/${case_name}.cmake")
    file(WRITE "${case_script}"
        "set(CMAKE_SYSTEM_PROCESSOR [==[${system_processor}]==])\n"
        "set(APPLE ${apple})\n"
        "set(CMAKE_OSX_ARCHITECTURES [==[${osx_architectures}]==])\n"
        "${target_architecture_model}\n"
        "if(NOT \"\${LEO2_TARGET_PROCESSOR}\" STREQUAL "
            "[==[${expected_processor}]==])\n"
        "  message(FATAL_ERROR \"${case_name}: target processor differs: "
            "'\${LEO2_TARGET_PROCESSOR}'\")\n"
        "endif()\n"
        "if(NOT \"\${LEO2_X86_TARGET}\" STREQUAL \"${expected_x86}\")\n"
        "  message(FATAL_ERROR \"${case_name}: x86 classification differs\")\n"
        "endif()\n"
        "if(NOT \"\${LEO2_X86_64_TARGET}\" STREQUAL "
            "\"${expected_x86_64}\")\n"
        "  message(FATAL_ERROR \"${case_name}: x86-64 classification differs\")\n"
        "endif()\n"
        "if(NOT \"\${LEO2_APPLE_UNIVERSAL_TARGET}\" STREQUAL "
            "\"${expected_universal}\")\n"
        "  message(FATAL_ERROR \"${case_name}: universal classification differs\")\n"
        "endif()\n")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -P "${case_script}"
        RESULT_VARIABLE case_result
        OUTPUT_VARIABLE case_stdout
        ERROR_VARIABLE case_stderr)
    if(NOT case_result EQUAL 0)
        message(FATAL_ERROR
            "Target-architecture case '${case_name}' failed (${case_result})\n"
            "${case_stdout}\n${case_stderr}")
    endif()
endfunction()

leopard2_check_target_architecture(
    nonapple-x86 AMD64 OFF "" AMD64 ON ON OFF)
leopard2_check_target_architecture(
    nonapple-arm aarch64 OFF "" aarch64 OFF OFF OFF)
leopard2_check_target_architecture(
    apple-cross-arm x86_64 ON arm64 arm64 OFF OFF OFF)
leopard2_check_target_architecture(
    apple-cross-x86 arm64 ON x86_64 x86_64 ON ON OFF)
leopard2_check_target_architecture(
    apple-x86-64h arm64 ON x86_64h x86_64h ON ON OFF)
leopard2_check_target_architecture(
    apple-universal x86_64 ON "arm64;x86_64" "" OFF OFF ON)
