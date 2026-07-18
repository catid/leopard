if(NOT DEFINED LEO2_CTEST_COMMAND OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_EXPECTED_TEST)
    message(FATAL_ERROR
        "Portable-ISA registration test is missing required inputs")
endif()
if(NOT LEO2_EXPECTED_TEST STREQUAL "archive" AND
   NOT LEO2_EXPECTED_TEST STREQUAL "controls")
    message(FATAL_ERROR
        "Unknown portable-ISA registration mode '${LEO2_EXPECTED_TEST}'")
endif()

execute_process(
    COMMAND "${LEO2_CTEST_COMMAND}" --test-dir "${LEO2_BINARY_DIR}" -N
    RESULT_VARIABLE inventory_result
    OUTPUT_VARIABLE inventory_stdout
    ERROR_VARIABLE inventory_stderr)
if(NOT inventory_result EQUAL 0)
    message(FATAL_ERROR
        "Could not inspect CTest inventory (${inventory_result})\n"
        "${inventory_stdout}\n${inventory_stderr}")
endif()

string(REGEX MATCHALL
    "Test +#[0-9]+: leopard2_portable_isa(\r?\n|$)"
    archive_tests "${inventory_stdout}")
string(REGEX MATCHALL
    "Test +#[0-9]+: leopard2_portable_isa_checker_self_test(\r?\n|$)"
    control_tests "${inventory_stdout}")
string(REGEX MATCHALL
    "Test +#[0-9]+: leopard2_field_options(\r?\n|$)"
    ordinary_tests "${inventory_stdout}")
list(LENGTH archive_tests archive_count)
list(LENGTH control_tests control_count)
list(LENGTH ordinary_tests ordinary_count)

if(LEO2_EXPECTED_TEST STREQUAL "archive")
    if(NOT archive_count EQUAL 1 OR NOT control_count EQUAL 0)
        message(FATAL_ERROR
            "Clean build did not register exactly one archive audit and no "
            "self-test-only substitute:\n${inventory_stdout}")
    endif()
else()
    if(NOT archive_count EQUAL 0 OR NOT control_count EQUAL 1)
        message(FATAL_ERROR
            "Sanitizer build did not replace only the archive audit with its "
            "checker self-test:\n${inventory_stdout}")
    endif()
endif()
if(NOT ordinary_count EQUAL 1)
    message(FATAL_ERROR
        "Portable-ISA policy unexpectedly removed ordinary CTests:\n"
        "${inventory_stdout}")
endif()
