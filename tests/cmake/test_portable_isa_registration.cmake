if(NOT DEFINED LEO2_CTEST_COMMAND OR NOT DEFINED LEO2_BINARY_DIR OR
   NOT DEFINED LEO2_MULTI_CONFIG OR
   NOT DEFINED LEO2_SANITIZED_CONFIGURATIONS OR
   NOT DEFINED LEO2_UNSANITIZED_CONFIGURATIONS)
    message(FATAL_ERROR
        "Portable-ISA registration test is missing required inputs")
endif()

string(REPLACE "," ";" sanitized_configurations
    "${LEO2_SANITIZED_CONFIGURATIONS}")
string(REPLACE "," ";" unsanitized_configurations
    "${LEO2_UNSANITIZED_CONFIGURATIONS}")

function(leopard2_check_portable_inventory configuration expected_mode)
    set(inventory_command
        "${LEO2_CTEST_COMMAND}" --test-dir "${LEO2_BINARY_DIR}" -N)
    if(LEO2_MULTI_CONFIG)
        list(APPEND inventory_command -C "${configuration}")
    endif()
    execute_process(
        COMMAND ${inventory_command}
        RESULT_VARIABLE inventory_result
        OUTPUT_VARIABLE inventory_stdout
        ERROR_VARIABLE inventory_stderr)
    if(NOT inventory_result EQUAL 0)
        message(FATAL_ERROR
            "Could not inspect CTest inventory for '${configuration}' "
            "(${inventory_result})\n${inventory_stdout}\n${inventory_stderr}")
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

    if(expected_mode STREQUAL "archive")
        if(NOT archive_count EQUAL 1 OR NOT control_count EQUAL 0)
            message(FATAL_ERROR
                "Unsanitized configuration '${configuration}' did not register "
                "exactly one archive audit and no self-test-only substitute:\n"
                "${inventory_stdout}")
        endif()
    elseif(expected_mode STREQUAL "controls")
        if(NOT archive_count EQUAL 0 OR NOT control_count EQUAL 1)
            message(FATAL_ERROR
                "Sanitized configuration '${configuration}' did not replace "
                "only the archive audit with its checker self-test:\n"
                "${inventory_stdout}")
        endif()
    else()
        message(FATAL_ERROR
            "Unknown portable-ISA registration mode '${expected_mode}'")
    endif()
    if(NOT ordinary_count EQUAL 1)
        message(FATAL_ERROR
            "Portable-ISA policy unexpectedly removed ordinary CTests for "
            "'${configuration}':\n${inventory_stdout}")
    endif()
endfunction()

if(LEO2_MULTI_CONFIG)
    list(LENGTH sanitized_configurations sanitized_count)
    list(LENGTH unsanitized_configurations unsanitized_count)
    math(EXPR classified_count "${sanitized_count} + ${unsanitized_count}")
    if(classified_count EQUAL 0)
        message(FATAL_ERROR
            "Multi-config portable-ISA registration has no configurations")
    endif()
    foreach(configuration IN LISTS sanitized_configurations)
        list(FIND unsanitized_configurations "${configuration}"
            conflicting_configuration_index)
        if(NOT conflicting_configuration_index EQUAL -1)
            message(FATAL_ERROR
                "Configuration '${configuration}' has conflicting sanitizer "
                "classifications")
        endif()
        leopard2_check_portable_inventory("${configuration}" controls)
    endforeach()
    foreach(configuration IN LISTS unsanitized_configurations)
        leopard2_check_portable_inventory("${configuration}" archive)
    endforeach()
else()
    list(LENGTH sanitized_configurations sanitized_count)
    list(LENGTH unsanitized_configurations unsanitized_count)
    math(EXPR classified_count "${sanitized_count} + ${unsanitized_count}")
    if(NOT classified_count EQUAL 1)
        message(FATAL_ERROR
            "Single-config portable-ISA registration must classify exactly one "
            "configuration")
    endif()
    if(sanitized_count EQUAL 1)
        list(GET sanitized_configurations 0 configuration)
        leopard2_check_portable_inventory("${configuration}" controls)
    else()
        list(GET unsanitized_configurations 0 configuration)
        leopard2_check_portable_inventory("${configuration}" archive)
    endif()
endif()
