if(NOT DEFINED LEO2_API_EXECUTABLE)
    message(FATAL_ERROR "API group contract requires the test executable")
endif()
execute_process(COMMAND "${LEO2_API_EXECUTABLE}" --list
    RESULT_VARIABLE result OUTPUT_VARIABLE groups ERROR_VARIABLE errors)
set(expected "dispatch\ncompat\ndecode\nrepair\nexpanded\nlarge-high\nlarge-low\n")
string(REPLACE "\r\n" "\n" groups "${groups}")
if(NOT result EQUAL 0 OR NOT groups STREQUAL expected OR NOT errors STREQUAL "")
    message(FATAL_ERROR "API group inventory mismatch: ${result}: ${groups} ${errors}")
endif()
foreach(bad_argument IN ITEMS unknown --lis "")
    execute_process(COMMAND "${LEO2_API_EXECUTABLE}" "${bad_argument}"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors)
    if(NOT result STREQUAL "1" OR NOT output STREQUAL "" OR
       NOT errors MATCHES "unknown API group")
        message(FATAL_ERROR "API selector accepted invalid argument: ${bad_argument}")
    endif()
endforeach()
execute_process(COMMAND "${LEO2_API_EXECUTABLE}" --list extra
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors)
if(NOT result STREQUAL "1" OR NOT output STREQUAL "" OR
   NOT errors MATCHES "expected zero arguments or one API group")
    message(FATAL_ERROR "API selector ignored extra arguments")
endif()
