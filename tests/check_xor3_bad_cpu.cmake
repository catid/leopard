if(NOT DEFINED PROGRAM OR NOT DEFINED CASE)
    message(FATAL_ERROR "PROGRAM and CASE are required")
endif()

if(CASE STREQUAL "empty")
    set(BAD_CPU "")
elseif(CASE STREQUAL "whitespace")
    set(BAD_CPU " ")
elseif(CASE STREQUAL "overflow")
    set(BAD_CPU "65536")
elseif(CASE STREQUAL "erange")
    set(BAD_CPU "9999999999999999999999999999999999999999")
else()
    message(FATAL_ERROR "Unknown bad-CPU case: ${CASE}")
endif()

execute_process(
    COMMAND "${PROGRAM}" --verify-only --cpu "${BAD_CPU}"
    RESULT_VARIABLE RESULT
    OUTPUT_VARIABLE STDOUT
    ERROR_VARIABLE STDERR
    TIMEOUT 10)

if(NOT "${RESULT}" STREQUAL "2")
    message(FATAL_ERROR
        "Bad --cpu case '${CASE}' returned ${RESULT}, expected 2\n"
        "stdout:\n${STDOUT}\n"
        "stderr:\n${STDERR}")
endif()
