# Classify each build configuration independently for the portable-ISA CTest
# policy.  A multi-configuration tree may legitimately contain a sanitized
# Debug/ASan archive and an unsanitized Release archive at the same time; one
# configure-time boolean would either scan sanitizer-generated instructions or
# weaken the clean Release audit.

set(LEO2_SANITIZER_FLAG_PATTERN
    "(^|[ \t;\"'])(-fsanitize=[^ \t;\"']+|/fsanitize(=[^ \t;\"']+)?)")
set(LEO2_SANITIZED_CONFIGURATIONS)
set(LEO2_UNSANITIZED_CONFIGURATIONS)
set(LEO2_MULTI_CONFIG OFF)
set(LEO2_AUDIT_USES_DEFAULT_FLAGS OFF)

if(CMAKE_GENERATOR MATCHES "Visual Studio|Xcode|Multi-Config")
    set(LEO2_MULTI_CONFIG ON)
    set(LEO2_AUDIT_CONFIGURATIONS ${CMAKE_CONFIGURATION_TYPES})
else()
    if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
        set(LEO2_AUDIT_CONFIGURATIONS "${CMAKE_BUILD_TYPE}")
    else()
        # An empty CMAKE_BUILD_TYPE is a valid single-config CMake model.  Keep
        # one explicit list element so the general flags are still classified.
        set(LEO2_AUDIT_CONFIGURATIONS "__LEO2_DEFAULT__")
        set(LEO2_AUDIT_USES_DEFAULT_FLAGS ON)
    endif()
endif()

list(LENGTH LEO2_AUDIT_CONFIGURATIONS
    LEO2_AUDIT_CONFIGURATION_COUNT)
if(LEO2_MULTI_CONFIG AND LEO2_AUDIT_CONFIGURATION_COUNT EQUAL 0)
    message(FATAL_ERROR
        "Leopard2 sanitizer classification has no active build configuration")
endif()

foreach(LEO2_AUDIT_CONFIGURATION IN LISTS LEO2_AUDIT_CONFIGURATIONS)
    if(LEO2_AUDIT_USES_DEFAULT_FLAGS)
        set(LEO2_CONFIGURATION_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        string(TOUPPER "${LEO2_AUDIT_CONFIGURATION}"
            LEO2_AUDIT_CONFIGURATION_UPPER)
        set(LEO2_CONFIGURATION_FLAGS_VARIABLE
            "CMAKE_CXX_FLAGS_${LEO2_AUDIT_CONFIGURATION_UPPER}")
        set(LEO2_CONFIGURATION_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} ${${LEO2_CONFIGURATION_FLAGS_VARIABLE}}")
    endif()
    if("${LEO2_CONFIGURATION_CXX_FLAGS}" MATCHES
       "${LEO2_SANITIZER_FLAG_PATTERN}")
        list(APPEND LEO2_SANITIZED_CONFIGURATIONS
            "${LEO2_AUDIT_CONFIGURATION}")
    else()
        list(APPEND LEO2_UNSANITIZED_CONFIGURATIONS
            "${LEO2_AUDIT_CONFIGURATION}")
    endif()
endforeach()

# Retain the historical scalar for the strict single-config release gate.
# Registration below uses the two exact configuration lists and therefore does
# not collapse a mixed multi-config tree to this conservative summary.
set(LEO2_PRODUCTION_ARCHIVE_SANITIZED OFF)
list(LENGTH LEO2_SANITIZED_CONFIGURATIONS
    LEO2_SANITIZED_CONFIGURATION_COUNT)
list(LENGTH LEO2_UNSANITIZED_CONFIGURATIONS
    LEO2_UNSANITIZED_CONFIGURATION_COUNT)
if(LEO2_SANITIZED_CONFIGURATION_COUNT GREATER 0 AND
   LEO2_UNSANITIZED_CONFIGURATION_COUNT EQUAL 0)
    set(LEO2_PRODUCTION_ARCHIVE_SANITIZED ON)
endif()
