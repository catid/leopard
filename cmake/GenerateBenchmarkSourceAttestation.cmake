if(NOT DEFINED LEO2_SOURCE_DIR OR LEO2_SOURCE_DIR STREQUAL "")
    message(FATAL_ERROR "LEO2_SOURCE_DIR is required")
endif()
if(NOT DEFINED LEO2_OUTPUT_FILE OR LEO2_OUTPUT_FILE STREQUAL "")
    message(FATAL_ERROR "LEO2_OUTPUT_FILE is required")
endif()

get_filename_component(LEO2_SOURCE_DIR
    "${LEO2_SOURCE_DIR}" REALPATH)
get_filename_component(LEO2_OUTPUT_DIRECTORY
    "${LEO2_OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${LEO2_OUTPUT_DIRECTORY}")

# Separate build-tool processes can execute the always-run refresh target at
# the same time.  Serialize capture and publication for this one output so a
# slower capture cannot overwrite a newer one after racing through the
# compare/rename sequence.
file(LOCK "${LEO2_OUTPUT_FILE}.lock"
    GUARD PROCESS
    TIMEOUT 120
    RESULT_VARIABLE LEO2_LOCK_RESULT)
if(NOT LEO2_LOCK_RESULT STREQUAL "0")
    message(FATAL_ERROR
        "Unable to lock Leopard2 benchmark source attestation: "
        "${LEO2_LOCK_RESULT}")
endif()

if(NOT DEFINED LEO2_GIT_EXECUTABLE OR
   LEO2_GIT_EXECUTABLE STREQUAL "" OR
   NOT EXISTS "${LEO2_GIT_EXECUTABLE}")
    find_program(LEO2_GIT_EXECUTABLE NAMES git)
endif()

set(LEO2_SOURCE_IS_GIT_WORKTREE FALSE)
if(LEO2_GIT_EXECUTABLE)
    execute_process(
        COMMAND "${LEO2_GIT_EXECUTABLE}" -C "${LEO2_SOURCE_DIR}"
            rev-parse --show-toplevel
        RESULT_VARIABLE LEO2_WORKTREE_RESULT
        OUTPUT_VARIABLE LEO2_WORKTREE_OUTPUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(LEO2_WORKTREE_RESULT EQUAL 0 AND
       NOT LEO2_WORKTREE_OUTPUT STREQUAL "")
        # rev-parse walks parent directories, so merely succeeding would make
        # an unpacked source archive inherit an unrelated enclosing repository.
        # Canonicalizing both sides still admits linked worktrees, whose
        # top-level is the linked checkout root even though .git is a file.
        get_filename_component(LEO2_WORKTREE_TOPLEVEL
            "${LEO2_WORKTREE_OUTPUT}" REALPATH)
        if(LEO2_WORKTREE_TOPLEVEL STREQUAL LEO2_SOURCE_DIR)
            set(LEO2_SOURCE_IS_GIT_WORKTREE TRUE)
        endif()
    endif()
endif()

function(leo2_capture_git_state prefix)
    execute_process(
        COMMAND "${LEO2_GIT_EXECUTABLE}" -C "${LEO2_SOURCE_DIR}"
            rev-parse --verify HEAD
        RESULT_VARIABLE commit_result
        OUTPUT_VARIABLE commit
        ERROR_VARIABLE commit_error
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND "${LEO2_GIT_EXECUTABLE}" -C "${LEO2_SOURCE_DIR}"
            rev-parse --verify HEAD^{tree}
        RESULT_VARIABLE tree_result
        OUTPUT_VARIABLE tree
        ERROR_VARIABLE tree_error
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Preserve the established schema-v5 contract: despite the historical
    # field name, relevant untracked files also make the source dirty.
    execute_process(
        COMMAND "${LEO2_GIT_EXECUTABLE}" -C "${LEO2_SOURCE_DIR}"
            status --porcelain=v1 --untracked-files=normal
            --ignore-submodules=none
        RESULT_VARIABLE status_result
        OUTPUT_VARIABLE status
        ERROR_VARIABLE status_error
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT commit_result EQUAL 0 OR
       NOT tree_result EQUAL 0 OR
       NOT status_result EQUAL 0)
        message(FATAL_ERROR
            "Unable to capture a complete Leopard2 source identity.\n"
            "commit: ${commit_error}\n"
            "tree: ${tree_error}\n"
            "status: ${status_error}")
    endif()
    if(NOT commit MATCHES "^[0-9a-f]+$" OR
       NOT tree MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR
            "Git returned a non-canonical Leopard2 source identity")
    endif()

    string(SHA256 status_digest "${status}")
    if(status STREQUAL "")
        set(dirty 0)
    else()
        set(dirty 1)
    endif()

    set(${prefix}_COMMIT "${commit}" PARENT_SCOPE)
    set(${prefix}_TREE "${tree}" PARENT_SCOPE)
    set(${prefix}_STATUS_DIGEST "${status_digest}" PARENT_SCOPE)
    set(${prefix}_DIRTY "${dirty}" PARENT_SCOPE)
endfunction()

if(LEO2_SOURCE_IS_GIT_WORKTREE)
    # Capture twice around generation.  This cannot make an actively edited
    # source tree atomic, but it prevents a clean attestation from being
    # published when HEAD, the index, or worktree status changed while the
    # identity was sampled.
    leo2_capture_git_state(LEO2_BEFORE)
    leo2_capture_git_state(LEO2_AFTER)
    if(NOT LEO2_BEFORE_COMMIT STREQUAL LEO2_AFTER_COMMIT OR
       NOT LEO2_BEFORE_TREE STREQUAL LEO2_AFTER_TREE OR
       NOT LEO2_BEFORE_STATUS_DIGEST STREQUAL LEO2_AFTER_STATUS_DIGEST OR
       NOT LEO2_BEFORE_DIRTY STREQUAL LEO2_AFTER_DIRTY)
        message(FATAL_ERROR
            "Leopard2 source identity changed while it was captured; "
            "rerun the build")
    endif()
    set(LEO2_SOURCE_COMMIT "${LEO2_AFTER_COMMIT}")
    set(LEO2_SOURCE_TREE "${LEO2_AFTER_TREE}")
    set(LEO2_SOURCE_DIRTY "${LEO2_AFTER_DIRTY}")
else()
    # Source archives remain buildable, but cannot claim a committed identity.
    set(LEO2_SOURCE_COMMIT "unknown")
    set(LEO2_SOURCE_TREE "unknown")
    set(LEO2_SOURCE_DIRTY 1)
endif()

set(LEO2_HEADER_CONTENT
"#ifndef LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H
#define LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H

#undef LEO2_BENCHMARK_SOURCE_COMMIT
#undef LEO2_BENCHMARK_SOURCE_TREE
#undef LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY
#define LEO2_BENCHMARK_SOURCE_COMMIT \"${LEO2_SOURCE_COMMIT}\"
#define LEO2_BENCHMARK_SOURCE_TREE \"${LEO2_SOURCE_TREE}\"
#define LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY ${LEO2_SOURCE_DIRTY}

#endif
")

# Publish through a same-directory temporary.  Comparing content first
# preserves the header timestamp for an unchanged identity; file(RENAME) then
# replaces a changed header atomically on the build filesystem.
string(RANDOM LENGTH 16 ALPHABET 0123456789abcdef LEO2_TEMP_SUFFIX)
set(LEO2_TEMP_FILE
    "${LEO2_OUTPUT_FILE}.${LEO2_TEMP_SUFFIX}.tmp")
file(WRITE "${LEO2_TEMP_FILE}" "${LEO2_HEADER_CONTENT}")
file(SHA256 "${LEO2_TEMP_FILE}" LEO2_NEW_HEADER_SHA256)
set(LEO2_HEADER_CHANGED TRUE)
if(EXISTS "${LEO2_OUTPUT_FILE}")
    file(SHA256 "${LEO2_OUTPUT_FILE}" LEO2_OLD_HEADER_SHA256)
    if(LEO2_NEW_HEADER_SHA256 STREQUAL LEO2_OLD_HEADER_SHA256)
        set(LEO2_HEADER_CHANGED FALSE)
    endif()
endif()
if(LEO2_HEADER_CHANGED)
    file(RENAME "${LEO2_TEMP_FILE}" "${LEO2_OUTPUT_FILE}")
else()
    file(REMOVE "${LEO2_TEMP_FILE}")
endif()
