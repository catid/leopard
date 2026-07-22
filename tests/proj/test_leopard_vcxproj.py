#!/usr/bin/env python3
"""Fail-closed structural contract for the Visual Studio library project.

The test runs on any host with Python.  It cannot replace a native MSVC build,
but it makes the hand-maintained project fail when the production CMake graph,
MSBuild configurations, ISA isolation, or optional-CUDA contract drifts.
"""

from collections import Counter
import copy
from pathlib import Path, PurePosixPath
import re
import sys
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
SOLUTION = ROOT / "proj" / "Leopard.sln"
PROJECT = ROOT / "proj" / "Leopard.vcxproj"
FILTERS = ROOT / "proj" / "Leopard.vcxproj.filters"
CMAKE = ROOT / "CMakeLists.txt"
NS = {"msb": "http://schemas.microsoft.com/developer/msbuild/2003"}
EXPECTED_TOOLS_VERSION = "14.0"
LEGACY_PROJECTS = (
    PROJECT,
    ROOT / "tests" / "proj" / "Benchmark.vcxproj",
    ROOT / "tests" / "proj" / "Experiments.vcxproj",
)

COMPILE_SUFFIXES = (".c", ".cc", ".cpp", ".cxx")
HEADER_SUFFIXES = (".h", ".hh", ".hpp", ".hxx", ".inl")
CUDA_SUFFIXES = (".cu", ".cuh")
KNOWN_SOURCE_SUFFIXES = COMPILE_SUFFIXES + HEADER_SUFFIXES + CUDA_SUFFIXES
AVX2_SOURCE_FILES = {
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
}
EXPECTED_CONFIGS = {
    "debug|win32": ("Debug", "Win32"),
    "debug|x64": ("Debug", "x64"),
    "release|win32": ("Release", "Win32"),
    "release|x64": ("Release", "x64"),
}
BACKEND_DEFINITIONS = {
    "LEO2_DISABLE_SSSE3_CODEGEN=1",
    "LEO2_DISABLE_AVX2_CODEGEN=1",
    "LEO2_HAVE_SSSE3_BACKEND=1",
    "LEO2_HAVE_AVX2_BACKEND=1",
}
PROTECTED_MACRO_NAMES = {
    definition.split("=")[0] for definition in BACKEND_DEFINITIONS
}
MSBUILD_PROPERTY = re.compile(r"\$\([^)]+\)")
MSBUILD_METADATA = re.compile(r"%\([^)]+\)")
BACKEND_OPTION_OVERRIDE = re.compile(
    r"(?:/|-)[du]\s*leo2_(?:backend_force|disable|have)_",
    re.IGNORECASE)
WPO_OPTION = re.compile(
    r"/(?:gl|ltcg)(?=$|[\s:,+-])|-flto(?:=\S+)?|/qipo(?=$|[\s:,+-])",
    re.IGNORECASE)

# MSBuild item, property, and metadata names are case-insensitive even though
# XML QName searches are not.  Keep every security-relevant name in one
# canonical spelling so a case variant cannot hide from a later exact scan.
MSBUILD_CANONICAL_NAMES = {
    name.lower(): name for name in (
        "AdditionalDependencies", "AdditionalLibraryDirectories",
        "AdditionalOptions",
        "BufferSecurityCheck", "CharacterSet", "Choose", "ClCompile",
        "ClInclude", "CLToolExe", "CLToolPath", "Command", "Configuration",
        "ConfigurationType", "CudaCompile", "CustomBuild",
        "CustomBuildStep", "EnableCOMDATFolding",
        "EnableEnhancedInstructionSet", "ExcludedFromBuild",
        "ExecutablePath", "FavorSizeOrSpeed", "Filter",
        "ForcedIncludeFiles", "ForcedUsingFiles", "FunctionLevelLinking",
        "GenerateDebugInformation", "Import", "ImportGroup", "IntDir",
        "InlineFunctionExpansion", "IntrinsicFunctions",
        "ItemDefinitionGroup", "ItemGroup", "Lib", "Link",
        "LinkTimeCodeGeneration", "OmitFramePointers", "OpenMPSupport",
        "Optimization", "OptimizeReferences", "Otherwise", "OutDir",
        "Platform", "PlatformToolset", "PostBuildEvent", "PreBuildEvent",
        "PreLinkEvent", "PreprocessorDefinitions",
        "ProfileGuidedOptimization", "Project", "ProjectConfiguration",
        "ProjectGuid", "ProjectName", "ProjectReference", "PropertyGroup",
        "RootNamespace", "RuntimeLibrary", "SDLCheck", "Target", "ToolExe",
        "ToolPath", "UndefinePreprocessorDefinitions", "UseDebugLibraries",
        "UseEnv", "UseLinkTimeCodeGeneration", "UsingTask", "VCInstallDir",
        "VCTargetsPath", "VCToolsInstallDir", "VCToolsPath", "WarningLevel",
        "When", "WholeProgramOptimization",
    )
}


class ContractError(AssertionError):
    pass


_CMAKE_BRACKET_OPEN = re.compile(r"\[(=*)\[")


def cmake_bracket_span(text, index, label):
    """Return (content, next index) for a CMake [=[...]=] construct."""
    match = _CMAKE_BRACKET_OPEN.match(text, index)
    if not match:
        return None
    closing = "]" + match.group(1) + "]"
    end = text.find(closing, match.end())
    if end < 0:
        raise ContractError("unterminated CMake bracket " + label)
    content = text[match.end():end]
    if content.startswith("\r\n"):
        content = content[2:]
    elif content.startswith("\n"):
        content = content[1:]
    return content, end + len(closing)


def cmake_comment_end(text, index):
    """Return the first index after a line or bracket comment at '#'."""
    bracket = cmake_bracket_span(text, index + 1, "comment")
    if bracket is not None:
        return bracket[1]
    newline = text.find("\n", index)
    return len(text) if newline < 0 else newline + 1


def normalized_msbuild_option_text(text):
    return (text or "").replace('"', "").replace("'", "")


def msbuild_values(text):
    return {
        value.strip().strip('"').strip("'")
        for value in (text or "").split(";")
        if value.strip()
    }


def reject_msbuild_expansion(text, label, allowed_metadata=()):
    if MSBUILD_PROPERTY.search(text or ""):
        raise ContractError(label + " contains unresolved MSBuild property expansion")
    allowed = set(allowed_metadata)
    for expansion in MSBUILD_METADATA.findall(text or ""):
        if expansion not in allowed:
            raise ContractError(label + " contains unresolved MSBuild metadata: " +
                                expansion)


def reject_backend_option_override(text, label):
    normalized = normalized_msbuild_option_text(text)
    if BACKEND_OPTION_OVERRIDE.search(normalized):
        raise ContractError(label + " overrides backend isolation macros")


def reject_wpo_options(text, label):
    normalized = normalized_msbuild_option_text(text)
    if WPO_OPTION.search(normalized):
        raise ContractError(label + " enables whole-program optimization")


def cmake_commands(text):
    """Return (lower-case command, argument text) with balanced parentheses."""
    commands = []
    index = 0
    length = len(text)
    while index < length:
        if text[index].isspace():
            index += 1
            continue
        if text[index] == "#":
            index = cmake_comment_end(text, index)
            continue
        if cmake_bracket_span(text, index, "argument") is not None:
            raise ContractError("stray CMake bracket argument")
        match = re.match(r"[A-Za-z_][A-Za-z0-9_]*", text[index:])
        if not match:
            index += 1
            continue
        name = match.group(0).lower()
        cursor = index + len(match.group(0))
        while cursor < length:
            if text[cursor].isspace():
                cursor += 1
                continue
            if text[cursor] == "#":
                cursor = cmake_comment_end(text, cursor)
                continue
            break
        if cursor >= length or text[cursor] != "(":
            index = cursor
            continue
        start = cursor + 1
        cursor = start
        depth = 1
        quoted = False
        escaped = False
        while cursor < length and depth:
            char = text[cursor]
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quoted = not quoted
            elif not quoted and char == "#":
                cursor = cmake_comment_end(text, cursor)
                continue
            elif not quoted:
                bracket = (cmake_bracket_span(
                    text, cursor, "argument") if char == "[" else None)
                if bracket is not None:
                    cursor = bracket[1]
                    continue
                if char == "(":
                    depth += 1
                elif char == ")":
                    depth -= 1
                    if depth == 0:
                        commands.append((name, text[start:cursor]))
                        cursor += 1
                        break
            cursor += 1
        if depth:
            raise ContractError("unterminated CMake command: " + name)
        index = cursor
    return commands


def cmake_tokens(body):
    tokens = []
    token = []
    token_started = False
    quoted = False
    escaped = False

    def flush():
        nonlocal token_started
        if token_started:
            tokens.append("".join(token))
            del token[:]
            token_started = False

    index = 0
    while index < len(body):
        char = body[index]
        if escaped:
            token.append(char)
            token_started = True
            escaped = False
        elif char == "\\":
            token_started = True
            escaped = True
        elif char == '"':
            token_started = True
            quoted = not quoted
        elif not quoted and char == "#":
            flush()
            index = cmake_comment_end(body, index)
            continue
        elif not quoted and char == "[":
            bracket = cmake_bracket_span(body, index, "argument")
            if bracket is not None:
                raise ContractError(
                    "CMake bracket arguments are unsupported by the "
                    "production graph proof")
            token.append(char)
            token_started = True
        elif not quoted and char.isspace():
            flush()
        elif not quoted and char == ";":
            if token_started:
                flush()
            else:
                tokens.append("")
        else:
            token.append(char)
            token_started = True
        index += 1
    if quoted or escaped:
        raise ContractError("unterminated quote or escape in CMake arguments")
    flush()
    return tokens


def cmake_condition_tokens(body):
    """Tokenize an if() body while retaining quotes and grouping."""
    tokens = []
    token = []
    token_started = False
    quoted = False
    token_was_quoted = False
    escaped = False

    def flush():
        nonlocal token_started
        if token_started:
            tokens.append(("".join(token), token_was_quoted))
            del token[:]
            token_started = False

    index = 0
    while index < len(body):
        char = body[index]
        if escaped:
            token.append(char)
            token_started = True
            escaped = False
        elif char == "\\":
            token_started = True
            escaped = True
        elif char == '"':
            token_started = True
            quoted = not quoted
            token_was_quoted = True
        elif not quoted and char == "#":
            flush()
            index = cmake_comment_end(body, index)
            continue
        elif not quoted and char == "[":
            bracket = cmake_bracket_span(body, index, "argument")
            if bracket is not None:
                raise ContractError(
                    "CMake bracket arguments are unsupported by the "
                    "production graph proof")
            token.append(char)
            token_started = True
        elif not quoted and char in "()":
            flush()
            tokens.append((char, False))
            token_was_quoted = False
        elif not quoted and char.isspace():
            flush()
            token_was_quoted = False
        elif not quoted and char == ";":
            if token_started:
                flush()
            else:
                tokens.append(("", False))
            token_was_quoted = False
        else:
            token.append(char)
            token_started = True
        index += 1
    if quoted or escaped:
        raise ContractError("unterminated quote or escape in CMake condition")
    flush()
    return tokens


# Boolean formula helpers used by the constrained CMake interpreter.  Formulae
# retain identities for compiler-probe results, so contradictions such as
# FLAG AND NOT FLAG cannot be mistaken for reachable configurations.
BOOL_TRUE = ("constant", True)
BOOL_FALSE = ("constant", False)
BOOL_SYMBOL_PREFIX = "__cmake_boolean_symbol__:"
CONDITIONAL_ASSIGNMENT_PREFIX = "conditional CMake source variable is ambiguous: "


def bool_atom(name):
    return ("atom", name)


def bool_not(formula):
    if formula == BOOL_TRUE:
        return BOOL_FALSE
    if formula == BOOL_FALSE:
        return BOOL_TRUE
    if formula[0] == "not":
        return formula[1]
    return ("not", formula)


def bool_and(*formulae):
    flattened = []
    for formula in formulae:
        if formula == BOOL_FALSE:
            return BOOL_FALSE
        if formula == BOOL_TRUE:
            continue
        flattened.extend(formula[1:] if formula[0] == "and" else (formula,))
    unique = []
    for formula in flattened:
        if formula in unique:
            continue
        if bool_not(formula) in unique:
            return BOOL_FALSE
        unique.append(formula)
    if not unique:
        return BOOL_TRUE
    if len(unique) == 1:
        return unique[0]
    return tuple(["and"] + sorted(unique, key=repr))


def bool_or(*formulae):
    flattened = []
    for formula in formulae:
        if formula == BOOL_TRUE:
            return BOOL_TRUE
        if formula == BOOL_FALSE:
            continue
        flattened.extend(formula[1:] if formula[0] == "or" else (formula,))
    unique = []
    for formula in flattened:
        if formula in unique:
            continue
        if bool_not(formula) in unique:
            return BOOL_TRUE
        unique.append(formula)
    if not unique:
        return BOOL_FALSE
    if len(unique) == 1:
        return unique[0]
    return tuple(["or"] + sorted(unique, key=repr))


def bool_formula_atoms(formula):
    if formula[0] == "atom":
        return {formula[1]}
    if formula[0] == "constant":
        return set()
    atoms = set()
    for child in formula[1:]:
        atoms.update(bool_formula_atoms(child))
    return atoms


def bool_formula_value(formula, assignment):
    operation = formula[0]
    if operation == "constant":
        return formula[1]
    if operation == "atom":
        return assignment[formula[1]]
    if operation == "not":
        return not bool_formula_value(formula[1], assignment)
    values = [bool_formula_value(child, assignment)
              for child in formula[1:]]
    return all(values) if operation == "and" else any(values)


def bool_satisfiable(formula):
    if formula == BOOL_TRUE:
        return True
    if formula == BOOL_FALSE:
        return False
    atoms = sorted(bool_formula_atoms(formula))
    if len(atoms) > 16:
        raise ContractError(
            "CMake condition exceeds fail-closed symbolic proof limit")
    for mask in range(1 << len(atoms)):
        assignment = {
            atom: bool(mask & (1 << index))
            for index, atom in enumerate(atoms)
        }
        if bool_formula_value(formula, assignment):
            return True
    return False


def bool_tautology(formula):
    return not bool_satisfiable(bool_not(formula))


class CMakeProductionGraph(object):
    """Constrained, fail-closed evaluator for the MSVC production graph.

    The checked-in project represents Windows x86/x64 MSVC builds.  This
    interpreter therefore fixes only authoritative CMake built-ins for that
    configuration and retains compiler-probe results as symbolic booleans.
    Variables and branch decisions are evaluated in command order.
    """

    _variable = re.compile(r"^\$\{([A-Za-z_][A-Za-z0-9_]*)\}$")
    _embedded_variable = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
    _target_objects = re.compile(
        r"^\$<TARGET_OBJECTS:([A-Za-z_][A-Za-z0-9_.+\-]*)>$")
    _cuda_graph_marker = re.compile(
        r"cuda|nvcc|nvrtc|nvidia|cudart|ptxas|nvlink|fatbinary|"
        r"cuobjdump|nvc\+\+",
        re.IGNORECASE)

    _target_build_mutation_commands = {
        "target_compile_definitions", "target_compile_features",
        "target_compile_options", "target_include_directories",
        "target_link_directories", "target_link_libraries",
        "target_link_interface_libraries", "target_link_options",
        "target_precompile_headers",
    }
    _directory_build_mutation_commands = {
        "add_compile_definitions", "add_compile_options", "add_definitions",
        "add_link_options", "include_directories", "link_directories",
        "link_libraries", "remove_definitions",
    }
    _build_extension_commands = {
        "add_custom_command", "add_custom_target", "add_dependencies",
        "cmake_parse_arguments", "configure_file",
        "configure_package_config_file", "execute_process", "file",
        "try_compile", "try_run",
    }
    _modeled_command_names = {
        "add_library", "add_subdirectory", "check_c_compiler_flag",
        "check_cxx_compiler_flag", "cmake_dependent_option",
        "cmake_language", "cmake_minimum_required", "cmake_parse_arguments",
        "cmake_path",
        "configure_package_config_file", "find_package",
        "get_cmake_property", "get_directory_property",
        "get_filename_component", "get_property",
        "get_source_file_property", "get_target_property", "include",
        "list", "load_command", "math", "option", "project",
        "separate_arguments", "set", "set_directory_properties",
        "set_property", "set_source_files_properties",
        "set_target_properties", "string", "subdirs", "target_sources",
        "unset", "enable_language",
    }
    _safe_root_command_names = {
        "add_executable", "add_test", "else", "elseif", "enable_testing",
        "endif", "find_program", "if", "install", "message",
        "set_tests_properties",
        "function", "endfunction", "macro", "endmacro", "foreach",
        "endforeach", "while", "endwhile", "block", "endblock",
    }
    _protected_variable = re.compile(
        r"^(?:CMAKE_.+|WIN32|MSVC|LEOPARD_INSTALL_CMAKEDIR|"
        r"LEOPARD_ENABLE_GF(?:8|16)|"
        r"ENABLE_OPENMP|LEO2_FLAG_ARCH_AVX2|"
        r"LEO2_(?:BACKEND_VARIANT(?:_NORMALIZED)?|BUILD_BENCHMARKS|BUILD_FUZZERS|"
        r"BUILD_TESTS|ENABLE_CUDA|PORTABLE_ISA_RELEASE_AUDIT|"
        r"EXPERIMENT_K1_DECODE_VALIDATOR_INLINE)|"
        r"CXX_FLAG_(?:O2|Oy|Zi|W4)|"
        r"(?:OpenMP|OPENMP|Threads|THREADS)_.+)$")
    _approved_production_mutations = {
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_DISABLE_SSSE3_CODEGEN=1",
            "LEO2_DISABLE_AVX2_CODEGEN=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_HAVE_SSSE3_BACKEND=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_HAVE_AVX2_BACKEND=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_SCALAR=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_SSSE3=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_AVX2=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_AVX512=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF8=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF16=1")),
        ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE=1")),
        ("leopard", "target_include_directories", (
            "PUBLIC", "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>",
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")),
        ("leopard", "target_link_libraries", (
            "PUBLIC", "Threads::Threads")),
        ("leopard", "target_link_libraries", (
            "PUBLIC", "OpenMP::OpenMP_CXX")),
        ("leopard", "target_link_libraries", (
            "PUBLIC", "${OpenMP_CXX_FLAGS}")),
        ("leopard2_backend_ssse3", "target_include_directories", (
            "PRIVATE", "${CMAKE_CURRENT_SOURCE_DIR}")),
        ("leopard2_backend_avx2", "target_include_directories", (
            "PRIVATE", "${CMAKE_CURRENT_SOURCE_DIR}")),
        ("leopard2_backend_avx2", "target_compile_options", (
            "PRIVATE", "/arch:AVX2")),
        ("leopard2_backend_ssse3", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF8=1")),
        ("leopard2_backend_ssse3", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF16=1")),
        ("leopard2_backend_avx2", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF8=1")),
        ("leopard2_backend_avx2", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF16=1")),
    }
    _approved_protected_assignments = {
        ("CMAKE_CONFIGURATION_TYPES", "Debug;Release"),
        ("CMAKE_CXX_STANDARD", "11"),
        ("CMAKE_BUILD_TYPE", "Release"),
        ("CMAKE_CXX_FLAGS_RELEASE", "${CMAKE_CXX_FLAGS_RELEASE} /O2"),
        ("CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Oy"),
        ("CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Zi"),
        ("CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} /W4"),
        ("CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}"),
        ("CMAKE_EXE_LINKER_FLAGS",
         "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}"),
        ("LEOPARD_INSTALL_CMAKEDIR",
         "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"),
        ("LEO2_BACKEND_VARIANT", "auto"),
        ("LEO2_BACKEND_VARIANT", "${LEO2_BACKEND_VARIANT_NORMALIZED}"),
    }
    _approved_protected_set_commands = {
        ("CMAKE_CONFIGURATION_TYPES", "Debug;Release", "CACHE", "STRING",
         "", "FORCE"),
        ("CMAKE_CXX_STANDARD", "11"),
        ("CMAKE_BUILD_TYPE", "Release"),
        ("CMAKE_CXX_FLAGS_RELEASE", "${CMAKE_CXX_FLAGS_RELEASE} /O2"),
        ("CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Oy"),
        ("CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Zi"),
        ("CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} /W4"),
        ("CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}"),
        ("CMAKE_EXE_LINKER_FLAGS",
         "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}"),
        ("LEOPARD_INSTALL_CMAKEDIR",
         "${CMAKE_INSTALL_LIBDIR}/cmake/leopard", "CACHE", "STRING",
         "Install directory for Leopard CMake package files"),
        ("LEO2_BACKEND_VARIANT", "auto", "CACHE", "STRING",
         "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512"),
        ("LEO2_BACKEND_VARIANT", "${LEO2_BACKEND_VARIANT_NORMALIZED}",
         "CACHE", "STRING",
         "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512",
         "FORCE"),
    }
    _approved_package_configure = (
        "cmake/leopardConfig.cmake.in",
        "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake",
        "INSTALL_DESTINATION", "${LEOPARD_INSTALL_CMAKEDIR}")
    _locator_git_find = (
        "LEO2_LOCATOR_GIT_EXECUTABLE", "NAMES", "git")
    _locator_git_revision = (
        "COMMAND", "${LEO2_LOCATOR_GIT_EXECUTABLE}", "rev-parse", "HEAD",
        "WORKING_DIRECTORY", "${CMAKE_CURRENT_SOURCE_DIR}",
        "RESULT_VARIABLE", "LEO2_LOCATOR_GIT_RESULT",
        "OUTPUT_VARIABLE", "LEO2_LOCATOR_GIT_OUTPUT",
        "OUTPUT_STRIP_TRAILING_WHITESPACE")
    _locator_git_tree = (
        "COMMAND", "${LEO2_LOCATOR_GIT_EXECUTABLE}", "rev-parse",
        "HEAD^{tree}",
        "WORKING_DIRECTORY", "${CMAKE_CURRENT_SOURCE_DIR}",
        "RESULT_VARIABLE", "LEO2_LOCATOR_TREE_RESULT",
        "OUTPUT_VARIABLE", "LEO2_LOCATOR_TREE_OUTPUT",
        "OUTPUT_STRIP_TRAILING_WHITESPACE")
    _locator_git_status = (
        "COMMAND", "${LEO2_LOCATOR_GIT_EXECUTABLE}", "status",
        "--porcelain", "--untracked-files=normal",
        "WORKING_DIRECTORY", "${CMAKE_CURRENT_SOURCE_DIR}",
        "RESULT_VARIABLE", "LEO2_LOCATOR_STATUS_RESULT",
        "OUTPUT_VARIABLE", "LEO2_LOCATOR_STATUS_OUTPUT",
        "OUTPUT_STRIP_TRAILING_WHITESPACE")
    _required_locator_provenance_commands = Counter({
        ("find_program", _locator_git_find): 1,
        ("execute_process", _locator_git_revision): 1,
        ("execute_process", _locator_git_tree): 1,
        ("execute_process", _locator_git_status): 1,
    })
    _sparse_sidecar_link_depends = (
        "TARGET", "bench_leopard2_sparse_encode", "APPEND", "PROPERTY",
        "LINK_DEPENDS",
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/"
        "WriteSparseEncodeEvidenceSidecar.cmake")
    _sparse_sidecar_post_build = (
        "TARGET", "bench_leopard2_sparse_encode", "POST_BUILD", "COMMAND",
        "${CMAKE_COMMAND}",
        "-DOUTPUT=$<TARGET_FILE:bench_leopard2_sparse_encode>."
        "leopard2-evidence",
        "-DEXECUTABLE=$<TARGET_FILE:bench_leopard2_sparse_encode>",
        "-DPRODUCTION_ARCHIVE=$<TARGET_FILE:leopard>",
        "-DBENCHMARK_OBJECT=$<TARGET_OBJECTS:"
        "leopard2_sparse_encode_benchmark_object>",
        "-DORACLE_OBJECT=$<TARGET_OBJECTS:"
        "leopard2_sparse_encode_oracle_object>",
        "-DBENCHMARK_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/"
        "bench_leopard2_sparse_encode.dir/link.txt",
        "-DPRODUCTION_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/"
        "leopard.dir/link.txt",
        "-DBUILD_PROGRAM=${CMAKE_MAKE_PROGRAM}",
        "-DBUILD_ROOT=${CMAKE_CURRENT_BINARY_DIR}",
        "-DBUILD_GENERATOR=${CMAKE_GENERATOR}",
        "-P",
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/"
        "WriteSparseEncodeEvidenceSidecar.cmake",
        "VERBATIM")
    _required_sparse_sidecar_commands = Counter({
        ("set_property", _sparse_sidecar_link_depends): 1,
        ("add_custom_command", _sparse_sidecar_post_build): 1,
    })
    _required_trusted_commands = Counter({
        ("cmake_minimum_required", ("VERSION", "3.7")): 1,
        ("project", ("leopard",)): 1,
        ("add_library", ("leopard", "STATIC", "${LIB_SOURCE_FILES}")): 1,
        ("add_library", ("libleopard", "ALIAS", "leopard")): 1,
        ("include", ("CMakeDependentOption",)): 1,
        ("include", ("CheckCXXCompilerFlag",)): 1,
        ("include", ("CMakePackageConfigHelpers",)): 1,
        ("include", ("GNUInstallDirs",)): 1,
        ("find_package", ("OpenMP",)): 1,
        ("find_package", ("Threads", "REQUIRED")): 1,
        ("option", (
            "LEO2_BUILD_TESTS", "Build Leopard2 correctness tests", "ON")): 1,
        ("option", (
            "LEO2_BUILD_BENCHMARKS", "Build Leopard benchmark programs",
            "ON")): 1,
        ("option", (
            "LEO2_BUILD_ALLK_DIAGNOSTIC",
            "Build the source-attested all-K Leopard1 comparison benchmark",
            "OFF")): 1,
        ("option", (
            "LEO2_BUILD_FUZZERS", "Build Leopard2 libFuzzer targets",
            "OFF")): 1,
        ("option", (
            "LEO2_ENABLE_CUDA", "Build the optional Leopard2 CUDA backend",
            "OFF")): 1,
        ("option", (
            "LEOPARD_ENABLE_GF8", "Include the GF(2^8) codec", "ON")): 1,
        ("option", (
            "LEOPARD_ENABLE_GF16", "Include the GF(2^16) codec", "ON")): 1,
        ("option", (
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT",
            "Require the strict x86-64 portable-ISA release audit", "OFF")): 1,
        ("option", (
            "LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE",
            "Diagnostic: force-inline the K=1 decode buffer validator",
            "OFF")): 1,
        ("string", (
            "TOLOWER", "${LEO2_BACKEND_VARIANT}",
            "LEO2_BACKEND_VARIANT_NORMALIZED")): 1,
        ("configure_package_config_file", _approved_package_configure): 1,
        ("check_cxx_compiler_flag", ("/O2", "CXX_FLAG_O2")): 1,
        ("check_cxx_compiler_flag", ("/Oy", "CXX_FLAG_Oy")): 1,
        ("check_cxx_compiler_flag", ("/Zi", "CXX_FLAG_Zi")): 1,
        ("check_cxx_compiler_flag", ("/W4", "CXX_FLAG_W4")): 1,
        ("check_cxx_compiler_flag",
         ("/arch:AVX2", "LEO2_FLAG_ARCH_AVX2")): 1,
        ("cmake_dependent_option", (
            "ENABLE_OPENMP", "Enable OpenMP support", "ON",
            "OPENMP_FOUND", "OFF")): 1,
        ("install", (
            "TARGETS", "leopard", "EXPORT", "leopardTargets", "ARCHIVE",
            "DESTINATION", "${CMAKE_INSTALL_LIBDIR}", "LIBRARY", "DESTINATION",
            "${CMAKE_INSTALL_LIBDIR}", "RUNTIME", "DESTINATION",
            "${CMAKE_INSTALL_BINDIR}")): 1,
        ("install", (
            "FILES", "leopard.h", "leopard2.h", "DESTINATION",
            "${CMAKE_INSTALL_INCLUDEDIR}")): 1,
        ("install", (
            "FILES", "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake",
            "DESTINATION", "${LEOPARD_INSTALL_CMAKEDIR}")): 1,
        ("install", (
            "EXPORT", "leopardTargets", "FILE", "leopardTargets.cmake",
            "NAMESPACE", "leopard::", "DESTINATION",
            "${LEOPARD_INSTALL_CMAKEDIR}")): 1,
    })
    # CMake is an imperative language: proving that each security-sensitive
    # command appears once is not enough when moving a package discovery,
    # compiler probe, or cache assignment changes what later commands see.
    # Keep one exact data-flow order for every trusted command, protected
    # assignment, and production-target mutation.  Membership and guard checks
    # below remain responsible for useful missing/duplicate diagnostics; this
    # sequence closes reorder-only attacks.
    _required_contract_event_order = (
        ("trusted", ("cmake_minimum_required", ("VERSION", "3.7"))),
        ("trusted", ("project", ("leopard",))),
        ("trusted", ("include", ("CMakeDependentOption",))),
        ("trusted", ("include", ("CheckCXXCompilerFlag",))),
        ("trusted", ("include", ("CMakePackageConfigHelpers",))),
        ("trusted", ("include", ("GNUInstallDirs",))),
        ("protected", (
            "CMAKE_CONFIGURATION_TYPES", "Debug;Release", "CACHE",
            "STRING", "", "FORCE")),
        ("protected", ("CMAKE_CXX_STANDARD", "11")),
        ("trusted", ("option", (
            "LEO2_BUILD_TESTS", "Build Leopard2 correctness tests", "ON"))),
        ("trusted", ("option", (
            "LEO2_BUILD_BENCHMARKS", "Build Leopard benchmark programs",
            "ON"))),
        ("trusted", ("option", (
            "LEO2_BUILD_ALLK_DIAGNOSTIC",
            "Build the source-attested all-K Leopard1 comparison benchmark",
            "OFF"))),
        ("trusted", ("option", (
            "LEO2_BUILD_FUZZERS", "Build Leopard2 libFuzzer targets",
            "OFF"))),
        ("trusted", ("option", (
            "LEO2_ENABLE_CUDA", "Build the optional Leopard2 CUDA backend",
            "OFF"))),
        ("trusted", ("option", (
            "LEOPARD_ENABLE_GF8", "Include the GF(2^8) codec", "ON"))),
        ("trusted", ("option", (
            "LEOPARD_ENABLE_GF16", "Include the GF(2^16) codec", "ON"))),
        ("trusted", ("option", (
            "LEO2_PORTABLE_ISA_RELEASE_AUDIT",
            "Require the strict x86-64 portable-ISA release audit", "OFF"))),
        ("trusted", ("option", (
            "LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE",
            "Diagnostic: force-inline the K=1 decode buffer validator",
            "OFF"))),
        ("protected", (
            "LEO2_BACKEND_VARIANT", "auto", "CACHE", "STRING",
            "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512")),
        ("trusted", ("string", (
            "TOLOWER", "${LEO2_BACKEND_VARIANT}",
            "LEO2_BACKEND_VARIANT_NORMALIZED"))),
        ("protected", (
            "LEO2_BACKEND_VARIANT", "${LEO2_BACKEND_VARIANT_NORMALIZED}",
            "CACHE", "STRING",
            "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512",
            "FORCE")),
        ("guarded", ("enable_language", ("CUDA",))),
        ("protected", ("CMAKE_BUILD_TYPE", "Release")),
        ("trusted", ("check_cxx_compiler_flag", ("/O2", "CXX_FLAG_O2"))),
        ("trusted", ("check_cxx_compiler_flag", ("/Oy", "CXX_FLAG_Oy"))),
        ("trusted", ("check_cxx_compiler_flag", ("/Zi", "CXX_FLAG_Zi"))),
        ("trusted", ("check_cxx_compiler_flag", ("/W4", "CXX_FLAG_W4"))),
        ("protected", (
            "CMAKE_CXX_FLAGS_RELEASE", "${CMAKE_CXX_FLAGS_RELEASE} /O2")),
        ("protected", (
            "CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Oy")),
        ("protected", (
            "CMAKE_CXX_FLAGS_DEBUG", "${CMAKE_CXX_FLAGS_DEBUG} /Zi")),
        ("protected", ("CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} /W4")),
        ("trusted", ("find_package", ("OpenMP",))),
        ("trusted", ("find_package", ("Threads", "REQUIRED"))),
        ("trusted", ("cmake_dependent_option", (
            "ENABLE_OPENMP", "Enable OpenMP support", "ON",
            "OPENMP_FOUND", "OFF"))),
        ("protected", (
            "CMAKE_CXX_FLAGS", "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")),
        ("protected", (
            "CMAKE_EXE_LINKER_FLAGS",
            "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")),
        ("trusted", ("add_library", (
            "leopard", "STATIC", "${LIB_SOURCE_FILES}"))),
        ("trusted", ("add_library", (
            "libleopard", "ALIAS", "leopard"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF8=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "NO_LEO_HAS_FF16=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_DISABLE_SSSE3_CODEGEN=1",
            "LEO2_DISABLE_AVX2_CODEGEN=1"))),
        ("mutation", ("leopard", "target_include_directories", (
            "PUBLIC", "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>",
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"))),
        ("mutation", ("leopard", "target_link_libraries", (
            "PUBLIC", "Threads::Threads"))),
        ("mutation", ("leopard2_backend_ssse3",
            "target_include_directories", (
                "PRIVATE", "${CMAKE_CURRENT_SOURCE_DIR}"))),
        ("trusted", ("check_cxx_compiler_flag", (
            "/arch:AVX2", "LEO2_FLAG_ARCH_AVX2"))),
        ("mutation", ("leopard2_backend_avx2",
            "target_include_directories", (
                "PRIVATE", "${CMAKE_CURRENT_SOURCE_DIR}"))),
        ("mutation", ("leopard2_backend_avx2", "target_compile_options", (
            "PRIVATE", "/arch:AVX2"))),
        ("mutation", ("leopard2_backend_ssse3",
            "target_compile_definitions", (
                "PRIVATE", "NO_LEO_HAS_FF8=1"))),
        ("mutation", ("leopard2_backend_ssse3",
            "target_compile_definitions", (
                "PRIVATE", "NO_LEO_HAS_FF16=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_HAVE_SSSE3_BACKEND=1"))),
        ("mutation", ("leopard2_backend_avx2",
            "target_compile_definitions", (
                "PRIVATE", "NO_LEO_HAS_FF8=1"))),
        ("mutation", ("leopard2_backend_avx2",
            "target_compile_definitions", (
                "PRIVATE", "NO_LEO_HAS_FF16=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_HAVE_AVX2_BACKEND=1"))),
        ("mutation", ("leopard", "target_link_libraries", (
            "PUBLIC", "OpenMP::OpenMP_CXX"))),
        ("mutation", ("leopard", "target_link_libraries", (
            "PUBLIC", "${OpenMP_CXX_FLAGS}"))),
        ("protected", (
            "LEOPARD_INSTALL_CMAKEDIR",
            "${CMAKE_INSTALL_LIBDIR}/cmake/leopard", "CACHE", "STRING",
            "Install directory for Leopard CMake package files")),
        ("trusted", ("configure_package_config_file",
            _approved_package_configure)),
        ("trusted", ("install", (
            "TARGETS", "leopard", "EXPORT", "leopardTargets", "ARCHIVE",
            "DESTINATION", "${CMAKE_INSTALL_LIBDIR}", "LIBRARY",
            "DESTINATION", "${CMAKE_INSTALL_LIBDIR}", "RUNTIME",
            "DESTINATION", "${CMAKE_INSTALL_BINDIR}"))),
        ("trusted", ("install", (
            "FILES", "leopard.h", "leopard2.h", "DESTINATION",
            "${CMAKE_INSTALL_INCLUDEDIR}"))),
        ("trusted", ("install", (
            "FILES", "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake",
            "DESTINATION", "${LEOPARD_INSTALL_CMAKEDIR}"))),
        ("trusted", ("install", (
            "EXPORT", "leopardTargets", "FILE", "leopardTargets.cmake",
            "NAMESPACE", "leopard::", "DESTINATION",
            "${LEOPARD_INSTALL_CMAKEDIR}"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_SCALAR=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_SSSE3=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_AVX2=1"))),
        ("mutation", ("leopard", "target_compile_definitions", (
            "PRIVATE", "LEO2_BACKEND_FORCE_AVX512=1"))),
        ("locator-provenance", ("find_program", _locator_git_find)),
        ("locator-provenance", (
            "execute_process", _locator_git_revision)),
        ("locator-provenance", (
            "execute_process", _locator_git_tree)),
        ("locator-provenance", (
            "execute_process", _locator_git_status)),
        ("sparse-sidecar", (
            "set_property", _sparse_sidecar_link_depends)),
        ("sparse-sidecar", (
            "add_custom_command", _sparse_sidecar_post_build)),
    )
    _dangerous_build_properties = {
        "COMPILE_DEFINITIONS", "COMPILE_FEATURES", "COMPILE_FLAGS",
        "COMPILE_OPTIONS", "CXX_COMPILER_LAUNCHER", "CXX_EXTENSIONS",
        "CXX_STANDARD", "CXX_STANDARD_REQUIRED", "IMPORTED_IMPLIB",
        "IMPORTED_LIBNAME", "IMPORTED_LINK_DEPENDENT_LIBRARIES",
        "IMPORTED_LINK_INTERFACE_LIBRARIES", "IMPORTED_LOCATION",
        "IMPORTED_OBJECTS", "INCLUDE_DIRECTORIES",
        "INTERFACE_COMPILE_DEFINITIONS", "INTERFACE_COMPILE_FEATURES",
        "INTERFACE_COMPILE_OPTIONS", "INTERFACE_INCLUDE_DIRECTORIES",
        "INTERFACE_SOURCES",
        "INTERFACE_LINK_DIRECTORIES", "INTERFACE_LINK_LIBRARIES",
        "INTERFACE_LINK_LIBRARIES_DIRECT",
        "INTERFACE_LINK_LIBRARIES_DIRECT_EXCLUDE",
        "INTERFACE_LINK_OPTIONS", "INTERFACE_PRECOMPILE_HEADERS",
        "INTERPROCEDURAL_OPTIMIZATION", "LINK_DIRECTORIES", "LINK_FLAGS",
        "LINK_INTERFACE_LIBRARIES", "LINK_LIBRARIES", "LINK_OPTIONS",
        "LINKER_TYPE", "MSVC_RUNTIME_LIBRARY", "POSITION_INDEPENDENT_CODE",
        "PRECOMPILE_HEADERS",
        "RULE_LAUNCH_COMPILE", "RULE_LAUNCH_CUSTOM", "RULE_LAUNCH_LINK",
        "STATIC_LIBRARY_OPTIONS", "UNITY_BUILD",
    }
    _dangerous_build_property_prefixes = (
        "COMPILE_DEFINITIONS_", "IMPORTED_IMPLIB_", "IMPORTED_LOCATION_",
        "IMPORTED_OBJECTS_", "INTERPROCEDURAL_OPTIMIZATION_", "LINK_FLAGS_",
        "VS_",
    )
    _approved_target_property_commands = {
        ("set_target_properties", (
            "leopard2_codec_options_abi_test", "PROPERTIES",
            "LINKER_LANGUAGE", "CXX")),
        ("set_property", (
            "TARGET", "leopard2_api_fuzzer", "APPEND_STRING", "PROPERTY",
            "LINK_FLAGS", " -fsanitize=fuzzer,address,undefined")),
        ("set_property", (
            "TARGET", "leopard2_pruned_plan_fuzzer", "APPEND_STRING",
            "PROPERTY", "LINK_FLAGS",
            " -fsanitize=fuzzer,address,undefined")),
        ("set_property", (
            "TARGET", "leopard2_code_identity_fuzzer", "APPEND_STRING",
            "PROPERTY", "LINK_FLAGS",
            " -fsanitize=fuzzer,address,undefined")),
    }
    _approved_nontarget_property_commands = {
        ("set_property", (
            "CACHE", "LEO2_BACKEND_VARIANT", "PROPERTY", "STRINGS",
            "auto", "scalar", "ssse3", "avx2", "avx512")),
    }

    def __init__(self, text, processor="AMD64", pointer_size="8",
                 platform_name="x64", require_mutation_contract=False):
        self.raw_commands = cmake_commands(text)
        self.commands = [(name, cmake_tokens(body))
                         for name, body in self.raw_commands]
        if (require_mutation_contract and
                self.commands[:2] != [
                    ("cmake_minimum_required", ["VERSION", "3.7"]),
                    ("project", ["leopard"])]):
            raise ContractError(
                "CMake bootstrap must begin with exact minimum and project "
                "commands")
        source_commands = {
            "add_library", "target_sources", "set_source_files_properties",
            "set_property", "set_target_properties"}
        self.source_variables = {
            match.group(1)
            for command, tokens in self.commands if command in source_commands
            for token in tokens
            for match in [self._variable.match(token)] if match
        }
        self.variables = {}
        self.targets = {}
        self.declared_targets = set()
        self.target_aliases = {}
        self.attachments = {}
        self.target_build_mutations = []
        self.trusted_command_counts = Counter()
        self.locator_provenance_counts = Counter()
        self.sparse_sidecar_counts = Counter()
        self.protected_assignments = []
        self.contract_events = []
        self.require_mutation_contract = require_mutation_contract
        self.source_property_references = []
        # Visual Studio 14 and v140 are the repository's declared legacy
        # configuration.  Compiler capability probes remain symbolic.
        for name, value in {
                "WIN32": "TRUE",
                "MSVC": "TRUE",
                "CMAKE_CXX_COMPILER_ID": "MSVC",
                "CMAKE_C_COMPILER_ID": "MSVC",
                "CMAKE_SYSTEM_PROCESSOR": processor,
                "CMAKE_SIZEOF_VOID_P": pointer_size,
                "CMAKE_VS_PLATFORM_NAME": platform_name,
                "CMAKE_GENERATOR_PLATFORM": platform_name}.items():
            self.variables[name] = [(BOOL_TRUE, (value,), ())]
        self.variables["CMAKE_CURRENT_SOURCE_DIR"] = [
            (BOOL_TRUE, (".",), ())]
        self._read_graph()

    @staticmethod
    def _merge_reasons(*reason_groups):
        merged = []
        for reasons in reason_groups:
            for reason in reasons:
                if reason not in merged:
                    merged.append(reason)
        return tuple(merged)

    def _assign(self, name, value, guard, reasons=()):
        if not bool_satisfiable(guard):
            return
        reasons = tuple(reasons)
        if bool_tautology(guard):
            self.variables[name] = [(BOOL_TRUE, value, reasons)]
            return
        previous = self.variables.get(name, [(BOOL_TRUE, None, ())])
        updated = []
        for old_guard, old_value, old_reasons in previous:
            replaced = bool_and(old_guard, guard)
            retained = bool_and(old_guard, bool_not(guard))
            if bool_satisfiable(replaced):
                updated.append((replaced, value, reasons))
            if bool_satisfiable(retained):
                updated.append((retained, old_value, old_reasons))
        self.variables[name] = self._merge_variants(updated)

    @staticmethod
    def _merge_variants(variants):
        merged = {}
        for guard, value, reasons in variants:
            key = (value, tuple(reasons))
            merged[key] = bool_or(merged.get(key, BOOL_FALSE), guard)
        return [(guard, key[0], key[1])
                for key, guard in merged.items() if bool_satisfiable(guard)]

    def _variable_variants(self, name, active_guard):
        if name not in self.variables:
            return []
        variants = []
        for guard, value, reasons in self.variables[name]:
            overlap = bool_and(active_guard, guard)
            if bool_satisfiable(overlap):
                variants.append((overlap, value, reasons))
        return variants

    def _unique_variable_value(self, name, active_guard):
        variants = self._variable_variants(name, active_guard)
        if not variants or any(value is None for unused, value, reasons in variants):
            raise ContractError("unresolved CMake source variable: " + name)
        for unused, unused_value, reasons in variants:
            if reasons:
                raise ContractError(reasons[0])
        values = {value for unused, value, unused_reasons in variants}
        if len(values) != 1:
            raise ContractError(
                "conditional CMake source variable is ambiguous: " + name)
        return next(iter(values))

    def _expand(self, tokens, active_guard, stack=None):
        stack = [] if stack is None else stack
        expanded = []
        for token in tokens:
            match = self._variable.match(token)
            if match:
                variable = match.group(1)
                if variable in stack:
                    raise ContractError(
                        "recursive CMake source variable: " + variable)
                values = self._unique_variable_value(variable, active_guard)
                expanded.extend(self._expand(
                    values, active_guard, stack + [variable]))
            elif "${" in token or "$ENV{" in token:
                raise ContractError(
                    "embedded or environment source variable is unsupported: " +
                    token)
            else:
                expanded.append(token)
        return expanded

    def _target_name(self, token, active_guard):
        names = self._expand([token], active_guard)
        if len(names) != 1 or not re.match(
                r"^[A-Za-z_][A-Za-z0-9_.+\-]*$", names[0]):
            raise ContractError("unresolved CMake target name: " + token)
        return names[0]

    def _canonical_target(self, name):
        visited = []
        while name in self.target_aliases:
            if name in visited:
                raise ContractError(
                    "cyclic CMake target ALIAS: " + " -> ".join(
                        visited + [name]))
            visited.append(name)
            name = self.target_aliases[name]
        return name

    def _record_target_build_mutation(self, command, tokens, guard, reasons):
        if not tokens:
            raise ContractError("CMake " + command + " has no target")
        target = self._target_name(tokens[0], guard)
        if target in self.target_aliases:
            raise ContractError(
                "CMake target ALIAS cannot be mutated: " + target)
        specification = tuple(tokens[1:])
        self.target_build_mutations.append(
            (target, command, specification, guard, tuple(reasons)))
        key = (target, command, specification)
        if key in self._approved_production_mutations:
            self.contract_events.append(("mutation", key))

    def _record_trusted_command(self, command, tokens, guard, reasons):
        key = (command, tuple(tokens))
        if key not in self._required_trusted_commands:
            return
        if reasons or not bool_tautology(guard):
            raise ContractError(
                "trusted CMake command guard drift: " + command + " " +
                repr(tuple(reasons)))
        self.trusted_command_counts[key] += 1
        self.contract_events.append(("trusted", key))

    def _record_locator_provenance_command(
            self, command, tokens, guard, reasons):
        key = (command, tuple(tokens))
        if key not in self._required_locator_provenance_commands:
            return False
        expected_guard = bool_atom("option:LEO2_BUILD_BENCHMARKS")
        if command == "execute_process":
            expected_guard = bool_and(
                expected_guard, bool_atom("locator-git-found"))
        if reasons or not self._formula_equivalent(guard, expected_guard):
            raise ContractError(
                "locator provenance command guard drift: " + command)
        self.locator_provenance_counts[key] += 1
        self.contract_events.append(("locator-provenance", key))
        return True

    def _record_sparse_sidecar_command(
            self, command, tokens, guard, reasons):
        key = (command, tuple(tokens))
        if key not in self._required_sparse_sidecar_commands:
            return False
        expected_guard = bool_and(
            bool_atom("option:LEO2_BUILD_BENCHMARKS"),
            bool_atom("option:LEOPARD_ENABLE_GF8"),
            bool_atom("option:LEOPARD_ENABLE_GF16"))
        if reasons or not self._formula_equivalent(guard, expected_guard):
            raise ContractError(
                "sparse evidence sidecar command guard drift: " + command)
        self.sparse_sidecar_counts[key] += 1
        self.contract_events.append(("sparse-sidecar", key))
        return True

    @classmethod
    def _is_protected_variable(cls, name):
        return cls._protected_variable.match(name) is not None

    @classmethod
    def _is_dangerous_build_property(cls, name):
        upper = name.upper()
        return (upper in cls._dangerous_build_properties or
                upper.startswith(cls._dangerous_build_property_prefixes))

    @staticmethod
    def _backend_variant_comparison(value):
        symbol = (
            BOOL_SYMBOL_PREFIX +
            "external-cache:leo2_backend_variant",)
        return bool_atom(
            "comparison:" + repr((symbol, "STREQUAL", (value,))))

    @classmethod
    def _expected_production_mutation_guard(cls, key):
        target, command, specification = key
        disabled_fields = {
            ("PRIVATE", "NO_LEO_HAS_FF8=1"): "LEOPARD_ENABLE_GF8",
            ("PRIVATE", "NO_LEO_HAS_FF16=1"): "LEOPARD_ENABLE_GF16",
        }
        if specification in disabled_fields:
            guard = bool_not(bool_atom(
                "option:" + disabled_fields[specification]))
            if target == "leopard2_backend_avx2":
                guard = bool_and(
                    guard, bool_atom("probe:LEO2_FLAG_ARCH_AVX2"))
            return guard
        forced_variants = {
            ("PRIVATE", "LEO2_BACKEND_FORCE_SCALAR=1"): "scalar",
            ("PRIVATE", "LEO2_BACKEND_FORCE_SSSE3=1"): "ssse3",
            ("PRIVATE", "LEO2_BACKEND_FORCE_AVX2=1"): "avx2",
            ("PRIVATE", "LEO2_BACKEND_FORCE_AVX512=1"): "avx512",
        }
        if specification in forced_variants:
            auto = cls._backend_variant_comparison("auto")
            scalar = cls._backend_variant_comparison("scalar")
            ssse3 = cls._backend_variant_comparison("ssse3")
            avx2 = cls._backend_variant_comparison("avx2")
            guard = bool_not(auto)
            variant = forced_variants[specification]
            if variant == "scalar":
                return bool_and(guard, scalar)
            guard = bool_and(guard, bool_not(scalar))
            if variant == "ssse3":
                return bool_and(guard, ssse3)
            guard = bool_and(guard, bool_not(ssse3))
            if variant == "avx2":
                return bool_and(guard, avx2)
            return bool_and(
                guard, bool_not(avx2),
                cls._backend_variant_comparison("avx512"))
        if (target == "leopard" and specification == (
                "PRIVATE",
                "LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE=1")):
            return bool_atom(
                "option:LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE")
        if (target == "leopard2_backend_avx2" or
                specification == ("PRIVATE", "LEO2_HAVE_AVX2_BACKEND=1")):
            return bool_atom("probe:LEO2_FLAG_ARCH_AVX2")
        if (command == "target_link_libraries" and specification ==
                ("PUBLIC", "OpenMP::OpenMP_CXX")):
            return bool_and(
                bool_atom("dependent-option:ENABLE_OPENMP"),
                bool_atom("predicate:TARGET:OpenMP::OpenMP_CXX"))
        if (command == "target_link_libraries" and specification ==
                ("PUBLIC", "${OpenMP_CXX_FLAGS}")):
            enabled = bool_atom("dependent-option:ENABLE_OPENMP")
            imported = bool_atom("predicate:TARGET:OpenMP::OpenMP_CXX")
            return bool_and(enabled, bool_not(bool_and(enabled, imported)))
        return BOOL_TRUE

    @staticmethod
    def _formula_equivalent(left, right):
        difference = bool_or(
            bool_and(left, bool_not(right)),
            bool_and(right, bool_not(left)))
        return not bool_satisfiable(difference)

    @classmethod
    def _validate_protected_assignment(cls, name, values, tokens):
        if (tuple(tokens) not in cls._approved_protected_set_commands or
                len(values) != 1 or (name, values[0]) not in
                cls._approved_protected_assignments):
            raise ContractError(
                "unapproved production compiler-control variable mutation: " +
                name)

    @staticmethod
    def _expected_protected_assignment_guard(name, values):
        value = values[0]
        probes = {
            "${CMAKE_CXX_FLAGS_RELEASE} /O2": "CXX_FLAG_O2",
            "${CMAKE_CXX_FLAGS_DEBUG} /Oy": "CXX_FLAG_Oy",
            "${CMAKE_CXX_FLAGS_DEBUG} /Zi": "CXX_FLAG_Zi",
            "${CMAKE_CXX_FLAGS} /W4": "CXX_FLAG_W4",
        }
        if value in probes:
            return bool_atom("probe:" + probes[value])
        if name in {"CMAKE_CXX_FLAGS", "CMAKE_EXE_LINKER_FLAGS"} and (
                "OpenMP" in value):
            return bool_atom("dependent-option:ENABLE_OPENMP")
        if name == "CMAKE_BUILD_TYPE":
            return bool_not(bool_atom("external:CMAKE_BUILD_TYPE"))
        return BOOL_TRUE

    def _validate_required_protected_assignments(self):
        observed = Counter(
            tokens
            for unused_name, unused_values, tokens, unused_guard,
            unused_reasons in
            self.protected_assignments)
        expected = Counter(self._approved_protected_set_commands)
        if observed != expected:
            missing = sorted((expected - observed).elements(), key=repr)
            extra = sorted((observed - expected).elements(), key=repr)
            raise ContractError(
                "missing or duplicate protected CMake assignment: missing=" +
                repr(missing) + " extra=" + repr(extra))
        for name, values, unused_tokens, guard, reasons in (
                self.protected_assignments):
            allowed_reasons = ()
            if name == "CMAKE_BUILD_TYPE":
                allowed_reasons = (
                    "unmodeled CMake conditional variable: CMAKE_BUILD_TYPE",)
            if reasons != allowed_reasons:
                raise ContractError(
                    "unsupported condition guards protected CMake assignment: " +
                    name + " " + repr(reasons))
            expected_guard = self._expected_protected_assignment_guard(
                name, values)
            if not self._formula_equivalent(guard, expected_guard):
                raise ContractError(
                    "protected CMake assignment guard drift: " + name)

    def _mutation_variable(self, token, active_guard, command):
        names = self._expand([token], active_guard)
        if len(names) != 1 or not re.match(
                r"^[A-Za-z_][A-Za-z0-9_]*$", names[0]):
            raise ContractError(
                "unsupported CMake " + command + " mutation destination: " +
                token)
        return names[0]

    def _expand_embedded_token(self, token, active_guard, stack=None):
        """Expand scalar ${VAR} fragments for mutation-destination discovery."""
        stack = [] if stack is None else stack
        if "$ENV{" in token:
            raise ContractError(
                "environment CMake variable is unsupported: " + token)
        match = self._embedded_variable.search(token)
        if not match:
            return token
        variable = match.group(1)
        if variable in stack:
            raise ContractError(
                "recursive CMake source variable: " + variable)
        values = self._unique_variable_value(variable, active_guard)
        replacement = ";".join(
            self._expand_embedded_token(
                value, active_guard, stack + [variable])
            for value in values)
        expanded = token[:match.start()] + replacement + token[match.end():]
        return self._expand_embedded_token(expanded, active_guard, stack)

    @staticmethod
    def _cmake_truth(value):
        if value is None:
            return BOOL_FALSE
        if len(value) == 1 and value[0].startswith(BOOL_SYMBOL_PREFIX):
            return bool_atom(value[0][len(BOOL_SYMBOL_PREFIX):])
        text = ";".join(value)
        upper = text.upper()
        if (not text or upper in {
                "0", "FALSE", "OFF", "NO", "N", "IGNORE", "NOTFOUND"} or
                upper.endswith("-NOTFOUND")):
            return BOOL_FALSE
        if upper in {"1", "TRUE", "ON", "YES", "Y"}:
            return BOOL_TRUE
        # CMake's numeric/truth parsing is its own language contract.  Values
        # outside the explicit constants above remain reachable rather than
        # being concretized with Python's numeric grammar.
        return bool_atom("truth:" + repr(text))

    def _condition_operand(self, token, active_guard, standalone=False):
        text, quoted = token
        explicit = self._variable.match(text)
        if explicit:
            name = explicit.group(1)
        elif quoted:
            return [(BOOL_TRUE, (text,), ())]
        elif text in self.variables:
            name = text
        elif standalone:
            if text.upper() in {
                    "0", "1", "FALSE", "TRUE", "OFF", "ON", "NO", "YES",
                    "N", "Y", "IGNORE", "NOTFOUND"}:
                return [(BOOL_TRUE, (text,), ())]
            reason = "unmodeled CMake conditional variable: " + text
            return [(BOOL_TRUE,
                     (BOOL_SYMBOL_PREFIX + "external:" + text,), (reason,))]
        else:
            constants = {
                "0", "1", "FALSE", "TRUE", "OFF", "ON", "NO", "YES",
                "N", "Y", "IGNORE", "NOTFOUND"}
            if (re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", text) and
                    text.upper() not in constants):
                reason = (
                    "unmodeled unquoted CMake conditional operand: " + text)
                return [(BOOL_TRUE,
                         (BOOL_SYMBOL_PREFIX + "external:" + text,),
                         (reason,))]
            return [(BOOL_TRUE, (text,), ())]
        variants = self._variable_variants(name, active_guard)
        if not variants:
            reason = "unmodeled CMake conditional variable: " + name
            return [(BOOL_TRUE,
                     (BOOL_SYMBOL_PREFIX + "external:" + name,), (reason,))]
        return [
            (guard, value, tuple(
                reason for reason in reasons
                if not reason.startswith(CONDITIONAL_ASSIGNMENT_PREFIX)))
            for guard, value, reasons in variants
        ]

    def _compare_condition_values(self, left, operation, right):
        left_symbol = (len(left) == 1 and
                       left[0].startswith(BOOL_SYMBOL_PREFIX))
        right_symbol = (len(right) == 1 and
                        right[0].startswith(BOOL_SYMBOL_PREFIX))
        if left_symbol or right_symbol:
            atom = "comparison:" + repr((left, operation, right))
            if operation == "STREQUAL" and left_symbol != right_symbol:
                symbol = left[0] if left_symbol else right[0]
                literal = right if left_symbol else left
                truth = self._cmake_truth(literal)
                if truth in (BOOL_TRUE, BOOL_FALSE):
                    formula = bool_atom(symbol[len(BOOL_SYMBOL_PREFIX):])
                    return formula if truth == BOOL_TRUE else bool_not(formula)
            return bool_atom(atom)
        left_text = ";".join(left)
        right_text = ";".join(right)
        if operation == "STREQUAL":
            return BOOL_TRUE if left_text == right_text else BOOL_FALSE
        if operation == "MATCHES":
            # CMake's regex engine and Python's are not dialect-equivalent.
            # Keep every MATCHES result symbolic so dialect differences cannot
            # make a real-reachable production mutation disappear from proof.
            return bool_atom(
                "matches:" + repr((left_text, right_text)))
        if operation == "EQUAL":
            return bool_atom(
                "numeric-equal:" + repr((left_text, right_text)))
        return bool_atom("unsupported-comparison:" +
                         repr((left_text, operation, right_text)))

    def _condition_comparison(self, left, operation, right, active_guard):
        formula = BOOL_FALSE
        reasons = []
        left_values = self._condition_operand(left, active_guard)
        right_values = self._condition_operand(right, active_guard)
        for left_guard, left_value, left_reasons in left_values:
            for right_guard, right_value, right_reasons in right_values:
                overlap = bool_and(active_guard, left_guard, right_guard)
                if not bool_satisfiable(overlap):
                    continue
                compared = self._compare_condition_values(
                    left_value or (), operation, right_value or ())
                formula = bool_or(formula, bool_and(overlap, compared))
                reasons.extend(left_reasons)
                reasons.extend(right_reasons)
                if any(len(value) == 1 and
                       value[0].startswith(BOOL_SYMBOL_PREFIX)
                       for value in (left_value or (), right_value or ())):
                    reasons.append(
                        "unsupported comparison of symbolic CMake boolean")
                if operation == "MATCHES":
                    reasons.append(
                        "CMake MATCHES regex dialect is intentionally symbolic")
                elif operation == "EQUAL":
                    reasons.append(
                        "CMake numeric comparison is intentionally symbolic")
                elif operation != "STREQUAL":
                    reasons.append(
                        "unsupported CMake comparison operator: " + operation)
        return formula, self._merge_reasons(reasons)

    def _condition_truth(self, token, active_guard):
        formula = BOOL_FALSE
        reasons = []
        for guard, value, value_reasons in self._condition_operand(
                token, active_guard, standalone=True):
            formula = bool_or(
                formula, bool_and(active_guard, guard,
                                  self._cmake_truth(value)))
            reasons.extend(value_reasons)
        return formula, self._merge_reasons(reasons)

    def _eval_condition(self, body, active_guard):
        tokens = cmake_condition_tokens(body)
        position = [0]
        comparisons = {
            "STREQUAL", "MATCHES", "EQUAL", "LESS", "GREATER",
            "LESS_EQUAL", "GREATER_EQUAL", "VERSION_LESS",
            "VERSION_GREATER", "VERSION_EQUAL", "VERSION_LESS_EQUAL",
            "VERSION_GREATER_EQUAL", "IN_LIST", "IS_NEWER_THAN"}
        unary_predicates = {
            "COMMAND", "POLICY", "TARGET", "TEST", "DEFINED", "EXISTS",
            "IS_READABLE", "IS_WRITABLE", "IS_EXECUTABLE", "IS_DIRECTORY",
            "IS_SYMLINK", "IS_ABSOLUTE"}

        def accept(operator):
            if position[0] >= len(tokens):
                return False
            text, quoted = tokens[position[0]]
            if not quoted and text.upper() == operator:
                position[0] += 1
                return True
            return False

        def parse_primary():
            if accept("("):
                result = parse_or()
                if not accept(")"):
                    raise ValueError("missing closing parenthesis")
                return result
            if position[0] >= len(tokens):
                raise ValueError("missing condition operand")
            token = tokens[position[0]]
            position[0] += 1
            upper = token[0].upper() if not token[1] else ""
            if upper in unary_predicates:
                if position[0] >= len(tokens):
                    raise ValueError("missing unary predicate operand")
                operand = tokens[position[0]]
                position[0] += 1
                atom = bool_atom("predicate:" + upper + ":" + operand[0])
                return atom, (
                    "unsupported CMake conditional predicate: " + upper,)
            if position[0] < len(tokens):
                operation, quoted = tokens[position[0]]
                operation = operation.upper() if not quoted else ""
                if operation in comparisons:
                    position[0] += 1
                    if position[0] >= len(tokens):
                        raise ValueError("missing comparison operand")
                    right = tokens[position[0]]
                    position[0] += 1
                    return self._condition_comparison(
                        token, operation, right, active_guard)
            return self._condition_truth(token, active_guard)

        def parse_not():
            if accept("NOT"):
                formula, reasons = parse_not()
                return bool_not(formula), reasons
            return parse_primary()

        def parse_and():
            formula, reasons = parse_not()
            while accept("AND"):
                right, right_reasons = parse_not()
                formula = bool_and(formula, right)
                reasons = self._merge_reasons(reasons, right_reasons)
            return formula, reasons

        def parse_or():
            formula, reasons = parse_and()
            while accept("OR"):
                right, right_reasons = parse_and()
                formula = bool_or(formula, right)
                reasons = self._merge_reasons(reasons, right_reasons)
            return formula, reasons

        try:
            formula, reasons = parse_or()
            if position[0] != len(tokens):
                raise ValueError("trailing condition tokens")
            return formula, reasons
        except ValueError as error:
            expression = " ".join(token for token, quoted in tokens)
            return (bool_atom("unsupported-condition:" + expression),
                    ("unsupported CMake conditional structure: " + str(error),))

    def _assignment_values(self, tokens, active_guard):
        values = []
        reasons = []
        for token in tokens:
            match = self._variable.match(token)
            if match:
                name = match.group(1)
                try:
                    values.extend(self._unique_variable_value(
                        name, active_guard))
                except ContractError as error:
                    reasons.append(str(error))
            elif "${" in token or "$ENV{" in token:
                reasons.append(
                    "embedded or environment CMake variable is unsupported: " +
                    token)
            else:
                values.append(token)
        return tuple(values), tuple(reasons)

    def _read_graph(self):
        starts = {
            "function": "endfunction", "macro": "endmacro",
            "foreach": "endforeach", "while": "endwhile",
            "block": "endblock"}
        ends = {value: key for key, value in starts.items()}
        stack = []
        guard = BOOL_TRUE
        reasons = ()
        unsupported_depth = 0
        conditional_depth = 0
        approved_includes = {
            "CMakeDependentOption", "CheckCXXCompilerFlag",
            "CMakePackageConfigHelpers", "GNUInstallDirs"}
        approved_packages = {
            ("OpenMP",),
            ("Threads", "REQUIRED"),
            ("PythonInterp", "QUIET"),
            ("Python3", "COMPONENTS", "Interpreter", "QUIET"),
        }
        for command, body in self.raw_commands:
            tokens = cmake_tokens(body)
            if any("$CACHE{" in token or "$ENV{" in token
                   for token in tokens):
                raise ContractError(
                    "external cache/environment expansion is unsupported in "
                    "the production graph proof")
            if command not in (
                    self._modeled_command_names |
                    self._target_build_mutation_commands |
                    self._directory_build_mutation_commands |
                    self._build_extension_commands |
                    self._safe_root_command_names):
                raise ContractError(
                    "unmodeled root CMake command: " + command)
            if command == "if":
                conditional_depth += 1
                if unsupported_depth:
                    stack.append({"type": "skipped-if"})
                    continue
                condition, condition_reasons = self._eval_condition(body, guard)
                branch = bool_and(guard, condition)
                stack.append({
                    "type": "if", "parent_guard": guard,
                    "parent_reasons": reasons, "taken": branch,
                    "decision_reasons": condition_reasons})
                guard = branch
                reasons = self._merge_reasons(reasons, condition_reasons)
                continue
            if command in ("elseif", "else"):
                if not stack or stack[-1]["type"] not in {"if", "skipped-if"}:
                    raise ContractError("unbalanced CMake " + command)
                if stack[-1]["type"] == "skipped-if":
                    continue
                frame = stack[-1]
                available = bool_and(
                    frame["parent_guard"], bool_not(frame["taken"]))
                if command == "elseif":
                    condition, condition_reasons = self._eval_condition(
                        body, available)
                    guard = bool_and(available, condition)
                    frame["taken"] = bool_or(frame["taken"], guard)
                    frame["decision_reasons"] = self._merge_reasons(
                        frame["decision_reasons"], condition_reasons)
                    reasons = self._merge_reasons(
                        frame["parent_reasons"],
                        frame["decision_reasons"])
                else:
                    guard = available
                    frame["taken"] = bool_or(frame["taken"], guard)
                    reasons = self._merge_reasons(
                        frame["parent_reasons"],
                        frame["decision_reasons"])
                continue
            if command == "endif":
                conditional_depth -= 1
                if conditional_depth < 0 or not stack or stack[-1]["type"] not in {
                        "if", "skipped-if"}:
                    raise ContractError("unbalanced CMake endif")
                frame = stack.pop()
                if frame["type"] == "if":
                    guard = frame["parent_guard"]
                    reasons = frame["parent_reasons"]
                continue
            if command in starts:
                if (command in {"function", "macro"} and tokens and
                        tokens[0].lower() in (
                            self._target_build_mutation_commands |
                            self._directory_build_mutation_commands |
                            self._build_extension_commands |
                            self._modeled_command_names)):
                    raise ContractError(
                        "CMake build-graph command shadowing is unsupported: " +
                        tokens[0])
                raise ContractError("unsupported CMake block: " + command)
            if command in ends:
                if not stack or stack[-1]["type"] != ends[command]:
                    raise ContractError("unbalanced CMake " + command)
                stack.pop()
                unsupported_depth -= 1
                continue

            if unsupported_depth:
                if (command == "add_library" and tokens and
                        tokens[0] == "leopard"):
                    raise ContractError(
                        "leopard add_library must be unconditional")
                if command in {
                        "add_library", "target_sources",
                        "set_source_files_properties", "set_property",
                        "set_target_properties", "set_directory_properties"}:
                    raise ContractError(
                        "source graph command in unsupported CMake block: " +
                        command)
                if command in {
                        "set", "unset", "list", "string", "file", "math",
                        "separate_arguments", "cmake_path", "get_property",
                        "get_cmake_property", "get_directory_property",
                        "get_filename_component", "get_source_file_property",
                        "get_target_property", "execute_process", "try_compile",
                        "try_run", "cmake_parse_arguments"}:
                    raise ContractError(
                        "variable mutation in unsupported CMake block: " +
                        command)
                if command in {
                        "include", "add_subdirectory", "subdirs",
                        "cmake_language", "load_command", "find_package",
                        "project", "enable_language"}:
                    raise ContractError(
                        "graph extension in unsupported CMake block: " + command)
                if command in self._target_build_mutation_commands:
                    raise ContractError(
                        "compile/link mutation in unsupported CMake block: " +
                        command)
                if command in self._directory_build_mutation_commands:
                    raise ContractError(
                        "directory compile/link mutation in unsupported CMake "
                        "block: " + command)
                if command in self._build_extension_commands:
                    raise ContractError(
                        "build extension in unsupported CMake block: " + command)
                continue

            cuda_guard = bool_atom("option:LEO2_ENABLE_CUDA")
            default_reachable = bool_satisfiable(
                bool_and(guard, bool_not(cuda_guard)))
            if default_reachable and command == "find_program":
                try:
                    program_tokens = self._expand(tokens, guard)
                except ContractError as error:
                    raise ContractError(
                        "unresolved tool discovery is reachable in the "
                        "CPU-only graph: " + str(error))
                approved_program_searches = {
                    self._locator_git_find,
                    ("LEO2_OBJDUMP_EXECUTABLE", "NAMES", "objdump",
                     "llvm-objdump"),
                    ("LEO2_SH_EXECUTABLE", "NAMES", "sh"),
                }
                if tuple(program_tokens) not in approved_program_searches:
                    raise ContractError(
                        "unapproved tool discovery is reachable in the CPU-only "
                        "graph")
            is_source_property = (
                command == "set_source_files_properties" or
                (command == "set_property" and tokens and
                 tokens[0].upper() == "SOURCE"))
            if default_reachable and is_source_property:
                try:
                    source_property_tokens = self._expand(tokens, guard)
                except ContractError as error:
                    raise ContractError(
                        "unresolved source property is reachable in the "
                        "CPU-only graph: " + str(error))
                upper_source_properties = {
                    token.upper() for token in source_property_tokens}
                if any("$<" in token for token in source_property_tokens):
                    raise ContractError(
                        "unmodeled source-property generator expression is "
                        "reachable in the CPU-only graph")
                if ("LANGUAGE" in upper_source_properties and
                        "CUDA" in upper_source_properties) or any(
                            self._cuda_graph_marker.search(token)
                            for token in source_property_tokens):
                    raise ContractError(
                        "CUDA source property is reachable in the CPU-only "
                        "graph")
            if default_reachable and command == "add_test":
                approved_cuda_optional_test = (
                    "NAME", "leopard2_cuda_optional", "COMMAND",
                    "${CMAKE_COMMAND}",
                    "-DLEO2_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}",
                    "-DLEO2_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/"
                    "cuda-optional-test",
                    "-DLEO2_GENERATOR=${CMAKE_GENERATOR}",
                    "-DLEO2_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}",
                    "-DLEO2_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}",
                    "-DLEO2_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}",
                    "-DLEO2_STATIC_LIBRARY_PREFIX="
                    "${CMAKE_STATIC_LIBRARY_PREFIX}",
                    "-DLEO2_STATIC_LIBRARY_SUFFIX="
                    "${CMAKE_STATIC_LIBRARY_SUFFIX}",
                    "-P", "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/"
                    "test_cuda_optional.cmake",
                )
                if tuple(tokens) == approved_cuda_optional_test:
                    test_tokens = []
                else:
                    try:
                        command_index = tokens.index("COMMAND")
                    except ValueError:
                        command_index = 0
                    test_tokens = tokens[command_index + 1:]
                expanded_test_tokens = []
                for token in test_tokens:
                    if self._variable.fullmatch(token):
                        try:
                            expanded_test_tokens.extend(
                                self._expand([token], guard))
                        except ContractError:
                            expanded_test_tokens.append(token)
                    elif "${" in token:
                        try:
                            expanded_test_tokens.append(
                                self._expand_embedded_token(token, guard))
                        except ContractError:
                            expanded_test_tokens.append(token)
                    else:
                        expanded_test_tokens.append(token)
                approved_test_generators = {"$<TARGET_FILE:bench_leopard2>"}
                if any("$<" in token and
                       token not in approved_test_generators
                       for token in expanded_test_tokens):
                    raise ContractError(
                        "unmodeled test generator expression is reachable in "
                        "the CPU-only graph")
                if any(self._cuda_graph_marker.search(token)
                       for token in expanded_test_tokens):
                    raise ContractError(
                        "CUDA test command is reachable in the CPU-only graph")
            inspected_tokens = tokens
            if default_reachable and command in {
                    "add_library", "add_executable", "target_sources"}:
                try:
                    inspected_tokens = self._expand(tokens, guard)
                except ContractError as error:
                    raise ContractError(
                        "unresolved source graph is reachable in the "
                        "CPU-only graph: " + str(error))
                if any("$<" in token and
                       self._target_objects.fullmatch(token) is None
                       for token in inspected_tokens):
                    raise ContractError(
                        "unmodeled source generator expression is reachable "
                        "in the CPU-only graph")
            if (default_reachable and
                    command in self._target_build_mutation_commands):
                approved_dependency_generators = {
                    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>",
                    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>",
                }
                if any("$<" in token and
                       token not in approved_dependency_generators
                       for token in tokens):
                    raise ContractError(
                        "unmodeled target dependency generator expression "
                        "is reachable in the CPU-only graph")
                dependency_tokens = []
                for token in tokens:
                    if token in approved_dependency_generators:
                        dependency_tokens.append(token)
                        continue
                    if self._variable.fullmatch(token):
                        try:
                            dependency_tokens.extend(
                                self._expand([token], guard))
                        except ContractError:
                            # Unknown package/toolchain variables are an
                            # external trust boundary.  Locally assigned
                            # indirection resolves above and is inspected.
                            dependency_tokens.append(token)
                            continue
                    elif "${" in token:
                        try:
                            dependency_tokens.append(
                                self._expand_embedded_token(token, guard))
                        except ContractError:
                            dependency_tokens.append(token)
                    else:
                        dependency_tokens.append(token)
                if any("$<" in token and
                       token not in approved_dependency_generators
                       for token in dependency_tokens):
                    raise ContractError(
                        "expanded target dependency generator expression is "
                        "reachable in the CPU-only graph")
                if any(re.search(
                        r"cuda|nvidia|cudart|nvrtc", token, re.IGNORECASE)
                       for token in dependency_tokens):
                    raise ContractError(
                        "CUDA target dependency is reachable in the CPU-only "
                        "graph")
            if default_reachable and any(
                    re.search(r"\.cuh?(?:$|[;>])", token, re.IGNORECASE) or
                    re.search(
                        r"(?:^|[/\\])(?:cuda(?:_[A-Za-z0-9_]+)?|nvrtc)\.h$",
                        token, re.IGNORECASE)
                    for token in inspected_tokens):
                raise ContractError(
                    "CUDA source/header is reachable in the CPU-only graph")

            if command == "cmake_minimum_required":
                if tokens != ["VERSION", "3.7"]:
                    raise ContractError(
                        "unapproved CMake minimum-version contract")
                if bool_satisfiable(guard):
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                continue
            if command == "include":
                if len(tokens) != 1 or tokens[0] not in approved_includes:
                    raise ContractError(
                        "unapproved CMake graph include: " + " ".join(tokens))
                if bool_satisfiable(guard):
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                continue
            if command == "find_package":
                if tuple(tokens) not in approved_packages:
                    raise ContractError(
                        "unapproved CMake package graph import: " +
                        " ".join(tokens))
                if (bool_satisfiable(guard) and tuple(tokens) in {
                        ("OpenMP",), ("Threads", "REQUIRED")}):
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                continue
            if command in {"add_subdirectory", "subdirs"}:
                raise ContractError(
                    "CMake " + command + " requires recursive graph proof")
            if command in {"cmake_language", "load_command"}:
                raise ContractError(
                    "CMake dynamic command execution is unsupported: " + command)
            if (command == "project" and bool_satisfiable(guard) and
                    tokens != ["leopard"]):
                raise ContractError(
                    "unapproved CMake project language/toolchain mutation")
            if command == "project" and bool_satisfiable(guard):
                self._record_trusted_command(command, tokens, guard, reasons)
                continue
            if command == "enable_language" and bool_satisfiable(guard):
                cuda_guard = bool_atom("option:LEO2_ENABLE_CUDA")
                if (tokens == ["CUDA"] and not reasons and
                        self._formula_equivalent(guard, cuda_guard)):
                    self.contract_events.append(
                        ("guarded", (command, tuple(tokens))))
                    continue
                raise ContractError(
                    "reachable CMake language enablement is unsupported")
            if command == "install" and bool_satisfiable(guard):
                key = (command, tuple(tokens))
                if key not in self._required_trusted_commands:
                    raise ContractError(
                        "unapproved CMake install/package mutation")
                self._record_trusted_command(command, tokens, guard, reasons)
                continue
            if (command == "configure_package_config_file" and
                    bool_satisfiable(guard)):
                if tuple(tokens) != self._approved_package_configure:
                    raise ContractError(
                        "unapproved generated package configuration command")
                self._record_trusted_command(command, tokens, guard, reasons)
                continue
            if (command in self._directory_build_mutation_commands and
                    bool_satisfiable(guard)):
                raise ContractError(
                    "CMake directory compile/link graph requires recursive "
                    "proof: " + command)
            if (command == "find_program" and bool_satisfiable(guard) and
                    self._record_locator_provenance_command(
                        command, tokens, guard, reasons)):
                self._assign(
                    "LEO2_LOCATOR_GIT_EXECUTABLE",
                    (BOOL_SYMBOL_PREFIX + "locator-git-found",), guard)
                continue
            if (command == "execute_process" and bool_satisfiable(guard) and
                    self._record_locator_provenance_command(
                        command, tokens, guard, reasons)):
                if tuple(tokens) == self._locator_git_revision:
                    result_variable = "LEO2_LOCATOR_GIT_RESULT"
                    output_variable = "LEO2_LOCATOR_GIT_OUTPUT"
                    symbol = "locator-git-revision"
                elif tuple(tokens) == self._locator_git_tree:
                    result_variable = "LEO2_LOCATOR_TREE_RESULT"
                    output_variable = "LEO2_LOCATOR_TREE_OUTPUT"
                    symbol = "locator-git-tree"
                else:
                    result_variable = "LEO2_LOCATOR_STATUS_RESULT"
                    output_variable = "LEO2_LOCATOR_STATUS_OUTPUT"
                    symbol = "locator-git-status"
                self._assign(
                    result_variable,
                    (BOOL_SYMBOL_PREFIX + symbol + "-result",), guard)
                self._assign(
                    output_variable,
                    (BOOL_SYMBOL_PREFIX + symbol + "-output",), guard)
                continue
            if (command in {"set_property", "add_custom_command"} and
                    bool_satisfiable(guard) and
                    self._record_sparse_sidecar_command(
                        command, tokens, guard, reasons)):
                continue
            if (command in self._build_extension_commands and
                    bool_satisfiable(guard)):
                raise ContractError(
                    "CMake generated/custom build extension requires recursive "
                    "proof: " + command)

            if command == "string" and bool_satisfiable(guard):
                key = (command, tuple(tokens))
                if key in self._required_trusted_commands:
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                    value = self._unique_variable_value(
                        "LEO2_BACKEND_VARIANT", guard)
                    if len(value) != 1:
                        raise ContractError(
                            "backend variant normalizer input is not scalar")
                    self._assign(
                        "LEO2_BACKEND_VARIANT_NORMALIZED",
                        (value[0].lower(),), guard, reasons)
                    continue

            if command == "set" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH can redirect approved graph includes")
                raw_values = list(tokens[1:])
                upper = [value.upper() for value in raw_values]
                has_cache = "CACHE" in upper
                if (bool_satisfiable(guard) and
                        "PARENT_SCOPE" in upper):
                    raise ContractError(
                        "CMake PARENT_SCOPE assignment is unsupported: " +
                        variable)
                if "CACHE" in upper:
                    raw_values = raw_values[:upper.index("CACHE")]
                if raw_values and raw_values[-1].upper() == "PARENT_SCOPE":
                    raw_values.pop()
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable)):
                    self._validate_protected_assignment(
                        variable, raw_values, tokens)
                    self.protected_assignments.append(
                        (variable, tuple(raw_values), tuple(tokens), guard,
                         tuple(reasons)))
                    self.contract_events.append(("protected", tuple(tokens)))
                elif bool_satisfiable(guard) and has_cache:
                    raise ContractError(
                        "unapproved CMake cache assignment: " + variable)
                if has_cache and "FORCE" not in upper:
                    # A non-FORCE cache initializer is only a default.  A
                    # caller's -D value survives and must keep every branch
                    # depending on it reachable in the proof.
                    value = (
                        BOOL_SYMBOL_PREFIX + "external-cache:" + variable,)
                    value_reasons = ()
                elif raw_values:
                    value, value_reasons = self._assignment_values(
                        raw_values, guard)
                else:
                    value = (BOOL_SYMBOL_PREFIX + "external:" + variable,)
                    value_reasons = (
                        "empty CMake set may reveal external cache variable: " +
                        variable,)
                assignment_reasons = self._merge_reasons(reasons, value_reasons)
                if conditional_depth:
                    assignment_reasons = self._merge_reasons(
                        assignment_reasons,
                        (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                self._assign(variable, value, guard, assignment_reasons)
                continue
            if command == "unset" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable)):
                    raise ContractError(
                        "unapproved production compiler-control variable "
                        "mutation: " + variable)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH mutation is unsupported")
                if bool_satisfiable(guard):
                    raise ContractError(
                        "CMake unset state/cache fallback is unsupported: " +
                        variable)
                continue
            if command == "option" and len(tokens) >= 2:
                variable = self._mutation_variable(tokens[0], guard, command)
                trusted_option = (
                    command, tuple(tokens)) in self._required_trusted_commands
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable) and
                        not trusted_option):
                    raise ContractError(
                        "unapproved production compiler-control variable "
                        "mutation: " + variable)
                if bool_satisfiable(guard) and not trusted_option:
                    raise ContractError(
                        "unapproved CMake cache option: " + variable)
                if bool_satisfiable(guard) and trusted_option:
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                option_reasons = reasons
                if conditional_depth:
                    option_reasons = self._merge_reasons(
                        option_reasons,
                        (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                # option() also preserves an existing cache value.  Its
                # authored default cannot prove the normal-command-line value.
                self._assign(
                    variable,
                    (BOOL_SYMBOL_PREFIX + "option:" + variable,),
                    guard, option_reasons)
                continue
            if command in {"check_cxx_compiler_flag",
                           "check_c_compiler_flag"} and len(tokens) >= 2:
                variable = self._mutation_variable(
                    tokens[-1], guard, command)
                trusted_probe = (
                    command, tuple(tokens)) in self._required_trusted_commands
                if bool_satisfiable(guard) and not trusted_probe:
                    raise ContractError(
                        "unapproved production compiler probe: " + command +
                        " " + " ".join(tokens))
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable) and
                        not trusted_probe):
                    raise ContractError(
                        "compiler probe may not overwrite production compiler "
                        "control variable: " + variable)
                if bool_satisfiable(guard):
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                symbol = BOOL_SYMBOL_PREFIX + "probe:" + variable
                self._assign(variable, (symbol,), guard, reasons)
                continue
            if command == "cmake_dependent_option" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                trusted_option = (
                    command, tuple(tokens)) in self._required_trusted_commands
                if bool_satisfiable(guard) and not trusted_option:
                    raise ContractError(
                        "unapproved production dependent option: " +
                        " ".join(tokens))
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable) and
                        not trusted_option):
                    raise ContractError(
                        "dependent option may not overwrite production compiler "
                        "control variable: " + variable)
                if bool_satisfiable(guard):
                    self._record_trusted_command(
                        command, tokens, guard, reasons)
                symbol = BOOL_SYMBOL_PREFIX + "dependent-option:" + variable
                self._assign(variable, (symbol,), guard, reasons)
                continue
            if command == "list" and len(tokens) >= 2:
                operation = tokens[0].upper()
                variable = self._mutation_variable(tokens[1], guard, command)
                if (bool_satisfiable(guard) and
                        self._is_protected_variable(variable)):
                    raise ContractError(
                        "unapproved production compiler-control variable "
                        "mutation: " + variable)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH mutation is unsupported")
                if (bool_satisfiable(guard) and
                        operation not in {"APPEND", "PREPEND"}):
                    raise ContractError(
                        "unsupported CMake list operation " + operation +
                        " touches " + variable)
                if operation in {"APPEND", "PREPEND"}:
                    additions, addition_reasons = self._assignment_values(
                        tokens[2:], guard)
                    try:
                        previous = self._unique_variable_value(variable, guard)
                    except ContractError as error:
                        previous = ()
                        addition_reasons += (str(error),)
                    value = (previous + additions if operation == "APPEND" else
                             additions + previous)
                    list_reasons = self._merge_reasons(
                        reasons, addition_reasons)
                    if conditional_depth:
                        list_reasons = self._merge_reasons(
                            list_reasons,
                            (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                    self._assign(variable, value, guard, list_reasons)
                continue

            if (command == "target_sources" and tokens and
                    tokens[0] == "leopard" and
                    not bool_satisfiable(guard)):
                # This backend is intentionally GNU/Clang-only: MSVC does not
                # provide the same per-source AVX-512VL compile contract.  The
                # exact unreachable attachment therefore does not belong in
                # the hand-maintained Visual Studio project.
                if tokens != [
                        "leopard", "PRIVATE",
                        "$<TARGET_OBJECTS:leopard2_backend_avx512>"]:
                    raise ContractError(
                        "leopard TARGET_OBJECTS has no MSVC-reachable "
                        "definition or attachment configuration")
            if (command == "add_library" and tokens and
                    tokens[0] == "leopard" and conditional_depth):
                raise ContractError(
                    "leopard add_library must be unconditional")
            if not bool_satisfiable(guard):
                if command == "add_library" and tokens:
                    self.declared_targets.add(tokens[0])
                continue
            if command in self._target_build_mutation_commands:
                self._record_target_build_mutation(
                    command, tokens, guard, reasons)
                continue
            variable_writers = {
                "file", "string", "math", "separate_arguments",
                "cmake_path", "get_property", "get_cmake_property",
                "get_directory_property", "get_filename_component",
                "get_source_file_property", "get_target_property",
                "execute_process", "find_program", "try_compile", "try_run",
                "cmake_parse_arguments"}
            expanded_identifiers = set()
            for token in tokens:
                try:
                    expanded_token = self._expand_embedded_token(token, guard)
                except ContractError:
                    if command in variable_writers:
                        # An unknown computed writer destination may name any
                        # production source variable.  Poison each until its
                        # next explicit assignment rather than accepting it.
                        expanded_identifiers.update(self.source_variables)
                else:
                    if re.match(
                            r"^[A-Za-z_][A-Za-z0-9_]*$", expanded_token):
                        expanded_identifiers.add(expanded_token)
                match = self._variable.match(token)
                if not match:
                    continue
                for unused_guard, values, value_reasons in (
                        self._variable_variants(match.group(1), guard)):
                    if values is None or value_reasons:
                        continue
                    expanded_identifiers.update(
                        value for value in values
                        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", value))
            if "CMAKE_MODULE_PATH" in expanded_identifiers:
                raise ContractError(
                    "CMAKE_MODULE_PATH mutation can redirect package graph "
                    "imports")
            protected_identifiers = sorted(
                identifier for identifier in expanded_identifiers
                if self._is_protected_variable(identifier))
            property_commands = {
                "set_property", "set_source_files_properties",
                "set_target_properties", "set_directory_properties"}
            if protected_identifiers and command not in property_commands:
                raise ContractError(
                    "unmodeled command may mutate production compiler-control "
                    "variable: " + command + " " +
                    ", ".join(protected_identifiers))
            if command in variable_writers:
                for variable in expanded_identifiers:
                    reason = (
                        "unmodeled CMake command may mutate source variable: " +
                        command + " " + variable)
                    self._assign(variable, (), guard,
                                 self._merge_reasons(reasons, (reason,)))
            touched_source_variables = (
                expanded_identifiers & self.source_variables)
            graph_commands = {
                "add_library", "target_sources",
                "set_source_files_properties", "set_property",
                "set_target_properties"}
            if touched_source_variables and command not in graph_commands:
                for variable in touched_source_variables:
                    reason = (
                        "unmodeled CMake command may mutate source variable: " +
                        command + " " + variable)
                    self._assign(variable, (), guard,
                                 self._merge_reasons(reasons, (reason,)))
            self._read_graph_command(
                command, tokens, guard, reasons, conditional_depth)

        if stack:
            raise ContractError("unbalanced CMake block: " + stack[-1]["type"])
        if (self.require_mutation_contract and
                self.trusted_command_counts != self._required_trusted_commands):
            missing = self._required_trusted_commands - self.trusted_command_counts
            extra = self.trusted_command_counts - self._required_trusted_commands
            raise ContractError(
                "missing or duplicate trusted CMake command: missing=" +
                repr(sorted(missing.elements(), key=repr)) + " extra=" +
                repr(sorted(extra.elements(), key=repr)))
        if self.require_mutation_contract:
            if (self.locator_provenance_counts !=
                    self._required_locator_provenance_commands):
                missing = (self._required_locator_provenance_commands -
                           self.locator_provenance_counts)
                extra = (self.locator_provenance_counts -
                         self._required_locator_provenance_commands)
                raise ContractError(
                    "missing or duplicate locator provenance command: "
                    "missing=" +
                    repr(sorted(missing.elements(), key=repr)) + " extra=" +
                    repr(sorted(extra.elements(), key=repr)))
            if (self.sparse_sidecar_counts !=
                    self._required_sparse_sidecar_commands):
                missing = (self._required_sparse_sidecar_commands -
                           self.sparse_sidecar_counts)
                extra = (self.sparse_sidecar_counts -
                         self._required_sparse_sidecar_commands)
                raise ContractError(
                    "missing or duplicate sparse evidence sidecar command: "
                    "missing=" +
                    repr(sorted(missing.elements(), key=repr)) + " extra=" +
                    repr(sorted(extra.elements(), key=repr)))
            self._validate_required_protected_assignments()
            expected_events = self._required_contract_event_order
            if (Counter(self.contract_events) == Counter(expected_events) and
                    tuple(self.contract_events) != expected_events):
                raise ContractError(
                    "security-sensitive CMake command order drift")

    def _read_graph_command(self, command, raw_tokens, guard, reasons,
                            conditional_depth):
        library_types = {"STATIC", "SHARED", "MODULE", "OBJECT", "INTERFACE"}
        source_commands = {
            "add_library", "target_sources", "set_source_files_properties"}
        if reasons and command in source_commands:
            raise ContractError(
                "unsupported conditional CMake source graph: " + reasons[0])
        if command == "add_library" and raw_tokens:
            if raw_tokens[0] in {"OpenMP::OpenMP_CXX", "Threads::Threads"}:
                raise ContractError(
                    "approved package target cannot be locally declared: " +
                    raw_tokens[0])
            target = self._target_name(raw_tokens[0], guard)
            self.declared_targets.add(target)
            if target == "leopard" and not bool_tautology(guard):
                raise ContractError(
                    "leopard add_library must be unconditional")
            tokens = self._expand(raw_tokens[1:], guard)
            upper_tokens = [token.upper() for token in tokens]
            if upper_tokens and upper_tokens[0] == "ALIAS":
                if len(tokens) != 2:
                    raise ContractError(
                        "unsupported CMake target ALIAS declaration: " + target)
                if target in self.target_aliases:
                    raise ContractError(
                        "duplicate CMake target ALIAS: " + target)
                if target == "libleopard":
                    if tuple(raw_tokens) != (
                            "libleopard", "ALIAS", "leopard"):
                        raise ContractError(
                            "libleopard must remain an exact alias of leopard")
                    self._record_trusted_command(
                        command, raw_tokens, guard, reasons)
                self.target_aliases[target] = tokens[1]
                return
            if target == "libleopard":
                raise ContractError(
                    "libleopard must remain an exact alias of leopard")
            if "IMPORTED" in upper_tokens:
                return
            kind = "DEFAULT"
            if tokens and tokens[0].upper() in library_types:
                kind = tokens.pop(0).upper()
            if tokens and tokens[0].upper() == "EXCLUDE_FROM_ALL":
                tokens.pop(0)
            if target == "leopard":
                if (kind != "STATIC" or
                        tuple(raw_tokens) != (
                            "leopard", "STATIC", "${LIB_SOURCE_FILES}")):
                    raise ContractError(
                        "leopard must be one exact STATIC library definition")
                self._record_trusted_command(
                    command, raw_tokens, guard, reasons)
            definitions = self.targets.setdefault(target, [])
            for definition in definitions:
                if (definition[0] != kind or definition[1] != tokens):
                    if bool_satisfiable(bool_and(definition[2], guard)):
                        raise ContractError("conflicting CMake target: " + target)
                elif definition[0] == kind and definition[1] == tokens:
                    definition[2] = bool_or(definition[2], guard)
                    return
            definitions.append([kind, tokens, guard])
        elif command == "target_sources" and raw_tokens:
            target = self._target_name(raw_tokens[0], guard)
            if target in self.target_aliases:
                raise ContractError(
                    "CMake target ALIAS cannot be mutated: " + target)
            tokens = self._expand(raw_tokens[1:], guard)
            if any(token.upper() == "FILE_SET" for token in tokens):
                raise ContractError(
                    "CMake FILE_SET source attachment requires parser support: " +
                    target)
            sources = [token for token in tokens if token.upper() not in {
                "PRIVATE", "PUBLIC", "INTERFACE", "SYSTEM", "BEFORE"}]
            if (conditional_depth and target == "leopard" and any(
                    not self._target_objects.match(token) for token in sources)):
                field_sources = {
                    ("PRIVATE", "LeopardFF8.cpp", "LeopardFF8.h"):
                        bool_atom("option:LEOPARD_ENABLE_GF8"),
                    ("PRIVATE", "LeopardFF16.cpp", "LeopardFF16.h"):
                        bool_atom("option:LEOPARD_ENABLE_GF16"),
                }
                specification = tuple(raw_tokens[1:])
                expected_guard = field_sources.get(specification)
                if (expected_guard is None or reasons or
                        not self._formula_equivalent(guard, expected_guard)):
                    raise ContractError(
                        "conditional direct leopard source attachment")
            self.attachments.setdefault(target, []).extend(
                (source, guard) for source in sources)
        elif command == "set_source_files_properties" and raw_tokens:
            self._record_source_properties(command, raw_tokens, guard, reasons)
        elif (command == "set_property" and raw_tokens and
              raw_tokens[0].upper() == "SOURCE"):
            self._record_source_properties(command, raw_tokens, guard, reasons)
        elif (command == "set_property" and raw_tokens and
              raw_tokens[0].upper() == "TARGET"):
            self._reject_target_graph_property(command, raw_tokens, guard)
        elif command == "set_property" and raw_tokens:
            self._reject_link_property(command, raw_tokens, guard)
        elif command == "set_target_properties" and raw_tokens:
            self._reject_target_graph_property(command, raw_tokens, guard)
        elif command == "set_directory_properties" and raw_tokens:
            self._reject_link_property(command, raw_tokens, guard)

    def _record_source_properties(self, command, raw_tokens, guard, reasons):
        if reasons:
            raise ContractError(
                "unsupported conditional CMake source properties: " + reasons[0])
        tokens = self._expand(raw_tokens, guard)
        upper = [token.upper() for token in tokens]
        if "PROPERTIES" in upper:
            property_index = upper.index("PROPERTIES")
        elif "PROPERTY" in upper:
            property_index = upper.index("PROPERTY")
        else:
            raise ContractError(
                "unsupported CMake source-property command")
        first_source = 1 if command == "set_property" else 0
        source_tokens = tokens[first_source:property_index]
        if any(token.upper() in {
                "APPEND", "APPEND_STRING", "DIRECTORY", "TARGET_DIRECTORY"}
                for token in source_tokens):
            raise ContractError(
                "unsupported scoped CMake source-property command")
        if not source_tokens:
            raise ContractError("CMake source-property command has no sources")
        properties = tokens[property_index + 1:]
        if not properties:
            raise ContractError("CMake source-property command has no properties")
        for token in source_tokens:
            reference = self._literal_source(token, reject_cuda=False)
            self.source_property_references.append(
                (reference, guard, list(properties)))

    def _reject_target_graph_property(self, command, raw_tokens, guard):
        tokens = self._expand(raw_tokens, guard)
        if (command, tuple(tokens)) in self._approved_target_property_commands:
            return
        upper = [token.upper() for token in tokens]
        marker = "PROPERTY" if command == "set_property" else "PROPERTIES"
        if marker not in upper:
            raise ContractError(
                "unsupported CMake target-property command")
        property_index = upper.index(marker)
        property_tokens = upper[property_index + 1:]
        if "SOURCES" in property_tokens:
            raise ContractError(
                "CMake target SOURCES property bypasses graph validation")
        if any(self._is_dangerous_build_property(property_name)
               for property_name in property_tokens):
            raise ContractError(
                "CMake target compile/link property bypasses graph validation")
        raise ContractError(
            "unapproved CMake target property mutation: " + command)

    def _reject_link_property(self, command, raw_tokens, guard):
        tokens = self._expand(raw_tokens, guard)
        upper = [token.upper() for token in tokens]
        if ((command, tuple(tokens)) in
                self._approved_nontarget_property_commands):
            return
        if tokens and upper[0] == "CACHE" and "PROPERTY" in upper:
            property_index = upper.index("PROPERTY")
            protected = [
                token for token in tokens[1:property_index]
                if (self._is_protected_variable(token) or
                    token == "OpenMP_CXX_FLAGS")
            ]
            if protected:
                raise ContractError(
                    "CMake cache property may mutate production "
                    "compiler-control variable: " + ", ".join(protected))
        if any(self._is_dangerous_build_property(token) for token in tokens):
            raise ContractError(
                "CMake compile/link property bypasses graph validation: " +
                command)
        raise ContractError(
            "unapproved CMake non-target property mutation: " + command)

    def _validate_production_mutations(self, production_presence):
        observed = Counter()
        for target, command, specification, mutation_guard, reasons in (
                self.target_build_mutations):
            canonical = self._canonical_target(target)
            presence = production_presence.get(canonical, BOOL_FALSE)
            overlap = bool_and(presence, mutation_guard)
            if not bool_satisfiable(overlap):
                continue
            key = (canonical, command, specification)
            if key not in self._approved_production_mutations:
                raise ContractError(
                    "unapproved production target compile/link mutation: " +
                    canonical + " " + command + " " +
                    " ".join(specification))
            allowed_reasons = ()
            if specification in {
                    ("PRIVATE", "LEO2_BACKEND_FORCE_SCALAR=1"),
                    ("PRIVATE", "LEO2_BACKEND_FORCE_SSSE3=1"),
                    ("PRIVATE", "LEO2_BACKEND_FORCE_AVX2=1"),
                    ("PRIVATE", "LEO2_BACKEND_FORCE_AVX512=1")}:
                allowed_reasons = (
                    "unsupported comparison of symbolic CMake boolean",)
            if (command == "target_link_libraries" and specification in {
                    ("PUBLIC", "OpenMP::OpenMP_CXX"),
                    ("PUBLIC", "${OpenMP_CXX_FLAGS}")}):
                allowed_reasons = (
                    "unsupported CMake conditional predicate: TARGET",)
            if reasons != allowed_reasons:
                raise ContractError(
                    "unsupported condition guards production target mutation: " +
                    canonical + " " + command + " " + repr(reasons))
            expected_guard = self._expected_production_mutation_guard(key)
            if not self._formula_equivalent(mutation_guard, expected_guard):
                raise ContractError(
                    "production target mutation guard drift: " + canonical +
                    " " + command + " actual=" + repr(mutation_guard) +
                    " expected=" + repr(expected_guard))
            observed[key] += 1
            if bool_satisfiable(
                    bool_and(mutation_guard, bool_not(presence))):
                raise ContractError(
                    "production target mutation is reachable without its "
                    "target: " + canonical + " " + command)
            if (key == ("leopard", "target_link_libraries",
                        ("PUBLIC", "${OpenMP_CXX_FLAGS}")) and
                    "OpenMP_CXX_FLAGS" in self.variables):
                raise ContractError(
                    "local OpenMP_CXX_FLAGS mutation can redirect leopard "
                    "link inputs")
        if self.require_mutation_contract:
            expected = Counter(self._approved_production_mutations)
            if observed != expected:
                missing = sorted((expected - observed).elements(), key=repr)
                extra = sorted((observed - expected).elements(), key=repr)
                raise ContractError(
                    "missing or duplicate production target mutation: " +
                    "missing=" + repr(missing) + " extra=" + repr(extra))

    @staticmethod
    def _literal_source(token, reject_cuda=True):
        if "$<" in token:
            raise ContractError(
                "unsupported CMake source generator expression: " + token)
        if "$" in token:
            raise ContractError("unresolved CMake source token: " + token)
        normalized = token.replace("\\", "/")
        if re.match(r"^[A-Za-z]:", normalized):
            raise ContractError(
                "production source must not use a drive path: " + token)
        path = PurePosixPath(normalized)
        if path.is_absolute() or ".." in path.parts:
            raise ContractError(
                "production source must be repository-relative: " + token)
        suffix = path.suffix.lower()
        if suffix not in KNOWN_SOURCE_SUFFIXES:
            raise ContractError(
                "unsupported production source path (fail closed): " + token)
        if reject_cuda and suffix in CUDA_SUFFIXES:
            raise ContractError(
                "CUDA source attached to ordinary CPU leopard: " + token)
        return path.as_posix()

    def resolve(self, target="leopard"):
        resolved = []
        attached_objects = set()
        object_sources = set()
        production_presence = {}
        visiting = []

        def visit(name, reach_guard, reached_as_object=False):
            name = self._canonical_target(name)
            for active_name, active_guard in visiting:
                if (name == active_name and
                        bool_satisfiable(bool_and(reach_guard, active_guard))):
                    raise ContractError(
                        "cyclic CMake TARGET_OBJECTS attachment: " + name)
            definitions = self.targets.get(name, [])
            reachable = [definition for definition in definitions
                         if bool_satisfiable(bool_and(
                             reach_guard, definition[2]))]
            if not reachable:
                if name in self.declared_targets:
                    raise ContractError(
                        "attached OBJECT target has no MSVC-reachable "
                        "definition: " + name)
                raise ContractError(
                    "attached CMake target is undefined: " + name)
            visiting.append((name, reach_guard))
            presence = BOOL_FALSE
            entries = []
            for kind, sources, definition_guard in reachable:
                active = bool_and(reach_guard, definition_guard)
                presence = bool_or(presence, active)
                entries.extend((token, active) for token in sources)
            entries.extend(
                (token, bool_and(presence, attachment_guard))
                for token, attachment_guard in self.attachments.get(name, []))
            production_presence[name] = bool_or(
                production_presence.get(name, BOOL_FALSE),
                bool_and(reach_guard, presence))
            for token, entry_guard in entries:
                if not bool_satisfiable(entry_guard):
                    continue
                object_match = self._target_objects.match(token)
                if object_match:
                    object_target = self._canonical_target(
                        object_match.group(1))
                    object_definitions = [
                        definition for definition in self.targets.get(
                            object_target, [])
                        if bool_satisfiable(bool_and(
                            entry_guard, definition[2]))]
                    if not object_definitions:
                        if object_target in self.declared_targets:
                            raise ContractError(
                                "attached OBJECT target has no MSVC-reachable "
                                "definition: " + object_target)
                        raise ContractError(
                            "TARGET_OBJECTS target is undefined: " + object_target)
                    if any(definition[0] != "OBJECT"
                           for definition in object_definitions):
                        raise ContractError(
                            "TARGET_OBJECTS does not name an OBJECT library: " +
                            object_target)
                    attached_objects.add(object_target)
                    visit(object_target, entry_guard, True)
                else:
                    literal = self._literal_source(token)
                    resolved.append((literal, entry_guard))
                    if reached_as_object:
                        object_sources.add(literal)
            visiting.pop()

        visit(target, BOOL_TRUE)
        self._validate_production_mutations(production_presence)
        if (self.require_mutation_contract and
                tuple(self.contract_events) !=
                self._required_contract_event_order):
            raise ContractError(
                "security-sensitive CMake command order drift")
        resolved_by_path = {}
        resolved_by_windows_key = {}
        for path, guard in resolved:
            resolved_by_path.setdefault(path, []).append(guard)
            resolved_by_windows_key.setdefault(path.casefold(), []).append(
                (path, guard))
        for reference, property_guard, properties in self.source_property_references:
            if any(bool_satisfiable(bool_and(property_guard, source_guard))
                   for unused_path, source_guard in
                   resolved_by_windows_key.get(reference.casefold(), [])):
                raise ContractError(
                    "CMake source properties affect production " + reference +
                    ": " + " ".join(properties))
        duplicates = []
        for entries in resolved_by_windows_key.values():
            guards = [guard for unused_path, guard in entries]
            if any(bool_satisfiable(bool_and(left, right))
                   for index, left in enumerate(guards)
                   for right in guards[index + 1:]):
                duplicates.append(" / ".join(sorted({
                    path for path, unused_guard in entries})))
        duplicates.sort()
        if duplicates:
            raise ContractError(
                "duplicate production source attachment: " + ", ".join(duplicates))
        return (sorted(resolved_by_path), sorted(attached_objects),
                sorted(object_sources))


def project_path(project_file, value):
    native = value.replace("\\", "/")
    return (project_file.parent / native).resolve().relative_to(ROOT).as_posix()


def item_paths(tree, kind, project_file):
    return [
        project_path(project_file, node.attrib["Include"])
        for node in tree.findall(".//msb:" + kind + "[@Include]", NS)
    ]


def xml_local_name(node):
    return node.tag.rsplit("}", 1)[-1]


def validate_msbuild_namespace(tree):
    namespace_prefix = "{" + NS["msb"] + "}"
    for node in tree.iter():
        if (not isinstance(node.tag, str) or
                not node.tag.startswith(namespace_prefix) or
                node.tag == namespace_prefix):
            raise ContractError(
                "foreign or missing MSBuild XML namespace: " +
                str(node.tag))


def validate_msbuild_element_casing(tree):
    for node in tree.iter():
        name = xml_local_name(node)
        canonical = MSBUILD_CANONICAL_NAMES.get(name.lower())
        if canonical is not None and name != canonical:
            raise ContractError(
                "noncanonical MSBuild element casing: " + name +
                " (expected " + canonical + ")")


def validate_source_item_structure(tree):
    allowed_compile_metadata = {"AdditionalOptions"}
    direct_item_groups = set(tree.getroot().findall("msb:ItemGroup", NS))
    seen = []
    for group in tree.findall(".//msb:ItemGroup", NS):
        source_items = [
            node for node in list(group)
            if xml_local_name(node) in ("ClCompile", "ClInclude")
        ]
        if not source_items:
            continue
        if group not in direct_item_groups:
            raise ContractError(
                "source ItemGroup must be a direct Project child")
        if "Condition" in group.attrib:
            raise ContractError("source ItemGroup must be unconditional")
        for node in source_items:
            kind = xml_local_name(node)
            if "Remove" in node.attrib or "Update" in node.attrib:
                raise ContractError(kind + " Remove/Update transforms are unsupported")
            if "Include" not in node.attrib:
                raise ContractError(kind + " source item lacks Include")
            if "Condition" in node.attrib:
                raise ContractError(kind + " source item must be unconditional")
            unexpected_attributes = set(node.attrib) - {"Include"}
            if unexpected_attributes:
                raise ContractError(
                    kind + " source item has unsupported attributes: " +
                    ", ".join(sorted(unexpected_attributes)))
            path = project_path(PROJECT, node.attrib["Include"])
            seen.append((kind, path))
            for metadata in list(node):
                name = xml_local_name(metadata)
                if name == "ExcludedFromBuild":
                    raise ContractError(
                        "ExcludedFromBuild can invalidate required source " +
                        "coverage: " + path)
                if "Condition" in metadata.attrib:
                    raise ContractError(
                        path + " source metadata must be unconditional")
                if metadata.attrib:
                    raise ContractError(
                        path + " source metadata has unsupported attributes")
                if kind == "ClInclude" or name not in allowed_compile_metadata:
                    raise ContractError(
                        path + " has unsupported source metadata: " + name)

    for node in tree.findall(".//msb:ExcludedFromBuild", NS):
        raise ContractError(
            "ExcludedFromBuild can invalidate required source coverage")

    duplicates = sorted(
        item for item, count in Counter(seen).items() if count != 1)
    if duplicates:
        raise ContractError("duplicate MSBuild source items: " + repr(duplicates))

    all_item_groups = tree.findall(".//msb:ItemGroup", NS)
    root_item_groups = tree.getroot().findall("msb:ItemGroup", NS)
    if set(all_item_groups) != set(root_item_groups):
        raise ContractError("all MSBuild ItemGroups must be direct Project children")
    if len(root_item_groups) != 3:
        raise ContractError("MSBuild ItemGroup topology is not exact")
    item_group_kinds = []
    for group in root_item_groups:
        child_names = {xml_local_name(child) for child in list(group)}
        if group.attrib == {"Label": "ProjectConfigurations"}:
            expected = {"ProjectConfiguration"}
            kind = "ProjectConfigurations"
        elif not group.attrib and child_names == {"ClInclude"}:
            expected = {"ClInclude"}
            kind = "ClInclude"
        elif not group.attrib and child_names == {"ClCompile"}:
            expected = {"ClCompile"}
            kind = "ClCompile"
        else:
            raise ContractError(
                "unapproved MSBuild ItemGroup contents: " +
                ", ".join(sorted(child_names)))
        if child_names != expected:
            raise ContractError("MSBuild " + kind + " ItemGroup is not exact")
        item_group_kinds.append(kind)
    if item_group_kinds != ["ProjectConfigurations", "ClInclude", "ClCompile"]:
        raise ContractError("MSBuild ItemGroup order is not exact")


def validate_exact_msbuild_children(parent, expected, label):
    children = list(parent)
    if [xml_local_name(child) for child in children] != [
            name for name, unused_text in expected]:
        raise ContractError(label + " children differ from the exact contract")
    for child, (name, expected_text) in zip(children, expected):
        if child.attrib or list(child):
            raise ContractError(label + " " + name + " must be a plain value")
        if (child.text or "").strip() != expected_text:
            raise ContractError(label + " " + name + " value differs")


def validate_msbuild_property_topology(tree):
    root = tree.getroot()
    all_groups = tree.findall(".//msb:PropertyGroup", NS)
    root_groups = root.findall("msb:PropertyGroup", NS)
    if set(all_groups) != set(root_groups):
        raise ContractError(
            "all MSBuild PropertyGroups must be direct Project children")
    if len(root_groups) != 10:
        raise ContractError("MSBuild PropertyGroup topology is not exact")

    globals_groups = [group for group in root_groups
                      if group.attrib == {"Label": "Globals"}]
    if len(globals_groups) != 1:
        raise ContractError("MSBuild Globals PropertyGroup is not exact")
    validate_exact_msbuild_children(globals_groups[0], (
        ("ProjectGuid", "{32176592-2F30-4BD5-B645-EB11C8D3453E}"),
        ("RootNamespace", "GF65536"),
        ("ProjectName", "Leopard"),
    ), "Globals PropertyGroup")

    user_groups = [group for group in root_groups
                   if group.attrib == {"Label": "UserMacros"}]
    if len(user_groups) != 1 or list(user_groups[0]):
        raise ContractError("MSBuild UserMacros PropertyGroup must be empty")

    configuration_groups = [
        group for group in root_groups
        if group.attrib.get("Label") == "Configuration"]
    configurations = index_conditions(
        configuration_groups, "Configuration PropertyGroup",
        {"Condition", "Label"})
    for key, group in configurations.items():
        configuration, unused_platform = EXPECTED_CONFIGS[key]
        del unused_platform
        validate_exact_msbuild_children(group, (
            ("ConfigurationType", "StaticLibrary"),
            ("UseDebugLibraries",
             "true" if configuration == "Debug" else "false"),
            ("WholeProgramOptimization", "false"),
            ("CharacterSet", "MultiByte"),
            ("PlatformToolset", "v140"),
        ), key + " Configuration PropertyGroup")

    output_groups = [group for group in root_groups if not group.attrib.get("Label")]
    outputs = index_conditions(output_groups, "output PropertyGroup")
    for key, group in outputs.items():
        validate_exact_msbuild_children(group, (
            ("OutDir", "Output/$(ProjectName)/$(Configuration)/$(Platform)/"),
            ("IntDir", "Obj/$(ProjectName)/$(Configuration)/$(Platform)/"),
        ), key + " output PropertyGroup")

    approved = set(globals_groups + user_groups + configuration_groups +
                   output_groups)
    if set(root_groups) != approved:
        raise ContractError("unapproved root MSBuild PropertyGroup")


def validate_msbuild_root_order(tree):
    """Pin evaluation phases, not merely the contents of each phase."""
    def descriptor(name, attributes=None):
        return (name, tuple(sorted((attributes or {}).items())))

    def condition(configuration, platform):
        return "'$(Configuration)|$(Platform)'=='" + \
            configuration + "|" + platform + "'"

    expected = [
        descriptor("ItemGroup", {"Label": "ProjectConfigurations"}),
        descriptor("ItemGroup"),
        descriptor("ItemGroup"),
        descriptor("PropertyGroup", {"Label": "Globals"}),
        descriptor("Import", {
            "Project": "$(VCTargetsPath)\\Microsoft.Cpp.Default.props"}),
    ]
    for configuration, platform in (
            ("Debug", "Win32"), ("Debug", "x64"),
            ("Release", "Win32"), ("Release", "x64")):
        expected.append(descriptor("PropertyGroup", {
            "Condition": condition(configuration, platform),
            "Label": "Configuration"}))
    expected.extend((
        descriptor("Import", {
            "Project": "$(VCTargetsPath)\\Microsoft.Cpp.props"}),
        descriptor("ImportGroup", {"Label": "ExtensionSettings"}),
    ))
    for configuration, platform in (
            ("Debug", "Win32"), ("Debug", "x64"),
            ("Release", "Win32"), ("Release", "x64")):
        expected.append(descriptor("ImportGroup", {
            "Condition": condition(configuration, platform),
            "Label": "PropertySheets"}))
    expected.append(descriptor("PropertyGroup", {"Label": "UserMacros"}))
    for configuration, platform in (
            ("Debug", "Win32"), ("Release", "Win32"),
            ("Debug", "x64"), ("Release", "x64")):
        expected.append(descriptor("PropertyGroup", {
            "Condition": condition(configuration, platform)}))
    for configuration, platform in (
            ("Debug", "Win32"), ("Debug", "x64"),
            ("Release", "Win32"), ("Release", "x64")):
        expected.append(descriptor("ItemDefinitionGroup", {
            "Condition": condition(configuration, platform)}))
    expected.extend((
        descriptor("Import", {
            "Project": "$(VCTargetsPath)\\Microsoft.Cpp.targets"}),
        descriptor("ImportGroup", {"Label": "ExtensionTargets"}),
    ))

    actual = [descriptor(xml_local_name(node), node.attrib)
              for node in list(tree.getroot())]
    if actual != expected:
        raise ContractError(
            "MSBuild root evaluation phase/order differs from the exact "
            "contract")


def validate_no_wpo_overrides(tree):
    allowed_additional_options = set()
    for path in AVX2_SOURCE_FILES:
        allowed_additional_options.update(tree.findall(
            (".//msb:ClCompile[@Include='..\\{}']/"
             "msb:AdditionalOptions").format(path), NS))
    for node in tree.iter():
        name = xml_local_name(node)
        text = node.text or ""
        if name == "AdditionalOptions":
            reject_msbuild_expansion(
                text, "AdditionalOptions", ("%(AdditionalOptions)",))
            reject_wpo_options(text, "AdditionalOptions")
            if node not in allowed_additional_options:
                raise ContractError(
                    "AdditionalOptions is allowed only on the AVX2 source")
        elif name == "WholeProgramOptimization":
            if text.strip().lower() != "false" or "Condition" in node.attrib:
                raise ContractError(
                    "WholeProgramOptimization must be unconditional false")
        elif name.lower() in {
                "linktimecodegeneration",
                "uselinktimecodegeneration",
                "profileguidedoptimization"}:
            raise ContractError(name + " can override WPO isolation")


def validate_msbuild_imports_and_toolchain(tree):
    """Keep imported build logic and compiler selection a closed contract."""
    root = tree.getroot()
    expected_root_attributes = {
        "DefaultTargets": "Build",
        "ToolsVersion": EXPECTED_TOOLS_VERSION,
    }
    if root.attrib != expected_root_attributes:
        raise ContractError(
            "MSBuild Project attributes differ from the VS2015 contract")
    allowed_root_children = {
        "Import", "ImportGroup", "ItemDefinitionGroup", "ItemGroup",
        "PropertyGroup"}
    unexpected_root_children = [
        xml_local_name(node) for node in list(root)
        if xml_local_name(node) not in allowed_root_children]
    if unexpected_root_children:
        raise ContractError(
            "unapproved root MSBuild elements: " +
            ", ".join(unexpected_root_children))

    parent = {
        child: node for node in root.iter() for child in list(node)
    }
    user_props = "$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props"
    expected_projects = [
        "$(VCTargetsPath)\\Microsoft.Cpp.Default.props",
        "$(VCTargetsPath)\\Microsoft.Cpp.props",
        user_props, user_props, user_props, user_props,
        "$(VCTargetsPath)\\Microsoft.Cpp.targets",
    ]
    imports = list(root.iter(
        "{http://schemas.microsoft.com/developer/msbuild/2003}Import"))
    if [node.attrib.get("Project") for node in imports] != expected_projects:
        raise ContractError("MSBuild imports differ from the approved toolchain")

    direct_imports = [node for node in imports if parent.get(node) is root]
    if [node.attrib.get("Project") for node in direct_imports] != [
            expected_projects[0], expected_projects[1], expected_projects[-1]]:
        raise ContractError("Microsoft C++ imports are not direct Project children")
    for node in direct_imports:
        if node.attrib != {"Project": node.attrib.get("Project")}:
            raise ContractError("Microsoft C++ import has unsupported attributes")

    property_sheets = index_conditions(
        root.findall("msb:ImportGroup[@Label='PropertySheets']", NS),
        "PropertySheets ImportGroup", allowed_attributes={"Condition", "Label"})
    user_condition = (
        "exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')")
    for key, group in property_sheets.items():
        if group.attrib.get("Label") != "PropertySheets":
            raise ContractError(key + " property-sheet group has the wrong label")
        children = list(group)
        if len(children) != 1 or xml_local_name(children[0]) != "Import":
            raise ContractError(key + " property-sheet import is not exact")
        imported = children[0]
        if imported.attrib != {
                "Project": user_props,
                "Condition": user_condition,
                "Label": "LocalAppDataPlatform"}:
            raise ContractError(key + " property-sheet import is not approved")

    import_groups = root.findall("msb:ImportGroup", NS)
    extension_groups = [group for group in import_groups
                        if group.attrib.get("Label") in {
                            "ExtensionSettings", "ExtensionTargets"}]
    if len(import_groups) != 6 or len(extension_groups) != 2:
        raise ContractError("MSBuild ImportGroup topology is not the exact contract")
    for label in ("ExtensionSettings", "ExtensionTargets"):
        groups = [group for group in extension_groups
                  if group.attrib.get("Label") == label]
        if (len(groups) != 1 or groups[0].attrib != {"Label": label} or
                list(groups[0])):
            raise ContractError(label + " must be an empty approved ImportGroup")

    forbidden_properties = {
        "cltoolexe", "cltoolpath", "vctoolspath", "vctoolsinstalldir",
        "vcinstalldir", "vctargetspath", "executablepath", "useenv",
        "toolexe", "toolpath", "userrootdir", "msbuilduserextensionspath",
        "msbuildextensionspath", "msbuildextensionspath32",
        "msbuildextensionspath64", "includepath", "librarypath",
        "referencepath", "sourcepath", "excludepath", "cltoolarchitecture",
    }
    for node in root.iter():
        name = xml_local_name(node).lower()
        parent_name = xml_local_name(parent[node]).lower() if node in parent else ""
        if name in {"configuration", "platform"}:
            if parent_name != "projectconfiguration":
                raise ContractError(
                    "MSBuild configuration/import input override is forbidden: " +
                    xml_local_name(node))
        if (name.endswith("dependson") or
                (name.endswith("targets") and name.startswith((
                    "forceimport", "custombefore", "customafter",
                    "importbywildcard")))):
            raise ContractError(
                "MSBuild import/target-chain hook is forbidden: " +
                xml_local_name(node))
        if name in forbidden_properties:
            raise ContractError(
                "MSBuild compiler tool override is forbidden: " +
                xml_local_name(node))
        if name in {"target", "usingtask", "exec"}:
            raise ContractError(
                "custom MSBuild execution logic is forbidden: " +
                xml_local_name(node))
        if name in {"custombuild", "custombuildstep"}:
            raise ContractError(
                "custom MSBuild build item is forbidden: " +
                xml_local_name(node))
        if name in {"lib", "additionaldependencies"}:
            raise ContractError(
                "unverified MSBuild librarian archive input is forbidden: " +
                xml_local_name(node))
        if name == "command" and ((node.text or "").strip() or node.attrib):
            raise ContractError("nonempty MSBuild build-event command is forbidden")
        for attribute in node.attrib:
            if attribute.rsplit("}", 1)[-1].lower() in (
                    forbidden_properties | {"toolexe", "toolpath"}):
                raise ContractError(
                    "MSBuild task tool override is forbidden: " + attribute)


def strip_c_family_comments(source_text):
    """Remove C/C++ comments while preserving literals and line numbers."""
    result = []
    index = 0
    length = len(source_text)
    while index < length:
        if source_text.startswith("//", index):
            result.extend((" ", " "))
            index += 2
            while index < length and source_text[index] not in "\r\n":
                result.append(" ")
                index += 1
            continue
        if source_text.startswith("/*", index):
            result.extend((" ", " "))
            index += 2
            while index < length and not source_text.startswith("*/", index):
                character = source_text[index]
                result.append(character if character in "\r\n" else " ")
                index += 1
            if index == length:
                raise ContractError("unterminated C/C++ block comment")
            result.extend((" ", " "))
            index += 2
            continue

        raw = None
        if source_text[index] in "RuUL":
            raw = re.match(
                r'(?:u8|u|U|L)?R"([^ ()\\\t\r\n]{0,16})\(',
                source_text[index:])
        if raw is not None:
            terminator = ")" + raw.group(1) + '"'
            end = source_text.find(
                terminator, index + raw.end())
            if end < 0:
                raise ContractError("unterminated C++ raw string literal")
            end += len(terminator)
            result.append(source_text[index:end])
            index = end
            continue

        character = source_text[index]
        if character in "\"'":
            quote = character
            result.append(character)
            index += 1
            while index < length:
                character = source_text[index]
                result.append(character)
                index += 1
                if character == "\\" and index < length:
                    result.append(source_text[index])
                    index += 1
                elif character == quote:
                    break
            continue
        result.append(character)
        index += 1
    return "".join(result)


def production_local_includes(source_text):
    source_text = strip_c_family_comments(
        re.sub(r"\\\r?\n", "", source_text))
    includes = []
    for operand in re.findall(
            r"^\s*#\s*include\b([^\n]*)", source_text, re.MULTILINE):
        literal = re.match(
            r'^\s*(?:<([^>\n]+)>|"([^"\n]+)")', operand)
        if literal is None:
            raise ContractError(
                "nonliteral production include cannot be proved CPU-only: " +
                operand.strip())
        includes.append(literal.group(1) or literal.group(2))
    return includes


def production_cuda_preprocessor_lines(source_text):
    source_text = strip_c_family_comments(
        re.sub(r"\\\r?\n", "", source_text))
    lines = [line for line in re.findall(
        r"^\s*#.*$", source_text, re.MULTILINE)
        if (re.search(r"\.cuh?(?:[>\"']|\s|$)", line, re.IGNORECASE) or
            re.search(
                r"(?:cuda(?:_[A-Za-z0-9_]+)?|nvrtc)\.h"
                r"(?:[>\"']|\s|$)",
                line, re.IGNORECASE) or
            re.search(
                r"cuda|nvrtc|nvidia|[<\"'/\\](?:thrust|cub)[/\\]",
                line, re.IGNORECASE))]
    lines.extend(
        line for line in source_text.splitlines()
        if (not re.match(r"^\s*#", line) and
            re.search(r"(?:__pragma|_Pragma)\s*\(.*comment", line,
                      re.IGNORECASE) and
            re.search(r"cuda|nvrtc|nvidia", line, re.IGNORECASE)))
    return lines


def production_cuda_source_lines(source_text):
    """Return any CPU-production lines that mention a CUDA runtime/toolchain."""
    marker = re.compile(
        r"cuda|nvcc|nvrtc|nvidia|cudart|ptxas|nvlink|fatbinary|"
        r"cuobjdump|nvc\+\+|__(?:global|device|host|shared|constant|managed)__|"
        r"__launch_bounds__|<<<|>>>|\.cuh?(?:[>\"']|\s|$)|"
        r"[<\"'/\\](?:thrust|cub)[/\\]",
        re.IGNORECASE)
    code = strip_c_family_comments(re.sub(r"\\\r?\n", "", source_text))
    return [line for line in code.splitlines() if marker.search(line)]


def windows_repository_include(base_directory, included):
    requested = PurePosixPath(included.replace("\\", "/"))
    if requested.is_absolute():
        return None
    try:
        base_parts = list(
            base_directory.resolve().relative_to(ROOT.resolve()).parts)
    except ValueError:
        raise ContractError("include base escapes the repository")
    parts = base_parts
    for part in requested.parts:
        if part in {"", "."}:
            continue
        if part == "..":
            if not parts:
                return None
            parts.pop()
        else:
            parts.append(part)
    current = ROOT.resolve()
    for part in parts:
        if not current.is_dir():
            return None
        matches = [child for child in current.iterdir()
                   if child.name.casefold() == part.casefold()]
        if len(matches) > 1:
            raise ContractError(
                "ambiguous Windows case-insensitive include component: " +
                part)
        if not matches:
            return None
        current = matches[0]
    try:
        current.resolve().relative_to(ROOT.resolve())
    except ValueError:
        raise ContractError("included file escapes the repository")
    return current


def production_graph(text=None, require_files=True,
                     require_mutation_contract=None):
    if require_mutation_contract is None:
        require_mutation_contract = text is None
    cmake = CMAKE.read_text(encoding="utf-8") if text is None else text
    processor_configurations = (
        ("x86", "4", "Win32"),
        ("X86", "4", "Win32"),
        ("i386", "4", "Win32"),
        ("i486", "4", "Win32"),
        ("i586", "4", "Win32"),
        ("i686", "4", "Win32"),
        ("x86_64", "8", "x64"),
        ("amd64", "8", "x64"),
        ("AMD64", "8", "x64"),
    )
    configurations = tuple(
        CMakeProductionGraph(
            cmake, processor=processor, pointer_size=pointer_size,
            platform_name=platform_name,
            require_mutation_contract=require_mutation_contract).resolve()
        for processor, pointer_size, platform_name in
        processor_configurations)
    if any(configuration != configurations[0]
           for configuration in configurations[1:]):
        raise ContractError(
            "CMake production graph differs between Win32 and x64")
    attached, object_targets, object_sources = configurations[0]
    compiled = [path for path in attached if path.endswith(COMPILE_SUFFIXES)]
    headers = {path for path in attached if path.endswith(HEADER_SUFFIXES)}

    if require_files:
        pending = list(compiled) + list(headers)
        visited = set()
        while pending:
            relative = pending.pop()
            if relative in visited:
                continue
            visited.add(relative)
            source = ROOT / relative
            if not source.is_file():
                raise ContractError(
                    "production graph names missing file: " + relative)
            source_text = source.read_text(encoding="utf-8")
            cuda_lines = production_cuda_source_lines(source_text)
            if cuda_lines:
                raise ContractError(
                    "CUDA source marker reachable from CPU leopard: " +
                    cuda_lines[0].strip())
            for included in production_local_includes(source_text):
                if re.match(
                        r"^(?:cuda(?:_[a-z0-9_]+)?|nvrtc)\.h$",
                        PurePosixPath(
                            included.replace("\\", "/")).name.lower()):
                    raise ContractError(
                        "CUDA header reachable from CPU leopard: " +
                        included)
                candidate = windows_repository_include(
                    source.parent, included)
                if candidate is None:
                    continue
                local = candidate.resolve().relative_to(
                    ROOT.resolve()).as_posix()
                if candidate.suffix.lower() in CUDA_SUFFIXES:
                    raise ContractError(
                        "CUDA header reachable from CPU leopard: " + local)
                if (candidate.is_file() and
                        candidate.suffix.lower() in HEADER_SUFFIXES):
                    headers.add(local)
                if candidate.is_file() and local not in visited:
                    pending.append(local)
    return (sorted(compiled), sorted(headers), object_targets,
            object_sources, cmake)


def normalized_condition(condition):
    match = re.match(
        r"^\s*'\$\(Configuration\)\|\$\(Platform\)'\s*==\s*"
        r"'([^'|]+)\|([^']+)'\s*$", condition or "")
    if not match:
        raise ContractError("unsupported MSBuild Condition: " + str(condition))
    configuration = match.group(1)
    platform = match.group(2)
    if configuration != configuration.strip() or platform != platform.strip():
        raise ContractError(
            "MSBuild Condition values contain semantic whitespace: " +
            str(condition))
    key = configuration.lower() + "|" + platform.lower()
    if key not in EXPECTED_CONFIGS:
        raise ContractError("unexpected MSBuild configuration: " + key)
    return key


def index_conditions(nodes, label, allowed_attributes=None):
    allowed_attributes = ({"Condition"} if allowed_attributes is None else
                          set(allowed_attributes))
    indexed = {}
    for node in nodes:
        unexpected = set(node.attrib) - allowed_attributes
        if unexpected:
            raise ContractError(
                label + " has unsupported attributes: " +
                ", ".join(sorted(unexpected)))
        key = normalized_condition(node.attrib.get("Condition"))
        if key in indexed:
            raise ContractError("duplicate " + label + " Condition: " + key)
        indexed[key] = node
    if set(indexed) != set(EXPECTED_CONFIGS):
        raise ContractError(
            label + " Conditions differ: " + ", ".join(sorted(indexed)))
    return indexed


def validate_msbuild_configurations(tree):
    root = tree.getroot()
    project_groups = tree.findall(
        ".//msb:ItemGroup[@Label='ProjectConfigurations']", NS)
    direct_project_groups = root.findall(
        "msb:ItemGroup[@Label='ProjectConfigurations']", NS)
    if (len(project_groups) != 1 or project_groups != direct_project_groups or
            project_groups[0].attrib != {"Label": "ProjectConfigurations"}):
        raise ContractError(
            "ProjectConfigurations must be one unconditional direct ItemGroup")
    project_group = project_groups[0]
    configurations = list(project_group)
    if any(xml_local_name(node) != "ProjectConfiguration"
           for node in configurations):
        raise ContractError(
            "ProjectConfigurations contains an unsupported child")
    all_configurations = [
        node for node in root.iter()
        if xml_local_name(node).lower() == "projectconfiguration"]
    if (len(all_configurations) != len(configurations) or
            set(all_configurations) != set(configurations)):
        raise ContractError(
            "ProjectConfiguration items outside the canonical group are "
            "forbidden")
    includes = [node.attrib.get("Include", "").lower()
                for node in configurations]
    if len(includes) != len(EXPECTED_CONFIGS) or set(includes) != set(
            EXPECTED_CONFIGS):
        raise ContractError("ProjectConfigurations are not the exact four configs")
    for node in configurations:
        key = node.attrib["Include"].lower()
        expected_configuration, expected_platform = EXPECTED_CONFIGS[key]
        expected_include = expected_configuration + "|" + expected_platform
        if node.attrib != {"Include": expected_include}:
            raise ContractError(
                key + " ProjectConfiguration attributes are not exact")
        children = list(node)
        if ([xml_local_name(child) for child in children] !=
                ["Configuration", "Platform"] or
                any(child.attrib for child in children)):
            raise ContractError(
                key + " ProjectConfiguration children are not exact")
        actual_configuration = node.findtext(
            "msb:Configuration", namespaces=NS)
        if actual_configuration != expected_configuration:
            raise ContractError(key + " has a mismatched Configuration value")
        if node.findtext("msb:Platform", namespaces=NS) != expected_platform:
            raise ContractError(key + " has a mismatched Platform value")

    property_nodes = tree.findall(
        ".//msb:PropertyGroup[@Label='Configuration']", NS)
    definition_nodes = tree.findall(".//msb:ItemDefinitionGroup", NS)
    if (set(property_nodes) != set(root.findall(
            "msb:PropertyGroup[@Label='Configuration']", NS)) or
            set(definition_nodes) != set(root.findall(
                "msb:ItemDefinitionGroup", NS))):
        raise ContractError(
            "configuration groups must be direct Project children")
    properties = index_conditions(property_nodes,
        "Configuration PropertyGroup", {"Condition", "Label"})
    definitions = index_conditions(
        definition_nodes, "ItemDefinitionGroup")

    root_children = list(root)
    default_props_index = next(
        index for index, node in enumerate(root_children)
        if (xml_local_name(node) == "Import" and
            node.attrib.get("Project", "").endswith(
                "Microsoft.Cpp.Default.props")))
    cpp_props_index = next(
        index for index, node in enumerate(root_children)
        if (xml_local_name(node) == "Import" and
            node.attrib.get("Project", "").endswith("Microsoft.Cpp.props")))
    if not all(default_props_index < root_children.index(node) < cpp_props_index
               for node in property_nodes):
        raise ContractError(
            "configuration groups must precede Microsoft.Cpp.props")

    approved_toolsets = set()
    approved_configuration_types = set()
    for key in sorted(EXPECTED_CONFIGS):
        configuration, unused_platform = EXPECTED_CONFIGS[key]
        del unused_platform
        prop = properties[key]
        if [xml_local_name(node) for node in list(definitions[key])] != [
                "ClCompile", "Link", "PostBuildEvent"]:
            raise ContractError(
                key + " ItemDefinitionGroup tool topology is not exact")

        def unique_text(parent, tag, label):
            nodes = parent.findall("msb:" + tag, NS)
            if len(nodes) != 1:
                raise ContractError(key + " must define exactly one " + label)
            if nodes[0].attrib:
                raise ContractError(
                    key + " has conditional/attributed " + label)
            return nodes[0].text or ""

        if unique_text(
                prop, "ConfigurationType", "ConfigurationType") != "StaticLibrary":
            raise ContractError(key + " is not a StaticLibrary")
        approved_configuration_types.update(
            prop.findall("msb:ConfigurationType", NS))
        if unique_text(prop, "PlatformToolset", "PlatformToolset") != "v140":
            raise ContractError(key + " does not preserve the v140 toolset")
        approved_toolsets.update(prop.findall("msb:PlatformToolset", NS))
        wpo = unique_text(
            prop, "WholeProgramOptimization", "WholeProgramOptimization")
        if (wpo or "").strip().lower() != "false":
            raise ContractError(key + " does not explicitly disable /GL")

        compile_nodes = definitions[key].findall("msb:ClCompile", NS)
        if len(compile_nodes) != 1:
            raise ContractError(key + " must define exactly one ClCompile")
        compile_node = compile_nodes[0]
        if compile_node.attrib:
            raise ContractError(
                key + " configuration ClCompile must be unconditional")
        expected_runtime = (
            "MultiThreadedDebug" if configuration == "Debug" else "MultiThreaded")
        if unique_text(
                compile_node, "RuntimeLibrary", "RuntimeLibrary") != expected_runtime:
            raise ContractError(key + " has the wrong static runtime library")
        raw_definitions = unique_text(
            compile_node, "PreprocessorDefinitions", "PreprocessorDefinitions")
        reject_msbuild_expansion(
            raw_definitions, key + " PreprocessorDefinitions",
            ("%(PreprocessorDefinitions)",))
        macro_set = msbuild_values(raw_definitions)
        required = set(BACKEND_DEFINITIONS)
        required.update(("_MBCS", "%(PreprocessorDefinitions)"))
        if macro_set != required:
            raise ContractError(
                key + " backend definitions are not the exact contract")
        if any(macro.startswith("LEO2_BACKEND_FORCE_") for macro in macro_set):
            raise ContractError(key + " forces a diagnostic backend")
        for forbidden in (
                "AdditionalOptions",
                "UndefinePreprocessorDefinitions",
                "EnableEnhancedInstructionSet",
                "ForcedIncludeFiles",
                "ForcedUsingFiles"):
            if compile_node.findall("msb:" + forbidden, NS):
                raise ContractError(
                    key + " must not define configuration-level " + forbidden)

        expected_compile = [
            ("WarningLevel", "Level3"),
            ("Optimization",
             "Disabled" if configuration == "Debug" else
             ("MaxSpeed" if key == "release|win32" else "Full")),
        ]
        if configuration == "Release":
            expected_compile.extend((
                ("FunctionLevelLinking", "true"),
                ("IntrinsicFunctions", "true"),
            ))
        expected_compile.append(("SDLCheck", "true"))
        if configuration == "Release":
            expected_compile.extend((
                ("InlineFunctionExpansion",
                 "AnySuitable" if key == "release|win32" else
                 "OnlyExplicitInline"),
                ("FavorSizeOrSpeed", "Speed"),
                ("OmitFramePointers",
                 "false" if key == "release|win32" else "true"),
            ))
        expected_compile.append(("RuntimeLibrary", expected_runtime))
        if configuration == "Release":
            expected_compile.append(("BufferSecurityCheck", "true"))
        expected_compile.append(("PreprocessorDefinitions", raw_definitions))
        if key == "release|x64":
            expected_compile.append(("OpenMPSupport", "true"))
        validate_exact_msbuild_children(
            compile_node, tuple(expected_compile),
            key + " configuration ClCompile")

        link_node = definitions[key].find("msb:Link", NS)
        expected_link = [("GenerateDebugInformation", "true")]
        if configuration == "Release":
            expected_link.extend((
                ("EnableCOMDATFolding", "true"),
                ("OptimizeReferences", "true"),
            ))
        expected_link.append(("AdditionalLibraryDirectories", ""))
        validate_exact_msbuild_children(
            link_node, tuple(expected_link), key + " configuration Link")
        validate_exact_msbuild_children(
            definitions[key].find("msb:PostBuildEvent", NS),
            (("Command", ""),), key + " PostBuildEvent")

    for node in tree.findall(".//msb:WholeProgramOptimization", NS):
        if (node.text or "").strip().lower() != "false":
            raise ContractError("WholeProgramOptimization must always be false")
    all_toolsets = {
        node for node in root.iter()
        if xml_local_name(node).lower() == "platformtoolset"}
    if all_toolsets != approved_toolsets:
        raise ContractError(
            "PlatformToolset is allowed only in the four configuration groups")
    all_configuration_types = {
        node for node in root.iter()
        if xml_local_name(node).lower() == "configurationtype"}
    if all_configuration_types != approved_configuration_types:
        raise ContractError(
            "ConfigurationType is allowed only in the four configuration "
            "groups")


def validate_per_file_isa(tree):
    compile_nodes = tree.findall(".//msb:ClCompile[@Include]", NS)
    options = {}
    for node in compile_nodes:
        path = project_path(PROJECT, node.attrib["Include"])
        if path in options:
            raise ContractError("duplicate per-file compile item: " + path)
        metadata = [xml_local_name(child) for child in list(node)]
        if path in AVX2_SOURCE_FILES:
            if metadata != ["AdditionalOptions"]:
                raise ContractError(
                    "AVX2 source metadata is not the exact contract")
        elif metadata:
            raise ContractError(
                "non-AVX2 source has per-file metadata: " + path)
        options[path] = " ".join(
            child.text or "" for child in node.findall(
                "msb:AdditionalOptions", NS))
        definitions = node.findall("msb:PreprocessorDefinitions", NS)
        if len(definitions) > 1:
            raise ContractError(
                "duplicate per-file PreprocessorDefinitions: " + path)
        if definitions:
            raw_definitions = definitions[0].text or ""
            reject_msbuild_expansion(
                raw_definitions, path + " PreprocessorDefinitions",
                ("%(PreprocessorDefinitions)",))
            macros = msbuild_values(raw_definitions)
            if "%(PreprocessorDefinitions)" not in macros:
                raise ContractError(
                    "per-file definitions do not inherit isolation macros: " +
                    path)
            if any(macro.startswith("LEO2_BACKEND_FORCE_")
                   for macro in macros):
                raise ContractError("per-file backend force: " + path)
        undefined = ";".join(
            child.text or "" for child in node.findall(
                "msb:UndefinePreprocessorDefinitions", NS))
        reject_msbuild_expansion(
            undefined, path + " UndefinePreprocessorDefinitions",
            ("%(UndefinePreprocessorDefinitions)",))
        undefined_macros = {
            macro.split("=")[0] for macro in msbuild_values(undefined)
        }
        if undefined_macros & PROTECTED_MACRO_NAMES:
            raise ContractError(
                "per-file metadata undefines isolation macros: " + path)
        reject_msbuild_expansion(
            options[path], path + " AdditionalOptions",
            ("%(AdditionalOptions)",))
        reject_backend_option_override(options[path], "per-file " + path)
        reject_wpo_options(options[path], "per-file " + path)
        enhanced = node.findall("msb:EnableEnhancedInstructionSet", NS)
        enhanced_text = " ".join(child.text or "" for child in enhanced)
        reject_msbuild_expansion(
            enhanced_text, path + " EnableEnhancedInstructionSet")
        if len(enhanced) > 1 or any(
                (child.text or "").strip() not in ("", "NotSet")
                for child in enhanced):
            raise ContractError(
                "per-file enhanced ISA metadata is unsupported: " + path)
    for path in AVX2_SOURCE_FILES:
        avx2 = options.get(path, "")
        if avx2.strip() != "/arch:AVX2 %(AdditionalOptions)":
            raise ContractError(
                "AVX2 backend options are not the exact contract: " + path)
    for path, flags in options.items():
        if path not in AVX2_SOURCE_FILES and re.search(
                r"/arch\s*:\s*avx", flags, re.IGNORECASE):
            raise ContractError("non-AVX2 source raises ISA: " + path)


def validate_visual_studio_project(tree):
    validate_msbuild_namespace(tree)
    validate_msbuild_element_casing(tree)
    validate_msbuild_imports_and_toolchain(tree)
    validate_source_item_structure(tree)
    validate_msbuild_configurations(tree)
    validate_msbuild_property_topology(tree)
    validate_per_file_isa(tree)
    validate_no_wpo_overrides(tree)
    validate_msbuild_root_order(tree)


def validate_legacy_solution_manifest(solution):
    actual = [line.strip() for line in solution.lstrip("\ufeff").splitlines()
              if line.strip()]
    project_type = "{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}"
    projects = (
        ("Leopard", "Leopard.vcxproj",
         "{32176592-2F30-4BD5-B645-EB11C8D3453E}"),
        ("LeopardBenchmark", "..\\tests\\proj\\Benchmark.vcxproj",
         "{97FCA15F-EAF3-4F1A-AFF8-83E693DA9D45}"),
        ("LeopardExperiments", "..\\tests\\proj\\Experiments.vcxproj",
         "{97FCA15F-EAF3-4F1A-AFF8-83E693DA9D65}"),
    )
    configurations = (
        "Debug|Win32", "Debug|x64", "Release|Win32", "Release|x64")
    expected = [
        "Microsoft Visual Studio Solution File, Format Version 12.00",
        "# Visual Studio 14",
        "VisualStudioVersion = 14.0.25420.1",
        "MinimumVisualStudioVersion = 10.0.40219.1",
    ]
    for name, path, guid in projects:
        expected.extend((
            'Project("' + project_type + '") = "' + name + '", "' +
            path + '", "' + guid + '"',
            "EndProject",
        ))
    expected.extend((
        "Global",
        "GlobalSection(SolutionConfigurationPlatforms) = preSolution",
    ))
    expected.extend(configuration + " = " + configuration
                    for configuration in configurations)
    expected.extend((
        "EndGlobalSection",
        "GlobalSection(ProjectConfigurationPlatforms) = postSolution",
    ))
    for unused_name, unused_path, guid in projects:
        for configuration in configurations:
            expected.extend((
                guid + "." + configuration + ".ActiveCfg = " +
                configuration,
                guid + "." + configuration + ".Build.0 = " +
                configuration,
            ))
    expected.extend((
        "EndGlobalSection",
        "GlobalSection(SolutionProperties) = preSolution",
        "HideSolutionNode = FALSE",
        "EndGlobalSection",
        "EndGlobal",
    ))
    if actual != expected:
        raise ContractError(
            "legacy Visual Studio solution graph/configuration manifest "
            "differs")


def validate_legacy_visual_studio_metadata():
    solution = SOLUTION.read_text(encoding="utf-8-sig")
    validate_legacy_solution_manifest(solution)
    if ("# Visual Studio 14" not in solution or
            "VisualStudioVersion = 14.0.25420.1" not in solution):
        raise ContractError(
            "legacy solution does not declare the repository's VS2015 version")
    for project in LEGACY_PROJECTS:
        tree = ET.parse(str(project))
        validate_msbuild_namespace(tree)
        validate_msbuild_element_casing(tree)
        validate_msbuild_imports_and_toolchain(tree)
        if tree.getroot().attrib.get("ToolsVersion") != EXPECTED_TOOLS_VERSION:
            raise ContractError(
                project.name + " does not use MSBuild ToolsVersion 14.0")
        toolsets = {
            (node.text or "").strip()
            for node in tree.findall(
                ".//msb:PropertyGroup[@Label='Configuration']/"
                "msb:PlatformToolset", NS)
        }
        if toolsets != {"v140"}:
            raise ContractError(
                project.name + " does not use the Visual Studio 2015 toolset")
        references = tree.findall(".//msb:ProjectReference", NS)
        if project == ROOT / "tests" / "proj" / "Benchmark.vcxproj":
            if (len(references) != 1 or references[0].attrib != {
                    "Include": "..\\..\\proj\\Leopard.vcxproj"} or
                    [xml_local_name(child)
                     for child in list(references[0])] != ["Project"] or
                    (references[0][0].text or "").strip().lower() !=
                    "{32176592-2f30-4bd5-b645-eb11c8d3453e}"):
                raise ContractError(
                    "Benchmark ProjectReference is not the exact contract")
        elif references:
            raise ContractError(
                project.name + " has an unapproved ProjectReference")


class LeopardVisualStudioProjectTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = ET.parse(str(PROJECT))
        cls.filters = ET.parse(str(FILTERS))
        (cls.expected_sources, cls.expected_headers,
         cls.object_targets, cls.object_sources,
         cls.cmake) = production_graph()
        cls.sources = item_paths(cls.project, "ClCompile", PROJECT)
        cls.headers = item_paths(cls.project, "ClInclude", PROJECT)

    def assert_unique_equal(self, actual, expected, label):
        duplicates = sorted(
            path for path, count in Counter(actual).items() if count != 1)
        self.assertEqual([], duplicates, label + " duplicated")
        self.assertEqual(sorted(expected), sorted(actual), label + " drifted")

    def test_compiled_sources_match_attached_cmake_graph(self):
        self.assert_unique_equal(
            self.sources, self.expected_sources, "compiled production sources")
        self.assertTrue(self.object_targets)
        self.assertTrue({
            "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
            "Leopard2BackendSSSE3.cpp"
        }.issubset(set(self.object_sources)))

    def test_reachable_production_headers_are_visible(self):
        self.assert_unique_equal(
            self.headers, self.expected_headers, "production headers")

    def test_filters_match_project_exactly(self):
        filter_sources = item_paths(self.filters, "ClCompile", FILTERS)
        filter_headers = item_paths(self.filters, "ClInclude", FILTERS)
        self.assert_unique_equal(filter_sources, self.sources, "source filters")
        self.assert_unique_equal(filter_headers, self.headers, "header filters")
        for node in self.filters.findall(".//msb:ClCompile", NS):
            self.assertEqual(
                "Source Files", node.findtext("msb:Filter", namespaces=NS))
        for node in self.filters.findall(".//msb:ClInclude", NS):
            self.assertEqual(
                "Header Files", node.findtext("msb:Filter", namespaces=NS))

    def test_exact_msbuild_configuration_contract(self):
        validate_visual_studio_project(self.project)

    def test_visual_studio_2015_metadata_is_internally_consistent(self):
        validate_legacy_visual_studio_metadata()

    def test_legacy_solution_graph_and_mappings_are_exact(self):
        solution = SOLUTION.read_text(encoding="utf-8-sig")
        leopard_guid = "{32176592-2F30-4BD5-B645-EB11C8D3453E}"
        mutations = (
            solution.replace(
                "\t" + leopard_guid +
                ".Debug|Win32.Build.0 = Debug|Win32\n", "", 1),
            solution.replace(
                "Global\n",
                'Project("{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}") = '
                '"Injected", "Injected.vcxproj", '
                '"{00000000-0000-0000-0000-000000000000}"\n'
                "EndProject\nGlobal\n", 1),
            solution.replace(
                leopard_guid +
                ".Release|x64.ActiveCfg = Release|x64",
                leopard_guid +
                ".Release|x64.ActiveCfg = Release|Win32", 1),
        )
        for mutated in mutations:
            with self.subTest(size=len(mutated)):
                self.assertNotEqual(solution, mutated)
                with self.assertRaisesRegex(
                        ContractError, "solution graph/configuration"):
                    validate_legacy_solution_manifest(mutated)

    def test_cuda_remains_absent_and_opt_in(self):
        project_text = PROJECT.read_text(encoding="utf-8-sig").lower()
        filters_text = FILTERS.read_text(encoding="utf-8-sig").lower()
        self.assertNotIn("cuda", project_text)
        self.assertNotIn("cuda", filters_text)
        self.assertEqual([], self.project.findall(".//msb:CudaCompile", NS))

        option = re.search(
            r"option\s*\(\s*LEO2_ENABLE_CUDA\s+\"[^\"]*\"\s+(ON|OFF)\s*\)",
            self.cmake, re.IGNORECASE)
        self.assertIsNotNone(option)
        self.assertEqual("OFF", option.group(1).upper())
        cuda_enable_count = 0
        condition_stack = []
        for name, body in cmake_commands(self.cmake):
            tokens = cmake_tokens(body)
            if name == "if":
                condition_stack.append(tokens)
            elif name == "endif":
                self.assertTrue(condition_stack)
                condition_stack.pop()
            elif name == "enable_language" and tokens == ["CUDA"]:
                cuda_enable_count += 1
                self.assertIn(["LEO2_ENABLE_CUDA"], condition_stack)
        self.assertEqual(1, cuda_enable_count)
        self.assertEqual([], condition_stack)

    def test_cuda_include_scanner_covers_angle_and_cu_forms(self):
        source = """
#include <Injected.cuh>
#include "Generated.cu"
# include <cuda_runtime.h>
#include \\
  "Spliced.cuh"
#define CUDA_HEADER "MacroInjected.cuh"
#pragma comment(lib, "cudart.lib")
__pragma(comment(lib, "nvcuda.lib"))
#include <cuda/std/atomic>
#include <thrust/device_vector.h>
#include <cub/cub.cuh>
"""
        self.assertEqual(
            ["Injected.cuh", "Generated.cu", "cuda_runtime.h",
             "Spliced.cuh", "cuda/std/atomic", "thrust/device_vector.h",
             "cub/cub.cuh"],
            production_local_includes(source))
        self.assertEqual(10, len(production_cuda_preprocessor_lines(source)))
        self.assertEqual(
            ["Hidden.inc", "Hidden.ipp", "Hidden"],
            production_local_includes(
                '#include "Hidden.inc"\n#include <Hidden.ipp>\n'
                '#include "Hidden"\n'))
        with self.assertRaisesRegex(ContractError, "nonliteral production include"):
            production_local_includes(
                '#define H "Hidden.inc"\n#include H\n')
        runtime_source = """
HMODULE module = LoadLibraryA("nvcuda.dll");
void *entry = GetProcAddress(module, "cuInit");
"""
        self.assertEqual(
            ['HMODULE module = LoadLibraryA("nvcuda.dll");'],
            production_cuda_source_lines(runtime_source))

        portability_header = (
            ROOT / "sse2neon" / "sse2neon.h").read_text(encoding="utf-8")
        author_comment = "//   Brandon Rowlett <browlett@nvidia.com>"
        self.assertIn(author_comment, portability_header)
        self.assertEqual(
            [], production_cuda_source_lines(portability_header))
        self.assertEqual(
            [], production_cuda_preprocessor_lines(portability_header))

        commented_examples = """
// #include <cuda_runtime.h>
// #define NVIDIA_BACKEND 1
/*
#include "kernel.cuh"
__global__ void documentation_only_kernel() {}
documentation_only_kernel<<<1, 1>>>();
const char *tool = "nvcc";
*/
#include "StillCpuOnly.h"
"""
        self.assertEqual(
            ["StillCpuOnly.h"],
            production_local_includes(commented_examples))
        self.assertEqual(
            [], production_cuda_source_lines(commented_examples))
        self.assertEqual(
            [], production_cuda_preprocessor_lines(commented_examples))

        cuda_mutations = (
            '#include <cuda_runtime.h>',
            '#include "kernel.cuh"',
            '#include <thrust/device_vector.h>',
            '#define NVIDIA_BACKEND 1',
            '#define __CUDACC__ 1',
            'const char *tool = "nvcc";',
            '#pragma comment(lib, "nvcuda.lib")',
            '__global__ void kernel() {}',
            'kernel<<<1, 1>>>();',
        )
        for mutation in cuda_mutations:
            with self.subTest(cuda_mutation=mutation):
                mutated = portability_header + "\n" + mutation + "\n"
                self.assertTrue(production_cuda_source_lines(mutated))

        self.assertTrue(production_cuda_source_lines(
            'const char *literal = R"tag(// nvidia runtime)tag";'))
        self.assertEqual(
            ROOT / "LeopardCommon.h",
            windows_repository_include(ROOT, "leopardcommon.h"))

    def test_c_family_comment_stripping_preserves_compiled_literals(self):
        source = (
            'const char *ordinary = "// nvidia /* cuda */"; '
            '// author@nvidia.com\n'
            'const char *raw = R"tag(/* cuda */ // nvcc)tag";\n'
            '/* #include <cuda_runtime.h>\n'
            'kernel<<<1, 1>>>(); */\n')
        stripped = strip_c_family_comments(source)
        self.assertIn('"// nvidia /* cuda */"', stripped)
        self.assertIn('R"tag(/* cuda */ // nvcc)tag"', stripped)
        self.assertNotIn("author@nvidia.com", stripped)
        self.assertNotIn("cuda_runtime.h", stripped)
        self.assertNotIn("kernel<<<1, 1>>>", stripped)
        self.assertEqual(source.count("\n"), stripped.count("\n"))

        with self.assertRaisesRegex(ContractError, "unterminated.*comment"):
            strip_c_family_comments("/* documentation")
        with self.assertRaisesRegex(ContractError, "unterminated.*raw string"):
            strip_c_family_comments('R"tag(documentation)wrong"')


class CMakeLexerTest(unittest.TestCase):
    def test_empty_quoted_arguments_and_list_elements_are_not_discarded(self):
        self.assertEqual(
            ["VALUE", "", "tail"], cmake_tokens('VALUE "" tail'))
        self.assertEqual(
            ["VALUE", "", "tail"], cmake_tokens("VALUE;;tail"))

    def test_bracket_comments_are_skipped_at_every_delimiter_depth(self):
        for equals in ("", "=", "==="):
            opening = "#[" + equals + "["
            closing = "]" + equals + "]"
            text = (
                opening + "\nmessage( hidden ( ) # \\\" )\n" + closing +
                "\nset(SAFE_VALUE retained)\n")
            with self.subTest(equals=equals):
                commands = cmake_commands(text)
                self.assertEqual(["set"], [name for name, unused in commands])
                self.assertEqual(
                    ["SAFE_VALUE", "retained"], cmake_tokens(commands[0][1]))

    def test_bracket_arguments_are_balanced_atomically_and_fail_closed(self):
        for equals in ("", "=", "==="):
            opening = "[" + equals + "["
            closing = "]" + equals + "]"
            text = (
                "message(" + opening + "( # \\\"\n" + closing + ")\n"
                "set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)\n"
                "message(" + opening + ") # \\\"\n" + closing + ")\n")
            with self.subTest(equals=equals):
                commands = cmake_commands(text)
                self.assertEqual(
                    ["message", "set", "message"],
                    [name for name, unused in commands])
                with self.assertRaisesRegex(
                        ContractError, "bracket arguments are unsupported"):
                    cmake_tokens(commands[0][1])

    def test_unterminated_bracket_constructs_are_rejected(self):
        for text in ("#[=[ unterminated", "message([==[ unterminated)"):
            with self.subTest(text=text):
                with self.assertRaisesRegex(
                        ContractError, "unterminated CMake bracket"):
                    cmake_commands(text)


class CMakeGraphMutationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cmake = CMAKE.read_text(encoding="utf-8")

    def resolve(self, mutation):
        return production_graph(
            self.cmake + "\n" + mutation + "\n", require_files=False,
            require_mutation_contract=True)

    def resolve_text(self, text, require_mutation_contract=False):
        return production_graph(
            text, require_files=False,
            require_mutation_contract=require_mutation_contract)

    def test_portable_release_audit_is_opt_in_and_platform_scoped(self):
        option = (
            'option(LEO2_PORTABLE_ISA_RELEASE_AUDIT\n'
            '    "Require the strict x86-64 portable-ISA release audit" OFF)')
        self.assertEqual(1, self.cmake.count(option))
        required_scope = (
            'if(LEO2_PORTABLE_ISA_RELEASE_AUDIT)\n'
            '    if(MSVC OR WIN32)')
        self.assertEqual(1, self.cmake.count(required_scope))
        self.assertEqual(
            1, self.cmake.count(
                'LEO2_PORTABLE_ISA_RELEASE_AUDIT currently supports x86-64 only'))
        self.assertEqual(
            1, self.cmake.count(
                'LEO2_PORTABLE_ISA_RELEASE_AUDIT requires a single-configuration '))
        self.assertEqual(
            1, self.cmake.count('if(CMAKE_GENERATOR MATCHES "Multi-Config")'))
        self.assertEqual(
            1, self.cmake.count('"CMAKE_EXPORT_COMPILE_COMMANDS=ON")'))

        default_on = self.cmake.replace(option, option[:-4] + 'ON)', 1)
        self.assertNotEqual(default_on, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "unapproved (?:CMake cache option|production "
                "compiler-control variable mutation)|missing or duplicate "
                "trusted CMake command"):
            self.resolve_text(default_on, require_mutation_contract=True)

    def test_legacy_target_references_are_explicitly_classified(self):
        # These four experiment families authenticate the historical CMake
        # target and archive names only in explicit versioned replay maps (and
        # their mutation fixtures). New evidence uses the canonical identity.
        # No production tool or new experiment family may silently add another
        # functional reference to the historical names.
        compatibility_contracts = {
            "CMakeLists.txt",
            "cmake/leopardConfig.cmake.in",
            "tests/cmake/test_cuda_optional.cmake",
            "tests/proj/test_leopard_vcxproj.py",
        }
        authenticated_replay_contracts = {
            "experiments/leopard2/backend_butterfly/run_abba.py",
            "experiments/leopard2/low_encode_copy/run_abba.py",
            "experiments/leopard2/main_compare/run_abba.py",
            "experiments/leopard2/main_compare/test_run_abba.py",
            "experiments/leopard2/non_power_of_two/c7/run_authoritative.py",
            "experiments/leopard2/non_power_of_two/c7/test_checkpoint.py",
            "experiments/leopard2/non_power_of_two/c7/validate_evidence.py",
        }
        legacy_reference = re.compile(
            r"liblibleopard|CMakeFiles/libleopard[.]dir|"
            r"leopard::libleopard|TARGET_FILE:libleopard|"
            r"(?:add_library|target_[A-Za-z_]+)[(][ \t]*libleopard\b|"
            r"install[(][ \t]*TARGETS[ \t]+libleopard\b|"
            r"(?:\"libleopard\"|'libleopard')|"
            r"--target[^\n]*\blibleopard\b")
        candidates = {CMAKE}
        for relative_root in ("cmake", "tests", "tools", "experiments/leopard2"):
            root = ROOT / relative_root
            for path in root.rglob("*"):
                if (not path.is_file() or "results" in path.parts or
                        path.suffix.lower() not in {
                            ".cmake", ".in", ".py", ".sh"}):
                    continue
                candidates.add(path)

        actual = set()
        for path in candidates:
            text = path.read_text(encoding="utf-8")
            if legacy_reference.search(text):
                actual.add(path.relative_to(ROOT).as_posix())

        self.assertEqual(
            compatibility_contracts | authenticated_replay_contracts,
            actual,
            "legacy CMake target/archive reference allowlist drifted")

        self.assertEqual(
            1, self.cmake.count("add_library(libleopard ALIAS leopard)"))
        self.assertNotIn("liblibleopard", self.cmake)

    def test_direct_target_sources_is_retained(self):
        (sources, unused_headers, unused_objects,
         unused_object_sources, unused_cmake) = self.resolve(
            "target_sources(leopard PRIVATE New.cpp)")
        del unused_headers, unused_objects, unused_object_sources, unused_cmake
        self.assertIn("New.cpp", sources)

    def test_future_attached_object_target_is_traversed(self):
        mutation = """
add_library(future_backend OBJECT FutureBackend.cpp)
target_sources(leopard PRIVATE $<TARGET_OBJECTS:future_backend>)
"""
        (sources, unused_headers, objects,
         object_sources, unused_cmake) = self.resolve(mutation)
        del unused_headers, unused_cmake
        self.assertIn("FutureBackend.cpp", sources)
        self.assertIn("future_backend", objects)
        self.assertIn("FutureBackend.cpp", object_sources)

    def test_resolved_target_objects_variable_is_traversed(self):
        mutation = """
add_library(variable_backend OBJECT VariableBackend.cpp)
set(OBJECT_EXPRESSION $<TARGET_OBJECTS:variable_backend>)
target_sources(leopard PRIVATE ${OBJECT_EXPRESSION})
"""
        (sources, unused_headers, objects,
         object_sources, unused_cmake) = self.resolve(mutation)
        del unused_headers, unused_cmake
        self.assertIn("VariableBackend.cpp", sources)
        self.assertIn("variable_backend", objects)
        self.assertIn("VariableBackend.cpp", object_sources)

    def test_unresolved_variable_is_rejected(self):
        with self.assertRaisesRegex(ContractError, "unresolved.*EXTRA_SOURCES"):
            self.resolve(
                "target_sources(leopard PRIVATE ${EXTRA_SOURCES})")

    def test_resolved_variable_is_retained(self):
        mutation = """
set(EXTRA_SOURCES ExtraOne.cpp ExtraTwo.cpp)
target_sources(leopard PRIVATE ${EXTRA_SOURCES})
"""
        (sources, unused_headers, unused_objects,
         unused_object_sources, unused) = self.resolve(mutation)
        del unused_headers, unused_objects, unused_object_sources, unused
        self.assertIn("ExtraOne.cpp", sources)
        self.assertIn("ExtraTwo.cpp", sources)

    def test_conditional_source_variable_is_rejected(self):
        mutation = """
if(WIN32)
    set(BRANCH_SOURCES BranchWin.cpp)
else()
    set(BRANCH_SOURCES BranchOther.cpp)
endif()
target_sources(leopard PRIVATE ${BRANCH_SOURCES})
"""
        with self.assertRaisesRegex(
                ContractError, "conditional.*BRANCH_SOURCES"):
            self.resolve(mutation)

    def test_unsupported_list_operation_on_source_variable_is_rejected(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        mutation = "list(REMOVE_ITEM LIB_SOURCE_FILES leopard2.cpp)"
        text = self.cmake.replace(marker, mutation + "\n" + marker, 1)
        with self.assertRaisesRegex(
                ContractError,
                "unsupported (?:CMake )?list operation.*LIB_SOURCE_FILES"):
            self.resolve_text(text)

    def test_source_variable_is_snapshotted_at_add_library_time(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        replacement = """set(SAVED_LIB_SOURCE_FILES ${LIB_SOURCE_FILES})
set(LIB_SOURCE_FILES Injected.cpp)
add_library(leopard STATIC ${LIB_SOURCE_FILES})
set(LIB_SOURCE_FILES ${SAVED_LIB_SOURCE_FILES})"""
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        sources = self.resolve_text(text)[0]
        self.assertIn("Injected.cpp", sources)
        self.assertNotIn("leopard2.cpp", sources)

    def test_indirect_mutation_destinations_are_resolved_at_command_time(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        mutations = (
            "set(${SOURCE_DEST} Injected.cpp)",
            "list(APPEND ${SOURCE_DEST} Injected.cpp)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                replacement = (
                    "set(SOURCE_DEST LIB_SOURCE_FILES)\n" + mutation +
                    "\n" + marker)
                text = self.cmake.replace(marker, replacement, 1)
                self.assertIn("Injected.cpp", self.resolve_text(text)[0])

    def test_unmodeled_source_variable_writers_are_rejected(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        mutations = (
            'file(GLOB LIB_SOURCE_FILES "*.cpp")',
            ('set(GLOB_OUTPUT LIB_SOURCE_FILES)\n'
             'file(GLOB ${GLOB_OUTPUT} "*.cpp")'),
            ('set(SOURCE_SUFFIX SOURCE_FILES)\n'
             'file(GLOB LIB_${SOURCE_SUFFIX} "*.cpp")'),
            'string(APPEND LIB_SOURCE_FILES ";Injected.cpp")',
            'custom_source_mutator(LIB_SOURCE_FILES)',
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                text = self.cmake.replace(
                    marker, mutation + "\n" + marker, 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "may mutate source variable|generated/custom build|"
                        "unmodeled root CMake command"):
                    self.resolve_text(text)

    def test_unproved_graph_imports_are_rejected(self):
        for mutation in (
                "add_subdirectory(cmake/injected_sources)",
                "subdirs(cmake/injected_sources)",
                "include(cmake/injected_sources.cmake)",
                "cmake_language(CALL target_sources leopard PRIVATE "
                "Injected.cpp)",
                "load_command(injected cmake/injected-command)",
                "find_package(Injected CONFIG REQUIRED PATHS "
                "cmake/injected NO_DEFAULT_PATH)"):
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "graph proof|unapproved CMake (?:graph|package)|"
                        "dynamic command"):
                    self.resolve(mutation)

    def test_approved_package_targets_cannot_be_locally_shadowed(self):
        mutations = (
            "add_library(Threads::Threads INTERFACE IMPORTED)",
            "add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)",
            "add_library(injected_runtime INTERFACE)\n"
            "add_library(Threads::Threads ALIAS injected_runtime)",
            "add_library(injected_openmp INTERFACE)\n"
            "add_library(OpenMP::OpenMP_CXX ALIAS injected_openmp)",
        )
        for mutation in mutations:
            with self.subTest(target=mutation.splitlines()[-1]):
                with self.assertRaisesRegex(
                        ContractError,
                        "approved package target cannot be locally declared"):
                    self.resolve(mutation)

    def test_find_package_in_unsupported_block_is_rejected(self):
        for block, closing in (
                ("function(inject_package)", "endfunction()"),
                ("macro(inject_package)", "endmacro()")):
            for package in (
                    "find_package(Injected CONFIG REQUIRED)",
                    "find_package(Threads REQUIRED)"):
                with self.subTest(block=block, package=package):
                    mutation = (block + "\n    " + package + "\n" + closing +
                                "\ninject_package()")
                    with self.assertRaisesRegex(
                            ContractError,
                            "unsupported CMake block"):
                        self.resolve(mutation)

    def test_object_library_link_attachment_is_rejected(self):
        mutations = (
            """
add_library(linked_backend OBJECT LinkedBackend.cpp)
target_link_libraries(leopard PRIVATE linked_backend)
""",
            """
add_library(linked_backend OBJECT LinkedBackend.cpp)
set(LINK_DESTINATION leopard)
target_link_libraries(${LINK_DESTINATION} PRIVATE linked_backend)
""",
            """
add_library(linked_backend OBJECT LinkedBackend.cpp)
function(inject_link)
    target_link_libraries(leopard PRIVATE linked_backend)
endfunction()
inject_link()
""",
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation.strip().splitlines()[-1]):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation|"
                        "compile/link mutation in unsupported CMake block|"
                        "unsupported CMake block"):
                    self.resolve(mutation)

    def test_equivalent_object_library_link_mutations_are_rejected(self):
        object_definition = (
            "add_library(linked_backend OBJECT LinkedBackend.cpp)")
        mutations = (
            "set_property(TARGET leopard APPEND PROPERTY "
            "LINK_LIBRARIES linked_backend)",
            "set_target_properties(leopard PROPERTIES "
            "LINK_LIBRARIES linked_backend)",
            "set_property(TARGET leopard APPEND PROPERTY "
            "INTERFACE_LINK_LIBRARIES linked_backend)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "link property bypasses graph"):
                    self.resolve(object_definition + "\n" + mutation)

        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        directory_mutations = (
            "link_libraries(linked_backend)",
            "set_property(DIRECTORY APPEND PROPERTY "
            "LINK_LIBRARIES linked_backend)",
            "set_directory_properties(PROPERTIES "
            "LINK_LIBRARIES linked_backend)",
        )
        for mutation in directory_mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                replacement = (
                    object_definition + "\n" + mutation + "\n" + marker)
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "directory compile/link graph|"
                        "compile/link property bypasses graph"):
                    self.resolve_text(text)

    def test_attached_object_compile_and_link_mutations_are_rejected(self):
        prelude = """
add_library(injected_backend OBJECT InjectedBackend.cpp)
target_sources(leopard PRIVATE $<TARGET_OBJECTS:injected_backend>)
"""
        mutations = (
            "target_link_libraries(injected_backend PRIVATE Injected::Injected)",
            "target_compile_options(injected_backend PRIVATE /arch:AVX2)",
            "target_compile_definitions(injected_backend PRIVATE "
            "LEO2_DISABLE_SSSE3_CODEGEN=0)",
            "target_compile_features(injected_backend PRIVATE cxx_std_20)",
            "target_include_directories(injected_backend PRIVATE injected)",
            "target_link_options(injected_backend PRIVATE /LTCG)",
            "target_link_directories(injected_backend PRIVATE injected)",
            "target_link_interface_libraries(injected_backend Injected::Injected)",
            "target_precompile_headers(injected_backend PRIVATE injected.h)",
            "set(MUTATION_TARGET injected_backend)\n"
            "target_compile_options(${MUTATION_TARGET} PRIVATE /arch:AVX2)",
            "TaRgEt_CoMpIlE_OpTiOnS(injected_backend PRIVATE /arch:AVX2)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.splitlines()[-1].split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation"):
                    self.resolve(prelude + mutation)

    def test_leopard_compile_mutations_are_rejected(self):
        mutations = (
            "target_compile_options(leopard PRIVATE /arch:AVX2)",
            "target_compile_definitions(leopard PRIVATE INJECTED=1)",
            "target_compile_definitions(leopard PRIVATE "
            "LEO2_BACKEND_FORCE_AVX2=1)",
            "target_compile_features(leopard PRIVATE cxx_std_20)",
            "target_include_directories(leopard PRIVATE injected)",
            "target_link_options(leopard PRIVATE /LTCG)",
            "target_link_directories(leopard PRIVATE injected)",
            "target_link_interface_libraries(leopard Injected::Injected)",
            "target_precompile_headers(leopard PRIVATE injected.h)",
            "set(MUTATION_TARGET leopard)\n"
            "target_compile_options(${MUTATION_TARGET} PRIVATE /arch:AVX2)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.splitlines()[-1].split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation|"
                        "unsupported condition guards production target mutation"):
                    self.resolve(mutation)

    def test_test_hooks_cannot_reenter_production_targets(self):
        mutations = (
            "target_compile_definitions(leopard PRIVATE "
            "LEO2_ENABLE_TEST_HOOKS=1)",
            "target_compile_definitions(leopard2_backend_ssse3 PRIVATE "
            "LEO2_ENABLE_TEST_HOOKS=1)",
            "target_compile_definitions(leopard2_backend_avx2 PRIVATE "
            "LEO2_ENABLE_TEST_HOOKS=1)",
        )
        for mutation in mutations:
            with self.subTest(target=mutation.split("(", 1)[1].split()[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation"):
                    self.resolve(mutation)

    def test_target_alias_cannot_hide_production_mutation(self):
        mutations = (
            """
add_library(leopard_alias ALIAS leopard)
target_compile_options(leopard_alias PRIVATE /arch:AVX2)
""",
            """
add_library(alias_backend OBJECT AliasBackend.cpp)
add_library(backend_alias ALIAS alias_backend)
target_sources(leopard PRIVATE $<TARGET_OBJECTS:alias_backend>)
target_compile_definitions(backend_alias PRIVATE INJECTED=1)
""",
            """
add_library(alias_backend OBJECT AliasBackend.cpp)
add_library(backend_alias ALIAS alias_backend)
target_sources(leopard PRIVATE $<TARGET_OBJECTS:backend_alias>)
target_compile_options(alias_backend PRIVATE /arch:AVX2)
""",
        )
        for mutation in mutations:
            with self.subTest(last=mutation.strip().splitlines()[-1]):
                with self.assertRaisesRegex(
                        ContractError,
                        "target ALIAS cannot be mutated|"
                        "unapproved production target compile/link mutation"):
                    self.resolve(mutation)

    def test_scoped_compile_and_link_mutations_are_rejected(self):
        commands = (
            "target_compile_options(leopard PRIVATE /arch:AVX2)",
            "target_compile_definitions(leopard PRIVATE INJECTED=1)",
            "target_include_directories(leopard PRIVATE injected)",
            "target_link_options(leopard PRIVATE /LTCG)",
        )
        for opening, closing in (
                ("function(inject)", "endfunction()"),
                ("macro(inject)", "endmacro()")):
            for command in commands:
                mutation = (opening + "\n" + command + "\n" + closing +
                            "\ninject()")
                with self.subTest(scope=opening.split("(", 1)[0],
                                  command=command.split("(", 1)[0]):
                    with self.assertRaisesRegex(
                            ContractError,
                            "compile/link mutation in unsupported CMake block|"
                            "unsupported CMake block"):
                        self.resolve(mutation)

    def test_directory_compile_and_link_mutations_are_rejected(self):
        mutations = (
            "add_compile_options(/arch:AVX2)",
            "add_compile_definitions(INJECTED=1)",
            "add_definitions(/DLEO2_DISABLE_SSSE3_CODEGEN=0)",
            "remove_definitions(/DLEO2_DISABLE_SSSE3_CODEGEN=1)",
            "include_directories(injected)",
            "add_link_options(/LTCG)",
            "link_directories(injected)",
            "link_libraries(injected)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "directory compile/link graph"):
                    self.resolve(mutation)

    def test_compile_link_properties_are_rejected(self):
        mutations = (
            "set_property(TARGET leopard APPEND PROPERTY "
            "COMPILE_OPTIONS /arch:AVX2)",
            "set_target_properties(leopard PROPERTIES "
            "COMPILE_DEFINITIONS INJECTED=1)",
            "set_property(TARGET leopard APPEND PROPERTY LINK_OPTIONS /LTCG)",
            "set_target_properties(leopard PROPERTIES "
            "INTERPROCEDURAL_OPTIMIZATION ON)",
            "set_target_properties(leopard PROPERTIES "
            "INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)",
            "set_target_properties(leopard PROPERTIES "
            "COMPILE_DEFINITIONS_RELEASE INJECTED=1)",
            "set_target_properties(leopard PROPERTIES "
            "MSVC_RUNTIME_LIBRARY MultiThreaded)",
            "set_property(TARGET Threads::Threads PROPERTY "
            "IMPORTED_LOCATION_RELEASE Injected.lib)",
            "set_property(TARGET OpenMP::OpenMP_CXX PROPERTY "
            "INTERFACE_COMPILE_OPTIONS /arch:AVX2)",
            "set_property(TARGET leopard PROPERTY "
            "INTERFACE_LINK_LIBRARIES_DIRECT Injected::Injected)",
            "set_property(TARGET leopard PROPERTY "
            "INTERFACE_SOURCES C:/Injected.cpp)",
            "set_property(TARGET leopard PROPERTY "
            "VS_USER_PROPS C:/Injected.props)",
            "set_property(TARGET leopard PROPERTY "
            "VS_PLATFORM_TOOLSET injected)",
            "set_property(TARGET leopard PROPERTY "
            "POSITION_INDEPENDENT_CODE ON)",
            "set_property(TARGET leopard APPEND PROPERTY "
            "STATIC_LIBRARY_OPTIONS /LTCG)",
            "set_property(DIRECTORY APPEND PROPERTY COMPILE_OPTIONS /arch:AVX2)",
            "set_directory_properties(PROPERTIES LINK_OPTIONS /LTCG)",
            "set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE injected)",
            "set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK injected)",
        )
        for mutation in mutations:
            with self.subTest(property=mutation):
                with self.assertRaisesRegex(
                        ContractError, "compile/link property bypasses graph"):
                    self.resolve(mutation)

    def test_compiler_control_variables_are_rejected(self):
        mutations = (
            "set(CMAKE_CXX_FLAGS /arch:AVX2)",
            "set(CMAKE_CXX_FLAGS_RELEASE \"${CMAKE_CXX_FLAGS_RELEASE} /GL\")",
            "string(APPEND CMAKE_CXX_FLAGS \" /arch:AVX2\")",
            "list(APPEND CMAKE_CXX_FLAGS /arch:AVX2)",
            "unset(CMAKE_CXX_FLAGS)",
            "set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)",
            "set(CMAKE_CXX_COMPILER_LAUNCHER injected)",
            "set(CMAKE_CXX_ARCHIVE_CREATE injected)",
            "set(CMAKE_STATIC_LINKER_FLAGS /LTCG)",
            "set(CMAKE_CXX_STANDARD 20)",
            "set(CMAKE_CXX_EXTENSIONS ON)",
            "set(CMAKE_UNITY_BUILD ON)",
            "set(CMAKE_VS_GLOBALS CLToolExe=Injected.exe)",
            "set(CMAKE_CURRENT_SOURCE_DIR injected)",
            "set(CMAKE_INSTALL_INCLUDEDIR injected)",
            "set(CMAKE_CXX_COMPILER_ID GNU)",
            "set(CMAKE_SYSTEM_PROCESSOR ARM64)",
            "set(CMAKE_USER_MAKE_RULES_OVERRIDE injected.cmake)",
            "set(CMAKE_PROJECT_INCLUDE injected.cmake)",
            "set(CMAKE_TOOLCHAIN_FILE injected.cmake)",
            "set(CMAKE_AR injected)",
            "set_property(CACHE CMAKE_CXX_FLAGS PROPERTY VALUE /arch:AVX2)",
            "set_property(CACHE OpenMP_CXX_FLAGS PROPERTY VALUE /arch:AVX2)",
            "set(OpenMP_CXX_LIB_NAMES evil CACHE STRING \"\" FORCE)",
            "set(OpenMP_evil_LIBRARY C:/evil.lib CACHE FILEPATH \"\" FORCE)",
            "set(OPENMP_FOUND TRUE)",
            "set(Threads_FOUND TRUE)",
            "set(THREADS_PTHREAD_ARG injected)",
            "set(LEOPARD_INSTALL_CMAKEDIR C:/evil)",
            "set(LEO2_FLAG_ARCH_AVX2 TRUE CACHE BOOL \"\" FORCE)",
            "option(LEO2_FLAG_ARCH_AVX2 \"poison\" ON)",
            "set(CXX_FLAG_O2 TRUE CACHE BOOL \"\" FORCE)",
            "set(ENABLE_OPENMP OFF CACHE BOOL \"\" FORCE)",
            "option(ENABLE_OPENMP \"poison\" OFF)",
            "set(LEO2_ENABLE_CUDA ON CACHE BOOL \"\" FORCE)",
            "set(LEO2_PORTABLE_ISA_RELEASE_AUDIT ON CACHE BOOL \"\" FORCE)",
            "set(LEO2_BUILD_TESTS OFF CACHE BOOL \"\" FORCE)",
            "set(LEO2_BUILD_FUZZERS ON CACHE BOOL \"\" FORCE)",
            "set(LEO2_BACKEND_VARIANT avx2 CACHE STRING \"\" FORCE)",
            "set_property(CACHE OpenMP_CXX_LIB_NAMES PROPERTY VALUE evil)",
            "set_property(CACHE Threads_FOUND PROPERTY VALUE TRUE)",
            "set_property(CACHE LEO2_FLAG_ARCH_AVX2 PROPERTY VALUE TRUE)",
            "set_property(CACHE ENABLE_OPENMP PROPERTY VALUE OFF)",
            "set(FLAG_DEST CMAKE_CXX_FLAGS)\n"
            "set(${FLAG_DEST} /arch:AVX2)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.splitlines()[-1]):
                with self.assertRaisesRegex(
                        ContractError, "compiler-control variable"):
                    self.resolve(mutation)

    def test_unapproved_compiler_probes_are_rejected(self):
        mutations = (
            'check_cxx_compiler_flag("/arch:AVX" LEO2_FLAG_ARCH_AVX2)',
            'check_cxx_compiler_flag("/GL" INJECTED_FLAG_GL)',
            'check_c_compiler_flag("/O2" INJECTED_C_FLAG_O2)',
        )
        for mutation in mutations:
            with self.subTest(probe=mutation):
                with self.assertRaisesRegex(
                        ContractError, "unapproved production compiler probe"):
                    self.resolve(mutation)

    def test_cache_state_mutation_properties_are_rejected(self):
        mutations = (
            "set_property(CACHE LEO2_BACKEND_VARIANT PROPERTY VALUE avx2)",
            "set_property(CACHE LEO2_ENABLE_CUDA PROPERTY VALUE ON)",
            "set_property(CACHE LEO2_PORTABLE_ISA_RELEASE_AUDIT PROPERTY VALUE ON)",
            "set_property(CACHE LEO2_BUILD_TESTS PROPERTY VALUE OFF)",
        )
        for mutation in mutations:
            with self.subTest(property=mutation):
                with self.assertRaisesRegex(
                        ContractError,
                        "non-target property mutation|compiler-control variable"):
                    self.resolve(mutation)

    def test_backend_variant_normalizer_is_exact(self):
        marker = ('string(TOLOWER "${LEO2_BACKEND_VARIANT}" '
                  'LEO2_BACKEND_VARIANT_NORMALIZED)')
        replacements = (
            "string(CONCAT LEO2_BACKEND_VARIANT_NORMALIZED avx2)",
            'string(TOLOWER "avx2" LEO2_BACKEND_VARIANT_NORMALIZED)',
            "string(REPLACE auto avx2 LEO2_BACKEND_VARIANT_NORMALIZED "
            "${LEO2_BACKEND_VARIANT})",
        )
        for replacement in replacements:
            with self.subTest(normalizer=replacement):
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "compiler-control variable|trusted CMake command"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

    def test_find_program_cannot_rewrite_production_control_state(self):
        cases = (
            ("add_library(leopard STATIC ${LIB_SOURCE_FILES})",
             "find_program(LEO2_X86_TARGET NAMES definitely_missing)"),
            ("if(LEO2_HAVE_SSSE3_BACKEND)",
             "find_program(LEO2_HAVE_SSSE3_BACKEND NAMES definitely_missing)"),
        )
        for marker, mutation in cases:
            with self.subTest(variable=mutation.split("(", 1)[1].split()[0]):
                text = self.cmake.replace(
                    marker, mutation + "\n" + marker, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaises(ContractError):
                    self.resolve_text(
                        text, require_mutation_contract=True)

    def test_locator_provenance_commands_are_exact_and_required(self):
        find_command = (
            "find_program(LEO2_LOCATOR_GIT_EXECUTABLE NAMES git)")
        revision_command = """execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_GIT_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_GIT_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)"""
        tree_command = """execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE} rev-parse HEAD^{tree}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_TREE_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_TREE_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)"""
        status_command = """execute_process(
            COMMAND ${LEO2_LOCATOR_GIT_EXECUTABLE}
                status --porcelain --untracked-files=normal
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE LEO2_LOCATOR_STATUS_RESULT
            OUTPUT_VARIABLE LEO2_LOCATOR_STATUS_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE)"""
        commands = (find_command, revision_command, tree_command, status_command)
        for command in commands:
            with self.subTest(required=command.split("(", 1)[0]):
                self.assertEqual(1, self.cmake.count(command))
                removed = self.cmake.replace(command, "", 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "locator provenance command (?:guard drift|"
                        "missing|duplicate)|missing or duplicate locator"):
                    self.resolve_text(
                        removed, require_mutation_contract=True)

            with self.subTest(duplicate=command.split("(", 1)[0]):
                duplicated = self.cmake.replace(
                    command, command + "\n" + command, 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate locator provenance command"):
                    self.resolve_text(
                        duplicated, require_mutation_contract=True)

        mutations = (
            ("NAMES git)", "NAMES git.exe)"),
            ("rev-parse HEAD", "rev-parse --verify HEAD"),
            ("RESULT_VARIABLE LEO2_LOCATOR_GIT_RESULT",
             "RESULT_VARIABLE LEO2_LOCATOR_GIT_RC"),
            ("OUTPUT_VARIABLE LEO2_LOCATOR_GIT_OUTPUT",
             "OUTPUT_VARIABLE LEO2_LOCATOR_GIT_SHA"),
            ("rev-parse HEAD^{tree}", "rev-parse HEAD^{commit}"),
            ("RESULT_VARIABLE LEO2_LOCATOR_TREE_RESULT",
             "RESULT_VARIABLE LEO2_LOCATOR_TREE_RC"),
            ("OUTPUT_VARIABLE LEO2_LOCATOR_TREE_OUTPUT",
             "OUTPUT_VARIABLE LEO2_LOCATOR_TREE_SHA"),
            ("status --porcelain --untracked-files=normal",
             "status --porcelain --untracked-files=no"),
            ("RESULT_VARIABLE LEO2_LOCATOR_STATUS_RESULT",
             "RESULT_VARIABLE LEO2_LOCATOR_STATUS_RC"),
            ("OUTPUT_VARIABLE LEO2_LOCATOR_STATUS_OUTPUT",
             "OUTPUT_VARIABLE LEO2_LOCATOR_STATUS_TEXT"),
            ("OUTPUT_VARIABLE LEO2_LOCATOR_GIT_OUTPUT\n"
             "            OUTPUT_STRIP_TRAILING_WHITESPACE)",
             "OUTPUT_VARIABLE LEO2_LOCATOR_GIT_OUTPUT)"),
            ("WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}\n"
             "            RESULT_VARIABLE LEO2_LOCATOR_GIT_RESULT",
             "WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}\n"
             "            RESULT_VARIABLE LEO2_LOCATOR_GIT_RESULT"),
        )
        for original, replacement in mutations:
            with self.subTest(mutation=replacement):
                mutated = self.cmake.replace(original, replacement, 1)
                self.assertNotEqual(mutated, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved tool discovery|generated/custom build "
                        "extension"):
                    self.resolve_text(
                        mutated, require_mutation_contract=True)

        benchmark_guard = (
            "if(LEO2_BUILD_BENCHMARKS)\n"
            "    add_executable(bench_leopard ${BENCH_SOURCE_FILES})")
        moved_find = self.cmake.replace(
            benchmark_guard,
            "if(NOT LEO2_BUILD_BENCHMARKS)\n"
            "    add_executable(bench_leopard ${BENCH_SOURCE_FILES})", 1)
        self.assertNotEqual(moved_find, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "locator provenance command guard drift"):
            self.resolve_text(moved_find, require_mutation_contract=True)

        inverted_probe = self.cmake.replace(
            "if(LEO2_LOCATOR_GIT_EXECUTABLE)\n        execute_process(",
            "if(NOT LEO2_LOCATOR_GIT_EXECUTABLE)\n        execute_process(",
            1)
        self.assertNotEqual(inverted_probe, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "locator provenance command guard drift"):
            self.resolve_text(inverted_probe, require_mutation_contract=True)

        sentinel = "__LEO2_LOCATOR_PROVENANCE_SENTINEL__"
        reordered = self.cmake.replace(revision_command, sentinel, 1)
        reordered = reordered.replace(tree_command, revision_command, 1)
        reordered = reordered.replace(sentinel, tree_command, 1)
        self.assertNotEqual(reordered, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "security-sensitive CMake command order drift"):
            self.resolve_text(reordered, require_mutation_contract=True)

    def test_sparse_evidence_sidecar_commands_are_exact_guarded_and_required(
            self):
        link_depends = """set_property(TARGET bench_leopard2_sparse_encode APPEND PROPERTY
            LINK_DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteSparseEncodeEvidenceSidecar.cmake)"""
        post_build = """add_custom_command(TARGET bench_leopard2_sparse_encode POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                -DOUTPUT=$<TARGET_FILE:bench_leopard2_sparse_encode>.leopard2-evidence
                -DEXECUTABLE=$<TARGET_FILE:bench_leopard2_sparse_encode>
                -DPRODUCTION_ARCHIVE=$<TARGET_FILE:leopard>
                -DBENCHMARK_OBJECT=$<TARGET_OBJECTS:leopard2_sparse_encode_benchmark_object>
                -DORACLE_OBJECT=$<TARGET_OBJECTS:leopard2_sparse_encode_oracle_object>
                -DBENCHMARK_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/bench_leopard2_sparse_encode.dir/link.txt
                -DPRODUCTION_LINK_RECIPE=${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/leopard.dir/link.txt
                -DBUILD_PROGRAM=${CMAKE_MAKE_PROGRAM}
                -DBUILD_ROOT=${CMAKE_CURRENT_BINARY_DIR}
                -DBUILD_GENERATOR=${CMAKE_GENERATOR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteSparseEncodeEvidenceSidecar.cmake
            VERBATIM)"""
        commands = (link_depends, post_build)
        for command in commands:
            with self.subTest(required=command.split("(", 1)[0]):
                self.assertEqual(1, self.cmake.count(command))
                removed = self.cmake.replace(command, "", 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate sparse evidence sidecar command"):
                    self.resolve_text(
                        removed, require_mutation_contract=True)

            with self.subTest(duplicate=command.split("(", 1)[0]):
                duplicated = self.cmake.replace(
                    command, command + "\n        " + command, 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate sparse evidence sidecar command"):
                    self.resolve_text(
                        duplicated, require_mutation_contract=True)

        adversarial_mutations = (
            (link_depends, link_depends.replace(
                "TARGET bench_leopard2_sparse_encode",
                "TARGET leopard", 1)),
            (link_depends, link_depends.replace(
                "APPEND PROPERTY", "APPEND_STRING PROPERTY", 1)),
            (link_depends, link_depends.replace(
                "LINK_DEPENDS", "LINK_OPTIONS", 1)),
            (link_depends, link_depends.replace(
                "WriteSparseEncodeEvidenceSidecar.cmake",
                "InjectedSparseSidecar.cmake", 1)),
            (post_build, post_build.replace(
                "TARGET bench_leopard2_sparse_encode POST_BUILD",
                "TARGET leopard POST_BUILD", 1)),
            (post_build, post_build.replace(
                "POST_BUILD", "PRE_LINK", 1)),
            (post_build, post_build.replace(
                "COMMAND ${CMAKE_COMMAND}", "COMMAND injected", 1)),
            (post_build, post_build.replace(
                "-DPRODUCTION_ARCHIVE=$<TARGET_FILE:leopard>",
                "-DPRODUCTION_ARCHIVE=$<TARGET_FILE:libleopard>", 1)),
            (post_build, post_build.replace(
                "-DBENCHMARK_OBJECT=$<TARGET_OBJECTS:"
                "leopard2_sparse_encode_benchmark_object>",
                "-DBENCHMARK_OBJECT=$<TARGET_OBJECTS:"
                "leopard2_sparse_encode_oracle_object>", 1)),
            (post_build, post_build.replace(
                "-DBUILD_ROOT=${CMAKE_CURRENT_BINARY_DIR}",
                "-DBUILD_ROOT=${CMAKE_CURRENT_SOURCE_DIR}", 1)),
            (post_build, post_build.replace(
                "WriteSparseEncodeEvidenceSidecar.cmake",
                "InjectedSparseSidecar.cmake", 1)),
            (post_build, post_build.replace("\n            VERBATIM)", ")", 1)),
        )
        for original, replacement in adversarial_mutations:
            with self.subTest(mutation=replacement.splitlines()[0]):
                self.assertNotEqual(original, replacement)
                mutated = self.cmake.replace(original, replacement, 1)
                self.assertNotEqual(mutated, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "generated/custom build extension|target property|"
                        "embedded or environment source variable|"
                        "compile/link property bypasses graph"):
                    self.resolve_text(
                        mutated, require_mutation_contract=True)

        field_guard = (
            "if(LEOPARD_ENABLE_GF8 AND LEOPARD_ENABLE_GF16)")
        self.assertEqual(2, self.cmake.count(field_guard))
        for replacement in (
                "if(LEOPARD_ENABLE_GF8)",
                "if(LEOPARD_ENABLE_GF8 OR LEOPARD_ENABLE_GF16)",
                "if(NOT LEOPARD_ENABLE_GF8 AND LEOPARD_ENABLE_GF16)"):
            with self.subTest(guard=replacement):
                mutated = self.cmake.replace(field_guard, replacement, 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "sparse evidence sidecar command guard drift"):
                    self.resolve_text(
                        mutated, require_mutation_contract=True)

        sentinel = "__LEO2_SPARSE_SIDECAR_COMMAND_SENTINEL__"
        reordered = self.cmake.replace(link_depends, sentinel, 1)
        reordered = reordered.replace(post_build, link_depends, 1)
        reordered = reordered.replace(sentinel, post_build, 1)
        self.assertNotEqual(reordered, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "security-sensitive CMake command order drift"):
            self.resolve_text(reordered, require_mutation_contract=True)

    def test_language_and_project_toolchain_mutations_are_rejected(self):
        mutations = (
            "enable_language(CXX)",
            "enable_language(CUDA)",
            "project(leopard LANGUAGES CXX CUDA)",
            "function(inject)\nenable_language(CUDA)\n"
            "endfunction()\ninject()",
            "macro(inject)\nproject(injected LANGUAGES CXX CUDA)\n"
            "endmacro()\ninject()",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "language|toolchain|unsupported CMake block"):
                    self.resolve(mutation)

    def test_generated_and_custom_build_extensions_are_rejected(self):
        mutations = (
            "add_custom_command(TARGET leopard POST_BUILD "
            "COMMAND injected)",
            "add_custom_command(OUTPUT Injected.obj COMMAND injected)",
            "add_custom_target(injected ALL COMMAND injected)",
            "add_dependencies(leopard injected)",
            "configure_file(Injected.cpp leopard2.cpp COPYONLY)",
            "file(GENERATE OUTPUT Injected.cpp CONTENT injected)",
            "execute_process(COMMAND injected)",
            "try_compile(RESULT build Injected.cpp)",
            "cmake_parse_arguments(LIB \"\" SOURCE_FILES \"\" "
            "SOURCE_FILES Injected.cpp)",
            "configure_package_config_file(Injected.cpp leopard2.cpp "
            "INSTALL_DESTINATION injected)",
            "function(inject)\n"
            "configure_package_config_file(Injected.cpp leopard2.cpp "
            "INSTALL_DESTINATION injected)\nendfunction()\ninject()",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "generated/custom build extension|"
                        "generated package configuration|build extension|"
                        "unsupported CMake block"):
                    self.resolve(mutation)

    def test_unmodeled_root_commands_are_rejected(self):
        mutations = (
            "exec_program(injected)",
            "write_file(leopard2.cpp injected)",
            "write_basic_package_version_file(injected.cmake VERSION 1.0)",
            "export(TARGETS leopard FILE injected.cmake)",
            "cmake_policy(SET CMP0001 OLD)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "unmodeled root CMake command"):
                    self.resolve(mutation)

    def test_bracket_comment_command_smuggling_is_rejected(self):
        mutation = """#[[
message(
]]
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
#[[
)
]]"""
        with self.assertRaisesRegex(
                ContractError, "compiler-control variable"):
            self.resolve(mutation)

        mutation = """set #[=[ ignored comment ]=]
(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)"""
        with self.assertRaisesRegex(
                ContractError, "compiler-control variable"):
            self.resolve(mutation)

    def test_bracket_argument_command_smuggling_is_rejected(self):
        mutation = """message([=[(]=])
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
message([=[)]=])"""
        self.assertEqual(
            ["message", "set", "message"],
            [name for name, unused in cmake_commands(mutation)])
        with self.assertRaisesRegex(
                ContractError, "bracket arguments are unsupported"):
            self.resolve(mutation)

    def test_regex_dialect_differences_cannot_make_mutations_unreachable(self):
        mutation = r'''if("a" MATCHES "\\a")
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
endif()'''
        with self.assertRaisesRegex(
                ContractError, "compiler-control variable"):
            self.resolve(mutation)

    def test_unknown_comparison_operands_remain_external_and_symbolic(self):
        conditions = (
            'USER_SWITCH STREQUAL "yes"',
            "USER_SWITCH EQUAL 7",
            'USER_SWITCH MATCHES "^yes$"',
        )
        for condition in conditions:
            with self.subTest(condition=condition):
                mutation = (
                    "if(" + condition + ")\n"
                    "set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)\n"
                    "endif()")
                with self.assertRaisesRegex(
                        ContractError, "compiler-control variable"):
                    self.resolve(mutation)

    def test_external_cache_and_environment_expansions_are_rejected(self):
        mutations = (
            'set(DEST CMAKE_VS_GLOBALS CACHE STRING "")\n'
            'string(CONCAT $CACHE{DEST} "VCTargetsPath=C:/evil")',
            'string(CONCAT $CACHE{DEST} "VCTargetsPath=C:/evil")',
            'get_filename_component($CACHE{DEST} injected ABSOLUTE)',
            'get_property($ENV{DEST} GLOBAL PROPERTY INJECTED)',
        )
        for mutation in mutations:
            with self.subTest(command=mutation.splitlines()[-1]):
                with self.assertRaisesRegex(
                        ContractError,
                        "cache assignment|cache/environment expansion"):
                    self.resolve(mutation)

    def test_cache_options_and_backend_enum_remain_external(self):
        mutations = (
            "if(NOT LEO2_BUILD_TESTS)\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(LEO2_ENABLE_CUDA)\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(LEO2_PORTABLE_ISA_RELEASE_AUDIT)\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(LEO2_BUILD_FUZZERS)\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(NOT LEO2_BUILD_BENCHMARKS)\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(LEO2_BACKEND_VARIANT STREQUAL \"scalar\")\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
            "if(LEOPARD_INSTALL_CMAKEDIR STREQUAL \"evil\")\n"
            "target_compile_options(leopard PRIVATE /GL)\nendif()",
        )
        for mutation in mutations:
            with self.subTest(condition=mutation.splitlines()[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation"):
                    self.resolve(mutation)

    def test_every_supported_x86_processor_spelling_is_proved(self):
        for processor in (
                "x86", "X86", "i386", "i486", "i586", "i686",
                "x86_64", "amd64", "AMD64"):
            with self.subTest(processor=processor):
                mutation = (
                    "if(CMAKE_SYSTEM_PROCESSOR STREQUAL \"" + processor +
                    "\")\ntarget_compile_options(leopard PRIVATE /GL)\n"
                    "endif()")
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved production target compile/link mutation"):
                    self.resolve(mutation)

    def test_cuda_tools_and_sources_cannot_enter_default_graph(self):
        mutations = (
            "find_program(CUDA_NVCC_EXECUTABLE nvcc REQUIRED)",
            "set(PROGRAM nvcc)\nfind_program(TOOL ${PROGRAM} REQUIRED)",
            "find_program(TOOL ptxas REQUIRED)",
            "find_program(TOOL nvlink REQUIRED)",
            "find_program(TOOL fatbinary REQUIRED)",
            "find_program(TOOL cuobjdump REQUIRED)",
            "add_library(cuda_always STATIC kernels.cu)",
            "add_library(cuda_always STATIC kernels.$<LOWER_CASE:CU>)",
            "set(CUDA_ALWAYS_SOURCE kernels.cuh)",
            "add_library(extra_cuda STATIC tests.cpp)\n"
            "target_link_libraries(extra_cuda PRIVATE "
            "C:/CUDA/lib/x64/cudart.lib)",
            "add_library(extra_cuda STATIC tests.cpp)\n"
            "target_include_directories(extra_cuda PRIVATE "
            "C:/Program_Files/NVIDIA_GPU_Computing_Toolkit/CUDA/include)",
            "add_library(extra_cuda STATIC tests.cpp)\n"
            "target_link_libraries(extra_cuda PRIVATE "
            "$<JOIN:cu,>dart.lib)",
            "add_library(extra_cuda STATIC tests.cpp)\nset(VENDOR CUDA)\n"
            "target_include_directories(extra_cuda PRIVATE "
            "C:/${VENDOR}/include)",
            "add_library(extra_cuda STATIC tests.cpp)\nset(VENDOR nvidia)\n"
            "target_link_libraries(extra_cuda PRIVATE "
            "C:/${VENDOR}/lib/foo.lib)",
            "add_library(extra_cuda STATIC tests.cpp)\n"
            "set(GEN \"$<JOIN:cu,>dart.lib\")\n"
            "target_link_libraries(extra_cuda PRIVATE ${GEN})",
            "add_library(extra_cuda STATIC tests/benchmark.cpp)\n"
            "set_source_files_properties(tests/benchmark.cpp "
            "PROPERTIES LANGUAGE CUDA)",
            "add_library(extra_cuda STATIC tests/benchmark.cpp)\n"
            "set_property(SOURCE tests/benchmark.cpp PROPERTY LANGUAGE CUDA)",
            "add_library(extra_cuda STATIC tests/benchmark.cpp)\n"
            "set_source_files_properties(tests/benchmark.cpp PROPERTIES "
            "VS_TOOL_OVERRIDE CudaCompile)",
            "add_library(extra_cuda STATIC tests/benchmark.cpp)\n"
            "set_source_files_properties(tests/benchmark.cpp PROPERTIES "
            "COMPILE_OPTIONS \"-x;cuda\")",
            "add_library(extra_cuda STATIC tests/benchmark.cpp)\n"
            "set_source_files_properties(tests/benchmark.cpp PROPERTIES "
            "COMPILE_OPTIONS $<JOIN:cu,da>)",
            "add_test(NAME cuda_probe COMMAND nvcc --version)",
            "set(TEST_TOOL nvcc)\n"
            "add_test(NAME cuda_probe COMMAND ${TEST_TOOL} --version)",
            "add_test(NAME cuda_probe COMMAND $<JOIN:nv,cc> --version)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "CUDA .*CPU-only graph|source generator expression|"
                        "dependency generator expression|"
                        "expanded target dependency generator expression|"
                        "unapproved tool discovery|unmodeled source-property|"
                        "unmodeled test generator expression"):
                    self.resolve(mutation)

    def test_unapproved_cache_declarations_cannot_hide_external_state(self):
        mutations = (
            'option(USER_SWITCH "external switch" OFF)',
            'set(USER_SWITCH OFF CACHE BOOL "external switch")',
            'set(USER_SWITCH)\nif(USER_SWITCH)\n'
            'set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)\nendif()',
            "unset(USER_SWITCH)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.splitlines()[0]):
                with self.assertRaisesRegex(
                        ContractError,
                        "cache option|cache assignment|compiler-control "
                        "variable|unset state"):
                    self.resolve(mutation)

    def test_install_execution_extensions_are_rejected(self):
        for mutation in (
                'install(CODE "execute_process(COMMAND injected)")',
                "install(SCRIPT injected.cmake)"):
            with self.subTest(command=mutation):
                with self.assertRaisesRegex(
                        ContractError,
                        "unapproved CMake install/package mutation"):
                    self.resolve(mutation)

    def test_arbitrary_target_properties_are_rejected(self):
        mutations = (
            "set_property(TARGET leopard PROPERTY LINKER_LANGUAGE CUDA)",
            "set_property(TARGET leopard PROPERTY EXCLUDE_FROM_ALL TRUE)",
            "set_property(TARGET leopard PROPERTY SOVERSION 99)",
            "set_target_properties(leopard2_codec_options_abi_test "
            "PROPERTIES LINKER_LANGUAGE CUDA)",
        )
        for mutation in mutations:
            with self.subTest(property=mutation):
                with self.assertRaisesRegex(ContractError, "property"):
                    self.resolve(mutation)

    def test_build_graph_command_shadowing_is_rejected(self):
        mutations = (
            "function(target_compile_options)\nendfunction()",
            "macro(target_link_libraries)\nendmacro()",
            "function(add_custom_command)\nendfunction()",
            "macro(target_sources)\nendmacro()",
            "macro(check_cxx_compiler_flag)\nendmacro()",
            "macro(find_package)\nendmacro()",
            "macro(option)\nendmacro()",
            "macro(cmake_dependent_option)\nendmacro()",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "command shadowing"):
                    self.resolve(mutation)

    def test_scoped_directory_properties_are_rejected(self):
        for mutation in (
                "function(inject)\nset_directory_properties(PROPERTIES "
                "COMPILE_OPTIONS /arch:AVX2)\nendfunction()\ninject()",
                "macro(inject)\nset_directory_properties(PROPERTIES "
                "LINK_OPTIONS /LTCG)\nendmacro()\ninject()"):
            with self.subTest(scope=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "unsupported CMake block"):
                    self.resolve(mutation)

    def test_required_production_mutations_cannot_be_removed_or_duplicated(self):
        commands = (
            """target_compile_definitions(leopard PRIVATE
        LEO2_DISABLE_SSSE3_CODEGEN=1
        LEO2_DISABLE_AVX2_CODEGEN=1)""",
            """target_include_directories(leopard PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)""",
            "target_link_libraries(leopard PUBLIC Threads::Threads)",
            """target_include_directories(leopard2_backend_ssse3 PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})""",
            """target_include_directories(leopard2_backend_avx2 PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR})""",
            "target_compile_options(leopard2_backend_avx2 PRIVATE /arch:AVX2)",
            "target_compile_definitions(leopard PRIVATE "
            "LEO2_HAVE_SSSE3_BACKEND=1)",
            "target_compile_definitions(leopard PRIVATE "
            "LEO2_HAVE_AVX2_BACKEND=1)",
            "target_link_libraries(leopard PUBLIC OpenMP::OpenMP_CXX)",
            "target_link_libraries(leopard PUBLIC \"${OpenMP_CXX_FLAGS}\")",
        )
        for command in commands:
            with self.subTest(command=command.split("(", 1)[0]):
                index = self.cmake.rfind(command)
                self.assertNotEqual(-1, index)
                text = self.cmake[:index] + self.cmake[index + len(command):]
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate production target mutation"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

        with self.assertRaisesRegex(
                ContractError,
                "missing or duplicate production target mutation"):
            self.resolve(
                "target_link_libraries(leopard PUBLIC Threads::Threads)")

    def test_required_trusted_cmake_commands_cannot_be_removed_or_duplicated(self):
        commands = (
            "cmake_minimum_required(VERSION 3.7)",
            "project(leopard)",
            "add_library(leopard STATIC ${LIB_SOURCE_FILES})",
            "add_library(libleopard ALIAS leopard)",
            "include(CMakeDependentOption)",
            "include(CheckCXXCompilerFlag)",
            "include(CMakePackageConfigHelpers)",
            "include(GNUInstallDirs)",
            "find_package(OpenMP)",
            "find_package(Threads REQUIRED)",
            'option(LEO2_BUILD_TESTS "Build Leopard2 correctness tests" ON)',
            'option(LEO2_BUILD_BENCHMARKS "Build Leopard benchmark programs" '
            'ON)',
            'option(LEO2_BUILD_ALLK_DIAGNOSTIC\n'
            '    "Build the source-attested all-K Leopard1 comparison '
            'benchmark" OFF)',
            'option(LEO2_BUILD_FUZZERS "Build Leopard2 libFuzzer targets" OFF)',
            'option(LEO2_ENABLE_CUDA "Build the optional Leopard2 CUDA backend" '
            'OFF)',
            'option(LEOPARD_ENABLE_GF8 "Include the GF(2^8) codec" ON)',
            'option(LEOPARD_ENABLE_GF16 "Include the GF(2^16) codec" ON)',
            'option(LEO2_PORTABLE_ISA_RELEASE_AUDIT\n'
            '    "Require the strict x86-64 portable-ISA release audit" OFF)',
            'option(LEO2_EXPERIMENT_K1_DECODE_VALIDATOR_INLINE\n'
            '    "Diagnostic: force-inline the K=1 decode buffer validator" '
            'OFF)',
            'string(TOLOWER "${LEO2_BACKEND_VARIANT}" '
            'LEO2_BACKEND_VARIANT_NORMALIZED)',
            'check_cxx_compiler_flag("/O2" CXX_FLAG_O2)',
            'check_cxx_compiler_flag("/Oy" CXX_FLAG_Oy)',
            'check_cxx_compiler_flag("/Zi" CXX_FLAG_Zi)',
            'check_cxx_compiler_flag("/W4" CXX_FLAG_W4)',
            'check_cxx_compiler_flag("/arch:AVX2" LEO2_FLAG_ARCH_AVX2)',
            'cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" '
            'ON "OPENMP_FOUND" OFF)',
            """configure_package_config_file(
    cmake/leopardConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake"
    INSTALL_DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")""",
            """install(TARGETS leopard
    EXPORT leopardTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")""",
            """install(FILES leopard.h leopard2.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")""",
            """install(FILES "${CMAKE_CURRENT_BINARY_DIR}/leopardConfig.cmake"
    DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")""",
            """install(EXPORT leopardTargets
    FILE leopardTargets.cmake
    NAMESPACE leopard::
    DESTINATION "${LEOPARD_INSTALL_CMAKEDIR}")""",
        )
        for command in commands:
            with self.subTest(command=command.split("(", 1)[0]):
                self.assertEqual(1, self.cmake.count(command))
                text = self.cmake.replace(command, "", 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate trusted CMake command|"
                        "locator provenance command guard drift|"
                        "unsupported conditional CMake source graph|"
                        "reachable CMake language enablement|"
                        "CMake bootstrap must begin"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

        with self.assertRaisesRegex(
                ContractError,
                "missing or duplicate trusted CMake command"):
            self.resolve("find_package(Threads REQUIRED)")

    def test_minimum_and_project_are_the_first_executable_commands(self):
        preambles = (
            """if(NOT MSVC)
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
endif()""",
            """if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
endif()""",
            """set(X yes)
if("X" STREQUAL "yes")
set(CMAKE_VS_GLOBALS VCTargetsPath=C:/evil)
endif()""",
        )
        for preamble in preambles:
            with self.subTest(first=preamble.splitlines()[0]):
                with self.assertRaisesRegex(
                        ContractError, "CMake bootstrap must begin"):
                    self.resolve_text(
                        preamble + "\n" + self.cmake,
                        require_mutation_contract=True)

    def test_trusted_cmake_command_arguments_and_guards_are_exact(self):
        marker = ('check_cxx_compiler_flag("/arch:AVX2" '
                  'LEO2_FLAG_ARCH_AVX2)')
        text = self.cmake.replace(
            marker,
            'check_cxx_compiler_flag("/arch:AVX" LEO2_FLAG_ARCH_AVX2)', 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError,
                "unapproved production compiler probe|"
                "missing or duplicate trusted CMake command"):
            self.resolve_text(text, require_mutation_contract=True)

        marker = "include(GNUInstallDirs)"
        replacement = (
            "if(COMMAND injected_gate)\n" + marker + "\nendif()")
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "trusted CMake command guard drift"):
            self.resolve_text(text, require_mutation_contract=True)

    def test_security_sensitive_cmake_command_order_is_exact(self):
        discovery = "find_package(OpenMP)\n"
        insertion = "endif(ENABLE_OPENMP)\n\nset(LEO2_X86_TARGET OFF)"
        self.assertEqual(1, self.cmake.count(discovery))
        self.assertEqual(1, self.cmake.count(insertion))
        text = self.cmake.replace(discovery, "", 1).replace(
            insertion,
            "endif(ENABLE_OPENMP)\nfind_package(OpenMP)\n\n"
            "set(LEO2_X86_TARGET OFF)", 1)
        with self.assertRaisesRegex(
                ContractError, "security-sensitive CMake command order"):
            self.resolve_text(text, require_mutation_contract=True)

    def test_required_protected_assignments_cannot_be_removed_or_duplicated(self):
        commands = (
            'set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" '
            'FORCE)',
            "set(CMAKE_CXX_STANDARD 11)",
            "set(CMAKE_BUILD_TYPE Release)",
            'set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")',
            'set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Oy")',
            'set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")',
            'set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")',
            'set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")',
            'set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} '
            '${OpenMP_EXE_LINKER_FLAGS}")',
            """set(LEOPARD_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"
    CACHE STRING "Install directory for Leopard CMake package files")""",
            """set(LEO2_BACKEND_VARIANT "auto" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512")""",
            """set(LEO2_BACKEND_VARIANT "${LEO2_BACKEND_VARIANT_NORMALIZED}" CACHE STRING
    "Diagnostic backend variant: auto, scalar, ssse3, avx2, or avx512" FORCE)""",
        )
        for command in commands:
            with self.subTest(variable=command.split("(", 1)[1].split()[0]):
                self.assertEqual(1, self.cmake.count(command))
                text = self.cmake.replace(command, "", 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "missing or duplicate protected CMake assignment|"
                        "may mutate source variable|"
                        "unresolved CMake source variable"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

        with self.assertRaisesRegex(
                ContractError,
                "missing or duplicate protected CMake assignment"):
            self.resolve("set(CMAKE_CXX_STANDARD 11)")

    def test_protected_assignment_guard_drift_is_rejected(self):
        marker = "set(CMAKE_CXX_STANDARD 11)"
        replacement = (
            "if(COMMAND injected_gate)\n" + marker + "\nendif()")
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError,
                "guards protected CMake assignment|"
                "protected CMake assignment guard drift"):
            self.resolve_text(text, require_mutation_contract=True)

    def test_protected_assignment_cache_and_scope_modifiers_are_exact(self):
        replacements = (
            ("set(CMAKE_CXX_STANDARD 11)",
             "set(CMAKE_CXX_STANDARD 11 PARENT_SCOPE)"),
            ('set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" '
             'FORCE)',
             'set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "")'),
            ("""set(LEOPARD_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"
    CACHE STRING "Install directory for Leopard CMake package files")""",
             """set(LEOPARD_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/leopard"
    CACHE STRING "Install directory for Leopard CMake package files" FORCE)"""),
        )
        for marker, replacement in replacements:
            with self.subTest(replacement=replacement.splitlines()[0]):
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "compiler-control variable mutation|PARENT_SCOPE"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

    def test_source_list_parent_scope_cannot_preserve_attacker_state(self):
        marker = "set(LIB_SOURCE_FILES\n"
        text = self.cmake.replace(
            marker, "set(LIB_SOURCE_FILES tests.cpp)\n" + marker, 1)
        text = text.replace(
            "        LeopardCommon.h)",
            "        LeopardCommon.h PARENT_SCOPE)", 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(ContractError, "PARENT_SCOPE"):
            self.resolve_text(text, require_mutation_contract=True)

    def test_computed_list_output_variables_cannot_rewrite_source_graph(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        operations = (
            "list(GET EVIL 0 LIB_SOURCE_FILES)",
            "list(LENGTH EVIL LIB_SOURCE_FILES)",
            'list(JOIN EVIL "" LIB_SOURCE_FILES)',
            "list(POP_BACK EVIL LIB_SOURCE_FILES)",
            "list(TRANSFORM EVIL TOUPPER OUTPUT_VARIABLE LIB_SOURCE_FILES)",
        )
        for operation in operations:
            with self.subTest(operation=operation.split("(", 1)[1].split()[0]):
                mutation = "set(EVIL tests.cpp)\n" + operation
                text = self.cmake.replace(
                    marker, mutation + "\n" + marker, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError, "unsupported CMake list operation"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

    def test_leopard_library_kind_and_definition_are_exact(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        replacements = (
            "add_library(leopard SHARED ${LIB_SOURCE_FILES})",
            "add_library(leopard MODULE ${LIB_SOURCE_FILES})",
            "add_library(leopard OBJECT ${LIB_SOURCE_FILES})",
            "add_library(leopard ${LIB_SOURCE_FILES})",
            "add_library(leopard STATIC EXCLUDE_FROM_ALL "
            "${LIB_SOURCE_FILES})",
        )
        for replacement in replacements:
            with self.subTest(definition=replacement):
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError, "exact STATIC library definition"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

        with self.assertRaisesRegex(
                ContractError,
                "missing or duplicate trusted CMake command"):
            self.resolve(marker)

    def test_legacy_target_is_alias_only(self):
        marker = "add_library(libleopard ALIAS leopard)"
        replacements = (
            "add_library(libleopard STATIC ${LIB_SOURCE_FILES})",
            "add_library(libleopard ALIAS leopard2_backend_ssse3)",
        )
        for replacement in replacements:
            with self.subTest(definition=replacement):
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "exact alias|missing or duplicate trusted"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

        with self.assertRaisesRegex(
                ContractError, "missing or duplicate trusted"):
            self.resolve_text(
                self.cmake.replace(marker, "", 1),
                require_mutation_contract=True)

    def test_production_mutation_guard_drift_is_rejected(self):
        marker = "target_link_libraries(leopard PUBLIC Threads::Threads)"
        replacements = (
            "if(COMMAND injected_gate)\n" + marker + "\nendif()",
            "if(CMAKE_VS_PLATFORM_NAME STREQUAL \"x64\")\n" + marker +
            "\nendif()",
            "if(FALSE)\n" + marker + "\nendif()",
        )
        for replacement in replacements:
            with self.subTest(guard=replacement.splitlines()[0]):
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "guards production target mutation|guard drift|"
                        "missing or duplicate production target mutation"):
                    self.resolve_text(
                        text, require_mutation_contract=True)

    def test_approved_link_flag_variable_cannot_be_locally_redirected(self):
        marker = ('target_link_libraries(leopard PUBLIC '
                  '"${OpenMP_CXX_FLAGS}")')
        mutation = "set(OpenMP_CXX_FLAGS linked_backend)\n" + marker
        text = self.cmake.replace(marker, mutation, 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError,
                "OpenMP_CXX_FLAGS.*mutation|mutation.*OpenMP_CXX_FLAGS"):
            self.resolve_text(text)

    def test_indirect_module_path_mutation_is_rejected(self):
        mutation = (
            "set(MODULE_DEST CMAKE_MODULE_PATH)\n"
            "string(APPEND ${MODULE_DEST} \";cmake/injected\")")
        with self.assertRaisesRegex(ContractError, "CMAKE_MODULE_PATH"):
            self.resolve(mutation)

    def test_variable_mutation_in_unsupported_block_is_rejected(self):
        mutations = (
            """set(SOURCE_DEST LIB_SOURCE_FILES)
macro(inject_sources)
    set(${SOURCE_DEST} Injected.cpp)
endmacro()
inject_sources()""",
            """set(SOURCE_DEST LIB_SOURCE_FILES)
function(inject_sources)
    set(${SOURCE_DEST} Injected.cpp PARENT_SCOPE)
endfunction()
inject_sources()""",
        )
        for mutation in mutations:
            with self.subTest(block=mutation.splitlines()[1]):
                with self.assertRaisesRegex(
                        ContractError,
                        "mutation in unsupported CMake block|"
                        "unsupported CMake block"):
                    self.resolve(mutation)

    def test_win32_only_source_variable_mutation_is_rejected(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        conditions = (
            'CMAKE_SYSTEM_PROCESSOR STREQUAL "x86"',
            'CMAKE_VS_PLATFORM_NAME STREQUAL "Win32"',
        )
        for condition in conditions:
            with self.subTest(condition=condition):
                mutation = (
                    "if(" + condition + ")\n"
                    "    list(APPEND LIB_SOURCE_FILES Win32Only.cpp)\n"
                    "endif()")
                text = self.cmake.replace(
                    marker, mutation + "\n" + marker, 1)
                with self.assertRaisesRegex(
                        ContractError, "conditional.*LIB_SOURCE_FILES"):
                    self.resolve_text(text)

    def test_conditional_leopard_definition_is_rejected(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        replacement = """if(NOT WIN32)
add_library(leopard STATIC ${LIB_SOURCE_FILES})
endif()"""
        text = self.cmake.replace(marker, replacement, 1)
        with self.assertRaisesRegex(
                ContractError, "add_library must be unconditional"):
            self.resolve_text(text)

    def test_contextual_leopard_definition_is_rejected(self):
        marker = "add_library(leopard STATIC ${LIB_SOURCE_FILES})"
        wrappers = (
            ("function(wrapper)", "endfunction()"),
            ("macro(wrapper)", "endmacro()"),
            ("foreach(item IN ITEMS one)", "endforeach()"),
            ("while(FALSE)", "endwhile()"),
        )
        for opening, closing in wrappers:
            with self.subTest(block=opening.split("(", 1)[0]):
                replacement = opening + "\n" + marker + "\n" + closing
                text = self.cmake.replace(marker, replacement, 1)
                with self.assertRaisesRegex(
                        ContractError,
                        "add_library must be unconditional|"
                        "unsupported CMake block"):
                    self.resolve_text(text)

    def test_contextual_target_property_source_mutation_is_rejected(self):
        mutations = (
            "set_property(TARGET leopard APPEND PROPERTY "
            "SOURCES Injected.cpp)",
            "set_target_properties(leopard PROPERTIES "
            "SOURCES Injected.cpp)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                text = (
                    "function(inject_sources)\n" + mutation +
                    "\nendfunction()\ninject_sources()")
                with self.assertRaisesRegex(
                        ContractError, "unsupported CMake block"):
                    self.resolve(text)

    def test_attached_object_requires_msvc_reachable_definition(self):
        texts = ("""
set(LIB_SOURCE_FILES Base.cpp Base.h)
add_library(leopard STATIC ${LIB_SOURCE_FILES})
if(NOT WIN32)
    add_library(hidden_backend OBJECT HiddenBackend.cpp)
endif()
target_sources(leopard PRIVATE $<TARGET_OBJECTS:hidden_backend>)
""", """
set(LIB_SOURCE_FILES Base.cpp Base.h)
add_library(leopard STATIC ${LIB_SOURCE_FILES})
if(MSVC)
    if(FALSE)
        add_library(hidden_backend OBJECT HiddenBackend.cpp)
    endif()
endif()
target_sources(leopard PRIVATE $<TARGET_OBJECTS:hidden_backend>)
""")
        for text in texts:
            with self.subTest(definition=text.splitlines()[3].strip()):
                with self.assertRaisesRegex(
                        ContractError, "no MSVC-reachable definition"):
                    self.resolve_text(text)

    def test_current_msvc_object_branch_is_required(self):
        text = self.cmake.replace("elseif(MSVC)", "elseif(FALSE)", 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "no MSVC-reachable definition"):
            self.resolve_text(text)

    def test_impossible_msvc_compiler_conjunction_is_rejected(self):
        impossible_branches = (
            'elseif(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")',
            'elseif(MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")',
            'elseif(MSVC AND NOT MSVC)',
        )
        for replacement in impossible_branches:
            with self.subTest(condition=replacement):
                text = self.cmake.replace("elseif(MSVC)", replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError,
                        "no MSVC-reachable|unsupported conditional"):
                    self.resolve_text(text)

    def test_target_objects_attachment_must_reach_msvc(self):
        marker = """    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_ssse3>)"""
        for condition in ("FALSE", "NOT WIN32"):
            with self.subTest(condition=condition):
                replacement = (
                    "    if(" + condition + ")\n" + marker +
                    "\n    endif()")
                text = self.cmake.replace(marker, replacement, 1)
                self.assertNotEqual(text, self.cmake)
                with self.assertRaisesRegex(
                        ContractError, "no MSVC-reachable"):
                    self.resolve_text(text)

    def test_target_objects_attachment_rejects_unmodeled_condition(self):
        marker = """    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_ssse3>)"""
        replacement = (
            "    if(COMMAND injected_attachment_gate)\n" + marker +
            "\n    endif()")
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "unsupported conditional"):
            self.resolve_text(text)

    def test_target_objects_else_preserves_unsupported_condition(self):
        marker = """    target_sources(leopard PRIVATE
        $<TARGET_OBJECTS:leopard2_backend_ssse3>)"""
        replacement = (
            "    if(COMMAND injected_attachment_gate)\n"
            "        message(FATAL_ERROR \"unreachable\")\n"
            "    else()\n" + marker + "\n    endif()")
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        with self.assertRaisesRegex(
                ContractError, "unsupported conditional"):
            self.resolve_text(text)

    def test_conditional_direct_source_attachment_is_rejected(self):
        mutation = """
if(WIN32)
    target_sources(leopard PRIVATE Conditional.cpp)
endif()
"""
        with self.assertRaisesRegex(
                ContractError, "conditional direct leopard"):
            self.resolve(mutation)

    def test_source_properties_on_production_source_are_rejected(self):
        mutations = (
            "set_source_files_properties(leopard2.cpp "
            "PROPERTIES HEADER_FILE_ONLY TRUE)",
            "set_property(SOURCE leopard2.cpp "
            "PROPERTY HEADER_FILE_ONLY TRUE)",
            "set_source_files_properties(leopard2.cpp "
            "PROPERTIES GENERATED TRUE)",
            "set_property(SOURCE leopard2.cpp "
            "PROPERTY COMPILE_OPTIONS /arch:AVX2)",
            "set_property(SOURCE Leopard2BackendSSSE3.cpp "
            "PROPERTY COMPILE_DEFINITIONS LEO2_DISABLE_SSSE3_CODEGEN=0)",
            "set_source_files_properties(leopardff8.cpp "
            "PROPERTIES COMPILE_OPTIONS /arch:AVX2)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "source properties affect production"):
                    self.resolve(mutation)

        with self.assertRaisesRegex(
                ContractError, "duplicate production source attachment"):
            self.resolve(
                "target_sources(leopard PRIVATE leopardff8.cpp)")

    def test_target_sources_property_mutation_is_rejected(self):
        mutations = (
            "set_property(TARGET leopard APPEND PROPERTY "
            "SOURCES Injected.cpp)",
            "set_target_properties(leopard PROPERTIES "
            "SOURCES Injected.cpp)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "SOURCES property"):
                    self.resolve(mutation)

    def test_conditional_generator_expression_is_rejected(self):
        with self.assertRaisesRegex(ContractError, "generator expression"):
            self.resolve(
                "target_sources(leopard PRIVATE "
                "\"$<$<BOOL:1>:Conditional.cpp>\")")

    def test_cuda_source_or_header_is_rejected(self):
        for path in ("accidental.cu", "accidental.cuh"):
            with self.subTest(path=path):
                with self.assertRaisesRegex(ContractError, "CUDA source"):
                    self.resolve(
                        "target_sources(leopard PRIVATE " + path + ")")

    def test_path_escape_is_rejected(self):
        for path in ("../Escape.cpp", "sub/../Escape.cpp", "C:/Escape.cpp"):
            with self.subTest(path=path):
                with self.assertRaisesRegex(
                        ContractError, "repository-relative|drive path"):
                    self.resolve(
                        "target_sources(leopard PRIVATE " + path + ")")


class MSBuildToolchainMutationTest(unittest.TestCase):
    @staticmethod
    def fresh_tree():
        return ET.parse(str(PROJECT))

    def test_unapproved_or_duplicate_import_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        projects = (
            "..\\Injected.props",
            "$(VCTargetsPath)\\Microsoft.Cpp.props",
        )
        for imported_project in projects:
            with self.subTest(project=imported_project):
                tree = self.fresh_tree()
                imported = ET.Element(namespace + "Import")
                imported.set("Project", imported_project)
                tree.getroot().append(imported)
                with self.assertRaisesRegex(ContractError, "imports differ"):
                    validate_visual_studio_project(tree)

    def test_compiler_tool_override_family_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        for property_name in (
                "CLToolExe", "cltoolpath", "VCToolsPath", "VCTargetsPath"):
            with self.subTest(property=property_name):
                tree = self.fresh_tree()
                group = ET.SubElement(tree.getroot(), namespace + "PropertyGroup")
                override = ET.SubElement(group, namespace + property_name)
                override.text = "injected-compiler.exe"
                with self.assertRaisesRegex(
                        ContractError, "tool override|noncanonical"):
                    validate_visual_studio_project(tree)

    def test_custom_execution_logic_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        target = ET.SubElement(tree.getroot(), namespace + "Target")
        target.set("BeforeTargets", "ClCompile")
        with self.assertRaisesRegex(
                ContractError, "execution logic|unapproved root MSBuild"):
            validate_visual_studio_project(tree)

    def test_build_event_and_custom_build_commands_are_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        for mutation in ("post-build", "custom-build"):
            with self.subTest(mutation=mutation):
                tree = self.fresh_tree()
                if mutation == "post-build":
                    command = tree.find(
                        ".//msb:ItemDefinitionGroup/msb:PostBuildEvent/"
                        "msb:Command", NS)
                    command.text = "injected-compiler.exe"
                else:
                    group = ET.SubElement(
                        tree.getroot(), namespace + "ItemGroup")
                    custom = ET.SubElement(group, namespace + "CustomBuild")
                    custom.set("Include", "injected.targets")
                    command = ET.SubElement(custom, namespace + "Command")
                    command.text = "injected-compiler.exe"
                with self.assertRaisesRegex(
                        ContractError, "build-event|custom MSBuild build"):
                    validate_visual_studio_project(tree)

    def test_configuration_compile_tool_contract_rejects_extra_inputs(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        compile_tool = tree.find(
            ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)
        include_dirs = ET.SubElement(
            compile_tool, namespace + "AdditionalIncludeDirectories")
        include_dirs.text = "..\\injected;%(AdditionalIncludeDirectories)"
        with self.assertRaisesRegex(
                ContractError, "ClCompile.*exact contract"):
            validate_visual_studio_project(tree)

    def test_foreign_namespace_cannot_satisfy_exact_tool_properties(self):
        tree = self.fresh_tree()
        sdl_check = tree.find(
            ".//msb:ItemDefinitionGroup/msb:ClCompile/msb:SDLCheck", NS)
        sdl_check.tag = "{urn:foreign}SDLCheck"
        with self.assertRaisesRegex(
                ContractError, "foreign or missing MSBuild XML namespace"):
            validate_visual_studio_project(tree)

    def test_unapproved_build_bearing_item_types_are_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        for item_name in (
                "ProjectReference", "ResourceCompile", "Object",
                "PackageReference", "Reference", "Midl"):
            with self.subTest(item=item_name):
                tree = self.fresh_tree()
                group = ET.SubElement(tree.getroot(), namespace + "ItemGroup")
                item = ET.SubElement(group, namespace + item_name)
                item.set("Include", "injected")
                with self.assertRaisesRegex(
                        ContractError, "ItemGroup topology|ItemGroup contents"):
                    validate_visual_studio_project(tree)

    def test_import_and_toolchain_input_properties_are_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        for property_name in (
                "UserRootDir", "MSBuildUserExtensionsPath", "Platform",
                "Configuration", "IncludePath", "CLToolArchitecture",
                "BuildDependsOn", "ForceImportAfterCppTargets",
                "CustomAfterMicrosoftCommonTargets"):
            with self.subTest(property=property_name):
                tree = self.fresh_tree()
                group = ET.SubElement(
                    tree.getroot(), namespace + "PropertyGroup")
                value = ET.SubElement(group, namespace + property_name)
                value.text = "injected"
                with self.assertRaisesRegex(
                        ContractError,
                        "override is forbidden|hook is forbidden|"
                        "PropertyGroup topology"):
                    validate_visual_studio_project(tree)

    def test_root_evaluation_phase_order_is_exact(self):
        mutations = ("definitions-before-props", "sheets-after-definitions")
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                tree = self.fresh_tree()
                root = tree.getroot()
                children = list(root)
                if mutation == "definitions-before-props":
                    definitions = [node for node in children
                                   if xml_local_name(node) ==
                                   "ItemDefinitionGroup"]
                    for node in definitions:
                        root.remove(node)
                    cpp_props = next(
                        index for index, node in enumerate(list(root))
                        if (xml_local_name(node) == "Import" and
                            node.attrib.get("Project", "").endswith(
                                "Microsoft.Cpp.props")))
                    for offset, node in enumerate(definitions):
                        root.insert(cpp_props + offset, node)
                else:
                    sheets = [node for node in children
                              if (xml_local_name(node) == "ImportGroup" and
                                  node.attrib.get("Label") ==
                                  "PropertySheets")]
                    for node in sheets:
                        root.remove(node)
                    targets = next(
                        index for index, node in enumerate(list(root))
                        if (xml_local_name(node) == "Import" and
                            node.attrib.get("Project", "").endswith(
                                "Microsoft.Cpp.targets")))
                    for offset, node in enumerate(sheets):
                        root.insert(targets + offset, node)
                with self.assertRaisesRegex(
                        ContractError, "root evaluation phase/order"):
                    validate_visual_studio_project(tree)

    def test_librarian_archive_inputs_are_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        definition = tree.findall(
            ".//msb:ItemDefinitionGroup", NS)[0]
        librarian = ET.SubElement(definition, namespace + "Lib")
        dependencies = ET.SubElement(
            librarian, namespace + "AdditionalDependencies")
        dependencies.text = (
            "..\\Injected.obj;%(AdditionalDependencies)")
        with self.assertRaisesRegex(
                ContractError, "librarian archive input"):
            validate_visual_studio_project(tree)

    def test_unscoped_platform_toolset_override_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        group = ET.Element(namespace + "PropertyGroup")
        toolset = ET.SubElement(group, namespace + "PlatformToolset")
        toolset.text = "LLVM-vs2014"
        root = tree.getroot()
        cpp_props = next(
            index for index, node in enumerate(list(root))
            if (xml_local_name(node) == "Import" and
                node.attrib.get("Project", "").endswith(
                    "Microsoft.Cpp.props")))
        root.insert(cpp_props, group)
        with self.assertRaisesRegex(ContractError, "PlatformToolset"):
            validate_visual_studio_project(tree)

    def test_case_insensitive_platform_toolset_override_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        group = ET.SubElement(tree.getroot(), namespace + "PropertyGroup")
        toolset = ET.SubElement(group, namespace + "platformtoolset")
        toolset.text = "LLVM-vs2014"
        with self.assertRaisesRegex(ContractError, "PlatformToolset"):
            validate_visual_studio_project(tree)

    def test_unscoped_configuration_type_override_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        for spelling in ("ConfigurationType", "configurationtype"):
            with self.subTest(spelling=spelling):
                tree = self.fresh_tree()
                group = ET.Element(namespace + "PropertyGroup")
                configuration_type = ET.SubElement(
                    group, namespace + spelling)
                configuration_type.text = "Application"
                root = tree.getroot()
                cpp_props = next(
                    index for index, node in enumerate(list(root))
                    if (xml_local_name(node) == "Import" and
                        node.attrib.get("Project", "").endswith(
                            "Microsoft.Cpp.props")))
                root.insert(cpp_props, group)
                with self.assertRaisesRegex(
                        ContractError, "ConfigurationType"):
                    validate_visual_studio_project(tree)

    def test_case_variants_cannot_hide_controlled_msbuild_elements(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"

        def hidden_source(tree):
            group = tree.findall(".//msb:ItemGroup", NS)[-1]
            node = ET.SubElement(group, namespace + "clcompile")
            node.set("Include", "..\\tests\\benchmark.cpp")

        def hidden_metadata(tree, name, value):
            compile_node = tree.findall(
                ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
            node = ET.SubElement(compile_node, namespace + name)
            node.text = value

        def hidden_item_definition(tree):
            definition = tree.findall(
                ".//msb:ItemDefinitionGroup", NS)[0]
            compile_node = ET.SubElement(
                definition, namespace + "clcompile")
            options = ET.SubElement(
                compile_node, namespace + "additionaloptions")
            options.text = "/GL %(AdditionalOptions)"

        mutations = (
            ("source", hidden_source),
            ("definitions", lambda tree: hidden_metadata(
                tree, "preprocessordefinitions",
                "LEO2_BACKEND_FORCE_AVX2=1;%(PreprocessorDefinitions)")),
            ("options", lambda tree: hidden_metadata(
                tree, "additionaloptions", "/GL %(AdditionalOptions)")),
            ("runtime", lambda tree: hidden_metadata(
                tree, "runtimeLibrary", "MultiThreadedDLL")),
            ("item-definition", hidden_item_definition),
        )
        for label, mutate in mutations:
            with self.subTest(mutation=label):
                tree = self.fresh_tree()
                mutate(tree)
                with self.assertRaisesRegex(
                        ContractError, "noncanonical MSBuild element casing"):
                    validate_visual_studio_project(tree)

    def test_configuration_groups_must_precede_cpp_props(self):
        tree = self.fresh_tree()
        root = tree.getroot()
        groups = root.findall(
            "msb:PropertyGroup[@Label='Configuration']", NS)
        for group in groups:
            root.remove(group)
            root.append(group)
        with self.assertRaisesRegex(ContractError, "precede Microsoft.Cpp.props"):
            validate_visual_studio_project(tree)

    def test_non_vs2015_tools_version_is_rejected(self):
        tree = self.fresh_tree()
        tree.getroot().set("ToolsVersion", "15.0")
        with self.assertRaisesRegex(ContractError, "VS2015 contract"):
            validate_visual_studio_project(tree)


class MSBuildConditionMutationTest(unittest.TestCase):
    @staticmethod
    def fresh_tree():
        return ET.parse(str(PROJECT))

    def test_duplicate_normalized_condition_is_rejected(self):
        for xpath, label in (
                (".//msb:PropertyGroup[@Label='Configuration']", "property"),
                (".//msb:ItemDefinitionGroup", "definition")):
            with self.subTest(group=label):
                tree = self.fresh_tree()
                node = tree.findall(xpath, NS)[0]
                duplicate = copy.deepcopy(node)
                duplicate.set(
                    "Condition",
                    " '$(Configuration)|$(Platform)' == 'Debug|Win32' ")
                tree.getroot().append(duplicate)
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_project_configuration_group_and_item_must_be_unconditional(self):
        mutations = ("group", "item")
        for mutation in mutations:
            with self.subTest(scope=mutation):
                tree = self.fresh_tree()
                group = tree.find(
                    ".//msb:ItemGroup[@Label='ProjectConfigurations']", NS)
                node = (group if mutation == "group" else
                        group.find("msb:ProjectConfiguration", NS))
                node.set("Condition", "false")
                with self.assertRaisesRegex(
                        ContractError,
                        "ProjectConfigurations|ProjectConfiguration attributes|"
                        "unapproved MSBuild ItemGroup"):
                    validate_visual_studio_project(tree)

    def test_project_configuration_transform_outside_group_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        group = ET.SubElement(tree.getroot(), namespace + "ItemGroup")
        transformed = ET.SubElement(
            group, namespace + "ProjectConfiguration")
        transformed.set("Remove", "Debug|Win32")
        with self.assertRaisesRegex(
                ContractError,
                "outside the canonical group|ItemGroup topology"):
            validate_visual_studio_project(tree)

    def test_missing_condition_is_rejected(self):
        for xpath, label in (
                (".//msb:PropertyGroup[@Label='Configuration']", "property"),
                (".//msb:ItemDefinitionGroup", "definition")):
            with self.subTest(group=label):
                tree = self.fresh_tree()
                node = tree.findall(xpath, NS)[0]
                tree.getroot().remove(node)
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_semantic_condition_whitespace_is_rejected(self):
        tree = self.fresh_tree()
        node = tree.findall(
            ".//msb:PropertyGroup[@Label='Configuration']", NS)[0]
        node.set(
            "Condition",
            "'$(Configuration)|$(Platform)'==' Debug | Win32 '")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_false_item_definition_or_compile_condition_is_rejected(self):
        mutations = ("definition", "compile")
        for mutation in mutations:
            with self.subTest(scope=mutation):
                tree = self.fresh_tree()
                definition = tree.findall(
                    ".//msb:ItemDefinitionGroup", NS)[0]
                if mutation == "definition":
                    definition.set("Condition", "false")
                else:
                    definition.find("msb:ClCompile", NS).set(
                        "Condition", "false")
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_nested_item_definition_group_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        definition = tree.findall(".//msb:ItemDefinitionGroup", NS)[0]
        tree.getroot().remove(definition)
        choose = ET.SubElement(tree.getroot(), namespace + "Choose")
        when = ET.SubElement(choose, namespace + "When")
        when.set("Condition", "false")
        when.append(definition)
        with self.assertRaisesRegex(
                ContractError,
                "direct Project children|unapproved root MSBuild"):
            validate_visual_studio_project(tree)


class MSBuildPerFileMutationTest(unittest.TestCase):
    @staticmethod
    def fresh_tree():
        return ET.parse(str(PROJECT))

    @staticmethod
    def compile_item(tree, filename):
        for node in tree.findall(".//msb:ClCompile[@Include]", NS):
            if project_path(PROJECT, node.attrib["Include"]) == filename:
                return node
        raise AssertionError("compile item not found: " + filename)

    @staticmethod
    def add_metadata(node, name, value):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        child = ET.SubElement(node, namespace + name)
        child.text = value

    @classmethod
    def source_group(cls, tree, filename="leopard2.cpp"):
        node = cls.compile_item(tree, filename)
        for group in tree.getroot().findall("msb:ItemGroup", NS):
            if node in list(group):
                return group
        raise AssertionError("source ItemGroup not found")

    def test_per_file_backend_force_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "leopard2.cpp")
        self.add_metadata(
            node, "PreprocessorDefinitions",
            "LEO2_BACKEND_FORCE_AVX2=1;%(PreprocessorDefinitions)")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_per_file_definition_replacement_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "leopard2.cpp")
        self.add_metadata(node, "PreprocessorDefinitions", "LOCAL_ONLY=1")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_per_file_enhanced_isa_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "leopard2.cpp")
        self.add_metadata(
            node, "EnableEnhancedInstructionSet",
            "AdvancedVectorExtensions2")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_per_file_undefine_isolation_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "leopard2.cpp")
        self.add_metadata(
            node, "UndefinePreprocessorDefinitions",
            "LEO2_DISABLE_AVX2_CODEGEN")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_duplicate_avx2_compile_item_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "Leopard2BackendAVX2.cpp")
        for group in tree.getroot().findall("msb:ItemGroup", NS):
            if node in list(group):
                group.append(copy.deepcopy(node))
                break
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_configuration_undefine_isolation_is_rejected(self):
        tree = self.fresh_tree()
        node = tree.findall(
            ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
        self.add_metadata(
            node, "UndefinePreprocessorDefinitions",
            "LEO2_DISABLE_AVX2_CODEGEN")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_configuration_option_backend_force_is_rejected(self):
        tree = self.fresh_tree()
        node = tree.findall(
            ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
        self.add_metadata(
            node, "AdditionalOptions",
            "/DLEO2_BACKEND_FORCE_AVX2=1 %(AdditionalOptions)")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_excluded_from_build_is_rejected_globally_and_conditionally(self):
        for condition in (None,
                          "'$(Configuration)|$(Platform)'=='Release|x64'"):
            with self.subTest(condition=condition or "global"):
                tree = self.fresh_tree()
                node = self.compile_item(tree, "leopard2.cpp")
                namespace = (
                    "{http://schemas.microsoft.com/developer/msbuild/2003}")
                excluded = ET.SubElement(node, namespace + "ExcludedFromBuild")
                excluded.text = "true"
                if condition:
                    excluded.set("Condition", condition)
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_configuration_excluded_from_build_is_rejected(self):
        tree = self.fresh_tree()
        node = tree.findall(
            ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
        self.add_metadata(node, "ExcludedFromBuild", "true")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_clcompile_remove_is_rejected(self):
        tree = self.fresh_tree()
        group = self.source_group(tree)
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        removed = ET.SubElement(group, namespace + "ClCompile")
        removed.set("Remove", "..\\leopard2.cpp")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_clcompile_update_is_rejected(self):
        tree = self.fresh_tree()
        group = self.source_group(tree)
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        updated = ET.SubElement(group, namespace + "ClCompile")
        updated.set("Update", "..\\leopard2.cpp")
        self.add_metadata(
            updated, "AdditionalOptions",
            "/arch:AVX2 %(AdditionalOptions)")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_source_group_item_and_metadata_conditions_are_rejected(self):
        mutations = ("group", "item", "metadata")
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                tree = self.fresh_tree()
                node = self.compile_item(tree, "leopard2.cpp")
                condition = "'$(Configuration)'=='Release'"
                if mutation == "group":
                    self.source_group(tree).set("Condition", condition)
                elif mutation == "item":
                    node.set("Condition", condition)
                else:
                    namespace = (
                        "{http://schemas.microsoft.com/developer/msbuild/2003}")
                    metadata = ET.SubElement(
                        node, namespace + "AdditionalOptions")
                    metadata.text = "%(AdditionalOptions)"
                    metadata.set("Condition", condition)
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_nested_source_item_group_is_rejected(self):
        tree = self.fresh_tree()
        group = self.source_group(tree)
        root = tree.getroot()
        root.remove(group)
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        choose = ET.SubElement(root, namespace + "Choose")
        when = ET.SubElement(choose, namespace + "When")
        when.set("Condition", "'$(Configuration)'=='Release'")
        when.append(group)
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_unmodeled_source_metadata_is_rejected(self):
        tree = self.fresh_tree()
        node = self.compile_item(tree, "leopard2.cpp")
        self.add_metadata(node, "CompileAs", "CompileAsC")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)

    def test_quoted_backend_flags_are_rejected(self):
        options = (
            '/D"LEO2_BACKEND_FORCE_AVX2=1"',
            '/U"LEO2_DISABLE_AVX2_CODEGEN"',
        )
        for scope in ("configuration", "per-file"):
            for option in options:
                with self.subTest(scope=scope, option=option):
                    tree = self.fresh_tree()
                    if scope == "configuration":
                        node = tree.findall(
                            ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
                    else:
                        node = self.compile_item(tree, "leopard2.cpp")
                    self.add_metadata(
                        node, "AdditionalOptions",
                        option + " %(AdditionalOptions)")
                    with self.assertRaises(ContractError):
                        validate_visual_studio_project(tree)

    def test_unresolved_msbuild_property_is_rejected_in_controlled_fields(self):
        mutations = ("config-option", "file-option",
                     "config-definition", "file-definition",
                     "config-isa", "file-isa")
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                tree = self.fresh_tree()
                if mutation.startswith("config"):
                    node = tree.findall(
                        ".//msb:ItemDefinitionGroup/msb:ClCompile", NS)[0]
                else:
                    node = self.compile_item(tree, "leopard2.cpp")
                if mutation.endswith("option"):
                    self.add_metadata(
                        node, "AdditionalOptions",
                        "$(HiddenISA) %(AdditionalOptions)")
                elif mutation.endswith("isa"):
                    self.add_metadata(
                        node, "EnableEnhancedInstructionSet",
                        "$(HiddenISA)")
                elif mutation == "config-definition":
                    definitions = node.find(
                        "msb:PreprocessorDefinitions", NS)
                    definitions.text += ";$(InjectedBackendMacro)"
                else:
                    self.add_metadata(
                        node, "PreprocessorDefinitions",
                        "$(InjectedBackendMacro);%(PreprocessorDefinitions)")
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_wpo_option_forms_are_rejected(self):
        for option in ("/GL", "/LTCG", "/LTCG:INCREMENTAL", "-flto=thin"):
            with self.subTest(option=option):
                tree = self.fresh_tree()
                link = tree.findall(
                    ".//msb:ItemDefinitionGroup/msb:Link", NS)[0]
                self.add_metadata(
                    link, "AdditionalOptions",
                    option + " %(AdditionalOptions)")
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_non_avx2_additional_options_are_rejected(self):
        for option in ("@hidden.rsp", "/FIhidden.h"):
            with self.subTest(option=option):
                tree = self.fresh_tree()
                link = tree.findall(
                    ".//msb:ItemDefinitionGroup/msb:Link", NS)[0]
                self.add_metadata(
                    link, "AdditionalOptions",
                    option + " %(AdditionalOptions)")
                with self.assertRaises(ContractError):
                    validate_visual_studio_project(tree)

    def test_wpo_property_override_is_rejected(self):
        tree = self.fresh_tree()
        link = tree.findall(".//msb:ItemDefinitionGroup/msb:Link", NS)[0]
        self.add_metadata(
            link, "LinkTimeCodeGeneration", "UseLinkTimeCodeGeneration")
        with self.assertRaises(ContractError):
            validate_visual_studio_project(tree)


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]))
    sys.exit(0 if result.wasSuccessful() else 1)
