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
            newline = text.find("\n", index)
            index = length if newline < 0 else newline + 1
            continue
        match = re.match(r"[A-Za-z_][A-Za-z0-9_]*", text[index:])
        if not match:
            index += 1
            continue
        name = match.group(0).lower()
        cursor = index + len(match.group(0))
        while cursor < length and text[cursor].isspace():
            cursor += 1
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
                newline = text.find("\n", cursor)
                cursor = length if newline < 0 else newline
                continue
            elif not quoted and char == "(":
                depth += 1
            elif not quoted and char == ")":
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
    quoted = False
    escaped = False

    def flush():
        if token:
            tokens.append("".join(token))
            del token[:]

    index = 0
    while index < len(body):
        char = body[index]
        if escaped:
            token.append(char)
            escaped = False
        elif char == "\\":
            escaped = True
        elif char == '"':
            quoted = not quoted
        elif not quoted and char == "#":
            flush()
            newline = body.find("\n", index)
            index = len(body) if newline < 0 else newline
            continue
        elif not quoted and (char.isspace() or char == ";"):
            flush()
        else:
            token.append(char)
        index += 1
    if quoted or escaped:
        raise ContractError("unterminated quote or escape in CMake arguments")
    flush()
    return tokens


def cmake_condition_tokens(body):
    """Tokenize an if() body while retaining quotes and grouping."""
    tokens = []
    token = []
    quoted = False
    token_was_quoted = False
    escaped = False

    def flush():
        if token:
            tokens.append(("".join(token), token_was_quoted))
            del token[:]

    index = 0
    while index < len(body):
        char = body[index]
        if escaped:
            token.append(char)
            escaped = False
        elif char == "\\":
            escaped = True
        elif char == '"':
            quoted = not quoted
            token_was_quoted = True
        elif not quoted and char == "#":
            flush()
            newline = body.find("\n", index)
            index = len(body) if newline < 0 else newline
            continue
        elif not quoted and char in "()":
            flush()
            tokens.append((char, False))
            token_was_quoted = False
        elif not quoted and (char.isspace() or char == ";"):
            flush()
            token_was_quoted = False
        else:
            token.append(char)
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

    def __init__(self, text, processor="AMD64", pointer_size="8",
                 platform_name="x64"):
        self.raw_commands = cmake_commands(text)
        self.commands = [(name, cmake_tokens(body))
                         for name, body in self.raw_commands]
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
        self.attachments = {}
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
        try:
            return BOOL_TRUE if float(text) != 0 else BOOL_FALSE
        except ValueError:
            return BOOL_TRUE

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
            reasons = ()
            if (re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", text) and
                    text.upper() not in constants):
                reasons = (
                    "unmodeled unquoted CMake conditional operand: " + text,)
            return [(BOOL_TRUE, (text,), reasons)]
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
            try:
                matched = re.search(right_text, left_text) is not None
            except re.error as error:
                raise ContractError(
                    "invalid CMake MATCHES expression: " + str(error))
            return BOOL_TRUE if matched else BOOL_FALSE
        if operation == "EQUAL":
            try:
                matched = float(left_text) == float(right_text)
            except ValueError:
                matched = False
            return BOOL_TRUE if matched else BOOL_FALSE
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
                if operation not in {"STREQUAL", "MATCHES", "EQUAL"}:
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
                stack.append({"type": command})
                unsupported_depth += 1
                continue
            if command in ends:
                if not stack or stack[-1]["type"] != ends[command]:
                    raise ContractError("unbalanced CMake " + command)
                stack.pop()
                unsupported_depth -= 1
                continue

            if unsupported_depth:
                if (command == "add_library" and tokens and
                        tokens[0] == "libleopard"):
                    raise ContractError(
                        "libleopard add_library must be unconditional")
                if command in {
                        "add_library", "target_sources",
                        "set_source_files_properties", "set_property",
                        "set_target_properties"}:
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
                        "cmake_language", "load_command"}:
                    raise ContractError(
                        "graph extension in unsupported CMake block: " + command)
                continue

            if command == "include":
                if len(tokens) != 1 or tokens[0] not in approved_includes:
                    raise ContractError(
                        "unapproved CMake graph include: " + " ".join(tokens))
                continue
            if command == "find_package":
                if tuple(tokens) not in approved_packages:
                    raise ContractError(
                        "unapproved CMake package graph import: " +
                        " ".join(tokens))
                continue
            if command in {"add_subdirectory", "subdirs"}:
                raise ContractError(
                    "CMake " + command + " requires recursive graph proof")
            if command in {"cmake_language", "load_command"}:
                raise ContractError(
                    "CMake dynamic command execution is unsupported: " + command)

            if command == "set" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH can redirect approved graph includes")
                raw_values = list(tokens[1:])
                upper = [value.upper() for value in raw_values]
                if "CACHE" in upper:
                    raw_values = raw_values[:upper.index("CACHE")]
                if raw_values and raw_values[-1].upper() == "PARENT_SCOPE":
                    raw_values.pop()
                value, value_reasons = self._assignment_values(
                    raw_values, guard)
                assignment_reasons = self._merge_reasons(reasons, value_reasons)
                if conditional_depth:
                    assignment_reasons = self._merge_reasons(
                        assignment_reasons,
                        (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                self._assign(variable, value, guard, assignment_reasons)
                continue
            if command == "unset" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH mutation is unsupported")
                unset_reasons = reasons
                if conditional_depth:
                    unset_reasons = self._merge_reasons(
                        unset_reasons,
                        (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                self._assign(variable, None, guard, unset_reasons)
                continue
            if command == "option" and len(tokens) >= 2:
                variable = self._mutation_variable(tokens[0], guard, command)
                default = tokens[-1] if len(tokens) >= 3 else "OFF"
                option_reasons = reasons
                if conditional_depth:
                    option_reasons = self._merge_reasons(
                        option_reasons,
                        (CONDITIONAL_ASSIGNMENT_PREFIX + variable,))
                self._assign(variable, (default,), guard, option_reasons)
                continue
            if command in {"check_cxx_compiler_flag",
                           "check_c_compiler_flag"} and len(tokens) >= 2:
                variable = self._mutation_variable(
                    tokens[-1], guard, command)
                symbol = BOOL_SYMBOL_PREFIX + "probe:" + variable
                self._assign(variable, (symbol,), guard, reasons)
                continue
            if command == "cmake_dependent_option" and tokens:
                variable = self._mutation_variable(tokens[0], guard, command)
                symbol = BOOL_SYMBOL_PREFIX + "dependent-option:" + variable
                self._assign(variable, (symbol,), guard, reasons)
                continue
            if command == "list" and len(tokens) >= 2:
                operation = tokens[0].upper()
                variable = self._mutation_variable(tokens[1], guard, command)
                if variable == "CMAKE_MODULE_PATH":
                    raise ContractError(
                        "CMAKE_MODULE_PATH mutation is unsupported")
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
                else:
                    reason = (
                        "unsupported list operation touches source variable: " +
                        variable)
                    self._assign(variable, (), guard,
                                 self._merge_reasons(reasons, (reason,)))
                continue

            if (command == "target_sources" and tokens and
                    tokens[0] == "libleopard" and
                    not bool_satisfiable(guard)):
                raise ContractError(
                    "libleopard TARGET_OBJECTS has no MSVC-reachable "
                    "definition or attachment configuration")
            if (command == "add_library" and tokens and
                    tokens[0] == "libleopard" and conditional_depth):
                raise ContractError(
                    "libleopard add_library must be unconditional")
            if not bool_satisfiable(guard):
                if command == "add_library" and tokens:
                    self.declared_targets.add(tokens[0])
                continue
            variable_writers = {
                "file", "string", "math", "separate_arguments",
                "cmake_path", "get_property", "get_cmake_property",
                "get_directory_property", "get_filename_component",
                "get_source_file_property", "get_target_property",
                "execute_process", "try_compile", "try_run",
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

    def _read_graph_command(self, command, raw_tokens, guard, reasons,
                            conditional_depth):
        library_types = {"STATIC", "SHARED", "MODULE", "OBJECT", "INTERFACE"}
        source_commands = {
            "add_library", "target_sources", "set_source_files_properties"}
        if reasons and command in source_commands:
            raise ContractError(
                "unsupported conditional CMake source graph: " + reasons[0])
        if command == "add_library" and raw_tokens:
            target = self._target_name(raw_tokens[0], guard)
            self.declared_targets.add(target)
            if target == "libleopard" and not bool_tautology(guard):
                raise ContractError(
                    "libleopard add_library must be unconditional")
            tokens = self._expand(raw_tokens[1:], guard)
            if "ALIAS" in tokens or "IMPORTED" in tokens:
                return
            kind = "DEFAULT"
            if tokens and tokens[0].upper() in library_types:
                kind = tokens.pop(0).upper()
            if tokens and tokens[0].upper() == "EXCLUDE_FROM_ALL":
                tokens.pop(0)
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
            tokens = self._expand(raw_tokens[1:], guard)
            if any(token.upper() == "FILE_SET" for token in tokens):
                raise ContractError(
                    "CMake FILE_SET source attachment requires parser support: " +
                    target)
            sources = [token for token in tokens if token.upper() not in {
                "PRIVATE", "PUBLIC", "INTERFACE", "SYSTEM", "BEFORE"}]
            if (conditional_depth and target == "libleopard" and any(
                    not self._target_objects.match(token) for token in sources)):
                raise ContractError(
                    "conditional direct libleopard source attachment")
            self.attachments.setdefault(target, []).extend(
                (source, guard) for source in sources)
        elif command == "set_source_files_properties" and raw_tokens:
            self._record_source_properties(command, raw_tokens, guard, reasons)
        elif (command == "set_property" and raw_tokens and
              raw_tokens[0].upper() == "SOURCE"):
            self._record_source_properties(command, raw_tokens, guard, reasons)
        elif (command == "set_property" and raw_tokens and
              raw_tokens[0].upper() == "TARGET"):
            self._reject_target_source_property(command, raw_tokens, guard)
        elif command == "set_target_properties" and raw_tokens:
            self._reject_target_source_property(command, raw_tokens, guard)

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

    def _reject_target_source_property(self, command, raw_tokens, guard):
        tokens = self._expand(raw_tokens, guard)
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
                "CUDA source attached to ordinary CPU libleopard: " + token)
        return path.as_posix()

    def resolve(self, target="libleopard"):
        resolved = []
        attached_objects = set()
        object_sources = set()
        visiting = []

        def visit(name, reach_guard, reached_as_object=False):
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
            for token, entry_guard in entries:
                if not bool_satisfiable(entry_guard):
                    continue
                object_match = self._target_objects.match(token)
                if object_match:
                    object_target = object_match.group(1)
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
        resolved_by_path = {}
        for path, guard in resolved:
            resolved_by_path.setdefault(path, []).append(guard)
        for reference, property_guard, properties in self.source_property_references:
            if any(bool_satisfiable(bool_and(property_guard, source_guard))
                   for source_guard in resolved_by_path.get(reference, [])):
                raise ContractError(
                    "CMake source properties affect production " + reference +
                    ": " + " ".join(properties))
        duplicates = []
        for path, guards in resolved_by_path.items():
            if any(bool_satisfiable(bool_and(left, right))
                   for index, left in enumerate(guards)
                   for right in guards[index + 1:]):
                duplicates.append(path)
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


def validate_no_wpo_overrides(tree):
    allowed_additional_options = set(tree.findall(
        ".//msb:ClCompile[@Include='..\\Leopard2BackendAVX2.cpp']/"
        "msb:AdditionalOptions", NS))
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
        "toolexe", "toolpath",
    }
    for node in root.iter():
        name = xml_local_name(node).lower()
        if name in forbidden_properties:
            raise ContractError(
                "MSBuild compiler tool override is forbidden: " +
                xml_local_name(node))
        if name in {"target", "usingtask"}:
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


def production_graph(text=None, require_files=True):
    cmake = CMAKE.read_text(encoding="utf-8") if text is None else text
    configurations = (
        CMakeProductionGraph(
            cmake, processor="x86", pointer_size="4",
            platform_name="Win32").resolve(),
        CMakeProductionGraph(
            cmake, processor="i686", pointer_size="4",
            platform_name="Win32").resolve(),
        CMakeProductionGraph(
            cmake, processor="AMD64", pointer_size="8",
            platform_name="x64").resolve(),
    )
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
            for included in re.findall(
                    r'^\s*#\s*include\s+"([^"\n]+\.(?:h|hh|hpp|hxx|inl|cuh))"',
                    source.read_text(encoding="utf-8"), re.MULTILINE):
                candidate = (source.parent / included).resolve()
                try:
                    local = candidate.relative_to(ROOT).as_posix()
                except ValueError:
                    continue
                if candidate.suffix.lower() in CUDA_SUFFIXES:
                    raise ContractError(
                        "CUDA header reachable from CPU libleopard: " + local)
                if candidate.is_file() and local not in headers:
                    headers.add(local)
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
        if path == "Leopard2BackendAVX2.cpp":
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
    avx2 = options.get("Leopard2BackendAVX2.cpp", "")
    if avx2.strip() != "/arch:AVX2 %(AdditionalOptions)":
        raise ContractError("AVX2 backend options are not the exact contract")
    for path, flags in options.items():
        if path != "Leopard2BackendAVX2.cpp" and re.search(
                r"/arch\s*:\s*avx", flags, re.IGNORECASE):
            raise ContractError("non-AVX2 source raises ISA: " + path)


def validate_visual_studio_project(tree):
    validate_msbuild_element_casing(tree)
    validate_msbuild_imports_and_toolchain(tree)
    validate_source_item_structure(tree)
    validate_msbuild_configurations(tree)
    validate_per_file_isa(tree)
    validate_no_wpo_overrides(tree)


def validate_legacy_visual_studio_metadata():
    solution = SOLUTION.read_text(encoding="utf-8-sig")
    if ("# Visual Studio 14" not in solution or
            "VisualStudioVersion = 14.0.25420.1" not in solution):
        raise ContractError(
            "legacy solution does not declare the repository's VS2015 version")
    for project in LEGACY_PROJECTS:
        tree = ET.parse(str(project))
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
            "Leopard2BackendAVX2.cpp", "Leopard2BackendSSSE3.cpp"
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


class CMakeGraphMutationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cmake = CMAKE.read_text(encoding="utf-8")

    def resolve(self, mutation):
        return production_graph(
            self.cmake + "\n" + mutation + "\n", require_files=False)

    def resolve_text(self, text):
        return production_graph(text, require_files=False)

    def test_direct_target_sources_is_retained(self):
        (sources, unused_headers, unused_objects,
         unused_object_sources, unused_cmake) = self.resolve(
            "target_sources(libleopard PRIVATE New.cpp)")
        del unused_headers, unused_objects, unused_object_sources, unused_cmake
        self.assertIn("New.cpp", sources)

    def test_future_attached_object_target_is_traversed(self):
        mutation = """
add_library(future_backend OBJECT FutureBackend.cpp)
target_sources(libleopard PRIVATE $<TARGET_OBJECTS:future_backend>)
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
target_sources(libleopard PRIVATE ${OBJECT_EXPRESSION})
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
                "target_sources(libleopard PRIVATE ${EXTRA_SOURCES})")

    def test_resolved_variable_is_retained(self):
        mutation = """
set(EXTRA_SOURCES ExtraOne.cpp ExtraTwo.cpp)
target_sources(libleopard PRIVATE ${EXTRA_SOURCES})
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
target_sources(libleopard PRIVATE ${BRANCH_SOURCES})
"""
        with self.assertRaisesRegex(
                ContractError, "conditional.*BRANCH_SOURCES"):
            self.resolve(mutation)

    def test_unsupported_list_operation_on_source_variable_is_rejected(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
        mutation = "list(REMOVE_ITEM LIB_SOURCE_FILES leopard2.cpp)"
        text = self.cmake.replace(marker, mutation + "\n" + marker, 1)
        with self.assertRaisesRegex(
                ContractError, "unsupported list operation.*LIB_SOURCE_FILES"):
            self.resolve_text(text)

    def test_source_variable_is_snapshotted_at_add_library_time(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
        replacement = """set(SAVED_LIB_SOURCE_FILES ${LIB_SOURCE_FILES})
set(LIB_SOURCE_FILES Injected.cpp)
add_library(libleopard STATIC ${LIB_SOURCE_FILES})
set(LIB_SOURCE_FILES ${SAVED_LIB_SOURCE_FILES})"""
        text = self.cmake.replace(marker, replacement, 1)
        self.assertNotEqual(text, self.cmake)
        sources = self.resolve_text(text)[0]
        self.assertIn("Injected.cpp", sources)
        self.assertNotIn("leopard2.cpp", sources)

    def test_indirect_mutation_destinations_are_resolved_at_command_time(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
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
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
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
                        ContractError, "may mutate source variable"):
                    self.resolve_text(text)

    def test_unproved_graph_imports_are_rejected(self):
        for mutation in (
                "add_subdirectory(cmake/injected_sources)",
                "subdirs(cmake/injected_sources)",
                "include(cmake/injected_sources.cmake)",
                "cmake_language(CALL target_sources libleopard PRIVATE "
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
                        ContractError, "mutation in unsupported CMake block"):
                    self.resolve(mutation)

    def test_win32_only_source_variable_mutation_is_rejected(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
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

    def test_conditional_libleopard_definition_is_rejected(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
        replacement = """if(NOT WIN32)
add_library(libleopard STATIC ${LIB_SOURCE_FILES})
endif()"""
        text = self.cmake.replace(marker, replacement, 1)
        with self.assertRaisesRegex(
                ContractError, "add_library must be unconditional"):
            self.resolve_text(text)

    def test_contextual_libleopard_definition_is_rejected(self):
        marker = "add_library(libleopard STATIC ${LIB_SOURCE_FILES})"
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
                        ContractError, "add_library must be unconditional"):
                    self.resolve_text(text)

    def test_contextual_target_property_source_mutation_is_rejected(self):
        mutations = (
            "set_property(TARGET libleopard APPEND PROPERTY "
            "SOURCES Injected.cpp)",
            "set_target_properties(libleopard PROPERTIES "
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
add_library(libleopard STATIC ${LIB_SOURCE_FILES})
if(NOT WIN32)
    add_library(hidden_backend OBJECT HiddenBackend.cpp)
endif()
target_sources(libleopard PRIVATE $<TARGET_OBJECTS:hidden_backend>)
""", """
set(LIB_SOURCE_FILES Base.cpp Base.h)
add_library(libleopard STATIC ${LIB_SOURCE_FILES})
if(MSVC)
    if(FALSE)
        add_library(hidden_backend OBJECT HiddenBackend.cpp)
    endif()
endif()
target_sources(libleopard PRIVATE $<TARGET_OBJECTS:hidden_backend>)
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
                        ContractError, "no MSVC-reachable"):
                    self.resolve_text(text)

    def test_target_objects_attachment_must_reach_msvc(self):
        marker = """    target_sources(libleopard PRIVATE
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
        marker = """    target_sources(libleopard PRIVATE
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
        marker = """    target_sources(libleopard PRIVATE
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
    target_sources(libleopard PRIVATE Conditional.cpp)
endif()
"""
        with self.assertRaisesRegex(
                ContractError, "conditional direct libleopard"):
            self.resolve(mutation)

    def test_source_properties_on_production_source_are_rejected(self):
        mutations = (
            "set_source_files_properties(leopard2.cpp "
            "PROPERTIES HEADER_FILE_ONLY TRUE)",
            "set_property(SOURCE leopard2.cpp "
            "PROPERTY HEADER_FILE_ONLY TRUE)",
        )
        for mutation in mutations:
            with self.subTest(command=mutation.split("(", 1)[0]):
                with self.assertRaisesRegex(
                        ContractError, "source properties affect production"):
                    self.resolve(mutation)

    def test_target_sources_property_mutation_is_rejected(self):
        mutations = (
            "set_property(TARGET libleopard APPEND PROPERTY "
            "SOURCES Injected.cpp)",
            "set_target_properties(libleopard PROPERTIES "
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
                "target_sources(libleopard PRIVATE "
                "\"$<$<BOOL:1>:Conditional.cpp>\")")

    def test_cuda_source_or_header_is_rejected(self):
        for path in ("accidental.cu", "accidental.cuh"):
            with self.subTest(path=path):
                with self.assertRaisesRegex(ContractError, "CUDA source"):
                    self.resolve(
                        "target_sources(libleopard PRIVATE " + path + ")")

    def test_path_escape_is_rejected(self):
        for path in ("../Escape.cpp", "sub/../Escape.cpp", "C:/Escape.cpp"):
            with self.subTest(path=path):
                with self.assertRaisesRegex(
                        ContractError, "repository-relative|drive path"):
                    self.resolve(
                        "target_sources(libleopard PRIVATE " + path + ")")


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
        with self.assertRaisesRegex(ContractError, "execution logic"):
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
                        "ProjectConfigurations|ProjectConfiguration attributes"):
                    validate_visual_studio_project(tree)

    def test_project_configuration_transform_outside_group_is_rejected(self):
        namespace = "{http://schemas.microsoft.com/developer/msbuild/2003}"
        tree = self.fresh_tree()
        group = ET.SubElement(tree.getroot(), namespace + "ItemGroup")
        transformed = ET.SubElement(
            group, namespace + "ProjectConfiguration")
        transformed.set("Remove", "Debug|Win32")
        with self.assertRaisesRegex(
                ContractError, "outside the canonical group"):
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
        with self.assertRaisesRegex(ContractError, "direct Project children"):
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
