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
PROJECT = ROOT / "proj" / "Leopard.vcxproj"
FILTERS = ROOT / "proj" / "Leopard.vcxproj.filters"
CMAKE = ROOT / "CMakeLists.txt"
NS = {"msb": "http://schemas.microsoft.com/developer/msbuild/2003"}

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


class CMakeProductionGraph(object):
    """Small fail-closed evaluator for source attachment commands."""

    _variable = re.compile(r"^\$\{([A-Za-z_][A-Za-z0-9_]*)\}$")
    _target_objects = re.compile(
        r"^\$<TARGET_OBJECTS:([A-Za-z_][A-Za-z0-9_.+\-]*)>$")

    def __init__(self, text):
        self.commands = [(name, cmake_tokens(body))
                         for name, body in cmake_commands(text)]
        self.annotated_commands = self._annotate_commands(self.commands)
        self.variables = {}
        self.conditional_variables = set()
        self.unsupported_variables = set()
        self.targets = {}
        self.target_definition_contexts = {}
        self.attachments = {}
        self.source_property_references = []
        self._read_variables()
        self._read_targets()

    @staticmethod
    def _annotate_commands(commands):
        annotated = []
        stack = []
        starts = {
            "if": "endif",
            "function": "endfunction",
            "macro": "endmacro",
            "foreach": "endforeach",
            "while": "endwhile",
            "block": "endblock",
        }
        ends = {value: key for key, value in starts.items()}
        for name, tokens in commands:
            if name in starts:
                stack.append((name, tuple(tokens)))
                continue
            if name in ("elseif", "else"):
                if not stack or stack[-1][0] != "if":
                    raise ContractError("unbalanced CMake " + name)
                branch = tuple(tokens) if name == "elseif" else ("ELSE",)
                stack[-1] = ("if", branch)
                continue
            if name in ends:
                if not stack or stack[-1][0] != ends[name]:
                    raise ContractError("unbalanced CMake " + name)
                stack.pop()
                continue
            annotated.append((name, tokens, tuple(stack)))
        if stack:
            raise ContractError("unbalanced CMake block: " + stack[-1][0])
        return annotated

    def _read_variables(self):
        for name, tokens, context in self.annotated_commands:
            conditional = bool(context)
            if name == "set" and tokens:
                variable = tokens[0]
                if conditional:
                    self.conditional_variables.add(variable)
                values = list(tokens[1:])
                if "CACHE" in values:
                    values = values[:values.index("CACHE")]
                if values and values[-1] == "PARENT_SCOPE":
                    values.pop()
                self.variables[variable] = values
            elif name == "unset" and tokens:
                if conditional:
                    self.conditional_variables.add(tokens[0])
                self.variables.pop(tokens[0], None)
            elif name == "list" and len(tokens) >= 2:
                operation = tokens[0].upper()
                variable = tokens[1]
                if conditional:
                    self.conditional_variables.add(variable)
                if operation == "APPEND":
                    self.variables.setdefault(variable, []).extend(tokens[2:])
                elif operation == "PREPEND":
                    self.variables[variable] = (
                        list(tokens[2:]) + self.variables.get(variable, []))
                else:
                    for token in tokens[1:]:
                        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", token):
                            self.unsupported_variables.add(token)

    def _expand(self, tokens, stack=None):
        stack = [] if stack is None else stack
        expanded = []
        for token in tokens:
            match = self._variable.match(token)
            if match:
                variable = match.group(1)
                if variable in stack:
                    raise ContractError(
                        "recursive CMake source variable: " + variable)
                if variable not in self.variables:
                    raise ContractError(
                        "unresolved CMake source variable: " + variable)
                if variable in self.conditional_variables:
                    raise ContractError(
                        "conditional CMake source variable is ambiguous: " +
                        variable)
                if variable in self.unsupported_variables:
                    raise ContractError(
                        "unsupported list operation touches source variable: " +
                        variable)
                expanded.extend(self._expand(
                    self.variables[variable], stack + [variable]))
            elif "${" in token or "$ENV{" in token:
                raise ContractError(
                    "embedded or environment source variable is unsupported: " +
                    token)
            else:
                expanded.append(token)
        return expanded

    def _target_name(self, token):
        names = self._expand([token])
        if len(names) != 1 or not re.match(
                r"^[A-Za-z_][A-Za-z0-9_.+\-]*$", names[0]):
            raise ContractError("unresolved CMake target name: " + token)
        return names[0]

    def _read_targets(self):
        library_types = {"STATIC", "SHARED", "MODULE", "OBJECT", "INTERFACE"}
        for command, raw_tokens, context in self.annotated_commands:
            conditional = bool(context)
            if command == "add_library" and raw_tokens:
                target = self._target_name(raw_tokens[0])
                if target == "libleopard" and conditional:
                    raise ContractError(
                        "libleopard add_library must be unconditional")
                tokens = self._expand(raw_tokens[1:])
                if "ALIAS" in tokens or "IMPORTED" in tokens:
                    continue
                kind = "DEFAULT"
                if tokens and tokens[0].upper() in library_types:
                    kind = tokens.pop(0).upper()
                if tokens and tokens[0].upper() == "EXCLUDE_FROM_ALL":
                    tokens.pop(0)
                definition = (kind, tokens, conditional)
                self.target_definition_contexts.setdefault(target, []).append(
                    context)
                if target in self.targets:
                    # The GNU/Clang and MSVC branches define the same isolated
                    # object target with different flags.  Coalesce only an
                    # identical source definition; any source drift is
                    # ambiguous and therefore rejected.
                    if self.targets[target][:2] != definition[:2]:
                        raise ContractError("conflicting CMake target: " + target)
                    previous = self.targets[target]
                    self.targets[target] = (
                        previous[0], previous[1],
                        previous[2] or conditional)
                else:
                    self.targets[target] = definition
            elif command == "target_sources" and raw_tokens:
                target = self._target_name(raw_tokens[0])
                tokens = self._expand(raw_tokens[1:])
                if any(token.upper() == "FILE_SET" for token in tokens):
                    raise ContractError(
                        "CMake FILE_SET source attachment requires parser support: " +
                        target)
                sources = [token for token in tokens if token.upper() not in {
                    "PRIVATE", "PUBLIC", "INTERFACE", "SYSTEM", "BEFORE"}]
                if conditional and target == "libleopard" and any(
                        not self._target_objects.match(token)
                        for token in sources):
                    raise ContractError(
                        "conditional direct libleopard source attachment")
                self.attachments.setdefault(target, []).extend(
                    (source, conditional) for source in sources)
            elif command == "set_source_files_properties" and raw_tokens:
                self._record_source_properties(
                    command, raw_tokens, conditional)
            elif (command == "set_property" and raw_tokens and
                  raw_tokens[0].upper() == "SOURCE"):
                self._record_source_properties(
                    command, raw_tokens, conditional)
            elif (command == "set_property" and raw_tokens and
                  raw_tokens[0].upper() == "TARGET"):
                self._reject_target_source_property(command, raw_tokens)
            elif command == "set_target_properties" and raw_tokens:
                self._reject_target_source_property(command, raw_tokens)

    def _record_source_properties(self, command, raw_tokens, conditional):
        tokens = self._expand(raw_tokens)
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
                (reference, conditional, list(properties)))

    def _reject_target_source_property(self, command, raw_tokens):
        tokens = self._expand(raw_tokens)
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

    @staticmethod
    def _context_is_msvc_reachable(context):
        if not context:
            return True
        positive_windows_or_msvc = False
        for block, tokens in context:
            if block != "if":
                return False
            expression = " ".join(tokens).upper()
            if (re.search(
                    r"(?:^|[\s(])(?:0|FALSE|OFF|NO|N|IGNORE|NOTFOUND)"
                    r"(?:$|[\s)])", expression) and
                    not re.search(r"\bOR\b", expression)):
                return False
            if re.search(r"\bNOT\s+(?:WIN32|MSVC)\b", expression):
                return False
            if ("CMAKE_CXX_COMPILER_ID" in expression and
                    re.search(r"\b(?:MATCHES|STREQUAL)\b", expression) and
                    "MSVC" not in expression):
                return False
            if re.search(r"\b(?:WIN32|MSVC)\b", expression):
                positive_windows_or_msvc = True
            if ("CMAKE_CXX_COMPILER_ID" in expression and
                    "MSVC" in expression):
                positive_windows_or_msvc = True
        return positive_windows_or_msvc

    def resolve(self, target="libleopard"):
        resolved = []
        attached_objects = []
        object_sources = []
        visiting = []

        def visit(name, reached_as_object=False):
            if name in visiting:
                raise ContractError(
                    "cyclic CMake TARGET_OBJECTS attachment: " + name)
            if name not in self.targets:
                raise ContractError("attached CMake target is undefined: " + name)
            visiting.append(name)
            kind, sources, definition_conditional = self.targets[name]
            entries = [
                (token, definition_conditional, True) for token in sources
            ] + [
                (token, conditional, False)
                for token, conditional in self.attachments.get(name, [])
            ]
            for token, conditional, from_definition in entries:
                object_match = self._target_objects.match(token)
                if object_match:
                    object_target = object_match.group(1)
                    if object_target not in self.targets:
                        raise ContractError(
                            "TARGET_OBJECTS target is undefined: " + object_target)
                    if self.targets[object_target][0] != "OBJECT":
                        raise ContractError(
                            "TARGET_OBJECTS does not name an OBJECT library: " +
                            object_target)
                    attached_objects.append(object_target)
                    visit(object_target, True)
                else:
                    if conditional and not from_definition:
                        raise ContractError(
                            "conditional CMake source attachment: " + token)
                    literal = self._literal_source(token)
                    resolved.append(literal)
                    if reached_as_object:
                        object_sources.append(literal)
            visiting.pop()

        visit(target)
        for object_target in sorted(set(attached_objects)):
            contexts = self.target_definition_contexts.get(object_target, [])
            if not any(self._context_is_msvc_reachable(context)
                       for context in contexts):
                raise ContractError(
                    "attached OBJECT target has no MSVC-reachable definition: " +
                    object_target)
        resolved_set = set(resolved)
        for reference, conditional, properties in self.source_property_references:
            if reference in resolved_set:
                qualifier = "conditional " if conditional else ""
                raise ContractError(
                    qualifier + "CMake source properties affect production " +
                    reference + ": " + " ".join(properties))
        duplicates = sorted(
            path for path, count in Counter(resolved).items() if count != 1)
        if duplicates:
            raise ContractError(
                "duplicate production source attachment: " + ", ".join(duplicates))
        return (sorted(resolved), sorted(set(attached_objects)),
                sorted(set(object_sources)))


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


def production_graph(text=None, require_files=True):
    cmake = CMAKE.read_text(encoding="utf-8") if text is None else text
    graph = CMakeProductionGraph(cmake)
    attached, object_targets, object_sources = graph.resolve()
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


def index_conditions(nodes, label):
    indexed = {}
    for node in nodes:
        key = normalized_condition(node.attrib.get("Condition"))
        if key in indexed:
            raise ContractError("duplicate " + label + " Condition: " + key)
        indexed[key] = node
    if set(indexed) != set(EXPECTED_CONFIGS):
        raise ContractError(
            label + " Conditions differ: " + ", ".join(sorted(indexed)))
    return indexed


def validate_msbuild_configurations(tree):
    configurations = tree.findall(
        ".//msb:ItemGroup[@Label='ProjectConfigurations']/"
        "msb:ProjectConfiguration", NS)
    includes = [node.attrib.get("Include", "").lower()
                for node in configurations]
    if len(includes) != len(EXPECTED_CONFIGS) or set(includes) != set(
            EXPECTED_CONFIGS):
        raise ContractError("ProjectConfigurations are not the exact four configs")
    for node in configurations:
        key = node.attrib["Include"].lower()
        expected_configuration, expected_platform = EXPECTED_CONFIGS[key]
        actual_configuration = node.findtext(
            "msb:Configuration", namespaces=NS)
        if actual_configuration != expected_configuration:
            raise ContractError(key + " has a mismatched Configuration value")
        if node.findtext("msb:Platform", namespaces=NS) != expected_platform:
            raise ContractError(key + " has a mismatched Platform value")

    properties = index_conditions(tree.findall(
        ".//msb:PropertyGroup[@Label='Configuration']", NS),
        "Configuration PropertyGroup")
    definitions = index_conditions(tree.findall(
        ".//msb:ItemDefinitionGroup", NS), "ItemDefinitionGroup")

    for key in sorted(EXPECTED_CONFIGS):
        configuration, unused_platform = EXPECTED_CONFIGS[key]
        del unused_platform
        prop = properties[key]

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
        if unique_text(prop, "PlatformToolset", "PlatformToolset") != "v140":
            raise ContractError(key + " does not preserve the v140 toolset")
        wpo = unique_text(
            prop, "WholeProgramOptimization", "WholeProgramOptimization")
        if (wpo or "").strip().lower() != "false":
            raise ContractError(key + " does not explicitly disable /GL")

        compile_nodes = definitions[key].findall("msb:ClCompile", NS)
        if len(compile_nodes) != 1:
            raise ContractError(key + " must define exactly one ClCompile")
        compile_node = compile_nodes[0]
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
    validate_source_item_structure(tree)
    validate_msbuild_configurations(tree)
    validate_per_file_isa(tree)
    validate_no_wpo_overrides(tree)


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
