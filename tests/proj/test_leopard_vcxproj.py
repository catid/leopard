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


class ContractError(AssertionError):
    pass


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
        self.variables = {}
        self.targets = {}
        self.attachments = {}
        self._read_variables()
        self._read_targets()

    def _read_variables(self):
        for name, tokens in self.commands:
            if name == "set" and tokens:
                variable = tokens[0]
                values = list(tokens[1:])
                if "CACHE" in values:
                    values = values[:values.index("CACHE")]
                if values and values[-1] == "PARENT_SCOPE":
                    values.pop()
                self.variables[variable] = values
            elif name == "unset" and tokens:
                self.variables.pop(tokens[0], None)
            elif name == "list" and len(tokens) >= 2:
                operation = tokens[0].upper()
                variable = tokens[1]
                if operation == "APPEND":
                    self.variables.setdefault(variable, []).extend(tokens[2:])
                elif operation == "PREPEND":
                    self.variables[variable] = (
                        list(tokens[2:]) + self.variables.get(variable, []))

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
        for command, raw_tokens in self.commands:
            if command == "add_library" and raw_tokens:
                target = self._target_name(raw_tokens[0])
                tokens = self._expand(raw_tokens[1:])
                if "ALIAS" in tokens or "IMPORTED" in tokens:
                    continue
                kind = "DEFAULT"
                if tokens and tokens[0].upper() in library_types:
                    kind = tokens.pop(0).upper()
                if tokens and tokens[0].upper() == "EXCLUDE_FROM_ALL":
                    tokens.pop(0)
                definition = (kind, tokens)
                if target in self.targets:
                    # The GNU/Clang and MSVC branches define the same isolated
                    # object target with different flags.  Coalesce only an
                    # identical source definition; any source drift is
                    # ambiguous and therefore rejected.
                    if self.targets[target] != definition:
                        raise ContractError("conflicting CMake target: " + target)
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
                self.attachments.setdefault(target, []).extend(sources)

    @staticmethod
    def _literal_source(token):
        if "$<" in token:
            raise ContractError(
                "unsupported CMake source generator expression: " + token)
        if "$" in token:
            raise ContractError("unresolved CMake source token: " + token)
        path = PurePosixPath(token.replace("\\", "/"))
        if path.is_absolute() or ".." in path.parts:
            raise ContractError(
                "production source must be repository-relative: " + token)
        suffix = path.suffix.lower()
        if suffix not in KNOWN_SOURCE_SUFFIXES:
            raise ContractError(
                "unsupported production source path (fail closed): " + token)
        if suffix in CUDA_SUFFIXES:
            raise ContractError(
                "CUDA source attached to ordinary CPU libleopard: " + token)
        return path.as_posix()

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
            kind, sources = self.targets[name]
            for token in list(sources) + self.attachments.get(name, []):
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
                    literal = self._literal_source(token)
                    resolved.append(literal)
                    if reached_as_object:
                        object_sources.append(literal)
            visiting.pop()

        visit(target)
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
    key = match.group(1).strip().lower() + "|" + match.group(2).strip().lower()
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
        macro_set = set(raw_definitions.split(";"))
        required = set(BACKEND_DEFINITIONS)
        required.update(("_MBCS", "%(PreprocessorDefinitions)"))
        if not required.issubset(macro_set):
            raise ContractError(key + " is missing backend isolation definitions")
        if any(macro.startswith("LEO2_BACKEND_FORCE_") for macro in macro_set):
            raise ContractError(key + " forces a diagnostic backend")
        options = " ".join(
            node.text or "" for node in compile_node.findall(
                "msb:AdditionalOptions", NS))
        if re.search(r"/arch\s*:\s*avx", options, re.IGNORECASE):
            raise ContractError(key + " raises the project-wide ISA floor")
        enhanced = " ".join(
            node.text or "" for node in compile_node.findall(
                "msb:EnableEnhancedInstructionSet", NS))
        if enhanced.strip() not in ("", "NotSet"):
            raise ContractError(key + " raises the project-wide ISA floor")

    for node in tree.findall(".//msb:WholeProgramOptimization", NS):
        if (node.text or "").strip().lower() != "false":
            raise ContractError("WholeProgramOptimization must always be false")


def validate_per_file_isa(tree):
    compile_nodes = tree.findall(".//msb:ClCompile[@Include]", NS)
    options = {}
    for node in compile_nodes:
        path = project_path(PROJECT, node.attrib["Include"])
        options[path] = " ".join(
            child.text or "" for child in node.findall(
                "msb:AdditionalOptions", NS))
    avx2 = options.get("Leopard2BackendAVX2.cpp", "")
    if "/arch:AVX2" not in avx2 or "%(AdditionalOptions)" not in avx2:
        raise ContractError("AVX2 backend lacks its inherited /arch:AVX2 option")
    for path, flags in options.items():
        if path != "Leopard2BackendAVX2.cpp" and re.search(
                r"/arch\s*:\s*avx", flags, re.IGNORECASE):
            raise ContractError("non-AVX2 source raises ISA: " + path)


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
        validate_msbuild_configurations(self.project)
        validate_per_file_isa(self.project)

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
                with self.assertRaisesRegex(ContractError, "duplicate"):
                    validate_msbuild_configurations(tree)

    def test_missing_condition_is_rejected(self):
        for xpath, label in (
                (".//msb:PropertyGroup[@Label='Configuration']", "property"),
                (".//msb:ItemDefinitionGroup", "definition")):
            with self.subTest(group=label):
                tree = self.fresh_tree()
                node = tree.findall(xpath, NS)[0]
                tree.getroot().remove(node)
                with self.assertRaisesRegex(ContractError, "Conditions differ"):
                    validate_msbuild_configurations(tree)


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]))
    sys.exit(0 if result.wasSuccessful() else 1)
