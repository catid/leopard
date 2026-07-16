#!/usr/bin/env python3
"""Structural contract for the checked-in Visual Studio library project.

This test intentionally runs on any host with Python.  It cannot replace a
native MSVC build, but it makes the hand-maintained project fail closed when
the production CMake source graph or the ISA/CUDA isolation contract changes.
"""

from collections import Counter
from pathlib import Path
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


def cmake_body(text, expression):
    match = re.search(expression, text, re.MULTILINE | re.DOTALL)
    if not match:
        raise AssertionError("CMake expression not found: " + expression)
    return match.group(1)


def cmake_paths(body):
    body = re.sub(r"#[^\n]*", "", body)
    return [
        token.strip('"')
        for token in re.split(r"\s+", body.strip())
        if token.strip('"').endswith((".c", ".cc", ".cpp", ".cxx", ".h"))
    ]


def project_path(project_file, value):
    native = value.replace("\\", "/")
    return (project_file.parent / native).resolve().relative_to(ROOT).as_posix()


def item_paths(tree, kind, project_file):
    return [
        project_path(project_file, node.attrib["Include"])
        for node in tree.findall(".//msb:" + kind + "[@Include]", NS)
    ]


def production_graph():
    text = CMAKE.read_text(encoding="utf-8")
    library = cmake_paths(cmake_body(
        text, r"set\s*\(\s*LIB_SOURCE_FILES\b(.*?)\)"))
    compiled = [path for path in library if path.endswith(COMPILE_SUFFIXES)]
    for target in ("leopard2_backend_ssse3", "leopard2_backend_avx2"):
        body = cmake_body(
            text,
            r"add_library\s*\(\s*" + re.escape(target) +
            r"\s+OBJECT\b(.*?)\)")
        compiled.extend(
            path for path in cmake_paths(body)
            if path.endswith(COMPILE_SUFFIXES))

    # Display every repository-local production header reachable from the
    # authoritative CMake graph.  This includes implementation-only headers
    # such as Leopard2Dispatch.h even when CMake need not list them explicitly.
    headers = {path for path in library if path.endswith(".h")}
    pending = list(compiled) + list(headers)
    visited = set()
    while pending:
        relative = pending.pop()
        if relative in visited:
            continue
        visited.add(relative)
        source = ROOT / relative
        if not source.is_file():
            raise AssertionError("production graph names missing file: " + relative)
        if source.suffix not in COMPILE_SUFFIXES + (".h",):
            continue
        for included in re.findall(
                r'^\s*#\s*include\s+"([^"/\\]+\.h)"',
                source.read_text(encoding="utf-8"), re.MULTILINE):
            if (ROOT / included).is_file() and included not in headers:
                headers.add(included)
                pending.append(included)
    return compiled, sorted(headers), text


class LeopardVisualStudioProjectTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = ET.parse(str(PROJECT))
        cls.filters = ET.parse(str(FILTERS))
        cls.expected_sources, cls.expected_headers, cls.cmake = production_graph()
        cls.sources = item_paths(cls.project, "ClCompile", PROJECT)
        cls.headers = item_paths(cls.project, "ClInclude", PROJECT)

    def assert_unique_equal(self, actual, expected, label):
        duplicates = sorted(path for path, count in Counter(actual).items() if count != 1)
        self.assertEqual([], duplicates, label + " duplicated")
        self.assertEqual(sorted(expected), sorted(actual), label + " drifted")

    def test_compiled_sources_match_production_cmake_graph(self):
        self.assert_unique_equal(
            self.sources, self.expected_sources, "compiled production sources")

    def test_reachable_production_headers_are_visible(self):
        self.assert_unique_equal(
            self.headers, self.expected_headers, "production headers")

    def test_filters_match_project_exactly(self):
        filter_sources = item_paths(self.filters, "ClCompile", FILTERS)
        filter_headers = item_paths(self.filters, "ClInclude", FILTERS)
        self.assert_unique_equal(filter_sources, self.sources, "source filters")
        self.assert_unique_equal(filter_headers, self.headers, "header filters")
        for node in self.filters.findall(".//msb:ClCompile", NS):
            self.assertEqual("Source Files", node.findtext("msb:Filter", namespaces=NS))
        for node in self.filters.findall(".//msb:ClInclude", NS):
            self.assertEqual("Header Files", node.findtext("msb:Filter", namespaces=NS))

    def test_baseline_and_backend_isa_are_isolated(self):
        required = {
            "LEO2_DISABLE_SSSE3_CODEGEN=1",
            "LEO2_DISABLE_AVX2_CODEGEN=1",
            "LEO2_HAVE_SSSE3_BACKEND=1",
            "LEO2_HAVE_AVX2_BACKEND=1",
            "%(PreprocessorDefinitions)",
        }
        definitions = self.project.findall(
            ".//msb:ItemDefinitionGroup/msb:ClCompile/msb:PreprocessorDefinitions", NS)
        self.assertEqual(4, len(definitions))
        for node in definitions:
            self.assertTrue(required.issubset(set((node.text or "").split(";"))))

        compile_nodes = self.project.findall(".//msb:ClCompile[@Include]", NS)
        options = {}
        for node in compile_nodes:
            path = project_path(PROJECT, node.attrib["Include"])
            options[path] = " ".join(
                child.text or "" for child in node.findall("msb:AdditionalOptions", NS))
        self.assertIn("/arch:AVX2", options["Leopard2BackendAVX2.cpp"])
        self.assertIn("%(AdditionalOptions)", options["Leopard2BackendAVX2.cpp"])
        for path, flags in options.items():
            if path != "Leopard2BackendAVX2.cpp":
                self.assertNotRegex(flags.lower(), r"/arch\s*:\s*avx")

        project_options = " ".join(
            node.text or "" for node in self.project.findall(
                ".//msb:ItemDefinitionGroup/msb:ClCompile/"
                "msb:AdditionalOptions", NS))
        self.assertNotRegex(project_options.lower(), r"/arch\s*:\s*avx")

        for node in self.project.findall(".//msb:EnableEnhancedInstructionSet", NS):
            self.assertIn((node.text or "").strip(), ("", "NotSet"))
        for node in self.project.findall(".//msb:WholeProgramOptimization", NS):
            self.assertEqual("false", (node.text or "").strip().lower())

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
        gate = cmake_body(
            self.cmake, r"if\s*\(\s*LEO2_ENABLE_CUDA\s*\)(.*?)endif\s*\(\s*\)")
        self.assertRegex(gate, r"enable_language\s*\(\s*CUDA\s*\)")
        self.assertNotIn(".cu", " ".join(self.expected_sources).lower())


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        LeopardVisualStudioProjectTest)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
