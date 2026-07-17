#!/usr/bin/env python3
"""Fail closed when current documentation uses pre-rename library names."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]

# These documents intentionally describe either the deprecated compatibility
# target or byte-exact evidence produced before the canonical-target rename.
# Every occurrence must also carry an explicit historical/compatibility label.
LEGACY_REFERENCE_DOCS = {
    "README_CUDA.md",
    "docs/leopard2_backend_butterfly_tier.md",
    "docs/leopard2_backend_isolation.md",
    "docs/leopard2_c6_gf256_rescue.md",
    "docs/leopard2_c7_exact_low.md",
    "docs/reproduction/leopard2_low_encode_copy.md",
    "experiments/leopard2/low_encode_copy/README.md",
    "experiments/leopard2/main_compare/README.md",
}

LEGACY_TARGET_NAME = "lib" + "leopard"
LEGACY_ARCHIVE_NAME = "lib" + LEGACY_TARGET_NAME + ".a"
LEGACY_NAMESPACE_NAME = "leopard::" + LEGACY_TARGET_NAME
LEGACY_REFERENCE = re.compile(
    re.escape(LEGACY_ARCHIVE_NAME) + "|" +
    re.escape(LEGACY_NAMESPACE_NAME) + "|" +
    r"(?<![:A-Za-z0-9_])" + re.escape(LEGACY_TARGET_NAME) +
    r"(?![.A-Za-z0-9_:])")
FUNCTIONAL_COMMAND = re.compile(
    r"--target(?:=|[ \t]+)[^\n]*\b" +
    re.escape(LEGACY_TARGET_NAME) + r"\b|" +
    r"(?:^|[/= \t'\"])(?:[^\n'\"]*/)?" +
    re.escape(LEGACY_ARCHIVE_NAME),
    re.MULTILINE)
CLASSIFICATION_MARKER = re.compile(
    r"historical|pre-rename|deprecated|compatib|retained|frozen|"
    r"original|not\s+named",
    re.IGNORECASE)


def documentation_files():
    paths = set(ROOT.glob("*.md"))
    for relative_root in ("docs", "experiments/leopard2"):
        directory = ROOT / relative_root
        for suffix in ("*.md", "*.txt"):
            paths.update(directory.rglob(suffix))
    return sorted(
        path for path in paths
        if "results" not in path.parts and "artifacts" not in path.parts)


def legacy_references(text):
    return list(LEGACY_REFERENCE.finditer(text))


class CanonicalLibraryDocumentationTest(unittest.TestCase):

    def test_current_commands_use_canonical_target_and_archive(self):
        stale = []
        for path in documentation_files():
            text = path.read_text(encoding="utf-8")
            for match in FUNCTIONAL_COMMAND.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                stale.append(
                    f"{path.relative_to(ROOT).as_posix()}:{line}: "
                    f"{match.group(0).strip()}")
        self.assertEqual([], stale)

    def test_legacy_prose_is_explicitly_classified(self):
        actual = set()
        unclassified = []
        for path in documentation_files():
            text = path.read_text(encoding="utf-8")
            matches = legacy_references(text)
            if not matches:
                continue
            relative = path.relative_to(ROOT).as_posix()
            actual.add(relative)
            for match in matches:
                context = text[max(0, match.start() - 240):match.end() + 240]
                if not CLASSIFICATION_MARKER.search(context):
                    line = text.count("\n", 0, match.start()) + 1
                    unclassified.append(f"{relative}:{line}")

        self.assertEqual([], unclassified)
        self.assertEqual(
            LEGACY_REFERENCE_DOCS,
            actual,
            "historical/compatibility documentation allowlist drifted")

    def test_scanner_rejects_stale_reproduction_examples(self):
        stale_examples = (
            "cmake --build build/release --target " + LEGACY_TARGET_NAME,
            "c++ probe.cpp build/release/" + LEGACY_ARCHIVE_NAME,
            "--library=build/release/" + LEGACY_ARCHIVE_NAME,
        )
        for example in stale_examples:
            with self.subTest(example=example):
                self.assertIsNotNone(FUNCTIONAL_COMMAND.search(example))

        historical = (
            "Retained historical v1 evidence binds `" +
            LEGACY_ARCHIVE_NAME + "`."
        )
        match = legacy_references(historical)[0]
        context = historical[max(0, match.start() - 240):match.end() + 240]
        self.assertIsNotNone(CLASSIFICATION_MARKER.search(context))


if __name__ == "__main__":
    unittest.main()
