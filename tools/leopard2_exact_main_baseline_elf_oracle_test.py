#!/usr/bin/env python3
"""Linux-only differential oracle for the exact-main ELF64 parser."""

from __future__ import annotations

import hashlib
import importlib.util
import pathlib
import platform
import re
import subprocess
import sys
import tempfile
import unittest


MODULE_PATH = pathlib.Path(__file__).with_name(
    "leopard2_exact_main_baseline.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_exact_main_baseline_oracle_subject", MODULE_PATH,
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load exact-main baseline contract")
contract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contract)

FIXTURE_PATH = pathlib.Path(__file__).with_name(
    "leopard2_exact_main_baseline_test.py")
FIXTURE_SPEC = importlib.util.spec_from_file_location(
    "leopard2_exact_main_baseline_oracle_fixtures", FIXTURE_PATH,
)
if FIXTURE_SPEC is None or FIXTURE_SPEC.loader is None:
    raise RuntimeError("cannot load exact-main ELF fixtures")
fixtures = importlib.util.module_from_spec(FIXTURE_SPEC)
FIXTURE_SPEC.loader.exec_module(fixtures)


READELF = pathlib.Path("/usr/bin/readelf")
ROW_PATTERN = re.compile(r"^\s*\[\s*(\d+)\]\s+(.*)$")
ORACLE_ROOTS = {
    "adapter_source_root": "/lib",
    "baseline_source_root": "/usr/share",
    "build_root": "/__leopard_oracle__/absent",
}

# Independent ELF64 constants.  Do not derive this table from the production
# parser: the point of this test is to compare its numeric interpretation with
# GNU readelf's textual interpretation.
READELF_SECTION_TYPES = {
    "NULL": 0,
    "PROGBITS": 1,
    "SYMTAB": 2,
    "STRTAB": 3,
    "RELA": 4,
    "HASH": 5,
    "DYNAMIC": 6,
    "NOTE": 7,
    "NOBITS": 8,
    "REL": 9,
    "SHLIB": 10,
    "DYNSYM": 11,
    "INIT_ARRAY": 14,
    "FINI_ARRAY": 15,
    "PREINIT_ARRAY": 16,
    "GROUP": 17,
    "SYMTAB_SHNDX": 18,
    "RELR": 19,
    "GNU_HASH": 0x6FFFFFF6,
    "VERDEF": 0x6FFFFFFD,
    "VERNEED": 0x6FFFFFFE,
    "VERSYM": 0x6FFFFFFF,
    "X86_64_UNWIND": 0x70000001,
}
EXCLUDED_NAMES = {
    ".comment", ".note.gnu.build-id", ".strtab", ".symtab",
}
EXCLUDED_PREFIXES = (".debug", ".zdebug")


def readelf_section_table(path: pathlib.Path) -> list[dict[str, object]]:
    result = subprocess.run(
        [str(READELF), "-SWt", str(path)],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={"LANG": "C", "LC_ALL": "C"},
        check=False,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(
            f"readelf failed for {path}: rc={result.returncode} "
            f"stderr={result.stderr!r}")
    try:
        lines = result.stdout.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise RuntimeError("readelf emitted non-ASCII section output") from error

    rows: list[dict[str, object]] = []
    for line_index, line in enumerate(lines):
        match = ROW_PATTERN.fullmatch(line)
        if match is None:
            continue
        index = int(match.group(1), 10)
        name = match.group(2).strip()
        if line_index + 2 >= len(lines):
            raise RuntimeError(f"readelf emitted a truncated row: {line!r}")
        tokens = lines[line_index + 1].split()
        if len(tokens) != 8:
            raise RuntimeError(
                f"readelf emitted unexpected section details: "
                f"{lines[line_index + 1]!r}")
        type_name = tokens.pop(0)
        if type_name not in READELF_SECTION_TYPES:
            raise RuntimeError(
                f"readelf emitted an unmapped section type {type_name!r}")
        if len(tokens) != 7:
            raise RuntimeError(
                f"readelf emitted unexpected section values: {line!r}")
        address, offset, size, entsize = (
            int(token, 16) for token in tokens[:4]
        )
        flags_match = re.fullmatch(
            r"\s*\[([0-9A-Fa-f]{16})\]:.*",
            lines[line_index + 2],
        )
        if flags_match is None:
            raise RuntimeError(
                f"readelf emitted unexpected numeric section flags: "
                f"{lines[line_index + 2]!r}")
        flags = int(flags_match.group(1), 16)
        link, info, alignment = (
            int(token, 10) for token in tokens[4:]
        )
        rows.append({
            "index": index,
            "name": name,
            "type": READELF_SECTION_TYPES[type_name],
            "type_name": type_name,
            "flags": flags,
            "address": address,
            "offset": offset,
            "size": size,
            "alignment": alignment,
            "entsize": entsize,
            "link": link,
            "info": info,
        })
    if not rows or [row["index"] for row in rows] != list(range(len(rows))):
        raise RuntimeError("readelf section rows are absent or non-contiguous")
    return rows


def real_elf_paths() -> list[pathlib.Path]:
    candidates = [READELF, pathlib.Path(sys.executable).resolve()]
    for candidate in (
        pathlib.Path("/usr/lib/x86_64-linux-gnu/libc.so.6"),
        pathlib.Path("/lib/x86_64-linux-gnu/libc.so.6"),
    ):
        if candidate.is_file():
            candidates.append(candidate.resolve())
            break
    unique: list[pathlib.Path] = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file() and resolved not in unique:
            unique.append(resolved)
    return unique


def comparable_section_table(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    return [
        {key: value for key, value in row.items() if key != "type_name"}
        for row in rows
    ]


def selected_oracle_rows(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    return [
        row for row in rows
        if int(row["flags"]) & 0x2 != 0 and
        str(row["name"]) not in EXCLUDED_NAMES and
        not any(str(row["name"]).startswith(prefix)
                for prefix in EXCLUDED_PREFIXES)
    ]


def normalized_sections_from_oracle(
    data: bytes,
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    sections = []
    for row in selected_oracle_rows(rows):
        section_type = str(row["type_name"])
        if section_type == "NOBITS":
            content_sha256 = None
        else:
            offset = int(row["offset"])
            size = int(row["size"])
            content_sha256 = hashlib.sha256(
                data[offset:offset + size]).hexdigest()
        sections.append({
            "index": int(row["index"]),
            "name": str(row["name"]),
            "type": section_type,
            "flags": int(row["flags"]),
            "address": int(row["address"]),
            "size": int(row["size"]),
            "alignment": int(row["alignment"]),
            "content_sha256": content_sha256,
        })
    return sections


def path_census_from_oracle(
    data: bytes,
    rows: list[dict[str, object]],
    roots_by_role: dict[str, str] = ORACLE_ROOTS,
) -> dict[str, object]:
    retained = selected_oracle_rows(rows)
    roots = []
    for role, path in roots_by_role.items():
        needle = path.encode("utf-8")
        section_counts = []
        for row in retained:
            if row["type_name"] == "NOBITS":
                occurrences = 0
            else:
                offset = int(row["offset"])
                size = int(row["size"])
                occurrences = data[offset:offset + size].count(needle)
            section_counts.append({
                "name": str(row["name"]),
                "occurrences": occurrences,
            })
        roots.append({
            "role": role,
            "path": path,
            "occurrences": sum(
                section["occurrences"] for section in section_counts),
            "sections": section_counts,
        })
    return {
        "match_rule": "non-overlapping-exact-utf8-byte-substring/v1",
        "roots": roots,
    }


HAS_READELF = sys.platform.startswith("linux") and READELF.is_file()
HOST_IS_X86_64 = platform.machine().lower() in ("x86_64", "amd64")


class ExactMainElfOracleTest(unittest.TestCase):
    @unittest.skipUnless(HAS_READELF, "GNU readelf is unavailable")
    def test_synthetic_conventional_and_extended_tables_match_readelf(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-elf-oracle-") as temporary:
            root = pathlib.Path(temporary)
            for extended in (False, True):
                with self.subTest(extended=extended):
                    data = fixtures.synthetic_elf(
                        extended_numbering=extended)
                    path = root / (
                        "extended.elf" if extended else "conventional.elf")
                    path.write_bytes(data)
                    oracle = readelf_section_table(path)
                    parsed = contract.parse_elf64_section_table(data)
                    oracle_rows = comparable_section_table(oracle)
                    self.assertEqual(len(parsed), len(oracle_rows))
                    for row_index, (actual, expected) in enumerate(
                            zip(parsed, oracle_rows)):
                        self.assertEqual(
                            actual, expected,
                            f"section row {row_index}")
                    self.assertEqual(
                        len(parsed), 0xFF01 if extended else 7)
                    identity = (
                        contract.normalized_code_identity_from_elf_bytes(
                            data, roots=fixtures.ELF_ROOTS))
                    self.assertEqual(
                        identity["sections"],
                        normalized_sections_from_oracle(data, oracle),
                    )
                    expected_census = path_census_from_oracle(
                        data, oracle, fixtures.ELF_ROOTS)
                    self.assertEqual(
                        identity["path_string_census"], expected_census)
                    self.assertEqual(
                        [root["occurrences"]
                         for root in expected_census["roots"]],
                        [1, 1, 1],
                    )

    @unittest.skipUnless(
        HAS_READELF and HOST_IS_X86_64,
        "host-ELF differential requires Linux x86-64",
    )
    def test_parser_and_projection_match_gnu_readelf(self) -> None:
        paths = real_elf_paths()
        self.assertGreaterEqual(len(paths), 2)
        positive_occurrences = 0
        for path in paths:
            with self.subTest(path=str(path)):
                data = path.read_bytes()
                oracle = readelf_section_table(path)
                parsed = contract.parse_elf64_section_table(data)
                self.assertEqual(parsed, comparable_section_table(oracle))

                identity = contract.normalized_code_identity_from_elf_bytes(
                    data, roots=ORACLE_ROOTS)
                self.assertEqual(
                    identity["sections"],
                    normalized_sections_from_oracle(data, oracle),
                )
                expected_census = path_census_from_oracle(data, oracle)
                self.assertEqual(
                    identity["path_string_census"], expected_census)
                positive_occurrences += sum(
                    int(root["occurrences"])
                    for root in expected_census["roots"])
                self.assertEqual(
                    contract.verify_normalized_code_identity_against_elf_bytes(
                        data, identity, roots=ORACLE_ROOTS),
                    identity,
                )
        self.assertGreater(positive_occurrences, 0)


if __name__ == "__main__":
    unittest.main()
