#!/usr/bin/env python3
"""Host-independent tests for the exact-main baseline identity contract."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import pathlib
import struct
from types import MappingProxyType
import unittest
from unittest import mock


MODULE_PATH = pathlib.Path(__file__).with_name(
    "leopard2_exact_main_baseline.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_exact_main_baseline_under_test", MODULE_PATH,
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load exact-main baseline contract")
contract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contract)


ELF_ROOTS = {
    "adapter_source_root": "/src/adapter",
    "baseline_source_root": "/src/baseline",
    "build_root": "/build/main",
}


def _align_up(value: int, alignment: int) -> int:
    return (value + alignment - 1) // alignment * alignment


def synthetic_elf(
    *,
    text: bytes = b"\x90\xc3",
    rodata: bytes | None = None,
    debug: bytes = b"synthetic debug payload",
    extended_numbering: bool = False,
) -> bytes:
    if rodata is None:
        rodata = b"|".join(path.encode("ascii") for path in ELF_ROOTS.values())
    specifications = [
        {
            "name": ".note.gnu.build-id", "type": 7, "flags": 0x2,
            "address": 0x400200, "alignment": 4,
            "content": b"synthetic-build-id",
        },
        {
            "name": ".text", "type": 1, "flags": 0x6,
            "address": 0x401000, "alignment": 1, "content": text,
        },
        {
            "name": ".rodata", "type": 1, "flags": 0x2,
            "address": 0x402000, "alignment": 1, "content": rodata,
        },
        {
            "name": ".bss", "type": 8, "flags": 0x3,
            "address": 0x404000, "alignment": 32,
            "content": None, "size": 64,
        },
        {
            "name": ".debug_info", "type": 1, "flags": 0,
            "address": 0, "alignment": 1, "content": debug,
        },
    ]
    section_names = [specification["name"] for specification in specifications]
    section_names.append(".shstrtab")
    name_offsets: dict[str, int] = {}
    string_table = bytearray(b"\0")
    for name in section_names:
        name_offsets[name] = len(string_table)
        string_table.extend(name.encode("ascii") + b"\0")
    specifications.append({
        "name": ".shstrtab", "type": 3, "flags": 0,
        "address": 0, "alignment": 1, "content": bytes(string_table),
    })

    payload = bytearray(b"\0" * 64)
    retained = []
    for specification in specifications:
        content = specification["content"]
        alignment = int(specification["alignment"])
        if content is None:
            offset = len(payload)
            size = int(specification["size"])
        else:
            aligned = _align_up(len(payload), alignment)
            payload.extend(b"\0" * (aligned - len(payload)))
            offset = len(payload)
            payload.extend(content)
            size = len(content)
        retained.append({
            **specification,
            "name_offset": name_offsets[str(specification["name"])],
            "offset": offset,
            "size": size,
        })

    if extended_numbering:
        filler_count = 0xFF00 - len(retained)
        header_rows = (
            retained[:-1] + ([None] * filler_count) + retained[-1:])
    else:
        header_rows = list(retained)
    section_table_offset = _align_up(len(payload), 8)
    payload.extend(b"\0" * (section_table_offset - len(payload)))
    section_count = len(header_rows) + 1
    string_table_index = len(header_rows)
    payload.extend(b"\0" * (section_count * 64))
    encoded_count = 0 if extended_numbering else section_count
    encoded_string_index = 0xFFFF if extended_numbering else string_table_index
    ident = b"\x7fELF\x02\x01\x01" + (b"\0" * 9)
    struct.pack_into(
        "<16sHHIQQQIHHHHHH", payload, 0,
        ident, 3, 62, 1, 0x401000, 0, section_table_offset, 0,
        64, 0, 0, 64, encoded_count, encoded_string_index,
    )
    struct.pack_into(
        "<IIQQQQIIQQ", payload, section_table_offset,
        0, 0, 0, 0, 0,
        section_count if extended_numbering else 0,
        string_table_index if extended_numbering else 0,
        0, 0, 0,
    )
    for index, specification in enumerate(header_rows, 1):
        if specification is None:
            continue
        struct.pack_into(
            "<IIQQQQIIQQ", payload,
            section_table_offset + index * 64,
            int(specification["name_offset"]), int(specification["type"]),
            int(specification["flags"]), int(specification["address"]),
            int(specification["offset"]), int(specification["size"]),
            0, 0, int(specification["alignment"]), 0,
        )
    return bytes(payload)


def rebuilt_identity(
    identity: dict,
    *,
    artifact: dict | None = None,
    sections: list[dict] | None = None,
    census: dict | None = None,
) -> dict:
    return contract.normalized_code_identity_record(
        artifact=copy.deepcopy(
            identity["artifact"] if artifact is None else artifact),
        sections=copy.deepcopy(
            identity["sections"] if sections is None else sections),
        path_string_census=copy.deepcopy(
            identity["path_string_census"] if census is None else census),
    )


def section_fixture() -> list[dict]:
    return [
        {
            "index": 1,
            "name": ".interp",
            "type": "PROGBITS",
            "flags": 0x2,
            "address": 0x318,
            "size": 28,
            "alignment": 1,
            "content_sha256": "1" * 64,
        },
        {
            "index": 4,
            "name": ".text",
            "type": "PROGBITS",
            "flags": 0x6,
            "address": 0x4000,
            "size": 4096,
            "alignment": 64,
            "content_sha256": "2" * 64,
        },
        {
            "index": 5,
            "name": ".rodata",
            "type": "PROGBITS",
            "flags": 0x2,
            "address": 0x6000,
            "size": 1024,
            "alignment": 32,
            "content_sha256": "3" * 64,
        },
        {
            "index": 9,
            "name": ".bss",
            "type": "NOBITS",
            "flags": 0x3,
            "address": 0x9000,
            "size": 128,
            "alignment": 32,
            "content_sha256": None,
        },
    ]


def census_fixture(
    *,
    build_root: str = "/tmp/leopard-build-primary",
    source_root: str = "/tmp/leopard-source",
) -> dict:
    root_values = (
        ("adapter_source_root", "/home/catid/leopard"),
        ("baseline_source_root", source_root),
        ("build_root", build_root),
    )
    return {
        "match_rule": "non-overlapping-exact-utf8-byte-substring/v1",
        "roots": [{
            "role": role,
            "path": path,
            "occurrences": 0,
            "sections": [
                {"name": section["name"], "occurrences": 0}
                for section in section_fixture()
            ],
        } for role, path in root_values],
    }


def identity_fixture(
    *,
    artifact_sha256: str = "a" * 64,
    build_root: str = "/tmp/leopard-build-primary",
) -> dict:
    return contract.normalized_code_identity_record(
        artifact={"size": 1_234_567, "sha256": artifact_sha256},
        sections=section_fixture(),
        path_string_census=census_fixture(build_root=build_root),
    )


class ExactMainBaselineContractTest(unittest.TestCase):
    def assertRejected(self, value: object) -> None:  # noqa: N802
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.validate_normalized_code_identity(value)

    def test_literal_golden_and_canonical_roundtrip(self) -> None:
        identity = identity_fixture()
        self.assertEqual(
            identity["schema"],
            "leopard2-gf8-exact-main-normalized-code-identity/v1",
        )
        self.assertEqual(
            identity["combined_sha256"],
            "c67c0f4432732ecb767697e528ac9a46fe51fadd8eb08c763251eda80f33ab20",
        )
        self.assertEqual(
            identity["record_sha256"],
            "894d586cdd932c183a747f275abee5643ed237df5126e80c4220578881a542c2",
        )
        encoded = contract.canonical_json_bytes(identity)
        self.assertEqual(
            contract.load_normalized_code_identity(encoded), identity)
        self.assertEqual(encoded, contract.canonical_json_bytes(
            contract.load_normalized_code_identity(encoded)))

    def test_raw_path_variant_is_separately_bound_but_normalizes_equal(self) -> None:
        primary = identity_fixture(artifact_sha256="a" * 64)
        path_variant = identity_fixture(
            artifact_sha256="b" * 64,
            build_root="/tmp/leopard-build-path-variant",
        )
        self.assertNotEqual(primary["artifact"], path_variant["artifact"])
        self.assertNotEqual(
            primary["path_string_census"]["roots"][2]["path"],
            path_variant["path_string_census"]["roots"][2]["path"],
        )
        self.assertEqual(
            primary["combined_sha256"], path_variant["combined_sha256"])
        self.assertEqual(primary["sections"], path_variant["sections"])

        shifted_sections = section_fixture()
        for section in shifted_sections:
            section["index"] += 10
        shifted_indices = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "d" * 64},
            sections=shifted_sections,
            path_string_census=census_fixture(),
        )
        self.assertNotEqual(primary["sections"], shifted_indices["sections"])
        self.assertEqual(
            primary["combined_sha256"], shifted_indices["combined_sha256"])

        altered_sections = section_fixture()
        altered_sections[2]["content_sha256"] = "4" * 64
        content_variant = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "c" * 64},
            sections=altered_sections,
            path_string_census=census_fixture(),
        )
        self.assertNotEqual(
            primary["combined_sha256"],
            content_variant["combined_sha256"],
        )

        census_variant = contract.normalized_code_identity_record(
            artifact=primary["artifact"],
            sections=section_fixture(),
            path_string_census=census_fixture(
                build_root="/tmp/leopard-build-census-only"),
        )
        self.assertEqual(
            primary["combined_sha256"], census_variant["combined_sha256"])
        self.assertNotEqual(
            primary["record_sha256"], census_variant["record_sha256"])
        self.assertNotEqual(
            primary["path_string_census"],
            census_variant["path_string_census"],
        )

    def test_elf_byte_derivation_golden_and_reverification(self) -> None:
        data = synthetic_elf()
        sections = contract.parse_elf64_section_table(data)
        self.assertEqual(
            [section["name"] for section in sections],
            [
                "", ".note.gnu.build-id", ".text", ".rodata", ".bss",
                ".debug_info", ".shstrtab",
            ],
        )
        self.assertEqual(
            set(contract.SECTION_TYPE_NAMES.values()),
            {
                "PROGBITS", "STRTAB", "RELA", "HASH", "DYNAMIC",
                "NOTE", "NOBITS", "REL", "DYNSYM", "INIT_ARRAY",
                "FINI_ARRAY", "PREINIT_ARRAY", "RELR", "GNU_HASH",
                "VERDEF", "VERNEED", "VERSYM", "X86_64_UNWIND",
            },
        )
        self.assertEqual(len(contract.SECTION_TYPE_NAMES), 18)
        record = contract.normalized_code_identity_from_elf_bytes(
            data, roots=ELF_ROOTS)
        self.assertEqual(
            [section["name"] for section in record["sections"]],
            [".text", ".rodata", ".bss"],
        )
        self.assertEqual(
            record["artifact"],
            {"size": len(data), "sha256": hashlib.sha256(data).hexdigest()},
        )
        for root in record["path_string_census"]["roots"]:
            self.assertEqual(root["occurrences"], 1)
            self.assertEqual(
                [row["occurrences"] for row in root["sections"]],
                [0, 1, 0],
            )
        self.assertEqual(
            record["combined_sha256"],
            "08f5ab4e48a323d9543500c2f1dc2fb424246de8e989c8191a205067a999768c",
        )
        self.assertEqual(
            record["record_sha256"],
            "d886689ffb6bfe7283281fb5c1387fdc58b4e4f625180b16ec4b7813a47587a5",
        )
        self.assertEqual(
            contract.verify_normalized_code_identity_against_elf_bytes(
                data, record, roots=ELF_ROOTS),
            record,
        )
        self.assertEqual(
            contract.path_string_census_from_elf(data, sections, ELF_ROOTS),
            record["path_string_census"],
        )
        encoded = contract.canonical_json_bytes(record)
        self.assertEqual(contract.load_normalized_code_identity(encoded), record)

    def test_elf_derivation_hashes_retained_bytes_once(self) -> None:
        data = synthetic_elf()
        parsed = contract.parse_elf64_section_table(data)
        retained_bytes = []
        for row in parsed:
            if row["name"] in (".text", ".rodata"):
                offset = row["offset"]
                retained_bytes.append(data[offset:offset + row["size"]])
        self.assertEqual(len(retained_bytes), 2)

        original_sha256 = contract.hashlib.sha256
        hashed_inputs = []

        def recording_sha256(value: object = b""):
            hashed_inputs.append(bytes(value))
            return original_sha256(value)

        with mock.patch.object(
                contract.hashlib, "sha256", side_effect=recording_sha256):
            contract.normalized_code_identity_from_elf_bytes(
                data, roots=ELF_ROOTS)
        self.assertEqual(
            [hashed_inputs.count(content) for content in retained_bytes],
            [1, 1],
        )

    def test_elf_reverification_rejects_self_consistent_forged_records(self) -> None:
        data = synthetic_elf()
        record = contract.normalized_code_identity_from_elf_bytes(
            data, roots=ELF_ROOTS)
        forgeries = []

        deleted_sections = copy.deepcopy(record["sections"])
        del deleted_sections[1]
        deleted_census = copy.deepcopy(record["path_string_census"])
        for root in deleted_census["roots"]:
            del root["sections"][1]
            root["occurrences"] = 0
        forgeries.append(rebuilt_identity(
            record, sections=deleted_sections, census=deleted_census))

        appended_sections = copy.deepcopy(record["sections"])
        appended_sections.append({
            "index": 5,
            "name": ".invented",
            "type": "PROGBITS",
            "flags": 0x2,
            "address": 0x405000,
            "size": 1,
            "alignment": 1,
            "content_sha256": "f" * 64,
        })
        appended_census = copy.deepcopy(record["path_string_census"])
        for root in appended_census["roots"]:
            root["sections"].append({
                "name": ".invented", "occurrences": 0})
        forgeries.append(rebuilt_identity(
            record, sections=appended_sections, census=appended_census))

        altered_sections = copy.deepcopy(record["sections"])
        altered_sections[1]["content_sha256"] = "4" * 64
        forgeries.append(rebuilt_identity(record, sections=altered_sections))
        forgeries.append(rebuilt_identity(
            record,
            artifact={"size": len(data), "sha256": "b" * 64},
        ))

        altered_census = copy.deepcopy(record["path_string_census"])
        altered_census["roots"][0]["occurrences"] = 0
        altered_census["roots"][0]["sections"][1]["occurrences"] = 0
        forgeries.append(rebuilt_identity(record, census=altered_census))

        altered_root_census = copy.deepcopy(record["path_string_census"])
        altered_root = altered_root_census["roots"][2]
        altered_root["path"] = "/build/other"
        altered_root["occurrences"] = 0
        for row in altered_root["sections"]:
            row["occurrences"] = 0
        forgeries.append(rebuilt_identity(record, census=altered_root_census))

        for index, forgery in enumerate(forgeries):
            with self.subTest(index=index):
                self.assertEqual(
                    contract.validate_normalized_code_identity(forgery),
                    forgery,
                )
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.verify_normalized_code_identity_against_elf_bytes(
                        data, forgery, roots=ELF_ROOTS)

        reordered = copy.deepcopy(record)
        reordered["sections"][0], reordered["sections"][1] = (
            reordered["sections"][1], reordered["sections"][0])
        for root in reordered["path_string_census"]["roots"]:
            root["sections"][0], root["sections"][1] = (
                root["sections"][1], root["sections"][0])
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.verify_normalized_code_identity_against_elf_bytes(
                data, reordered, roots=ELF_ROOTS)

        tampered = bytearray(data)
        tamper_offset = data.index(b"/build/main")
        tampered[tamper_offset + 1] ^= 1
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.verify_normalized_code_identity_against_elf_bytes(
                bytes(tampered), record, roots=ELF_ROOTS)

    def test_extended_numbering_and_path_variant_semantics(self) -> None:
        conventional_data = synthetic_elf()
        extended_data = synthetic_elf(extended_numbering=True)
        conventional = contract.normalized_code_identity_from_elf_bytes(
            conventional_data, roots=ELF_ROOTS)
        extended = contract.normalized_code_identity_from_elf_bytes(
            extended_data, roots=ELF_ROOTS)
        self.assertNotEqual(
            conventional["artifact"]["sha256"],
            extended["artifact"]["sha256"],
        )
        self.assertEqual(
            conventional["combined_sha256"], extended["combined_sha256"])

        selected_a = synthetic_elf(rodata=b"selected:/build/path-a")
        selected_b = synthetic_elf(rodata=b"selected:/build/path-b")
        roots_a = {**ELF_ROOTS, "build_root": "/build/path-a"}
        roots_b = {**ELF_ROOTS, "build_root": "/build/path-b"}
        selected_record_a = contract.normalized_code_identity_from_elf_bytes(
            selected_a, roots=roots_a)
        selected_record_b = contract.normalized_code_identity_from_elf_bytes(
            selected_b, roots=roots_b)
        self.assertNotEqual(
            selected_record_a["artifact"]["sha256"],
            selected_record_b["artifact"]["sha256"],
        )
        self.assertNotEqual(
            selected_record_a["combined_sha256"],
            selected_record_b["combined_sha256"],
        )

        excluded_a = synthetic_elf(
            rodata=b"stable", debug=b"excluded:/build/path-a")
        excluded_b = synthetic_elf(
            rodata=b"stable", debug=b"excluded:/build/path-b")
        excluded_record_a = contract.normalized_code_identity_from_elf_bytes(
            excluded_a, roots=roots_a)
        excluded_record_b = contract.normalized_code_identity_from_elf_bytes(
            excluded_b, roots=roots_b)
        self.assertNotEqual(
            excluded_record_a["artifact"]["sha256"],
            excluded_record_b["artifact"]["sha256"],
        )
        self.assertEqual(
            excluded_record_a["combined_sha256"],
            excluded_record_b["combined_sha256"],
        )
        self.assertEqual(
            excluded_record_a["path_string_census"]["roots"][2][
                "occurrences"],
            0,
        )

    def test_elf_path_census_nonoverlap_straddle_and_nobits(self) -> None:
        repeated_roots = {
            "adapter_source_root": "/a/a",
            "baseline_source_root": "/baseline",
            "build_root": "/build",
        }
        repeated = contract.normalized_code_identity_from_elf_bytes(
            synthetic_elf(rodata=b"/a/a/a/a"), roots=repeated_roots)
        self.assertEqual(
            repeated["path_string_census"]["roots"][0]["occurrences"],
            2,
        )
        self.assertEqual(
            [root["sections"][-1]["occurrences"]
             for root in repeated["path_string_census"]["roots"]],
            [0, 0, 0],
        )

        straddled_data = synthetic_elf(
            text=b"\x90/src/adap", rodata=b"ter")
        self.assertGreaterEqual(straddled_data.find(b"/src/adapter"), 0)
        straddled = contract.normalized_code_identity_from_elf_bytes(
            straddled_data, roots=ELF_ROOTS)
        self.assertEqual(
            straddled["path_string_census"]["roots"][0]["occurrences"],
            0,
        )

    def test_malformed_elf_section_tables_fail_closed(self) -> None:
        data = synthetic_elf()
        section_table_offset = struct.unpack_from("<Q", data, 40)[0]

        def changed(offset: int, field_format: str, value: int) -> bytes:
            payload = bytearray(data)
            struct.pack_into(field_format, payload, offset, value)
            return bytes(payload)

        malformed = []
        malformed.append(b"short")
        bad_magic = bytearray(data)
        bad_magic[:4] = b"NOPE"
        malformed.append(bytes(bad_magic))
        malformed.append(changed(4, "<B", 1))
        malformed.append(changed(5, "<B", 2))
        malformed.append(changed(6, "<B", 2))
        malformed.append(changed(16, "<H", 1))
        malformed.append(changed(18, "<H", 3))
        malformed.append(changed(20, "<I", 2))
        malformed.append(changed(52, "<H", 63))
        malformed.append(changed(58, "<H", 63))
        malformed.append(changed(40, "<Q", len(data)))
        malformed.append(changed(60, "<H", 4_097))
        malformed.append(changed(60, "<H", 8))
        malformed.append(changed(62, "<H", 0))
        malformed.append(changed(62, "<H", 7))
        malformed.append(changed(
            section_table_offset + 6 * 64 + 4, "<I", 1))
        malformed.append(changed(
            section_table_offset + 2 * 64, "<I", 0xFFFFFFFF))
        malformed.append(changed(60, "<H", 0))
        malformed.append(changed(60, "<H", 0xFF00))
        malformed.append(changed(
            section_table_offset + 32, "<Q", 7))
        malformed.append(changed(62, "<H", 0xFFFF))
        malformed.append(changed(62, "<H", 0xFF00))
        malformed.append(changed(
            section_table_offset + 40, "<I", 6))
        malformed.append(changed(
            section_table_offset + 3 * 64 + 32, "<Q", 1 << 63))
        malformed.append(changed(
            section_table_offset + 3 * 64 + 48, "<Q", (1 << 63) + 1))

        string_offset = struct.unpack_from(
            "<Q", data, section_table_offset + 6 * 64 + 24)[0]
        string_size = struct.unpack_from(
            "<Q", data, section_table_offset + 6 * 64 + 32)[0]
        nonempty_zero_name = bytearray(data)
        nonempty_zero_name[string_offset] = ord("x")
        malformed.append(bytes(nonempty_zero_name))

        unterminated = bytearray(data)
        unterminated[string_offset + string_size - 1] = ord("x")
        malformed.append(bytes(unterminated))

        non_ascii_name = bytearray(data)
        text_name = data.find(b".text\0", string_offset,
                              string_offset + string_size)
        self.assertGreaterEqual(text_name, 0)
        non_ascii_name[text_name] = 0xFF
        malformed.append(bytes(non_ascii_name))
        malformed.append(changed(
            section_table_offset + 3 * 64 + 32, "<Q", len(data)))
        text_offset = struct.unpack_from(
            "<Q", data, section_table_offset + 2 * 64 + 24)[0]
        malformed.append(changed(section_table_offset + 8, "<Q", 1))

        for index, payload in enumerate(malformed):
            with self.subTest(index=index):
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.parse_elf64_section_table(payload)

        extended_data = synthetic_elf(extended_numbering=True)
        extended_section_table_offset = struct.unpack_from(
            "<Q", extended_data, 40)[0]
        oversized_extended = bytearray(extended_data)
        struct.pack_into(
            "<Q", oversized_extended,
            extended_section_table_offset + 32, 65_537)
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.parse_elf64_section_table(bytes(oversized_extended))

        normalized_malformed = (
            changed(section_table_offset + 3 * 64 + 24, "<Q", text_offset),
            changed(section_table_offset + 3 * 64 + 4, "<I", 0x70000002),
        )
        for index, payload in enumerate(normalized_malformed):
            with self.subTest(normalized_index=index):
                with mock.patch.object(
                        contract.hashlib, "sha256",
                        wraps=contract.hashlib.sha256) as sha256:
                    with self.assertRaises(contract.ExactMainBaselineError):
                        contract.normalized_code_identity_from_elf_bytes(
                            payload, roots=ELF_ROOTS)
                    self.assertEqual(sha256.call_count, 0)

        control_name_unknown_type = bytearray(data)
        rodata_name = data.find(
            b".rodata\0", string_offset, string_offset + string_size)
        self.assertGreaterEqual(rodata_name, 0)
        control_name_unknown_type[rodata_name:rodata_name + 7] = \
            b"\x1b[2J.ev"
        struct.pack_into(
            "<I", control_name_unknown_type,
            section_table_offset + 3 * 64 + 4, 0x70000002)
        with self.assertRaises(contract.ExactMainBaselineError) as error:
            contract.normalized_code_identity_from_elf_bytes(
                bytes(control_name_unknown_type), roots=ELF_ROOTS)
        self.assertNotIn("\x1b", str(error.exception))
        self.assertIn("\\x1b", str(error.exception))

        old_maximum = contract.MAX_ELF_INPUT_BYTES
        try:
            contract.MAX_ELF_INPUT_BYTES = len(data) - 1
            with self.assertRaises(contract.ExactMainBaselineError):
                contract.parse_elf64_section_table(data)
        finally:
            contract.MAX_ELF_INPUT_BYTES = old_maximum

        for format_name in (
                "ELF64_HEADER_FORMAT", "ELF64_SECTION_HEADER_FORMAT"):
            with self.subTest(format_name=format_name):
                original_format = getattr(contract, format_name)
                try:
                    setattr(contract, format_name, "<Q")
                    with self.assertRaises(contract.ExactMainBaselineError):
                        contract.parse_elf64_section_table(data)
                finally:
                    setattr(contract, format_name, original_format)

        original_text_length = contract.MAX_TEXT_LENGTH
        try:
            contract.MAX_TEXT_LENGTH = 3
            with self.assertRaises(contract.ExactMainBaselineError):
                contract.parse_elf64_section_table(data)
        finally:
            contract.MAX_TEXT_LENGTH = original_text_length

        original_maximum_sections = contract.MAX_SECTION_COUNT
        try:
            contract.MAX_SECTION_COUNT = 2
            with mock.patch.object(
                    contract.hashlib, "sha256",
                    wraps=contract.hashlib.sha256) as sha256:
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.normalized_code_identity_from_elf_bytes(
                        data, roots=ELF_ROOTS)
                self.assertEqual(sha256.call_count, 0)
        finally:
            contract.MAX_SECTION_COUNT = original_maximum_sections

        original_type_names = contract.SECTION_TYPE_NAMES
        try:
            contract.SECTION_TYPE_NAMES = MappingProxyType({
                number: name
                for number, name in original_type_names.items()
                if number != 1
            })
            with self.assertRaises(contract.ExactMainBaselineError):
                contract.normalized_code_identity_from_elf_bytes(
                    data, roots=ELF_ROOTS)
        finally:
            contract.SECTION_TYPE_NAMES = original_type_names

        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_from_elf_bytes(
                data, roots={"unexpected": "/x"})
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_from_elf_bytes(
                data,
                roots={
                    **ELF_ROOTS,
                    "build_root": ELF_ROOTS["adapter_source_root"] + "/sub",
                },
            )

        nobits_beyond_eof = changed(
            section_table_offset + 4 * 64 + 24, "<Q", len(data) + 4_096)
        nobits_record = contract.normalized_code_identity_from_elf_bytes(
            nobits_beyond_eof, roots=ELF_ROOTS)
        self.assertEqual(nobits_record["sections"][-1]["type"], "NOBITS")

        parsed = contract.parse_elf64_section_table(data)
        forged_table = copy.deepcopy(parsed)
        del forged_table[1]
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.path_string_census_from_elf(
                data, forged_table, ELF_ROOTS)

    def test_constructor_and_validator_detach_inputs(self) -> None:
        artifact = {"size": 1_234_567, "sha256": "a" * 64}
        sections = section_fixture()
        census = census_fixture()
        identity = contract.normalized_code_identity_record(
            artifact=artifact,
            sections=sections,
            path_string_census=census,
        )
        artifact["size"] = 7
        sections[1]["name"] = ".changed"
        census["roots"][0]["sections"][0]["name"] = ".changed"
        self.assertEqual(identity, contract.validate_normalized_code_identity(
            identity))

        validated = contract.validate_normalized_code_identity(identity)
        validated["sections"][0]["name"] = ".changed"
        self.assertEqual(identity["sections"][0]["name"], ".interp")

    def test_exact_keys_types_and_cross_bindings(self) -> None:
        identity = identity_fixture()
        mutations = []

        extra = copy.deepcopy(identity)
        extra["extra"] = None
        mutations.append(extra)

        missing = copy.deepcopy(identity)
        del missing["artifact"]
        mutations.append(missing)

        bool_size = copy.deepcopy(identity)
        bool_size["artifact"]["size"] = True
        mutations.append(bool_size)

        changed_rule = copy.deepcopy(identity)
        changed_rule["selection_rule"]["required_flag"] = "SHF_EXECINSTR"
        mutations.append(changed_rule)

        changed_digest = copy.deepcopy(identity)
        changed_digest["combined_sha256"] = "f" * 64
        mutations.append(changed_digest)

        inconsistent_total = copy.deepcopy(identity)
        inconsistent_total["path_string_census"]["roots"][2][
            "occurrences"] = 2
        mutations.append(inconsistent_total)

        stale_artifact_binding = copy.deepcopy(identity)
        stale_artifact_binding["artifact"]["sha256"] = "b" * 64
        mutations.append(stale_artifact_binding)

        stale_census_binding = copy.deepcopy(identity)
        stale_census_binding["path_string_census"]["roots"][2]["path"] = \
            "/tmp/leopard-build-other"
        mutations.append(stale_census_binding)

        slash_root = copy.deepcopy(identity)
        slash_root["path_string_census"]["roots"][0]["path"] = "/"
        mutations.append(slash_root)

        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(mutation)

    def test_census_root_paths_are_canonical_display_safe_and_bounded(
        self,
    ) -> None:
        structurally_invalid_paths = (
            "",
            "/",
            "tmp/leopard-build-primary",
            "/tmp//leopard-build-primary",
            "/tmp/./leopard-build-primary",
            "/tmp/zz/../leopard-build-primary",
            "/" + ("x" * 4_097),
        )
        for path in structurally_invalid_paths:
            with self.subTest(path=repr(path)):
                census = census_fixture()
                census["roots"][2]["path"] = path
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.normalized_code_identity_record(
                        artifact={"size": 1_234_567, "sha256": "a" * 64},
                        sections=section_fixture(),
                        path_string_census=census,
                    )

        for path in (
            "/tmp/caf\u00e9",
            "/tmp/leo\u200bpard",
            "/tmp/leo\u00a0pard",
            "/tmp/leo\ufe0fpard",
        ):
            with self.subTest(non_ascii_path=repr(path)):
                census = census_fixture()
                census["roots"][2]["path"] = path
                with self.assertRaisesRegex(
                        contract.ExactMainBaselineError,
                        "portable display-safe path"):
                    contract.normalized_code_identity_record(
                        artifact={"size": 1_234_567, "sha256": "a" * 64},
                        sections=section_fixture(),
                        path_string_census=census,
                    )

    def test_hash_section_name_and_artifact_size_canonical_forms(self) -> None:
        for digest in (("2" * 63) + "A", "2" * 63):
            with self.subTest(digest=digest):
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.normalized_code_identity_record(
                        artifact={"size": 1_234_567, "sha256": digest},
                        sections=section_fixture(),
                        path_string_census=census_fixture(),
                    )

        invalid_name_sections = section_fixture()
        invalid_name_sections[0]["name"] = "interp/alias"
        invalid_name_census = census_fixture()
        for root in invalid_name_census["roots"]:
            root["sections"][0]["name"] = "interp/alias"
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 1_234_567, "sha256": "a" * 64},
                sections=invalid_name_sections,
                path_string_census=invalid_name_census,
            )

        zero_file_backed_sizes = section_fixture()
        for section in zero_file_backed_sizes:
            if section["type"] != "NOBITS":
                section["size"] = 0
                section["content_sha256"] = (
                    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                )
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 63, "sha256": "a" * 64},
                sections=zero_file_backed_sizes,
                path_string_census=census_fixture(),
            )

    def test_section_selection_order_and_content_rules(self) -> None:
        mutations: list[tuple[list[dict], dict]] = []

        reordered = section_fixture()
        reordered[0], reordered[1] = reordered[1], reordered[0]
        reordered_census = census_fixture()
        for root in reordered_census["roots"]:
            root["sections"][0], root["sections"][1] = (
                root["sections"][1], root["sections"][0],
            )
        mutations.append((reordered, reordered_census))

        duplicate_index = section_fixture()
        duplicate_index[1]["index"] = duplicate_index[0]["index"]
        mutations.append((duplicate_index, census_fixture()))

        duplicate_name = section_fixture()
        duplicate_name[0]["name"] = ".text"
        duplicate_census = census_fixture()
        for root in duplicate_census["roots"]:
            root["sections"][0]["name"] = ".text"
        mutations.append((duplicate_name, duplicate_census))

        non_alloc = section_fixture()
        non_alloc[0]["flags"] = 0
        mutations.append((non_alloc, census_fixture()))

        build_id = section_fixture()
        build_id[0]["name"] = ".note.gnu.build-id"
        build_id_census = census_fixture()
        for root in build_id_census["roots"]:
            root["sections"][0]["name"] = ".note.gnu.build-id"
        mutations.append((build_id, build_id_census))

        debug = section_fixture()
        debug[0]["name"] = ".debug_info"
        debug_census = census_fixture()
        for root in debug_census["roots"]:
            root["sections"][0]["name"] = ".debug_info"
        mutations.append((debug, debug_census))

        missing_text = section_fixture()
        missing_text[1]["name"] = ".init"
        missing_text_census = census_fixture()
        for root in missing_text_census["roots"]:
            root["sections"][1]["name"] = ".init"
        mutations.append((missing_text, missing_text_census))

        inert_text = section_fixture()
        inert_text[1]["flags"] = 0x2
        mutations.append((inert_text, census_fixture()))

        nobits_digest = section_fixture()
        nobits_digest[-1]["content_sha256"] = "4" * 64
        mutations.append((nobits_digest, census_fixture()))

        file_backed_null = section_fixture()
        file_backed_null[0]["content_sha256"] = None
        mutations.append((file_backed_null, census_fixture()))

        empty_with_wrong_digest = section_fixture()
        empty_with_wrong_digest[0]["size"] = 0
        mutations.append((empty_with_wrong_digest, census_fixture()))

        nonempty_with_empty_digest = section_fixture()
        nonempty_with_empty_digest[0]["content_sha256"] = (
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        )
        mutations.append((nonempty_with_empty_digest, census_fixture()))

        prefixed_nobits = section_fixture()
        prefixed_nobits[-1]["type"] = "SHT_NOBITS"
        prefixed_nobits[-1]["content_sha256"] = "4" * 64
        mutations.append((prefixed_nobits, census_fixture()))

        executable_nobits = section_fixture()
        executable_nobits[1]["type"] = "NOBITS"
        executable_nobits[1]["flags"] = 0x7
        executable_nobits[1]["content_sha256"] = None
        mutations.append((executable_nobits, census_fixture()))

        readonly_nobits = section_fixture()
        readonly_nobits[2]["type"] = "NOBITS"
        readonly_nobits[2]["content_sha256"] = None
        mutations.append((readonly_nobits, census_fixture()))

        bad_alignment = section_fixture()
        bad_alignment[1]["alignment"] = 3
        mutations.append((bad_alignment, census_fixture()))

        misaligned_address = section_fixture()
        misaligned_address[1]["address"] += 1
        mutations.append((misaligned_address, census_fixture()))

        nobits_census = census_fixture()
        nobits_census["roots"][2]["sections"][-1]["occurrences"] = 1
        nobits_census["roots"][2]["occurrences"] += 1
        mutations.append((section_fixture(), nobits_census))

        overflow = section_fixture()
        overflow[-1]["address"] = (1 << 64) - 64
        overflow[-1]["size"] = 128
        mutations.append((overflow, census_fixture()))

        impossible_count = census_fixture()
        impossible_count["roots"][2]["occurrences"] = 100
        impossible_count["roots"][2]["sections"][0]["occurrences"] = 100
        mutations.append((section_fixture(), impossible_count))

        for index, (sections, census) in enumerate(mutations):
            with self.subTest(index=index):
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.normalized_code_identity_record(
                        artifact={"size": 1_234_567, "sha256": "a" * 64},
                        sections=sections,
                        path_string_census=census,
                    )

        empty_file_backed = section_fixture()
        empty_file_backed[0]["size"] = 0
        empty_file_backed[0]["content_sha256"] = (
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        )
        empty_record = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "a" * 64},
            sections=empty_file_backed,
            path_string_census=census_fixture(),
        )
        self.assertEqual(empty_record["sections"][0]["size"], 0)

        zero_alignment_sections = section_fixture()
        zero_alignment_sections[0]["alignment"] = 0
        zero_alignment = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "a" * 64},
            sections=zero_alignment_sections,
            path_string_census=census_fixture(),
        )
        self.assertEqual(zero_alignment["sections"][0]["alignment"], 0)

        non_dot_sections = section_fixture()
        non_dot_sections[2]["name"] = "rodata.cst32"
        non_dot_census = census_fixture()
        for root in non_dot_census["roots"]:
            root["sections"][2]["name"] = "rodata.cst32"
        non_dot = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "a" * 64},
            sections=non_dot_sections,
            path_string_census=non_dot_census,
        )
        self.assertEqual(non_dot["sections"][2]["name"], "rodata.cst32")

        relr_sections = section_fixture()
        relr_sections[0]["name"] = ".relr.dyn"
        relr_sections[0]["type"] = "RELR"
        relr_census = census_fixture()
        for root in relr_census["roots"]:
            root["sections"][0]["name"] = ".relr.dyn"
        relr = contract.normalized_code_identity_record(
            artifact={"size": 1_234_567, "sha256": "a" * 64},
            sections=relr_sections,
            path_string_census=relr_census,
        )
        self.assertEqual(relr["sections"][0]["type"], "RELR")

        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 64, "sha256": "a" * 64},
                sections=section_fixture(),
                path_string_census=census_fixture(),
            )

    def test_strict_json_rejects_noncanonical_inputs(self) -> None:
        identity = identity_fixture()
        encoded = contract.canonical_json_bytes(identity)
        self.assertEqual(contract.strict_json_loads(encoded), identity)
        duplicate_digest = b'{"combined_sha256":"' + (b"f" * 64) + b'",' + encoded[1:]
        cases = (
            b'{"x":1,"x":2}',
            duplicate_digest,
            encoded + encoded,
            b'{"x":NaN}',
            b'{"x":Infinity}',
            b'{"x":1e309}',
            b'"bad\\ud800"',
            b'"bad\xed\xa0\x80"',
            (b"[" * 10_000) + (b"]" * 10_000),
            b"",
            b"\xff",
        )
        for index, payload in enumerate(cases):
            with self.subTest(index=index):
                with self.assertRaises(contract.ExactMainBaselineError):
                    contract.strict_json_loads(payload)

    def test_canonical_encoder_rejects_coercions_and_unsafe_values(self) -> None:
        invalid_values = (
            "bad\ud800",
            float("nan"),
            float("inf"),
            (1, 2),
            {"k": (1, 2)},
            {1: "a"},
            {True: "a"},
            {1.0: "a"},
            {None: "a"},
        )
        for function in (
            contract.canonical_json_bytes,
            contract.canonical_sha256,
        ):
            for value in invalid_values:
                with self.subTest(function=function.__name__, value=repr(value)):
                    with self.assertRaises(contract.ExactMainBaselineError):
                        function(value)
        self.assertEqual(
            contract.canonical_json_bytes(
                {"bool": True, "integer": 1, "null": None}),
            b'{"bool":true,"integer":1,"null":null}\n',
        )

    def test_record_builder_requires_exact_mapping_inputs(self) -> None:
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact=[],
                sections=section_fixture(),
                path_string_census=census_fixture(),
            )
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 1, "sha256": "a" * 64},
                sections="not-sections",
                path_string_census=census_fixture(),
            )
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 1_234_567, "sha256": "a" * 64},
                sections=[1],
                path_string_census=census_fixture(),
            )
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact=MappingProxyType({
                    "size": 1_234_567, "sha256": "a" * 64,
                }),
                sections=section_fixture(),
                path_string_census=census_fixture(),
            )

        nested: object = 0
        for unused in range(2_000):
            nested = {"nested": nested}
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": nested, "sha256": "a" * 64},
                sections=section_fixture(),
                path_string_census=census_fixture(),
            )
        overlapping_roots = census_fixture(
            source_root="/home/catid/leopard-baseline")
        with self.assertRaises(contract.ExactMainBaselineError):
            contract.normalized_code_identity_record(
                artifact={"size": 1_234_567, "sha256": "a" * 64},
                sections=section_fixture(),
                path_string_census=overlapping_roots,
            )


if __name__ == "__main__":
    unittest.main()
