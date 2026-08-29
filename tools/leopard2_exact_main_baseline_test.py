#!/usr/bin/env python3
"""Host-independent tests for the exact-main baseline identity contract."""

from __future__ import annotations

import copy
import importlib.util
import pathlib
from types import MappingProxyType
import unittest


MODULE_PATH = pathlib.Path(__file__).with_name(
    "leopard2_exact_main_baseline.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_exact_main_baseline_under_test", MODULE_PATH,
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load exact-main baseline contract")
contract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contract)


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
        invisible_format_code_points = (
            0x00AD, 0x0600, 0x0605, 0x061C, 0x06DD, 0x070F,
            0x0890, 0x0891, 0x08E2, 0x180E, 0x200B, 0x200F,
            0x2028, 0x202E, 0x2060, 0x2064, 0x2066, 0x206F,
            0xFEFF, 0xFFF9, 0xFFFB, 0x110BD, 0x110CD, 0x13430,
            0x1343F, 0x1BCA0, 0x1BCA3, 0x1D173, 0x1D17A,
            0xE0001, 0xE0020, 0xE007F,
            0x034F, 0x115F, 0x1160, 0x17B4, 0x180B, 0x180F,
            0x2065, 0x2800, 0x3164, 0xFE00, 0xFE0F, 0xFFA0,
            0x00A0, 0x2007, 0xE0000, 0xE0100,
        )
        invalid_paths = (
            "",
            "/",
            "tmp/leopard-build-primary",
            "/tmp//leopard-build-primary",
            "/tmp/./leopard-build-primary",
            "/tmp/zz/../leopard-build-primary",
            "/tmp/caf\u00e9",
            "/" + ("x" * 4_097),
        ) + tuple(
            "/tmp/leo" + chr(code_point) + "pard"
            for code_point in invisible_format_code_points
        )
        for path in invalid_paths:
            with self.subTest(path=repr(path)):
                census = census_fixture()
                census["roots"][2]["path"] = path
                with self.assertRaises(contract.ExactMainBaselineError):
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
