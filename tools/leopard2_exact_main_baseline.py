#!/usr/bin/env python3
"""Pure contracts for reproducible exact-main baseline evidence.

This module deliberately performs no process execution, host discovery, file I/O,
or build acquisition.  It turns an already-retained ELF section projection into a
closed normalized-code identity and independently validates that identity.  The
build/acquisition envelope is layered on this contract separately so path-sensitive
raw ELF hashes are never mistaken for the portable qualification identity.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import struct
from pathlib import PurePosixPath
from types import MappingProxyType
from typing import Any, Mapping, NoReturn, Sequence


NORMALIZED_CODE_IDENTITY_SCHEMA = \
    "leopard2-gf8-exact-main-normalized-code-identity/v1"
NORMALIZED_CODE_ALGORITHM = \
    "sha256-canonical-elf64-shf-alloc-section-projection/v1"

ELF_CLASS = "ELF64"
ELF_ENDIANNESS = "little"
ELF_MACHINE = "x86-64"
SHF_ALLOC = 0x2
SHF_EXECINSTR = 0x4
SHF_WRITE = 0x1

# These are parser/resource bounds, not evidence thresholds.
MAX_JSON_BYTES = 64 << 20
MAX_ELF_INPUT_BYTES = 256 << 20
MAX_ARTIFACT_BYTES = (1 << 63) - 1
MAX_SECTION_COUNT = 4096
MAX_ELF_SECTION_TABLE_COUNT = 1 << 16
MAX_SECTION_INDEX = (1 << 31) - 1
MAX_SECTION_ADDRESS = (1 << 64) - 1
MAX_SECTION_SIZE = (1 << 63) - 1
MAX_SECTION_FLAGS = (1 << 64) - 1
MAX_SECTION_ALIGNMENT = 1 << 63
MAX_PATH_OCCURRENCES = (1 << 31) - 1
MAX_TEXT_LENGTH = 4096

EMPTY_CONTENT_SHA256 = \
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

IDENTITY_KEYS = frozenset((
    "schema", "artifact", "elf_class", "endianness", "machine",
    "selection_rule", "sections", "path_string_census",
    "combined_sha256", "record_sha256",
))
ARTIFACT_KEYS = frozenset(("size", "sha256"))
SELECTION_RULE_KEYS = frozenset((
    "algorithm", "required_flag", "excluded_names",
    "excluded_name_prefixes", "order", "nobits",
))
SECTION_KEYS = frozenset((
    "index", "name", "type", "flags", "address", "size", "alignment",
    "content_sha256",
))
CENSUS_KEYS = frozenset(("match_rule", "roots"))
CENSUS_ROOT_KEYS = frozenset(("role", "path", "occurrences", "sections"))
CENSUS_SECTION_KEYS = frozenset(("name", "occurrences"))

EXCLUDED_SECTION_NAMES = (
    ".comment", ".note.gnu.build-id", ".strtab", ".symtab",
)
EXCLUDED_SECTION_PREFIXES = (".debug", ".zdebug")
PATH_CENSUS_MATCH_RULE = "non-overlapping-exact-utf8-byte-substring/v1"
PATH_CENSUS_ROOT_ROLES = (
    "adapter_source_root", "baseline_source_root", "build_root",
)
_LOWER_HEX_64 = re.compile(r"[0-9a-f]{64}")
_SECTION_NAME = re.compile(r"[A-Za-z0-9_.$+-]+")

# GNU readelf and llvm-readelf spell these without an SHT_ prefix.  Keep the
# domain closed because NOBITS changes content, size, and census semantics.
ALLOC_SECTION_TYPES = frozenset((
    "PROGBITS", "NOBITS", "NOTE", "RELA", "REL", "RELR", "DYNSYM", "STRTAB",
    "GNU_HASH", "HASH", "VERSYM", "VERNEED", "VERDEF", "INIT_ARRAY",
    "FINI_ARRAY", "PREINIT_ARRAY", "DYNAMIC", "X86_64_UNWIND",
))

ELF64_HEADER_SIZE = 64
ELF64_SECTION_HEADER_SIZE = 64
ELF64_HEADER_FORMAT = "<16sHHIQQQIHHHHHH"
ELF64_SECTION_HEADER_FORMAT = "<IIQQQQIIQQ"
ET_EXEC = 2
ET_DYN = 3
EM_X86_64 = 62
SHN_XINDEX = 0xFFFF
SHN_LORESERVE = 0xFF00
SHT_NULL = 0
SHT_STRTAB = 3
SHT_NOBITS = 8
SECTION_TYPE_NAMES = MappingProxyType({
    1: "PROGBITS",
    3: "STRTAB",
    4: "RELA",
    5: "HASH",
    6: "DYNAMIC",
    7: "NOTE",
    8: "NOBITS",
    9: "REL",
    11: "DYNSYM",
    14: "INIT_ARRAY",
    15: "FINI_ARRAY",
    16: "PREINIT_ARRAY",
    19: "RELR",
    0x6FFFFFF6: "GNU_HASH",
    0x6FFFFFFD: "VERDEF",
    0x6FFFFFFE: "VERNEED",
    0x6FFFFFFF: "VERSYM",
    0x70000001: "X86_64_UNWIND",
})
ELF_SECTION_KEYS = frozenset((
    "index", "name", "type", "flags", "address", "offset", "size",
    "alignment", "entsize", "link", "info",
))


class ExactMainBaselineError(ValueError):
    """An exact-main baseline object violates the pure contract."""


def _fail(message: str) -> NoReturn:
    raise ExactMainBaselineError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _require_dict(
    value: Any,
    expected_keys: frozenset[str],
    label: str,
) -> dict[str, Any]:
    _require(type(value) is dict, f"{label} is not an exact JSON object")
    _require(set(value) == expected_keys, f"{label} has an unexpected key set")
    return value


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int, f"{label} is not an exact integer")
    _require(minimum <= value <= maximum,
             f"{label} is outside its structural bound")
    return value


def _bounded_text(value: Any, label: str) -> str:
    _require(type(value) is str and 0 < len(value) <= MAX_TEXT_LENGTH,
             f"{label} is not a bounded non-empty string")
    _require(all(
        ord(character) >= 0x20 and
        not (0x7F <= ord(character) <= 0x9F) and
        not (0xD800 <= ord(character) <= 0xDFFF)
        for character in value
    ), f"{label} contains a non-display-safe character")
    return value


def _sha256(value: Any, label: str) -> str:
    _require(type(value) is str and _LOWER_HEX_64.fullmatch(value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _absolute_posix_path(value: Any, label: str) -> str:
    path = _bounded_text(value, label)
    # Acquisition must stage under a portable printable-ASCII root.  The path is
    # both a byte-search needle and a human-reviewed evidence label, so accepting
    # visually ignorable or Unicode-confusable spellings would be unsafe.
    _require(all(0x21 <= ord(character) <= 0x7E for character in path),
             f"{label} is not a portable display-safe path")
    _require(path != "/" and path.startswith("/") and "//" not in path,
             f"{label} is not a canonical absolute POSIX path")
    parsed = PurePosixPath(path)
    _require(str(parsed) == path and all(
        part not in (".", "..") for part in parsed.parts
    ), f"{label} is not a canonical absolute POSIX path")
    return path


def exact_json_equal(left: Any, right: Any) -> bool:
    """Compare JSON values without Python bool/int or int/float coercions."""
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            exact_json_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            exact_json_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def strict_json_loads(data: bytes, label: str = "exact-main baseline JSON") -> Any:
    """Load one canonicalizable finite JSON value with duplicate-key rejection."""
    _require(type(data) is bytes, f"{label} is not a byte string")
    _require(0 < len(data) <= MAX_JSON_BYTES,
             f"{label} has an invalid byte length")

    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ExactMainBaselineError(
                    f"{label} contains duplicate object key {key!r}")
            result[key] = value
        return result

    def reject_constant(token: str) -> NoReturn:
        raise ExactMainBaselineError(
            f"{label} contains non-finite number {token}")

    def finite_float(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            raise ExactMainBaselineError(
                f"{label} contains a non-finite float")
        return value

    try:
        text = data.decode("utf-8")
        value = json.loads(
            text,
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
            parse_float=finite_float,
        )
        canonical_json_bytes(value)
        return value
    except ExactMainBaselineError:
        raise
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError,
            RecursionError) as error:
        raise ExactMainBaselineError(
            f"{label} is not strict singleton JSON") from error


def _require_exact_json_value(value: Any) -> None:
    value_type = type(value)
    if value_type is dict:
        for key, child in value.items():
            _require(type(key) is str,
                     "canonical JSON object key is not an exact string")
            _require_exact_json_value(child)
        return
    if value_type is list:
        for child in value:
            _require_exact_json_value(child)
        return
    if value_type in (str, int, bool, type(None)):
        return
    if value_type is float:
        _require(math.isfinite(value),
                 "canonical JSON contains a non-finite float")
        return
    _fail("canonical JSON contains a non-JSON value type")


def canonical_json_bytes(value: Any) -> bytes:
    """Encode finite JSON in the repository's canonical newline form."""
    try:
        _require_exact_json_value(value)
        return (json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ) + "\n").encode("utf-8")
    except (TypeError, ValueError, UnicodeError, RecursionError) as error:
        raise ExactMainBaselineError(
            "value is not canonical finite JSON") from error


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _selection_rule() -> dict[str, Any]:
    """Return a detached policy so callers cannot mutate validation semantics."""
    return {
        "algorithm": NORMALIZED_CODE_ALGORITHM,
        "required_flag": "SHF_ALLOC",
        "excluded_names": list(EXCLUDED_SECTION_NAMES),
        "excluded_name_prefixes": list(EXCLUDED_SECTION_PREFIXES),
        "order": "relative-elf-section-header-order-index-not-digested",
        "nobits": "structural-only-content-sha256-null",
    }


def _artifact_identity(value: Any, label: str) -> dict[str, Any]:
    artifact = _require_dict(value, ARTIFACT_KEYS, label)
    result = {
        "size": _bounded_int(
            artifact["size"], 64, MAX_ARTIFACT_BYTES, f"{label} size"),
        "sha256": _sha256(artifact["sha256"], f"{label} SHA-256"),
    }
    return result


def _section_record(value: Any, label: str) -> dict[str, Any]:
    section = _require_dict(value, SECTION_KEYS, label)
    index = _bounded_int(
        section["index"], 1, MAX_SECTION_INDEX, f"{label} index")
    name = _bounded_text(section["name"], f"{label} name")
    _require(_SECTION_NAME.fullmatch(name) is not None,
             f"{label} name is not a canonical ELF section name")
    section_type = _bounded_text(section["type"], f"{label} type")
    _require(section_type in ALLOC_SECTION_TYPES,
             f"{label} type is not an allowed allocatable ELF section type")
    flags = _bounded_int(
        section["flags"], 0, MAX_SECTION_FLAGS, f"{label} flags")
    _require(flags & SHF_ALLOC != 0, f"{label} is not SHF_ALLOC")
    _require(section_type != "NOBITS" or flags & SHF_WRITE != 0,
             f"{label} NOBITS section is not writable")
    _require(section_type != "NOBITS" or flags & SHF_EXECINSTR == 0,
             f"{label} NOBITS section is executable")
    _require(name not in EXCLUDED_SECTION_NAMES and not any(
        name.startswith(prefix) for prefix in EXCLUDED_SECTION_PREFIXES
    ), f"{label} is excluded by the normalized selection rule")
    address = _bounded_int(
        section["address"], 0, MAX_SECTION_ADDRESS, f"{label} address")
    size = _bounded_int(
        section["size"], 0, MAX_SECTION_SIZE, f"{label} size")
    alignment = _bounded_int(
        section["alignment"], 0, MAX_SECTION_ALIGNMENT,
        f"{label} alignment")
    _require(alignment in (0, 1) or alignment & (alignment - 1) == 0,
             f"{label} alignment is not a power of two")
    _require(alignment in (0, 1) or address % alignment == 0,
             f"{label} address violates its alignment")
    _require(address <= MAX_SECTION_ADDRESS - size,
             f"{label} address range overflows ELF64")
    content_sha256 = section["content_sha256"]
    if section_type == "NOBITS":
        _require(content_sha256 is None,
                 f"{label} NOBITS section has a content SHA-256")
    else:
        content_sha256 = _sha256(
            content_sha256, f"{label} content SHA-256")
        _require((size == 0) == (content_sha256 == EMPTY_CONTENT_SHA256),
                 f"{label} empty content SHA-256 is inconsistent")
    return {
        "index": index,
        "name": name,
        "type": section_type,
        "flags": flags,
        "address": address,
        "size": size,
        "alignment": alignment,
        "content_sha256": content_sha256,
    }


def _census_record(
    value: Any,
    sections: Sequence[Mapping[str, Any]],
    section_names: Sequence[str],
) -> dict[str, Any]:
    census = _require_dict(value, CENSUS_KEYS, "path-string census")
    _require(census["match_rule"] == PATH_CENSUS_MATCH_RULE and
             type(census["match_rule"]) is str,
             "path-string census match rule changed")
    roots_value = census["roots"]
    _require(type(roots_value) is list and
             len(roots_value) == len(PATH_CENSUS_ROOT_ROLES),
             "path-string census has the wrong root count")
    roots: list[dict[str, Any]] = []
    paths: list[str] = []
    for root_index, (root_value, expected_role) in enumerate(
            zip(roots_value, PATH_CENSUS_ROOT_ROLES)):
        root_label = f"path-string census root {root_index}"
        root = _require_dict(root_value, CENSUS_ROOT_KEYS, root_label)
        _require(root["role"] == expected_role and
                 type(root["role"]) is str,
                 f"{root_label} role changed")
        path = _absolute_posix_path(root["path"], f"{root_label} path")
        paths.append(path)
        rows_value = root["sections"]
        _require(type(rows_value) is list and
                 len(rows_value) == len(section_names),
                 f"{root_label} has the wrong section count")
        rows: list[dict[str, Any]] = []
        for section_index, (row_value, section_name, section) in enumerate(
                zip(rows_value, section_names, sections)):
            label = f"{root_label} section {section_index}"
            row = _require_dict(row_value, CENSUS_SECTION_KEYS, label)
            _require(row["name"] == section_name and
                     type(row["name"]) is str,
                     f"{label} does not align with the normalized section")
            occurrences = _bounded_int(
                row["occurrences"], 0, MAX_PATH_OCCURRENCES,
                f"{label} occurrence count")
            if section["type"] == "NOBITS":
                _require(occurrences == 0,
                         "NOBITS section claims file-backed path-string bytes")
            else:
                maximum_occurrences = (
                    section["size"] // len(path.encode("utf-8")))
                _require(occurrences <= maximum_occurrences,
                         f"{label} occurrence count exceeds section bytes")
            rows.append({"name": section_name, "occurrences": occurrences})
        total = _bounded_int(
            root["occurrences"], 0, MAX_PATH_OCCURRENCES,
            f"{root_label} total")
        _require(total == sum(row["occurrences"] for row in rows),
                 f"{root_label} total is inconsistent")
        roots.append({
            "role": expected_role,
            "path": path,
            "occurrences": total,
            "sections": rows,
        })
    _require(all(
        left not in right and right not in left
        for index, left in enumerate(paths)
        for right in paths[index + 1:]
    ), "path-string census roots overlap by exact byte substring")
    return {
        "match_rule": PATH_CENSUS_MATCH_RULE,
        "roots": roots,
    }


def _normalized_section_projection(
    section: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        key: copy.deepcopy(value)
        for key, value in section.items()
        if key != "index"
    }


def _combined_identity_sha256(sections: Sequence[Mapping[str, Any]]) -> str:
    projection = {
        "schema": NORMALIZED_CODE_IDENTITY_SCHEMA,
        "elf_class": ELF_CLASS,
        "endianness": ELF_ENDIANNESS,
        "machine": ELF_MACHINE,
        "selection_rule": _selection_rule(),
        "sections": [
            _normalized_section_projection(section) for section in sections
        ],
    }
    return canonical_sha256(projection)


def _record_sha256(value: Mapping[str, Any]) -> str:
    projection = {
        key: copy.deepcopy(field)
        for key, field in value.items()
        if key != "record_sha256"
    }
    return canonical_sha256(projection)


def validate_normalized_code_identity(value: Any) -> dict[str, Any]:
    """Validate and detach one normalized exact-main ELF identity."""
    identity = _require_dict(
        value, IDENTITY_KEYS, "normalized code identity")
    _require(identity["schema"] == NORMALIZED_CODE_IDENTITY_SCHEMA and
             type(identity["schema"]) is str,
             "normalized code identity schema changed")
    artifact = _artifact_identity(
        identity["artifact"], "normalized code identity artifact")
    _require(identity["elf_class"] == ELF_CLASS and
             type(identity["elf_class"]) is str,
             "normalized code identity ELF class changed")
    _require(identity["endianness"] == ELF_ENDIANNESS and
             type(identity["endianness"]) is str,
             "normalized code identity endianness changed")
    _require(identity["machine"] == ELF_MACHINE and
             type(identity["machine"]) is str,
             "normalized code identity machine changed")
    selection_rule = _selection_rule()
    _require_dict(
        identity["selection_rule"], SELECTION_RULE_KEYS,
        "normalized code identity selection rule")
    _require(exact_json_equal(identity["selection_rule"], selection_rule),
             "normalized code identity selection rule changed")

    sections_value = identity["sections"]
    _require(type(sections_value) is list and
             0 < len(sections_value) <= MAX_SECTION_COUNT,
             "normalized code identity has an invalid section count")
    sections = [
        _section_record(section, f"normalized section {index}")
        for index, section in enumerate(sections_value)
    ]
    indices = [section["index"] for section in sections]
    names = [section["name"] for section in sections]
    _require(indices == sorted(indices) and len(set(indices)) == len(indices),
             "normalized sections are not in unique header-index order")
    _require(len(set(names)) == len(names),
             "normalized sections contain a duplicate name")
    text_sections = [
        section for section in sections if section["name"] == ".text"
    ]
    _require(len(text_sections) == 1 and
             text_sections[0]["flags"] & SHF_EXECINSTR != 0,
             "normalized identity lacks one executable .text section")
    file_backed_size = sum(
        section["size"] for section in sections
        if section["type"] != "NOBITS"
    )
    _require(file_backed_size <= artifact["size"],
             "normalized file-backed sections exceed the raw artifact size")

    census = _census_record(
        identity["path_string_census"], sections, names)
    combined = _sha256(
        identity["combined_sha256"],
        "normalized code identity combined SHA-256")
    _require(combined == _combined_identity_sha256(sections),
             "normalized code identity combined SHA-256 is inconsistent")
    canonical = {
        "schema": NORMALIZED_CODE_IDENTITY_SCHEMA,
        "artifact": artifact,
        "elf_class": ELF_CLASS,
        "endianness": ELF_ENDIANNESS,
        "machine": ELF_MACHINE,
        "selection_rule": selection_rule,
        "sections": sections,
        "path_string_census": census,
        "combined_sha256": combined,
        "record_sha256": _sha256(
            identity["record_sha256"],
            "normalized code identity record SHA-256"),
    }
    _require(canonical["record_sha256"] == _record_sha256(canonical),
             "normalized code identity record SHA-256 is inconsistent")
    return copy.deepcopy(canonical)


def normalized_code_identity_record(
    *,
    artifact: dict[str, Any],
    sections: list[dict[str, Any]] | tuple[dict[str, Any], ...],
    path_string_census: dict[str, Any],
) -> dict[str, Any]:
    """Construct one closed normalized identity from retained section facts."""
    _require(type(artifact) is dict,
             "normalized artifact is not an exact JSON object")
    _require(type(sections) in (list, tuple),
             "normalized sections are not a sequence")
    _require(all(type(section) is dict for section in sections),
             "normalized sections contain a non-object entry")
    _require(type(path_string_census) is dict,
             "path-string census is not an exact JSON object")
    try:
        detached_artifact = copy.deepcopy(artifact)
        detached_sections = copy.deepcopy(list(sections))
        detached_census = copy.deepcopy(path_string_census)
    except (TypeError, ValueError, RecursionError) as error:
        raise ExactMainBaselineError(
            "normalized identity input is not bounded JSON data") from error
    value = {
        "schema": NORMALIZED_CODE_IDENTITY_SCHEMA,
        "artifact": detached_artifact,
        "elf_class": ELF_CLASS,
        "endianness": ELF_ENDIANNESS,
        "machine": ELF_MACHINE,
        "selection_rule": _selection_rule(),
        "sections": detached_sections,
        "path_string_census": detached_census,
        "combined_sha256": _combined_identity_sha256(detached_sections),
        "record_sha256": "0" * 64,
    }
    value["record_sha256"] = _record_sha256(value)
    return validate_normalized_code_identity(value)


def load_normalized_code_identity(data: bytes) -> dict[str, Any]:
    """Strictly load and validate one canonicalizable identity record."""
    return validate_normalized_code_identity(
        strict_json_loads(data, "normalized exact-main identity JSON"))


def _raw_elf64_section_header(
    data: bytes,
    table_offset: int,
    index: int,
) -> dict[str, int]:
    offset = table_offset + index * ELF64_SECTION_HEADER_SIZE
    fields = struct.unpack_from(ELF64_SECTION_HEADER_FORMAT, data, offset)
    return {
        "name_offset": fields[0],
        "type": fields[1],
        "flags": fields[2],
        "address": fields[3],
        "offset": fields[4],
        "size": fields[5],
        "link": fields[6],
        "info": fields[7],
        "alignment": fields[8],
        "entsize": fields[9],
    }


def parse_elf64_section_table(data: bytes) -> list[dict[str, Any]]:
    """Parse one complete little-endian x86-64 ELF section table."""
    _require(
        struct.calcsize(ELF64_HEADER_FORMAT) == ELF64_HEADER_SIZE and
        struct.calcsize(ELF64_SECTION_HEADER_FORMAT) ==
        ELF64_SECTION_HEADER_SIZE,
        "ELF struct layout diverges from the fixed ELF64 sizes",
    )
    _require(type(data) is bytes, "ELF artifact is not an exact byte string")
    _require(ELF64_HEADER_SIZE <= len(data) <= MAX_ELF_INPUT_BYTES,
             "ELF artifact has an invalid byte length")
    fields = struct.unpack_from(ELF64_HEADER_FORMAT, data, 0)
    ident = fields[0]
    _require(ident[:4] == b"\x7fELF", "ELF artifact has invalid magic")
    _require(ident[4] == 2, "ELF artifact is not ELF64")
    _require(ident[5] == 1, "ELF artifact is not little-endian")
    _require(ident[6] == 1, "ELF artifact has an invalid ident version")
    elf_type = fields[1]
    machine = fields[2]
    version = fields[3]
    section_table_offset = fields[6]
    elf_header_size = fields[8]
    section_header_size = fields[11]
    encoded_section_count = fields[12]
    encoded_string_table_index = fields[13]
    _require(elf_type in (ET_EXEC, ET_DYN),
             "ELF artifact is not an executable image")
    _require(machine == EM_X86_64, "ELF artifact is not x86-64")
    _require(version == 1, "ELF artifact has an invalid object version")
    _require(elf_header_size == ELF64_HEADER_SIZE,
             "ELF artifact has an invalid header size")
    _require(section_header_size == ELF64_SECTION_HEADER_SIZE,
             "ELF artifact has an invalid section-header size")
    _require(section_table_offset != 0 and
             section_table_offset <= len(data) - ELF64_SECTION_HEADER_SIZE,
             "ELF section table starts outside the artifact")

    section_zero = _raw_elf64_section_header(
        data, section_table_offset, 0)
    _require(
        section_zero["name_offset"] == 0 and
        section_zero["type"] == SHT_NULL and
        section_zero["flags"] == 0 and
        section_zero["address"] == 0 and
        section_zero["offset"] == 0,
        "ELF section zero is not the null section",
    )
    if encoded_section_count == 0:
        section_count = section_zero["size"]
        _require(section_count >= SHN_LORESERVE,
                 "ELF extended section count is not canonical")
    else:
        _require(encoded_section_count < SHN_LORESERVE,
                 "ELF direct section count uses a reserved encoding")
        _require(section_zero["size"] == 0,
                 "ELF section zero has an unexpected extended count")
        section_count = encoded_section_count
    if encoded_string_table_index == SHN_XINDEX:
        string_table_index = section_zero["link"]
        _require(string_table_index >= SHN_LORESERVE,
                 "ELF extended string-table index is not canonical")
    else:
        _require(encoded_string_table_index < SHN_LORESERVE,
                 "ELF direct string-table index uses a reserved encoding")
        _require(section_zero["link"] == 0,
                 "ELF section zero has an unexpected extended string index")
        string_table_index = encoded_string_table_index
    _require(type(section_count) is int and
             1 <= section_count <= MAX_ELF_SECTION_TABLE_COUNT,
             "ELF artifact has an invalid section count")
    table_size = section_count * ELF64_SECTION_HEADER_SIZE
    _require(section_table_offset <= len(data) - table_size,
             "ELF section table extends outside the artifact")
    _require(type(string_table_index) is int and
             0 < string_table_index < section_count,
             "ELF section-name string table index is invalid")

    raw_sections = [
        _raw_elf64_section_header(data, section_table_offset, index)
        for index in range(section_count)
    ]
    for index, section in enumerate(raw_sections):
        _require(section["size"] <= MAX_SECTION_SIZE,
                 f"ELF section {index} size exceeds the structural bound")
        _require(section["alignment"] <= MAX_SECTION_ALIGNMENT,
                 f"ELF section {index} alignment exceeds the structural bound")
        if section["type"] != SHT_NOBITS:
            _require(section["offset"] <= len(data) and
                     section["size"] <= len(data) - section["offset"],
                     f"ELF section {index} extends outside the artifact")

    string_table = raw_sections[string_table_index]
    _require(string_table["type"] == SHT_STRTAB and
             string_table["size"] > 0,
             "ELF section-name table is not a file-backed string table")
    string_table_bytes = data[
        string_table["offset"]:
        string_table["offset"] + string_table["size"]
    ]
    _require(string_table_bytes.startswith(b"\0"),
             "ELF section-name table lacks the empty name")

    sections: list[dict[str, Any]] = []
    for index, section in enumerate(raw_sections):
        name_offset = section["name_offset"]
        _require(name_offset < len(string_table_bytes),
                 f"ELF section {index} name offset is outside the string table")
        terminator = string_table_bytes.find(
            b"\0", name_offset, name_offset + MAX_TEXT_LENGTH + 1)
        _require(terminator >= 0,
                 f"ELF section {index} name is unterminated or too long")
        try:
            name = string_table_bytes[name_offset:terminator].decode("ascii")
        except UnicodeDecodeError as error:
            raise ExactMainBaselineError(
                f"ELF section {index} name is not ASCII") from error
        sections.append({
            "index": index,
            "name": name,
            "type": section["type"],
            "flags": section["flags"],
            "address": section["address"],
            "offset": section["offset"],
            "size": section["size"],
            "alignment": section["alignment"],
            "entsize": section["entsize"],
            "link": section["link"],
            "info": section["info"],
        })
    _require(all(set(section) == ELF_SECTION_KEYS for section in sections),
             "ELF section parser emitted an unexpected key set")
    return sections


def _normalized_sections_from_elf(
    data: bytes,
    sections: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[Mapping[str, Any]]]:
    _require(len(SECTION_TYPE_NAMES) == len(ALLOC_SECTION_TYPES) and
             set(SECTION_TYPE_NAMES.values()) == ALLOC_SECTION_TYPES,
             "ELF section-type map diverges from the normalized domain")
    selected: list[tuple[Mapping[str, Any], str]] = []
    for row in sections:
        if int(row["flags"]) & SHF_ALLOC == 0:
            continue
        name = str(row["name"])
        if name in EXCLUDED_SECTION_NAMES or any(
                name.startswith(prefix)
                for prefix in EXCLUDED_SECTION_PREFIXES):
            continue
        section_type = SECTION_TYPE_NAMES.get(int(row["type"]))
        _require(section_type is not None,
                 f"allocatable ELF section has an unknown type: {name!r}")
        selected.append((row, section_type))
        _require(len(selected) <= MAX_SECTION_COUNT,
                 "ELF artifact has an invalid retained section count")
    _require(selected,
             "ELF artifact has an invalid retained section count")

    file_ranges = sorted(
        (int(row["offset"]),
         int(row["offset"]) + int(row["size"]), str(row["name"]))
        for row, section_type in selected
        if section_type != "NOBITS" and int(row["size"]) != 0
    )
    _require(all(
        previous[1] <= current[0]
        for previous, current in zip(file_ranges, file_ranges[1:])
    ), "retained file-backed ELF sections overlap")

    normalized: list[dict[str, Any]] = []
    retained_rows: list[Mapping[str, Any]] = []
    for row, section_type in selected:
        offset = int(row["offset"])
        size = int(row["size"])
        if section_type == "NOBITS":
            content_sha256: str | None = None
        else:
            content_sha256 = hashlib.sha256(
                memoryview(data)[offset:offset + size]).hexdigest()
        normalized.append({
            "index": int(row["index"]),
            "name": str(row["name"]),
            "type": section_type,
            "flags": int(row["flags"]),
            "address": int(row["address"]),
            "size": size,
            "alignment": int(row["alignment"]),
            "content_sha256": content_sha256,
        })
        retained_rows.append(row)
    return normalized, retained_rows


def _validated_census_roots(roots: Any) -> dict[str, str]:
    _require(type(roots) is dict and
             set(roots) == set(PATH_CENSUS_ROOT_ROLES),
             "ELF path census roots have an unexpected key set")
    canonical = {
        role: _absolute_posix_path(roots[role], f"ELF {role}")
        for role in PATH_CENSUS_ROOT_ROLES
    }
    paths = [canonical[role] for role in PATH_CENSUS_ROOT_ROLES]
    _require(all(
        left not in right and right not in left
        for index, left in enumerate(paths)
        for right in paths[index + 1:]
    ), "ELF path census roots overlap by exact byte substring")
    return canonical


def _nonoverlapping_occurrence_count(
    data: bytes,
    needle: bytes,
    start: int,
    size: int,
) -> int:
    end = start + size
    count = 0
    cursor = start
    while cursor <= end - len(needle):
        match = data.find(needle, cursor, end)
        if match < 0:
            break
        count += 1
        _require(count <= MAX_PATH_OCCURRENCES,
                 "ELF path occurrence count exceeds its structural bound")
        cursor = match + len(needle)
    return count


def _path_string_census_from_normalized(
    data: bytes,
    normalized: Sequence[Mapping[str, Any]],
    retained_rows: Sequence[Mapping[str, Any]],
    roots: Any,
) -> dict[str, Any]:
    canonical_roots = _validated_census_roots(roots)
    census_roots = []
    for role in PATH_CENSUS_ROOT_ROLES:
        path = canonical_roots[role]
        needle = path.encode("utf-8")
        rows = []
        for normalized_section, retained_row in zip(
                normalized, retained_rows):
            if normalized_section["type"] == "NOBITS":
                occurrences = 0
            else:
                occurrences = _nonoverlapping_occurrence_count(
                    data,
                    needle,
                    int(retained_row["offset"]),
                    int(retained_row["size"]),
                )
            rows.append({
                "name": normalized_section["name"],
                "occurrences": occurrences,
            })
        census_roots.append({
            "role": role,
            "path": path,
            "occurrences": sum(row["occurrences"] for row in rows),
            "sections": rows,
        })
    census = {
        "match_rule": PATH_CENSUS_MATCH_RULE,
        "roots": census_roots,
    }
    return _census_record(
        census, normalized,
        [section["name"] for section in normalized])


def _path_string_census_from_parsed_elf(
    data: bytes,
    sections: Sequence[Mapping[str, Any]],
    roots: Any,
) -> dict[str, Any]:
    normalized, retained_rows = _normalized_sections_from_elf(data, sections)
    return _path_string_census_from_normalized(
        data, normalized, retained_rows, roots)


def path_string_census_from_elf(
    data: bytes,
    sections: list[dict[str, Any]],
    roots: dict[str, str],
) -> dict[str, Any]:
    """Recompute the selected-section path census from exact ELF bytes."""
    parsed = parse_elf64_section_table(data)
    _require(type(sections) is list and exact_json_equal(sections, parsed),
             "ELF path census section table was not derived from the artifact")
    return _path_string_census_from_parsed_elf(data, parsed, roots)


def normalized_code_identity_from_elf_bytes(
    data: bytes,
    *,
    roots: dict[str, str],
) -> dict[str, Any]:
    """Derive a closed normalized identity directly from retained ELF bytes."""
    sections = parse_elf64_section_table(data)
    normalized, retained_rows = _normalized_sections_from_elf(data, sections)
    census = _path_string_census_from_normalized(
        data, normalized, retained_rows, roots)
    return normalized_code_identity_record(
        artifact={"size": len(data), "sha256": hashlib.sha256(data).hexdigest()},
        sections=normalized,
        path_string_census=census,
    )


def verify_normalized_code_identity_against_elf_bytes(
    data: bytes,
    identity: Any,
    *,
    roots: dict[str, str],
) -> dict[str, Any]:
    """Bind an identity to exact ELF bytes and acquisition-provided roots."""
    validated = validate_normalized_code_identity(identity)
    canonical_roots = _validated_census_roots(roots)
    recorded_roots = {
        root["role"]: root["path"]
        for root in validated["path_string_census"]["roots"]
    }
    _require(exact_json_equal(recorded_roots, canonical_roots),
             "normalized identity path roots do not match acquisition")
    recomputed = normalized_code_identity_from_elf_bytes(
        data, roots=canonical_roots)
    _require(exact_json_equal(recomputed, validated),
             "normalized identity does not match the retained ELF bytes")
    return recomputed


__all__ = (
    "ExactMainBaselineError",
    "MAX_ELF_INPUT_BYTES",
    "NORMALIZED_CODE_IDENTITY_SCHEMA",
    "SECTION_TYPE_NAMES",
    "canonical_json_bytes",
    "canonical_sha256",
    "exact_json_equal",
    "load_normalized_code_identity",
    "normalized_code_identity_from_elf_bytes",
    "normalized_code_identity_record",
    "parse_elf64_section_table",
    "path_string_census_from_elf",
    "strict_json_loads",
    "validate_normalized_code_identity",
    "verify_normalized_code_identity_against_elf_bytes",
)
