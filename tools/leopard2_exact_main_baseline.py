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
from pathlib import PurePosixPath
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
MAX_ARTIFACT_BYTES = (1 << 63) - 1
MAX_SECTION_COUNT = 4096
MAX_SECTION_INDEX = (1 << 31) - 1
MAX_SECTION_ADDRESS = (1 << 64) - 1
MAX_SECTION_SIZE = (1 << 63) - 1
MAX_SECTION_FLAGS = (1 << 64) - 1
MAX_SECTION_ALIGNMENT = 1 << 63
MAX_PATH_OCCURRENCES = (1 << 31) - 1
MAX_TEXT_LENGTH = 4096

# Fixed Unicode-15 format/control ranges.  Do not use unicodedata.category here:
# record validity must not drift with the Python runtime's Unicode database.
_DISPLAY_EXCLUDED_CODE_POINT_RANGES = (
    (0x00AD, 0x00AD), (0x0600, 0x0605), (0x061C, 0x061C),
    (0x06DD, 0x06DD), (0x070F, 0x070F), (0x0890, 0x0891),
    (0x08E2, 0x08E2), (0x180E, 0x180E), (0x200B, 0x200F),
    (0x2028, 0x202E), (0x2060, 0x2064), (0x2066, 0x206F),
    (0xFEFF, 0xFEFF), (0xFFF9, 0xFFFB), (0x110BD, 0x110BD),
    (0x110CD, 0x110CD), (0x13430, 0x1343F), (0x1BCA0, 0x1BCA3),
    (0x1D173, 0x1D17A), (0xE0001, 0xE0001), (0xE0020, 0xE007F),
)
DISPLAY_EXCLUDED_CODE_POINTS = frozenset(
    code_point
    for first, last in _DISPLAY_EXCLUDED_CODE_POINT_RANGES
    for code_point in range(first, last + 1)
)
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
        not (0xD800 <= ord(character) <= 0xDFFF) and
        ord(character) not in DISPLAY_EXCLUDED_CODE_POINTS
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


__all__ = (
    "ExactMainBaselineError",
    "NORMALIZED_CODE_IDENTITY_SCHEMA",
    "canonical_json_bytes",
    "canonical_sha256",
    "exact_json_equal",
    "load_normalized_code_identity",
    "normalized_code_identity_record",
    "strict_json_loads",
    "validate_normalized_code_identity",
)
