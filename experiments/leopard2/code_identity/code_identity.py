#!/usr/bin/env python3
"""Experimental Leopard2 wire-code identity serializer.

This is deliberately not part of the public API.  It is a small independent
reference used to freeze and test the Experiment W binary proposal.
"""

from __future__ import annotations

import dataclasses
import struct
from typing import Iterable, Tuple


MAGIC = b"L2ID"
FORMAT_VERSION = 1
HEADER_BYTES = 36
MAX_ID_BYTES = 65535
MAX_METADATA = 64
MAX_METADATA_VALUE = 4096

PROFILE_LEGACY_HIGH = 1
PROFILE_LOW = 2
# This is the persistent identity assigned to the public enum value
# LEO2_PROFILE_EXACT_EXPERIMENTAL_V1.  C7 freezes its v1 coordinate semantics
# as exact-prefix low-rate; production codec construction remains unsupported.
PROFILE_EXACT_LOW = 3
# Family 4 is deliberately not accepted yet.  Reserving the number here keeps
# C8's future exact-high profile from reusing family 3 or silently changing the
# public enum-3 exact-low mathematics.
PROFILE_EXACT_HIGH_RESERVED = 4
FIELD_GF8 = 1
FIELD_GF16 = 2

# Critical metadata known by this reader.  The first four values are SHA-256
# digests of canonical mathematical objects; shard layout is a one-byte wire
# framing selector.  None of these values names an implementation choice.
META_PROFILE_PARAMETERS_SHA256 = 0x8001
META_COORDINATE_SET_SHA256 = 0x8002
META_SHORTENING_SET_SHA256 = 0x8003
META_PUNCTURING_SET_SHA256 = 0x8004
META_SHARD_LAYOUT = 0x8005
SHARD_LAYOUT_NATIVE_V1 = 0
SHARD_LAYOUT_GF16_PADDED_ODD_V1 = 1
_KNOWN_DIGESTS = {
    META_PROFILE_PARAMETERS_SHA256,
    META_COORDINATE_SET_SHA256,
    META_SHORTENING_SET_SHA256,
    META_PUNCTURING_SET_SHA256,
}
_KNOWN_CRITICAL = _KNOWN_DIGESTS | {META_SHARD_LAYOUT}


class IdentityError(ValueError):
    """The bytes are malformed, non-canonical, or semantically unsupported."""


@dataclasses.dataclass(frozen=True, order=True)
class Metadata:
    type: int
    value: bytes


@dataclasses.dataclass(frozen=True)
class CodeIdentity:
    profile: int
    profile_version: int
    field: int
    field_version: int
    original_count: int
    recovery_count: int
    parent_count: int
    padded_side: int
    coordinate_map_version: int
    metadata: Tuple[Metadata, ...] = ()


_HEADER = struct.Struct(">4sBBHI4B4I2H")
_TLV = struct.Struct(">HH")
if _HEADER.size != HEADER_BYTES:
    raise RuntimeError("code-identity header size constant is inconsistent")


def _ceil_pow2(value: int) -> int:
    if value <= 0:
        raise IdentityError("counts must be positive")
    return 1 << (value - 1).bit_length()


def _validate(identity: CodeIdentity) -> None:
    scalar_u8 = (
        ("profile", identity.profile),
        ("profile_version", identity.profile_version),
        ("field", identity.field),
        ("field_version", identity.field_version),
    )
    for name, value in scalar_u8:
        if not isinstance(value, int) or not 0 <= value <= 0xFF:
            raise IdentityError(f"{name} does not fit uint8")
    scalar_u32 = (
        ("original_count", identity.original_count),
        ("recovery_count", identity.recovery_count),
        ("parent_count", identity.parent_count),
        ("padded_side", identity.padded_side),
    )
    for name, value in scalar_u32:
        if not isinstance(value, int) or not 0 <= value <= 0xFFFFFFFF:
            raise IdentityError(f"{name} does not fit uint32")
    if not isinstance(identity.coordinate_map_version, int) or not (
        0 <= identity.coordinate_map_version <= 0xFFFF
    ):
        raise IdentityError("coordinate_map_version does not fit uint16")

    if identity.profile not in (
        PROFILE_LEGACY_HIGH, PROFILE_LOW, PROFILE_EXACT_LOW
    ):
        raise IdentityError("unknown profile family")
    if identity.profile_version != 1:
        raise IdentityError("unsupported profile version")
    if identity.field not in (FIELD_GF8, FIELD_GF16):
        raise IdentityError("unknown field")
    if identity.field_version != 1:
        raise IdentityError("unsupported field representation version")
    if identity.coordinate_map_version != 1:
        raise IdentityError("unsupported coordinate-map version")
    if identity.original_count == 0 or identity.recovery_count == 0:
        raise IdentityError("K and R must be positive")

    if identity.profile == PROFILE_LEGACY_HIGH:
        side = _ceil_pow2(identity.recovery_count)
        parent = _ceil_pow2(identity.original_count + side)
    elif identity.profile == PROFILE_LOW:
        side = _ceil_pow2(identity.original_count)
        parent = _ceil_pow2(side + identity.recovery_count)
    else:
        # Exact-low V1 has no dyadic parent or padded side.  These redundant
        # envelope fields carry its exact transmitted length and K-side size.
        side = identity.original_count
        parent = identity.original_count + identity.recovery_count
        if parent > 0xFFFFFFFF:
            raise IdentityError("exact code length overflows uint32")
    if identity.padded_side != side or identity.parent_count != parent:
        raise IdentityError("parent or padded side is inconsistent with profile")
    if parent > 65536:
        raise IdentityError("parent exceeds supported field order")
    if identity.field == FIELD_GF8 and parent > 256:
        raise IdentityError("GF8 parent exceeds 256 coordinates")

    if len(identity.metadata) > MAX_METADATA:
        raise IdentityError("too many metadata entries")
    previous = -1
    seen_base_types = set()
    for item in identity.metadata:
        if not isinstance(item, Metadata):
            raise IdentityError("metadata entries must be Metadata")
        if not 0 < item.type <= 0xFFFF:
            raise IdentityError("metadata type zero is reserved")
        if item.type <= previous:
            raise IdentityError("metadata is not in canonical type order")
        previous = item.type
        base_type = item.type & 0x7FFF
        if base_type in seen_base_types:
            raise IdentityError("duplicate metadata base type")
        seen_base_types.add(base_type)
        if not item.type & 0x8000 and (item.type | 0x8000) in _KNOWN_CRITICAL:
            raise IdentityError("noncritical alias of known critical metadata")
        if item.type & 0x8000 and item.type not in _KNOWN_CRITICAL:
            raise IdentityError("unknown critical metadata")
        if not isinstance(item.value, bytes):
            raise IdentityError("metadata values must be bytes")
        if len(item.value) > MAX_METADATA_VALUE:
            raise IdentityError("metadata value is too long")
        if item.type in _KNOWN_DIGESTS and len(item.value) != 32:
            raise IdentityError("known digest metadata must be 32 bytes")
        if identity.profile == PROFILE_EXACT_LOW and item.type in (
            META_COORDINATE_SET_SHA256,
            META_SHORTENING_SET_SHA256,
            META_PUNCTURING_SET_SHA256,
        ):
            raise IdentityError(
                "exact-low V1 map is fixed and has no coordinate/shortening/puncturing TLV"
            )
        if item.type == META_SHARD_LAYOUT:
            if len(item.value) != 1:
                raise IdentityError("shard-layout metadata must be one byte")
            layout = item.value[0]
            if layout == SHARD_LAYOUT_NATIVE_V1:
                raise IdentityError("explicit native shard layout is non-canonical")
            if layout != SHARD_LAYOUT_GF16_PADDED_ODD_V1:
                raise IdentityError("unsupported shard layout")
            if identity.field != FIELD_GF16:
                raise IdentityError("padded-odd shard layout requires GF16")


def serialize(identity: CodeIdentity) -> bytes:
    """Return the one canonical v1 representation of identity."""
    _validate(identity)
    total = HEADER_BYTES + sum(_TLV.size + len(item.value) for item in identity.metadata)
    if total > MAX_ID_BYTES:
        raise IdentityError("serialized identity is too long")
    result = bytearray(
        _HEADER.pack(
            MAGIC,
            FORMAT_VERSION,
            0,  # reserved flags: must stay zero in v1
            HEADER_BYTES,
            total,
            identity.profile,
            identity.profile_version,
            identity.field,
            identity.field_version,
            identity.original_count,
            identity.recovery_count,
            identity.parent_count,
            identity.padded_side,
            identity.coordinate_map_version,
            len(identity.metadata),
        )
    )
    for item in identity.metadata:
        result += _TLV.pack(item.type, len(item.value))
        result += item.value
    return bytes(result)


def deserialize(data: bytes) -> CodeIdentity:
    """Parse a canonical v1 identity, failing closed on semantic extensions."""
    if not isinstance(data, bytes):
        raise IdentityError("serialized identity must be bytes")
    if len(data) < HEADER_BYTES:
        raise IdentityError("truncated header")
    (
        magic,
        format_version,
        flags,
        header_bytes,
        total_bytes,
        profile,
        profile_version,
        field,
        field_version,
        original_count,
        recovery_count,
        parent_count,
        padded_side,
        coordinate_map_version,
        metadata_count,
    ) = _HEADER.unpack_from(data)
    if magic != MAGIC:
        raise IdentityError("bad magic")
    if format_version != FORMAT_VERSION:
        raise IdentityError("unsupported envelope version")
    if flags != 0:
        raise IdentityError("reserved flags are nonzero")
    if header_bytes != HEADER_BYTES:
        raise IdentityError("non-canonical header size")
    if total_bytes != len(data) or total_bytes > MAX_ID_BYTES:
        raise IdentityError("declared length does not match input")
    if metadata_count > MAX_METADATA:
        raise IdentityError("too many metadata entries")

    metadata = []
    offset = HEADER_BYTES
    for _ in range(metadata_count):
        if len(data) - offset < _TLV.size:
            raise IdentityError("truncated metadata header")
        type_id, length = _TLV.unpack_from(data, offset)
        offset += _TLV.size
        if length > MAX_METADATA_VALUE or length > len(data) - offset:
            raise IdentityError("truncated or oversized metadata value")
        metadata.append(Metadata(type_id, data[offset : offset + length]))
        offset += length
    if offset != len(data):
        raise IdentityError("unaccounted trailing bytes")

    identity = CodeIdentity(
        profile,
        profile_version,
        field,
        field_version,
        original_count,
        recovery_count,
        parent_count,
        padded_side,
        coordinate_map_version,
        tuple(metadata),
    )
    _validate(identity)
    if serialize(identity) != data:
        raise IdentityError("non-canonical representation")
    return identity


def make_identity(
    profile: int,
    field: int,
    original_count: int,
    recovery_count: int,
    metadata: Iterable[Metadata] = (),
) -> CodeIdentity:
    """Construct a current-profile identity and derive its redundant fields."""
    if profile == PROFILE_LEGACY_HIGH:
        side = _ceil_pow2(recovery_count)
        parent = _ceil_pow2(original_count + side)
    elif profile == PROFILE_LOW:
        side = _ceil_pow2(original_count)
        parent = _ceil_pow2(side + recovery_count)
    elif profile == PROFILE_EXACT_LOW:
        side = original_count
        parent = original_count + recovery_count
    else:
        raise IdentityError("unknown profile family")
    return CodeIdentity(
        profile=profile,
        profile_version=1,
        field=field,
        field_version=1,
        original_count=original_count,
        recovery_count=recovery_count,
        parent_count=parent,
        padded_side=side,
        coordinate_map_version=1,
        metadata=tuple(sorted(metadata)),
    )
