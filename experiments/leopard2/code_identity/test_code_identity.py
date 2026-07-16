#!/usr/bin/env python3
"""Dependency-free validation for the Experiment W reference format."""

from __future__ import annotations

import hashlib
import json
import pathlib
import random
import struct
import unittest

from code_identity import (
    CodeIdentity,
    FIELD_GF16,
    FIELD_GF8,
    HEADER_BYTES,
    IdentityError,
    MAX_METADATA,
    MAX_METADATA_VALUE,
    META_COORDINATE_SET_SHA256,
    META_SHARD_LAYOUT,
    Metadata,
    META_PUNCTURING_SET_SHA256,
    META_SHORTENING_SET_SHA256,
    PROFILE_EXACT_LOW,
    PROFILE_EXACT_HIGH_RESERVED,
    PROFILE_LEGACY_HIGH,
    PROFILE_LOW,
    SHARD_LAYOUT_GF16_PADDED_ODD_V1,
    SHARD_LAYOUT_NATIVE_V1,
    deserialize,
    make_identity,
    serialize,
)


HERE = pathlib.Path(__file__).resolve().parent


class CodeIdentityTests(unittest.TestCase):
    def assert_rejected(self, data: bytes) -> None:
        with self.assertRaises(IdentityError):
            deserialize(data)

    def test_golden_vectors(self) -> None:
        vectors = json.loads((HERE / "golden_vectors.json").read_text())
        for vector in vectors:
            identity = make_identity(
                vector["profile"],
                vector["field"],
                vector["K"],
                vector["R"],
                tuple(
                    Metadata(int(item["type"], 0), bytes.fromhex(item["hex"]))
                    for item in vector["metadata"]
                ),
            )
            encoded = bytes.fromhex(vector["serialized_hex"])
            self.assertEqual(serialize(identity), encoded, vector["name"])
            self.assertEqual(deserialize(encoded), identity, vector["name"])
            self.assertEqual(hashlib.sha256(encoded).hexdigest(), vector["sha256"])

    def test_profile_boundary_round_trips(self) -> None:
        counts = (1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63,
                  64, 65, 127, 128, 129, 191, 192, 193, 223, 224, 225,
                  239, 240, 241, 247, 248, 249, 255, 256, 257, 511, 512,
                  513, 1000, 4096)
        tested = 0
        for profile in (PROFILE_LEGACY_HIGH, PROFILE_LOW, PROFILE_EXACT_LOW):
            for k in counts:
                for r in counts:
                    for field in (FIELD_GF8, FIELD_GF16):
                        try:
                            identity = make_identity(profile, field, k, r)
                            encoded = serialize(identity)
                        except IdentityError:
                            continue
                        self.assertEqual(deserialize(encoded), identity)
                        tested += 1
        self.assertEqual(tested, 6085)

    def test_exact_low_v1_semantics_and_enum_reconciliation(self) -> None:
        self.assertEqual(PROFILE_EXACT_LOW, 3)
        self.assertEqual(PROFILE_EXACT_HIGH_RESERVED, 4)
        gf8 = make_identity(PROFILE_EXACT_LOW, FIELD_GF8, 3, 253)
        self.assertEqual(gf8.parent_count, 256)
        self.assertEqual(gf8.padded_side, 3)
        self.assertEqual(gf8.coordinate_map_version, 1)
        self.assertEqual(deserialize(serialize(gf8)), gf8)

        gf16 = make_identity(PROFILE_EXACT_LOW, FIELD_GF16, 1000, 4096)
        self.assertEqual(gf16.parent_count, 5096)
        self.assertEqual(gf16.padded_side, 1000)
        self.assertEqual(deserialize(serialize(gf16)), gf16)

        with self.assertRaises(IdentityError):
            serialize(make_identity(PROFILE_EXACT_LOW, FIELD_GF8, 4, 253))
        with self.assertRaises(IdentityError):
            make_identity(PROFILE_EXACT_HIGH_RESERVED, FIELD_GF16, 253, 3)
        for metadata_type in (
            META_COORDINATE_SET_SHA256,
            META_SHORTENING_SET_SHA256,
            META_PUNCTURING_SET_SHA256,
        ):
            for digest in (bytes(32), bytes((0xff,)) * 32):
                with self.assertRaises(IdentityError):
                    serialize(make_identity(
                        PROFILE_EXACT_LOW,
                        FIELD_GF16,
                        129,
                        1000,
                        (Metadata(metadata_type, digest),),
                    ))

    def test_exhaustive_practical_metadata_round_trips(self) -> None:
        tested = 0
        for count in range(MAX_METADATA + 1):
            metadata = tuple(
                Metadata(type_id, bytes(((type_id * 17 + i) & 0xFF) for i in range(type_id % 257)))
                for type_id in range(6, 6 + count)
            )
            identity = make_identity(PROFILE_LOW, FIELD_GF16, 129, 100, metadata)
            self.assertEqual(deserialize(serialize(identity)), identity)
            tested += 1
        for length in range(0, MAX_METADATA_VALUE + 1):
            item = Metadata(0x1234, bytes((i * 29 + length) & 0xFF for i in range(length)))
            identity = make_identity(PROFILE_LEGACY_HIGH, FIELD_GF8, 240, 16, (item,))
            self.assertEqual(deserialize(serialize(identity)), identity)
            tested += 1
        self.assertEqual(tested, 4162)

    def test_unknown_optional_round_trips_and_critical_fails_closed(self) -> None:
        optional = make_identity(
            PROFILE_LOW, FIELD_GF16, 129, 100, (Metadata(0x1234, b"future"),)
        )
        self.assertEqual(deserialize(serialize(optional)), optional)
        with self.assertRaises(IdentityError):
            serialize(make_identity(
                PROFILE_LOW, FIELD_GF16, 129, 100, (Metadata(0x9234, b"future"),)
            ))

    def test_known_critical_digest(self) -> None:
        digest = hashlib.sha256(b"canonical-coordinate-set").digest()
        identity = make_identity(
            PROFILE_LOW,
            FIELD_GF16,
            129,
            100,
            (Metadata(META_COORDINATE_SET_SHA256, digest),),
        )
        self.assertEqual(deserialize(serialize(identity)), identity)

    def test_padded_odd_shard_layout_identity(self) -> None:
        native = make_identity(PROFILE_LOW, FIELD_GF16, 2, 256)
        native_bytes = serialize(native)
        self.assertEqual(len(native_bytes), HEADER_BYTES)
        self.assertEqual(deserialize(native_bytes), native)

        padded = make_identity(
            PROFILE_LOW,
            FIELD_GF16,
            2,
            256,
            (Metadata(
                META_SHARD_LAYOUT,
                bytes((SHARD_LAYOUT_GF16_PADDED_ODD_V1,)),
            ),),
        )
        encoded = serialize(padded)
        self.assertEqual(encoded[:8], native_bytes[:8])
        self.assertEqual(encoded[8:12], struct.pack(">I", HEADER_BYTES + 5))
        self.assertEqual(encoded[12:34], native_bytes[12:34])
        self.assertEqual(encoded[34:36], b"\x00\x01")
        self.assertEqual(encoded[HEADER_BYTES:], b"\x80\x05\x00\x01\x01")
        self.assertEqual(deserialize(encoded), padded)

        rejected = (
            Metadata(META_SHARD_LAYOUT, b""),
            Metadata(META_SHARD_LAYOUT, bytes((SHARD_LAYOUT_NATIVE_V1,))),
            Metadata(META_SHARD_LAYOUT, b"\x02"),
            Metadata(META_SHARD_LAYOUT, b"\x01\x00"),
        )
        for metadata in rejected:
            with self.assertRaises(IdentityError):
                serialize(make_identity(
                    PROFILE_LOW, FIELD_GF16, 2, 256, (metadata,)
                ))
        for base_type in range(1, 6):
            with self.assertRaises(IdentityError):
                serialize(make_identity(
                    PROFILE_LOW, FIELD_GF16, 2, 256,
                    (Metadata(base_type, b"downgrade"),)
                ))
        with self.assertRaises(IdentityError):
            serialize(make_identity(
                PROFILE_LOW,
                FIELD_GF8,
                2,
                5,
                (Metadata(META_SHARD_LAYOUT, b"\x01"),),
            ))
        with self.assertRaises(IdentityError):
            serialize(CodeIdentity(**{
                **padded.__dict__,
                "metadata": (
                    Metadata(META_SHARD_LAYOUT & 0x7FFF, b"duplicate"),
                    Metadata(META_SHARD_LAYOUT, b"\x01"),
                ),
            }))

    def test_every_truncation_and_trailing_byte_rejected(self) -> None:
        identity = make_identity(
            PROFILE_LOW, FIELD_GF16, 129, 100, (Metadata(0x1234, b"future"),)
        )
        encoded = serialize(identity)
        for length in range(len(encoded)):
            self.assert_rejected(encoded[:length])
        self.assert_rejected(encoded + b"\x00")

    def test_malformed_fixed_header(self) -> None:
        base = bytearray(serialize(make_identity(PROFILE_LEGACY_HIGH, FIELD_GF8, 240, 16)))
        mutations = (
            (0, 0),        # magic
            (4, 2),        # future envelope version
            (5, 1),        # reserved flags
            (7, 35),       # non-canonical header length
            (12, 99),      # profile
            (13, 2),       # profile version
            (14, 99),      # field
            (15, 2),       # field representation version
            (31, 15),      # padded side disagrees
            (33, 2),       # coordinate-map version
            (35, 65),      # excessive metadata count
        )
        for offset, value in mutations:
            candidate = bytearray(base)
            candidate[offset] = value
            self.assert_rejected(bytes(candidate))
        for declared in (0, HEADER_BYTES - 1, HEADER_BYTES + 1, 0xFFFFFFFF):
            candidate = bytearray(base)
            struct.pack_into(">I", candidate, 8, declared)
            self.assert_rejected(bytes(candidate))

    def test_noncanonical_and_invalid_metadata(self) -> None:
        base = make_identity(PROFILE_LOW, FIELD_GF16, 129, 100)
        for metadata in (
            (Metadata(0, b""),),
            (Metadata(7, b"b"), Metadata(6, b"a")),
            (Metadata(6, b"a"), Metadata(0x8006, bytes(32))),
            (Metadata(META_COORDINATE_SET_SHA256, bytes(31)),),
            (Metadata(6, b"x" * (MAX_METADATA_VALUE + 1)),),
        ):
            with self.assertRaises(IdentityError):
                serialize(CodeIdentity(**{**base.__dict__, "metadata": metadata}))

    def test_overflow_and_semantic_inconsistency(self) -> None:
        base = make_identity(PROFILE_LOW, FIELD_GF16, 129, 100)
        replacements = (
            {"original_count": 0},
            {"recovery_count": 0},
            {"original_count": 0x100000000},
            {"parent_count": 255},
            {"padded_side": 7},
            {"field": FIELD_GF8},
            {"coordinate_map_version": 0x10000},
        )
        for replacement in replacements:
            with self.assertRaises(IdentityError):
                serialize(dataclasses_replace(base, **replacement))

    def test_deterministic_invalid_fuzz(self) -> None:
        rng = random.Random(0x4C324944)
        valid = serialize(make_identity(
            PROFILE_LOW, FIELD_GF16, 129, 100, (Metadata(0x1234, b"future"),)
        ))
        rejected = 0
        accepted_canonical = 0
        for _ in range(20000):
            candidate = bytearray(valid)
            for _ in range(rng.randrange(1, 5)):
                candidate[rng.randrange(len(candidate))] ^= 1 << rng.randrange(8)
            try:
                decoded = deserialize(bytes(candidate))
                self.assertEqual(serialize(decoded), bytes(candidate))
                accepted_canonical += 1
            except IdentityError:
                rejected += 1
        self.assertEqual(rejected + accepted_canonical, 20000)
        self.assertGreater(rejected, 15000)


def dataclasses_replace(identity: CodeIdentity, **changes: int) -> CodeIdentity:
    values = dict(identity.__dict__)
    values.update(changes)
    return CodeIdentity(**values)


if __name__ == "__main__":
    unittest.main(verbosity=2)
