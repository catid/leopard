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
    Metadata,
    PROFILE_LEGACY_HIGH,
    PROFILE_LOW,
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
        for profile in (PROFILE_LEGACY_HIGH, PROFILE_LOW):
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
        self.assertEqual(tested, 3952)

    def test_exhaustive_practical_metadata_round_trips(self) -> None:
        tested = 0
        for count in range(MAX_METADATA + 1):
            metadata = tuple(
                Metadata(type_id, bytes(((type_id * 17 + i) & 0xFF) for i in range(type_id % 257)))
                for type_id in range(1, count + 1)
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
            (Metadata(2, b"b"), Metadata(1, b"a")),
            (Metadata(1, b"a"), Metadata(0x8001, bytes(32))),
            (Metadata(META_COORDINATE_SET_SHA256, bytes(31)),),
            (Metadata(1, b"x" * (MAX_METADATA_VALUE + 1)),),
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
