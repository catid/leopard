#!/usr/bin/env python3
"""Generate a deterministic libFuzzer corpus for Experiment W."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import struct
import sys

sys.dont_write_bytecode = True

from code_identity import (  # noqa: E402
    FIELD_GF16,
    FIELD_GF8,
    PROFILE_EXACT_LOW,
    PROFILE_LEGACY_HIGH,
    PROFILE_LOW,
    serialize,
    make_identity,
)
from test_code_identity_c import malformed_corpus  # noqa: E402


HERE = pathlib.Path(__file__).resolve().parent


def builder_seed(
    profile: int,
    field: int,
    original_count: int,
    recovery_count: int,
    metadata: tuple[tuple[int, bytes], ...] = (),
) -> bytes:
    result = bytearray((profile, field))
    result += struct.pack(">II", original_count, recovery_count)
    result.append(len(metadata))
    for type_id, value in metadata:
        if len(value) > 255:
            raise ValueError("builder seed metadata exceeds one-byte length")
        result += struct.pack(">H", type_id)
        result.append(len(value))
        result += value
    return bytes(result)


def corpus() -> tuple[bytes, ...]:
    seeds: set[bytes] = set()
    vectors = json.loads((HERE / "golden_vectors.json").read_text())
    seeds.update(bytes.fromhex(vector["serialized_hex"]) for vector in vectors)
    seeds.update(malformed_corpus())

    for profile in (PROFILE_LEGACY_HIGH, PROFILE_LOW, PROFILE_EXACT_LOW):
        for field in (FIELD_GF8, FIELD_GF16):
            for original_count, recovery_count in (
                (1, 1), (3, 5), (129, 1), (129, 100), (240, 16),
                (255, 1), (256, 256), (4096, 512),
            ):
                seeds.add(builder_seed(
                    profile, field, original_count, recovery_count,
                    ((0x1234, b"future"),),
                ))
                try:
                    seeds.add(serialize(make_identity(
                        profile, field, original_count, recovery_count
                    )))
                except (TypeError, ValueError):
                    pass

    # Exercise builder parsing for critical, duplicate-base, zero-length, and
    # maximum one-byte value lengths even when the complete record is invalid.
    seeds.add(builder_seed(
        PROFILE_LOW, FIELD_GF16, 129, 100,
        ((0x8001, bytes(32)), (0x0001, b"downgrade")),
    ))
    seeds.add(builder_seed(
        PROFILE_LOW, FIELD_GF16, 2, 256,
        ((0x8005, b"\x01"), (0x1234, bytes(range(255)))),
    ))
    seeds.add(builder_seed(
        PROFILE_EXACT_LOW, FIELD_GF8, 3, 253,
        ((0x8002, bytes(32)), (0x1234, b"")),
    ))
    return tuple(sorted(seeds, key=lambda value: (len(value), value)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=pathlib.Path)
    arguments = parser.parse_args()

    arguments.output.mkdir(parents=True, exist_ok=True)
    written: list[dict[str, object]] = []
    for value in corpus():
        digest = hashlib.sha256(value).hexdigest()
        path = arguments.output / f"{digest}.bin"
        path.write_bytes(value)
        written.append({"sha256": digest, "bytes": len(value)})
    print(json.dumps({
        "schema": "leopard2-code-identity-fuzz-corpus-v1",
        "count": len(written),
        "entries": written,
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
