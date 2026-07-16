#!/usr/bin/env python3
"""Compile and differentially test the independent C99 code-identity oracle."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import random
import struct
import subprocess
import sys
import tempfile
from typing import Iterable, Sequence

sys.dont_write_bytecode = True

from code_identity import (  # noqa: E402
    FIELD_GF16,
    FIELD_GF8,
    HEADER_BYTES,
    IdentityError,
    MAX_METADATA,
    MAX_METADATA_VALUE,
    META_COORDINATE_SET_SHA256,
    META_SHARD_LAYOUT,
    Metadata,
    PROFILE_EXACT_LOW,
    PROFILE_LEGACY_HIGH,
    PROFILE_LOW,
    deserialize,
    make_identity,
    serialize,
)


HERE = pathlib.Path(__file__).resolve().parent
COUNTS = (
    1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
    127, 128, 129, 191, 192, 193, 223, 224, 225, 239, 240, 241,
    247, 248, 249, 255, 256, 257, 511, 512, 513, 1000, 4096,
)
COUNT_EDGE_PAIRS = (
    (0, 1),
    (1, 0),
    (32767, 32768),
    (32768, 32768),
    (32769, 32767),
    (65535, 1),
    (65536, 1),
    (65537, 1),
    (0xFFFFFFFF, 1),
    (1, 0xFFFFFFFF),
    (0xFFFFFFFF, 0xFFFFFFFF),
)


class COracle:
    def __init__(self, executable: pathlib.Path) -> None:
        self.process = subprocess.Popen(
            [str(executable), "--protocol"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="ascii",
        )

    def query(self, command: str) -> str:
        if self.process.stdin is None or self.process.stdout is None:
            raise RuntimeError("C oracle pipes are unavailable")
        self.process.stdin.write(command + "\n")
        self.process.stdin.flush()
        response = self.process.stdout.readline()
        if not response:
            stderr = ""
            if self.process.stderr is not None:
                stderr = self.process.stderr.read()
            raise RuntimeError(
                f"C oracle exited {self.process.poll()}: {stderr.strip()}"
            )
        return response.rstrip("\r\n")

    def close(self) -> None:
        if self.process.poll() is None and self.process.stdin is not None:
            self.process.stdin.write("Q\n")
            self.process.stdin.flush()
        stdout, stderr = self.process.communicate(timeout=10)
        if self.process.returncode != 0 or stdout or stderr:
            raise RuntimeError(
                f"C oracle shutdown failed: rc={self.process.returncode} "
                f"stdout={stdout!r} stderr={stderr!r}"
            )


def compile_oracle(cc: str, output: pathlib.Path, sanitizer: bool) -> None:
    command = [
        cc,
        "-std=c99",
        "-O1" if sanitizer else "-O2",
        "-g" if sanitizer else "-DNDEBUG",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wconversion",
        "-Wsign-conversion",
        "-Wshadow",
        "-Werror",
    ]
    if sanitizer:
        command += [
            "-fno-omit-frame-pointer",
            "-fsanitize=address,undefined",
            "-fno-sanitize-recover=all",
        ]
    command += [
        str(HERE / "code_identity.c"),
        str(HERE / "test_code_identity_c.c"),
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True)


def encode_command(
    profile: int,
    field: int,
    original_count: int,
    recovery_count: int,
    metadata: Sequence[Metadata],
) -> str:
    fields = [
        "E", str(profile), str(field), str(original_count),
        str(recovery_count), str(len(metadata)),
    ]
    for item in metadata:
        fields += [str(item.type), item.value.hex() or "-"]
    return " ".join(fields)


def expected_decode(data: bytes) -> str:
    try:
        identity = deserialize(data)
    except IdentityError:
        return "ERR"
    return "OK " + serialize(identity).hex()


def check_encode(
    oracle: COracle,
    profile: int,
    field: int,
    original_count: int,
    recovery_count: int,
    metadata: Sequence[Metadata] = (),
) -> bytes:
    try:
        identity = make_identity(
            profile, field, original_count, recovery_count, metadata
        )
        expected = "OK " + serialize(identity).hex()
    except IdentityError:
        expected = "ERR"
    actual = oracle.query(
        encode_command(
            profile, field, original_count, recovery_count, metadata
        )
    )
    if actual != expected:
        raise AssertionError(
            f"encode mismatch: {profile=} {field=} {original_count=} "
            f"{recovery_count=} metadata={len(metadata)}\n"
            f"expected {expected[:200]}\nactual   {actual[:200]}"
        )
    return bytes.fromhex(expected[3:]) if expected.startswith("OK ") else b""


def check_decode(oracle: COracle, data: bytes) -> bool:
    expected = expected_decode(data)
    actual = oracle.query("D " + data.hex())
    if actual != expected:
        raise AssertionError(
            f"decode mismatch for {data[:64].hex()} ({len(data)} bytes): "
            f"expected {expected[:200]}, actual {actual[:200]}"
        )
    return expected.startswith("OK ")


def raw_identity_with_tlvs(
    tlvs: Iterable[tuple[int, int, bytes]],
    *,
    profile: int = PROFILE_LOW,
    field: int = FIELD_GF16,
    original_count: int = 129,
    recovery_count: int = 100,
) -> bytes:
    result = bytearray(serialize(make_identity(
        profile, field, original_count, recovery_count
    )))
    items = tuple(tlvs)
    struct.pack_into(">H", result, 34, len(items))
    for type_id, declared_bytes, value in items:
        result += struct.pack(">HH", type_id, declared_bytes)
        result += value
    struct.pack_into(">I", result, 8, len(result))
    return bytes(result)


def malformed_corpus() -> list[bytes]:
    valid = serialize(make_identity(
        PROFILE_LOW,
        FIELD_GF16,
        129,
        100,
        (Metadata(0x1234, b"future"),),
    ))
    corpus = [valid[:length] for length in range(len(valid))]
    corpus += [valid + b"\x00"]

    base = bytearray(serialize(make_identity(
        PROFILE_LEGACY_HIGH, FIELD_GF8, 240, 16
    )))
    mutations = (
        (0, 0), (4, 2), (5, 1), (7, 35), (12, 99), (13, 2),
        (14, 99), (15, 2), (31, 15), (33, 2), (35, 65),
    )
    for offset, value in mutations:
        candidate = bytearray(base)
        candidate[offset] = value
        corpus.append(bytes(candidate))
    for declared in (0, HEADER_BYTES - 1, HEADER_BYTES + 1, 0xFFFFFFFF):
        candidate = bytearray(base)
        struct.pack_into(">I", candidate, 8, declared)
        corpus.append(bytes(candidate))

    corpus += [
        raw_identity_with_tlvs(((0, 0, b""),)),
        raw_identity_with_tlvs(((7, 1, b"b"), (6, 1, b"a"))),
        raw_identity_with_tlvs(((6, 1, b"a"), (0x8006, 32, bytes(32)))),
        raw_identity_with_tlvs(((META_COORDINATE_SET_SHA256, 31, bytes(31)),)),
        raw_identity_with_tlvs(((0x9234, 6, b"future"),)),
        raw_identity_with_tlvs(((0x1234, 8, b"short"),)),
        raw_identity_with_tlvs(((0x1234, MAX_METADATA_VALUE + 1, b""),)),
        raw_identity_with_tlvs(((META_SHARD_LAYOUT, 0, b""),)),
        raw_identity_with_tlvs(((META_SHARD_LAYOUT, 1, b"\x00"),)),
        raw_identity_with_tlvs(((META_SHARD_LAYOUT, 1, b"\x02"),)),
        raw_identity_with_tlvs(((META_SHARD_LAYOUT, 2, b"\x01\x00"),)),
        *(raw_identity_with_tlvs(((base_type, 9, b"downgrade"),))
          for base_type in range(1, 6)),
        raw_identity_with_tlvs(
            ((META_SHARD_LAYOUT, 1, b"\x01"),),
            field=FIELD_GF8,
            original_count=2,
            recovery_count=5,
        ),
        raw_identity_with_tlvs((
            (META_SHARD_LAYOUT & 0x7FFF, 9, b"duplicate"),
            (META_SHARD_LAYOUT, 1, b"\x01"),
        )),
    ]
    excessive_count = bytearray(base)
    struct.pack_into(">H", excessive_count, 34, MAX_METADATA + 1)
    corpus.append(bytes(excessive_count))
    return corpus


def run_differential(oracle: COracle) -> dict[str, int]:
    counts: dict[str, int] = {
        "golden": 0,
        "profile": 0,
        "count_edges": 0,
        "metadata": 0,
        "size_limit_valid": 0,
        "oversize_rejected": 0,
        "malformed": 0,
        "mutations": 0,
        "mutation_accepted": 0,
    }

    vectors = json.loads((HERE / "golden_vectors.json").read_text())
    for vector in vectors:
        metadata = tuple(
            Metadata(int(item["type"], 0), bytes.fromhex(item["hex"]))
            for item in vector["metadata"]
        )
        encoded = check_encode(
            oracle, vector["profile"], vector["field"], vector["K"],
            vector["R"], metadata,
        )
        if encoded.hex() != vector["serialized_hex"]:
            raise AssertionError(f"golden mismatch: {vector['name']}")
        if hashlib.sha256(encoded).hexdigest() != vector["sha256"]:
            raise AssertionError(f"golden digest mismatch: {vector['name']}")
        if not check_decode(oracle, encoded):
            raise AssertionError(f"C rejected golden: {vector['name']}")
        counts["golden"] += 1

    for profile in (PROFILE_LEGACY_HIGH, PROFILE_LOW, PROFILE_EXACT_LOW):
        for original_count in COUNTS:
            for recovery_count in COUNTS:
                for field in (FIELD_GF8, FIELD_GF16):
                    encoded = check_encode(
                        oracle, profile, field, original_count, recovery_count
                    )
                    if encoded:
                        if not check_decode(oracle, encoded):
                            raise AssertionError("canonical profile rejected")
                        counts["profile"] += 1
    if counts["profile"] != 6085:
        raise AssertionError(f"unexpected profile count {counts['profile']}")
    for profile in (PROFILE_LEGACY_HIGH, PROFILE_LOW, PROFILE_EXACT_LOW):
        for original_count, recovery_count in COUNT_EDGE_PAIRS:
            for field in (FIELD_GF8, FIELD_GF16):
                encoded = check_encode(
                    oracle, profile, field, original_count, recovery_count
                )
                if encoded and not check_decode(oracle, encoded):
                    raise AssertionError("canonical edge-count profile rejected")
                counts["count_edges"] += 1
    if counts["count_edges"] != 66:
        raise AssertionError(
            f"unexpected count-edge cases {counts['count_edges']}"
        )

    for count in range(MAX_METADATA + 1):
        metadata = tuple(
            Metadata(
                type_id,
                bytes((type_id * 17 + index) & 0xFF
                      for index in range(type_id % 257)),
            )
            for type_id in range(6, 6 + count)
        )
        # Reverse the command order periodically: both make_identity and the C
        # builder must independently produce canonical type ordering.
        command_metadata = tuple(reversed(metadata)) if count & 1 else metadata
        encoded = check_encode(
            oracle, PROFILE_LOW, FIELD_GF16, 129, 100, command_metadata
        )
        if not check_decode(oracle, encoded):
            raise AssertionError("canonical metadata-count case rejected")
        counts["metadata"] += 1
    for length in range(MAX_METADATA_VALUE + 1):
        item = Metadata(
            0x1234,
            bytes((index * 29 + length) & 0xFF for index in range(length)),
        )
        encoded = check_encode(
            oracle, PROFILE_LEGACY_HIGH, FIELD_GF8, 240, 16, (item,)
        )
        if not check_decode(oracle, encoded):
            raise AssertionError("canonical metadata-length case rejected")
        counts["metadata"] += 1
    if counts["metadata"] != 4162:
        raise AssertionError(f"unexpected metadata count {counts['metadata']}")

    exact_limit = tuple(
        Metadata(type_id, bytes([type_id]) * 4096)
        for type_id in range(6, 21)
    ) + (Metadata(21, bytes([21]) * 3995),)
    encoded = check_encode(
        oracle, PROFILE_LOW, FIELD_GF16, 129, 100, exact_limit
    )
    if len(encoded) != 65535:
        raise AssertionError(f"maximum identity has {len(encoded)} bytes")
    if not check_decode(oracle, encoded):
        raise AssertionError("maximum-size canonical identity rejected")
    counts["size_limit_valid"] += 1
    oversized = exact_limit[:-1] + (Metadata(21, bytes([21]) * 3996),)
    if check_encode(
        oracle, PROFILE_LOW, FIELD_GF16, 129, 100, oversized
    ):
        raise AssertionError("oversized identity unexpectedly serialized")
    counts["oversize_rejected"] += 1

    for candidate in malformed_corpus():
        if check_decode(oracle, candidate):
            raise AssertionError("malformed corpus unexpectedly accepted")
        counts["malformed"] += 1

    rng = random.Random(0x4C324944)
    valid = serialize(make_identity(
        PROFILE_LOW, FIELD_GF16, 129, 100,
        (Metadata(0x1234, b"future"),),
    ))
    for _ in range(20000):
        candidate = bytearray(valid)
        for _ in range(rng.randrange(1, 5)):
            candidate[rng.randrange(len(candidate))] ^= 1 << rng.randrange(8)
        if check_decode(oracle, bytes(candidate)):
            counts["mutation_accepted"] += 1
        counts["mutations"] += 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cc", default=os.environ.get("CC", "gcc"))
    parser.add_argument(
        "--sanitizers", action="store_true",
        help="compile the C oracle with AddressSanitizer and UBSan",
    )
    arguments = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="leo2-code-identity-c-") as temp:
        executable = pathlib.Path(temp) / "test_code_identity_c"
        compile_oracle(arguments.cc, executable, arguments.sanitizers)
        subprocess.run([str(executable)], check=True)
        oracle = COracle(executable)
        try:
            counts = run_differential(oracle)
        finally:
            oracle.close()

    print(
        "C/Python code-identity differential passed: "
        + " ".join(f"{key}={value}" for key, value in counts.items())
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
