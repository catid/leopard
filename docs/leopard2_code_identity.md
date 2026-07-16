# Leopard2 serialized code-identity experiment

Status: Experiment W reference proposal, not public API or frozen wire format.
Bead: `leopard-79h.18.23`.  Validation date: 2026-07-16 UTC.

The reference implementation is
`experiments/leopard2/code_identity/code_identity.py`.  Its identifier describes
the code mathematics and transmitted layout.  It deliberately has no backend,
SIMD ISA, thread count, tile size, dispatcher threshold, calibration, or
autotuning field.

## Version 1 binary proposal

All multibyte integers are unsigned big-endian (network order).  The fixed
header is 36 bytes:

| Offset | Bytes | Meaning |
| ---: | ---: | --- |
| 0 | 4 | ASCII magic `L2ID` |
| 4 | 1 | envelope format version, `1` |
| 5 | 1 | reserved flags, zero |
| 6 | 2 | header bytes, `36` |
| 8 | 4 | total serialized bytes |
| 12 | 1 | profile family |
| 13 | 1 | profile version |
| 14 | 1 | field family |
| 15 | 1 | field-representation version |
| 16 | 4 | systematic count `K` |
| 20 | 4 | transmitted recovery count `R` |
| 24 | 4 | parent count `N` |
| 28 | 4 | padded side (`T` for high, `P` for low) |
| 32 | 2 | coordinate-map version |
| 34 | 2 | metadata TLV count |

Profile family 1 is legacy high and family 2 is low; both currently have profile
version 1.  Field family 1 is legacy GF8 and family 2 is legacy GF16, both with
field-representation version 1.  Coordinate-map version 1 means the exact maps
in `docs/leopard2_math_and_sources.md`.  Redundant `N` and padded-side values are
intentional: readers recompute them and reject an identity whose claimed values
do not match its profile.

Each metadata record is a 16-bit type, 16-bit value length, then the value.  TLVs
must be strictly increasing by type and base types must be unique.  Type zero is
reserved.  Bit 15 marks critical metadata.  Known critical types 1 through 4 are
32-byte SHA-256 digests of canonical profile parameters, evaluation coordinates,
shortening coordinates, and puncturing coordinates respectively.  Unknown
critical types fail closed.  Unknown noncritical types are retained byte for
byte, so an old reader can parse and reserialize informational extensions without
changing the identity.  Metadata is limited to 64 values of at most 4096 bytes;
the complete identifier is limited to 65535 bytes.

## Canonical and forward-compatible behavior

There is one serialization for an identity: zero reserved bits, exact lengths,
minimal fixed-width fields, sorted unique TLVs, and no trailing bytes.  A reader
rejects noncanonical input instead of normalizing it silently.  Unknown envelope,
profile, field-representation, coordinate-map versions, and critical metadata are
semantic changes and therefore fail closed.  A future implementation must assign
a new version rather than interpreting old bytes according to a new coordinate
map.  Only unknown noncritical metadata may be skipped or preserved by an older
reader.

The experimental exact profile is intentionally absent.  Its coordinate set is
not frozen; assigning it a persistent identifier now would create an accidental
compatibility promise.  When frozen, it needs a new profile family/version and a
critical coordinate-set digest or an equally complete versioned map.

## Migration rules

An application migrating from the legacy API creates family 1/version 1,
field/version 1, map version 1, and derives `T` and `N` from `K,R`.  It should
store the full identifier beside parity data and compare the parsed mathematical
fields before decoding.  A bare `(K,R)` tuple is insufficient because it omits
profile, field representation, parent, and coordinate map.  Existing data is not
rewritten.  Low-profile data uses family 2/version 1 and must never be decoded as
legacy high even if `K,R,N` happen to agree.

Metadata that only names a CPU backend or measured execution policy is forbidden.
Deployment labels and checksums of shard contents also belong outside this
identifier.  Critical hashes are appropriate only for immutable algebraic sets
whose canonical byte representation is itself documented.

## Validation and disposition

Run:

    cd experiments/leopard2/code_identity
    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v test_code_identity.py

The suite checks golden bytes and hashes, 3,952 profile/field/count boundary
round trips, 4,162 metadata count/length round trips, every truncation, trailing
bytes, reserved bits, bad lengths, integer overflow, inconsistent derived fields,
unknown critical extensions, malformed digest TLVs, ordering/duplicate rules,
and 20,000 deterministic mutation cases.  Unknown noncritical data is explicitly
round-tripped.  The test has no third-party dependencies.

Promotion requires independent C and Python implementations to match every
golden vector, public malformed-input fuzzing under sanitizers, a frozen exact
profile decision, and an API/storage use case that justifies committing a public
wire contract.  Kill or revise this proposal if the fixed header cannot represent
a frozen profile without profile-specific interpretation, forward-compatible
parsing becomes ambiguous, or the envelope adds material operational complexity
without preventing real profile mismatches.  Until those gates pass, this code
remains isolated experimental evidence.
