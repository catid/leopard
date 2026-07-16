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

Profile family 1 is legacy high, family 2 is padded low, and family 3 is exact
prefix low; all currently have profile version 1.  Family 3 is the persistent
identity corresponding to public enum value
`LEO2_PROFILE_EXACT_EXPERIMENTAL_V1 == 3`; the stable codec constructor still
rejects it, so this is a research-format freeze rather than production ABI.
Family 4 is reserved for a possible future exact-high profile and version-1
readers reject it.  Exact-high must never reinterpret family 3.

For family 3/map 1, `N=K+R`, the padded-side field equals `K`, and the ordered
evaluation coordinates are the first `K+R` Leopard Cantor coordinates.  The
first `K` are systematic and the next `R` are transmitted parity.  There are
no shortened or punctured positions, so shortening- or puncturing-set metadata
is rejected.  Field family 1 is legacy GF8 and field family 2 is legacy GF16,
both with field-representation version 1.  Coordinate-map version 1 means the
exact maps in `docs/leopard2_math_and_sources.md` and the C7 addendum.  Redundant
`N` and padded-side values are intentional: readers recompute them and reject
an identity whose claimed values do not match its profile.

Each metadata record is a 16-bit type, 16-bit value length, then the value.  TLVs
must be strictly increasing by type and base types must be unique.  Type zero is
reserved.  Bit 15 marks critical metadata.  Base types `0x0001` through
`0x0005` are reserved for their critical forms and are rejected, so clearing
bit 15 cannot silently downgrade their semantics.  Critical types `0x8001` through
`0x8004` are 32-byte SHA-256 digests of canonical profile parameters, evaluation
coordinates, shortening coordinates, and puncturing coordinates respectively.
Critical type `0x8005` is the one-byte shard-layout selector.  Its only current
value is `1`, `GF16_PADDED_ODD_V1`; it is valid only with GF16.  Absence of the
TLV is the single canonical spelling of `NATIVE_V1`.  An explicit value zero is
therefore noncanonical, and other lengths or values fail closed.  This preserves
every previously generated native identifier byte for byte while making the
padded systematic framing unambiguous.  Unknown critical types fail closed.
Unknown noncritical types are retained byte for byte, so an old reader can parse
and reserialize informational extensions without changing the identity.
Metadata is limited to 64 values of at most 4096 bytes; the complete identifier
is limited to 65535 bytes.

## Canonical and forward-compatible behavior

There is one serialization for an identity: zero reserved bits, exact lengths,
minimal fixed-width fields, sorted unique TLVs, and no trailing bytes.  A reader
rejects noncanonical input instead of normalizing it silently.  Unknown envelope,
profile, field-representation, coordinate-map versions, and critical metadata are
semantic changes and therefore fail closed.  A future implementation must assign
a new version rather than interpreting old bytes according to a new coordinate
map.  Only unknown noncritical metadata may be skipped or preserved by an older
reader.

The C7 exact-low coordinate set is frozen by family 3/version 1/map 1 rather
than by an optional digest.  Changing its point order, field representation, or
systematic/parity partition requires a new profile or coordinate-map version.
Family 4 remains an unimplemented reservation, not a decodable identity.

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
    PYTHONDONTWRITEBYTECODE=1 python3 test_code_identity_c.py --cc gcc
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      PYTHONDONTWRITEBYTECODE=1 python3 test_code_identity_c.py \
        --cc gcc --sanitizers

The suite checks golden bytes and hashes, 6,085 profile/field/count boundary
round trips, 4,162 metadata count/length round trips, every truncation, trailing
bytes, reserved bits, bad lengths, integer overflow, inconsistent derived fields,
unknown critical extensions, malformed digest TLVs, ordering/duplicate rules,
and 20,000 deterministic mutation cases.  Unknown noncritical data is explicitly
round-tripped.  The test has no third-party dependencies.

The independent C99 implementation is `code_identity.c` with its deliberately
experimental header `code_identity.h`.  It is not included by the root CMake
build, installed, or exposed through the Leopard public ABI.  It parses metadata
as borrowed views into caller-owned input and serializes into caller-owned
storage; the implementation performs no allocation.  Borrowed input and metadata
storage must not overlap the destination objects.  The size-query operation
reports the exact required buffer size and undersized output fails before any
write.  Counts, TLV values, TLV count, and complete identifiers retain the same
4096-byte, 64-entry, and 65535-byte bounds as the Python reference.

`test_code_identity_c.c` contains 211 direct C checks and also provides a test-
only line protocol.  `test_code_identity_c.py` compiles it with strict C99 GCC
warnings and compares both C serialization and C deserialization with the
Python implementation.  The deterministic differential matrix covers all 5
golden vectors, 6,085 valid profile/field/count combinations, 66 zero, field-
boundary, and `UINT32_MAX` count cases, 4,162 metadata count/length cases, one
exact 65,535-byte valid identifier, one 65,536-byte rejection, 81 explicitly
malformed identifiers, and 20,000 deterministic mutated inputs.
Of the mutations, 1,384 remain valid canonical identities; both implementations
accept and reserialize those exact bytes, while both reject the other 18,616.
The same complete matrix passes under combined GCC AddressSanitizer and
UndefinedBehaviorSanitizer.  Sanitized truncation tests use exact-size backing
allocations, so an out-of-bounds parser read cannot hide inside a larger test
buffer.  The harness compiles into an automatically removed temporary directory
and leaves no generated binary or Python bytecode in the source tree.

The independent-implementation and deterministic sanitized differential gates
are now satisfied, including the frozen experimental exact-low identity.
Promotion still requires a sustained coverage-guided C API fuzzing campaign,
cross-endian or emulated-endian validation, and a concrete API/storage use case
that justifies committing a public wire contract.  Kill or revise this proposal
if the fixed header cannot represent a frozen profile without profile-specific
interpretation, forward-compatible parsing becomes ambiguous, or the envelope
adds material operational complexity without preventing real profile
mismatches.  Until those gates pass, this code remains isolated experimental
evidence rather than public ABI.
