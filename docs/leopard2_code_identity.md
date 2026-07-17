# Leopard2 serialized code-identity experiment

Status: Experiment W reference proposal, not public API or a frozen production
wire contract.  Family 3/version 1/map 1 nevertheless freezes the experimental
C7 exact-low mathematics so retained research artifacts cannot reinterpret it.
Bead: `leopard-79h.18.23`.  Validation date: 2026-07-17 UTC.

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
Family 4 is reserved for any future serialization of C8's default-off
exact-high candidate and version-1 readers reject it.  The current C8
checkpoint is deliberately unserialized; exact-high must never reinterpret
family 3.

For family 3/map 1, `N=K+R`, the padded-side field equals `K`, and the ordered
evaluation coordinates are the first `K+R` Leopard Cantor coordinates.  The
first `K` are systematic and the next `R` are transmitted parity.  There are
no shortened or punctured positions, so shortening- or puncturing-set metadata
is rejected.  Field family 1 is legacy GF8 and field family 2 is legacy GF16,
both with field-representation version 1.  Coordinate-map version 1 means the
exact maps in `docs/leopard2_math_and_sources.md` and the C7 addendum.  Redundant
`N` and padded-side values are intentional: readers recompute them and reject
an identity whose claimed values do not match its profile.  Because map 1
already determines every coordinate and has no shortened or punctured set,
family 3 rejects coordinate-, shortening-, and puncturing-set digest TLVs as
redundant noncanonical spellings, regardless of digest contents.

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
Family 4 remains an unimplemented serialization reservation, not a decodable
identity.  C8's executable candidate and evidence remain independently named
`exact_high_prefix_v1_candidate` until a later C9/C10 decision freezes a
decoder and serialized profile.

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

Build the isolated coverage-guided target and generate its deterministic seed
corpus with:

    cmake -S ../../.. -B ../../../build/code-id-fuzz \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=OFF \
      -DLEO2_BUILD_FUZZERS=ON -DLEO2_ENABLE_CUDA=OFF \
      -DCMAKE_C_COMPILER=clang-18 -DCMAKE_CXX_COMPILER=clang++-18
    JOBS="$(nproc)"; if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    cmake --build ../../../build/code-id-fuzz \
      --target leopard2_code_identity_fuzzer -j "$JOBS"
    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v \
      test_fuzz_campaign.py
    test ! -e ../../../build/code-id-corpus-independent
    PYTHONDONTWRITEBYTECODE=1 python3 make_fuzz_corpus.py \
      --output ../../../build/code-id-corpus-independent
    test ! -e ../../../build/code-id-campaign-independent
    PYTHONDONTWRITEBYTECODE=1 python3 run_fuzz_campaign.py \
      --fuzzer ../../../build/code-id-fuzz/leopard2_code_identity_fuzzer \
      --corpus ../../../build/code-id-corpus-independent \
      --results ../../../build/code-id-campaign-independent \
      --runs 50000 --max-len 65535 --rss-limit-mb 8192 \
      --external-rss-limit-mb 10240

Use the allowed CPU count rather than a literal 128 when fewer CPUs are in the
process affinity mask.  The explicit libFuzzer RSS limit avoids a Clang 18
high-water-accounting false positive seen in long-lived launch environments;
an experiment runner should still impose an external per-job memory cap.

The same complete C/Python differential can execute as a real big-endian s390x
binary under QEMU.  On Debian/Ubuntu the optional validation-only packages are
`gcc-s390x-linux-gnu`, `libc6-dev-s390x-cross`, and `qemu-user-static`; they are
not build or runtime dependencies of Leopard:

    PYTHONDONTWRITEBYTECODE=1 python3 test_code_identity_c.py \
      --cc s390x-linux-gnu-gcc --static \
      --runner /usr/bin/qemu-s390x-static

The suite checks golden bytes and hashes, 6,085 profile/field/count boundary
round trips, 4,162 metadata count/length round trips, every truncation, trailing
bytes, reserved bits, bad lengths, integer overflow, inconsistent derived fields,
unknown critical extensions, malformed digest TLVs, ordering/duplicate rules,
and 20,000 deterministic mutation cases.  Unknown noncritical data is explicitly
round-tripped.  The test has no third-party dependencies.

The independent C99 implementation is `code_identity.c` with its deliberately
experimental header `code_identity.h`.  It is not linked into `libleopard`,
installed, or exposed through the Leopard public ABI.  Root CMake sees it only
inside the default-off `leopard2_code_identity_fuzzer` target when
`LEO2_BUILD_FUZZERS=ON`; a tests/benchmarks/fuzzers-off archive contains neither
identity nor fuzzer symbols.  It parses metadata as borrowed views into
caller-owned input and serializes into caller-owned storage; the implementation
performs no allocation.  Borrowed input and metadata storage must not overlap the
destination objects.  The size-query operation reports the exact required buffer
size and undersized output fails before any write.  Counts, TLV values, TLV count,
and complete identifiers retain the same 4096-byte, 64-entry, and 65535-byte
bounds as the Python reference.

`test_code_identity_c.c` contains 217 direct C checks and also provides a test-
only line protocol.  `test_code_identity_c.py` compiles it with strict C99 GCC
warnings and compares both C serialization and C deserialization with the
Python implementation.  The deterministic differential matrix covers all 5
golden vectors, 6,085 valid profile/field/count combinations, 66 zero, field-
boundary, and `UINT32_MAX` count cases, 4,162 metadata count/length cases, one
exact 65,535-byte valid identifier, one 65,536-byte rejection, 87 explicitly
malformed identifiers (including six exact-low redundant-digest encodings),
and 20,000 deterministic mutated inputs.
Of the mutations, 1,384 remain valid canonical identities; both implementations
accept and reserialize those exact bytes, while both reject the other 18,616.
The same complete matrix passes under combined GCC AddressSanitizer and
UndefinedBehaviorSanitizer.  Sanitized truncation tests use exact-size backing
allocations, so an out-of-bounds parser read cannot hide inside a larger test
buffer.  The harness compiles into an automatically removed temporary directory
and leaves no generated binary or Python bytecode in the source tree.

The initial campaign at implementation commit `d5d605b` is rejected as evidence:
all 28 workers reused seed `1279609156` while mutating one shared corpus.  It
found no crash, but it was neither a set of reproducible independent campaigns
nor compliant with the logged-distinct-seed requirement.  The replacement
runner derives a stable nonzero seed per worker, gives every affinity-pinned
worker an isolated input corpus and artifact directory, validates exact logged
seed and execution counts, rejects sanitizer markers and timeouts, retains each
log, enforces an externally observed resident-memory limit in addition to
libFuzzer's own accounting, and merges successful corpora by content hash.  It
refuses dirty source or nonempty result directories.  `test_fuzz_campaign.py`
covers stable distinct seeds, isolated corpora, order-independent merging, stale
directories, malformed logs, sanitizer markers, a bounded timeout, and the
external RSS cap.  A replacement 28-worker campaign must update
`fuzz_checkpoint.json` before the fuzz gate is considered satisfied.

The same 217 direct C checks and complete differential matrix passed in a
statically linked s390x binary under QEMU 8.2.2 with zero mismatch, providing
executed big-endian evidence rather than an inference from byte-wise loads.  Raw
corpora and worker logs remain ignored build artifacts.

The independent implementation, deterministic sanitized differential, and
emulated big-endian research gates are satisfied, including the frozen
experimental exact-low identity.  The coverage-guided target is implemented,
but its independent-seed campaign gate remains pending the replacement run
described above.  Public promotion still requires a concrete API/storage use
case and a deliberate
decision that the remaining exact-profile mathematics are stable enough to
justify a long-term wire promise.  Kill or revise this proposal if the fixed
header cannot represent a frozen profile without profile-specific
interpretation, forward-compatible parsing becomes ambiguous, or the envelope
adds material operational complexity without preventing real profile
mismatches.  Until that decision, this code remains isolated experimental
evidence rather than public ABI.
