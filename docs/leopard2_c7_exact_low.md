# Leopard2 C7 exact low-rate profile checkpoint

## Outcome and status

C7 selects and freezes a separately versioned **experimental** exact-low
mathematical profile.  It does not promote a production codec or dispatcher
rule.  The stable constructor continues to return `LEO2_UNSUPPORTED` for public
enum value `LEO2_PROFILE_EXACT_EXPERIMENTAL_V1 == 3`, and AUTO never selects
it.  The default CMake graph does not compile the C7 implementation.

The family-3 additions first recorded in commit `d154cc3` were a provisional
Experiment-W envelope checkpoint.  The independent coordinate study in this
checkpoint is the decision that justifies the experimental freeze: prefix
coordinates retain the strong C6 field-boundary evidence, use one simple map,
and avoid a generated coordinate table.  The aligned dyadic-union alternative
is frequently unavailable at the GF8 boundary and is generally not an affine
image of the prefix.  The offline search finds small-field improvements in
some symbolic cells, but those sets change the generator, lack C6 execution
evidence, and would require another coordinate-map version.

This checkpoint includes a scalar-setup/SIMD-execution encoder and original
repair plan for both legacy GF8 and GF16 representations, an independent direct
field oracle, arbitrary valid physical tails, requested parity subsets,
maximum-loss and missing-parity tests, re-encoding after recovery, allocation
interposition, strict builds, and ASan+UBSan.  Its timing harness records raw
setup and execution samples for exact and padded low profiles.  The retained
affinity-selected smoke is explicitly non-authoritative, is rebuilt from the
same committed source as the correctness matrix, and cannot support a promotion
claim.
Authoritative crossover timing must be rebuilt from the final integrated
production commit.

The v4 evidence tooling and post-R1 sanitizer contract are prepared, but final
A/B artifacts are intentionally not regenerated in this change.  Generation
must use the final integrated core containing the forced-scalar warning fix and
must first reconfirm the fail-closed sanitizer attribution for that exact
22-file closure.

## Exact-low family 3/version 1/map 1

For positive public counts `K,R` with `K+R` no greater than the selected field
order, define the ordered evaluation coordinates in Leopard Cantor-coordinate
notation as:

    systematic: omega_0, omega_1, ..., omega_(K-1)
    parity:     omega_K, omega_(K+1), ..., omega_(K+R-1)
    degree:     less than K

The first `K` transmitted shards are the systematic values.  The following `R`
transmitted shards are the polynomial evaluations at the parity coordinates.
There is no dyadic parent, no shortening, and no puncturing.  In the existing
36-byte `L2ID` Experiment-W envelope:

- profile family is 3, matching public enum value 3;
- profile version is 1;
- field/version is legacy GF8/1 or legacy GF16/1;
- `parent_count` is exactly `K+R`;
- `padded_side` is exactly `K`;
- coordinate-map version is 1.

Profile family 4 is reserved for any future serialization of C8's completed,
default-off exact-high candidate and is rejected by the version-1 readers.
C8 remains an unserialized research profile; exact-high must not reinterpret
family 3.  Coordinate-,
shortening-, and puncturing-set digest TLVs are invalid for exact low V1 because
the versioned map already fixes the coordinates and the latter two sets do not
exist.  Rejecting all three values, including all-zero and all-one digests,
leaves one canonical identity spelling.  The stable codec API remains unchanged
and unsupported.

### Canonical coordinate bytes

The map is fixed by the version fields, so a coordinate-set digest is redundant
and is rejected in a family-3 `L2ID`.  An external research artifact may still
fingerprint the set.  Its one canonical byte representation is:

1. ASCII `LEO2-EXACT-LOW-PREFIX-V1` followed by a zero byte;
2. one byte containing field width, 8 or 16;
3. unsigned big-endian 32-bit `K`, then unsigned big-endian 32-bit `R`;
4. role byte 1, unsigned big-endian 32-bit `K`, then the systematic coordinate
   indices in order;
5. role byte 2, unsigned big-endian 32-bit `R`, then the parity coordinate
   indices in order.

A coordinate index is one byte for GF8 and two unsigned big-endian bytes for
GF16.  It is the coefficient bit mask in Leopard's declared Cantor basis, not
the corresponding polynomial-basis integer.  For GF8 that basis is
`1,214,152,146,86,200,88,230` over polynomial `0x11d`; for GF16 it is the 16
values declared at the top of `LeopardFF16.cpp` over polynomial `0x1002d`.
Changing any basis, point order, role partition, or coordinate-byte encoding
requires a new field/profile/map version.

## Generator derivation and density

Write `x_i=omega_i` for systematic point `i`, `q_j=omega_(K+j)`, and

    Z(X) = product over 0 <= s < K of (X + x_s).

Characteristic two makes subtraction and addition identical.  The systematic
parity row is

    G[j,i] = Z(q_j) / ((q_j + x_i) Z'(x_i))

where

    Z'(x_i) = product over s != i of (x_i + x_s).

Every evaluation point is distinct.  Therefore every factor in the numerator
and denominator is nonzero, and every one of the `K*R` parity coefficients is
nonzero.  This is an intrinsically dense direct map; coordinate search cannot
turn it into a sparse generator while retaining distinct evaluation points.

The Python oracle derives the same rows independently by constructing the
monomial Vandermonde matrix on the systematic points, inverting it with field
Gaussian elimination, and multiplying parity power vectors by that inverse.
The C++ byte executor separately derives barycentric rows with production field
helpers and checks every coefficient against carryless polynomial arithmetic
and the declared basis maps.  A declared GF16 `K=3,R=500` case checks all 500
Vandermonde rows (1,500 coefficients) in both independent oracles rather than
sampling a few rows.

## Affine invariance

For one global affine coordinate change

    phi(X) = a X + b, with a != 0,

each difference acquires one factor `a`.  `Z(phi(q))` has `K` factors of `a`;
the explicit `(phi(q)+phi(x_i))` term and the `K-1` factors in the derivative
together have the same `K` factors.  They cancel, so every systematic generator
coefficient is unchanged.  This proves that a global scale/translation cannot
improve density, coefficient-one count, wire bytes, or direct execution cost.

The test exhausts all 240 invertible affine maps of GF(16) for all 120 valid
prefix `(K,R)` geometries: 28,800 mapped generators and 734,400 coefficient
comparisons.  GF8/GF16 sampling adds 177,600 complete transformed-row
comparisons.  This theorem applies only to one global affine image of the whole
ordered set; it is not used to claim equivalence for a non-affine union.

## Prefix, aligned union, and searched sets

The real aligned low alternative studied is:

    P = ceil_pow2(K)
    systematic = omega_0 .. omega_(K-1)
    parity     = omega_P .. omega_(P+R-1)

It can be viewed as removing the padded-low shortened interval from
transmission while retaining an aligned parity start.  It is unavailable when
`P+R` exceeds the field order.  When `K` is not a power of two, the ordered
systematic-plus-parity set is generally not one global affine image of the
prefix, and its generator must be treated as different.

Across GF(16), the aligned union exists for 85 of 120 geometries and is globally
affine in only 49.  The non-affine cases change 758 coefficients.  Both maps
remain fully dense across the 1,812 compared coefficients.  Complete aligned
dyadic-block fragments decrease from 358 to 328, but coefficient-one
specializations also decrease from 297 to 198.  In the large-field witnesses:

- GF8 `(3,253)`, `(7,249)`, `(127,129)`, and `(248,8)` cannot use the aligned
  union because `P+R` crosses 256;
- GF16 `(129,100)` reduces symbolic fragments from 11 to 5 but changes all
  12,900 coefficients and loses 152 coefficient-one specializations;
- GF16 `(1000,17)` reduces fragments from 9 to 8 and changes all 17,000
  coefficients.

The transform-aware GF(16) search keeps systematic points `0..K-1` and
enumerates every nonempty parity subset of the remaining field: 65,519
candidates.  It minimizes complete dyadic fragments, then maximizes coefficients
equal to one, then breaks ties lexicographically.  It chooses a non-prefix set
in 57 of 120 cells.  A separate full-field partition search examines all 65,534
nontrivial systematic sets.  Every candidate remains dense.  These symbolic
wins are not sufficient for promotion: they change parity bytes, require a map
table or another deterministic rule, and have no end-to-end evidence meeting
the 10% threshold.  Prefix map 1 is consequently the bounded C7 choice.

## Experimental encoder and decoder

`experiments/leopard2/non_power_of_two/c7/c7_exact_low.cpp` is not part of the
default build.  Immutable codec setup stores `K*R` nonzero multiplier logs.
Encoding evaluates exactly the requested transmitted parity outputs; a null
parity output is skipped.  Execution has zero transform scratch and allocates
nothing.  Setup is quadratic in the transmitted generator size, so this direct
implementation is a crossover candidate, not a balanced-code production
answer.

The immutable original-repair plan accepts a sorted missing-original set and a
parity-presence bitmap, deterministically selects the lowest available parity
equations, inverts only the missing-original minor, and folds surviving-original
terms into fixed execution coefficients.  Execution restores missing originals
only.  A no-loss plan returns without inspecting byte count or any pointer.
Before the first output write, execution validates every restored destination,
every selected parity term, and every surviving-original term.  A null required
term or any unsupported overlap therefore rejects the whole call without
partial output.  Parity rebuild is explicit re-encoding and is checked byte for
byte.

GF8 accepts every positive physical byte count.  GF16 accepts positive even
physical counts in Leopard's established full-tile/compact-tail layout.  C7
explicitly rejects odd GF16 physical sizes; applications must use the separately
versioned padded-odd shard layout to represent an odd payload.  Read-only input
shards may alias each other.  Every non-null encode output must be disjoint from
all inputs and other outputs.  Every restored decode output must be disjoint
from every received input and other restored output.  Invalid overlap is
rejected before a fixed-multiplier helper runs.  The decode alias test uses a
valid nonzero constant-polynomial codeword, points every surviving-original and
parity input at one shared read-only buffer, snapshots the complete shared
input, and restores one original into separate guarded storage for every
retained GF8 and GF16 byte length.

## Correctness checkpoint

The independent algebra result records:

| Gate | Result |
| --- | ---: |
| GF(16) public geometries | 120 |
| GF(16) MDS coordinate subsets | 131,038 |
| Dense coefficients matched to Vandermonde oracle | 3,060 |
| Exhaustive affine maps / coefficients | 28,800 / 734,400 |
| Fixed-systematic transform-search candidates | 65,519 |
| Full-field systematic partitions | 65,534 |
| Large GF8/GF16 geometries / dense coefficients | 10 / 84,781 |
| Declared GF16 Vandermonde rows / coefficients | 500 / 1,500 |

Each scalar, SSSE3, AVX2, and AUTO C++ artifact records the same:

| Gate | Per backend |
| --- | ---: |
| GF8/GF16 geometries | 9 / 5 |
| Independently checked coefficients | 118,717 |
| Independent GF16 Vandermonde coefficients | 1,500 |
| Encode executions / symbol comparisons | 117 / 1,030,423 |
| Requested-subset encodes | 117 |
| Decode plans/executions | 403 / 403 |
| Recovered-symbol comparisons | 272,487 |
| Maximum-loss plans | 117 |
| Plans with an unavailable low parity | 175 |
| Re-encode parity checks | 403 |
| Odd-GF16 / all overlap rejections | 10 / 59 |
| Parity-output / restored-output / restored-input overlap rejects | 13 / 12 / 20 |
| Null selected-parity / surviving-original rejects | 14 / 6 |
| Bytes checked unchanged after atomic rejection | 61,570 |
| Read-only input-alias calls / symbols checked | 13 / 2,139 |
| Decode read-only input-alias calls / symbols checked | 117 / 6,025 |
| Hot-path allocations | 0 |
| Deterministic digest | `0xec4179e9f2776a58` |

Lengths are GF8 `1,2,3,7,31,64,65,257` bytes and GF16
`2,4,6,14,62,64,66,130,514` physical bytes.  Buffers are unaligned and guarded.
Combined ASan+UBSan with leak detection passes the same complete matrix.  The
sanitized standalone source has a compile-time feature gate that fails unless
both instruments are active; the retained build/run manifest also records the
exact CMake-selected compiler, make, archiver, ranlib, and linker roles plus the
standalone link driver and compiler-selected linker.  It binds each tool's
bytes and complete `--version` output, the path-normalized retained CMake
cache, verbose configure, core-build, and standalone-link logs, linked archive,
executable, source closure, run environment, child affinity, and artifact
hashes.  Each build records the same exact 22-file core closure (all 11
translation units, their complete project-header dependency set, top-level
CMake input, and package-config template).  It also retains all 11 compiler
dependency files after exact checkout-root normalization; the validator
reparses those records and requires that they reconstruct precisely that
closure.  Each raw dependency record must also use the exact CMake target path
under its declared build directory, name the same retained dependency artifact,
keep that artifact in the exact dependency directory for the declared build
role, and declare that exact object target on the dependency file's left-hand
side; preserving a basename while relabeling or swapping the file within or
across builds is rejected.  The
build job count is a typed integer in `1..8`,
and the normalized build command must end in the corresponding exact `-jN`
argument.  The checker parses those logs rather than accepting a success
label.  In particular, the forced-scalar build must have empty stderr; its
compile-time-only backend selection explicitly consumes the otherwise unused
feature record without changing the selected scalar backend.  Manifest v4's
post-R1 probe freezes 320 ASan and 54 UBSan references in
the standalone harness and 348 ASan and 86 UBSan references across all 11 named
core-archive members; normal builds contain none.  The exact sanitized archive
attribution is:

| Archive member | ASan | UBSan |
| --- | ---: | ---: |
| `leopard.cpp.o` | 13 | 7 |
| `leopard2.cpp.o` | 141 | 15 |
| `Leopard2Backend.cpp.o` | 40 | 9 |
| `Leopard2BackendScalar.cpp.o` | 16 | 6 |
| `Leopard2CpuFeatures.cpp.o` | 9 | 5 |
| `Leopard2Plan.cpp.o` | 7 | 5 |
| `LeopardCommon.cpp.o` | 13 | 5 |
| `LeopardFF16.cpp.o` | 27 | 10 |
| `LeopardFF8.cpp.o` | 28 | 9 |
| `Leopard2BackendSSSE3.cpp.o` | 26 | 8 |
| `Leopard2BackendAVX2.cpp.o` | 28 | 7 |

The delta from retained v3's 329/87 archive attribution is expected after the
integrated backend work: active-parent locator setup, context-selected and
native fused XOR operations, SSSE3/AVX2/NEON additions, and GF8/GF16 kernel
refactoring changed the compiler's instrumented sites.  These totals are not a
range or minimum.  The v4 runner still compares the exact totals and every
named member, aborting before correctness execution if any count changes.  The
unchanged 320/54 standalone count independently confirms that the experimental
harness instrumentation did not drift.
The counts first came from a fresh five-variant compile at post-R1 commit
`808ca28`; the runner stopped on the old fail-closed constants before executing
any child.  That closure remained byte-identical through `219831f`.  The later
one-line forced-scalar warning fix changes `Leopard2Backend.cpp`, however, so
the final selected core must receive a fresh fail-closed sanitizer scan before
A/B generation.  The v4 constants remain exact expectations, not an inference
that a changed source closure preserves instrumentation.  No timing or
correctness result is taken from the stopped run.

### Dual Git identity and historical evidence

Manifest v4 separates two immutable Git identities:

- `tooling_git_sha` is the clean checked-out `HEAD` containing the exact
  committed C7 harness, runner, and validator.  Generation rejects staged,
  modified, or untracked non-ignored files.
- `core_git_sha` is that commit or an ancestor.  Every byte of the exact
  22-file `EXPECTED_SOURCE_CLOSURE` in the checkout and every build record must
  equal the corresponding blob at this commit.

This lets a tooling-only evidence repair describe an earlier frozen core
without pretending that older commit contained the repaired tools.  The runner
checks ancestry and both byte closures before configuring.  The validator
binds all three tooling records to `tooling_git_sha`, every core record to
`core_git_sha`, and checks ancestry again.  Trusted live validation also
requires a clean checkout at `tooling_git_sha`.  A/B comparison requires both
identities, tooling records, all five binary fingerprints, and all program
records to agree.  Peer-attestation schema v3 records both identities and the
clean-tooling-head/core-ancestor check.

The committed checkpoint under `results/` remains historical manifest v3 and
peer-attestation v2 evidence for its pre-R1 core.  It is never rewritten or
accepted by the v4 generator.  The validator retains a read-only v3 path that
authenticates its tooling and core bytes directly from the recorded Git commit
and applies its original exact 329 ASan / 87 UBSan archive proof.  Relabeling
those bytes as v4 fails both the 348/86 member proof and the required separate
tooling identity.

All nested identity and accounting comparisons are recursively type-strict.
JSON `false` does not satisfy an expected integer zero, `true` does not satisfy
field or profile identifier one, and Boolean values are never accepted as
numeric measurements.  Finite timing samples and their raw-derived summaries
may use numerically equivalent JSON integer or floating-point syntax; this
includes historical v3 evidence whose exact-zero MAD is serialized as integer
`0` rather than floating-point `0.0`.  Artifact sizes, normalization-token
counts, build jobs, CPU identifiers, scratch sizes, and all non-measurement
fields retain explicit type and range checks.  Mutation tests cover
ordinary-build sanitizer zeros,
per-member counts, exact-field and profile identifiers, correctness zeros,
affinity, sample summaries, A/B tooling records, and comparison scan counts.

Manifest v4 normalizes only the exact checkout-root prefix to literal
`${LEO2_SOURCE_ROOT}` in retained text and argv.  All core and standalone
compiles use supported `-ffile-prefix-map` and `-fdebug-prefix-map`, plus
`-fmacro-prefix-map` when both the paired C and C++ drivers accept it.  Their
stable binary mapping target is literal `LEO2_SOURCE_ROOT`.  The hash-bound
runner also launches core compilers with checkout inputs spelled relative to
the stable build directory, because Clang ASan global metadata preserves an
absolute source argument independently of those maps; the standalone harness
uses checkout-relative source, archive, and output arguments for the same
reason.  A two-checkout
proof requires every backend's archive and executable hashes to agree and
scans both retained text and binary bytes for either checkout root.  The
scan authenticates and examines the same one-read byte snapshot for each
artifact, so a concurrent path replacement cannot alter what was scanned.  The
final manifest retains a machine-checkable `comparison` attestation: the exact
peer-manifest artifact and SHA-256, a deterministic peer-evidence bundle, the
canonical matching-fingerprint SHA-256, the five build names, and the exact
normalized-text/archive/executable scan counts.  It stores no checkout path.
The bundle contains only the normalized configure/build/compile logs, CMake
caches, symbol scans, dependency records, and child result/stdout/stderr needed
for portable semantic replay; it deliberately does not copy generated archives
or executables.  Its gzip and tar metadata are canonical, exactly one gzip
member is permitted with no trailing padding or data, every tar member is
indexed by size and SHA-256, and validation rejects traversal, links, devices,
duplicates, noncanonical metadata, and bounded-size/decompression violations.
It also retains a canonical peer-attestation JSON artifact that binds both
retained peer artifacts, the peer's committed tooling and core closure, exact
programs, binary records, normalized records, run records, clean tooling Git
HEAD, ancestor core identity, one-time live validation, and root-byte scan.

The runner reads the original peer manifest once, captures all peer evidence
into a private immutable snapshot, and uses those exact bytes for portable and
live validation, equality comparison, hashing, scanning, bundling, and
attestation.  Subsequent changes to the peer manifest, logs, dependency files,
archives, or executables cannot change the result.  The trusted local validator
first checks the snapshot portably, then requires exact program-record
equality and complete peer/current archive and executable record equality
(path, byte count, and SHA-256), and only then performs one-time live tool and
output replay;
peer-controlled tool redirection is never executed.  Every non-tool artifact
path is checkout-relative, contains no
parent traversal, and resolves within the declared Git worktree even through
symlinks.  A `not-run` record is distinguishable from `pass` and is not
accepted by the retained-checkpoint test.  The
portable validator checks retained bytes, semantic logs, exact program records,
sanitizer equality, and member attribution without executing or requiring the
recorded tools or unretained peer build outputs.  `--live` explicitly requires
the current checkout's exact tools and independently rebuilt outputs and
byte-replays both current `nm` scans.  Each archive and executable is read and
authenticated once, then `nm` scans a private file made from that exact byte
snapshot rather than reopening the mutable build-output path.  Because all five current output hashes
and all program records must equal the peer records, this rechecks equivalent
bytes; it does not claim to recover or re-open the historical peer inodes after
the private generation snapshot has been discarded.
The retained validator likewise reads the peer manifest, semantic bundle, and
attestation only once each, authenticates the size and SHA-256 of that exact
in-memory byte string, and parses those same bytes.  A concurrent replacement
of any backing path after its read cannot redirect subsequent parsing.
Every normalized log, dependency file, symbol scan, and child result is also
interpreted from the exact authenticated byte snapshot rather than reopening
its backing path after hashing.  Live archives and executables receive the same
snapshot treatment before symbol replay.
All retained JSON parsing rejects duplicate object keys, non-standard numeric
constants, and finite-syntax numbers that overflow to non-finite values.  The
exact peer manifest and peer attestation must also
match the runner's canonical sorted, indented serialization, so alternative
whitespace or ambiguous-but-hash-consistent encodings cannot be substituted.

The provenance boundary is explicit.  Dependency reconstruction proves the
complete in-checkout project-header closure; absolute compiler and system
dependency roots are intentionally not a hermetic sysroot inventory.  Live
validation hashes and versions the recorded system tools and requires exact
peer/current program records, but it assumes those trusted system-tool paths
remain stable while they are invoked.  Concurrent privileged replacement of a
compiler, linker, `nm`, or `taskset` pathname is outside this evidence model.

Run-result validation fixes the exact build role, correctness/smoke kind,
sanitizer leak/undefined-behavior environment, selected runtime backend,
profile identity, allocation mode, OpenMP settings, and complete child schema.
The smoke cell must retain exactly seven positive finite samples per metric and
median/MAD summaries recomputed from those samples.

Final current-core attestation requires a v4 manifest whose A/B comparison
status is `pass`; a single-checkout v4 manifest is only an intermediate A-side
artifact.  The historical v3 checkpoint remains independently readable but is
not current-core evidence.

The Experiment-W Python suite passes 12 tests.  The strict C99 implementation
passes 217 direct checks and the differential matrix covers 5 golden vectors,
6,085 valid profile/count/field identities, 66 edge cases, 4,162 metadata cases,
87 malformed identifiers (including the six redundant-digest encodings), and
20,000 deterministic mutations.  Existing family-1/family-2 golden bytes remain
unchanged.

## Performance scope

The C++ harness reports exact and padded-low setup separately, then records raw
alternating encode and decode samples.  Its 12-cell matrix covers low,
balanced, and high transmitted rates; GF8/GF16; 64 B, 1 KiB, and 64 KiB; batch
one/eight; and loss counts three/eight.  Exact scratch is zero; padded scratch
is queried from the public API.

`results/smoke-nonauthoritative.json` is only a harness smoke pinned to the
first deterministically selected allowed CPU unless an explicit valid CPU is
provided.
It labels itself `non-authoritative-smoke` and retains seven raw samples for
each of four setup and four execution measurements.  No number from it is a
promotion result.  The final matrix must be rebuilt after the two-way SIMD and
comparison-adapter work is integrated, on a separately isolated physical core,
before any crossover conclusion.

## Reproduction

Rebuild the four normal backend archives, the Clang ASan+UBSan archive, all five
strict standalone executables, five independently pinned correctness runs, and
the affinity-selected smoke from a clean committed tooling revision.  Keep
builds and results under ignored paths so final clean-HEAD validation remains
meaningful.  Select the frozen integrated core explicitly; it must be an
ancestor of tooling `HEAD`, and all 22 core-closure paths must still match it
byte for byte.  Omitting
`--cpus` chooses the lowest five IDs from `sched_getaffinity`; the smoke uses the
first selected CPU.  An explicit `--cpus` still requires exactly five distinct
allowed IDs:

    test -z "$(git status --porcelain=v1 --untracked-files=all)"
    CORE_SHA=<frozen-integrated-core-commit>
    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/run_matrix.py \
      --core-git-sha "$CORE_SHA" --jobs-per-build 4 \
      --build-dir .research/leopard2/c7-final/build \
      --results-dir .research/leopard2/c7-final/results \
      --manifest .research/leopard2/c7-final/results/build-run-manifest.json

The separate authoritative runner consumes the resulting non-sanitized AVX2
archive and standalone executable.  Its schedule is not configurable: all 12
declared cells run in their frozen order, with seven retained samples for each
setup and execution metric.  The runner supports distinct clean tooling and
ancestor core commits, records exact Git trees plus runner, harness, archive,
executable, Python, and `taskset` hashes, proves the same fixed 22-file core
source closure used by the build runner, and checks the harness's embedded
source, archive, and core identities.  A coordinator must first reserve one allowed
physical core and its SMT sibling using canonical
`leopard2-cpu-reservation/v1` JSON.  The runner takes a nonblocking exclusive
lock for the whole campaign, moves itself to housekeeping CPUs, pins the child
to the timing CPU, and rejects evidence if `/proc/stat` reports any non-idle
work on the sibling.  Do not run other performance work during this window.

Run the synthetic validator and mutation suite before reserving a core:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/run_authoritative.py self-test
    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev -m unittest -v \
      experiments.leopard2.non_power_of_two.c7.test_run_authoritative

Then run one authoritative campaign, substituting an allowed SMT pair and the
exact ignored build paths produced above:

    TOOLING_SHA="$(git rev-parse HEAD)"
    CPU=<isolated-logical-cpu>
    SIBLING=<that-cpu-smt-sibling>
    RESERVATION=.research/leopard2/c7-final/cpu-reservation.json
    OUTPUT=.research/leopard2/c7-final/authoritative
    # The coordinator writes RESERVATION as canonical JSON without a newline:
    # {"benchmark_cpu":CPU,"nonce":"...","owner":"...",
    #  "reserved_sibling":SIBLING,"schema":"leopard2-cpu-reservation/v1",
    #  "status":"held"}
    python3 experiments/leopard2/non_power_of_two/c7/run_authoritative.py run \
      --source-root . --expected-tooling-commit "$TOOLING_SHA" \
      --expected-core-commit "$CORE_SHA" \
      --archive .research/leopard2/c7-final/build/core-avx2/liblibleopard.a \
      --executable .research/leopard2/c7-final/build/c7-avx2 \
      --reservation-file "$RESERVATION" --output "$OUTPUT" \
      --cpu "$CPU" --sibling "$SIBLING"

Success writes canonical `raw.json` and `manifest.json` plus the exact child
result and stdout/stderr bytes.  Any failure after output-directory creation
retains a signed `failure.json` and whatever child artifacts were available;
an existing output directory is never replaced.  Verification is portable and
does not reopen the original checkout, build, or reservation paths:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/run_authoritative.py verify \
      --manifest "$OUTPUT/manifest.json"

Regenerate independent algebra and validate retained evidence:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/algebra.py \
      --output experiments/leopard2/non_power_of_two/c7/results/algebra.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/validate_evidence.py \
      experiments/leopard2/non_power_of_two/c7/results/build-run-manifest.json

    # Optional, host-specific exact tool/output and nm replay:
    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/validate_evidence.py --live \
      experiments/leopard2/non_power_of_two/c7/results/build-run-manifest.json

The A/B gate runs the same clean tooling commit and core commit in a second
checkout, then gives
that checkout's manifest and root to the other runner with
`--compare-reproducibility-manifest` and
`--compare-reproducibility-root`.  Comparison fails when either Git identity,
any tooling record, any backend hash, any program record, or a checkout-root
byte differs, and the resulting local manifest records
the path-free comparison attestation.  It writes the exact captured bytes to
`results/peer-manifest.json`, the portable semantic inputs to
`results/peer-evidence-bundle.tar.gz`, and their binding report to
`results/peer-reproducibility-attestation.json` (with `results` interpreted as
the selected `--results-dir`).  These generated proof artifacts belong in the
ignored research results area unless repository policy explicitly approves a
small retained checkpoint; generated archives and copied executables are never
placed in this bundle.  The peer manifest may be the initial
single-checkout input; only the comparison-producing manifest is retained as
final evidence.

Run the Experiment-W identity gates from its directory:

    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v test_code_identity.py
    PYTHONDONTWRITEBYTECODE=1 python3 test_code_identity_c.py --cc gcc
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
      PYTHONDONTWRITEBYTECODE=1 python3 test_code_identity_c.py \
      --cc gcc --sanitizers

The runner performs that strict source/core/archive binding and records it in
`results/build-run-manifest.json`.  Sanitizer builds additionally define the
allocation-interposition opt-out, the sanitizer mode, and the compile-time
ASan+UBSan requirement because custom global allocation interposition would
mask sanitizer allocation diagnostics.  The standalone C7 harness and runner
are currently POSIX/Linux-only: they use `posix_memalign`, `sched_getaffinity`,
and `taskset`.  The field mathematics and `L2ID` byte format are not
Linux-specific.  Do not use historical C6 archives for final performance
evidence.

## Disposition

Keep family 3 and its implementation experimental/default-off.  C7 resolves
the coordinate-map and correctness questions, but production promotion remains
blocked on final integrated-core timing, a region dispatcher (C10), broader API
review, and exact-specialized decoder work (C9).  Family 4 remains reserved for
C8 exact-high work.
