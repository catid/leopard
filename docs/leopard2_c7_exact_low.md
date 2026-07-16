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
CPU-0 smoke is explicitly non-authoritative, is rebuilt from the same committed
source as the correctness matrix, and cannot support a promotion claim.
Authoritative crossover timing must be rebuilt from the final integrated
production commit.

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

Profile family 4 is reserved for future exact-high work and is rejected by the
version-1 readers.  Exact-high must not reinterpret family 3.  Coordinate-,
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
rejected before a fixed-multiplier helper runs.

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
| Hot-path allocations | 0 |
| Deterministic digest | `0xec4179e9f2776a58` |

Lengths are GF8 `1,2,3,7,31,64,65,257` bytes and GF16
`2,4,6,14,62,64,66,130,514` physical bytes.  Buffers are unaligned and guarded.
Combined ASan+UBSan with leak detection passes the same complete matrix.  The
sanitized standalone source has a compile-time feature gate that fails unless
both instruments are active; the retained build/run manifest also records the
exact compiler, flags, linked archive, executable, source closure, sanitizer
symbol scan, run environment, child affinity, and artifact hashes.

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

`results/smoke-nonauthoritative.json` is only a harness smoke pinned to CPU 0.
It labels itself `non-authoritative-smoke` and retains seven raw samples for
each of four setup and four execution measurements.  No number from it is a
promotion result.  The final matrix must be rebuilt after the two-way SIMD and
comparison-adapter work is integrated, on a separately isolated physical core,
before any crossover conclusion.

## Reproduction

Rebuild the four normal backend archives, the Clang ASan+UBSan archive, all five
strict standalone executables, five independently pinned correctness runs, and
the CPU-0 smoke from a committed source revision:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/run_matrix.py \
      --core-git-sha "$(git rev-parse HEAD)" --cpus 0,1,2,3,4 \
      --smoke-cpu 0 --jobs-per-build 4

Regenerate independent algebra and validate retained evidence:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/algebra.py \
      --output experiments/leopard2/non_power_of_two/c7/results/algebra.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
      experiments/leopard2/non_power_of_two/c7/test_checkpoint.py -v

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
