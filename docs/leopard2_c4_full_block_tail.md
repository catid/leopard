# Leopard2 C4: full dyadic block plus exact or dense tail

## Result

C4 does not promote production code or a wire profile.

The isolated scalar study validates a full `2^a` head block followed by a
`d`-symbol tail for `q = 2^a + d`, with `0 < d < 2^a`.  It keeps two different
maps separate:

- The parent-compatible map fixes the remaining `2^a-d` evaluations to zero
  and computes the existing `2^(a+1)`-point Leopard inverse.  Direct dense,
  smaller-padded, and recursive-prefix tail executors match that parent.
- The exact-q map interpolates the unique degree-less-than-q polynomial from
  exactly q prefix points.  Its dense Schur-complement and Newton-tail
  executors agree, but 8,711 tested vectors distinguish it from the padded
  parent.  It is new-profile research, not a legacy optimization.

The implementation and retained evidence are under
`experiments/leopard2/non_power_of_two/c4/`.  Nothing there is referenced by
CMake, installed headers, the public API, or runtime dispatch.

The deterministic operation model predicts useful regions, particularly one
or a few symbols just beyond a large power of two.  Pinned scalar Python also
shows execution savings for several tiny-tail cases.  Neither establishes an
end-to-end C++/SIMD codec win, memory traffic, instruction-cache cost, or a
neighbor-safe dispatcher.  The disposition is therefore:

- **Promoted:** nothing.
- **Smaller padded plus dyadic lifts:** production candidate, inconclusive
  until a C++ kernel calls the existing fused transforms.
- **Recursive prefix tail:** production candidate, inconclusive until its
  boundary operations are fused and measured.  C2's general dictionary-based
  executor remains a negative performance result; this prefix specialization
  does not reverse that result without a C++ checkpoint.
- **Direct dense parent tail:** inconclusive for very small d and high plan
  reuse.  Its `2^a*d` table and setup cost kill broad use.
- **Exact direct/Newton tail:** mathematically correct new-profile research.
  It cannot be selected under a legacy code identifier.

## Parent-wire-compatible factorization

Let `b = 2^a`, `N = 2b`, and split the supplied evaluations into `y_H` of
length b and `y_T` of length d.  Let `B_b(s)` be the b-point LCH evaluation map
on the aligned coset beginning at shift s.  The padded parent inverse first
forms

    u = inverse(B_b(s)) * y_H
    v = inverse(B_b(s+b)) * [y_T, 0_(b-d)]

and then applies the ordinary N-point inverse root butterfly to `(u,v)`.  For
the root skew factor m, each coefficient pair is

    c_0 = u + m * (u + v)
    c_1 = u + v

with the existing zero/one specializations.  The experiment computes v in
three ways:

1. **Direct dense tail:** retain the first d columns of
   `inverse(B_b(s+b))` and multiply that b-by-d matrix by `y_T`.
2. **Smaller padded tail:** use `e = ceil_pow2(d)`, run one regular e-point
   inverse on `[y_T,0_(e-d)]`, then lift the result through the remaining
   dyadic inverse layers while the right child is known zero.
3. **Recursive prefix tail:** recursively execute only the child containing
   active prefix values, call a complete transform for every full child, and
   use a known-zero unary lift at ragged ancestors.

All three produce the same b coefficients v.  The direct matrix is constructed
from independently evaluated normalized LCH polynomials and Gaussian inversion
through b=32.  Larger correctness cells build unit columns from the padded
reference and then compare all three executors with the complete parent; those
larger direct-column cells are differential evidence, not an independent
algebra oracle.

This factorization removes computation only.  The zero suffix and all N parent
coordinates remain part of the code definition.

## Exact-q factorization

For the new exact map, partition the q-by-q prefix evaluation matrix as

    A = [ A00  A01 ]
        [ A10  A11 ]

where A00 is b-by-b.  Define

    E = inverse(A00) * A01
    S = A11 + A10 * E

where addition is subtraction in characteristic two.  Execution is

    h = inverse(A00) * y_H
    r = y_T + A10 * h
    t = inverse(S) * r
    c_H = h + E * t

and the exact coefficients are `(c_H,t)`.  The b-point computation of h is a
normal Leopard inverse.  Dense tail execution stores `inverse(S)`.

On the zero coset, A01 and E are zero, and S is exactly the d-point prefix
evaluation map of the lower LCH basis on the tail coset beginning at b.  The
Newton-tail executor computes divided differences on those d points, expands
the Newton polynomials, and converts the result to the same lower LCH basis.
The implementation asserts the matrix identity before use and compares every
Newton result with the dense Schur solve.

This transparent Newton implementation is not Coxon's fast recursive basis
conversion and is not Tang-Han's epsilon inverse.  It is retained as an
independent exact-tail oracle.  The mathematical distinction and source scope
follow the C2/C3 notes and the bibliography in
`docs/leopard2_math_and_sources.md`, especially:

- Nicholas Coxon, *Fast Transforms over Finite Fields of Characteristic Two*,
  https://arxiv.org/pdf/1807.07785
- Nianqi Tang and Yunghsiang S. Han, *New Decoding of Reed-Solomon Codes Based
  on FFT and Modular Approach*, https://arxiv.org/pdf/2207.11079

## Correctness and MDS evidence

The retained 32-worker run used every CPU allowed to this checkout, CPUs 0-31.
A one-worker deterministic recomputation matched the complete analysis object.
The artifact carries the exact C4 source hash and the C3 oracle hash; `verify`
checks both, validates an integrity digest, recomputes all evidence, and
compares canonical JSON.

| Evidence | Count |
| --- | ---: |
| Correctness jobs | 65 |
| Input vectors | 9,672 |
| Field/result comparisons | 435,714 |
| Parent-versus-exact difference witnesses | 8,711 |
| Exhaustive GF(2^4) received-coordinate subsets | 50,708 |
| Cost geometries | 402 |
| Byte/batch/reuse/profile cells | 38,592 |

Correctness coverage includes:

- GF(2^4) at every non-power q from 3 through 15, on the zero and final
  aligned parent cosets where distinct;
- legacy-coordinate GF8 at q=3,5,6,7,9,15,17,31,33,63,65,127,129,255;
- legacy-coordinate GF16 at q=3,5,7,9,15,17,31,33,63,65,127,129;
- zero, basis, structured, and SHA-256-seeded random vectors, plus every
  `16^3` input at q=3 in GF(2^4);
- direct normalized-LCH evaluation matrices, Gaussian inverse oracles,
  padded butterfly inverses, dense Schur solves, and Newton-tail conversion;
- all three parent tail executors at shifted cosets; and
- systematic/MDS rank of the exact prefix-point code for every q-coordinate
  subset of all 16 GF(2^4) evaluation points, for every non-power q=3..15.

Every listed GF8/GF16 geometry covers all three parent-compatible tail
executors.  The independently inverted dense exact split is additionally run
only when b<=32 and d<=17; the Newton exact-tail comparison is the zero-shift
subset of those cells.  The 8,711 difference count and exact-oracle claims
refer to that explicitly gated subset, not every large parent geometry.

The MDS check enumerates coordinate subsets, not all field-valued messages.
For a linear code, full rank of every q-row generator submatrix is the direct
MDS condition.

## Cost sweep and offline table

The sweep covers every d for GF8 head blocks b=2,4,8,16,32,64,128: all 247
non-power q values through 255.  GF16 covers every d through b=64 and the
representative tails 1,2,3,b/8,b/4,b/2,b-1 at b=128,256,512,1024,2048.  Each
geometry is evaluated at shard sizes 64, 1,024, 65,536, and 1,048,576 bytes;
batches 1,8,64; and plan reuse 1,8,64,1024.

The score is the same deliberately simple heuristic used by C3:

    3 * fixed_multiplications + XORs + loads + 2 * stores
      + 12 * inversions

Setup is counted separately and divided by reuse.  Execution is multiplied by
the symbol count and batch.  Complete transform counts use the actual active
coset skew factors; the dense model conservatively treats every matrix entry
as a nontrivial multiplier.  Shared field/basis initialization is excluded.
This is an operation and logical-access model, not a prediction of cache or
SIMD behavior.

Direct-parent table setup assumes d unit-column parent transforms rather than
the slower Gaussian correctness oracle.  Smaller-padded setup has no private
table.  Recursive setup uses a compact descriptor-store proxy, not a measured
C++ scheduler build.  Those assumptions make setup/reuse comparisons useful
for screening but insufficient for promotion.

The compact offline choice table is `analysis.modeled_choices` in
`result.json`.  Every cell is

    [bytes, batch, reuse, parent_method_id, parent_gain,
     exact_method_id, exact_gain]

and method IDs index `analysis.choice_method_ids`.  Including the padded parent
as a parent candidate makes this table deterministic and non-regressing within
its own score.  It is research evidence, not a production dispatch table.

| Result | Cell |
| --- | --- |
| Best parent-compatible modeled gain | 1.758844225x, GF16 q=2049/N=4096, d=1, 64 KiB, batch 64, reuse 1024, direct dense tail |
| Best exact-profile modeled gain versus padded computation | 1.983767781x, same geometry, batch 8, exact dense tail |
| Parent cells at or above 1.10x | 13,536 / 19,296 |
| Exact cells at or above 1.10x | 2,446 / 19,296 |

Parent winner counts are 13,536 recursive-prefix, 4,320 smaller-padded, 672
direct-dense, and 768 full-padded cells.  Exact winner counts are 19,112 dense
and 184 Newton cells.  Fewer modeled operations do not satisfy the promotion
gate without a regular C++ implementation and whole-codec measurement.

## Pinned scalar timing

Timing was intentionally isolated on allowed CPU 16 with
`OMP_NUM_THREADS=1`.  The host is an AMD Ryzen 9 9950X3D with 32 logical/16
physical CPUs and one NUMA node.  Python 3.12.3 ran on Linux 6.8.0-134.  Values
are median nanoseconds for one field-symbol position, with median absolute
deviation in parentheses.  They do not process shard byte arrays and are not
throughput claims.

| Field / q | Padded parent | Direct parent tail | Smaller padded | Recursive prefix | Exact dense | Exact Newton |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GF8 / 3 | 2,276 (5) | 2,129 (6) | 2,240 (4) | 2,154 (5) | 2,230 (13) | 2,814 (5) |
| GF8 / 5 | 5,142 (20) | 3,936 (10) | 4,243 (12) | 4,221 (25) | 3,784 (9) | 4,322 (44) |
| GF8 / 9 | 11,386 (54) | 7,716 (24) | 7,980 (49) | 8,056 (34) | 6,998 (54) | 7,423 (29) |
| GF8 / 17 | 24,431 (103) | 15,652 (39) | 15,410 (30) | 15,545 (52) | 13,839 (23) | 13,964 (41) |
| GF8 / 33 | 52,217 (120) | 32,053 (57) | 30,685 (63) | 30,796 (230) | 28,288 (69) | 27,911 (104) |
| GF16 / 3 | 4,398 (17) | 4,251 (29) | 4,377 (12) | 4,295 (19) | 4,185 (21) | 6,736 (17) |
| GF16 / 5 | 14,686 (66) | 12,564 (45) | 13,123 (68) | 13,100 (31) | 12,034 (41) | 14,512 (121) |
| GF16 / 9 | 39,087 (132) | 32,357 (30) | 32,928 (87) | 32,844 (114) | 31,952 (92) | 34,296 (95) |
| GF16 / 17 | 100,065 (168) | 82,547 (106) | 80,168 (124) | 80,110 (139) | 77,748 (284) | 80,438 (109) |

Plan setup is material.  For direct parent tails it rises from 7,611 ns at
GF8 q=3 to 2,814,235 ns at q=33, and from 59,400 ns at GF16 q=3 to 7,539,507
ns at q=17.  Exact dense/Newton setup ranges are retained per case in the
timing JSON.  Smaller-padded and recursive executors use no separately built
dense plan in this scalar model; production recursion would still need a
precompiled schedule whose setup and storage must be measured.

These scalar results identify a C++ checkpoint worth trying.  They cannot
promote a kernel: Python list and field-operation overhead does not model
Leopard's fused byte-wide SIMD transforms, and only q=b+1 was timed.  Neighbor
tails, shard byte loops, setup amortization, scratch, and end-to-end encode and
decode remain required.

## Artifact validation and reproduction

From the repository root:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py self-test

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py analyze \
      --workers 32 \
      --output experiments/leopard2/non_power_of_two/c4/result.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py verify \
      experiments/leopard2/non_power_of_two/c4/result.json --workers 1

    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
      PYTHONMALLOC=debug taskset -c 16 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py benchmark \
      --field gf8 --counts 3 5 9 17 33 --samples 11 --iterations 64 \
      --output experiments/leopard2/non_power_of_two/c4/timing_gf8.json

    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
      PYTHONMALLOC=debug taskset -c 16 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py benchmark \
      --field gf16 --counts 3 5 9 17 --samples 9 --iterations 24 \
      --output experiments/leopard2/non_power_of_two/c4/timing_gf16.json

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py verify-timing \
      experiments/leopard2/non_power_of_two/c4/timing_gf8.json

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c4/study.py verify-timing \
      experiments/leopard2/non_power_of_two/c4/timing_gf16.json

The verifier rejects a modified integrity payload before recomputation and
rejects an artifact whose C4 source or C3 dependency hash no longer matches.
The self-test also mutates a sealed in-memory object and requires rejection.

Final SHA-256 values:

| File | SHA-256 |
| --- | --- |
| `study.py` | `721d5ec70a3d8b8e54b0ab5a11062fe52efc398164e4d2ca0550db1ec41c96bc` |
| `result.json` | `be5acb2861c2b54fd0efcd9c76fb1911d3b19e8d947777ddad1edda22e5b1cfd` |
| `timing_gf8.json` | `47a5484e3cb600bbae81db176ca9ad8680ff7238a038d9198db2cce752d6a3d3` |
| `timing_gf16.json` | `67522cb743b5796a4179f6dbc5b03afab66db9eaffc0328c7bb7c81663be7ee4` |
| C3 dependency `study.py` | `c0030150633f133d0b04dda762096e3d6164173e99a60888b214476cb2d1b795` |

## Highest-value next step

Build one standalone C++ prefix-tail checkpoint that calls the existing fused
GF8/GF16 kernels for complete blocks, emits compact unary-lift runs for the
ragged boundary, and compares byte-for-byte with the padded production
transform.  Measure q=b+d for d=1,2,3,b/8,b/4,b/2,b-1 across shard bytes and
plan reuse, with direct dense enabled only under a table-size cap.  Promote
only if the whole-codec gain remains at least 10%, setup is amortizable, and a
deterministic dispatcher avoids neighboring regressions.  Exact-q work should
continue only under the separately versioned C7/C8 profile gates.
