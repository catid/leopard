# Leopard2 C5 binary dyadic-block decomposition

Status: bounded scalar checkpoint complete with a non-promotion disposition for
the naive local-to-global join-plus-parent route studied here.  Do not close the
parent C5 Bead:
its broader acceptance still requires representative GF16 and measured
byte/batch/reuse evidence for any surviving optimized factorization.  The
coordinator can retain this result as a completed scalar-checkpoint child.

This checkpoint does not justify adding a separate binary-block encoder or
wire profile.  The parent-preserving identity is correct but is already a
scheduling option within C1/C2.  The explicit exact-profile join adds
cross-block work; a different factored or pruned-output construction remains
inconclusive rather than disproved.

The executable model is
`experiments/leopard2/non_power_of_two/c5/dyadic_blocks.py`.  It is an algebra,
correctness, and logical-operation-count program.  It does not benchmark,
modify production code, register a CMake target, or select a runtime path.

## Two code identities, kept separate

Let `omega_i` be the field element at additive coordinate `i`; addition of
coordinates is XOR.  For any positive `q`, the canonical binary decomposition
is

    q = b_0 + b_1 + ... + b_s,  b_j = 2^a_j,  a_0 > ... > a_s.

Block `j` starts at `o_j = sum_{h<j} b_h`.  Every `o_j` is divisible by `b_j`,
so `[o_j,o_j+b_j)` is the aligned additive coset
`omega_o_j + V_log2(b_j)`.

The experiment studies two different meanings of these blocks:

1. **Parent-wire-compatible computation.**  The first `q` coefficients of an
   existing power-of-two LCH parent are nonzero.  The missing coefficient
   suffix remains mathematically zero and all outputs remain the existing
   parent outputs.  Splitting the work does not define a new code.
2. **Exact prefix-evaluation code.**  Exactly `q` source values at points
   `omega_0,...,omega_(q-1)` are interpolated, then evaluated at later points.
   This is an exact systematic RS code, but it is a new profile unless its
   generator matrix is separately proved equal to an existing profile.  C5
   makes no legacy-compatibility claim for it.

Successful recovery is not wire equivalence.  The result artifact carries
these two identities explicitly.

## Parent-preserving LCH block identity

Write the normalized LCH basis as

    X_i(x) = product over set bits k of i: s_k(x) / s_k(omega_(2^k)),

where `s_k(x) = product over alpha in V_k of (x + alpha)`.  If a block has
offset `o` and size `b`, the low `log2(b)` bits of `o` are zero.  Therefore the
set bits of `o` and any `t < b` are disjoint, and normalization factors also
separate:

    X_(o+t)(x) = X_o(x) X_t(x).

`X_o` is constant on every aligned coset of `V_log2(b)`.  A block contribution
can consequently be computed by a normal shifted `b`-point LCH transform,
scaled once per output by `X_o`, then XOR-accumulated with other blocks.  A
coset whose factor is zero is skipped.  The scalar oracle constructs every
`X_i` independently as a monomial polynomial and checks this identity entry by
entry against the padded parent evaluation matrix.

The retained operation counts are exact for that abstract schedule:

- `kernel_butterfly_equivalents` counts complete radix-2 butterflies in active
  shifted block transforms;
- `block_factor_fixed_multiplications` excludes factors zero and one;
- `cross_block_accumulation_xors` assigns the first contribution and counts
  XOR only for later contributions;
- loads and stores count two logical operands/results per butterfly; and
- scratch is the largest local block, excluding caller-owned input and output.

They do not include packing, scatter, schedule setup, cache traffic, or the
zero/one skew specialization inside a butterfly.

Representative legacy-GF8 counts are:

| q | Blocks | Parent | Block butterflies / padded | Factor muls | Join XORs | Local scratch / padded |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 33 | 32+1 | 64 | 160 / 192 | 0 | 32 | 32 / 64 |
| 65 | 64+1 | 128 | 384 / 448 | 0 | 64 | 64 / 128 |
| 129 | 128+1 | 256 | 896 / 1,024 | 0 | 128 | 128 / 256 |
| 191 | 128+32+16+8+4+2+1 | 256 | 1,856 / 1,024 | 640 | 768 | 128 / 256 |
| 255 | 128+64+32+16+8+4+2+1 | 256 | 2,240 / 1,024 | 704 | 896 | 128 / 256 |

Across all 255 GF8 prefix lengths, the block schedule uses fewer raw
butterflies in only 19 cases, the same number in 18, and more in 218.  Even the
favorable `2^a+1` cases add an output accumulation pass.  C1/C2 dependency
pruning represents the same parent map without committing to this coarse block
schedule, so the identity does not justify another production implementation.

## Exact-code joins in three bases

For a basis `B`, let `E_B` be its `q`-point source evaluation matrix.  Let
`D_B` be the block diagonal matrix that maps independently computed local block
coefficients to source values.  The unique local-to-global join into basis `B` is

    J_B = inverse(E_B) D_B.

For LCH, each diagonal block is the existing shifted-coset evaluation of the
global `X_0,...,X_(b-1)` basis.  For Newton it uses roots from that block (a
translated local Newton basis), and for Lagrange it is the identity on values.

For parity evaluation matrix `O_B`, the complete local-coordinate encoder is

    O_B J_B.

The independent oracle constructs the direct Lagrange evaluation matrix `L`
from products of point differences and verifies

    O_B J_B = L D_B

entry by entry.  LCH bases are built from monomial subspace polynomials;
Newton bases are built from their defining prefix products; Lagrange uses
source values directly.  Thus the candidate does not verify itself with a
second execution of the same transform.

The same-basis join is not uniformly dense, but it does contain cross-block
work for LCH and Newton.  Lagrange has an identity join because its coordinates
are the source values, after which its parity map is dense.  Representative
GF8 results below use 17 parity outputs:

| q | Basis | Join nonzeros / entries | Cross-block nonzeros / entries | Parity nonzeros / entries |
| ---: | --- | ---: | ---: | ---: |
| 15 | LCH | 49 / 225 | 30 / 140 | 255 / 255 |
| 15 | Newton | 51 / 225 | 36 / 140 | 255 / 255 |
| 15 | Lagrange | 15 / 225 | 0 / 140 | 255 / 255 |
| 33 | LCH | 65 / 1,089 | 32 / 64 | 561 / 561 |
| 33 | Newton | 65 / 1,089 | 32 / 64 | 561 / 561 |
| 63 | LCH | 321 / 3,969 | 232 / 2,604 | 1,071 / 1,071 |
| 63 | Newton | 421 / 3,969 | 358 / 2,604 | 1,071 / 1,071 |
| 65 | LCH | 129 / 4,225 | 64 / 128 | 1,105 / 1,105 |
| 65 | Newton | 129 / 4,225 | 64 / 128 | 1,105 / 1,105 |

The JSON reports `q` slots for a straightforward immutable-input,
out-of-place join and the largest local block separately.  That is an
implementation accounting point, not a scratch lower bound: Lagrange's
identity join needs no coefficient copy, and a triangular join may permit
additional overwrite/reuse.

Partitioning the requested parity interval into maximal aligned dyadic cosets
does not decouple the input blocks.  Every one of 189 GF(2^4) input/output block
pairs and all 291 pairs in the selected GF8 matrix cases has at least one
nonzero coefficient.  The effective local-coordinate maps are already close
to dense in the larger sampled cases:

| q | Basis | Effective nonzeros / entries | Coupled block pairs / pairs |
| ---: | --- | ---: | ---: |
| 15 | LCH / Newton / Lagrange | 244 / 255; 244 / 255; 255 / 255 | 8 / 8 each |
| 33 | LCH / Newton / Lagrange | 485 / 561; 495 / 561; 561 / 561 | 10 / 10 each |
| 63 | LCH / Newton / Lagrange | 1,012 / 1,071; 1,014 / 1,071; 1,071 / 1,071 | 12 / 12 each |
| 65 | LCH / Newton / Lagrange | 965 / 1,105; 969 / 1,105; 1,105 / 1,105 | 10 / 10 each |

The parity density is structural, not a sample accident.  For an output point
`beta` outside the source prefix:

- every Lagrange coefficient is a product of nonzero point differences;
- every Newton basis value is a product of `beta + omega_i`, all nonzero; and
- a set bit `k` of an LCH basis index `j < q` satisfies `2^k <= j < q`, while
  every root of `s_k` has coordinate below `2^k`; hence `s_k(beta)` is nonzero.

Therefore every global-basis parity row is fully dense.  Joining local blocks
and then running a padded transform adds the join to existing padded work.
Fusing all the way from source values to parity reproduces the ordinary dense
direct evaluation matrix, which is already a separate small-code dispatcher
candidate rather than a C5 transform.

## Correctness evidence

The deterministic artifact is
`experiments/leopard2/non_power_of_two/c5/results/self_test.json`.

| Check | Count |
| --- | ---: |
| GF(2^4) prefix geometries | 15 |
| GF(2^4) LCH/Newton/Lagrange join cases | 45 |
| GF(2^4) effective join entries compared with direct algebra | 2,040 |
| GF(2^4) coupled input/output dyadic-block pairs | 189 / 189 |
| GF(2^4) parent block-map entries compared | 1,585 |
| GF(2^4) systematic generators | 15 |
| GF(2^4) full-length K-coordinate subsets rank-tested | 65,534 |
| Public GF(2^4) `(K,R)` geometries implied by those full-length tests | 120 |
| Legacy-GF8 prefix geometries | 255 |
| Legacy-GF8 parent block-map entries compared | 7,146,545 |
| Legacy-GF8 detailed join cases, at ten boundary-focused q values | 30 |
| Legacy-GF8 effective join entries compared with direct algebra | 12,648 |
| Legacy-GF8 coupled input/output dyadic-block pairs | 291 / 291 |
| Legacy-GF8 LCH parity coefficients checked nonzero | 2,796,160 |
| Parity coefficients proved nonzero across all three bases | 8,388,480 |

For GF(2^4), the test constructs the full 16-coordinate systematic generator
for every `1 <= K < 16` and Gaussian-eliminates every subset of `K` transmitted
coordinates.  Any shorter recovery prefix is a row subset of that full code,
so this exhausts all 120 valid `(K,R)` geometries without repeating equivalent
rank tests.  It does not enumerate all `16^K` messages; exhaustive linear-map
and rank checks prove the corresponding claims for every message.

GF8 uses Leopard's polynomial `0x11d`, legacy Cantor basis, symbol values, and
coordinate order.  Detailed matrix inversions are deliberately bounded to
`q = 3,5,7,9,15,17,31,33,63,65`; the parent identity and parity-density checks
cover every `q` from 1 through 255.

## Disposition

Do not promote the naive explicit local-to-global join followed by parent
evaluation.  It adds cross-block setup/execution before the parent work, and
fusing all the way from source values yields the direct evaluator already
covered by direct/exact-profile work.  Matrix density and all-to-all block
coupling are not runtime lower bounds: they do not exclude a different fast
factorization or a pruned-output schedule.  Such a construction was not
implemented here and remains an open parent-C5/C7/C8 question.

Do not promote a separate parent-wire block executor.  Its algebra is correct,
but 218 of 255 GF8 prefixes already have more block butterflies before factor
multiplications and accumulation are charged.  The few favorable geometries
are better represented by the more general C1/C2 scheduled pruning work.

No timing was run, intentionally: this scalar checkpoint ran while other host
validation was active, and the studied explicit route did not justify an
optimized candidate.  The density proof uses only distinct
prefix points and active-subspace root sets and therefore extends structurally
to GF16, but that proof does not satisfy the parent Bead's executable GF16 or
measured byte/batch/reuse gates.  Those gates remain open.

Only promote a surviving block factorization if a C1/C2 production schedule
subsuming it passes all wire/correctness gates and demonstrates at least a 10%
end-to-end improvement after setup and memory traffic in the parent C5's GF16,
byte, batch, and reuse cells.  C10 can then incorporate its crossover.  An
exact coordinate profile, if pursued for field-boundary rescue, belongs to
C6/C7/C8 and needs a new serialized code identity and its own compatibility
evidence.

## Reproduction

From the repository root:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c5/dyadic_blocks.py \
        --output /tmp/c5-self-test.json

    cmp /tmp/c5-self-test.json \
        experiments/leopard2/non_power_of_two/c5/results/self_test.json

    sha256sum \
        experiments/leopard2/non_power_of_two/c5/dyadic_blocks.py \
        experiments/leopard2/non_power_of_two/c5/results/self_test.json

The program uses only the Python standard library.  Retained hashes:

    de4191028948b1dc2b53032dd22a978282347a62711c3bfb29a56c5d3b4808b4  dyadic_blocks.py
    204230d8008f398c7a9369e8ae2910f0852356d9bc16f4418288ac5e59fd7301  self_test.json

Mathematical lineage and the distinction between a truncated parent transform
and a new exact code follow R15 (Coxon's truncated additive transforms) and the
profile rules in `docs/leopard2_math_and_sources.md`.  C5 is a clean scalar
reimplementation and does not copy code from a research repository.
