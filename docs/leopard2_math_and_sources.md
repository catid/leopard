# Leopard2 mathematics, coordinates, and sources

Status: implementation specification, 2026-07-16.  A formula in this document is
not enabled in a production decoder until the validation gates in "Independent
oracles" pass.  In particular, the full-field decoder formulas in R10 require the
active-parent factors derived below when a parent of size `2^n` is embedded in a
larger field `GF(2^m)`.

## Terms and code identity

Leopard2 uses the following names consistently:

- `K`: transmitted systematic/original shard count.
- `R`: transmitted recovery/parity shard count.
- `N`: power-of-two parent-code length.
- `P = ceil_pow2(K)`: low-profile transform block length.
- `T = ceil_pow2(R)`: high-profile transform block length.
- shortening: a systematic parent coordinate fixed to zero and not transmitted.
- puncturing: a parity parent coordinate not transmitted.
- erasure: a missing coordinate whose position is known.

A wire profile fixes the field, field representation, evaluation-coordinate order,
parent construction, shortening set, puncturing set, and parity order.  A SIMD
backend, thread count, transform schedule, or calibrated dispatch threshold is not
part of the wire profile and may not alter parity bytes.

## Field and coordinate representation

The legacy profiles retain the repository's exact representation.

| Property | GF8 | GF16 |
| --- | --- | --- |
| irreducible polynomial | `0x11d` | `0x1002d` |
| stored symbol | Cantor-basis coordinate bits | Cantor-basis coordinate bits |
| addition | bitwise XOR | bitwise XOR |
| coordinate `omega_i` | stored value `i` | stored value `i` |

The GF8 Cantor basis, expressed in the polynomial basis used while constructing the
log table, is:

    1, 214, 152, 146, 86, 200, 88, 230

The GF16 basis is:

    0001, acca, 3c0e, 163e, c582, ed2e, 914c, 4012,
    6c98, 10d8, 6a72, b900, fdb8, fb34, ff38, 991e

For binary `i = sum(i_j 2^j)`, the mathematical evaluation point is
`omega_i = sum(i_j v_j)`.  Because the API representation consists of those same
Cantor coordinates, `omega_i + omega_j` has stored value `i XOR j`.  This is the
coordinate identity used by sparse and incremental locator experiments; it does
not license changing the polynomial-basis constants or byte representation.

## LCH basis on an active parent

Let the ambient field be `F = GF(2^m)` and let `0 <= j <= n <= m`.  Define the
active additive subspaces

    V_j = { omega_0, ..., omega_(2^j - 1) }
    s_j(x) = product over a in V_j of (x - a)
    c_j = s_j'(x) = product over nonzero a in V_j of a

The `s_j` are linearized polynomials, so `s_j(x+y)=s_j(x)+s_j(y)`, and their
formal derivatives are the nonzero constants `c_j`.  The normalized LCH basis is

    Xbar_i(x) = product_j s_j(x)^(i_j) / p_i
    p_i       = product_j s_j(v_j)^(i_j)

for the binary digits `i_j` of `i`.  Transform length `2^q`, `q <= n`, evaluates
a polynomial of degree below `2^q` on the coset `V_q + beta`.  A shift/coset
argument is therefore a field value (stored in Cantor coordinates), not a block
ordinal that may be renumbered.

R16 and R17 define this basis and the additive FFT.  R10 equations (1)-(16) give
the normalization and the block FFT/IFFT identities used here.  Existing Leopard
tables encode the quotients involving `s_j(beta)/s_j(v_j)` as logarithmic skew
factors.  New scalar tests independently generate `s_j`, `c_j`, and `p_i` from
field multiplication before comparing those tables.

### Production skew, normalizer, and derivative constants

The production transform stores one logarithmic multiplier for each logical
butterfly index `ell`.  Let `j` be the number of trailing one bits of `ell`, let
`h=2^j`, and let `beta=ell-(h-1)`.  The butterfly joins two `h`-coefficient
halves on the coset beginning at `beta`, so its field multiplier is

    M_ell = s_j(beta) / s_j(v_j)
          = product_(a in V_j) (beta-a)
            / product_(a in V_j) (v_j-a).

This identifies every `FFTSkew` entry without reproducing its table generator.
The oracle evaluates both products with direct polynomial-basis field
multiplication and compares all 255 GF8 and 65,535 GF16 entries.

The basis normalizers used by active-parent setup are generated from

    b_j = s_j(v_j),
    p_i = product_(j: i_j=1) b_j.

Production caches the `log(b_j)` constants and forms a requested `log(p_i)` by
summing only the entries selected by `i`.  The low and high prepared factors are

    coefL(block) = c_k / s_k(omega_(block*P)),

    coefH(x) = c_n / [p_(N-T) s_t(x)].

Every `p_i` in GF8 and GF16, every active power-of-two subspace/coset constant,
and every legal power-of-two low/high prepared-factor array are compared with
direct field products.  This evidence is independent of codec recovery
succeeding.

The product rule gives the normalized-basis derivative identity

    derivative(Xbar_i)
      = sum_(j: i_j=1) [c_j / s_j(v_j)] Xbar_(i-2^j).

For Leopard's declared Cantor bases these ratios permit the existing triangular
XOR implementation.  The byte kernel accumulates `f'` into `f` rather than
clearing the old coefficients.  Decoder output is evaluated only at locator
roots, where `f=0`, so `(f+f')(e)=f'(e)`.  Tests compare that exact accumulation
with the direct product-rule formula.

Shift-zero kernels historically formed a pointer one element before the
`FFTSkew` array.  Logical indexing is now preserved with a real sentinel storage
element, obtaining the same arithmetic and parity bytes without an out-of-bounds
C++ pointer.

### Why ambient `m` cannot simply become active `n`

R10 specializes its decoders to the full field `V_m = GF(2^m)`.  In that case

    s_m(x) = x^(2^m) - x
    c_m = s_m'(x) = 1

in characteristic two.  Leopard frequently uses an active parent `V_n` inside
GF16 with `n < 16`.  The correct congruence is then modulo `s_n(x)`, not modulo
`x^(2^m)-x`, and in general `c_n != 1`.  The following must remain indexed by
the active dimension `n` while field operations and the available basis remain
indexed by ambient `m`:

- the parent vanishing polynomial `s_n` and derivative `c_n`;
- the highest LCH basis index `N-1`;
- the block count and block shifts `omega_(i*P)` or `omega_(i*T)`;
- the high-rate normalization `p_(N-T)`;
- the erasure universe and locator evaluations, which are restricted to `V_n`;
- the sets used by the low-rate derivative/interpolation combination.

The basis polynomials `s_j` for `j<n` are the ambient field's existing prefix
basis polynomials; they are not regenerated in a new field.  Thus an active
parent is an evaluation restriction, not a new GF(2^n) representation.

## Wire profiles

### `LEO2_PROFILE_LEGACY_HIGH_V1`

Define

    T = ceil_pow2(R)
    N = ceil_pow2(K + T)
    D = N - T
    S = D - K

The parent is `[N,D]`, with coordinates:

    0 .. R-1       transmitted parity, in public parity order
    R .. T-1       punctured parent parity
    T .. T+K-1     transmitted systematic message
    T+K .. N-1     shortened known-zero message

The high encoder interpolates the degree-below-`D` polynomial from the systematic
region.  The existing implementation uses R10 Lemma 1 property 2: sum shifted
`T`-point IFFTs of message blocks into the last coefficient block, followed by a
`T`-point FFT at shift zero.  Only parity coordinates `0..R-1` are emitted.  The
profile claims compatibility only for cells that match the legacy golden vectors
and the large differential gate byte for byte.

### `LEO2_PROFILE_LOW_V1`

Define

    P = ceil_pow2(K)
    N = ceil_pow2(P + R)

The parent is `[N,P]`, with coordinates:

    0 .. K-1       transmitted systematic message
    K .. P-1       shortened known-zero message
    P .. P+R-1     transmitted parity, in public parity order
    P+R .. N-1     punctured parent parity

The encoder interpolates the padded message with one `P`-point IFFT at shift zero
and evaluates the coefficients in shifted `P`-point blocks beginning at `P`.
Completely punctured blocks are skipped and the last block is output-truncated.
This profile is new and makes no legacy parity claim.

For both profiles, the transmitted code is `[K+R,K]`.  Shortening is represented
as known zero input to parent interpolation; it is not included in the erasure
locator.  Puncturing is a permanent erasure.  MDS behavior is checked against
direct interpolation for every tested subset of `K` transmitted coordinates.

`AUTO` is deterministic: it selects legacy high when `R <= K` and low when
`R > K`.  Runtime ISA selection may change kernels only.  Field auto-selection
uses GF8 when the selected `N <= 256`, otherwise GF16 when `N <= 65536`.
Explicit field selection is part of the code identity and impossible combinations
are rejected.

## Locator and generic erasure identity

For erasure set `E` in active coordinates, define

    Lambda(x)  = product over e in E of (x - e)
    Lambda'(e) = product over u in E, u != e of (e - u)

At a surviving point `a`, setup stores `log Lambda(a)`.  At an erased output
`e`, it stores `log Lambda'(e)`; reveal applies the additive inverse of that
logarithm, which is multiplication by `1/Lambda'(e)`.  These values depend on
the erasure pattern but not shard bytes and therefore belong in an immutable
decode plan.

For a code polynomial `f` of degree below `D`, if `|E| <= N-D`, then
`fhat=f*Lambda` has degree below `N`.  Its evaluations are zero at erasures and
`f(a)*Lambda(a)` at survivors.  Interpolate `fhat`, formally differentiate it,
and evaluate only requested erased coordinates.  Since `Lambda(e)=0`,

    fhat'(e) = f(e) Lambda'(e)
    f(e)     = fhat'(e) / Lambda'(e)

This direct identity is the production fallback.  It is independent of the
specialized low/high complexity reductions and works with fewer than the maximum
number of erasures.

### Direct repair rows from an existing locator

The following is a Leopard2 derivation used only when a reusable plan already
owns both the transform locator and a direct-repair execution view.  Let `V` be
the active parent coordinate set, `E` the plan's exact erasure set, and
`A = V minus E`.  Thus `A` includes known-zero shortened systematic coordinates
and exactly the deterministic `K` received public coordinates selected by the
plan.  Write `Z_V = Lambda Z_A`.  The derivative `c_N = Z_V'` is constant and
nonzero on the additive parent subspace, hence for a survivor `s` and erased
original `x`:

    c_N = Lambda(s) Z_A'(s)
    c_N = Lambda'(x) Z_A(x)

The Lagrange coefficient carrying received value `f(s)` into `f(x)` is

    Z_A(x) / ((x + s) Z_A'(s))
      = Lambda(s) / ((x + s) Lambda'(x)).

Characteristic two turns subtraction into addition, and Cantor coordinate
addition is index xor.  With the locator stored in logarithmic form, the fixed
multiplier is therefore

    term_log = locator[s] - locator[x] - log(s xor x)  (mod 2^m - 1).

Known-zero shortened members of `A` affect the canceled parent derivative but
need no byte-execution term.  Native-high and xor-`P` translated-low views use
their own execution coordinates; translation preserves every pairwise xor.
Tests compare every constructed row against the independent systematic
generator-matrix solver before the locator shortcut is allowed in production.

R20 and the generic decoder in R16 provide related locator constructions.  Legacy
Leopard evaluates locator logarithms with a full-field FWHT.  Leopard2 plan setup
first preserves that calculation as a compatibility oracle, then compares an
active-parent direct/product construction.  No full-field setup optimization is
accepted unless its `Lambda` and derivative values match direct products.

For a power-of-two active parent `N=2^n` embedded in `GF(2^m)`, dense locator
setup can remain entirely inside the additive group `V_n`.  Let `q=2^m-1`, let
`a[e]` be the erasure indicator, and define `ell[i]=log(omega_i)` with
`ell[0]=0`.  Cantor-coordinate addition gives `omega_x+omega_e=omega_(x xor e)`,
so the logarithmic locator evaluation, including the derivative convention at
a root, is the XOR convolution

    lambda[x] = sum over e in V_n of a[e] ell[x xor e]  (mod q).

This is a new active-subspace derivation of the legacy convolution, not a
truncation of its transformed full-field table.  If `H_N` is the unnormalized
Walsh transform, then `H_N H_N = N I (mod q)`.  The full-field implementation
needs no visible scale because `2^m = 1 (mod q)`.  A proper active parent does:

    lambda = H_N( (N^-1 H_N(ell)) pointwise-multiplied-by H_N(a) ),
    N^-1 = 2^(m-n) = 2^m / N  (mod q).

Leopard2 precomputes the scaled fixed kernel `N^-1 H_N(ell)` once for every
supported proper parent and transforms the caller-owned `N`-entry output buffer
in place.  Dense setup therefore uses `O(N log N)` work and `O(N)` output storage
without a full-field temporary.  Omitting the `N^-1` factor is correct only at
full field size and is an active-parent normalization bug.  The derivation is
checked exhaustively for every erasure subset of every parent in test-only
`GF(2^4)`, exhaustively through `N=16` in production GF8 coordinates, and
differentially against the retained full-field oracle for all GF8/GF16 parent
sizes.  Dense large-GF16 tests additionally sample coordinates with the direct
`O(E)` sum, so correctness does not rely only on two paths sharing Walsh and
modular-add primitives.  The samples explicitly exercise sums that wrap modulo
65,535 and canonicalize the legacy zero-log sentinel.

The permanent-locator API accepts either the union bitmap used by current plans
or a dynamic-only bitmap.  Its sparse path adds dynamic products to the cached
permanent base; its dense path forms the explicit union before the active Walsh
transform.  Thus the setup result does not depend on which calling convention
crosses the measured direct/Walsh threshold.

Linearity gives `lambda_(P union D) = lambda_P + lambda_D (mod q)` for
disjoint permanent and dynamic sets.  This does not make a dense cache useful:
both a cached coordinate-domain output and a cached post-kernel transform-domain
vector retain the two dynamic Walsh transforms and require one additional
`N`-entry modular-add pass.  The explicit union path needs neither that pass nor
an `N`-symbol persistent cache.  Current deterministic survivor selection keeps
exactly `K` of the `K+R` transmitted coordinates, so the real-plus-virtual set
outside the permanent mask has exactly `R` members.  The codec's sparse-cache
decision from `recovery_count` therefore exactly matches the later dynamic
direct/Walsh decision.  The derivation, candidate operation counts, exhaustive
decomposition tests, and negative dense promotion result are recorded in
`docs/leopard2_permanent_locator_cache.md`.

When a specialized decoder requires exactly `N-D` erasures but more than `D`
coordinates are available, the plan marks deterministic surplus received
coordinates as virtual erasures.  The baseline rule prefers surviving systematic
coordinates and then the lowest coordinate index.  A received-subset experiment
may reduce work, but may not change decoded bytes.

## Active-parent low-rate decoder

For `K_parent=P=2^k`, split the active `N=2^n` evaluations into `P`-point cosets.
After weighting received values by `Lambda`, interpolate each block to
`fhat^(i)(x)`.  Corollary 1 of R10 applies for any `k <= n <= m`; the active form
uses `V_n`, not the full ambient field:

    g(x) = derivative(fhat^(0)(x))
         + (sum over a in V_n \\ V_k of 1/a) fhat^(0)(x)
         + sum from i=1 to N/P-1
             [ c_k / s_k(omega_(i*P)) ] fhat^(i)(x)
             modulo s_k(x)

Evaluate `g` only on the message subspace `V_k`.  For each missing original
coordinate `e`, recover

    f(e) = g(e) / Lambda'(e)

This is R10 Algorithm 4 with its outer set changed from `V_m` to the active
`V_n`.  No `c_n` appears in this formula because `fhat=f*Lambda` has degree below
the active parent and the recovery uses its ordinary derivative directly.  Tests
must nevertheless regenerate every `c_k/s_k(shift)` coefficient independently.

Parity recovery is deliberately absent from this hot path; explicit parity rebuild
re-encodes after original recovery.

## Active-parent high-rate decoder

Let `T=2^t`, `D=N-T`, and select exactly `T` parent erasures (real, permanent
punctures, plus deterministic virtual erasures).  Let `ftilde` equal a received
value at selected survivors and zero at selected erasures.  Equality of evaluations
on `V_n` gives

    f(x)Lambda(x) = ftilde(x)Lambda(x) + q(x)s_n(x)

for `deg(q)<T`.  Divide `ftilde` at LCH coefficient boundary `D` to obtain the
`T`-coefficient polynomial `h`.  R10 equation (10), valid for `t<=n<=m`, yields

    z(x) = h(x)Lambda(x) + p_(N-T) q(x)s_t(x)

with `deg(z)<T`.  Algorithm 5 obtains `h` by summing shifted `T`-point IFFTs,
evaluates it on `V_t`, multiplies by locator values, and inverse-transforms those
values to `z`.  Only message blocks containing requested missing originals need
the final shifted evaluation.

The shift-zero block admits one exact execution cancellation.  Write its raw
evaluation vector as `F_0`, and let `H_later` be the coefficient-vector sum of
all nonzero-shift blocks.  Since the forward and inverse transforms on the same
active coset are mutual inverses and the transform is linear,

    FFT_0(IFFT_0(F_0) + H_later)
        = F_0 + FFT_0(H_later).

Production therefore does not inverse-transform block zero.  If
`H_later != 0`, the first later inverse transform seeds the coefficient
accumulator; after the shift-zero forward transform, the live rows of `F_0`
are XORed into the result.  If `H_later = 0`, `F_0` is already the required
evaluation vector and both transforms are omitted.  The equality uses the
same active length, normalized LCH basis, and shift on the inverse/forward
pair, so it introduces no new subspace, normalization, locator, shortening,
puncturing, or coordinate-map factor.  It is a direct linear simplification of
R10 Algorithm 5 lines 3-4.

The promoted GF8 AVX2 execution composes, but does not reorder, the locator
diagonal with the first two inverse layers.  For four rows, let
`u_i = h_i Lambda_i` when row `i` is live and zero otherwise.  With the three
ordinary inverse-butterfly skew factors `a`, `b`, and `c`, the stored boundary
is exactly

    v1 = u1 + u0              v0 = u0 + a v1
    v3 = u3 + u2              v2 = u2 + b v3
    o2 = v2 + v0              o3 = v3 + v1
    o0 = v0 + c o2            o1 = v1 + c o3

where addition is XOR.  This is a labeled execution derivation obtained by
substitution into Leopard's existing two inverse radix-two layers, not a new
decoder identity.  Locator weight logs zero and the GF8 modulus log 255 both
mean multiplication by one.  In the skew positions 255 instead means the
zero-skew sentinel, so the corresponding product term is absent while the XOR
edge remains.  Scalar-reference known-answer tests enumerate all sixteen live
masks, identity and ordinary locator weights, each skew sentinel independently,
disjoint and exact in-place storage, vector widths, and arbitrary byte tails.
The full scale-then-IFFT path remains the oracle outside the measured dispatch
region documented in `docs/leopard2_high_decode_no_copy.md`.

For a missing original `e` outside `V_t`, differentiating the active congruence
introduces `c_n=s_n'`:

    f(e) Lambda'(e) = q(e) c_n
    q(e)             = z(e) / [p_(N-T) s_t(e)]
    f(e)             = c_n z(e)
                       / [p_(N-T) s_t(e) Lambda'(e)]

R10 equation (63) omits `c_n` only because its full-field `c_m=1`.  For a missing
parity `e` in `V_t`, the active form is

    f(e) = c_n [z'(e) - h(e)Lambda'(e)]
           / [p_(N-T) c_t Lambda'(e)]

The stable API defaults to original recovery, so this parity-only derivative work
is omitted unless parity rebuilding is explicitly requested (and re-encoding is
normally simpler).  Both active formulas are new derivations from R10 equations
(55)-(66) and are gated by exhaustive small-field and direct GF8/GF16 tests.

## C6 exact-prefix GF8 boundary code (experimental only)

C6 studies a separately defined exact code for the cells where a frozen dyadic
profile has a 512-coordinate GF16 parent even though `K+R <= 256`.  It uses the
legacy GF8 field representation but not a legacy parent coordinate map.
Integers `0..255` are symbol bit patterns in Leopard's Cantor-coordinate
representation, so `i xor s` denotes `omega_i + omega_s`; it is not arithmetic
on raw polynomial-basis integer encodings:

    source points = {0,...,K-1}
    parity points = {K,...,K+R-1}.

Let `S={0,...,K-1}` and `Z_S(x)=product_(s in S)(x+s)`.  Characteristic-two
addition is coordinate XOR.  Direct Lagrange evaluation gives

    w_i    = 1 / product_(s in S, s != i)(i+s)
    G(q,i) = Z_S(q) w_i / (q+i),       q in {K,...,K+R-1}.

Every denominator is nonzero because all public evaluation coordinates are
distinct.  The systematic generator is identity rows followed by `G`; any `K`
rows are an evaluation map at `K` distinct field elements and are nonsingular.
This is the generic exact-parameter RS construction described in the context of
R13, R14, R20, and R21, specialized by a new derivation to Leopard's public
Cantor-coordinate elements.  `experiments/leopard2/non_power_of_two/c6/algebra.py`
checks the derivation using independent carryless polynomial arithmetic and
exhaustive GF(2^4) rank tests.

For an ordered missing-original set `M={m_0,...,m_(L-1)}`, the current C6
checkpoint requires the parity prefix to be received and deterministically
selects rows `0,...,L-1`.  It does not yet plan around parity erasures.  Let

    A[e,u] = G(e,m_u)
    b_e    = p_e + sum_(i not in M) G(e,i) d_i.

Then `A d_M = b`, so plan setup computes `A^-1` and folds both kinds of input
into one fixed coefficient list per missing output:

    d_(m_u) = sum_e A^-1[u,e] p_e
              + sum_(i not in M) (
                    sum_e A^-1[u,e] G(e,i)) d_i.

All signs are additions because the field has characteristic two.  Distinct
evaluation points make `A` nonsingular: otherwise a nonzero degree-`<K`
polynomial could vanish at the `K-L` surviving systematic points and the `L`
selected parity points.  Byte-heavy execution consequently needs only fixed
multiply-adds and no inverse, matrix operation, allocation, or scratch.  The
Python oracle independently constructs and inverts `A`, folds the same terms,
and evaluates them with carryless polynomial arithmetic.  This direct repair
derivation follows the generic interpolation context of R20/R21 but is new for
the exact-prefix coordinate convention above.

No active-parent normalization, shortening factor, or puncturing factor can
make this code wire-equivalent to the frozen GF16 parent.  That parent has 512
distinct coordinates, while GF8 has only 256 elements; a GF8 transform cannot
evaluate the same coordinate set.  C2-style truncation can preserve the parent
only while retaining GF16.  C6 therefore recorded explicit unequal generator
coefficients and treated the exact prefix map as a then-unversioned new-profile
candidate.  Its coordinate digests remain external research fingerprints, not
serialized code IDs.  C7 subsequently froze that same prefix mathematics as
experimental
family 3/version 1/map 1 after the independent comparison below.  See
`docs/leopard2_c6_gf256_rescue.md` for the C6 executable and promotion gates.

## C7 exact-low coordinate selection (experimental only)

C7 generalizes the C6 prefix construction to either declared legacy field when
`K+R` fits the field.  The parity-row derivation above is equivalently

    G[j,i] = Z(q_j) / ((q_j + x_i) Z'(x_i)),

with `x_i=omega_i`, `q_j=omega_(K+j)`, and
`Z(X)=product_s(X+x_s)`.  Distinct points make every factor nonzero, proving
that all `K*R` entries are nonzero.  A separately implemented monomial
Vandermonde inverse produces every row again; this is the independent oracle
rather than a restatement of the barycentric formula.  The declared GF16
`K=3,R=500` witness compares all 500 rows and 1,500 coefficients in both the
Python algebra and standalone C++ gates.

For one global affine map `phi(X)=aX+b`, `a != 0`, the numerator gains `a^K`.
The explicit difference and the `K-1` derivative factors together gain the same
`a^K`, so the systematic generator is invariant.  The exhaustive GF(2^4) test
checks all 240 affine maps for all 120 valid `(K,R)` geometries.  This identity
does not apply to the separately studied non-affine aligned union

    systematic = omega_0 .. omega_(K-1)
    parity     = omega_P .. omega_(P+R-1), P=ceil_pow2(K).

C7 constructs that union explicitly, rejects it when `P+R` exceeds the field,
classifies whether the whole ordered set is one affine image, and otherwise
compares its changed coefficients, coefficient-one specializations, and exact
aligned-dyadic fragment count.  It also exhausts every GF(16) parity subset with
the systematic prefix fixed.  Those searches remain dense and can define
different parity bytes.  They do not inherit prefix identity or C6 timings.

Family 3/version 1/map 1 consequently fixes the ordered prefix coordinates,
has `parent_count=K+R`, `padded_side=K`, and has no shortening or puncturing.
The Experiment-W serializer rejects coordinate-, shortening-, and
puncturing-set digest TLVs for this family as redundant noncanonical spellings;
changing any of those sets requires a new version instead of a digest override.
Family 4 is reserved and rejected for future exact-high work.  The canonical
coordinate-byte encoding, buffer layout, executable evidence, and remaining
promotion gates are defined self-contained in
`docs/leopard2_c7_exact_low.md`.  The production constructor and AUTO dispatcher
continue to reject this experimental profile.

## C8 exact high-rate parity solve (experimental only)

C8 freezes a separate parity-first candidate solely for research:

    parity points = {0,...,R-1}
    message points = {R,...,R+K-1}
    polynomial degree < K.

Let the complete transmitted point set be `X={0,...,K+R-1}` in Leopard's
Cantor-coordinate representation and define the dual weight

    u_i = 1 / product_(ell != i) (omega_i + omega_ell).

For each `0 <= a < R`, the Lagrange leading-coefficient identity gives

    sum_i u_i omega_i^a f(omega_i) = 0

for every polynomial `f` of degree less than `K`.  Partitioning the first `R`
columns from the final `K` columns yields the exact parity constraints

    A p + B m = 0,
    p = inverse(A) B m.

The sign disappears in characteristic two.  This is a new specialization of
the generic generalized-RS dual and interpolation identities in R20/R21 to the
C8 point order; no transform implementation is copied.  The independent C8
oracle constructs `A` and `B` from the full-set vanishing derivative, then
requires `inverse(A)B` to equal cardinal Lagrange rows built only from the
message set.  Newton divided differences and a dyadic Schur-complement solve
provide two additional algebraically independent execution forms.

For non-power `R`, with `b` the largest lower power of two, partition

    A = [ A11 A12 ; A21 A22 ]

and use

    S = A22 + A21 inverse(A11) A12.

The scalar Schur executor must reproduce the full solve exactly.  Its dense
cross-block work is measured rather than assumed sparse.

The C8 and legacy-high generators coincide only if legacy padding disappears:

    R = T = ceil_pow2(R),
    K + R = N = ceil_pow2(K + T).

Thus both `R` and `K+R` must be powers of two.  In that case both definitions
interpolate the same `K` values at `R..R+K-1` and evaluate `0..R-1`.  Outside
that gate legacy interpolation additionally fixes shortened parent values, or
uses different message coordinates, so C8 records explicit coefficient and
byte witnesses and claims a new profile.  Exhaustive GF(2^4) MDS/rank checks,
GF8/GF16 direct algebra, all available SIMD backends, and ASan/UBSan gate the
derivation.  The candidate remains unserialized and default-off; details and
the no-promotion result are in `docs/leopard2_c8_exact_high.md`.

R13 Appendix-A Algorithm 8 remains a prefix/coset epsilon transform.  C8's
known message points are the suffix `R..R+K-1`, which is not generally an
aligned `shift XOR (0..K-1)` set.  Applying Algorithm 8 by renumbering would
change this coordinate profile.  C8 therefore labels its Tang--Han row a
control/lower model until a suffix-set derivation is proved.

## Independent oracles and promotion gates

The optimized LCH code is never its own only oracle.

1. Scalar GF arithmetic is generated from the declared irreducible polynomial and
   Cantor basis.  All nonzero products and inverses are cross-checked against
   direct polynomial multiplication/reduction.
2. A direct Lagrange encoder interpolates the parent systematic coordinates and
   evaluates requested parity coordinates without calling an LCH transform.
3. A direct Gaussian/Lagrange decoder reconstructs from any selected `K` public
   coordinates.
4. A test-only GF(2^4) enumerates parent/profile choices and erasure subsets.  It
   verifies the active `c_n` factors, locator derivatives, all systematic basis
   messages, and MDS behavior.
5. GF8/GF16 compare direct, generic, specialized, and pruned paths.  Legacy-high
   parity is additionally compared with old Leopard.
6. Full GF8 `N=256` cases are compared with XDRS only after its compile-time
   selection is changed from default Algorithm 1 to low Algorithm 2 or high
   Algorithm 3 and its `CHECK` path is enabled.

The active-parent gates now cover production constants and bare transforms in
addition to end-to-end recovery.  Production may therefore select the
specialized low/high decoder in its measured region; the generic decoder remains
a forced oracle and fallback.  Unknown error correction from R11-R13 is outside
the erasure hot path.

The committed direct-oracle test currently reports:

    field_checks=727117
    subspace_checks=1375
    normalization_checks=10
    high_basis_cases=127250
    high_recovered_symbols=916190
    missing_cn_detected=31
    profiles=170
    mds_subsets=110868

The `missing_cn_detected` counter consists of active-parent high-rate cases where
deliberately omitting `c_n` produces the wrong result.  This makes the active
normalization test sensitive to the exact full-field-to-parent error described
above, rather than merely exercising a formula that would also pass without the
new factor.

The production active-LCH differential additionally reports:

    skew_factors=65790
    normalizers=65792
    subspace_factors=131608
    prepared_factors=1837940
    forward_symbols=269115
    inverse_symbols=666336
    independent_inverse_symbols=196544
    derivative_symbols=198652
    lane_symbols=532480
    transforms=2280

It covers every aligned coset of every GF8 dyadic transform length and the first
and last aligned cosets of every GF16 power-of-two length through 65,536.  Small
inverse transforms are compared with independent Lagrange interpolation; larger
GF16 inverse transforms recover coefficient vectors from independently generated
full-coset values.  Those values use the subspace identity
`s_(j+1)(x)=s_j(x)(s_j(x)+s_j(v_j))`, with every `s_j(v_j)` first generated by
direct products; strategic outputs are also checked with the direct-product
evaluator.  Lane-varying GF8 and GF16 cases verify every byte lane and every
low/high ALTMAP pair rather than repeating one symbol across a SIMD tile.  Forced
scalar and available SIMD builds execute the same test-only hooks into the
production byte kernels.  Production branch-free pruning is deliberately not
claimed here; it remains a separate task.

## GF16 byte-tail and physical-storage profiles

Legacy GF16 uses an ALTMAP tile containing 32 field symbols as 32 low bytes
followed by 32 high bytes.  A partial 64-byte tile cannot simply be zero-padded
and its parity truncated: multiplying a present low byte by a GF16 coefficient
can produce a nonzero high byte, and discarding that byte makes later erasure
recovery impossible.

Native Leopard2 GF16 safely supports every positive even byte length.  A final
`2q` application bytes are compactly scattered as `q` low and `q` high bytes of
complete GF16 symbols; the inverse gather restores the application buffer.
Native odd lengths remain unsupported.

The separately identified `LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1` embeds an odd
`B`-byte application payload into a physical `W=B+1` systematic shard as
`payload || 0`.  Every systematic and parity coordinate transmits all `W`
bytes, and the ordinary compact-even GF16 code is applied at physical length
`W`.  If `G_S` is the generator submatrix for any `K` public coordinates, its
GF16 inverse recovers every padded symbol lane and therefore the injected
payload.  The final parity byte is generally nonzero and may not be punctured.
This changes application framing and persistent shard-layout identity, but not
the parent-code coordinates, GF16 representation, transform, locator, or
decoder normalization mathematics.  See `docs/leopard2_gf16_tails.md` for the
alphabet impossibility proof, alternative-construction disposition, and full
wire semantics.

## Formula-to-source map

| Item | Source or derivation |
| --- | --- |
| Cantor coordinates, `s_j`, normalized `Xbar_i`, additive FFT | R16, R17; R10 eqs. (1)-(12) |
| `FFTSkew` quotient and `p_i` constant generation | R16, R17; R10 eqs. (1)-(12), direct specialization above |
| normalized-basis derivative accumulation | R16 decoder lineage; product-rule derivation above |
| shifted-block encode and IFFT summation | R10 Lemma 1, eqs. (15)-(21) |
| full-parent direct systematic generator row `Z_S(q) / ((q+x_i) Z'_S(x_i))` | shared new Lagrange derivation in `docs/leopard2_direct_encode.md`; independently checked by `tests/leopard2/direct_oracle.cpp` and `test_direct_encode.cpp` for every bounded high/low GF8/GF16 profile |
| generic locator/derivative erasure identity and active-parent Walsh normalization | R10 eqs. (22)-(29); R16 decoder; new XOR-convolution derivation above |
| low weighted block derivative combination | R10 Theorem 1, Corollary 1, Algorithm 4 |
| high `h`, `z`, and message-only evaluation | R10 Corollary 2, eqs. (52)-(72), Algorithm 5 |
| active-parent factor `c_n` in high recovery | new derivation above from congruence modulo `s_n`; direct tests required |
| exact/truncated parent-preserving candidates | R15; R22 for general TFT design principles |
| arbitrary-parameter new-profile candidates | R13, R14, R18-R21; never assumed wire-equivalent |
| C6 exact-prefix GF8 generator `Z_S(q) w_i / (q+i)` | new specialization of direct Lagrange evaluation, with arbitrary-parameter context from R13, R14, R20, R21; exhaustive/direct evidence in the C6 algebra experiment |
| C8 exact-high dual solve `p=inverse(A)B m` and Schur split | new specialization of generalized-RS dual/Lagrange identities from R20/R21; exhaustive/direct evidence in the C8 algebra and executable experiments |
| exact-prefix Lagrange/Newton-to-LCH conversion | R15 Algorithms 1 and 3; executable derivation and direct-algebra checks in `experiments/leopard2/non_power_of_two/c3b/fast_inverse.py` |
| arbitrary-epsilon inverse and completed evaluations | R13 Appendix A, Lemmas 8-9 and Algorithm 8; executable derivation and direct-algebra checks in the C3b experiment |

## Literature refresh

The required arXiv, IEEE Xplore, and DBLP searches were run on 2026-07-16 with
queries `Reed-Solomon LCH-FFT` and `Reed-Solomon LCH-FFT erasure decoding`, and
with a title search for R10.  R10 is now indexed as IEEE Transactions on
Information Theory 72(6), 3784-3798 (2026), DOI
`10.1109/TIT.2026.3685291`.  Candidates reviewed included the earlier ISIT
low-rate precursor, 2024 work on Vandermonde decomposition, the 2025 generalized
RS/alternant decoder, and error/error-erasure or different-code-family results.
No later pure-erasure LCH result was found that supersedes R10 or changes the
production plan.  This is a scoped search result, not a universal state-of-the-art
claim.  The same endpoints must be checked again before release reporting.

## Retrieval record

Research files are ignored under `.research/leopard2`; PDFs and reference
repositories are not committed.  The following primary PDFs were retrieved on
2026-07-16:

| Ref | SHA-256 |
| --- | --- |
| R10 | `8b23e555dc985a5cf51336e6f90a0a5f18e8cc20c5c74b1961e2aa0970abab50` |
| R11 | `51e0381a59ad6f896d17d879b51bfc22fcd6d1e578691f459ba1f390fca72cfe` |
| R12 | `6e8ff625c93bb0935c69510a85c745727ef031e9e2a38d60f4d34bc9d89102f0` |
| R13 | `7c25b702f693a4e8a4cd27dfa9d8ec53a8501b0f61cb1ece952d1e65dfd90ba0` |
| R15 | `971e4ec04c5c066df28105a17952100486cdbb206b643ddd84bcea2ce08d12fd` |
| R16 | `8aade67a885b8fb83eb7cc93e68650d5f6771a8d161042fa201ea5a6a8652617` |
| R17 | `052c1e7667a96f2865927ce8f22630a3afabd888313b786f1be46b88cfc830f2` |
| R19 | `e5dafd9ee609a9ebbf00202febb83f42d71389146f3d6c0a5da836944ee0e523` |
| R20 | `264e2ba4ce166278048c482e496b95a0861654d1e81a5b79ace3e6421cfecce5` |
| R21 | `e04a017716ba021629ebe84a10c3a4b591d4c82b4f2f76121e5da447998d6802` |
| R22 | `8b1748fbcd326771daeb60f7ffc6ef6267696835bfba75ab2b45c8ab663dd8f4` |

XDRS was inspected at commit
`ae05a779e7f44be13c3d34e14d15b08b4bc02404` (Apache-2.0).  Its GF256
polynomial `0x11d`, Cantor basis, and high/low coordinate conventions match the
corresponding mathematical layouts, but its default benchmark selects Algorithm 1
and disables checking.  No XDRS source is copied into Leopard2.

Comparison repositories were refreshed into the ignored research cache on
2026-07-16. No third-party source is committed:

| Ref | Commit | License observed in checkout |
| --- | --- | --- |
| R04 ISA-L | `e8cc5e87fc64b4da434f32bc1fa18184622a4998` | BSD-3-Clause |
| R05 Jerasure | `de1739cc8483696506829b52e7fda4f6bb195e6a` | BSD-3-Clause |
| R06 FastECC | `b8ca7db6bf5556185c96009b161e8aec82af734e` | Apache-2.0 |
| R07 ECC-Benchmark | `c43d4290f8525351821f7b04791cee3bdfbaccdd` | MIT |

`tools/leopard2_external_comparison.py cache` verifies these identities when
the optional cache is present. The current x86-64 host has no system NASM, so
the default-off ISA-L checkpoint builds official NASM 2.16.03 from
https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz after
verifying archive SHA-256
`1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148`.
It remains inside the ignored research cache. The host has no installed
GF-Complete package, which Jerasure 2.0 requires. Those host facts are not
mathematical exclusions and do not turn an unmeasured cell into a benchmark
result. The adapter methodology and ISA-L BSD notice are recorded in
`docs/leopard2_isal_comparison.md` and
`experiments/leopard2/isal_compare/NOTICE` respectively.

R14 and R18 were available only through their publisher landing pages in this
pass.  They are inputs to later exact-size experiments, not to the initial
production implementation.

## Complete URL bibliography

### Projects and comparison implementations

- R01 Leopard-RS: https://github.com/catid/leopard
- R02 XDRS: https://github.com/fastecc/xdrs
- R03 Beads: https://github.com/gastownhall/beads and https://beads.gascity.com/
- R04 Intel ISA-L: https://github.com/intel/isa-l
- R05 Jerasure: https://github.com/ceph/jerasure
- R06 FastECC: https://github.com/Bulat-Ziganshin/FastECC
- R07 ECC-Benchmark: https://github.com/Bulat-Ziganshin/ECC-Benchmark
- R08 SSE2NEON: https://github.com/DLTcollab/sse2neon

### Papers and literature endpoints

- R10 Chen et al., *Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes Based on LCH-FFT*: https://i4ai.org/hanyunghsiang/IT2026.pdf
- R11 Chen et al., *Parallel Welch-Berlekamp Algorithm*: https://i4ai.org/hanyunghsiang/IT2025.pdf
- R12 Tang et al., *Fast Error and Erasure Decoding Algorithm for Reed-Solomon Codes*: https://i4ai.org/hanyunghsiang/CL2024.pdf
- R13 Tang and Han, *New Decoding of Reed-Solomon Codes Based on FFT and Modular Approach*: https://arxiv.org/abs/2207.11079 and https://arxiv.org/pdf/2207.11079
- R14 Tang and Lin, *Fast Encoding and Decoding Algorithms for Arbitrary (n,k) Reed-Solomon Codes Over GF(2^m)*: https://ieeexplore.ieee.org/document/8955804/
- R15 Coxon, *Fast Transforms over Finite Fields of Characteristic Two*: https://arxiv.org/abs/1807.07785 and https://arxiv.org/pdf/1807.07785
- R16 Lin, Chung, and Han, *Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes*: https://arxiv.org/abs/1404.3458 and https://arxiv.org/pdf/1404.3458
- R17 Lin, Al-Naffouri, and Han, *FFT Algorithm for Binary Extension Finite Fields and its Application to Reed-Solomon Codes*: https://arxiv.org/abs/1503.05761 and https://arxiv.org/pdf/1503.05761
- R18 Yu et al., *Reed-Solomon Coding Algorithms Based on Reed-Muller Transform for Any Number of Parities*: https://ieeexplore.ieee.org/document/10086680/
- R19 Tang et al., *A Fast Decoding Algorithm for Generalized Reed-Solomon Codes and Alternant Codes*: https://arxiv.org/abs/2502.02356 and https://arxiv.org/pdf/2502.02356
- R20 Didier, *Efficient Erasure Decoding of Reed-Solomon Codes*: https://arxiv.org/abs/0901.1886 and https://arxiv.org/pdf/0901.1886
- R21 Bostan and Schost, *Polynomial Evaluation and Interpolation on Special Sets of Points*: https://mathexp.eu/bostan/publications/BoSc05.pdf
- R22 van der Hoeven, *Notes on the Truncated Fourier Transform*: https://www.texmacs.org/joris/tft/tft.pdf
- R23 refresh endpoints: https://arxiv.org/search/?query=Reed-Solomon+LCH-FFT&searchtype=all , https://ieeexplore.ieee.org/search/searchresult.jsp?newsearch=true&queryText=Reed-Solomon%20LCH-FFT , and https://dblp.org/search?q=Reed-Solomon%20LCH-FFT

### Processor, profiling, and verification references

- R30 Intel Intrinsics Guide: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
- R31 Intel 64 and IA-32 manuals: https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html
- R32 Arm Neon intrinsics: https://arm-software.github.io/acle/neon_intrinsics/advsimd.html
- R33 Arm C Language Extensions: https://arm-software.github.io/acle/main/acle.html
- R34 Linux perf documentation: https://perfwiki.github.io/
- R35 Agner Fog optimization manuals: https://www.agner.org/optimize/
- R36 RISC-V Vector specification: https://github.com/riscv/riscv-v-spec
- R37 RISC-V vector intrinsics: https://github.com/riscv-non-isa/rvv-intrinsic-doc
- R38 NVIDIA CUDA documentation: https://docs.nvidia.com/cuda/
- R39 AMD HIP documentation: https://rocm.docs.amd.com/projects/HIP/
- R40 LLVM libFuzzer: https://llvm.org/docs/LibFuzzer.html
- R41 Clang AddressSanitizer: https://clang.llvm.org/docs/AddressSanitizer.html
- R42 Clang UndefinedBehaviorSanitizer: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
- R43 Clang ThreadSanitizer: https://clang.llvm.org/docs/ThreadSanitizer.html
