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

At a surviving point `a`, setup stores `Lambda(a)`.  At an erased output `e`, it
stores `1/Lambda'(e)`.  These values depend on the erasure pattern but not shard
bytes and therefore belong in an immutable decode plan.

For a code polynomial `f` of degree below `D`, if `|E| <= N-D`, then
`fhat=f*Lambda` has degree below `N`.  Its evaluations are zero at erasures and
`f(a)*Lambda(a)` at survivors.  Interpolate `fhat`, formally differentiate it,
and evaluate only requested erased coordinates.  Since `Lambda(e)=0`,

    fhat'(e) = f(e) Lambda'(e)
    f(e)     = fhat'(e) / Lambda'(e)

This direct identity is the production fallback.  It is independent of the
specialized low/high complexity reductions and works with fewer than the maximum
number of erasures.

R20 and the generic decoder in R16 provide related locator constructions.  Legacy
Leopard evaluates locator logarithms with a full-field FWHT.  Leopard2 plan setup
first preserves that calculation as a compatibility oracle, then compares an
active-parent direct/product construction.  No full-field setup optimization is
accepted unless its `Lambda` and derivative values match direct products.

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
| generic locator/derivative erasure identity | R10 eqs. (22)-(29); R16 decoder; direct derivation above |
| low weighted block derivative combination | R10 Theorem 1, Corollary 1, Algorithm 4 |
| high `h`, `z`, and message-only evaluation | R10 Corollary 2, eqs. (52)-(72), Algorithm 5 |
| active-parent factor `c_n` in high recovery | new derivation above from congruence modulo `s_n`; direct tests required |
| exact/truncated parent-preserving candidates | R15; R22 for general TFT design principles |
| arbitrary-parameter new-profile candidates | R13, R14, R18-R21; never assumed wire-equivalent |
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
the optional cache is present. The current x86-64 host has neither NASM (which
current ISA-L requires for its optimized build) nor an installed GF-Complete
package (which Jerasure 2.0 requires). Those host facts are not mathematical
exclusions and do not turn an unmeasured cell into a benchmark result.

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
