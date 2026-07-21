# IT2026 Algorithm 4/5 implementation audit

This audit compares Leopard2 at commit `10c8a3992f7983e5530d149274364bb11da63bbc`
with Chen et al., *Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes
Based on LCH-FFT* (R10):

<https://i4ai.org/hanyunghsiang/IT2026.pdf>

The conclusion is narrow: no correctness-level Algorithm 4 or message-only
Algorithm 5 formula gap was found in GF8 or GF16.  Several execution and
dispatch opportunities remain, including a newly identified translation that
allows Algorithm 4 to decode some legacy-high codewords without changing their
wire bytes.

## Source integrity and extraction limits

The locally cached PDF has 16 pages, is 349,188 bytes, and has SHA-256
`8b23e555dc985a5cf51336e6f90a0a5f18e8cc20c5c74b1961e2aa0970abab50`.
The `pdftotext -layout` extraction has 1,421 lines, 14,103 words, 16 form-feed
page boundaries, and SHA-256
`ed4819bbc12633266bb06ab32f472b16d6c7d4005ceba8b9c269aa1140243b8d`.
A fresh layout extraction was byte-identical to the cache.

This proves complete reproducible layout text, not semantic Markdown.  The
two-column ordering and some stacked fractions, products, superscripts, and
overbarred basis symbols remain ambiguous in plain text.  Algorithms 4 and 5
and their surrounding equations were checked visually against PDF pages 7 and
9 rather than inferred from a research summary.

The main paper anchors used here are:

- equations (1)-(12), pages 2-3: field points, additive subspaces, subspace
  polynomials, and the normalized LCH basis;
- equations (18)-(21), page 4: low/high systematic encoding;
- equations (22)-(29), pages 4-5: generic locator/derivative decoding;
- Theorem 1, Corollary 1, and equations (31)-(36), pages 5-6;
- Corollary 2 and equations (37)-(40), page 6;
- equations (41)-(51) and Algorithm 4, pages 6-7;
- equations (52)-(72) and Algorithm 5, pages 7-9; and
- the complexity and simulation protocol in Sections V-VI, pages 9-10.

## Algorithm 4 mapping

The paper's `N=2^m` and `K=2^k` are parent parameters.  Leopard2 uses an
active parent `N=2^n` embedded in the ambient field and padded low side
`P=2^k`; public K may be smaller than P through shortening.

| Paper step | Leopard2 source | Assessment |
| --- | --- | --- |
| 1. Locator values and inverse derivatives | plan setup `leopard2.cpp:5385-5410`; active Walsh/direct setup `LeopardFF8.cpp:2702-2847`, `LeopardFF16.cpp:2728-2875` | Faithful; dense setup is active-N, not full-field. |
| 2. Zero erasures, weight survivors by `Lambda` | GF8 `LeopardFF8.cpp:3487-3494,3781-3791`; GF16 `LeopardFF16.cpp:3501-3509,3802-3813` | Faithful. |
| 3. Shifted P-point inverse transforms | GF8 `LeopardFF8.cpp:3496-3506,3584-3616,3820-3858`; GF16 `LeopardFF16.cpp:3511-3522,3603-3634,3814-3885` | Faithful, with valid exact input pruning. |
| 4. Differentiate block zero | GF8 `LeopardFF8.cpp:3508-3511,3618,3818`; GF16 `LeopardFF16.cpp:3524-3527,3637,3840` | Correct at requested roots. |
| 5. Reduce later blocks with `c_P/s_P(shift)` | constants GF8 `LeopardFF8.cpp:3021-3068`, GF16 `LeopardFF16.cpp:3058-3105`; reductions GF8 `3513-3522,3620-3629,3859-3861`, GF16 `3529-3539,3639-3649,3886-3890` | Faithful. |
| 6. P-point forward transform | GF8 `LeopardFF8.cpp:3524-3529,3631-3649,3865-3883`; GF16 `3541-3546,3651-3669,3894-3913` | Faithful, with requested-output pruning. |
| 7. Reveal requested messages with inverse locator derivative | GF8 `LeopardFF8.cpp:3531-3534,3652-3659,3886-3894`; GF16 `3548-3552,3671-3680,3916-3924` | Faithful and message-only. |

`AddFormalDerivative` forms `f+f'` rather than materializing `f'` alone
(`LeopardFF8.cpp:2869-2883`, `LeopardFF16.cpp:2897-2916`).  Algorithm 4
evaluates the result only at requested locator roots, where the extra `f` term
is zero.  The same argument makes the active-parent scalar multiple of block
zero irrelevant at requested roots.  This is a root-only optimization, not a
general polynomial equality.

Production factors are compared with direct subspace-polynomial products in
`tests/leopard2/test_active_lch.cpp:403-481`.  End-to-end low recovery is
compared with direct interpolation in
`tests/leopard2/test_decode_low_acceptance.cpp:537-565`.

## Algorithm 5 mapping

The paper's `N-K=2^t` becomes the active parent redundancy `T=2^t`.  Public R
may be smaller than T through puncturing.

| Paper step | Leopard2 source | Assessment |
| --- | --- | --- |
| 1. Locator setup | same immutable active-N setup as Algorithm 4 | Faithful. |
| 2. Form `ftilde` from selected survivors | mapping/virtual erasures `leopard2.cpp:5279-5311`; GF8 staging/fusion `LeopardFF8.cpp:1673-1777,4114-4194,4352-4450`; GF16 copy-first staging `LeopardFF16.cpp:1540-1566,4138,4389-4425` | Faithful; selection supplies exactly T erasures. |
| 3. Sum shifted T-point inverse transforms into h | GF8 `LeopardFF8.cpp:4114-4194,4352-4450`; GF16 `LeopardFF16.cpp:4140-4221,4391-4486` | Faithful, with exact pruning/fused accumulation where qualified. |
| 4. Evaluate h on `V_t` | GF8 `LeopardFF8.cpp:4196,4452`; GF16 `LeopardFF16.cpp:4224,4489` | Faithful, but algebraically reducible below. |
| 5-6. Form `z=h*Lambda` and inverse-transform z | GF8 weighted boundary `LeopardFF8.cpp:1585-1670,4197-4199,4453-4455`; GF16 separate pass `LeopardFF16.cpp:4225-4233,4490-4500` | Faithful. |
| 7. Differentiate z | omitted | Faithful message-only specialization. |
| 8. Evaluate z on requested message blocks | GF8 `LeopardFF8.cpp:4201-4266,4457-4520`; GF16 `LeopardFF16.cpp:4235-4299,4502-4566` | Faithful and output-pruned. |
| 9. Evaluate z' on parity | omitted | Faithful message-only specialization. |
| 10. Reveal missing originals | GF8 `LeopardFF8.cpp:4267-4275,4521-4532`; GF16 `LeopardFF16.cpp:4300-4311,4567-4579` | Faithful with active-parent normalization. |
| 11. Recover parity erasures | omitted; explicit rebuild re-encodes | Faithful message-only specialization. |

The paper explicitly permits deleting lines 7, 9, and 11 when only messages
must be recovered.  The behavior is checked in
`tests/leopard2/test_decode_high_acceptance.cpp:565-600`.

For a proper active parent, the high reveal factor is

    f(e) = c_n z(e) / [p_(N-T) s_t(e) Lambda'(e)].

R10 equation (63) has no visible `c_n=s_n'` because the paper uses the full
field, where `c_m=1`.  Retaining `c_n` is a required generalization, not a
deviation.  `tests/leopard2/test_direct_oracle.cpp:301-422` exhausts the
small-field high identity and proves that omitting `c_n` fails on a proper
active parent.

## Active-parent result

The implementation does not blindly replace ambient m with active n.  It
retains the ambient GF8/GF16 Cantor representation and LCH basis prefix while:

- restricting transforms and locator evaluation to additive cosets in `V_n`;
- applying the active Walsh normalization `N^-1=2^(m-n)`;
- deriving low `c_P/s_P(shift)` factors on the active blocks;
- retaining high `c_n` and normalized `p_(N-T)`;
- treating shortened coordinates as known interpolation zeros;
- treating punctures as permanent erasures; and
- adding virtual erasures to select exactly the parent dimension of received
  symbols.

The derivation is in `docs/leopard2_math_and_sources.md:145-375`; coordinate
construction is in `leopard2.cpp:773-785,4527-4590`; immutable plan metadata is
built in `leopard2.cpp:836-1055,5207-5421`.

## Half-parent high/low translation identity

Let

    P = ceil_pow2(public K) = ceil_pow2(public R) = T.

Both profiles then have `N=2P` and parent dimension P.  Define
`tau(i)=i xor P`.  R10 equation (1) defines field points by binary basis
coordinates, hence

    omega_tau(i) = omega_i + omega_P.

Translation swaps the two P-point cosets.  If `f_H` is the high-profile
polynomial, `f_L(x)=f_H(x+omega_P)` has:

- the same K public systematic values;
- translated shortened zeros;
- low parity at coordinate `P+j` equal to high parity at coordinate j; and
- translated punctures.

The public parity ordinal and bytes are therefore identical, although the two
declared internal coordinate descriptions differ.  This follows from R10
equations (18)-(21) and the additive point definition.

The locator translates too:

    Lambda_L(x)  = Lambda_H(x + omega_P),
    Lambda_L'(x) = Lambda_H'(x + omega_P).

Consequently, a legacy-high codec can execute Algorithm 4 through translated
coordinate, locator, and output views without changing its encoder, profile,
field, or wire bytes.  This applies to every pair whose K and R round to the
same P, not only exact K=R; `(127,65)` is an example.  The candidate must remain
forced until exhaustive equivalence and pinned crossover evidence exist.

## Algorithm 5 zero-shift cancellation

Paper Algorithm 5 line 3 includes the shift-zero block and line 4 immediately
evaluates accumulated h on the same `V_t`.  Linearity gives

    FFT_0(IFFT_0(F_0) + H_other) = F_0 + FFT_0(H_other).

GF8 and GF16 now keep raw block zero out of the coefficient accumulator.  When
a later block contributes, its first inverse transform seeds the accumulator,
the remaining live blocks accumulate normally, and the transform at shift zero
is followed by XOR of the raw live block-zero rows.  When no later block is
live, block zero is staged directly and neither transform executes.  This
removes one T-point inverse transform for every Algorithm 5 execution and also
removes the T-point forward transform in the block-zero-only case.  Balanced
full-original-loss recovery is exactly such a case under the current
deterministic survivor selection.

`FinishHighDecodeSyndrome` implements the shared boundary in both fields for
prepared, materialized, tiled, and pruned execution.  Test-only counters prove
which transforms were elided, and synthetic hook tests cover an empty block
zero, a sparse block zero, a block-zero-only input, and a first later
exact-pruned contribution followed by accumulation.  Forced-specialized
GF8/GF16 end-to-end tests compare the result with the generic decoder.  This is
an execution simplification of paper lines 3-4, not a new decoder identity.

Other performance gaps are extra GF16 copy/diagonal passes, low tiled handling
of an empty block zero, a selector that is mostly workspace-driven outside two
narrow measured exceptions, and correctness-first rather than transform-first
surplus survivor selection.  None changes decoded bytes.

## Benchmark interpretation

The paper's throughput figures are important evidence, but they are not a
per-cell Leopard2-versus-Leopard1 guarantee:

- they use full-length GF8 N=256 parents, while arbitrary public Leopard2 cells
  can be shortened and punctured;
- the paper's K is parent dimension, not necessarily public K;
- the XDRS harness selects the low implementation at K<=128 and high above it,
  whereas Leopard2 AUTO currently selects legacy high whenever public R<=K;
- the reported groups contain 1,024 codewords with locator setup amortized once
  and random erasure patterns; and
- the executable comparison baseline is not necessarily the exact current
  Leopard1 compiler/backend/API path.

Plan-cold, one-stripe, maximum-loss, and adversarial-pattern measurements answer
different questions.  Production promotion therefore still requires
same-source byte equivalence, setup/execution separation, pinned timing, and
neighbor-cell regression gates.
