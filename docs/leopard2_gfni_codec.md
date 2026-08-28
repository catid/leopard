# GFNI affine multiplication inside the Leopard2 codec

## Status

**Promoted: production explicit backend member plus one exact AUTO encode
policy.**  Default builds compile
`Leopard2BackendGFNI.cpp` into the archive whenever the toolchain accepts
`-mgfni`; a GFNI-capable host selects it with `LEO2_BACKEND_GFNI` at context
creation.  On calibrated AMD family 1Ah/model 08h, an AUTO context may also use
it for the exact single-thread legacy-high GF16 `K=1000`, `R=200`, `T=256`,
64-KiB, full-output one-shot or ordinary one-item-batch encode.  AUTO still
reports AVX2, all other operations and cells retain their context table, and
parity is byte-identical to the AVX2 member.  Requirements 1-3 of the
production integration list at the end of this document are satisfied;
requirements 4-6 (isolated exact-main evidence, a second microarchitecture,
and the 512-bit follow-up) remain open.

The historical evaluation path — `LEO2_GFNI_VARIANT` via compiler flags,
running under the AVX2 identity — remains available for A/B experiments and
deliberately still fails the ISA-isolation audit.

`docs/leopard2_gfni_affine.md` closed with "a future production bead should
integrate the ZMM affine kernel behind dispatch, then repeat complete
encode/decode benchmarks".  This is that integration, at 256-bit width.

## What changed

`Leopard2BackendAVX2.cpp` gained a `LEO2_GFNI_VARIANT` compile-time variant
alongside the existing `LEO2_AVX512_VARIANT`.  A fixed multiplication by one
field element is a GF(2)-linear map on the symbol bits, so each byte-wide
component of the product is an 8-by-8 bit matrix and `VGF2P8AFFINEQB` evaluates
it in one instruction.

- **GF8.**  `product = vgf2p8affineqb(data, matrix)` replaces
  `shuffle(low_table, data & 15) ^ shuffle(high_table, (data >> 4) & 15)`:
  one instruction instead of six.
- **GF16.**  Leopard's GF16 block layout is 32 low bytes followed by 32 high
  bytes of the 16-bit symbols, so the 16-by-16 multiplication matrix decomposes
  into four 8-by-8 blocks and the product is
  `product_low  = A00 * low ^ A01 * high`,
  `product_high = A10 * low ^ A11 * high`.
  That is four affine operations and two XORs per 64-byte tile, replacing four
  ANDs, two shifts, eight shuffles and eight XORs.
- **Storage.**  GF8 reuses its 32-byte nibble-table storage shape: each 16-byte
  row holds one affine matrix duplicated, so the existing 128-bit broadcast
  fills all four 64-bit lanes.  GF16 stores only the four required 64-bit
  matrices in `FF16AffineTable`, reducing its table from 8 MiB to 2 MiB, and
  broadcasts each packed matrix at its use site.  The matrices are derived
  from the same `ff8_multiply_log` / `ff16_multiply_log` callbacks that build
  the nibble tables, so both backends come from one field definition.
- **Tails.**  Sub-vector scalar tails evaluate the stored matrix bitwise
  (`GFNIApplyMatrix`), so the variant needs no second table and no extra memory.
- **Unrolling.**  The GF8 multiply / multiply-add loop advances 64 bytes per
  iteration in this variant.  The affine form leaves one arithmetic operation
  per vector, so loop bookkeeping would otherwise be a large share of the issue
  slots; this also matches Leopard1's `mul_mem` / `muladd_mem` stride.
- **GF16 radix-four fusion at every size.**  `UseFusedButterfly4`
  (`LeopardFF16.cpp:114-130`) caps fusion at 128 bytes, so above that every GF16
  radix-four group runs as four separate radix-two sweeps over the same four
  shards.  The cap exists because a nibble multiplier needs eight broadcast
  table vectors, and three multipliers plus eight data vectors do not fit in
  sixteen `ymm`.  An affine multiplier needs four, and they are re-broadcast at
  each use rather than held, so the live set is eight data vectors plus four
  matrices plus two products.  In the GFNI variant `AVX2FF16Butterfly4` fuses at
  every size and `AVX2FF16Butterfly4Range` ignores the caller's `prefer_fused`
  hint, which exists only to bound nibble-table pressure.  The four shards are
  loaded and stored once per group instead of four times.

The vector data path stays 256 bits wide.  Backend identity, dispatch,
selectors, plan construction, scheduling policies, and the wire profile are
untouched, and none of the GF8/GF16 AVX2 nibble code paths changed.

The empirical operand order is the one the earlier experiment established:

    operand_bit = 8 * (7 - output_bit) + input_bit

`experiments/leopard2/gfni_codec/kernel_screen.cpp` re-derives the nibble tables
and the affine operand from one shared GF(2)-linear map and requires the two
kernel families to agree, so running it also re-pins that order on the host.

## Correctness evidence

- The full test suite passes on the GFNI build: 99 of 100 tests, including
  `leopard2_legacy_golden`, `leopard2_transform_differential`,
  `leopard2_boundaries`, `leopard2_gf16_tails`, `leopard2_gf16_padded_odd`,
  `leopard2_direct_encode`, `leopard2_direct_repair`, `leopard2_r1_xor`,
  `leopard2_backend_ops`, and both fuzz smoke targets.
- (Evaluation build only.)  The single failure was `leopard2_portable_isa`,
  which reports `vgf2p8affineqb` inside the object labelled as the AVX2 member.
  That audit was correct and drove the promotion: the production member owns a
  `gfni` checker class and the full suite is 103/103, with the audit passing on
  the GFNI-bearing archive.
- A randomized cross-build differential of 400 shapes
  (`experiments/leopard2/gfni_codec/cross_build_differential.py`, seed
  20260724) over random `K`, `R`, byte counts including ragged and sub-vector
  sizes, loss counts, forced fields, all five backend selections, one and two
  threads, and batch sizes one and two found **zero** digest differences and
  zero accept/reject differences against the stock build.  123 of the 400 shapes
  were rejected identically by both builds.
- Directed spot checks additionally matched on GF8 high, GF8 XOR `R=1`, GF16
  high, GF16 forced inflation, the `R > K` low profile, the `K=65 R=65` direct
  repair shape, `K=2 R=2` at 1 KiB, 64-byte shards, and a ragged 1000-byte GF16
  shard, for each of the `auto`, `avx2`, and `avx512` backend selections.

## Kernel screen

`kernel_screen.cpp`, one AMD Ryzen 9 9950X3D logical CPU, 32 shards, five
butterfly layers, best of 40 repetitions.  Times are microseconds for the whole
sweep; lower is better.

| GF8 butterfly, bytes/shard | avx2-256 | avx512-256 | avx512-512 | gfni-256 | gfni-512 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 256 | 0.61 | 0.62 | 0.22 | 0.30 | 0.18 |
| 1024 | 1.37 | 1.14 | 0.77 | 1.05 | 0.66 |
| 4096 | 7.25 | 6.19 | 4.61 | 5.41 | 4.11 |
| 16384 | 27.29 | 22.49 | 17.69 | 18.13 | 14.91 |
| 65536 | 103.06 | 77.80 | 71.44 | 74.62 | 71.38 |
| 262144 | 364.09 | 345.89 | 278.95 | 290.20 | 276.36 |

| GF16 multiply-add, bytes/shard | avx2-nibble | gfni-256 | gfni-512 |
| ---: | ---: | ---: | ---: |
| 256 | 0.39 | 0.29 | 0.22 |
| 1024 | 1.48 | 0.93 | 0.60 |
| 4096 | 5.83 | 3.62 | 3.53 |
| 16384 | 23.06 | 14.01 | 13.47 |
| 65536 | 94.58 | 72.57 | 70.99 |
| 262144 | 388.23 | 284.41 | 254.28 |

The 256-byte GF8 row is a sub-microsecond total and moved by nearly 2x between
two runs of this screen; treat it as directional only.  The rest reproduced
across runs.

Two separable effects are visible.  Widening the existing nibble kernel from
256 to 512 bits is worth about 1.09x to 1.57x on its own, which the shipped
AVX-512VL variant does not take because it keeps a 256-bit data path.  Replacing
the nibble tables with the affine instruction is worth a further 1.04x to 1.38x
at 256 bits.  Both together reach 1.32x to 2.08x for GF8 and 1.33x to 2.47x for
GF16 above the noisy smallest cell.  Only the 256-bit affine change is
implemented in this candidate.

## Whole-codec effect versus exact Leopard main

The table below is a clustered-ABBA screen against the exact Leopard main
adapter (baseline commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, built with
its canonical `-march=native` Release flags), on one pinned logical CPU of an
AMD Ryzen 9 9950X3D, three rounds per cell, five timing samples per child,
reuse 8, one thread.  Ratios are exact-main time divided by Leopard2 time, so
above one favours Leopard2.  All three workload digests matched exact main in
every cell for both builds.

This is a **screen, not promotion evidence**: it does not hold the CPU-pair
lease, does not verify build closure, and does not gate on sibling idleness.
The isolated `experiments/leopard2/main_compare/run_abba.py` campaign is still
required before any of this enters `docs/leopard2_vs_main_benchmark.md`.

Three configurations are compared: `stock` is `HEAD`, `GFNI` is the affine
variant with the GF16 fusion cap still in place, and `fused` adds the GF16
radix-four fusion.  Against stock Leopard2 itself, `fused` is a median 1.41x on
encode and 1.30x on decode, reaching 1.63x on the transform-bound cells.  The
kernel-level gains are larger — see the radix-four screen below — and the
difference is Amdahl: XOR-only zero-skew layers, staging, gather/scatter and
fixed per-call cost do not speed up.

| Cell | stock enc | GFNI enc | fused enc | stock dec | fused dec |
| --- | ---: | ---: | ---: | ---: | ---: |
| GF8 K=2 R=2, 1 KiB | 0.574 | 0.577 | **0.603** | 5.98 | **7.26** |
| GF8 K=4 R=4, 1 KiB | 0.682 | 0.853 | **0.859** | 5.27 | **6.51** |
| GF8 K=16 R=4, 1 KiB | 0.823 | 1.028 | **1.077** | 4.16 | **5.46** |
| GF8 K=8 R=8, 1 KiB | 0.885 | 1.146 | **1.121** | 4.42 | **5.72** |
| GF8 K=16 R=16, 4 KiB | 1.019 | 1.307 | **1.318** | 1.27 | **1.61** |
| GF8 K=129 R=1, 4 KiB | 1.008 | 1.045 | **1.071** | 1.02 | **1.05** |
| GF8 K=129 R=1, 64 KiB | 1.152 | 1.169 | **1.141** | 1.15 | **1.12** |
| GF8 K=128 R=128, 4 KiB | 0.972 | 1.549 | **1.551** | 1.24 | **1.54** |
| GF8 K=128 R=128, 64 KiB | 1.070 | 1.598 | **1.625** | 1.25 | **1.45** |
| GF8 K=100 R=100, 64 KiB | 1.070 | 1.507 | **1.610** | 1.27 | **1.45** |
| GF8 K=64 R=64, 64 KiB | 1.141 | 1.506 | **1.542** | 1.35 | **1.49** |
| GF8 K=192 R=64, 64 KiB | 1.093 | 1.515 | **1.609** | 1.39 | **1.82** |
| GF8 K=224 R=32, 64 KiB | 1.134 | 1.338 | **1.361** | 1.46 | **1.67** |
| GF8 K=240 R=16, 4 KiB | 1.168 | 1.540 | **1.529** | 3.16 | **4.10** |
| GF8 K=240 R=16, 64 KiB | 1.226 | 1.423 | **1.414** | 3.03 | **3.56** |
| GF16 K=255 R=129, 4 KiB | 1.003 | 1.385 | **1.616** | 2.04 | **2.81** |
| GF16 K=255 R=129, 64 KiB | 0.946 | 1.240 | **1.464** | 1.67 | **2.22** |
| GF16 K=200 R=50, 64 KiB | 1.194 | 1.577 | **1.764** | 3.10 | **4.65** |
| GF16 K=300 R=100, 64 KiB | 1.077 | 1.374 | **1.554** | 2.01 | **2.97** |
| GF16 K=1000 R=200, 4 KiB | 1.091 | 1.555 | **1.781** | 2.28 | **3.64** |
| GF16 K=1000 R=200, 64 KiB | 1.128 | 1.382 | **1.547** | 4.16 | **6.13** |
| GF16 K=1000 R=999, 4 KiB | 1.034 | 1.560 | **1.689** | 1.67 | **2.38** |
| GF16 K=2000 R=500, 4 KiB | 1.097 | 1.520 | **1.729** | 1.97 | **2.91** |
| GF16 K=4096 R=512, 4 KiB | 1.119 | 1.462 | **1.695** | 3.42 | **4.92** |

Every cell where the shipped build still loses to exact main flips except the
two smallest: `K=2 R=2` and `K=4 R=4` at 1 KiB.  Those are not transform-bound.
The whole encode call there is about 70 ns against exact main's 40 ns, so the
residual is Leopard2's fixed per-call cost (buffer validation walks over `K` and
`R`, plan/dispatch selection) rather than field arithmetic, and no kernel change
can close it.

Retained data: `experiments/leopard2/gfni_codec/results/`.

## Radix-four kernel screen

`experiments/leopard2/gfni_codec/fused_radix4_screen.cpp` times one radix-four
inverse group — four shards through the four butterflies
`(v0,v1,log01)`, `(v2,v3,log23)`, `(v0,v2,log02)`, `(v1,v3,log02)` — for each
schedule, and requires all schedules to agree byte for byte before timing.  The
baseline column is the schedule production actually runs: GF8 already fuses
radix-four with nibble tables held in registers
(`AVX2FF8Butterfly4RangePreparedImpl`), while GF16 above 128 bytes uses the
split sweeps.  Microseconds, best of 400, one pinned CPU.

| GF8, bytes/shard | fused nibble (production) | GFNI fused | gain |
| ---: | ---: | ---: | ---: |
| 1024 | 0.060 | 0.030 | 2.00x |
| 4096 | 0.210 | 0.100 | 2.10x |
| 16384 | 0.860 | 0.550 | 1.56x |
| 65536 | 3.460 | 2.090 | 1.66x |
| 262144 | 14.220 | 7.160 | 1.99x |

| GF16, bytes/shard | split nibble (production) | GFNI fused | gain |
| ---: | ---: | ---: | ---: |
| 1024 | 0.090 | 0.050 | 1.80x |
| 4096 | 0.360 | 0.190 | 1.89x |
| 16384 | 1.380 | 0.800 | 1.72x |
| 65536 | 5.610 | 3.040 | 1.85x |
| 262144 | 22.060 | 12.470 | 1.77x |

Three GF16 fused forms were compared.  Re-broadcasting the once-used `log01` and
`log23` matrices from the existing 16-byte table rows inside the tile loop beat
both holding them in registers and using 32-byte replicated memory operands:
`vbroadcasti128` from a loop-invariant L1 address is a single cheap load, while
holding twelve matrices forces data vectors to spill and full memory operands
add a load per affine.  The chosen form therefore needs no table layout change
and no extra memory.

So the kernels themselves are at or near 2x, and end-to-end is 1.3x-1.63x.

## GF8 out-of-place radix-eight staging

A stubbed-operation profile of legacy-high GF8 encode attributes about 40% of
the call to `AVX2FF8Butterfly4Out`, the staging kernel that reads caller source
shards into the work slots, and a further 24% to the in-place radix-four range.
An isolation screen
(`experiments/leopard2/gfni_codec/butterfly4out_floor_screen.cpp`) then showed
that the staging kernel is within 4% of a kernel that only moves the bytes at
64 KiB per shard: it is bound by load/store throughput, not arithmetic.  No
instruction change can help a movement floor.  Moving less can.

An out-of-place radix-eight group carries three transform layers per
load/store round instead of two — eight loads and eight stores for
twenty-four shard-layers against sixteen for two radix-four rounds, a third
less traffic per unit of transform work.  Register pressure is not the obstacle
it appears to be: out-of-place means the eight inputs are consumed into
registers before any output is written, so eight data vectors are live rather
than sixteen, and every multiplier matrix is a memory operand, which costs
nothing on a memory-bound kernel.  The same idea does not pay on the in-place
ALU-bound kernels, which is why it is applied only here.

The isolated screen measured 1.15x at 4 KiB, 1.18x at 16 KiB and 1.26x at
64 KiB against radix-four staging plus one radix-two sweep.  In the codec, with
the selector at `m == 32`:

| Cell | encode | decode |
| --- | ---: | ---: |
| GF8 K=224 R=32, 64 KiB | **1.332x** | 1.129x |
| GF8 K=224 R=32, 16 KiB | 1.094x | 1.001x |
| GF8 K=224 R=32, 4 KiB | 1.109x | 0.994x |
| GF8 K=224 R=32, 1 KiB | 1.084x | 0.989x |
| GF8 K=200 R=32, 64 KiB | 1.133x | 1.032x |
| GF8 K=160 R=32, 64 KiB | 1.146x | 1.133x |

All other measured cells are neutral within run-to-run spread.

The selector is exactly `m == 32` for measured reasons, not caution.  At
`m == 8` a radix-eight group consumes every layer of the transform, so the
accumulating final layer that folds a message block into `xor_result` never
runs and every block after the first is dropped — a randomized cross-build
differential caught precisely that, and it is why the differential runs before
any timing.  At `m == 16` three layers leave only one, which displaces the
fused `dist4 == m` accumulating radix-four with a plain radix-two sweep plus a
separate XOR, measured at 0.88x to 0.94x.  At `m == 128` the staging is a
smaller share of the call and the result is marginal, 1.05x to 1.07x with one
neutral cell.

The operation is optional and published only by the GFNI table; backends that
leave `ff8_ifft_butterfly8_out` null keep the radix-four staging schedule
unchanged.

### Forward radix-eight for the low-rate transform

A stubbed-operation profile of a LOW_V1 encode at `K=32 R=192`, 64 KiB
(573.7 us) attributes 34% to `AVX2FF8Butterfly4Out` and **21% to
`AVX2FF8Butterfly2`** — the in-place radix-two butterfly, which at 2.0 shard
touches per shard-layer is the most expensive form there is, against
radix-four's 1.0 and radix-eight's 0.67.

It is used because of layer parity.  `FFT_DIT_FromCoefficients`
(`LeopardFF8.cpp:2293-2417`) covers `L = log2(m)` layers as a radix-four first
stage plus radix-four rounds plus, when `L` is odd, one trailing radix-two
sweep.  At `m == 32`, `L = 5`: two radix-four rounds leave exactly one layer,
and that sweep runs once per parity block.  A radix-eight first stage covers
three layers, so five become one radix-eight round plus one radix-four round and
**the radix-two disappears entirely**.

`ff8_fft_butterfly8_out` is the forward mirror of the inverse staging
operation.  The group order reverses to largest distance first and the butterfly
becomes `x ^= mul(y, log); y ^= x`.  The skew for a layer at distance `D` inside
a group is `skew[subgroup base + D]`, the same rule the radix-four stage uses,
giving `skew[4d]`, then `skew[2d]` and `skew[6d]`, then `skew[d]`, `skew[3d]`,
`skew[5d]` and `skew[7d]` with `d = m >> 3`.

Measured against the same build with both radix-eight operations absent:

| Cell | encode | decode |
| --- | ---: | ---: |
| GF8 low K=32 R=192, 64 KiB | **1.214x** | 1.003x |
| GF8 low K=32 R=96, 64 KiB | **1.209x** | 0.997x |
| GF8 low K=32 R=224, 64 KiB | **1.188x** | 1.000x |
| GF8 low K=17 R=200, 64 KiB | 1.147x | 1.003x |
| GF8 low K=32 R=192, 16 KiB | 1.144x | 1.000x |
| GF8 low K=32 R=192, 4 KiB | 1.032x | 0.992x |
| GF8 high K=240 R=16, 64 KiB | 1.007x | 1.003x |
| GF8 high K=128 R=128, 64 KiB | 1.008x | 1.016x |

Gated to `m == 32` for the same reason the inverse is: that is the profiled
shape and the one whose layer count leaves no radix-two.  `m == 16` would save
zero touches — `ceil(4/2) == ceil(4/3) == 2` — while still displacing a fused
radix-four, which is exactly how the inverse side regressed to 0.88x.  The
`m == 8` correctness hazard does **not** transfer here, because `FFT_DIT` and
`FFT_DIT_FromCoefficients` have no `xor_result` parameter to strand.

**Status: production.  The GFNI codec ships in default builds as the explicit
`LEO2_BACKEND_GFNI` member; `LEO2_GFNI_VARIANT` without `LEO2_GFNI_MEMBER`
remains the in-place A/B evaluation configuration under the AVX2 identity.**

**Kernel status: measured, kernel published, call site WIRED.**  The operation
`ff8_fft_butterfly8_out` exists, is in the AVX2 GFNI table, and
`FFT_DIT_FromCoefficients` now dispatches to it (`LeopardFF8.cpp:2374`, gated on
`m == 32` and a non-NULL op).  The two contract tests this section predicted
would break were generalized rather than relaxed —
`tests/leopard2/test_direct_encode.cpp:701-715` restates the no-copy contract as
the schedule-independent invariant `2*b2 + 4*b4 + 8*b8 == executed_blocks * p`,
and `tests/leopard2/test_context_backends.cpp:900-913` restates amortization as
`forward_coarse_calls < callsites.fft_dit4 + eight_out*2`.

The original prediction, retained because it is the reason those tests were
restated rather than weakened:

- `leopard2_direct_encode` — "GF8 fused out-of-place call count mismatch".
  `tests/leopard2/test_direct_encode.cpp:701-704` pins
  `fft_butterfly4_out_of_place == executed_blocks * (p / 4)`.  A radix-eight
  first stage issues `m >> 3` calls of a different operation instead of
  `m >> 2` calls of this one, so the count changes.  The no-copy property the
  test exists to protect is *preserved* — the radix-eight operation is also
  out-of-place — but `TestOnlyLowEncodeCounts` has no field for it.
- `leopard2_context_backends` — "low-profile encode did not reduce GF8 forward
  indirect dispatches".  `tests/leopard2/test_context_backends.cpp:977-979`
  requires `ff8_fft_four_range_calls() != 0`.  At `m == 32` the radix-eight
  stage leaves `dist == 1`, which routes to `ff8_fft_butterfly4_out` rather than
  `FFT_DIT4_Range`, so the forward range operation is never called and the
  assertion fails on a true statement about the new schedule.

Landing it therefore requires extending `TestOnlyLowEncodeCounts` with a
radix-eight field, restating the direct-encode expectation in terms of the
schedule actually selected, and restating the context-backends dispatch-reduction
invariant so it measures indirect dispatches rather than assuming the range
operation specifically.  That is a deliberate change to two documented contracts
and should not be done in the same pass as the kernel; it is left unlanded
rather than papered over.

## Production integration requirements

Requirements 1-3 are now satisfied; requirements 4-6 remain open.

1. **Own qualified member — SATISFIED.**  `Leopard2BackendGFNI.cpp` re-includes
   the AVX2 source with `LEO2_GFNI_MEMBER` + `LEO2_GFNI_VARIANT`, compiled
   `-mavx2 -mgfni -mno-avx512f`, publishing `LEO2_BACKEND_GFNI = 6`
   (name `"avx2-gfni"`).  Qualification is the AVX2 gate plus CPUID.(7,0):ECX
   bit 8 (`X86Features.gfni`).  Explicit selection remains exact; the bounded
   AUTO policy described in requirement 3 can qualify and borrow the same
   immutable table without changing the context's reported AVX2 baseline.  The
   startup KAT gained `TestFF8Butterfly8Out`, which checks both fused radix-eight
   kernels against a reference composed from the audited two-point butterfly —
   including the sentinel-skew contract, where the fused ops absorb the skip
   branch that two-point callers perform themselves.  The portable-ISA checker
   carries a `gfni` member class (AVX2 VEX set + `vgf2p8affineqb` +
   four plain-AVX2 mnemonics emitted only by the GFNI-only kernels), with
   mutation-tested leak fixtures; `vgf2p8mulb` is explicitly rejected.
   Evidence: parity byte-identical to AVX2 on randomized shapes; same-binary
   speedup over the AVX2 member 1.331x (GF8 K=224 R=32, matching the
   experimental 1.33x) to 1.732x (GF16 K=200 R=50); full suite
   103/103 including `leopard2_portable_isa` on the GFNI-bearing archive and
   three GFNI fault-injection tests.  Every AVX2-gated policy predicate
   (tiling, R=1 XOR, fused-final, direct encode/repair, sink metadata,
   syndrome fusion) was extended to the new identity explicitly — thirty
   sites — so the member inherits the calibrated AVX2 policies it was
   measured under; the AUTO-AVX512 encode exception and the AVX512-only
   `ff8_high_encode_one_block` gate deliberately remain AVX2/AVX512-scoped.
2. **Table ownership — SATISFIED.**  The GFNI variant now stores
   `FF16AffineTable { uint64_t block[4]; }`, 32 bytes per logarithm (2 MiB);
   vector kernels broadcast each matrix with `vpbroadcastq`, which is
   register-identical to the former duplicated-row broadcast.  Performance is
   unchanged within noise across the measured cells and wire output is
   byte-identical (see
   `experiments/leopard2/optimization_log/24-gf16-affine-table-packing.md`).
   GF8 keeps its 32-byte rows deliberately: 8 KB total, and the radix-eight
   kernels load the duplicated form directly.
3. **Selector policy — SATISFIED FOR ONE EXACT CELL.**  The production selector
   admits only AMD family 1Ah/model 08h, an AUTO AVX2 baseline, native-layout
   flags-zero legacy-high GF16 `K=1000`, `R=200`, `T=256`, one context thread,
   exactly 65,536 bytes, all 200 outputs, and either `leo2_encode` or the
   ordinary one-item batch path without scalable preflight scratch.  Explicit
   backends, decode, scalable-preflight, multi-item/reusable batches, and every
   neighboring identity remain on their context table.  The sealed v2
   same-binary campaign retained 25 target ABBA rounds and 12 inactive cells:
   GFNI/AVX2 was 1.449623x for ordinary encode (95% CI
   `[1.445842, 1.453414]`) and 1.452917x for one-shot encode
   (`[1.448359, 1.457489]`).  All route, call-count, digest, timer-floor, and
   inactive-cell gates passed.  The bound source was commit
   `dbe26daac245739a9692e481622af0a0a077462c`, tree
   `167fd0034673901f015ea5d812604b6c2cbee5f6`; the campaign report SHA-256 is
   `680b42c53e90dad3ec3529fbc048fcb273763d410eb84dd8421100014b8e2233`.
   This qualifies the bounded selector; it is not an exact-Leopard1
   performance claim.
4. **Isolated exact-main evidence.**  The counterbalanced
   `experiments/leopard2/main_compare/run_abba.py` lineage with the CPU-pair
   lease, source/build closure, and the zero-sibling-jiffy gate is still
   required before any claim enters `docs/leopard2_vs_main_benchmark.md`.
   Its frozen schema v16 predates the production-default operation selector
   and cannot distinguish this exact GFNI route from its reported AVX2 context
   baseline.  It must not be reused for current source.  A successor contract
   must attest an actual production-default GFNI call in every Leopard2 child
   (or explicitly pin and attest the selector off for a separately named AVX2
   control) before the exact-main comparison is valid.
5. **Second microarchitecture.**  Everything here is one Zen 5 part.  Intel
   parts with GFNI have different affine throughput and different 512-bit
   frequency behaviour.
6. **512-bit follow-up.**  The kernel screen says the remaining 1.04x to 1.66x
   from widening to `zmm` is separable from the affine change.  That is a much
   larger edit and should be its own candidate.
