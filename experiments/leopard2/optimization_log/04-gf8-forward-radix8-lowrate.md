# 04 — GF8 forward radix-eight for the low-rate transform

**Disposition: LANDED** (GFNI variant, gated `m == 32`)

## Idea
A stubbed-operation profile of a `LOW_V1` encode at `K=32 R=192`, 64 KiB
(573.7 us) attributes 34% to `AVX2FF8Butterfly4Out` and **21% to
`AVX2FF8Butterfly2`** — the in-place radix-two butterfly, the most expensive
form at 2.0 shard touches per shard-layer against radix-four's 1.0 and
radix-eight's 0.67.

It is used because of layer parity. `FFT_DIT_FromCoefficients` covers
`L = log2(m)` layers as a radix-four first stage plus radix-four rounds plus,
when `L` is odd, one trailing radix-two sweep. At `m == 32`, `L = 5`: two
radix-four rounds leave exactly one layer, and that sweep runs once per parity
block. A radix-eight first stage covers three layers, so five become one
radix-eight round plus one radix-four round and **the radix-two disappears**.

## Result, final — quiet machine, versus the same build without radix-eight

| Cell | encode | decode |
| --- | ---: | ---: |
| GF8 low K=32 R=192, 64 KiB | **1.233x** | 1.001x |
| GF8 low K=32 R=96, 64 KiB | **1.215x** | 0.997x |
| GF8 low K=32 R=224, 64 KiB | **1.212x** | 0.998x |
| GF8 low K=17 R=200, 64 KiB | 1.158x | 1.005x |
| GF8 low K=32 R=192, 16 KiB | 1.139x | 1.001x |
| GF8 low K=32 R=192, 4 KiB | 1.033x | 1.001x |
| GF8 high K=240 R=16, 64 KiB | 1.000x | 1.004x |
| GF8 high K=128 R=128, 64 KiB | 1.004x | 0.998x |

**No leopard1 column is possible.** leopard1 rejects `R > K`, so it has no
low-rate profile to compare against. These are Leopard2-versus-Leopard2.

Two earlier timing runs were discarded: both started while a ctest suite was
running and showed base times inflated 2-4x (K=32 R=192 base 2095 us against
536 us here). The lesson is recorded in the README — gate a measurement on the
test log printing its summary line, not on a bounded process poll, which times
out and proceeds.

## Earlier result (pre-landing, same conclusion)
| Cell | encode | decode |
| --- | ---: | ---: |
| GF8 low K=32 R=192, 64 KiB | **1.214x** | 1.003x |
| GF8 low K=32 R=96, 64 KiB | **1.209x** | 0.997x |
| GF8 low K=32 R=224, 64 KiB | **1.188x** | 1.000x |
| GF8 low K=17 R=200, 64 KiB | 1.147x | 1.003x |
| GF8 low K=32 R=192, 16 KiB | 1.144x | 1.000x |
| GF8 low K=32 R=192, 4 KiB | 1.032x | 0.992x |
| GF8 high K=240 R=16, 64 KiB | 1.007x | 1.003x |

**No leopard1 column is possible.** leopard1 rejects `R > K`, so it has no
low-rate profile to compare against. These are Leopard2-versus-Leopard2.

Wire-identical over 250 randomized shapes plus 8 directed low-profile shapes
including ragged and 64-byte cases.

## The test contracts it changed
It changes the low encoder's operation mix, and two contract tests correctly
detected that. Both are now restated. **All three restatements reduce exactly to
the original assertions when no radix-eight operation is published**, so nothing
is weakened for existing backends.

1. `tests/leopard2/test_direct_encode.cpp` — **RESOLVED.** It pinned
   `fft_butterfly4_out_of_place == blocks * (p/4)`. The property being protected
   is the no-copy contract: every coefficient enters the transform exactly once
   through an out-of-place first layer. That is now asserted directly, in
   coefficients, as
   `2*b2 + 4*b4 + 8*b8 == executed_blocks * p`, which is schedule-independent
   and strictly stronger. A `fft_butterfly8_out_of_place` counter was added to
   `TestOnlyLowEncodeCounts`. This change is kept.
2. `tests/leopard2/test_context_backends.cpp` — **RESOLVED on the fourth
   attempt**, by deriving the counter semantics instead of guessing them.
   `TestFFTDIT4Calls` does `fetch_add(dist)` immediately before dispatching, so
   `callsites.fft_dit4` counts **coordinate groups** while
   `ff8_fft_four_range_calls` counts **dispatches**. The property is that one
   dispatch amortizes over many groups. Three consequences:
   - radix-eight is not an `fft_dit4` callsite, so it needs its own counter and
     must not be added to `ff8_fft_four_calls` (that equality proves the
     context's own table is called);
   - the denominator must count radix-four-**equivalent** groups covered,
     `callsites.fft_dit4 + 2 * eight_out_calls`, because one radix-eight
     dispatch covers eight coordinates;
   - the separate no-copy check must assert "the first layer is out of place",
     whatever its radix, rather than naming the radix-four operation.

   Historical record of what did not work: It requires
   `ff8_fft_four_range_calls() != 0 && < callsites.fft_dit4`. Two attempts
   failed:
   - Counting radix-eight into `ff8_fft_four_calls` broke a separate equality
     (`ff8_fft_four_calls() == callsites.fft_dit4`) that exists to prove the
     context's own table is being called. Radix-eight is not an `fft_dit4`
     callsite, so it must not be counted as one.
   - Adding a dedicated `ff8_fft_eight_out_calls` to the numerator still fails,
     because the **denominator** `callsites.fft_dit4` is itself reduced by
     radix-eight — the first stage no longer issues those callsites at all.

   The fix is a denominator expressed in **coordinate groups covered** rather
   than radix-four callsites, so that "coarse dispatches < work units" stays
   meaningful when a coarser operation exists. That is a deliberate redesign of
   what the invariant means and should be done on its own, not wedged in beside
   a kernel change — a test that passes because it stopped checking is worse
   than no test.

## To revisit
Kernel `AVX2FF8FFTButterfly8Out` and op `ff8_fft_butterfly8_out` are in the tree
and in the AVX2 GFNI table. Skews, verified against the radix-four convention
(`skew[subgroup base + D]`), are `skew[4d]`; `skew[2d]`, `skew[6d]`; then
`skew[d]`, `skew[3d]`, `skew[5d]`, `skew[7d]` with `d = m >> 3`. Re-wiring is
~30 lines in `FFT_DIT_FromCoefficients`. Only the context-backends invariant
stands in the way.

`m == 16` must stay excluded for a different reason than the inverse side:
`ceil(4/2) == ceil(4/3) == 2`, so radix-eight saves **zero** touches there while
still displacing a fused radix-four. The `m == 8` correctness hazard does not
transfer, because the forward path has no `xor_result` to strand.
