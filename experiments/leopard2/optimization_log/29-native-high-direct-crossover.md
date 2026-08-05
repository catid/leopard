# Native-high whole-direct crossover extension

## Result

**Promoted through 7168 bytes.** Commit
`322cdf04d088357218e106456d57a6ea6e3e2f96` widens the allocation-free
whole-direct one-shot route from 256 to 7168 bytes for legacy-high GF8,
pure-AVX2 `N=32,T=8`, `K=9..16`, `R=5..8`, and `L=5..R`.

The previous 256-byte ceiling was a conservative first measurement, not an
arithmetic or safety constraint. The widened route continues to build the
active locator and `K` by `L` direct coefficient matrix in caller scratch and
executes the exact requested byte count without a heap plan. Because direct
execution no longer consumes Algorithm 5 work slots, stale transform-tile and
work-slot predicates were removed from this route.

At 7169 bytes production deliberately returns to the ordinary heap-owned
direct plan. The profile, field, backend, count, loss, flag, diagnostic, alias,
and failure fallbacks are unchanged.

## Crossover discovery

The pinned directional program proceeded in three stages:

1. Seventy cells from 257 through 4095 bytes covered five representative
   count/loss shapes. Whole-direct improved over the preceding production
   route by 1.3556x--7.6111x and over exact Leopard main by
   1.3952x--4.7522x.
2. All 80 valid `(K,R,L)` shapes were swept at 4096, 4097, 4608, 5120, 6144,
   and 7168 bytes. Across 480 selected cells, the minimum same-source speedup
   was 1.1063x and the geometric mean was 1.3325x.
3. The 13 weakest count/loss shapes were swept at every byte from 7105 through
   7168, covering all 64 tail residues. All 832 cells matched digests and
   recovered data; the minimum speedup was 1.1070x.

The byte ceiling is therefore not inferred from one aligned benchmark.

## Frozen same-source gate

Candidate and prior-cap control binaries were built from clean committed
sources, copied to a read-only lane, and hashed before and after timing. Three
ABBA rounds used CPU 28, reserved SMT sibling 12, 31 retained samples, four
warmups, reuse 200, one codec thread, and the canonical campaign lock. The
reserved sibling accumulated zero non-idle jiffies.

Selected cells improved by 1.1240x--4.9723x; the weakest lower 95-percent
confidence bound was 1.1106x. Representative results are:

| Cell | Prior us | Candidate us | Speedup | 95% CI |
| --- | ---: | ---: | ---: | ---: |
| K=9 R=5 L=5, 257 B | 1.3937 | 0.2803 | 4.9723x | 4.8490--5.0987 |
| K=12 R=7 L=5, 4096 B | 3.2469 | 2.4996 | 1.2989x | 1.2831--1.3150 |
| K=16 R=8 L=5, 7136 B | 5.9949 | 5.3288 | 1.1250x | 1.1139--1.1362 |
| K=16 R=8 L=5, 7168 B | 6.0312 | 5.3658 | 1.1240x | 1.1106--1.1375 |
| K=16 R=8 L=8, 7168 B | 10.7954 | 8.0655 | 1.3385x | 1.3364--1.3405 |

The 256-byte lower inert control was 1.0002x with a confidence interval
spanning one. The 7169-byte upper inert control was 0.9877x, a 1.23-percent
point difference with a confidence interval of 0.9506--1.0264; it remains
below the two-percent point-regression gate and the new route is not selected
there.

## Exact Leopard-main gate

The exact-main replay used the same frozen candidate and Leopard main commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. A clean repeat recorded zero
reserved-sibling non-idle jiffies. Leopard main over Leopard2 ratios ranged
from 1.8856x to 5.1333x, and the weakest lower 95-percent bound was 1.8811x.
Exact main physically rounds nonmultiple byte counts to 64 bytes with a zero
suffix; both programs fingerprinted the same requested logical prefix.

## Correctness and memory-safety gates

The focused Release tests now cover:

- 12,544 allocation-audited public one-shot cases;
- 11,886 independent-generator and reusable-Algorithm-5 differentials;
- every valid count/loss shape at selected power, tail, crossover, and maximum
  boundaries through 7168 bytes;
- 7169-byte fallback, mixed parity availability, unaligned buffers, exact
  scratch, failure atomicity at the maximum, and zero eligible allocations;
- Clang 18 ASan+UBSan; and
- GF8-only and GF16-only library builds.

Builds and tests were serialized with one compiler job and explicit memory
caps to avoid the host OOM failures seen earlier in the campaign. ASan runtime
was not virtual-address capped because its shadow map requires a large address
space.

## Source-major continuation

A second prototype kept the cheap raw locator setup but switched to the mature
source-major direct loop above the output-major crossover. It is promising,
but not promoted here. Across all 80 shapes its minimum same-source speedup
fell below five percent at 24 KiB, and seed/tail cells approached parity near
8 KiB. At 32 KiB, six of 80 cells were below the five-percent gate; by 64 KiB,
75 of 80 were. A padded 32-byte tail reduced some ragged overhead but did not
remove the heterogeneous transition.

That prototype needs its own pattern-aware or shared-helper design and frozen
neighbor campaign. No source-major prototype code remains in production.

Complete identities, hashes, all frozen ABBA cells, discovery summaries, test
counts, and the source-major negative/continuation data are in
`experiments/leopard2/direct_repair/results/raw_native_high_direct_crossover_20260805.json`.
