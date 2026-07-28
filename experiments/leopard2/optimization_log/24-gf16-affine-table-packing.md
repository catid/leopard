# 24 — GF16 affine table packing: 8 MiB -> 2 MiB (GFNI requirement 2)

**Disposition: SHIPPED — 4x table footprint reduction, performance unchanged
to slightly better, wire output identical**

## The waste
The GFNI evaluation deliberately reused the nibble-table storage shape so no
vector call site changed: `FF16NibbleTable` is 128 bytes per logarithm
(`low[4][16]` + `high[4][16]`), but the GFNI variant stores each of the four
8-byte affine matrix blocks *duplicated* within a 16-byte row and writes the
`high` half as a byte-for-byte copy of `low`.  8 MiB allocated to hold 2 MiB
of matrices — flagged as production requirement 2 in the GFNI doc from day
one, and tracked as `leopard-79h.45`.

## The change
A GFNI-variant-only table type, selected where the variant macro already
splits the code:

    struct FF16AffineTable { uint64_t block[4]; };   // 32 B/log, 2 MiB total

The equivalence argument is short and exact: `GFNIStoreMatrix` wrote byte `b`
of the matrix to row byte `b` (then duplicated), and the old preamble
broadcast the 16-byte row with `vbroadcasti128` — matrix repeated four times
across the 64-bit lanes.  `_mm256_set1_epi64x(block[i])` (vpbroadcastq)
produces the identical register image from 8 bytes, so `VGF2P8AFFINEQB` sees
bit-identical operands.  Scalar tails get `GFNIApplyMatrix64`, the same
parity loop reading matrix bytes by shift instead of by row index.

Implementation was type-driven: splitting the table type by variant let the
compiler enumerate every consumer.  Eight sites converted (scalar product,
four multiply/multiply-add kernel preambles, the fused butterfly4 helper, the
zero-skew butterfly2 preamble, and the publish loop); the stock nibble build
inside the shared initializer is compiled out of GFNI variants (it was
unreachable after the preprocessor-selected early return but still had to
type-check).  `high_tables` disappears from the GFNI data path entirely — it
was already unused by `AVX2FF16ProductVectors` in this variant.

## Measured
Same-binary, packed GFNI member vs AVX2 member (fresh build, pinned core):

| Cell | AVX2 | GFNI packed | speedup | previous (16B rows) |
| --- | ---: | ---: | ---: | ---: |
| GF16 K=200 R=50, 64 KiB | 1337.1 us | 771.0 us | **1.734x** | 1.732x |
| GF16 K=1000 R=200, 64 KiB | 9112.0 us | 5925.3 us | **1.538x** | 1.514x |
| GF16 K=2000 R=500, 256 KiB | 91903.8 us | 64852.3 us | **1.417x** | 1.457x |
| GF8 K=192 R=64, 64 KiB | 706.3 us | 500.6 us | **1.411x** | 1.433x |
| GF8 K=224 R=32, 64 KiB | 809.4 us | 605.9 us | **1.336x** | 1.331x |

Within noise of the unpacked form throughout — as expected, since a transform
touches only a few hundred distinct logarithms and the win is footprint, not
throughput.  Table state now reports `ff16_bytes = 2097152`; the
backend-failure accounting gained a GFNI-specific expectation (the packed
member deliberately no longer matches the SIMD nibble pin).

GF8 tables keep the 32-byte duplicated-row shape: they are 8 KB total, and
the fused radix-eight kernels load 32-byte rows relying on that duplication —
repacking them would touch measured kernels for a three-orders-of-magnitude
smaller prize.

## Validation
- Full suite at the packed build: 102/103 with one real failure that was the
  drift contract doing its job — `leopard2_balanced_promotion_plan_self_test`
  mirrors `run_abba.py`'s translation-unit list and detected that the earlier
  campaign-model fix updated the producer but not this consumer.  Mirror
  updated (AVX2Xor + GFNI, count 20 -> 22); test green; every other test
  passed including the GFNI KAT against the packed tables at startup.
- **Randomized backend differential: 200 shapes, 0 mismatches** (avx2 vs
  packed-gfni in one binary; 42 shapes rejected identically).
- Parity byte-identical to AVX2 on the API smoke; backend 7 still
  `LEO2_INVALID_ARGUMENT`.
- `TestGetGFNITableState` reports 2,097,152 GF16 bytes; the backend-failure
  accounting asserts it GFNI-specifically.
