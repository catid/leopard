# Default-off T32/B256 generated AVX2 encoder

Status: static candidate only.  It has not been built, disassembled, tested, or
timed.  The CMake option `LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED` is OFF by
default.

## Gap and attribution target

The source-attested all-K diagnostic at Leopard2 `97ee49a` measured exact
legacy-high GF8 K=R=T=32, B=256 at exact-main/Leopard2 ratios 0.8035 (loss-1
row) and 0.8327 (loss-32 row).  Loss count does not change the encode workload;
the two rows expose campaign noise.  The later broad packed-validator prototype
recovered 1.0906x same-source but remained 0.9676x versus exact main.  Therefore
this arithmetic candidate must be timed on top of that packed route, not be
credited with validation savings it did not create.

## Independent stage derivation

For T=32, the inverse DIT order is distance 1, 2, 4, 8, 16 and the forward DIT
order is distance 16, 8, 4, 2, 1.  Thus one complete encode is exactly:

1. Four contiguous inverse radix-8 groups (distances 1/2/4).
2. Eight four-row columns at offsets c,c+8,c+16,c+24, each receiving the
   remaining inverse radix-4 (distances 8/16) immediately followed by the
   leading forward radix-4 (distances 16/8).
3. Four contiguous forward radix-8 groups (distances 4/2/1).

The mature decomposition stores after both outer radix-4 groups.  Keeping one
four-row column live across their boundary removes one full 32-row load/store
round without commuting any butterflies within either transform.  At B=256
the removed traffic is 32 * 256 bytes read plus 32 * 256 bytes written, or
16,384 bytes per encode.  The middle live set is four data vectors, not eight.

The canonical inverse skew slice `FFTSkewStorage + 32` is:

    255,196,219,76,153,54,7,99,
    17,19,111,49,102,67,28,52,
    85,137,183,226,51,26,224,108,
    34,27,131,134,187,38,222,139

The canonical forward slice `FFTSkewStorage` is:

    0,255,255,85,255,17,85,34,
    255,153,17,102,85,51,34,187,
    255,219,153,7,17,111,102,28,
    85,183,51,224,34,131,187,222

For an inverse radix-8 group at base b, constants are selected in order
`b+1,b+3,b+5,b+7,b+2,b+6,b+4`.  This yields:

    b=0:  196,76,54,99,219,7,153
    b=8:  19,49,67,52,111,28,102
    b=16: 137,226,26,108,183,224,51
    b=24: 27,134,38,139,131,222,187

The fused outer inverse radix-4 uses indices 8,24,16 and therefore constants
`17,34,85`.  The outer forward radix-4 uses the same indices in the forward
slice and therefore `255,85,255`.

For a forward radix-8 group the same seven index values are supplied, but the
implementation executes them in distance order 4,2,1.  The four tuples are:

    b=0:  255,85,17,34,255,85,255
    b=8:  153,102,51,187,17,34,85
    b=16: 219,7,111,28,153,102,17
    b=24: 183,224,131,222,51,187,34

The candidate reads immutable nibble tables through `GetAVX2FF8Tables()` from
a separate AVX2 object.  It neither extends `Ops` nor changes the mature AVX2
translation unit.  Its exact-shape route follows the shared aggregate packed
validator and transform-backend selection, before the mature scratch pointer
schedule is materialized.  The public routing branch is preprocessed away when
the experiment is OFF.  ON builds use a volatile selector in the isolated
object so candidate and same-source control retain identical executable text;
the control changes only the selector's initialized data word and falls
through to the mature transform.

## Required qualification before timing

Correctness kill gates:

- Compare all 32 systematic-coordinate basis messages and deterministic random
  stripes against both the independent direct GF8 systematic generator and
  legacy Leopard parity.
- Compare with the mature same-binary transform, including unaligned packed
  slabs and ordinary plus one-item batch entry points.
- Prove B255/B257, K31/R32, K32/R31, sparse/detached/reordered outputs, low
  profile, GF16, scalar, SSSE3, AVX-512, and GFNI never enter the candidate.
- Exercise scratch, metadata, input/output overlap, duplicate output, null
  output, guard, and failure-atomicity cases before any candidate write.
- Build with GCC and Clang warnings, ASan+UBSan, and the portable field-reduced
  configuration.

Static/assembly kill gates:

- With the option OFF, compare compile commands, preprocessed production
  sources, mature object hashes/sections, and hot-symbol bytes against the
  parent commit.  The only new source must be absent from the build graph.
- With it ON, reject EVEX/ZMM/opmask instructions and any vector spill/reload in
  the generated entry.  Reject generated text above 8 KiB unless a tighter
  measured budget is established.

Predeclared timing cells:

- Primary: exact K32/R32/B256/q1 at losses 1,2,4,8,16,32, candidate versus
  packed-validator mature control and exact main, nine-round mirrored ABBA on
  a pinned physical core with its SMT sibling idle and frozen executable
  hashes.
- Batch control: B256/q8 is a neighbor, because the generated route is
  deliberately inert for the prevalidated multi-item batch entry point.
- Byte controls: B255, B257, B64, B1024.  Exact main receives B255/B257 via
  its source-attested zero-padding adapter (`--bytes 256/320` plus
  `--logical-bytes 255/257`), and validation requires the v2 schema and common
  logical-prefix digests.
- Shape controls at B256: K31/R32, K32/R31, K31/R31, K33/R32, K32/R33, plus
  existing K16/R8.  The two R>K shapes are same-source-only because Leopard1's
  public API rejects R>K; every other shape also retains exact-main evidence.
- Side controls: K16/R16/B256 and K64/R64/B256.

Promotion requires lower-95% speedup at least 1.05x against both packed mature
control and exact main, with no neighboring lower-95% ratio below 0.98.  Any
parity mismatch, wrong route, spill, EVEX instruction, or larger regression
kills the candidate.
