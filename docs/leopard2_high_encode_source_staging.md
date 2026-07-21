# High-rate encoder source staging

## Result

The legacy-high and low-profile encoders normally no longer copy a complete active
input prefix into scratch before beginning a block inverse transform.  For a
transform side greater than four, each complete group of four caller shards is
read directly by the selected scalar, SSSE3, AVX2, or explicit AVX512 inverse
radix-four backend and the two-layer result is written to scratch.  A final
group of one to three shards retains the established copy-and-zero path.  Sides of four or
less retain the old path because their only inverse stage can already be fused
with high-rate XOR accumulation.  A measured GF16 AVX2 crossover also retains
copy-first for every legacy-high block when the transform side is at least 256
and the current aligned byte pass exceeds 16 KiB.  Low-profile interpolation
still uses direct source staging in that region.

This is an execution-kernel change only.  It does not change the field,
coordinates, parent size, shortening, puncturing, parity order, or arithmetic
order within a four-way butterfly, so it cannot select a different wire
profile.

## Transform identity

For a distance-one inverse group beginning at coordinate `r`, the existing
radix-four schedule applies the two radix-two layers with multiplier logs

    log01 = FFTSkew[r + 1]
    log23 = FFTSkew[r + 3]
    log02 = FFTSkew[r + 2]

The out-of-place backend consumes the four immutable source rows using exactly
those logs and writes the same four rows that the in-place
`IFFT_DIT4_Range(..., dist = 1)` call produced after the old copy.  Subsequent
layers therefore resume at distance four (radix-four distance sixteen).  When
the transform has an odd number of radix-two layers, the established final
distance-half layer and its XOR-accumulating form are unchanged.

For each complete group, the old sequence moved the active bytes once during
the copy and once during the first two inverse layers.  The new sequence reads
the caller bytes once and writes the transformed scratch rows once.  Thus it
removes one active-input read/write memory pass without changing the number of
field butterflies.  The input pointers remain read-only.

The inactive suffix is still initialized to zero because later regular layers
may combine it with the active prefix.  Eliminating those writes requires a
true sparse-input inverse schedule, not an assumption that caller scratch was
previously zero.  That follow-up is tracked as
`leopard-79h.18.1.12.2` and is not claimed by this change.

## Boundary behavior

- A complete four-shard group performs no input-shard copy.
- A ragged final group copies only its one to three live rows; its remaining
  rows were already initialized as shortened zeros.
- A completely inactive group is not transformed.
- The GF16 AVX2 large-block crossover is decided once per pass from backend,
  profile family, transform side, and byte count.  It is not an online
  benchmark and cannot alter the wire profile.
- GF16 uses the same padded physical-byte contract as the rest of the GF16
  transform path; no odd-byte projection is introduced here.
- The caller's selected context backend is used.  The process-global default
  table is not consulted by the new calls.
- Test-only counters report the number of direct inverse radix-four groups and
  the number of live input shards that still required a copy.

For `(K,R)=(32,16)`, encoding executes eight direct groups and copies no input
shards.  For `(33,16)`, it executes the same eight complete groups and copies
only the one-shard final block.  Both identities are checked independently in
GF8 and GF16.

## Correctness and safety evidence

The focused tests compare every requested parity byte with the independent
direct systematic generator, retain a snapshot of every caller source, and
trace the selected backend table.  The integrated validation at the candidate
checkpoint included:

- the complete Release suite, 73 of 73 tests;
- focused strict GCC 13 and Clang 18 suites, 8 of 8 each;
- Clang 18 ASan plus UBSan, 8 of 8 focused tests;
- no-OpenMP GCC ThreadSanitizer under the host's required `setarch -R`
  invocation, 4 of 4 concurrency/API tests;
- legacy golden parity, GF16 legacy matrices and padded odd payloads;
- scalar, SSSE3, and AVX2 context routing; and
- arbitrary vector tails, unaligned guarded buffers, repeated encoding, and
  concurrent immutable-codec execution.

## Performance disposition

The production candidate and its immediate clean control were measured as
separate test-hook-free binaries in counterbalanced process order on CPU 14,
with SMT sibling 30 reserved.  Six process medians per provider were retained
for the main matrix.  Sixteen invocations with observed sibling activity were
rejected rather than folded into the result; 464 accepted JSON documents remain
under `.research/leopard2/high-source-staging/`.  Every accepted pair had
identical original, parity, and recovered-output digests.

Representative control-time/candidate-time geometric means before applying
the deterministic fallback were:

| Field/backend and shape | 64 B | 4 KiB | 64 KiB |
| --- | ---: | ---: | ---: |
| GF8 AVX2, K=240 R=16 | 1.086x | 1.251x | 1.288x |
| GF16 AVX2, K=4096 R=16 | 1.170x | 1.117x | 1.007x |
| GF16 scalar, K=4096 R=16 | 1.073-1.132x across the range | | |
| GF16 SSSE3, K=4096 R=16 | 1.098-1.148x across the range | | |

GF16 K=1000/R=200 at 4 KiB improved by 1.029x scalar, 1.042x SSSE3,
and 1.047x AVX2.  The AVX2 crossover for the same T=256 shape was 1.054x
at 8 KiB, 1.031x at 16 KiB, 0.993x at 32 KiB, and 0.917x at 48 KiB; at
64 KiB it was 0.881x.  A two-block K=512/R=256 case similarly fell to
0.913x at 64 KiB, while the one-block K=256/R=256 case remained approximately
neutral.

An initial deterministic fallback applied copy-first only to later accumulated
blocks.  A fresh test-hook-free ABBA gate rejected that policy: K=1000/R=200
was still 0.956x at 48 KiB and 0.963x at 64 KiB, and K=512/R=256 was 0.912x at
64 KiB.  All digests matched.  The final policy therefore returns every
legacy-high block, including the first, to copy-first in the qualified region.
Scalar, SSSE3, GF8, low-profile interpolation, smaller sides, and passes no
larger than 16 KiB retain the measured source-staging path.

A later GF8 follow-up replaces the legacy copy-first accumulator specifically
for qualified `T=2` and `T=4` AVX2 cells.  It preserves the same wire profile
while feeding complete message blocks directly into inverse-transform XOR
accumulation.  Its conservative selector, exact-main results, and explicit
negative cells are recorded under
`experiments/leopard2/gf8_high_encode/results/8c54abd-small-t-amd9950x3d/`.

The rebuilt final policy then passed 56 of 56 isolated ABBA launches with no
rejections and exact original, parity, and recovery digests. Control-time over
candidate-time geometric means were 1.023x for K=1000/R=200 at 16 KiB, 0.997x
at 48 KiB, and 1.002x at 64 KiB; K=512/R=256 at 64 KiB was 1.005x, and
K=4096/R=512 at 64 KiB was 0.999x. Thus every qualified large fallback cell
was within 0.5% of the mature path by geometric mean. Retained source-staging
cells improved by 1.122x for GF16 K=4096/R=16 at 4 KiB and 1.268x for GF8
K=240/R=16 at 4 KiB. The bound manifest and summary are under
`.research/leopard2/high-source-staging-final-policy/` in the ignored research
cache.

## Exact Leopard-main gate

Commit `a1fa5bf` was rebuilt from a clean detached, test-hooks-off tree and
compared with the independently linked exact Leopard-main commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  To preserve the runner's strict
zero-work rule on the reserved SMT sibling, each representative cell was a
separate three-round ABBA bundle.  All eight bundles passed source, build,
runtime, parity, recovered-data, affinity, and retained-result verification;
every sibling non-idle delta was zero.  Ratios are Leopard main time divided by
Leopard2 time, so values above one favor Leopard2:

| Cell | Encode speedup | 95% confidence interval |
| --- | ---: | ---: |
| GF8 XOR, K=129 R=1 | 0.975x | 0.958-0.992 |
| GF8 high, K=240 R=16, one loss | 1.172x | 1.164-1.181 |
| GF8 high, K=240 R=16, full loss | 1.194x | 1.172-1.217 |
| GF8 balanced, K=128 R=128 | 1.037x | 1.035-1.039 |
| GF16 inflation, K=200 R=50 | 1.023x | 1.020-1.027 |
| GF16 high, K=1000 R=200, one loss | 0.894x | 0.888-0.900 |
| GF16 high, K=1000 R=200, full loss | 0.895x | 0.890-0.900 |
| GF16 large, K=4096 R=512 | 0.949x | 0.939-0.959 |

The source-staging change closes the prior GF8 high-rate deficit and makes the
balanced and field-inflation cells modestly faster.  It does not close the
remaining GF16 encode gap: the measured large-shard AVX2 fallback is necessary
to avoid an even larger regression, but its mature copy-first path remains
about 5-11% behind exact main in these cells.  That residual work remains under
the parent high-encode optimization bead rather than being hidden by this
milestone.  The verified manifests are retained under the ignored research
cache at
`.research/leopard2/exact-main-a1fa5bf-cpu4-cells-20260718/`.
