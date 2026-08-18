# 37 — GF8 T=32 Q3 packed B64 boundary

## Hypothesis

The post-fix Leopard1 atlas showed a small but credible GF8 encode deficit at
`K=77..79`, `R=32`, and 64-byte shards.  A frozen 25-round diagnostic put
Leopard1/Leopard2 at 0.974x for K79.  The transform itself was not the obvious
problem: these cells already use Leopard's redundancy-side `T=32` high-rate
encoder.  At only 64 bytes per shard, the public Leopard2 path's pointer,
range, and sparse-plan staging was large relative to the byte arithmetic.

The existing packed B64 terminal already removed that overhead for
`K=33..64`.  Extending its aggregate range proof through the exact three-block
`K=65..96,R=32` band should recover the fixed overhead without introducing a
new transform, changing the wire profile, or perturbing neighboring shapes.

## Implementation

The selected range is exactly:

- legacy-high GF8 with native shard layout and the effective AVX2 backend;
- `K=65..96`, `R=32`, and exactly 64 bytes per shard;
- packed, contiguous source and output slabs with ordinary AUTO routing;
- public one-shot encode and the one-item batch wrapper.

The wire construction remains parent `[N,D]=[128,96]`, with `T=32`, parity
coordinates `0..31`, message coordinates beginning at 32, and shortened zero
coordinates ending at 127.  The terminal validates aggregate input, output,
scratch, pointer-array, and optional batch-descriptor spans, derives its 64
work pointers from the already queried `2T` scratch rows, and then calls the
unchanged mature `ReedSolomonEncode` transform.  It does not stage the generic
per-shard address metadata or build the dense sparse-plan descriptor.

Sparse, detached, overlapping, non-64-byte, non-AVX2, GF16, low-profile,
`R!=32`, and `K>96` calls retain the mature path.  Multi-item batch and reusable
batch bindings are also intentionally unchanged.  A separately derived fused
T32/Q3 circuit was therefore not promoted: the simpler public-boundary change
already cleared the performance threshold while keeping the arithmetic shared.

## Correctness and safety

The diagnostic suite sweeps every `K=65..96` at `R=32,B=64` twice: once with
dense deterministic data and once with only the final live source nonzero.
Every parity byte matches both the independent direct systematic generator and
legacy `leo_encode`.  The production-linked suite checks K65, K79, and K96,
including unaligned packed buffers, short/misaligned scratch, scratch/input,
scratch/output, source/output, and metadata overlap.  It poisons the scratch
metadata prefix and proves the production terminal leaves it untouched, which
distinguishes the selected route from a correct generic fallback.

Additional gates passed:

- release focused tests: 2/2 in 105.09 seconds;
- GF8-only focused test: 1/1 in 91.43 seconds;
- GCC 13 and Clang 18 strict builds for both focused targets;
- Clang 18 ASan+UBSan+LSan: 2/2 in 141.00 seconds;
- operation-count/source audit: 489/489 checks;
- project graph audit: 171/171 normal and 171/171 under `python -O`;
- hardened runner contract tests under normal Python and `python -O`.

## Performance

The accepted candidate is commit `1ef58ad17fa8f30cd8868a323fd8374929bdcaf8`,
executable SHA-256
`4305ce3a1025f0e565d37fc3a16478eeaa10800125e187ac1a2222bb9d8190ca`.
The exact Leopard1 reference is commit `6e5725eb`, executable SHA-256
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.
Ratios above 1.0 favor the candidate.

The runner used 25 independent balanced rounds per cell, 31 inner samples,
64 warmups, reuse 8192, CPU 13 with sibling 29 reserved, one shared immutable
candidate/control inode, and a runner-owned reproducible-build proof.  It
retained 1,850 accepted processes and discarded 36 processes from six
contaminated round attempts before retrying them.

| K | one-item batch: mature/candidate | one-item batch: Leopard1/candidate | one-shot: mature/candidate | one-shot: Leopard1/candidate |
| ---: | ---: | ---: | ---: | ---: |
| 65 | 1.146 [1.133,1.159] | 1.142 [1.133,1.152] | 1.115 [1.104,1.126] | 1.103 [1.095,1.112] |
| 76 | 1.136 [1.129,1.143] | 1.137 [1.129,1.146] | 1.112 [1.106,1.118] | 1.106 [1.098,1.114] |
| 77 | 1.131 [1.123,1.139] | 1.126 [1.116,1.136] | 1.110 [1.103,1.117] | 1.094 [1.085,1.104] |
| 78 | 1.128 [1.121,1.135] | 1.125 [1.118,1.131] | 1.108 [1.100,1.116] | 1.095 [1.087,1.102] |
| 79 | 1.132 [1.123,1.142] | 1.115 [1.106,1.125] | 1.113 [1.105,1.121] | 1.084 [1.075,1.093] |
| 80 | 1.125 [1.118,1.131] | 1.132 [1.124,1.139] | 1.106 [1.101,1.110] | 1.098 [1.091,1.106] |
| 81 | 1.126 [1.121,1.131] | 1.126 [1.118,1.135] | 1.105 [1.100,1.110] | 1.097 [1.090,1.105] |
| 95 | 1.103 [1.097,1.109] | 1.098 [1.092,1.104] | 1.084 [1.079,1.090] | 1.067 [1.062,1.072] |
| 96 | 1.104 [1.098,1.111] | 1.121 [1.113,1.130] | 1.081 [1.075,1.086] | 1.092 [1.083,1.100] |

Every inactive candidate/control interval was wholly inside the predeclared
`[1/1.02,1.02]` equivalence band:

| Inactive cell | one-item batch CI95 | one-shot CI95 |
| --- | ---: | ---: |
| K79/R31/B64 (Q3) | [0.996,1.006] | [0.995,1.007] |
| K79/R33/B64 (Q2) | [0.990,1.003] | [0.991,1.003] |
| K79/R32/B63 (Q3) | [0.998,1.003] | [0.996,1.002] |
| K79/R32/B65 (Q3) | [0.996,1.004] | [0.998,1.005] |
| K97/R32/B64 (Q4) | [0.993,1.003] | [0.993,1.002] |

The complete machine-readable projection is
`results/t32_q3_b64_boundary_checkpoint_20260818.json`.  Its ignored raw and
summary inputs have SHA-256
`9b8ac5d89934868ed18b2a6d41aff2388a299e6b285ba5dbaeb4aa45c7268a65`
and `7e082601154820f8bb635d78b6ac764f4d1d42dc58a4eb7872875afb8dfae83e`.

Two launches were rejected before timing because their CMake cache types did
not match the canonical production contract: first the compiler was recorded
as `STRING` rather than `FILEPATH`, then AVX-512 probe values were recorded as
`INTERNAL` rather than explicit `UNINITIALIZED=FALSE`.  No timing conclusion
was retained from either attempt.

## Result

Landed for the exact packed GF8/AVX2 B64 region.  The original K77–79 deficit
is eliminated: Leopard2 is now 1.115–1.126x faster than exact Leopard1 for
one-item batch execution and 1.084–1.095x faster for one-shot encode in those
three focus cells.  The win comes from removing fixed public-boundary staging;
the paper-derived redundancy-side arithmetic and legacy wire bytes remain
unchanged.
