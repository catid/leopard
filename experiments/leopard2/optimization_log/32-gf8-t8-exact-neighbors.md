# GF8 T8 exact-neighbor encoders

## Result

**Promoted for five adjacent cells.** Commit
`96cda006fadfb16ad66bcadd525f76f4724bd508` closes the measured AVX2 gaps
around the exact `K=12,R=8,T=8` encoder without changing the legacy wire
profile:

- `K=11,R=8` at 256 and 1024 bytes;
- `K=12,R=7` at 256 and 1024 bytes; and
- `K=13,R=8` at 256 bytes.

Public loss-free encoding is 1.111x to 1.282x faster than exact Leopard main
commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. Reusable binding execution is
1.129x to 1.373x faster. `K=13,R=8,B=1024` did not clear the gate and remains
on the mature fallback.

| K | R | bytes | main | Leopard2 one-shot | main / Leopard2, 95% CI |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 11 | 8 | 256 | 0.13835 us | 0.11007 us | **1.2570x** `[1.2426,1.2715]` |
| 11 | 8 | 1024 | 0.45489 us | 0.40927 us | **1.1115x** `[1.1043,1.1187]` |
| 12 | 7 | 256 | 0.13954 us | 0.11007 us | **1.2677x** `[1.2469,1.2888]` |
| 12 | 7 | 1024 | 0.45750 us | 0.40955 us | **1.1171x** `[1.1017,1.1327]` |
| 13 | 8 | 256 | 0.15484 us | 0.12076 us | **1.2822x** `[1.2657,1.2990]` |

## Exact legacy constructions

All three shapes use the existing parent `[N=32,D=24]` with `T=8`. They are
shortening or puncturing identities, not new codes.

For `K=11,R=8`, parent message coordinate 19 is shortened. The executor binds
the twelfth input of the exact `K=12` circuit to an immutable zero shard. Since
the encoder is linear, this is exactly the `K=11` legacy code and retains all
eight transmitted parity rows.

For `K=12,R=7`, parent parity coordinate 7 is punctured. The exact `R=8`
circuit runs unchanged and scatters that one unused row into stack-owned
storage. Coordinates 0 through 6 are therefore byte-identical to the legacy
encoder.

For `K=13,R=8`, parent message coordinate 20 becomes active. An independent
direct systematic-generator oracle derived its contribution at parity
coordinates 0 through 7 as the legacy GF8 logarithms:

    229, 94, 57, 147, 121, 151, 78, 228

The AVX2 specialization keeps that source live through the exact `K=12`
schedule and multiply-adds each coefficient immediately before its destination
row is scattered. The compile-time-false `K=12` specialization is unchanged at
1505 bytes. The `K=13` specialization is 1841 bytes and uses three YMM spill
slots rather than two.

The independent tests do not merely compare the new transform with itself.
They compare all selected shapes with the direct systematic-generator oracle,
and targeted GF(2) basis-bit cases verify every `K=13` source-coordinate-20
coefficient.

## Public and reusable paths

The packed public terminals were widened only for the qualified shapes. The
1024-byte reusable shape mask gained `K=12,R=7`; all selector decisions remain
codec-setup work rather than branches in the byte loop. Layouts that fail the
packed validation, sparse requested parity, forced diagnostics, and non-AVX2
backends retain the mature path.

| K | R | bytes | main | Leopard2 binding | main / binding, 95% CI |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 11 | 8 | 256 | 0.13835 us | 0.10337 us | **1.3384x** `[1.3216,1.3555]` |
| 11 | 8 | 1024 | 0.45489 us | 0.40285 us | **1.1292x** `[1.1236,1.1348]` |
| 12 | 7 | 256 | 0.13954 us | 0.10220 us | **1.3653x** `[1.3401,1.3910]` |
| 12 | 7 | 1024 | 0.45750 us | 0.40247 us | **1.1367x** `[1.1192,1.1545]` |
| 13 | 8 | 256 | 0.15484 us | 0.11281 us | **1.3726x** `[1.3543,1.3910]` |

The same frozen candidate was also run with only the terminal selector toggled.
This isolates the new route from compiler or binary-layout differences:

| K | R | bytes | control / selected, 95% CI |
| ---: | ---: | ---: | ---: |
| 11 | 8 | 256 | **1.5844x** `[1.5709,1.5980]` |
| 11 | 8 | 1024 | **1.1567x** `[1.1446,1.1689]` |
| 12 | 7 | 256 | **1.5635x** `[1.5443,1.5830]` |
| 12 | 7 | 1024 | **1.1325x** `[1.1187,1.1464]` |
| 13 | 8 | 256 | **1.5464x** `[1.5380,1.5547]` |

The public scratch contract is unchanged: 4672 bytes per stripe at 256 bytes
and 16960 at 1024 bytes, aligned to 64 bytes. These exact byte kernels do not
consume shard-data scratch, but the query still covers retained fallbacks.

## Frozen campaign

The clean candidate was copied read-only and hashed as
`d42b7729af8c047af661c89731f3763b832e9daa7c728979af25f4d42abc9302`.
The independently linked exact-main binary remained
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.
Both hashes were verified after timing.

Each cell used nine alternating ABBA/BAAB rounds, two process medians per side
per round, 31 retained samples, 16 warmups, reuse 8192, batch one, and one
codec thread. CPU 30 was isolated for timing and SMT sibling 14 was reserved.
All 117 rounds were accepted on the first attempt, all accumulated zero
non-idle sibling jiffies, and all 468 process outputs agreed in their 117
seed-scoped original/parity digest groups. Confidence intervals are Student-t
intervals over the nine independent mean log contrasts.

## Correctness and rejection

Release focused tests and the 42-shape production sweep passed. GCC 13 and
Clang 18 strict-warning builds passed. The focused and production tests passed
under Clang 18 ASan+UBSan+LSan at peak test RSS of 53,504 and 127,680 KiB.
GF8-only and GF16-only archive builds also passed. Work stayed serial and
memory-bounded because this host had previously OOM-killed parallel builds.

The `K=13,R=8` 1024-byte circuit was not promoted. Its prototype failed the
five-percent criterion, and the retained fallback measures 0.9857x exact main
for public one-shot encoding (`[0.9787,0.9927]`). Reusable execution is neutral
at 0.9995x (`[0.9950,1.0040]`). The negative cell remains visible rather than
being hidden by the successful 256-byte result.

The compact checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_exact_neighbors_checkpoint_20260810.json`.
