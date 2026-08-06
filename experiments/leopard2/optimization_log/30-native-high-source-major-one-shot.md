# Native-high source-major one-shot execution

## Result

**Promoted for 12,288 through 16,384-byte shards.** Commit
`9498e93485b46c747d6440b056f893facadea66f` removes the heap-owned transient
decode plan in this interval for legacy-high GF8, pure-AVX2 `N=32,T=8`,
`K=9..16`, `R=5..8`, and `L=5..R` missing originals.

The route keeps the caller-scratch locator setup from report 29, derives the
direct coefficients and exact deterministic source order into bounded stack
metadata, and calls the same pointer-level source-major executor used by
reusable plans. It therefore avoids both plan allocation and a second
maintenance copy of the byte-heavy loop. Shards through 7168 bytes retain the
output-major raw route; 7169 through 12,287 bytes and 16,385 bytes upward retain
the mature plan fallback.

The source-major executor handles a remainder of zero through 15 bytes exactly.
For a remainder of 16 through 31 bytes it runs one padded, aligned 32-byte AVX2
tile and copies only the requested prefix to the destination. This tail is
bounded stack storage and does not change the public scratch requirement,
aliasing contract, arithmetic, or wire profile.

## Correctness gate

The focused Release suite covers 20,992 allocation-audited public one-shot
cases and 12,432 independent-generator/reusable-Algorithm-5 differentials.
Coverage includes every valid `(K,R,L)` shape, both route transitions, every
32- and 64-byte residue, logarithmic sizes through 1 MiB, mixed parity
availability, unaligned buffers, alias rejection, exact scratch, and failure
atomicity. Every selected public call performed zero allocation.

The complete dual-field Release suite passed 137/137. The three focused tests
passed under Clang 18 ASan+UBSan. GF8-only and GF16-only archives built, and
the GF16-only suite passed 9/9. The production tail can also be compiled out
independently; that build succeeded. Full Release validation was serialized or
restricted to two tests at a time after earlier high-parallelism runs exhausted
host memory.

## Frozen same-source gate

Candidate and control came from the same committed source and compile the same
route. They differ only in the initialized volatile diagnostic word that
selects the ordinary plan fallback. Their `.text` sections have the same
SHA-256. Frozen binaries were copied to a lane-owned directory and verified
unchanged before and after timing.

All 80 valid `(K,R,L)` shapes were measured at both aligned interval endpoints,
for 160 cells total. Three ABBA rounds used CPU 28, reserved SMT sibling 12, 31
retained samples, four warmups, reuse 200, one codec thread, and the canonical
campaign lock. Speedup versus the plan fallback was:

| Statistic | Result |
| --- | ---: |
| Minimum cell geomean | 1.0906x |
| Maximum cell geomean | 1.2144x |
| Overall geomean | 1.1283x |
| Weakest 95% lower bound | 1.0699x |

The reserved sibling accumulated 4430 idle jiffies and no non-idle jiffies.
Every workload digest matched.

## Ragged-tail and neighbor gates

A 320-cell screen crossed five representative shapes, both interval regions,
and every residue. Its overall geomean was 1.1127x, and every residue aggregate
was at least 1.0637x. Twenty-two noisy single-round cells fell below 1.05x, so
the two weakest tail shapes were rerun with 15 independent ABBA seeds:

| Cell | Speedup | 95% CI |
| --- | ---: | ---: |
| K=12 R=8 L=5, 16,382 B | 1.1065x | 1.0863--1.1270x |
| K=14 R=5 L=5, 12,298 B | 1.1127x | 1.0906--1.1353x |

Six inert neighbors at 12,287 and 16,385 bytes were also rerun with 15 seeds.
Their point ratios ranged from 0.9899x to 1.0099x, every confidence interval
spanned one, and no point regression exceeded 1.01 percent. The reserved
sibling was fully idle in that campaign.

## Exact Leopard-main gate

Fourteen representative selected cells were replayed against a separately
built exact Leopard main at commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. Leopard main over Leopard2 ranged
from 2.1216x to 5.3870x, with an overall geomean of 3.4837x. The weakest cell's
95-percent interval was 2.1000--2.1434x. All recovered-data digests matched.

Complete identities, hashes, aggregate statistics, robust tail confirmations,
neighbor results, and test counts are in
`experiments/leopard2/direct_repair/results/raw_native_high_source_major_20260805.json`.
