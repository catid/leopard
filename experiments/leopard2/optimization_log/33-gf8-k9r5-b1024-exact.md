# GF8 K9/R5 1 KiB exact encoder

## Result

**Promoted.** The GF8 AVX2 legacy-high encoder now routes exact
`K=9,R=5,B=1024` packed layouts through the existing exact `K9/T8` circuit.
The qualified ordinary one-item batch path is **1.4501x faster than exact
Leopard main** and the one-shot path is **1.4511x faster**. The same frozen
Leopard2 binary with only this terminal disabled is 1.4737x and 1.4842x slower
respectively.

| Metric | candidate | baseline | baseline / candidate, 95% CI |
| --- | ---: | ---: | ---: |
| one-item batch vs disabled terminal | 0.30017 us | 0.44236 us | **1.4737x** `[1.4645,1.4829]` |
| one-item batch vs exact main | 0.30017 us | 0.43527 us | **1.4501x** `[1.4422,1.4579]` |
| one-shot vs disabled terminal | 0.29997 us | 0.44520 us | **1.4842x** `[1.4743,1.4941]` |
| one-shot vs exact main | 0.29997 us | 0.43527 us | **1.4511x** `[1.4440,1.4582]` |

The retained `K=9,R=5,B=256` terminal remains healthy: one-item batch
execution is 1.5439x exact main (`[1.5335,1.5543]`) and one-shot is 1.5757x
(`[1.5647,1.5869]`).

## Exact legacy construction

This is not a new code or wire profile. `T=8`, `N=32`, and `D=24`; transmitted
parity coordinates are 0 through 4, actual message coordinates are 8 through
16, parity 5 through 7 is punctured, and message 17 through 31 is shortened to
zero. The ninth source is parent coordinate 16. Its independent direct
systematic-generator contribution to the five transmitted parity rows has the
legacy GF8 logarithms:

    121, 151, 78, 228, 229

The AVX2 route transforms the first eight sources with the existing exact T8
circuit, multiply-adds that ninth source column, and writes only the five
transmitted rows. Direct algebra and legacy byte-oracle tests cover every
output byte. Re-encoding after source mutation proves reusable bindings read
live application buffers rather than cached values.

## Routing and contract

The route is restricted to AVX2, GF8, legacy-high, exact `K=9,R=5`, exactly
1024 bytes, packed complete inputs, and five packed outputs. Ordinary encode,
one-item batch encode, and immutable reusable bindings select it. The binding
snapshots its route, so later changes to the diagnostic selector do not alter
an existing binding. This campaign times ordinary one-item batch and one-shot
encoding; reusable bindings are correctness-covered but are not a performance
claim in this report.

Sparse parity requests, detached or non-packed layouts, byte tails, other
counts and sizes, lower backends, and forced diagnostics retain the mature
path. Scratch remains caller-owned and allocation-free during execution. The
public query is 16,832 bytes per stripe at 64-byte alignment; validation covers
undersized and misaligned scratch, aggregate overlap and overflow, unaligned
buffers, changed sources, and multi-item batches.

## Neighbor proof

Candidate and control were the same executable; control differed only by
disabling this one terminal at runtime. Every 95% CI below is wholly contained
in the predeclared equivalence band `[1/1.02,1.02]`.

| inactive cell | control / candidate execution | control / candidate one-shot |
| --- | ---: | ---: |
| K9/R5/B1023 | 1.0032x `[0.9974,1.0091]` | 1.0023x `[0.9948,1.0099]` |
| K9/R5/B1025 | 1.0037x `[0.9980,1.0095]` | 1.0027x `[0.9967,1.0086]` |
| K9/R6/B1024 | 1.0035x `[0.9963,1.0108]` | 1.0014x `[0.9940,1.0088]` |
| K8/R5/B1024 | 1.0004x `[0.9965,1.0043]` | 1.0016x `[0.9978,1.0054]` |
| K10/R5/B1024 | 1.0000x `[0.9964,1.0037]` | 1.0008x `[0.9980,1.0036]` |

## Frozen campaign and validation

The clean pure-AVX2 candidate at commit `95bae0613a14ec9c0b15c220cfc77d0a30e566a4`
was copied read-only and hashed as
`d2604bfc93baa39340062831cfb1222f43a4017abd2e7b2fbb580d638f3bc0ce`.
The independently linked exact-main binary at commit `6e5725eb` was
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.
The runner proved full source/object/archive/link closure and reproduced the
candidate byte-for-byte with a serial clean rebuild before timing, then
rechecked inputs, frozen binaries, runner dependencies, and closure afterward.

Each of seven cells used 25 independent alternating paired rounds, 31 retained
samples per process, 64 warmups, reuse 8192, batch one, and one codec thread.
CPU 29 was pinned and SMT sibling 13 reserved. There were 175 accepted rounds
and 800 accepted process outputs; all workload/output digests matched and all
accepted rounds recorded zero sibling non-idle jiffies. One four-process round
attempt observed one sibling jiffy and was discarded and repeated before
analysis. Confidence intervals use Student-t over mean log contrasts with
24 degrees of freedom.

Before the confirmatory run, the full serial Release suite passed 143/143 and
the pure-AVX2 focused production suite passed 3/3. After the campaign, GCC 13
and Clang 18 strict library builds passed with
`-Wall -Wextra -Wpedantic -Werror`; the three focused production tests passed
under Clang 18 ASan+UBSan+LSan at a peak test RSS of 347,648 KiB. An initial
VA-limited launch was invalid because ASan could not reserve shadow memory and
failed before test code executed; the passing rerun removed the VA cap and used
a 1 GiB RSS cgroup with swap disabled. The benchmark wrapper passed under
normal Python and `python -O`. Independent production, math, API,
alias/overflow, compile-variant, evidence-schema, and final-diff reviews found
no remaining concrete defect.

Two earlier nine-round runs were sample-size pilots: each failed the strict
full-CI equivalence gate on a different near-1.0 neighbor. They were not pooled
or retried until favorable. The confirmatory 25-round design was predeclared,
versioned as evidence schema v2, and run once.

The compact checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_k9r5_b1024_exact_checkpoint_20260810.json`. Its raw and summary hashes bind
the uncommitted large artifacts retained under `/tmp`.
