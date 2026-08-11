# GF8 K13/R8 1 KiB AVX2 candidates

## Result

**Rejected; no production code retained.** The existing legacy-high encoder is
about 1.45% slower than exact Leopard main in public one-shot timing at
`K=13,R=8,B=1024`; reusable execution is statistically tied. Five materially
different AVX2 implementations were correct, but each improved the current
fallback by only 1.4-2.4%, below the predeclared 5% promotion threshold.

| Candidate | batch execution control / candidate | one-shot control / candidate |
| --- | ---: | ---: |
| direct generator column | 1.0212x `[1.0174,1.0250]` | 1.0199x `[1.0089,1.0310]` |
| pairwise tail reload | 1.0197x `[1.0146,1.0248]` | **1.0221x** `[1.0165,1.0277]` |
| separate tail pass | **1.0239x** `[1.0169,1.0310]` | 1.0175x `[1.0100,1.0250]` |
| integrated sparse IFFT | 1.0213x `[1.0139,1.0287]` | 1.0143x `[1.0092,1.0196]` |
| materialized sparse block | 1.0210x `[1.0147,1.0273]` | 1.0175x `[1.0096,1.0254]` |

Ratios above one favor the candidate. These are bounded diagnostic screens,
not promotion evidence: they held the canonical benchmark lock and pinned CPU
6 while reserving sibling 22, but the CPU pair was not externally sealed. No
candidate approached the promotion gate, so a larger authoritative neighbor
campaign would not change the kill decision.

## Algorithms tested

The legacy parent is `T=8,N=32,D=24`: parity coordinates 0-7, messages 8-20,
and shortened zeros 21-31. Source coordinate 20 contributes generator logs
`229,94,57,147,121,151,78,228` to the eight parity outputs.

The direct form reused the exact K12 circuit and multiply-added that source
column. The reload form shortened the extra source's live range. The separate
pass reduced register pressure at the cost of another output pass. The
integrated sparse form evaluated the five-input shortened second block while
the first block remained live, which caused spills. Finally, the materialized
sparse form performed the paper-faithful shortened T8 inverse transform with
only eight live vectors and no zero-shard loads; its 2,237-byte function had
one YMM spill slot and no EVEX instructions. Their nearly identical timings
show that the unavoidable fifth-source arithmetic, rather than unused-zero
work, is the practical floor on this AVX2 host.

## Correctness and method

Every timed process passed the Leopard2 round trip and candidate/control
digests matched across nine independent seed groups. A maximum-loss seed also
matched the control for original, parity, and recovered digests:
`cbba74dca053af8c`, `cd850a4e81074898`, and `ec7775d85639644e`.
The direct systematic generator and legacy parity oracle independently checked
the wire bytes.

Each candidate used nine alternating ABBA/BAAB rounds, four invocations per
round, 15 retained samples, eight warmups, reuse 8192, one codec thread, and a
256 MiB process limit. Intervals are Student-t 95% intervals over the nine
paired mean-log contrasts. Candidate binaries and raw-directory manifest
hashes are recorded in
`experiments/leopard2/gf8_high_encode/results/`
`t8_k13r8_b1024_avx2_negative_20260811.json`.

The mature production path remains selected. This is preferable to adding a
large exact kernel for a gain that is too small to satisfy the project's
maintenance and neighboring-regression policy.

After reverting every prototype, a serial memory-capped rebuild matched all
specialized executable sections byte-for-byte. The whole benchmark hash differs
because it embeds the newer Git commit/tree identity; only one `.text` byte
changed, the displacement of that reporting string inside the benchmark's
`Run()` function. No codec or terminal instruction differs from the frozen
control.
