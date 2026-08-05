# Whole-shard direct repair for tiny native-high decode

## Result

**Promoted.** Commit `edeb558e95edda97268fb93d4a4cbb3634bbd903`
replaces the scratch-resident Algorithm 5 executor in the bounded native-high
one-shot region with locator-derived direct repair over the complete shard.
Against exact Leopard `main`, all 25 tested cells win by 1.2886x--4.5130x.
The former `K=16,R=8,L=8,B=65` gap moves from 0.829x to **2.2282x**.

## Why the previous route lost

The raw setup introduced in report 27 removed heap-plan construction, but its
byte execution still ran a padded Algorithm 5 transform tile. A 65-byte shard
therefore performed one full 64-byte tile plus another complete transform for
the one-byte suffix, then gathered that suffix. The paper's asymptotic
`O(N log T)` advantage cannot repay fixed transform scheduling and staging at
these tiny dimensions.

The original hypothesis was to repair only the ragged suffix directly. A
standalone screen showed that the locator-derived direct row was cheaper over
the entire 1--256-byte shard, so the final candidate removes every transform
tile in this bounded region instead of adding a prefix/tail branch.

## Arithmetic and implementation

Raw setup already computes the erasure locator and deterministically chooses
exactly `K` received public coordinates. For a missing original coordinate
`x` and a selected survivor `s`, the direct coefficient is

    Lambda(s) / ((x xor s) Lambda'(x)).

In Leopard's logarithmic representation this is

    locator[s] - locator[x] - log(s xor x).

The derivation is recorded in `docs/leopard2_math_and_sources.md`. Setup emits
an output-major `K` by `L` log-coefficient matrix on the stack, then
`AVX2FF8LinearCombinationTiny` computes all `L` outputs over the requested
positive byte count. The kernel supports one through eight outputs and handles
the tail without reading or writing outside the caller's buffers.

Production selection remains deliberately narrow:

- legacy-high GF8 with the effective pure-AVX2 backend;
- native `N=32,T=8`, `K=9..16`, `R=5..8`, and `L=5..R`;
- public one-shot decode and shard sizes 1..256 bytes.

There is no hot-path allocation. Reusable plans, other fields/backends,
profiles, forced diagnostics, and the 257-byte boundary keep their established
routes.

## Correctness and safety

The focused Release matrix covers every `K=9..16`, `R=5..8`, `L=5..R`, and
every byte count 1..256:

- 9,088 allocation-audited public one-shot cases;
- 10,752 comparisons with both the independent systematic generator oracle
  and reusable Algorithm 5;
- zero eligible hot-path allocations;
- GF8-only and GF16-only library builds; and
- Clang 18 ASan+UBSan runs of both focused executables.

The Release tests peaked at 22 and 30 MiB RSS; the sanitizer tests peaked at
78 and 106 MiB. A RelWithDebInfo sanitizer compile exceeded an intentionally
tight 512-MiB virtual cap inside LLVM debug-value analysis. A one-job Release
sanitizer build succeeded under 768 MiB at 296,824 KiB peak RSS, distinguishing
the compiler cap from a codec defect. ASan execution was not virtual-address
capped because its shadow mapping requires a large address range.

## Pinned comparison with exact Leopard main

Candidate and exact-main executables were copied to a read-only lane and
hashed before and after timing. Exact main is commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`; both builds have a pure-AVX2 ISA
ceiling. Three ABBA rounds per cell ran on CPU 28 with SMT sibling 12 reserved,
one codec thread, 31 samples, four warmups, and reuse 400. One non-idle sibling
jiffy occurred over the complete 300-process campaign, so this is controlled
pinned evidence without an external affinity-supervisor seal.

Representative exact-main-over-Leopard2 one-shot ratios are:

| Cell | Leopard main us | Leopard2 us | Speedup | 95% CI |
| --- | ---: | ---: | ---: | ---: |
| K=9 R=5 L=5, 1 B | 0.8160 | 0.1808 | 4.5130x | 4.4500--4.5770 |
| K=9 R=5 L=5, 65 B | 0.9209 | 0.2057 | 4.4779x | 4.4743--4.4815 |
| K=12 R=7 L=7, 256 B | 1.1478 | 0.4174 | 2.7498x | 2.7114--2.7887 |
| K=16 R=8 L=8, 63 B | 0.8336 | 0.6469 | 1.2886x | 1.2830--1.2943 |
| K=16 R=8 L=8, 65 B | 0.9474 | 0.4252 | 2.2282x | 2.2029--2.2538 |
| K=16 R=8 L=8, 256 B | 1.1756 | 0.5999 | 1.9599x | 1.9395--1.9804 |

Leopard main rounds nonmultiple byte counts to physical 64-byte storage while
both programs fingerprint the same logical prefix. This gives application
equivalence, not identical physical byte work, and is recorded explicitly in
the result artifact.

## Same-source attribution and boundary

The whole-direct candidate was also compared with the pre-change raw
Algorithm 5 control. All seven affected boundary sizes improved by
1.1987x--2.2812x. At the deliberately ineligible 257-byte boundary, the ratio
was 0.9846x with a 95% interval of 0.9708--0.9986. That is a statistically
credible 1.54-percent binary-layout/process difference, below the predeclared
two-percent regression gate; the new code path is not selected there. It would
be inaccurate to call this control neutral.

The direct executor is likely useful beyond 256 bytes, but that requires a new
direct-versus-Algorithm-5 crossover sweep rather than an unmeasured widening.

Complete identities, hashes, all 25 exact-main cells, all eight same-source
cells, test counts, confidence intervals, and limitations are in
`experiments/leopard2/direct_repair/results/raw_native_high_whole_direct_avx2_20260805.json`.
