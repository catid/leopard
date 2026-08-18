# Leopard2 versus exact Leopard main

## R10 RS(256,K) decoder reproduction: 2026-08-11

The clean pure-AVX2 build at commit `a0d781c` reproduces the central performance
result of Chen et al. R10: the specialized low-rate Algorithm 4 and high-rate
Algorithm 5 decoder is faster at every paper K value. This screen uses
1024-byte shards, corresponding to the paper's 1024 codewords sharing one
erasure pattern, and charges plan construction on every public one-shot decode.
Each row is the geometric mean over five independently shuffled patterns with
exactly `R=256-K` erasures. Ratios above one favor Leopard2.

| K,R | Profile / algorithm | Retained generic / Leopard2 | Exact Leopard main / Leopard2 |
| --- | --- | ---: | ---: |
| 8,248 | low / Algorithm 4 | **8.000x** | unavailable (`R>K`) |
| 16,240 | low / Algorithm 4 | **3.254x** | unavailable (`R>K`) |
| 32,224 | low / Algorithm 4 | **1.876x** | unavailable (`R>K`) |
| 64,192 | low / Algorithm 4 | **1.597x** | unavailable (`R>K`) |
| 128,128 | low / Algorithm 4 | **1.337x** | **1.130x** (different wire profile) |
| 192,64 | legacy high / Algorithm 5 | **1.359x** | **1.149x** |
| 224,32 | legacy high / Algorithm 5 | **1.584x** | **1.398x** |
| 240,16 | legacy high / Algorithm 5 | **1.604x** | **1.469x** |
| 248,8 | legacy high / Algorithm 5 | **2.791x** | **2.489x** |

The initial K=64/128 and K=224 deficits were not extra Algorithm 4/5 byte
work. Reusable specialized execution already won. The public one-shot path was
compiling exact pruned schedules whose construction cost could not be amortized
at one use. Exact-byte one-shot plans now execute the mature regular transforms
through 1 KiB for native low rate and for the measured GF8 T=32 region;
caller-created reusable plans still retain and use exact schedules.

The T=32 promotion screen covered 81 `(K,R,loss,bytes)` cells. At 1 KiB its
one-shot geometric-mean improvement over the preceding Leopard2 build is
1.462x, with a 1.229x weakest cell; the 1025-byte controls are neutral at
1.000x geometric mean and remain on the exact-schedule path. At K=224/R=32,
Leopard2 one-shot decode is 1.799x, 1.436x, 1.385x, and 1.387x faster than the
in-tree legacy Leopard1 path for 2, 16, 31, and 32 missing originals
respectively. All workload,
parity, and recovered-output comparisons match.

The source tree was clean, candidate and exact-main executables were frozen and
rehashed after timing, and disassembly found no ZMM, opmask, `vpternlog`, or
GFNI instructions. Release correctness tests passed 6/6 and focused Clang
ASan+UBSan tests passed 2/2. This is a pinned promotion screen, not a claim that
this host reproduces the paper's absolute throughput or a substitute for a
cross-machine ABBA campaign. Full numbers and raw JSON are in
`experiments/leopard2/r10_results/a0d781c-avx2-summary.json` and its attested
raw bundle.

## GF8 K=12/R=8 exact high-encode update: 2026-08-10

The exact legacy-high GF8/AVX2 circuit at `K=12,R=8` now skips the fully
shortened third message block, inverse-transforms only four active rows in the
second block, folds the fixed shift, and removes the zero multipliers from the
T=8 forward schedule. A packed public terminal also avoids the general range
sort and transform-geometry setup at the two measured byte counts.

| bytes/shard | exact main | Leopard2 | exact main / Leopard2, 95% CI |
| ---: | ---: | ---: | ---: |
| 256 | 0.13720 us | 0.11033 us | **1.2436x** `[1.2207,1.2668]` |
| 1024 | 0.45386 us | 0.40666 us | **1.1161x** `[1.1114,1.1208]` |

The clean candidate is commit `64a5dd9`; exact Leopard main is commit
`6e5725e`. Nine balanced ABBA/BAAB rounds per cell used 31 retained samples,
reuse 8192, one thread, CPU 30, and reserved sibling 14. All accepted rounds
recorded zero sibling work and all parity digests matched. A same-executable
selector comparison measured 1.8763x at 256 bytes and 1.2061x at 1024 bytes
over the ordinary-validation control. Eight K/R/byte selector neighbors had
confidence intervals spanning parity.

This is deliberately not generalized to every nearby T=8 shape. Exact-main
measurements still show 1024-byte gaps at K11/R8, K12/R7, and K13/R8, plus
larger 256-byte gaps in the first two. Those are tracked separately. Complete
binary identities, neighbor results, validation, and the raw-bundle hash are
in `experiments/leopard2/gf8_high_encode/results/`
`t8_k12r8_exact_checkpoint_20260810.json`.

## GF8 R=1 small-shard update: 2026-08-02

The latest isolated campaign compares a frozen pure-AVX2 Leopard2 candidate
with exact Leopard main commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.
It completed 111 cells and 2,664 fresh timed processes on CPU 14 while sibling
30 accumulated zero non-idle jiffies.  Original, parity, and recovered-output
digests matched in every process.  Candidate and control were equal-length
hard links to one immutable executable, and both executables contained zero
EVEX instructions.

Four exact selector cells cleared a strict five-percent lower-confidence gate
in batch encode, reused decode, one-shot encode, and one-shot decode:

| K, bytes | Internal batch encode | Internal reused decode | Exact main / Leopard2 one-shot encode | Exact main / Leopard2 one-shot decode |
| --- | ---: | ---: | ---: | ---: |
| 2, 64 | 1.856x | 2.126x | 0.494x | 0.571x |
| 2, 256 | 1.784x | 1.910x | 0.623x | 0.661x |
| 2, 1024 | 1.616x | 1.745x | 0.954x | 0.857x |
| 4, 1024 | 1.262x | 1.156x | 0.922x | 0.892x |

The internal columns compare the selected kernel with Leopard2's preceding
pairwise route.  The exact-main columns include the real public one-shot APIs.
Thus the new kernels materially reduce work but do not yet erase Leopard2's
validation, dispatch, and setup overhead at these tiny cells.  Across the full
111-cell sample, Leopard2 was credibly slower than exact main in 76 of 111
one-shot-encode cells and 92 of 111 one-shot-decode cells.  This is scoped to
GF8, legacy-high R=1, one loss, one thread, batch one, and cache-hot repeated
execution; it is not a statement about the high/low transform decoders.

K=3 fused-final at 64/256 bytes, K=4 dense at 64/256 bytes, and broad 64/256
byte coarse thresholds were rejected.  Their crossovers were losing,
inconclusive, or nonmonotonic.  The machine-readable checkpoint is
`experiments/leopard2/r1_xor/results/`
`8ff4ed9-small-reduction-abba-20260802.json`.

## GF8 T=8 ragged 65--928-byte update: 2026-07-31

Dense legacy-high GF8/AVX2 T=8 profiles now reuse the prepared binding over
measured ragged byte tiers from 65 through 928 bytes. The deterministic
selector covers 1,213 qualified `(K,R,bytes)` target cells and leaves five
measured exclusions plus 840 other route neighbors on the prior path.

| population | comparison | geometric-mean speedup | minimum 95% lower bound |
| --- | --- | ---: | ---: |
| 1,213-cell discovery/holdout selector | selector-off control / Leopard2 | 2.0296x | 1.0504x |
| 1,213-cell discovery/holdout selector | padded exact Leopard main / Leopard2 | 1.2851x | 1.000026x |
| 14-cell final-source confirmation | selector-off control / Leopard2 | 1.7715x | 1.0834x |
| 14-cell final-source confirmation | padded exact Leopard main / Leopard2 | 1.0598x | 1.00091x |

The discovery and predeclared holdout contain 2,058 cells and 32,004 accepted
processes. The separate final-source campaign contains 14 targets, five
neighbors, and 816 accepted processes. Every logical input, parity, and
recovery digest matched; all accepted rounds passed the isolation gate; and
the frozen binary hashes were unchanged after timing.

Exact Leopard main requires a physical shard length divisible by 64, so the
main comparison processes a zero-padded
`ceil(logical_shard_bytes/64)*64`-byte shard and verifies the common logical
prefix. The same-source comparison processes the exact logical byte count in
both paths. Five cells remain excluded:
`(K,R,bytes)=(5,5,191),(6,5,319),(6,6,319),(7,5,319),(7,6,319)`.
The compact checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_ragged_checkpoint_20260731.json`.

## GF8 T=8 arbitrary-tail update: 2026-07-31

Every dense legacy-high GF8/AVX2 T=8 profile with `K=5..16` and
`R=5..min(K,8)` now uses the fused prepared binding at byte counts
`1,2,3,7,8,15,16,17,31,32,33,63`. The combined frozen campaign contains 504
target cells, 84 64/65-byte route neighbors, and a predeclared nine-round
holdout for the only two initially inconclusive cells.

| comparison | geometric-mean speedup | minimum point | minimum 95% lower bound |
| --- | ---: | ---: | ---: |
| Same-source selector-off control / Leopard2 | 2.9609x | 1.8349x | 1.6586x |
| Padded exact Leopard main / Leopard2 | 1.5669x | 1.1947x | 1.0263x |

All original, parity, and recovery digests matched, every accepted round had
an idle reserved SMT sibling, and no neighbor failed the two-percent
regression gate. Candidate and control executable instruction sections were
identical. The final production-default binary has that same executable-section
digest.

Leopard main requires a physical shard length divisible by 64. For the
sub-64-byte row it therefore processes a zero-padded 64-byte shard; input,
parity, and recovery digests cover only the matching logical prefix. This is
an application-equivalent comparison, whereas the same-source ratio times the
exact requested byte count in both paths. The compact checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_tiny_checkpoint_20260731.json`.

## GF8 T=8 one-kibibyte update: 2026-07-31

Four additional dense legacy-high GF8/AVX2 encode profiles now use the
wire-identical direct-input T=8 binding at exactly 1024 bytes. The final
current-source comparison used batch/reuse 64, nine target rounds, an exact
Leopard main binary from commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, matching workload/parity/recovery
digests, and a reserved idle SMT sibling.

| K,R | exact main / Leopard2, 95% CI |
| --- | ---: |
| 6,5 | 1.1396x `[1.1300,1.1492]` |
| 6,6 | 1.1404x `[1.1302,1.1508]` |
| 10,8 | 1.1582x `[1.1503,1.1661]` |
| 16,5 | 1.3028x `[1.2934,1.3123]` |

All four same-source improvements also have lower confidence bounds above
1.05, and no one of the 38 unmodified K/R neighbors had a credible regression
over two percent. `K=11,R=5` and `K=11,R=6` remain on the prior path because
their nine-round lower same-source bounds were 1.0475x and 1.0477x. The
machine-readable checkpoint is
`experiments/leopard2/gf8_high_encode/results/`
`t8_one_kib_checkpoint_20260731.json`.

## Optimization update: 2026-07-21

The production stack after the paper-faithfulness audit has materially moved
the remaining balanced and small-loss cells.  The compact checkpoints below
compare against exact Leopard default-branch commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  Ratios remain exact-main time
divided by Leopard2 time.  Each result used counterbalanced pinned processes,
matching original/parity/recovery digests, and separately reported Leopard2
plan setup.  These are scoped optimization checkpoints rather than a rerun of
the entire uniform matrix in the next section.

| Decode cell | Leopard2 path | Plan reuse | Amortized speedup |
| --- | --- | ---: | ---: |
| GF8 K=240 R=16, 4 KiB, 16 losses | Algorithm 5 with zero-shift cancellation | 64 | 1.974x |
| GF8 K=127 R=65, 4 KiB, 65 losses | translated Algorithm 4 | 64 | 1.195x |
| GF8 K=100 R=100, 4 KiB, 50 losses | translated Algorithm 4 | 64 | 1.299x |
| GF8 K=128 R=128, 64 KiB, 128 losses | translated Algorithm 4 | 8 | 1.261x |
| GF16 K=255 R=129, 4 KiB, 129 losses | translated Algorithm 4 | 32 | 1.653x |
| GF8 K=65 R=65, 64 B, 8 losses | direct eight-loss repair | 512 | 1.572x |
| GF8 K=65 R=65, 1 KiB, 8 losses | direct eight-loss repair | 256 | 3.068x |
| GF8 K=65 R=65, 64 KiB, 8 losses | direct eight-loss repair | 16 | 3.870x |
| GF8 K=65 R=65, 1 MiB, 8 losses | direct eight-loss repair | 4 | 7.561x |

The 64-byte direct cell is reuse-sensitive: charging its 11.215-us plan setup
once gives 13.563 us versus Leopard1's 3.725 us, and the aggregate crossover
is reuse nine.  At 1 KiB and above it already wins with setup charged once.
The direct selector is deliberately limited to AVX2 GF8 legacy-high K=65,
R=65..128, and at most eight missing originals; neighboring shapes retain the
paper transform paths.

Balanced GF8 encoding also now uses a runtime-qualified, wire-identical
AVX-512VL whole-transform callback on the measured AMD processor family.
Exact-main speedups are 1.091x/1.191x for T=16 at 4/64 KiB;
1.095x/1.085x/1.158x for T=32 at 2/4/64 KiB; and
1.102x/1.081x/1.117x for T=64 at 2/4/64 KiB.  A noisy T=16, 2-KiB cell was
rejected and remains on AVX2.  A later GF16 terminal-range callback experiment
was also rejected at only 1.016x and 1.000x in its two targets, so none of that
candidate code is present.

The dedicated R=1 coarse XOR path closes the earlier GF8 XOR regression.  Five
separate exact-main campaigns around the K=128 power boundary all passed the
zero-sibling-work and byte-identity gates.  At K=129, Leopard2 encode is
1.104x at 4 KiB, 1.329x at 64 KiB, and 1.098x at 1 MiB; first-use decode is
1.013x, 1.183x, and 1.103x respectively.  The K=128 and K=130 4-KiB neighbors
show no credible regression over two percent, and K=130 encode is credibly
1.081x.  The apparent K=130 decode loss in an earlier contaminated campaign
did not reproduce under the accepted isolation gate, so no narrower selector
is needed.

The later GF16 K=1000/R=200/64-KiB byte-tiling checkpoint halves Leopard2
scratch from 33,585,728 to 16,808,512 bytes and improves the prior Leopard2
implementation by 1.135x for full parity (1.108x--1.147x for measured parity
subsets).  A clean exact-main ABBA gate nevertheless measures only 0.97678x,
95% CI [0.97349, 0.98008], so Leopard2 remains 2.32% slower than exact Leopard
main in that full-output cell.  The optimization is retained as a composable
scratch and same-source improvement, but the GF16 exact-main encoder gap is not
closed.  A distinct K=1000/R=199 neighbor measures 1.01049x, 95% CI
[1.00875, 1.01224]; it is not evidence for the R=200 prefix-output branch,
because Leopard1 has no parity-subset API.

Machine-readable evidence, including binary/source identities, raw-bundle
hashes, isolation limits, selector bounds, and validation, is retained in:

- `experiments/leopard2/high_decode_copy/results/algorithm5_zero_shift_cancel_checkpoint.json`
- `experiments/leopard2/high_decode_copy/results/high_low_translation_checkpoint.json`
- `experiments/leopard2/direct_repair/results/k65_exact_main_checkpoint.json`
- `experiments/leopard2/gf8_balanced_encode/results/checkpoint.json`
- `experiments/leopard2/gf16_high_encode/results/byte_tiling_positive_20260721.json`
- `experiments/leopard2/gf16_high_encode/results/terminal_range_negative_20260721.json`
- `experiments/leopard2/r1_xor/results/e621109-exact-main-amd9950x3d/summary.json`

## Earlier full-matrix checkpoint: 2026-07-18

The current authoritative comparison measures detached, test-hooks-off
Leopard2 commit `a1fa5bf02641f0af8c344b21c879b9da0bbd133a` against the
independently linked exact Leopard default-branch commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  Ratios are exact-main time
divided by Leopard2 time, so values above one favor Leopard2.  Each cell is an
independent three-round ABBA bundle.  All eight bundles passed source/build
closure, original/parity/recovery digest, retained-stream, affinity, and
statistics verification, with zero non-idle jiffies on the reserved SMT
sibling.

| Cell | Encode, 95% CI | Decode first use, 95% CI | Decode reuse 8, 95% CI |
| --- | ---: | ---: | ---: |
| GF8 XOR, K=129 R=1 | 0.975 [0.958,0.992] | 0.977 [0.965,0.989] | 0.978 [0.966,0.991] |
| GF8 high, K=240 R=16, one loss | 1.172 [1.164,1.181] | 2.212 [2.196,2.228] | 2.218 [2.202,2.234] |
| GF8 high, K=240 R=16, full loss | 1.194 [1.172,1.217] | 1.666 [1.655,1.676] | 1.679 [1.669,1.689] |
| GF8 balanced, K=128 R=128 | 1.037 [1.035,1.039] | 1.049 [1.049,1.049] | 1.052 [1.052,1.052] |
| GF16 inflation, K=200 R=50 | 1.023 [1.020,1.027] | 1.932 [1.932,1.932] | 1.946 [1.945,1.946] |
| GF16 high, K=1000 R=200, one loss | 0.894 [0.888,0.900] | 3.383 [3.373,3.393] | 3.396 [3.387,3.406] |
| GF16 high, K=1000 R=200, full loss | 0.895 [0.890,0.900] | 2.283 [2.271,2.295] | 2.295 [2.283,2.307] |
| GF16 large, K=4096 R=512 | 0.949 [0.939,0.959] | 2.022 [1.962,2.085] | 2.212 [2.140,2.288] |

The planned high-rate decoder now provides credible speedups at every
nontrivial transform cell: about 1.05x at balanced GF8, 1.67-2.22x at GF8
high rate, and 1.93-3.40x at the measured GF16 first-use cells.  This includes
plan setup; reuse is slightly better.  GF8 high-rate encoding is now 1.17-1.19x
faster than main, and balanced/field-inflation encoding is modestly faster.
The remaining production priorities are explicit rather than hidden: the GF8
XOR cell is about 2.2-2.5% slower, GF16 K=1000/R=200 encode is about 10.5%
slower, and the larger GF16 encode cell is about 5.1% slower.  Those residuals
are tracked under the exact-main encoder gap work.

Verified current manifests are retained in the ignored research cache at
`.research/leopard2/exact-main-a1fa5bf-cpu4-cells-20260718/`.  Their individual
SHA-256 identities are recorded in the corresponding Beads evidence and can be
recomputed with the verifier while the detached source/build closure exists.

## Earlier checkpoint: 2026-07-16

The historical comparison below measures Leopard2 commit
`620667c93fad4eb9d97c8c5ccf7dcf6994914346` against the same exact detached
baseline.  It is retained to show the measured progression; its slower encode
and decode results must not be quoted as the current implementation.

## Method

The standalone adapter builds only the four codec translation units from the
detached baseline. GCC 13.3 compiled baseline code with its canonical effective
Release flags, including whole-build `-march=native`. Leopard2 used its shipped
portable-core policy with SSSE3 and AVX2 isolated in dedicated translation
units. Both used `-Wall -Wextra -fopenmp -g -O0 -O3`; the Leopard2 build had
tests, fuzzers, and CUDA disabled so test hooks could not affect dispatch.

The campaign used eight legacy-high wire-compatible cells, three independent
ABBA rounds per cell, order `main,Leopard2,Leopard2,main`, nine retained timing
samples per child, eight calls per sample, two warmups, one stripe, and one
thread. This produced 96 fresh child processes and 864 retained samples per
operation. Within each cell, every invocation used the same deterministic data
and loss set. Original, transmitted-parity, and recovered-original FNV-1a-64
digests matched between implementations in every cell.

Timing ran on logical CPU 15 of an AMD Ryzen 9 9950X3D. Its SMT sibling 31 was
reserved. The host has 16 physical cores, 32 logical CPUs, one socket, and one
NUMA node. `amd-pstate-epp` reported governor `powersave` with EPP
`performance`; starting package temperature was 38.1 C. During the 137.86
second window, CPU 31 accumulated 13,786 idle jiffies and zero non-idle
jiffies. Twenty-six user-owned threads were moved off the pair and restored
without error; one additional user-owned `systemd` thread could not be moved
because the kernel returned `EPERM`. The unprivileged coordinator also could
not repin root-owned or per-CPU kernel threads; clustered ABBA ordering limits
that residual noise.
Hardware counters were unavailable because `perf_event_paranoid=4`.

Speedup is `main time / Leopard2 time`: values above 1 favor Leopard2.
Confidence intervals use one mean log contrast per independent ABBA round,
Student-t with 2 degrees of freedom. Absolute times below are the median across
six process-level retained-sample medians; small differences from the paired
geometric ratio are expected.

## Encoding

All cells use the exact legacy-high coordinate profile and matched parity
bytes. Leopard2 encode throughput was 7.2% to 23.4% lower than main, with every
95% interval below parity.

| Cell | K/R | Bytes/shard | main GB/s | Leopard2 GB/s | Speedup, 95% CI |
|---|---:|---:|---:|---:|---:|
| XOR GF8 | 129/1 | 64 KiB | 101.74 | 82.88 | 0.806 [0.784, 0.829] |
| GF8 high, one loss | 240/16 | 64 KiB | 22.68 | 17.06 | 0.766 [0.695, 0.845] |
| GF8 high, full loss | 240/16 | 64 KiB | 22.56 | 18.12 | 0.785 [0.737, 0.836] |
| GF8 balanced, full loss | 128/128 | 64 KiB | 9.78 | 9.09 | 0.928 [0.921, 0.936] |
| GF16 field inflation, 8 loss | 200/50 | 64 KiB | 9.46 | 8.54 | 0.903 [0.898, 0.907] |
| GF16 high, one loss | 1000/200 | 64 KiB | 7.02 | 6.20 | 0.884 [0.865, 0.904] |
| GF16 high, full loss | 1000/200 | 64 KiB | 6.98 | 6.19 | 0.887 [0.885, 0.889] |
| GF16 large parent, 8 loss | 4096/512 | 4 KiB | 7.76 | 6.41 | 0.825 [0.824, 0.827] |

The comparison does not isolate a single cause. The measured policies include
main's whole-build `-march=native`, Leopard2's portable core plus runtime AVX2
fixed kernels, and Leopard2 API/schedule overhead. Profiling and the open
pruned/tiled encoder work are required before attributing the regression.
Candidate `leopard-79h.26.3` makes ragged source staging (and Low V1 parity-tail
staging) fixed at one 64-byte tile per shard, but the table above used aligned
shards and predates that candidate; it is not evidence of a throughput change.

## Decoding

Main has no reusable erasure plan, so each decode call includes locator and
pattern setup. Leopard2 first use below is a derived median plan creation plus
one execution with the codec already created. Reuse-8 is execution plus one
eighth of plan creation. Plan and execution were timed in separate loops; these
are not presented as jointly timed calls. Codec creation is reported in the raw
evidence and excluded from both implementations' operation comparison.

| Cell, losses | main setup-included | L2 first use | L2 reuse-8 | First-use speedup, 95% CI | Reuse-8 speedup, 95% CI |
|---|---:|---:|---:|---:|---:|
| XOR GF8, 1 | 84.1 us | 121.4 us | 121.3 us | 0.693 [0.690, 0.696] | 0.694 [0.691, 0.697] |
| GF8 240/16, 1 | 1.976 ms | 1.752 ms | 1.751 ms | 1.123 [1.116, 1.130] | 1.124 [1.117, 1.131] |
| GF8 240/16, 16 | 2.286 ms | 2.318 ms | 2.317 ms | 0.986 [0.982, 0.990] | 0.986 [0.983, 0.990] |
| GF8 128/128, 128 | 2.185 ms | 3.526 ms | 3.524 ms | 0.617 [0.581, 0.654] | 0.617 [0.581, 0.654] |
| GF16 200/50, 8 | 6.637 ms | 4.037 ms | 4.028 ms | 1.648 [1.621, 1.675] | 1.652 [1.625, 1.679] |
| GF16 1000/200, 1 | 51.141 ms | 25.574 ms | 25.402 ms | 2.005 [1.989, 2.022] | 2.019 [2.002, 2.035] |
| GF16 1000/200, 200 | 61.417 ms | 36.695 ms | 36.521 ms | 1.675 [1.664, 1.685] | 1.683 [1.672, 1.694] |
| GF16 4096/512, 8 | 10.307 ms | 5.891 ms | 5.705 ms | 1.749 [1.698, 1.801] | 1.805 [1.752, 1.861] |

The strongest result is one-loss GF16 1000/200: Leopard2 cuts the derived
first-use time almost in half. All four GF16 cells show credible 1.65x to 2.01x
first-use gains. GF8 is regime-dependent: one-loss high rate improves 12.3%,
full-loss high rate regresses about 1.4%, XOR regresses about 44% in time, and
balanced full-loss Leopard2 takes about 1.61x as long as main.

## Current equal-rounded GF8 multi-loss result

The table above predates the source-major AVX2 direct-repair promotion.  On
final candidate source `00191af`, a frozen 47-target/five-neighbor campaign
plus one nine-round holdout compared reusable-plan execution with exact
Leopard1 `6e5725e`.  The combined exact-main-over-Leopard2 execution geomean is
5.0310x, with a 1.1041x weakest lower confidence bound.  The 64-KiB subset has
a 7.5633x geomean and 2.8748x weakest lower bound; the 1-KiB subset has a
6.4154x geomean and 2.8550x weakest lower bound.  All logical original,
parity, and recovered-output digests match.

This result applies to the measured equal-rounded legacy-high GF8/AVX2 region:
`17 <= K <= 128`, `ceil_pow2(K) == ceil_pow2(R)`, `K != 65`, and two through
eight missing originals.  It does not claim that a new plan is always faster.
At one through 65 bytes, execution wins all 14 cells but first use has a
0.6228x geomean versus Leopard1 because plan construction dominates; the
observed reuse crossover is two through fifteen calls.  That remaining gap is
`leopard-79h.38.5.10.44`.

The compact authoritative evidence, binary identities, isolation bounds,
confidence intervals, and reproduction commands are in
`experiments/leopard2/direct_repair/results/equal_rounded_avx2_authoritative_20260731.json`.

## What Leopard2 adds

The benchmark covers only the bit-identical legacy-high profile, but Leopard2
also adds architecture absent from main:

- Versioned high and low wire profiles. High uses `T=ceil_pow2(R)`,
  `N=ceil_pow2(K+T)`, parity coordinates `0..R-1`, originals `T..T+K-1`, a
  shortened systematic tail, and punctured parity `R..T-1`. Low uses
  `P=ceil_pow2(K)`, `N=ceil_pow2(P+R)`, originals `0..K-1`, a shortened tail to
  `P`, and parity starting at `P`; it supports `R>K`.
- Active-parent Algorithm 5 high-rate message-only decoding and Algorithm 4
  low-rate decoding, with the generic `O(N log N)` decoder retained as a
  fallback and oracle.
- Immutable reusable codecs and erasure plans, allocation-free encode and
  reusable-plan execution with caller scratch, no-loss no-op behavior, explicit
  scratch queries, arbitrary GF8 byte tails, complete-symbol native GF16
  lengths, a versioned padded-odd GF16 layout, parity subsets, batch APIs, and a
  persistent worker pool.  The one-shot decode wrapper allocates its plan.
- Direct paths for XOR parity, `K=1`, bounded tiny systematic encoding, and up
  to four missing originals in small codes.
- Runtime-isolated scalar, SSSE3, AVX2, and explicit AVX-512 fixed-operation
  backends with startup self-tests instead of raising the entire library's ISA
  floor; native NEON remains open production work.

The low profile is not included in an exact-main performance ratio: main
rejects `R>K`, and comparing a different generator/wire profile would not be a
like-for-like compatibility benchmark.

## Current GF8 R=32 regression closure (2026-08-18)

The original performance atlas exposed a finite set of GF8 legacy-high AVX2
regressions against exact Leopard main.  After the targeted fixes, one fresh
predeclared confirmation reran 97 former-loss and boundary-control workloads
at `R=32`, 64- and 1024-byte shards, and deterministic one, two, ten-percent,
and maximum-loss patterns.  Each cell used 75 independent balanced rounds;
one-shot decode kept setup charged.

All 194 exact-main comparisons passed with a lower 95% confidence bound above
1.0.  The weakest encode result was 1.0332x, CI95 [1.0314,1.0350], and the
weakest setup-inclusive decode result was 1.0560x, CI95
[1.0544,1.0577].  All 194 identical-binary control intervals were wholly
inside the predeclared `[1/1.02,1.02]` equivalence band.  Five contaminated
attempts were retained as rejected evidence and retried under the unchanged
external-isolation rule; all 7,275 accepted rounds recorded zero sibling
activity.

This is a targeted current-source regression closure, not a replacement for
the entire older atlas and not a claim about GF16, low-rate, large-shard, batch,
or multicore regimes.  The compact checkpoint and exact reproduction command
are in
`experiments/leopard2/gf8_high_encode/results/current_atlas_final97_confirmation_checkpoint_20260818.json`
and
`experiments/leopard2/optimization_log/38-current-atlas-final97-regression-closure.md`.

## Historical findings and current gaps

The 2026-07-16 audit found no high/low coordinate-map or Algorithm 4/5
arithmetic mismatch.  It did find the following four production gaps, but all
four have since been implemented; they are retained here to prevent the older
benchmark snapshot from being mistaken for current source:

1. Low encode formerly copied all `P` coefficient shards for every nonempty
   parity block.  Production now evaluates directly from immutable coefficient
   pointers through an out-of-place first layer.  `leopard-79h.26.1` remains
   open only for its authoritative isolated crossover/nonregression evidence.
2. Dense locator setup formerly used the ambient field order.  Production now
   restricts Walsh setup to active parent `N`; product-tree and epsilon locator
   construction remain optional comparison experiments under
   `leopard-79h.29` rather than missing production correctness work.
3. Specialized decode formerly staged an `N`-shard workspace.  Production now
   normally uses `min(N,2P)` for low and `min(N,2T+L)` for high where side
   tiling is profitable, while a narrow measured high-rate region deliberately
   retains `N`.  Ragged public-coordinate staging is bounded to one 64-byte
   tile per coordinate.  Some generic/fallback configurations still exceed the
   plan's aspirational `O(min(P,T))` scratch target and output pruning remains
   under `leopard-79h.26.4`.
4. Backend selection formerly required a process-global forced diagnostic
   build.  Contexts in one production binary can now independently select
   immutable qualified scalar, SSSE3, AVX2, or explicit AVX-512 operation
   tables.  The context's persistent worker pool remains lazily mutable.  Native
   NEON remains production work under `leopard-79h.13.7`.

Current core gaps are the exact-main encoder and dispatcher crossovers under
`leopard-79h.38`, sparse encoder output pruning under `leopard-79h.26.4`, and
allocation-free prewarm, static scheduling, affinity, 128-core scaling, and
NUMA policy under `leopard-79h.14`, `.23`, and `.24`.  Exact-size stable
profiles, streaming/incremental encode, serialized code identity, and CUDA are
separate experiments.  CUDA remains a correctly isolated, default-off build
scaffold and is intentionally sequenced after the production CPU codec.

## Evidence and reproduction

The committed evidence bundle is under
`experiments/leopard2/main_compare/results/620667c-amd9950x3d-cpu15/`:

- `manifest.json`: SHA-256
  `485dc78d2924d5e138dcc96d1b01ec59a3ed60014f1b5308d81c8ef8d9ad9967`
- `raw.json`: SHA-256
  `3b8586a1ddd964037c3d56e2dbda9909d54b4db576da0adf24d6d919fe286bc2`
- `summary.json`: SHA-256
  `f0c99d16f34418c6f13efdb5bb960dc695724a99cf208c5fb16bd7d0fdff0cc6`
- `isolation.json`: SHA-256
  `5ca4cec26af2c3394aaf9ef017ef6fbbb57a5b6638ed0945e0b69510013308e8`

`experiments/leopard2/main_compare/README.md` contains the exact fresh-build,
reservation, run, and verification commands. From a different machine, replay
the retained streams and analysis structurally with:

    python3 experiments/leopard2/main_compare/run_abba.py verify \
      --manifest experiments/leopard2/main_compare/results/620667c-amd9950x3d-cpu15/manifest.json \
      --no-current-input-check

Omitting `--no-current-input-check` additionally requires the original exact
source and build paths and re-hashes their full closure.

Before opening the replaceable reservation file, the current collector holds a
nonblocking exclusive flock on the owned `/run/user/UID` directory inode for
the complete campaign. C7 and the backend-butterfly collector hold the same
stable anchor, so renaming or recreating a reservation file or its containing
directory cannot split current cooperating runners onto independent file
locks. This stable layer conservatively serializes all current Leopard2
evidence campaigns for the UID, including campaigns on disjoint pairs. The
retained evidence schema and historical portable replay are unchanged; this
anchor is runtime authority rather than a new wire claim.
