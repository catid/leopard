# Leopard2 versus exact Leopard main

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
