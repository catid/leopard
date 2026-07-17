# Leopard2 versus exact Leopard main

Date: 2026-07-16

This comparison measures Leopard2 commit
`620667c93fad4eb9d97c8c5ccf7dcf6994914346` against the exact detached
Leopard default-branch baseline
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. The result is mixed: reusable
Leopard2 decoding is materially faster in the intended high-rate GF16 region,
but the legacy-compatible Leopard2 encoder is slower in every measured cell.
XOR and balanced full-loss decode also regress.

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
- Immutable reusable codecs and erasure plans, allocation-free execution,
  no-loss no-op behavior, explicit scratch queries, arbitrary byte tails,
  parity subsets, batch APIs, and a persistent worker pool.
- Direct paths for XOR parity, `K=1`, bounded tiny systematic encoding, and up
  to four missing originals in small codes.
- Runtime-isolated scalar, SSSE3, and AVX2 fixed-operation backends with startup
  self-tests instead of raising the entire library's ISA floor.

The low profile is not included in an exact-main performance ratio: main
rejects `R>K`, and comparing a different generator/wire profile would not be a
like-for-like compatibility benchmark.

## Confirmed gaps

The algorithm audit found no high/low coordinate-map or Algorithm 4/5 arithmetic
mismatch, but it confirmed these production gaps:

1. Low encode copies all `P` coefficient shards into a second `P` workspace for
   each nonempty parity block. Tracked as `leopard-79h.26.1`.
2. Dense locator setup still runs ambient-field 256/65,536-entry Walsh
   transforms rather than always scaling with active parent `N`. Tracked as
   `leopard-79h.29.1`.
3. Follow-up candidate `leopard-79h.26.2.2` removes the `N` full-shard
   specialized workspace where smaller: low uses `min(N,2P)` and high uses
   `min(N,2T+L)`. Forced generic retains `N`; the materialized specialized path
   also remains when tiling would use at least `N`. Follow-up candidate
   `leopard-79h.26.2.3` splits ragged execution so its `K+R` staging term is
   fixed at one 64-byte tile per public coordinate instead of full-shard slots.
   Correctness and scratch-slope tests pass; isolated promotion timing remains
   open, so the benchmark numbers above do not yet claim these candidates.
4. Runtime backend choice is process-global; lower scalar/SSSE3 contexts in the
   same production binary are not supported. Tracked as `leopard-79h.13.1`.

General dependency pruning, active topology/NUMA scheduling, exact-size stable
profiles, native NEON, streaming/incremental encode, serialized code identity,
and CUDA kernels also remain open. CUDA is currently only a correctly isolated,
default-off build scaffold.

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
