# GF8 K12/R8 exact AVX2 high encoder

## Result

**Promoted at 256 and 1024 bytes.** Commit `64a5dd933e1901fe50d4d18ef1ed5c35d0d3e1ee`
adds an exact legacy-high GF8 circuit for `K=12,R=8,T=8` and a bounded
packed-call terminal. Against exact Leopard main commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, public loss-free parity encoding
is 1.2436x faster at 256 bytes and 1.1161x faster at 1024 bytes. Both
nine-round confidence intervals clear the five-percent promotion threshold.

| shard bytes | Leopard main | Leopard2 | main / Leopard2, 95% CI |
| ---: | ---: | ---: | ---: |
| 256 | 0.13720 us | 0.11033 us | **1.2436x** `[1.2207,1.2668]` |
| 1024 | 0.45386 us | 0.40666 us | **1.1161x** `[1.1114,1.1208]` |

The same executable contains diagnostic words that turn only the packed
terminal off. The selected route is 1.8763x faster at 256 bytes and 1.2061x
faster at 1024 bytes than that ordinary-validation control. This distinction
matters: the arithmetic circuit alone does not pay for sorting `K+R` public
ranges, reconstructing transform geometry, and staging padded pointer arrays.

## Exact circuit

This is the existing legacy wire profile, not a new code. For `K=12,R=8`,
`T=8`, `N=32`, `D=24`, and `S=12`. Coordinates are parity `0..7`, actual
message `8..19`, and shortened zeros `20..31`. Only two of the three message
blocks contain data.

For each 32-byte AVX2 tile, the circuit:

1. inverse-transforms the first eight-message block with the legacy shift;
2. inverse-transforms only the four active rows of the second block;
3. folds the known-zero upper radix-four group into the first accumulator with
   the derived legacy log-68 factor;
4. executes the exact T=8 forward schedule, replacing its three zero
   multiplications with XORs; and
5. scatters completed parity pairs as soon as their live ranges end.

The caller's first four parity rows are used as controlled spill storage after
the public alias checks have completed. The kernel itself accepts arbitrary
shard pointers and uses unaligned loads and stores. Packed layout is required
only for the fast public validation terminal; detached and sparse layouts
retain fully validated fallbacks. The compiler emits a 1505-byte kernel with
50 `vpshufb`, 28 table broadcasts, and two YMM spill slots.

The public scratch query remains 4672 bytes at B=256 and 16960 bytes at
B=1024. The exact byte kernel does not consume shard-data scratch, but keeping
the established query preserves fallback and binding geometry. This milestone
therefore claims an execution-time win, not a scratch-contract reduction.

## Frozen campaign

The candidate was rebuilt from clean commit `64a5dd9`, copied read-only, and
hashed before and after timing. Its SHA-256 is
`da9118650725bcd4dd3a39bae90b3b659fa3ab66d4c9cfdf77f54ef9ef01176b`.
The independently linked main binary is
`78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93`.

Each cell used nine alternating ABBA/BAAB rounds, two process medians per side
per round, 31 retained samples, 16 warmups, reuse 8192, batch one, and one
codec thread. Timing ran on CPU 30 while sibling 14 was reserved. Normal agent
affinity excluded both CPUs; every one of the 162 accepted rounds accumulated
zero non-idle sibling jiffies. There were no discarded attempts. All 648
original-data and parity digests matched. Confidence intervals are Student-t
intervals over the nine independent mean log contrasts.

The selector was also toggled in the same binary at the three K/R neighbors
and five byte neighbors. Every confidence interval spans parity:

| inert control | mode 0 / mode 1, 95% CI |
| --- | ---: |
| K11/R8/B1024 | 0.9984x `[0.9950,1.0018]` |
| K12/R7/B1024 | 1.0006x `[0.9974,1.0038]` |
| K13/R8/B1024 | 0.9977x `[0.9927,1.0027]` |
| K12/R8/B255 | 1.0024x `[0.9943,1.0105]` |
| K12/R8/B257 | 1.0042x `[0.9978,1.0106]` |
| K12/R8/B512 | 1.0040x `[0.9960,1.0120]` |
| K12/R8/B1023 | 0.9915x `[0.9769,1.0064]` |
| K12/R8/B1025 | 1.0144x `[0.9924,1.0369]` |

The exact-main neighbor map remains useful rather than flattering. K11/R8,
K12/R7, and K13/R8 are still slower than main at 1024 bytes; K11/R8 and
K12/R7 are substantially slower at 256 bytes. Those gaps are now tracked by
`leopard-79h.38.5.10.56` instead of widening a shape-specific selector without
evidence.

## Correctness and safety

The focused test uses an independent systematic generator-matrix oracle and
covers packed and unaligned slabs, sparse parity sets, detached source and
output rows, one-item batches, reusable bindings, forced mature fallbacks,
overlap errors, exact scratch, and failure atomicity at both selected byte
counts. It passed in Release, GCC 13 and Clang 18 strict-warning builds, and a
Clang 18 ASan+UBSan+LSan build. The sanitizer test peaked at 53 MiB RSS; all
builds were serialized after this host's earlier OOM events.

## Rejected variants

- A packed-address-only fixed circuit lost to main at 1024 bytes and could not
  serve arbitrary caller pointers. It was removed.
- A two-pass tail-first circuit reduced register pressure but was 2.55 percent
  slower than the selected one-pass controlled-spill circuit in the four-round
  preliminary comparison. It was removed.
- A high-level stagewise helper did not improve the complete public call. It
  was removed rather than left as dead production code.

The compact machine-readable checkpoint, including all identities, ratios,
neighbors, validation, and the retained raw-bundle hash, is
`experiments/leopard2/gf8_high_encode/results/`
`t8_k12r8_exact_checkpoint_20260810.json`.
