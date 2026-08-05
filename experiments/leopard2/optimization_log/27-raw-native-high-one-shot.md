# Scratch-resident tiny native-high one-shot decode

## Problem and candidate

The public `leo2_decode` wrapper historically created and destroyed a complete
heap-owned decode plan before executing it.  In the small native legacy-high
region, that setup dominated the useful byte work and left the representative
`K=16`, `R=8`, eight-loss, 64-byte cell behind exact Leopard1 even though the
reusable Algorithm 5 executor was faster.

Commit `a3ab73d8dcc415106c31a1507ef660cd61de0d5e` adds a bounded raw path for
GF8/AVX2 legacy-high parents with `N=32`, `T=8`, `K=9..16`, `R=5..8`, losses
`5..R`, and shard sizes through 256 bytes.  It constructs the deterministic
received set, erasure mask, active-parent locator, block prefixes, output
descriptors, and scatter metadata inside caller-provided scratch, then invokes
the existing tiled Algorithm 5 kernel.  It performs no allocation.  Complete
64-byte prefixes reveal directly into application outputs; only a ragged tail
uses the already reserved output slots and a final gather.  Every other
profile, field, backend, forced diagnostic mode, shape, and byte size keeps the
ordinary plan path.

## Correctness and safety

The allocation-audited public test covers 768 eligible combinations over every
`K=9..16`, `R=5..8`, loss count `5..R`, bytes
`1,63,64,65,255,256`, and varied original/parity loss shapes.  A second 168-case
suite compares recovered values against both the independent systematic
generator oracle and a reusable forced-Algorithm-5 plan.  The 257-byte boundary
proves fallback.  Focused Release, Clang ASan+UBSan, TSan, and strict GCC
warning builds pass, including short/misaligned scratch, pointer-presence
mismatches, overlap rejection, overflow, insufficient data, and failure
atomicity.  Eligible calls record zero hot-path allocations.

## Same-source attribution

The final pure-AVX2 executable was measured against itself with diagnostic mode
2 selecting the heap transient plan and mode 3 selecting the scratch-resident
path.  Three `control,candidate,candidate,control` rounds used CPU 28, reserved
sibling 12, 31 retained samples, four warmups, reuse 200, and the canonical
campaign lock.  The sibling accumulated zero non-idle jiffies.  Ratios below
are control time divided by candidate time.

| Bytes | Speedup | 95% CI | Control us | Candidate us |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1.3248x | 1.2985--1.3517 | 1.0389 | 0.7841 |
| 63 | 1.3395x | 1.3099--1.3697 | 1.0021 | 0.7518 |
| 64 | 1.3417x | 1.2609--1.4276 | 0.9543 | 0.6992 |
| 65 | 1.2514x | 1.2078--1.2966 | 1.4102 | 1.1248 |
| 256 | 1.2237x | 1.2009--1.2470 | 1.2498 | 1.0257 |
| 257, inert | 1.0014x | 0.9732--1.0304 | 1.7188 | 1.7191 |

The eligible boundary cells all clear the five-percent promotion threshold;
the first ineligible byte is neutral and its confidence interval spans one.

## Exact Leopard1 comparison

Both implementations were rebuilt from clean committed sources with an
explicit pure-AVX2 ceiling and frozen read-only before timing.  Exact Leopard1
is commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  Three ABBA rounds on CPU
29 reserved sibling 13 also recorded zero sibling non-idle jiffies.  Ratios are
exact-main time divided by Leopard2 one-shot time.

| Bytes | Speedup | 95% CI | Leopard1 us | Leopard2 us |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1.0629x | 1.0482--1.0779 | 0.8289 | 0.7791 |
| 63 | 1.1238x | 1.1114--1.1363 | 0.8283 | 0.7399 |
| 64 | **1.1885x** | **1.1752--1.2019** | 0.8264 | 0.6980 |
| 65 | 0.8289x | 0.8152--0.8429 | 0.9503 | 1.1455 |

The intended 64-byte regression is therefore closed.  The sharp reversal at
65 bytes is a separate algorithmic boundary: Leopard2 runs a second padded
transform tile for a one-byte tail.  It is the next optimization target; a
direct repair of only the ragged tail is the leading candidate.  Exact main's
nonmultiple cells process a zero-padded physical multiple of 64, so those rows
are application-equivalent rather than equal physical work.

Machine-readable identities, validation counts, results, and limitations are
in
`experiments/leopard2/direct_repair/results/raw_native_high_one_shot_20260805.json`.
