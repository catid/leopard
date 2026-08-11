# 35 — R10 one-shot schedule elision

## Hypothesis

The Algorithm 4/5 kernels already reduce the transform depth to the smaller
side, but an exact-byte one-shot decode cannot amortize compilation of an exact
input/output dependency graph. At the R10 1024-codeword size, executing mature
regular transforms may cost less than compiling and executing the pruned graph.
Caller-created reusable plans should keep the graph because their setup is
amortized.

## Result

Landed for pure-AVX2 GF8 exact-byte one-shot plans through 1024 bytes:

- all native `LOW_V1` Algorithm 4 profiles;
- native-high T=32/N=256 profiles with K=97..224, R=17..32, and 2..R-1
  missing originals.

The existing maximum-loss and T=64 policies are unchanged. The 1025-byte
neighbor retains schedules, as do all public reusable plans and attribution
mode 2.

## Evidence

The clean commit is `a0d781c`. A five-pattern RS(256,K) screen used 1024-byte
shards, exactly R erasures, 25 timing samples, reuse 512, one thread, and a
frozen pure-AVX2 executable. Setup-inclusive specialized/generic geometric-mean
speedups were:

| K | 8 | 16 | 32 | 64 | 128 | 192 | 224 | 240 | 248 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| speedup | 8.000x | 3.254x | 1.876x | 1.597x | 1.337x | 1.359x | 1.584x | 1.604x | 2.791x |

Every exact Leopard-main comparison also wins: 1.130x at K128 (different low
wire profile), then 1.149x, 1.398x, 1.469x, and 2.489x at K192, K224, K240,
and K248. The wire-identical high-profile digests match exact main.

The T=32 boundary screen has 81 candidate/control pairs across
K={97,160,224}, R={17,24,32}, three loss counts, and 576/1024/1025 bytes.
One-shot speedup is 1.734x geometric mean at 576 bytes and 1.462x at 1024;
the weakest target is 1.229x. The 1025-byte control is 1.000x geometric mean,
and reusable execution is neutral at every size. All 81 digest pairs match.

Release tests passed 6/6 and focused ASan+UBSan tests passed 2/2. The complete
summary and compressed raw JSON are in `experiments/leopard2/r10_results/`.

## Interpretation

This closes the apparent contradiction with R10. The specialized arithmetic
was faster; transient schedule construction hid that advantage in first-use
measurements. The fix changes setup policy only, not formulas, coordinate maps,
parity bytes, or reusable execution.
