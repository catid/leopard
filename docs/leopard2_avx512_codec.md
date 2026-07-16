# Leopard2 AVX-512BW/VBMI codec experiment

## Disposition

AVX-512BW nibble kernels are a strong production-integration candidate, but
this checkpoint does **not** enable them in the library.  On the tested AMD
Ryzen 9 9950X3D, the complete four-parity codec-like chain was 37.9% to 119.2%
faster than its AVX2 counterpart over 64 B through 16 MiB shards.  Radix-4 and
radix-8 chains also favored AVX-512BW at every tested size.  Promotion still
requires runtime dispatch, full Leopard2 encode/decode measurements, another
AVX-512 microarchitecture, and a cleaner repeat of the noisy 16 MiB codec
cell.  In other words: advance AVX-512BW to production integration testing,
but do not make it the default from this standalone result alone.

The tested AVX-512VBMI four-quadrant 256-byte lookup is rejected.  It was
slower than the AVX-512BW two-nibble lookup in every multiplication,
butterfly, and codec cell.  It also regressed against AVX2 for a 64 B radix-4
(-40.4%), a 64 B radix-8 (-22.6%), and a 16 MiB radix-4 (-2.7%).  The extra
three masked byte permutes, high-bit comparisons, table loads, and larger code
outweighed removal of the nibble split.

The 512-bit copy/XOR loops are useful supporting kernels, but their production
decision should be made with the AVX-512BW codec backend rather than from
standalone bandwidth cells.

## Scope and arithmetic oracle

The experiment is isolated under `experiments/leopard2/avx512_codec/`; it does
not include or call optimized Leopard arithmetic and does not change CMake or
runtime dispatch.  The scalar oracle reconstructs the legacy GF8
representation from:

- irreducible polynomial `0x11d`;
- Cantor basis bytes `1, 214, 152, 146, 86, 200, 88, 230`;
- polynomial-basis multiplication followed by an independently constructed
  polynomial-to-Cantor inverse map.

Its 256-by-256 multiplication table has FNV-1a hash
`0x5a9c078c013b3e83`, matching the separately implemented GFNI experiment.
The AVX2 and AVX-512BW paths use Leopard-style low/high nibble tables.  The
VBMI path splits each multiplier's complete 256-byte row into four 64-byte
tables, selects by the input's high two bits, and uses the low six bits as
`VPERMB` indices.

The fused radix-4 and radix-8 kernels implement the legacy inverse-butterfly
shape:

    y ^= x
    x ^= coefficient * y

The codec-like chain reads eight independent source shards, evaluates four
fixed GF8 linear combinations (including zero and one coefficient
specializations), performs a fused radix-4 output transform, and writes four
parity-like outputs.  Its logical traffic model counts 32 source reads and
four output writes per shard byte.  This is deliberately more representative
than `mul_mem`, but is not claimed to be a complete Leopard encoder.

## Correctness and safety

GCC 13.3.0 compiled the source with `-O3 -std=c++17 -march=x86-64-v2
-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Werror`.
The same warnings passed in the ASan+UBSan build.  Both ordinary and sanitized
validation reported:

- 65,536 independent scalar coefficient/input products;
- 197,376 SIMD product bytes: all 256 coefficients, every input byte, and a
  257-byte unaligned tail through AVX2, AVX-512BW, and AVX-512VBMI;
- 66 backend/length cases covering copy, XOR, fused radix-4, fused radix-8,
  and the complete chain;
- lengths `0,1,2,3,7,8,15,16,17,31,32,33,63,64,65,127,128,129,255,256,257,1023`;
- deliberately different unaligned input and output offsets;
- no ASan, leak, or UBSan finding.

## Pinned measurements

The retained measurements used physical CPU 0 (its SMT sibling is CPU 16)
from the allowed `0-31` CPU set.  All other project benchmark and validation
jobs were explicitly held during this phase.  The host has 16 physical cores,
one NUMA node, AVX2, AVX-512F/BW/VBMI, and a 128 MiB aggregate L3 split across
two cache instances.  This machine cannot provide the requested 128-core or
multi-NUMA evidence.

Each of two complete runs measured 240 cells: six operations, four backends,
and ten sizes from 64 B through 16 MiB.  Each cell has nine samples and median
and MAD, with iterations chosen to execute at least 32 MiB of logical traffic
per sample unless one invocation was already larger.  The table below uses
the median of the two retained run medians; positive percentages mean lower
time than AVX2.

| Operation | Bytes | AVX-512BW vs AVX2 | AVX-512VBMI vs AVX2 |
|---|---:|---:|---:|
| fixed multiply | 64 | +18.1% | +8.2% |
| fixed multiply | 16 KiB | +174.8% | +112.1% |
| fixed multiply | 1 MiB | +88.3% | +63.2% |
| fixed multiply | 16 MiB | +100.3% | +54.8% |
| radix-4 | 64 | +12.5% | -40.4% |
| radix-4 | 16 KiB | +84.7% | +31.2% |
| radix-4 | 1 MiB | +64.5% | +12.4% |
| radix-4 | 16 MiB | +39.1% | -2.7% |
| radix-8 | 64 | +21.1% | -22.6% |
| radix-8 | 16 KiB | +112.0% | +33.0% |
| radix-8 | 1 MiB | +69.5% | +7.8% |
| radix-8 | 16 MiB | +36.6% | +4.8% |
| codec chain | 64 | +77.8% | +15.2% |
| codec chain | 16 KiB | +119.2% | +33.1% |
| codec chain | 1 MiB | +72.3% | +7.7% |
| codec chain | 4 MiB | +70.8% | +5.1% |
| codec chain | 16 MiB | +37.9% | +5.4% |

Median within-run MAD was 0.28% in run 1 and 0.34% in run 2.  Most key cells
were stable, but the 16 MiB AVX-512BW codec medians differed by 29.8% between
runs.  That cell is reported rather than discarded and is a mandatory rerun
in full-codec integration.  Even its slower retained run still beat the
corresponding AVX2 run, but the combined percentage is not promotion-quality
evidence by itself.

## Frequency and register pressure

The readable `amd-pstate-epp` interface reported the `powersave` governor and
a 600 MHz to 5.752 GHz range.  Median post-cell `scaling_cur_freq` snapshots
were 5.511 GHz for AVX2, 5.286 GHz for AVX-512BW, and 5.235 GHz for VBMI across
both runs.  These snapshots suggest roughly a 4-5% wide-vector frequency
effect, already included in elapsed times, but they are not APERF/MPERF
measurements and should not be interpreted as a precise downclock ratio.
Hardware counters were unavailable (`perf_event_paranoid=4`).

Static GCC output illustrates the register-pressure tradeoff.  The AVX2
codec function is 1,453 bytes with 14 vector stack references, while
AVX-512BW is 1,497 bytes with four and VBMI is 1,594 bytes with four.  The
radix-8 functions are 946/1,045/1,255 bytes with 4/1/1 vector stack
references respectively (AVX2/BW/VBMI).  Counts conservatively include tail
paths, so they diagnose generated-code pressure rather than proving every
reference is a hot-loop spill.  The retained assembly summary makes these
counts reproducible.

## Reproduction and retained evidence

From the repository root:

    experiments/leopard2/avx512_codec/run.sh

Set `LEO2_AVX512_CPU=N` to select another allowed physical CPU and `CXX` to
select the compiler.  The runner uses temporary binaries under `/tmp`, removes
them on exit, runs strict and sanitized validation, performs two pinned runs,
merges results deterministically, summarizes generated assembly, and writes
`SHA256SUMS`.

Important retained artifacts are `results_run1.csv`, `results_run2.csv`,
`comparison.csv`, `assembly_summary.csv`, validation outputs, compiler/CPU
metadata, and hashes.  The tested repository baseline was
`50e7858ed9945d108b8035d7f2cfe7cc563582f1`.

References:

- Intel Intrinsics Guide:
  https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
- Intel architecture manuals:
  https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html
- Leopard source and legacy GF8 representation:
  https://github.com/catid/leopard
