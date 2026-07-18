# Leopard2 AVX-512BW/VBMI codec experiment

## 2026-07 full-operation register-file follow-up

The earlier full-codec candidate below changed only the fused radix-four
range operation.  Assembly attribution subsequently found that the dominant
GF16 high-rate encoder operations were instead the two-way and final
accumulating butterflies: the portable AVX2 translation unit has 16 YMM
registers and repeatedly spills fixed-multiplier table vectors, while exact
Leopard main built with `-march=native` uses the AVX-512VL register file and
keeps those vectors in registers.

A distinct full-operation candidate now recompiles the reviewed AVX2 operation
table in `Leopard2BackendAVX512.cpp` with AVX2 plus AVX-512F/BW/VL and
`-mprefer-vector-width=256`.  It retains 256-bit loads, stores, and arithmetic;
the additional ISA contract is used for registers YMM16 through YMM31.  Runtime
qualification requires the corresponding CPUID bits and XCR0 state.  The
candidate shares the immutable AVX2 GF8/GF16 nibble tables, so it does not add
a second roughly 8 MiB GF16 table.  `AUTO` remains AVX2 and the candidate must
be selected explicitly while promotion evidence is collected.

The backend startup KAT, explicit-context concurrency tests, odd-byte GF16,
R=1, direct-encode, pruned-transform, decode-schedule, failure-injection, and
portable-ISA gates pass.  The archive audit treats the new object as a distinct
ISA member, rejects EVEX in the AVX2 member, permits only reviewed 256-bit EVEX
forms in the AVX-512VL member, and rejects ZMM and GFNI instructions.  Assembly
inspection of the accumulating GF16 butterfly confirms that its vector loop
uses the expanded register file without fixed-table stack spills.

The following same-binary measurements are directional only: the machine was
busy with independent sanitizer/fuzzer work, so they are not promotion
evidence.  They compare explicit AVX-512VL with explicit AVX2 using identical
Leopard2 code and data.

| Cell | AVX2 encode | AVX-512VL encode | Time reduction |
| --- | ---: | ---: | ---: |
| K=1000 R=200, 64 KiB | 8407 us | 7486 us | 11.0% |
| K=4096 R=512, 4 KiB | 2255 us | 1961 us | 13.0% |

Exact-Leopard-main ABBA evidence on an isolated core remains the promotion
gate.  Until it is recorded, this full-operation candidate is experimental and
does not alter `AUTO`.

Two alternative AVX2 table layouts were rejected during attribution.  Direct
32-byte row loads duplicated the fixed table and caused 3.5% to 4.5%
neighboring regressions.  A four-table packed-symbol design was 51% to 58%
slower.  Neither variant remains in the production source.

## Earlier narrow-candidate disposition

The full-codec GF16 AVX-512VL/BW integration candidate is **rejected**.  It
did not meet the experiment's promotion rule: a credible improvement of at
least 5% in the target regime, no unexplained neighboring regression above
2%, and closure of the corresponding exact-Leopard-main encode gap.  No
AVX-512 backend, public backend selector, CMake source, or runtime-dispatch
change from that narrow follow-up was promoted.  The later full-operation
candidate described above is a distinct implementation and decision.

The candidate isolated a 256-bit AVX-512VL/BW GF16 radix-four kernel in its
own ISA translation unit.  It used the expanded vector register file while
retaining YMM-width memory operations and the legacy wire representation.
`AUTO` remained AVX2; the experimental backend required explicit selection.
Exhaustive kernel tests covered all 65,535 ordinary GF16 logarithm values in
each multiplier position, forward and inverse radix-four forms, and scalar
equality.  Backend KAT, context, public API, direct encode, portable-ISA,
tests-off archive, and benchmark-parser gates passed.  Parity and recovered
data digests matched exact main in every retained timing invocation.

Pinned full-codec results used an AMD Ryzen 9 9950X3D, CPU 4, reserved SMT
sibling 20, three independent ABBA rounds, nine samples per child, two
warmups, loss count 8, and plan reuse 8.  Ratios are exact-main time divided
by candidate time, so values above one favor the candidate.  Exact main was
commit `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`; the experiment and hardened
evidence-runner commits were `66701c69074385a62bae8edf66038c67a76684a7`
and `eba41f7c47d860dfcb03d9199473780a3aeac0e7` respectively.

| Cell | Encode, 95% CI | Decode first use, 95% CI | Decode reuse 8, 95% CI |
| --- | ---: | ---: | ---: |
| K=1000 R=200, 64 KiB | 0.917 [0.914,0.919] | 2.981 [2.970,2.992] | 3.000 [2.989,3.012] |
| K=4096 R=512, 4 KiB | 1.018 [1.012,1.024] | 2.026 [2.006,2.045] | 2.156 [2.136,2.177] |
| K=1000 R=200, 32 KiB neighbor | 0.973 [0.910,1.039] | 1.939 [1.935,1.944] | 1.963 [1.958,1.968] |

Thus the candidate still encoded the primary K=1000/R=200 cell 8.35% slower
than exact main.  It closed the K=4096/R=512 gap, but only by 1.8%, whose
lower confidence bound was far below 5%.  The neighboring encode result was
inconclusive and did not exclude a regression.  The 1.94x to 3.00x decode
speedups confirm that the experiment did not regress Leopard2's specialized
decoder, but they do not justify a backend whose purpose was to close the
remaining encode gap.  A preliminary same-binary comparison also showed only
about 2% over AVX2 at K=1000/R=200, consistent with the exact-main result and
below the kill threshold.

Two initial combined campaigns were deliberately rejected rather than quoted:
the reserved sibling accumulated 14 non-idle jiffies during an overlapping
parallel build and 3 during the coordinated retry.  Splitting the cells into
shorter independent campaigns produced zero-sibling-work manifests.  The
ignored research cache retains the rejected `failure.json` files and the
three verified manifests:

- `.research/leopard2/gf16-avx512-eba41f7-exact-main-cpu4-20260718/`
  (`failure.json` SHA-256
  `885f97c799c86d23f1221adce09425e0cf4c6150bcb1f30ad0e642aa683da9d1`);
- `.research/leopard2/gf16-avx512-eba41f7-exact-main-cpu4-retry1-20260718/`
  (`failure.json` SHA-256
  `9b56551d251fcc5029744c0b0b3d8d1bc0d24f5bbfe0756077387f25dd38f40e`);
- `.research/leopard2/gf16-avx512-eba41f7-k1000-r200-64k-cpu4-20260718/`
  (`manifest.json` SHA-256
  `2c8bae43f6b2bd9c5323f5152b484a6c049bcd8a3dcd912be846e42b72ae8a2a`);
- `.research/leopard2/gf16-avx512-eba41f7-k4096-r512-4k-cpu4-20260718/`
  (`manifest.json` SHA-256
  `dd121b9d38f8b6a01d6f33a87721b7218e33c65bcc6ecfe1245f4dce78d92907`);
- `.research/leopard2/gf16-avx512-eba41f7-k1000-r200-32k-cpu4-20260718/`
  (`manifest.json` SHA-256
  `45b8d46d4146502e9d840d4b1eff85b6ffe5330c71abec18aabaa14607e29d34`).

## Earlier standalone checkpoint

The earlier isolated GF8 microkernel experiment remains useful as a warning
against promoting from a favorable synthetic chain.  Its AVX-512BW
four-parity codec-like chain was 37.9% to 119.2% faster than AVX2 over 64 B
through 16 MiB shards, and its radix-4 and radix-8 chains favored AVX-512BW at
every tested size.  Those results justified the bounded full-codec follow-up
above; they did not predict its end-to-end gain.

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
