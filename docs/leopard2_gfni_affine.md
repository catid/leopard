# Leopard2 GFNI affine experiment

## Decision

This is a positive but still experimental result.  Direct
`VGF2P8AFFINEQB` fixed multiplication is correct in Leopard's actual GF8
Cantor representation and materially faster than AVX-512 nibble tables in a
mid-cache region.  It is not promoted to production by this experiment alone:
the codec-like gain is not monotonic, 16 MiB multiplication is about 2.2%
slower, and this standalone chain is not a complete Leopard transform.

The production hypothesis was that a fixed GF8 multiplier is an 8-by-8 linear
map over GF(2), so one affine instruction can replace two nibble shuffles and an
XOR.  The promotion threshold was at least 5% in the target end-to-end region,
with neighboring regressions no worse than 2%.  The kill criterion was no such
region after including matrix setup, zero/one specialization, complete chains,
and wide-vector behavior.  The kill criterion was not met: the 4 KiB and 16 KiB
chain cells improved by 23.2% and 18.5%, with adjacent 1 KiB and 64 KiB cells at
+7.7% and +0.3%.  A future production bead should integrate the ZMM affine
kernel behind dispatch, then repeat complete encode/decode benchmarks on more
than one microarchitecture before promotion.

## Representation and empirical matrix order

The scalar oracle reconstructs Leopard's field without calling optimized
Leopard code:

- polynomial `0x11d`;
- Cantor basis bytes `1,214,152,146,86,200,88,230`;
- polynomial multiplication followed by the inverse Cantor linear map.

For every one of the 64 bits in the GFNI matrix operand, the experiment applies
that bit to each of the eight input basis vectors and observes the hardware
output.  The discovered operand rule is:

    operand_bit = 8 * (7 - output_bit) + input_bit

Thus Leopard's identity map encodes as `0x0102040810204080`.  All 256 fixed
multiplier matrices are generated from the eight scalar products of each
multiplier with the input basis vectors.  Their little-endian FNV-1a hash is
`0xc5ae814681e06043`; the full 256-by-256 scalar multiplication table hash is
`0x5a9c078c013b3e83`.

The alternate `VGF2P8MULB` experiment constructs a field isomorphism rather
than incorrectly treating Leopard bytes as AES-field bytes.  In the AES
`0x11b` field, `0x03` is the selected root of Leopard's `0x11d` polynomial.
The empirically encoded maps are:

    Leopard Cantor -> AES field: 0xefd0aaca822e5cae
    AES field -> Leopard Cantor: 0xffb08466dc727ea0

This route requires affine conversion, `VGF2P8MULB`, and inverse conversion.
For standalone ZMM multiplication it was 7.9-30.1% slower than direct affine
for 64-16,384 bytes, converging in bandwidth-bound cells.  It did not provide a
consistent chain advantage over direct affine.

## Correctness

The following all passed:

- every 256 multiplier by 256 input-byte pair, 65,536 cases total;
- raw affine XMM, YMM, and ZMM forms;
- basis-converted `GF2P8MULB` XMM, YMM, and ZMM forms;
- AVX2 and AVX-512 nibble-table references;
- non-vector-multiple tails through 257 bytes;
- explicit coefficient zero (`memset`) and one (`memcpy`) specializations;
- a 12-stage multiply/XOR chain containing general and one factors, with zero
  specialization tested separately so it cannot erase earlier chain work.

Zero specialization measured about 224-227 GiB/s and one specialization about
112-114 GiB/s at 64 KiB; all backends reach the same libc paths by design.

## Measurements

Host: AMD Ryzen 9 9950X3D, CPU 0, one NUMA node, GFNI, AVX2, AVX-512F/BW/VL.
GCC 13.3.0 flags were `-O3 -g -std=c++17 -march=native -Wall -Wextra
-Werror`.  The process affinity allowed CPUs 0-31.  Authoritative runs were
pinned to CPU 0 because concurrent work would invalidate cache-sensitive
latency and frequency results.  The `amd-pstate-epp` driver reported the
`powersave` governor and a readable range of 600 MHz to 5.752 GHz.  Hardware
performance counters were unavailable because `perf_event_paranoid=4`.

Each CSV cell contains nine timed samples and reports median and MAD.  Two
complete pinned runs were retained; the table below uses the median of their
medians.  Improvement is lower time relative to AVX-512 nibble tables.

| Work | Bytes | Direct GFNI ZMM improvement |
|---|---:|---:|
| fixed multiply | 64 | +7.2% |
| fixed multiply | 256 | +32.6% |
| fixed multiply | 1 KiB | +56.7% |
| fixed multiply | 16 KiB | +52.2% |
| fixed multiply | 64 KiB | +25.5% |
| fixed multiply | 1 MiB | +0.7% |
| fixed multiply | 16 MiB | -2.2% |
| 12-stage chain | 256 | +14.3% |
| 12-stage chain | 1 KiB | +7.7% |
| 12-stage chain | 4 KiB | +23.2% |
| 12-stage chain | 16 KiB | +18.5% |
| 12-stage chain | 64 KiB | +0.3% |
| 12-stage chain | 256 KiB | +7.0% |
| 12-stage chain | 1 MiB | +0.3% |

Run-to-run variation was visible in the smallest cells, which is another reason
not to promote directly from this standalone result.  Full raw samples,
repetition counts, per-run MAD, cycles per input byte, effective throughput,
and internal-pass throughput are retained under
`experiments/leopard2/gfni_affine/`.

## Reproduction and sources

From the repository root:

    experiments/leopard2/gfni_affine/run.sh

Override the pinned allowed CPU with `GFNI_CPU=N` and the compiler with `CXX`.
The runner compiles, executes exhaustive validation, performs two pinned
benchmark runs, merges them deterministically, and writes `SHA256SUMS`.
Experiment source hash: `7dafaf7ef4851beebabae05ef62af3b1ca1116a46520506d2d9c4ac60a0afb48`.
Repository baseline when the retained run was produced:
`570888f87406ac492cac0fe6a4d6bd600104ce37`.

Primary references used:

- Intel Intrinsics Guide:
  https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
- Intel GFNI technology guide:
  https://builders.intel.com/docs/networkbuilders/galois-field-new-instructions-gfni-technology-guide.pdf
- Intel architecture manuals:
  https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html
- Leopard source representation and nibble kernels:
  https://github.com/catid/leopard
