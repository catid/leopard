# Experimental FF8 XOR-circuit backend

## Status and goal

This branch adds a separate, experimental CPU backend for Leopard's FF8 path.
It preserves the existing FFT-based encoder and decoder, but replaces fixed
GF(256) payload multiplication with generated binary XOR circuits. It is not
the default backend and does not change the data output of valid
`leo_encode()` or `leo_decode()` calls. This branch also hardens the packed
APIs' documented count/work-count validation and fixes a `k=1` no-loss decode
null-recovery bug; those changes affect invalid or formerly crashing calls,
not valid packed codewords. The experimental backend does not implement FF16,
CUDA, or AVX-512 four-buffer fusion.

The payload loops use loads, stores, XORs, copies, zeroing, and address/loop
arithmetic. They do not use payload-indexed multiplication tables, byte
shuffles, GFNI, CLMUL, integer field multiplication, AND-based general field
multiplication, or a dynamically indexed gate interpreter. These circuits
implement fixed coefficients; they are not a variable-by-variable GF(256)
multiplier. Scalar field operations remain in initialization and decoder
metadata calculations.

## Native eight-plane shard layout

An experimental B-byte shard is one ordinary contiguous allocation divided
into eight equal planes:

```text
plane_bytes = B / 8
plane j     = shard + j * plane_bytes, j = 0..7
```

For group byte `g`, bit lane `lane`, and coordinate bit `j`:

```text
packed[8*g + lane] bit j = plane[j][g] bit lane
```

Thus one byte position across the eight planes represents eight parallel
GF(256) elements. `B` retains Leopard's positive multiple-of-64 requirement,
so every plane is a multiple of eight bytes and the portable `uint64_t` tail is
always exact.

The API does not transpose at its boundary. Original shards remain systematic:
their bytes are not modified or rearranged by encoding. Recovery bytes in the
native plane layout are not directly byte-compatible with packed Leopard
recovery bytes. Applying the documented 8x8 bit transpose converts between the
two layouts. Tests and the optional packed-boundary benchmark perform this
transpose outside the native hot path.

Both representations use Leopard's existing Cantor-basis coordinate bits, so
no field-basis conversion is involved.

## Field model and generated circuits

The generator reproduces `LeopardFF8.cpp` exactly:

- polynomial: `0x11D`;
- Cantor basis: `1, 214, 152, 146, 86, 200, 88, 230`;
- 256-entry logarithm and exponent tables, including Leopard's redundant log
  representation of exponent zero.

For a fixed nonzero coefficient, GF(256) multiplication is an 8x8 binary
linear map. Column `j` is the coordinate byte produced by multiplying `1 << j`
by the coefficient. The generator reduces each invertible matrix to identity
with reversible row additions and reverses those additions to obtain literal
CNOT operations:

```text
wire[destination] ^= wire[source]
```

It also constructs the complete 16x16 forward and inverse butterfly maps. For
each butterfly it compares the synthesized full-map circuit with the direct
construction, then chooses by gate count, dependency depth, and lexical order.
Adjacent duplicate gates cancel. Every emitted specialization names its live
wires `x0` through `x7` and `y0` through `y7`; there is no runtime gate array.

The checked-in generated file is
`generated/LeopardFF8XorCircuits.inl`. Its deterministic SHA-256 circuit
checksum is:

```text
e6947bc7968024941c1ed95dc1b7f607e2c44d9c1a55abaeb7cb528c4ce42ff7
```

Generated gate statistics:

| Circuit family | Minimum | Maximum | Average | Maximum depth | Average depth |
|---|---:|---:|---:|---:|---:|
| Multiply | 0 | 43 | 31.671875 | 35 | 24.2578125 |
| Forward FFT | 8 | 51 | 40.000000 | 16 | 12.2734375 |
| Inverse FFT | 8 | 51 | 40.000000 | 16 | 12.2734375 |

### The two meanings of 255

The dispatch tables deliberately remain separate:

- multiplication log `255` is the same nonzero multiplier as log `0`, namely
  multiplication by one; its circuit has zero gates but still copies when the
  source and destination differ;
- butterfly skew `255` is a sentinel that omits the multiply-add term, while
  retaining the eight `y ^= x` gates in the appropriate forward or inverse
  order.

Conflating these meanings would corrupt both decoding and sentinel
butterflies.

## Kernel and transform architecture

There are three coefficient-specialized whole-buffer dispatch tables: multiply,
forward butterfly, and inverse butterfly. The selected function pointer is
called once for a whole-buffer operation. Inside that function, the loop loads
eight named plane values for multiplication or sixteen named plane values for
a butterfly, executes literal XOR expressions, and stores them.

`Auto` deliberately retains the established baseline selection while the
AVX-512 modes remain opt-in. For a selected mode, the available chunk order is:

1. AVX-512 ZMM (64 bytes) or AVX-512VL YMM (32 bytes), when explicitly selected;
2. AVX2, 32 bytes from each plane;
3. 128-bit SIMD, 16 bytes from each plane;
4. portable `uint64_t`, 8 bytes from each plane.

Each smaller path handles only the exact remainder of one plane. Multiplication
loads all eight source planes before any store, so `source == destination` is
supported.

The backend copies Leopard's FF8 metadata initialization and uses a padded skew
array so the transform's historical one-before-logical-begin view remains a
valid pointer. The packed FF8 and FF16 tables now use the same safety pattern.
Encoder chunking, zero padding, transform truncation, work ordering, formal
derivative, decoder ErrorBitfield selection, and the `k=1`/`r=1` special cases
remain intact.

This first experiment executes two-way butterflies. It intentionally does not
use the packed backend's fused two-layer four-buffer path: four plane-sliced
FF8 buffers would require 32 live vector registers, while AVX2 provides 16.
This adds intermediate load/store passes and is an important performance
disadvantage, not a benchmark artifact.

### Choosing lower-cost locator constants

The decoder's evaluated error locator is defined only up to one common nonzero
scale. The backend scans all 255 possible common logarithmic shifts after
locator evaluation and minimizes:

```text
(total multiplier gate count, total multiplier depth, numeric shift)
```

The score includes present recovery inputs, present original inputs, and final
inverse-locator multipliers for missing originals. Missing and padded zero
buffers are excluded. Applying the same field factor to every received input
passes unchanged through the linear IFFT, formal derivative, and truncated
FFT; the inverse shifted locator cancels it at reveal. Transform skew constants
are fixed and are never changed by this optimization. This reduces generated
XOR gate counts and dependency depth where the algebra leaves a free choice;
it does not reduce the fixed eight-plane loads and stores. Even an identity
multiplier copies eight planes when source and destination differ.

## Experimental C API

Include `leopard_ff8xor.h` and call `leo_init()` before using the backend:

```c
unsigned work_count = leo_ff8xor_encode_work_count(k, r);
LeopardResult result = leo_ff8xor_encode(
    buffer_bytes, k, r, work_count, original_data, work_data);
```

The corresponding decode functions are
`leo_ff8xor_decode_work_count()` and `leo_ff8xor_decode()`. Work counts and
recovered-output indices match the existing Leopard API rules. The caller
provides an array of `work_count` work pointers, and every work allocation is
`buffer_bytes` bytes. Encode recovery shards occupy `work_data[0..r-1]`.
Decode marks unavailable original or recovery inputs with `NULL`; recovered
original `i` is returned in `work_data[i]`. All original, recovery, and work
pointers refer to native eight-plane shards.

The first experimental encode or decode call after `leo_init()` initializes
ff8xor metadata through a C++11 thread-safe lazy boundary. This keeps the
generated circuit object out of packed-only static-library clients; it does
not relax the requirement to call `leo_init()` first. Work-count helpers
return zero for invalid counts or transforms outside the FF8 range.

The experimental API accepts only parameter sets whose padded FF8 transform
size is at most 256 and returns `Leopard_TooMuchData` otherwise. It retains the
existing count, pointer-array, insufficient-recovery, initialization, and
multiple-of-64 size validation conventions.

## Regeneration and tests

The normal library build consumes checked-in generated code and does not
require Python. With Python 3 available:

```sh
python3 tools/generate_ff8_xor_circuits.py
python3 tools/generate_ff8_xor_circuits.py --check
cmake --build build --target generate_ff8_xor_circuits
cmake --build build --target check_ff8_xor_circuits
```

`--check` returns nonzero when the checked-in output differs from regenerated
output.

Build and run the finite automated test with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The test exhaustively checks all 256 logarithmic multipliers against all 256
coordinate bytes. It validates all butterfly skews on basis, edge, and
deterministic random states in both directions and proves the two operations
are mutual inverses. It verifies transpose helpers, API errors, every forced
common locator shift from 0 through 254, Portable/SIMD128/AVX2 modes and tails,
in-place multiplication, packed parity and decode equivalence, native round
trips, and packed-backend regression cases. End-to-end coverage includes
`(k,r)` values `(4,2)`, `(8,2)`, `(16,4)`, `(32,8)`, `(64,16)`, `(128,32)`,
and `(128,128)` at 64, 1024, 4096, and 65536 bytes, plus the `k=1`, `r=1`,
non-power-of-two recovery counts, no-loss, mixed-loss, and
insufficient-recovery edge cases.

## Benchmarking

The finite benchmark compares the packed backend and native ff8xor in the same
process with identical deterministic data. Allocation and correctness checks
are outside encode/decode timing. It reports median and best time plus input
and output MB/s.

```sh
./build/bench_leopard_ff8xor --quick
./build/bench_leopard_ff8xor
./build/bench_leopard_ff8xor --include-transpose
./build/bench_leopard_ff8xor --csv
./build/bench_leopard_ff8xor --quick --json --abba --counters
./build/bench_leopard_ff8xor --quick --json --abba --cache-color
```

Native rows exclude transpose time. Rows named `ff8xor_packed_boundary`
include packed-to-plane input transposes and plane-to-packed outputs. The full
mode covers `(k,r)` values `(8,2)`, `(16,4)`, `(32,8)`, `(64,16)`, `(128,32)`,
and `(128,128)`; buffer sizes 1 KiB, 4 KiB, 64 KiB, and 1 MiB; and each unique
loss count in `1`, `min(4,r)`, `max(1,r/2)`, and `r`.
Microbenchmarks report packed multiply-add, XOR-circuit multiply, forward and
inverse two-way butterflies, and both transpose directions. The reported
speed ratio is `packed time / ff8xor time`, so a value greater than one means
the XOR backend is faster.

For repeatable optimization work, the benchmark attempts to pin itself to the
first CPU in its allowed affinity set. `--cpu N` selects a specific allowed
logical CPU and `--no-pin` disables pinning. `--abba` measures matched packed
and native end-to-end calls in A-B-B-A order; standalone microbenchmarks and
optional packed-boundary rows remain sequential and are labeled as such. ABBA
rows report twice the requested sample count because each round contains two
samples of each backend. Allocation and equivalence checks remain outside the
timed regions.

### Optional transform-buffer cache coloring

`--cache-color` tests an allocation-only mitigation for the native layout's
worst plane-stride conflict. When `plane_bytes` is a multiple of 4096, all
eight plane starts in one shard have the same conventional 64-byte-line L1D
set index. A two-buffer butterfly can therefore present sixteen live lines to
one set; this exceeds the 12-way L1D associativity of the implementation host.
The benchmark leaves other sizes on the normal compact allocator because
extra coloring can add page/TLB cost without addressing a conflict.

[`tests/FF8XorCacheColoredBuffers.h`](../tests/FF8XorCacheColoredBuffers.h)
over-allocates a selected shard and exposes a 64-byte-aligned, still-contiguous
`buffer_bytes` region at a deterministic offset within a 4 KiB period. The
rank-six linear color map uses masks `1,2,4,8,16,32,21,42`; every index-bit
column is nonzero, so every radix-2 partner through transform size 256 has a
different color. Decode work uses salt 5. Exhaustive startup checks cover every
transform-path recovery padding (`m=2..128`) and valid original count,
source-to-work copies, final
recovery scaling pairs, all bounded radix-2 partners, address bounds, and all
eight plane-start colors. The previously considered salt 63 is not used: it
collides for `m=64`, original index 64.

This option changes benchmark allocation only. It does not alter
`leo_ff8xor_encode()`, `leo_ff8xor_decode()`, caller pointers, or the native
eight-plane format. Each exposed shard remains one contiguous allocation view;
the owner can reserve up to 8127 extra bytes, and allocation stays outside the
timed region. JSON metadata records `cache_coloring_requested`, while both JSON
and CSV rows expose the per-result `cache_coloring_applied` field. In
particular, microbenchmarks and non-aliasing sizes are never labeled colored.

On the implementation host, 12 outer rounds alternated uncolored-colored-
colored-uncolored and the reverse order, yielding 24 colored and 24 uncolored
process observations for every quick-mode row. Each process was pinned to CPU
0 and used the benchmark's internal packed/native ABBA order, one warm-up,
three iterations, and a 250 us minimum sample. The following values are the
ratio of aggregate median native `median_us` (uncolored / colored), exclude
transposes, and use the corrected decode salt 5:

| k,r | Encode | Decode loss 1 | Decode middle loss | Decode loss r |
|---|---:|---:|---:|---:|
| 8,2 | 1.220x | 1.193x | - | 1.220x |
| 32,8 | 1.356x | 1.171x | 1.204x (loss 4) | 1.205x |

Packed-normalized double ratios were 1.171x to 1.360x, supporting that the
gain was not merely process-to-process frequency drift. PMU events were
unavailable because this host has `perf_event_paranoid=4`, so this evidence is
elapsed-time only. The page-relative model and benefit must be rechecked on
CPUs with different L1 geometry. Reverse-store ordering and narrower chunk
schedules did not help in isolated trials. An internally padded layout helped
microkernels, but conversion-inclusive full-codec evidence was not strong
enough to retain it; persistent or tiled layouts remain separate experiments.

`--json` emits newline-delimited JSON with one environment record, three
circuit records, and one record per result. CSV and JSON include the generated
CMake configuration and a configured-flag summary (global, configuration,
C++ standard, FF8 XOR, and isolated AVX-512 flags), OS, affinity, SIMD selection,
deterministic schedule ID, modeled payload traffic, and selected decoder
locator shift. A verbose build log remains authoritative for implicit compiler
driver flags. The traffic model counts expected payload loads and stores,
including the three-buffer-byte skew-sentinel butterfly fast path; it does not
claim to measure cache-line transfers.

On Linux, `--counters` requests cycles, instructions, reference cycles, cache
references/misses, front- and back-end stalled cycles, L1D load misses, DTLB
load misses, and ITLB load misses through `perf_event_open`. PMU samples use
separate repetitions after wall-clock timing, so counter ioctls do not
contaminate reported times. Events are collected in PMU groups of at most three
so each group can schedule as a unit even on CPUs with a small programmable
counter file. Cycles, instructions, and reference cycles share one exact
payload window for IPC and frequency calculations; enabled/running time scales
multiplexed groups. IPC and effective GHz are derived only when their source
counters are available. Kernel policy, unsupported events, and
permission failures are represented as JSON `null`/empty CSV fields plus the
actual error string; the benchmark never substitutes invented zeroes.

### Deterministic schedule-frequency corpus

[`generated/FF8XorScheduleCorpus.json`](../generated/FF8XorScheduleCorpus.json)
is a checked-in 104-record corpus covering the full benchmark matrix. An
independent Python model reproduces the Cantor field tables, FFT skew table,
encoder chunking, deterministic erasures, locator construction, locator-shift
selection, natural FFT/IFFT calls, and decoder ErrorBitfield pruning. Records
contain input-scaling and recovery-scaling multiplier-log frequencies, FFT and
IFFT skew frequencies, direction-separated complete two-layer four-buffer
tuple frequencies, exact erasure indices, and modeled payload traffic. Combined
multiplier and tuple histograms are retained for aggregate weighting. Benchmark
`schedule_id` values join directly to corpus record IDs.

```sh
python3 tools/generate_ff8xor_schedule_corpus.py
python3 tools/generate_ff8xor_schedule_corpus.py --check
cmake --build build --target check_ff8xor_schedule_corpus
```

The corpus contains checksums of both its ordered schedule records and the
generated circuits. Fixed seeds and stable JSON key ordering make `--check`
detect changes to either schedule semantics or workload order.

### Assembly regression census

`tools/inspect_ff8xor_assembly.py` accepts the FF8 XOR object or Leopard static
archive. It inventories all 256 specializations in each available base,
AVX-512VL, and AVX-512 ZMM family, including code size, XOR2, XOR3
(`vpternlog[dq]` immediate `0x96`), loads/stores, calls, stack references, and
mnemonics. Direct-branch control-flow analysis distinguishes one-per-buffer
dispatch calls from calls inside cyclic payload loops. Strict mode fails on
missing mandatory families or specializations, payload-loop calls, possible
scaled stack indexing, direct or indexed static-table reads,
narrow dynamically indexed loads, shuffle/permute/blend/shift/gather/GFNI/
CLMUL/integer-multiply/vector-mask instructions, or non-XOR3 ternary truth
tables. Possible compiler-generated vector spills are reported separately;
`--fail-on-spills` makes those warnings fatal for a code-generation quality
gate. AVX-512-enabled CMake builds require the complete VL and ZMM family set.

```sh
python3 tools/inspect_ff8xor_assembly.py build/liblibleopard.a \
    --strict --fail-on-spills
cmake --build build --target inspect_ff8xor_assembly
cmake --build build --target check_ff8xor_assembly
```

The non-strict target writes `ff8xor_assembly_census.json` in the build
directory. Supported x86 Linux Release builds run the structural strict check
under CTest. The developer `check_ff8xor_assembly` target additionally requires
zero spills. On this machine, GCC's Release kernels pass that stronger gate.
Clang 18 produced genuine base-path loop spills in 104 of 256 FFT and 251 of
256 IFFT specializations, while its generated AVX-512VL and ZMM families had
none. Disabling Clang's loop and SLP vectorizers reduced code size and some
spills but did not eliminate those FFT/IFFT spills, so the experimental branch
reports them rather than treating correct code from that compiler as a failed
implementation.

### Standalone AVX-512 XOR3 experiment

`bench_leopard_ff8xor_xor3` isolates the question of whether one AVX-512
`vpternlogd` can profitably replace the two instructions needed for
`a XOR b XOR c`. It tests XMM, YMM, and ZMM widths and keeps three forms
separate: a forced control containing exactly two XOR instructions, ordinary
two-XOR intrinsic source whose lowering is left to the compiler, and an
explicit ternary-logic intrinsic with truth-table immediate `0x96`. The
compiler-auto form is important because an optimizing compiler may already
select `vpternlogd`; the forced form measures the actual two-instruction
counterfactual.

```sh
./build/bench_leopard_ff8xor_xor3 --verify-only
./build/bench_leopard_ff8xor_xor3 --quick
./build/bench_leopard_ff8xor_xor3
cmake --build build --target check_ff8xor_xor3_assembly
```

Runtime entry requires AVX2, AVX-512F, AVX-512VL, and matching OS-enabled
extended state. Unsupported builds or hosts report a CTest skip rather than
executing the isolated AVX-512 translation unit. `--verify-only` checks all
eight truth-table inputs, 10,000 deterministic random vectors for every form
and width, and streaming results without reporting timings.

Timed mode covers a dependent latency chain, four independent chains, and
four-array streaming working sets. Each explicit-versus-control comparison is
warmed first and sampled in A-B-B-A order on a pinned logical CPU. Wall-clock
time comes from `steady_clock`; fenced TSC samples provide a reference-clock
cross-check, and RDPRU APERF samples are printed when supported. TSC GHz is not
core frequency, and even APERF cannot remove thermal, scheduler, or wide-vector
frequency effects, so repeat runs and compare medians, best times, and the
reported clock samples rather than treating one best sample as authoritative.
Streaming `logical_GiB/s` counts output bytes; the printed traffic estimate is
four times that value for three loads and one store.

The strict assembly target inspects stable noinline probes. It requires two
XORs for the forced control, one `vpternlog[dq] 0x96` for the explicit form,
the requested vector width, and no calls or stack references. On the GCC 13.3
implementation host the complete probe sizes, including `endbr64` and `ret`,
were 13/13/17 bytes for forced XMM/YMM/ZMM and 12 bytes for each explicit
probe. These sizes are compiler-, ABI-, and hardening-dependent; they are a
regression signal, not a portable encoding-size guarantee.

The corrected full Release run on the implementation host found explicit
VPTERNLOG 1.903x--2.203x faster than forced XOR2 and 1.497x--1.552x faster than
compiler auto on the dependent chains. With four independent chains the
corresponding ranges were 2.094x--2.103x and 1.361x--1.451x. In contrast, the
15 streaming ratios against compiler auto ranged from 0.994x to 1.021x; the
largest stream gain against forced XOR2 was 1.114x for YMM at a 4 KiB aggregate
working set. Pairwise APERF was generally unchanged, so the streaming result
was not explained by systematic downclocking. Thus the register-only
microbenchmark confirms that ternary XOR can reduce dependency depth, but it
does not establish an end-to-end ff8xor codec gain. Generated-circuit changes
must still pass the full codec benchmark and assembly census before retention.

## Measurements from the implementation host

These results were collected on an AMD Ryzen Threadripper PRO 9985WX
(64 cores, 128 threads), Linux 6.8, GCC 13.3.0, and CMake 3.28.3.
Historical measurements below used `-march=native -Wall -Wextra -fopenmp -g
-O0 -O3 -std=gnu++11`; the experimental translation unit additionally used
`-mno-avx512f -fno-tree-reassoc`. A CMake typo had copied Debug flags into the
Release flags before appending `-O3`, so `NDEBUG` was unset. The guardrail work
now preserves the toolchain's standard Release flags (`-O3 -DNDEBUG`) and
embeds the configured flag summary in every CSV/JSON run. Treat the table below
as the pre-guardrail historical baseline rather than a result from the fixed
Release configuration.

Full-mode end-to-end rows used two warm-ups, seven measured samples, and enough
repeated calls per sample to reach at least 1 ms; microbenchmarks used three
warm-ups and 31 samples. The table reports median/best. Allocation and
correctness checks were outside timed regions. These 1 MiB rows use the native
plane layout and exclude transposes. Ratio is packed median time divided by
ff8xor median time.

| k/r/op/loss | Packed us | ff8xor us | ff8xor input MB/s | ff8xor output MB/s | Ratio |
|---|---:|---:|---:|---:|---:|
| 8/2/encode/0 | 262.272/257.813 | 710.877/677.106 | 11800.365/12388.914 | 2950.091/3097.228 | 0.369 |
| 8/2/decode/1 | 1237.686/1154.329 | 2851.991/2657.015 | 2941.316/3157.155 | 367.665/394.644 | 0.434 |
| 16/4/encode/0 | 591.707/552.136 | 2158.510/2144.689 | 7772.591/7822.680 | 1943.148/1955.670 | 0.274 |
| 16/4/decode/1 | 5519.352/5327.431 | 8056.096/7792.927 | 2082.549/2152.877 | 130.159/134.555 | 0.685 |
| 32/8/encode/0 | 2341.047/2069.966 | 5689.580/5516.548 | 5897.524/6082.505 | 1474.381/1520.626 | 0.411 |
| 32/8/decode/1 | 14807.675/14198.862 | 20530.905/20124.109 | 1634.338/1667.375 | 51.073/52.105 | 0.721 |
| 64/16/encode/0 | 6711.959/6606.179 | 13837.755/13727.067 | 4849.693/4888.798 | 1212.423/1222.200 | 0.485 |
| 64/16/decode/1 | 36589.695/36320.658 | 50921.519/49906.401 | 1317.888/1344.695 | 20.592/21.011 | 0.719 |
| 128/32/encode/0 | 20211.210/19872.827 | 35336.324/34073.130 | 3798.293/3939.108 | 949.573/984.777 | 0.572 |
| 128/32/decode/1 | 85031.532/83478.469 | 123994.538/118312.521 | 1082.449/1134.434 | 8.457/8.863 | 0.686 |
| 128/128/encode/0 | 36628.827/36064.503 | 72663.957/71723.351 | 1847.102/1871.325 | 1847.102/1871.325 | 0.504 |
| 128/128/decode/1 | 92633.997/92134.921 | 141694.073/139725.591 | 947.236/960.581 | 7.400/7.505 | 0.654 |

Across all 104 matched end-to-end cases, native encode ratios ranged from
0.239 to 0.887 (mean 0.585), or 1.13x to 4.18x slower. Native decode ratios
ranged from 0.376 to 0.721 (mean 0.557), or 1.39x to 2.66x slower. With both
packed-boundary transposes included, encode ranged from 0.052 to 0.291 (mean
0.129) and decode from 0.191 to 0.525 (mean 0.332). No measured end-to-end
ff8xor case was faster than packed Leopard.

Representative 1 MiB microbenchmarks were:

| Operation | Coefficient/circuit | Median/best us | Median input/output MB/s |
|---|---|---:|---:|
| Packed multiply-add | log 22; 32 gates in corresponding XOR circuit | 15.849/15.668 | 132320.115/66160.057 |
| ff8xor multiply | log 22; 32 gates, depth 26 | 19.463/19.359 | 53875.842/53875.842 |
| ff8xor FFT butterfly | skew 1; 40 gates, depth 11 | 98.409/97.988 | 21310.552/21310.552 |
| ff8xor IFFT butterfly | skew 1; 40 gates, depth 11 | 98.846/98.485 | 21216.337/21216.337 |
| Packed to plane transpose | standalone | 260.158/258.899 | 4030.540/4030.540 |
| Plane to packed transpose | standalone | 42.797/42.638 | 24501.204/24501.204 |

The packed microkernel is multiply-add while the circuit kernel is
multiply-only, so their input-byte rates are not directly equivalent. The
sampled low/average/high gate circuits also show that these large-buffer
microbenchmarks are dominated by memory behavior and instruction scheduling,
not monotonically by source gate count.

## Representative optimized code inspection

`objdump -dr -C` inspection used strict AVX2 (`-march=haswell`) and the final
native flags above. Representative multiplier/skew 42 has 34 gates at depth 26
for multiply and 43 gates at depth 13 for each butterfly. With GCC's default
reassociation, 217 of 256 FFT and 241 of 256 IFFT specializations spilled;
`-fno-tree-reassoc` eliminated vector stack spills in all 256 multiply, FFT,
and IFFT specializations inspected.

The final pure-AVX2 loops for coefficient/skew 42 contained 23 `vpxor` plus
eight stores for multiply, and 42 `vpxor` plus sixteen stores for each
butterfly, along with vector loads and scalar pointer/loop control. There were
no vector stack references, calls inside the chunk loop, `pshufb`, integer
multiply, `pand`, gathers, payload-table accesses, or dynamic gate indexing.
Dispatch compiled to one indexed function-pointer load and tail jump before
the whole-buffer loop. `-mno-avx512f` prevents the Zen 5 compiler from replacing
the experiment with AVX-512VL `vpternlog` and high vector registers.

## Known disadvantages and next experiment

- Applications must retain the plane layout to avoid boundary transpose cost.
- Each whole-buffer operation incurs one coefficient function-pointer dispatch.
- Two-way butterflies lose the packed backend's existing fused two-layer
  memory pass.
- Sixteen live AVX2 vectors leave no spare vector register, so compiler spills
  remain possible for some generated circuits.
- Static specialization of all multiplier and butterfly coefficients increases
  source, object, and instruction-cache footprint.
- The implementation-host static archive was 62 MiB and the experimental
  benchmark executable was 30 MiB. Lazy experimental initialization kept the
  packed-only benchmark at 762 KiB, so packed-only static clients do not pull
  the generated object.
- Packed compatibility requires explicit 8x8 transposes.
- The inherited CMake build still applies `-march=native` globally. The
  `uint64_t` implementation is a functional fallback on the build host, but a
  binary built this way is not yet guaranteed to run on an older x86 host;
  baseline-ISA object isolation is tracked as a required follow-up.
- This backend is FF8-only and CPU-only.

The highest-value next optimization is AVX-512 four-buffer fused generation:
32 architectural vector registers can hold four plane-sliced FF8 buffers and
recover the fused two-layer memory behavior. Other follow-ups include
tuple-specialized four-way circuits, CUDA plane kernels, persistent plane
pipelines, packed-boundary transpose/locator fusion, and synthesis optimized
for depth and register pressure rather than gate count alone.
