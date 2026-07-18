# Experimental FF8 XOR-circuit backend

## Status and goal

This branch adds a separate, experimental CPU backend for Leopard's FF8 path.
It preserves the existing FFT-based encoder and decoder, but replaces fixed
GF(256) payload multiplication with generated binary XOR circuits. It is not
the default backend and does not change the data output of valid
`leo_encode()` or `leo_decode()` calls. This branch also hardens the packed
APIs' documented count/work-count validation and fixes a `k=1` no-loss decode
null-recovery bug; those changes affect invalid or formerly crashing calls,
not valid packed codewords. The experimental backend does not implement FF16
or CUDA. Its generated AVX-512 four-buffer circuit corpus is now integrated
into the production transform path behind an explicit experimental selector.
The selector defaults to disabled; normal ff8xor and packed-Leopard behavior
therefore remains unchanged unless a caller deliberately enables the fused
path.

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

It also constructs the complete 16x16 forward and inverse butterfly maps. The
deterministic synthesis portfolio varies column order, pivot choice,
elimination direction, inverse/transpose orientation, commuting cancellation,
and exact shortest rewrites for windows touching at most four wires. An
in-place Boyar-Paar-style bidirectional greedy search uses legal CNOTs on both
sides of the matrix plus bounded two-gate lookahead. Multiplier candidates also
use a complete radius-three GL(8,2) breadth-first ball (52,277 maps) in a
meet-in-the-middle search that proves optimality through six CNOTs and rewrites
full eight-wire windows. A deterministic SSA linear-form cost rejects shorter
source programs that compile to more XOR work. Applying the greedy search to
16-wire butterflies reduced source XORs but made GCC spill 380
portable/SIMD128/AVX2 specializations, so those wider circuits deliberately
retain the verified zero-spill portfolio. Selection is by gate count,
dependency depth, then lexical order, subject to the compiled-work and
register-pressure gates. Every emitted specialization names its live wires
`x0` through `x7` and `y0` through `y7`; there is no runtime gate array.

The checked-in generated file is
`generated/LeopardFF8XorCircuits.inl`. Its deterministic SHA-256 circuit
checksum is:

```text
f03cc4d1230d128c6900745f3372b6ce946ede037904b887bc760abd557b8028
```

Generated gate statistics:

| Circuit family | Minimum | Maximum | Average | Maximum depth | Average depth |
|---|---:|---:|---:|---:|---:|
| Multiply | 0 | 23 | 19.152344 | 17 | 11.375000 |
| Forward FFT | 8 | 51 | 40.000000 | 14 | 10.859375 |
| Inverse FFT | 8 | 51 | 40.000000 | 14 | 10.859375 |

There are 4,903 multiplier gates, 10,240 forward gates, and 10,240 inverse
gates in total. Checked-in regeneration takes about 43.7 seconds and 59 MiB on
the implementation host. Compared with the original single-row-reduction
generator, multiplier gates fell from 8,108 to 4,903 (39.53%) and multiplier
depth from 6,210 to 2,912 (53.11%). Including its current cost-profile
provenance, the generated file is 1,280,484 bytes, 6.77% smaller than the
original 1,373,440-byte file. Relative to the preceding
greedy portfolio, bounded exact search reduced source multiplier gates by
1.13% and depth by 1.52%. A representative GCC AVX2 corpus compile reduced
`vpxor` instructions from 4,060 to 4,027 and text from 43,758 to 43,566 bytes,
with no stack spills; the guard rejects exact rewrites that do not improve its
deterministic compiled-work proxy.

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

Payload `KernelMode::Auto` deliberately retains the established baseline
selection while the AVX-512 circuit modes remain opt-in. This is distinct from
the boundary-transpose `Mode::Auto`, which selects the retained BITALG/VBMI
paths when their narrower feature contracts are available. For a selected
payload mode, the available chunk order is:

1. AVX-512 ZMM (64 bytes) or AVX-512VL YMM (32 bytes), when explicitly selected;
2. AVX2, 32 bytes from each plane;
3. 128-bit SIMD, 16 bytes from each plane;
4. portable `uint64_t`, 8 bytes from each plane.

On x86-64 GCC/Clang builds, the public API gateways, runtime dispatcher,
portable path, and SSE2 128-bit tail are compiled explicitly with
`-march=x86-64`. AVX2 and AVX-512 live in separate translation units with exact
source flags. MSVC x64 likewise keeps AVX2 in a separate `/arch:AVX2`
translation unit while its public dispatcher retains the compiler's SSE2
baseline. The baseline dispatcher checks CPUID, OSXSAVE, and XCR0 before
entering either targeted object; an unavailable forced mode resolves to a safe
available baseline (AVX2, 128-bit, or portable). `LEO_FF8XOR_ENABLE_AVX2=OFF`
omits both the isolated AVX2 circuit path and the AVX-512 circuit path that
depends on its tails.

Each smaller path handles only the exact remainder of one plane. Multiplication
loads all eight source planes before any store, so `source == destination` is
supported.

Encoder accumulation and formal-derivative additions resolve the SIMD mode
once for the complete batch, then group independent buffers four, two, and one
at a time. For buffers through 1 KiB, named two- and four-stream loops expose
independent load/XOR/store chains in portable, SIMD128, AVX2, AVX-512VL, and
AVX-512-ZMM forms. Larger buffers retain the one-time mode lookup but use the
better measured sequential vector-loop shape; interleaving at 4 KiB and above
was neutral to slower on the implementation host. Each individual
`source == destination` pair is supported, including odd batch tails. The
internal transform helper accepts independent stream pairs; ranges belonging
to different entries must not overlap, as at both production call sites.

The backend copies Leopard's FF8 metadata initialization and uses a padded skew
array so the transform's historical one-before-logical-begin view remains a
valid pointer. The packed FF8 and FF16 tables now use the same safety pattern.
Encoder chunking, zero padding, transform truncation, work ordering, formal
derivative, decoder ErrorBitfield selection, and the `k=1`/`r=1` special cases
remain intact.

The production transform has an opt-in AVX-512 ZMM radix-4 path for the 64
generated reachable skew tuples. `FourBufferMode::Disabled`, the default,
always retains the established two-way schedule. `FourBufferMode::Xor2`
selects literal two-input XOR circuits and `FourBufferMode::Xor3` selects the
generated mix of two-input XOR and `VPTERNLOGD 0x96`. This selector is
independent of `KernelMode`; fusion is attempted only when the resolved payload
mode is `KernelMode::Avx512Zmm`. AVX2 still has only the 16 vector registers
needed by one plane-sliced two-buffer butterfly and therefore remains two-way.
`FourBufferMode` is currently a developer-facing C++ control used by the test
and benchmark programs; it does not add a selector to the experimental C ABI.
Ordinary `leo_ff8xor_encode()` and `leo_ff8xor_decode()` calls therefore retain
the disabled default unless an embedding C++ experiment explicitly changes
that internal control.

The transform resolves one reachable tuple and calls one tuple-specialized
whole-buffer function before entering its offset loop. The loop itself has no
coefficient or gate dispatch: it loads the 32 named plane registers, executes
the literal generated circuit, and stores those 32 registers. A fused unit is
used only when both layers are complete, equivalently
`range + 2 * low_distance < count_truncated`, and every plane consists of one
or more complete 64-byte ZMM chunks. A partial truncated unit, an unknown
tuple, a sub-ZMM or non-64-byte plane tail, an unavailable build/CPU/OS state,
an unsafe tracked-buffer state, or a decoder layer omitted by its
`ErrorBitfield` falls back for the whole unit to the exact existing two-way
schedule. It never applies a four-buffer map to only part of a plane.

Forward fusion composes the natural higher-distance layer followed by its two
lower-distance pairs; inverse fusion uses the reverse natural ordering. The
tracked deferred-zero planner advances the same logical buffer states only
after a fused call succeeds. Decoder `ErrorBitfield` decisions remain distinct
for the higher layer and each lower pair, and fusion requires all applicable
decisions to select the complete unit. These constraints preserve transform
truncation, zero materialization, state identities, work ordering, and the
decoder's per-layer pruning. When fusion does not apply, its extra intermediate
load/store pass remains the same performance disadvantage documented for the
two-way implementation.

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

### Rejected locator-boundary fusion

A bounded experiment fused identity-scaled or structural-zero inputs directly
into the first decoder IFFT layer, retaining the staged scale/copy path for
every nonidentity pair. The bounded corpus required 768 generated cases (three
nonempty presence masks times 256 skews); a dense version covering both input
coefficients would require 16,646,400 maps (`255 * 255 * 256`). The candidate
was exact at the map level, but was rejected and removed from the runtime.

On the GCC 13.3 implementation host, text in
`LeopardFF8Xor.cpp.o + LeopardFF8XorAVX2.cpp.o` grew from 1,285,847 to
2,563,407 bytes: an increase of 1,277,560 bytes, or 1.993555x baseline. A
pinned 1 MiB AVX2 decode benchmark used four independent repetitions of 21
paired A-B-B-A rounds. Across the three rows containing one eligible fused
pair, equal-row/equal-repetition geometric mean staged/fused time was
0.988577x, so the candidate was about 1.16% slower. The small-transform result
was a repeatable 0.966218x, while the larger transforms were effectively
neutral. Transposes and allocation were excluded.

The checked negative evidence has checksum
`3d948e4159b338058a66d6fe5adc91d7d51035575641302a5633f2dae46d223d`.
Its independent checker reconstructs all 256 skews and three presence masks,
validates 16 basis vectors per map, checks the schedule-corpus cardinalities,
and recomputes the reported code-size and timing aggregates. The evidence also
records the removed C++ candidate's 11,520-case map census and a temporary
end-to-end harness over all 255 forced locator shifts. There is an important
limitation: the direct map proof defines and stores every destination plane,
but that temporary end-to-end harness started with zeroed work buffers. The
requested rerun with distinct nonzero poison and canaries was not completed
before the candidate was rejected. This is not presented as poisoned-work
end-to-end coverage; such coverage is required before reconsidering the
runtime hook.

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
require Python. The primary generator writes both the production XOR2 file and
the separately checked experimental XOR3 file. With Python 3 available:

```sh
python3 tools/generate_ff8_xor_circuits.py
python3 tools/generate_ff8_xor_circuits.py --check
python3 tests/test_ff8_xor3_schedule.py
python3 tests/test_ff8xor_generated_xor3.py
cmake --build build --target generate_ff8_xor_circuits
cmake --build build --target check_ff8_xor_circuits
```

`--check` returns nonzero when either checked-in output differs from regenerated
output. The four-buffer corpus and rejected locator-boundary evidence have
their own finite generation and checking commands:

```sh
python3 tools/generate_ff8_xor_four_buffer_circuits.py
python3 tools/generate_ff8_xor_four_buffer_circuits.py --check
python3 tests/test_ff8xor_four_buffer_circuits.py
python3 tools/inspect_ff8xor_four_buffer_avx512.py \
    --compiler g++ --compiler clang++-18 --strict
python3 tools/inspect_ff8xor_four_buffer_production.py \
    --archive build/liblibleopard.a --strict
python3 tools/evaluate_ff8xor_locator_boundary.py --check
cmake --build build --target check_ff8_xor_four_buffer_circuits
cmake --build build --target check_ff8xor_four_buffer_production_assembly
cmake --build build --target check_ff8xor_locator_boundary_evidence
ctest --test-dir build --output-on-failure \
    -R '(leopard_ff8xor_four_buffer_runtime|ff8xor_(generated_xor3|four_buffer_circuits|four_buffer_production|locator_boundary_evidence))'
```

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
and `(128,128)`; across these cases, selected buffer sizes include 64, 1024,
4096, and 65536 bytes. Additional cases cover `k=1`, `r=1`, non-power-of-two
recovery counts, no-loss, mixed-loss, and
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
./build/bench_leopard_ff8xor --quick --json --abba \
    --ff8xor-mode avx512zmm
./build/bench_leopard_ff8xor --quick --json --abba \
    --ff8xor-mode avx512zmm --four-buffer-mode xor2
./build/bench_leopard_ff8xor --quick --json --abba \
    --ff8xor-mode avx512zmm --four-buffer-mode xor3
```

`--ff8xor-mode` accepts `auto`, `portable`, `simd128`, `avx2`, `avx512vl`,
or `avx512zmm`. Selection occurs before timing and the requested value is
printed in human, CSV, and JSON metadata. A forced mode that is unavailable in
the build, on the CPU, or in the OS-enabled extended state exits with status
77; it never silently benchmarks a fallback. `auto` retains normal runtime
selection. `--four-buffer-mode` accepts `disabled`, `xor2`, or `xor3` and
defaults to `disabled`. Either active fusion mode requires an explicit
`--ff8xor-mode avx512zmm`; the benchmark rejects a non-ZMM combination instead
of printing a misleading fusion label. A requested ZMM mode that is not
available retains the status-77 behavior.

Human, CSV, and JSON results report the requested four-buffer mode, successful
fused-unit count, and its estimated payload bytes elided. The latter models
the intermediate loads and stores removed by successful complete radix-4
units; it is combined with, but remains separately visible from, the existing
deferred-materialization traffic estimate. Public ff8xor entry points reset
both diagnostic structures on every call, including validation errors and
special-case/no-loss returns, so benchmark rows cannot inherit stale fusion
counts.

Native rows exclude transpose time. Rows named `ff8xor_packed_boundary`
include packed-to-plane input transposes and plane-to-packed outputs. The full
mode covers `(k,r)` values `(8,2)`, `(16,4)`, `(32,8)`, `(64,16)`, `(128,32)`,
and `(128,128)`; buffer sizes 1 KiB, 4 KiB, 64 KiB, and 1 MiB; and each unique
loss count in `1`, `min(4,r)`, `max(1,r/2)`, and `r`.
Microbenchmarks report packed multiply-add, XOR-circuit multiply, forward and
inverse two-way butterflies, both transpose directions, and sequential versus
two/four-stream batched plain XOR. The reported speed ratio is
`packed time / ff8xor time`, so a value greater than one means the XOR backend
is faster.

For repeatable optimization work, the benchmark attempts to pin itself to the
first CPU in its allowed affinity set. `--cpu N` selects a specific allowed
logical CPU and `--no-pin` disables pinning. `--abba` measures matched packed
and native end-to-end calls in A-B-B-A order. The paired sequential/batched
XOR microbenchmarks and portable/automatic transpose pairs always use A-B-B-A
to control drift; other standalone microbenchmarks and optional packed-boundary
rows remain sequential. Every row
labels the actual order. ABBA rows report twice the requested sample count
because each round contains two samples of each implementation. Allocation and
equivalence checks remain outside the timed regions.

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
archive. It inventories all 256 specializations in each available base, AVX2,
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

`tools/inspect_ff8xor_xor_batch.py` separately requires the two- and
four-stream plain-XOR loops, the selected YMM/ZMM width, at least two or four
independent vector XORs in the cyclic component, and no calls, vector spills,
scaled stack indexing, RIP-relative tables, shuffles, or multiplication. It is
registered as `ff8xor_xor_batch_assembly` for supported Release builds and is
also available as `check_ff8xor_xor_batch_assembly`.

The non-strict target writes `ff8xor_assembly_census.json` in the build
directory. Supported x86 Linux Release builds run the structural strict check
under CTest. The developer `check_ff8xor_assembly` target additionally requires
zero spills. On this machine, GCC's Release kernels pass that stronger gate.
Clang 18 produced genuine loop spills in 167 of 256 base FFT, 254 of 256 base
IFFT, 70 of 256 AVX2 FFT, and 248 of 256 AVX2 IFFT specializations (739 total),
while its multiplier, AVX-512VL, and AVX-512 ZMM families had none. Disabling
Clang's loop and SLP vectorizers reduced code size and some spills but did not
eliminate those FFT/IFFT spills, so the experimental branch reports them
rather than treating correct code from that compiler as a failed
implementation.

`tools/inspect_ff8xor_baseline_isa.py` is a separate x86-64-v1 gate. It checks
the public gateways and the complete experimental dispatcher object, plus the
CPU detection and packed FF8/FF16 metadata initialization functions reached by
`leo_init()`. The latter remain in packed translation units, so their symbol
filters also include compiler-created FWHT/OpenMP clones. The check rejects
SSE3, SSSE3, SSE4, AVX/AVX-512, GFNI, CLMUL, BMI, and other post-v1 opcodes.
The only deliberate exception is `xgetbv` inside a named CPUID/OSXSAVE-guarded
feature probe. The M=1, encoder accumulation, and formal-derivative additions
call the experimental XOR dispatcher directly, so forced Portable does not
reach the packed `-march=native` XOR helpers.

```sh
cmake --build build --target check_ff8xor_baseline_isa
ctest --test-dir build --output-on-failure \
    -R 'ff8xor_baseline_isa_(inspector|census)'
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release \
    -DLEO_FF8XOR_ENABLE_AVX2=OFF -DLEO_FF8XOR_ENABLE_AVX512=OFF
```

The unsanitized suite also builds `check_ff8xor_baseline_isa` as a CTest case,
which protects the shell quoting of its demangled-symbol regular expressions
in addition to exercising the inspector directly.

The baseline sources append individually compiler-probed `-mno-*` feature
flags after `-march=x86-64`, so a caller's global `-mssse3` or host-native flag
cannot silently contaminate the dispatcher. Clang 18 Release with an
adversarial global `-mssse3` passed the strict artifact census with zero
violations. A separate Clang 18 Debug ASan+UBSan build passed all 18 registered
functional and parser tests at that checkpoint; later additions expanded the
registered suite. Sanitizers intentionally inject helper opcodes
and calls (Clang ASan inserts `lahf`), so CMake omits all production-shape
artifact assembly tests and custom targets when any single- or multi-config
C++ flag set enables a sanitizer. It prints that decision; parser unit tests
and all functional coverage remain enabled. A configure-only regression test
covers both Release single-config and Debug Ninja Multi-Config sanitizer flags.

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

### Generated explicit-XOR3 experiment

The generator also emits the separate checked artifact
`generated/LeopardFF8XorCircuitsXor3.inl`. It schedules literal named-wire
`XorValue` and three-input `Xor3Value` operations without changing the base
binary maps. `Xor3Value` is exact parity (`VPTERNLOG` immediate `0x96`), and
the deterministic schedule checksum is:

```text
a054d879212b80bdc03e418cdc3a6fe65b569a082784f7c4227b1146d4bce26e
```

The explicit schedules reduce source operation count as follows:

| Family | Source XOR2 gates | Scheduled XOR2 | Scheduled XOR3 | Total instructions | Reduction |
|---|---:|---:|---:|---:|---:|
| Multiply | 4,903 | 3,113 | 895 | 4,008 | 18.25% |
| FFT | 10,240 | 3,072 | 3,584 | 6,656 | 35.00% |
| IFFT | 10,240 | 3,072 | 3,584 | 6,656 | 35.00% |

Fewer source operations did not produce a codec win. On the GCC 13.3 host,
nine alternating pinned quick-mode processes per variant and width used
native buffers, internal A-B-B-A order, and packed-normalized end-to-end
times. Full explicit scheduling measured 0.981208x at AVX-512VL and 0.976871x
at AVX-512 ZMM, where values above one favor explicit scheduling. An FFT-only
variant measured 1.003649x and 0.992819x respectively, but its normalized FFT
microbenchmarks measured 0.994077x and 0.992180x. The inconsistent width-level
result and contrary microbenchmark do not justify retaining it.

The compiler already forms legal XOR3 instructions from the historical XOR2
source. Explicit scheduling reduced FFT code and instruction count, but grew
the multiplier and IFFT production census; all control and candidate assembly
variants nevertheless passed strict inspection with zero spills or inner-loop
calls. Two-buffer production therefore continues to include only the
historical XOR2 artifact. This negative result does not apply to the separately
generated radix-4 XOR3 schedules documented below. Checked evaluation checksum
`03412e947d515369eae7ca9929687e8f54fb89362e272f872e38fcee866420ba`
binds the artifacts, build, assembly totals, timings, and decision. Its raw
temporary JSONL and assembly reports were not retained, so their manifest
hashes authenticate the summary but cannot make those inputs reparsable.

### Generated and integrated AVX-512 four-buffer corpus

`tools/generate_ff8_xor_four_buffer_circuits.py` exhaustively enumerates all
21,845 valid FF8 `(k,r)` parameter pairs and finds 64 reachable three-skew
radix-4 tuples spanning 128 skew values. An independent cross-check against
the 104-case schedule corpus finds the same 64 tuples. For each direction and
tuple, the generator verifies direct butterfly composition, a composition of
the selected two-buffer portfolio, and synthesized whole-map candidates. It
then emits the better verified direct or composed circuit; no whole-map
candidate won. The checked metadata and named-register C++ share checksum:

```text
0301b24993ab889caa17aa00e64e969523ace698f63466833af500710244142b
```

The exhaustive valid-shape census contains 4,747,804 complete FFT radix-4
calls and 4,945,408 complete IFFT calls. A unit is counted only when both
lower-layer pairs execute (`range_start + 2 * distance < count_truncated`);
partial truncated units must retain the natural two-way schedule.

| Family | Lowering | Minimum ops | Maximum ops | Average ops | Total ops | Selected direct/composed |
|---|---|---:|---:|---:|---:|---:|
| FFT | XOR2 | 53 | 183 | 155.000 | 9,920 | 48 / 16 |
| IFFT | XOR2 | 53 | 183 | 155.000 | 9,920 | 53 / 11 |
| FFT | XOR2+XOR3 | 43 | 115 | 101.484 | 6,495 | 51 / 13 |
| IFFT | XOR2+XOR3 | 41 | 110 | 97.500 | 6,240 | 57 / 7 |

The selected FFT XOR3 schedules contain 3,425 ternary operations; IFFT
contains 3,680. Every generated body uses the 32 statically named wires
`a0..a7`, `b0..b7`, `c0..c7`, and `d0..d7`, with tuple dispatch intended
outside the vector-chunk loop.

All 32 ZMM registers are live, so the lowering contract matters. Ordinary
intrinsic wrappers let GCC and Clang reassociate the circuit DAG and were
observed to spill. The production translation unit therefore uses destructive
read/write `+v` inline-assembly wrappers for `vpxord` and
`vpternlogd $0x96`. The same wrappers were first proved on synthetic
frequency-representative and maximum-operation stress kernels, then the strict
production inspector was extended to all 256 compiled specializations: 64
tuples, both FFT directions, and both XOR2/XOR3 lowerings.

That full-object inspection found a subtler residency bug which the synthetic
stress sample had missed. Although all 32 values were initially loaded, GCC
kept their equivalence to the backing payload addresses and rematerialized
late circuit sources, producing 32 through 48 payload loads in a
specialization. Four zero-instruction barriers, each tying one group of eight
named registers with `+v`, make every loaded value opaque while staying below
GCC's inline-assembly operand limit. After that fix, both GCC 13 and Clang 18
compile all 256 production specializations with exactly 32 payload vector
loads and 32 stores, zero vector-stack references, zero calls, zero excess
reloads, the exact generated `vpxord`/`vpternlogd` counts, and no RIP-relative
tables, byte shuffles, GFNI, or CLMUL field operations. Tuple 63 could
legitimately let the compiler elide unchanged stores before the opacity
barrier; the production traffic model deliberately and conservatively charges
all four output buffers.

The GCC 13 object is 417,744 bytes (SHA-256
`72c9aa10a3ad41db4dccf6ef1474d47a41ca4408f12953b0c962241e34b224cf`).
Its per-specialization code sizes are 897--1,681 bytes (average 1,509.5) for
FFT and IFFT XOR2, 849--1,329 bytes (average 1,241.75) for FFT XOR3, and
833--1,305 bytes (average 1,221.625) for IFFT XOR3. The independently checked
Clang 18 object is 559,288 bytes (SHA-256
`f41bc4b79a47e9ff9738342f3798f9b6c7f68f99c9d46a5e0c16797dfe825311`)
and has the same exact memory/call/operation census. Only GCC
per-specialization code ranges were retained at this checkpoint; no Clang
range is inferred.

Production runtime tests compare every tuple/direction/lowering against the
two-way reference, cover unavailable/narrow/unknown/disabled fallbacks, and
exercise native encode/decode equivalence through truncation, tracked-buffer,
mixed-loss, and maximum-loss cases. These checks convert the earlier generated
corpus from synthetic evidence into an opt-in production implementation; they
do not make it the default.

#### Four-buffer performance checkpoint

The following is an integration checkpoint, not a broad promotion result.
GCC 13 Release runs used native plane buffers, pinned logical CPU 24, two
warm-ups and seven measured iterations per row, with adaptive inner calls
reaching at least 1 ms. Each size aggregates the full six encode `(k,r)` rows
or all 20 distinct decode `(k,r,loss)` rows by geometric mean. Ratios are
candidate time divided by the disabled two-way control time, so values below
1 are faster. Allocation, checking, and transposes were excluded.

| Shard bytes | Repeated XOR2 encode | Repeated XOR2 decode | Repeated XOR3 encode | Repeated XOR3 decode |
|---:|---:|---:|---:|---:|
| 1 KiB | 0.755x | 0.907x | 0.739x | 0.936x |
| 4 KiB | 0.784x | 0.901x | 0.793x | 0.907x |
| 64 KiB | 1.090x | 0.953x | 1.177x | 0.993x |
| 1 MiB | 1.086x | 0.988x | 1.136x | 1.001x |

For each lowering, the cells are the geometric mean of matching row times
from two disabled and two candidate reports before the cross-row aggregation.
These reports were regenerated after the final default-disabled fast-reject
ordering, so the control and candidate describe the checked-in code. The
result is useful but not uniformly favorable: XOR3 is best for 1 KiB encode,
while XOR2 is slightly better in the other aggregated small cells and in the
large decode cells. Both lowerings clearly regress large encode; large decode
ranges from a modest XOR2 win to an approximately neutral XOR3 result.
Consequently all builds still default to `FourBufferMode::Disabled`. A
possible policy restricted to payloads at or below 4 KiB is a separate tuning
experiment, not enabled by this checkpoint.

### Offline ISA-aware circuit costing

Circuit selection is reproducible and does not calibrate at runtime. The normal
generator reads `generated/FF8XorCostProfiles.json` and uses
`portable-default-v1`. That profile deliberately preserves the historical
gate-count, dependency-depth, and lexical ordering, so adding the cost model
does not silently change the checked-in circuits. Its modeled score is
provenance-only and comes after the lexical program in the key; it does not
break portable lexical ties. The generated C++ records the selected profile ID
and checksum; the current portable profile checksum is
`96a9fb14ba04b508831daa564dd1ecec7454eac234f4d12e80a75f4455070112`.
JSON, CSV, and human benchmark metadata print that provenance beside the
circuit checksum.

An explicit `amd-zen5-gcc13-avx2-avx512-v1` profile scores literal XOR2 count,
a deterministic distinct-linear-form XOR2 estimate, dependency depth,
live-range events, a fixed width-versus-register-budget overflow estimate,
fixed plane loads/stores, estimated code bytes, and estimated I-cache lines.
Its calibrated XOR3 weight is retained as provenance, but source generation
assigns every candidate in this machine-profile selection stage an XOR3 count
of zero: it does not predict which XOR pairs the compiler will combine into
`vpternlog`, nor does it simulate the compiler's register allocator. The
separate explicit-XOR3 artifact is a post-selection rescheduling experiment,
not an input to this cost-profile choice. Actual XOR2/XOR3 and spill counts
enter only the post-compile assembly evidence. These are known limitations of
the current model and reasons not to infer an AVX-512 benefit from its source
score.
Checked evidence includes the complete 256-coefficient AVX2
calibration, the standalone XOR3 calibration, the 104-case transform-schedule
corpus (including coefficient, skew, and four-buffer-tuple frequencies), and
all 362 equivalent-circuit timing records under `tools/profiles/`. The candidate
artifact retains summarized production-assembly totals and codec statistics,
plus SHA-256 hashes for all 14 codec JSONL inputs; those raw codec JSONLs are
not checked into the repository. Every retained JSON evidence file carries a
deterministic checksum. The evidence test also binds the pair timings to the
current portable generated-file hash and to the exact rejected machine file
summarized by the candidate evaluation.

The coefficient calibration improved Spearman correlation with measured
multiplier time from 0.323 for source gate count to 0.470 for the source-only
ISA model; all five deterministic coefficient-modulo-five slices improved.
This is only a correlation across different multiplication maps. It is not
evidence that the model can choose between equivalent schedules for one map.
The latter question was tested directly rather than inferred from that result.

The explicit machine profile changed 51 multiplier, 176 FFT, and 135 IFFT
schedules. A GCC 13.3 production census found no vector spills and the following
schedule-corpus-weighted changes versus the portable corpus:

| Family | Source gates | Modeled cost | AVX2 XOR2 | AVX2 code bytes | AVX-512VL XOR2 | AVX-512VL XOR3 |
|---|---:|---:|---:|---:|---:|---:|
| Multiply | +1.069% | -0.765% | -1.353% | -0.302% | -0.615% | -2.389% |
| FFT | 0.000% | -1.536% | -2.036% | -0.772% | +2.936% | -5.321% |
| IFFT | 0.000% | -1.451% | -0.480% | -0.011% | -0.271% | -0.691% |

Those compiled reductions did not translate into a retained optimization. A
same-process, pinned, 15-round ABBA harness first compared every changed pair's
output and then timed its literal named-register AVX2 kernels. The candidate
won 164 of 362 median comparisons, measured 0.986x by unweighted geometric
mean and 0.968x when weighted by the schedule corpus. Predictor-saving versus
measured-speedup Spearman was -0.099 (raw gate saving was -0.113), so even the
slightly less negative rank statistic has no useful selection meaning. The
family workload-weighted ratios were 1.052x multiply, 0.926x FFT, and 0.986x
IFFT.

Seven alternating quick codec repetitions then provided 98 observations over
14 cases. After normalizing each native result by the packed control from its
own build, the geometric mean of case medians was 0.993945x (7 of 14 case
medians won), a 0.605% slowdown. The run used native buffers and excluded
transposes. This is a small noisy result, but it is not a reproducible benefit;
combined with the pair test, it rejects the current machine-specific corpus.
The portable circuits remain the default. The checked codec evidence records
historical benchmark schema v1 inputs captured before traffic-accounting schema
v2. Those raw historical JSONLs are not retained, so current tests can verify
the checked summary and its 14 exact input hashes but cannot reparse those old
rows. For new inputs the evaluator supports both known schemas, requires exact
complete row sets, and cross-links current v2 circuit/profile provenance to the
generated files. It validates and binds the scheduled/elided/adjusted traffic
fields in v2 row descriptors, but performance normalization remains based only
on measured `median_us`; it does not traffic-adjust elapsed time.

Regenerate or check the deterministic profile and portable circuits with:

```sh
python3 tools/generate_ff8xor_cost_profiles.py --check
python3 tools/generate_ff8_xor_circuits.py --check
python3 tools/generate_ff8_xor_circuits.py \
    --cost-profile amd-zen5-gcc13-avx2-avx512-v1 \
    --output /tmp/LeopardFF8XorCircuits-machine.inl
ctest --test-dir build --output-on-failure -R ff8xor_cost
```

The calibration and candidate tools are offline developer tools; use their
`--help` output to record compiler, CPU pinning, repetitions, input artifacts,
and output path. Normal library builds consume only checked-in generated files
and do not require Python or execute calibration.

## Measurements from the implementation host

These results were collected on an AMD Ryzen Threadripper PRO 9985WX
(64 cores, 128 threads), Linux 6.8, GCC 13.3.0, and CMake 3.28.3.
The full reference build used GCC Release flags `-O3 -DNDEBUG`, C++11, global
`-march=native -Wall -Wextra -fopenmp`, baseline ff8xor flags
`-march=x86-64 -fno-tree-reassoc`, isolated AVX2 flags
`-march=x86-64 -mavx2 -fno-tree-reassoc`, and isolated AVX-512 flags
`-march=x86-64 -mavx2 -mavx512f -mavx512vl -fno-tree-reassoc`. Runtime auto
dispatch selected packed AVX2 and ff8xor AVX2. The benchmark was pinned to
logical CPU 0; PMU counters and cache coloring were disabled.

Full-mode end-to-end rows used two ABBA warm-up rounds, seven measured ABBA
rounds (14 observations per backend), and enough repeated calls per observation
to reach at least 1 ms; microbenchmarks used three
warm-ups and 31 samples. The table reports median/best. Allocation and
correctness checks were outside timed regions. These 1 MiB rows use the native
plane layout and exclude transposes. Ratio is packed median time divided by
ff8xor median time.

| k/r/op/loss | Packed us | ff8xor us | ff8xor input MB/s | ff8xor output MB/s | Ratio |
|---|---:|---:|---:|---:|---:|
| 8/2/encode/0 | 265.218/256.485 | 585.748/559.135 | 14321.2 | 3580.3 | 0.453 |
| 8/2/decode/1 | 1677.50/1361.24 | 2548.33/2402.69 | 3291.8 | 411.5 | 0.658 |
| 16/4/encode/0 | 770.209/554.854 | 1953.57/1862.43 | 8588.0 | 2147.0 | 0.394 |
| 16/4/decode/1 | 5727.86/4567.96 | 7292.75/6560.88 | 2300.5 | 143.8 | 0.785 |
| 32/8/encode/0 | 2356.80/1974.29 | 5927.90/5691.51 | 5660.4 | 1415.1 | 0.398 |
| 32/8/decode/1 | 13493.1/13121.6 | 19818.4/19156.7 | 1693.1 | 52.9 | 0.681 |
| 64/16/encode/0 | 6762.15/6244.27 | 13466.2/13034.0 | 4983.5 | 1245.9 | 0.502 |
| 64/16/decode/1 | 35639.7/34592.2 | 46095.2/45566.9 | 1455.9 | 22.7 | 0.773 |
| 128/32/encode/0 | 19089.8/18525.8 | 34326.9/33359.1 | 3910.0 | 977.5 | 0.556 |
| 128/32/decode/1 | 79560.2/78756.6 | 107555/106579 | 1247.9 | 9.7 | 0.740 |
| 128/128/encode/0 | 34499.1/33759.6 | 66088.1/65209.9 | 2030.9 | 2030.9 | 0.522 |
| 128/128/decode/1 | 86477.5/84880.3 | 128470/127872 | 1044.7 | 8.2 | 0.673 |

Across all 104 matched end-to-end cases, median-time encode ratios ranged from
0.247 to 0.868 (geometric mean 0.532), while decode ratios ranged from 0.346 to
0.813 (geometric mean 0.570). Thus ff8xor took 1.88x as long for encode and
1.76x as long for decode on the equal-case geometric mean; the overall 104-case
factor was 1.78x. No measured end-to-end ff8xor case was faster. These are
summaries from one run rather than raw-sample confidence intervals, so they do
not establish a causal improvement versus an older generated-circuit build.

The exact-synthesis/batched-XOR checkpoint was then measured in quick mode
against a frozen executable using checksum `4300215e...`, alternating old/new
process order while each executable retained its internal packed/native ABBA
order. Across nine outer repetitions, packed-drift-normalized median speedup
was 1.017x for the four encode cases and 1.014x for the ten decode cases.
Shorter repeats, including best-time comparisons, moved to either side of
1.0, so this is not claimed as a reproducible codec-level speedup. Five
same-binary microbenchmark repetitions did establish the intended retained
region:

| Batched kernel | Bytes per stream | Median speedup | Best-time speedup |
|---|---:|---:|---:|
| XOR2 | 64 | 1.355x | 1.367x |
| XOR2 | 1,024 | 1.057x | 1.046x |
| XOR4 | 64 | 1.793x | 1.794x |
| XOR4 | 1,024 | 1.117x | 1.096x |

At 4 KiB and 64 KiB the selected implementation was effectively neutral
because it deliberately returns to sequential vector loops after one hoisted
mode lookup. True interleaving and a manual four-vector AVX2 unroll were both
slower there and were rejected. Allocation and correctness checks were outside
timing, modeled traffic counted two loads plus one store, and transposes were
excluded. The current quick old-versus-packed ratios remained well below one,
so these local gains do not change the experiment's overall slower result.

### Deferred-zero and redundant-materialization checkpoint

The transform planner now tracks a small logical state for each work buffer:
logical zero, or an opaque nonzero identity number whose payload is physically
materialized. A nonzero identity number proves symbolic equality, not that the
payload contains a nonzero bit. The planner never compares payload bytes.
Copies preserve an identity; an arbitrary input or a general butterfly
receives a fresh one. This permits
exact reductions for zero/zero, zero/input, equal-input, sentinel, and
coefficient-one butterflies, and lets the formal-derivative accumulation skip
zero sources, cancel equal identities, or replace `zero ^= source` with a
copy. A deferred zero is physically cleared immediately before any general
kernel could read it, and every public output is materialized before return.

Tracking is deliberately selective. Decode uses it at shard sizes of 64 KiB
or larger. Encode uses it at those sizes only when its final `m`-sized input
chunk is partial; full chunks contain no structural padding to eliminate.
Smaller payloads and full-chunk encodes retain the exact prior schedule because
the state-machine overhead did not amortize on the implementation host. This
threshold is a measured policy for this experimental branch, not a universal
CPU constant.

The benchmark reports the static scheduled traffic, a signed estimated byte
elision, the adjusted total, and seven operation counters. These are
deterministic load/store models, not PMU measurements. For `k=7`, `r=3`, and
64 KiB shards, encode deferred one clear and reduced one butterfly, modeling
65,536 bytes elided. A deterministic one-original-loss decode deferred seven
clears and reduced seven butterflies without adding a clear, modeling 458,752
bytes elided. Public entry points reset the thread-local diagnostics even on
validation errors, no-loss, `k=1`, and `r=1` fast paths, so a row cannot inherit
stale counters from an earlier transform.

A focused same-process A-B-B-A comparison against the untracked schedule for
`k=7`, `r=3` measured approximately 8--9% lower encode time and 23--27% lower
one-loss decode time at 64 and 128 KiB. The broader 64 KiB decode matrix was
13.6% faster by geometric mean. Host noise moved some sequential individual
cells in both directions, so retention is based on counterbalanced focused
comparisons, the explicit 64 KiB threshold, and leaving non-benefiting regions
on the old path. The two production objects grew by 6,744 text bytes (0.738%)
and 64 BSS bytes on GCC 13.3.

Packed-equivalence, tracked-path, and native-round-trip tests poison every work
shard they allocate with distinct nonzero data before encode and before each
decode. Forced-locator-shift decoding likewise refills its work with a nonzero
pattern on every iteration. Native round trips include arbitrary native bytes
at `k=7,r=3` and 64 KiB, so the tracked path is exercised without relying on a
packed transpose. Forced-backend coverage includes `k=5,r=5` (`k<m`) and
`k=13,r=3` (four chunks `4/4/4/1`), plus no loss, one, multiple, maximum, and
mixed losses. The poisoned tracked cases specifically prevent a missing
deferred clear from passing merely because a fresh allocation happened to
contain zeros. The Release build and a Debug ASan+UBSan build pass this matrix.

The highest-value follow-up for this checkpoint is a generated fused
copy-plus-`(1+c)` multiplier for the inverse `(x,0)` case. The current reduced
path copies `x` to `y` and then multiplies `x` in place; except when `c=1`, it
still moves the same four payload-buffer equivalents as a general butterfly.

### Retained blocked AVX2 boundary-transpose checkpoint

The first retained packed-compatibility acceleration was a separately compiled,
runtime-gated AVX2 8x8 transpose in both directions. It remains the fallback
and tail path beneath the later AVX-512 kernels. It processes 32 groups
(256 shard bytes) per block and hands every shorter remainder to the portable
word-transpose implementation. The native codec still neither calls nor times
this helper. Source and destination are caller-owned, non-overlapping buffers;
conversion
allocates no memory. `--portable-transpose` selects the portable control for
the transpose-inclusive codec rows, while standalone microbenchmarks always
compare portable and automatic dispatch directly.

The following full-mode results used the same implementation host and GCC 13.3
Release build described above, with OpenMP disabled and pinning disabled.
Each portable/AVX2 microbenchmark pair used ABBA ordering, three warm-ups, and
31 rounds (62 observations per implementation); correctness and allocation
were outside timing. Throughput is one input shard byte per reported byte.

| Direction | Shard bytes | Portable MB/s | AVX2 MB/s | Median speedup |
|---|---:|---:|---:|---:|
| Packed to plane | 1,024 | 3,937 | 21,121 | 5.42x |
| Packed to plane | 4,096 | 3,957 | 21,098 | 5.34x |
| Packed to plane | 65,536 | 4,021 | 21,134 | 5.26x |
| Packed to plane | 1,048,576 | 3,948 | 20,586 | 5.21x |
| Plane to packed | 1,024 | 7,781 | 35,072 | 4.55x |
| Plane to packed | 4,096 | 7,791 | 35,787 | 4.61x |
| Plane to packed | 65,536 | 7,868 | 36,025 | 4.58x |
| Plane to packed | 1,048,576 | 7,749 | 34,003 | 4.39x |

Two complete `--include-transpose` runs then compared automatic AVX2 dispatch
against `--portable-transpose` in the same binary and source revision. These
packed-boundary rows were sequential rather than ABBA; the table summarizes
geometric means across all `(k,r,loss)` cases at each size. Native rows remain
separately available and exclude both directions of transpose.

| Operation | 1 KiB | 4 KiB | 64 KiB | 1 MiB | All cases |
|---|---:|---:|---:|---:|---:|
| Encode (6 cases/size) | 3.11x | 3.25x | 2.19x | 1.64x | 2.46x |
| Decode (20 cases/size) | 1.47x | 1.75x | 1.55x | 1.22x | 1.48x |

The automatic packed-boundary backend remained slower than packed Leopard:
its packed-time/ff8xor-time geometric means were 0.271 for the 24 encode rows
and 0.498 for the 80 decode rows. The transpose gain is therefore real but
does not reverse the experiment's end-to-end result.

GCC lowered the forward AVX2 cyclic payload component to 266 instructions with
64 `vpmovmskb`, 56 `vpsllw`, 16 vector loads/stores, no calls, and no stack
references. The inverse cyclic component used a bit-swap/unpack network: 162
instructions and no calls, but three YMM spill slots (six stack references).
Rewriting its unpack stages
in place to shorten live ranges did not remove those spills, so they are
reported rather than hidden; the retained loop still measured 34--36 GB/s.
Exact-size ASan/UBSan tests cover 1 through 67 groups (including the 31/32/33
block boundary), portable/auto/forced-AVX2 modes, source preservation,
destination canaries, the AVX2 portions of all 4,096 one-hot states in a
512-byte block, dense byte patterns, and deterministic random round trips.

### Retained AVX-512 boundary-transpose checkpoint

Automatic packed-boundary conversion now has two additional, independently
compiled ZMM paths. Packed-to-plane uses `VPSHUFBITQMB` eight times per
64-byte packed block and stores the eight masks into their planes. Its runtime
contract is OS-enabled XMM/YMM/opmask/ZMM state plus AVX-512F, AVX-512BW, and
AVX-512BITALG. Plane-to-packed uses a 512-byte hierarchy: eight plane loads,
an 8x8 bit-swap network, then `VPERMT2B` pair/quad/octet interleaving. Its
contract substitutes AVX-512VBMI for BITALG. Neither path requires AVX2,
AVX-512VL, or the other path's special feature. The normal build keeps both in
isolated translation units and runtime-gates entry; setting
`LEO_FF8XOR_ENABLE_AVX512=OFF` builds their baseline-safe stubs instead.

The exact automatic order is directional. Every legal 64-byte shard can use a
complete BITALG forward block. The VBMI inverse is selected only for each
complete 512-byte block; an available AVX2 256-byte block and then portable
groups handle its remainder. Thus the retained inverse threshold is exactly
512 shard bytes, not a claim that ZMM is best below that size. Native ff8xor
still performs no boundary transpose.

The retained choice followed isolated same-process ABBA prototypes rather
than ISA preference. The 256-bit BITALG forward path beat AVX2 by roughly
1.51--1.66x from 256 bytes through 1 MiB, but the 512-bit form won at every
measured size, so the redundant VL implementation was rejected. Two first
inverse attempts were also rejected: mask extraction plus scalar plane loads
ran at about 0.52--0.53x AVX2 from 256 bytes upward, while a gather form fell
to 0.125--0.189x. The retained VBMI permutation hierarchy instead won from
its first full 512-byte block onward.

The following production numbers used GCC 13.3.0 Release on the AMD Ryzen
Threadripper PRO 9985WX implementation host, pinned to logical CPU 24. OpenMP
was enabled in the build but these helpers are single-threaded. Each row is a
same-binary forced-AVX2/forced-ZMM ABBA pair with three warm-ups and 31 rounds
(62 observations per implementation); allocation and correctness were
outside timing. Throughput counts one input shard byte.

| Direction | Shard bytes | AVX2 MB/s | Retained ZMM MB/s | Median speedup |
|---|---:|---:|---:|---:|
| Packed to plane (BITALG) | 512 | 20,031 | 47,476 | 2.37x |
| Plane to packed (VBMI) | 512 | 31,990 | 58,216 | 1.82x |
| Packed to plane (BITALG) | 1,024 | 20,078 | 48,058 | 2.39x |
| Plane to packed (VBMI) | 1,024 | 33,051 | 73,608 | 2.23x |
| Packed to plane (BITALG) | 4,096 | 20,091 | 55,608 | 2.77x |
| Plane to packed (VBMI) | 4,096 | 33,828 | 77,960 | 2.30x |
| Packed to plane (BITALG) | 65,536 | 20,077 | 49,501 | 2.47x |
| Plane to packed (VBMI) | 65,536 | 32,562 | 74,317 | 2.28x |
| Packed to plane (BITALG) | 1,048,576 | 19,765 | 47,685 | 2.41x |
| Plane to packed (VBMI) | 1,048,576 | 32,772 | 54,234 | 1.65x |

The 1 MiB inverse rate falls relative to its cache-resident peak, but its
AVX2-relative gain remains 1.65x; no tested large-buffer result suggests a ZMM
frequency effect large enough to erase the benefit. Linux denied PMU access
(`perf_event_paranoid=4`), so this run cannot separate frequency, memory, and
cache effects and does not claim a measured clock rate.

Strict archive inspection reports 321 code bytes and 66 instructions for the
BITALG function: eight `vpshufbitqmb`, eight `kmovq`, one payload ZMM load,
no calls, and no explicit stack references or vector spill slots. The VBMI
function is 971 bytes and 158 instructions: GCC folded the bit-swap ternary
expressions into 12 `vpternlogd`, followed by 24 `vpermt2b` and eight ZMM
stores, again with no calls or stack references. A checked CMake target,
`check_ff8xor_transpose_avx512_assembly`, preserves these shape checks. The
baseline archive census also covers the CPUID/XCR0 dispatcher and rejects any
VEX/EVEX leakage there. Each targeted source appends probed exclusions for
non-contract ISA extensions and then its exact required feature set; this
ordering preserves Clang's AVX-512 prerequisite closure. A GCC build with
adversarial global AVX2, AVX-512VL/DQ/VBMI2/VNNI, GFNI, VAES, VPCLMUL, BMI2,
and POPCNT flags passed both the baseline census and the direction-specific
instruction allowlists, and the same sources passed Clang 18 ASan+UBSan.

Correctness now covers groups 1 through 67 across portable, Auto, forced AVX2,
forced BITALG, and forced VBMI modes; this crosses 8-group BITALG, 32-group
AVX2, and 64-group VBMI boundaries with exact-size ASan canaries. All 4,096
one-hot input bits in a full VBMI block and all dense byte values round-trip.
Pure feature tests exhaust all 256 relevant/irrelevant CPUID combinations, all
32 required XCR0-state combinations, and basic-leaf availability.

Representative 1 MiB microbenchmarks from the earlier full reference run were:

| Operation | Coefficient/circuit | Median/best us | Median input/output MB/s |
|---|---|---:|---:|
| Packed multiply-add | log 1; corresponding circuit 20 gates | 15.163/14.904 | 138304/69152 |
| ff8xor multiply | log 1; 20 gates, depth 7 | 18.727/18.341 | 55992/55992 |
| ff8xor FFT butterfly | skew 1; 40 gates, depth 10 | 97.402/97.032 | 21531/21531 |
| ff8xor IFFT butterfly | skew 1; 40 gates, depth 10 | 98.297/97.730 | 21335/21335 |
| Packed to plane transpose | standalone | 258.939/258.200 | 4050/4050 |
| Plane to packed transpose | standalone | 43.218/42.729 | 24263/24263 |

The packed microkernel is multiply-add while the circuit kernel is
multiply-only, so their input-byte rates are not directly equivalent. The
sampled low/average/high circuits also show that performance is not monotonic
in source gate count, which is consistent with memory behavior and instruction
scheduling effects. The two transpose rows in this historical table predate
the retained AVX2 and AVX-512 boundary kernels; the directional AVX2 and ZMM
tables above supersede them for current transpose performance.

## Representative optimized code inspection

`objdump -drC` inspection covers all 256 coefficients in every built base,
AVX2, AVX-512VL, and AVX-512-ZMM family. GCC 13.3 Release with
`-fno-tree-reassoc` passed the strict gate with zero vector-stack references in
all specialized functions. Representative coefficient/skew 42 has 20 source
gates at depth 10 for multiply and 43 source gates at depth 11 for each
butterfly. Its AVX2 loops compiled to 16 `vpxor` for multiply, 41 for FFT, and
40 for IFFT, with the expected plane loads/stores and scalar loop control.
There were no calls inside vector loops, `pshufb`, integer multiply, `pand`,
gathers, payload-table accesses, or dynamic gate indexing. The separate batch
census reported AVX2 XOR2/XOR4 kernels of 116/188 bytes with 2/4 loop XORs;
AVX-512VL used 99/159 bytes and AVX-512 ZMM used 99/175 bytes. All six had zero
calls, vector-stack references, static-table reads, and forbidden operations.

The wider greedy candidate was rejected after the same gate found spills in
380 GCC butterfly functions. With the retained circuits, Clang 18 passes all
strict structural checks but still spills in 739 base/AVX2 butterfly
specializations; multiplier and AVX-512 families remain spill-free. This is a
reported compiler/code-generation limitation, not concealed as a pass of the
optional `--fail-on-spills` gate. GCC AVX-512 output also combines some XORs
into `vpternlog[dq] 0x96`; the separately generated explicit-XOR3 schedules
were measured and rejected for production as documented above.

`tools/inspect_ff8xor_four_buffer_production.py` complements that coefficient
census by extracting the actual `LeopardFF8XorAVX512Four.cpp.o` archive member
and checking all 256 radix-4 functions. Its strict contract requires exact
generated XOR2/XOR3 counts, exactly 32 loads, and at most 32 stores, while
rejecting payload vector spills, calls, non-move vector memory operands,
non-`0x96` ternary immediates, excess reloads, RIP-relative tables, and field
shuffle/GFNI/CLMUL instructions. Both current compilers emit all 32 stores.
The production-object sizes, residency-barrier fix, and current GCC code ranges
are recorded in the four-buffer section above.

## Known disadvantages and next experiment

- Applications must retain the plane layout to avoid boundary transpose cost.
- Each two-way whole-buffer operation incurs one coefficient function-pointer
  dispatch; each successful radix-4 unit incurs one tuple-specialized
  whole-buffer dispatch.
- The fused path requires AVX-512 ZMM and a complete 64-byte-chunked plane.
  AVX2, partial/truncated radix-4 units, tails, and pruned decoder layers retain
  the extra two-way intermediate memory pass.
- The integrated four-buffer path is opt-in and default-disabled. Its current
  checkpoint improves the broad 1--4 KiB matrix but regresses large encode and
  does not provide a uniform large-decode win, so enabling it unconditionally
  would make the experiment worse on the implementation host.
- Encoder chunk accumulation and the decoder formal-derivative sweep expose
  two/four-stream XOR instruction-level parallelism only through 1 KiB. Larger
  buffers use the measured faster sequential-loop shape and still lack fusion
  into adjacent transform passes.
- Sixteen live AVX2 vectors leave no spare vector register, so compiler spills
  remain possible for some generated circuits.
- Static specialization of all multiplier, two-way butterfly, and 256 radix-4
  functions increases source, object, and instruction-cache footprint. The
  current GCC Release/OpenMP checkpoint is about 3.70 MiB for the static archive
  and 2.95 MiB for the experimental benchmark executable, versus about
  93.4 KiB for the packed benchmark. The radix-4 object alone is 417,744 bytes
  with GCC and 559,288 bytes with Clang. Packed-only static clients do not pull
  the generated circuit objects.
- Packed compatibility requires explicit 8x8 transposes.
- The inherited packed payload objects still use the project's global
  `-march=native`; this experiment does not make the existing packed codec a
  cross-host binary. The ff8xor public gateways, initialization call graph,
  dispatcher, portable kernels, and SSE2 tails are separately constrained and
  inspected for x86-64-v1, while its AVX2/AVX-512 entry points are feature
  gated. Applications needing a portable packed backend still require a
  broader packed-code multiversioning change.
- This backend is FF8-only and CPU-only.

The highest-value next optimization for this path is a separately validated
small-payload policy and schedule tuned for the region at or below 4 KiB where
fusion currently wins. That experiment must repeat correctness, production
assembly inspection, and counterbalanced full-matrix timing before changing
the default; this checkpoint deliberately does not enable such a policy.
Other follow-ups include fusing the formal derivative with its adjacent FFT
layer, CUDA plane kernels, persistent plane pipelines, packed-boundary
transpose/scaling fusion, and synthesis optimized for depth and register
pressure rather than gate count alone.
