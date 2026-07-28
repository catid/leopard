# Open correctness and optimization findings, 2026-07-24

These are the findings from a source audit of the production codec that survived
an adversarial verification pass.  Each was re-derived from the source by a
second reader that was instructed to refute it; seven of the fourteen verified
candidates were refuted and are not listed.  Magnitudes below are the corrected
ones, not the original claims.

None of these is fixed.  They should become Beads issues; the issue tracker was
not writable during this session (see the note at the end).

## 1. Dead block-0 pruned input schedule on every native-high decode plan — RESOLVED

**Shipped 2026-07-28.**  The high-profile planner loop now starts at block 1;
measured decode plan setup 1.06x-1.15x faster across eight native-high cells
with unchanged execution and digests.  The executor discard guards are
conditional (`.block == 0`), so absent entries flow through unchanged.  See
`experiments/leopard2/optimization_log/22-dead-block0-schedule.md`.

Historical statement:

`PreparePlanExecutionMetadata` compiles a pruned inverse schedule for block 0
(`leopard2.cpp:1264-1286`) on essentially every native-high transform decode
plan.  All four Algorithm 5 executors then discard that entry with a bare
`++input_plan_index` before starting at block 1 (`LeopardFF8.cpp:4409-4416`,
`:4665-4672`; `LeopardFF16.cpp:4252-4259`, `:4526-4533`).

Block 0's inverse is elided outright by the
`FFT_0(IFFT_0(F_0) + H_later) = F_0 + FFT_0(H_later)` cancellation
(`docs/leopard2_math_and_sources.md:337-353`), so the compiled schedule has no
consumer at all.  Root cause is a missed planner update in commit `bd13ef7`
"perf: cancel Algorithm 5 zero-shift transforms", which moved the executor loops
to block 1 and added the discard guard without touching `leopard2.cpp`.

Cost: one `Theta((T/2) log2 T)` sparsity-independent plan compile (2,304 raw
butterflies at `T=512`) plus up to `(T/2) log2 T * 12` bytes of retained
operation storage held for the plan's lifetime (~27 KB at `T=512`).  It is a
first-use setup and plan-memory reduction, not a byte-proportional win, so it
matters most at small shard sizes.

Fix: start the high branch at block 1.  The low branch at
`leopard2.cpp:1185-1207` must stay, because `LeopardFF8.cpp:3834-3846` and
`:4041-4050` genuinely execute block 0.  The executors' `block == 0` guard
becomes dead but harmless and the `input_plan_index == input_plan_count`
debug assertions still hold.

## 2. Batch preflight sort is 2x-16x slower than it needs to be — PARTIALLY RESOLVED, SHARES CORRECTED

**2026-07-28:** the writable-partitioned scheme landed
(`experiments/leopard2/optimization_log/23-writable-partitioned-preflight.md`)
on top of the earlier heapsort->introsort step, and the end-to-end batch entry
point moves only **1.01x-1.06x** — the component-level shares below
substantially overstate the sort's weight in the production path.  The LSD
radix option is declined: not worth its scratch-size ABI change at the true
share.  Finding 3's Amdahl ceilings inherit the same correction.

Historical statement:

`SortBatchRangesWithoutAllocation` / `SiftBatchRangeHeap`
(`leopard2.cpp:5123-5160`, called at `:6037` and `:6133`) is a hand-rolled
heapsort.  Its documented justification — "without relying on an allocator",
`docs/leopard2_batch_preflight.md:32` — does not hold: `<algorithm>` is already
included at `leopard2.cpp:39` and `PrepareEncodeBatchRanges` already calls
`std::sort` unconditionally at `:5480`/`:5489` inside the path the same document
calls allocation-free.

Measured on this host over the exact record set the batch benchmark produces for
`K=100 R=28 B=1024` (134,145 records, 3.2 MB): heapsort 8.66 ms, `std::sort`
4.35 ms (2.0x).  A decomposition of the scalable preflight gives flatten
0.057 us/item (0.6%), heapsort 8.75 us/item (98.6%), sweep 0.069 us/item (0.8%).
Two better options were measured: a 3x16-bit LSD radix on `range.begin` into a
second buffer is 0.52 us/item (16.6x), needing an additive doubling of
`leo2_encode_batch_preflight_scratch_size` and a check that the sweep is
insensitive to tie order at equal `begin`; alternatively, since
`ValidateSortedBatchRanges` permits readable/readable overlap, only the
`B*(R+1)` writable records need ordering, and sorting just those plus binary
search measured 1.39 ms (6.2x).

## 3. Batch preflight is a serial section inside the thread pool

`leo2_thread_pool::Run` executes `preflight` on the calling thread while holding
`run_mutex_` (`leopard2.cpp:485-491`), before `EnsureStarted` and before the
generation bump that wakes workers, while every worker blocks at
`leopard2.cpp:607`.  All four batch entry points route through it (`:7595`,
`:7657`, `:8390`, `:8452`).

Measured (`K=100 R=28 B=1024`, GF8 legacy-high, pinned, median of 7), the
preflight is a byte-independent 8.5-9.5 us/item serial section.  With the
AUTO/AVX2 backend it is 91.3% of the batch at 64-byte shards — an Amdahl ceiling
of 1.09x at any thread count — 59.2% at 1 KiB (1.69x) and still 13.9% at 8 KiB.

Finding 2 is the cheap fix for most of this: the serial section is almost
entirely the sort.  Parallelising the flatten pass would attack 0.6% of it.

## 4. One-loss direct repair is structurally unreachable at scale

`CanPrepareDirectRepair` (`leopard2.cpp:1491-1502`) caps `K` at
`kDirectLegacyMaxRepairOriginals = 16` and parent dimension at 256 outside the
narrow measured equal-rounded window.  Consequently GF8 high `K=240 R=16`
(parent dimension 240), low GF8 `K=32 R=192`, and low GF16 `K=128 R=896` are
blocked *only* by the `K<=16` cap.

`tools/leopard2_operation_counts.py` at one loss gives shard-multiply ratios of
2.37x (GF8 240/16), 16.45x (low GF8 32/192), 8.56x (low GF16 128/896), 3.23x
(GF16 1000/200) and 2.90x (GF16 4096/512) in favour of direct repair, and the
butterfly count is loss-independent, so a single lost original pays full
Algorithm 4/5 cost today.

Narrow this to **exactly one loss**: the retained directional screen
(`experiments/leopard2/direct_repair/results/one_loss_equal_rounded_exact_main.json`)
shows `L=8` won only 245/390 cells with a worst cell of 0.133x, which is why the
`K=65` exception was kept separate.  The three `K<=256` cells need only a
selector relaxation plus a pinned exact-main gate; `K=1000`/`4096` additionally
need heap-backed generator rows and a coset-decomposed weight product, because
`kDirectMaxParentDimension = 256` is a stack-array bound
(`leopard2.cpp:2142`, `:4811-4815`).

## 5. Sparse encode output pruning has a hard ceiling; the real win is a direct systematic row

Legacy-high encode runs its `ceil(K/T)` inverse `T`-transforms unconditionally
(`LeopardFF16.cpp:2278-2333`, `LeopardFF8.cpp:2500-2550`) and lets the requested
parity mask affect only the final forward transform.  One evaluation point
depends on all `T` coefficients, so exact output pruning can never reach the
accumulation, which is `ceil(K/T)/(ceil(K/T)+1)` of the graph: 81.5% for GF16
`K=1000 R=200` and 75% for GF8 `K=192 R=64`.  The ceiling on exact-mask sparse
pruning for those cells is therefore 18.5% and 25%, of which prefix truncation
already captures 11.7% for a parity-0 mask.

The unexploited item is a direct systematic generator row: `K` shard
multiply-adds versus 4,396 butterflies plus 1,000 copies plus 768 XOR shard-ops
at GF16 `K=1000 R=200`, with roughly 10x less modelled memory traffic.  That
path exists and is validated for the high profile —
`experiments/leopard2/direct_encode/results/checkpoint.json` region
`excluded_high_profile` records six pinned high-profile cells with median gain
+541.3% and one regression — but is unreachable at scale because
`leopard2.cpp:337-338` and `:1428-1452` cap the table at `K,R <= 16` and
`CanAutoDirectEncodeCodec` restricts AUTO to `LEO2_PROFILE_LOW_V1`.

**2026-07-28 correction (report 20):** the bounded K,R <= 16 re-screen is
invalid.  Its extractor reported codec-setup medians as encode-execution
medians, omitted the historical 17.88x cell, lacked frozen-binary ABBA
provenance, and made a vacuous outer digest comparison.  The reported
0.875x-1.017x “ceiling” therefore neither confirms nor refutes the small-table
opportunity.  The bounded and at-scale variants both remain unmeasured on the
current tree; `leopard-79h.42.1` tracks a corrected explicit-AVX2 execution
campaign.

## 6. GF16 decode keeps a full-prefix workspace — RESOLVED

**Shipped.**  `AVX2DecodeExecutionTileBytes` now carries a GF16 legacy-high
entry for `padded_side` 256 and 512 at a 16 KiB tile, floored at a 64 KiB
aligned prefix.  A 13-cell sweep over 4/8/16/32 KiB measured **1.15x to 1.91x**
decode with matching digests, the largest single win of the AVX2 campaign;
scratch falls 67.7 MB -> 17 MB at `K=2000 R=500` and 33.6 MB -> 8 MB at side
256.  Evidence and the rejected per-regime table are in
`experiments/leopard2/optimization_log/13-gf16-decode-tiling.md`.

**Both follow-ups below are now resolved by the live-set rule** documented in
`experiments/leopard2/optimization_log/15-liveset-tiling-rule.md`. The
`padded_side` key was a proxy for `2 * padded_side * aligned_prefix` against
L3, and it only held at the 64 KiB shards it was swept at. Keying on the live
set directly — tile to a 16 MB target once the untiled set passes 64 MB, while
retaining the calibrated 64 KiB entries — extended the win to
`padded_side` 64 and 128 and to larger shards, reaching **3.226x** decode and
**2.735x** encode. The same treatment was applied to the GF16 encode path,
where these cells had been at parity with leopard1 and are now 2.44x-2.54x.

Historical statement of the two follow-ups:

- The sweep found two regimes — single message-block counts
  (`original_count <= padded_side`) peak at 4-8 KiB, multi-block counts peak at
  32 KiB for side 256 and 8-16 KiB for side 512.  The shipped uniform 16 KiB is
  within ~5% of the best measured tile everywhere but is nobody's optimum.
  Promoting the per-regime table needs the same offline
  `padded_side` x shard size x block count sweep that produced the GF8 table.
- Sides other than 256/512 and shard sizes below 64 KiB remain untiled for GF16
  and are unmeasured, not rejected.

## 7. Pruned-transform compilation materializes the DAG three times

`CompilePrunedTransformPlan` materializes the complete unpruned DAG three times
(`Leopard2Plan.cpp:1174-1181` at 12 B/op, `:1195` at 2 B/op, `:1231` at 12 B/op)
and makes about 11 heap allocations per plan, even though
`CompileSparseNodeReverse` (`:529-581`) proves the operands are recomputable
from the deterministic traversal with zero heap.

The regime is small-side high plans, where per-plan allocator traffic dominates:
`leopard2.cpp:1053-1062` records that omitting these schedules "reduced plan
setup by 3.0x at T=8 and 5.8--7.9x at T=64", and the shipped remedy is a
four-condition policy skip covering only GF8/AVX2/high/side in {8,64}.  GF16
high `K=60000 R=1000` still pays roughly 1,400 allocations of pure setup.

The minimal rewrite is not "emit inside the backward recursion" — `operations`
must stay in forward order for `MatchFusedFour` and
`ExecutePrunedTransformPlan`, and the retained-count pass at `:1313-1318` is a
deliberate persistent-footprint decision.  Three recursions over the existing
traversal (forward liveness into 2 packed bits/op, backward need/write into
6 packed bits/op while counting, then forward emission) replace 26 B/op of
transient state (6.4 MB at side 32768) with 1 B/op (~245 KB) and cut allocations
from ~11 to ~4.

## 8. GF8 range butterflies are never known-answer-tested with `prefer_fused == false`

`Ops::ff8_ifft_butterfly4_range` and `Ops::ff8_fft_butterfly4_range` take a
`prefer_fused` flag that selects materially different schedules.
`TestFF8ButterflyRanges` passes only `true` (`Leopard2Backend.cpp:794`, `:813`).
The GF16 companion does loop over both (`Leopard2Backend.cpp:1677`, `:1712`).
The non-fused GF8 range schedules are therefore not covered by the startup
qualification known-answer test.  This is a test-coverage gap, not a known
defect; it was reported by the survey and is included here because it is cheap
to close.

## Measured negative results, 2026-07-24 (do not re-propose without new evidence)

### N1. Cache-blocking the GF16 AVX2 radix-four schedules does nothing

Two independent restructurings were implemented, verified wire-identical over
150-200 randomized shapes each, and measured against `HEAD` with clustered ABBA
on one pinned CPU over twelve GF8/GF16 cells:

1. Byte-blocking `AVX2FF16Butterfly4Split` into 8 KiB chunks so the four-shard
   group stays inside L1 across all four radix-two layers.  Encode ratios
   0.980-1.015, i.e. neutral to slightly negative.
2. Chunking `AVX2FF16Butterfly4Range`'s coordinate range so all four layers run
   over a group of coordinates before advancing, targeting a 256 KiB live set,
   instead of four full sweeps.  Encode ratios 0.985-1.005.

The first result also revealed why the survey's mechanism was misleading:
`UseFusedButterfly4` (`LeopardFF16.cpp:114-130`) returns false above 128 bytes,
so `AVX2FF16Butterfly4Split` is not even on the hot path for 4 KiB or 64 KiB
shards; the range sweeps are.  The second candidate targeted the real path and
was still neutral.

Conclusion: at these shard sizes the GF16 AVX2 transform is ALU/port bound, not
L1/L2-miss bound.  The hardware prefetchers and the out-of-order window already
absorb the layer-major reuse distance.  Op-count reduction is what moves this
kernel — which is exactly the shape of the GFNI result — and cache restructuring
is not.  Do not re-propose blocking here without a profile showing stalls.

### N1b. Lifting the GF8 `kFusedByteLimit` under GFNI does not pay

`AVX2FF8FFTButterfly4` and `AVX2FF8FFTButterfly4Range` fall back to streaming
radix-two sweeps above 1 KiB per shard (`Leopard2BackendAVX2.cpp`, constant
`kFusedByteLimit = 1024`).  A port model argued the constant is stale under
GFNI: with nibble tables a radix-two butterfly costs eight vector-ALU ops per
32-byte tile so the streaming and fused forms are both ALU bound at the same
rate, while an affine butterfly costs three, which should leave the streaming
form dispatch bound and predict about 1.50x on the GF8 forward transform above
1 KiB from changing one constant.

Implemented, verified wire-identical over 180 randomized shapes, and measured
against the same build with the constant unchanged, ten GF8 cells, three rounds:
**encode 0.975-1.048 (median about 1.00, with one 2.5% regression at
K=192 R=64), decode 1.001-1.082 (median about 1.05)**.  The predicted 1.50x did
not appear at all.  Reverted: a decode-only ~5% that overturns an explicitly
measured streaming policy is not worth the risk, and the model that motivated
it was wrong by more than an order of magnitude.

This is the third static model in this session to predict a large win that
measured neutral, after the two cache-blocking restructurings in N1.  Treat
port and traffic models for these kernels as hypothesis generators only; on
this part they have not once predicted a real effect.  The one change in this
area that did pay — GF16 radix-four fusion — was validated by a standalone
kernel screen before being integrated.

### N1c. Radix-eight pays only on the cold streaming round, not on in-place rounds

After the out-of-place radix-eight staging landed (1.332x at K=224 R=32, 64 KiB),
the obvious extension was an in-place radix-eight range, which by round counting
should unlock the wider window: covering log2(m) layers with rounds of three,
two and one layer gives 3 rounds today versus 2 with radix-eight at m = 32 and
m = 64, and 4 versus 3 at m = 128 and m = 256 — a modelled 1.50x and 1.33x.

Touch counts support the model: a radix-four round is 1.00 touches per
shard-layer, a radix-two round is 2.00, and a radix-eight round is 0.67.

It does not pay.  An isolated in-place screen measured 1.119x at 1 KiB and
0.98x-1.01x at 4 KiB, 16 KiB and 64 KiB — neutral.  The real codec agrees: with
staging radix-eight alone, m = 64 (K=192 R=64) measured 1.019x and m = 128
(K=128 R=128) 1.05x-1.07x, against a model predicting 1.50x and 1.33x.

The reason is which data is cold.  The out-of-place staging round reads the
caller's source shards, 224 of them, 14.7 MB at 64 KiB, streamed once and never
resident.  Every in-place round works on the shared work array, 32 shards and
2.1 MB at 64 KiB, which stays in L2 or L3 across rounds.  Reducing touches only
buys time when the touches actually go to memory, so the same restructuring that
is worth 1.26x on the streaming round is worth nothing on the resident ones.

Two cautions for whoever picks this up.  First, the isolated in-place screen is
*warmer* than the real encoder, which interleaves work-array rounds with cold
source reads, so treat its exact numbers as a lower bound rather than a
measurement of the production access pattern — the real-codec m = 64 and m = 128
figures above are the trustworthy ones and they agree.  Second, the round-count
model is not predictive here and should not be used to justify further radix
changes on its own; only the cold round matters.

The m == 16 case cannot be rescued this way either.  At m == 16 a radix-eight
staging round leaves exactly one layer, and that layer is the accumulating fold
into xor_result.  An accumulating radix-eight cannot absorb it, because at
m == 16 the final layer pairs coordinates i and i+8, which crosses the two
radix-eight groups rather than living inside one.

### N1d. GF8 encode byte tiling at 64 KiB: a sharp 4-5%, not worth shipping on this evidence

The GF8 encode tiling policy (`leopard2.cpp`, the `padded_side` switch setting
`high_gf8_tile_bytes` and `high_gf8_min_shard_bytes`) does not engage at 64 KiB
shards: for `padded_side == 32` it requires at least 512 KiB.  The work array is
then 32 x 64 KiB = 2.1 MB, which exceeds the 1 MB L2, so every in-place round
goes to L3.  Since that policy was "generated offline from pinned exact-backend
measurements" taken with nibble kernels, and GFNI plus radix-eight have since
made the arithmetic much cheaper, the calibration was a good candidate for being
stale — the same reasoning that correctly predicted the GF16 radix-four cap was
stale.

It is only slightly stale.  Lowering the threshold and sweeping the tile size,
against the current GFNI + radix-eight build, all digests matching:

| Cell | untiled | 8 KiB tile | 16 KiB tile | 32 KiB tile |
| --- | ---: | ---: | ---: | ---: |
| K=224 R=32, 64 KiB | 446.4 us | **1.043x** | 0.900x | 0.970x |
| K=224 R=32, 16 KiB | 114.3 us | **1.040x** | 0.989x | 0.990x |
| K=200 R=32, 64 KiB | 426.3 us | **1.053x** | 0.919x | 0.838x |
| K=160 R=32, 64 KiB | 339.1 us | **1.044x** | 0.910x | 0.840x |

An 8 KiB tile is consistently worth 4-5%, but 16 KiB and 32 KiB lose 10-16%.
That is a sharp optimum measured on four cells, against a policy the project
derived from a pinned offline sweep.  Shipping it would be exactly the kind of
overfit the existing side-specific table exists to prevent, so it is recorded
here rather than landed.  If it is picked up, it needs the same offline sweep
across `padded_side`, shard size and block count that produced the current
table, not four cells.

### N2. Removing the third encode-validation pass is not measurable

`ValidateEncodeBuffers` walks `recovery[]`, then `original[]`, then calls
`ValidateMonotonicEncodePointers`, which re-walks both.  A fast path was added
that proves output pairwise disjointness during the first pass and input/output
disjointness from the two address spans, skipping the third walk entirely for
the ordinary packed-slab layout.  Measured against `HEAD` at 64/512/4096 bytes
for K/R in {2/2, 8/8, 16/16, 32/32, 128/128, 240/16}: every ratio moved by less
than the run-to-run spread.

The third pass is pointer arithmetic only; the cost is in the first two passes,
which do the mandatory `MakeRange` overflow checks, scratch-overlap checks and
systematic-pad checks.  Those are the API contract, not overhead.

What the same experiment did establish, from an attribution build that returns
success from `ValidateEncodeBuffers` immediately: encode buffer validation costs
1.1-1.8 ns per shard pointer, and removing it moves GF8 K=32 R=32 at 64 bytes
from 0.67 to 0.94 against exact main.  Separately, a size sweep against exact
main shows Leopard2's byte-proportional encode cost is already equal or better
in every shape measured (K=2 R=2: 0.0137 vs 0.0282 us/KiB; K=32 R=32: 1.89 vs
2.23 us/KiB).  The entire remaining small-shard deficit is fixed per-call cost,
and validation is its largest single remediable component.  Any future attempt
should therefore attack the per-element work in the first two passes, or hoist
validation across calls, not the third walk.

## The AVX2 ceiling, measured

A stubbed-kernel attribution build gives the fraction of an encode call that the
vector kernels own, and therefore the hard Amdahl bound on any kernel work.  The
build makes the chosen dispatch-table operations return immediately, bypasses
the startup known-answer test and the benchmark round-trip check under one
measurement-only macro, and is never a shipping configuration.

**Stubbing every operation in the table** gives the true floor:

| Cell | share in the ops | ceiling with infinitely fast ops |
| --- | ---: | ---: |
| GF8 K=240 R=16, 4 KiB | 98.4% | 63x |
| GF8 K=240 R=16, 64 KiB | 99.9% | 1087x |
| GF8 K=224 R=32, 64 KiB | 99.9% | 1281x |
| GF8 K=192 R=64, 64 KiB | 99.9% | 1222x |
| GF8 K=128 R=128, 64 KiB | 99.9% | 801x |
| GF16 K=1000 R=200, 4 KiB | 97.8% | 46x |
| GF16 K=1000 R=200, 64 KiB | 66.5% | 2.98x |
| GF16 K=255 R=129, 64 KiB | 79.5% | 4.87x |

So for GF8 there is **no meaningful Amdahl limit**: essentially the whole encode
call is vector kernel time.  Only the large GF16 cells have real non-kernel
cost — 20% to 33%, which is the scratch and staging traffic of a 33 MB working
set — and even there the ceiling is 3x to 5x.

An earlier revision of this section reported ceilings of 1.40x to 3.68x from a
*partial* stub set and concluded that high-rate GF8 was mostly non-kernel.  That
was wrong: the stub set omitted the out-of-place, xor-output, weighted and
multiply-add-outputs kernels, which are exactly the ones the high-rate encoder
uses most.  The corrected numbers above supersede it.  Pure XOR operations
(`xor_mem`, `xor_mem_2to1`, `xor_memory_sources`, `xor_memory_dense`) were
measured separately and are 0% to 4.6% of encode, so XOR bandwidth is not a
factor either.

### Per-operation profile, GF8 K=224 R=32, 64 KiB

Stubbing one operation at a time against a 616 us full call:

| Operation | owns | share |
| --- | ---: | ---: |
| `AVX2FF8Butterfly4Out` | 246 us | **40%** |
| `AVX2FF8Butterfly4RangePreparedImpl` | 145 us | 24% |
| `AVX2FF8WeightedIFFTButterfly4` | 4 us | 0.6% |
| `AVX2FF8IFFTButterfly4XorKernel` | ~0 us | ~0% |

`AVX2FF8Butterfly4Out` is the out-of-place radix-four that reads caller source
shards straight into the work slots, and it is the single dominant kernel of
legacy-high GF8 encoding.  It is already fused radix-four and already picks up
the affine multiply through `AVX2FF8ProductVector`.

Its combined read plus write throughput is **flat at 116-188 GB/s from 1 KiB to
64 KiB per shard**.  Flatness across a range that spans cache-resident to
far-beyond-cache rules out a bandwidth limit, and the rate is well under the
port and ALU bounds for its instruction mix, so what limits it is not yet
identified.  A two-tile, single-offset unrolling of it — on the hypothesis that
eight pointer increments per 32-byte tile made it dispatch bound — was
implemented, verified wire-identical, and measured at about 1.01x: rejected.

### Why it is that fast and no faster: it is at the memory-movement floor

`experiments/leopard2/gfni_codec/butterfly4out_floor_screen.cpp` reproduces the
kernel's exact access shape — 7 message blocks of 8 radix-four groups, 224
source shards read out-of-place into 32 work shards — and compares it against
two floors: `memfloor`, which does the four loads and four stores and **no
arithmetic whatsoever**, and `xorfloor`, which adds only the eight XORs the
algorithm requires.  Combined read plus write GB/s, one pinned CPU, best of 60:

| bytes/shard | memfloor (no arithmetic) | xorfloor | production | production / memfloor |
| ---: | ---: | ---: | ---: | ---: |
| 1024 | 211.4 | 188.0 | 153.4 | 1.38x |
| 4096 | 185.7 | 191.1 | 157.6 | 1.18x |
| 16384 | 170.3 | 168.4 | 156.7 | 1.09x |
| 65536 | 135.8 | 136.0 | 130.4 | **1.04x** |

At 64 KiB per shard the production kernel is within **4%** of a kernel that only
moves the bytes.  Deleting every affine multiply and every XOR would buy 4%.
The multiplies are already free, the XORs are already free, and the kernel is
doing nothing but load/store.  The earlier reading of "2.5x under its bound" was
wrong: it compared against an ALU and port bound, but the binding constraint is
load/store throughput, and the kernel is essentially at it.

Two consequences.  First, **a 2x AVX2 target is unreachable for high-rate GF8**,
and not for an Amdahl reason — it is because the operation that owns 40% of the
call cannot go more than 4% faster at production shard sizes no matter what
instructions it uses.  Wider vectors do not help a load/store floor either, so
this conclusion is not lifted by the deferred AVX-512 work.  Second, the only
lever left against a memory floor is **moving fewer bytes**: a radix-eight
out-of-place pass would carry three transform layers per load/store round
instead of two, cutting traffic by a third on the 64% of the call that these two
streaming kernels own, for a modelled ~1.27x.  Note that the radix-eight
proposal was rejected earlier in this session on the grounds that it does not
change vector-ALU operations per byte-layer — which is correct and now
irrelevant, because this kernel is not ALU bound.  The real obstacle is
register pressure: eight live inputs plus eight live outputs exceeds sixteen
`ymm` before any matrix or temporary.

Non-temporal stores were also measured and are 2x to 3x *worse* (39-61 GB/s):
the work slots are re-read by the following transform layers, so bypassing the
cache is the wrong trade here.

### A note on method

Five separate static models in this session predicted large wins that measured
flat: two cache-blocking restructurings, a validation-pass elimination, a port
model predicting 1.50x from one constant, and this dispatch-bound unrolling.
The only two changes that paid — GFNI affine multiplication and GF16
radix-four fusion — were both validated by a standalone kernel screen before
integration.  For these kernels, measure first.

## A shard-touch model that predicts these results

Every transform stage costs exactly `2m` shard-touches regardless of radix: a
radix-four stage is `m/4` groups of 8 touches, radix-eight is `m/8` groups of 16,
radix-two is `m/2` pairs of 4.  So a transform of `L = log2(m)` layers costs
`2m * ceil(L/2)` with radix-four rounds and `2m * ceil(L/3)` with radix-eight.

For legacy-high GF8 encode at `K=224 R=32` (`m=32`, `L=5`, 7 message blocks):

| | inverse | accumulate | forward | total |
| --- | ---: | ---: | ---: | ---: |
| radix-four staging | 1344 | 192 | 192 | 1728 |
| radix-eight staging | 896 | 192 | 192 | 1280 |

Predicted 1.350x against a measured 1.332x — **1.4% error**.  The model is worth
using for sizing candidates.

It is not sufficient on its own, because per-touch cost varies by roughly 2x with
temperature.  Dividing the stubbed-operation profile: `AVX2FF8Butterfly4Out` is
246 us for 448 touches, 0.55 us/touch, reading 14 MB of caller shards that are
never resident; `AVX2FF8Butterfly4RangePreparedImpl` is 145 us for 576 touches,
0.25 us/touch, on the 2.1 MB work array it just wrote.  That ratio is the whole
explanation for N1c: saving touches on the cold staging round pays, saving the
same touches on warm in-place rounds does not.

Use the model to bound a candidate and the temperature split to discount it.

### Where the untouched forward surface sits

`FFT_DIT` (`LeopardFF8.cpp:2233-2280`) is strictly in-place and both ends are
pinned — its input is the work array the last inverse block wrote, and for
legacy-high its output slots *are* the caller's recovery buffers
(`leopard2.cpp:7416-7424`), so there is no mandatory copy to fuse a radix-eight
into, unlike the inverse staging case.  The mirror shape does exist one function
over, in `FFT_DIT_FromCoefficients` (`LeopardFF8.cpp:2293-2417`), whose first and
last stages are already out-of-place through `ff8_fft_butterfly4_out`; that
serves the low-rate encoder and high-rate decode.

The forward is 15.0% of touches at `K=224 R=32`, rising to 40% at `K=64 R=32`,
60% at `K=32 R=32`, and 80-90% of a `LOW_V1` encode.  So a forward radix-eight is
worth about 1.05x at the profiled high-rate shape — under the project's 5% gate,
and less than that in time because the forward's touches are the warmest in the
call — but it is the largest remaining surface for low-rate and small-`K/m`
codecs, which this session never profiled.

Two constraints if it is picked up.  `m == 16` must be excluded for a *different*
reason than the inverse regression: at `L=4`, `ceil(4/2) == ceil(4/3) == 2`, so
radix-eight saves exactly zero touches while still displacing a fused
radix-four.  The `m == 8` correctness hazard does *not* transfer, because
`FFT_DIT` has no `xor_result` parameter at all.  And it must be gated on the
encode-only `prefer_fused` flag so decode's `FFT_DIT` call sites stay
bit-for-bit unaffected until separately measured.

## Deferred AVX-512 work

Parked on 2026-07-24 in favour of plain-AVX2 work, because the AVX2 backend is
what every x86-64-v3 host runs and it is the AUTO default.  Cross-referenced
from `docs/leopard2_avx512_codec.md` ("DEFERRED AVX-512 WORK — start here") so
it is findable from the AVX-512 side too.

### D1. Widen the AVX-512VL member's data path to 512 bits

`Leopard2BackendAVX512.cpp` recompiles the AVX2 kernels with
`-mprefer-vector-width=256`, so it uses AVX-512 only for `ymm16..ymm31`.  The
kernel screen in `experiments/leopard2/gfni_codec/results/kernel_screen.txt`
measures the width change separately from the instruction change, on one pinned
CPU of an AMD Ryzen 9 9950X3D (microseconds for a 32-shard, 5-layer sweep):

| GF8 butterfly, bytes/shard | avx512-256 (shipped) | avx512-512 | width gain |
| ---: | ---: | ---: | ---: |
| 1024 | 1.14 | 0.77 | 1.48x |
| 4096 | 6.19 | 4.61 | 1.34x |
| 16384 | 22.49 | 17.69 | 1.27x |
| 65536 | 77.80 | 71.44 | 1.09x |
| 262144 | 345.89 | 278.95 | 1.24x |

The gain is largest at small and mid shard sizes and decays into the
bandwidth-bound region, which is the opposite of the GFNI change's profile, so
the two compose rather than overlap.

Why it is deferred: unlike the affine change, this cannot reuse the nibble-table
storage shape or leave call sites alone.  Every kernel's loads, stores,
accumulators, tail thresholds, and the 64-byte-tile assumptions in the split and
range schedules have to move together, and the AVX-512VL member's measured
schedules (`Leopard2BackendAVX2.cpp` guards on `LEO2_AVX512_VARIANT`) were tuned
for 256-bit tiles.  It also needs its own frequency-behaviour check on Intel
parts, where 512-bit operation can down-clock.

### D2. Promote GFNI to its own qualified backend member

Tracked in `docs/leopard2_gfni_codec.md` under "Production integration
requirements".  GFNI itself is not AVX-512 — the VEX form needs only AVX plus
GFNI — so this belongs with the AVX2-family work, but it shares the "new
qualified member with its own allowed-instruction list" requirement with D1 and
is listed here so both are visible from one place.

## Also fixed in this session

`leopard2_visual_studio_project_self_test` was already failing at `HEAD`.  Commit
`ed3fd8a` added the `LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN` CMake option and its
`target_compile_definitions` mutation on the production `leopard` target without
registering either in the CMake-graph contract, so the guard fired on every one
of its 391 mutation subtests.  Registering the option, the approved production
mutation, its guard formula, and its position in the required command order
restores the suite to green (134 tests).

## Issue-tracker note

These findings were not filed as Beads issues.  `.beads/` is the 1.x
embedded-Dolt format (`bd-v1`, `embeddeddolt` present) and reports
"upgraded from v1.1.0 to v0.47.0 since last use", but the only `bd` on `PATH` is
`0.47.0` and it prints `LEGACY DATABASE DETECTED` on every invocation.
`CLAUDE.md` forbids running a legacy 0.x binary against this repository, so no
mutating `bd` command was run.  Install the 1.x binary and file these before
starting on them.
