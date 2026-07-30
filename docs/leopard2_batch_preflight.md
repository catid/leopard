# Scalable batch alias preflight

Leopard2 API version 6 exposes an opt-in caller-owned alias-preflight workspace
without changing the existing batch ABI:

- `leo2_encode_batch_preflight_scratch_size` and
  `leo2_encode_batch_with_preflight_scratch`;
- `leo2_decode_plan_batch_preflight_scratch_size` and
  `leo2_decode_plan_execute_batch_with_preflight_scratch`.

Version 6 also adds an encode-only reusable binding for workloads whose shard,
parity, and scratch addresses remain stable:

- `leo2_encode_batch_binding_create`;
- `leo2_encode_batch_binding_execute`;
- `leo2_encode_batch_binding_item_count`;
- `leo2_encode_batch_binding_destroy`.

The original `leo2_encode_batch` and `leo2_decode_plan_execute_batch` remain
correct, allocation-free compatibility entry points. They deliberately retain
their pairwise `O(B^2 * (K+R))` cross-item sweep. The new entry points use the
same item descriptors and per-item scratch. A query returns the additional
batch-wide workspace needed for a maximum parity subset; callers can allocate
or reserve it during setup and reuse it for later calls with the same or a
smaller batch.

## Algorithm and contract

The fast preflight first performs a read-only pass over every item. It checks
item-count arithmetic, layouts, physical byte lengths, pointer-array spans,
shard spans, GF16 padding, presence maps, per-item scratch size/alignment, and
the complete supplied batch-workspace span. It does not modify the workspace,
any item scratch, metadata, or output during this pass. In particular, a
workspace overlap with the item array, a pointer array, a shard, or per-item
scratch is rejected before the workspace is written.

After that pass, the implementation flattens all intervals into caller storage.
Each record holds a half-open address range and one writable bit. Item and
pointer arrays plus received shards are immutable; output shards and every
supplied scratch span are writable. An allocation-free partition moves writable
records to the front, and `std::sort` orders only that smaller set. A linear
sweep rejects writable/writable overlap. Each readable record then performs one
predecessor probe in the disjoint writable set; readable/readable overlap
remains permitted. The comparator is a strict total order on
`(begin, end, writable)`.

This evolved from a hand-rolled in-place heapsort to full introsort and then to
the writable-partitioned form. Sorting preflight metadata can dominate
whole-batch time for small shards, so every change is measured at the complete
public call rather than by comparison count alone. Historical full-sort
measurements on an AMD Ryzen 9 9950X3D, one pinned logical CPU, best of three
runs of `bench_leopard2_batch_preflight`, were:

| Cell | Encode speedup | Decode speedup |
| --- | ---: | ---: |
| K=3 R=2, B=64 | 1.76x | 2.58x |
| K=3 R=2, B=1024 | 1.73x | 1.35x |
| K=100 R=28, B=64 | 1.36x | 1.30x |
| K=100 R=28, B=1024 | 1.72x | 1.52x |

The writable-partitioned implementation made that original nine-item cutoff
stale. API version 5 began selecting scalable preflight from two items for every
codec except legacy-high `K=1,R=1`, whose specialized compatibility validator
remains faster through eight items. `B=1` is unchanged. A three-pass LSD radix
on `range.begin` measured faster still, but it needs a second equal-size buffer
and therefore a larger size query; that remains open work. The algorithm
deliberately permits input/input, metadata/metadata, and input/metadata sharing.

For `M = O(B * (K+R))` intervals, setup is `O(M log M)` time and `O(M)`
caller-owned bytes. Encoding, decoding, sorting, and the interval sweep perform
no allocation. A context pool may still start or grow worker threads lazily as
documented by `leo2_context_options`; this change does not hide that separate
scheduler setup inside the preflight workspace.

The size query returns zero for one item. Legacy-high `K=1,R=1` also returns
zero through eight items; all other codecs return a nonzero aligned span from
two items onward. No-loss decode always returns zero and remains a true no-op:
it does not inspect per-item pointers, byte sizes, aliases, or workspace.

The general `bench_leopard2` harness follows this recommended caller pattern:
it queries and allocates both batch-wide workspaces before timing, uses the
scalable entry points whenever the query is nonzero, and otherwise retains the
compatibility calls. Its reported batch scratch totals include the reusable
batch-wide workspace as well as every per-stripe codec scratch span. Benchmark
batch counts therefore measure the best safe public batch API rather than
accidentally timing the deliberately quadratic compatibility fallback at 64 or
1,024 items.

## Reusable encode binding

The scalable scratch API still validates and sorts `O(B * (K+R))` address
ranges on every call. That is the correct safe default when application
descriptors or buffers may change, but it can dominate 64-byte GF8 stripes.
The binding API moves that work into explicit setup. Creation deep-copies the
batch descriptors and all original/recovery pointer-array entries, validates
the full alias contract once, and captures each per-item scratch address.
Execution invokes only the already-validated codec items and allocates no
memory.

The caller may release or modify its descriptor and pointer arrays after
creation. The codec, context, captured shard buffers, and captured scratch
buffers must outlive the binding and remain at the same addresses and sizes.
The source bytes at captured addresses may change between calls. Because a
binding captures writable parity and scratch ranges, the same binding is not a
concurrent-execution object; callers use separate disjoint bindings for
concurrent batches. Context-worker startup remains lazy and should be paid
before a latency-sensitive execution if the context uses a thread pool.

A frozen-binary, one-core diagnostic on 2026-07-30 compared the exact
`main`-branch Leopard codec with the reusable Leopard2 binding. Each batch item
had independent output and scratch storage; shard size was 64 bytes. Times are
median microseconds, so smaller is better:

| K,R | Batch | Leopard main | Leopard2 binding | main / binding |
|---:|---:|---:|---:|---:|
| 5,5 | 1 | 0.049145 | 0.042465 | 1.157x |
| 5,5 | 8 | 0.413207 | 0.336215 | 1.229x |
| 5,5 | 64 | 2.939488 | 2.560965 | 1.148x |
| 5,5 | 1024 | 56.314828 | 45.698453 | 1.232x |
| 8,8 | 1 | 0.046328 | 0.067578 | 0.686x |
| 8,8 | 8 | 0.377621 | 0.533207 | 0.708x |
| 8,8 | 64 | 3.013469 | 4.310516 | 0.699x |
| 8,8 | 1024 | 60.508781 | 78.259438 | 0.773x |

This establishes the setup/execution separation and a production-useful K=5
win; it does not close the K=8 AVX2 kernel gap. An AUTO-disabled transform
comparison was within roughly 2–9% of the K=8 direct path across these batch
counts, so the remaining deficit is not explained by one dispatch threshold.
It remains tracked as codec-kernel work rather than being hidden by API timing.

Unsupported overlap and any ordinary later-item validation error reject the
whole batch before an earlier item executes. The preflight workspace contents
are unspecified on return, but item scratch, shard outputs, and immutable
metadata are untouched on rejection. The existing `UINT32_MAX` item limit and
checked `size_t`/address arithmetic are unchanged.

## Correctness evidence

`leopard2_batch_aliasing` covers both scalar and four-thread contexts and now
adds:

- exact workspace queries at the 1/2 cutoff, the legacy `K=1,R=1` exception,
  and overflow/error behavior;
- valid encode and decode at batches 64 and 1024 with alternating 17- and
  33-byte shard tails;
- shared immutable inputs and pointer arrays;
- partial and duplicate output overlaps, output/input overlaps, metadata
  overlaps, scratch/scratch overlaps, and workspace overlap with metadata,
  shards, and per-item scratch;
- the non-adjacent containing-interval case described above;
- 512 randomized encode/decode result-and-byte comparisons against the
  quadratic compatibility oracle across disjoint and colliding layouts;
- read-only rejection of a later invalid item before item scratch, output, or
  preflight scratch changes;
- allocation-audited repeated execution after setup;
- concurrent calls sharing immutable codecs and decode plans while using
  independent item and batch scratch;
- binding deep-copy independence, live source-byte reuse, setup-time conflict
  rejection, scalar/four-thread execution, and allocation-free repeated
  execution;
- a 1024-item no-loss plan with deliberately invalid per-item state, proving
  true no-op behavior.

The focused Release command is:

    cmake -S . -B build/release -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$(nproc)"
    ctest --test-dir build/release -R leopard2_batch_aliasing \
      --output-on-failure

ASan/UBSan and TSan builds should run the same test independently; the C++
allocation counter is disabled under TSan because that runtime supplies its own
replacement allocation functions.

## Diagnostic crossover evidence

The original crossover benchmark used scalar GF8, one-byte shards, one thread,
eight warmups, and the median of seven alternating-order measurements. It
compared complete valid batch calls, so each number included both alias
preflight and the same codec execution. On 2026-07-18 this 30-CPU host was
pinned to allowed CPU 0. Other project workers were active, so these are
historical diagnostic data, not release-authoritative latency claims:

| K,R | Batch | Encode compatibility / scalable ns | Encode gain | Decode compatibility / scalable ns | Decode gain | Encode / decode workspace bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 3,2 | 1 | 78.123 / 78.060 | 1.001x | 43.467 / 43.515 | 0.999x | 0 / 0 |
| 3,2 | 8 | 1204.062 / 1196.264 | 1.007x | 724.605 / 724.823 | 1.000x | 0 / 0 |
| 3,2 | 64 | 37253.100 / 13315.615 | 2.798x | 33214.085 / 9805.110 | 3.387x | 12352 / 12352 |
| 3,2 | 1024 | 7405664.500 / 401051.500 | 18.466x | 7786925.500 / 274079.500 | 28.411x | 196672 / 196672 |
| 100,28 | 1 | 5482.819 / 5489.383 | 0.999x | 7513.231 / 7521.176 | 0.999x | 0 / 0 |
| 100,28 | 8 | 56560.743 / 56582.913 | 1.000x | 63052.817 / 63043.820 | 1.000x | 0 / 0 |
| 100,28 | 64 | 801693.400 / 664546.750 | 1.206x | 948273.300 / 724241.150 | 1.309x | 201280 / 161344 |
| 100,28 | 1024 | 141596200.000 / 14166648.000 | 9.995x | 138352107.000 / 14577945.000 | 9.491x | 3219520 / 2580544 |

At that time the 1- and 8-item cells exercised the exact compatibility
implementation; their sub-one-percent differences measured
run noise/call-wrapper cost. The large-batch cells established that caller
scratch removes the quadratic setup term.

After writable partitioning, a pinned complete-call comparison of the old
nine-item dispatch and the new general two-item dispatch produced these encode
speedups at 64-byte shards:

| K,R | Batch 2 | Batch 4 | Batch 8 |
|---:|---:|---:|---:|
| 3,2 | 1.186x | 1.294x | 1.651x |
| 5,5 | 1.256x | 1.445x | 1.865x |
| 8,8 | 1.111x | 1.324x | 1.920x |
| 16,8 | 1.309x | 1.503x | 1.944x |
| 100,28 | 1.287x | 1.331x | 1.704x |
| 240,16 | 1.283x | 1.321x | 1.629x |

The 4-KiB neighbors also improved; the weakest tested general neighbor was
1.004x at K=100,R=28,batch=4. Legacy-high K=1,R=1 was the sole measured
exception, with new/old ratios of 0.725x, 0.783x, and 0.974x at batches 2, 4,
and 8 for 64-byte shards. It therefore retains the old nine-item dispatch
without penalizing other K=1 profiles.

Reproduce the dedicated preflight raw JSONL with:

    taskset -c 0 env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
      ./build/release/bench_leopard2_batch_preflight
