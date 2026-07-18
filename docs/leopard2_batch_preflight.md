# Scalable batch alias preflight

Leopard2 API version 4 adds an opt-in caller-owned alias-preflight workspace
without changing the existing batch ABI:

- `leo2_encode_batch_preflight_scratch_size` and
  `leo2_encode_batch_with_preflight_scratch`;
- `leo2_decode_plan_batch_preflight_scratch_size` and
  `leo2_decode_plan_execute_batch_with_preflight_scratch`.

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
supplied scratch span are writable. A custom in-place heapsort orders the
records without relying on an allocator. One max-end sweep then rejects an
overlap whenever either interval is writable. It deliberately permits
input/input, metadata/metadata, and input/metadata sharing. Tracking the
furthest prior end, rather than comparing only adjacent intervals, handles a
long immutable interval that contains another immutable interval before a
later output overlaps only the long interval.

For `M = O(B * (K+R))` intervals, setup is `O(M log M)` time and `O(M)`
caller-owned bytes. Encoding, decoding, sorting, and the interval sweep perform
no allocation. A context pool may still start or grow worker threads lazily as
documented by `leo2_context_options`; this change does not hide that separate
scheduler setup inside the preflight workspace.

The size query returns zero below nine items, and the new execution entry point
routes those calls through the existing implementation. This fixed cutoff
makes batches 1 and 8 identical to the compatibility path and avoids a small
batch regression. No-loss decode likewise returns zero and remains a true
no-op: it does not inspect per-item pointers, byte sizes, aliases, or workspace.

Unsupported overlap and any ordinary later-item validation error reject the
whole batch before an earlier item executes. The preflight workspace contents
are unspecified on return, but item scratch, shard outputs, and immutable
metadata are untouched on rejection. The existing `UINT32_MAX` item limit and
checked `size_t`/address arithmetic are unchanged.

## Correctness evidence

`leopard2_batch_aliasing` covers both scalar and four-thread contexts and now
adds:

- exact workspace queries at the 8/9 cutoff and overflow/error behavior;
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

The benchmark uses scalar GF8, one-byte shards, one thread, eight warmups, and
the median of seven alternating-order measurements. It compares complete valid
batch calls, so each number includes both alias preflight and the same codec
execution. When the size query returns zero, the harness invokes the ordinary
compatibility entry point, matching the intended caller dispatch. On
2026-07-18 this 30-CPU host was pinned to allowed CPU 0. Other
project workers were active, so these are diagnostic crossover data, not
release-authoritative latency claims:

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

The 1- and 8-item cells exercise the exact compatibility implementation; their
sub-one-percent differences measure run noise/call-wrapper cost and remain
inside the 2 percent gate. The large-batch cells show the intended removal of
quadratic setup work. Reproduce the raw JSONL with:

    taskset -c 0 env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
      ./build/release/bench_leopard2_batch_preflight
