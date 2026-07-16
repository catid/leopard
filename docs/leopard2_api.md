# Leopard2 C API

`leopard2.h` adds an object-based API without changing `leopard.h` or the old
`leo_*` behavior.  The old entry points retain their count and 64-byte-size rules.

## Object lifetimes

Create one `leo2_context`, then any number of immutable codecs.  A codec fixes
`K`, `R`, profile, field, parent length, and coordinate map.  A decode plan copies
one presence/erasure pattern and precomputes its locator values, profile-specific
normalization factors, and deterministic received-coordinate selection.  The
context and codec must outlive their codecs and plans respectively.

Codec and plan execution is read-only and may be called concurrently when every
call has distinct outputs and scratch.  Setup may allocate; encode and reusable
plan execution do not.

## Encoding

Call `leo2_encode_scratch_size`, allocate at least that many bytes aligned to
`leo2_scratch_alignment()`, then call `leo2_encode`.  Every original pointer is
required.  A null recovery output requests that parity shard be skipped.  Inputs
may alias other inputs, but outputs and scratch must be disjoint from all shard
ranges; unsupported overlap returns `LEO2_OVERLAP`.

The high profile calls the existing block-IFFT accumulation and truncated final
FFT and matches old Leopard in the compatibility test.  The low profile performs
one padded-`P` IFFT followed by only the shifted parity blocks needed by the
requested output prefix/subset.

## Decoding

Create a plan from byte presence arrays with `leo2_decode_plan_create`.  Execution
must provide pointers that match that pattern.  Missing original inputs are null,
and their corresponding `restored_original` entries are writable.  Present
original output entries are ignored.  Missing parity is never rebuilt implicitly.

`leo2_decode_plan_execute` restores missing originals only.  A pattern with no
missing originals has zero scratch and is a true no-op, including when parity is
missing.  `R=1` high-profile repair and `K=1` low-profile repair use direct paths.
Other patterns use the specialized IT2026 low- or high-rate decoder selected by
the codec profile.  Locator construction, profile normalization, pruning inputs,
and output selection depend only on the erasure pattern and are performed during
plan creation; execution is the byte-heavy reusable step.

The plan uses exactly `K` received public coordinates.  It keeps every surviving
systematic shard, then keeps the lowest-index received parity shards needed to
reach `K`; any surplus received parity is treated as a deterministic virtual
erasure.  This selection affects only the work schedule, not the decoded message.
Applications still pass the original presence pattern and may leave pointers for
surplus received parity populated.

For differential testing and diagnosis, set
`LEO2_CODEC_FORCE_GENERIC_DECODE` in `leo2_codec_options.flags`.  Non-direct
transform repairs then use the retained full `O(N log N)` active-parent decoder
instead of the profile-specialized decoder.  The flag is not a distinct wire
profile and does not change encoded data.

`leo2_decode` is a convenience wrapper that allocates and destroys a plan.  Use a
reusable plan when setup amortization matters.

## Profiles

Legacy high V1:

    T = ceil_pow2(R)
    N = ceil_pow2(K + T)
    D = N - T
    [ parity 0..T-1 ][ message T..T+K-1 ][ shortened zero tail ]

Only parity `0..R-1` is transmitted; `R..T-1` is punctured.

Low V1:

    P = ceil_pow2(K)
    N = ceil_pow2(P + R)
    [ message 0..K-1 ][ shortened zeros to P ]
    [ parity P..P+R-1 ][ punctured parity tail ]

`AUTO` deterministically uses high when `R <= K`, otherwise low.  Field auto uses
GF8 for a parent of at most 256 coordinates and GF16 otherwise.  Profile and field
are mathematical code identity; backend selection never changes them.

## Byte lengths and batch calls

GF8 supports every positive byte length and internally handles SIMD tails through
caller scratch.  GF16 supports every positive even byte length.  Complete
64-byte ALTMAP tiles remain unchanged and use the zero-copy encoding path.  A
partial final GF16 tile contains `q` complete symbols in `2q` application bytes:
the first `q` bytes are scattered to ALTMAP low lanes `0..q-1`, the next `q`
bytes to high lanes `0..q-1`, and unused lanes are zero.  Outputs use the inverse
gather.  Odd GF16 byte lengths return `LEO2_UNSUPPORTED` because an unpaired byte
is not a complete GF16 symbol.  A no-loss decode plan remains a true no-op with
zero scratch even for an otherwise unsupported byte length.  See
`leopard2_gf16_tails.md` for the construction and the odd-length limitation.

The batch entry points execute independent items using each item's own buffers and
scratch.  A context owns a persistent worker pool when its effective
`thread_count` is greater than one, so batch calls do not create threads in the
hot path.  The calling thread participates in the batch; the pool contains
`thread_count - 1` workers.  A count of one executes batch items serially.

`thread_count = 0` selects `std::thread::hardware_concurrency()`, falls back to one
when it is unavailable, and caps the result at 128.  An explicit count from 1
through 128 is accepted; a larger count returns `LEO2_INVALID_ARGUMENT`.
`leo2_context_thread_count` reports the effective value.  This option controls
parallelism across batch items, not the wire profile or the result bytes.

Only one batch call at a time uses a given context pool.  Concurrent batch calls
sharing that context are serialized at the scheduler, while calls using separate
contexts may run independently.  Ordinary encode or plan-execute calls remain
safe to invoke concurrently with distinct outputs and scratch.  Within a batch,
the implementation reports the result from the lowest-index failing item
deterministically; completion order is otherwise unspecified.
